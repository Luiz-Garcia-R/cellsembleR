# =====================================
# Check gene set coverage
# =====================================
#' @keywords internal
.check_gene_set_coverage <- function(expr_matrix,
                                     signatures,
                                     min_genes = 3) {

  coverage <- lapply(names(signatures), function(gs_name) {

    genes <- signatures[[gs_name]]

    detected_genes <- intersect(genes, rownames(expr_matrix))

    detected <- length(detected_genes)
    total <- length(genes)

    list(
      detected = detected,
      total = total,
      proportion = detected / total,
      missing_genes = setdiff(genes, rownames(expr_matrix)),
      detected_genes = detected_genes
    )
  })

  names(coverage) <- names(signatures)

  # -------------------------
  # Warning
  # -------------------------
  low_cov <- sapply(coverage, function(x) x$detected < min_genes)

  if (any(low_cov)) {

    warning(
      "Some gene sets have low coverage (< ",
      min_genes,
      " detected genes): ",
      paste(names(coverage)[low_cov], collapse = ", ")
    )
  }

  return(coverage)
}

# ============================
# Auxiliary print functions
# ============================
#' @keywords internal
.print_header <- function(title) {
  cat("\n")
  cat(strrep("=", 50), "\n")
  cat(title, "\n")
  cat(strrep("=", 50), "\n")
}

#' @keywords internal
.print_block <- function(title, content, width = 40) {
  cat("\n", title, "\n", sep = "")
  cat(strrep("-", width), "\n")
  content()
  cat(strrep("-", width), "\n")
}

# ============================
# Bootstrap
# ============================
#' @keywords internal
.boot_two_sample <- function(x, y, stat_fun, B = 2000, conf = 0.95) {

  x <- x[!is.na(x)]
  y <- y[!is.na(y)]

  boot_vals <- replicate(B, {

    xb <- sample(x, length(x), replace = TRUE)
    yb <- sample(y, length(y), replace = TRUE)

    stat_fun(xb, yb)

  })

  alpha <- (1 - conf) / 2

  ci <- stats::quantile(
    boot_vals,
    probs = c(alpha, 1 - alpha),
    na.rm = TRUE
  )

  list(
    boot = boot_vals,
    ci_low = ci[1],
    ci_high = ci[2]
  )
}

# ============================
# Calculate p value for corr_mat
# ============================
#' @keywords internal
.cor_pmat <- function(mat, method = "pearson") {

  n <- ncol(mat)

  p_mat <- matrix(
    NA,
    nrow = n,
    ncol = n
  )

  colnames(p_mat) <- colnames(mat)
  rownames(p_mat) <- colnames(mat)

  for (i in seq_len(n)) {

    for (j in seq_len(n)) {

      test <- suppressWarnings(
        stats::cor.test(
          mat[, i],
          mat[, j],
          method = method
        )
      )

      p_mat[i, j] <- test$p.value
    }
  }

  return(p_mat)
}

# =====================================
# Get features
# =====================================
#' @keywords internal
.cs_get_features <- function(x) {

  if (!is.null(x$signatures)) {
    return(names(x$signatures))
  }

  if (!is.null(x$signatures)) {
    return(x$signatures)
  }

  stop("No features/signatures found.")
}

# ============================
# Check dependencies
# ============================
#' @keywords internal
.cs_check_dependencies <- function(
    pkgs,
    source = c("cran", "bioc", "github")
) {

  source <- match.arg(source)

  missing <- pkgs[
    !vapply(
      pkgs,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]

  if (length(missing) == 0) {
    return(invisible(TRUE))
  }

  installer <- switch(

    source,

    cran = paste0(
      "install.packages(c(",
      paste(sprintf('"%s"', missing),
            collapse = ", "),
      "))"
    ),

    bioc = paste0(
      "BiocManager::install(c(",
      paste(sprintf('"%s"', missing),
            collapse = ", "),
      "))"
    ),

    github = paste0(
      'remotes::install_github("',
      missing,
      '")'
    )
  )

  stop(
    paste0(
      "Missing required packages: ",
      paste(missing, collapse = ", "),
      "\n\nInstall using:\n",
      installer
    ),
    call. = FALSE
  )
}

# ============================
# CS Prepare
# ============================
#' @keywords internal
.cs_prepare_expression_for_method <- function(project, id = NULL, method = c("mcp_counter", "epic", "xcell")) {

  method <- match.arg(method)

  expr_matrix <- .cs_get_expr(project, id)
  gene_annot  <- .cs_get_gene_annot(project, id)

  if (is.null(expr_matrix) || is.null(gene_annot)) {
    stop("Missing expression or annotation.")
  }

  expr_matrix <- as.matrix(expr_matrix)

  # -------------------------
  # align annotation
  # -------------------------
  gene_annot <- .cs_align_gene_annot(gene_annot, expr_matrix)

  # -------------------------
  # resolve gene symbols (non destructive)
  # -------------------------
  gene_symbol <- gene_annot$symbol
  gene_id     <- gene_annot$gene_id

  # fallback if symbol is unusable
  if (all(is.na(gene_symbol)) || method == "mcp_counter") {

    gene_symbol <- AnnotationDbi::mapIds(
      org.Hs.eg.db::org.Hs.eg.db,
      keys = gene_id,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  }

  # -------------------------
  # method-specific filtering rules
  # -------------------------
  valid <- !is.na(gene_symbol) & gene_symbol != ""

  if (method == "mcp_counter") {
    # MCP is strict
    valid <- valid & !grepl("^LOC", gene_symbol)
  }

  # -------------------------
  # build safe matrix (DO NOT mutate project)
  # -------------------------
  expr_matrix <- expr_matrix[valid, , drop = FALSE]
  gene_symbol <- gene_symbol[valid]

  rownames(expr_matrix) <- gene_symbol

  # collapse duplicates safely
  expr_matrix <- rowsum(expr_matrix, group = rownames(expr_matrix))

  return(expr_matrix)
}

# ============================
# CS attach
# ============================
#' @keywords internal
.cs_attach_if_needed <- function(pkg) {

  search_name <- paste0("package:", pkg)

  if (!search_name %in% search()) {

    suppressPackageStartupMessages(
      library(pkg,
              character.only = TRUE)
    )
  }
}
