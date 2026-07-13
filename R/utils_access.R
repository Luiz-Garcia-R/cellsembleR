# ============================
# Internal helpers for rna_project access
# ============================

# -------------------------------------
# Generic accessor
# -------------------------------------
#' @keywords internal

.cs_get_from_container <- function(container,
                                   id = NULL,
                                   name = NULL) {

  if (is.null(container)) {
    stop("Container not found.")
  }

  # Resolve id
  if (is.null(id) || identical(id, "last")) {

    if (is.null(container$last)) {
      stop("No 'last' entry available in container.")
    }

    obj <- container[[container$last]]

  } else if (is.numeric(id)) {

    if (length(id) != 1) {
      stop("'id' must be a single numeric index.")
    }

    ids <- setdiff(names(container), "last")

    if (id < 1 || id > length(ids)) {
      stop("Index out of bounds.")
    }

    obj <- container[[ids[id]]]

  } else if (is.character(id)) {

    if (!id %in% names(container)) {
      stop("Invalid ID.")
    }

    obj <- container[[id]]

  } else {

    stop("Invalid 'id' argument.")

  }

  if (!is.null(name)) {

    if (!name %in% names(obj)) {
      stop(paste0("Field '", name, "' not found in object."))
    }

    return(obj[[name]])
  }

  obj
}

# -------------------------------------
# Expression matrix
# -------------------------------------
#' @keywords internal

.cs_get_expr <- function(project, id = NULL) {

  .cs_get_from_container(
    project$data$normalized_data,
    id,
    "expr_matrix"
  )
}

# -------------------------------------
# Metadata
# -------------------------------------
#' @keywords internal

.cs_get_meta <- function(project, id = NULL) {

  .cs_get_from_container(
    project$data$normalized_data,
    id,
    "metadata"
  )
}

# -------------------------------------
# Gene annotation
# -------------------------------------
#' @keywords internal

.cs_get_gene_annot <- function(project, id = NULL) {

  .cs_get_from_container(
    project$input$imp_data,
    id,
    "gene_annotation"
  )
}

# -------------------------------------
# Get gene annotation
# -------------------------------------
#' @keywords internal

.cs_align_gene_annot <- function(gene_annotation, expr_mat) {

  if (is.null(gene_annotation)) {
    stop("Gene annotation not found.")
  }

  genes_expr <- rownames(expr_mat)

  if (is.null(genes_expr)) {
    stop("Expression matrix must have rownames.")
  }

  genes_expr_clean <- sub("\\..*$", "", genes_expr)
  gene_annotation$gene_id <- sub("\\..*$", "", gene_annotation$gene_id)

  idx <- match(genes_expr_clean, gene_annotation$gene_id)

  aligned <- gene_annotation[idx, , drop = FALSE]
  aligned$gene_id <- genes_expr_clean

  return(aligned)
}
