.onAttach <- function(libname, pkgname) {

  version_text <- paste0(
    "cellsembleR v.",
    utils::packageVersion("cellsembleR")
  )

  if (requireNamespace("crayon", quietly = TRUE)) {
    version_text <- crayon::green(version_text)
  }

  packageStartupMessage(
    "\n",
    version_text, " loaded successfully!\n",
    "--------------------------------------------------\n",
    "A package for streamlined RNA-seq data analysis.\n",
    "GitHub: https://github.com/Luiz-Garcia-R/cellsembleR\n"
  )
}

