#!/usr/bin/env Rscript

# Run all notes and assignment R Makrdown notebooks.

NOTES_DIR <- "notes"
ASSIGNMENTS_DIR <- "assignments"

render_all_rmds <- function(dir) {
  cat(glue::glue("Rendering all R Markdowns in '{dir}':"), "\n")
  rmds <- list.files(dir, full.names = TRUE, pattern = "Rmd$", recursive = FALSE)
  rmds <- sort(rmds)
  for (rmd in rmds) {
    if (stringr::str_detect(rmd, "template")) {
      next
    }
    cat("  ", rmd, "\n")
    rmarkdown::render(rmd)
  }
}

render_all_rmds(NOTES_DIR)
render_all_rmds(ASSIGNMENTS_DIR)
