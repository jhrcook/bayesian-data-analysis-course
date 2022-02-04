#!/usr/bin/env Rscript

library(glue)

# ---- Configuration ----

# Directory contains the Stan models.
MODELS_DIR <- file.path("models")
# Name of the notebook to generate.
OUTPUT_RMD_FILE <- file.path("stan_models.Rmd")

# ----

stan_model_files <- list.files(
  MODELS_DIR,
  pattern = "stan$",
  full.names = TRUE,
  recursive = FALSE,
  ignore.case = TRUE
)
stan_model_files <- sort(stan_model_files)
cat(glue("number Stan models: {length(stan_model_files)}"), "\n")

stan_model_to_markdown <- function(model_file) {
  model_txt <- readLines(model_file)
  model_txt <- paste(model_txt, collapse = "\n")
  md_txt <- "
  ## Model: `{basename(model_file)}`

  ```stan
  {model_txt}
  ```

  "
  return(glue(md_txt))
}

rmd_text <- "# (PART) Models {-}

# Stan models

Below are the Stan models built as a part of this course.
The original files are available in the GitHub repo in the \"models\" directory.

"

for (mdl_file in stan_model_files) {
  model_md <- stan_model_to_markdown(mdl_file)
  rmd_text <- paste(rmd_text, model_md, sep = "\n")
}

rmd_text <- paste(rmd_text, "\n")
writeLines(rmd_text, OUTPUT_RMD_FILE)
cat(glue("output: {OUTPUT_RMD_FILE}"), "\n")
