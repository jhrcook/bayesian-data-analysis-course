
# Get the path to a data file.
data_file_path <- function(file_name) {
  return(here::here("data", file_name))
}

# Read data from the "data/" directory.
read_data <- function(file_name, .f = readLines, ...) {
  return(.f(data_file_path(file_name), ...))
}

# This file has weird formatting so requires special processing.
read_bioassay_data <- function(...) {
  read_data("bioassay.txt", read_tsv, ...) %>%
    separate(`x n y`, into = c("x", "n", "y"), sep = " +") %>%
    mutate_all(as.numeric)
}
