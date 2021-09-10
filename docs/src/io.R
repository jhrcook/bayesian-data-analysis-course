
# Get the path to a data file.
data_file_path <- function(file_name) {
  return(here::here("data", file_name))
}

# Read data from the "data/" directory.
read_data <- function(file_name, .f = readLines) {
  return(.f(data_file_path(file_name)))
}
