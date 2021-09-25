
# Check is two values (or vectors) are close to each other.
close_to <- function(x, y, epsilon = 0.01) {
  check_lengths <- length(x) == length(y)
  check_diff <- abs(x - y) < epsilon
  return(check_lengths & check_diff)
}

# Stop if the observed value is not close to the expected value.
# All inputs are passed to `close_to()`.
stop_if_not_close_to <- function(...) {
  stopifnot(close_to(...))
}

# Stop if not all observed values are close to expected values.
# All inputs are passed to `close_to()`.
stop_if_not_all_close_to <- function(...) {
  stopifnot(all(close_to(...)))
}
