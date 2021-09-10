
# Check is two values (or vectors) are close to each other.
close_to <- function(x, y, epsilon = 0.01) {
  return(abs(x - y) < epsilon)
}
