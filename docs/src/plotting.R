

plot_dist <- function(x, y, xlab, ylab) {
  plot(
    x,
    y,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    ylim = c(0, max(y) * 1.02),
    yaxs = "i"
  )
}
