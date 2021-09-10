
# Plot a distribution using base plotting methods.
plot_dist <- function(x, y, xlab, ylab, ...) {
  plot(
    x,
    y,
    type = "l",
    xlab = xlab,
    ylab = ylab,
    ylim = c(0, max(y) * 1.02),
    yaxs = "i",
    ...
  )
}


# Plot two distributions using ggplot.
plot_two_distributions <- function(x1, name0, x2, name1) {
  bind_rows(
    tibble(x = x1, name = name0),
    tibble(x = x2, name = name1)
  ) %>%
    ggplot(aes(x = x)) +
    geom_density(aes(color = name, fill = name), alpha = 0.2) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    theme(legend.position = c(0.85, 0.7))
}


# Plot a single distribution using ggplot.
plot_single_distribution <- function(x) {
  tibble(x = x) %>%
    ggplot(aes(x = x)) +
    geom_density() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02)))
}
