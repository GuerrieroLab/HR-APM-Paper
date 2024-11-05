## values are averaged by bin_values (binning by 20 pixels, then moving average by 10 running bins)
moving_average <- function(y, n = 3) {
  # Create a filter that averages 'n' points
  filter <- rep(1/n, n)
  # Apply the filter to the y values
  stats::filter(y, filter, sides = 2)
}

bin_values <- function(x, width = 20) {
  bin_start <- floor((x + width/2) / width) * width
  return(bin_start)
}

