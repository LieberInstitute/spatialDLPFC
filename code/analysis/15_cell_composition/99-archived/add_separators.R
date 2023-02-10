# https://stackoverflow.com/questions/69534248/how-can-i-make-a-discontinuous-axis-in-r-with-ggplot2

add_separators <- function(x, y = 0, angle = 45, length = .1) {
    add_y <- length * sin(angle * pi / 180)
    add_x <- length * cos(angle * pi / 180)
    ## making the list for your segments
    myseg <- list(
        x = x - add_x, xend = x + add_x,
        y = rep(y - add_y, length(x)), yend = rep(y + add_y, length(x))
    )
    ## this function returns an annotate layer with your segment coordinates
    annotate("segment",
        x = myseg$x, xend = myseg$xend,
        y = myseg$y, yend = myseg$yend
    )
}
