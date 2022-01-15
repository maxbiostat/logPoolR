coord_radar <- function (theta = "x", start = 0, direction = 1) 
{ ## from http://web-r.org/board_ISVF22/8271
  theta <- match.arg(theta, c("x", "y"))
  r <- if (theta == "x") 
    "y"
  else "x"
  ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
          direction = sign(direction),
          is_linear = function(coord) TRUE)
}
number_ticks <- function(n) {function(limits) pretty(limits, n)}