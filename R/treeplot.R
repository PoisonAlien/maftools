#This code is a direct translation of python package 'squarify' by Uri Laserson (https://github.com/laserson/squarify)
#See https://github.com/laserson/squarify for more details

# Squarified Treemap Layout
# Implements algorithm from Bruls, Huizing, van Wijk, "Squarified Treemaps"
#   (but not using their pseudocode)

pad_rectangle <- function(rect) {
  if (rect$dx > 2) {
    rect$x <- rect$x + 1
    rect$dx <- rect$dx - 2
  }
  if (rect$dy > 2) {
    rect$y <- rect$y + 1
    rect$dy <- rect$dy - 2
  }
}

layoutrow <- function(sizes, x, y, dx, dy) {
  # generate rects for each size in sizes
  # dx >= dy
  # they will fill up height dy, and width will be determined by their area
  # sizes should be pre-normalized wrt dx * dy (i.e., they should be same units)
  covered_area <- sum(sizes)
  width <- covered_area / dy
  rects <- list()
  for (size in sizes) {
    rects <- c(rects, list(list(x = x, y = y, dx = width, dy = size / width)))
    y <- y + size / width
  }
  return(rects)
}

layoutcol <- function(sizes, x, y, dx, dy) {
  # generate rects for each size in sizes
  # dx < dy
  # they will fill up width dx, and height will be determined by their area
  # sizes should be pre-normalized wrt dx * dy (i.e., they should be same units)
  covered_area <- sum(sizes)
  height <- covered_area / dx
  rects <- list()
  for (size in sizes) {
    rects <- c(rects, list(list(x = x, y = y, dx = size / height, dy = height)))
    x <- x + size / height
  }
  return(rects)
}

layoutMain <- function(sizes, x, y, dx, dy) {
  if (dx >= dy) {
    return(layoutrow(sizes, x, y, dx, dy))
  } else {
    return(layoutcol(sizes, x, y, dx, dy))
  }
}

leftoverrow <- function(sizes, x, y, dx, dy) {
  # compute remaining area when dx >= dy
  covered_area <- sum(sizes)
  width <- covered_area / dy
  leftover_x <- x + width
  leftover_y <- y
  leftover_dx <- dx - width
  leftover_dy <- dy
  return(list(leftover_x, leftover_y, leftover_dx, leftover_dy))
}

leftovercol <- function(sizes, x, y, dx, dy) {
  # compute remaining area when dx >= dy
  covered_area <- sum(sizes)
  height <- covered_area / dx
  leftover_x <- x
  leftover_y <- y + height
  leftover_dx <- dx
  leftover_dy <- dy - height
  return(list(leftover_x, leftover_y, leftover_dx, leftover_dy))
}

leftover <- function(sizes, x, y, dx, dy) {
  if (dx >= dy) {
    return(leftoverrow(sizes, x, y, dx, dy))
  } else {
    return(leftovercol(sizes, x, y, dx, dy))
  }
}

worst_ratio <- function(sizes, x, y, dx, dy) {
  rects <- layoutMain(sizes, x, y, dx, dy)
  return(max(sapply(rects, function(rect) {
    max(rect$dx / rect$dy, rect$dy / rect$dx)
  })))
}

squarify <- function(sizes, x, y, dx, dy) {
  sizes <- as.numeric(sizes)

  if (length(sizes) == 0) {
    return(list())
  }

  if (length(sizes) == 1) {
    return(layoutMain(sizes, x, y, dx, dy))
  }

  # figure out where 'split' should be
  i <- 1
  while (i < length(sizes) && worst_ratio(sizes[1:i], x, y, dx, dy) >= worst_ratio(sizes[1:(i + 1)], x, y, dx, dy)) {
    i <- i + 1
  }
  current <- sizes[1:i]
  remaining <- sizes[(i + 1):length(sizes)]

  leftover_info <- leftover(current, x, y, dx, dy)
  leftover_x <- leftover_info[[1]]
  leftover_y <- leftover_info[[2]]
  leftover_dx <- leftover_info[[3]]
  leftover_dy <- leftover_info[[4]]

  return(c(layoutMain(current, x, y, dx, dy), squarify(remaining, leftover_x, leftover_y, leftover_dx, leftover_dy)))
}

padded_squarify <- function(sizes, x, y, dx, dy) {
  rects <- squarify(sizes, x, y, dx, dy)
  for (rect in rects) {
    pad_rectangle(rect)
  }
  return(rects)
}

normalize_sizes <- function(sizes, dx, dy) {
  total_size <- sum(sizes)
  total_area <- dx * dy
  sizes <- sizes * (total_area / total_size)
  return(sizes)
}

