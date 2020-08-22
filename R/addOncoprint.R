add_oncoprint <- function(type, x, y, width, height, variant.classes, type_col) {
  for (i in 1:length(variant.classes)) {
    if (any(type %in% variant.classes[i])) {
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
        grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = type_col[variant.classes[i]]))
    } else if (any(type %in% "Amp")) {
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
        grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = bg))
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
        unit(8, "mm"), gp = grid::gpar(col = NA, fill = type_col["Amp"]))
    } else if (any(type %in% "Del")) {
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
        grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = bg))
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height - grid::unit(8, "mm"),
        gp = grid::gpar(col = NA, fill = type_col["Del"])
      )
    }
  }

  if (any(type %in% "")) {
    grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
      grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = bg))
  }
}


add_oncoprint2 <- function(type, x, y, width, height, variant.classes, type_col) {
  for (i in 1:length(variant.classes)) {
    if (any(type %in% variant.classes[i])) {
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
        grid::unit(1, "mm"), gp = grid::gpar(col = NA, fill = type_col[variant.classes[i]]))
    } else if (any(type %in% "Amp")) {
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height -
        unit(8, "mm"), gp = grid::gpar(col = NA, fill = type_col["Amp"]))
    } else if (any(type %in% "Del")) {
      grid::grid.rect(x, y, width - unit(0.5, "mm"), height - grid::unit(8, "mm"),
        gp = grid::gpar(col = NA, fill = type_col["Del"])
      )
    }
  }
}
