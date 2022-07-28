as.aaplot <- function(object) {
  if (inherits(object, "aaplot"))
    return(object)

  if (!inherits(object, 'aplot')) {
    stop("input should be a 'aplot' object.")
  }

  attr(object,"class") <- "aaplot"

  return(object)
}


##' @method print aaplot
##' @importFrom patchwork plot_layout
##' @importFrom patchwork plot_spacer
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 theme_void
##' @export
print.aaplot <- function(x, ...) {
  grid.draw(x)
}


##' @importFrom grid grid.draw
##' @importFrom aplot as.patchwork
##' @method grid.draw aaplot
##' @export
grid.draw.aaplot <- function(x, recoding = TRUE) {
  attr(x,"class") <- "aplot"
  grid::grid.draw(as.patchwork(x,align="none"))
}
