is.categorical <- function(x) {
  is.character(x) || is.factor(x)
}