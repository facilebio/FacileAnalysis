is.categorical <- function(x) {
  is.character(x) || is.factor(x)
}

is.integerish <- function(x) {
  test_integerish(x)
}
