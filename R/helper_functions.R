# Helper function to format LOLE measures
lole_bind_rows <- function(x) {
  name_ <- substitute(x)
  class(x) <- "list"
  o <- bind_rows(x)
  colnames(o) <- paste(name_, colnames(o), sep = ".")
  return(o)
}