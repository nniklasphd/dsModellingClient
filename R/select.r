#' Required by the 'datashield.glm' function
#' @title Selects specific elements in a list
#' @param l the list to select elements from 
#' @params field the location of the elements to select
#' @return the selected elements
.select <- function(l, field) {
  lapply(l, function(obj) {obj[[field]]})
}
