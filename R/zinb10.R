#' Example dataset 'zinb10'
#'
#' An example dataset generated from the proposed model with a zero-inflated negative binomial mediator (K=1). The mediator contains 10% zero values in which half are false zeros.
#'
#' @docType data
#' @usage data(zinb10)
#' @format An object of class \code{'data.frame'} with 100 rows and 3 variables:
#' \describe{
#'   \item{X}{independent variable, continuous data type}
#'   \item{Y}{outcome, continuous data type}
#'   \item{Mobs}{observed mediator values with possibly false zeros, count data type}
#' }
#' @keywords datasets
#' @examples
#' data(zinb10)
#' head(zinb10)
#'
"zinb10"
