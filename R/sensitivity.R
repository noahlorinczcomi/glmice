#' Sensitivity Estimate for Logistic Regression Models.
#'
#' Returns a sensitivity estimate when comparing your predicted classifications to observed classifications.
#' @param predicted this parameter can either take (i) a 2x2 table/confusion matrix of predictions and observed classes (with predictions displayed on the vertical, left side), or (ii) your predicted classifications.
#' @param actual if param predicted is not a table, this parameter will take your actual (observed) values of your outcome.
#' @return sensitivity estimate
#' @export


sensitivity <- function(predicted, actual) {
  if(is.table(predicted) == TRUE) {
    a <- predicted[1]
    b <- predicted[3]
    c <- predicted[2]
    d <- predicted[4]

    sensitivity <- a/(a+c)
    return(sensitivity)
  } else {

    table1 <- table(predicted, actual)

    a <- table1[1]
    b <- table1[3]
    c <- table1[2]
    d <- table1[4]

    sensitivity <- a / (a+c)
    return(sensitivity)
  }
}

