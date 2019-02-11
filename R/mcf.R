#' Calculates McFadden's Pseudo R-Squared
#'
#' Returns McFadden's pseudo r-squared for logistic regression models performed on 'mice'-imputed data sets.
#' @param model a logit model from which you would like to return McFadden's pseudo r-squared. This can be a model created either with 'glm.mids()' or 'with()'
#' @return mcfs2: McFadden's pseudo r-squared
#' @export

mcf <- function( model ) {

  iterations <- model$analyses[[1]]$iter

  for(i in 1:iterations) {
    null_ds <-  model$analyses[[i]]$null.deviance
    res_ds <-  model$analyses[[i]]$deviance
  }

  ds <- cbind( as.numeric(null_ds), as.numeric(res_ds))

  mcfs <- 1 - (res_ds/null_ds)
  ds <- cbind(mcfs, ds)

  mcfs2 <- (sum(mcfs)/nrow(ds))

  return(paste("McFadden's pseudo R-squared is",
               round(mcfs2, digits=2)))

}
