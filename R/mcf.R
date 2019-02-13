#' Calculates McFadden's Pseudo R-Squared
#'
#' Returns McFadden's pseudo r-squared for logistic regression models performed on 'mice'-imputed data sets.
#' @param model a logit model from which you would like to return McFadden's pseudo r-squared. This can be a model created either with 'glm.mids()' or 'with()'
#' @return mcfs2: McFadden's pseudo r-squared
#' @export

mcf <- function (model) 
{
  iterations <- model$call1$m
  null_ds <- as.numeric()
  res_ds <- as.numeric()
  for (i in 1:iterations) {
    null_ds[i] <- model$analyses[[i]]$null.deviance
    res_ds[i] <- model$analyses[[i]]$deviance
  }
  ds <- cbind(as.numeric(null_ds), as.numeric(res_ds))
  m_null <- mean(null_ds)
  m_res <- mean(res_ds)
  mcfs <- round(((1 - (m_res / m_null)) * 100), 4)
  mcfs <- paste0(mcfs, "%")
  # end
  return(mcfs)
}
