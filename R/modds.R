#' Calculates an OR and associated probability for an 'average' case
#'
#' Returns the odds ratio and associated probability for a user-defined outcome (using logistic regression; 'glm()') for a case whose values are all at their median, estimating an 'average' case.
#' @param x a logit model performed using 'glm(...,family=binomial(link="logit"))'
#' @return odds: an odds ratio; prob: the associated probability analogous to this odds ratio
#' @export

modds <- function( x ) {

  # identifying if there is an interaction included in the model
  p <- t(x$coefficients)
  pnames <- colnames(p)
  p <- data.frame(p)
  colnames(p) <- pnames
  int_cols <- select(p, contains(":"))

  # (if no interaction)
  if(ncol(int_cols) == 0) {

    mf <- data.frame(model.frame( x ))
    mf <- data.frame(lapply(mf, function(x) as.numeric(as.character(x)))) # this turns all
    # factor variables to numeric variables

    g <- as.matrix(as.vector(apply(mf[,-1], 2, median)))

    x_values <- as.matrix(rbind( 1, g )) # added 1 for intercept

    coefs <- as.numeric(x$coefficients)

    bind <- as.matrix( cbind( coefs, x_values )) # can't bind if there are interaction terms
    log_odds <- crossprod( x_values,
                           coefs )

    odds <- exp( log_odds )
    prob <- odds / ( odds + 1 )

    # baseline marginal distribution of outcome
    oc <- mf[ , 1]
    bmd <- (length(oc[oc==1]) / length(oc)) * 100

  } else {

    mf <- data.frame(model.frame( x ))
    mf <- data.frame(lapply(mf, function(x) as.numeric(as.character(x)))) # this turns all
    # factor variables to numeric variables

    g <- as.matrix(as.vector(apply(mf[,-1], 2, median)))

    x_values <- as.matrix(rbind( 1, g )) # added 1 for intercept

    coefs <- as.numeric(x$coefficients)

    # working on getting names of interaction terms
    int_cols_name <- colnames(int_cols)

    int1 <- gsub(":.*","",int_cols_name)
    # gives just character string before the ":"

    # now you can use "s" to identify the other variable (i.e., the part of
    # int_cols_name that does not contain "int1")

    # here is the col_name without the ":"
    s2 <- gsub(":","",int_cols_name)

    int2 <- gsub(int1, "", s2)

    # now have terms 1 and 2 of the interaction, need to make an x_value for them

    x_int1 <- mf[,which(colnames(mf) == paste(int1))]
    x_int2 <- mf[,which(colnames(mf) == paste(int2))]

    x_int <- median((x_int1*x_int2))
    x_values2 <- rbind( x_values, as.numeric(x_int ))

    bind <- as.matrix( cbind( coefs, x_values2 )) # can't bind if there are interaction terms
    log_odds <- crossprod( x_values2,
                           coefs )

    odds <- exp( log_odds )
    prob <- odds / ( odds + 1 )

    # Baseline marginal distribution
    oc <- mf[ , 1]
    bmd <- (length(oc[oc==1]) / length(oc)) * 100

  }


  return(c( "Odds" = round(odds, digits=2),
            "Probability" = round(prob, digits=2),
            "% Baseline Marginal Distribution" = round(bmd, digits=2)))
}
