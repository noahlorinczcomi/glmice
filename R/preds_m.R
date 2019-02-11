#' Prediction Accuracy and Diagnostic Statistics for Logistic Regression Models Performed on 'mice'-Imputed Data Sets
#'
#' Returns accuracy of predictions, a chi-square test of model vs. null model, sensitivity, specificity, and % accuracy of predictions vs. null model for logistic regression models performed on 'mice'-imputed data sets. This will use pooled coefficients from all models to calculate the accuracy of correct classifications using your model.
#' @param model a logistic regression model performed using either 'glm.mids()' or 'with()', from package 'mice'
#' @param imp the object to which your imputed data sets are assigned (i.e., the object that took the output of the mice() function)
#' @param plevel the level at which you would like to classify cases vs. non-cases (default is p=0.5)
#' @return accuracy of predictions, chi-square test of model vs. null, sensitivity, specificity, % advantage in classifications using model vs. null model.
#' @export


preds_m <- function(model, imp, plevel = 0.5) {

  # the function needs to identify if it is dealing with an interaction or not

  # getting all pooled coefficients
  pld <- pool(model)
  p_e_big <- pld$pooled$estimate[1:50]
  p_e_sml <- p_e_big[!is.na(p_e_big)]

  p <- t(pld$pooled[1])
  pnames <- colnames(p)
  p <- data.frame(p)
  colnames(p) <- pnames
  int_cols <- select(p, contains(":"))

  # this is if there IS a single interaction going on
  if(ncol(int_cols) != 0) {

    if(missing(plevel)) {


      # getting first complete data frame from imputations
      df <- data.frame(complete(imp, 1))

      # identifying variables used in interaction
      p <- t(pld$pooled[1])
      pnames <- colnames(p)
      p <- data.frame(p)
      colnames(p) <- pnames
      int_cols <- select(p, contains(":"))
      rownames(int_cols) <- NULL

      # working on getting names of interaction terms
      int_cols_name <- colnames(int_cols)
      int1 <- gsub(":.*","",int_cols_name)
      s2<-gsub(":","",int_cols_name)
      int2 <- gsub(int1, "", s2)

      # now have terms 1 and 2 of the interaction, need to make an x_value for them
      # for the later calculation of the odds/prob when all values are at median
      x_int1 <- df[,which(colnames(mf) == paste(int1))]
      x_int2 <- df[,which(colnames(mf) == paste(int2))]

      x_int <- median((x_int1*x_int2))

      # getting coeffcients and coefficient names to (1) make sure the df has the same
      # column ordering as the term ordering in the call, and (2) to identify variables
      # in the interaction to calculate odds/probs
      p_coefs <- data.frame(t(pld$pooled[1]))
      colnames(p)[1] <- colnames(df)[1]

      # names of coefficients (will include intercept and interactions)
      pnames2 <- colnames(p_coefs)

      # making new column to make room for the interaction, and renaming it to the name
      # of the interaction in the call
      df <- df %>%
        mutate(new = 99)

      colnames(df)[which(names(df) == "new")] <- paste(pnames[length(pnames)])

      # changing "intercept" name in the coefficients to the name of the outcome in the df
      # (for ordering the columns of the df)
      pnames[1] <- colnames(df)[1]

      df <- df %>%
        select(pnames)

      # making the actual interaction values
      df[,ncol(df)] <- x_int1 * x_int2

      # identifying outcome and coefficients (numeric only, no names)
      actual_outcome <- df[,1]
      coefs <- as.numeric(p)

      # making the outcome = 1 because when computing odds/probs the intercept takes
      # the value of 1 as its partner
      df[,1] <- 1

      # calculating log odds, odds, and probs
      for(i in 1:nrow(df)) {
        logs <- sum(as.numeric(df[i,]) * coefs)
      }

      for(i in 1:nrow(df)) {
        logs[i] <- sum(as.numeric(df[i,]) * coefs)
      }

      odds <- exp(logs)
      pred_probs <- odds / (1+odds)

      # predicted values
      predictions <- ifelse(pred_probs > 0.5, 1, 0)

      # table of predicted an actual values. Associated accuracy of predictions
      table1 <- table(Predicted_Values = predictions, Actual_Values = actual_outcome)
      acc <- sum(diag(table1)) / sum(table1) * 100

      # calculating a null model to test the improvement in accuract of the model
      null_model <- glm(actual_outcome~1,
                        data=df, family=binomial(link="logit"))
      null_preds <- ifelse(fitted(null_model) > 0.5, 1, 0)
      null_bind <- data.frame(cbind(actual_outcome, null_preds))
      null_acc1 <- nrow(null_bind[null_bind$actual_outcome == null_bind$null_preds,])
      null_acc2 <- null_acc1 / nrow(df) * 100

      # % advantage of model over null model in making predictions
      mod_acc_adv <- acc - null_acc2

      # odds and prob when variable are at their median
      df <- data.frame(lapply(df, function(x) as.numeric(as.character(x)))) # turns all factors to numeric for median calculation
      g <- as.matrix(as.vector(apply(df[,-1], 2, median)))
      x_values <- as.matrix(rbind( 1, g )) # added 1 for intercept

      modd <- exp(crossprod(coefs, x_values))
      prob <- modd / (1+modd)

      # chi-square test
      cs <- chisq.test(table1)

      # sensitivity and specificity
      sens <- sensitivity(table1)
      spec <- specificity(table1)


      # if a plevel is specified
    } else if(!missing(plevel)) {

      # getting all pooled coefficients
      pld <- pool(model)
      p_e_big <- pld$pooled$estimate[1:50]
      p_e_sml <- p_e_big[!is.na(p_e_big)]

      # getting first complete data frame from imputations
      df <- data.frame(complete(imp, 1))

      # identifying variables used in interaction
      p <- t(pld$pooled[1])
      pnames <- colnames(p)
      p <- data.frame(p)
      colnames(p) <- pnames
      int_cols <- select(p, contains(":"))
      rownames(int_cols) <- NULL

      # working on getting names of interaction terms
      int_cols_name <- colnames(int_cols)
      int1 <- gsub(":.*","",int_cols_name)
      s2<-gsub(":","",int_cols_name)
      int2 <- gsub(int1, "", s2)

      # now have terms 1 and 2 of the interaction, need to make an x_value for them
      # for the later calculation of the odds/prob when all values are at median
      x_int1 <- df[,which(colnames(mf) == paste(int1))]
      x_int2 <- df[,which(colnames(mf) == paste(int2))]

      x_int <- median((x_int1*x_int2))

      # getting coeffcients and coefficient names to (1) make sure the df has the same
      # column ordering as the term ordering in the call, and (2) to identify variables
      # in the interaction to calculate odds/probs
      p_coefs <- data.frame(t(pld$pooled[1]))
      colnames(p)[1] <- colnames(df)[1]

      # names of coefficients (will include intercept and interactions)
      pnames2 <- colnames(p_coefs)

      # making new column to make room for the interaction, and renaming it to the name
      # of the interaction in the call
      df <- df %>%
        mutate(new = 99)

      colnames(df)[which(names(df) == "new")] <- paste(pnames[length(pnames)])

      # changing "intercept" name in the coefficients to the name of the outcome in the df
      # (for ordering the columns of the df)
      pnames[1] <- colnames(df)[1]

      df <- df %>%
        select(pnames)

      # making the actual interaction values
      df[,ncol(df)] <- x_int1 * x_int2

      # identifying outcome and coefficients (numeric only, no names)
      actual_outcome <- df[,1]
      coefs <- as.numeric(p)

      # making the outcome = 1 because when computing odds/probs the intercept takes
      # the value of 1 as its partner
      df[,1] <- 1

      # calculating log odds, odds, and probs
      for(i in 1:nrow(df)) {
        logs <- sum(as.numeric(df[i,]) * coefs)
      }

      for(i in 1:nrow(df)) {
        logs[i] <- sum(as.numeric(df[i,]) * coefs)
      }

      odds <- exp(logs)
      pred_probs <- odds / (1+odds)

      # predicted values
      predictions <- ifelse(pred_probs > plevel, 1, 0)

      # table of predicted an actual values. Associated accuracy of predictions
      table1 <- table(Predicted_Values = predictions, Actual_Values = actual_outcome)
      acc <- sum(diag(table1)) / sum(table1) * 100

      # calculating a null model to test the improvement in accuract of the model
      null_model <- glm(actual_outcome~1,
                        data=df, family=binomial(link="logit"))
      null_preds <- ifelse(fitted(null_model) > plevel, 1, 0)
      null_bind <- data.frame(cbind(actual_outcome, null_preds))
      null_acc1 <- nrow(null_bind[null_bind$actual_outcome == null_bind$null_preds,])
      null_acc2 <- null_acc1 / nrow(df) * 100

      # % advantage of model over null model in making predictions
      mod_acc_adv <- acc - null_acc2

      # odds and prob when variable are at their median
      df <- data.frame(lapply(df, function(x) as.numeric(as.character(x)))) # turns all factors to numeric for median calculation
      g <- as.matrix(as.vector(apply(df[,-1], 2, median)))
      x_values <- as.matrix(rbind( 1, g )) # added 1 for intercept

      modd <- exp(crossprod(coefs, x_values))
      prob <- modd / (1+modd)

      # chi-square test
      cs <- chisq.test(table1)

      # sensitivity and specificity
      sens <- sensitivity(table1)
      spec <- specificity(table1)

    }


    # if there is no interaction going on
  } else if(ncol(int_cols) == 0) {

    # and there is no specified level
    if(missing(plevel)) {

      pld <- pool(model)
      p_e_big <- pld$pooled$estimate[1:50]
      p_e_sml <- p_e_big[!is.na(p_e_big)]


      df <- data.frame(complete(imp, 1))

      p <- data.frame(t(pld$pooled[1]))
      colnames(p)[1] <- colnames(df)[1]

      pnames <- colnames(p)

      df <- df %>%
        select(pnames)

      actual_outcome <- df[,1]
      coefs <- as.numeric(p)

      df[,1] <- 1

      for(i in 1:nrow(df)) {
        logs <- sum(as.numeric(df[i,]) * coefs)
      }

      for(i in 1:nrow(df)) {
        logs[i] <- sum(as.numeric(df[i,]) * coefs)
      }

      odds <- exp(logs)
      pred_probs <- odds / (1+odds)

      predictions <- ifelse(pred_probs > 0.5, 1, 0)

      table1 <- table(Predicted_Values = predictions, Actual_Values = actual_outcome)
      acc <- sum(diag(table1)) / sum(table1) * 100

      null_model <- glm(actual_outcome~1,
                        data=df, family=binomial(link="logit"))
      null_preds <- ifelse(fitted(null_model) > 0.5, 1, 0)
      null_bind <- data.frame(cbind(actual_outcome, null_preds))
      null_acc1 <- nrow(null_bind[null_bind$actual_outcome == null_bind$null_preds,])
      null_acc2 <- null_acc1 / nrow(df) * 100

      mod_acc_adv <- acc - null_acc2

      cs <- chisq.test(table1)

      # sensitivity and specificity
      sens <- sensitivity(table1)
      spec <- specificity(table1)

      # if a plevel is specified
    } else if(!missing(plevel)){

      pld <- pool(model)
      p_e_big <- pld$pooled$estimate[1:50]
      p_e_sml <- p_e_big[!is.na(p_e_big)]


      df <- data.frame(complete(imp, 1))

      p <- data.frame(t(pld$pooled[1]))
      colnames(p)[1] <- colnames(df)[1]

      pnames <- colnames(p)

      df <- df %>%
        select(pnames)

      actual_outcome <- df[,1]
      coefs <- as.numeric(p)

      df[,1] <- 1

      for(i in 1:nrow(df)) {
        logs <- sum(as.numeric(df[i,]) * coefs)
      }

      for(i in 1:nrow(df)) {
        logs[i] <- sum(as.numeric(df[i,]) * coefs)
      }

      odds <- exp(logs)
      pred_probs <- odds / (1+odds)

      predictions <- ifelse(pred_probs > plevel, 1, 0)

      table1 <- table(Predicted_Values = predictions, Actual_Values = actual_outcome)
      acc <- sum(diag(table1)) / sum(table1) * 100

      null_model <- glm(actual_outcome~1,
                        data=df, family=binomial(link="logit"))
      null_preds <- ifelse(fitted(null_model) > plevel, 1, 0)
      null_bind <- data.frame(cbind(actual_outcome, null_preds))
      null_acc1 <- nrow(null_bind[null_bind$actual_outcome == null_bind$null_preds,])
      null_acc2 <- null_acc1 / nrow(df) * 100

      mod_acc_adv <- acc - null_acc2

      cs <- chisq.test(table1)

      # sensitivity and specificity
      sens <- sensitivity(table1)
      spec <- specificity(table1)


    }


  }
  output <- list(table1,
                 paste("The accuracy of Predictions is",round(acc, digits = 2), "%"),
                 paste("The model is", round(mod_acc_adv, digits = 2), "% better than the null model"),
                 cs,
                 paste("Sensitivity:", round(sens,digits=4)*100, "%;"),
                 paste("Specificity:", round(spec,digits=4)*100, "%"))
  return(output)
}
