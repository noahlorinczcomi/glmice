#' Prediction Accuracy and Diagnostic Statistics for Logistic Regression Models
#'
#' Returns accuracy of predictions, a chi-square test of model vs. null model, sensitivity, specificity, and % accuracy of predictions vs. null model for logistic regression models.
#' @param x a logistic regression model performed using 'glm(...family=binomial(link="logit"))'
#' @param plevel the level at which you would like to classify cases vs. non-cases (default is p=0.5)
#' @return accuracy of predictions, chi-square test of model vs. null, sensitivity, specificity, % advantage in classifications using model vs. null model.
#' @export

preds <- function(x, plevel = 0.5) {

  if(missing(plevel)) {

    pred_probs <- x$fitted.values
    predictions <- ifelse(pred_probs > 0.5, 1, 0)

    df <- data.frame(model.frame(x))
    y <- df[ , 1]

    table1 <- table(Predicted_Values = predictions, Actual_Values = y) # actual predictions
    acc <- sum(diag(table1)) / sum(table1) * 100 # accuracy of predictions

    null_model <- glm(y~1,
                      data=df, family=binomial(link="logit"))
    null_preds <- ifelse(fitted(null_model) > 0.5, 1, 0)
    null_bind <- data.frame(cbind(y, null_preds))
    null_acc1 <- nrow(null_bind[null_bind$y == null_bind$null_preds,])
    null_acc2 <- null_acc1 / nrow(df) * 100

    mod_acc_adv <- acc - null_acc2

    cs <- chisq.test(table1)

    sens <- sensitivity(table1)
    spec <- specificity(table1)

  } else if(!missing(plevel)){

    pred_probs <- x$fitted.values
    predictions <- ifelse(pred_probs > plevel, 1, 0)

    df <- data.frame(model.frame(x))
    y <- df[ , 1]

    table1 <- table(Predicted_Values = predictions, Actual_Values = y) # actual predictions
    acc <- sum(diag(table1)) / sum(table1) * 100 # accuracy of predictions

    null_model <- glm(y~1,
                      data=df, family=binomial(link="logit"))
    null_preds <- ifelse(fitted(null_model) > 0.5, 1, 0)
    null_bind <- data.frame(cbind(y, null_preds))
    null_acc1 <- nrow(null_bind[null_bind$y == null_bind$null_preds,])
    null_acc2 <- null_acc1 / nrow(df) * 100

    mod_acc_adv <- acc - null_acc2

    cs <- chisq.test(table1)

    sens <- sensitivity(table1)
    spec <- specificity(table1)


  }

  output <- list(table1,
                 paste("Accuracy of Predictions is",round(acc, digits = 2), "%"),
                 paste("Model is", round(mod_acc_adv, digits = 2), "% better than the median"),
                 cs,
                 paste("Sensitivity:", round(sens, digits=4)*100, "%;",
                       paste("Specificity:", round(spec, digits=4)*100, "%")))

  return(output)

}
