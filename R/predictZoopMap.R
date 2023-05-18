#' ZOOP MAP prediction
#' 
#' @description This function predicts ZOOP in the given day.
#'   The only covariates used are whale and regime (change if more are included).
#'   The spatial fields and probabilities of presence of whales must be 
#'   computed in advance and being introduced as arguments.
#' @param model model fitted with the function \code{\link{GibbsZoop}}
#' @param date new date "yyyy-mm-dd"
#' @param predictWhales probabilities of whale presence, rows are simulations,
#'   cols are days 181 x 17
#' @param gamma0s spatial field in temperature
#' @param etats spatial field in ZOOP. 
#'   Choose the eta from the year to predict.
#' @return matrix: rows are simulations columns are ZOOP
#' @export 
predictZoopMap <- function(model, date = "2011-04-29", predictWhales, gamma0s, etats) {
  
  B <- nrow(model)
  m <- ncol(gamma0s)
  
  day_year <- lubridate::yday(date)
  year <- lubridate::year(date)
  
  if (day_year < 33) {
    regime <- 0
  } else if (day_year < 90) {
    regime <- model[,"beta2"]
  } else {
    regime <- model[,"beta3"]
  }
  
  # temp <- matrix(nrow = B, ncol = m)
  temp <- model[,"sine"] * sin(2 * pi * day_year / 365) + 
    model[,"cosine"] * cos(2 * pi * day_year / 365) + 
    model[,paste0("psit", year-2002)] + gamma0s +
    rnorm(B * m, sd = 1 / sqrt(model[,"precTemp"]))
  
  # ZOOP <- matrix(nrow = B, ncol = m)
  p <- predictWhales[, 181 * (year - 2003) + day_year]
  ZOOP0 <- regime + model[,"betaTemp"] * temp + etats
  ZOOP1 <- ZOOP0 + model[,"beta1"]
  ZOOP  <- p * exp(ZOOP1) + (1 - p) * exp(ZOOP0)
  
  return(ZOOP)
}