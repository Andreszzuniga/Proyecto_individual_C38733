# scripts/03_score.R
# Scoring a 10 años por sexo; usa modelo combinado si falta el específico

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

mods   <- readRDS("data/processed/cox_models_ALL.rds")
fit_m  <- mods$fit_m
fit_f  <- mods$fit_f
fit_all<- mods$fit_all

get_S0_at <- function(fit, t0 = 3652.5){
  bh <- basehaz(fit, centered = TRUE)
  idx <- max(which(bh$time <= t0)); if(is.infinite(idx)) idx <- nrow(bh)
  exp(-bh$hazard[idx])
}

S0_m_10y   <- if(!is.null(fit_m))  get_S0_at(fit_m)  else NA_real_
S0_f_10y   <- if(!is.null(fit_f))  get_S0_at(fit_f)  else NA_real_
S0_all_10y <- if(!is.null(fit_all))get_S0_at(fit_all)else NA_real_

risk_side <- function(fit, S0_10y, newx){
  coefs <- coef(fit); vars <- names(coefs)[!is.na(coefs)]
  means <- fit$means; means <- means[intersect(names(means), vars)]
  x <- newx[vars]
  for(v in intersect(names(means), names(x))) x[v] <- x[v] - means[[v]]
  L <- sum(coefs[vars] * x)
  1 - (S0_10y)^(exp(L))
}

# sex: 1=hombre, 0=mujer
risk_10y <- function(sex, age, sbp, chol, hdl, bmi, smoker = NA, diabetes = NA, htn_treated = NA){
  sex <- as.integer(sex)
  if(sex == 1L && !is.null(fit_m)){
    x <- c(age=age, sbp=sbp, chol=chol, hdl=hdl, bmi=bmi, smoker=smoker, diabetes=diabetes, htn_treated=htn_treated)
    return(risk_side(fit_m, S0_m_10y, x))
  }
  if(sex == 0L && !is.null(fit_f)){
    x <- c(age=age, sbp=sbp, chol=chol, hdl=hdl, bmi=bmi, smoker=smoker, diabetes=diabetes, htn_treated=htn_treated)
    return(risk_side(fit_f, S0_f_10y, x))
  }
  # respaldo: modelo combinado (incluye 'sex' en la LP)
  if(!is.null(fit_all)){
    x <- c(sex=sex, age=age, sbp=sbp, chol=chol, hdl=hdl, bmi=bmi, smoker=smoker, diabetes=diabetes, htn_treated=htn_treated)
    return(risk_side(fit_all, S0_all_10y, x))
  }
  stop("No hay modelo disponible para calcular el riesgo.")
}

predict_risk_df <- function(df){
  stopifnot(all(c("sex","age","sbp","chol","hdl","bmi") %in% names(df)))
  mapply(risk_10y,
         sex  = df$sex,
         age  = df$age,
         sbp  = df$sbp,
         chol = df$chol,
         hdl  = df$hdl,
         bmi  = df$bmi,
         smoker = if("smoker" %in% names(df)) df$smoker else NA,
         diabetes = if("diabetes" %in% names(df)) df$diabetes else NA,
         htn_treated = if("htn_treated" %in% names(df)) df$htn_treated else NA)
}
