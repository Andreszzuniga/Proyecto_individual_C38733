# scripts/04_calibration.R  (PARCHE: unir predictores al núcleo)
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
})

t0_days <- 3652.5
out_dir <- "output"
dir.create(out_dir, showWarnings = FALSE)

# Carga funciones de scoring y modelos
source("scripts/03_score.R")   # define: risk_10y, predict_risk_df
core <- readRDS("data/processed/nhanes_surv_core_9914.rds") %>%
  mutate(time_days = time_mo * (365.25/12))

# --- NUEVO: traer predictores completos y unir por seqn ---
pred <- readRDS("data/processed/nhanes_pred_all_full_9914.rds") %>%
  select(seqn, sbp, chol, hdl, bmi, smoker, diabetes, htn_treated)

# función para asegurar columnas (si faltan en 'pred', crearlas como NA)
ensure_cols <- function(df, cols_types = list(
  sbp=NA_real_, chol=NA_real_, hdl=NA_real_, bmi=NA_real_,
  smoker=NA_integer_, diabetes=NA_integer_, htn_treated=NA_integer_
)){
  for(nm in names(cols_types)){
    if(!nm %in% names(df)) df[[nm]] <- cols_types[[nm]]
  }
  df
}
pred <- ensure_cols(pred)

# unir al núcleo (no filtra filas; solo añade columnas)
core_pred <- core %>% left_join(pred, by = "seqn")

# construir el data.frame mínimo que espera predict_risk_df()
need <- c("sex","age","sbp","chol","hdl","bmi","smoker","diabetes","htn_treated")
have <- intersect(need, names(core_pred))
# si falta alguna, la añadimos como NA
missing <- setdiff(need, have)
if(length(missing)){
  for(m in missing) core_pred[[m]] <- if(m %in% c("smoker","diabetes","htn_treated")) NA_integer_ else NA_real_
}
core_pred <- core_pred[, need]

# calcular riesgo individual a 10 años
pred_risk <- predict_risk_df(core_pred)

dat <- readRDS("data/processed/nhanes_surv_core_9914.rds") %>%   # reabrimos para recuperar solo time/status/sex
  mutate(time_days = time_mo * (365.25/12),
         status = as.integer(status),
         sex_lab = ifelse(sex==1L, "Hombres", "Mujeres"),
         risk10 = pred_risk)
