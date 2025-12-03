suppressPackageStartupMessages({ library(dplyr) })

# Usa las funciones de scoring (risk_10y, predict_risk_df)
source("scripts/03_score.R")

core <- readRDS("data/processed/nhanes_surv_core_9914.rds")
pred <- readRDS("data/processed/nhanes_pred_all_full_9914.rds") %>%
  dplyr::select(seqn, sbp, chol, hdl, bmi, smoker, diabetes, htn_treated)

# Asegurar columnas opcionales si faltan
need <- c("sbp","chol","hdl","bmi","smoker","diabetes","htn_treated")
for (nm in need) if (!nm %in% names(pred)) pred[[nm]] <- NA

core_pred <- core %>% dplyr::left_join(pred, by = "seqn")

scores <- core_pred %>%
  dplyr::transmute(
    seqn, sex, age, sbp, chol, hdl, bmi, smoker, diabetes, htn_treated,
    risk10 = predict_risk_df(
      data.frame(
        sex = sex, age = age, sbp = sbp, chol = chol, hdl = hdl, bmi = bmi,
        smoker = smoker, diabetes = diabetes, htn_treated = htn_treated
      )
    )
  )

dir.create("output", showWarnings = FALSE)
write.csv(scores, "output/risk10_scores.csv", row.names = FALSE)
cat("OK -> output/risk10_scores.csv\n")
