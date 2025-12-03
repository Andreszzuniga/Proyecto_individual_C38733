# scripts/02_fit_cox_ALL.R
# Ajusta modelos de Cox por sexo y un modelo combinado de respaldo

suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

core <- readRDS("data/processed/nhanes_surv_core_9914.rds") %>%
  mutate(time_days = time_mo * (365.25/12))
pred <- readRDS("data/processed/nhanes_pred_all_full_9914.rds")

print(with(core, table(sex, status)))

get_S0_at <- function(fit, t0 = 3652.5){
  bh <- basehaz(fit, centered = TRUE)
  idx <- max(which(bh$time <= t0)); if(is.infinite(idx)) idx <- nrow(bh)
  exp(-bh$hazard[idx])
}
cindex_safe <- function(fit, df){
  if (is.null(fit)) return(NA_real_)
  used_idx <- if (!is.null(fit$na.action)) setdiff(seq_len(nrow(df)), fit$na.action) else seq_len(nrow(df))
  df_used <- df[used_idx, , drop = FALSE]
  lp <- predict(fit, type = "lp")
  survival::concordance(Surv(time_days, status) ~ lp, data = df_used)$concordance
}

optional <- c("sbp","chol","hdl","bmi","smoker","diabetes","htn_treated")
pick_vars <- function(df, min_prop = 0.40){
  avail <- intersect(optional, names(df))
  keep <- vapply(avail, function(v){
    nn <- sum(!is.na(df[[v]])); prop <- nn/nrow(df)
    uniq <- length(unique(na.omit(df[[v]]))) > 1
    prop >= min_prop && uniq
  }, logical(1))
  avail[keep]
}

fit_by_sex <- function(sex_flag, label){
  df <- core %>% filter(sex == sex_flag) %>%
    select(seqn, sex, age, time_days, status) %>%
    left_join(pred %>% select(seqn, all_of(optional)), by = "seqn")
  if(nrow(df) == 0 || sum(df$status, na.rm = TRUE) == 0){
    message("⚠ Sin datos o sin eventos para ", label)
    return(NULL)
  }
  vars <- c("age", pick_vars(df))
  message("Variables usadas en ", label, ": ", paste(vars, collapse = ", "))
  form <- as.formula(paste("Surv(time_days, status) ~", paste(vars, collapse = " + ")))
  fit  <- coxph(form, data = df, x = TRUE, y = TRUE)
  print(summary(fit))
  s0  <- get_S0_at(fit, 3652.5); cix <- cindex_safe(fit, df)
  cat(sprintf("S0(10y) %s: %.5f | C-index %s: %.3f\n", label, s0, label, cix))
  fit
}

fit_m <- fit_by_sex(1, "Hombres")
fit_f <- fit_by_sex(0, "Mujeres")

# Modelo combinado (respaldo)
df_all <- core %>%
  select(seqn, sex, age, time_days, status) %>%
  left_join(pred %>% select(seqn, all_of(optional)), by = "seqn")
vars_all <- c("sex", "age", pick_vars(df_all))
form_all <- as.formula(paste("Surv(time_days, status) ~", paste(vars_all, collapse = " + ")))
fit_all  <- coxph(form_all, data = df_all, x = TRUE, y = TRUE)
print(summary(fit_all))
cat(sprintf("S0(10y) Combinado: %.5f | C-index Combinado: %.3f\n",
            get_S0_at(fit_all,3652.5), cindex_safe(fit_all, df_all)))

saveRDS(list(fit_m = fit_m, fit_f = fit_f, fit_all = fit_all), "data/processed/cox_models_ALL.rds")
cat("Modelos guardados en data/processed/cox_models_ALL.rds\n")

# Exportables
out_dir <- "output"; dir.create(out_dir, showWarnings = FALSE)
save_cox <- function(fit, tag){
  if (is.null(fit)) return(invisible(NULL))
  s <- summary(fit)
  tab <- data.frame(
    var     = rownames(s$coef),
    beta    = s$coef[, "coef"],
    HR      = exp(s$coef[, "coef"]),
    HR_low  = exp(s$conf.int[, "lower .95"]),
    HR_high = exp(s$conf.int[, "upper .95"]),
    pvalue  = s$coef[, "Pr(>|z|)"],
    row.names = NULL
  )
  write.csv(tab, file.path(out_dir, paste0("cox_", tag, "_coeffs.csv")), row.names = FALSE)
  write.csv(data.frame(t_days = 3652.5, S0_10y = get_S0_at(fit, 3652.5)),
            file.path(out_dir, paste0("cox_", tag, "_S0_10y.csv")), row.names = FALSE)
}
save_cox(fit_m,  "H")
save_cox(fit_f,  "M")
save_cox(fit_all,"ALL")
cat("Coeficientes y S0(10y) guardados en 'output/'\n")
