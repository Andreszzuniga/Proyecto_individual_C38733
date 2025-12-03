# scripts/09_risk_plots.R
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(tidyr); library(survival)
})

dir.create("output", showWarnings = FALSE)

# ---------- 0) Carga datos + funciones de scoring ----------
core <- readRDS("data/processed/nhanes_surv_core_9914.rds") %>%
  mutate(time_days = time_mo*(365.25/12))

pred <- readRDS("data/processed/nhanes_pred_all_full_9914.rds")

dat <- core %>%
  select(seqn, sex, age, time_days, status) %>%
  left_join(pred %>% select(seqn, sbp, chol, hdl, bmi, smoker, diabetes, htn_treated), by="seqn")

# Tipos seguros
numfix <- c("age","sbp","chol","hdl","bmi")
for(v in numfix) if(v %in% names(dat)) dat[[v]] <- suppressWarnings(as.numeric(dat[[v]]))
for(v in c("sex","smoker","diabetes","htn_treated")){
  if(v %in% names(dat)) dat[[v]] <- suppressWarnings(as.integer(dat[[v]]))
}

# Carga funciones de riesgo (usa tus modelos entrenados)
source("scripts/03_score.R")  # define risk_10y() y predict_risk_df()

# Aplica riesgo 10 años solo a filas con columnas requeridas
need <- c("sex","age","sbp","chol","hdl","bmi")
have <- intersect(need, names(dat))
if(length(have) < length(need)){
  stop("Faltan columnas para riesgo: necesito ", paste(need, collapse=", "))
}

df_sc <- dat %>% filter(!is.na(sex), !is.na(age), !is.na(sbp), !is.na(chol), !is.na(hdl), !is.na(bmi))
df_sc$risk10 <- predict_risk_df(df_sc)

# ---------- 1) Barras: riesgo medio por bandas de edad y sexo ----------
brks <- seq(20, 90, by=10)
df_sc <- df_sc %>%
  mutate(sex_fac = factor(sex, levels=c(0,1), labels=c("Mujeres","Hombres")),
         age_band = cut(age, breaks = brks, right = FALSE, include.lowest = TRUE))

risk_age_sex <- df_sc %>%
  group_by(age_band, sex_fac) %>%
  summarise(n = n(),
            risk_mean = mean(risk10, na.rm=TRUE),
            .groups="drop")

ggplot(risk_age_sex, aes(age_band, risk_mean, fill=sex_fac)) +
  geom_col(position=position_dodge(width=.8)) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  labs(title="Riesgo a 10 años por bandas de edad y sexo",
       x="Rango de edad", y="Riesgo promedio (10 años)", fill="Sexo") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=25, hjust=1))
ggsave("output/fig10_risk_bar_age_sex.png", width=9, height=5, dpi=300)

# ---------- 2) Líneas: riesgo vs edad a distintos niveles de SBP ----------
# Medias por sexo para fijar covariables (chol/hdl/bmi; y opcional: smoker/diabetes/htn_treated)
means_by_sex <- df_sc %>%
  group_by(sex) %>%
  summarise(across(c(chol, hdl, bmi, smoker, diabetes, htn_treated), ~mean(.x, na.rm=TRUE)),
            .groups="drop")

ages_grid <- seq(30, 85, by=1)
sbp_levels <- c(120, 140, 160)

grid <- tidyr::expand_grid(
  sex = c(0L,1L),
  age = ages_grid,
  sbp = sbp_levels
) %>%
  left_join(means_by_sex, by="sex") %>%
  mutate(
    sex_fac = factor(sex, levels=c(0,1), labels=c("Mujeres","Hombres"))
  )

grid$risk10 <- predict_risk_df(grid %>%
                                 transmute(sex, age, sbp, chol, hdl, bmi,
                                           smoker = ifelse(is.na(smoker), 0L, round(smoker)),
                                           diabetes = ifelse(is.na(diabetes), 0L, round(diabetes)),
                                           htn_treated = ifelse(is.na(htn_treated), 0L, round(htn_treated)))
)

ggplot(grid, aes(age, risk10, color=factor(sbp))) +
  geom_line(size=1) +
  facet_wrap(~sex_fac) +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  labs(title="Riesgo a 10 años vs Edad a distintos niveles de SBP",
       x="Edad (años)", y="Riesgo (10 años)", color="SBP (mmHg)") +
  theme_minimal()
ggsave("output/fig11_risk_vs_age_lines_by_sbp.png", width=9, height=5.5, dpi=300)

# ---------- 3) Serie “temporal”: S0(t) (supervivencia base) vs tiempo ----------
mods <- readRDS("data/processed/cox_models_ALL.rds")
fits <- list(
  Mujeres = mods$fit_f,
  Hombres = mods$fit_m
)
# Si falta hombres, ignora; si ambos faltan, intenta combinado leyendo del 02_fit (opcional)

s0_list <- list()
for(lbl in names(fits)){
  fit <- fits[[lbl]]
  if(!is.null(fit)){
    bh <- basehaz(fit, centered = TRUE)  # tiempo en días
    s0 <- dplyr::tibble(
      sex_lbl = lbl,
      years = bh$time / 365.25,
      S0 = exp(-bh$hazard)
    )
    s0_list[[lbl]] <- s0
  }
}
if(length(s0_list)){
  s0_all <- bind_rows(s0_list)
  ggplot(s0_all, aes(years, S0, color=sex_lbl)) +
    geom_line(size=1) +
    scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
    labs(title="Supervivencia base S0(t) por sexo",
         x="Tiempo (años)", y="S0(t)", color="Modelo") +
    theme_minimal()
  ggsave("output/fig12_baseline_survival_curves.png", width=8, height=5, dpi=300)
}

# ---------- 4) Línea (Kaplan–Meier) por sexo ----------
km <- survfit(Surv(time_days, status) ~ sex, data = core)
km_df <- do.call(rbind, lapply(seq_along(km$strata), function(i){
  idx <- sum(km$strata[seq_len(i)]) - km$strata[i] + seq_len(km$strata[i])
  data.frame(
    sex = names(km$strata)[i],
    time_years = km$time[idx]/365.25,
    surv = km$surv[idx]
  )
}))
km_df$sex <- factor(km_df$sex, labels=c("Mujeres","Hombres"))

ggplot(km_df, aes(time_years, surv, color=sex)) +
  geom_step() +
  scale_y_continuous(labels=scales::percent_format(accuracy = 1)) +
  labs(title="Curvas de Kaplan–Meier por sexo",
       x="Tiempo (años)", y="Supervivencia observada", color="Sexo") +
  theme_minimal()
ggsave("output/fig13_km_by_sex.png", width=8, height=5, dpi=300)

message("Listo: fig10_risk_bar_age_sex.png, fig11_risk_vs_age_lines_by_sbp.png, fig12_baseline_survival_curves.png, fig13_km_by_sex.png en /output")
