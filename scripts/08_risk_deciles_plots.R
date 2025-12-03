# scripts/10_risk_deciles_plots.R
suppressPackageStartupMessages({
  library(dplyr); library(ggplot2); library(readr); library(tidyr)
})

dir.create("output", showWarnings = FALSE)

# 1) Cargar puntajes y outcome a 10 años
scores <- readr::read_csv("output/risk10_scores.csv", show_col_types = FALSE)
core   <- readRDS("data/processed/nhanes_surv_core_9914.rds") %>%
  mutate(time_days = time_mo*(365.25/12),
         y10 = as.integer(time_days <= 3652.5 & status == 1)) %>%
  select(seqn, sex, y10)

dat <- scores %>%
  inner_join(core, by = c("seqn","sex")) %>%
  filter(is.finite(risk10), risk10 >= 0, risk10 <= 1)

# 2) Deciles globales de riesgo
dat <- dat %>%
  mutate(decile = ntile(risk10, 10L),
         sex_fac = factor(sex, levels=c(0,1), labels=c("Mujeres","Hombres")))

# 3) Tabla de calibración por decil: Global y por sexo
mk_tab <- function(df, tag){
  df %>% group_by(decile) %>%
    summarise(
      n = n(),
      risk_pred_mean = mean(risk10, na.rm=TRUE),
      events_obs = sum(y10, na.rm=TRUE),
      rate_obs = events_obs / n,
      .groups = "drop"
    ) %>% mutate(group = tag)
}

tab_all <- mk_tab(dat, "Todos")
tab_m   <- mk_tab(filter(dat, sex == 1L), "Hombres")
tab_f   <- mk_tab(filter(dat, sex == 0L), "Mujeres")
tab_all3 <- bind_rows(tab_all, tab_m, tab_f)

readr::write_csv(tab_all3, "output/deciles_calibration_table.csv")

# 4) Barras lado a lado (observado vs predicho) por decil – Global
tab_all_long <- tab_all %>%
  select(decile, rate_obs, risk_pred_mean) %>%
  pivot_longer(cols = c(rate_obs, risk_pred_mean),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         rate_obs = "Tasa observada",
                         risk_pred_mean = "Riesgo predicho"))

p1 <- ggplot(tab_all_long, aes(factor(decile), value, fill = metric)) +
  geom_col(position = position_dodge(width=.8)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Calibración por deciles (Global)",
       x = "Decil de riesgo (1=bajo … 10=alto)",
       y = "Proporción en 10 años",
       fill = NULL) +
  theme_minimal()
ggsave("output/fig14_deciles_global_bars.png", p1, width=9, height=5, dpi=300)

# 5) Barras lado a lado por decil – Por sexo (facet)
tab_sex_long <- bind_rows(tab_m %>% mutate(group="Hombres"),
                          tab_f %>% mutate(group="Mujeres")) %>%
  select(group, decile, rate_obs, risk_pred_mean) %>%
  pivot_longer(cols = c(rate_obs, risk_pred_mean),
               names_to = "metric", values_to = "value") %>%
  mutate(metric = recode(metric,
                         rate_obs = "Tasa observada",
                         risk_pred_mean = "Riesgo predicho"))

p2 <- ggplot(tab_sex_long, aes(factor(decile), value, fill = metric)) +
  geom_col(position = position_dodge(width=.8)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~group) +
  labs(title = "Calibración por deciles (por sexo)",
       x = "Decil de riesgo (1=bajo … 10=alto)",
       y = "Proporción en 10 años",
       fill = NULL) +
  theme_minimal()
ggsave("output/fig15_deciles_by_sex_bars.png", p2, width=10, height=5.5, dpi=300)

# 6) Línea Observado vs Predicho por decil – Global
p3 <- ggplot(tab_all, aes(decile)) +
  geom_line(aes(y = rate_obs), linewidth=1) +
  geom_point(aes(y = rate_obs), size=2) +
  geom_line(aes(y = risk_pred_mean), linetype=2) +
  geom_point(aes(y = risk_pred_mean), shape=21, fill="white", size=2) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(title = "Observado y predicho por deciles (Global)",
       x = "Decil de riesgo",
       y = "Proporción en 10 años") +
  theme_minimal()
ggsave("output/fig16_deciles_global_lines.png", p3, width=8, height=5, dpi=300)

# 7) Línea Observado vs Predicho – Por sexo (facet)
tab_sex_lines <- bind_rows(
  tab_m %>% mutate(group="Hombres"),
  tab_f %>% mutate(group="Mujeres")
)

p4 <- ggplot(tab_sex_lines, aes(decile)) +
  geom_line(aes(y = rate_obs, color="Observado"), linewidth=1) +
  geom_point(aes(y = rate_obs, color="Observado"), size=2) +
  geom_line(aes(y = risk_pred_mean, color="Predicho"), linetype=2) +
  geom_point(aes(y = risk_pred_mean, color="Predicho"), shape=21, fill="white", size=2) +
  scale_color_manual(values = c("Observado" = "#1f77b4", "Predicho" = "#d62728")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  facet_wrap(~group) +
  labs(title = "Observado vs Predicho por deciles (por sexo)",
       x = "Decil de riesgo", y = "Proporción en 10 años", color = NULL) +
  theme_minimal()
ggsave("output/fig17_deciles_by_sex_lines.png", p4, width=9, height=5.5, dpi=300)

message("Listo:",
        "\n - output/deciles_calibration_table.csv",
        "\n - output/fig14_deciles_global_bars.png",
        "\n - output/fig15_deciles_by_sex_bars.png",
        "\n - output/fig16_deciles_global_lines.png",
        "\n - output/fig17_deciles_by_sex_lines.png")
