# scripts/06
#_eda_relations.R
# EDA con ggplot: relaciones entre edad/sexo y biomarcadores (sbp, chol, hdl, bmi)

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  ggally_ok <- requireNamespace("GGally", quietly = TRUE)
})

dir.create("output", showWarnings = FALSE)

# -------------------------------------------------------------------
# 1) Cargar bases ya generadas en pasos previos
# -------------------------------------------------------------------
core <- readRDS("data/processed/nhanes_surv_core_9914.rds")  # contiene: seqn, sex, age, time_mo, status
pred <- readRDS("data/processed/nhanes_pred_all_full_9914.rds")  # contiene predictores ricos por seqn

# -------------------------------------------------------------------
# 2) Armonizar tipos y coalesce sin choques de factor/integer
# -------------------------------------------------------------------
# Variables que usaremos
want_cont <- c("age","sbp","chol","hdl","bmi")
want_bin  <- c("sex","smoker","diabetes","htn_treated")
want_all  <- c("seqn", want_cont, want_bin)

# Join base (puede crear .x/.y)
dat_eda <- core %>% left_join(pred, by = "seqn")

# Helpers
na_like <- function(kind) if (kind == "int") NA_integer_ else NA_real_

coerce_numeric <- function(z) {
  # convierte factor/char a numeric de forma segura
  suppressWarnings(as.numeric(as.character(z)))
}
coerce_integer <- function(z) {
  # integer consistente desde factor/char/numeric
  suppressWarnings(as.integer(as.numeric(as.character(z))))
}

getcol <- function(d, nm) if (nm %in% names(d)) d[[nm]] else NULL

# Coalesce tipado (tolerante a NULL)
coalesce3_typed <- function(a, b, c, kind = c("num","int")){
  kind <- match.arg(kind)
  cast <- if (kind == "int") coerce_integer else coerce_numeric
  if (is.null(a)) a <- na_like(if (kind=="int") "int" else "num")
  if (is.null(b)) b <- na_like(if (kind=="int") "int" else "num")
  if (is.null(c)) c <- na_like(if (kind=="int") "int" else "num")
  dplyr::coalesce(cast(a), cast(b), cast(c))
}

# Construye columna canónica para cada variable deseada
for (v in want_all) {
  x <- paste0(v, ".x")
  y <- paste0(v, ".y")
  kind <- if (v %in% want_bin) "int" else "num"
  dat_eda[[v]] <- coalesce3_typed(getcol(dat_eda, x),
                                  getcol(dat_eda, y),
                                  getcol(dat_eda, v),
                                  kind = if (kind == "int") "int" else "num")
}

# Nos quedamos solo con las columnas canónicas
dat_eda <- dat_eda %>% select(any_of(want_all))

# Tipos finales y etiquetas para graficar
dat_eda <- dat_eda %>%
  mutate(
    across(all_of(want_cont), as.numeric),
    smoker      = as.integer(smoker),
    diabetes    = as.integer(diabetes),
    htn_treated = as.integer(htn_treated),
    sex         = factor(sex, levels = c(0,1), labels = c("Mujeres","Hombres"))
  )

message("Columns after fix: ", paste(names(dat_eda), collapse = ", "))

# -------------------------------------------------------------------
# 3) Utilidad para guardar solo si existen columnas requeridas
# -------------------------------------------------------------------
save_plot_if <- function(cols, plot_fun, filename, width=7, height=5, dpi=300){
  if (all(cols %in% names(dat_eda))) {
    p <- plot_fun()
    ggsave(file.path("output", filename), p, width=width, height=height, dpi=dpi)
    message("OK -> output/", filename)
  } else {
    message("Omitido (faltan columnas): ", filename,
            "  | requiere: ", paste(cols, collapse=", "),
            "  | tengo: ", paste(names(dat_eda), collapse=", "))
  }
}

# -------------------------------------------------------------------
# 4) Figuras
# -------------------------------------------------------------------

# 4.1 Edad por sexo (violin + box)
save_plot_if(c("sex","age"), function(){
  ggplot(dat_eda, aes(sex, age, fill=sex)) +
    geom_violin(trim=TRUE, alpha=.40) +
    geom_boxplot(width=.15, outlier.alpha=.25) +
    labs(title="Distribución de Edad por Sexo", x=NULL, y="Edad (años)") +
    theme_minimal() + theme(legend.position="none")
}, "fig05_age_by_sex.png")

# 4.2 Colesterol vs Edad (hexbin si hay muchos puntos; sino scatter + smooth), por sexo
save_plot_if(c("age","chol","sex"), function(){
  df <- dat_eda %>% filter(!is.na(age), !is.na(chol))
  if (nrow(df) > 30000) {
    ggplot(df, aes(age, chol)) +
      geom_bin2d() +
      facet_wrap(~sex) +
      scale_fill_continuous(name="Conteo") +
      labs(title="Colesterol (LBXTC) vs Edad (hexbin)", x="Edad", y="Colesterol") +
      theme_minimal()
  } else {
    ggplot(df, aes(age, chol, color=sex)) +
      geom_point(alpha=.25, size=1) +
      geom_smooth(se=FALSE) +
      labs(title="Colesterol (LBXTC) vs Edad por Sexo", x="Edad", y="Colesterol", color=NULL) +
      theme_minimal()
  }
}, "fig06_age_vs_chol.png", width=8, height=5.5)

# 4.3 Edad vs {SBP, HDL, BMI} con suavizado por sexo (formato largo)
save_plot_if(c("age","sex"), function(){
  vars <- intersect(c("sbp","hdl","bmi"), names(dat_eda))
  stopifnot(length(vars) >= 1)
  long <- dat_eda %>%
    select(age, sex, all_of(vars)) %>%
    pivot_longer(cols = all_of(vars), names_to="var", values_to="value") %>%
    filter(!is.na(age), !is.na(value))
  labs_y <- c(sbp="Presión Sistólica (SBP)", hdl="HDL", bmi="IMC")
  ggplot(long, aes(age, value, color=sex)) +
    geom_point(alpha=.12, size=.6) +
    geom_smooth(se=FALSE) +
    facet_wrap(~var, scales="free_y", labeller = as_labeller(labs_y)) +
    labs(title="Relaciones con la Edad por Sexo", x="Edad (años)", y=NULL, color=NULL) +
    theme_minimal()
}, "fig07_age_vs_sbp_hdl_bmi.png", width=10, height=5.5)

# 4.4 Densidades por sexo (chol, hdl, bmi, sbp)
save_plot_if(c("sex"), function(){
  dens_vars <- intersect(c("chol","hdl","bmi","sbp"), names(dat_eda))
  long <- dat_eda %>%
    select(sex, all_of(dens_vars)) %>%
    pivot_longer(-sex, names_to="var", values_to="value") %>%
    filter(!is.na(value))
  nice <- c(chol="Colesterol", hdl="HDL", bmi="IMC", sbp="SBP")
  ggplot(long, aes(value, fill=sex)) +
    geom_density(alpha=.35) +
    facet_wrap(~var, scales="free", labeller = as_labeller(nice)) +
    labs(title="Distribuciones por Sexo", x=NULL, y="Densidad", fill=NULL) +
    theme_minimal()
}, "fig08_densities_by_sex.png", width=10, height=6)

# 4.5 Boxplots por sexo (chol, hdl, bmi, sbp)
save_plot_if(c("sex"), function(){
  bx_vars <- intersect(c("chol","hdl","bmi","sbp"), names(dat_eda))
  long <- dat_eda %>%
    select(sex, all_of(bx_vars)) %>%
    pivot_longer(-sex, names_to="var", values_to="value")
  nice <- c(chol="Colesterol", hdl="HDL", bmi="IMC", sbp="SBP")
  ggplot(long, aes(sex, value, fill=sex)) +
    geom_boxplot(outlier.alpha=.2, width=.5) +
    facet_wrap(~var, scales="free_y", labeller=as_labeller(nice)) +
    labs(title="Comparación por Sexo", x=NULL, y=NULL) +
    theme_minimal() + theme(legend.position="none")
}, "fig09_box_by_sex.png", width=10, height=6)

# 4.6 Matriz de correlaciones (numéricas)
num_vars <- intersect(c("age","sbp","chol","hdl","bmi"), names(dat_eda))
if (length(num_vars) >= 2) {
  df_num <- dat_eda %>%
    select(all_of(num_vars)) %>%
    mutate(across(everything(), as.numeric))
  if (ggally_ok) {
    p_corr <- GGally::ggcorr(df_num, label=TRUE, label_size=3, hjust=.9, layout.exp=1) +
      ggtitle("Correlaciones (numéricas)")
  } else {
    cors <- cor(df_num, use="pairwise.complete.obs")
    cor_df <- as.data.frame(as.table(cors)); names(cor_df) <- c("Var1","Var2","Corr")
    p_corr <- ggplot(cor_df, aes(Var1, Var2, fill=Corr)) +
      geom_tile() + geom_text(aes(label=sprintf("%.2f", Corr)), size=3) +
      scale_fill_gradient2(limits=c(-1,1)) +
      labs(title="Correlaciones (numéricas)", x=NULL, y=NULL) +
      theme_minimal()
  }
  ggsave("output/eda_correlation_matrix.png", p_corr, width=6, height=5, dpi=300)
  message("OK -> output/eda_correlation_matrix.png")
} else {
  message("Omitido matriz de correlaciones (num_vars < 2).")
}

message("Listo. Revisa los PNG en /output y la consola por cualquier 'Omitido'.")
