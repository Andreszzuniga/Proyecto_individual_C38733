# scripts/01_build_predictors_ALL.R
# Construye PREDICTORES (ricos) y NÚCLEO de supervivencia (robusto) para NHANES 1999–2014

suppressPackageStartupMessages({
  library(dplyr)
  library(nhanesA)
  library(purrr)
})

options(timeout = 600)

# --- Mortalidad enlazada (hecha previamente con 00_get_mortality_NHANES.R) ---
lmf <- readRDS("data/processed/nhanes_lmf_1999_2014.rds")

nh_get <- function(code) tryCatch(nhanes(code), error = function(e) NULL)

cycles <- list(
  A = list(suf = "",    years = "1999-2000"),
  B = list(suf = "_B",  years = "2001-2002"),
  C = list(suf = "_C",  years = "2003-2004"),
  D = list(suf = "_D",  years = "2005-2006"),
  E = list(suf = "_E",  years = "2007-2008"),
  F = list(suf = "_F",  years = "2009-2010"),
  G = list(suf = "_G",  years = "2011-2012"),
  H = list(suf = "_H",  years = "2013-2014")
)

# -------- Predictores “ricos” por ciclo (pueden tener NA) --------------------
pull_cycle <- function(suf){
  demo <- nh_get(paste0("DEMO", suf)); if (is.null(demo)) return(NULL)
  demo <- demo %>% select(SEQN, RIAGENDR, RIDAGEYR)
  
  bp   <- nh_get(paste0("BPX", suf));  if(!is.null(bp))  bp  <- bp  %>% select(SEQN, BPXSY1, BPXSY2, BPXSY3)
  smk  <- nh_get(paste0("SMQ", suf));  if(!is.null(smk)) smk <- smk %>% select(SEQN, SMQ020)
  bmx  <- nh_get(paste0("BMX", suf));  if(!is.null(bmx)) bmx <- bmx %>% select(SEQN, BMXBMI)
  bpq  <- nh_get(paste0("BPQ", suf));  if(!is.null(bpq)) bpq <- bpq %>% select(SEQN, BPQ050A)
  
  # Colesterol total (LBXTC) – probamos varios códigos de tabla
  lip <- NULL
  for(tb in c(paste0("L13",suf), paste0("TCHOL",suf), paste0("LAB13",suf))){
    tmp <- nh_get(tb)
    if(!is.null(tmp) && "LBXTC" %in% names(tmp)){
      lip <- tmp[, c("SEQN","LBXTC")]
      break
    }
  }
  if (is.null(lip)) return(NULL)  # sin colesterol, no seguimos
  
  # HDL – si no existe en el ciclo, lo dejamos como NA (no abortar)
  hdl <- NULL
  for(tb in c(paste0("HDL",suf), paste0("L13",suf), paste0("LAB13",suf))){
    tmp <- nh_get(tb); if (is.null(tmp)) next
    hit <- intersect(c("LBDHDD","LBXHDD","LBXHDL","LBDHDL"), names(tmp))
    if(length(hit) >= 1){
      v <- hit[1]
      hdl <- tmp[, c("SEQN", v)]
      names(hdl) <- c("SEQN","LBDHDD")
      break
    }
  }
  
  # Diabetes (opcional)
  diab_raw <- nh_get(paste0("DIQ", suf))
  diab <- if(!is.null(diab_raw) && "DIQ010" %in% names(diab_raw)){
    diab_raw %>% select(SEQN, DIQ010) %>%
      mutate(diabetes = case_when(DIQ010 == 1 ~ 1L,
                                  DIQ010 %in% c(2,3) ~ 0L,
                                  TRUE ~ NA_integer_)) %>%
      select(SEQN, diabetes)
  } else NULL
  
  pred <- demo %>%
    { if(!is.null(bp))  left_join(., bp,  by="SEQN") else . } %>%
    left_join(lip, by = "SEQN") %>%
    { if(!is.null(hdl)) left_join(., hdl, by="SEQN") else mutate(., LBDHDD = NA_real_) } %>%
    { if(!is.null(smk)) left_join(., smk, by="SEQN") else . } %>%
    { if(!is.null(bmx)) left_join(., bmx, by="SEQN") else . } %>%
    { if(!is.null(bpq)) left_join(., bpq, by="SEQN") else . }
  
  if(!is.null(diab)) pred <- left_join(pred, diab, by="SEQN") else pred$diabetes <- NA_integer_
  
  pred %>% transmute(
    seqn = SEQN,
    sex  = dplyr::case_when(
      RIAGENDR %in% c(1, "1", "Male", "MALE", "Hombre", "M") ~ 1L,
      RIAGENDR %in% c(2, "2", "Female", "FEMALE", "Mujer", "F") ~ 0L,
      TRUE ~ NA_integer_
    ),
    age  = RIDAGEYR,
    sbp  = rowMeans(cbind(BPXSY1,BPXSY2,BPXSY3), na.rm=TRUE),
    chol = LBXTC,
    hdl  = LBDHDD,
    bmi  = BMXBMI,
    smoker = ifelse(SMQ020==1, 1L, 0L),
    diabetes = diabetes,
    htn_treated = ifelse(BPQ050A==1, 1L, 0L)
  )
}

pred_all_full <- purrr::map(cycles, ~ pull_cycle(.x$suf)) %>% bind_rows(.id = "cycle")

# -------- Núcleo DEMO para no perder hombres (sexo/edad) ---------------------
demo_all <- purrr::map(cycles, function(x){
  dm <- nh_get(paste0("DEMO", x$suf)); if (is.null(dm)) return(NULL)
  dm %>%
    dplyr::select(SEQN, RIAGENDR, RIDAGEYR) %>%
    dplyr::transmute(
      seqn = SEQN,
      sex  = dplyr::case_when(
        RIAGENDR %in% c(1, "1", "Male", "MALE", "Hombre", "M") ~ 1L,
        RIAGENDR %in% c(2, "2", "Female", "FEMALE", "Mujer", "F") ~ 0L,
        TRUE ~ NA_integer_
      ),
      age  = RIDAGEYR
    )
}) %>% dplyr::bind_rows()

# -------- Núcleo de supervivencia: LMF + DEMO (mínimo filtro) ---------------
dat_core <- lmf %>%
  dplyr::inner_join(demo_all, by = "seqn") %>%
  dplyr::filter(!is.na(sex), !is.na(age), !is.na(time_mo), !is.na(status))

dir.create("data/processed", showWarnings = FALSE)
saveRDS(pred_all_full, "data/processed/nhanes_pred_all_full_9914.rds")
saveRDS(dat_core,      "data/processed/nhanes_surv_core_9914.rds")

cat("OK -> guardados:\n",
    "  data/processed/nhanes_pred_all_full_9914.rds  (predictores)\n",
    "  data/processed/nhanes_surv_core_9914.rds      (núcleo supervivencia)\n")


