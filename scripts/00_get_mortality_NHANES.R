# scripts/00_get_mortality_NHANES.R
# Descarga y lectura de LMF (NHANES 1999–2014) + creación de Surv(10 años, ictus)
suppressPackageStartupMessages({ library(readr); library(dplyr); library(fs) })

# Carpetas
dir_create("data/raw/nhanes_lmf")
dir_create("data/processed")

# Archivos a usar (NHANES, públicos con mortalidad 2019)
base <- "https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/"
files <- c(
  "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
  "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
  "NHANES_2003_2004_MORT_2019_PUBLIC.dat",
  "NHANES_2005_2006_MORT_2019_PUBLIC.dat",
  "NHANES_2007_2008_MORT_2019_PUBLIC.dat",
  "NHANES_2009_2010_MORT_2019_PUBLIC.dat",
  "NHANES_2011_2012_MORT_2019_PUBLIC.dat",
  "NHANES_2013_2014_MORT_2019_PUBLIC.dat"
)
targets <- file.path("data/raw/nhanes_lmf", files)

# Descarga (si faltan)
for(i in seq_along(files)){
  url <- paste0(base, files[i]); dest <- targets[i]
  if(!fs::file_exists(dest)){
    message("Descargando: ", files[i])
    try(utils::download.file(url, destfile = dest, mode = "wb", quiet = TRUE), silent = TRUE)
  } else {
    message("Ya existe: ", files[i])
  }
}

# Lector con el layout NHANES del programa oficial (fwf)
read_lmf_nhanes <- function(path){
  read_fwf(file = path,
           col_types = "iiiiiiii",
           fwf_cols(
             seqn         = c(1,6),
             eligstat     = c(15,15),
             mortstat     = c(16,16),
             ucod_leading = c(17,19),
             diabetes     = c(20,20),
             hyperten     = c(21,21),
             permth_int   = c(43,45),
             permth_exm   = c(46,48)
           ),
           na = c("", "."))
}

# Leer y unir
dfs <- lapply(targets, function(p){
  if(!fs::file_exists(p)) return(NULL)
  tryCatch(read_lmf_nhanes(p),
           error=function(e){ message("Error leyendo ", p, ": ", e$message); NULL })
})
dfs <- Filter(Negate(is.null), dfs)
stopifnot(length(dfs) > 0)

lmf <- bind_rows(dfs, .id = "source_file")

# Elegibles adultos y construcción de Surv(10 años) para ictus (UCOD_LEADING==5)
lmf <- lmf %>%
  filter(eligstat == 1L) %>%
  mutate(
    time_mo = pmin(permth_exm, 120L),          # meses desde examen, truncado a 10 años
    stroke_cause = (ucod_leading == 5L),       # 5 = cerebrovascular (I60–I69)
    status = as.integer(mortstat == 1L & stroke_cause & permth_exm <= 120L)
  )

# Guardar
out_path <- "data/processed/nhanes_lmf_1999_2014.rds"
saveRDS(lmf, out_path)

# Chequeos rápidos
cat("Filas:", nrow(lmf), "\n")
cat("Defunciones (cualquier causa):", sum(lmf$mortstat == 1L, na.rm = TRUE), "\n")
cat("Defunciones por ictus <=10a:", sum(lmf$status == 1L, na.rm = TRUE), "\n")
cat("SEQN únicos:", length(unique(lmf$seqn)), "\n")
cat("Guardado en:", out_path, "\n")

