# ============================================================
# Dépendances
# ============================================================
library(INLA)        # install.packages("INLA", repos=c(getOption("repos"),
#                     INLA="https://inla.r-inla-download.org/R/stable"))
library(readxl)
library(dplyr)
library(sf)
library(ggplot2)
library(viridis)
library(lubridate)
library(geodata)

# ============================================================
# Paramètres utilisateur
# ============================================================
workdir      <- "C:/Users/Ousmane Diao/Desktop/Training/R-Inla_PFC1"
excel_path   <- file.path(workdir, "df_mal.xlsx")
country_iso3 <- "SEN"      # ex. "CMR" pour Cameroun
crs_proj     <- 3857       # Web Mercator (mètres). UTM local possible.

setwd(workdir)

# ============================================================
# 1) Données : lecture & prétraitement (ordre clé)
# ============================================================
df <- read_excel(excel_path) %>%
  rename(
    cases      = `Nb.total.cas.paludisme.confirmés`,
    population = pop,
    lon        = Longitude,
    lat        = Latitude,
    annee      = Annee,
    mois       = Mois
  ) %>%
  mutate(
    annee   = as.integer(annee),
    mois    = as.integer(mois),
    yearmon = as.Date(sprintf("%d-%02d-01", annee, mois))
  ) %>%
  arrange(yearmon) %>%                              # <-- important : ordre temporel
  mutate(
    timeID = match(yearmon, sort(unique(yearmon)))  # index 1..n_time SANS trous
  )

# ------------------------------------------------------------
# Covariables : standardisation (optionnel)
# ------------------------------------------------------------
covariables <- c("ndvi", "Humidite", "Precipitation", "Temperature", "population")
df <- df %>%
  mutate(across(all_of(covariables), ~ suppressWarnings(as.numeric(.)))) %>%
  mutate(
    ndvi_s      = scale(ndvi)[,1],
    humidite_s  = scale(Humidite)[,1],
    precip_s    = scale(Precipitation)[,1],
    temp_s      = scale(Temperature)[,1]
  )

# Remplacer les éventuels NA/Inf post-standardisation par 0 (moyenne centrée)
df <- df %>%
  mutate(
    ndvi_s     = ifelse(is.finite(ndvi_s), ndvi_s, 0),
    humidite_s = ifelse(is.finite(humidite_s), humidite_s, 0),
    precip_s   = ifelse(is.finite(precip_s), precip_s, 0),
    temp_s     = ifelse(is.finite(temp_s), temp_s, 0)
  )

# Offset d'exposition (comptes Poisson) : log(population)
df$offset_log_pop <- log(pmax(df$population, 1e-9))

# ============================================================
# 2) Géométries & projection
#    (Créer les points APRES l'ordonnancement pour garder l'alignement)
# ============================================================
pts_ll <- st_as_sf(df, coords = c("lon", "lat"), crs = 4326)   # points lon/lat
pts    <- st_transform(pts_ll,  crs_proj)                      # en mètres

# Frontière pays (niveau 0)
adm0_ll <- gadm(country = country_iso3, level = 0, path = workdir) |> st_as_sf()
adm0    <- st_transform(adm0_ll, crs_proj)

# Pour INLA: segment de frontière (objet 'sp')
adm0_sp          <- as(adm0, "Spatial")
boundary_segment <- inla.sp2segment(adm0_sp)
boundary_segment$loc <- inla.mesh.map(boundary_segment$loc)

# Matrice de coordonnées (X,Y) en mètres correspondant EXACTEMENT à df
xy         <- st_coordinates(pts)
coords_mat <- as.matrix(xy[, 1:2])

# ============================================================
# 3) Maillage SPDE & modèle Matérn
# ============================================================
mesh <- inla.mesh.2d(
  loc      = coords_mat,
  boundary = boundary_segment,
  max.edge = c(50e3, 150e3),  # ~50 km & 150 km (adapte au pays)
  cutoff   = 10e3,            # ~10 km (évite triangles trop serrés)
  offset   = c(30e3, 100e3)   # marge hors frontière
)
cat("Nb. de sommets du mesh:", mesh$n, "\n")

# ============================================================
# (Optionnel) Visualisation rapide du mesh et des points
# ============================================================
plot(mesh, main = "Maillage SPDE")
points(coords_mat[,1], coords_mat[,2], col = "red", pch = 16, cex = 0.5)


spde <- inla.spde2.matern(mesh = mesh, alpha = 2)
# (Option plus stable : spde <- inla.spde2.pcmatern(mesh, alpha=2,
#        prior.range=c(100e3,0.5), prior.sigma=c(1,0.01)) )

# (Optionnel) Visualisation rapide
# plot(mesh, main = "Maillage SPDE"); points(coords_mat[,1], coords_mat[,2], col="red", pch=16, cex=0.5)

# ============================================================
# 4) Index spatio-temporel & matrice A
# ============================================================
n_time <- length(unique(df$timeID))

spde_index <- inla.spde.make.index(
  name    = "spatial",
  n.spde  = spde$n.spde,
  n.group = n_time
)

A_obs <- inla.spde.make.A(
  mesh   = mesh,
  loc    = coords_mat,
  group  = df$timeID,
  n.group = n_time
)

# Garde-fous
stopifnot(nrow(A_obs) == nrow(df))
stopifnot(all(df$timeID >= 1 & df$timeID <= n_time))
stopifnot(ncol(A_obs) == spde$n.spde * n_time)

# ============================================================
# 5) Stack INLA (estimation)
# ============================================================
stk_obs <- inla.stack(
  data = list(y = as.integer(df$cases)),
  A    = list(A_obs, 1),
  effects = list(
    spatial = spde_index,   # fournit 'spatial' et 'spatial.group'
    data.frame(
      Intercept      = 1,
      timeID         = df$timeID,
      ndvi_s         = df$ndvi_s,
      humidite_s     = df$humidite_s,
      precip_s       = df$precip_s,
      temp_s         = df$temp_s,
      offset_log_pop = df$offset_log_pop
    )
  ),
  tag = "dat"
)

# ============================================================
# 6) Modèle & ajustement INLA
# ============================================================
formula <- y ~ -1 + Intercept +
  ndvi_s + humidite_s + precip_s + temp_s +
  offset(offset_log_pop) +
  f(spatial, model = spde, group = spatial.group,
    control.group = list(model = "ar1")) +
  f(timeID, model = "rw1")

fit <- INLA::inla(
  formula, family = "poisson",
  data = inla.stack.data(stk_obs),
  control.predictor = list(A = inla.stack.A(stk_obs), compute = TRUE),
  control.compute   = list(dic = TRUE, waic = TRUE, cpo = TRUE, config = TRUE),
  # control.inla    = list(strategy = "simplified.laplace", int.strategy = "eb"), # option robustesse
  verbose = FALSE
)

# Résumés
print(fit$summary.fixed)
print(fit$summary.hyperpar)
cat("DIC:", fit$dic$dic, "  WAIC:", fit$waic$waic, "\n")

