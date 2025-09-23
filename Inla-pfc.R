# Dépendances
library(INLA)
library(inlabru)
library(sf)
library(sp)
library(dplyr)
library(lubridate)
library(ggplot2)
library(viridis)

# --- Chargement données et shapefile ---
df <- readxl::read_excel("df_inla1.xlsx")

# Renommer colonnes si besoin (adapter si noms différents)
df <- df %>%
  rename(cases = `Total.cas.confirmés`,
         population = pop,
         lon = Longitude,
         lat = Latitude,
         annee = Annee,    # suppose colonnes numériques existantes
         mois = Mois)

# lire le shapefile couvrant la zone (remplacer par ton fichier)
shp <- st_read("Districts_sanitaires_Senegal_correct.shp", quiet = TRUE)

# --- Prétraitement temporel à partir de annee/mois numériques ---
# s'assurer que annee et mois sont numériques
df <- df %>%
  mutate(annee = as.integer(annee),
         mois  = as.integer(mois))

# Créer une date "année-mois" et un index mensuel consécutif (timeID)
df <- df %>%
  mutate(yearmon = as.Date(paste0(annee, "-", sprintf("%02d", mois), "-01")))

# Option 1 (simple, garanti consécutif même s'il manque mois) :
min_year <- min(df$annee, na.rm = TRUE)
df <- df %>%
  mutate(timeID = (annee - min_year) * 12 + mois)

# --- Géométrie et projection cohérente ---
# transformer shapefile et données dans la même CRS métrique (ex: UTM zone 28N EPSG:32628)
target_crs <- 32628
shp <- st_transform(shp, crs = target_crs)

df_sf <- st_as_sf(df, coords = c("lon","lat"), crs = 4326, remove = FALSE)
df_sf <- st_transform(df_sf, crs = target_crs)
coords <- st_coordinates(df_sf)
df_sf$x <- coords[,1]; df_sf$y <- coords[,2]

# --- Maillage SPDE contraint au shapefile ---
# Convertir sf -> Spatial pour INLA si nécessaire
shp_sp <- as(shp, "Spatial")   # require sp compatibility
boundary_segment <- inla.sp2segment(shp_sp)
boundary_segment
# Construire un mesh qui respecte le polygone du shapefile
mesh <- inla.mesh.2d(boundary = boundary_segment,
                     loc = cbind(df_sf$x, df_sf$y),
                     max.edge = c(15000, 40000),
                     cutoff = 5000,
                     max.n = 3000)

cat("Nombre de vertices du mesh :", mesh$n, "\n")
plot(mesh); points(df_sf$x, df_sf$y, col='red', pch=16, cex=0.6)

spde <- inla.spde2.pcmatern(mesh=mesh, alpha=2,
                            prior.range = c(20000, 0.5),
                            prior.sigma = c(1, 0.01))

# --- Gestion du temps (observations + prédiction) ---
n_time_obs <- length(unique(df_sf$timeID))
next_n_months <- 24
n_time <- n_time_obs + next_n_months
n_time
spde_index <- inla.spde.make.index(name = "spatial",
                                   n.spde = spde$n.spde,
                                   n.group = n_time)

# Matrice A pour observations (groupe = timeID)
A_obs <- inla.spde.make.A(mesh = mesh,
                          loc = cbind(df_sf$x, df_sf$y),
                          group = df_sf$timeID,
                          n.group = n_time)

# --- Covariables : s'assurer qu'elles sont dans df_sf ---
# Exemple : ndvi, humidite, Precipitation, Temperature
# Si noms différents, adapter ci-dessous.
# Remplir les NA ou vérifier distribution avant modélisation.
covariates <- c("ndvi", "Humidite", "Precipitation", "Temperature")
missing_covs <- setdiff(covariates, names(df_sf))
if(length(missing_covs)>0) stop("Colonnes covariables manquantes : ", paste(missing_covs, collapse=", "))

# Standardiser optionnellement (recommandé) pour stabilité numérique
df_sf <- df_sf %>%
  mutate(ndvi_s = scale(ndvi)[,1],
         humidite_s = scale(Humidite)[,1],
         precip_s = scale(Precipitation)[,1],
         temp_s = scale(Temperature)[,1])

# --- Stack estimation : inclure covariables comme effets fixes ---
offset_obs <- log(df_sf$population / 1000 + 1e-9)

stk_obs <- inla.stack(
  data = list(y = df_sf$cases),
  A = list(A_obs, 1),
  effects = list(spatial = spde_index,
                 data.frame(Intercept = 1,
                            timeID = df_sf$timeID,
                            population = df_sf$population,
                            ndvi_s = df_sf$ndvi_s,
                            humidite_s = df_sf$humidite_s,
                            precip_s = df_sf$precip_s,
                            temp_s = df_sf$temp_s)),
  tag = "estimation"
)

# --- Formule : covariables + effet spatial + effet temporel RW1 ---
formula <- y ~ -1 + Intercept +
  ndvi_s + humidite_s + precip_s + temp_s +
  f(spatial, model = spde, group = spatial.group,
    control.group = list(model = "ar1")) +
  f(timeID, model = "rw1")

# --- Options INLA et exécution ---
control_predictor <- list(A = inla.stack.A(stk_obs), compute = TRUE, link = 1)
control_compute <- list(dic = TRUE, waic = TRUE, cpo = F)

res <- inla(
  formula,
  family = "nbinomial",     # <-- on change juste ici
  data = inla.stack.data(stk_obs),
  control.predictor = control_predictor,
  control.compute   = control_compute,
  control.inla      = list(strategy = "simplified.laplace"),
  offset = offset_obs,
  verbose = TRUE
)


qaQQprint(res$summary.fixed)
print(res$summary.hyperpar)

# --- Préparations prédiction future ---
# Grille de prédiction sur l'emprise du shapefile (en CRS métrique)
bb <- st_bbox(st_transform(shp, st_crs(df_sf)))
grid_spacing <- 5000
x_seq <- seq(bb$xmin, bb$xmax, by = grid_spacing)
y_seq <- seq(bb$ymin, bb$ymax, by = grid_spacing)
pp <- expand.grid(x = x_seq, y = y_seq)
pp_sf <- st_as_sf(pp, coords = c("x","y"), crs = st_crs(df_sf))

# conserver uniquement points dans le shapefile (si souhaité)
pp_sf <- pp_sf[st_intersects(pp_sf, shp, sparse = FALSE), , drop = FALSE]

# pour les covariables futures : ici on met la moyenne observée (OPTION)
mean_ndvi <- mean(df_sf$ndvi_s, na.rm = TRUE)
mean_hum  <- mean(df_sf$humidite_s, na.rm = TRUE)
mean_prec <- mean(df_sf$precip_s, na.rm = TRUE)
mean_temp <- mean(df_sf$temp_s, na.rm = TRUE)
mean_pop  <- mean(df_sf$population, na.rm = TRUE)

future_timeIDs <- seq(n_time_obs + 1, n_time_obs + next_n_months)
# pour chaque mois futur on peut créer une pile de prédictions ; ici on prédit sur le dernier mois
pred_timeID <- max(future_timeIDs)

A_pred <- inla.spde.make.A(mesh = mesh,
                           loc = st_coordinates(pp_sf),
                           group = pred_timeID,
                           n.group = n_time)

stk_pred <- inla.stack(
  data = list(y = NA),
  A = list(A_pred, 1),
  effects = list(spatial = spde_index,
                 data.frame(Intercept = 1,
                            timeID = rep(pred_timeID, nrow(pp_sf)),
                            population = rep(mean_pop, nrow(pp_sf)),
                            ndvi_s = rep(mean_ndvi, nrow(pp_sf)),
                            humidite_s = rep(mean_hum, nrow(pp_sf)),
                            precip_s = rep(mean_prec, nrow(pp_sf)),
                            temp_s = rep(mean_temp, nrow(pp_sf)))),
  tag = "prediction"
)

stk_full <- inla.stack(stk_obs, stk_pred)

res_pred <- inla(formula,
                 family = "poisson",
                 data = inla.stack.data(stk_full),
                 control.predictor = list(A = inla.stack.A(stk_full), compute = TRUE),
                 control.compute = control_compute,
                 offset = c(offset_obs, rep(log(mean_pop/1000 + 1e-9), nrow(pp_sf))),
                 verbose = TRUE)

# --- Extraire prédictions ---
index_pred <- inla.stack.index(stk_full, tag = "prediction")$data
pred_mean_counts <- res_pred$summary.fitted.values[index_pred, "mean"]
pred_incidence_per1000 <- pred_mean_counts * 1000 / mean_pop
pp_sf$pred_incidence_per1000 <- pred_incidence_per1000

pp_wgs84 <- st_transform(pp_sf, crs = 4326)

# --- Carte ---
ggplot() +
  geom_sf(data = st_transform(df_sf, 4326), aes(color = cases), size = 0.6) +
  geom_sf(data = pp_wgs84, aes(fill = pred_incidence_per1000), shape = 15, size = 2) +
  scale_fill_viridis_c() +
  labs(title = paste("Prédiction incidence/1000 pour timeID =", pred_timeID)) +
  theme_minimal()
