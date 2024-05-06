## %######################################################%##
#                                                          #
####                   3-Pre-modeling                   ####
#                                                          #
## %######################################################%##
{
  require(dplyr)
  require(terra)
  require(flexsdm)
  require(here)
  require(progress)
}


## %######################################################%##
#                                                          #
####             Create directory structure             ####
#                                                          #
## %######################################################%##

dir <- flexsdm::sdm_directory(
  main_dir = file.path(getwd(), "3-Models"),
  projections = list.dirs("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/Predictors/Future_CFP", recursive = FALSE) %>% basename(),
  calibration_area = TRUE,
  algorithm = c("gam", "gau", "gbm", "glm", "max", "net", "raf", "svm", "meanthr"),
  ensemble = NULL,
  threshold = FALSE,
  return_vector = TRUE
)
# dir[1] %>% fs::dir_tree(., recurse = TRUE)
dir %>% head()

## %######################################################%##
#                                                          #
####                  Calibration area                  ####
#                                                          #
## %######################################################%##

# Unfiltered occurrence database
unfilt_occ <- data.table::fread(here("2-AllOccurrences/2_all_cleaned_data_final.gz")) %>% tibble()
unfilt_occ <- unfilt_occ %>% select(id = IDr, species = search_name, x = longitude_m, y = latitude_m, year)

# Filtered occurrence database (used for modeling)
occ <- readr::read_tsv(here("2-AllOccurrences/2_all_cleaned_data_final_filtered.gz"))
occ <- occ %>% select(id = IDr, species = search_name, x = longitude_m, y = latitude_m)

# Jepson ecoregion

cfp <- terra::vect("D:/83-NSF_spatial_and_species_traits/3-Variables/JepsonEcoregion/JepsonRegions.shp")
cfp$Region <- NULL

# Process species
sp <- unfilt_occ$species %>%
  unique() %>%
  sort()

for (i in 1:length(sp)) {
  x2 <-
    flexsdm::calib_area(
      data = unfilt_occ[unfilt_occ$species == sp[i], ],
      x = "x",
      y = "y",
      method = c("mask", cfp, "RegionCode")
    )
  x2$Region <- NULL
  terra::writeVector(x2, file.path(grep("Calibration", dir, value = TRUE), paste0(sp[i], ".shp")), overwrite = TRUE)
}
rm(unfilt_occ)

## %######################################################%##
#                                                          #
####                 Data partitioning                  ####
#                                                          #
## %######################################################%##
# Raterize CFP database

pred <-
  file.path(getwd(), "3-Models/1_Inputs/2_Predictors/1_Current") %>%
  list.files(., patter = ".tif$", full.names = TRUE) %>%
  terra::rast()
# pred <- homogenize_na(pred)

cfp_r <- rasterize(cfp, pred$aet, field = "RegionCode")
cfp_r[] <- as.numeric(cfp_r[])

# Extract ecoregion value
occ$region <- terra::extract(cfp_r, occ[, c("x", "y")])[, "RegionCode"]
occ <- occ %>% filter(!is.na(region))
occ$pr_ab <- 1

# Block partition for species with > 30 presences
nr <- occ %>% count(species)
sp <- nr$species[nr$n > 30]

pred <- pred[[!names(pred) %in% "category"]]
names(pred)
part <- list()
for (i in 1:length(sp)) {
  print(i)
  part[[i]] <-
    flexsdm::part_sblock(
      env_layer = pred,
      data = occ %>% filter(species == sp[i]),
      x = "x",
      y = "y",
      pr_ab = "pr_ab",
      n_part = 4, # four blocks
      min_res_mult = 50,
      max_res_mult = 300,
      num_grids = 30,
      min_occ = 10,
      prop = 0.5
    )
}
names(part) <- sp

# Species without good partition
sp <- names(which(sapply(part, function(x) is.null(names(x)))))
sp # Partition were found for all species

part <- part[sapply(part, length) == 3]

# saver block partition details and rasters
best_part_info <- bind_rows(lapply(part, function(x) x$best_part_info), .id = "species")
occ_part <- bind_rows(lapply(part, function(x) x$part), .id = "species")
occ_part %>%
  group_by(species, .part) %>%
  count() %>%
  filter(n == min(n))

# data.table::fwrite(best_part_info, file.path(dir[3], "best_part_info_block.gz"))
# data.table::fwrite(occ_part, file.path(dir[3], "occ_part_block.gz"))


blocks <- lapply(part, function(x) {
  get_block(pred$aet, x$grid)
})
blocks <- rast(blocks)
names(blocks) <- names(part)
# writeRaster(blocks, file.path(dir[3], paste0('part_',names(blocks), ".tif")))


## %######################################################%##
#                                                          #
####    Create pseudo_absences and background points    ####
####      for species with geographical partition       ####
#                                                          #
## %######################################################%##
# Partition layer
part <- list.files(here("3-Models/1_Inputs/1_Occurrences"), pattern = ".tif", full.names = TRUE)
names(part) <- gsub(".tif$", "", gsub("part_", "", basename(part)))
# Calibration areas
clibarea <- list.files(here("3-Models/1_Inputs/3_Calibration_area"), pattern = ".shp$", full.names = TRUE)
names(clibarea) <- gsub(".shp$", "", gsub("part_", "", basename(clibarea)))
# Occurrences
occ_part <- list.files(here("3-Models/1_Inputs/1_Occurrences"), patter = "occ_part", full.names = TRUE)
occ_part <- lapply(occ_part, function(x) data.table::fread(x) %>% tibble())
occ_part <- bind_rows(occ_part)
occ_part$species %>% unique()

sp <- names(part)
db_bg <- db_pa <- list() # object to store ps-abs and bgp


for (i in 1:length(sp)) {
  print(i)
  coord <- occ_part %>%
    dplyr::filter(species == sp[i]) %>%
    dplyr::select(x, y, .part)
  r <- terra::rast(part[sp[i]])
  v <- terra::vect(clibarea[[sp[i]]])
  r <- r %>%
    crop(., v) %>%
    mask(., v)
  part_n <- coord$.part %>%
    table() %>%
    c()
  abs <- lapply(names(part_n), function(x) {
    flexsdm::sample_pseudoabs(
      data = coord,
      x = "x",
      y = "y",
      n = part_n[x],
      method = c("geo_const", width = "5000"),
      rlayer = r,
      maskval = as.numeric(x)
    )
  })
  names(abs) <- names(part_n)
  abs <- bind_rows(abs, .id = ".part") %>% dplyr::mutate(.part = as.numeric(.part))

  db_pa[[i]] <- abs %>% select(x, y, pr_ab, .part)

  abs <- lapply(names(part_n), function(x) {
    flexsdm::sample_background(
      data = coord,
      x = "x",
      y = "y",
      n = round(10000 / length(part_n)),
      method = "random",
      rlayer = r,
      maskval = as.numeric(x)
    )
  })
  names(abs) <- names(part_n)
  abs <- bind_rows(abs, .id = ".part") %>% dplyr::mutate(.part = as.numeric(.part))

  db_bg[[i]] <- abs %>% select(x, y, pr_ab, .part)
}
names(db_pa) <- names(db_bg) <- sp
bind_rows(db_pa, .id = "species")
db_bg <- bind_rows(db_bg, .id = "species")
db_pa <- bind_rows(db_pa, .id = "species")

occ_part <- bind_rows(occ_part, db_pa) %>% arrange(species, desc(pr_ab))

# data.table::fwrite(occ_part, file.path("3-Models/1_Inputs/1_Occurrences", "1_occ_presabs_geopart.gz"))
# data.table::fwrite(db_bg, file.path("3-Models/1_Inputs/1_Occurrences", "1_occ_bkground_geopart.gz"))


## %######################################################%##
#                                                          #
####                   Merge all data                   ####
#                                                          #
## %######################################################%##

# IT IS NOT NECESSARY

# occ_presabs <- bind_rows(
#   data.table::fread(
#     file.path(
#       "3-Models/1_Inputs/1_Occurrences",
#       "1_occ_presabs_geopart.gz"
#     )
#   ) %>% tibble(),
#   data.table::fread(
#     file.path(
#       "3-Models/1_Inputs/1_Occurrences",
#       "1_occ_presabs_randompart.gz"
#     )
#   ) %>% tibble()
# )
#
# occ_bkground <- bind_rows(
#   data.table::fread(
#     file.path(
#       "3-Models/1_Inputs/1_Occurrences",
#       "1_occ_bkground_geopart.gz"
#     )
#   ) %>% tibble(),
#   data.table::fread(
#     file.path(
#       "3-Models/1_Inputs/1_Occurrences",
#       "1_occ_bkground_randompart.gz"
#     )
#   ) %>% tibble()
# )

# data.table::fwrite(occ_presabs, file.path("3-Models/1_Inputs/1_Occurrences", "2_occ_presabs_part_allsp.gz"))
# data.table::fwrite(occ_bkground, file.path("3-Models/1_Inputs/1_Occurrences", "2_occ_bkground_part_allsp.gz"))
