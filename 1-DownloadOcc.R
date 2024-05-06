##%######################################################%##
#                                                          #
####       California flora occurrences database        ####
#                                                          #
##%######################################################%##
# Written by: Santiago J.E. Velazco
setwd("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/6-PresencesOnlySDMs_Tribal")
# packages and functions
require(here)
require(readxl)
require(dplyr)
require(janitor)
require(stringr)
require(textclean)
require(BIEN)
require(ridigbio)
require(rgbif)
require(terra)
require(bdc)
require(rnaturalearth)

sp <- readxl::read_excel("Tribal_presence_only_sdm.xlsx", 1)
sp <- sp$species

# studya <- terra::rast("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/StudyAreaRaster.tif")
# plot(studya)
#
# wld <- rnaturalearth::ne_coastline()
#
# us <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/US/gadm_US.shp")
# # plot(us, col='red')
#
# cal <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/US/gadm_US_CA.shp")
#
# cfp <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/JepsonEcoregion/JepsonRegions_CFP.shp")
#
# cfp2 <- cfp
# # crs(cfp2) <- "+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs"
# cfp2 <- terra::project(cfp2, crs(studya))
#
# plot(cfp2, add=F)


coord_to_m <- function(database, x, y, from_crs="+proj=longlat +datum=WGS84", to_crs){
  coord <- database %>% dplyr::select(x=dplyr::all_of(x), y=dplyr::all_of(y)) %>% as_tibble()
  coord <- coord %>% terra::vect(., geom = c('x', 'y'), from_crs)
  coord <- terra::project(coord, to_crs)
  coord <- terra::geom(coord)[,c('x', 'y')] %>% data.frame %>% dplyr::rename(longitude_degree=x, latitude_degree=y)
  return(coord)
}


##%######################################################%##
#                                                          #
####                Prepare plant codes                 ####
#                                                          #
##%######################################################%##
# vroom::vroom_write(plant_codes, 'California USDA Plants codes.gz')

##%######################################################%##
#                                                          #
####                       CalFlora                     ####
#                                                          #
##%######################################################%##

# Merge all databasets
require(readr)
occdb <- list.files(("./1-CalFlora/"), full.names = T, pattern = 'txt') #%>% as.list
names(occdb) <- gsub(".txt$", "", gsub("calflora-", "", basename(occdb)))
occdb <- as.list(occdb)
# cc <- readr::cols(
#   ID = col_character(),
#   Taxon = col_character(),
#   Observer = col_character(),
#   Source = col_character(),
#   Date = col_date(format = ""),
#   Phenology = col_character(),
#   InfestedArea = col_character(),
#   GrossArea = col_character(),
#   `Percent Cover` = col_logical(),
#   Distribution = col_character(),
#   `National Ownership` = col_character(),
#   County = col_character(),
#   Latitude = col_double(),
#   Longitude = col_double(),
#   `Location Description` = col_character(),
#   `!Comments` = col_logical()
# )

for(i in 1:length(occdb)){
  occdb[[i]] <- readr::read_tsv(occdb[[i]])
}
sapply(occdb,function(x) class(x$Latitude))
sapply(occdb,function(x) class(x$Longitude))
occdb <- bind_rows(occdb, .id = 'search_name') %>%
  janitor::clean_names()
names(occdb)
occdb$date_12 <- NULL
occdb$county_9 <- NULL
occdb$location_description_16 <- NULL

occdb <- rename(occdb, date_=date_5,  county=county_7, location_description=location_description_8)
names(occdb)
# vroom::vroom_write(occdb, here('1-CalFlora', 'occ_CalFlora_Tribal_2022.gz'))


##%######################################################%##
#                                                          #
####               California_Consortium                ####
#                                                          #
##%######################################################%##
require(dplyr)
require(readr)

occdb <- list.files(here('1-CalCons'), full.names = T, pattern = '.csv')
names(occdb) <- gsub(".csv$", "", basename(occdb))
occdb <- as.list(occdb)
occdb

# Pinus californianum
# p <- readr::read_csv("./1-CalCons/Pinus californiarum.csv")
# names(p)
# p$specificEpithet
# p <- p %>% dplyr::filter(scientificName %>% grepl("californiarum", .))
# p <- p[!is.na(p$decimalLongitude),]
# p <- p[!is.na(p$decimalLatitude),]
# p$scientificName %>% table
# readr::write_csv(p, "./1-CalCons/Pinus californiarum.csv")

# n <- read_csv(occdb[[1]]) %>% names

for(i in 1:length(occdb)) {
  occdb[[i]] <-
    readr::read_csv(occdb[[i]]) %>%
    janitor::clean_names()

    # col_types = readr::cols(
    #   specimen_number = col_character(),
    #   `taxon name` = col_character(),
    #   collector = col_character(),
    #   coll_number_prefix = col_character(),
    #   coll_number = col_character(),
    #   coll_number_suffix = col_character(),
    #   collection_date = col_character(),
    #   late_collection_date = col_character(),
    #   verbatim_date = col_character(),
    #   county = col_character(),
    #   elevation = col_character(),
    #   locality = col_character(),
    #   latitude = col_double(),
    #   longitude = col_double(),
    #   datum = col_character(),
    #   source = col_character(),
    #   error_distance = col_character(),
    #   units = col_character(),
    #   citation = col_character()
    # )

# %>%dplyr::select({{coln}})
}

dim(occdb[[1]])
dim(occdb[[2]])

for(i in 1:length(occdb)){
  occdb[[i]]$date_identified <- occdb[[i]]$date_identified %>% as.numeric()
}
occdb <- bind_rows(occdb, .id = 'search_name')
occdb <- occdb %>% rename(scientific_name = scientific_name)
occdb <- data.table::fread(here('1-CalCons', 'CalCons_Tribal_2022.gz')) %>% tibble()
occdb <- occdb[!is.na(occdb$decimal_latitude)| !is.na(occdb$decimal_longitude), ]
# vroom::vroom_write(occdb, here('1-CalCons', 'CalCons_Tribal_2022.gz'))



##%######################################################%##
#                                                          #
####                        BIEN                        ####
#                                                          #
##%######################################################%##
ext_gen_sp <- function(x){
  textclean::replace_non_ascii(x) %>%
  word(1:2) %>%
    paste(., collapse = ' ')
}

sp <- readxl::read_excel("Tribal_presence_only_sdm.xlsx", 1) %>%
  dplyr::pull(species)

occdb <- as.list(rep(NA, length(sp)))
names(occdb) <- sp

for(i in 1:length(sp)){
  message(paste(i, 'from', length(occdb), sp[i]))
  occdb[[i]] <-
    BIEN::BIEN_occurrence_species(
      species = names(occdb[i]) %>% ext_gen_sp,
      observation.type = TRUE,
      all.taxonomy = TRUE,
      natives.only = TRUE,
      political.boundaries = TRUE,
      collection.info = TRUE
    )

  message('Found: ', nrow(occdb[[i]]))


  if(nrow(occdb[[i]])>0) {
    occdb[[i]] <- occdb[[i]][!duplicated(names(occdb[[i]]))] %>%
      as_tibble() %>%
      dplyr::mutate(search_name=names(occdb[i])) %>% #search_name is a column that indicat the names used for record downloading
      dplyr::relocate(search_name)
  }
}

# Delete all those element of the list occdb without occurrences
filt <- sapply(occdb, function(x) nrow(x))
filt <- unlist(filt)>0
filt <- names(filt[filt])
occdb <- occdb[filt]

# Merge all data.frames allocated within the occdb list
occdb <- bind_rows(occdb)

# clean columns names
occdb <- occdb %>% janitor::clean_names()
dim(occdb)

dir.create('./1-BIEN')
# vroom::vroom_write(occdb, here('1-BIEN', 'occ_bien_Tribal_2022.gz'))


##%######################################################%##
#                                                          #
####                        GBIF                        ####
#                                                          #
##%######################################################%##
sp <- readxl::read_excel("Tribal_presence_only_sdm.xlsx", 1) %>%
  dplyr::pull(species)
occdb <- as.list(rep(NA, length(sp)))
names(occdb) <- sp



for(i in 1:length(occdb)){
  message(paste(i, 'from', length(occdb), sp[i]))
  occdb[[i]] <-
    rgbif::occ_data(
      scientificName = names(occdb[i]),
      country = "US",
      hasCoordinate = TRUE,
      limit = 200000
    )[[2]]

  message('Found: ', nrow(occdb[[i]]))

  if(nrow(occdb[[i]])>0) { #E
    occdb[[i]] <- occdb[[i]] %>%
      dplyr::select(-17) %>%
      as_tibble() %>%
      dplyr::mutate(search_name=names(occdb[i])) %>% #search_name is a column that indicate the names used for record dowloading
      dplyr::relocate(search_name)
  }
}

# Delete all those element of the list occdb without occurrences
filt <- sapply(occdb, function(x) nrow(x))
filt <- unlist(filt)>0
filt <- names(filt[filt])
occdb <- occdb[filt]

# Merge all data.frames allocated within the occdb list
occdb <- bind_rows(occdb)

# clean columns names
occdb <- occdb %>% janitor::clean_names()
dim(occdb)
occdb$search_name %>% table
dir.create('./1-GBIF')
# vroom::vroom_write(occdb, here('1-GBIF', 'occ_GBIF_Tribal_2022.gz'))

##%######################################################%##
#                                                          #
####                      iDigBio                       ####
#                                                          #
##%######################################################%##
sp <- readxl::read_excel("Tribal_presence_only_sdm.xlsx", 1) %>%
  dplyr::pull(species)

occdb <- as.list(rep(NA, length(sp)))
names(occdb) <- sp

for(i in 1:length(occdb)) {
  message(paste(i, 'from', length(occdb), sp[i]))
  occdb[[i]] <-
    ridigbio::idig_search_records(rq = list(scientificname = names(occdb[i])),
                                  limit = 100000) %>% as_tibble()

  message('Found: ', nrow(occdb[[i]]))

  if(nrow(occdb[[i]])>0) { #E
    occdb[[i]] <- occdb[[i]] %>%
      dplyr::mutate(search_name=names(occdb[i])) %>% #search_name is a column that indicat the names used for record dowloading
      dplyr::relocate(search_name)
  }
}

# Delete all those element of the list occdb without occurrences
filt <- sapply(occdb, function(x) nrow(x))
filt <- unlist(filt)>0
filt <- names(filt[filt])
occdb <- occdb[filt]

# Merge all data.frames allocated within the occdb list
occdb <- bind_rows(occdb)

# clean columns names
occdb <- occdb %>% janitor::clean_names()


dir.create('./1-IDIGBIO')
# vroom::vroom_write(occdb, here('1-IDIGBIO', 'occ_IDIGBIO_Tribal_2022.gz'))




##%######################################################%##
#                                                          #
####               Cal Fish and Wildlife                ####
#                                                          #
##%######################################################%##
sp <- readxl::read_excel("Tribal_presence_only_sdm.xlsx", 1) %>%
  dplyr::pull(species)

raw_data <- data.table::fread("./1-CalFishAndWildlife/cfw_final_presences.gz") %>% tibble()
raw_data <- raw_data %>% dplyr::filter(species%in%sp)
raw_data <- raw_data %>% mutate(search_name=species)
# vroom::vroom_write(raw_data,
                   # './1-CalFishAndWildlife/cfw_final_presences_Tribal_2022.gz')



raw_data <- data.table::fread('./1-CalFishAndWildlife/cfw_final_presences_Tribal_2022.gz')

ca_crs <-   ("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

raw_data <- tibble(raw_data, coord_to_m(raw_data, x='x_coords', y='y_coords', from_crs = ca_crs,  to_crs = "+proj=longlat +datum=WGS84")
)

# vroom::vroom_write(raw_data,
# './1-CalFishAndWildlife/cfw_final_presences_Tribal_2022.gz')


##%######################################################%##
#                                                          #
####                        Ryan                        ####
#                                                          #
##%######################################################%##

#### Pinus californiarum
db <- readr::read_csv("./1-Ryan/Pinus californiarum.csv")
db$species <- "Pinus californiarum"
db <- unique(db)

# Code to recover date
db_dbif <- rgbif::occ_data(
  scientificName = "Pinus monophylla",
  country = "US",
  hasCoordinate = TRUE,
  limit = 200000
)[[2]]
db_dbif$year
db_dbif <- rename(db_dbif, latitude=decimalLatitude, longitude=decimalLongitude)

db2 <- left_join(db, db_dbif %>% dplyr::select(longitude, latitude, year))
plot(db2[2:3])
a <- rnaturalearth::ne_countries(continent = "North america", scale = "medium")
plot(a, add=T)


db2$search_name <- "Pinus californiarum"
hist(db2$year)
# vroom::vroom_write(db2, './1-Ryan/Ryan_Tribal_2022.gz')
# db2 <- data.table::fread('./1-Ryan/Ryan_Tribal_2022.gz')
# vroom::vroom_write(db2, './1-Ryan/Pinus californiarum_year.csv')

# #### Pinus monophylla
# db <- readr::read_csv("./1-Ryan/Pinus monophylla.csv")
# db$species <- "Pinus monophylla"
# db <- unique(db)
#
# # Code to recover date
# db_dbif <- rgbif::occ_data(
#   scientificName = "Pinus monophylla",
#   country = "US",
#   hasCoordinate = TRUE,
#   limit = 200000
# )[[2]]
# db_dbif$year
# db_dbif <- rename(db_dbif, latitude=decimalLatitude, longitude=decimalLongitude)
#
# db2 <- left_join(db, db_dbif %>% dplyr::select(longitude, latitude, year))
# plot(db2[2:3])
# db2 <- unique(db2)
#
# db2$search_name <- "Pinus monophylla"
# hist(db2$year)
# # vroom::vroom_write(db2, './1-Ryan/Pinus monophylla_year.csv')
#
# #### Pinus quadrifolia
# db <- readr::read_csv("./1-Ryan/Pinus quadrifolia.csv")
# db$species <- "Pinus quadrifolia"
# db <- unique(db)
#
# # Code to recover date
# db_dbif <- rgbif::occ_data(
#   scientificName = "Pinus quadrifolia",
#   country = "US",
#   hasCoordinate = TRUE,
#   limit = 200000
# )[[2]]
# db_dbif$year
# db_dbif <- rename(db_dbif, latitude=decimalLatitude, longitude=decimalLongitude)
#
# db2 <- left_join(db, db_dbif %>% dplyr::select(longitude, latitude, year))
# plot(db2[2:3])
# db2 <- unique(db2)
#
# db2$search_name <- "Pinus quadrifolia"
#
# # vroom::vroom_write(db2, './1-Ryan/Pinus quadrifolia_year.csv')


#########################################################################################

##%######################################################%##
#                                                          #
####               Integrate all databases               ####
#                                                          #
##%######################################################%##
require(data.table)
require(dplyr)
require(ggplot2)
require(janitor)
require(stringr)
require(terra)
require(ggplot2)

# Species to be modeled
sp <- readxl::read_excel("Tribal_presence_only_sdm.xlsx", 1)
sp$species %>% length

### Basic layers for filtering data -----
# Raster layer of study area
# studya <- terra::rast("aet.tif)
# names(studya) <- 'study_area'


# World
wld <- rnaturalearth::ne_coastline(returnclass = "sp")
crs(wld)

# US WGS84
us <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/US/gadm_US_disolved.shp")
# plot(us, col='red')
us$NAME_1 <- "US"
us <- us[,"NAME_1"]

# California WGS84
studya <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/US/gadm_US_CA.shp")
studya <- studya[,"NAME_1"]
# plot(studya, add=T)







##%######################################################%##
#                                                          #
####         Prepare each dataset separately           ####
#                                                          #
##%######################################################%##

### Read occurrences databases -----
{
  dr <- list.dirs() %>% grep('1-', ., value = T)
  dr <-
    sapply(dr, function(x)
      list.files(x, pattern = '_Tribal_2022.gz', full.names = TRUE))
  raw_data <- list()
  for (i in 1:length(dr)) {
    raw_data[[i]] <- data.table::fread(dr[[i]]) %>% tibble()
  }
  names(raw_data) <- gsub('./1-', '', names(dr))
}


##%######################################################%##
#                                                          #
####                       BIEN                         ####
#                                                          #
##%######################################################%##

names(raw_data)
db <- 'BIEN'
names(raw_data[[db]])

raw_data[[db]] <- raw_data[[db]] %>% mutate(datetime = as.Date(date_collected)) %>%
  dplyr::mutate(year = lubridate::year(datetime))

raw_data[[db]] %>% head
raw_data[[db]] %>% names

raw_data[[db]] <- raw_data[[db]] %>% dplyr::select(
  search_name,
  occ_id = catalog_number,
  scientific_name = scrubbed_species_binomial,
  scientific_family = scrubbed_family,
  longitude,
  latitude,
  locality,
  county,
  state_province,
  country,
  recorded_by,
  basis_of_record=observation_type,
  institution_code=dataset,
  datasource,
  year,
  date=datetime,
) %>%
  mutate(occ_id=as.character(occ_id), date=as.character(date)) %>%
  as_tibble()

# remove records without coordinate
dim(raw_data[[db]])
raw_data[[db]] <- dplyr::filter(raw_data[[db]], !is.na(longitude), !is.na(latitude))




##%######################################################%##
#                                                          #
####                CalFishAndWildlife                  ####
#                                                          #
##%######################################################%##

names(raw_data)
db <- 'CalFishAndWildlife'
names(raw_data[[db]])

raw_data[[db]]$survey_date
raw_data[[db]] <- raw_data[[db]] %>% mutate(datetime = as.Date(survey_date)) %>% dplyr::mutate(
  year = lubridate::year(datetime))
names(raw_data[[db]])
raw_data[[db]] <- raw_data[[db]]

raw_data[[db]] %>% head
raw_data[[db]] %>% names %>% sort

raw_data[[db]] <- raw_data[[db]] %>% dplyr::select(
  search_name,
  occ_id = wypt_id,
  scientific_name = species,
  longitude_m=x_coords,
  latitude_m=y_coords,
  longitude=longitude_degree,
  latitude=latitude_degree,
  basis_of_record=survey_type,
  year,
  date=datetime,
) %>%
  mutate(occ_id=as.character(occ_id), date=as.character(date)) %>%
  as_tibble()

# remove records without coordinate
dim(raw_data[[db]])
raw_data[[db]] <- dplyr::filter(raw_data[[db]], !is.na(longitude), !is.na(latitude))
names(raw_data[[db]])

raw_data[[db]] <- raw_data[[db]] # %>% dplyr::rename(longitude_m=longitude, latitude_m=latitude)



##%######################################################%##
#                                                          #
####                     CalFlora                       ####
#                                                          #
##%######################################################%##

names(raw_data)
db <- 'CalFlora'
names(raw_data[[db]])

raw_data[[db]] <- raw_data[[db]] %>% mutate(datetime = as.Date(date_)) %>% dplyr::mutate(
  year = lubridate::year(datetime))

raw_data[[db]] %>% head
raw_data[[db]] %>% names %>% sort
raw_data[[db]]$plant
names(raw_data[[db]])
raw_data[[db]] <-
  raw_data[[db]] %>% dplyr::select(
  search_name,
  occ_id = id,
  scientific_name = plant,
  longitude,
  latitude,
  datasource=source,
  locality=location_description,
  county,
  recorded_by=observer,
  year,
  date=datetime,
  location_quality
) %>%
  mutate(occ_id=as.character(occ_id), date=as.character(date)) %>%
  as_tibble()

# remove records without coordinate
dim(raw_data[[db]])
raw_data[[db]] <- dplyr::filter(raw_data[[db]], !is.na(longitude), !is.na(latitude))
dim(raw_data[[db]])


##%######################################################%##
#                                                          #
####         Consortium of California herbaria          ####
#                                                          #
##%######################################################%##

names(raw_data)
db <- 'CalCons'
names(raw_data[[db]])

# Rescue year collection data
raw_data[[db]]$event_date
raw_data[[db]]$event_date <- gsub('[.]', '', raw_data[[db]]$event_date)

raw_data[[db]]$year <- lubridate::date(raw_data[[db]]$event_date) %>% lubridate::year()
raw_data[[db]] <- raw_data[[db]][which(raw_data[[db]]$year>=1950|is.na(raw_data[[db]]$year)), ]

# filter coordinates precision

fitl <- (((raw_data[[db]]$coordinate_uncertainty_in_meters)<=13500)) |
  is.na(raw_data[[db]]$coordinate_uncertainty_in_meters)
raw_data[[db]] <- raw_data[[db]][fitl, ]

raw_data[[db]]$associated_collectors
raw_data[[db]] <-
  raw_data[[db]] %>% dplyr::select(
  search_name,
  occ_id = id,
  scientific_name = scientific_name,
  scientific_family =  scientific_name,
  longitude=decimal_longitude,
  latitude=decimal_latitude,
  locality,
  county,
  recorded_by=associated_collectors,
  year,
) %>%
  mutate(occ_id=as.character(occ_id), date=NA) %>%
  as_tibble()

# remove records without coordinate
dim(raw_data[[db]])
raw_data[[db]] <- dplyr::filter(raw_data[[db]], !is.na(longitude), !is.na(latitude))
dim(raw_data[[db]])



##%######################################################%##
#                                                          #
####                       GBIF                         ####
#                                                          #
##%######################################################%##

names(raw_data)
db <- 'GBIF'

raw_data[[db]] <-
  raw_data[[db]] %>% select(-starts_with('http'), -ends_with('_key')) %>%
  mutate(event_date = as.Date(event_date))

names(raw_data[[db]]) %>% sort

raw_data[[db]] <-
  raw_data[[db]] %>% dplyr::select(
  search_name,
  occ_id = key,
  scientific_name = species,
  scientific_family =  family,
  longitude = decimal_longitude,
  latitude = decimal_latitude,
  coord_unc_m = coordinate_uncertainty_in_meters,
  locality,
  municipality,
  county,
  state_province,
  country,
  country_code,
  issues,
  taxon_rank,
  basis_of_record,
  institution_code,
  recorded_by,
  year,
  date=event_date,
) %>%
  mutate(occ_id=as.character(occ_id), date=as.character(date)) %>%
  as_tibble()
raw_data[[db]]$year <- as.numeric(raw_data[[db]]$year)

# remove records without coordinate
dim(raw_data[[db]])
raw_data[[db]] <- dplyr::filter(raw_data[[db]], !is.na(longitude), !is.na(latitude))
dim(raw_data[[db]])

raw_data[[db]] <- raw_data[[db]] %>% mutate(coord_unc_m=ifelse(is.na(coord_unc_m), 0, coord_unc_m))

raw_data[[db]] <- raw_data[[db]] %>% dplyr::filter(coord_unc_m<(270*5))

##%######################################################%##
#                                                          #
####                      IDIGBIO                       ####
#                                                          #
##%######################################################%##

names(raw_data)
db <- 'IDIGBIO'

raw_data[[db]] %>% names
raw_data[[db]] %>% head
raw_data[[db]] <-
  raw_data[[db]] %>% mutate(date = as.Date(datecollected)) %>%
  dplyr::mutate(year = lubridate::year(date))

names(raw_data[[db]]) %>% sort

raw_data[[db]] <-
  raw_data[[db]] %>% dplyr::select(
  search_name,
  occ_id = occurrenceid,
  scientific_name = scientificname,
  scientific_family =  family,
  longitude = geopoint_lon,
  latitude = geopoint_lat,
  state_province=stateprovince,
  country,
  recorded_by=collector,
  year,
  date,
) %>%
  mutate(occ_id=as.character(occ_id), date=as.character(date)) %>%
  as_tibble()

# remove records without coordinate
dim(raw_data[[db]])
raw_data[[db]] <- dplyr::filter(raw_data[[db]], !is.na(longitude), !is.na(latitude))
dim(raw_data[[db]])



##%######################################################%##
#                                                          #
####                 MERGE ALL DATASETS                 ####
#                                                          #
##%######################################################%##
require(ggplot2)

all <- raw_data %>% dplyr::bind_rows(.id='data_base')

# Filter by US and California
filt <- raster::extract(us,
                        vect(dplyr::select(all, longitude, latitude),
                             geom=c('longitude', 'latitude')),
                        df = TRUE,sp = TRUE,  cellnumbers = TRUE)[,-1]

all <- mutate(all, within_us=!is.na(filt))

filt <- raster::extract(studya,
                        vect(dplyr::select(all, longitude, latitude),
                             geom=c('longitude', 'latitude')),
                        df = TRUE,sp = TRUE,  cellnumbers = TRUE)[,-1]
all <- mutate(all, within_ca=!is.na(filt))


all <- all %>% dplyr::filter(within_us)

ggplot(all, aes(longitude, latitude)) + geom_hex(bins=150) + coord_equal() +
  scale_fill_gradientn(colours = pals::jet(10)) + theme_minimal()

all %>% count(data_base) %>%
  ggplot(aes(data_base, n, fill=data_base)) +
  geom_bar(stat = 'identity') +
  coord_flip() +
  theme_minimal() +
  theme()

# Filtering by precision
# 1° = 111 km
# 0.1° = 11.1 km
# 0.01° = 1.11 km
# 0.001° = 111 m
# 0.0001° = 11.1 m
# 0.00001° = 1.11 m
# 0.000001° = 0.11 m



# Filtering by year and coordinate precision
all <- bdc::bdc_coordinates_precision(all, lat = "latitude", lon = "longitude", ndec = 3)
all <- bdc::bdc_year_outOfRange(all, eventDate = "year", year_threshold = 1950)
all <- tibble(all)
all <- all %>% dplyr::filter(.rou, .year_outOfRange)



# Transform coordinates
ca_crs <- ("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")

coord_m <-
  coord_to_m(all, x='longitude', y='latitude', from_crs = "+proj=longlat +datum=WGS84",  to_crs = ca_crs)
names(coord_m) <- c("longitude_m", "latitude_m" )
all$latitude_m <- NULL
all$longitude_m <- NULL
all <- tibble(all, coord_m)
str(all)
all <- bdc_filter_out_flags(all)
data.table::fwrite(all, './2-AllOccurrences/0_all_raw_data.gz')



##%######################################################%##
#                                                          #
####               Save all raw figures                 ####
#                                                          #
##%######################################################%##

dir.create('./2-AllOccurrences/Figures_raw')
us <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/US/gadm_US.shp")
us <- us[,"NAME_1"]
# us2 <- rnaturalearth::ne_states(country = "united states of america") %>% vect()
# us2 <- crop(us2, us)
plot(us)

coord <-  data.table::fread('./2-AllOccurrences/0_all_raw_data.gz') %>% tibble
coord <- coord %>% dplyr::select(search_name, longitude, latitude, data_base)
sp <- coord$search_name %>% unique %>% sort
require(sf)
require(ggspatial)

i=2
for(i in 1:length(sp)){
  xy <- coord %>% dplyr::filter(search_name==sp[i]) %>% dplyr::select(-1)
  xy <- unique(xy)
  xy <- terra::vect(xy, geom=c("longitude", "latitude"))
  us2 <- crop(us, ext(xy))
  us2 <- sf::st_as_sf(us2)
  xy <- coord %>% dplyr::filter(search_name==sp[i]) %>% dplyr::select(-1) %>% unique()

  p <- ggplot(us2) +
    geom_sf() +
    geom_point(data=xy, aes(longitude, latitude, color=data_base)) +
    theme_bw() +
    labs(subtitle = sp[i]) +
    theme(legend.position = "bottom")

  ggsave(filename = file.path('./2-AllOccurrences/Figures_raw', paste0(sp[i], ".png")),
         plot = p, dpi = 300, units = "cm", width = 25, height = 25)
}



