##%######################################################%##
#                                                          #
####       Occurrence data cleaning and filtering       ####
#                                                          #
##%######################################################%##

# Written by: Santiago S.E. Velazco

# Packages and functions
{
  require(terra)
  require(dplyr)
  require(data.table)
  require(ggplot2)
  require(ggspatial)
  require(patchwork)
  require(tidyr)
  require(here)
}

# Useful terra and polygons
# Study area
studya <- terra::rast("aet.tif")

env_variables <-
  './bcm_v65_current_env/' %>%
  list.files(.,patter='.tif$', full.names = TRUE) %>%
  terra::rast()
# env_variables <- flexsdm::homogenize_na(env_variables)

env_variables$aet %>% plot

# CFP
cfp <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/JepsonEcoregion/JepsonRegions.shp")

# Databases:
db0 <- data.table::fread('./2-AllOccurrences/0_all_raw_data.gz') %>% tibble
unique(db0$search_name) %>% sort

# Salix and Prosopis species will be removed
db0 <- db0 %>% dplyr::filter(!search_name%in%c("Salix nigra", "Prosopis pubescens"))
unique(db0$search_name) %>% sort

dim(db0)
##%######################################################%##
#                                                          #
####                     Filtering                      ####
#                                                          #
##%######################################################%##

# Sort columns
ncell <-
  terra::cellFromXY(studya,
                    db0 %>% dplyr::select(longitude_m, latitude_m) %>% as.matrix)
db0 <- tibble(db0, ncell)

db0 <- db0 %>% dplyr::arrange(search_name, ncell, desc(year), data_base)
db0 %>% dplyr::filter(search_name=='Asclepias erosa') %>%  dplyr::select(search_name, ncell, year)

##%######################################################%##
#                                                          #
####                  Filter by cell                    ####
#                                                          #
##%######################################################%##
db0 <- db0 %>% group_by(search_name) %>% dplyr::filter(!duplicated(ncell))
db0 <- db0 %>% dplyr::select(-ncell)
db0 <- db0 %>% group_by()
dim(db0)

##%######################################################%##
#                                                          #
####        Filter occurrence outside study area        ####
#                                                          #
##%######################################################%##

filt <- raster::extract(studya,
                        db0 %>% dplyr::select(longitude_m, latitude_m) %>% data.frame)
filt <- filt$aet
db0 <- db0 %>% dplyr::filter(!is.na(filt))
dim(db0)

##%######################################################%##
#                                                          #
####             Filtering by institution               ####
#                                                          #
##%######################################################%##
cal_inst <- data.table::fread("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/AuxiliaryData/CalInstitutions.txt") %>% tibble()


ins_buf <- terra::vect(cal_inst, geom=c('x_m', 'y_m'), crs=crs(studya))
ins_buf <- terra::buffer(ins_buf, 270*2)
ins_buf$v <- 1
ins_buf <- terra::aggregate(ins_buf, dissolve=TRUE, vars='v')

filt <- terra::extract(ins_buf, terra::vect(db0, geom=c('longitude_m', 'latitude_m')))
filt <- data.frame(filt)
unique(filt$id.x)
filt <- is.na(filt$id.x)
sum(!filt)
db0 <- db0[filt,]
dim(db0)

##%######################################################%##
#                                                          #
####                  Filter by year                    ####
#                                                          #
##%######################################################%##
table(is.na(db0$year))
dbNA <- db0 %>% dplyr::filter(is.na(year))
dbNA %>% count(search_name) %>% arrange(desc(n))

db0 <- db0 %>% dplyr::filter(year>=1950)
db0 %>% count(search_name) %>% arrange(desc(n))
plot(studya)
points(db0$longitude_m, db0$latitude_m)

##%######################################################%##
#                                                          #
####              Filtering by occurrence               ####
####            geographical precision (for             ####
####         this filter NA will no be removed)         ####
#                                                          #
##%######################################################%##
db0$location_quality
db0 <- db0 %>% dplyr::filter(location_quality %in% c("medium", "high", ""))


##%######################################################%##
#                                                          #
####     detecting potential not wild occurrences       ####
#                                                          #
##%######################################################%##

non_wild <-
  c('botanic',
    'botanical',
    'zoo',
    'campus',
    'cultivated',
    'nursery',
    'garden',
    'campus')

filt <- grepl(paste(non_wild,collapse="|"),
              tolower(db0$locality))
table(filt) #38 occurrences have this words

db0 <- db0 %>% dplyr::mutate(wild=!filt)

# individualize each row
db0 <- db0 %>% group_by(data_base) %>% mutate(IDr=paste(data_base, 1:length(data_base)))
db0 <- db0 %>% group_by()


plot(cfp)
points(db0 %>%
         dplyr::select(longitude_m, latitude_m), col='red', pch=19)


##%######################################################%##
#                                                          #
####                  Detect outliers                   ####
#                                                          #
##%######################################################%##
require(flexsdm)
# pseudo absences database
cfp %>% plot
spp <- db0$search_name %>% unique
absences <- as.list(spp)
names(absences) <- db0$search_name %>% unique

xy <- db0 %>% dplyr::select(longitude_m, latitude_m) %>% data.frame %>% vect(geom=c("longitude_m", "latitude_m"))

filt <- raster::extract(cfp, xy)[, "RegionCode"]
head(filt)
db0 <- mutate(db0, ecofilt=filt)
absences_eco <- db0 %>% dplyr::select(search_name, ecofilt) %>% distinct() %>% na.omit()
absences_nr <- db0 %>% count(search_name)
i=1

studya
for (i in 1:length(absences)) {
  print(i)
  f <-
    absences_eco %>% dplyr::filter(search_name == names(absences[i])) %>% pull(ecofilt)
  cfp_r2 <- terra::mask(studya, cfp[cfp$RegionCode %in% f, ])
  pres <- db0 %>% dplyr::filter(search_name== names(absences[i]))
  set.seed(1)
  absences[[names(absences[i])]] <- flexsdm::sample_pseudoabs(
    data = pres,
    x = 'longitude_m',
    y = 'latitude_m',
    n = 10000,
    method = "random",
    rlayer = cfp_r2
  )
}
absences <- lapply(absences, data.frame)
absences <- dplyr::bind_rows(absences, .id='search_name')
absences <- absences %>% mutate(IDr=as.character(1:nrow(absences)))

# data.table::fwrite(absences, './Pseudo_absences_for_outliers.gz')

db2 <- db0 %>% select(search_name, longitude_m, latitude_m, IDr )
db2$pr_ab <- 1
db2 <- bind_rows(db2, absences)
# save(db2, file = "C:/Users/santi/Documents/GitHub/spatial_sp_traits/Data/spp_pres_psabs.RData")

require(flexsdm)
db2 <- split(db2, db2$search_name)

i <- 1
out <- list()
for(i in 1:length(db2)){
  print(i)
  out[[i]] <- env_outliers(da = db2[[i]],
                           x = 'longitude_m',
                           y = 'latitude_m',
                           pr_ab = 'pr_ab',
                           env_layer = env_variables[[!is.factor(env_variables)]],
                           id = 'IDr')
}

out <- bind_rows(out)
db0 <- bind_rows(db0)
names(db0)
names(out)
db0 <-
  left_join(db0, out, by=c('IDr', 'search_name'))
nrow(db0)


db0$latitude_m.y <- NULL
db0$longitude_m.y <- NULL
names(db0)[colnames(db0)=='latitude_m.x'] <- 'latitude_m'
names(db0)[colnames(db0)=='longitude_m.x'] <- 'longitude_m'

# data.table::fwrite(db0, './2-AllOccurrences/1_all_cleaned_data.gz')


##%######################################################%##
#                                                          #
####              Some exploratory figures              ####
#                                                          #
##%######################################################%##
dir.create('./Figures')
db0 <- data.table::fread('./2-AllOccurrences/1_all_cleaned_data.gz') %>% tibble

# number of records by species
rspecies <- db0 %>% group_by(search_name, data_base) %>% count %>% arrange(1)

f1 <- ggplot(db0, aes(year)) + geom_histogram() + theme_minimal()

f2 <- ggplot(rspecies, aes(data_base, sort(n))) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  coord_flip() +
  labs(y='n records', x=element_blank(), fill=element_blank()) +
  theme(legend.position = 'bottom')

f3 <- ggplot(rspecies, aes(search_name, sort(n))) +
  geom_bar(stat = 'identity', aes(fill=data_base)) +
  theme_classic() +
  coord_flip() +
  labs(y='n records', x=element_blank(), fill=element_blank()) +
  theme(legend.position = 'bottom')

plot(cfp)
f4 <- ggplot() +
  ggspatial::layer_spatial(cfp, col='black') +
  geom_hex(data = db0, aes(longitude_m, latitude_m), bins=50) +
  scale_fill_gradientn(colours = pals::jet(15)) +
  theme_minimal()

ggsave(
  plot = f3 |
    f1 / f2,
  filename = "./Figures/Number_of_records.png",
  scale = 1.2,
  width = 18,
  height = 18,
  units = 'cm',
  dpi=500
)

ggsave(
  plot = f4,
  filename = "./Figures/Records_density.png",
  scale = 1.2,
  width = 18,
  height = 18,
  units = 'cm',
  dpi=500
)

##%######################################################%##
#                                                          #
####                 individual figures                 ####
#                                                          #
##%######################################################%##
require(rasterVis)
sp <- unique(db0$search_name)

cfp_r <- terra::rasterize(cfp, studya, field="RegionCode")
plot(cfp_r)

for(i in 1:length(sp)){
  db <- db0 %>% dplyr::filter(search_name==sp[i]) %>% group_by()

  xy <- db %>% dplyr::select(longitude_m, latitude_m ) %>% apply(., 2, range)

  # base <- ggplot() +
  #   ggspatial::layer_spatial(cfp2, col='black', fill='gray90') +
  #   coord_sf(xlim = xy[,1], ylim =xy[,2])

  base <- terra::crop(cfp_r, ext(xy[,1], xy[,2]))
  base0 <-
  gplot(cfp_r) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradientn(colours=pals::parula(6), na.value = 'transparent') +
    # scale_fill_gradient(low = 'red', high = 'gray90', na.value = 'transparent') +
    coord_equal()
  base <- gplot(base) +
    geom_tile(aes(fill = value), show.legend = FALSE) +
    scale_fill_gradient(low = 'gray90', high = 'gray90', na.value = 'transparent') +
    coord_equal()

  a <- base0 +
    # ggplot() +
    # ggspatial::layer_spatial(cfp2, col='black', fill='gray90') +
    geom_point(data = db, aes(longitude_m, latitude_m), col='black', alpha=0.5) +
    theme_minimal() + theme(axis.title = element_blank(), legend.position = 'none')

  b <- base +
    geom_point(data = db, aes(longitude_m, latitude_m, col=year)) +
    scale_color_gradientn(colours = pals::parula(15)[1:12]) +
    theme_minimal() + theme(legend.position = 'bottom') + labs(col=element_blank()) +
    theme(axis.title = element_blank())+labs(subtitle = 'Year')

  c <- base +
    geom_point(data = db, aes(longitude_m, latitude_m, col=data_base)) +
    theme_minimal() + theme(legend.position = 'bottom') + labs(col=element_blank()) +
    theme(axis.title = element_blank())+labs(subtitle = 'Source') +
    guides(col = guide_legend(nrow = 2))

  d <- base +
    geom_point(data = db, aes(longitude_m, latitude_m, col=factor(.out_sum))) +
    theme_minimal() + theme(legend.position = 'bottom') + labs(col=element_blank()) +
    theme(axis.title = element_blank())+labs(subtitle = 'Outliers') +
    guides(col = guide_legend(nrow = 2))

  ggsave(
    plot = ((a + b + c + d)) +
      plot_annotation(title = paste(sp[i], ':', nrow(db))) +
      plot_layout(nrow = 1),
    filename = paste0("./Figures/points ", sp[i], '.png'),
    scale = 1,
    width = 30,
    height = 17,
    units = 'cm',
    dpi=500
  )
}


##### Save shapefiles
# All cleaned records
names(db0)
coord <-
coord <- vect(db0 %>% dplyr::select(-date), geom=c("longitude_m", "latitude_m"), crs=crs(studya))
plot(coord)
terra::writeVector(coord, "./2-AllOccurrences/1_all_cleaned_data.shp", overwrite=TRUE)


# Cleaned records for each species
for(i in 1:length(sp)) {
  terra::writeVector(coord[coord$search_name == sp[i], ],
                    paste0("./2-AllOccurrences/", sp[i], '.shp'),
                    overwrite = TRUE)
}


source("C:/Users/santi/Dropbox/R functions/Scraping_Jepson_Herbarium.R")
s <- c("Asclepias erosa","Bursera microphylla","Pinus californiarum","Washingtonia filifera")
s2 <- jepson_url(s)
s3 <- jepson_info(s2)
# readr::write_tsv(s3, './2-AllOccurrences/0_species_info.txt')




##%######################################################%##
#                                                          #
####    Database and Figure for cleaned occurrences     ####
#                                                          #
##%######################################################%##
require(terra)
require(dplyr)
require(data.table)
require(ggplot2)
require(patchwork)
require(tidyr)
require(here)

# After the check of made by Janet I removed those records she pointed out
###### Process manually cleaned database

# Occurrence database
db <- data.table::fread('./2-AllOccurrences/1_all_cleaned_data.gz') %>% tibble
# Occurrence database from shapefile
coord <- terra::vect("./2-AllOccurrences/1_all_cleaned_data_manually_revised.shp")
dim(db)
dim(coord)
db <- db %>% dplyr::filter(IDr%in%coord$IDr)

db %>% group_by(search_name) %>% count %>% arrange(n) %>% data.frame()
# data.table::fwrite(db,'./2-AllOccurrences/2_all_cleaned_data_final.gz')



##%######################################################%##
#                                                          #
####            number of records by species            ####
#                                                          #
##%######################################################%##
require(rasterVis)
db0 <- data.table::fread('./2-AllOccurrences/2_all_cleaned_data_final.gz') %>% tibble()
dim(db0)

studya <- terra::rast("aet.tif")
cfp_r <- terra::vect("C:/Users/santi/OneDrive/Documentos/FORESTAL/1-Trabajos/83-NSF_spatial_and_species_traits/3-Variables/JepsonEcoregion/JepsonRegions.shp") %>%
  terra::rasterize(., studya, field="RegionCode")
plot(cfp_r)

base <-
  gplot(cfp_r) +
  geom_tile(aes(fill = value)) +
  scale_fill_gradientn(colours=pals::parula(6), na.value = 'transparent') +
  # scale_fill_gradient(low = 'red', high = 'gray90', na.value = 'transparent') +
  coord_equal()

rspecies <- db0 %>% group_by(search_name, data_base) %>% count %>% arrange(1)

f1 <- ggplot(db0, aes(year)) + geom_histogram() + theme_minimal()

f2 <- ggplot(rspecies, aes(data_base, sort(n))) +
  geom_bar(stat = 'identity') +
  theme_classic() +
  coord_flip() +
  labs(y='n records', x=element_blank(), fill=element_blank()) +
  theme(legend.position = 'bottom')

f3 <- ggplot(rspecies, aes(search_name, sort(n))) +
  geom_bar(stat = 'identity', aes(fill=data_base)) +
  theme_classic() +
  coord_flip() +
  labs(y='n records', x=element_blank(), fill=element_blank()) +
  theme(legend.position = 'bottom')


f4 <- base +
  geom_point(data = db0, aes(longitude_m, latitude_m)) +
  # scale_fill_gradientn(colours = pals::jet(15)) +
  theme_minimal()
f4
dir.create("./FiguresCleanedOcc/")
ggsave(
  plot = f3 |
    f1 / f2,
  filename = "./FiguresCleanedOcc/Number_of_records.png",
  scale = 1.2,
  width = 18,
  height = 18,
  units = 'cm',
  dpi=500
)

ggsave(
  plot = f4,
  filename = "./FiguresCleanedOcc/Records_density.png",
  scale = 1.2,
  width = 18,
  height = 18,
  units = 'cm',
  dpi=500
)

sp <- db0$search_name %>% sort %>% unique
for(i in 1:length(sp)){
  db <- db0 %>% dplyr::filter(search_name==sp[i]) %>% group_by()

  xy <- db %>% dplyr::select(longitude_m, latitude_m ) %>%
    apply(., 2, range)


  # base <- ggplot() +
  #   ggspatial::layer_spatial(cfp2, col='black', fill='gray90') +
  #   coord_sf(xlim = xy[,1], ylim =xy[,2])
  base <- terra::crop(studya, ext(xy[,1], xy[,2]))
  base <- gplot(base) +
    geom_tile(aes(fill = value)) +
    scale_fill_gradient(low = 'gray90', high = 'gray90', na.value = 'transparent') +
    coord_equal()

  c <- base +
    geom_point(data = db, aes(longitude_m, latitude_m), col='blue') +
    theme_minimal() + theme(legend.position = 'bottom') + labs(col=element_blank()) +
    theme(axis.title = element_blank())+
    guides(col = guide_legend(nrow = 2))

  ggsave(
    plot = c +
      plot_annotation(title = paste(sp[i], ':', nrow(db))),
    filename = file.path("./FiguresCleanedOcc", paste0("points ", sp[i], '.png')),
    scale = 1.2,
    width = 16,
    height = 16,
    units = 'cm',
    dpi=400
  )
}


# TODO CONTINU HERE
##%######################################################%##
#                                                          #
####           Correct spatial sampling bias            ####
#                                                          #
##%######################################################%##
{
  require(flexsdm)
  require(terra)
  require(dplyr)
  require(data.table)
  require(ggplot2)
  require(ggspatial)
  require(patchwork)
  require(tidyr)
  require(here)
}


db0 <-
  data.table::fread('./2-AllOccurrences/2_all_cleaned_data_final.gz') %>%
  tibble

env_variables <-
  'bcm_v65_current_env/' %>%
  list.files(.,patter='.tif$', full.names = TRUE) %>%
  terra::rast()
env_variables$category <- NULL
names(env_variables)
env_variables$terrain <- NULL # remove categorical variable
env_variables <- flexsdm::homogenize_na(env_variables)

sp_table <-
  db0 %>%
  group_by(search_name) %>%
  count %>%
  arrange(n) %>%
  data.frame()
sp_table <- sp_table %>% mutate(need_bias_corr = FALSE)
sp_table[sp_table$n > 100, 3] <- TRUE

sp <- sp_table %>% dplyr::filter(need_bias_corr==T) %>% pull(search_name)

bin_2 <- list()
for(i in 1:length(sp)){
  bin_2[[i]] <- flexsdm::occfilt_env(
    data = db0 %>% dplyr::filter(search_name == sp[i]),
    x = "longitude_m",
    y = "latitude_m",
    id = "IDr",
    env_layer = env_variables,
    nbins = 2
  )
}
names(bin_2) <- sp

# data.table::fwrite(bind_rows(bin_2, .id="search_name"), './2-AllOccurrences/2_all_cleaned_data_final_sp_filtered2bin.gz')

bin_3 <- list()
for(i in 1:length(sp)) {
  message(paste(sp[i], i))
  bin_3[[i]] <- occfilt_env(
    data = db0 %>% dplyr::filter(search_name == sp[i]),
    x = "longitude_m",
    y = "latitude_m",
    id = "IDr",
    env_layer = env_variables,
    nbins = 3
  )
}

names(bin_3) <- sp
# data.table::fwrite(bind_rows(bin_3, .id="search_name"), './2-AllOccurrences/2_all_cleaned_data_final_sp_filtered3bin.gz')

bin_4 <-list()
for(i in 1:length(sp)) {
  message(paste(sp[i], i))
  try(bin_4[[i]] <- occfilt_env(
    data = db0 %>% dplyr::filter(search_name == sp[i]),
    x = "longitude_m",
    y = "latitude_m",
    id = "IDr",
    env_layer = env_variables,
    nbins = 4
  ))
}

names(bin_4) <- sp
# data.table::fwrite(bind_rows(bin_4, .id="search_name"), './2-AllOccurrences/2_all_cleaned_data_final_sp_filtered4bin.gz')


bin_6 <-list()
for(i in 1:length(sp)) {
  message(paste(sp[i], i))
  try(bin_6[[i]] <- occfilt_env(
    data = db0 %>% dplyr::filter(search_name == sp[i]),
    x = "longitude_m",
    y = "latitude_m",
    id = "IDr",
    env_layer = env_variables,
    nbins = 6
  ))
}

# names(bin_6) <- sp
# data.table::fwrite(bind_rows(bin_6, .id="search_name"), './2-AllOccurrences/2_all_cleaned_data_final_sp_filtered6bin.gz')

bin_8 <-list()
for(i in 1:length(sp)) {
  message(paste(sp[i], i))
  try(bin_8[[i]] <- occfilt_env(
    data = db0 %>% dplyr::filter(search_name == sp[i]),
    x = "longitude_m",
    y = "latitude_m",
    id = "IDr",
    env_layer = env_variables,
    nbins = 8
  ))
}

# names(bin_8) <- sp
# data.table::fwrite(bind_rows(bin_8, .id="search_name"), './2-AllOccurrences/2_all_cleaned_data_final_sp_filtered8bin.gz')

##%######################################################%##
#                                                          #
####             Autocorrelation analysis               ####
#                                                          #
##%######################################################%##
require(ape)

# Environmental variables
pred <-
  './bcm_v65_current_env/' %>%
  list.files(.,patter='.tif$', full.names = TRUE) %>%
  terra::rast()

pred$category <- NULL # remove factors
pred <- flexsdm::homogenize_na(pred)


# Autocorrelation with 2 bins
bin2 <- data.table::fread('./2-AllOccurrences/2_all_cleaned_data_final_sp_filtered2bin.gz') %>% tibble
bin2.1 <- bin2 %>% group_by(search_name) %>%
  group_split()

imoran <- list()
for(i in 1:length(bin2.1)){
  print(i)
  coord <- bin2.1[[i]] %>% dplyr::select(longitude_m, latitude_m)
  data <- data.frame(terra::extract(pred, coord))[-1]
  distm <- dist(coord)
  distm <- as.matrix(distm)
  distm <- 1/distm
  diag(distm) <- 0
  try(imoran[[i]] <-
        apply(data, 2, function(x)
          Moran.I(x, distm, na.rm = T)[c(1, 4)] %>% unlist) %>% data.frame() %>% as_tibble()
  )
  try(imoran[[i]] <- imoran[[i]][1,])
}

names(imoran) <- unique(bin2$search_name)
imorandf <- dplyr::bind_rows(imoran, .id='search_name') %>% mutate(nbins='bin2') %>% dplyr::relocate(nbins)
imorandf <- imorandf %>% mutate(mean_bin=apply(imorandf[, names(pred)], 1, mean))
imorandf$nrecords <- sapply(bin2.1, nrow)
# readr::write_tsv(imorandf,here('2-AllOccurrences/2_filtered2bin_peformance.txt'))


# Autocorrelation with 3 bins
bin2 <- data.table::fread(here('2-AllOccurrences/2_all_cleaned_data_final_sp_filtered3bin.gz')) %>% tibble
bin2.1 <- bin2 %>% group_by(search_name) %>%
  group_split()

length(bin2.1)

imoran <- list()
for(i in 1:length(bin2.1)){
  print(i)
  coord <- bin2.1[[i]] %>% dplyr::select(longitude_m, latitude_m)
  data <- data.frame(terra::extract(pred, coord))[-1]
  distm <- dist(coord)
  distm <- as.matrix(distm)
  distm <- 1/distm
  diag(distm) <- 0
  try(imoran[[i]] <-
        apply(data, 2, function(x)
          Moran.I(x, distm, na.rm = T)[c(1, 4)] %>% unlist) %>% data.frame() %>% as_tibble()
  )
  try(imoran[[i]] <- imoran[[i]][1,])
}

names(imoran) <- unique(bin2$search_name)
imorandf <- dplyr::bind_rows(imoran, .id='search_name') %>% mutate(nbins='bin3') %>% dplyr::relocate(nbins)
imorandf <- imorandf %>% mutate(mean_bin=apply(imorandf[, names(pred)], 1, mean))
imorandf$nrecords <- sapply(bin2.1, nrow)
# readr::write_tsv(imorandf,here('2-AllOccurrences/2_filtered3bin_peformance.txt'))


# Autocorrelation with 4 bins
bin2 <- data.table::fread(here('2-AllOccurrences/2_all_cleaned_data_final_sp_filtered4bin.gz')) %>% tibble
bin2.1 <- bin2 %>% group_by(search_name) %>%
  group_split()

length(bin2.1)

imoran <- list()
for(i in 1:length(bin2.1)){
  print(i)
  coord <- bin2.1[[i]] %>% dplyr::select(longitude_m, latitude_m)
  data <- data.frame(terra::extract(pred, coord))[-1]
  distm <- dist(coord)
  distm <- as.matrix(distm)
  distm <- 1/distm
  diag(distm) <- 0
  try(imoran[[i]] <-
        apply(data, 2, function(x)
          Moran.I(x, distm, na.rm = T)[c(1, 4)] %>% unlist) %>% data.frame() %>% as_tibble()
  )
  try(imoran[[i]] <- imoran[[i]][1,])
}

names(imoran) <- unique(bin2$search_name)
imorandf <- dplyr::bind_rows(imoran, .id='search_name') %>% mutate(nbins='bin3') %>% dplyr::relocate(nbins)
imorandf <- imorandf %>% mutate(mean_bin=apply(imorandf[, names(pred)], 1, mean))
imorandf$nrecords <- sapply(bin2.1, nrow)
# readr::write_tsv(imorandf,here('2-AllOccurrences/2_filtered4bin_peformance.txt'))


# Autocorrelation with 6 bins
bin2 <- data.table::fread(here('2-AllOccurrences/2_all_cleaned_data_final_sp_filtered6bin.gz')) %>% tibble
bin2.1 <- bin2 %>% group_by(search_name) %>%
  group_split()

length(bin2.1)

imoran <- list()
for(i in 1:length(bin2.1)){
  print(i)
  coord <- bin2.1[[i]] %>% dplyr::select(longitude_m, latitude_m)
  data <- data.frame(terra::extract(pred, coord))[-1]
  distm <- dist(coord)
  distm <- as.matrix(distm)
  distm <- 1/distm
  diag(distm) <- 0
  try(imoran[[i]] <-
        apply(data, 2, function(x)
          Moran.I(x, distm, na.rm = T)[c(1, 4)] %>% unlist) %>% data.frame() %>% as_tibble()
  )
  try(imoran[[i]] <- imoran[[i]][1,])
}

names(imoran) <- unique(bin2$search_name)
imorandf <- dplyr::bind_rows(imoran, .id='search_name') %>% mutate(nbins='bin3') %>% dplyr::relocate(nbins)
imorandf <- imorandf %>% mutate(mean_bin=apply(imorandf[, names(pred)], 1, mean))
imorandf$nrecords <- sapply(bin2.1, nrow)
# readr::write_tsv(imorandf,here('2-AllOccurrences/2_filtered6bin_peformance.txt'))

# Autocorrelation with 8 bins
bin2 <- data.table::fread(here('2-AllOccurrences/2_all_cleaned_data_final_sp_filtered8bin.gz')) %>% tibble
bin2.1 <- bin2 %>% group_by(search_name) %>%
  group_split()

length(bin2.1)

imoran <- list()
for(i in 1:length(bin2.1)){
  print(i)
  coord <- bin2.1[[i]] %>% dplyr::select(longitude_m, latitude_m)
  data <- data.frame(terra::extract(pred, coord))[-1]
  distm <- dist(coord)
  distm <- as.matrix(distm)
  distm <- 1/distm
  diag(distm) <- 0
  try(imoran[[i]] <-
        apply(data, 2, function(x)
          Moran.I(x, distm, na.rm = T)[c(1, 4)] %>% unlist) %>% data.frame() %>% as_tibble()
  )
  try(imoran[[i]] <- imoran[[i]][1,])
}

names(imoran) <- unique(bin2$search_name)
imorandf <- dplyr::bind_rows(imoran, .id='search_name') %>% mutate(nbins='bin3') %>% dplyr::relocate(nbins)
imorandf <- imorandf %>% mutate(mean_bin=apply(imorandf[, names(pred)], 1, mean))
imorandf$nrecords <- sapply(bin2.1, nrow)
# readr::write_tsv(imorandf,here('2-AllOccurrences/2_filtered8bin_peformance.txt'))




bindsperformance <- here('2-AllOccurrences/') %>%
  list.files(pattern = 'bin_peformance', full.names = TRUE) %>%
  lapply(., readr::read_tsv)
names(bindsperformance) <- paste0(c(2:4, 6, 8), "bin")

bindsperformance <- bind_rows(bindsperformance, .id='nbins')

dim(bindsperformance)
bindsperformance <- bindsperformance %>% arrange(search_name, nbins) %>% select(nbins, search_name, mean_bin, nrecords)


selected <- bindsperformance %>%
  group_by(search_name) %>%
  filter(mean_bin <= mean(mean_bin)) %>%
  filter(nrecords == max(nrecords)) %>%
  group_by()

selected$selected <- TRUE
bindsperformance <- left_join(bindsperformance, selected)
bindsperformance$selected[is.na(bindsperformance$selected)] <- FALSE
# readr::write_tsv(bindsperformance, here('2-AllOccurrences/2_bins_selected.txt'))



##%######################################################%##
#                                                          #
####        Final database with filtered species        ####
#                                                          #
##%######################################################%##
binsperf <- readr::read_tsv(here('2-AllOccurrences/2_bins_selected.txt'))  %>%
  filter(selected) %>% arrange(nbins)

selected # Only will be used 3 and 4 bind

# Read occurrences databases
db0 <- data.table::fread('./2-AllOccurrences/2_all_cleaned_data_final.gz') %>%
  tibble # entire database


# Filter bin3
sp <- binsperf %>% filter(nbins=='3bin') %>% pull(search_name)
db_bin3 <- data.table::fread('./2-AllOccurrences/2_all_cleaned_data_final_sp_filtered3bin.gz') %>%
  tibble %>% filter(search_name %in% sp)

# Filter bin4
sp <- binsperf %>% filter(nbins=='4bin') %>% pull(search_name)
db_bin4 <- data.table::fread('./2-AllOccurrences/2_all_cleaned_data_final_sp_filtered4bin.gz') %>%
  tibble %>% filter(search_name %in% sp)

# # Remove filtered species from the unfiltered databases
# db0 <- db0 %>% filter(!search_name%in%binsperf$search_name)

# Merge all databases
# db0 <- bind_rows(db0, db_bin3, db_bin4)
db0 <- bind_rows(db_bin3, db_bin4)
db0$search_name %>% unique %>% length()


# readr::write_tsv(db0, here('2-AllOccurrences/2_all_cleaned_data_final_filtered.gz'))

db_raw <- data.table::fread('./2-AllOccurrences/2_all_cleaned_data_final.gz') %>%
  tibble # entire database

db_filt <- readr::read_tsv(here('2-AllOccurrences/2_all_cleaned_data_final_filtered.gz'))

db_raw %>% group_by(search_name) %>% count()
db_filt %>% group_by(search_name) %>% count()
