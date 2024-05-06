## %######################################################%##
#                                                          #
####                     4-Modeling                     ####
#                                                          #
## %######################################################%##
# devtools::install_github("sjevelazco/flexsdm")
{
  require(dplyr)
  require(terra)
  require(flexsdm)
  require(here)
  require(progress)
  require(raster)
  require(ggplot2)
}
memory.limit(1000000)


# wmean ensemble
meanw_ens <- function(m, w, thr) {
  m <- terra::weighted.mean(m, w)
  m_thr <- m * (m >= thr)
  m <- c(m, m_thr)
  names(m) <- c("wmean", "max_sorensen")
  return(m)
}

meanthr_ens <- function(m, thr_m, thr) {
  for (r in 1:terra::nlyr(m)) {
    m[[r]][m[[r]] < thr_m[r]] <- 0
  }
  m <- terra::mean(m)
  m_thr <- m * (m >= thr)
  m <- c(m, m_thr)
  names(m) <- c("meanthr", "max_sorensen")
  return(m)
}



## %######################################################%##
#                                                          #
####             Prepare data for modeling              ####
####           species with > 50 occurrences            ####
#                                                          #
## %######################################################%##
# Read occurrences databases
occ <- data.table::fread(here("3-Models/1_Inputs/1_Occurrences", "1_occ_presabs_geopart.gz")) %>% tibble()
bkg <- data.table::fread(here("3-Models/1_Inputs/1_Occurrences", "1_occ_bkground_geopart.gz")) %>% tibble()

# Environmental variables - Current conditions
here()
env <-
  here("3-Models/1_Inputs/2_Predictors/1_Current") %>%
  list.files(patter = ".tif$", full.names = TRUE) %>%
  terra::rast()
env_names <- names(env)


n_occ <- occ %>%
  dplyr::filter(pr_ab == 1) %>%
  pull(species) %>%
  table() %>%
  sort()
sp <- names(n_occ)[n_occ > 50]

env_fut <- here("3-Models/1_Inputs/2_Predictors/2_Projection") %>% list.dirs(recursive = F)
names(env_fut) <- env_fut
env_fut <- as.list(env_fut)

for (i in 1:length(env_fut)) {
  env_fut[[i]] <- env_fut[[i]] %>%
    list.files(., patter = "tif$", full.names = TRUE) %>%
    terra::rast()
}
names(env_fut) <- basename(names(env_fut))

# Extract environmental variables
occ <- sdm_extract(occ, x = "x", y = "y", env_layer = env, filter_na = TRUE)
bkg <- sdm_extract(bkg, x = "x", y = "y", env_layer = env, filter_na = TRUE)

# Count number of presences
n_occ <- occ %>%
  dplyr::filter(pr_ab == 1) %>%
  pull(species) %>%
  table() %>%
  sort()
sp <- names(n_occ)[n_occ > 50]

factor <- "terrain"

# CODES BELOW ARE NOT NECESSARY
# ## %######################################################%##
# #                                                          #
# ####             Prepare  data for modeling             ####
# ####           species with <=50 occurrences            ####
# #                                                          #
# ## %######################################################%##
# # Read occurrences databases
# occ <- data.table::fread(here("3-Models/1_Inputs/1_Occurrences", "2_occ_presabs_part_allsp.gz")) %>% tibble()
# bkg <- data.table::fread(here("3-Models/1_Inputs/1_Occurrences", "2_occ_bkground_part_allsp.gz")) %>% tibble()
#
# # Environmental variables - Current conditions
# topo <- "C:/Users/SVelazco/Documents/2-FireModels/SelectedVariables/Current_CFP/bcm_v65_current_env/topo_hetero.tif" %>% terra::rast()
#
# env <-
#   here("3-Models/1_Inputs/2_Predictors/1_Current/env_var.tif") %>%
#   terra::rast()
# env$terrain <- NULL
# env <- rast(list(env, topo))
# # plot(env)
# env_names <- names(env)
#
# # Future predictors
# env_fut <- here("3-Models/1_Inputs/2_Predictors/2_Projection") %>% list.dirs(recursive = F)
# names(env_fut) <- env_fut
# env_fut <- as.list(env_fut)
#
# for (i in 1:length(env_fut)) {
#   env_fut[[i]] <- env_fut[[i]] %>%
#     list.files(., patter = "env_var.tif$", full.names = TRUE) %>%
#     terra::rast()
#   env_fut[[i]]$terrain <- NULL
#   env_fut[[i]] <- rast(list(env_fut[[i]], topo))
# }
#
# names(env_fut) <- basename(names(env_fut))
#
#
# # Extract env conditions
# occ <- sdm_extract(occ, x = "x", y = "y", env_layer = env, filter_na = TRUE)
# bkg <- sdm_extract(bkg, x = "x", y = "y", env_layer = env, filter_na = TRUE)
#
# # Count number of presences
# n_occ <- occ %>%
#   dplyr::filter(pr_ab == 1) %>%
#   pull(species) %>%
#   table() %>%
#   sort()
# sp <- names(n_occ)[n_occ <= 50]
#
# factor <- NULL
#
#


##%######################################################%##
#                                                          #
####          Transform factor levels names it          ####
####           cause error in some algorithms           ####
#                                                          #
##%######################################################%##
# env <- "./3-Models/1_Inputs/2_Predictors/1_Current/terrain.tif"  %>% terra::rast()
# levels(env) <-
# as.character(as.numeric(as.factor(levels(env)[[1]])))
# env <- flexsdm::homogenize_na(env)
# names(env) <- "terrain"
# terra::writeRaster(env, here("3-Models/1_Inputs/2_Predictors/1_Current/terrain2.tif"))
# file.remove("./3-Models/1_Inputs/2_Predictors/1_Current/terrain.tif")

# # Factor variables - Future conditions
# env <- "./3-Models/1_Inputs/2_Predictors/1_Current/terrain2.tif"  %>% terra::rast()
# env_fut <- here("3-Models/1_Inputs/2_Predictors/2_Projection") %>% list.dirs(recursive = F)
#
# for(i in 1:length(env_fut)) {
#   terra::writeRaster(env, here(env_fut[[i]], "terrain2.tif"))
#   file.remove(here(env_fut[[i]], "terrain.tif"))
# }



## %######################################################%##
#                                                          #
####                 Loop for modeling                  ####
#                                                          #
## %######################################################%##
perf_dir <- here("3-Models/2_Outputs/0_Model_performance")
i=3
TUNE=FALSE
# RUN GAIN species 2 and 3
for (i in 3) {
  message(paste("Modeling sp", i, sp[i]))
  pa <- occ %>% dplyr::filter(species == sp[i])
  b <- bkg %>% dplyr::filter(species == sp[i])

  #### Boosted regression trees ####

    m_gbm <- flexsdm::tune_gbm(
      data = pa,
      response = "pr_ab",
      predictors = env_names[!env_names %in% factor],
      predictors_f = factor,
      partition = ".part",
      grid = expand.grid(
        n.trees = seq(10, 200, 10),
        shrinkage = seq(0.1, 2, 0.1),
        n.minobsinnode = seq(1, 15, 2)
      ),
      thr = "max_sorensen",
      metric = "SORENSEN",
      n_cores = 8 # length(pa$.part %>% unique())
    )

    if (length(m_gbm) > 1) {
      h <- m_gbm$hyper_performance
      p1 <- h %>% ggplot(aes(
        x = n.trees, y = SORENSEN_mean,
        col = as.factor(n.minobsinnode)
      ), group = shrinkage) +
        geom_line() +
        facet_wrap(. ~ shrinkage) +
        theme(legend.position = "bottom") +
        labs(x = "n.minobsinnode")
      ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_gbm_sorensen", ".png")), dpi = 200)

      p1 <- h %>% ggplot(aes(
        x = n.trees, y = AUC_mean,
        col = as.factor(n.minobsinnode)
      ), group = shrinkage) +
        geom_line() +
        facet_wrap(. ~ shrinkage) +
        theme(legend.position = "bottom") +
        labs(x = "n.minobsinnode")
      ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_gbm_auc", ".png")), dpi = 200)

      p1 <- h %>% ggplot(aes(
        x = n.trees, y = TSS_mean,
        col = as.factor(n.minobsinnode)
      ), group = shrinkage) +
        geom_line() +
        facet_wrap(. ~ shrinkage) +
        theme(legend.position = "bottom") +
        labs(x = "n.minobsinnode")
      ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_gbm_tss", ".png")), dpi = 200)

      readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_gbm.txt")))
    }


  #### Maximum entropy ####
  m_max <- flexsdm::tune_max(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    background = b,
    partition = ".part",
    clamp = FALSE,
    grid = expand.grid(
      regmult = seq(0.1, 5, 0.2),
      classes = c("l", "lq", "lqh", "lqhp", "lqhpt")
    ),
    thr = "max_sorensen",
    metric = "SORENSEN",
    n_cores = 8
  )

  if (length(m_max) > 1) {
    h <- m_max$hyper_performance
    p1 <- h %>% ggplot(aes(x = regmult, y = SORENSEN_mean, col = classes)) +
      geom_line()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_max_sorensen", ".png")), dpi = 200)
    p1 <- h %>% ggplot(aes(x = regmult, y = AUC_mean, col = classes)) +
      geom_line()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_max_auc", ".png")), dpi = 200)
    p1 <- h %>% ggplot(aes(x = regmult, y = TSS_mean, col = classes)) +
      geom_line()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_max_tss", ".png")), dpi = 200)
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_max.txt")))
  }


  #### Neural Network ####
  m_net <- tune_net(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    grid = expand.grid(
      size = (2:length(env_names)),
      decay = c(seq(0.01, 1, 0.05), 1, 3, 4, 5, 10)
    ),
    thr = "max_sorensen",
    metric = "SORENSEN",
    n_cores = 8
  )

  if (length(m_net) > 1) {
    h <- m_net$hyper_performance
    p1 <- h %>% ggplot(aes(x = decay, y = SORENSEN_mean, col = as.factor(size))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_net_sorensen", ".png")), dpi = 200)

    p1 <- h %>% ggplot(aes(x = decay, y = AUC_mean, col = as.factor(size))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_net_auc", ".png")), dpi = 200)

    p1 <- h %>% ggplot(aes(x = decay, y = TSS_mean, col = as.factor(size))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_net_tss", ".png")), dpi = 200)

    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_net.txt")))
  }


  #### Random forest ####
  m_raf <- tune_raf(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    grid = expand.grid(mtry = seq(1, length(env_names), 1)),
    thr = "max_sorensen",
    metric = "SORENSEN",
    n_cores = 8
  )

  if (length(m_raf) > 1) {
    h <- m_raf$hyper_performance
    h %>% ggplot(aes(x = mtry, y = SORENSEN_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_raf_sorensen", ".png")), dpi = 200)
    h %>% ggplot(aes(x = mtry, y = AUC_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_raf_auc", ".png")), dpi = 200)
    h %>% ggplot(aes(x = mtry, y = TSS_mean)) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_raf_tss", ".png")), dpi = 200)
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_raf.txt")))
  }


  #### Support Vector Machine ####
  m_svm <- tune_svm(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    grid = expand.grid(
      C = seq(2, 80, 4),
      sigma = c(seq(0.001, 0.2, 0.002))
    ),
    thr = "max_sorensen",
    metric = "SORENSEN",
    n_cores = 8
  )

  if (length(m_svm) > 1) {
    h <- m_svm$hyper_performance
    h %>% ggplot(aes(x = sigma, y = SORENSEN_mean, col = as.factor(C))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_svm_sorensen", ".png")), dpi = 200)
    h %>% ggplot(aes(x = sigma, y = AUC_mean, col = as.factor(C))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_svm_auc", ".png")), dpi = 200)
    h %>% ggplot(aes(x = sigma, y = TSS_mean, col = as.factor(C))) +
      geom_line() +
      theme_classic()
    ggsave(filename = here(perf_dir, paste0(sp[i], " hyp_svm_tss", ".png")), dpi = 200)
    readr::write_tsv(x = h, file = here(perf_dir, paste0(sp[i], " hyp_svm.txt")))
  }

  #### Generalized Additive Model ####
  n_t <- flexsdm:::n_training(data = pa, partition = ".part")

  candidate_k <- 20
  while (any(n_t < flexsdm:::n_coefficients(data = pa, predictors = env_names[!env_names %in% factor], predictors_f = factor, k = candidate_k))) {
    candidate_k <- candidate_k - 3
  }

  m_gam <- fit_gam(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    thr = "max_sorensen",
    k = candidate_k
  )

  # Error in gam.fit3(x = X, y = y, sp = L %*% lsp3 + lsp0, Eb = Eb, UrS = UrS,  :
  #                     inner loop 3; can't correct step size

  #### Generalized Linear Models ####
  if (sum(pa$pr_ab == 1) >= 100) {
    m_glm <- fit_glm(
      data = pa,
      response = "pr_ab",
      # predictors = env_names[c(1:8, 10)],
      predictors = env_names[!env_names %in% factor],
      predictors_f = factor,
      partition = ".part",
      thr = "max_sorensen",
      poly = 2
    )
  }


  #### Gaussian Process ####
  m_gau <- fit_gau(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    # background = b,
    thr = "max_sorensen"
  )


  # compile models objects
  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]

  ##### Ensemble ####
  filt <- flexsdm::sdm_summarize(lapply(models, get))
  filt_perf <- filt$BOYCE_mean >= 0.5 &
    filt$FPB_mean >= 1 &
    filt$SORENSEN_mean >= 0.7 &
    round(filt$thr_value, 2) != 0 &
    round(filt$thr_value, 2) != 1
  models_perf <- models[filt_perf]

  m_ensemble <-
    flexsdm::fit_ensemble(
      lapply(models_perf, get),
      ens_method = c("meanthr"),
      thr_model = "max_sorensen",
      metric = "SORENSEN",
      thr = "max_sorensen"
    )

  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]

  # Model performance
  performance <- flexsdm::sdm_summarize(lapply(models, function(x) {
    if (length(get(x)) > 0) {
      get(x)
    }
  }))

  readr::write_tsv(x = performance, file = here(perf_dir, paste0(sp[i], "_models_performance.txt")))


  ## %######################################################%##
  #                                                          #
  ####             Predict individual models              ####
  #                                                          #
  ## %######################################################%##
  models <- models[models != "m_ensemble"]

  message("Predicting models for species ", sp[i], " ", i)

  models_object <- lapply(models, function(x) {
    get(x)
  })

  prd <-
    sdm_predict(
      models = models_object,
      pred = env,
      thr = c("max_sorensen"),
      con_thr = TRUE,
      clamp = TRUE,
      pred_type = "cloglog",
      predict_area = NULL
    )

  for (mm in 1:length(prd)) {
    terra::writeRaster(prd[[mm]],
      here(
        "3-Models/2_Outputs/1_Current",
        "Algorithm",
        gsub("subsp.", "subsp", names(prd[mm])),
        paste0(sp[i], ".tif")
      ),
      overwrite = TRUE
    )
  }


  ##### Model projection #####
  # f=6
  for (f in 1:length(env_fut)) {
    print(f)
    prd <-
      sdm_predict(
        models = models_object,
        pred = env_fut[[f]],
        thr = c("max_sorensen"),
        con_thr = TRUE,
        clamp = TRUE,
        pred_type = "cloglog",
        predict_area = NULL
      )

    for (mm in 1:length(prd)) {
      terra::writeRaster(prd[[mm]],
        here(
          "3-Models/2_Outputs/2_Projection",
          names(env_fut[f]),
          "Algorithm",
          gsub("subsp.", "subsp", names(prd[mm])),
          paste0(sp[i], ".tif")
        ),
        overwrite = TRUE
      )
    }
    rm(prd)
  }

  rm(list = grep("m_", ls(), value = TRUE))
}




## %######################################################%##
#                                                          #
####               Thresholded ensemble                 ####
#                                                          #
## %######################################################%##
# modeled <- "3-Models/2_Outputs/1_Current/Algorithm/meanw" %>% list.files() %>% gsub(".tif$", "",.)
perf_dir <- here("3-Models/2_Outputs/0_Model_performance")
memory.limit(1000000)
sp <- c("Asclepias erosa", "Bursera microphylla", "Pinus californiarum")
sp="Washingtonia filifera"
i=1
for (i in 1:length(sp)) {
  message(paste("Modeling sp", i, sp[i]))
  pa <- occ %>% dplyr::filter(species == sp[i])
  b <- bkg %>% dplyr::filter(species == sp[i])
  perf <- readr::read_tsv(
    file = here(perf_dir, paste0(sp[i], "_models_performance.txt")),
    col_types = readr::cols()
  )

  #### Boosted regression trees ####
  perf_2 <- perf[perf$model == "gbm", ] %>% data.frame()
  m_gbm <- flexsdm::fit_gbm(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part", n_trees = perf_2$n.trees,
    n_minobsinnode = perf_2$n.minobsinnode,
    shrinkage = perf_2$shrinkage,
    thr = "max_sorensen"
  )


  #### Maximum entropy ####
  perf_2 <- perf[perf$model == "max", ] %>% data.frame()
  m_max <- flexsdm::tune_max(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    background = b,
    partition = ".part",
    grid = expand.grid(
      regmult = perf_2$regmult,
      classes = perf_2$classes
    ),
    thr = "max_sorensen",
    metric = "SORENSEN",
    n_cores = 8
  )


  #### Neural Network ####
  perf_2 <- perf[perf$model == "net", ] %>% data.frame()
  m_net <- fit_net(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    thr = "max_sorensen", size = perf_2$size, decay = perf_2$decay
  )


  #### Random forest ####
  perf_2 <- perf[perf$model == "raf", ] %>% data.frame()
  m_raf <- fit_raf(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    thr = "max_sorensen", mtry = perf_2$mtry
  )


  #### Support Vector Machine ####
  perf_2 <- perf[perf$model == "svm", ] %>% data.frame()
  m_svm <- fit_svm(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    thr = "max_sorensen", sigma = perf_2$sigma, C = perf_2$C
  )

  #### Generalized Additive Model ####
  n_t <- flexsdm:::n_training(data = pa, partition = ".part")

  candidate_k <- 20
  while (any(n_t < flexsdm:::n_coefficients(data = pa, predictors = env_names[!env_names %in% factor], predictors_f = factor, k = candidate_k))) {
    candidate_k <- candidate_k - 3
  }

  m_gam <- fit_gam(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    thr = "max_sorensen",
    k = candidate_k
  )

  # Error in gam.fit3(x = X, y = y, sp = L %*% lsp3 + lsp0, Eb = Eb, UrS = UrS,  :
  #                     inner loop 3; can't correct step size

  #### Generalized Linear Models ####
  if (sum(pa$pr_ab == 1) >= 100) {
    m_glm <- fit_glm(
      data = pa,
      response = "pr_ab",
      # predictors = env_names[c(1:8, 10)],
      predictors = env_names[!env_names %in% factor],
      predictors_f = factor,
      partition = ".part",
      thr = "max_sorensen",
      poly = 2
    )
  }


  #### Gaussian Process ####
  m_gau <- fit_gau(
    data = pa,
    response = "pr_ab",
    predictors = env_names[!env_names %in% factor],
    predictors_f = factor,
    partition = ".part",
    # background = b,
    thr = "max_sorensen"
  )


  # Compile models objects
  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]

  # Filter by performance
  filt <- flexsdm::sdm_summarize(lapply(models, get))
  filt_perf <- filt$BOYCE_mean >= 0.5 &
    filt$FPB_mean >= 1 &
    filt$SORENSEN_mean >= 0.7 &
    round(filt$thr_value, 2) != 0 &
    round(filt$thr_value, 2) != 1
  models_perf <- models[filt_perf]

  ##### Ensemble ####
  m_ensemble <-
    flexsdm::fit_ensemble(
      lapply(models_perf, get),
      ens_method = c("meanthr"),
      thr_model = "max_sorensen",
      metric = "SORENSEN",
      thr = "max_sorensen"
    )

  models <- grep("m_", ls(), value = TRUE)
  filt <- sapply(models, function(x) {
    length(get(x))
  })
  models <- models[filt > 0]

  # Model performance
  performance <- flexsdm::sdm_summarize(lapply(models, function(x) {
    if (length(get(x)) > 0) {
      get(x)
    }
  }))

  readr::write_tsv(x = performance, file = here(perf_dir, paste0(sp[i], "_models_performance_2.txt")))


  ## %######################################################%##
  #                                                          #
  ####             Predict individual models              ####
  #                                                          #
  ## %######################################################%##
  models <- models[models != "m_ensemble"]

  message("Predicting models for species ", sp[i], " ", i)

  models_object <- lapply(models, function(x) {
    get(x)
  })

  prd <-
    sdm_predict(
      models = models_object,
      pred = env,
      thr = c("max_sorensen"),
      con_thr = TRUE,
      clamp = TRUE,
      pred_type = "cloglog",
      predict_area = NULL
    )

  for (mm in 1:length(prd)) {
    terra::writeRaster(prd[[mm]],
      here(
        "3-Models/2_Outputs/1_Current",
        "Algorithm",
        gsub("subsp.", "subsp", names(prd[mm])),
        paste0(sp[i], ".tif")
      ),
      overwrite = TRUE
    )
  }


  ##### Model projection #####
  for (f in 1:length(env_fut)) {
    print(f)
    prd <-
      sdm_predict(
        models = models_object,
        pred = env_fut[[f]],
        thr = c("max_sorensen"),
        con_thr = TRUE,
        clamp = TRUE,
        pred_type = "cloglog",
        predict_area = NULL
      )

    for (mm in 1:length(prd)) {
      terra::writeRaster(prd[[mm]],
        here(
          "3-Models/2_Outputs/2_Projection",
          names(env_fut[f]),
          "Algorithm",
          gsub("subsp.", "subsp", names(prd[mm])),
          paste0(sp[i], ".tif")
        ),
        overwrite = TRUE
      )
    }
    rm(prd)
  }

  rm(list = grep("m_", ls(), value = TRUE))
}




## %######################################################%##
#                                                          #
####                 Predict ensemble                   ####
#                                                          #
## %######################################################%##
# Read occurrences databases
occ <- data.table::fread(here("3-Models/1_Inputs/1_Occurrences/1_occ_presabs_geopart.gz")) %>% tibble()
sp <- occ$species %>% unique()
perf_dir <- here("3-Models/2_Outputs/0_Model_performance")

env_fut <- here("3-Models/1_Inputs/2_Predictors/2_Projection") %>% list.dirs(recursive = F)
names(env_fut) <- env_fut
env_fut <- as.list(env_fut)
names(env_fut) <- basename(names(env_fut))
sp <- c("Asclepias erosa", "Bursera microphylla", "Pinus californiarum")

for (i in 1:length(sp)) {
  print(i)
  ens_perf <- here(perf_dir, paste0(sp[i], "_models_performance.txt")) %>%
    readr::read_tsv()

  # models_perf <- gsub("m_", "", models_perf) %>% sort()


  filt_perf <- ens_perf$BOYCE_mean >= 0.5 &
    ens_perf$FPB_mean >= 1 &
    ens_perf$SORENSEN_mean >= 0.7 &
    round(ens_perf$thr_value, 2) != 0 &
    round(ens_perf$thr_value != 1, 2) &
    !ens_perf$model %in% c("meanthr", "meanw")
  models_perf <- ens_perf$model[filt_perf]



  thr <- ens_perf$thr_value[ens_perf$model == "meanthr"]
  thr_mod <- ens_perf$thr_value[ens_perf$model %in% models_perf]

  dd <- here("3-Models/2_Outputs/1_Current/Algorithm") %>%
    list.dirs(., recursive = FALSE)
  dd <- dd[basename(dd) %in% models_perf]
  dd <-
    lapply(dd, function(x) {
      list.files(
        x,
        pattern =  sp[i],
        recursive = FALSE,
        full.names = TRUE
      )
    }) %>%
    unlist() %>%
    terra::rast()
  dd <- dd[[grep("max_sorensen", names(dd), invert = TRUE)]]
  dd <- dd[[names(dd) %>% sort()]]

  # Predict ensemble
  prd <- meanthr_ens(m = dd, thr_m = thr_mod, thr = thr)

  terra::writeRaster(prd, here("3-Models/2_Outputs/1_Current/", "Algorithm/meanthr", paste0(gsub("subsp.", "subsp", sp[i]), ".tif")), overwrite = TRUE)
  rm(prd)


  ##### Model projection #####
  for (f in 1:length(env_fut)) {
    print(f)

    dd <- file.path("3-Models/2_Outputs/2_Projection", names(env_fut[f]), "Algorithm") %>%
      list.dirs(., recursive = FALSE)
    dd <- dd[basename(dd) %in% models_perf]
    dd <-
      lapply(dd, function(x) {
        list.files(
          x,
          pattern =  sp[i],
          recursive = FALSE,
          full.names = TRUE
        )
      }) %>%
      unlist() %>%
      terra::rast()
    dd <- dd[[grep("max_sorensen", names(dd), invert = TRUE)]]
    dd <- dd[[names(dd) %>% sort()]]
    prd <- meanthr_ens(m = dd, thr_m = thr_mod, thr = thr)

    terra::writeRaster(prd, here(
      "3-Models/2_Outputs/2_Projection",
      names(env_fut[f]), "Algorithm/meanthr", paste0(sp[i], ".tif")
    ), overwrite = TRUE)
    rm(prd)
  }
}


## %######################################################%##
#                                                          #
####            Save final performance table             ####
#                                                          #
## %######################################################%##

thr <- "./3-Models/2_Outputs/0_Model_performance/" %>% list.files(pattern = "_models_performance.txt", full.names = TRUE)
thr2 <- lapply(thr, readr::read_tsv)
names(thr2) <- gsub("_models_performance.txt", "", basename(thr))
thr2 <- bind_rows(thr2, .id = "species")
readr::write_tsv(thr2,  "./3-Models/2_Outputs/0_Model_performance/00_model_performance.txt")
