library(dplyr)
library(ggplot2)
library(lubridate)
library(sf)
library(spatstat)
library(spocc)
library(terra)

##################################################
# Helper functions
##################################################

# Return (x, y) which do not give any NA values on the covariates
remove_NA_on_covariates <- function(x, y, covariates) {
  keep <- rep(TRUE, length(x))
  for(i in seq_along(covariates)) {
    keep <- keep & inside.owin(x, y, as.owin(covariates[[i]]))
  }
  list(x = x[keep], y = y[keep])
}

# Predict at locations, given covariates.
# Note: I can't easily use predict.ppm since I had to adjust the intercepts by hand.
predict_at <- function(x, y, year, month, coef, environmental_covariates, no_bias_values) {
  bias_values <- Reduce('+', lapply(names(no_bias_values), function(nm) coef[nm] * no_bias_values[nm]))
  if(is.null(bias_values)) {
    bias_values <- 0
  }
  bias_values <- rep(bias_values, length(x))
  intercept <- rep(coef["(Intercept)"], length(x))
  month_coef <- rep(ifelse(month == "01", 0., coef[paste0("month", month)]), length(x))
  year_coef <- rep(ifelse(year == first_year, 0., coef[paste0("year", year)]), length(x))
  exp(intercept +
        month_coef +
        year_coef +
        Reduce('+', lapply(names(environmental_covariates), function(nm) coef[nm] * as.function(environmental_covariates[[nm]](year = year, month = as.numeric(month)))(x, y))) +
        bias_values)
}

integrate_model <- function(initial_fit,
                            new_fit,
                            month,
                            year,
                            quadrature_points,
                            weight,
                            distance,
                            environmental_covariates,
                            no_bias_values) {
  # Absolute difference between old and new predictions over a set of quadrature points
  pred_diff <- abs(predict_at(x = quadrature_points$x,
                              y = quadrature_points$y,
                              month = month,
                              year = year,
                              coef = initial_fit,
                              environmental_covariates = environmental_covariates,
                              no_bias_values = no_bias_values) -
                     predict_at(x = quadrature_points$x,
                                y = quadrature_points$y,
                                month = month,
                                year = year,
                                coef = new_fit,
                                environmental_covariates = environmental_covariates,
                                no_bias_values = no_bias_values))
  if(distance == "Linf") {
    # In this case, we want the max of weighted predictions
    value <- as.function(weight)(quadrature_points$x, quadrature_points$y) * pred_diff
    max(value)
  } else {
    # Otherwise, we want basically (1 / N \sum_{i = 1}^N w(x_i) * f(x_i)^k)^{1 / k} for quadrature points x_i
    value <- as.function(weight)(quadrature_points$x, quadrature_points$y) * power(pred_diff, distance)
    inverse_power(mean(value), distance)
  }
}

# Power that each loss function corresponds to
power <- function(x, distance) {
  if(distance == "L1") {
    x
  } else if(distance == "L2") {
    x^2
  } else if(distance == "L4") {
    x^4
  }
}

# Inverse of the above powers
inverse_power <- function(x, distance) {
  if(distance == "L1") {
    x
  } else if(distance == "L2") {
    sqrt(x)
  } else if(distance == "L4") {
    x^(1 / 4)
  }
}

# Convert to spatstat im
.raster_to_im <- function(x){
  r <- as.data.frame(x, xy = TRUE)
  im <- spatstat.geom::as.im(r)
  return(im)
}

# Note: to generalise the spatstat.geom::as.im to SpatRaster and RasterLayer
as.im.SpatRaster <- function(x){ .raster_to_im(x) }
as.im.RasterLayer <- function(x){ .raster_to_im(x) }

##################################################
# End of helper functions
##################################################

##################################################
# Format rasters
##################################################

dir.create("data-raw/australia", showWarnings = FALSE)
dir.create("data/australia", showWarnings = FALSE)

states <- c("Victoria")

eastern_states <- rnaturalearth::ne_states(country = "Australia",
                                           returnclass = "sf") %>%
  dplyr::filter(name %in% states) %>%
  sf::st_union() %>%
  terra::vect()

# GDA 2020 Lambert CRS
output_crs <- 7845

elevation <- geodata::elevation_30s(country = "AUS", path = "data-raw/australia")
elevation <- elevation %>%
  terra::crop(eastern_states) %>%
  terra::mask(eastern_states) %>%
  terra::project(paste0("epsg:", output_crs))

writeRaster(elevation, filename = "data/australia/elevation.tif", overwrite = TRUE)

ttime <- geodata::travel_time("city", size = 2, path = "data-raw/australia")
ttime <- ttime %>%
  terra::crop(eastern_states) %>%
  terra::mask(eastern_states) %>%
  terra::project(paste0("epsg:", output_crs))

writeRaster(ttime, filename = "data/australia/travel_time.tif", overwrite = TRUE)

# These were downloaded by hand from https://www.longpaddock.qld.gov.au/silo/gridded-data/ and put in data-raw/australia
raster_files <- c("2015.daily_rain",
                  "2016.daily_rain",
                  "2017.daily_rain",
                  "2018.daily_rain",
                  "2019.daily_rain",
                  "2020.daily_rain",
                  "2021.daily_rain",
                  "2022.daily_rain",
                  "2015.max_temp",
                  "2016.max_temp",
                  "2017.max_temp",
                  "2018.max_temp",
                  "2019.max_temp",
                  "2020.max_temp",
                  "2021.max_temp",
                  "2022.max_temp",
                  "2015.radiation",
                  "2016.radiation",
                  "2017.radiation",
                  "2018.radiation",
                  "2019.radiation",
                  "2020.radiation",
                  "2021.radiation",
                  "2022.radiation")

for(file in raster_files) {
  print(file)
  r <- rast(paste0("data-raw/australia/", file, ".nc"))
  crs(r) <- crs("+init=EPSG:4326")
  mn <- month(time(r))
  rr <- tapp(r, index = mn, fun = "mean")

  rr <- rr %>%
    terra::crop(eastern_states) %>%
    terra::mask(eastern_states) %>%
    terra::project(paste0("epsg:", output_crs))

  rr <- resample(rr, elevation)
  writeRaster(rr, filename = paste0("data/australia/", file, "_monthly.tif"), overwrite = TRUE)
}

##################################################
# End of format rasters
##################################################

# Load data from iNat
occ <- spocc::occ(query = "Phascolarctos cinereus",
                  from = "inat",
                  date = c("2015-01-01", "2022-12-31"),
                  has_coords = TRUE,
                  limit = 1e4)

# Drop unused columns
cleaned_occ <- as.data.frame(occ$inat$data$Phascolarctos_cinereus)[, c("quality_grade", "time_observed_at", "longitude", "latitude")]

# Drop missing values
cleaned_occ <- cleaned_occ[!is.na(cleaned_occ$time_observed_at) &
                             !is.na(cleaned_occ$longitude) &
                             !is.na(cleaned_occ$latitude), ]

# Format the date
cleaned_occ$year_month <- format(as.Date(cleaned_occ$time_observed_at, format = "%Y-%m-%d"),"%Y/%m")
cleaned_occ$time_observed_at <- NULL

# Drop non-research grade
cleaned_occ <- cleaned_occ[cleaned_occ$quality_grade == "research", ]
cleaned_occ$quality_grade <- NULL

# Drop observations not in the target state(s)
cleaned_occ <- sf::st_as_sf(cleaned_occ, agr = "constant", coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_crop(eastern_states) %>%
  as.data.frame()

cleaned_occ$longitude <- sapply(cleaned_occ$geometry, function(p) p[1])
cleaned_occ$latitude <- sapply(cleaned_occ$geometry, function(p) p[2])
cleaned_occ$geometry <- NULL

# Drop duplicated points
cleaned_occ <- cleaned_occ[!duplicated(cleaned_occ), ]

# Load rasters
first_year <- 2015
considered_years <- first_year:2022
considered_months <- 1:12
elevation_im <- as.im(rast(paste0("data/australia/elevation.tif")))
travel_time_im <- as.im(rast(paste0("data/australia/travel_time.tif")))
rainfall_ims <- lapply(considered_years, function(y) lapply(considered_months, function(m) as.im(rast(paste0("data/australia/", y, ".daily_rain_monthly.tif"))[[m]])))
max_temp_ims <- lapply(considered_years, function(y) lapply(considered_months, function(m) as.im(rast(paste0("data/australia/", y, ".max_temp_monthly.tif"))[[m]])))
radiation_ims <- lapply(considered_years, function(y) lapply(considered_months, function(m) as.im(rast(paste0("data/australia/", y, ".radiation_monthly.tif"))[[m]])))

covariates <- list(elevation_im)
covariates <- append(covariates, list(travel_time_im))
covariates <- append(covariates, unlist(rainfall_ims, recursive = FALSE))
covariates <- append(covariates, unlist(max_temp_ims, recursive = FALSE))
covariates <- append(covariates, unlist(radiation_ims, recursive = FALSE))

names(covariates) <- 1:length(covariates)

elevation <- function(year = 2015, month = 1) elevation_im
travel_time <- function(year = 2015, month = 1) travel_time_im
rainfall <- function(year = 2015, month = 1) rainfall_ims[[which(considered_years == year)[1]]][[month]]
max_temp <- function(year = 2015, month = 1) max_temp_ims[[which(considered_years == year)[1]]][[month]]
radiation <- function(year = 2015, month = 1) radiation_ims[[which(considered_years == year)[1]]][[month]]

xy <- as.data.frame(cleaned_occ)

occ_transformed <- sf::st_as_sf(xy, coords = c("longitude", "latitude"), crs = 4326) %>%
  sf::st_transform(crs = output_crs) %>%
  sf::st_coordinates() %>%
  as.data.frame()

occ_transformed$year_month <- cleaned_occ$year_month

window <- spatstat.geom::as.owin(elevation_im)
plot(window)

year_months <- sort(unique(occ_transformed$year_month))
configurations <- lapply(year_months, function(y) ppp(x = occ_transformed$X[occ_transformed$year_month == y],
                                                      y = occ_transformed$Y[occ_transformed$year_month == y],
                                                      window = window))
covariate_functions <- list(elevation = elevation,
                            rainfall = rainfall,
                            max_temp = max_temp)#,
#radiation = radiation)

ones <- elevation_im
ones$v[] <- 1

####################
cost_covariate = travel_time()
observer_bias_covariates = list(travel_time = travel_time)
environmental_covariates = covariate_functions
no_bias_values = c(travel_time = 0.)
weight = ones
distance = "L1"
n_quad = 1000L
n_predict = 400L
detection_probability = 1
W = window
only_benefit = FALSE
trnd = Y ~ month + year + elevation + rainfall + max_temp + travel_time

prediction_points <- as.ppp(gridcenters(window, nx = sqrt(n_predict), ny = sqrt(n_predict)), W = window)
which_nonna <- remove_NA_on_covariates(x = prediction_points$x,
                                       y = prediction_points$y,
                                       covariates = covariates)
xs = which_nonna$x
ys = which_nonna$y
months = rep("02", length(xs))
years = rep(2022, length(xs))

# Check basic preconditions
stopifnot(length(xs) == length(ys),
          length(xs) > 0, n_quad > 0,
          length(observer_bias_covariates) == length(no_bias_values),
          detection_probability >= 0 & detection_probability <= 1)

# Construct quadrature points, and remove points which give NA values
quadrature_points <- as.ppp(gridcenters(window, nx = sqrt(n_quad), ny = sqrt(n_quad)), W = window)
which_nonna <- remove_NA_on_covariates(x = quadrature_points$x,
                                       y = quadrature_points$y,
                                       covariates = covariates)
quadrature_points <- ppp(x = which_nonna$x,
                         y = which_nonna$y,
                         window = window)

# Construct data.frame for spatstat fit
H <- hyperframe(Y = configurations)
for(i in seq_len(length(covariate_functions))) {
  cov <- covariate_functions[[i]]
  nm <- names(covariate_functions)[i]
  H[, nm] <- lapply(seq_len(length(year_months)), function(k) {
    y <- as.numeric(gsub("^([0-9]+)/([0-9]+)$", "\\1", year_months[k]))
    m <- as.numeric(gsub("^([0-9]+)/([0-9]+)$", "\\2", year_months[k]))
    cov(year = y, month = m)
  })
}

for(i in seq_len(length(observer_bias_covariates))) {
  cov <- observer_bias_covariates[[i]]
  nm <- names(observer_bias_covariates)[i]
  H[, nm] <- H[, nm] <- lapply(seq_len(length(year_months)), function(k) {
    y <- as.numeric(gsub("^([0-9]+)/([0-9]+)$", "\\1", year_months[k]))
    m <- as.numeric(gsub("^([0-9]+)/([0-9]+)$", "\\2", year_months[k]))
    observer_bias_covariates[[i]](year = y, month = m)
  })
}

H$month <- gsub("^([0-9]+)/([0-9]+)$", "\\2", year_months)
H$year <- gsub("^([0-9]+)/([0-9]+)$", "\\1", year_months)

# Initial fit
initial_fit <- mppm(trnd, data = H)
initial_fit <- coef(initial_fit)

# Initial average predicted intensity
month <- months[1]
year <- years[1]
p <- predict_at(x = quadrature_points$x,
                y = quadrature_points$y,
                month = month,
                year = year,
                coef = initial_fit,
                environmental_covariates = environmental_covariates,
                no_bias_values = no_bias_values)
pred_initial <- mean(p)

# Construct additional fits when points in (x, y) are added.
fits <- lapply(seq_len(length(xs)), function(i) {
  # Index of latest month/year of the target month
  mo <- months[i]
  ye <- years[i]
  indices <- which(mo == H$month & ye == H$year)
  last_occurrence <- indices[length(indices)]

  # Configuration with an additional point
  new_conf <- ppp(x = c(xs[i], H$Y[[last_occurrence]]$x), y = c(ys[i], H$Y[[last_occurrence]]$y), window = window)
  H_plus <- H
  H_plus$Y[[last_occurrence]] <- new_conf

  # Fit the enlarged configuration
  fit <- mppm(trnd, data = H_plus)
  coefs <- coef(fit)

  # New average predicted intensity
  pred <- mean(predict_at(x = quadrature_points$x,
                          y = quadrature_points$y,
                          month = mo,
                          year = ye,
                          coef = coefs,
                          environmental_covariates = environmental_covariates,
                          no_bias_values = no_bias_values))

  # Adjust the intercept so that both this and the initial one give the same average value
  coefs["(Intercept)"] <- coefs["(Intercept)"] - log(pred) + log(pred_initial)

  coefs
})

# Construct the benefit at the points in (x, y)
result <- sapply(seq_len(length(xs)), function(i) {
  r <- as.function(cost_covariate)(xs[i], ys[i])
  if(!only_benefit || distance == "constant") {
    # Compute the predicted detection probability at (x, y)
    detection_probability <- predict_at(xs[i], ys[i], coef = initial_fit,
                                        month = months[i],
                                        year = years[i],
                                        environmental_covariates = environmental_covariates,
                                        no_bias_values = no_bias_values)
    r <- r * detection_probability
  }
  if(r > 0. && distance != "constant") {
    r <- r * integrate_model(initial_fit = initial_fit,
                             new_fit = fits[[i]],
                             month = months[i],
                             year = years[i],
                             quadrature_points = quadrature_points,
                             weight = weight,
                             distance = distance,
                             environmental_covariates = environmental_covariates,
                             no_bias_values = no_bias_values)
  }
  r
})

# Not strictly necessary, but it makes sense to divive by the length of the window
result <- list(result / volume(window))

make_into_factor <- function(df) {
  dat <- df
  bins <- seq(from = 0, to = 1, length.out = 11)
  bin_names <- paste0(bins[-1] - 0.1, " - ", bins[-1])
  for(i in seq_len(ncol(dat))) {
    if(!colnames(dat)[i] %in% c('X', 'Y')) {
      max_i <- max(dat[, i])
      dat[, i] <- cut(dat[, i] / max_i, breaks = bins, labels = bin_names)
    }
  }
  dat
}

dat <- data.frame(X = xs,
                  Y = ys)

for(i in seq_len(length(distance))) {
  name <- paste0(distance[i], "_no_weight")
  dat <- cbind(dat, result[[i]])
  colnames(dat)[ncol(dat)] <- name
}

library(tidyterra)

eastern_states <- rnaturalearth::ne_states(country = "Australia",
                                           returnclass = "sf") %>%
  dplyr::filter(name %in% states) %>%
  sf::st_union() %>%
  terra::vect() %>%
  terra::project(paste0("epsg:", output_crs))

theme_set(theme_minimal(base_size = 23))
plots <- lapply(3:ncol(dat), function(i) {
  benefit <- make_into_factor(dat)[, i]

  pred_df <- data.frame(x = xs,
                        y = ys,
                        utility = benefit)
  #log_pred = pmax(-2, log2(benefit) - max(log2(benefit))))
  pred_df <- sf::st_as_sf(pred_df, coords = c("x", "y"), crs = output_crs)

  m <- months[1]
  y <- years[1]
  indices <- which(m == H$month & y == H$year)
  last_occurrence <- indices[length(indices)]

  presence_df <- data.frame(x = H$Y[[last_occurrence]]$x,
                            y = H$Y[[last_occurrence]]$y)
  presence_df <- sf::st_as_sf(presence_df, coords = c("x", "y"), crs = output_crs)

  if(grepl("no_weight", colnames(dat)[i], fixed = TRUE)) {
    g <- ggplot(eastern_states) +
      geom_spatvector(alpha = 0.2)
  } else {
    g <- rasterVis::gplot(raster(weight)) +
      geom_tile(aes(fill = value), alpha = 0.4) +
      # Below corresponds to `palette = "Purple-Yellow" with c2 = 0
      colorspace::scale_fill_continuous_sequential(name = "Weight",
                                                   h1 = 320,
                                                   h2 = 80,
                                                   c1 = 60,
                                                   c2 = 0,
                                                   l1 = 30,
                                                   l2 = 95,
                                                   p1 = 0.7,
                                                   p2 = 1.3,
                                                   cmax = 65)
  }
  g <- g +
    geom_sf(data = pred_df, aes(colour = utility), size = 1.5, alpha = 0.5, inherit.aes = FALSE) +
    geom_sf(data = presence_df, colour = 'red', size = 0.8, alpha = 0.8, inherit.aes = FALSE) +
    #scale_color_manual(values = scales::seq_gradient_pal("blue", "red", "Lab")(seq(0, 1, length.out = 10))) +
    scale_color_viridis_d(name = "Relative utility", option = "A", direction = -1) +
    #ggtitle(paste0("Loss function ", colnames(dat)[i])) +
    xlab("") +
    ylab(element_blank()) +
    theme(strip.background = element_blank(),
          strip.text.x = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    guides(colour = guide_legend(override.aes = list(size = 5, alpha = 1)))
  g
})

plots[[1]]

