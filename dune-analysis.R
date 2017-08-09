library(dplyr)
library(data.table)
library(mgcv)
library(nlme)
library(ggplot2)
#library(foreign)

source("dune-analysis_funs.R")

# dune_data <- read.dbf("B:/0_scratchProcessing/Phantom3/desert_dunes/transects_dunes_points_dtm_slope.dbf", as.is=T)
# save(dune_data, file="dune_data.RData")
load("dune_data.RData")

# shrub_data <- read.dbf("transect_shapefiles/dunes_shrub_dens.dbf", as.is = T)
# save(shrub_data, file = "shrub_data.RData")
load("shrub_data.RData")

# add transect length/height -----------------------------------------------------

dune_data_seq <- dune_data %>%
  filter(h_asl > 0) %>%
  filter(!dune %in% c("dune5", "dune7a")) %>%
  group_by(Id) %>% #transect number
  mutate(distance = seq(0, by=0.1, length.out = n()),
         height = sapply(h_asl, function(x) x - min(h_asl)),
         height_slope = height * slope)



# calculate metrics -------------------------------------------------------

dune_data_metrics <- dune_data_seq %>%
  group_by(dune, fence, Id) %>%
  summarise(dune_area = sum(height) * 0.1,
            height_mean = mean(height),
            height_max = max(height),
            height_95 = quantile(height, probs = 0.95),
            slope_mean = mean(slope),
            slope_med = median(slope),
            slope_max = max(slope),
            slope_50 = quantile(slope, probs = 0.5),
            slope_95 = quantile(slope, probs = 0.95),
            slopeheight = mean(height_slope),
            slopeheight_max = max(height_slope),
            compactness = dune_area / height_max,
            east = mean(east), north = mean(north)) %>%
  mutate(fence_dune = paste0(fence,"_",dune))


shrub_metrics <- shrub_data %>%
  mutate(mulga = !is.na(Species)) %>%
  group_by(fence_dune) %>%
  summarise(count = n(),
            transect_area = mean(area),
            shrub_dens = count/transect_area,
            mulga_dens = sum(mulga)/transect_area)


dune_shrub_metrics <- inner_join(dune_data_metrics, shrub_metrics, by = "fence_dune")

shrub_metrics$fence <- unlist(lapply(strsplit(shrub_metrics$fence_dune, "_"), '[[', 1))
shrub_metrics$dune <- unlist(lapply(strsplit(shrub_metrics$fence_dune, "_"), '[[', 2))

# dune transect plots -----------------------------------------------------

ggplot(data = dune_data_seq, aes(x = distance, y = height, group = Id)) +
  geom_line(aes(colour = fence), size = 0.1, alpha = 0.5) +
  facet_wrap(~dune, scales = "free")

ggplot(data = dune_data_metrics, aes(x = Id, y = slopeheight)) +
  geom_point(aes(colour = fence), size = 2, alpha = 1) +
  facet_wrap(~dune, scales = "free")



# dune metric plots -------------------------------------------------------

ggplot(data = dune_data_metrics, aes(y = dune_area)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = height_mean)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = height_max)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = slope_95)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = slope_50)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = slope_max)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = slopeheight)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = slopeheight_max)) +
  geom_boxplot(aes(x = dune, fill = fence))

ggplot(data = dune_data_metrics, aes(y = compactness)) +
  geom_boxplot(aes(x = dune, fill = fence))



# shrub metric plots ------------------------------------------------------

boxplot(shrub_dens ~ fence, data = shrub_metrics)

boxplot(mulga_dens ~ fence, data = shrub_metrics)



# dune and shrub models ---------------------------------------------------

dune_shrub_metrics$dune <- factor(dune_shrub_metrics$dune) # need factors for s(by = ) input
dune_shrub_metrics$fence <- factor(dune_shrub_metrics$fence)
dune_shrub_metrics$fence_dune <- factor(dune_shrub_metrics$fence_dune)

# shrubiness across fence - YES
summary(lm(shrub_dens ~ fence, data = shrub_metrics))
summary(lm(mulga_dens ~ fence, data = shrub_metrics))

# shrubiness across space - YES
spatial_smooth <- gam(formula = shrub_dens ~ te(east,north), data = dune_shrub_metrics)
summary(spatial_smooth); plot(spatial_smooth, residuals = T)



# explore
metric_cors <- rbindlist(lapply(X = names(dune_shrub_metrics)[4:15], FUN = shrub_smooth, dune_shrub_metrics, F, F, T))

## models

# height
shrub_smooth("height_mean", dune_shrub_metrics)
shrub_smooth_by_dune("height_mean", dune_shrub_metrics)
shrub_smooth_dune_re("height_mean", dune_shrub_metrics)

shrub_smooth("height_max", dune_shrub_metrics)
shrub_smooth_by_dune("height_max", dune_shrub_metrics)
shrub_smooth_dune_re("height_max", dune_shrub_metrics)

# slope
shrub_smooth("slope_mean", dune_shrub_metrics)
shrub_smooth_by_dune("slope_mean", dune_shrub_metrics)
shrub_smooth_dune_re("slope_mean", dune_shrub_metrics)

shrub_smooth("slope_max", dune_shrub_metrics)
shrub_smooth_by_dune("slope_max", dune_shrub_metrics)
shrub_smooth_dune_re("slope_max", dune_shrub_metrics)

# shape
shrub_smooth("compactness", dune_shrub_metrics)
shrub_smooth_by_dune("compactness", dune_shrub_metrics)
shrub_smooth_dune_re("compactness", dune_shrub_metrics)

shrub_smooth("slopeheight", dune_shrub_metrics)
shrub_smooth_by_dune("slopeheight", dune_shrub_metrics)
shrub_smooth_dune_re("slopeheight", dune_shrub_metrics)

# area
shrub_smooth("dune_area", dune_shrub_metrics)
shrub_smooth_by_dune("dune_area", dune_shrub_metrics)
shrub_smooth_dune_re("dune_area", dune_shrub_metrics)




