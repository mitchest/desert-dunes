library(dplyr)
library(data.table)
library(mgcv)
library(nlme)
library(ggplot2)
library(itsadug)
#library(moments)
#library(foreign)

source("dune-analysis_funs.R")

#dune_data <- read.dbf("transects_dunes_points_dtm_slope_dens.dbf", as.is=T)
#save(dune_data, file="dune_data.RData")
load("dune_data.RData")



# add transect length/height -----------------------------------------------------

dune_data_seq <- dune_data %>%
  filter(h_asl > 0) %>%
  filter(!dune %in% c("dune5", "dune7a", "dune2")) %>%
  rename(shrub_dens = RASTERVALU) %>%
  group_by(Id) %>% #transect number
  mutate(distance = seq(0, by=0.1, length.out = n()),
         height = sapply(h_asl, function(x) x - min(h_asl)),
         r_j = height - mean(height),
         q_j = slope - mean(slope)) %>%
  ungroup() %>%
  filter(height > 0.2,
         !(distance > 175 & dune == "dune2"))

# dune 9/7 transects were switched somehow
#dune_data_seq$distance[dune_data_seq$dune == "dune7a" & dune_data_seq$fence == "inside"] <- abs(dune_data_seq$distance[dune_data_seq$dune == "dune7a" & dune_data_seq$fence == "inside"] - max(dune_data_seq$distance[dune_data_seq$dune == "dune7a" & dune_data_seq$fence == "inside"]))
dune_data_seq$distance[dune_data_seq$dune == "dune7b" & dune_data_seq$fence == "outside"] <- abs(dune_data_seq$distance[dune_data_seq$dune == "dune7b" & dune_data_seq$fence == "outside"] - max(dune_data_seq$distance[dune_data_seq$dune == "dune7b" & dune_data_seq$fence == "outside"]))
dune_data_seq$distance[dune_data_seq$dune == "dune9" & dune_data_seq$fence == "outside"] <- abs(dune_data_seq$distance[dune_data_seq$dune == "dune9" & dune_data_seq$fence == "outside"] - max(dune_data_seq$distance[dune_data_seq$dune == "dune9" & dune_data_seq$fence == "outside"]))


# calculate metrics -------------------------------------------------------

dune_shrub_metrics <- dune_data_seq %>%
  group_by(dune, fence, Id) %>%
  summarise(shrub_dens = mean(shrub_dens),
            dune_area = sum(height) * 0.1,
            dune_width = n() * 0.1,
            height_max = max(height),
            height_kurt = sum(r_j^4) / (n() * (sqrt(mean(r_j^2)))^4),
            rough_rms = sqrt(mean(r_j^2)),
            slope_cv = sd(slope) / mean(slope),
            flatness = dune_area / height_max,
            flatness2 = dune_width / height_max,
            sphericity = dune_area / (pi*height_max*height_max/2),
            rectangularity = dune_area / (dune_width * height_max),
            east = mean(east), north = mean(north))

# as factor for modelling
dune_shrub_metrics$dune <- factor(dune_shrub_metrics$dune) # need factors for s(by = ) input
dune_shrub_metrics$fence <- factor(dune_shrub_metrics$fence)


###
### I dont think the width metrics are a good idea - the dunes have such a long low slope into the swale,
### it seems attemtps to characterise their widths are just indicative of sampling length
###


# dune transect plots -----------------------------------------------------

transect_plot <- ggplot(data = dune_data_seq, aes(x = distance, y = height, group = Id)) +
  geom_line(aes(colour = fence), size = 1, alpha = 0.5) +
  ylab("Dune height (m)") + xlab("Distance along transect (m)") + theme_bw() +
  facet_grid(dune ~ ., scales = "fixed")
ggsave("plots/dune_transect_plot.png", transect_plot, width = 30, height = 20, units = "cm")

dune1 <- ggplot(data = dune_data_seq[dune_data_seq$dune == "dune1",], aes(x = distance, y = height, group = Id)) +
  geom_line(aes(colour = fence), size = 1, alpha = 0.5) +
  ylab("Dune height (m)") + xlab("Distance along transect (m)") + theme_bw()
ggsave("plots/dune1_plot.png", dune1, width = 30, height = 15, units = "cm")



# shrub metric plots ------------------------------------------------------

par(mfrow = c(3,2), mar = c(2,4,2,2), cex.lab = 1.3)
boxplot(shrub_dens ~ fence, data = dune_shrub_metrics, ylab = "Shrub density")
boxplot(height_99 ~ fence, data = dune_shrub_metrics, ylab = "Maximum height")
boxplot(rough_mean ~ fence, data = dune_shrub_metrics, ylab = "Roughness")
boxplot(slope_rms ~ fence, data = dune_shrub_metrics, ylab = "Slope-iness")
boxplot(slope_99 ~ fence, data = dune_shrub_metrics, ylab = "Slope99th")
boxplot(shape ~ fence, data = dune_shrub_metrics, ylab = "Shape")



# dune and shrub models ---------------------------------------------------

# shrubiness across fence - YES
shrubiness <- gam(formula = shrub_dens ~ fence + s(dune, bs = 're'),
                  correlation = corExp(form = ~ east + north, nugget = T),
                  data = dune_shrub_metrics, method = "REML")
summary(shrubiness)
plot(shrubiness, all.terms = T)

# shrubiness across space - YES
spatial_smooth <- gam(formula = shrub_dens ~ te(east,north), data = dune_shrub_metrics)
summary(spatial_smooth); plot(spatial_smooth, residuals = T)


# explore
metric_cors <- rbindlist(lapply(X = names(dune_shrub_metrics)[4:15], FUN = shrub_smooth, dune_shrub_metrics, F, F, T))


## models

# height
shrub_smooth("height_max", dune_shrub_metrics)
shrub_smooth_by_dune("height_max", dune_shrub_metrics)
shrub_smooth_dune_re("height_max", dune_shrub_metrics, plot = T)

# slope
shrub_smooth_by_dune("slope_rms", dune_shrub_metrics)
shrub_smooth_dune_re("slope_cv", dune_shrub_metrics, plot = T)

# shape
shrub_smooth("sphericity", dune_shrub_metrics)
shrub_smooth_by_dune("sphericity", dune_shrub_metrics)
shrub_smooth_dune_re("sphericity", dune_shrub_metrics, plot = T)

# shrub_smooth("rectangularity", dune_shrub_metrics)
# shrub_smooth_by_dune("rectangularity", dune_shrub_metrics)
# shrub_smooth_dune_re("rectangularity", dune_shrub_metrics, plot = T)

shrub_smooth("flatness", dune_shrub_metrics)
shrub_smooth_by_dune("flatness", dune_shrub_metrics)
shrub_smooth_dune_re("flatness", dune_shrub_metrics, plot = T)

# shrub_smooth("flatness2", dune_shrub_metrics)
# shrub_smooth_by_dune("flatness2", dune_shrub_metrics)
# shrub_smooth_dune_re("flatness2", dune_shrub_metrics, plot = T)

# roughness
shrub_smooth_by_dune("rough_rms", dune_shrub_metrics)
shrub_smooth_dune_re("rough_rms", dune_shrub_metrics, plot = T)


# width??
shrub_smooth_dune_re("dune_width", dune_shrub_metrics, plot = T)



