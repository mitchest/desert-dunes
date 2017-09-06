library(dplyr)
library(data.table)
library(mgcv)
library(nlme)
library(ggplot2)
library(gridExtra)
#library(moments)
#library(foreign)

source("dune-analysis_funs.R")

#dune_data <- read.dbf("transects_dunes_points_dtm_slope_dens.dbf", as.is=T)
#save(dune_data, file="dune_data.RData")
load("dune_data.RData")



# add transect length/height -----------------------------------------------------

dune_data_seq <- dune_data %>%
  filter(h_asl > 0) %>%
  filter(!dune %in% c("dune5", "dune7a", "dune2"),
         !Id %in% c(157,158,159,161)) %>% # few odd length transects
  rename(shrub_dens = RASTERVALU) %>%
  group_by(Id) %>% #transect number
  mutate(distance = seq(0, by=0.1, length.out = n()),
         height = sapply(h_asl, function(x) x - min(h_asl)),
         r_j = height - mean(height),
         q_j = slope - mean(slope)) %>%
  ungroup() %>%
  filter(height > 0.20,
         !(distance > 175 & dune == "dune2")) # lengths too uneven

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
            heightdev_rms = sqrt(mean(r_j^2)),
            slope_rms = sqrt(mean(q_j^2)),
            slope_cv = sd(slope) / mean(slope),
            flatness = log(dune_area / height_max),
            flatness2 = dune_width / height_max,
            sphericity = dune_area / (pi*height_max*height_max/2),
            rectangularity = dune_area / (dune_width * height_max),
            east = median(east), north = median(north)) %>%
  filter(height_max > 2) # any odd little transects

# as factor for modelling
dune_shrub_metrics$dune <- factor(dune_shrub_metrics$dune) # need factors for some mgcv args
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



# shrub metric plots ------------------------------------------------------

boxplot(shrub_dens ~ fence, data = dune_shrub_metrics, ylab = "Shrub density")
boxplot(height_max ~ fence, data = dune_shrub_metrics, ylab = "Maximum height")
boxplot(heightdev_rms ~ fence, data = dune_shrub_metrics, ylab = "Height dev.")
boxplot(slope_rms ~ fence, data = dune_shrub_metrics, ylab = "Slope-iness")
boxplot(slope_cv ~ fence, data = dune_shrub_metrics, ylab = "Slope CV")
boxplot(flatness ~ fence, data = dune_shrub_metrics, ylab = "Flatness")
boxplot(rectangularity ~ fence, data = dune_shrub_metrics, ylab = "Rectangularity")


# dune and shrub models ---------------------------------------------------

# shrubiness across fence - YES
dune_shrub_metrics$shrub_dens_log <- log1p(dune_shrub_metrics$shrub_dens)
shrubiness <- gamm(formula = shrub_dens ~ fence + s(dune, bs = 're'),
                  correlation = corGaus(form = ~ east + north, nugget = F),
                  data = dune_shrub_metrics, method = "REML")
summary(shrubiness$gam)
plot(shrubiness$lme, all.terms = T)
plot(predict(shrubiness$gam, type="response") ~ dune_shrub_metrics$fence)

# shrubiness across space - YES
spatial_smooth <- gam(formula = shrub_dens ~ te(east,north), data = dune_shrub_metrics)
summary(spatial_smooth); plot(spatial_smooth, residuals = T)



## models
make_plots <- F
print_summary <- F

# height
## height max is probably not that good an idea??
#fence_dune_re("height_max", dune_shrub_metrics)
height <- shrub_smooth_dune_re2("height_max", dune_shrub_metrics, plot = make_plots, summary = print_summary)
my.gam.check(height)
get_paper_values(height)

#fence_dune_re("heightdev_rms", dune_shrub_metrics)
heightdev <- shrub_smooth_dune_re2("heightdev_rms", dune_shrub_metrics, plot = make_plots, summary = print_summary)
my.gam.check(heightdev)
get_paper_values(heightdev)

# slope/roughness
#fence_dune_re("slope_rms", dune_shrub_metrics)
sloperms <- shrub_smooth_dune_re2("slope_rms", dune_shrub_metrics, plot = make_plots, summary = print_summary)
my.gam.check(sloperms)
get_paper_values(sloperms)

#fence_dune_re("slope_cv", dune_shrub_metrics)
slopecv <- shrub_smooth_dune_re2("slope_cv", dune_shrub_metrics, plot = make_plots, summary = print_summary)
my.gam.check(slopecv)
get_paper_values(slopecv)

# shape
#fence_dune_re("flatness", dune_shrub_metrics)
flat <- shrub_smooth_dune_re2("flatness", dune_shrub_metrics, plot = make_plots, summary = print_summary)
my.gam.check(flat)
get_paper_values(flat)

#fence_dune_re("rectangularity", dune_shrub_metrics)
rect <- shrub_smooth_dune_re2("rectangularity", dune_shrub_metrics, plot = make_plots, summary = print_summary)
my.gam.check(rect)
get_paper_values(rect)

# shrub_smooth_dune_re2("flatness2", dune_shrub_metrics, plot = T)
# shrub_smooth_dune_re2("sphericity", dune_shrub_metrics, plot = T)

heightplt <- plot_shrub_by_fence(height, dune_shrub_metrics, ylabel = "Maximum height")
heightdevplt <- plot_shrub_by_fence(heightdev, dune_shrub_metrics, ylabel = "Height deviation")
slopermsplt <- plot_shrub_by_fence(sloperms, dune_shrub_metrics, ylabel = "Roughness (slope rms)")
slopecvplt <- plot_shrub_by_fence(slopecv, dune_shrub_metrics, ylabel = "Roughness (slope CV)")
flatplt <- plot_shrub_by_fence(flat, dune_shrub_metrics, ylabel = "Flatness")
rectplt <- plot_shrub_by_fence(rect, dune_shrub_metrics, ylabel = "Rectangularity")

shrub_partials <- grid.arrange(heightplt, heightdevplt, slopermsplt, slopecvplt, flatplt, rectplt)


