library(ggplot2)
library(gridExtra)
library(nlme)
library(mgcv)
library(data.table)
library(dplyr)
library(grid)
#library(moments)
#library(foreign)

source("dune-analysis_funs.R")

#dune_data <- read.dbf("transects_dunes_points_dtm_slope_dens.dbf", as.is=T)
#save(dune_data, file="dune_data.RData")
load("dune_data.RData")

shrub_ts_data <- read.csv("historical_shrub_data.csv")



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
  summarise(shrub_dens = mean(shrub_dens) * 10000, # transform to shrubs/hectare
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


# secondary data set with more balanced shrub density insdie/outside
dune_shrub_metrics2 <- dune_shrub_metrics %>%
  filter(shrub_dens < max(dune_shrub_metrics$shrub_dens[dune_shrub_metrics$fence=="outside"]))



# dune transect plots -----------------------------------------------------

transect_plot <- ggplot(data = dune_data_seq, aes(x = distance, y = height, group = Id)) +
  geom_line(aes(colour = fence), size = 1, alpha = 0.2) +
  scale_color_manual(values = c("red", "blue")) +
  ylab("Dune height (m)") + xlab("Distance along transect (m)") + theme_bw() +
  facet_grid(dune ~ ., scales = "fixed")
ggsave("plots/dune_transect_plot.png", transect_plot, width = 30, height = 20, units = "cm")



# shrub time series -------------------------------------------------------

shrub_timeseries <- shrub_ts_data %>%
  mutate(#year = factor(year),
         shrub_density = count / (area/10000)) %>% # shrubs/hectare
  group_by(fence, year) %>%
  summarise(mean_shrub_density = mean(shrub_density),
            std_dev = sd(shrub_density),
            nsamps = n(),
            std_err = std_dev / sqrt(nsamps))

shrub_diff_ts <- rbindlist(lapply(X = unique(shrub_timeseries$year), FUN = shrub_dens_diff, shrub_timeseries))

ggplot(shrub_timeseries, aes(x = year, y = mean_shrub_density, colour = fence)) +
  geom_line(aes(group = fence), size = 1) +
  geom_errorbar(aes(ymax = mean_shrub_density + std_err, ymin = mean_shrub_density - std_err), size = 0.25) +
  scale_x_continuous(breaks = unique(shrub_timeseries$year)) +
  theme_classic()

history_plt <- ggplot(shrub_diff_ts, aes(x = year, y = diff)) +
  geom_line(size = 1) +
  geom_errorbar(aes(ymax = diff + std_err, ymin = diff - std_err), size = 0.25, width = 0.4) +
  scale_x_continuous(breaks = unique(shrub_diff_ts$year)) +
  ylim(-10, 80) +
  geom_hline(yintercept = 0, linetype = 2) +
  ylab("Shrub density (shrubs/Ha)") + xlab("Image date") +
  #ggtitle("Difference in shrub density inside and outside the Dingo fence") +
  # annotate("text", x = 1948, y = 8, 
  #          label = "Difference in shrub density inside and outside the Dingo fence", size = 6, hjust = 0) +
  #annotate("text", x = 1992, y = -5, label = sprintf('\u2193'), size = 9) +
  geom_segment(aes(x = 1992, y = -2, xend = 1992, yend = -8), arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", x = 1994, y = -5, label = "more shrubs outside", size = 4, hjust = 0) +
  #annotate("text", x = 1992, y = 5, label = sprintf('\u2191'), size = 9) +
  geom_segment(aes(x = 1992, y = 2, xend = 1992, yend = 8), arrow = arrow(length = unit(0.03, "npc"))) +
  annotate("text", x = 1994, y = 5, label = "more shrubs inside", size = 4, hjust = 0) +
  theme_classic()
ggsave(filename = "plots/history_plot.pdf", device = "pdf", plot = history_plt, width = 9, height = 6)


# shrub metric plots ------------------------------------------------------

boxplot(shrub_dens ~ fence, data = dune_shrub_metrics, ylab = "Shrub density")
my_violin("shrub_dens", dune_shrub_metrics, ylab = "Shrub density")

shrub_plt <- my_hist("shrub_dens", dune_shrub_metrics, ylab = "Shrub density", legend_off = F)
ggsave("plots/shrubiness.png", device = "png", plot = shrub_plt, width = 6, height = 4)

shrub_model <- lme(fixed = shrub_dens ~ fence, random = ~ 1 | dune,
                   correlation = corGaus(form = ~ east + north | dune, nugget = F),
                   data = dune_shrub_metrics)
summary(shrub_model)

par(mfrow=c(3,2))
boxplot(height_max ~ fence, data = dune_shrub_metrics, ylab = "Maximum height")
boxplot(heightdev_rms ~ fence, data = dune_shrub_metrics, ylab = "Height dev.")
boxplot(slope_rms ~ fence, data = dune_shrub_metrics, ylab = "Slope-iness")
boxplot(slope_cv ~ fence, data = dune_shrub_metrics, ylab = "Slope CV")
boxplot(flatness ~ fence, data = dune_shrub_metrics, ylab = "Flatness")
boxplot(rectangularity ~ fence, data = dune_shrub_metrics, ylab = "Rectangularity")

metric_plt <- grid.arrange(
  my_violin("height_max", data = dune_shrub_metrics, ylab = "Maximum height"),
  my_violin("heightdev_rms", data = dune_shrub_metrics, ylab = "Deviation in height profile"),
  my_violin("slope_rms", data = dune_shrub_metrics, ylab = "Absolute range in slope"),
  my_violin("slope_cv", data = dune_shrub_metrics, ylab = "Slope variability"),
  my_violin("flatness", data = dune_shrub_metrics, ylab = "Flatness"),
  my_violin("rectangularity", data = dune_shrub_metrics, ylab = "Rectangularity"),
  ncol = 3)
ggsave("plots/dune_metric_plots.png", device = "png", plot = metric_plt, width = 12, height = 8)

metric_dens_plt <- grid.arrange(
  my_hist("height_max", data = dune_shrub_metrics, ylab = "Maximum height"),
  my_hist("heightdev_rms", data = dune_shrub_metrics, ylab = "Deviation in height profile"),
  my_hist("slope_rms", data = dune_shrub_metrics, ylab = "Absolute range in slope"),
  my_hist("slope_cv", data = dune_shrub_metrics, ylab = "Slope variability", legend_off = F),
  my_hist("flatness", data = dune_shrub_metrics, ylab = "Flatness"),
  my_hist("rectangularity", data = dune_shrub_metrics, ylab = "Rectangularity"),
  ncol = 3)
ggsave("plots/dune_metric_dens-plots.pdf", device = "pdf", plot = metric_dens_plt, width = 12, height = 8)

# geomorph clustering -----------------------------------------------------
library(ggbiplot)

pca_dat <- dune_shrub_metrics %>%
  ungroup() %>%
  select(height_max, heightdev_rms, slope_rms, slope_cv, flatness, rectangularity)

pca <- prcomp(pca_dat, scale. = T)

# pca_plot_data <- data.frame(pca$x, shrub_density = dune_shrub_metrics$shrub_dens,
#                             dune = dune_shrub_metrics$dune, fence = dune_shrub_metrics$fence)
# ggplot(data = pca_plot_data, mapping = aes(x = PC1, y = PC2, col = fence)) +
#   geom_point(size = 2, alpha = 0.5) +
#   theme_classic()

pca_plt <- ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = dune_shrub_metrics$dune,
         ellipse = F, circle = F, alpha = 0) +
  geom_point(aes(shape = dune_shrub_metrics$fence, colour = dune_shrub_metrics$dune), size = 3) + 
  scale_color_brewer(palette = "Set1", name = 'Dune number') +
  scale_shape_discrete(name = 'Dingo fence') +
  theme_classic()
ggsave("plots/dune_metric_pca.png", device = "png", plot = pca_plt, width = 12, height = 10)

ggbiplot(pca, obs.scale = 1, var.scale = 1,
         groups = dune_shrub_metrics$shrub_dens,
         ellipse = F, circle = F, alpha = 0) +
  geom_point(aes(colour = dune_shrub_metrics$shrub_dens, shape = dune_shrub_metrics$fence), size = 3) + 
  scale_color_gradient(low = "coral", high = "green", name = "Shrub density") +
  scale_shape_discrete(name = 'Dingo fence') +
  theme_classic()



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
# make_plots <- F
# print_summary <- F

height <- shrub_smooth_dune_re("height_max", dune_shrub_metrics, bs = "'tp'")
heightdev <- shrub_smooth_dune_re("heightdev_rms", dune_shrub_metrics, bs = "'tp'")
sloperms <- shrub_smooth_dune_re("slope_rms", dune_shrub_metrics, bs = "'tp'")
slopecv <- shrub_smooth_dune_re("slope_cv", dune_shrub_metrics, bs = "'tp'")
flat <- shrub_smooth_dune_re("flatness", dune_shrub_metrics, bs = "'tp'")
rect <- shrub_smooth_dune_re("rectangularity", dune_shrub_metrics, bs = "'tp'")

my.gam.check(heightdev)
my.gam.check(height)
my.gam.check(sloperms)
my.gam.check(slopecv)
my.gam.check(flat)
my.gam.check(rect)

get_paper_values(height)
get_paper_values(heightdev)
get_paper_values(sloperms)
get_paper_values(slopecv)
get_paper_values(flat)
get_paper_values(rect)

# parial plots
heightplt <- plot_shrub_by_fence(height, dune_shrub_metrics, ylabel = "Maximum height", panel_lab = "a.")
heightdevplt <- plot_shrub_by_fence(heightdev, dune_shrub_metrics, ylabel = "Height deviation", panel_lab = "b.")
slopermsplt <- plot_shrub_by_fence(sloperms, dune_shrub_metrics, ylabel = "Roughness (slope rms)", panel_lab = "x.")
slopecvplt <- plot_shrub_by_fence(slopecv, dune_shrub_metrics, ylabel = "Roughness (slope CV)", panel_lab = "c.")
flatplt <- plot_shrub_by_fence(flat, dune_shrub_metrics, ylabel = "Flatness", panel_lab = "y.")
rectplt <- plot_shrub_by_fence(rect, dune_shrub_metrics, ylabel = "Rectangularity", panel_lab = "d.")

shrub_partials_all <- grid.arrange(heightplt, heightdevplt, slopermsplt, slopecvplt, flatplt, rectplt)
shrub_partials <- grid.arrange(heightplt, heightdevplt, slopecvplt, rectplt)
ggsave("plots/shrub_partials_all.png", device = "png", plot = shrub_partials_all, width = 12, height = 10)
ggsave("plots/shrub_partials.pdf", device = "pdf", plot = shrub_partials, width = 12, height = 7)

