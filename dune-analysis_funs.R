my_violin <- function(y, data, ylab = y) {
  ggplot(data = data, mapping = aes_string(x = "fence", y = y)) +
    ylab(ylab) + xlab("Dingo fence") +
    geom_violin(draw_quantiles = c(0.25,0.5,0.75)) + theme_bw()
}

my_hist <- function(x, data, ylab = y, legend_off = T) {
  if (legend_off) {
    plt <- ggplot(data = data, mapping = aes_string(x = x, fill = "fence")) +
      ylab("density") + xlab(ylab) +
      geom_density(alpha = 0.3) + scale_fill_manual(values = c("red", "blue")) +
      guides(fill = FALSE) + theme_classic()
  } else {
    plt <- ggplot(data = data, mapping = aes_string(x = x, fill = "fence")) +
      ylab("density") + xlab(ylab) +
      geom_density(alpha = 0.3) + scale_fill_manual(values = c("red", "blue")) +
      theme_classic() +
      theme(
        legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6))
  }
  plt
}

shrub_dens_diff <- function(x, data) {
  dat_in <- filter(data, year == x, fence == "inside")
  dat_out <- filter(data, year == x, fence == "outside")
  data.frame(year = x,
             diff = dat_in$mean_shrub_density - dat_out$mean_shrub_density,
             std_err = sqrt((dat_in$std_dev^2 / dat_in$nsamps) + (dat_out$std_dev^2 / dat_out$nsamps)))
}

fence_dune_re <- function(y, data, summary = T, plot = F) {
  fm <- gamm(formula = as.formula(paste0(y," ~ fence + s(dune, bs = 're')")),
            correlation = corGaus(form = ~ east + north | dune, nugget = F),
            data = data, method = "REML")
  if(summary) {print(anova(fm$gam)); print(summary(fm$gam))}
  if (plot) {plot(fm$gam)}
  invisible(fm)
}

shrub_only <- function(y, data, bs = "'tp'") {
  fm <- gamm(formula = as.formula(paste0(y," ~ s(shrub_dens, bs = ",bs,") + s(dune, bs = 're')")),
             correlation = corGaus(form = ~ east + north | dune, nugget = F),
             data = data, method = "REML")
  print(plot_shrub_by_fence(fm = fm, orig_data = data))
  print(summary(fm$gam))
  invisible(fm)
}

fence_shrub_dune <- function(y, data) {
  require(effects)
  fm <- lme(fixed = as.formula(paste0(y," ~ fence + shrub_dens")),
            random = ~ 1 | dune,
            correlation = corGaus(form = ~ east + north | dune, nugget = F),
            data = data, method = "REML")
  print(summary(fm))
  invisible(fm)
}

shrub_smooth_dune_re <- function(y, data, summary = F, plot = F, bs = "'tp'") {
  fm <- gamm(formula = as.formula(paste0(y," ~ s(shrub_dens, bs = ",bs,") + fence + s(dune, bs = 're')")),
            correlation = corGaus(form = ~ east + north | dune, nugget = F),
            data = data, method = "REML")
  if(summary) {print(summary(fm$gam)); print(anova(fm$gam))}
  if (plot) {
    plot_shrub_by_fence(fm = fm, orig_data = data)
    #plot(fm$gam, all.terms = T)
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="inside"), col = "red", rm.ranef = T, legend_plot_all = "bottomright")
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="outside"), col = "blue", rm.ranef = T, add = T, legend_plot_all = "bottomright")
    #plot_smooth(fm$gam, view = "shrub_dens", plot_all = "fence", rm.ranef = T, xlim = c(0,0.015), print.summary = F)
  }
  invisible(fm)
}

shrub_smooth_dune_re2 <- function(y, data, summary = T, plot = F, bs = "'tp'") {
  fm <- gamm(formula = as.formula(paste0(y," ~ fence + s(shrub_dens, bs = ",bs,", by = fence) + s(dune, bs = 're')")),
             correlation = corGaus(form = ~ east + north | dune, nugget = F),
             data = data, method = "REML")
  if(summary) {print(summary(fm$gam)); print(anova(fm$gam))}
  if (plot) {
    plot_shrub_by_fence(fm = fm, orig_data = data, print = T)
    #plot(fm$gam, all.terms = T)
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="inside"), col = "red", rm.ranef = T, legend_plot_all = "bottomright")
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="outside"), col = "blue", rm.ranef = T, add = T, legend_plot_all = "bottomright")
    #plot_smooth(fm$gam, view = "shrub_dens", plot_all = "fence", rm.ranef = T, xlim = c(0,0.015), print.summary = F)
  }
  invisible(fm)
}

plot_shrub_by_fence <- function(fm, orig_data, xmin = 0, xmax = 150, ylabel = NULL, panel_lab = "a.", print = F) { #no.legend = F
  library(grid)
  fm <- fm$gam
  max_out <- max(orig_data$shrub_dens[orig_data$fence == "outside"])
  # have to give the dune value, even though it'll be excluded....
  newdat_in <- data.frame(shrub_dens = seq(xmin, xmax, length.out = 100), fence = rep("inside", 100), dune = rep("dune1", 100))
  newdat_out <- data.frame(shrub_dens = seq(xmin, max_out, length.out = 100), fence = rep("outside", 100), dune = rep("dune1", 100))
  pred_in <- predict(fm, newdata = newdat_in, se.fit = T, type = "response", exclude = "s(dune)")
  pred_out <- predict(fm, newdata = newdat_out, se.fit = T, type = "response", exclude = "s(dune)")
  preds <- rbind(
    data.frame(fit = as.numeric(pred_in$fit), se = as.numeric(pred_in$se.fit), 
               fence = "inside", shrub_dens = newdat_in$shrub_dens),
    data.frame(fit = as.numeric(pred_out$fit), se = as.numeric(pred_out$se.fit), 
               fence = "outside", shrub_dens = newdat_out$shrub_dens)
  )
  #preds <- filter(preds, !(fence == "outside" & shrub_dens > max_out))
  #preds <- preds[preds$fence == "outside" & preds$shrub_dens > max(orig_data$shrub_dens[orig_data$fence == "outside"]),] <- NA
  if (is.null(ylabel)) {ylabel <- as.character(fm$formula)[2]}
  predplt <- ggplot(data = preds, aes(x = shrub_dens)) +
    geom_ribbon(aes(x = shrub_dens, ymax = fit + se, ymin = fit - se, fill = fence), alpha = "0.1") +
    geom_line(aes(y = fit, colour = fence)) +
    geom_rug(mapping = aes(x = shrub_dens, colour = fence), data = orig_data) +
    xlim(c(xmin, xmax)) +
    scale_color_manual(values = c("red","blue")) + scale_fill_manual(values = c("red","blue")) + 
    #geom_vline(xintercept = max(orig_data$shrub_dens[orig_data$fence == "outside"]), colour = "blue", linetype = "dotted") +
    labs(x = "Shrub Density (/Ha)", y = ylabel) + theme_classic()
  # if (no.legend) {
  #   predplt <- predplt + guides(fill = FALSE, colour = FALSE)
  # }
  predplt <- arrangeGrob(predplt, 
                         top = grid::textGrob(panel_lab,
                                        x = unit(0.1, "npc"), y = unit(0.9, "npc"),
                                        just = c("left","top"), gp = gpar(col = "black", fontsize = 14)))
  if(print){print(predplt)}
  invisible(predplt)
}

get_paper_values <- function(fm) {
  fm_sum <- summary(fm$gam)
  print(round(fm_sum$p.table, 3))
  print(round(fm_sum$pTerms.table, 3))
  print(round(fm_sum$s.table, 3))
}

my.gam.check <- function(fm) {
  par(mfrow=c(2,2))
  try(gam.check(fm$gam))
  par(mfrow=c(1,1))
}







