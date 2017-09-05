fence_dune_re <- function(y, data, summary = T, plot = F) {
  fm <- gamm(formula = as.formula(paste0(y," ~ fence + s(dune, bs = 're')")),
            correlation = corGaus(form = ~ east + north | dune, nugget = F),
            data = data, method = "REML")
  if(summary) {print(anova(fm$gam)); print(summary(fm$gam))}
  if (plot) {plot(fm$gam)}
  invisible(fm)
}

shrub_smooth_dune_re <- function(y, data, summary = T, plot = F, bs = "'ts'") {
  fm <- gamm(formula = as.formula(paste0(y," ~ fence + s(shrub_dens, bs = ",bs,") + s(dune, bs = 're')")),
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
    plot_shrub_by_fence(fm = fm, orig_data = data)
    #plot(fm$gam, all.terms = T)
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="inside"), col = "red", rm.ranef = T, legend_plot_all = "bottomright")
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="outside"), col = "blue", rm.ranef = T, add = T, legend_plot_all = "bottomright")
    #plot_smooth(fm$gam, view = "shrub_dens", plot_all = "fence", rm.ranef = T, xlim = c(0,0.015), print.summary = F)
  }
  invisible(fm)
}

plot_shrub_by_fence <- function(fm, orig_data, xmin = 0, xmax = 0.0155, ylabel = NULL) {
  fm <- fm$gam
  # have to give the dune value, even though it'll be excluded....
  newdat_in <- data.frame(shrub_dens = seq(xmin, xmax, length.out = 100), fence = rep("inside", 100), dune = rep("dune1", 100))
  newdat_out <- data.frame(shrub_dens = seq(xmin, xmax, length.out = 100), fence = rep("outside", 100), dune = rep("dune1", 100))
  pred_in <- predict(fm, newdata = newdat_in, se.fit = T, type = "response", exclude = "s(dune)")
  pred_out <- predict(fm, newdata = newdat_out, se.fit = T, type = "response", exclude = "s(dune)")
  preds <- rbind(
    data.frame(fit = as.numeric(pred_in$fit), se = as.numeric(pred_in$se.fit), 
               fence = "inside", shrub_dens = newdat_in$shrub_dens),
    data.frame(fit = as.numeric(pred_out$fit), se = as.numeric(pred_out$se.fit), 
               fence = "outside", shrub_dens = newdat_out$shrub_dens)
  )
  if (is.null(ylabel)) {ylabel <- as.character(fm$formula)[2]}
  predplt <- ggplot(data = preds, aes(x = shrub_dens)) +
    geom_ribbon(aes(x = shrub_dens, ymax = fit + se, ymin = fit - se, fill = fence), alpha = "0.1") +
    geom_line(aes(y = fit, colour = fence)) +
    geom_rug(mapping = aes(x = shrub_dens), data = orig_data) +
    xlim(c(xmin, xmax)) +
    scale_color_manual(values = c("red","blue")) + scale_fill_manual(values = c("red","blue")) + 
    xlab("Shrub density (shrubs/m2)") + ylab(ylabel) + theme_bw()
  print(predplt)
  invisible(predplt)
}

