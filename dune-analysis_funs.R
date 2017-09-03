fence_dune_re <- function(y, data, summary = T, plot = F) {
  fm <- gam(formula = as.formula(paste0(y," ~ fence + s(dune, bs = 're')")),
            correlation = corExp(form = ~ east + north, nugget = T),
            data = data, method = "REML")
  if(summary) {print(summary(fm))}
  if (plot) {plot(fm)}
  invisible(fm)
}

shrub_smooth <- function(y, data, summary = T, plot = F, test.p = F) {
  fm <- gam(formula = as.formula(paste0(y," ~ fence + s(shrub_dens, bs = 'tp')")),
            correlation = corExp(form = ~ east + north, nugget = T),
            data = data, method = "REML")
  if(test.p) {
    fm_sum <- summary(fm)
    return(data.frame(var = y, fence = round(fm_sum$p.pv[2], 4), shrub = round(fm_sum$s.pv,4)))
    }
  if(summary) {print(summary(fm))}
  if (plot) {plot(fm)}
  invisible(fm)
}

shrub_smooth_by_dune <- function(y, data, summary = T, plot = F) {
  fm <- gam(formula = as.formula(paste0(y," ~ fence + s(shrub_dens, bs = 'tp', by = dune)")),
            correlation = corExp(form = ~ east + north, nugget = T),
            data = data, method = "REML")
  if(summary) {print(summary(fm))}
  if (plot) {plot(fm)}
  invisible(fm)
}

shrub_smooth_dune_re <- function(y, data, summary = T, plot = F, bs = "'tp'") {
  fm <- gam(formula = as.formula(paste0(y," ~ fence + s(shrub_dens, bs = ",bs,") + s(dune, bs = 're')")),
            correlation = corExp(form = ~ east + north, nugget = T),
            data = data, method = "REML")
  if(summary) {print(summary(fm)); print(anova(fm))}
  if (plot) {
    #plot(fm, all.terms = T)
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="inside"), col = "red", rm.ranef = T, legend_plot_all = "bottomright")
    #plot_smooth(fm, view = "shrub_dens", cond = list(fence="outside"), col = "blue", rm.ranef = T, add = T, legend_plot_all = "bottomright")
    plot_smooth(fm, view = "shrub_dens", plot_all = "fence", rm.ranef = T)
  }
  invisible(fm)
}

mulga_smooth_by_dune <- function(y, data, summary = T, plot = F) {
  fm <- gam(formula = as.formula(paste0(y," ~ fence + s(mulga_dens, bs = 'tp', by = dune)")),
            correlation = corExp(form = ~ east + north, nugget = T),
            data = data, method = "REML")
  if(summary) {print(summary(fm))}
  if (plot) {plot(fm)}
  invisible(fm)
}