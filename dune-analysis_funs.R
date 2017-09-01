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

shrub_smooth_dune_re <- function(y, data, summary = T, plot = F) {
  fm <- gam(formula = as.formula(paste0(y," ~ fence + s(shrub_dens, bs = 'tp') + s(dune, bs = 're')")),
            correlation = corExp(form = ~ east + north, nugget = T),
            data = data, method = "REML")
  if(summary) {print(summary(fm))}
  if (plot) {plot(fm, select = 1)}
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