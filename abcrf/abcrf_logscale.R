# Control: abcrf regAbcrf on log10 of the scale parameters (Raynal et al. 2019 allow any h(theta)),
# to separate the effect of the response scale from the change of implementation.
suppressPackageStartupMessages({ library(abcrf); library(jsonlite) })
NT <- 1000; NC <- 16
d <- read.csv("reftable_sim9.csv"); o <- read.csv("observed_sim9.csv", row.names = 1)
st <- setdiff(names(d), c("model", grep("^p_", names(d), value = TRUE)))
two <- which(d$model %in% c("H2", "H3")); con <- which(d$model == "H3")
est <- function(rows, y, thr) {
  set.seed(7); dd <- data.frame(y = log10(y), d[rows, st])
  rf <- regAbcrf(y ~ ., data = dd, ntree = NT, paral = TRUE, ncores = NC)
  p <- predict(rf, obs = o["primary", st], training = dd, paral = TRUE, ncores = NC, rf.weights = TRUE, post.err.med = TRUE)
  w <- as.numeric(p$weights)
  c(list(median = 10^p$med, lo = 10^p$quantiles[1], hi = 10^p$quantiles[2], r2 = rf$model.rf$r.squared,
         prior_NMAE_med_log10 = p$prior.NMAE.med),
    setNames(lapply(thr, function(b) sum(w[y > b])), paste0("P_gt_", thr)),
    setNames(lapply(thr, function(b) mean(y > b)), paste0("prior_P_gt_", thr)))
}
R <- list(T_years = est(two, d$p_T_years[two], c(10e3, 19e3, 26.5e3, 115e3)),
          Nm_contact = est(con, d$p_Ne[con] * d$p_mcross[con], c(1, 10)),
          mcross_contact = est(con, d$p_mcross[con], c(0.001, 0.01)))
write(toJSON(R, auto_unbox = TRUE, digits = 6, pretty = TRUE), "../figdata/abcrf_logscale.json"); print(toJSON(R, auto_unbox = TRUE, digits = 4, pretty = TRUE))
