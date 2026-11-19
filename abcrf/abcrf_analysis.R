# ABC-RF model choice and parameter estimation with the reference implementation, abcrf 2.0
# (Pudlo et al. 2016; Raynal et al. 2019), replacing the scikit-learn reimplementation (18 Nov 2026).
#
# Inputs: reference tables exported from McCartney by scripts/abc/export_reftables.py
#   reftable_sim9.csv  literature mutation-rate prior (primary), 1,392 sims per history, 15 statistics
#   reftable_sim8.csv  broad mutation-rate prior (sensitivity)
# Output: ../figdata/abcrf_results.json, ../figdata/abcrf_weights_*.csv (posterior weights for Figure 4C)
#
# Standard abcrf steps: abcrf() with LDA axes added (package default), grouped scenarios via `group`,
# predict.abcrf() for votes and posterior probability, prior error rate from the OOB forest;
# regAbcrf() on the raw parameter, predict.regAbcrf() for median, 95% interval, OOB-median NMAE
# (global = prior, local = posterior, post.err.med = TRUE) and the RF weights.
# BESPOKE parts (no abcrf function exists): held-out interval coverage, posterior tail areas and mode
# from the RF weights, the contact-detection power curve, and the zero-migration calibration.
suppressPackageStartupMessages({ library(abcrf); library(jsonlite) })
setwd(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))))
NT <- 1000; NC <- 16

load_tab <- function(tag) {
  d <- read.csv(sprintf("reftable_%s.csv", tag))
  o <- read.csv(sprintf("observed_%s.csv", tag), row.names = 1)
  stats <- setdiff(names(d), c("model", grep("^p_", names(d), value = TRUE)))
  d$model <- factor(d$model, levels = c("H1", "H1b", "H2", "H3"))
  list(d = d, obs = o[, stats], stats = stats)
}
S9 <- load_tab("sim9"); S8 <- load_tab("sim8")

# ---------------- model choice ----------------
choice <- function(S, rows, groups = list(), seed = 1) {
  set.seed(seed)
  dd <- data.frame(model = droplevels(S$d$model[rows]), S$d[rows, S$stats])
  rf <- abcrf(model ~ ., data = dd, group = groups, ntree = NT, paral = TRUE, ncores = NC)
  p <- predict(rf, obs = S$obs["primary", ], training = dd, ntree = NT, paral = TRUE, ncores = NC)
  list(rf = rf, votes = as.numeric(p$vote) / NT, labels = colnames(p$vote),
       selected = as.character(p$allocation), posterior = as.numeric(p$post.prob),
       prior_error = rf$prior.err)
}
replicate10 <- function(S, rows, groups = list()) {
  R <- lapply(1:10, function(s) choice(S, rows, groups, seed = 100 + s))
  V <- do.call(rbind, lapply(R, `[[`, "votes"))
  list(labels = R[[1]]$labels, votes_mean = colMeans(V), votes_sd = apply(V, 2, sd),
       selected = sapply(R, `[[`, "selected"),
       posterior_mean = mean(sapply(R, `[[`, "posterior")), posterior_sd = sd(sapply(R, `[[`, "posterior")),
       prior_error_mean = mean(sapply(R, `[[`, "prior_error")), prior_error_sd = sd(sapply(R, `[[`, "prior_error")),
       confusion = unclass(R[[1]]$rf$model.rf$confusion.matrix), first = R[[1]])
}
all9 <- seq_len(nrow(S9$d)); two9 <- which(S9$d$model %in% c("H2", "H3"))
out <- list(abcrf_version = as.character(packageVersion("abcrf")), ranger_version = as.character(packageVersion("ranger")),
            ntree = NT, n_sims = nrow(S9$d), stats = S9$stats)
cat("model choice, sim9\n")
m4 <- replicate10(S9, all9)
m1 <- replicate10(S9, all9, groups = list(c("H1", "H1b"), c("H2", "H3")))
m2 <- replicate10(S9, two9)
strip <- function(m) m[setdiff(names(m), "first")]
out$four <- strip(m4); out$stage1 <- strip(m1); out$stage2 <- strip(m2)
# per-history recovery (OOB) and variable importance from the first four-way forest
cm <- m4$confusion[, 1:4]; out$four$recall <- diag(cm) / rowSums(cm)
imp <- m4$first$rf$model.rf$variable.importance; out$four$importance <- as.list(sort(imp, decreasing = TRUE))
# power: share of contact sims called contact, by true migrants per generation (BESPOKE summary of OOB calls)
rf2 <- m2$first$rf$model.rf
oobcall <- rf2$predictions; truth <- droplevels(S9$d$model[two9])
Nm <- S9$d$p_Ne[two9] * S9$d$p_mcross[two9]
h3 <- truth == "H3"
bins <- list(c(0, .2), c(.2, 1), c(1, 5), c(5, Inf)); lab <- c("< 0.2", "0.2-1", "1-5", "> 5")
out$power <- lapply(seq_along(bins), function(i) {
  k <- h3 & Nm >= bins[[i]][1] & Nm < bins[[i]][2]
  list(bin = lab[i], n = sum(k), called_contact = mean(as.character(oobcall[k]) == "H3")) })
# broad-prior sensitivity (single forest, as before)
all8 <- seq_len(nrow(S8$d)); two8 <- which(S8$d$model %in% c("H2", "H3"))
c81 <- choice(S8, all8, groups = list(c("H1", "H1b"), c("H2", "H3"))); c82 <- choice(S8, two8)
out$broad_stage1 <- c81[c("labels", "votes", "selected", "posterior", "prior_error")]
out$broad_stage2 <- c82[c("labels", "votes", "selected", "posterior", "prior_error")]

# ---------------- parameter estimation ----------------
QS <- seq(0.0005, 0.9995, by = 0.0005)
estimate <- function(S, rows, resp, obsrow = "primary", tag = NULL, seed = 7) {
  set.seed(seed)
  dd <- data.frame(y = resp, S$d[rows, S$stats])
  rf <- regAbcrf(y ~ ., data = dd, ntree = NT, paral = TRUE, ncores = NC)
  p <- predict(rf, obs = S$obs[obsrow, ], training = dd, quantiles = c(0.025, 0.975), paral = TRUE, ncores = NC,
               rf.weights = TRUE, post.err.med = TRUE)
  w <- as.numeric(p$weights)
  if (!is.null(tag)) write.csv(data.frame(value = resp, weight = w), sprintf("../figdata/abcrf_weights_%s.csv", tag), row.names = FALSE)
  list(median = as.numeric(p$med), lo = as.numeric(p$quantiles[1]), hi = as.numeric(p$quantiles[2]),
       mean = as.numeric(p$expectation), r2 = rf$model.rf$r.squared,
       prior_NMAE_med = as.numeric(p$prior.NMAE.med), post_NMAE_med = as.numeric(p$post.NMAE.med),
       prior_NMAE_mean = as.numeric(p$prior.NMAE.mean), post_NMAE_mean = as.numeric(p$post.NMAE.mean),
       w = w, resp = resp, rf = rf, dd = dd)
}
tail_p <- function(e, b) sum(e$w[e$resp > b])           # BESPOKE: posterior tail area from RF weights
mode_log <- function(e) {                                # BESPOKE: mode of the weighted density on log10 scale
  g <- seq(3, 7, length.out = 2000); x <- log10(e$resp)
  d <- sapply(g, function(z) sum(e$w * dnorm(z, x, 0.12))); 10^g[which.max(d)] }
summ <- function(e) {
  r <- e[c("median", "lo", "hi", "mean", "r2", "prior_NMAE_med", "post_NMAE_med", "prior_NMAE_mean", "post_NMAE_mean")]
  r$mode <- mode_log(e)
  for (b in c(10e3, 19e3, 26.5e3, 115e3)) r[[sprintf("P_older_than_%gkyr", b / 1000)]] <- tail_p(e, b)
  r }
cat("parameters, sim9\n")
Ty <- S9$d$p_T_years
eA <- estimate(S9, two9, Ty[two9], tag = "avg"); out$split_model_averaged <- summ(eA)
iso <- which(S9$d$model == "H2"); con <- which(S9$d$model == "H3")
out$split_isolation <- summ(estimate(S9, iso, Ty[iso], tag = "iso"))
out$split_contact   <- summ(estimate(S9, con, Ty[con], tag = "con"))
out$split_all_loci  <- summ(estimate(S9, two9, Ty[two9], obsrow = "all_loci"))
fl <- two9[Ty[two9] >= 1e4]; out$split_floor_1e4 <- summ(estimate(S9, fl, Ty[fl]))
gS <- two9[S9$d$p_g[two9] < 180]; gL <- two9[S9$d$p_g[two9] >= 180]
out$split_g_lt_180 <- summ(estimate(S9, gS, Ty[gS])); out$split_g_ge_180 <- summ(estimate(S9, gL, Ty[gL]))
for (k in c("mu_yr", "g", "Ne", "Nanc")) {
  e <- estimate(S9, two9, S9$d[[paste0("p_", k)]][two9])
  out[[paste0("param_", k)]] <- e[c("median", "lo", "hi", "r2")] }
# held-out coverage of the 95% interval (BESPOKE check; abcrf has no coverage function)
set.seed(11); hold <- sample(two9, 150); keep <- setdiff(two9, hold)
ddk <- data.frame(y = Ty[keep], S9$d[keep, S9$stats])
rk <- regAbcrf(y ~ ., data = ddk, ntree = NT, paral = TRUE, ncores = NC)
pk <- predict(rk, obs = S9$d[hold, S9$stats], training = ddk, quantiles = c(0.025, 0.975), paral = TRUE, ncores = NC)
out$coverage95 <- mean(Ty[hold] >= pk$quantiles[, 1] & Ty[hold] <= pk$quantiles[, 2])
# migration under contact
NmC <- S9$d$p_Ne[con] * S9$d$p_mcross[con]; eN <- estimate(S9, con, NmC)
out$Nm_contact <- c(eN[c("median", "lo", "hi", "r2", "prior_NMAE_med", "post_NMAE_med")],
                    list(P_gt_1 = tail_p(eN, 1), P_gt_1_prior = mean(NmC > 1), P_gt_10 = tail_p(eN, 10), P_gt_10_prior = mean(NmC > 10)))
eM <- estimate(S9, con, S9$d$p_mcross[con])
out$mcross_contact <- c(eM[c("median", "lo", "hi", "r2")],
                        list(P_gt_0.01 = tail_p(eM, 0.01), P_gt_0.01_prior = mean(S9$d$p_mcross[con] > 0.01),
                             P_gt_0.001 = tail_p(eM, 0.001), P_gt_0.001_prior = mean(S9$d$p_mcross[con] > 0.001)))
# zero-migration calibration (BESPOKE): Nm estimator applied to isolation simulations
pz <- predict(eN$rf, obs = S9$d[iso, S9$stats], training = eN$dd, quantiles = 0.5, paral = TRUE, ncores = NC)
out$Nm_isolation_sims_at_or_above_obs <- mean(pz$med >= eN$median)
# broad-prior sensitivity
eB <- estimate(S8, two8, S8$d$p_T_years[two8], tag = "broad"); out$split_broad_prior <- summ(eB)

write(toJSON(out, auto_unbox = TRUE, digits = 6, pretty = TRUE), "../figdata/abcrf_results.json")
cat("done\n")
