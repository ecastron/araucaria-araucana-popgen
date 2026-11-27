#!/usr/bin/env Rscript
# Revision (25 Nov 2026): diagnostics for the ABC-RF analysis.
#
# Input : abcrf/reftable_sim9.csv, observed_sim9.csv (repository; same tables as abcrf_analysis.R)
# Steps (abcrf 2.0; same settings as abcrf_analysis.R: 1,000 trees, LDA axes added)
#   1. Variable importance of the forests for one vs two refugia and isolation vs contact (seed 101,
#      the first of the ten replicate forests in abcrf_analysis.R)
#   2. Posterior calibration (BESPOKE; abcrf has no calibration function): hold out 20% of the
#      simulations (seed 21, stratified by scenario), train on the rest, compute the posterior
#      probability of the selected model for each held-out simulation with predict.abcrf, bin it,
#      and compare the mean posterior with the share of correct selections in each bin
#   3. Prior probabilities of the split-age events, from the H2 + H3 simulations, next to the
#      posterior tail areas already in figdata/abcrf_results.json; Bayes factors
#   4. Zero-migration calibration in full: distribution of the migrants-per-generation estimate
#      (regAbcrf trained on the contact simulations, seed 7 as abcrf_analysis.R) applied to all
#      isolation simulations (true migration zero)
# Output: figdata/revision/abcrf_diagnostics.json, abcrf_importance.tsv, abcrf_calibration.tsv (same folder)
suppressPackageStartupMessages({library(abcrf); library(jsonlite)})
here <- dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE))))
R <- dirname(here)   # repository root
NT <- 1000; NC <- max(1, parallel::detectCores() - 2)
d <- read.csv(file.path(R, "abcrf/reftable_sim9.csv"))
o <- read.csv(file.path(R, "abcrf/observed_sim9.csv"), row.names = 1)
stats <- setdiff(names(d), c("model", grep("^p_", names(d), value = TRUE)))
d$model <- factor(d$model, levels = c("H1", "H1b", "H2", "H3"))
res <- list(); out <- file.path(R, "figdata/revision"); dir.create(out, showWarnings = FALSE, recursive = TRUE)
G1 <- list(c("H1", "H1b"), c("H2", "H3"))
two <- which(d$model %in% c("H2", "H3"))

# ---- 1. variable importance ----
set.seed(101)
dd1 <- data.frame(model = d$model, d[, stats])
rf1 <- abcrf(model ~ ., data = dd1, group = G1, ntree = NT, paral = TRUE, ncores = NC)
set.seed(101)
dd2 <- data.frame(model = droplevels(d$model[two]), d[two, stats])
rf2 <- abcrf(model ~ ., data = dd2, ntree = NT, paral = TRUE, ncores = NC)
imp <- function(rf) { v <- rf$model.rf$variable.importance; v <- v / sum(v); sort(v, decreasing = TRUE) }
i1 <- imp(rf1); i2 <- imp(rf2)
write.table(rbind(data.frame(comparison = "one_vs_two_refugia", statistic = names(i1), share = round(i1, 4)),
                  data.frame(comparison = "isolation_vs_contact", statistic = names(i2), share = round(i2, 4))),
            file.path(out, "abcrf_importance.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
res$importance_one_vs_two <- as.list(round(i1, 4)); res$importance_iso_vs_contact <- as.list(round(i2, 4))

# ---- 2. calibration on held-out simulations ----
calib <- function(rows, groups, label) {
  set.seed(21)
  hold <- unlist(lapply(split(rows, d$model[rows]), function(r) sample(r, round(0.2 * length(r)))))
  trn <- setdiff(rows, hold)
  ddt <- data.frame(model = droplevels(d$model[trn]), d[trn, stats])
  rf <- abcrf(model ~ ., data = ddt, group = groups, ntree = NT, paral = TRUE, ncores = NC)
  p <- predict(rf, obs = d[hold, stats], training = ddt, ntree = NT, paral = TRUE, ncores = NC)
  truth <- as.character(d$model[hold])
  if (length(groups)) truth <- ifelse(truth %in% groups[[1]], "g1", "g2")
  sel <- as.character(p$allocation); pp <- as.numeric(p$post.prob)
  br <- c(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0001)
  b <- cut(pp, br, right = FALSE)
  tab <- data.frame(comparison = label, bin = levels(b),
                    n = as.integer(table(b)),
                    mean_posterior = round(tapply(pp, b, mean), 3),
                    share_correct = round(tapply(sel == truth, b, mean), 3))
  list(table = tab, n_heldout = length(hold), overall_mean_posterior = mean(pp), overall_correct = mean(sel == truth),
       heldout_ge_0.95 = list(n = sum(pp >= 0.95), share_correct = mean((sel == truth)[pp >= 0.95])))
}
c1 <- calib(seq_len(nrow(d)), G1, "one_vs_two_refugia")
c2 <- calib(two, list(), "isolation_vs_contact")
write.table(rbind(c1$table, c2$table), file.path(out, "abcrf_calibration.tsv"), sep = "\t", row.names = FALSE, quote = FALSE)
res$calibration_one_vs_two <- c1[-1]; res$calibration_iso_vs_contact <- c2[-1]
res$calibration_tables <- rbind(c1$table, c2$table)

# ---- 3. prior vs posterior for split-age events ----
J <- fromJSON(file.path(R, "figdata/abcrf_results.json"))
T <- d$p_T_years[two]
ev <- list(lt_10kyr = list(prior = mean(T < 1e4), post = 1 - J$split_model_averaged$P_older_than_10kyr),
           gt_26.5kyr = list(prior = mean(T > 26.5e3), post = J$split_model_averaged$`P_older_than_26.5kyr`),
           gt_115kyr = list(prior = mean(T > 115e3), post = J$split_model_averaged$P_older_than_115kyr))
odds <- function(p) p / (1 - p)
res$split_events <- lapply(ev, function(e) c(e, bayes_factor = odds(e$post) / odds(e$prior)))
res$split_prior_median_kyr <- median(T) / 1000

# ---- 4. migration estimate when true migration is zero ----
con <- which(d$model == "H3"); iso <- which(d$model == "H2")
NmC <- d$p_Ne[con] * d$p_mcross[con]
set.seed(7)
ddn <- data.frame(y = NmC, d[con, stats])
rn <- regAbcrf(y ~ ., data = ddn, ntree = NT, paral = TRUE, ncores = NC)
pn <- predict(rn, obs = o["primary", stats], training = ddn, quantiles = c(0.025, 0.975), paral = TRUE, ncores = NC)
pz <- predict(rn, obs = d[iso, stats], training = ddn, quantiles = 0.5, paral = TRUE, ncores = NC)
z <- as.numeric(pz$med)
res$Nm_zero_migration <- list(observed_median = as.numeric(pn$med),
  isolation_estimates_quantiles = as.list(setNames(quantile(z, c(0.05, 0.25, 0.5, 0.75, 0.95)), c("q05", "q25", "q50", "q75", "q95"))),
  share_isolation_at_or_above_observed = mean(z >= as.numeric(pn$med)),
  prior_median_contact = median(NmC))
write(toJSON(res, auto_unbox = TRUE, digits = 5, pretty = TRUE), file.path(out, "abcrf_diagnostics.json"))
print(res$calibration_tables); str(res[c("split_events", "Nm_zero_migration")])
cat("\ntop importance, one vs two:\n"); print(head(round(i1, 3), 8))
cat("\ntop importance, iso vs contact:\n"); print(head(round(i2, 3), 8))
