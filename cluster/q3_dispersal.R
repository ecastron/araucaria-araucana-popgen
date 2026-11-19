#!/usr/bin/env Rscript
# Q3: which dispersal process generates the spatial genetic structure found in Q2?
#
# Three candidate agents predict kin structure at distinguishable spatial scales:
#
#   Austral parakeet (Enicognathus ferrugineus)  ~100 m   long-distance, measured 10-153 m
#                                                         (mean 61 m) in a congener
#   Scatter-hoarding rodents (Abrothrix spp.)    ~10 m    caches, some beyond 40 m
#   Failed dispersal                             ~0 m     seedlings directly beneath mothers
#
# Q2 established that kinship declines with distance (Sp = 0.0071, p = 1.6e-10). Q3 asks a
# harder question: at WHAT distance does kinship reach background? That is a property of curve
# SHAPE, not slope, and shape is what degrades first under noise.
#
# PRE-DECLARED POWER CHECK, evaluated before any interpretation:
#   To separate ~10 m from ~100 m we need enough pairs in each of the fine bins. If fewer than
#   30 within-stand pairs fall below 25 m, the fine-scale bins cannot support a mean and Q3 is
#   NOT ANSWERABLE. That verdict is an acceptable outcome, declared in advance so it cannot be
#   rationalised away afterwards.
#
# The extent of detectable kin structure is estimated as the first distance bin whose 95%
# confidence interval includes the background kinship (the mean over pairs beyond 1 km).

set.seed(1)
P <- "/home/ecastron/araucaria"; OUT <- file.path(P, "results/Q3")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
sink(file.path(OUT, "q3_summary.txt"), split = TRUE)

x <- readRDS(file.path(P, "results/Q2/q2_sgs.rds"))
d <- x$d; k <- x$k; st <- x$stand
cat(sprintf("Q3 dispersal mechanism\n%d within-stand pairs, %d stands, %.1f-%.1f m\n\n",
            length(d), length(unique(st)), min(d), max(d)))

# ---- power check ---------------------------------------------------------------------------
cat("--- pair counts at the scales that separate the mechanisms ---\n")
brk <- c(0, 10, 25, 50, 100, 250, 500, 1000, Inf)
lab <- c("0-10 m", "10-25 m", "25-50 m", "50-100 m", "100-250 m", "250-500 m",
         "500-1000 m", ">1000 m")
cl <- cut(d, brk, labels = lab, include.lowest = TRUE)
cnt <- table(cl)
for (i in seq_along(lab)) cat(sprintf("  %-11s %5d pairs\n", lab[i], cnt[i]))
fine <- sum(d < 25)
cat(sprintf("\npairs below 25 m: %d\n", fine))
answerable <- fine >= 30
if (!answerable) {
  cat("\n>>> PRE-DECLARED VERDICT: Q3 IS NOT ANSWERABLE.\n")
  cat(">>> Fewer than 30 pairs below 25 m, so the bins that would separate rodent caching\n")
  cat(">>> (~10 m) from parakeet dispersal (~100 m) cannot support a mean. This is a\n")
  cat(">>> limitation of the sampling design, NOT evidence about the dispersal agent, and it\n")
  cat(">>> must not be reported as though it favoured any mechanism.\n\n")
}

# ---- kinship by fine distance bin ----------------------------------------------------------
cat("--- mean kinship by distance bin ---\n")
cat(sprintf("  %-11s %6s %9s %9s %20s\n", "bin", "n", "mean", "se", "95% CI"))
res <- data.frame()
for (i in seq_along(lab)) {
  s <- k[cl == lab[i]]
  s <- s[is.finite(s)]
  if (length(s) < 2) { cat(sprintf("  %-11s %6d        --\n", lab[i], length(s))); next }
  m <- mean(s); se <- sd(s) / sqrt(length(s))
  cat(sprintf("  %-11s %6d %9.4f %9.4f   [%7.4f, %7.4f]\n",
              lab[i], length(s), m, se, m - 1.96 * se, m + 1.96 * se))
  res <- rbind(res, data.frame(bin = lab[i], n = length(s), mean = m, se = se,
                               lo = m - 1.96 * se, hi = m + 1.96 * se))
}

bg <- k[d > 1000]; bg <- bg[is.finite(bg)]
bgm <- mean(bg)
cat(sprintf("\nbackground kinship (pairs > 1 km): %.4f  (n = %d)\n", bgm, length(bg)))

# ---- extent of detectable kin structure ----------------------------------------------------
above <- res[res$lo > bgm, ]
cat("\n--- extent of detectable kin structure ---\n")
if (nrow(above) == 0) {
  cat("  NO distance bin has kinship significantly above background.\n")
  cat("  The Q2 slope is detectable across the whole range without any bin standing out,\n")
  cat("  which is consistent with weak, diffuse structure rather than a sharp kin cluster.\n")
} else {
  cat(sprintf("  bins above background: %s\n", paste(above$bin, collapse = ", ")))
  cat(sprintf("  structure detectable out to roughly %s\n", above$bin[nrow(above)]))
}

# ---- comparison against the three predictions ----------------------------------------------
cat("\n--- can the three mechanisms be separated? ---\n")
pred <- data.frame(agent = c("failed dispersal", "rodent caching", "parakeet"),
                   scale_m = c(2, 10, 100))
for (i in seq_len(nrow(pred))) {
  n_at <- sum(d <= pred$scale_m[i])
  cat(sprintf("  %-18s predicted ~%3d m : %4d pairs at or below that distance\n",
              pred$agent[i], pred$scale_m[i], n_at))
}
cat("\n")
if (!answerable) {
  cat("VERDICT: cannot distinguish. The design samples 10-20 trees per stand over hundreds of\n")
  cat("metres, so the short-distance bins where these mechanisms differ are nearly empty.\n")
  cat("Q3 requires a sampling design targeting near neighbours, not a reanalysis of these data.\n")
} else {
  cat("VERDICT: fine-scale bins are populated; read the extent above against the predictions.\n")
}

write.csv(res, file.path(OUT, "q3_kinship_bins.csv"), row.names = FALSE)
saveRDS(list(res = res, background = bgm, answerable = answerable, n_below_25 = fine),
        file.path(OUT, "q3.rds"))
cat(sprintf("\nwritten: %s\n", file.path(OUT, "q3_kinship_bins.csv")))
sink()
