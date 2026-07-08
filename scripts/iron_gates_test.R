#!/usr/bin/env Rscript
# Compare Loschbour vs Iron Gates as the HG proxy in the 3-source qpAdm model.

library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

outgroup <- MODELS$european$outgroup

models <- list(
  loschbour    = c("Russia_Samara_EBA_Yamnaya", "Luxembourg_Loschbour_Mesolithic", "Turkey_N"),
  iron_gates   = c("Russia_Samara_EBA_Yamnaya", "Serbia_IronGates_Mesolithic",     "Turkey_N"),
  ehg          = c("Russia_Samara_EBA_Yamnaya", "Russia_Karelia_Mesolithic_HG",    "Turkey_N")
)

all_pops <- c(
  TARGET,
  MODELS$european$outgroup,
  "Russia_Samara_EBA_Yamnaya",
  "Luxembourg_Loschbour_Mesolithic",
  "Russia_Karelia_Mesolithic_HG",
  "Serbia_IronGates_Mesolithic",
  "Turkey_N"
)

IG_F2_DIR <- sub("f2/?$", "f2_ig/", F2_DIR)

if (length(list.files(IG_F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks for Iron Gates test (first run, ~2 min)...\n")
  extract_f2(MERGED_PREFIX, outdir = IG_F2_DIR, pops = all_pops, overwrite = TRUE)
} else {
  cat("Loading cached f2 blocks...\n")
}
f2_data <- f2_from_precomp(IG_F2_DIR)

results <- imap(models, function(sources, name) {
  res <- qpadm(f2_data, sources, outgroup, TARGET, verbose = FALSE)
  w   <- res$weights
  p <- res$popdrop %>% filter(!grepl("1", pat)) %>% pull(p)
  list(name=name, sources=sources, weights=w, p=if (length(p)) p[1] else NA)
})

cat("\n=== qpAdm model comparison: HG proxy ===\n\n")

for (r in results) {
  cat(sprintf("Model: %s\n", r$name))
  cat(sprintf("  Sources: %s\n", paste(r$sources, collapse=" + ")))
  cat(sprintf("  p-value: %.4f  %s\n", r$p,
              ifelse(r$p > 0.05, "(passes)", "(FAILS)")))
  wt <- r$weights %>% select(any_of(c("left","source","target","weight","se")))
  for (i in seq_len(nrow(wt))) {
    src <- if ("left" %in% names(wt)) wt$left[i] else wt$source[i]
    cat(sprintf("  %s: %.1f%% (SE %.1f%%)\n",
                src, wt$weight[i]*100, wt$se[i]*100))
  }
  cat("\n")
}
