library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# Outgroup-power test: does adding well-covered deep Steppe/WHG populations
# to the outgroup set sharpen (lower SE / raise power) the Slavic vs.
# pre-Slavic Balkan-IA decomposition, without the missingness pollution seen
# when Russia_Samara_EBA_Yamnaya + Luxembourg_Loschbour_Mesolithic (N=2) were
# added directly to SLAVIC_OUTGROUPS in scripts/config.R?
#
# Fix vs. that earlier attempt: swap the N=2 Loschbour genome for
# Serbia_IronGates_Mesolithic (N=73, already this project's established WHG
# proxy in the deep-ancestry model). Russia_Samara_EBA_Yamnaya (N=46) is kept
# — it was never actually low-coverage, Loschbour was the weak link.
#
# Runs in its own dedicated f2 cache (me/f2_slavic_deepout/) so this can't
# contaminate me/f2_slavic/ used by slavic_model.R.

BASE_OUTGROUPS <- SLAVIC_OUTGROUPS
DEEP_OUTGROUPS <- c(SLAVIC_OUTGROUPS, "Russia_Samara_EBA_Yamnaya", "Serbia_IronGates_Mesolithic")

SOURCES <- SLAVIC_BEST_SOURCES

TEST_F2_DIR <- sub("f2/?$", "f2_slavic_deepout/", F2_DIR)

all_pops <- unique(c(TARGET, SOURCES, DEEP_OUTGROUPS))

cat("=== Outgroup-power test: baseline vs. deep-outgroup Slavic model ===\n")
cat("Target:", TARGET, "\n")
cat("Sources:", paste(SOURCES, collapse = " + "), "\n\n")

cached <- list.dirs(TEST_F2_DIR, full.names = FALSE, recursive = FALSE)
if (!all(all_pops %in% cached)) {
  cat("Building dedicated f2 cache (may take a few minutes)...\n\n")
  dir.create(TEST_F2_DIR, showWarnings = FALSE, recursive = TRUE)
  extract_f2(MERGED_PREFIX, outdir = TEST_F2_DIR, pops = all_pops, overwrite = TRUE)
} else {
  cat("f2 blocks cached, loading...\n\n")
}
f2 <- f2_from_precomp(TEST_F2_DIR, pops = all_pops)

run_model <- function(outgroups, label) {
  r <- tryCatch(
    qpadm(f2, target = TARGET, left = SOURCES, right = outgroups),
    error = function(e) { cat("  [error:", conditionMessage(e), "]\n\n"); NULL }
  )
  if (is.null(r)) return(invisible(NULL))

  full <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  p    <- full$p[1]

  w_col <- intersect(c("weight", "w", "est"), names(r$weights))
  se_col <- intersect(c("se", "se_weight"), names(r$weights))
  z_col  <- intersect(c("z", "z_weight"), names(r$weights))
  wts <- setNames(r$weights[[w_col[1]]], r$weights$left)
  ses <- if (length(se_col) > 0) setNames(r$weights[[se_col[1]]], r$weights$left) else rep(NA, length(wts))
  zs  <- if (length(z_col)  > 0) setNames(r$weights[[z_col[1]]],  r$weights$left) else wts / ses
  feasible <- all(wts >= -0.01 & wts <= 1.01)

  cat(sprintf("--- %s (%d outgroups) ---\n", label, length(outgroups)))
  cat(sprintf("  p = %.4f  |  %s  |  weights %s\n\n",
              p,
              if (p > 0.05) "PASS" else "fail",
              if (feasible) "feasible" else "INFEASIBLE"))
  cat(sprintf("  %-30s  %8s  %8s  %6s\n", "source", "weight", "se", "z"))
  for (nm in names(wts)) {
    cat(sprintf("  %-30s  %+7.1f%%  %7.1f%%  %+6.2f\n",
                nm, wts[nm] * 100, ses[nm] * 100, zs[nm]))
  }
  cat("\n")
  invisible(list(p = p, wts = wts, ses = ses, feasible = feasible))
}

base_res <- run_model(BASE_OUTGROUPS, "Baseline outgroups")
deep_res <- run_model(DEEP_OUTGROUPS, "Baseline + Yamnaya(N=46) + IronGates_Mesolithic(N=73)")

if (!is.null(base_res) && !is.null(deep_res)) {
  cat("=== Comparison ===\n\n")
  cat(sprintf("  %-30s  %12s  %12s\n", "", "baseline SE", "deep-outgroup SE"))
  for (nm in names(base_res$ses)) {
    cat(sprintf("  %-30s  %10.1f%%  %14.1f%%\n", nm, base_res$ses[nm] * 100, deep_res$ses[nm] * 100))
  }
  cat(sprintf("\n  p-value: baseline = %.4f, deep-outgroup = %.4f\n", base_res$p, deep_res$p))
}

cat("\nDone.\n")
