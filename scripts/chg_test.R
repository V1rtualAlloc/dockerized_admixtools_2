library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# Direct CHG test: does the target have Caucasus HG ancestry beyond what
# Yamnaya (which is ~50% CHG) already accounts for?
#
# Compares three qpAdm models:
#   A — Yamnaya + WHG + EEF         (reference, bundles EHG+CHG)
#   B — EHG + CHG + WHG + EEF       (decomposed, 4-source)
#   C — EHG + CHG + EEF             (decomposed, no separate WHG)

OUTGROUPS <- MODELS[[BACKGROUND]]$outgroup

all_pops <- unique(c(
  TARGET,
  "Russia_Samara_EBA_Yamnaya",
  "Russia_Karelia_Mesolithic_HG",
  "Georgia_KotiasKlde_Mesolithic",
  "Luxembourg_Loschbour_Mesolithic",
  "Turkey_N",
  OUTGROUPS
))

cat("=== Direct CHG test ===\n")
cat("Target:", TARGET, "\n\n")
cat("Question: is CHG ancestry entirely mediated through Yamnaya,\n")
cat("or does the target carry additional direct CHG signal?\n\n")

f2 <- f2_from_precomp(F2_DIR, pops = all_pops)

run_model <- function(sources, label) {
  r <- tryCatch(
    qpadm(f2, target = TARGET, left = sources, right = OUTGROUPS),
    error = function(e) { cat("  [error:", conditionMessage(e), "]\n\n"); NULL }
  )
  if (is.null(r)) return(invisible(NULL))

  full <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  p    <- full$p[1]

  w_col <- intersect(c("weight", "w", "est"), names(r$weights))
  wts   <- setNames(r$weights[[w_col[1]]], r$weights$left)
  feasible <- all(wts >= -0.01 & wts <= 1.01)

  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  p = %.4f  |  %s  |  weights %s\n\n",
              p,
              if (p > 0.05) "PASS" else "fail",
              if (feasible) "feasible" else "INFEASIBLE"))
  for (nm in names(wts)) {
    cat(sprintf("  %-42s  %+.1f%%\n", nm, wts[nm] * 100))
  }
  cat("\n")
  invisible(list(p = p, wts = wts, feasible = feasible))
}

rA <- run_model(
  c("Russia_Samara_EBA_Yamnaya", "Luxembourg_Loschbour_Mesolithic", "Turkey_N"),
  "Model A — Yamnaya + WHG + EEF  [reference]"
)

rB <- run_model(
  c("Russia_Karelia_Mesolithic_HG", "Georgia_KotiasKlde_Mesolithic",
    "Luxembourg_Loschbour_Mesolithic", "Turkey_N"),
  "Model B — EHG + CHG + WHG + EEF  [4-source decomposed]"
)

rC <- run_model(
  c("Russia_Karelia_Mesolithic_HG", "Georgia_KotiasKlde_Mesolithic", "Turkey_N"),
  "Model C — EHG + CHG + EEF  [no separate WHG]"
)

cat("=== Interpretation guide ===\n\n")
cat("Model A: standard Yamnaya-based model. Yamnaya bundles EHG (~50%) + CHG (~50%),\n")
cat("  so this model cannot distinguish a direct CHG signal from a Yamnaya one.\n\n")
cat("Model B: if PASS + feasible, the data can be explained with EHG and CHG as\n")
cat("  separate sources — check whether both weights are substantially positive.\n")
cat("  Infeasible weights (one going negative) means EHG and CHG are too collinear\n")
cat("  to separate given the data, so Yamnaya remains the better proxy.\n\n")
cat("Model C: if PASS, WHG ancestry is already captured by EHG (which is itself\n")
cat("  ~50% WHG + ~50% ANE), so a separate WHG source is not needed.\n\n")
cat("Done.\n")
