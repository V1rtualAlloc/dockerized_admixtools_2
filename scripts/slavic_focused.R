library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

OUTGROUPS <- SLAVIC_OUTGROUPS

SLAVIC <- c("Poland_EarlyMedieval_Slav", "Poland_EarlySlav")

BALKAN_IA <- c(
  "NorthMacedonia_IA",          # N=15, central Balkans
  "Croatia_EIA",                # N=24, western Balkans
  "Bulgaria_KapitanAndreevo_EIA", # N=11, eastern Balkans
  "Serbia_ImperialRoman"        # N=30, Serbian Roman (bundled substrate)
)

ROMAN <- c("Turkey_Medieval_Byzantine", "Turkey_LateAntiquity_Byzantine",
           "Turkey_EarlyMedieval_Byzantine", "Italy_Lazio_ImperialRoman_Roman")

pops_needed <- unique(c(TARGET, SLAVIC, BALKAN_IA, ROMAN, OUTGROUPS))
cat("Populations in cache:", length(pops_needed), "\n")
cat(paste(" ", pops_needed, collapse="\n"), "\n\n")

FOCUSED_F2_DIR <- sub("f2/?$", "f2_focused/", F2_DIR)

cached <- list.dirs(FOCUSED_F2_DIR, full.names = FALSE, recursive = FALSE)
if (!all(pops_needed %in% cached)) {
  cat("Building f2 cache...\n")
  dir.create(FOCUSED_F2_DIR, showWarnings = FALSE, recursive = TRUE)
  extract_f2(MERGED_PREFIX, outdir = FOCUSED_F2_DIR, pops = pops_needed, overwrite = TRUE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(FOCUSED_F2_DIR, pops = pops_needed)

run <- function(sources, label) {
  r <- tryCatch(qpadm(f2, target=TARGET, left=sources, right=OUTGROUPS), error=function(e) NULL)
  if (is.null(r)) { cat(label, "ERROR\n\n"); return(invisible(NULL)) }
  full <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  w_col  <- intersect(c("weight","w","est"), names(r$weights))
  se_col <- intersect(c("se","se_weight"),   names(r$weights))
  wts <- setNames(r$weights[[w_col[1]]], r$weights$left)
  ses <- setNames(r$weights[[se_col[1]]], r$weights$left)
  feasible <- all(wts >= -0.01 & wts <= 1.01)
  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  p=%.3f  %s  %s\n", full$p[1],
              if (full$p[1] > 0.05) "PASS" else "fail",
              if (feasible) "feasible" else "INFEASIBLE"))
  for (nm in names(wts))
    cat(sprintf("  %-46s  %+6.1f%%  ±%5.1f%%  z=%+.2f\n", nm, wts[nm]*100, ses[nm]*100, wts[nm]/ses[nm]))
  cat("\n")
}

cat("=== 2-source: Slavic + Balkan IA ===\n\n")
for (sl in SLAVIC) for (ia in BALKAN_IA)
  run(c(sl, ia), sprintf("%s + %s", sl, ia))

for (ro in ROMAN) {
  cat(sprintf("=== 3-source: Slavic + Balkan IA + %s ===\n\n", ro))
  for (sl in SLAVIC) for (ia in BALKAN_IA)
    run(c(sl, ia, ro), sprintf("%s + %s", sl, ia))
}
