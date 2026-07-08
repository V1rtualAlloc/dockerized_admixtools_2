library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# qpAdm rotating: systematically test all 2-source and 3-source combinations
# from a pool of ancient populations. Reports which models pass (p > 0.05)
# with feasible weights (all between 0 and 1).

SOURCE_POOL <- c(
  "Russia_Samara_EBA_Yamnaya",       # Steppe
  "Serbia_IronGates_Mesolithic",      # Balkan HG (replaces WHG/Loschbour)
  "Luxembourg_Loschbour_Mesolithic",  # WHG (kept for comparison)
  "Turkey_N",                         # EEF
  "Russia_Karelia_Mesolithic_HG",     # EHG
  "Georgia_KotiasKlde_Mesolithic"     # CHG
  # Iran_GanjDareh_N and Israel_Natufian excluded: they are in the outgroup set
)

OUTGROUPS <- MODELS[[BACKGROUND]]$outgroup
all_pops  <- unique(c(TARGET, SOURCE_POOL, OUTGROUPS))

cat("=== qpAdm rotating ===\n")
cat("Target:", TARGET, "\n")
cat("Source pool:", paste(SOURCE_POOL, collapse = ", "), "\n\n")

if (length(list.files(F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_DIR, pops = all_pops, overwrite = FALSE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(F2_DIR, pops = all_pops)

combos <- c(
  combn(SOURCE_POOL, 2, simplify = FALSE),
  combn(SOURCE_POOL, 3, simplify = FALSE)
)
cat(sprintf("Testing %d models (%d 2-source, %d 3-source)...\n\n",
            length(combos),
            choose(length(SOURCE_POOL), 2),
            choose(length(SOURCE_POOL), 3)))

run_one <- function(sources) {
  r <- tryCatch(
    qpadm(f2, target = TARGET, left = sources, right = OUTGROUPS),
    error = function(e) NULL
  )
  if (is.null(r) || is.null(r$popdrop) || is.null(r$weights)) return(NULL)

  full <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  if (nrow(full) == 0) return(NULL)

  w_col <- intersect(c("weight", "w", "est"), names(r$weights))
  if (length(w_col) == 0) return(NULL)
  wts      <- setNames(r$weights[[w_col[1]]], r$weights$left)
  feasible <- all(wts >= -0.01 & wts <= 1.01)
  w_str    <- paste(sprintf("%s=%.1f%%", names(wts), wts * 100), collapse=", ")

  tibble(n = length(sources), model = paste(sources, collapse=" + "),
         p = full$p, feasible = feasible, weights = w_str)
}

results <- map_dfr(combos, run_one) %>% arrange(desc(p))

passing <- results %>% filter(p > 0.05,  feasible)
failing <- results %>% filter(p <= 0.05 | !feasible)

cat(sprintf("=== PASSING (p > 0.05, feasible weights): %d / %d ===\n\n",
            nrow(passing), nrow(results)))
if (nrow(passing) > 0) {
  cat(sprintf("  %3s  %7s  %-60s  %s\n", "src", "p", "model", "weights"))
  cat(strrep("-", 130), "\n")
  for (i in seq_len(nrow(passing))) {
    r <- passing[i, ]
    cat(sprintf("  %3d  %7.4f  %-60s  %s\n", r$n, r$p, r$model, r$weights))
  }
}

cat(sprintf("\n=== FAILING (p <= 0.05 or infeasible): %d / %d ===\n\n",
            nrow(failing), nrow(results)))
cat(sprintf("  %3s  %7s  %s\n", "src", "p", "model"))
cat(strrep("-", 100), "\n")
for (i in seq_len(nrow(failing))) {
  r <- failing[i, ]
  flag <- if (!r$feasible) " [infeasible weights]" else ""
  cat(sprintf("  %3d  %7.4f  %s%s\n", r$n, r$p, r$model, flag))
}

cat("\nDone.\n")
