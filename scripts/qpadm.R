library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

model   <- MODELS[[BACKGROUND]]
sources <- model$sources
right   <- model$outgroup
ROTATE_POOL <- c(
  "Russia_Samara_EBA_Yamnaya", "Serbia_IronGates_Mesolithic",
  "Luxembourg_Loschbour_Mesolithic", "Turkey_N",
  "Russia_Karelia_Mesolithic_HG", "Georgia_KotiasKlde_Mesolithic"
)
all_pops <- c(TARGET, REFERENCES, sources, right, ROTATE_POOL)

cat("=== qpAdm ancestry analysis ===\n")
cat("Target:    ", TARGET, "\n")
cat("Background:", BACKGROUND, "\n")
cat("Sources:   ", paste(sources, collapse = ", "), "\n\n")


# Load or compute f2 blocks
if (length(list.files(F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks (first run — may take several minutes)...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_DIR, pops = all_pops, overwrite = FALSE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(F2_DIR)
cat("f2 blocks loaded.\n\n")


# 3-source model for target
cat("=== 3-source model:", TARGET, "===\n\n")
result <- qpadm(f2, left = sources, right = right, target = TARGET)

cat("--- Ancestry proportions ---\n")
print(result$weights %>% select(target, left, weight, se, z))

cat("\n--- Model fit ---\n")
print(result$popdrop %>% select(pat, p, feasible) %>% head(10))


# 2-source sub-models
cat("\n=== 2-source sub-models ===\n\n")
for (pair in combn(sources, 2, simplify = FALSE)) {
  r <- tryCatch(qpadm(f2, left = pair, right = right, target = TARGET), error = function(e) NULL)
  if (is.null(r)) next
  full <- r$popdrop %>% filter(!grepl("1", pat))
  pv   <- if (nrow(full) > 0) full$p[1] else NA
  w    <- r$weights %>% filter(left == pair[1]) %>% pull(weight)
  cat(sprintf("  %-45s | weight1 = %5.3f | p = %s | %s\n",
              paste(pair, collapse = " + "),
              ifelse(length(w) > 0, w[1], NA),
              ifelse(!is.na(pv), sprintf("%.4f", pv), "NA"),
              ifelse(!is.na(pv) && pv > 0.05, "PASS", "fail")))
}


# Comparison table across reference populations
cat("\n=== Comparison table ===\n\n")

run_qpadm <- function(target) {
  r <- tryCatch(qpadm(f2, left = sources, right = right, target = target), error = function(e) NULL)
  if (is.null(r)) return(NULL)
  w  <- r$weights %>% select(left, weight) %>% deframe()
  pv <- r$popdrop %>% filter(!grepl("1", pat)) %>% pull(p)
  pv <- if (length(pv) > 0) pv[1] else NA
  tibble(
    population = target,
    !!sources[1] := w[sources[1]],
    !!sources[2] := w[sources[2]],
    !!sources[3] := w[sources[3]],
    p          = pv,
    fit        = ifelse(!is.na(pv) & pv > 0.05, "PASS", "fail")
  )
}

comparison <- bind_rows(run_qpadm(TARGET), map(REFERENCES, run_qpadm)) %>%
  arrange(desc(.data[[sources[1]]]))

cat(sprintf("  %-38s  %13s  %13s  %13s  %8s  %s\n",
            "Population", sources[1], sources[2], sources[3], "p-value", "fit"))
cat(strrep("-", 100), "\n")
for (i in seq_len(nrow(comparison))) {
  r <- comparison[i, ]
  marker <- if (r$population == TARGET) " <-- you" else ""
  cat(sprintf("  %-38s  %12.1f%%  %12.1f%%  %12.1f%%  %8s  %s%s\n",
              r$population,
              r[[sources[1]]] * 100, r[[sources[2]]] * 100, r[[sources[3]]] * 100,
              ifelse(is.na(r$p), "NA", sprintf("%.4f", r$p)),
              r$fit, marker))
}

cat("\nDone.\n")
