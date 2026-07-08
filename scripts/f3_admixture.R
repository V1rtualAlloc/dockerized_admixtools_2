library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# Admixture f3 test: f3(target; A, B) < 0 is evidence of admixture.
# If target is a mixture of A and B, its allele frequencies are intermediate
# between A and B, making the product (p_target - p_A)(p_target - p_B)
# negative on average. z < -3 = significant admixture signal.
# Tests all pairs from a pool of ancient and modern source populations.

SOURCES <- unique(c(
  "Russia_Samara_EBA_Yamnaya",
  "Luxembourg_Loschbour_Mesolithic",
  "Turkey_N",
  "Russia_Karelia_Mesolithic_HG",
  "Georgia_KotiasKlde_Mesolithic",
  "Iran_GanjDareh_N",
  "Israel_Natufian",
  "French", "Sardinian", "Spanish", "Russian",
  REFERENCES
))
SOURCES <- SOURCES[SOURCES != TARGET]

all_pops <- unique(c(TARGET, SOURCES))

cat(sprintf("=== Admixture f3: f3(%s; A, B) ===\n\n", TARGET))
cat("Significantly negative f3 (z < -3) = target is admixed between A and B.\n\n")

if (length(list.files(F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_DIR, pops = all_pops, overwrite = FALSE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(F2_DIR, pops = all_pops)

pairs_mat <- combn(SOURCES, 2)
cat(sprintf("Testing %d population pairs...\n\n", ncol(pairs_mat)))

results <- map_dfr(seq_len(ncol(pairs_mat)), function(i) {
  f3(f2, pop1 = TARGET, pop2 = pairs_mat[1, i], pop3 = pairs_mat[2, i])
}) %>%
  select(source1 = pop2, source2 = pop3, est, se, z) %>%
  arrange(est)

cat(sprintf("  %-35s  %-35s  %10s  %8s  %6s\n",
            "Source 1", "Source 2", "f3", "se", "z"))
cat(strrep("-", 100), "\n")

for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  marker <- if (r$z < -3) " *" else ""
  cat(sprintf("  %-35s  %-35s  %10.6f  %8.6f  %6.2f%s\n",
              r$source1, r$source2, r$est, r$se, r$z, marker))
}

n_sig <- sum(results$z < -3)
cat(sprintf("\n%d of %d pairs show significant admixture signal (z < -3)\n",
            n_sig, nrow(results)))
cat("\nDone.\n")
