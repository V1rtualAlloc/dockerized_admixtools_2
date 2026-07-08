library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# Outgroup f3: f3(Mbuti; TARGET, X) for a wide set of populations X.
# Higher f3 = more shared drift between TARGET and X = closer relationship.
# This is an unbiased way to rank all populations by affinity to the target.

OUTGROUP <- "Mbuti"

# Populations to test against — broad set covering ancient and modern diversity.
# Draws on the background's own source/outgroup/pool populations (already
# verified for BACKGROUND) plus REFERENCES (region-specific, see config.R)
# and WORLD_REFS (background-agnostic modern diversity panel).
model <- MODELS[[BACKGROUND]]
TEST_POPS <- unique(c(
  REFERENCES,
  if (!is.null(model$pool)) model$pool else model$sources,
  model$outgroup,
  WORLD_REFS
))
TEST_POPS <- setdiff(TEST_POPS, c(TARGET, OUTGROUP))

all_pops <- c(TARGET, OUTGROUP, TEST_POPS)

cat("=== Outgroup f3 analysis ===\n")
cat("f3(", OUTGROUP, ";", TARGET, ", X)\n\n")

# Load or compute f2 blocks
if (length(list.files(F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks (first run — may take several minutes)...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_DIR, pops = all_pops, overwrite = FALSE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(F2_DIR)

results <- f3(f2, pop1 = OUTGROUP, pop2 = TARGET, pop3 = TEST_POPS) %>%
  rename(pop = pop3) %>%
  arrange(desc(est))

cat(sprintf("  %-38s  %10s  %8s  %6s\n", "Population", "f3", "se", "z"))
cat(strrep("-", 70), "\n")
for (i in seq_len(nrow(results))) {
  r <- results[i, ]
  cat(sprintf("  %-38s  %10.6f  %8.6f  %6.2f\n", r$pop, r$est, r$se, r$z))
}

cat("\nDone.\n")
