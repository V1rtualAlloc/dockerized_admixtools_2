library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

# Pools three Iron Age Balkan populations into a single "Balkans_IA" label to
# increase effective N and reduce SEs. The three populations share Illyrian/
# Hallstatt heritage and are genetically close; pooling is a defensible shortcut
# when individual proxy SEs are too wide.
#
#   Croatia_EIA (N=24) + NorthMacedonia_IA (N=15) + Bulgaria_KapitanAndreevo_EIA (N=11)
#   → Balkans_IA (N=50) — optimal pool; adding smaller populations increases missingness and raises SE
#
# Creates me/merged_pooled.ind with remapped labels; symlinks .geno/.snp.

OUTGROUPS <- SLAVIC_OUTGROUPS

POOL_SOURCES <- SLAVIC_BALKAN_IA_CORE
POOLED_LABEL <- "Balkans_IA"

ROMAN <- SLAVIC_ROMAN_MINIMAL

# ── Step 1: build pooled .ind + symlinks ──────────────────────────────────────

ind_path    <- paste0(MERGED_PREFIX, ".ind")
pooled_base <- file.path(dirname(MERGED_PREFIX), "merged_pooled")

ind <- read.table(ind_path, header=FALSE, col.names=c("id","sex","pop"),
                  stringsAsFactors=FALSE)

n_before <- sum(ind$pop %in% POOL_SOURCES)
ind$pop[ind$pop %in% POOL_SOURCES] <- POOLED_LABEL
n_after  <- sum(ind$pop == POOLED_LABEL)

cat(sprintf("Pooled %d individuals into '%s' (from: %s)\n\n",
            n_after, POOLED_LABEL, paste(POOL_SOURCES, collapse=", ")))

write.table(ind, paste0(pooled_base, ".ind"),
            quote=FALSE, row.names=FALSE, col.names=FALSE, sep="\t")

for (ext in c(".geno", ".snp")) {
  src <- paste0(MERGED_PREFIX, ext)
  dst <- paste0(pooled_base, ext)
  if (!file.exists(dst)) file.symlink(src, dst)
}

# ── Step 2: focused f2 cache ─────────────────────────────────────────────────

all_pops <- unique(c(TARGET, "Poland_IA", POOLED_LABEL, ROMAN, OUTGROUPS))

F2_POOLED_DIR <- file.path(dirname(F2_DIR), "f2_pooled/")

cached <- list.dirs(F2_POOLED_DIR, full.names=FALSE, recursive=FALSE)
if (!all(all_pops %in% cached)) {
  cat("Building pooled f2 cache...\n\n")
  dir.create(F2_POOLED_DIR, showWarnings=FALSE, recursive=TRUE)
  extract_f2(pooled_base, outdir=F2_POOLED_DIR, pops=all_pops, overwrite=TRUE)
} else {
  cat("Pooled f2 cache loaded.\n\n")
}
f2 <- f2_from_precomp(F2_POOLED_DIR, pops=all_pops)

# ── Step 3: models ────────────────────────────────────────────────────────────

run_model <- function(sources, label) {
  r <- tryCatch(
    qpadm(f2, target=TARGET, left=sources, right=OUTGROUPS),
    error=function(e) { cat("  [error:", conditionMessage(e), "]\n\n"); NULL }
  )
  if (is.null(r)) return(invisible(NULL))

  full     <- r$popdrop %>% filter(!grepl("1", pat)) %>% slice(1)
  p        <- full$p[1]
  w_col    <- intersect(c("weight","w","est"), names(r$weights))
  se_col   <- intersect(c("se","se_weight"),   names(r$weights))
  wts      <- setNames(r$weights[[w_col[1]]], r$weights$left)
  ses      <- setNames(r$weights[[se_col[1]]], r$weights$left)
  feasible <- all(wts >= -0.01 & wts <= 1.01)

  cat(sprintf("--- %s ---\n", label))
  cat(sprintf("  p = %.4f  |  %s  |  %s\n\n",
              p,
              if (p > 0.05) "PASS" else "fail",
              if (feasible) "feasible" else "INFEASIBLE"))
  cat(sprintf("  %-46s  %8s  %7s  %6s\n", "source", "weight", "SE", "z"))
  for (nm in names(wts))
    cat(sprintf("  %-46s  %+7.1f%%  %6.1f%%  %+6.2f\n",
                nm, wts[nm]*100, ses[nm]*100, wts[nm]/ses[nm]))
  cat("\n")
  invisible(list(p=p, wts=wts, ses=ses, feasible=feasible))
}

cat("=== Pooled Balkans IA (N=50) — SE comparison ===\n\n")
cat(sprintf("Target: %s\n", TARGET))
cat(sprintf("Pooled: %s → %s (N=%d)\n\n",
            paste(POOL_SOURCES, collapse=" + "), POOLED_LABEL, n_after))

run_model(
  c("Poland_IA", POOLED_LABEL),
  "2-source: Slavic + Balkans_IA [pooled]"
)
run_model(
  c("Poland_IA", POOLED_LABEL, "Italy_Lazio_ImperialRoman_Roman"),
  "3-source: Slavic + Balkans_IA [pooled] + Roman Italian"
)
run_model(
  c("Poland_IA", POOLED_LABEL, "Turkey_Medieval_Byzantine"),
  "3-source: Slavic + Balkans_IA [pooled] + Byzantine"
)

cat("Done.\n")
