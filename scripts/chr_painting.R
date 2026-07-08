library(tidyverse)
source("/data/scripts/config.R")

# Chromosome ancestry painting via allele-frequency maximum likelihood.
# Requires the text EIGENSTRAT subset produced by the prep step:
#   docker run --rm -v $PROJECT_ROOT:/data eigensoft \
#     -p /data/me/subset_extract.par

SUBSET_GENO <- file.path(dirname(MERGED_PREFIX), "subset.geno")
SUBSET_SNP  <- file.path(dirname(MERGED_PREFIX), "subset.snp")
SUBSET_IND  <- file.path(dirname(MERGED_PREFIX), "subset.ind")

if (!file.exists(SUBSET_GENO))
  stop("subset.geno not found. Run convertf first:\n",
       "  docker run --rm -v $PROJECT_ROOT:/data eigensoft ",
       "-p /data/me/subset_extract.par")

model   <- MODELS[[BACKGROUND]]
sources <- model$sources          # exactly 3 sources expected; used as-is for labels
labels  <- sources

snp <- read.table(SUBSET_SNP, col.names = c("rsid","chr","gpos","pos","ref","alt"),
                  colClasses = c("character","integer","numeric","integer","character","character"))
ind <- read.table(SUBSET_IND, col.names = c("id","sex","pop"))

cat("=== Chromosome ancestry painting (allele-frequency MLE) ===\n")
cat("Target:", TARGET, "\n")
cat("Sources:", paste(sources, collapse=", "), "\n")
cat(sprintf("SNPs: %d | Individuals: %d\n\n", nrow(snp), nrow(ind)))

me_idx  <- which(ind$pop == TARGET)
src_idx <- lapply(sources, function(s) which(ind$pop == s))

# Read text EIGENSTRAT genotype matrix [n_ind x n_snp], values 0/1/2/9
cat("Reading genotypes...\n")
lines <- readLines(SUBSET_GENO)
# EIGENSTRAT is SNP-major: one row per SNP, one character per individual.
# nrow = n_snp, ncol = n_ind, byrow = TRUE (read left-to-right across individuals).
geno_mat <- matrix(as.integer(strsplit(paste(lines, collapse=""), "")[[1]]),
                   nrow = nrow(snp), ncol = nrow(ind), byrow = TRUE)
cat(sprintf("Matrix: %d SNPs x %d individuals (%s)\n\n",
            nrow(geno_mat), ncol(geno_mat),
            format(object.size(geno_mat), units="MB")))

# Allele frequency per source pop per SNP (Laplace pseudocount for stability)
# geno_mat is [n_snp x n_ind]: columns = individuals, rows = SNPs
allele_freq <- function(col_idx) {
  g     <- geno_mat[, col_idx, drop = FALSE]  # [n_snp x n_pop_ind]
  valid <- g != 9L
  n_obs <- rowSums(valid) * 2L
  n_alt <- rowSums(g * valid)
  (n_alt + 0.5) / (n_obs + 1.0)
}

cat("Computing source allele frequencies...\n")
p_src <- lapply(src_idx, allele_freq)

g_me <- geno_mat[, me_idx[1]]  # [n_snp] vector, SNP-indexed

# Unconstrained -> simplex via softmax (ensures alpha >= 0, sum = 1)
softmax3 <- function(x) { e <- exp(c(x, 0)); e / sum(e) }

neg_log_lik <- function(x, idx) {
  a    <- softmax3(x)
  pmix <- a[1]*p_src[[1]][idx] + a[2]*p_src[[2]][idx] + a[3]*p_src[[3]][idx]
  pmix <- pmax(pmin(pmix, 1 - 1e-9), 1e-9)
  g    <- g_me[idx]
  -sum(g * log(pmix) + (2L - g) * log(1 - pmix))
}

cat("Running per-chromosome MLE...\n\n")
cat(sprintf("  %-3s  %8s", "chr", "n_SNPs"))
for (lbl in labels) cat(sprintf("  %14s", lbl))
cat("\n")
cat(strrep("-", 20 + 16*length(labels)), "\n")

chr_results <- map_dfr(1:22, function(chr) {
  idx <- which(snp$chr == chr & g_me != 9L)
  if (length(idx) < 500) return(NULL)

  opt <- tryCatch(
    optim(c(0, 0), neg_log_lik, idx = idx,
          method = "BFGS", control = list(maxit = 1000, reltol = 1e-10)),
    error = function(e) NULL
  )
  if (is.null(opt)) return(NULL)

  a <- softmax3(opt$par)
  cat(sprintf("  %-3d  %8d", chr, length(idx)))
  for (v in a) cat(sprintf("  %13.1f%%", v*100))
  cat("\n")
  as_tibble(c(list(chr=chr, n_snps=length(idx)), setNames(as.list(a), sources)))
})

wm <- function(s) weighted.mean(chr_results[[s]], chr_results$n_snps) * 100
cat(strrep("-", 20 + 16*length(labels)), "\n")
cat(sprintf("  %-3s  %8s", "wt.m", ""))
for (s in sources) cat(sprintf("  %13.1f%%", wm(s)))
cat("\n")

for (i in seq_along(sources)) {
  cat(sprintf("Range (%s): %.0f%% - %.0f%%\n", labels[i],
              min(chr_results[[sources[i]]])*100, max(chr_results[[sources[i]]])*100))
}
cat("\n")

# Plot
PALETTE_FILL  <- setNames(c("#E69F00", "#56B4E9", "#009E73"), labels)
PALETTE_COLOR <- setNames(c("#B87D00", "#2E7CA8", "#007050"), labels)

plot_dat <- chr_results %>%
  pivot_longer(all_of(sources), names_to="ancestry", values_to="prop") %>%
  mutate(ancestry = factor(ancestry, levels=labels),
         chr = factor(chr, levels=1:22))

means_df <- tibble(
  ancestry = factor(labels, levels=labels),
  mean = sapply(sources, wm) / 100
)

p <- ggplot(plot_dat, aes(chr, prop*100, fill=ancestry)) +
  geom_col(width=0.8) +
  geom_hline(data=means_df, aes(yintercept=mean*100, color=ancestry),
             linetype="dashed", linewidth=0.7, show.legend=FALSE) +
  scale_fill_manual(values=PALETTE_FILL) +
  scale_color_manual(values=PALETTE_COLOR) +
  scale_y_continuous(limits=c(0,105), breaks=seq(0,100,20),
                     labels=function(x) paste0(x,"%")) +
  labs(title   = paste("Chromosome ancestry painting -", TARGET),
       x       = "Chromosome",
       y       = "Ancestry proportion",
       fill    = "Source",
       caption = paste0("Sources: ", paste(sources, collapse=" + "),
                        "\nDashed lines = weighted genome-wide mean | Method: diploid binomial MLE")) +
  theme_bw(base_size=12) +
  theme(legend.position="bottom", panel.grid.major.x=element_blank())

outfile <- file.path(dirname(F2_DIR), "chr_painting.pdf")
ggsave(outfile, p, width=14, height=7)

csvfile <- file.path(dirname(F2_DIR), "chr_painting.csv")
write_csv(chr_results, csvfile)

cat("Plot saved:", outfile, "\nResults saved:", csvfile, "\nDone.\n")
