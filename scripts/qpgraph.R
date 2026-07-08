library(admixtools)
library(tidyverse)
source("/data/scripts/config.R")

model <- MODELS[[BACKGROUND]]
GRAPH_POPS <- unique(c(
  TARGET,
  if (!is.null(model$pool)) model$pool else model$sources,
  "Mbuti"
))

all_pops <- unique(c(GRAPH_POPS, model$outgroup))

cat("=== qpGraph admixture graph ===\n")
cat("Target:", TARGET, "\n")
cat("Populations:", paste(GRAPH_POPS, collapse = ", "), "\n\n")

if (length(list.files(F2_DIR, pattern = "\\.rds$")) == 0) {
  cat("Computing f2 blocks (first run — may take several minutes)...\n")
  extract_f2(MERGED_PREFIX, outdir = F2_DIR, pops = all_pops, overwrite = FALSE)
} else {
  cat("f2 blocks cached, loading...\n")
}
f2 <- f2_from_precomp(F2_DIR, pops = GRAPH_POPS)

cat("Searching for best-fitting admixture graph...\n\n")
graphs <- find_graphs(
  f2,
  outpop   = "Mbuti",
  numadmix = 2,
  stop_gen = 300
)

best  <- graphs %>% arrange(score) %>% slice(1)
cat("Best graph score:", best$score[1], "\n\n")

graph <- best$graph[[1]]
fit   <- qpgraph(f2, graph)

cat("--- Graph fit ---\n")
cat(sprintf("  Score (lower is better): %.4f\n", fit$score))

# Print admixture edges (weight column name varies by version)
if (!is.null(fit$edges) && nrow(fit$edges) > 0) {
  cat("\n  Edge columns available:", paste(names(fit$edges), collapse = ", "), "\n")
  if ("type" %in% names(fit$edges)) {
    adm <- fit$edges %>% filter(type == "admix") %>% select(from, to, weight)
    if (nrow(adm) > 0) {
      cat("\n  Admixture events (dashed arrows in graph):\n")
      for (i in seq_len(nrow(adm))) {
        cat(sprintf("    %s → %s: %.1f%%\n", adm$from[i], adm$to[i], adm$weight[i] * 100))
      }
    }
  }
}

outfile <- file.path(dirname(F2_DIR), "qpgraph_best.pdf")
pdf(outfile, width = 12, height = 9)
plot_graph(fit$edges)
dev.off()
cat("\nGraph plot saved to:", outfile, "\n\nDone.\n")
