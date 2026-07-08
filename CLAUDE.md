# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project

Comparing personal 23andMe genotype data against the Allen Ancient DNA Resource (AADR) using ADMIXTOOLS2 (R package). The goal is to model ancestry as mixtures of ancient populations using qpAdm and related f-statistic methods.

## Toolchain

Two distinct tools are involved — do not confuse them:

- **AdmixTools** (C binaries) — original Reich Lab command-line suite: `convertf`, `mergeit`, `qpAdm`, etc. Used for data preparation (format conversion, merging). Built via the `Dockerfile`.
- **ADMIXTOOLS2** (`admixtools` R package, `uqrmaie1/admixtools`) — R reimplementation used for the actual analyses. Installed separately in R.

## Docker images

Every command below mounts the repo root into a container at `/data`. Set `PROJECT_ROOT` once per shell to the absolute path where this repo lives — copy `.env.example` to `.env`, set `PROJECT_ROOT` there, then run `set -a && source .env && set +a` before running the commands below.

Always use an absolute path — tilde expansion is unreliable with Docker's `-v` flag.

### `eigensoft` — AdmixTools + plink 1.9

```bash
docker build -t eigensoft $PROJECT_ROOT/docker/admixtools/

# convertf (default entrypoint)
docker run --rm -v $PROJECT_ROOT/aadr:/data eigensoft -p /data/convert_docker.par

# other binaries
docker run --rm --entrypoint mergeit -v $PROJECT_ROOT:/data eigensoft -p /data/merge.par
docker run --rm --entrypoint plink   -v $PROJECT_ROOT:/data eigensoft --23file /data/me/me_sorted.txt ...
```

### `admixtools2` — R + admixtools package

```bash
docker build -t admixtools2 $PROJECT_ROOT/docker/r-analysis/

docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/qpadm.R
```

### `dates` — DATES v753 (Moorjani lab, LD dating)

```bash
docker build -t dates $PROJECT_ROOT/docker/dates/

# Step 1: run DATES (writes output to rolloff/output/)
docker run --rm -v $PROJECT_ROOT:/data dates -p /data/rolloff/dates.par

# Step 2: fit exponential and plot (uses admixtools2 image)
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/rolloff_plot.R
```

Parameter file: `rolloff/dates.par`. Admixlist (source1, source2, target, outdir): `rolloff/admixlist.txt`.
Output in `rolloff/output/`. Single-individual estimates have very wide SE — interpret qualitatively only.

## Scripts

All scripts in `scripts/` source `scripts/config.R` for shared configuration. Edit config.R to change the target sample.

**Generic vs. case study:** `qpadm.R`, `qpadm_rotate.R`, `f3_outgroup.R`, `f3_admixture.R`, `f4.R`, `qpgraph.R`, `pca.R`, `rolloff_plot.R`, and `chr_painting.R` are driven entirely by `config.R` and work for any `TARGET`/`BACKGROUND`. Everything below the `chr_painting.R` row hardcodes specific historical population names (and, for `report.R`, prose about one specific genome) — they're a worked example of the author's own Balkan/Slavic ancestry investigation, not templates for a different target.

**REFERENCES caveat:** `REFERENCES` (comparison populations shown alongside TARGET) is inherently regional, so it's keyed by BACKGROUND as `REFERENCES_BY_BACKGROUND` in config.R. Only `european` is populated — any other BACKGROUND gets an empty comparison set until you add real, AADR-verified populations for that region. Everything else these scripts pull from `MODELS[[BACKGROUND]]` (`sources`, `outgroup`, `pool`) or `WORLD_REFS` is already background-keyed or background-agnostic.

| Script | Purpose |
|--------|---------|
| `config.R` | TARGET, BACKGROUND, MERGED_PREFIX, F2_DIR, MODELS, REFERENCES_BY_BACKGROUND, WORLD_REFS |
| `qpadm.R` | 3-source admixture model + comparison table (columns labeled by the background's actual source names) |
| `qpadm_rotate.R` | Systematic 2- and 3-source model search over `MODELS[[BACKGROUND]]$pool` (falls back to `sources` if no pool defined) |
| `f3_outgroup.R` | Outgroup f3 ranking by affinity to target |
| `f3_admixture.R` | Admixture f3: all source pairs tested for negative f3 (admixture signal) |
| `f4.R` | D-statistics: balance test for every pair of the background's 3 qpAdm sources |
| `qpgraph.R` | Admixture graph fitting |
| `pca.R` | MDS from pairwise f2 distances (separate f2 cache); Balkan-specific extra populations only added when `BACKGROUND == "european"` |
| `rolloff_plot.R` | Fits exponential to DATES output; reports date estimate with jackknife SE |
| `chr_painting.R` | Per-chromosome ancestry proportions via allele-frequency MLE (requires prep step below) |
| `chg_test.R` | *(case study)* 3-model comparison to test for direct CHG ancestry beyond Yamnaya |
| `iron_gates_test.R` | *(case study)* Compares Loschbour vs Iron Gates vs EHG as the HG proxy in the 3-source deep ancestry model; uses `me/f2_ig/` cache |
| `report.R` | *(case study)* 5-page PDF report combining all key results; requires `me/chr_painting.csv` |
| `slavic_model.R` | *(case study)* Historical Slavic/pre-Slavic model; structured rotating search (Slavic x Balkan-IA x Roman/Byzantine) |
| `slavic_modern.R` | *(case study)* Tests best Slavic models on modern Balkan/European reference populations; uses clean f2 cache |
| `slavic_pooled.R` | *(case study)* Pools Croatia_EIA + NorthMacedonia_IA + Bulgaria_KapitanAndreevo_EIA → Balkans_IA (N=50); builds `me/f2_pooled/` cache |
| `slavic_focused.R` | *(case study)* Focused Slavic model: 2 Slavic proxies × 4 Balkan IA proxies × 4 Roman proxies; uses `me/f2_focused/` cache |
| `slavic_outgroup_test.R` | *(case study)* Tests alternate outgroup sets for the Slavic model |
| `slavic_pie.R` | *(case study)* Pie chart for the best Slavic model (Poland_EarlyMedieval_Slav + Croatia_EIA); re-runs qpAdm and saves `me/slavic_pie.pdf`; uses `me/f2_focused/` cache |

### f2 cache pollution warning

`extract_f2` computes f2 over the SNP intersection across **all** populations in the call. Small ancient populations with high missingness (N<10, especially N=4) silently veto thousands of SNPs from every pairwise comparison in the cache, inflating SEs and distorting weights. Keep separate caches:

- `me/f2/` — main cache (deep ancestry: Yamnaya, IronGates, Turkey_N)
- `me/f2_modern/` — Slavic model cache (slavic_model.R + slavic_modern.R); excludes small ancient populations like Bulgaria_Kazanlak_LIA and Albania_Cinamak_IA
- `me/f2_pooled/` — pooled Balkan IA cache (slavic_pooled.R); uses `me/merged_pooled.ind` with remapped population labels
- `me/f2_focused/` — focused Slavic model cache (slavic_focused.R); 20 populations, excludes small pops
- `me/f2_ig/` — HG proxy comparison cache (iron_gates_test.R); Loschbour vs Iron Gates vs EHG

Do not add `Turkey_ImperialRoman` (N=4) or other N<6 populations to any shared cache — they silently veto SNPs across all pairs.

### SE floor for single-individual targets

The jackknife SE in qpAdm has a hard floor of ~±10% when the target is N=1, regardless of source population size. This comes from ~1700 genome blocks used in the jackknife — the signal varies that much block-to-block across one person's genome. Pooling sources (N=50) reduced SE from ±12% → ±10.8%, but cannot go lower. To reach ±5%, the target needs N≥16 individuals — either family members added to the merged dataset, or a modern reference population with sufficient sample size used as a proxy target.

### Chromosome painting prep step

`chr_painting.R` needs a text EIGENSTRAT subset (target + its 3 `MODELS[[BACKGROUND]]$sources`) that no earlier step produces. For `TARGET="me"` it's already checked out as `me/subset_extract.par` / `me/subset_poplist.txt`. For any other target, create the equivalent `<dir>/subset_extract.par` (genotypename/snpname/indivname pointing at `MERGED_PREFIX`, outputformat EIGENSTRAT, poplistname pointing at a poplist file) and `<dir>/subset_poplist.txt` (target + its 3 sources, one per line) yourself — see README.md's "Chromosome painting prep step" for the full template.

```bash
docker run --rm --entrypoint convertf -v $PROJECT_ROOT:/data eigensoft \
  -p /data/me/subset_extract.par
# Output: me/subset.geno / .snp / .ind  (~84 MB text)
```

Then run the analysis:

```bash
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/chr_painting.R
```

Output: `me/chr_painting.pdf` and `me/chr_painting.csv`. EIGENSTRAT format is SNP-major (rows = SNPs, cols = individuals). Steppe proportions vary 32–47% across chromosomes; Balkan HG 4–10%; Anatolian 44–63%.

### Report

Requires `me/chr_painting.csv` (run `chr_painting.R` first). Output: `me/ancestry_report.pdf`.

```bash
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/report.R
```

## Data

### AADR v66.p1 (1240k SNP panel)

Located in `aadr/`. Source: [Harvard Dataverse doi:10.7910/DVN/FFIDCW](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW).

| File | Format | Notes |
|------|--------|-------|
| `v66.p1_1240K.aadr.patch.PUB.geno` | TGENO | Original download, 7.1 GB |
| `v66.p1_1240K.aadr.patch.PUB.snp` | — | SNP info |
| `v66.p1_1240K.aadr.patch.PUB.ind` | — | 23,089 individuals |
| `v66.p1_1240K.aadr.PUB.anno` | TSV | Sample metadata |
| `v66.p1_1240K.geno/.snp/.ind` | PACKEDANCESTRYMAP | Converted, use these for analysis |

This conversion is not shipped in the repo (`aadr/` is gitignored) — a fresh clone needs to run it (see README.md's AADR section). Once produced locally, there's no need to re-run it unless the source files change.
