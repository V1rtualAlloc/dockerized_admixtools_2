# ADMIXTOOLS2 Ancient DNA Analysis

Comparing 23andMe personal genotype data against the Allen Ancient DNA Resource (AADR) using ADMIXTOOLS2.

## Prerequisites

All analysis runs inside Docker containers, so the only host requirement is Docker itself plus a couple of standard CLI tools.

### Linux

- **Docker Engine** — [official install guide](https://docs.docker.com/engine/install/). On Ubuntu/Debian:
  ```bash
  curl -fsSL https://get.docker.com | sh
  sudo usermod -aG docker "$USER"   # log out/in afterwards so docker works without sudo
  ```
- **curl, awk, sed, grep** — preinstalled on virtually every distro. If missing: `sudo apt install curl gawk sed grep`.
- **bash 4+** (for `scripts/download_aadr.sh`, which uses associative arrays) — preinstalled on all modern distros.

### Windows

- **Docker Desktop** — [official install guide](https://docs.docker.com/desktop/install/windows-install/). Requires WSL2, which Docker Desktop's installer offers to enable for you; if not, enable it first with `wsl --install` in an elevated PowerShell.
- **PowerShell 5.1+** — ships with Windows 10/11 by default, used for `scripts/download_aadr.ps1`.
- The `mkdir -p`, `grep`, `awk`, `sed` steps under [How to add a new 23andMe sample](#how-to-add-a-new-23andme-sample) are written for a POSIX shell. Run them in **WSL** (`wsl --install`, then use its Ubuntu shell) or **Git Bash** (bundled with [Git for Windows](https://git-scm.com/download/win)). Everything else (`docker build`, `docker run`, the PowerShell download script) works natively in PowerShell.

## Setup

Every command below mounts the repo root into a container at `/data`. Set `PROJECT_ROOT` once per shell to the absolute path where you cloned this repo.

Copy the example env file and edit it with your own path:

```bash
cp .env.example .env
# edit .env — set PROJECT_ROOT to the absolute path of this repo
```

Then, in each new shell before running the commands below:

```bash
set -a && source .env && set +a
```

Always use an absolute path — tilde expansion is unreliable with Docker's `-v` flag.

## Data

### AADR v66.p1 (1240k SNP panel)
Downloaded from [Harvard Dataverse](https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/FFIDCW).

Files in `aadr/`:
- `v66.p1_1240K.aadr.patch.PUB.geno` — genotype data (TGENO format, 7.1 GB)
- `v66.p1_1240K.aadr.patch.PUB.snp` — SNP info
- `v66.p1_1240K.aadr.patch.PUB.ind` — 23,089 individuals
- `v66.p1_1240K.aadr.PUB.anno` — sample metadata/annotations
- `v66.p1_1240K.geno/.snp/.ind` — converted PACKEDANCESTRYMAP format (used for analysis)

Download the four raw files above with the provided script (requires `PROJECT_ROOT` set, see [Setup](#setup)):

```bash
# Linux/macOS
./scripts/download_aadr.sh

# Windows
.\scripts\download_aadr.ps1
```

The PACKEDANCESTRYMAP conversion (`v66.p1_1240K.geno/.snp/.ind`) is not downloaded — produce it yourself with `convertf` from the `eigensoft` image (build it first, see [Docker images](#docker-images) below).

Create `aadr/convert_docker.par`:

```
genotypename:    /data/v66.p1_1240K.aadr.patch.PUB.geno
snpname:         /data/v66.p1_1240K.aadr.patch.PUB.snp
indivname:       /data/v66.p1_1240K.aadr.patch.PUB.ind
outputformat:    PACKEDANCESTRYMAP
genooutfilename: /data/v66.p1_1240K.geno
snpoutfilename:  /data/v66.p1_1240K.snp
indoutfilename:  /data/v66.p1_1240K.ind
```

Then run:

```bash
docker run --rm -v $PROJECT_ROOT/aadr:/data eigensoft -p /data/convert_docker.par
```

This reads the whole 7.1 GB `.geno` file, so expect it to take a few minutes and to need roughly another 7 GB of free disk for the output.

## Docker images

### `docker/admixtools/` — AdmixTools + plink 1.9

Compiles AdmixTools from source (convertf, mergeit, qpAdm, etc.) and installs plink 1.9. Used for all data preparation steps.

```bash
docker build -t eigensoft $PROJECT_ROOT/docker/admixtools/
```

### `docker/r-analysis/` — R + ADMIXTOOLS2

Runs ancestry analyses using the `admixtools` R package (v2.0.10).

```bash
docker build -t admixtools2 $PROJECT_ROOT/docker/r-analysis/
```

### `docker/dates/` — DATES v753 (Moorjani lab, LD dating)

```bash
docker build -t dates $PROJECT_ROOT/docker/dates/

# Step 1: run DATES
docker run --rm -v $PROJECT_ROOT:/data dates -p /data/rolloff/dates.par

# Step 2: fit exponential and plot
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/rolloff_plot.R
```

Parameter file: `rolloff/dates.par`. Admixlist file: `rolloff/admixlist.txt` (format: `source1\tsource2\ttarget\toutdir` on one line). Output in `rolloff/output/`.

## Scripts

All analysis scripts live in `scripts/`. Edit `scripts/config.R` to set the target sample and paths — all other scripts source it automatically.

### Generic — work for any target

Driven entirely by `config.R` (`TARGET`, `BACKGROUND`, `MODELS`). Point `config.R` at your own merged dataset and these run as-is, for any ancestry background — the qpAdm sources/outgroups and every population pool the scripts pull from (`MODELS[[BACKGROUND]]$pool`, `WORLD_REFS`) are keyed by `BACKGROUND` or background-agnostic.

**Exception:** `REFERENCES` (in `config.R`) is the one thing that can't be generic — it's the set of modern/historical populations shown alongside the target for comparison, and "relevant comparison population" is inherently regional. Only `REFERENCES_BY_BACKGROUND$european` is populated (this project's original Balkan/Slavic use case). If you set `BACKGROUND` to anything else, `REFERENCES` is empty by default — `qpadm.R`'s comparison table, `pca.R`'s "Reference" group, and the `f3`/`f4` scripts' comparison rows will just show the target with no peers until you add your own region's populations to `REFERENCES_BY_BACKGROUND` (verify names against the AADR `.ind` file first — see the f2 cache gotcha below for why small/high-missingness populations are risky to add).

| Script | What it does |
|--------|-------------|
| `config.R` | Shared config: TARGET, BACKGROUND, MERGED_PREFIX, F2_DIR, MODELS, REFERENCES_BY_BACKGROUND, WORLD_REFS |
| `qpadm.R` | 3-source admixture model + comparison table across reference populations |
| `qpadm_rotate.R` | Tests all 2- and 3-source combinations from `MODELS[[BACKGROUND]]$pool` (falls back to the 3 default sources if no pool is defined); reports passing models |
| `f3_outgroup.R` | Ranks all populations by genetic affinity to the target (outgroup f3) |
| `f3_admixture.R` | Tests all source pairs for admixture signal (negative f3) |
| `f4.R` | D-statistics: balance test for every pair of the background's qpAdm sources vs reference populations |
| `qpgraph.R` | Fits an admixture graph (automated search, numadmix=2) |
| `pca.R` | MDS from pairwise f2 distances; "all populations" and "regional focus" versions. Extra Balkan/medieval populations for finer resolution are only added when `BACKGROUND == "european"` (`REGIONAL_EXTRA` in the script) |
| `rolloff_plot.R` | Fits exponential to DATES output; reports date estimate with jackknife SE |
| `chr_painting.R` | Per-chromosome ancestry via allele-frequency MLE (see prep step below) |

### Case study — the author's own Balkan/Slavic ancestry analysis

These hardcode specific ancient/historical population names (e.g. `Poland_EarlyMedieval_Slav`, `Croatia_EIA`, `Serbia_IronGates_Mesolithic`) and, in `report.R`'s case, pre-written interpretive text about one specific genome. They're a worked example of the workflow, not templates — running them against a different target requires editing the population lists (and, for `report.R`, the prose) by hand, not just changing `config.R`.

| Script | What it does |
|--------|-------------|
| `chg_test.R` | 3-model comparison to test for direct CHG ancestry beyond Yamnaya |
| `iron_gates_test.R` | Compares Loschbour vs Iron Gates vs EHG as the HG proxy in the 3-source model |
| `report.R` | 5-page PDF report: title summary, qpAdm, rotating models, f3 outgroup, chromosome painting, Slavic model |
| `slavic_model.R` | Historical Slavic/pre-Slavic model for the target; structured rotating search (Slavic x Balkan-IA x Roman/Byzantine pool); uses `me/f2_slavic/` cache |
| `slavic_modern.R` | Tests the best Slavic models on modern Balkan/European reference populations (Bulgarian, Albanian, Greek, Hungarian, Czech, Polish, Russian); uses `me/f2_modern/` cache |
| `slavic_pooled.R` | Pools Croatia_EIA + NorthMacedonia_IA + Bulgaria_KapitanAndreevo_EIA → Balkans_IA (N=50) to reduce SE; creates `me/merged_pooled.ind` and `me/f2_pooled/` cache |
| `slavic_focused.R` | Focused Slavic model: 2 Slavic proxies × 4 Balkan IA proxies × 4 Roman proxies |
| `slavic_outgroup_test.R` | Tests alternate outgroup sets for the Slavic model |
| `slavic_pie.R` | Pie chart for the best Slavic model (Poland_EarlyMedieval_Slav + Croatia_EIA); saves `me/slavic_pie.pdf`; uses `me/f2_focused/` cache |

Run any script with:

```bash
docker run --rm \
  -v $PROJECT_ROOT:/data \
  admixtools2 /data/scripts/<script>.R
```

f2 blocks are cached in `F2_DIR` after the first run. If you change the population set in config.R, delete the cache directory and rerun.

> **f2 cache gotcha:** `extract_f2` computes f2 over the SNP intersection across all populations in the call. Including small ancient populations with high missingness (N<10) silently drops thousands of SNPs from every pairwise comparison, inflating SEs and distorting weights throughout the cache. Keep small or high-missingness ancient populations in a separate dedicated cache.

### Chromosome painting prep step

`chr_painting.R` needs a text EIGENSTRAT subset containing just the target and its 3 qpAdm sources (`MODELS[[BACKGROUND]]$sources`) — this isn't produced by any earlier step, so create it once per target/background.

For the author's own sample (`TARGET="me"`, `BACKGROUND="european"`) this is already checked out at `me/subset_extract.par` / `me/subset_poplist.txt`. For a new target — say `TARGET="NAME"` with `MERGED_PREFIX="/data/samples/NAME/merged"` as set up in [How to add a new 23andMe sample](#how-to-add-a-new-23andme-sample) — create both files yourself.

`samples/NAME/subset_poplist.txt` — one population per line: the target, then its 3 sources from `MODELS[[BACKGROUND]]$sources` in `config.R` (e.g. for `european`: `Russia_Samara_EBA_Yamnaya`, `Serbia_IronGates_Mesolithic`, `Turkey_N`):

```
NAME
Russia_Samara_EBA_Yamnaya
Serbia_IronGates_Mesolithic
Turkey_N
```

`samples/NAME/subset_extract.par`:

```
genotypename:    /data/samples/NAME/merged.geno
snpname:         /data/samples/NAME/merged.snp
indivname:       /data/samples/NAME/merged.ind
outputformat:    EIGENSTRAT
genooutfilename: /data/samples/NAME/subset.geno
snpoutfilename:  /data/samples/NAME/subset.snp
indoutfilename:  /data/samples/NAME/subset.ind
poplistname:     /data/samples/NAME/subset_poplist.txt
```

Then run:

```bash
docker run --rm --entrypoint convertf \
  -v $PROJECT_ROOT:/data \
  eigensoft -p /data/samples/NAME/subset_extract.par
```

Then run the analysis normally with `admixtools2 /data/scripts/chr_painting.R` — it reads the subset from `dirname(MERGED_PREFIX)`, so as long as `MERGED_PREFIX` in `config.R` points at `samples/NAME/merged`, it finds the files created above automatically.

### Report

Generates `me/ancestry_report.pdf` (5 pages). Requires `me/chr_painting.csv` — run `chr_painting.R` first if it doesn't exist.

```bash
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/report.R
```

---

## How to add a new 23andMe sample

Replace `NAME` with the sample label (no spaces) and `input.txt` with the path to the raw 23andMe file.

### Step 1 — Convert 23andMe → PLINK bed

23andMe files have MT before X and Y, but plink requires X, Y, MT order. Sort first, then convert.

**Sort (run in terminal — standard Unix tools, no Docker needed):**

```bash
mkdir -p $PROJECT_ROOT/samples/NAME

grep '^#' input.txt > samples/NAME/sorted.txt

grep -v '^#' input.txt | awk 'BEGIN{OFS="\t"} {
  c=$2
  if (c=="X")  n=23
  else if (c=="Y")  n=24
  else if (c=="XY") n=25
  else if (c=="MT") n=26
  else n=int(c)
  print n, $0
}' | sort -k1,1n -k4,4n | cut -f2- >> samples/NAME/sorted.txt
```

**Convert to PLINK (via Docker):**

```bash
docker run --rm --entrypoint plink \
  -v $PROJECT_ROOT:/data \
  eigensoft \
  --23file /data/samples/NAME/sorted.txt NAME NAME \
  --make-bed --out /data/samples/NAME/NAME
```

Output: `samples/NAME/NAME.bed/.bim/.fam`

### Step 2 — Convert PLINK → EIGENSTRAT

23andMe data often has duplicate rsIDs which `convertf` cannot handle. Plink's `--list-duplicate-vars` misses cases where the same rsID appears at two slightly different positions, so use awk instead.

**Remove duplicate rsIDs (run in terminal):**

```bash
awk 'seen[$2]++ {print $2}' samples/NAME/NAME.bim > samples/NAME/NAME_duprsids.txt
```

**Write clean bed without duplicates (via Docker):**

```bash
docker run --rm --entrypoint plink \
  -v $PROJECT_ROOT:/data \
  eigensoft \
  --bfile /data/samples/NAME/NAME \
  --exclude /data/samples/NAME/NAME_duprsids.txt \
  --make-bed \
  --out /data/samples/NAME/NAME_nodup
```

**Convert to EIGENSTRAT:**

Create a parameter file `samples/NAME/plink2eigenstrat.par`:

```
genotypename:    /data/samples/NAME/NAME_nodup.bed
snpname:         /data/samples/NAME/NAME_nodup.bim
indivname:       /data/samples/NAME/NAME_nodup.fam
outputformat:    EIGENSTRAT
genooutfilename: /data/samples/NAME/NAME.geno
snpoutfilename:  /data/samples/NAME/NAME.snp
indoutfilename:  /data/samples/NAME/NAME.ind
```

Then run:

```bash
docker run --rm --entrypoint convertf \
  -v $PROJECT_ROOT:/data \
  eigensoft -p /data/samples/NAME/plink2eigenstrat.par
```

**Fix population label (run in terminal):**

`convertf` sets the population to `???` when converting from PLINK (which has no population column). Fix it before merging:

```bash
sed -i "s/NAME:NAME M        ???/NAME:NAME M             NAME/" samples/NAME/NAME.ind
```

Output: `samples/NAME/NAME.geno/.snp/.ind`

### Step 3 — Merge with AADR

Create a parameter file `samples/NAME/merge.par`:

```
geno1:           /data/aadr/v66.p1_1240K.geno
snp1:            /data/aadr/v66.p1_1240K.snp
ind1:            /data/aadr/v66.p1_1240K.ind
geno2:           /data/samples/NAME/NAME.geno
snp2:            /data/samples/NAME/NAME.snp
ind2:            /data/samples/NAME/NAME.ind
genooutfilename: /data/samples/NAME/merged.geno
snpoutfilename:  /data/samples/NAME/merged.snp
indoutfilename:  /data/samples/NAME/merged.ind
outputformat:    PACKEDANCESTRYMAP
testmismatch:    NO
```

Then run:

```bash
docker run --rm --entrypoint mergeit \
  -v $PROJECT_ROOT:/data \
  eigensoft -p /data/samples/NAME/merge.par
```

Output: `samples/NAME/merged.geno/.snp/.ind` — ready for analysis.

### Step 4 — Run analysis

Edit `scripts/config.R`:

```r
TARGET        <- "NAME"
BACKGROUND    <- "european"   # european / south_asian / east_asian /
#                               middle_eastern / central_asian / native_american
MERGED_PREFIX <- "/data/samples/NAME/merged"
F2_DIR        <- "/data/samples/NAME/f2/"
```

Then run whichever analyses you want:

```bash
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/qpadm.R
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/f3_outgroup.R
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/f4.R
docker run --rm -v $PROJECT_ROOT:/data admixtools2 /data/scripts/qpgraph.R
```

**Outgroup selection note (european background):** outgroups must distinguish the three sources from each other. Include Near Eastern populations (Iran_N, Natufian) for resolution; exclude EHG and CHG which are direct components of Yamnaya. The HG source is `Serbia_IronGates_Mesolithic` (Balkan HG, ~8000 BCE) rather than Loschbour — Iron Gates fits marginally better for a Macedonian/Serbian individual and is geographically more appropriate.

**Reading qpAdm output:**
- `weight` — ancestry proportion (should be 0–1, sum to 1 for a good model)
- `p` in model fit table — p > 0.05 means the model fits; full model is always the first row (pat = `000...`)
