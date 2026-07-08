#!/usr/bin/env bash
# Downloads the AADR v66.p1 1240K files used by this project from Harvard Dataverse
# (doi:10.7910/DVN/FFIDCW) into $PROJECT_ROOT/aadr/.
#
# Usage: PROJECT_ROOT=/path/to/ADMIXTOOLS2 ./scripts/download_aadr.sh

set -euo pipefail

: "${PROJECT_ROOT:?Set PROJECT_ROOT first (see .env.example)}"

OUT_DIR="$PROJECT_ROOT/aadr"
mkdir -p "$OUT_DIR"

# label -> Dataverse file id
declare -A FILES=(
  [v66.p1_1240K.aadr.patch.PUB.geno]=13994829
  [v66.p1_1240K.aadr.patch.PUB.snp]=13994514
  [v66.p1_1240K.aadr.patch.PUB.ind]=13994513
  [v66.p1_1240K.aadr.PUB.anno]=13994515
)

for name in "${!FILES[@]}"; do
  id="${FILES[$name]}"
  dest="$OUT_DIR/$name"
  if [ -f "$dest" ]; then
    echo "Skipping $name (already exists)"
    continue
  fi
  echo "Downloading $name..."
  curl -L --fail -o "$dest" "https://dataverse.harvard.edu/api/access/datafile/$id"
done

echo "Done. Files written to $OUT_DIR"
