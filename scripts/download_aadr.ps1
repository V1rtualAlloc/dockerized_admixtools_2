# Downloads the AADR v66.p1 1240K files used by this project from Harvard Dataverse
# (doi:10.7910/DVN/FFIDCW) into $env:PROJECT_ROOT\aadr\.
#
# Usage: $env:PROJECT_ROOT = "C:\path\to\ADMIXTOOLS2"; .\scripts\download_aadr.ps1

$ErrorActionPreference = "Stop"

if (-not $env:PROJECT_ROOT) {
    throw "Set `$env:PROJECT_ROOT first (see .env.example)"
}

$OutDir = Join-Path $env:PROJECT_ROOT "aadr"
New-Item -ItemType Directory -Force -Path $OutDir | Out-Null

$Files = @{
    "v66.p1_1240K.aadr.patch.PUB.geno" = 13994829
    "v66.p1_1240K.aadr.patch.PUB.snp"  = 13994514
    "v66.p1_1240K.aadr.patch.PUB.ind"  = 13994513
    "v66.p1_1240K.aadr.PUB.anno"       = 13994515
}

foreach ($name in $Files.Keys) {
    $id = $Files[$name]
    $dest = Join-Path $OutDir $name
    if (Test-Path $dest) {
        Write-Host "Skipping $name (already exists)"
        continue
    }
    Write-Host "Downloading $name..."
    Invoke-WebRequest -Uri "https://dataverse.harvard.edu/api/access/datafile/$id" -OutFile $dest
}

Write-Host "Done. Files written to $OutDir"
