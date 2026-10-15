#!/usr/bin/env bash
# Module 1 data step: download every open input into data/raw/mr/<source>/ and
# write data/raw/mr/MANIFEST.tsv (path, url, bytes, sha256, downloaded_utc, purpose).
#
# Usage: bash src/mr/download_inputs.sh            (from the repository root)
#
# - Resumable: curl -C - continues a partial file; a file whose size already
#   equals the server's Content-Length is not fetched again.
# - Open sources only (GWAS Catalog FTP over HTTPS, FinnGen R10 public bucket,
#   eQTLGen, UCSC, MAGMA's 1000 Genomes EUR reference). No login, no token.
# - Alzheimer's / dementia OUTCOME files go to data/raw/mr/quarantine_outcome/.
#   This script only downloads and checksums them; it never opens them. They
#   must not be read until prereg/mr.md is frozen (only the header line may be
#   inspected, for the column map).
# - TLS verification is never disabled.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
RAW="$ROOT/data/raw/mr"
MANIFEST="$RAW/MANIFEST.tsv"
PROGRESS="$RAW/download_progress.log"
mkdir -p "$RAW"

GC=https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics
FG=https://storage.googleapis.com/finngen-public-data-r10/summary_stats
EQ=https://molgenis26.gcc.rug.nl/downloads/eqtlgen/cis-eqtl

# relative path under data/raw/mr | URL | purpose
FILES=(
  "eqtlgen/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz|$EQ/2019-12-11-cis-eQTLsFDR-ProbeLevel-CohortInfoRemoved-BonferroniAdded.txt.gz|exposure: eQTLGen full cis-eQTL results (blood, GRCh37, Z-scores)"
  "eqtlgen/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz|$EQ/2018-07-18_SNP_AF_for_AlleleB_combined_allele_counts_and_MAF_pos_added.txt.gz|exposure: eQTLGen allele frequencies of the assessed allele (hg19)"
  "eqtlgen/README_cis|$EQ/README_cis|exposure: eQTLGen cis README (column definitions)"
  "eqtlgen/README_allele_frequency|$EQ/README_allele_frequency|exposure: eQTLGen AF README (column definitions)"
  "gwascatalog/GCST007609/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz|$GC/GCST007001-GCST008000/GCST007609/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.rhrmean.assoc.gz|positive control: resting heart rate (SLC6A2, ADRB1)"
  "gwascatalog/GCST007609/readme.txt|$GC/GCST007001-GCST008000/GCST007609/readme.txt|positive control: resting heart rate README"
  "gwascatalog/GCST006630/Evangelou_30224653_DBP.txt.gz|$GC/GCST006001-GCST007000/GCST006630/Evangelou_30224653_DBP.txt.gz|positive control: diastolic blood pressure (ADRA1A)"
  "gwascatalog/GCST006630/README_Evangelou_30224653.txt|$GC/GCST006001-GCST007000/GCST006630/README_Evangelou_30224653.txt|positive control: DBP README"
  "gwascatalog/GCST006624/Evangelou_30224653_SBP.txt.gz|$GC/GCST006001-GCST007000/GCST006624/Evangelou_30224653_SBP.txt.gz|positive control: systolic blood pressure (DBH, CACNA1D)"
  "gwascatalog/GCST006624/README_Evangelou_30224653.txt|$GC/GCST006001-GCST007000/GCST006624/README_Evangelou_30224653.txt|positive control: SBP README"
  "gwascatalog/GCST90002232/GCST90002232_buildGRCh37.tsv.gz|$GC/GCST90002001-GCST90003000/GCST90002232/GCST90002232_buildGRCh37.tsv.gz|positive control: fasting glucose (ADRA2A)"
  "gwascatalog/GCST90002232/GCST90002232_buildGRCh37.tsv.gz-meta.yaml|$GC/GCST90002001-GCST90003000/GCST90002232/GCST90002232_buildGRCh37.tsv.gz-meta.yaml|positive control: fasting glucose metadata"
  "gwascatalog/GCST90002232/README.txt|$GC/GCST90002001-GCST90003000/GCST90002232/README.txt|positive control: fasting glucose README"
  "gwascatalog/GCST90002232/md5sum.txt|$GC/GCST90002001-GCST90003000/GCST90002232/md5sum.txt|positive control: fasting glucose md5"
  "gwascatalog/GCST007561/Dashti_30846698_sleepdurationsumstats.txt.zip|$GC/GCST007001-GCST008000/GCST007561/Dashti_30846698_sleepdurationsumstats.txt.zip|positive control: sleep duration (HCRTR2)"
  "gwascatalog/GCST007561/Dashti_30846698_fullUKBB_Sleepduration_summary_stats_README.txt|$GC/GCST007001-GCST008000/GCST007561/Dashti_30846698_fullUKBB_Sleepduration_summary_stats_README.txt|positive control: sleep duration README"
  "gwascatalog/GCST002222/jointGwasMc_LDL.txt.gz|$GC/GCST002001-GCST003000/GCST002222/jointGwasMc_LDL%20(1).txt.gz|validation mr_ldl_hmgcr_cad: LDL-C exposure GWAS (GLGC 2013)"
  "gwascatalog/GCST003116/cad.add.160614.website.txt|$GC/GCST003001-GCST004000/GCST003116/cad.add.160614.website.txt|validation mr_ldl_hmgcr_cad: CAD outcome GWAS (CARDIoGRAMplusC4D 1000G)"
  "gwascatalog/GCST003116/cad.add.readme|$GC/GCST003001-GCST004000/GCST003116/cad.add.readme|validation mr_ldl_hmgcr_cad: CAD README"
  "liftover/hg38ToHg19.over.chain.gz|https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz|liftover: UCSC GRCh38 -> GRCh37 chain"
  "ld_1kg_eur/g1000_eur.zip|https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download|LD reference: 1000 Genomes Phase 3 EUR PLINK bfile, GRCh37 (MAGMA reference data, 'European')"
  "finngen_r10/R10_manifest.tsv|$FG/R10_manifest.tsv|FinnGen R10 manifest (endpoint names, case/control counts)"
  "quarantine_outcome/GCST90027158/GCST90027158_buildGRCh38.tsv.gz|$GC/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz|AD OUTCOME (QUARANTINED, not read): Bellenguez 2022 stage 1"
  "quarantine_outcome/GCST90027158/GCST90027158_buildGRCh38.tsv.gz-meta.yaml|$GC/GCST90027001-GCST90028000/GCST90027158/GCST90027158_buildGRCh38.tsv.gz-meta.yaml|AD outcome metadata (study-level, no data rows)"
  "quarantine_outcome/GCST90027158/README.txt|$GC/GCST90027001-GCST90028000/GCST90027158/README.txt|AD outcome README (study-level, no data rows)"
  "quarantine_outcome/GCST90027158/md5sum.txt|$GC/GCST90027001-GCST90028000/GCST90027158/md5sum.txt|AD outcome md5 (study-level, no data rows)"
  "quarantine_outcome/GCST007511/Kunkle_etal_Stage1_results.txt|$GC/GCST007001-GCST008000/GCST007511/Kunkle_etal_Stage1_results.txt|AD OUTCOME (QUARANTINED, not read): Kunkle 2019 stage 1, proxy-free sensitivity"
  "quarantine_outcome/finngen_R10_G6_ALZHEIMER/finngen_R10_G6_ALZHEIMER.gz|$FG/finngen_R10_G6_ALZHEIMER.gz|AD OUTCOME (QUARANTINED, not read): FinnGen R10 G6_ALZHEIMER replication"
)

log() { echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) $*" | tee -a "$PROGRESS"; }

remote_size() {
  curl -sS -I -L --max-time 60 "$1" | tr -d '\r' | awk 'tolower($1)=="content-length:"{n=$2} END{print n+0}'
}

for entry in "${FILES[@]}"; do
  IFS='|' read -r rel url purpose <<< "$entry"
  dest="$RAW/$rel"
  mkdir -p "$(dirname "$dest")"
  want=$(remote_size "$url" || echo 0)
  have=$( [ -f "$dest" ] && stat -c %s "$dest" || echo 0 )
  if [ "$want" -gt 0 ] && [ "$have" -eq "$want" ]; then
    log "present  $rel ($have bytes)"
  else
    log "fetching $rel ($have of $want bytes)"
    curl -sS -L --fail --retry 5 --retry-delay 10 -C - -o "$dest" "$url" || {
      # a server that ignores Range for a complete file answers 416; anything else is a real failure
      log "FAILED   $rel (curl exit $?)"; continue; }
    log "done     $rel ($(stat -c %s "$dest") bytes)"
  fi
  [ -f "$dest.downloaded_utc" ] || date -u +%Y-%m-%dT%H:%M:%SZ > "$dest.downloaded_utc"
done

# Manifest: rebuilt from the files on disk every run.
printf 'path\turl\tbytes\tsha256\tdownloaded_utc\tpurpose\n' > "$MANIFEST.tmp"
for entry in "${FILES[@]}"; do
  IFS='|' read -r rel url purpose <<< "$entry"
  dest="$RAW/$rel"
  if [ -f "$dest" ]; then
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "data/raw/mr/$rel" "$url" "$(stat -c %s "$dest")" \
      "$(sha256sum "$dest" | cut -d' ' -f1)" "$(cat "$dest.downloaded_utc" 2>/dev/null || echo NA)" "$purpose" >> "$MANIFEST.tmp"
  else
    printf '%s\t%s\tNA\tNA\tNA\tMISSING: %s\n' "data/raw/mr/$rel" "$url" "$purpose" >> "$MANIFEST.tmp"
  fi
done
mv "$MANIFEST.tmp" "$MANIFEST"
log "manifest written: $MANIFEST"
