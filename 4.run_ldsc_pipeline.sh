#!/usr/bin/env bash
set -ueo pipefail

############################################################
# 0. USER-CONFIGURABLE PATHS (must be edited manually)
############################################################

# Path to ldsc.py
LDSC="${LDSC:-/home/mza/ldsc/ldsc-master/ldsc.py}"

# PLINK reference panel prefix (without chromosome numbers or extensions)
REF="${REF:-/io/mza/ldsc/1000G_Phase3_plink/1000G_EUR_Phase3_plink/1000G.EUR.QC.}"

# HapMap3 SNP list (optional)
HM3="${HM3:-/io/mza/ldsc/w_hm3.snplist}"

# Baseline model prefix (e.g., baseline_plus_EC)
BASE="${BASE:-/io/mza/ldsc/work_ATAC_v2/ldscores_merged/baseline_plus_EC/baseline_plus_EC.hm3pruned.}"

# LD weights prefix
WEIGHTS="${WEIGHTS:-/io/mza/ldsc/1000G_Phase3_weights/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC.}"

############################################################
# 1. INPUT/OUTPUT DIRECTORIES
############################################################

# Directory containing *.specific_peaks.bed files
BEDDIR="${BEDDIR:-./ldsc_peaksets}"

# GWAS directories
declare -a GWAS_ROOTS=(
  ./ldsc/gwas
)

# All LDSC outputs will be stored here
OUT_ROOT="${OUT_ROOT:-./ldsc}"

WORKROOT="${WORKROOT:-${OUT_ROOT}/work}"
L2ROOT="${L2ROOT:-${WORKROOT}/ldscores}"
SCATAC="${SCATAC:-${OUT_ROOT}/scatac}"
LOGDIR="${LOGDIR:-${OUT_ROOT}/_logs}"
BINDIR="${BINDIR:-${OUT_ROOT}/_bin}"
CACHEDIR="${CACHEDIR:-${OUT_ROOT}/_cache}"

mkdir -p "$LOGDIR" "$BINDIR" "$L2ROOT" "$SCATAC" "$CACHEDIR"

############################################################
# 2. RUNTIME PARAMETERS
############################################################

PAR="${PAR:-8}"         # parallel jobs for h2 step
ONLY_CELL="${ONLY_CELL:-}"   # run only one cell (optional)
AGG_ONLY="${AGG_ONLY:-0}"    # skip compute; only aggregate summaries

if [[ "${DEBUG:-0}" == "1" ]]; then set -x; fi

ts(){ date "+%F %T"; }

############################################################
# 3. PATH CHECKING
############################################################

check_path() {
    if [[ ! -e "$1" ]]; then
        echo "[FATAL] Missing required path: $1"
        echo "Please edit the top section of this script to set LDSC/REF/BASE/WEIGHTS correctly."
        exit 1
    fi
}

echo "=== [CHECK] Verifying required paths ==="
check_path "$LDSC"
check_path "${REF}1.bim"
check_path "${BASE}1.l2.ldscore.gz"
check_path "${WEIGHTS}1.l2.ldscore.gz"
if [[ -n "$HM3" ]]; then check_path "$HM3"; fi
echo "[OK] All required paths found"

############################################################
# 4. BUILD A GLOBAL SNP LIST FROM BASELINE (keep exact order)
############################################################

BASELINE_SNPLIST_GLOBAL="${CACHEDIR}/baseline.hm3pruned.snplist"

if [[ ! -s "$BASELINE_SNPLIST_GLOBAL" ]]; then
  echo "[$(ts)] Building unified baseline SNP list..."
  : > "$BASELINE_SNPLIST_GLOBAL"
  for c in {1..22}; do
    f="${BASE}${c}.l2.ldscore.gz"
    [[ -f "$f" ]] || { echo "[FATAL] Missing baseline l2: $f"; exit 2; }
    zcat "$f" \
      | awk -F'\t' '
          NR==1 && ($2=="SNP" || $1=="SNP"){next}
          { if ($2!="") print $2; else print $1 }
        ' >> "$BASELINE_SNPLIST_GLOBAL"
  done
fi

[[ -s "$BASELINE_SNPLIST_GLOBAL" ]] || { echo "[FATAL] Failed to build baseline SNP list"; exit 2; }
echo "[$(ts)] Using baseline SNP list: $BASELINE_SNPLIST_GLOBAL"

############################################################
# 5. CREATE h2 WORKER SCRIPT
############################################################

H2WORKER="${BINDIR}/h2_worker.sh"

cat > "$H2WORKER" <<'EOS'
#!/usr/bin/env bash
set -ueo pipefail

GWAS_FULL="$1"
CELL="$2"

LDSC="${LDSC:?}"
BASE="${BASE:?}"
WEIGHTS="${WEIGHTS:?}"
L2ROOT="${L2ROOT:?}"
SCATAC="${SCATAC:?}"

GWAS_BASE=$(basename "$GWAS_FULL")
SUM="$GWAS_FULL"

OUTROOT="${SCATAC}/${GWAS_BASE%.*}"
OUTDIR="${OUTROOT}/${CELL}"
mkdir -p "$OUTDIR"

LOG="${OUTDIR}/h2_${CELL}_vs_baseline.log"

if [[ -f "$LOG" ]] && grep -q "^Lambda GC:" "$LOG"; then
  echo "[SKIP] ${CELL} × ${GWAS_BASE} already completed"
  exit 0
fi

# Check completeness of l2 files
missing=0
for c in {1..22}; do
  [[ -f "${L2ROOT}/${CELL}/${CELL}.${c}.l2.ldscore.gz" ]] || missing=1
done
if [[ $missing -eq 1 ]]; then
  echo "[ERR] Incomplete l2 for ${CELL}"
  exit 1
fi

python "$LDSC" \
  --h2 "$SUM" \
  --ref-ld-chr "${BASE},${L2ROOT}/${CELL}/${CELL}." \
  --w-ld-chr  "${WEIGHTS}" \
  --not-M-5-50 \
  --out "${OUTDIR}/h2_${CELL}_vs_baseline" \
  --print-coefficients \
  > >(tee "$LOG") 2>&1

coef=$(sed -n '/^Coefficients:/,/^Coefficient SE:/p' "$LOG" \
        | sed 's/^Coefficients://; /^Coefficient SE:/d' \
        | tr ' ' '\n' | grep -E '^-?[0-9.eE+-]+$' | tail -n1 || true)

se=$(sed -n '/^Coefficient SE:/,/^Lambda GC:/p' "$LOG" \
        | sed 's/^Coefficient SE://; /^Lambda GC:/d' \
        | tr ' ' '\n' | grep -E '^-?[0-9.eE+-]+$' | tail -n1 || true)

z="NA"
if [[ -n "${coef:-}" && -n "${se:-}" && "$se" != "0" ]]; then
  z=$(python - <<PY
c=float(${coef})
s=float(${se})
print(c/s)
PY
)
fi

echo -e "cell\tcoef\tse\tz\n${CELL}\t${coef}\t${se}\t${z}" \
  > "${OUTDIR}/summary_${GWAS_BASE%.*}_${CELL}.tsv"
EOS

chmod +x "$H2WORKER"

############################################################
# 6. DISCOVER CELLS AND GWAS FILES
############################################################

mapfile -t CELLS < <(
  ls -1 "${BEDDIR}"/*.specific_peaks.bed 2>/dev/null \
    | xargs -I{} basename {} \
    | sed -E 's/\.specific_peaks\.bed$//' \
    | sort
)

declare -a GWAS_FILES=()
for root in "${GWAS_ROOTS[@]}"; do
  if [[ -d "$root" ]]; then
    while IFS= read -r f; do GWAS_FILES+=("$f"); done \
      < <(ls -1 "$root"/*_ldsc.sumstats.gz "$root"/*.sumstats.gz 2>/dev/null || true)
  fi
done
mapfile -t GWAS_FILES < <(printf "%s\n" "${GWAS_FILES[@]}" | sort -u)

echo "[INFO] Cells found: ${#CELLS[@]}"
echo "[INFO] GWAS files:  ${#GWAS_FILES[@]}"

############################################################
# 7. FUNCTION: Build annotation (.annot.gz) per cell
############################################################
build_annot_for_cell() {
  local cell="$1"
  local bed="${BEDDIR}/${cell}.specific_peaks.bed"
  local outdir="${L2ROOT}/${cell}"
  mkdir -p "$outdir"

  local need=0
  for c in {1..22}; do
    local f="${outdir}/${cell}.${c}.annot.gz"
    if [[ ! -f "$f" ]]; then need=1; break; fi
    local hdr
    hdr=$(zcat "$f" | head -n1 || true)
    [[ "$hdr" == $'CHR\tBP\tSNP\tCM\tANNOT' ]] || { need=1; break; }
  done

  if [[ $need -eq 0 ]]; then
    echo "[ANNOT] ${cell}: annotation already exists"
    return 0
  fi

  echo "[ANNOT] Rebuilding annotation for ${cell}"

  for c in {1..22}; do
    awk 'OFS="\t"{print "chr"$1, $4-1, $4, $2, $4}' \
      ${REF}${c}.bim > ${outdir}/tmp.${c}.snpbed

    bedtools intersect -c \
      -a ${outdir}/tmp.${c}.snpbed \
      -b "$bed" \
      | awk 'BEGIN{OFS="\t"}{ann=($6>0?1:0); print $1,$5,$4,0,ann}' \
      > ${outdir}/${cell}.${c}.annot

    (echo -e "CHR\tBP\tSNP\tCM\tANNOT"; cat ${outdir}/${cell}.${c}.annot) \
      | gzip -c > ${outdir}/${cell}.${c}.annot.gz

    rm -f ${outdir}/${cell}.${c}.annot ${outdir}/tmp.${c}.snpbed
  done

  echo "[ANNOT] ${cell} annotation rebuilt"
}

############################################################
# 8. FUNCTION: Run LD Score computation per cell
############################################################
run_l2_for_cell() {
  local cell="$1"
  local outdir="${L2ROOT}/${cell}"
  mkdir -p "$outdir"

  echo "[L2] Checking baseline SNP counts for ${cell}"

  declare -A BASE_M
  for c in {1..22}; do
    local bf="${BASE}${c}.l2.ldscore.gz"
    local m
    m=$(zcat "$bf" | awk -F'\t' 'NR==1 && ($2=="SNP" || $1=="SNP"){next} {n++} END{print (n+0)}')
    BASE_M[$c]=$m
  done

  for c in {1..22}; do
    local f="${outdir}/${cell}.${c}.l2.ldscore.gz"
    if [[ -f "$f" ]]; then
      if ! gzip -t "$f" 2>/dev/null; then
        echo "[WARN] Corrupted: removing $f"
        rm -f "$f"
      else
        local m2
        m2=$(zcat "$f" | awk -F'\t' 'NR==1 && ($2=="SNP"||$1=="SNP"){next}{n++}END{print (n+0)}')
        if [[ "${m2:-0}" -ne "${BASE_M[$c]}" ]]; then
          echo "[WARN] ${cell}.chr${c}: SNP mismatch ${m2} != ${BASE_M[$c]} → forcing recompute"
          rm -f "$f"
        fi
      fi
    fi
  done

  for c in {1..22}; do
    local f="${outdir}/${cell}.${c}.l2.ldscore.gz"
    [[ -f "$f" ]] && continue

    echo "[L2] Computing LD scores for ${cell}, chr ${c}"

    python "$LDSC" \
      --l2 \
      --bfile     "${REF}${c}" \
      --ld-wind-cm 1 \
      --annot     "${outdir}/${cell}.${c}.annot.gz" \
      --out       "${outdir}/${cell}.${c}" \
      --print-snps "$BASELINE_SNPLIST_GLOBAL" \
      > "${outdir}/${cell}.${c}.l2.log" 2>&1
  done

  local miss=0
  local bad=0
  for c in {1..22}; do
    local f="${outdir}/${cell}.${c}.l2.ldscore.gz"
    [[ -f "$f" ]] || { miss=$((miss+1)); continue; }
    local m2
    m2=$(zcat "$f" | awk -F'\t' 'NR==1 && ($2=="SNP"||$1=="SNP"){next}{n++}END{print (n+0)}')
    [[ "${m2:-0}" -eq "${BASE_M[$c]}" ]] || bad=$((bad+1))
  done

  if [[ $miss -eq 0 && $bad -eq 0 ]]; then
    echo "[L2] ${cell}: all l2 files OK"
  else
    echo "[ERR] ${cell}: l2 incomplete/mismatched (missing=$miss, mismatched=$bad)"
  fi
}

############################################################
# 9. AGGREGATE: per-disease summary tables
############################################################
rebuild_per_disease_summaries() {
  echo "[$(ts)] Rebuilding per-disease summary tables"

  while IFS= read -r -d '' s; do
    disease=$(basename "$(dirname "$(dirname "$s")")")
    cell=$(basename "$(dirname "$s")")
    outdir=$(dirname "$(dirname "$s")")
    summary_file=$(ls "${outdir}"/summary_*.tsv 2>/dev/null | head -n1 || true)

    [[ -z "$summary_file" ]] && summary_file="${outdir}/summary_${disease}.tsv"
    [[ -f "$summary_file" ]] || echo -e "cell\tcoef\tse\tz" > "$summary_file"

    tmp=$(mktemp)
    awk -F'\t' -v OFS='\t' -v CELL="$cell" '
      NR==1 {print; next}
      $1==CELL {next}
      {print}
    ' "$summary_file" > "$tmp"
    mv "$tmp" "$summary_file"

    tail -n +2 "$s" >> "$summary_file"
  done < <(find "$SCATAC" -mindepth 3 -maxdepth 3 -type f -name 'summary_*_*.tsv' -print0)
}

############################################################
# 10. AGGREGATE: build long & wide matrices
############################################################
build_global_tables() {
  echo "[$(ts)] Building global long/wide tables"

  local OUTDIR="$SCATAC/_combined"
  mkdir -p "$OUTDIR"

  local tmp_long; tmp_long=$(mktemp)
  echo -e "disease\tcell\tcoef\tse\tz\tsource" > "$tmp_long"

  while IFS= read -r -d '' S; do
    local disease; disease=$(basename "$(dirname "$S")")
    awk -v d="$disease" -v src="$S" -F'\t' '
      BEGIN {OFS="\t"}
      NR==1 {next}
      {print d,$1,$2,$3,$4,src}
    ' "$S" >> "$tmp_long"
  done < <(find "$SCATAC" -mindepth 1 -maxdepth 2 -type f -name 'summary_*.tsv' -print0)

  local LONG="$OUTDIR/combined_long.tsv"
  mv "$tmp_long" "$LONG"

  python - <<'PY'
import os, pandas as pd, numpy as np

root = os.environ.get("SCATAC_COMBINED_ROOT", "./ldsc/scatac/_combined")
longf = os.path.join(root, "combined_long.tsv")
zwide = os.path.join(root, "combined_z.tsv")
cwide = os.path.join(root, "combined_coef.tsv")
swide = os.path.join(root, "combined_se.tsv")

df = pd.read_csv(longf, sep="\t")
df = df[[c for c in df.columns if not c.startswith("Unnamed")]]

for col in ("coef","se","z"):
    df[col] = pd.to_numeric(df[col], errors="coerce")

z = df.pivot(index="cell", columns="disease", values="z")
c = df.pivot(index="cell", columns="disease", values="coef")
s = df.pivot(index="cell", columns="disease", values="se")

z_order_cols = z.abs().max().sort_values(ascending=False).index
z_order_rows = z.abs().max(axis=1).sort_values(ascending=False).index

z = z.loc[z_order_rows, z_order_cols]
c = c.reindex_like(z)
s = s.reindex_like(z)

z.to_csv(zwide, sep="\t")
c.to_csv(cwide, sep="\t")
s.to_csv(swide, sep="\t")

print("=== Summary ===")
print("Long: ", longf)
print("Wide z: ", zwide)
print("Wide coef: ", cwide)
print("Wide se: ", swide)
PY
}

############################################################
# 11. MAIN LOOP: per cell → annot → l2 → h2
############################################################

export L2ROOT SCATAC BASELINE_SNPLIST_GLOBAL \
       LDSC REF BASE WEIGHTS

for CELL in "${CELLS[@]}"; do
  [[ -n "$ONLY_CELL" && "$CELL" != "$ONLY_CELL" ]] && continue

  echo "[$(ts)] === CELL START: ${CELL} ==="

  build_annot_for_cell "$CELL"
  run_l2_for_cell "$CELL"

  ok=1
  for c in {1..22}; do
    [[ -f ${L2ROOT}/${CELL}/${CELL}.${c}.l2.ldscore.gz ]] || ok=0
  done

  if [[ $ok -eq 0 ]]; then
    echo "[WARN] ${CELL}: skipping h2 due to incomplete l2"
    continue
  fi

  echo "[H2] Launching parallel h2 jobs for ${CELL} (P=${PAR})"
  printf "%s\n" "${GWAS_FILES[@]}" \
    | xargs -P "$PAR" -I{} bash "$H2WORKER" "{}" "$CELL"

  echo "[$(ts)] === CELL END: ${CELL} ==="
done

############################################################
# FINAL: aggregate all results
############################################################

rebuild_per_disease_summaries
build_global_tables

echo "[$(ts)] DONE"
