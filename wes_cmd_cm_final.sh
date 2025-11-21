#!/bin/bash
#SBATCH --job-name=multi_wes_cmd_cm4
#SBATCH --output=/home/shamita.s/RAO/logs/multi_%j.out
#SBATCH --error=/home/shamita.s/RAO/logs/multi_%j.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=36:00:00
#SBATCH --partition=compute

set -eo pipefail
source ~/miniconda3/etc/profile.d/conda.sh
conda activate gatk_env

# =========================
# USER CONFIG
# =========================
BASE_DIR="/home/shamita.s/RAO"
FASTQ_DIR="$BASE_DIR/fastq"
RUNS_FILE="$BASE_DIR/runs.txt"
REF="/home/shamita.s/miniconda3/Grch38/hg38.fa"
KNOWN="/home/shamita.s/miniconda3/gatk/known_sites/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"

GATK="gatk"
PICARD_JAR="/home/shamita.s/miniconda3/envs/picard-env/share/picard-*/picard.jar"

# ANNOVAR: adjust to where you unpacked annovar/annovar
ANNOVAR_DIR="$HOME/tools/annovar/annovar"
HUMANDB="$ANNOVAR_DIR/humandb"

# Existing DB references (kept for other uses)
CLINVAR="/home/shamita.s/databases/clinvar/clinvar.vcf.gz"
GNOMAD_SITE="/home/shamita.s/databases/gnomad/gnomad.genomes.v3.1.2.sites.chr1.vcf.bgz"
DBNSFP="/home/shamita.s/databases/dbNSFP/dbNSFP4.4a_grch38.txt.gz"

HIGH_COMPLEXITY_BED="/home/shamita.s/databases/filters/high_complexity_regions.hg38.bed"
TARGET_BED=""

OUT_COHORT="$BASE_DIR/OUTPUT"
mkdir -p "$OUT_COHORT" "$FASTQ_DIR" "$BASE_DIR/logs"

# Threads: prefer SLURM allocation if present
THREADS="${SLURM_CPUS_PER_TASK:-16}"
export THREADS

# SAMPLES
mapfile -t SAMPLES < <(cut -f1 "$RUNS_FILE" | sed '/^$/d')

# =========================
# STEP 0: DOWNLOAD (integrated)
# =========================
function download_fastqs() {
  echo "[DOWNLOAD] Checking FASTQs from $RUNS_FILE"
  tail -n +2 "$RUNS_FILE" | while IFS=$'\t' read -r RUN_ID URLS; do
    [[ -z "$RUN_ID" || -z "$URLS" || "$RUN_ID" == "run_accession" || "$URLS" == "fastq_ftp" || "$URLS" == "-" ]] && continue
    SAMPLE_DIR="$FASTQ_DIR/$RUN_ID"
    mkdir -p "$SAMPLE_DIR"
    if ls "$SAMPLE_DIR/${RUN_ID}_1.fastq.gz" "$SAMPLE_DIR/${RUN_ID}_2.fastq.gz" >/dev/null 2>&1; then
      echo "  - $RUN_ID: already present"
      continue
    fi
    echo "  - $RUN_ID: downloading"
    pushd "$SAMPLE_DIR" >/dev/null
      for url in ${URLS//;/ }; do
        if [[ "$url" =~ ^ftp\.sra\.ebi\.ac\.uk/.*\.fastq\.gz$ ]]; then
          wget -c "ftp://$url" &
        elif [[ "$url" =~ ^https?://.*\.fastq\.gz$ ]]; then
          wget -c "$url" &
        else
          echo "    [WARN] Skipping malformed URL: $url"
        fi
      done
      # wait for parallel wget jobs for this sample
      wait

      # rename logic unchanged
      if compgen -G "*.fastq.gz" >/dev/null; then
        for f in *.fastq.gz; do
          if [[ "$f" == *"_1.fastq.gz" ]]; then
            target="${RUN_ID}_1.fastq.gz"
            if [[ "$f" -ef "$target" ]]; then continue; fi
            if [[ -f "$target" ]]; then
              src_size=$(stat -c%s "$f"); tgt_size=$(stat -c%s "$target")
              if (( src_size > tgt_size )); then mv -f "$f" "$target"; else rm -f "$f"; fi
            else
              mv -f "$f" "$target"
            fi
          fi
          if [[ "$f" == *"_2.fastq.gz" ]]; then
            target="${RUN_ID}_2.fastq.gz"
            if [[ "$f" -ef "$target" ]]; then continue; fi
            if [[ -f "$target" ]]; then
              src_size=$(stat -c%s "$f"); tgt_size=$(stat -c%s "$target")
              if (( src_size > tgt_size )); then mv -f "$f" "$target"; else rm -f "$f"; fi
            else
              mv -f "$f" "$target"
            fi
          fi
        done
      fi

    popd >/dev/null
  done
  echo "[DOWNLOAD] Done"
}

# =========================
# STEP 1–7: Per-sample processing ? GVCF (threads used where possible)
# =========================
function process_sample() {
  local SAMPLE_ID="$1"
  echo "========== Processing $SAMPLE_ID (threads=$THREADS) =========="
  local OUTDIR="$OUT_COHORT/OUTPUT/$SAMPLE_ID"
  local VCFBAM="$OUTDIR/VCF_and_BAM"
  mkdir -p "$OUTDIR/logs" "$VCFBAM/vcf_plots"
  cd "$OUTDIR"
  

 # Skip if GVCF already exists
  if [[ -f "$VCFBAM/${SAMPLE_ID}.g.vcf.gz" ]]; then
    echo "[SKIP] $SAMPLE_ID: GVCF already exists at $VCFBAM/"
    return 0
  fi


  local FQ1="$FASTQ_DIR/$SAMPLE_ID/${SAMPLE_ID}_1.fastq.gz"
  local FQ2="$FASTQ_DIR/$SAMPLE_ID/${SAMPLE_ID}_2.fastq.gz"
  if [[ ! -f "$FQ1" || ! -f "$FQ2" ]]; then
    echo "[ERROR] FASTQs missing for $SAMPLE_ID at $FASTQ_DIR/$SAMPLE_ID/"
    exit 1
  fi

  # 1) QC & trimming (fastp threaded)
  conda activate base
  fastp -i "$FQ1" -I "$FQ2" \
        -o ${SAMPLE_ID}_trimmed_1.fastq.gz -O ${SAMPLE_ID}_trimmed_2.fastq.gz \
        -h ${SAMPLE_ID}_fastp.html -j ${SAMPLE_ID}_fastp.json -w $THREADS

  # 2) Alignment (bwa mem threaded)
  conda activate bwa_env
  bwa mem -t "$THREADS" "$REF" ${SAMPLE_ID}_trimmed_1.fastq.gz ${SAMPLE_ID}_trimmed_2.fastq.gz > ${SAMPLE_ID}.sam

  # 3) BAM processing + RG (samtools threaded)
  conda activate samtools_env
  samtools view -Sb ${SAMPLE_ID}.sam | samtools sort -@ "$THREADS" -o ${SAMPLE_ID}_sorted.bam
  samtools addreplacerg -r "ID:${SAMPLE_ID}\tLB:lib1\tPL:ILLUMINA\tPU:unit1\tSM:${SAMPLE_ID}" \
      -o ${SAMPLE_ID}_rg.bam ${SAMPLE_ID}_sorted.bam
  samtools index ${SAMPLE_ID}_rg.bam

  # 4) MarkDuplicates (picard; single-threaded by default, use MarkDuplicatesSpark if available)
  conda activate picard-env
  # try MarkDuplicatesSpark if available (faster, uses multiple cores); fall back to MarkDuplicates
  if java -jar $PICARD_JAR MarkDuplicatesSpark --help >/dev/null 2>&1; then
    java -Xmx32g -jar $PICARD_JAR MarkDuplicatesSpark \
      -I ${SAMPLE_ID}_rg.bam -O $VCFBAM/${SAMPLE_ID}_dedup.bam \
      -M $VCFBAM/${SAMPLE_ID}_metrics.txt --conf 'spark.executor.cores='"$THREADS" || true
  else
    java -Xmx32g -jar $PICARD_JAR MarkDuplicates \
      -I ${SAMPLE_ID}_rg.bam \
      -O $VCFBAM/${SAMPLE_ID}_dedup.bam \
      -M $VCFBAM/${SAMPLE_ID}_metrics.txt \
      --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
  fi

  # 5) Base Recalibration (GATK) - GATK's BaseRecalibrator is single threaded; we give java options
  conda activate gatk_env
  export _JAVA_OPTIONS="-Xmx48g -XX:+UseParallelGC"
  $GATK BaseRecalibrator -R "$REF" -I $VCFBAM/${SAMPLE_ID}_dedup.bam --known-sites "$KNOWN" -O ${SAMPLE_ID}_recal.table
  $GATK ApplyBQSR       -R "$REF" -I $VCFBAM/${SAMPLE_ID}_dedup.bam --bqsr-recal-file ${SAMPLE_ID}_recal.table -O ${SAMPLE_ID}.recal.bam

  # Optional QC: Qualimap (unchanged)
  if [[ -n "$TARGET_BED" && -f "$TARGET_BED" ]]; then
    qualimap bamqc -bam ${SAMPLE_ID}.recal.bam -outdir ${SAMPLE_ID}_qualimap -outfile ${SAMPLE_ID}_qualimap.pdf -outformat PDF --paint-chromosome-limits --java-mem-size=8G || true
  fi

  # 6) HaplotypeCaller ? GVCF (use native-pair-hmm threads)
  $GATK HaplotypeCaller -R "$REF" -I ${SAMPLE_ID}.recal.bam \
        -O $VCFBAM/${SAMPLE_ID}.g.vcf.gz -ERC GVCF \
        --native-pair-hmm-threads "$THREADS"

  echo "========== ${SAMPLE_ID}: GVCF done =========="
}

# =========================
# STEP 8: Joint Genotyping
# =========================
function joint_genotyping() {
  echo "[JOINT] Combining GVCFs"
  local GVCF_ARGS=()
  for S in "${SAMPLES[@]}"; do
    GVCF_ARGS+=("-V" "$OUT_COHORT/OUTPUT/$S/VCF_and_BAM/${S}.g.vcf.gz")
  done
  $GATK CombineGVCFs -R "$REF" "${GVCF_ARGS[@]}" -O $OUT_COHORT/combined.g.vcf.gz

  echo "[JOINT] GenotypeGVCFs"
  $GATK GenotypeGVCFs -R "$REF" -V $OUT_COHORT/combined.g.vcf.gz -O $OUT_COHORT/cohort.raw.vcf.gz
}



# =========================
# STEP 9: Artifact & Hard Filters
# =========================
function filter_artifacts_and_hard() {
    cd "$OUT_COHORT"
    local INVCF="cohort.raw.vcf.gz"

   echo "[INFO] Step-9: Starting artifact & hard filters..."

    # ---- 1. Collect artifact metrics (optional) ----
    $GATK CollectSequencingArtifactMetrics \
        -I OUTPUT/${SAMPLES[0]}/VCF_and_BAM/${SAMPLES[0]}_dedup.bam \
        -O cohort_artifacts \
        -R "$REF" || true

    # ---- 2. GATK hard filters (flag only, records remain) ----
    $GATK VariantFiltration \
        -R "$REF" \
        -V "$INVCF" \
        --filter-expression "QD < 1.5" --filter-name "LowQD" \
        --filter-expression "FS > 100.0" --filter-name "HighFS" \
        --filter-expression "MQ < 30.0" --filter-name "LowMQ" \
        --filter-expression "SOR > 3.0" --filter-name "StrandBias" \
        --filter-expression "MQRankSum < -12.5"  --filter-name "LowMQRank" \
        --filter-expression "ReadPosRankSum < -8.0" --filter-name "LowReadPos" \
        -O cohort.step1.filtered.vcf.gz

    # ---- 3. Keep only PASS variants ----
    bcftools view -f PASS cohort.step1.filtered.vcf.gz -Oz -o cohort.step1.PASSonly.vcf.gz
    tabix -f cohort.step1.PASSonly.vcf.gz

    # ---- 4. Remove high-complexity regions if BED given ----
    if [[ -f "$HIGH_COMPLEXITY_BED" ]]; then
        bedtools intersect -header -v \
            -a cohort.step1.PASSonly.vcf.gz \
            -b "$HIGH_COMPLEXITY_BED" > cohort.step2.noHCR.vcf
        bgzip -f cohort.step2.noHCR.vcf && tabix -f cohort.step2.noHCR.vcf.gz
    else
        ln -sf cohort.step1.PASSonly.vcf.gz cohort.step2.noHCR.vcf.gz
        ln -sf cohort.step1.PASSonly.vcf.gz.tbi cohort.step2.noHCR.vcf.gz.tbi || true
    fi

    # ---- 5. Depth + allele-count + VAF filter using GATK (safe) ----
   $GATK VariantFiltration \
        -R "$REF" \
        -V cohort.step2.noHCR.vcf.gz \
        --genotype-filter-expression "DP < 3 || (AD != null && (AD[1] < 2 || (AD[0]+AD[1]) > 0 && AD[1]/(AD[0]+AD[1]) < 0.05))" \
        --genotype-filter-name "LowDP_AD_VAF" \
        -O cohort.step4.dp_ac_vaf.vcf.gz
    tabix -f cohort.step4.dp_ac_vaf.vcf.gz

    echo "[INFO] Step-9 filtering completed successfully!"
    echo "[INFO] Step-4 output file:"
    ls -lh cohort.step4.dp_ac_vaf.vcf.gz

    # Optional: quick sanity check
    echo "[INFO] Total variants after Step-4:"
    bcftools view -H cohort.step4.dp_ac_vaf.vcf.gz | wc -l
 
 # ---- 6. Remove sites with no ALT allele after genotype filtering ----
 bcftools +fill-tags cohort.step4.dp_ac_vaf.vcf.gz -- -t AC,AN \
    | bcftools view -i 'AC>0' -Oz -o cohort.step5.final.vcf.gz
tabix -f cohort.step5.final.vcf.gz

echo "[INFO] Final strict VCF:"
ls -lh cohort.step5.final.vcf.gz
}


# =========================       
# STEP 10: ANNOVAR Annotation
# =========================
function annovar_annotate() {
  cd "$OUT_COHORT"

  #input: cohort.step4.dp_ac_vaf.vcf.gz
  local INVCF="cohort.step5.final.vcf.gz"
  local OUT_PREFIX="cohort.annovar"

   #run table_annovar.pl directly on VCF
 perl "$ANNOVAR_DIR/table_annovar.pl" "$INVCF" "$HUMANDB" \
       -buildver hg38 \
       -out "$OUT_PREFIX" \
       -remove \
       -protocol refGene,clinvar,gnomad211_genome,avsnp150,dbnsfp42a \
       -operation g,f,f,f,f \
       -nastring . \
       -vcfinput

#Run ANNOVAR annotation with multiple databases


   #compress & index VCF output if present
  if [[ -f "${OUT_PREFIX}.hg38_multianno.vcf" ]]; then
    bcftools sort -Oz -o "${OUT_PREFIX}.hg38_multianno.sorted.vcf.gz" \
                  "${OUT_PREFIX}.hg38_multianno.vcf"
    tabix -f "${OUT_PREFIX}.hg38_multianno.sorted.vcf.gz"
    ln -sf "${OUT_PREFIX}.hg38_multianno.sorted.vcf.gz" cohort.annovar.vcf.gz
  fi
}

#============================================================
#  STEP 11: MERGE SAMPLE I'D WITH ANNOTATION FILE + memory-optimized + formatting + Excel-generation
#============================================================

# ---------- User-editable paths (already set to your locations) ---------- 
VCF="$OUT_COHORT/cohort.step5.final.vcf.gz"
ANNOVAR_TXT="$OUT_COHORT/cohort.annovar.hg38_multianno.txt"

TMPDIR="$OUT_COHORT/tmp_merge"
mkdir -p "$TMPDIR" "$OUT_COHORT/logs"
export TMPDIR

# ---------- derived paths ----------
VCF_MATRIX="$OUT_COHORT/vcf_GT_AD_matrix.tsv"
SAMPLE_LIST="$OUT_COHORT/vcf_sample_list.txt"

KEYED_ANNO="$TMPDIR/annovar_keyed.tsv"
KEYED_VCF="$TMPDIR/vcf_keyed.tsv"
SORTED_ANNO="$TMPDIR/annovar_keyed.sorted.tsv"
SORTED_VCF="$TMPDIR/vcf_keyed.sorted.tsv"

OUT_HUMAN="$OUT_COHORT/cohort.final_annotated_with_samples2.tsv"
OUT_HUMAN_GZ="$OUT_HUMAN.gz"
OUT_TECH="$OUT_COHORT/cohort.final_annotated_with_samples.technical2.tsv"
OUT_XLSX="$OUT_COHORT/cohort.final_annotated_with_samples2.xlsx"

# ---------- quick checks ----------
echo "[INFO] Checking inputs..."
if [[ ! -f "$VCF" ]]; then
  echo "[FATAL] VCF not found: $VCF" >&2
  exit 1
fi
if [[ ! -f "$ANNOVAR_TXT" ]]; then
  echo "[FATAL] ANNOVAR file not found: $ANNOVAR_TXT" >&2
  exit 1
fi
if [[ ! -x "$ANNOVAR/table_annovar.pl" && ! -f "$ANNOVAR/table_annovar.pl" ]]; then
  # make sure ANNOVAR path exists (not required to run annovar here, but good to warn)
  echo "[WARN] ANNOVAR directory or table_annovar.pl not found at $ANNOVAR (we don't run it here, but please verify)." >&2
fi

# ---------- ensure sample list exists ----------
if [[ ! -f "$SAMPLE_LIST" ]]; then
  echo "[INFO] Creating sample list from VCF..."
  bcftools query -l "$VCF" > "$SAMPLE_LIST"
fi
echo "[INFO] sample list: $(wc -l < "$SAMPLE_LIST") samples"

# ---------- ensure VCF_MATRIX exists ----------
if [[ ! -f "$VCF_MATRIX" ]]; then
  echo "[INFO] Extracting GT+AD matrix from VCF (this can take some minutes)..."
  # produce CHROM POS REF ALT [GT AD]...
  echo -e "CHROM\tPOS\tREF\tALT\t$(paste -sd '\t' $SAMPLE_LIST)" > "$VCF_MATRIX"

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT[\t%GT\t%AD]\n' "$VCF" >> "$VCF_MATRIX"

fi
echo "[INFO] VCF matrix present: $VCF_MATRIX  (lines: $(wc -l < "$VCF_MATRIX"))"

# ---------- Step: build keyed ANNOVAR stream (KEY\tcols...) ----------
echo "[STEP] Building keyed ANNOVAR file..."
python3 - <<PY
import csv, os, sys
in_annovar = os.path.expanduser("$ANNOVAR_TXT")
out_keyed = os.path.expanduser("$KEYED_ANNO")
with open(in_annovar, 'r', newline='') as inf, open(out_keyed, 'w', newline='') as outf:
    r = csv.reader(inf, delimiter='\t')
    w = csv.writer(outf, delimiter='\t', lineterminator='\n')
    header = next(r)
    header_l = [h.lower() for h in header]
    # find columns
    def find_any(names):
        for cand in names:
            if cand in header_l:
                return header_l.index(cand)
        return None
    c_idx = find_any(['chr','#chrom','chrom'])
    p_idx = find_any(['start','pos','position'])
    r_idx = find_any(['ref'])
    a_idx = find_any(['alt'])
    if None in (c_idx,p_idx,r_idx,a_idx):
        print("[ERROR] Could not find CHR/POS/REF/ALT columns in ANNOVAR header", file=sys.stderr)
        sys.exit(1)
    # write header: KEY + original header
    w.writerow(['KEY'] + header)
    for row in r:
        try:
            chrom = row[c_idx].strip()
            pos = row[p_idx].strip()
            ref = row[r_idx].strip()
            alt = row[a_idx].strip()
        except Exception:
            continue
        key = f"{chrom}:{pos}:{ref}:{alt}"
        w.writerow([key] + row)
print("[INFO] keyed ANNOVAR written to", out_keyed)
PY

# ---------- Step: build keyed VCF matrix ----------
echo "[STEP] Building keyed VCF matrix file..."
awk 'BEGIN{OFS="\t"} NR==1{print "KEY",$0; next} { key=$1":"$2":"$3":"$4; $1=""; $2=""; $3=""; $4=""; sub(/^\t/,""); print key,$0 }' "$VCF_MATRIX" > "$KEYED_VCF"
echo "[INFO] keyed VCF matrix lines: $(wc -l < "$KEYED_VCF")"

# ---------- Step: sort keyed files (disk-based) ----------
echo "[STEP] Sorting keyed files (disk-backed GNU sort)..."
# tune -S to available memory; default set to 40G which is safe on this job (64G total)
LC_ALL=C sort -T "$TMPDIR" -k1,1 "$KEYED_ANNO" -S 40G -o "$SORTED_ANNO"
LC_ALL=C sort -T "$TMPDIR" -k1,1 "$KEYED_VCF"  -S 40G -o "$SORTED_VCF"
echo "[INFO] sorted files ready."

# ---------- Step: streaming merge - produce two TSVs (human + technical) ----------
echo "[STEP] Streaming merge join (writing human & technical TSV)..."
python3 - <<PY
import csv, gzip, os, sys

sorted_anno = os.path.expanduser("$SORTED_ANNO")
sorted_vcf  = os.path.expanduser("$SORTED_VCF")
sample_list = os.path.expanduser("$SAMPLE_LIST")
out_human = os.path.expanduser("$OUT_HUMAN")
out_tech  = os.path.expanduser("$OUT_TECH")

# read sample names
with open(sample_list,'r') as sf:
    samples = [s.strip() for s in sf if s.strip()]

# open files
r_anno = open(sorted_anno, 'rt', newline='')
r_vcf  = open(sorted_vcf,  'rt', newline='')

reader_anno = csv.reader(r_anno, delimiter='\t')
reader_vcf  = csv.reader(r_vcf,  delimiter='\t')

# headers
anno_header = next(reader_anno)  # KEY + annovar columns
vcf_header  = next(reader_vcf)   # KEY + vcf matrix columns

ann_cols = anno_header[1:]  # drop KEY

# prepare writers
out_h = open(out_human, 'wt', newline='')
out_t = open(out_tech,  'wt', newline='')
w_h = csv.writer(out_h, delimiter='\t', lineterminator='\n')
w_t = csv.writer(out_t, delimiter='\t', lineterminator='\n')

# write headers: annovar columns + sample names
w_h.writerow(ann_cols + samples)
w_t.writerow(ann_cols + samples)

def parse_pairs(rest):
    # rest: list of alternating GT,AD,GT,AD...
    pairs = []
    for i in range(0, len(rest), 2):
        gt = rest[i] if i < len(rest) else '.'
        ad = rest[i+1] if (i+1) < len(rest) else '.'
        pairs.append((gt, ad))
    return pairs

def format_human(gt, ad):
    # return "." for 0/0 or missing; else GT (vaf; P%)
    if gt in (None, ".", "./."): return '.'
    if gt == "0/0" or gt == "0|0": return '.'
    if ad is None or ad == ".":
        return f"{gt} (.; .)"
    parts = ad.split(',')
    try:
        ref = int(parts[0]) if parts[0] != "." else 0
        alt = int(parts[1]) if len(parts)>1 and parts[1] != "." else 0
    except:
        return f"{gt} (.; .)"
    denom = ref + alt
    if denom == 0:
        return f"{gt} (.; .)"
    vaf = round(alt/denom,2)
    pct = f"{int(round(vaf*100))}%"
    return f"{gt} ({vaf}; {pct})"

def format_tech(gt, ad):
    if gt in (None, ".", "./."): return '.'
    if ad is None or ad == ".":
        return f"{gt}|.|."
    parts = ad.split(',')
    try:
        ref = int(parts[0]) if parts[0] != "." else 0
        alt = int(parts[1]) if len(parts)>1 and parts[1] != "." else 0
    except:
        return f"{gt}|{ad}|."
    denom = ref + alt
    vaf = round(alt/denom,3) if denom>0 else 0.0
    return f"{gt}|{ad}|{vaf}"

# prime rows
try:
    a_row = next(reader_anno)
except StopIteration:
    a_row = None
try:
    v_row = next(reader_vcf)
except StopIteration:
    v_row = None

processed = 0
while a_row is not None:
    a_key = a_row[0]
    # advance v_row until >= a_key
    while v_row is not None and v_row[0] < a_key:
        try:
            v_row = next(reader_vcf)
        except StopIteration:
            v_row = None
            break
    if v_row is not None and v_row[0] == a_key:
        # parse sample pairs
        rest = v_row[1:]
        pairs = parse_pairs(rest)
        sample_fields_h = []
        sample_fields_t = []
        for i,s in enumerate(samples):
            if i < len(pairs):
                gt,ad = pairs[i]
            else:
                gt,ad = '.', '.'
            sample_fields_h.append(format_human(gt,ad))
            sample_fields_t.append(format_tech(gt,ad))
        # write merged row
        out_ann = a_row[1:]  # annovar cols
        w_h.writerow(out_ann + sample_fields_h)
        w_t.writerow(out_ann + sample_fields_t)
        # advance both
        try:
            a_row = next(reader_anno)
        except StopIteration:
            a_row = None
        try:
            v_row = next(reader_vcf)
        except StopIteration:
            v_row = None
    else:
        # no vcf match: write annovar row with dots
        sample_dots_h = ['.'] * len(samples)
        sample_dots_t = ['.'] * len(samples)
        out_ann = a_row[1:]
        w_h.writerow(out_ann + sample_dots_h)
        w_t.writerow(out_ann + sample_dots_t)
        try:
            a_row = next(reader_anno)
        except StopIteration:
            a_row = None
    processed += 1
    if processed % 500000 == 0:
        print(f"[INFO] processed {processed} rows...", flush=True)

# close
out_h.close()
out_t.close()
r_anno.close()
r_vcf.close()
print("[INFO] streaming merge done. Wrote human and technical TSVs.")
PY

# ---------- gzip the human-readable TSV to save space ----------
echo "[INFO] Compressing human-readable TSV..."
gzip -f "$OUT_HUMAN"   # produces $OUT_HUMAN.gz

# ---------- Create Excel-friendly summary ----------
echo "[STEP] Creating Excel summary (compact when large)..."
python3 - <<PY
import pandas as pd, subprocess, os
OUT_TSV = os.path.expanduser("$OUT_TECH")   # use technical TSV for reliable parsing
OUT_XLSX = os.path.expanduser("$OUT_XLSX")

# count rows
rows = int(subprocess.check_output(["wc","-l", OUT_TSV]).split()[0])
print("[INFO] Rows in technical TSV:", rows)

if rows <= 400000:
    # safe to read whole file
    df = pd.read_csv(OUT_TSV, sep='\t', low_memory=False)
    df.to_excel(OUT_XLSX, index=False)
else:
    # create compact excel: pick core annotation columns + SAMPLES_PRESENT
    chunk_iter = pd.read_csv(OUT_TSV, sep='\t', chunksize=200000, low_memory=False)
    first = True
    compact = []
    for chunk in chunk_iter:
        if first:
            cols = list(chunk.columns)
            pick = []
            for cand in ['Chr','Start','End','Ref','Alt','Gene.refGene','ExonicFunc.refGene','Func.refGene','clinvar','gnomad211_genome_AF','avsnp150']:
                if cand in cols:
                    pick.append(cand)
            if len(pick) < 5:
                pick = cols[:10]
            sample_cols = [c for c in cols if c not in pick]
            first = False
        # SAMPLES_PRESENT: sample names with non-dot in any sample field (technical format uses GT|AD|VAF or .)
        def make_present(row):
            present = [c for c in sample_cols if row[c] != '.']
            return ','.join(present)
        chunk['SAMPLES_PRESENT'] = chunk.apply(make_present, axis=1)
        compact.append(chunk[pick + ['SAMPLES_PRESENT']])
        if sum(len(x) for x in compact) >= 500000:
            break
    compact_df = pd.concat(compact, ignore_index=True)
    compact_df.to_excel(OUT_XLSX, index=False)
print("[INFO] Excel summary written to", OUT_XLSX)
PY

echo "[END] $(date)"
echo "[OUTPUTS]"
ls -lh "$OUT_TECH" "$OUT_HUMAN_GZ" "$OUT_XLSX"

# =========================
# STEP 11: Post-ANNOVAR biological filters (cleaned VCF integrated)
# =========================

conda activate R_env

Rscript - << 'EOF'


library(data.table)
library(openxlsx)

#fUNCTIONAL FILTERING FOR NO INTRONIC REGION
 
input_file <- "cohort.final_annotated_with_samples.technical2.tsv"
variants <- fread(input_file)
 
cat("Input variants:", nrow(variants), "\n")
 
# Allowed functional classes
func_keep <- c("exonic", "splicing", "exonic;splicing", "splicing;exonic", "UTR5", "UTR3")
 
filtered_variants <- variants[ Func.refGene %in% func_keep ]
 
cat("Remaining after functional filter:", nrow(filtered_variants), "\n")
cat("Removed:", nrow(variants) - nrow(filtered_variants), "\n")
 
output_file <- "annotation.step1_functional.txt"
fwrite(filtered_variants, output_file, sep="\t")
 
cat("Output written:", output_file, "\n")


#REMOVE SYNONYMOUS KEEP ONLY NON-SYNONYMOUS

input_file <- "annotation.step1_functional.txt"
variants <- fread(input_file)
 
cat("Input variants:", nrow(variants), "\n")
 
# Keep only NON-synonymous coding variants
nonsyn_keep <- c(
  "nonsynonymous SNV",
  "stopgain",
  "stoploss",
  "frameshift insertion",
  "frameshift deletion",
  "frameshift block substitution",
  "nonframeshift insertion",
  "nonframeshift deletion",
  "nonframeshift block substitution"
)
 
filtered_variants <- variants[ ExonicFunc.refGene %in% nonsyn_keep ]
 
cat("Remaining after nonsynonymous filter:", nrow(filtered_variants), "\n")
cat("Removed:", nrow(variants) - nrow(filtered_variants), "\n")
 
output_file <- "annotation.step2_nonsynonymous.txt"
fwrite(filtered_variants, output_file, sep="\t")
cat("Output written:", output_file, "\n")


# ALLELE FREQUENCY AND POPULATION FREQUENCY FILTER

input_file <- "annotation.step2_nonsynonymous.txt"
variants <- fread(input_file)
 
cat("Input variants:", nrow(variants), "\n")
 
# SAFE Frequency columns
af_cols <- c(
  "AF", "AF_popmax", "AF_sas", "AF_afr", "AF_eas",
  "AF_nfe", "AF_amr", "AF_asj", "AF_oth"
)
 
# Replace missing with 0
for (col in af_cols) {
  if (col %in% names(variants)) {
    variants[[col]][variants[[col]] == "."] <- NA
	variants[[col]] <- as.numeric(variants[[col]])
  }
}
 
# Create MAX AF column
variants[, MAX_AF := do.call(pmax, c(.SD, list(na.rm = TRUE))), .SDcols = af_cols]
 
# Filter variants with max frequency < 1%
filtered_variants <- variants[MAX_AF < 0.01 | is.na(MAX_AF)]
 
cat("Remaining after AF < 0.01 filter:", nrow(filtered_variants), "\n")
cat("Removed:", nrow(variants) - nrow(filtered_variants), "\n")
 
output_file <- "annotation.step3_AF_filtered.txt"
fwrite(filtered_variants, output_file, sep="\t")
 
cat("Output written:", output_file, "\n")

#MUTATION TESTER SIFT-POLYPHEN FILTER 

input_file <- "annotation.step3_AF_filtered.txt"
variants <- fread(input_file)
 
cat("Input variants:", nrow(variants), "\n")
 
# Prediction columns
predictors <- c("SIFT_pred", "Polyphen2_HDIV_pred", "Polyphen2_HVAR_pred",
                "PROVEAN_pred", "MutationTaster_pred", "LRT_pred")
 
# Only use columns that actually exist
available_predictors <- predictors[predictors %in% colnames(variants)]
cat("Using predictors:", paste(available_predictors, collapse = ", "), "\n")
 
# Count how many predictors marked the variant as Deleterious ("D")
damaging_count <- apply(
  variants[, ..available_predictors],
  1,
  function(x) sum(toupper(trimws(x)) == "D", na.rm = TRUE)
)
 
# Important genes to ALWAYS keep
important_genes <- c("LMNA", "POMT2", "NEB", "TPM3")
 
# Gene column
gene_col <- if ("Gene.refGene" %in% colnames(variants)) "Gene.refGene" else stop("Gene.refGene not found.")
 
# Final filtering logic
filtered_variants <- variants[
  damaging_count >= 3 |
  grepl("splicing", Func.refGene, ignore.case = TRUE) |
  get(gene_col) %in% important_genes
]
 
# Output
output_file <- "annotation.step4_Mutation.txt"
fwrite(filtered_variants, output_file, sep="\t")
 
cat("Applied STRICT pathogenicity filter (≥3 damaging predictors)\n")
cat("Splicing and important genes retained\n")
cat("Input:", nrow(variants), "\n")
cat("Remaining:", nrow(filtered_variants), "\n")
cat("Output written:", output_file, "\n")



# REPORTED GENES FILTER 

input_file <- "annotation.step4_Mutation.txt"
cmd_gene_file <- "CMD_genes2.txt"
cm_gene_file <- "CM_genes2.txt"
 
# Load data
variants <- fread(input_file)
cmd_genes <- unique(trimws(readLines(cmd_gene_file)))
cm_genes <- unique(trimws(readLines(cm_gene_file)))
reported_genes <- unique(c(cmd_genes, cm_genes))
 
# Important genes
important_genes <- c("LMNA", "POMT2", "NEB", "TPM3")
 
# Check gene column
gene_col <- if ("Gene.refGene" %in% colnames(variants)) "Gene.refGene" else "genename"
 
# Keep variants in reported or important genes or containing splicing
filtered_variants <- variants[
  get(gene_col) %in% reported_genes |
  get(gene_col) %in% important_genes |
  grepl("splic", variants$Func.refGene, ignore.case = TRUE)
]
 
output_file <- "annotation.step5_reported.txt"
fwrite(filtered_variants, output_file)
 
cat("✅ Kept reported, important, and splicing gene variants.\n")
cat("Remaining variants:", nrow(filtered_variants), "\n")
cat("Filtered data written to:", output_file, "\n")

data <- fread("annotation.step5_reported.txt")
write.xlsx(data, "annotation.step5_reported.xlsx")

EOF


# =========================
# MAIN
# =========================
download_fastqs

for S in "${SAMPLES[@]}"; do
  process_sample "$S"
done

joint_genotyping
filter_artifacts_and_hard
annovar_annotate


echo "[DONE] Full multi-sample WES cohort pipeline complete."
echo "Outputs:"
echo " - $OUT_COHORT/cohort.raw.vcf.gz (raw joint calls)"
echo " - $OUT_COHORT/cohort.annovar.filtered.clinvar.vcf.gz (final annotated - clinvar annotated)"
echo " - $OUT_COHORT/cohort.CMD_prioritized.vcf.gz  (CMD gene-panel)"
echo " - $OUT_COHORT/cohort.CM_prioritized.vcf.gz   (CM gene-panel)"


