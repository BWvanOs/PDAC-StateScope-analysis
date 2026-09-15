#!/usr/bin/env bash
set -euo pipefail

# rnaseq_star_align.sh
# STAR RNA-seq alignment pipeline (2-pass) for paired-end FASTQs.
# Defaults input directory to the folder containing this script—so you can drop
# the script right into your FASTQ folder and run it with just -g and -a.

# Defaults
threads=8
sjdb_overhang=100
outdir=""
index_dir=""
sample_glob="*R1*.f*q*"
single_end=0
dry_run=0

# Compute script directory (default FASTQ dir)
script_path="$(realpath "$0")"
script_dir="$(dirname "$script_path")"
fastq_dir_default="$script_dir"
fastq_dir="$fastq_dir_default"

# Parse args
ARGS=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    -i|--input) fastq_dir="$2"; shift 2;;
    -g|--genome) genome_fa="$2"; shift 2;;
    -a|--gtf) gtf="$2"; shift 2;;
    -o|--outdir) outdir="$2"; shift 2;;
    -t|--threads) threads="$2"; shift 2;;
    --index-dir) index_dir="$2"; shift 2;;
    --sjdb-overhang) sjdb_overhang="$2"; shift 2;;
    --sample-glob) sample_glob="$2"; shift 2;;
    --single-end) single_end=1; shift 1;;
    --dry-run) dry_run=1; shift 1;;
    -h|--help) usage; exit 0;;
    *) ARGS+=("$1"); shift 1;;
  esac
done
set -- "${ARGS[@]:-}"

# Validate required args
if [[ -z "${genome_fa:-}" || -z "${gtf:-}" ]]; then
  usage; exit 1
fi

# Defaults dependent on fastq_dir
if [[ -z "${outdir}" ]]; then
  outdir="${fastq_dir%/}/STAR_out"
fi
if [[ -z "${index_dir}" ]]; then
  index_dir="${outdir%/}/STAR_index"
fi

# Resolve absolute paths
fastq_dir="$(realpath -m "$fastq_dir")"
outdir="$(realpath -m "$outdir")"
index_dir="$(realpath -m "$index_dir")"
genome_fa="$(realpath -m "$genome_fa")"
gtf="$(realpath -m "$gtf")"

echo "FASTQs : $fastq_dir (defaulted to script folder if -i not provided)"
echo "Genome : $genome_fa"
echo "GTF    : $gtf"
echo "Outdir : $outdir"
echo "Index  : $index_dir"
echo "Threads: $threads"
echo "Overhang: $sjdb_overhang"
echo "Single-end: $single_end"
echo

# Check tools
need() { command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found in PATH."; exit 2; }; }
need STAR
need samtools
zcat_cmd="zcat"; command -v pigz >/dev/null 2>&1 && zcat_cmd="pigz -dc"

mkdir -p "$outdir"

# Build STAR index if missing
maybe_build_index() {
  if [[ -f "${index_dir}/Genome" && -s "${index_dir}/Genome" ]]; then
    echo "STAR index found at ${index_dir}"
    return
  fi
  echo "Building STAR index at ${index_dir} ..."
  mkdir -p "$index_dir"
  set -x
  STAR \
    --runThreadN "$threads" \
    --runMode genomeGenerate \
    --genomeDir "$index_dir" \
    --genomeFastaFiles "$genome_fa" \
    --sjdbGTFfile "$gtf" \
    --sjdbOverhang "$sjdb_overhang"
  { set +x; } 2>/dev/null
}

# Derive sample name from R1 path
sample_name_from_r1() {
  local r1="$1" base s
  base="$(basename "$r1")"
  s="${base%.gz}"; s="${s%.gzip}"; s="${s%.fastq}"; s="${s%.fq}"
  s="${s/_R1_/_}"; s="${s/_R1/-}"; s="${s/.R1/-}"
  s="${s/_1_/_}";  s="${s/_1/-}";  s="${s/.1/-}"
  s="${s%%[-_ ]}"
  echo "$s"
}

# Find matching R2 given R1
guess_r2_for_r1() {
  local r1="$1" r2=""
  if   [[ "$r1" =~ \.R1\.f(ast)?q(\.gz)?$ ]]; then r2="${r1/.R1./.R2.}"
  elif [[ "$r1" =~ _R1 ]]; then r2="${r1/_R1/_R2}"
  elif [[ "$r1" =~ _1\.f ]]; then r2="${r1/_1./_2.}"
  elif [[ "$r1" =~ \.1\.f ]]; then r2="${r1/.1./.2.}"
  fi
  [[ -n "$r2" && -f "$r2" ]] && echo "$r2" || echo ""
}

# Collect input files safely (supports spaces)
mapfile -d '' r1_files < <(find "$fastq_dir" -maxdepth 1 -type f -iname "$sample_glob" -print0 || true)
if (( ${#r1_files[@]} == 0 )); then
  echo "No R1 files found under: $fastq_dir (glob: $sample_glob)"
  exit 3
fi

maybe_build_index

# Process each sample
for r1 in "${r1_files[@]}"; do
  r1="${r1%$'\n'}"
  sample="$(sample_name_from_r1 "$r1")"
  samp_out="${outdir}/${sample}"
  mkdir -p "$samp_out"

  r2=""
  if (( single_end == 0 )); then
    r2="$(guess_r2_for_r1 "$r1")"
    if [[ -z "$r2" ]]; then
      echo "WARN: could not find R2 for $r1 — skipping (or use --single-end)"; continue
    fi
  fi

  echo "=== Sample: $sample ==="
  echo "R1: $r1"
  [[ $single_end -eq 0 ]] && echo "R2: $r2"

  # Decide decompression
  in_cmd=()
  if [[ "$r1" =~ \.gz(ip)?$ ]]; then
    in_cmd=(--readFilesCommand $zcat_cmd)
  fi

  # Run STAR (2-pass)
  set -x
  if (( single_end == 0 )); then
    STAR \
      --runThreadN "$threads" \
      --twopassMode Basic \
      --genomeDir "$index_dir" \
      --readFilesIn "$r1" "$r2" \
      "${in_cmd[@]}" \
      --outFileNamePrefix "${samp_out}/" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM GeneCounts \
      --outSAMattributes NH HI AS nM XS
  else
    STAR \
      --runThreadN "$threads" \
      --twopassMode Basic \
      --genomeDir "$index_dir" \
      --readFilesIn "$r1" \
      "${in_cmd[@]}" \
      --outFileNamePrefix "${samp_out}/" \
      --outSAMtype BAM SortedByCoordinate \
      --quantMode TranscriptomeSAM GeneCounts \
      --outSAMattributes NH HI AS nM XS
  fi
  { set +x; } 2>/dev/null

  # Index BAM
  bam="${samp_out}/Aligned.sortedByCoord.out.bam"
  if [[ -f "$bam" ]]; then
    samtools index -@ "$threads" "$bam"
  fi

  echo "Done: $sample → $samp_out"
done

echo "All done. Outputs in: $outdir"

