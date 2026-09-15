#!/usr/bin/env bash
# Gene-level read counting of RNA-seq BAMs with featureCounts (exon features,
# unstranded, fragment counting). Column headers are reduced to sample names.
#
# Usage: rnaseq/03_count.sh <out.tsv> <sample> [<sample> ...]
# Input: ${MP_DATA}/rnaseq/star/<sample>/<sample>.Aligned.sortedByCoord.out.bam
set -euo pipefail
: "${MP_DATA:?source config.env first}"
: "${MP_REFS:?source config.env first}"

OUT="$1"; shift
GTF="$MP_REFS/gencode.vM25.annotation.gtf"
THREADS="${MP_THREADS:-16}"

bams=()
for s in "$@"; do
    bams+=("$MP_DATA/rnaseq/star/$s/${s}.Aligned.sortedByCoord.out.bam")
done
mkdir -p "$(dirname "$OUT")"

featureCounts \
    -T "$THREADS" \
    -a "$GTF" -F GTF \
    -t exon -g gene_name \
    -p --countReadPairs \
    -s 0 \
    -B -C --primary \
    -o "$OUT" \
    "${bams[@]}"

python3 - "$OUT" <<'PY'
import os, re, sys
path = sys.argv[1]
with open(path) as fh:
    lines = fh.readlines()
hdr = lines[1].rstrip("\n").split("\t")
hdr[6:] = [re.sub(r"\.Aligned\.sortedByCoord\.out\.bam$", "", os.path.basename(h)) for h in hdr[6:]]
lines[1] = "\t".join(hdr) + "\n"
with open(path, "w") as fh:
    fh.writelines(lines)
PY
