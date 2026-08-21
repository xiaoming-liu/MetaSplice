#!/bin/bash
# Quick test for MetaSplice scoring pipeline
# Run from the MetaSplice directory

set -e

echo "=== MetaSplice Test ==="

# Create test VCF
cat > /tmp/metasplice_test.vcf <<'EOF'
##fileformat=VCFv4.2
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
1	861407	.	C	T	.	.	.
7	117199644	.	T	G	.	.	.
17	41246481	.	G	A	.	.	.
EOF

echo "Test VCF: 3 variants"

# Run with minimal features (requires FASTA + GTF)
if [ -z "$1" ]; then
    echo "Usage: bash test.sh /path/to/data"
    echo "  Data directory should contain:"
    echo "    hg19.fa (indexed)"
    echo "    gencode.v19.annotation.gtf.gz"
    echo "    dbNSFP5.3_gene.gz (optional)"
    exit 1
fi

DATA=$1

python score.py \
    --input /tmp/metasplice_test.vcf \
    --output /tmp/metasplice_test_scores.tsv \
    --fasta "${DATA}/hg19.fa" \
    --gtf "${DATA}/gencode.v19.annotation.gtf.gz" \
    --dbnsfp-gene "${DATA}/dbNSFP5.3_gene.gz" \
    --dann-dir "${DATA}/DANN" \
    --eigen-dir "${DATA}/Eigen"

echo ""
echo "=== Results ==="
cat /tmp/metasplice_test_scores.tsv
echo ""
echo "=== Test complete ==="

echo ""
echo "=== Eigen two-block regression tests ==="
python "$(dirname "$0")/test_eigen_regression.py" "${DATA}"
