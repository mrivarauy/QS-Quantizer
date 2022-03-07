# Set base dir
QUANT=original #change for quantized data. i.e. 2bin or 8bin depending on current quantization
COV=20 #set for working coverage. i.e. for 20X coverage "20" for 60X coverage "60"
BASE=./${QUANT}

# Set up input data
REF_DIR=./ref_dir
INPUT_DIR="${BASE}/input"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
TRUTH_VCF="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# Set the number of CPUs to use
THREADS="15"

# Set up output directory
OUTPUT_DIR="${BASE}/${COV}"
PREFIX="HG003_${QUANT}_${COV}"

# Run the pipeline

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${REF_DIR}":"${REF_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${PREFIX}.sorted.bam" \
-f "${REF_DIR}"/"${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}" \
-t "${THREADS}" \
--ont


#Evaluate performance
# Run hap.py
sudo docker run --ipc=host \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${REF_DIR}":"${REF_DIR}" \
pkrusche/hap.py:latest \
/opt/hap.py/bin/hap.py \
${REF_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${PREFIX}.vcf \
-f "${REF_DIR}/${TRUTH_BED}" \
-r "${REF_DIR}/${REF}" \
-o "${OUTPUT_DIR}/happy.output" \
--pass-only \
--engine=vcfeval \
--threads="${THREADS}"
