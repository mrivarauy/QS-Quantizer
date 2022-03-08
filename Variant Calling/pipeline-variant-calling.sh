#Parameters:
#   $1: quantizer used in non-run regions: Q0 (no quantization), Q2, Q4, Q8, F10
#   $2: quantizer used in run regions: Q0 (no quantization), Q2, Q4, Q8, F10


# Quantizer
QUANT="$1$2"

# Set base dir
BASE="./bins/${QUANT}"
# Set up input data
REF_DIR="./ref_dir"
INPUT_DIR="${BASE}/input"
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"
TRUTH_VCF="HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
TRUTH_BED="HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# Set up output directory
OUTPUT_DIR="${BASE}/output"
PREFIX="HG003_${QUANT}_"

#names of quantized files
READS1="${INPUT_DIR}/HG003_1_${QUANT}.fastq"
READS2="${INPUT_DIR}/HG003_2_${QUANT}.fastq"
READS3="${INPUT_DIR}/HG003_3_${QUANT}.fastq"

LOGFILE="./log${QUANT}.txt"
printf "creating quantized files\n" > "${LOGFILE}"

#create quantized files
cd ./
mkdir -p "${INPUT_DIR}"
python ../quantizer.py ./original-fastq/GM24149_1_Guppy_4.2.2_prom.fastq.gz $1 $2 "${READS1}" &
python ../quantizer.py ./original-fastq/GM24149_2_Guppy_4.2.2_prom.fastq.gz $1 $2 "${READS2}" &
python ../quantizer.py ./original-fastq/GM24149_3_Guppy_4.2.2_prom.fastq.gz $1 $2 "${READS3}" &
wait

# Set the number of CPUs to use
THREADS="30"

mkdir -p "${OUTPUT_DIR}"


printf "creating BAM file\n" >> "${LOGFILE}"

# Make the bam file within a conda environment
conda activate human-env
minimap2 -ax map-ont -t "${THREADS}" "${REF_DIR}"/"${REF}" "${READS1}" "${READS2}" "${READS3}" | \
  samtools view -hb -F 0x904 | \
  samtools sort -@"${THREADS}" -o "${INPUT_DIR}"/"${PREFIX}".sorted.bam -
samtools index -@"${THREADS}" "${INPUT_DIR}"/"${PREFIX}".sorted.bam
conda deactivate

printf "executing variant calling pipeline 90x\n" >> "${LOGFILE}"

# Run the variant calling pipeline
singularity exec \
-B "${INPUT_DIR}" \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
pepper_deepvariant-r0.4.simg \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${PREFIX}.sorted.bam" \
-f "${REF_DIR}"/"${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}" \
-t "${THREADS}" \
--ont


# Subsampling
SAMPLING="90x"
# Set up report directory
REPORT_DIR="${BASE}/report/${SAMPLING}"

mkdir -p "${REPORT_DIR}"

printf "generating report 90x\n" >> "${LOGFILE}"

# Run hap.py
singularity exec \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
-B "${REPORT_DIR}" \
hap.py-latest.simg \
/opt/hap.py/bin/hap.py \
${REF_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${PREFIX}.vcf.gz \
-f "${REF_DIR}/${TRUTH_BED}" \
-r "${REF_DIR}/${REF}" \
-o "${REPORT_DIR}/_happy.output" \
--pass-only \
--engine=vcfeval \
--threads=10


#######################################################################
# 50x
#######################################################################

# Subsampling
SAMPLING="50x"
DS_BAM="${INPUT_DIR}"/"${PREFIX}"_"${SAMPLING}".sorted.bam
printf "downsampling to 50x\n" >> "${LOGFILE}"
conda activate human-env
samtools view -s 0.55 -@"${THREADS}" -b "${INPUT_DIR}"/"${PREFIX}".sorted.bam > "${DS_BAM}"
samtools index -@"${THREADS}" "${DS_BAM}"
conda deactivate

printf "executing variant calling pipeline 50x\n" >> "${LOGFILE}"

# Run the variant calling pipeline
singularity exec \
-B "${INPUT_DIR}" \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
pepper_deepvariant-r0.4.simg \
run_pepper_margin_deepvariant call_variant \
-b "${DS_BAM}" \
-f "${REF_DIR}"/"${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}_${SAMPLING}" \
-t "${THREADS}" \
--ont

# Set up report directory
REPORT_DIR="${BASE}/report/${SAMPLING}"
mkdir -p "${REPORT_DIR}"


printf "generating report 50x\n" >> "${LOGFILE}"
# Run hap.py
singularity exec \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
-B "${REPORT_DIR}" \
hap.py-latest.simg \
/opt/hap.py/bin/hap.py \
${REF_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${PREFIX}_${SAMPLING}.vcf.gz \
-f "${REF_DIR}/${TRUTH_BED}" \
-r "${REF_DIR}/${REF}" \
-o "${REPORT_DIR}/_happy.output" \
--pass-only \
--engine=vcfeval \
--threads=10


#######################################################################
# 40x
#######################################################################

# Subsampling
SAMPLING="40x"
DS_BAM="${INPUT_DIR}"/"${PREFIX}"_"${SAMPLING}".sorted.bam
printf "downsampling to ${SAMPLING}\n" >> "${LOGFILE}"
conda activate human-env
samtools view -s 0.44 -@"${THREADS}" -b "${INPUT_DIR}"/"${PREFIX}".sorted.bam > "${DS_BAM}"
samtools index -@"${THREADS}" "${DS_BAM}"
conda deactivate

printf "executing variant calling pipeline ${SAMPLING}\n" >> "${LOGFILE}"

# Run the variant calling pipeline
singularity exec \
-B "${INPUT_DIR}" \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
pepper_deepvariant-r0.4.simg \
run_pepper_margin_deepvariant call_variant \
-b "${DS_BAM}" \
-f "${REF_DIR}"/"${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}_${SAMPLING}" \
-t "${THREADS}" \
--ont

# Set up report directory
REPORT_DIR="${BASE}/report/${SAMPLING}"
mkdir -p "${REPORT_DIR}"


printf "generating report ${SAMPLING}\n" >> "${LOGFILE}"
# Run hap.py
singularity exec \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
-B "${REPORT_DIR}" \
hap.py-latest.simg \
/opt/hap.py/bin/hap.py \
${REF_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${PREFIX}_${SAMPLING}.vcf.gz \
-f "${REF_DIR}/${TRUTH_BED}" \
-r "${REF_DIR}/${REF}" \
-o "${REPORT_DIR}/_happy.output" \
--pass-only \
--engine=vcfeval \
--threads=10


#######################################################################
# 30x
#######################################################################

# Subsampling
SAMPLING="30x"
DS_BAM="${INPUT_DIR}"/"${PREFIX}"_"${SAMPLING}".sorted.bam
printf "downsampling to ${SAMPLING}\n" >> "${LOGFILE}"
conda activate human-env
samtools view -s 0.33 -@"${THREADS}" -b "${INPUT_DIR}"/"${PREFIX}".sorted.bam > "${DS_BAM}"
samtools index -@"${THREADS}" "${DS_BAM}"
conda deactivate

printf "executing variant calling pipeline ${SAMPLING}\n" >> "${LOGFILE}"

# Run the variant calling pipeline
singularity exec \
-B "${INPUT_DIR}" \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
pepper_deepvariant-r0.4.simg \
run_pepper_margin_deepvariant call_variant \
-b "${DS_BAM}" \
-f "${REF_DIR}"/"${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}_${SAMPLING}" \
-t "${THREADS}" \
--ont

# Set up report directory
REPORT_DIR="${BASE}/report/${SAMPLING}"
mkdir -p "${REPORT_DIR}"


printf "generating report ${SAMPLING}\n" >> "${LOGFILE}"
# Run hap.py
singularity exec \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
-B "${REPORT_DIR}" \
hap.py-latest.simg \
/opt/hap.py/bin/hap.py \
${REF_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${PREFIX}_${SAMPLING}.vcf.gz \
-f "${REF_DIR}/${TRUTH_BED}" \
-r "${REF_DIR}/${REF}" \
-o "${REPORT_DIR}/_happy.output" \
--pass-only \
--engine=vcfeval \
--threads=10


#######################################################################
# 20x
#######################################################################

# Subsampling
SAMPLING="20x"
DS_BAM="${INPUT_DIR}"/"${PREFIX}"_"${SAMPLING}".sorted.bam
printf "downsampling to ${SAMPLING}\n" >> "${LOGFILE}"
conda activate human-env
samtools view -s 0.22 -@"${THREADS}" -b "${INPUT_DIR}"/"${PREFIX}".sorted.bam > "${DS_BAM}"
samtools index -@"${THREADS}" "${DS_BAM}"
conda deactivate

printf "executing variant calling pipeline ${SAMPLING}\n" >> "${LOGFILE}"

# Run the variant calling pipeline
singularity exec \
-B "${INPUT_DIR}" \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
pepper_deepvariant-r0.4.simg \
run_pepper_margin_deepvariant call_variant \
-b "${DS_BAM}" \
-f "${REF_DIR}"/"${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}_${SAMPLING}" \
-t "${THREADS}" \
--ont

# Set up report directory
REPORT_DIR="${BASE}/report/${SAMPLING}"
mkdir -p "${REPORT_DIR}"


printf "generating report ${SAMPLING}\n" >> "${LOGFILE}"
# Run hap.py
singularity exec \
-B "${OUTPUT_DIR}" \
-B "${REF_DIR}" \
-B "${REPORT_DIR}" \
hap.py-latest.simg \
/opt/hap.py/bin/hap.py \
${REF_DIR}/${TRUTH_VCF} \
${OUTPUT_DIR}/${PREFIX}_${SAMPLING}.vcf.gz \
-f "${REF_DIR}/${TRUTH_BED}" \
-r "${REF_DIR}/${REF}" \
-o "${REPORT_DIR}/_happy.output" \
--pass-only \
--engine=vcfeval \
--threads=10

###################################################################################
printf "compressing files\n" >> "${LOGFILE}"
gzip "${READS1}" &
gzip "${READS2}" &
gzip "${READS3}" &
wait
ls -l "${INPUT_DIR}"/*.gz > "${BASE}/report/compressed_sizes.txt"
