# Set base dir
BASE="/path/to/base/directory/variant-calling" #path to base drirectory

# Set up input data
REF_DIR="/media/Disco9/datos_DeepVariant/ref_dir" #reference directory (it contains human genome GRCh38 reference file and reference vcf file)
INPUT_DIR="${BASE}/input" # input directory (it contains fastq files and bam files of reads mapped to reference)
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna" #human genome GRCh38 reference file

# Set the number of CPUs to use
THREADS="15"

# Set up output directory
OUTPUT_DIR="${BASE}/output_directory"
PREFIX="files_prefix"

# Run the pipeline

sudo docker run --ipc=host \
-v "${INPUT_DIR}":"${INPUT_DIR}" \
-v "${OUTPUT_DIR}":"${OUTPUT_DIR}" \
-v "${REF_DIR}":"${REF_DIR}" \
kishwars/pepper_deepvariant:r0.4 \
run_pepper_margin_deepvariant call_variant \
-b "${INPUT_DIR}/${PREFIX}.bam" \
-f "${REF_DIR}"/"${REF}" \
-o "${OUTPUT_DIR}" \
-p "${PREFIX}" \
-t "${THREADS}" \
--ont
