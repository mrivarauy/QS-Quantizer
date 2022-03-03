
# Set base dir
BASE="/path/to/base/directory" #base directory

# Set up input data
REF_DIR="/path/to/reference/directory/ref_dir" #reference directory (it contains human genome GRCh38 reference)
INPUT_DIR="${BASE}/input" #input directory (it contains fastq files )
REF="GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz" #human genome GRCh38 reference file
READS1="${INPUT_DIR}/GM24149_1_Guppy_4.2.2_prom.fastq.gz" #fastq file 1
READS2="${INPUT_DIR}/GM24149_2_Guppy_4.2.2_prom.fastq.gz" #fastq file 2
READS3="${INPUT_DIR}/GM24149_3_Guppy_4.2.2_prom.fastq.gz" #fastq file 3

# Set the number of CPUs to use
THREADS="28" #number of threads to use

# Set up output directory
OUTPUT_DIR="${BASE}/output" #output directory
PREFIX="files_prefix" #files prefix

mkdir -p "${OUTPUT_DIR}"

# Make the bam file

minimap2 -ax map-ont -t "${THREADS}" "${REF_DIR}"/"${REF}" "${READS1}" "${READS2}" "${READS3}" | \
  samtools view -hb -F 0x904 | \
  samtools sort -@"${THREADS}" -o "${INPUT_DIR}"/"${PREFIX}".sorted.bam -
samtools index -@"${THREADS}" "${INPUT_DIR}"/"${PREFIX}".sorted.bam
