# Set base dir
QUANT=original #change for quantized data. i.e. 2bin or 8bin depending on current quantization
BASE=./${QUANT}

# Set up input data
REF_DIR=./ref_dir
INPUT_DIR=${BASE}/input
REF=GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
READS1=${INPUT_DIR}/GM24149_1_Guppy_4.2.2_prom.fastq.gz
READS2=${INPUT_DIR}/GM24149_2_Guppy_4.2.2_prom.fastq.gz
READS3=${INPUT_DIR}/GM24149_3_Guppy_4.2.2_prom.fastq.gz

# Set the number of CPUs to use
THREADS=28

# Set up output directory
OUTPUT_DIR=${BASE}/output
PREFIX=HG003_${QUANT}

mkdir -p ${OUTPUT_DIR}

# Make the bam file

minimap2 -ax map-ont -t ${THREADS} ${REF_DIR}/${REF} ${READS1} ${READS2} ${READS3} | \
  samtools view -hb -F 0x904 | \
  samtools sort -@${THREADS} -o ${INPUT_DIR}/${PREFIX}.sorted.bam -
samtools index -@${THREADS} ${INPUT_DIR}/${PREFIX}.sorted.bam
