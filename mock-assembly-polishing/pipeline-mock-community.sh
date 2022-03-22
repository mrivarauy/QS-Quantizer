#!/usr/bin/env bash

#Parameters:
#   $1: quantizer used in non-run regions: Q0 (no quantization), Q2, Q4, Q8, F10
#   $2: quantizer used in run regions: Q0 (no quantization), Q2, Q4, Q8, F10

# Quantizer
PREFIX="$1$2"

#Set up input data
BASE="./bins/${PREFIX}" 
REF_DIR_ILL="./ref_dir" #assembly reference, used for performance evaluation
INPUT_DIR="${BASE}/input"

READS_QUAN="${INPUT_DIR}/Zymo_${PREFIX}.fastq"
READS_ORI="./original-fastq/Zymo-GridION-EVEN-3Peaks-R103-merged.fq.gz" #path to fastq.gz files
MEDAKA_MODEL="./medaka-model/r103_min_high_g345" #used medaka model
MP_PARAMS="../marginPolish/params/allParams.np.microbial.r103g324.json" #marginPolish parameters for microbial data and specific ONT pore, available in https://github.com/UCSC-nanopore-cgl/MarginPolish/tree/master/params
HELEN_MODEL="../helen-models/HELEN_r103_guppy_microbial.pkl" # HELEN trained model for used ONT pore

#Set up output directory
OUT_DIR="${BASE}/output"
METAQUAST_OUTDIR="${OUT_DIR}/mQ_report" #output directory name for metaQuast evaluation

NPROC=32 #number of threads to use

ASSEMBLY=${OUT_DIR}/assembly-${PREFIX}/raw-${PREFIX}.fasta #assembly location

## flags, change to false to not run
flyeflag=true
raconflag=true
medakaflag=true
marginflag=true
helenflag=true
metaquastflag=true

# -----------------------------


#create quantized files

mkdir -p "${INPUT_DIR}"
python ../quantizer.py ${READS_ORI} $1 $2 "${READS_QUAN}" &
wait


##assembly using flye
if [ $flyeflag == true ]; then
	conda activate flye-env
	flye --nano-raw ${READS_QUAN} --meta --iterations 0 --threads ${NPROC} --out-dir ${OUT_DIR}/assembly-${PREFIX}
	conda deactivate
	infoseq --only --length --name ${OUT_DIR}/assembly-${PREFIX}/assembly.fasta | awk '{print $2, $1}' | 
		sort -nr | head -n8 | cut -d' ' -f2 > ${OUT_DIR}/assembly-${PREFIX}/best_contigs
	fastaUtils.pl -u ${OUT_DIR}/assembly-${PREFIX}/assembly.fasta | 
		grep --no-group-separator -A 1 -w -f ${OUT_DIR}/assembly-${PREFIX}/best_contigs > ${OUT_DIR}/assembly-${PREFIX}/raw-${PREFIX}.fasta
	ASSEMBLY=${OUT_DIR}/assembly-${PREFIX}/raw-${PREFIX}.fasta
fi

conda activate mock-env
##evaluation using metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --no-icarus --fragmented --min-identity 90 --min-contig 5000 \
        --threads ${NPROC} -r ${REF_DIR_ILL} -o ${METAQUAST_OUTDIR}_raw ${ASSEMBLY}
fi

##polishing using racon
if [ $raconflag == true ]; then
    minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_QUAN} > ${OUT_DIR}/map1-${PREFIX}.sam
    racon -t ${NPROC} ${READS_QUAN} ${OUT_DIR}/map1-${PREFIX}.sam ${ASSEMBLY} > ${OUT_DIR}/racon_r1-${PREFIX}.fasta
    minimap2 -ax map-ont -t ${NPROC} ${OUT_DIR}/racon_r1-${PREFIX}.fasta ${READS_QUAN} > ${OUT_DIR}/map2-${PREFIX}.sam
    racon -t ${NPROC} ${READS_QUAN} ${OUT_DIR}/map2-${PREFIX}.sam ${OUT_DIR}/racon_r1-${PREFIX}.fasta > ${OUT_DIR}/racon_r2-${PREFIX}.fasta
fi

##evaluation using metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --no-icarus --fragmented --min-identity 90 --min-contig 5000 \
        --threads ${NPROC} -r ${REF_DIR_ILL} -o ${METAQUAST_OUTDIR}_racon ${OUT_DIR}/racon_r2-${PREFIX}.fasta
fi

## polishing using medaka (requires previous racon polishing)
if [ $medakaflag == true ]; then
    medaka_consensus -i ${READS_QUAN} -d ${OUT_DIR}/racon_r1-${PREFIX}.fasta -o ${OUT_DIR}/r1-medaka-${PREFIX} -m ${MEDAKA_MODEL} -t ${NPROC}
    mv ${OUT_DIR}/r1-medaka-${PREFIX}/consensus.fasta ${OUT_DIR}/r1-medaka-${PREFIX}.fasta && rm -rf ${OUT_DIR}/r1-medaka-${PREFIX}
    medaka_consensus -i ${READS_QUAN} -d ${OUT_DIR}/racon_r2-${PREFIX}.fasta -o ${OUT_DIR}/r2-medaka-${PREFIX} -m ${MEDAKA_MODEL} -t ${NPROC}
    mv ${OUT_DIR}/r2-medaka-${PREFIX}/consensus.fasta ${OUT_DIR}/r2-medaka-${PREFIX}.fasta && rm -rf ${OUT_DIR}/r2-medaka-${PREFIX}
fi

##evaluation using metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --no-icarus --fragmented --min-identity 90 --min-contig 5000 \
        --threads ${NPROC} -r ${REF_DIR_ILL} -o ${METAQUAST_OUTDIR}_medaka ${OUT_DIR}/r2-medaka-${PREFIX}.fasta
fi

##polishing using marginPolish
if [ $marginflag == true ]; then
	if [ -f ${OUT_DIR}/map1-${PREFIX}.sam ]; then
	    samtools view -T ${ASSEMBLY} -F 2308 -Sb ${OUT_DIR}/map1-${PREFIX}.sam | samtools sort -@ ${NPROC} - -o ${OUT_DIR}/${ASSEMBLY}.bam
	else
		minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_QUAN} | samtools view -T ${ASSEMBLY} -F 2308 -b - |
			samtools sort -@ ${NPROC} -o ${OUT_DIR}/${ASSEMBLY}.bam
	fi
	samtools index -@ ${NPROC} ${OUT_DIR}/${ASSEMBLY}.bam
	mkdir ${OUT_DIR}/MP-${PREFIX}
	marginPolish ${OUT_DIR}/${ASSEMBLY}.bam ${ASSEMBLY} ${MP_PARAMS} -t ${NPROC} -o ${OUT_DIR}/MP-${PREFIX} -f
	mv ${OUT_DIR}/MP-${PREFIX}/output.fa ${OUT_DIR}/MP-${PREFIX}/MP-${PREFIX}.fasta
fi

##evaluation using metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --no-icarus --fragmented --min-identity 90 --min-contig 5000 \
        --threads ${NPROC} -r ${REF_DIR_ILL} -o ${METAQUAST_OUTDIR}_MP ${OUT_DIR}/MP-${PREFIX}/MP-${PREFIX}.fasta
fi

##polishing using helen (requires previous marginPolish polishing)
if [ $helenflag == true ]; then
	cp ${HELEN_MODEL} $(pwd)
	model=${HELEN_MODEL##*/}
	docker run --rm -it --ipc=host -v $(pwd):/helen/ kishwars/helen:latest \
		helen polish -i /helen/${OUT_DIR}/MP-${PREFIX}/ -m /helen/${model} -b 256 -w 4 -t 8 -o /helen/${OUT_DIR}/HELEN-${PREFIX} -p MP-helen-${PREFIX}
	mv ${OUT_DIR}/HELEN-${PREFIX}/MP-helen-${PREFIX}.fa ${OUT_DIR}/HELEN-${PREFIX}/MP-helen-${PREFIX}.fasta
fi

##evaluation using metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --no-icarus --fragmented --min-identity 90 --min-contig 5000 \
        --threads ${NPROC} -r ${REF_DIR_ILL} -o ${METAQUAST_OUTDIR}_HELEN ${OUT_DIR}/HELEN-${PREFIX}/MP-helen-${PREFIX}.fasta
fi
