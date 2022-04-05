#!/usr/bin/env bash

FRACTION=20 #define percentage of original fastq file used for this run (coverage variation)
QUANT=4bin
PREFIX=${QUANT}-${FRACTION}p-1 # files prefix (change for original or 4bin run, and change last number for replicates)
#RUN=${FRACTION} #run name
NPROC=35 #number of threads to use
READS_NANO=./input/HG00733_2-${QUANT}-${FRACTION}p.fastq # fastq file 
MP_PARAMS=../marginPolish/params/allParams.np.human.r94-g235.json # marginPolish parameters for human data and specific ONT pore, available in https://github.com/UCSC-nanopore-cgl/MarginPolish/tree/master/params
MP_PROGRAM=../marginPolish/build # path to marginPolish software 
HELEN_MODEL=../helen-models/r941_flip231_v001.pkl # HELEN trained model for used ONT pore
REF_DIR=./truth-assemblies/hg00733_truth_assembly.fa #assembly reference, used for performance evaluation


## flags
wtdbgflag=true
marginflag=true
quastmargin=true
helenflag=true
quasthelen=true
# -----------------------------

source activate human-env


##assembly using wtdbg2
if [ $wtdbgflag == true ]; then
	mkdir wtdbg2-${FRACTION}
	cd wtdbg2-${FRACTION}
	wtdbg2 -t ${NPROC} -x ont -L 10000 -g 3.3g -i ../${READS_NANO} -o wtdbg2-assembly-${FRACTION}
	source activate human-env
	wtpoa-cns -t ${NPROC} -i wtdbg2-assembly-${FRACTION}.ctg.lay.gz -f -o raw-wtdbg2-assembly-${FRACTION}.fasta
	ASSEMBLY=raw-wtdbg2-assembly-${FRACTION}.fasta
	cd ../
fi


#polishing using marginPolish
if [ $marginflag == true ]; then
	ASSEMBLY=raw-wtdbg2-assembly-${FRACTION}.fasta
	source activate human-env
	mkdir ./wtdbg2-${FRACTION}/MP-${PREFIX}
	minimap2 -ax map-ont -t ${NPROC} ./wtdbg2-${FRACTION}/${ASSEMBLY} ${READS_NANO} | samtools sort -@ ${NPROC} | samtools view -hb -F 0x104 > ./wtdbg2-${FRACTION}/${ASSEMBLY}.bam #mapping reads to assembly with minimap2
	samtools index -@ ${NPROC} ./wtdbg2-${FRACTION}/${ASSEMBLY}.bam
	${MP_PROGRAM}/marginPolish ./wtdbg2-${FRACTION}/${ASSEMBLY}.bam ./wtdbg2-${FRACTION}/${ASSEMBLY} ${MP_PARAMS} -t ${NPROC} -o ./wtdbg2-${FRACTION}/MP-${PREFIX} -f
	mv ./wtdbg2-${FRACTION}/MP-${PREFIX}/output.fa ./wtdbg2-${FRACTION}/MP-${PREFIX}/MP-${PREFIX}.fasta
fi

#marginPolish performance evaluation using QUAST
if [ $quastmargin == true ]; then
	quast-lg.py --threads ${NPROC} -r ${REF_DIR} --large --no-icarus --space-efficient --min-identity 80 --fragmented ./wtdbg2-${FRACTION}/MP-${PREFIX}/MP-${PREFIX}.fasta -o ./wtdbg2-${FRACTION}/MP-${PREFIX}/QUAST

fi


#polishing using HELEN (it uses marginPolish images as input)
if [ $helenflag == true ]; then
	cp ${HELEN_MODEL} $(pwd)
	model=${HELEN_MODEL##*/}
	docker run -v $(pwd):/helen/ --rm -i --ipc=host --user="$(id -u):$(id -g)" kishwars/helen:latest helen polish --image_dir /helen/wtdbg2-${FRACTION}/MP-${PREFIX} --model_path /helen/${model} --threads ${NPROC} --output_dir /helen/wtdbg2-${FRACTION}/HELEN-${PREFIX} --output_prefix MP-helen-${PREFIX}

fi


#HELEN performance evaluation using QUAST
if [ $quasthelen == true ]; then
	quast-lg.py --threads ${NPROC} -r ${REF_DIR} --large --no-icarus --space-efficient --min-identity 80 --fragmented ./wtdbg2-${FRACTION}/HELEN-${PREFIX}/MP-helen-${PREFIX}.fa -o ./wtdbg2-${FRACTION}/HELEN-${PREFIX}/QUAST
fi
