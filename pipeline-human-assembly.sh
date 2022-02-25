#!/usr/bin/env bash

PREFIX=r94-4bin-60p # files prefix 
RUN=r94-60p #run name
ASSEMBLY=raw-wtdbg2-assembly-r94.fasta #wtdbg2 assembly file in fasta format
NPROC=35 #number of threads to use
READS_NANO=/path/to/fastq/_file/HG00733_2-4bin-60p.fastq # fastq file 
MP_PARAMS=/path/to/marginPolish/params/allParams.np.human.r94-g235.json # marginPolish parameters for used ONT pore, available in https://github.com/UCSC-nanopore-cgl/MarginPolish/tree/master/params
MP_PROGRAM=/path/to/marginPolish/software/marginPolish/build # path to marginPolish software 
HELEN_MODEL=/path/to/helen-models/helen_trained_models_v0.0.1_r941_flip231_v001.pkl # HELEN trained model for used ONT pore
REF_DIR=/path/to/truth-assemblies/truth_assemblies_HG00733_hg00733_truth_assembly.fa #assembly reference, used for performance evaluation


## flags
wtdbgflag=true
marginflag=true
helenflag=true
quastmargin=true
quasthelen=true
# -----------------------------

source activate ambiente


##assembly using wtdbg2
if [ $wtdbgflag == true ]; then
	mkdir wtdbg2-${PREFIX}
	cd wtdbg2-${PREFIX}
	wtdbg2 -t ${NPROC} -x ont -L 10000 -g 3.3g -i ${READS_NANO} -o wtdbg2-assembly-${PREFIX}
	source activate cuantinano
	wtpoa-cns -t ${NPROC} -i wtdbg2-assembly-${PREFIX}.ctg.lay.gz -f -o raw-wtdbg2-assembly-${PREFIX}.fasta
	ASSEMBLY=raw-wtdbg2-assembly-${PREFIX}.fasta
fi


#polishing using marginPolish
if [ $marginflag == true ]; then
	ASSEMBLY=raw-wtdbg2-assembly-${RUN}.fasta
	source activate cuantinano
	mkdir /path/to/output/directory/MP-${PREFIX}
	minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_NANO} | samtools sort -@ ${NPROC} | samtools view -hb -F 0x104 > ${ASSEMBLY}.bam #mapping reads to assembly with minimap2
	samtools index -@ ${NPROC} ${ASSEMBLY}.bam
	${MP_PROGRAM}/marginPolish ${ASSEMBLY}.bam ${ASSEMBLY} ${MP_PARAMS} -t ${NPROC} -o /path/to/output/directory/MP-${PREFIX} -f
	mv /path/to/output/directory/MP-${PREFIX}/output.fa /path/to/output/directory/MP-${PREFIX}/MP-${PREFIX2}.fasta
fi


#polishing using HELEN (it uses marginPolish images as input)
if [ $helenflag == true ]; then
	cp ${HELEN_MODEL} $(pwd)
	model=${HELEN_MODEL##*/}
	docker run -v $(pwd):/helen/ --rm -i --ipc=host --user="$(id -u):$(id -g)" kishwars/helen:latest helen polish --image_dir /helen/MP-${PREFIX} --model_path /helen/${model} --threads ${NPROC} --output_dir /helen/HELEN-${PREFIX2} --output_prefix MP-helen-${PREFIX}
	mv MP-helen-${PREFIX}.fa MP-helen-${PREFIX}.fasta
fi


#marginPolish performance evaluation using QUAST
if [ $quastmargin == true ]; then
	quast-lg.py --threads ${NPROC} -r ${REF_DIR} --large --no-icarus --space-efficient --min-identity 80 --fragmented /path/to/output/directory/MP-${PREFIX}/MP-${PREFIX}.fasta -o /path/to/output/directory/MP-${PREFIX}/QUAST

fi

#HELEN performance evaluation using QUAST
if [ $quasthelen == true ]; then
	quast-lg.py --threads ${NPROC} -r ${REF_DIR} --large --no-icarus --space-efficient --min-identity 80 --fragmented HELEN-${PREFIX}/MP-helen-${PREFIX}.fasta -o HELEN-${PREFIX}/QUAST
fi




