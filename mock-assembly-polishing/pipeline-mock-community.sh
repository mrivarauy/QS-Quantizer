#!/usr/bin/env bash

PREFIX=r10-3 #files prefix 
ASSEMBLY=raw-${PREFIX}.fasta #assembly file name
NPROC=32 #number of threads to use
READS_NANO=/path/to/zymo_mock_community/fastq.gz #path to fastq.gz files
MEDAKA_MODEL=r103_min_high_g345 #used medaka model
MP_PARAMS=/path/to/marginPolish/params/allParams.np.microbial.r103g324.json #marginPolish parameters for microbial data and specific ONT pore, available in https://github.com/UCSC-nanopore-cgl/MarginPolish/tree/master/params
HELEN_MODEL=/path/to/helen-models/HELEN_r103_guppy_microbial.pkl # HELEN trained model for used ONT pore
REF_DIR_ILL=/path/to/ref_dir/ #assembly reference, used for performance evaluation
METAQUAST_OUTDIR=outdir #output directory name for metaQuast evaluation
## flags, change to false to not run
flyeflag=true
raconflag=true
medakaflag=true
marginflag=true
helenflag=true
metaquastflag=true

# -----------------------------

##assembly using flye
if [ $flyeflag == true ]; then
	flye --nano-raw ${READS_NANO} --meta --iterations 0 --threads ${NPROC} --out-dir assembly-${PREFIX}
	infoseq --only --length --name assembly-${PREFIX}/assembly.fasta | awk '{print $2, $1}' | 
		sort -nr | head -n8 | cut -d' ' -f2 > best_contigs
	fastaUtils.pl -u assembly-${PREFIX}/assembly.fasta | 
		grep --no-group-separator -A 1 -w -f best_contigs > raw-${PREFIX}.fasta
	ASSEMBLY=raw-${PREFIX}.fasta
fi

##polishing using racon
if [ $raconflag == true ]; then
    minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_NANO} > map1-${PREFIX}.sam
    racon -t ${NPROC} ${READS_NANO} map1-${PREFIX}.sam ${ASSEMBLY} > racon_r1-${PREFIX}.fasta
    minimap2 -ax map-ont -t ${NPROC} racon_r1-${PREFIX}.fasta ${READS_NANO} > map2-${PREFIX}.sam
    racon -t ${NPROC} ${READS_NANO} map2-${PREFIX}.sam racon_r1-${PREFIX}.fasta > racon_r2-${PREFIX}.fasta
fi

## polishing using medaka (requires previous racon polishing)
if [ $medakaflag == true ]; then
    medaka_consensus -i ${READS_NANO} -d racon_r1-${PREFIX}.fasta -o r1-medaka-${PREFIX} -m ${MEDAKA_MODEL} -t ${NPROC}
    mv r1-medaka-${PREFIX}/consensus.fasta r1-medaka-${PREFIX}.fasta && rm -rf r1-medaka-${PREFIX}
    medaka_consensus -i ${READS_NANO} -d racon_r2-${PREFIX}.fasta -o r2-medaka-${PREFIX} -m ${MEDAKA_MODEL} -t ${NPROC}
    mv r2-medaka-${PREFIX}/consensus.fasta r2-medaka-${PREFIX}.fasta && rm -rf r2-medaka-${PREFIX}
fi

##polishing using marginPolish
if [ $marginflag == true ]; then
	if [ -f map1-${PREFIX}.sam ]; then
	    samtools view -T ${ASSEMBLY} -F 2308 -Sb map1-${PREFIX}.sam | samtools sort -@ ${NPROC} - -o ${ASSEMBLY}.bam
	else
		minimap2 -ax map-ont -t ${NPROC} ${ASSEMBLY} ${READS_NANO} | samtools view -T ${ASSEMBLY} -F 2308 -b - |
			samtools sort -@ ${NPROC} -o ${ASSEMBLY}.bam
	fi
	samtools index -@ ${NPROC} ${ASSEMBLY}.bam
	mkdir MP-${PREFIX}
	marginPolish ${ASSEMBLY}.bam ${ASSEMBLY} ${MP_PARAMS} -t ${NPROC} -o MP-${PREFIX} -f
	mv MP-${PREFIX}/output.fa MP-${PREFIX}.fasta
fi

##polishing using helen (requires previous marginPolish polishing)
if [ $helenflag == true ]; then
	cp ${HELEN_MODEL} $(pwd)
	model=${HELEN_MODEL##*/}
	docker run --rm -it --ipc=host -v $(pwd):/helen/ kishwars/helen:latest \
		helen polish -i /helen/MP-${PREFIX}/ -m /helen/${model} -b 256 -w 4 -t 8 -o /helen/ -p MP-helen-${PREFIX}
	mv MP-helen-${PREFIX}.fa MP-helen-${PREFIX}.fasta
fi

##evaluacion using metaquast
if [ $metaquastflag == true ]; then
    metaquast.py --no-icarus --fragmented --min-identity 90 --min-contig 5000 \
        --threads ${NPROC} -r ${REF_DIR_ILL} -o metaquast-${METAQUAST_OUTDIR} *.fasta
fi
