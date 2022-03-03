# QS-Quantizer

We investigate the effect of quality score information loss on downstream analysis from nanopore sequencing FASTQ files. We polished *denovo* assemblies for a mock microbial community and a human genome, and we called variants on a human genome. We repeated these experiments using various pipelines, under various coverage level scenarios,  and various quality score quantizers.  In all cases we found  that  the  quantization  of  quality  scores  cause  little  difference  on (or even improves) the results obtained with the original (non-quantized) data. This suggests that the precision that is currently used for nanopore quality scores is unnecessarily high, and motivates the use of lossy compression algorithms for this kind of data. Moreover, we show that even a non-specialized compressor, like gzip, yields large storage space savings after quantization of quality scores.


Softwares used in this work are listed in the table below. 


Software versions:

| Software        | Version   |  | Software        | Version   |
| -------------   |:---------:|  | -------------   |:---------:|
| Longshot        | 0.4.1     |  | bedtools        | 2.27.1    |
| deepVariant     | r.0.4     |  | deepVariant     | r.0.4     |
| minimap2 (H.A.) | 2.15      |  | flye            | 2.8.2     |
| minimap2 (Zymo) | 2.14      |  | HELEN           | 0.0.23    |
| racon           | 1.3.2     |  | Longshot        | 0.4.1     |
| marginPolish    | 1.3.0     |  | marginPolish    | 1.3.0     |
| HELEN           | 0.0.23    |  | medaka          | 0.5.0     |
| QUAST           | 5.0.2     |  | metaQUAST       | 5.0.2     |
| metaQUAST       | 5.0.2     |  | minimap2 (H.A.) | 2.15      |
| bedtools        | 2.27.1    |  | minimap2 (Zymo) | 2.14      |
| samtools        | 1.9       |  | QUAST           | 5.0.2     |
| wtdbg2          | 2.3       |  | racon           | 1.3.2     |
| medaka          | 0.5.0     |  | samtools        | 1.9       |
| flye            | 2.8.2     |  | wtdbg2          | 2.3       |


## Assembly and Polishing of mock community

We evaluated the impact of quality score quantization on the genome assembly polishing for a Zymo-BIOMICS Microbial Community Standard. We used data generated in this [article](https://pubmed.ncbi.nlm.nih.gov/31089679/). We assembled the genomes with Flye, and we polished this raw assembly using various polishing pipelines (Racon, Medaka, MarginPolish and HELEN). We tested each polishing pipeline on the original (non-quantized) data, and in various quantized versions. "pipeline-mock-community.sh" is the script used for these experiments.

## Assembly and Polishing of human genome

We evaluated the impact of quality score quantization on human genome assembly polishing for sample HG00733 with the polishing pipelines MP and Helen. We used data generated in this [article](https://pubmed.ncbi.nlm.nih.gov/32686750/). We used wtdbg2 for human genome assembly and Margin Polish/HELEN pipeline for polishing. These polishing pipelines were executed both for the orginal FASTQ files and for 4 bin quantized data. We carried on this comparison for several coverage scenarios, which we obtained by randomly selecting a fraction of the dataset reads. "pipeline-human-assembly.sh" is the script used for these experiments.
