# QS-Quantizer

We investigate the effect of quality score information loss on downstream analysis from nanopore sequencing FASTQ files. We polished *denovo* assemblies for a mock microbial community and a human genome, and we called variants on a human genome. We repeated these experiments using various pipelines, under various coverage level scenarios,  and various quality score quantizers.  In all cases we found  that  the  quantization  of  quality  scores  cause  little  difference  on (or even improves) the results obtained with the original (non-quantized) data. This suggests that the precision that is currently used for nanopore quality scores is unnecessarily high, and motivates the use of lossy compression algorithms for this kind of data. Moreover, we show that even a non-specialized compressor, like gzip, yields large storage space savings after quantization of quality scores.


Softwares used in this work are listed in the table below. 


Software versions:

| Software        | Version   |
| -------------   |:---------:|
| bedtools        | 2.27.1    |
| deepVariant     | r.0.4     |
| flye            | 2.8.2     |
| HELEN           | 0.0.23    |
| Longshot        | 0.4.1     |
| marginPolish    | 1.3.0     |
| medaka          | 0.5.0     |
| metaQUAST       | 5.0.2     |
| minimap2        | 2.15      |
| minimap2 (Zymo) | 2.14      |
| QUAST           | 5.0.2     |
| racon           | 1.3.2     |
| samtools        | 1.9       |
| wtdbg2          | 2.3       |


## Quantization of quality scores
We tested various quantizers, each mapping quality scores to values on a (small) set referred to as the quantization alphabet. We denote by Qi a quantizer for a quantization alphabet of size i, i > 1, where the quantization Qi (x) of a quality score x depends solely on x. The specific definitions of Q2, Q4, and Q8 are presented in tables 1, 2, and 3, respectively. All these quantizers collapse a large set of high quality scores into a single value, and define a finer partition for lower scores, which occur more frequently. We also tested constant quantizers, which map every quality score to a fixed prescribed value.


### Table 1 (Q2)
| Quality scores | Quantized quality score | 
| -------------  |:-----------------------:|
| 0....7         | 5                       |
| 8...93         | 15                      |


### Table 2 (Q4)
| Quality scores | Quantized quality score | 
| -------------  |:-----------------------:|
| 0....7         | 5                       |
| 8...13         | 12                      |
| 14....19       | 18                      |
| 20...93        | 24                      |

### Table 3 (Q8)
| Quality scores | Quantized quality score | 
| -------------  |:-----------------------:|
| 0....6         | 5                       |
| 7...11         | 10                      |
| 12....16       | 15                      |
| 17...21        | 20                      |
| 22....26       | 25                      |
| 27...31        | 30                      |
| 32....36       | 35                      |
| 37...93        | 40                      |

### Usage
``quantizer.py [optional arguments] input_file Qprimary Qrun output_file``

| Positional arguments                                                                                                                            | 
| ------------------   | ------------------------------------------------------------------------------------------------------------------------ |
| input_file           | Input fastq file, either in plain text formtat or gzip compressed                                                        |
| Qprimary             | Quantizer used in non-repetitive regions. One of Q2, Q4, Q8, or F10                                                      |
| Qrun                 | Quantizer used near repetitive regions. One of Q2, Q4, Q8, or F10. Specify the same as Qprimary to use a fixed quantizer |
| output_file          | Output fastq file in plain text format                                                                                   |
| optional arguments                                                                                                                              |
| ------------------   | ------------------------------------------------------------------------------------------------------------------------ |
| -h, --help           | show this help message and exit                                                                                          |


## Assembly and Polishing of mock community

We evaluated the impact of quality score quantization on the genome assembly polishing for a Zymo-BIOMICS Microbial Community Standard. We used data generated in this [article](https://pubmed.ncbi.nlm.nih.gov/31089679/). We assembled the genomes with Flye, and we polished this raw assembly using various polishing pipelines (Racon, Medaka, MarginPolish and HELEN). We tested each polishing pipeline on the original (non-quantized) data, and in various quantized versions. "pipeline-mock-community.sh" is the script used for these experiments.

## Assembly and Polishing of human genome

We evaluated the impact of quality score quantization on human genome assembly polishing for sample HG00733 with the polishing pipelines MP and Helen. We used data generated in this [article](https://pubmed.ncbi.nlm.nih.gov/32686750/). We used wtdbg2 for human genome assembly and Margin Polish/HELEN pipeline for polishing. These polishing pipelines were executed both for the orginal FASTQ files and for 4 bin quantized data. We carried on this comparison for several coverage scenarios, which we obtained by randomly selecting a fraction of the dataset reads. "pipeline-human-assembly.sh" is the script used for these experiments.

## Variant Calling

We compared the nanopore variant calling performance of PEPPER-Margin-DeepVariant on human sample HG003, against variant calling on quantized versions of the same data. We performed this comparison at various coverages, ranging from 20X to 90X, and for various quantizers. We used data generated in this [article](https://pubmed.ncbi.nlm.nih.gov/34725481/). "bam-script.sh" and "pipeline-variant-calling.sh" are the scripts used for these experiments.
