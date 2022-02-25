# QS-Quantizer

We  investigate  the  effect  of  quality  score  information  loss  on  downstream  analysis  from  nanoporesequencing FASTQ files.  We polished denovo assemblies for a mock microbial community and a humangenome,  and  we  called  variants  on  a  human  genome.   We  repeated  these  experiments  using  variouspipelines,  under various coverage level scenarios,  and various quality score quantizers.  In all cases wefound  that  the  quantization  of  quality  scores  cause  little  difference  on  (or  even  improves)  the  resultsobtained with the original (non-quantized) data.  This suggests that the precision that is currently usedfor nanopore quality scores is unnecessarily high, and motivates the use of lossy compression algorithmsfor this kind of data.  Moreover, we show that even a non-specialized compressor, like gzip, yields largestorage space savings after quantization of quality scores.


Softwares used in this work are listed in the table below. 


Software versions:

| Software      | Version   |
| ------------- |:---------:|
| Longshot      | 0.4.1     |
| deepVariant   | r.0.4     |
| minimap2      | 2.15-r905 |
| racon         | 1.3.2     |
| marginPolish  | 1.3.0     |
| HELEN         | 0.0.23    |
| QUAST         | 5.0.2     |
| metaQUAST     | 5.0.2     |
| bedtools      | 2.27.1    |
