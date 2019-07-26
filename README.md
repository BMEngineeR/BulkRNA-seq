# BulkRNA-seq
* BulkRNAseq pipeline
1. __workflow__

  - fastq -> bam -> gene(in row)\*sample(in column) matrix -> R script (DE, pathway)

2.  Prequirement: 
  - aligner: [hisat2]()
  - reference:
    - download from hisat2 offcial website
    - download from ensemble then use `hisat2-built` to index the reference
