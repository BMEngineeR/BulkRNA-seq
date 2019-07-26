# BulkRNA-seq
***BulkRNAseq pipeline (Based on Linux system)***
## 1. __workflow__

  - fastq ---> bam ---> gene (in row) \* sample (in column) matrix ---> R script (DE, pathway)

## 2.  Prequirement:
### software:
  - aligner: [hisat2](http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip)
  - counter: [HTseq]()
  - reference:
    - download from hisat2 offcial website: 
      - Mouse: [mm10](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz) 
      - human: [GRCh38](ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz)
    - download from ensemble, remove alternative chromosome, then use `hisat2-built` to index the reference:
      - Mouse: [GRCm38](ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa.gz)(GRCm38 = mm10)
      - human: [GRCh38](ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz)
  - annotation(*please make sure the consistent version between reference and annotation. *)
    - human: [gtf (GRCh38)](ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.chr.gtf.gz) 
    - mouse: [gtf (GRCm38)](ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.chr.gtf.gz)
