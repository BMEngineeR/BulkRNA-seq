# BulkRNA-seq
***BulkRNAseq pipeline (Based on Linux system)***
## 1. __workflow__

  - fastq ---> bam ---> gene (in row) \* sample (in column) matrix ---> R script (DE, pathway)
## 2. notice:
  - if you only have two samples per treatment(or control), please use edgeR to avoid outlier gene effect. 

## 2.  Pre-requirement:
### software:
  - aligner: [hisat2](https://ccb.jhu.edu/software/hisat2/index.shtml) `wget http://ccb.jhu.edu/software/hisat2/dl/hisat2-2.1.0-Linux_x86_64.zip`
  - gene-counter: please follow the install instruction [HTseq](https://htseq.readthedocs.io/en/release_0.11.1/install.html#installation-on-linux)
### reference:
  - download from hisat2 offcial website: 
    - Mouse: mm10 `wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz` 
    - human: GRCh38 `wget ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/grch38.tar.gz`
  - download from ensemble, remove alternative chromosome, then use `hisat2-built` to index the reference:
    - Mouse: GRCm38(GRCm38 = mm10) `wget ftp://ftp.ensembl.org/pub/release-97/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna_rm.primary_assembly.fa.gz`
    - human: GRCh38 `wget ftp://ftp.ensembl.org/pub/release-97/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz`
### annotation:
> (*please make sure the consistent version between reference and annotation. *)
  - human: gtf (GRCh38) `wget ftp://ftp.ensembl.org/pub/release-97/gtf/homo_sapiens/Homo_sapiens.GRCh38.97.chr.gtf.gz` 
  - mouse: gtf (GRCm38) `wget ftp://ftp.ensembl.org/pub/release-97/gtf/mus_musculus/Mus_musculus.GRCm38.97.chr.gtf.gz`
