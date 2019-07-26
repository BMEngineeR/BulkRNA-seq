HISAT=/home/hisat2-2.1.0/hisat2
SAMTOOLS=/usr/local/bin/samtools
MOUSE_SEQ=/home/reference/genome/mm10/genome_tran/genome_tran
MOUSE_ANNO=/home/reference/annotation/mm10/Mus_musculus.GRCm38.83.gtf
HUMAN_SEQ=/home/reference/genome/GRCh38_hisat/grch38_tran/genome_tran
HUMAN_ANNO=/home/reference/annotation/GRCh38.91/Homo_sapiens.GRCh38.91.gtf

READ_SUFFIX=_001.fastq.gz
TRIM3=0

rand=$RANDOM

#generate list of sample names
ls *R1_001.fastq*|cut -d_ -f1 > mylist




for i in `cat mylist`
do
echo $i
#run hisat
$HISAT -p 2 -3 $TRIM3 --summary-file ${i}.hisat_summary --new-summary -x $MOUSE_SEQ -1 ${i}*R1${READ_SUFFIX} -2 ${i}*R2${READ_SUFFIX}  | $SAMTOOLS view -Shb - | $SAMTOOLS sort -T ${i}_${rand} - > ${i}.sort.bam
$SAMTOOLS index ${i}.sort.bam
# new version use mapQ to filter >= 5 seems good
$SAMTOOLS view -q 5 -h ${i}.sort.bam | $SAMTOOLS view -Shb - | $SAMTOOLS sort -T ${i}_${rand} -n - > ${i}.uniq.nsort.bam
#run HTSeq, you may change parameter --stranded= reverse, no, forward, based on library preparation protocol.
$SAMTOOLS view ${i}.uniq.nsort.bam | python -m HTSeq.scripts.count --mode=union --stranded=reverse - $MOUSE_ANNO > ${i}.count.union.reverse
#$SAMTOOLS view ${i}.uniq.nsort.bam | python -m HTSeq.scripts.count --mode=union --stranded=no - $HUMAN_ANNO > ${i}.count.union.no

done


