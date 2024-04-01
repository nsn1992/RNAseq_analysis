
# Non stranded en htseq y alineados como tal
mkdir -p htseq/non-strandedboth

htseq-count --format=bam --stranded=no --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/non-stranded/SRR479052.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/non-strandedboth/SRR479052.htseq

htseq-count --format=bam --stranded=no --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/non-stranded/SRR479054.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/non-strandedboth/SRR479054.htseq

# Non stranded en htseq y alineados como stranded=R
mkdir -p htseq/non-strandedvsstrandedR

htseq-count --format=bam --stranded=no --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/SRR479052.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/non-strandedvsstrandedR/SRR479052.htseq

htseq-count --format=bam --stranded=no --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/SRR479054.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/non-strandedvsstrandedR/SRR479054.htseq

# Stranded=reverse en htseq sobre archivos alineados como non-stranded
mkdir -p htseq/strandedRvsnon-stranded
htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/non-stranded/SRR479052.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/strandedRvsnon-stranded/SRR479052.htseq

htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/non-stranded/SRR479054.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/strandedRvsnon-stranded/SRR479054.htseq


# Stranded=reverse sobre archivos alineados como tal
mkdir -p htseq/strandedRboth
htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/SRR479052.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/strandedRboth/SRR479052.htseq

htseq-count --format=bam --stranded=reverse --mode=intersection-nonempty --minaqual=10 --type=exon --idattr=gene_id \
 --additional-attr=gene_name hisat2/SRR479054.refined.bam input/Homo_sapiens.GRCh38.109.chr21.gtf > htseq/strandedRboth/SRR479054.htseq