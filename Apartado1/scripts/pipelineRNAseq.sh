# EJECUTAR DESDE EL DIRECTORIO Apartado1 
# CONTROL DE CALIDAD
mkdir -p QCrawfiles
fastqc input/*.fastq -o QCrawfiles
# TRIMMING DE SECUENCIAS
# Prueba 1. Cutadapt.
# Ejecuto cutadapt.sh que va a generar los archivos trimmed que se usan después.
bash scripts/cutadapt.sh
# Prueba 2. Cutadapt.
# Muestra SRR479052. Elimino adaptadores, sobrerrepresentadas y filtro por calidad.
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG \
-g CTTTTACTTCCTCTAGATAGTCAAGTTCGACCGTCTTCTCAGCGCTCCGC -G CTAACACGTGCGCGAGTCGGGGGCTCGCACGAAAGCCGCCGTGGCGCAAT \
 --quality-cutoff=20 \
 -o cutadapt/SRR479052.chr21_1_trimmed2.fastq -p cutadapt/SRR479052.chr21_2_trimmed2.fastq \
 input/SRR479052.chr21_1.fastq input/SRR479052.chr21_2.fastq
# Muestra SRR479054 Elimino adaptadores y filtro por calidad.
cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG --quality-cutoff=20 \
 -o cutadapt/SRR479054.chr21_1.trimmed2.fastq -p cutadapt/SRR479054.chr21_2.trimmed2.fastq \
 input/SRR479054.chr21_1.fastq input/SRR479054.chr21_2.fastq

# Prueba 3. Trimmomatic
# Muestra SRR479052. Elimino sobrerrepresentadas y recorto las últimas 21pb.
mkdir -p trimmomatic
trimmomatic PE -threads 8 input/SRR479052.chr21_1.fastq input/SRR479052.chr21_2.fastq \
 trimmomatic/SRR479052.chr21_1.trimmed.fastq trimmomatic/SRR479052.chr21_1.trimmedunpaired.fastq \
 trimmomatic/SRR479052.chr21_2.trimmed.fastq trimmomatic/SRR479052.chr21_2.trimmedunpaired.fastq \
 CROP:80 ILLUMINACLIP:sobrerrepresentadas.fa:1:54:10 MINLEN:36
# Muestra SRR479054. Recorto las últimas 21pb.
trimmomatic PE -threads 8 input/SRR479054.chr21_1.fastq input/SRR479054.chr21_2.fastq \
 trimmomatic/SRR479054.chr21_1.trimmed.fastq trimmomatic/SRR479054.chr21_1.trimmedunpaired.fastq \
 trimmomatic/SRR479054.chr21_2.trimmed.fastq trimmomatic/SRR479054.chr21_2.trimmedunpaired.fastq \
 CROP:80   

# ALINEAMIENTO DE SECUENCIAS CON HISAT2

#Índice del genoma de referencia
mkdir -p input/REF
hisat2-build --seed 123 -p 2 input/Homo_sapiens.GRCh38.dna.chromosome.21.fa input/REF/chr21_GRCh38

# Alineamiento asumiendo orientación específica del tipo stranded=reverse
mkdir -p hisat2

hisat2 --new-summary --summary-file hisat2/SRR479052.hisat2.summary \
       -1 cutadapt/SRR479052.chr21_1_trimmed.fastq \
       -2 cutadapt/SRR479052.chr21_2_trimmed.fastq \
       --rna-strandness FR --seed 123 --phred33 -p 2 -k 1 \
       -x input/REF/chr21_GRCh38 -S hisat2/SRR479052.sam

hisat2 --new-summary --summary-file hisat2/SRR479054.hisat2.summary \
       -1 cutadapt/SRR479054.chr21_1_trimmed.fastq \
       -2 cutadapt/SRR479054.chr21_2_trimmed.fastq \
       --rna-strandness FR --seed 123 --phred33 -p 2 -k 1 \
       -x input/REF/chr21_GRCh38 -S hisat2/SRR479054.sam

# REFINAMIENTO DEL ALINEAMIENTO
# 1.Transformo el archivo del alineamiento al formato binario .bam.
samtools view -bS hisat2/SRR479052.sam > hisat2/SRR479052.bam
samtools view -bS hisat2/SRR479054.sam > hisat2/SRR479054.bam
# 2.Ordeno el archivo .bam del alineamiento
samtools sort hisat2/SRR479052.bam -o hisat2/SRR479052.sorted.bam 
samtools sort hisat2/SRR479054.bam -o hisat2/SRR479054.sorted.bam
# 3.Elimino duplicados con picard.jar Markduplicates. Contar con picard.jar dentro del directorio de ejecución.
java -jar picard.jar MarkDuplicates I=hisat2/SRR479052.sorted.bam O=hisat2/SRR479052.refined.bam \
 M=hisat2/redupSRR479052.txt ASSUME_SORTED=true

java -jar picard.jar MarkDuplicates I=hisat2/SRR479054.sorted.bam O=hisat2/SRR479054.refined.bam \
 M=hisat2/redupSRR479054.txt ASSUME_SORTED=true

# 4.Creo un índice para cada alineamiento. 
samtools index hisat2/SRR479052.refined.bam
samtools index hisat2/SRR479054.refined.bam

# Alineamiento non-stranded tras sugerirlo resultados de visualización de bam en igv.
mkdir -p hisat2/non-stranded
hisat2 --new-summary --summary-file hisat2/non-stranded/SRR479052.hisat2.summary \
       -1 cutadapt/SRR479052.chr21_1_trimmed.fastq \
       -2 cutadapt/SRR479052.chr21_2_trimmed.fastq \
        --seed 123 --phred33 -p 2 -k 1 \
       -x input/REF/chr21_GRCh38 -S hisat2/non-stranded/SRR479052.sam

hisat2 --new-summary --summary-file hisat2/non-stranded/SRR479054.hisat2.summary \
       -1 cutadapt/SRR479054.chr21_1_trimmed.fastq \
       -2 cutadapt/SRR479054.chr21_2_trimmed.fastq \
       --seed 123 --phred33 -p 2 -k 1 \
       -x input/REF/chr21_GRCh38 -S hisat2/non-stranded/SRR479054.sam
# 1.Transformo el archivo del alineamiento al formato binario .bam. 
samtools view -bS hisat2/non-stranded/SRR479052.sam > hisat2/non-stranded/SRR479052.bam
samtools view -bS hisat2/non-stranded/SRR479054.sam > hisat2/non-stranded/SRR479054.bam
# 2.Ordeno el archivo .bam del alineamiento
samtools sort hisat2/non-stranded/SRR479052.bam -o hisat2/non-stranded/SRR479052.sorted.bam 
samtools sort hisat2/non-stranded/SRR479054.bam -o hisat2/non-stranded/SRR479054.sorted.bam
# 3.Elimino duplicados con picard.jar Markduplicates.
java -jar picard.jar MarkDuplicates I=hisat2/non-stranded/SRR479052.sorted.bam O=hisat2/non-stranded/SRR479052.refined.bam \
 M=hisat2/non-stranded/redupSRR479052.txt ASSUME_SORTED=true

java -jar picard.jar MarkDuplicates I=hisat2/non-stranded/SRR479054.sorted.bam O=hisat2/non-stranded/SRR479054.refined.bam \
 M=hisat2/non-stranded/redupSRR479054.txt ASSUME_SORTED=true

# 4.Creo un índice para cada alineamiento. 
samtools index hisat2/non-stranded/SRR479052.refined.bam
samtools index hisat2/non-stranded/SRR479054.refined.bam

# CUANTIFICACIÓN CON HTSEQ.

bash scripts/counts.sh

# NORMALIZACIÓN CON bamCoverage

bamCoverage -b hisat2/non-stranded/SRR479052.refined.bam -o htseq/non-strandedboth/SRR479052.bw --normalizeUsing BPM
bamCoverage -b hisat2/non-stranded/SRR479054.refined.bam -o htseq/non-strandedboth/SRR479054.bw --normalizeUsing BPM

