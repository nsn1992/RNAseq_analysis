mkdir -p cutadapt
for file in  input/*.fastq; do
	new_file_name=$(basename "$file" .fastq)
	cutadapt -u -21 -o "cutadapt/${new_file_name}_trimmed.fastq" "$file" 
done
