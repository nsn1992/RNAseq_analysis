#!/bin/bash

# Directorio de entrada y salida
mkdir -p trimmed
input_dir="input"
output_dir="trimmed"

# Loop sobre todos los archivos en el directorio de entrada
for forward_read in "$input_dir"/*_1.fastq; do
    # Obtener el nombre base del archivo forward
    base_name=$(basename "$forward_read" _1.fastq)
    # Definir los nombres de archivo de entrada y salida
    forward_read_input="$input_dir/$base_name"_1.fastq
    reverse_read_input="$input_dir/$base_name"_2.fastq
    output_base_name="$output_dir/$base_name"trimmed
    # Ejecutar Trimmomatic
    trimmomatic PE -threads 8 "$forward_read_input" "$reverse_read_input" \
        "$output_base_name"_1.fastq "$output_base_name"_unpaired_1.fastq \
        "$output_base_name"_2.fastq "$output_base_name"_unpaired_2.fastq \
        CROP:80 ILLUMINACLIP:sobrerrepresentadas.fa:1:54:10.true MINLEN:36
done

