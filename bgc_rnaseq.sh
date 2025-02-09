#!/bin/bash

# Script Parameters:
# $1: eggnog database
# $2: Input fasta file 
# $3: Output directory  
# $4: py script
# $5: input RNAseq .fastq dir


# Step 1: 
mkdir -p  "$3/eggnog/"
mkdir -p  "$3/log_files/"

# Step 2: Run eggnog mapper
echo -e "Running eggnog mapper...\n"
emapper.py --cpu 48  --data_dir $1  --decorate_gff yes   --tax_scope eukaryota   --itype genome  -m mmseqs  -i $2 --override  --output_dir "$3/eggnog/"   -o eggnog  --dbmem --override  > "$3/log_files/chr_test.log"

# Step 3: Process eggnog annotations
echo -e "Processing eggnog annotations...\n"
sed '/^##/d' "$3/eggnog/eggnog.emapper.annotations" | sed 's/#//' > "$3/eggnog/eggnog.anno"

# Step 4: Process GFF decorated file
echo -e "Processing GFF decorated file...\n"
grep 'ID' "$3/eggnog/eggnog.emapper.decorated.gff" | cut -d ';' -f1 | tr -d ' ' | sed 's/ID=//g' | \
  awk -v OFS='\t' '{print $9,$1,$4,$5,$7}' > "$3/eggnog/eggnog.emapper.decorated.gff.saf"

# Step 5: Match gene IDs
echo -e "Matching gene IDs...\n"
seqkit fx2tab -n "$3/eggnog/eggnog.emapper.genepred.fasta" | awk '{print $1,$9}' | \
  awk -F';' '{print $1}' | sed 's/ID=//g' > "$3/eggnog/geneid_match.tab"

# Step 6: Generate GenBank file
echo -e "Creating GenBank file...\n"
python  $4  "$3/eggnog/eggnog.emapper.decorated.gff"  $2  "$3/generated.gbk"


# Step 7: Index genome with Hisat2
echo -e "Building mapping index with Hisat2...\n"
mkdir -p "$3/hisat2_index/"
hisat2-build -p 8 "$2" "$3/hisat2_index/geno_index"

# Step 8: Grep the list of fastq files
echo -e "Grep the list of fastq files...\n"
find $5 -type f \( -iname "*.fastq" -o -iname "*.fq" -o -iname "*.fq.gz" -o -iname "*.clean.fq" -o -iname "*.clean.fq.gz" \) | sort | paste - - > "$3/fq_list.txt"

# Step 9: Run Fastp and Hisat2
echo -e "Running Fastp and Hisat2...\n"
cat "$3/fq_list.txt" | while read file1 file2
do
  filename1=$(basename "$file1")
  filename2=$(basename "$file2")
  dirfile1=$(dirname $file1)
  dirname=$(basename $dirfile1)

  echo DIR: "$dirname"
  mkdir -p "$3/$dirname"

  fastp -i "$file1"  -I "$file2" -o "${3}/${dirname}/clean_${filename1}" -O "${3}/${dirname}/clean_${filename2}" -j "${3}/${dirname}/${dirname}.json" -h "${3}/${dirname}/${dirname}.html" -w 8
  hisat2 -p 18 -x "$3/hisat2_index/geno_index" -1 "${3}/${dirname}/clean_${filename1}" -2 "${3}/${dirname}/clean_${filename2}" | samtools sort -@ 18 -o "${3}/${dirname}_sorted.bam"
  samtools index "${3}/${dirname}_sorted.bam"

done

# Step 10: Run featureCounts
echo -e "Running featureCounts...\n"
featureCounts -T 18 -p -B -t CDS -F SAF -a "$3/eggnog/eggnog.emapper.decorated.gff.saf" -o "$3/all_counts.tsv" "$3/"*".bam"


echo -e "\nEND\n${1}\n${2}"
