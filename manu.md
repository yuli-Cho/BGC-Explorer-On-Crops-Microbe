**1.Running emapper.py**  
performs the gene annotation and mapping task, with results stored in the specified output directory. Log output is redirected to a log file for tracking the process.  
```
nohup emapper.py --cpu 48  --data_dir /your/path/of/eggnog-mapper-software/    --decorate_gff yes   --tax_scope eukaryota   --itype genome  -m mmseqs  -i /your/path/of/fasta/all.chrs.fasta --override  --output_dir /your/output/dir/output_name   -o eggnog > ./your/output/log_files/output_name.log  2>&1  &
```


**2.Clean Annotation File**   
removes lines starting with ## from the annotation file, which are typically comments. It then removes the # character from the remaining lines, creating a cleaned annotation file (eggnog.anno) for further use.  
```
sed '/^##/d' "/your/output/dir/output_name/eggnog.emapper.annotations" |sed 's/#//' > "/your/output/dir/output_name/eggnog.anno"
```


**3.Extract Specific Fields from GFF File**  
extracts lines containing ID from the GFF file, processes the relevant fields to remove unnecessary characters, and formats the output into a tab-separated file (eggnog.emapper.decorated.gff.saf) with specific columns.  
```
grep 'ID' "/your/output/dir/output_name/eggnog.emapper.decorated.gff" | cut -d ';' -f1 | tr -d ' ' | sed 's/ID=//g' | awk -v OFS='\t' '{print $9,$1,$4,$5,$7}' > "/your/output/dir/output_name/eggnog.emapper.decorated.gff.saf"
```


 **4.Extract Gene IDs from Fasta File**  
extracts gene identifiers from the eggnog.emapper.genepred.fasta file by converting it to a tab-delimited format, stripping the ID= prefix, and saving the processed gene IDs into a separate file (geneid_match.tab).  
```
seqkit fx2tab -n  "/your/output/dir/output_name/eggnog.emapper.genepred.fasta" | awk '{print $1,$9}' | awk -F';' '{print $1}' | sed 's/ID=//g' > "/your/output/dir/output_name/geneid_match.tab"
```


**5.run RNAseq analysis**  
Please download and save the egg_rice_rnaseq_part.sh script for RNAseq analysis.  
input_dir is the folder path where the reads of the transcriptome double-end sequencing are saved.  
reference_file is the save path of the fasta file of the species.  
annotation_file is the "eggnog.emapper.decorated.gff.saf" file generated in step 3.  
Note: please create output_dir in advance.  

Usage:
```
bash egg_rice_rnaseq_part.sh <input_dir> <output_dir> <reference_file> <annotation_file>

example:
nohup bash  /your/path/of/scripts/egg_rice_rnaseq_part.sh  /your/path/of/clean_data/  /your/path/of/output/  /your/path/of/ref/all.chrs.con  /your/path/of/ref/eggnog.emapper.decorated.gff.saf &

```


**6.Prepare Genebank files**  
**script.py** converts genome data in FASTA format and gene annotation data in GFF3 format into a GenBank file format.   
**script.py** reads the FASTA file for chromosome sequences and the GFF3 file for gene features (mRNA and CDS), then associates these features with the chromosomes in the GenBank format.  

Purpose:  
To standardize chromosome IDs.  
To parse the GFF3 file for mRNA and CDS annotations.  
To generate a GenBank file containing mRNA and CDS features for each chromosome in the FASTA file.  
It is designed for plant genome data processing but can be used for any genome with appropriate input data.  

Input:  
FASTA file: Contains genome sequence data for multiple chromosomes. The file should have the sequence data in FASTA format with chromosome IDs.  
GFF3 file: Contains annotation data for the genome, including mRNA and CDS features. The file should be in GFF3 format.  

Usage:
```
python script.py <input.fasta> <input.gff3> <output.gbk>

example:
python scripts/create_chr_genbank.py  ./ref/all.chrs.fa  ./ref/all.gff3  ./ref/all.gbk
```

Output:  
The script produces a GenBank file (.gbk), where each chromosome has its corresponding mRNA and CDS features annotated.  
The output file will be saved with the name specified as <output.gbk>.  


**6.Run deepbgc**  
Next, we use DeepBGC to predict biosynthetic gene clusters (BGCs) and use deep learning models to identify BGCs from genomic data. The input is the GenBank format file (all.gbk) generated in the previous step, which contains the genomic data to be used for DeepBGC pipeline operations.  

```
deepbgc pipeline <input.gbk> -o <output_directory> > <log_file> 2>&1

example:
nohup deepbgc pipeline ./ref/all.gbk  -o your/output/path/ > .your/output/log_files/deepbgc.log 2>&1 &

```









