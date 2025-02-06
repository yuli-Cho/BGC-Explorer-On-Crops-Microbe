## 
This script converts genome data in FASTA format and gene annotation data in GFF3 format into a GenBank file format.   
The script reads the FASTA file for chromosome sequences and the GFF3 file for gene features (mRNA and CDS), then associates these features with the chromosomes in the GenBank format.  

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

