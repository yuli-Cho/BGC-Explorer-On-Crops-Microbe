#!/usr/bin/env python3

import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from datetime import datetime
import re

def translate_cds(sequence, strand, phase):
    """Translate CDS sequence to amino acids considering strand and phase"""
    # Adjust sequence based on phase
    if phase not in [0, 1, 2]:
        phase = 0
    seq = sequence[phase:]
    
    # Ensure sequence length is multiple of 3
    seq = seq[:len(seq) - (len(seq) % 3)]
    
    # Reverse complement if on minus strand
    if strand == -1:
        seq = str(Seq(seq).reverse_complement())
    
    try:
        # Translate to amino acids
        protein = str(Seq(seq).translate(table=1, to_stop=True))
        return protein
    except:
        return None

def parse_gff(gff_file):
    """Parse GFF file and return features grouped by chromosome"""
    features_by_chr = {}
    
    with open(gff_file) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
                
            parts = line.strip().split('\t')
            if len(parts) != 9 or parts[2].lower() != 'cds':
                continue
                
            chr_id = parts[0]
            if chr_id not in features_by_chr:
                features_by_chr[chr_id] = []
            
            start = int(parts[3]) - 1  # Convert to 0-based
            end = int(parts[4])
            strand = 1 if parts[6] == '+' else -1
            phase = int(parts[7]) if parts[7].isdigit() else 0
            
            # Parse attributes
            attrs = {}
            for attr in parts[8].split(';'):
                if '=' in attr:
                    key, value = attr.split('=', 1)
                    attrs[key] = value
            
            feature = {
                'start': start,
                'end': end,
                'strand': strand,
                'phase': phase,
                'attrs': attrs
            }
            features_by_chr[chr_id].append(feature)
    
    return features_by_chr

def create_gbk_record(chr_id, sequence, features, user_info):
    """Create a GenBank record for a chromosome"""
    # Create the basic record
    record = SeqRecord(
        Seq(sequence),
        id=chr_id,
        name=chr_id,
        description=f"Converted from eggnog GFF for chromosome {chr_id}"
    )
    
    # Add required annotations
    record.annotations["molecule_type"] = "DNA"
    record.annotations["topology"] = "linear"
    record.annotations["date"] = datetime.now().strftime("%d-%b-%Y").upper()
    record.annotations["data_file_division"] = "BCT"
    record.annotations["source"] = "Converted from eggnog GFF"
    record.annotations["comment"] = f"Converted by {user_info}"
    
    # Add features
    for feat in features:
        qualifiers = {}
        
        # Add standard qualifiers
        if 'ID' in feat['attrs']:
            qualifiers['locus_tag'] = feat['attrs']['ID']
        if 'Name' in feat['attrs']:
            qualifiers['gene'] = feat['attrs']['Name']
        if 'protein_id' in feat['attrs']:
            qualifiers['protein_id'] = feat['attrs']['protein_id']
            
        # Add eggnog specific qualifiers
        if 'eggNOG_OGs' in feat['attrs']:
            qualifiers['note'] = f"eggNOG_OGs: {feat['attrs']['eggNOG_OGs']}"
        
        # Extract and translate the CDS sequence
        start, end = feat['start'], feat['end']
        cds_seq = sequence[start:end]
        protein_seq = translate_cds(cds_seq, feat['strand'], feat['phase'])
        
        if protein_seq:
            qualifiers['translation'] = protein_seq
        
        # Set codon start based on phase
        qualifiers['codon_start'] = feat['phase'] + 1
        
        # Create the feature
        feature = SeqFeature(
            FeatureLocation(feat['start'], feat['end'], feat['strand']),
            type='CDS',
            qualifiers=qualifiers
        )
        
        record.features.append(feature)
    
    return record

def main():
    parser = argparse.ArgumentParser(description='Convert eggnog GFF to GenBank format')
    parser.add_argument('gff_file', help='Input GFF file from eggnog')
    parser.add_argument('fasta_file', help='Input FASTA file with sequences')
    parser.add_argument('output_file', help='Output GenBank file')
    args = parser.parse_args()
    
    # Get current user and time information
    current_time = "2025-02-09 05:43:20"  # Using the provided time
    user_info = "HHeng-bioinfo"  # Using the provided user login
    user_string = f"Conversion performed by {user_info} at {current_time} UTC"
    
    # Parse GFF file
    features_by_chr = parse_gff(args.gff_file)
    
    # Read sequences from FASTA
    sequences = SeqIO.to_dict(SeqIO.parse(args.fasta_file, "fasta"))
    
    # Create GenBank records
    records = []
    for chr_id in features_by_chr:
        if chr_id in sequences:
            record = create_gbk_record(
                chr_id,
                str(sequences[chr_id].seq),
                features_by_chr[chr_id],
                user_string
            )
            records.append(record)
    
    # Write output
    SeqIO.write(records, args.output_file, "genbank")
    print(f"Conversion completed. Processed {len(records)} sequences.")
    print(f"Output written to: {args.output_file}")

if __name__ == "__main__":
    main()