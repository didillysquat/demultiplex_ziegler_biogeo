"""
A quick script to demultiplex the .sff files for the ziegler ITS2 biogeography paper.
"""
import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

base_dir = "/home/humebc/projects/20210413_ziegler"
dir_list = [_[0] for _ in os.walk(base_dir) if "rawdata" in _[0] and "OM1" not in _[0] and "AD" not in _[0] and "RS1" not in _[0]]
for target_dir in dir_list:
    print(f"\n\nProcessing {target_dir}")
    # Do the demultiplexing for each of the regions
    base_name = target_dir.split("/")[-1].split("_")[0]
    fastq_path = os.path.join(target_dir, f"{base_name}.fastq")
    barcode_df = pd.read_csv(os.path.join(target_dir, f"{base_name}.barcodes.csv"))
    code_to_sample_dict = {barcode_df.loc[ind, "Barcodes"]: barcode_df.loc[ind, "Sample names"] for ind in barcode_df.index.values}
    foo = "bar"

    # Dictionary to hold a list of the SeqIO records per sample
    demultiplexed_dict = defaultdict(list)

    # Go barcode by barcode
    # keep the sequence if none of the other barcodes are found in it
    # Let's try for an exact match to start with and we can see how many of the sequences we lose
    seq_record_list = list(SeqIO.parse(fastq_path, "fastq"))
    
    conflict_count = 0
    barcode_count = 0
    for barcode, sample in code_to_sample_dict.items():
        barcode_count += 1
        conflict_count_barcode = 0
        found = []
        to_search = []
        other_barcodes = [_ for _ in code_to_sample_dict.keys() if _ != barcode]
        for record in seq_record_list:
            first_20 = record.seq[:20]
            # If the barcode in found in the seq but not any of the other samples
            conflict = False
            if barcode in first_20:
                for other_barcode in other_barcodes:
                    if other_barcode in first_20:
                        conflict = True
                        break
                if not conflict:
                    demultiplexed_dict[sample].append(record)
                    found.append(record)
                else:
                    conflict_count += 1
                    conflict_count_barcode += 1
                    to_search.append(record)
            else:
                to_search.append(record)
        # remove the found records from the seq_record_list
        print(f"Found {len(found)} sequences for {sample}; ({conflict_count_barcode} conflict seqs); {len(to_search)} remaining to search; {barcode_count}/{len(code_to_sample_dict)}")
        seq_record_list = to_search
        # Write out the found seqs as fastq
        os.makedirs(os.path.join(target_dir, "sample_fastqs"), exist_ok=True)
        SeqIO.write(found, os.path.join(target_dir, "sample_fastqs", f"{sample}.fastq"), "fastq")
        

    # Here we have the individual fastq contents in the demultiplexed_dict
    # and we can output how many seqs were not found and how many had a conflict
    foo = "bar"
    
                    

            

    
foo = "bar"