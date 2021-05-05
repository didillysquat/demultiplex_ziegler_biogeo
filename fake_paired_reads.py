"""
Now that we have demultiplexed the sff files, we need to create fwd and rev files to
feed into SymPortal.

We should get the average length of the seqs and then aim to have a 100 bp overlap.
The quality scores should all be high
"""

import os
import pandas as pd
from Bio import SeqIO
from collections import defaultdict

base_dir = "/home/humebc/projects/20210413_ziegler"
dir_list = [_[0] for _ in os.walk(base_dir) if "rawdata" in _[0] and "sample" not in _[0]]
for target_dir in dir_list:
    print(f"\n\nProcessing {target_dir}")
    # Do the demultiplexing for each of the regions
    base_name = target_dir.split("/")[-1].split("_")[0]
    fastq_paths = [_ for _ in list(os.walk(os.path.join(target_dir, "sample_fastqs")))[0][2] if "R1" not in _ and "R2" not in _]
    # Do the generation
    
    for fastq_to_split in fastq_paths:
        print(f"splitting {fastq_to_split}")
        lengths = []
        fwd_records = []
        rev_records = []
        seq_record_list = list(SeqIO.parse(os.path.join(target_dir, "sample_fastqs", fastq_to_split), "fastq"))
        # For each record generate a fwd and rev seq if above 150bp in length
        for seq_record in seq_record_list:
            record_len = len(seq_record)
            lengths.append(record_len)
            if record_len > 150:
                fwd_records.append(seq_record[:int(record_len/2) + 25])
                rc_seq = seq_record.reverse_complement()[:int(record_len/2) + 25]
                rc_seq.id = seq_record.id
                rc_seq.name = seq_record.name
                rc_seq.description = seq_record.description
                rev_records.append(rc_seq)
        
        
        # Now write out the fwd and rev
        if fwd_records:
            print(f"average length was {sum(lengths)/len(lengths)}")
            SeqIO.write(fwd_records, os.path.join(target_dir, "sample_fastqs", fastq_to_split.replace(".fastq", ".R1.fastq")), "fastq")
            SeqIO.write(rev_records, os.path.join(target_dir, "sample_fastqs", fastq_to_split.replace(".fastq", ".R2.fastq")), "fastq")
        else:
            print(f"No seqs for {fastq_to_split}")

    foo = "bar"