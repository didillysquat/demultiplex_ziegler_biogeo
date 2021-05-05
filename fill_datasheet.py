"""
We have all of the fastq file in the one directory, now we want to fill out the .csv with the correct file
names and remove those samples that are missing files.
"""

import pandas as pd
import os

df = pd.read_csv("/home/humebc/projects/20210413_ziegler/fastqs_for_symportal_submission/ziegler_datasheet.short.csv")
df = df.astype({"fastq_fwd_file_name": str, "fastq_rev_file_name": str})
file_dir = "/home/humebc/projects/20210413_ziegler/fastqs_for_symportal_submission/"

foo = "bar"

# for index in the df, check for the file based on the sample name.
# Check that only two files found
# delete if not found
keep = []
for ind in df.index:
    sample_name = df.at[ind, 'sample_name']
    fwd_read_path = os.path.join(file_dir, f"{sample_name}.R1.fastq.gz")
    rev_read_path = os.path.join(file_dir, f"{sample_name}.R2.fastq.gz")
    if os.path.exists(fwd_read_path) and os.path.exists(rev_read_path):
        df.at[ind, "fastq_fwd_file_name"] = f"{sample_name}.R1.fastq.gz"
        df.at[ind, "fastq_rev_file_name"] = f"{sample_name}.R2.fastq.gz"
        keep.append(ind)
remove = [_ for _ in df.index if _ not in keep]
# remove the rows to remove
df = df.loc[keep,]
df.to_csv("/home/humebc/projects/20210413_ziegler/fastqs_for_symportal_submission/ziegler_datasheet.short.complete.csv", index=False)
print(f"The following {len(remove)} samples had no valid files:")
for ind in remove:
    print(f"\t{ind}")
