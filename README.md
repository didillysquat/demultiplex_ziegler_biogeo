# demultiplex_ziegler_biogeo
This is the code used to demultiplex the .sff files that were generated by 454 sequencing for the Ziegler et al project
https://onlinelibrary.wiley.com/doi/full/10.1111/jbi.12913.
There are three scripts, one to demultiplex, on to create fake fastq R1 R2 files and one to populated the symportal datasheet.
I will aim to upload the single read fastq files to NCBI to replace the old reads.
To do the conversion from .sff to fastq, prior to running the demultiplexing python script we used: https://github.com/indraniel/sff2fastq
