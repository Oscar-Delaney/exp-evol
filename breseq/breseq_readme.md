# Bioinformatic analysis of genomic data using Breseq

We used two strains of A. baylyi: AB3, which is the ancestral strain capable of recombination, and AB13 which has a recombination knockout. The genomes of each are available in the 'genomes' folder.
As explained in the methods section of our paper, we performed metagenomic sequencing on the endpoint populations of the evolution experiment, and the resulting sequencing fastq files should be placed in the 'data' folder. However they are many gigabytes, and so are not stored in Github, and are instead on NCBI.
Running Breseq to determine the likely differences of the evolved strain from their ancestor takes a long time, and the outputs also include very large files, so these are also ommitted from the 'output' folder.
This data is then processed to isolate the most relevant results, which are saved in the 'analyses' folder.
The figures used in the paper are also shown in the 'figures' folder.
Within the 'code' folder, the workflow is as follows:
- run_Breseq.sh compares the sequence information in 'data' to the original sequences in 'genomes' and finds any likely differences, storing the results in 'outputs'. This code took about 5 days to run on our single machine.
- process_outputs.sh takes the files in 'outputs' and extracts and summarises the most useful information, saving it in 'analyses'
- mapping_success.sh is a helper file that is used to remove sequencing runs that do not map well onto the A. baylyi genome.
- make_Breseq_figures.R takes the simplified results in 'analyses' and further processes the data to be nicely graphable, and generates the figures.