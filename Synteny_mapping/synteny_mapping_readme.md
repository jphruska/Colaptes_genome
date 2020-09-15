## Synteny mapping workflow

1. Aligned Caur_TTU_1.0 assembly to Chicken genome (NCBI GCF_000002315.6), using nucmer module of MUMMER (01_nucmer.sh). 
2. Filtered alignments using delta-filter module while setting the minimum alignment identity to 70% and allowing many-to-many alignments (02_filter.sh). Rename scaffolds according to the Chicken chromosome they are syntenous with.
3. Produce tab-delimited text file including information on the position, percent identity, and length of each alignment using MUMMER's show-coords module (03_show_coords.sh). 
4. Create synteny plot with OmicCircos (04_create_circos_synteny_plot.R). 
