Hi-C scaffolding of existing colaptes reference genome. 

Workflow: 

1. Run juicer.sh, using reference genome c05_colaptes_arcs.fasta found in /lustre/work/johruska/juicer/references/. 
2. Run r script to take split* files and create merged_nodups.txt file. This process does not complete. Use merged_no_dups.R script to create no_dup.txt files. 
3. Merge no_dups.txt files into merged_nodups.txt with the following command: for i in {1..238}; do cat no_dupes_${i}.txt >> merged_nodups.txt; done. 
4. Run 3d-dna.sh on merged_nodups.txt file, with variable --editor-repeat-coverage flags (3,4,5). Output each to a directory called working${1,2,3}. 
5. Import appropriate .assembly and .hic files into Juicebox for chromosome/scaffold manipulation. Use heat map to correct for misassamblies and to construct putative "chromosomes". 
6. Export modified .assembly file from Juicebox, and re-input into 3d-dna. Run the asm pipeline post review tool to input the manually edited assembly file and create the final
fasta file and hic file for visualization. 
7. Calculate physical distance coverage using the steps outlined in physical_distance_coverage_workflow.txt. 

