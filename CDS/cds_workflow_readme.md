# Workflow describing gene discovery and annotation, followed by mutation rate analysis. 

## Gene discovery and annotation

1. Run MAKER with protein Apaloderma, Buceros, Merops and Picoides protein datasets. 
2. Train SNAP with round 1 MAKER output. 
3. Run Augustus with the round 1 MAKER output. 
4. Run second round of MAKER with Augustus and SNAP. 
5. Prepare gff files for round 2 of MAKER. 
6. Rename MAKER round 2 output. 
7. Extract CDS (colaptes_cds.fasta). 

## Mutation rate analysis

8. Download CDS for species representing closely related orders and perform reciprocal BLAST of all species vs Colaptes. 
9. Download phylogenetic tree comprising all orders of Neoaves (Jarvis et al. 2014) and prune tree to the four representative orders (Bucerotiformes, Piciformes, Trogoniformes and Coraciiformes). 
10. Use T-Coffee to align the putative homologues betweeen the four passerine species. 
11. Prior to back-translating, remove any gaps in the protein alignment with trimAl. 
12. With alignments, test for selection using the gene-wide and branch-specific tests for selection in CODEML. Remove alignments with gene-wide or branch-specific selection, after correction using Benjamini and Hochberg. 
13. Extract four-fold degenerate sites from gene alignments that passed selection filter (step 12) with rphast, Biostrings, and seqinr. 
14. Concatenate four-fold sites and use jModelTest2 to determine an appropriate model of sequence evolution. 
15. Use appropriate model of sequence evolution and user-specified tree (Jarvis et al. 2014) to estimate branch lengths based on four-fold degenerate sites. 
16. Use Colaptes-specific branch length of resulting tree to estimate a mean and 95% HPD distribution of Colaptes-specific mutation rates. 

## Window analysis

17. Reduce colaptes_cds.fasta file so it only includes information on the the start and end points of cds content on the scaffolds of interest (in this case the scaffolds that showed strong synteny to Gallus assembly). 
18. Use the subsequent reduced cds content file, run colaptes_gene_window.R to obtain cds content (colaptes_cds_props.txt) along non-overlapping 100 kb windows. 
19. Use colaptes_cds_props.txt as input for plot_window_CDS.R to produce plot graph of CDS distribution across assembly. 
