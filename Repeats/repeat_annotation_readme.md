## Repeat Annotation and Window Analysis Workflow

# Custom Repeat Library

1. Run repeatmodeler on assembly lacking HiC scaffolding (c05_colaptes_arcs3.FINAL.reordered_2.fasta) to identify repeats de novo. Prior to running repeatmodeler, build database of reference. 
2. Run repeatmodeler. 
3. Build database of vertebrate repbase. 
4. Rmblast the repeatmodeler output to vertebrate repbase database. 
5. Remove any repeatmodeler sequences that are >= 98% to repbase sequences. 
6. Make a blast database of the colaptes assembly. 
7. Blast repeatmodeler output to colaptes assembly and take top 100 hits. 
8. Filter top 50 hits for each sequence and create a bed file to use for extraction of genomic sequences and flanks, using colaptes_filter_repeatmodeler.blast.R. 
9. Use bedtools to extract fasta sequences from the reference genome. 
10. Take all the fasta files, input into geneious and make alignments with mafft to develop consensus sequences. 

# Repeat Masking

11. With custom repeat library, use RepeatMasker to assess transposable element content in genome. 
12. From RepeatMasker output, create repeat landscape using caldivergencefromalign.pl and createrepeatlandscape.pl. 
13. Prior to window analysis, remove TE overlaps using summarize_repeatmasker_remove_overlaps.R. 

# Repeat Window Analysis

14. With input from summarize_repeatmasker_remove_overlaps.R, run colaptes_te_window_analysis.R to get TE content along non-overlapping 100 kb windows. 
15. Using output from colaptes_te_window_analysis.R run plot_window_te.R to produce plot of TE content along scaffolds. 
