Pilon polishing workflow. 

1. Cleaned raw chromium-10x reads and removed adapters with bbduk. 
2. Index input assembly (PacBio) and create sequence dictionary. 
3. Align cleaned chromium-10x reads to PacBio assembly with bwa-mem. Pipe output to samtools, and sort the resulting bam file. Bam file was also indexed with samtools using samtools index. 
4. Run pilon with PacBio assembly and bam file from Step 3. 
