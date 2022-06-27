# rwc-and-wcHUH-seq

These two python scripts are variants of the original HUH-seq script (https://academic.oup.com/nar/article/49/2/1046/6067397#226585372), which was written in R.

These scripts use the Biopython package to parse next-generation sequencing (NGS) data for subsequent manipulation and quantitation.

These scripts were created to analyze data from experiments that are designed 

Briefly, wcHUH-seq takes in forward NGS reads, trims them down to the portion of interest, and counts the number of reads for each of sixteen potential sequences. This script takes in files for a reference, a magnesium reaction condition, and a manganese reaction condition -- each in triplicate. It compares the reads in reaction conditions to the reference and generates percent reductions through subtractive sequence (i.e., if a sequence occurs in the reference but not in the reaction condition it has been cleaved by the nuclease).

rwcHUH-seq takes in forward NGS reads, trims them down to the portion of interest, and counts the number of reads for each of sixteen potential sequences. This script takes in files for magnesium and manganese reaction conditions only -- both in triplicate. It calculates the percentage that each potential sequences makes up of the total number of reads to quantify reaction efficiency and preference in Rep-mediated reunion.

We have provided the Python code, raw data, processed data, and percent reductions/enrichments for each of these experiments.
