## Pacific cod, Genomic Comparisons of Korean Peninsula and US/CAN coast populations
<br>
<br>

**Purpose:** This repository contains analyses conducted to compare RAD-sequencing data sets of east (US/CAN west coast) and west (Korean peninsula) populations of Pacific cod. 
<br>
In short, raw RAD-sequencing data was re-analyzed using the stacks pipeline and a *de novo* reference assembled using individuals from the Korean peninsula data set. After filtering for polymorphic SNPs, minor allele frequency, missing data, and Hardy-Weinberg equilibrium, the consensus sequences of the remaining loci were aligned to the Atlantic cod genome.
<br>
<br>
<br>

### Directory Structure: 

**analyses**: Analysis conducted after data processing and filtering. Includes: 
1. steps for alignment to the Atlantic cod genome (gadMor2) and corresponding annotation file
2. outlier analysis
3. sliding window analysis of Fst across the genome, with visualization 
4. exploratory Manhattan Plots 


**notebooks**: Jupyter notebooks. Data processing in stacks and subsequent filtering steps are held in `pipeline` folders. 


**stacks_b8_wgenome_r05**: Genepop files from the eight batch of stacks, run with a low stack depth (-m 3) and a low proportion of individuals per population required to process a locus (-r 0.5)


**stacks_pipeline_analyses**: This analysis encountered a major roadblock from missing data in 1/2 of the samples used for the study. An extensive exploration of the stacks pipeline was completed to understand why this was occurring. Files from those exploratory analyses can be found in this folder. 



<br>
<br>
*This repo had to be recreated on 3/24/2018 because a large file was pushed to github.*