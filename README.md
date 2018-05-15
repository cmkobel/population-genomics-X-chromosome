This repo is the hand in of a semester end project in the course: Population Genomics spring 2018, at BiRC, Aarhus University.


***

## Population Genetics on X-chromosome 


The data consists of a vcf file of 150 male full X chromosomes, a bed file with callable regions, a gif gene annotation file, a metafile with information about the samples and a set of files for use with REHH.

Gene annotation:

Gene annotation (gtf format) for Hg19 can be found in the following website
https://www.gencodegenes.org/releases/17.html
It was also uploaded to the dropbox as **gencode.v17.annotation.gtf**

Fst Calculation:

You will do the analysis from scratch by reading the genotype file of each population into different tables (rememeber rows are SNP positions and columns are individuals), the information about the snps are in the .snp file (ancestral and derived alleles). The data is haploid (n) therefore calculating Fst consists of estimating the allele frequencies for each position and calculating the expected heterozygosity within population Hs and contrasting Expected Heterozygosity across populations Ht. 

Fst = (Ht - Hs)/Ht. 

This can be done by averaging Fst values for a set of consecutive markers in a given window size (100 SNPs).


### Investigate the following

A. Perform an Fst scan between sets of populations in a sliding window of 100 SNP positions, including at least the contrast between Africa and Europe, between Europe and East Asia, and between East Asia and Africa. Identify the 10 strongest Fst outlier regions in each case. Identify their genomic position and the genes covered by these Fst peaks. Discuss potential adaptive explanations. 

B. Perform an iHS scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A.

C. Perform an XP-EHH scan of the whole X chromosome for at least three populations. Identify the 10 most significant regions and associated with genes as in A. 

D. Intersect the analysis of Fst and XP-EHH

E. Perform any additional analysis of your own choice, such as (diversity along the C X chromosome)

### Data

A DropBox link uploaded on BlackBoard (./Materials/Week 12: Projects Materials)
