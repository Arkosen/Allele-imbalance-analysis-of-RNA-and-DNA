# Allele imbalance analysis for RNA and DNA sequencing.

# Motivation:

ASE quantifies the difference in expression of two alleles of a gene and can be measured using RNA-seq reads that align to heterozygous sites. Compared to standard differential gene expression analysis, ASE is insensitive to environmental or trans-acting factors, which generally affect both alleles equally. This makes ASE a powerful tool for revealing genes that are affected by cis-acting mutations, including noncoding regulatory mutations that affect sequences such as promoters, enhancers, and insulators as well as protein-coding or splicing mutations that result in nonsense-mediated decay. Another advantage of ASE is that it is detectable even when the identity of the pathogenic variants causing dysregulation are unknown and it can reveal the effects of rare germline or somatic mutation. Thus, ASE is a powerful tool for the identification of genes with altered gene dosage due to cis-acting genome alterations.

# Method:

# Reference bias correction, read counts and filtering
We filtered RNA-seq reads for mapping bias and obtained allele specific read counts at heterozygous positions using WASP (https://github.com/bmvdgeijn/WASP). WASP uses random sampling to ensure that allele counts at nearby heterozygous sites are independent. Specifically, when a read overlaps multiple sites, the allelic count is incremented at only one of the sites, which is selected randomly. We observed that at most heterozygous sites overlapping RNA-seq reads only match the reference or alternate alleles. However, some sites also have reads that match neither allele, which we refer to as ‘other’ reads. ‘Other’ reads may reflect sequencing errors or mis-mapped reads so we removed all heterozygous sites where the ‘other’ read count was greater than two prior to ASE analysis. We only considered heterozygous sites with atleast 10 reads.

# Description
To estimate ASE in these samples, we implemented a statistical model that utilizes allele-specific read counts at heterozygous sites mapping to exons of genes, while accounting for genotyping errors, sequencing errors, overdispersion of RNA-seq read counts and phasing of alleles on the same gene. This model estimates allele imbalance for each gene, which is how far the reference allele proportion differs from the expected value of 0.5. For more details on modelling read counts please refer to Sen et al, 2022, Cancer Research (DOI: https://doi.org/10.1158/0008-5472.CAN-21-0810) and Sen et al, 2022, Genome Biology (DOI:pending). This GitHub contains scripts from Genome Biology paper.

# Instructions for running code.

# For running rna_allele_imbalance.R:

R version 4.1.2 (2021-11-01), Platform: x86_64-conda-linux-gnu (64-bit), Running under: CentOS Linux 7 (Core), Dependencies : rmutil_1.1.5, dplyr_1.0.7, plyr_1.8.6, data.table_1.14.2

~/Rscript --vanilla rna_allele_imbalance.R allele_count.csv result.csv

# Format of allele_count.csv input file for estimating RNA allele imbalance. Requires a comma seperated file with following columns. Column names MUST be included:

seqnames,start,end,ref,alt,snp_id,gene_id,ref.matches,alt.matches,N,errors,GQ,genotype.error;
1,889188,889188,G,T,1:889188_G/T,ENSG00000188976.10_2,88,0,88,0,69,1.26E-07
1,976242,976242,C,A,1:976242_C/A,ENSG00000188157.14_4,15,0,15,0,25,0.003162278
1,983546,983546,G,T,1:983546_G/T,ENSG00000188157.14_4,51,0,51,0,69,1.26E-07

# Expected output is a Comma seperated file with following columns:

gene_id,pval,a,fdr;
ENSG00000188157.14_4,3.55E-11,0.499933893,1.88E-07
ENSG00000034152.18_3,6.60E-10,0.499933893,1.75E-06
ENSG00000162591.15_2,6.24E-08,0.499933893,0.00011016

# For running dna_allele_imbalance.R:















