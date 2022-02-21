# Allele specific expression analysis.

# Motivation:

ASE quantifies the difference in expression of two alleles of a gene and can be measured using RNA-seq reads that align to heterozygous sites. Compared to standard differential gene expression analysis, ASE is insensitive to environmental or trans-acting factors, which generally affect both alleles equally. This makes ASE a powerful tool for revealing genes that are affected by cis-acting mutations, including noncoding regulatory mutations that affect sequences such as promoters, enhancers, and insulators as well as protein-coding or splicing mutations that result in nonsense-mediated decay. Another advantage of ASE is that it is detectable even when the identity of the pathogenic variants causing dysregulation are unknown and it can reveal the effects of rare germline or somatic mutation. Thus, ASE is a powerful tool for the identification of genes with altered gene dosage due to cis-acting genome alterations.

# Pre-processing
We filtered RNA-seq reads for mapping bias and obtained allele specific read counts at heterozygous positions using WASP (https://github.com/bmvdgeijn/WASP). WASP uses random sampling to ensure that allele counts at nearby heterozygous sites are independent. Specifically, when a read overlaps multiple sites, the allelic count is incremented at only one of the sites, which is selected randomly. We observed that at most heterozygous sites overlapping RNA-seq reads only match the reference or alternate alleles. However, some sites also have reads that match neither allele, which we refer to as ‘other’ reads. ‘Other’ reads may reflect sequencing errors or mis-mapped reads so we removed all heterozygous sites where the ‘other’ read count was greater than two prior to ASE analysis. We only considered heterozygous sites with atleast 10 reads.

# Method description
To estimate ASE in these samples, we implemented a statistical model that utilizes allele-specific read counts at heterozygous sites mapping to exons of genes, while accounting for genotyping errors, sequencing errors, overdispersion of RNA-seq read counts and phasing of alleles on the same gene. This model estimates allele imbalance for each gene, which is how far the reference allele proportion differs from the expected value of 0.5. For more details on modelling read counts please refer to Sen et al, 2022, Cancer Research (DOI: https://doi.org/10.1158/0008-5472.CAN-21-0810) and Sen et al, 2022, Genome Biology (DOI:pending). This GitHub contains scripts from Genome Biology paper.

# For running rna_allele_imbalance.R:

R version 4.1.2 (2021-11-01), Platform: x86_64-conda-linux-gnu (64-bit), Running under: CentOS Linux 7 (Core), Dependencies : rmutil_1.1.5, dplyr_1.0.7, plyr_1.8.6, data.table_1.14.2

~/Rscript --vanilla rna_allele_imbalance.R rna_allele_count.csv ase.csv

Example of rna_allele_count.csv and ase.csv files are provided in the Github.

# DNA allele imbalance analysis.

# Motivation
While several existing tools leverage read depth to predict SCNAs, these methods have limited precision and report many false positive focal SCNAs. In Sen et al, 2022 Genome Biology we demonstrated that quantifying difference in allele-imbalance in DNA-seq reads between tumor and normal samples can be a powerful method to detect SCNAs.

# Pre-processing
We realigned exome-seq reads from tumor and normal samples to the reference genome using BWA-ALN and filtered the sequencing reads to remove mapping bias using WASP. We obtained allele-specific read counts at heterozygous sites (excluding multiallelic sites) for normal tissues using the CollectAllelicCounts from GATK (4.1.1). We assumed that most heterozygous sites in normal tissues are germline polymorphisms and obtained allele-specific read counts at the shared positions for matched tumor samples. We analyzed shared heterozygous positions because this facilitates direct comparison of reference allele proportions between tumor and paired normal tissues. To model DNA allelic imbalance over large genomic segments, we sorted exons by their genomic coordinates and grouped 20 consecutive exons into genomic bins. Next, we assigned the heterozygous sites which were covered by at least 10 reads to the genomic bins. To ensure robust regional DNA allele imbalance estimates, we retained genomic bins with at least 10 heterozygous sites for DNA allele imbalance analysis

# Method description
We estimate allele imbalance for tumor (a_tumor) and normal (a_normal) using method adapted from ASE analysis. We then calculate the difference in a between tumor and normal samples, δ_a as:

δ_a=|a_tumor|-|a_normal|

Finally, to create contiguous segments of allelic imbalance, we performed Circular Binary Segmentation (CBS) on δ_a. For detailed description of the method refer to Sen et al, 2022, Genome Biology (DOI:pending) 

# For running dna_allele_imbalance.R:

R version 4.1.2 (2021-11-01), Platform: x86_64-conda-linux-gnu (64-bit), Running under: CentOS Linux 7 (Core), Dependencies : rmutil_1.1.5, dplyr_1.0.7, plyr_1.8.6, data.table_1.14.2, rtracklayer_1.54.0 

~/Rscript --vanilla ~/dna_allele_imbalance.R dna_allele_count_tumor.tsv dna_test_regions.bed tumor.csv

~/Rscript --vanilla ~/dna_allele_imbalance.R dna_allele_count_normal.tsv dna_test_regions.bed normal.csv

Run dna_allele_imbalance.R script seperately for tumor and normal DNA-seq data. Example of dna_allele_count_tumor.tsv, dna_allele_count_normal.tsv, dna_test_region.bed, tumor.csv and normal.csv are provided in this Github. With the results file (tumor.csv and normal.csv) run the segmentation.R to perform circular binary segmentation (CBS).

# For running segmentation.R:

R version 4.1.2 (2021-11-01), Platform: x86_64-conda-linux-gnu (64-bit), Running under: CentOS Linux 7 (Core), Dependencies : rmutil_1.1.5, dplyr_1.0.7, plyr_1.8.6, data.table_1.14.2, DNAcopy_1.64.0

~/Rscript --vanilla ~/segmentation.R tumor.csv normal.csv segmented.csv

Example of the expected output is provided in the Github.

















