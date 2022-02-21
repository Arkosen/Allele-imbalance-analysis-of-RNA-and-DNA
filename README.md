# Allele Specific Expression (ASE)

Motivation:

ASE quantifies the difference in expression of two alleles of a gene and can be measured using RNA-seq reads that align to heterozygous sites. Compared to standard differential gene expression analysis, ASE is insensitive to environmental or trans-acting factors, which generally affect both alleles equally. This makes ASE a powerful tool for revealing genes that are affected by cis-acting mutations, including noncoding regulatory mutations that affect sequences such as promoters, enhancers, and insulators as well as protein-coding or splicing mutations that result in nonsense-mediated decay. Another advantage of ASE is that it is detectable even when the identity of the pathogenic variants causing dysregulation are unknown and it can reveal the effects of rare germline or somatic mutation. Thus, ASE is a powerful tool for the identification of genes with altered gene dosage due to cis-acting genome alterations.

Method:

Reference bias correction, read counts and filtering
We filtered RNA-seq reads for mapping bias and obtained allele specific read counts at heterozygous positions using WASP. WASP uses random sampling to ensure that allele counts at nearby heterozygous sites are independent. Specifically, when a read overlaps multiple sites, the allelic count is incremented at only one of the sites, which is selected randomly. We observed that at most heterozygous sites overlapping RNA-seq reads only match the reference or alternate alleles. However, some sites also have reads that match neither allele, which we refer to as ‘other’ reads. ‘Other’ reads may reflect sequencing errors or mis-mapped reads so we removed all heterozygous sites where the ‘other’ read count was greater than two prior to ASE analysis. 

Modelling read count data.
To estimate ASE in these samples, we implemented a statistical model that utilizes allele-specific read counts at heterozygous sites mapping to exons of genes, while accounting for genotyping errors, sequencing errors and overdispersion of RNA-seq read counts. This model estimates allele imbalance (aRNA) for each gene, which is how far the reference allele proportion differs from the expected value of 0.5. For more details on modelling read counts please refer to Sen et al, 2022, Cancer Research (DOI: https://doi.org/10.1158/0008-5472.CAN-21-0810) and Sen et al, 2022, Genome Biology (DOI:pending)
