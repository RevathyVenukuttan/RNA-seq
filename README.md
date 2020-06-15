## Introduction:

Alternative Splicing or Differential splicing is a regulated process during gene expression that results in single gene coding for multiple proteins. In this process, multiple messenger RNAs (mRNAs) can be generated by joining exons together in different combinations. This process increases both protein diversity and regulate gene expression. Alternative splicing is generally controlled by proteins that bind directly to regulatory sequence elements and either activate or repress splicing of adjacent splice sites in a target pre-mRNA. In humans, NOVA1 and NOVA2 are the genes which regulate the alternative splicing activity. To understand better which other genes are involved/affected, a knockdown experiment is performed where Pasilla, a drosophila gene that regulates alternative splicing is knockdown as Brook et al., showed that the RNA regulatory map of Pasilla and NOVA1, NOVA2 are highly conserved between insects and humans. 

## Goal:

The project aims at identifying differential expression in genes and the associated gene ontology affected by knockdown of Pasilla in Drosophila, the gene that regulates the process of alternate splicing. 

## Data and Preprocessing:

Raw data files contains 4 untreated and 3 treated RNAi samples where 2 among untreated and 2 among treated are paired end sequenced and 2 among untreated and 1 among treated are single end sequenced samples. 
Reads were aligned to Drosophila reference genome using TopHat aligner and reads were summarized at gene-level using htseq-count from HTSeq.

## Steps of Analysis:

1. Filtering out of lowly expressed genes using a threshold value based on read counts (2 or 3 cpm in atleast 3 samples)
2. Filtered genes are converted into a DGEList object
3. Hierarchical clustering is performed to visualize variable genes in the list
4. Normalization performed using Trimmed Mean of Means (TMM) method performed on the DGEList
5. A linear additive model is used for Differential Expression Analysis using the function lmFit() and eBayes() correction used for optimizing the fit
6. Results from the Differential Expression analysis is annotated to their gene names and ENTREZ IDs using FLYBASE data from org.Dm.eg.db database.
7. goana() is used to identify the GO term associated with the differentially expressed genes
