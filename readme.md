RNA-seq data analysis
=======================

This project is an analysis of RNAseq data of mouse Gcn1 and Atf4 knockdown cells and WT control cells.


Data source
-----------

The raw data for this project comes from Jonanthan Hernandez Lab at NIH

- RawData/filteredCounts has the filtered gene counts generated using [CCBR Pipeliner](https://ccbr.github.io/Pipeliner/). Genes with >= 1 CPM in >= 2 samples are included.
- RawData/result_webScrap.txt is gene Entrez IDs scrapped from NCBI using /analysis/spath_forDa.py
- RawData/filteredCounts/sampletable.txt also has metadata including sampleName, fileName, condition and label.


Objectives
--------

- Generate differential expression using Limma Voom.
- Perform over-representation analysis(ORA) and Gene Set Enrichment Analysis (GSEA) using Clusterprofiler.
- Investigate the genes and pathways regulated by Gcn1 and Atf4 either independently or together (overlap).


Analysis code is written for R version 3.6.1

---

Summary
--------
Sinha and Hernandez et al. identified the general control of amino acid synthesis 1-like 1 gene (Gcn1) in an in vivo genetic screen to isolate drivers of metastatic outgrowth. Preliminary in vitro and in vivo experiments suggest that Gcn1 drives tumor growth and metastasis during glutamine deprivation by activating the master regulator of the integrated stress response pathway ATF4. 

This dataset has 6 samples: two control sample (SCR) and four Gcn1 knockdown samples (SH1 and SH2) under 0.5 mM treatment. First the knock-down of Gcn1 is verified by importing bigwig/bam files into IGV and checking the exons of Gcn1 gene. SH1 and SH2 has lower counts on the first few exons. RSEM counts were normalized using the voom algorithm from the Limma R package (version 3.40.6). These normalized counts were used for clustering and data visualization. PCA shows SH1 knockdown samples have larger variation compared to SCR controls than SH2 compared to SCR. Heatmap generated using the top 500 most variable genes across samples visually confirmed this. 

The Limma package was used to test for differential gene expression between experimental conditions. Significant differentially expressed genes were identified with a false-discovery rate ≤ 0.05 and a log(fold change) > 1. The number of DEGs in SH1_0pt5mM-SCR_0pt5mM: 3300; SH2_0pt5mM-SCR_0pt5mM: 1028; SCR_4mM-SCR_0pt5mM: 244. Since SCR_4mM-SCR_0pt5mM serves as negative control, the DEGs in this contrast that overlap with two other contrasts were excluded for subsequent analysis. 

Over Representation Analysis (ORA) is used to identify biological functions or processes are over-represented (enriched) in differentially expressed genes (DEGs).Pathway enrichment was performed using pre-ranked Gene Set Enrichment Analysis (GSEA) with the KEGG and Reactome genesets from the Molecular Signatures Database. Differentially expressed genes that led to the perturbation of aminoacyl tRNA biosynthesis and serine family amino acid biosynthetic process were identified.
