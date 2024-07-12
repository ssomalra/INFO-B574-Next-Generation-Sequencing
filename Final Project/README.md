# Genomic Insights into ATF3 Regulation of Cell Proliferation of HepG2 Cells
A collaborative project by: Sahiti Somalraju, Chandrima Modak, Rameesha Syed

## Abstract
Previous research has indicated that increased expression of the Activating Transcription Factor 3 (ATF3) protein contributes to the inhibition of cell proliferation, migration, and growth in HepG2 hepatocarcinoma cells. In our investigation, we conducted ChIP-seq analysis to uncover the genomic loci to which ATF3 binds, aiming to elucidate the genes involved in the tumor-suppressive cascade. Our results revealed a specific ATF3 binding motif, TGACATCA, within the promoter regions of the genome. Furthermore, our study identified novel genes, including MAP2K3, which is involved in the cardiovascular fibrosis pathway response.

## Introduction
Hepatocellular carcinoma (HCC) stands as a significant subclass of cancer, constituting a substantial portion of the health-cost burden, estimated at 55% (Sung et al., 2020). However, the formidable challenge lies in the drug resistance exhibited by hepatocarcinoma cells, hindering the attainment of cancer-free states post-treatment. One of the strategies proposed to regulate tumor cell suppression is Activating Transcription Factor (ATF3) {cite Li here too}. ATF3 plays a crucial role in cellular stress response, significantly influencing metabolism, immunity, and tumor progression. Recognized as a master modulator in maintaining metabolic homeostasis, ATF3 functions as a central hub in coordinating cellular adaptive responses (Gilchirst et al., 2010). Research conducted by  Li et al (2019),. revealed that overexpression of activating transcription factor 3 (ATF3) resulted in decreased cell proliferation, elevated cell apoptosis rates, and inhibition of cell cycle progression. In this study, we employed ChIP-seq (Chromatin Immunoprecipitation followed by Next-Generation Sequencing) to comprehensively analyze the genome-wide binding activity of ATF3 in HepG2 cells. ChIP-seq enabled high-resolution mapping of ATF3 binding sites, providing a detailed understanding of its genomic targets. By characterizing the regulatory landscape of ATF3 in HepG2 cells, our genomic data analysis uncovered novel genes interacting with ATF3, and elucidated novel pathways, thereby contributing to the broader understanding of the HCC. 

## Methods
### Data Retrieval
The study utilizes ChIP-seq data gathered from two sets of HepG2 cells subjected to ATF3 antibody pulldown, sourced from ENCODE under the accession numbers [`ENCFF522PUA`](https://www.encodeproject.org/files/ENCFF522PUA/) and [`ENCFF094LXX`](https://www.encodeproject.org/files/ENCFF094LXX). In parallel, two sets of control datasets lacking a specific target, [`ENCFF247MLF`](https://www.encodeproject.org/files/ENCFF247MLF) and [`ENCFF086CF`](https://www.encodeproject.org/files/ENCFF086CFB), were collected to establish a baseline. Sequencing was conducted on the Illumina HiSeq 2000 platform, generating 50nt single-end short-reads.

### ChIP-seq Analysis
The computation workflow for the ChIP-seq experiment is highlighted in Figure 1 below and described in detail.

<p align="center">
  <img width="750" alt="pipeline" src="https://github.com/user-attachments/assets/c19ddc0d-b442-4c9a-b1a6-d4ef1c9dc431">
</p>
