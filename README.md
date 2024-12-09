# Proteomics_IME
Protein-based classification reveals an immune-hot subtype in IDH mutant astrocytoma with worse prognosis

## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk/)
* **Repository Development**: Jihong Tang

## Status 
Active Development 

## Last update time
TIME_UPDATE

## Introduction
This repository contains the code for the multi-omics and spatial single cell investigation of IDH-mutant astrocytoma. We analyzed MS-based Proteomics, bulk DNA sequencing, RNA sequencing, DNA methylation array, whole side imaging and single cell RNA-seq data, as well as 10x Visium HD spatial transcriptomics and CODEX multiplex imaging data. The data were used for protein-based clustering, survival analysis, multi-omics characterization, cell type annotation, and spatial analysis. We also developed an AI-aided multiomics  classifier based on multi-omics features. 


Code for these analysis were included in each corresponding folder. A small dataset is also provided to demo the code.

## System requirements and dependencies
Preprocessing of the raw sequencing data requires high performance clusters and workstations. We used a cluster equipped with Linux CentOS 7 (kernel version 3.10.0-1062.el7.x86_64) with 40 cores, 256 GB RAM and at least 100TB storage, and another workstation equipped with Windows with .

The other parts of the code can be run on a desktop. Required dependencies include R 4.3.2 - R 5.1.0 (for VisuimHD data) and python 3.9.6. The required R and python packages are included in each of the code snippets. We also used QuPath 0.5.0 which is downloadable at https://qupath.github.io/.

## Installation

To install, the user can either download a zip file from https://github.com/Jihong-Tang/Proteomics_IME, or by cloning the repository via
```
git clone https://github.com/Jihong-Tang/Proteomics_IME.git
```
To install the dependent R packages, please refer to the manual of each individual package.

## Demo
For each analysis, the user can enter the directory and run the R code inside the folder. The code will read the demo data provided in the folder (will be public upon publication of manuscript), run the analysis, and generate result files within the folder.

## Raw data availability
The mass spectrometry proteomics data have been deposited to the ProteomeXchange Consortium via the iProX partner repository (https://www.iprox.org/), under accession number Project ID: IPX0010405000 and will be publicly available as of the date of publication. Raw sequencing data of the newly sequenced samples will be deposited to Genome Sequence Archive (GSA, https://ngdc.cncb.ac.cn/gsa/) and publicly available as of the date of publication. Data from TCGA were downloaded from NCI Genomics Data Commons (GDC) data portal (https://portal.gdc.cancer.gov). Previously published CGGA data have been uploaded to the Genome Sequence Archive in BIG Data Center, Beijing Institute of Genomics (BIG), Chinese Academy of Sciences, under accession number BioProject ID: PRJCA001636 (https://bigd.big.ac.cn/bioproject/browse/PRJCA001636) and PRJCA001747 (https://bigd.big.ac.cn/bioproject/browse/PRJCA001747). All the other data supporting the findings of this study are available upon reasonable request.

## Contact
For any questions, please contact Professor Jiguang Wang via email: jgwang AT ust DOT hk
