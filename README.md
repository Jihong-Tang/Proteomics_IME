# Proteomics_IME
Protein-based classification reveals an immune-hot subtype in IDH-mutant astrocytoma with worse prognosis

## Ownership
[Wang Lab at HKUST](http://wang-lab.ust.hk/)
* **Repository Development**: Jihong Tang

## Status 
Active Development 

## Latest update
Last updated: Thu Sep 11 00:23:15 HKT 2025

## Introduction
This repository contains the code for the multi-omics and spatial single cell investigation of IDH-mutant astrocytoma. We analyzed MS-based Proteomics, bulk DNA sequencing, RNA sequencing, DNA methylation array, whole side imaging and single cell RNA-seq data, as well as 10x Visium HD spatial transcriptomics and CODEX multiplex imaging data. The data were used for protein-based clustering, survival analysis, multi-omics characterization, cell type annotation, and spatial analysis. We also developed an AI-aided multiomics  classifier based on multi-omics features. 


Code for these analysis were included in each corresponding folder. A small dataset is also provided to demo the code.

## System requirements and dependencies
Preprocessing of the raw sequencing data requires high performance clusters and workstations. We used a cluster equipped with Linux CentOS 7 (kernel version 3.10.0-1062.el7.x86_64) with 40 cores, 256 GB RAM and at least 100TB storage, and another workstation equipped with Windows for mass spectrametry data library searching.

The other parts of the code can be run on a desktop. Required dependencies include R 4.3.2 - R 5.1.0 (for VisuimHD data) and python 3.9.6. The required R and python packages are included in each of the code snippets. We also used QuPath 0.5.0 which is downloadable at https://qupath.github.io/.

## Installation

To install, the user can either download a zip file from https://github.com/Jihong-Tang/Proteomics_IME, or by cloning the repository via
```
git clone https://github.com/Jihong-Tang/Proteomics_IME.git
```
To install the dependent R packages, please refer to the manual of each individual package.

## Raw data availability
Raw sequencing data of the profiled samples have been deposited in the Genome Sequence Archive (GSA, https://ngdc.cncb.ac.cn/gsa) at the National Genomics Data Center, China National Center for Bioinformation / Beiling Institute of Genomics, Chinese Academy of Sciences, under accession numbers: HRA012868, HRA012869, HRA012870 and HRA012903. Mass spectrometry proteomics, phosphoproteomics and DNA methylation data have been deposited in OMIX (https://ngdc.cncb.ac.cn/omix) at the same institute under accession numbers: OMIX011526, OMIX011528 and OMIX011529, respectively. Spatial omics data downloading and visualizations of regions of interest are available at https://wang-lab.hkust.edu.hk/software/STP. Data from TCGA were downloaded from NCI Genomics Data Commons data portal (https://portal.gdc.cancer.gov). Previously published CGGA data have been uploaded to GSA, under BioProject ID: PRJCA001636 and PRJCA001747. All the other data supporting the findings of this study are available from the lead contact upon reasonable request.

## Processed data availability
Processed level 3 data could be downloaded from https://www.cgga.org.cn/download.jsp under the DataSet ID CCell_4083.  

## Contact
For any questions, please contact Professor Jiguang Wang via email: jgwang AT ust DOT hk
