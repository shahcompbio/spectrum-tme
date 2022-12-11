# Ovarian cancer mutational processes drive site-specific immune evasion

This repository contains supporting material for the manuscript:

> I. Vázquez-García, F. Uhlitz et al, ["Ovarian cancer mutational processes drive site-specific immune evasion"](https://www.nature.com/articles/s41586-022-05496-1), *Nature* (2022)

## Table of Contents

- [Code](#code)
- [Data](#data)
  - [Overview](#overview)
  - [Data availability](#data-availability)
    - [scRNA-seq](#scrna-seq)
    - [H&E](#hne)
    - [mpIF](#mpif)
    - [Bulk WGS](#bulk-wgs)
    - [MSK-IMPACT](#msk-impact)
              
## Code

The source code contains RMarkdown notebooks to reproduce the manuscript figures and tables.

## Data

### Overview

The MSK SPECTRUM study is registered on dbGaP under accession number [phs002857.v1.p1](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002857.v1.p1).

A Synapse page for the MSK SPECTRUM study is available to provide access to multi-modal datasets from one central location. This page can be accessed under accession number [syn25569736](https://www.synapse.org/msk_spectrum).

### Data availability

#### scRNA-seq
  - Expression counts are available from [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE180661).
  - Processed objects are available from [Synapse](https://www.synapse.org/#!Synapse:syn33521743/datasets/).
  - Interactive data visualizations are available via [cellxgene](https://cellxgene.cziscience.com/collections/4796c91c-9d8f-4692-be43-347b1727f9d8).

#### Tumor-normal bulk WGS
  - Raw sequencing reads are available for controlled access from the NCBI Sequence Read Archive via [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gapprev/gap/cgi-bin/study.cgi?study_id=phs002857.v1.p1).
  - Somatic mutations and copy number data can be accessed from [Synapse](https://www.synapse.org/#!Synapse:syn33521770/datasets/).
  - Somatic mutations, copy number and structural variants can be visualized through [cBioPortal](https://cbioportal.org/study/summary?id=msk_spectrum_tme_2022).

#### Tumor-normal targeted panel sequencing (MSK-IMPACT)
  - Somatic mutations, copy number and structural variants can be visualized through [cBioPortal](https://cbioportal.org/study/summary?id=msk_spectrum_tme_2022).

#### H&E
  - Deidentified images are available via [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002857.v1.p1).
  - Cell segmentation, tissue segmentation and lymphocyte classification are available from [Synapse](https://www.synapse.org/#!Synapse:syn33521762/datasets/).

#### mpIF
  - Deidentified images are available via [dbGaP](https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs002857.v1.p1).
  - Cell segmentation, tissue segmentation and cell phenotyping are available from [Synapse](https://www.synapse.org/#!Synapse:syn33520881/datasets/).