# Project Title

This repo contains the code and data for https://retinomics.org, an interactive website, that was developed to help visualize the GWAS and metabolomics results in the paper. 

The website has interactive pages for each type of biomarker of interest, such as GWAS, metabolites, hematology, infections, polygenic risk scores and PheCodes. For each type of biomarker, there are two views. First, there are heatmap plots to show the effect and pvalue of selected biomarker of interest. Second, there are macula location plots that shows the link between all significant biomarkers for that macula location. The location plots are manhattan/region plots for the GWAS and volcano and bar plots for the other types of biomarkers. Users can interact with the plots by clicking or selecting a biomarker of interest to trigger the loading of its heatmaps or click a macula location in the heatmaps to load its corresponding location plots. The website was implemented in HTML, native Javascript, and d3js.

## Table of Contents
- html folder contains the html and JS code that renders the interactive visualisations
- html/data contains the various data files by biomarkers of interest, such as GWAS, immunomics, metabolomics, phenomics etc.
- html/data/GWAS is organized by Chromosomes, chr1 to chr22 and chrX, which contain the spatial effects and p-values for every SNP on that chromosome that was significant at least one spatial location
- html/data/GWAS also contains location_summary, which contains information about SNPs that were signicant at the spatial location for each of the approximately 30,000 spatial locations across the macula.
- html/data/GWAS/manhattan_imgs contain pre-processed Manhattan plots for each spatial location
- most files are compressed to minimize bandwidth requirements
- Total data is approximately 100GB.  

## Installation

No installation required.