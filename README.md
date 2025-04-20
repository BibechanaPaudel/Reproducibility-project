# Reproducibility-project

**Title:** Impact of Plant Growth Promoting Rhizobacteria (PGPR) Treatments on Potato Virus Y (PVY) Accumulation in *Nicotiana benthamiana*

This repository contains the full R-based analytical workflow assessing the effect of PGPR on the accumulation of PVY in *Nicotiana benthamiana* across three independent biological replicates.

### Project Overview

**Objective:**  
To evaluate whether PGPR treatments mitigate PVY accumulation, based on qPCR-derived Cq values, analyzed over time (Dpi) in both inoculated and systemic leaf tissues.

**Treatments:**
- Control  
- *Pseudomonas fluorescens*  
- *Serratia marcescens*  
- *Bacillus amyloliquefaciens*  
- *Bacillus subtilis*

### Data Summary

- **Inputs:** Cq values from RT-qPCR for PVY detection  
- **Factors:** Treatment, Days post-inoculation (Dpi), Replication (technical), and Tissue type  
- **Format:** Three CSV files (one per biological replicate), located in the `Data/` folder  
- **Preprocessing:** Calculation of viral load using the regression equation from Feng et al. (2006), followed by log10 transformation

### Analysis Highlights

- All three biological replicates are analyzed within a **single RMarkdown (`.Rmd`) file**  
- Data is merged and labeled by replicate for downstream statistical modeling  
- Linear model assess treatment and time interactions  
- Tukey-adjusted comparisons are performed using `emmeans` and `multcompView`  
- Visualizations include bar plots for inoculated and systemic leaves with error bars and significance letters  
- Outputs include combined plots with common legends for effective comparison

```
├── Nb_PGPR+PVY.html  
├── Nb_PGPR+PVY.md
├── Nb_PGPR+PVY.Rmd
├── Nb_PGPR+PVY_1st_Reproducibility.csv  #1st trial data
├── Nb_PGPR+PVY_2nd_Reproducibility.csv  #2nd trial data
├── Nb_PGPR+PVY_3rd_Reproducibility.csv  #3rd trial data
├── Nb_PGPR+PVY_files
│   └── figure-gfm
│       ├── Combine fig 1st trial-1.png
│       ├── Combine fig 2nd trial-1.png
│       ├── Combine fig 3rd trial-1.png
│       ├── PVY on IL_1st trial_NB-1.png
│       ├── PVY on IL_2nd trial_NB-1.png
│       ├── PVY on IL_3rd trial_NB-1.png
│       ├── PVY on SL_1st trial_NB-1.png
│       ├── PVY on SL_2ndtrial_NB-1.png
│       └── PVY on SL_3rdtrial_NB-1.png
├── README.md
└── Reproducibility-project.Rproj
```