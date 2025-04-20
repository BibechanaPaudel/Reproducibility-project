# Impact of PGPR on PVY Accumulation in *Nicotiana benthamiana*

This repository contains the full R-based analytical workflow assessing the effect of PGPR on the accumulation of PVY in *Nicotiana benthamiana* across three independent biological replicates.The complete dataset used in the project is included in the repository with the folder name data to ensure full transparency and reproducibility. 

**Note:** PGPR treatments were applied one week prior to PVY inoculation, allowing time for induced resistance responses to develop.

### Folder and File Structure
- [Data](https://github.com/BibechanaPaudel/Reproducibility-project/tree/main/Nb_PGPR%2BPVY_files/figure-gfm): Contains .csv files for each of the three biological replicates.

- [Nb_PGPR+PVY.Rmd](https://github.com/BibechanaPaudel/Reproducibility-project/blob/main/Nb_PGPR%2BPVY.Rmd): The complete R script with code, output, and narrative.

- [Nb_PGPR+PVY.html](https://github.com/BibechanaPaudel/Reproducibility-project/blob/main/Nb_PGPR%2BPVY.html): Rendered HTML version of the RMarkdown.

- [Nb_PGPR+PVY.md](https://github.com/BibechanaPaudel/Reproducibility-project/blob/main/Nb_PGPR%2BPVY.md): Markdown-exported version for GitHub-friendly viewing.

- [Nb_PGPR+PVY_files/figure-gfm](https://github.com/BibechanaPaudel/Reproducibility-project/tree/main/Nb_PGPR%2BPVY_files/figure-gfm): Output plots in high-quality vector format (.png). 


### Project Overview

**Objective:**  
To evaluate whether PGPR treatments mitigate PVY accumulation, based on qPCR-derived Cq values, analyzed over time (Dpi) in both inoculated and systemic leaf tissues.

**Treatments:**
- Control  
- *Pseudomonas fluorescens*  
- *Serratia marcescens*  
- *Bacillus amyloliquefaciens*  
- *Bacillus subtilis*

### Script workflow

#### **1. Data Import and Preparation**

- Three biological replicate datasets (`.csv`) are stored in the `Data/` folder.
- All data are processed in a single RMarkdown (`.Rmd`) script.
- Each file contains Cq values from RT-qPCR for PVY quantification, along with metadata (Treatment, Dpi, Replicate).
- Viral load is calculated from Cq values using the regression equation from [Feng, J.L et al., 2006](https://academic.oup.com/abbs/article/38/10/669/217), then log-transformed for normality.

#### **2. Data Grouping and Statistical Analysis**

- Data are grouped by Treatment, Days post-inoculation (Dpi), and Replicate.
- A linear model is fit to the log-transformed viral load to analyze interaction effects.
- Estimated marginal means are calculated using `emmeans`, followed by Tukey-adjusted pairwise comparisons.
- Significance groupings (letters) are extracted for visualization.

#### **3. Visualisation**

- Bar plots are generated for:
  - **Inoculated leaves** at 1, 4, 7, and 10 Dpi
  - **Systemic leaves** at 7 and 10 Dpi
- Visuals include:
  - Grouped bar charts with error bars (standard error)
  - Statistical significance labels over bars
  - Combined plots with shared legends
- Color-blind friendly palettes are used.

#### **4. Statistical Testing**

- ANOVA and type-II tests (via `car::Anova`) assess Treatment*Dpi interactions.
- Post hoc multiple comparisons are performed using `emmeans` and `multcompView`.
- Model fit is evaluated using residuals and diagnostics.

#### **5. Helper Functions (within the RMarkdown)**

- **Viral Load Transformation:** Applies regression formula for converting Cq to copy number.
- **calculateMeansAndSE:** Aggregates replicate data to compute average and standard error.
- **addSigLetters:** Merges statistical letters for use in plots.
- **LeafTypeClassifier:** Categorizes sample type (Inoculated vs. Systemic) based on Dpi.
- **PlotTemplates:** Modular ggplot2 code to streamline visual consistency across figures.
- **CombinePlots:** Uses `ggpubr::ggarrange()` to create unified multi-panel figures.

```
├── Data
│   ├── Nb_PGPR+PVY_1st_Reproducibility.csv   ## Data from first biological replication
│   ├── Nb_PGPR+PVY_2nd_Reproducibility.csv   ## Data from second biological replication
│   └── Nb_PGPR+PVY_3rd_Reproducibility.csv   ## Data from third biological replication
├── Nb_PGPR+PVY.html   ##script in html format
├── Nb_PGPR+PVY.md     ##script in github favoured markdowm format
├── Nb_PGPR+PVY.Rmd    ##script in rmd format
├── Nb_PGPR+PVY_files
│   └── figure-gfm
│       ├── Combine fig 1st trial-1.png  ## combine fig from 1st biological replication
│       ├── Combine fig 2nd trial-1.png  ## combine fig from 2nd biological replication
│       ├── Combine fig 3rd trial-1.png  ## combine fig from 3rd biological replication
│       ├── PVY on IL_1st trial_NB-1.png 
│       ├── PVY on IL_2nd trial_NB-1.png
│       ├── PVY on IL_3rd trial_NB-1.png
│       ├── PVY on SL_1st trial_NB-1.png
│       ├── PVY on SL_2ndtrial_NB-1.png
│       └── PVY on SL_3rdtrial_NB-1.png
├── README.md 
└── Reproducibility-project.Rproj
```