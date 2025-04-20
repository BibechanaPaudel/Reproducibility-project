# **First trial**

\###Libraries used

``` r
#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("agricolae")
#install.packages("tidyr")
#install.packages("ggpubr")
#install.packages("multcomp")
#install.packages("emmeans")
#install.packages("multcompView")
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.2     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(ggplot2)
library(agricolae)
library(tidyr)
library(ggpubr)
library(multcomp)
```

    ## Loading required package: mvtnorm
    ## Loading required package: survival
    ## Loading required package: TH.data
    ## Loading required package: MASS
    ## 
    ## Attaching package: 'MASS'
    ## 
    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select
    ## 
    ## 
    ## Attaching package: 'TH.data'
    ## 
    ## The following object is masked from 'package:MASS':
    ## 
    ##     geyser

``` r
library(emmeans)
```

    ## Welcome to emmeans.
    ## Caution: You lose important information if you filter this package's results.
    ## See '? untidy'

``` r
library(multcompView)
```

\###Load and verify data of Rep 1

``` r
rep1 <- read.csv("Data/Nb_PGPR+PVY_1st_Reproducibility.csv", header=T)
head(rep1)  # View the first few rows
```

    ##   Treatment  Dpi Replication    Cq
    ## 1   Control IL-1           1 19.72
    ## 2   Control IL-1           2 19.56
    ## 3   Control IL-1           3 19.98
    ## 4   Control IL-4           1 17.55
    ## 5   Control IL-4           2 17.48
    ## 6   Control IL-4           3 17.58

``` r
str(rep1)   # Check the structure of the data
```

    ## 'data.frame':    90 obs. of  4 variables:
    ##  $ Treatment  : chr  "Control" "Control" "Control" "Control" ...
    ##  $ Dpi        : chr  "IL-1" "IL-1" "IL-1" "IL-4" ...
    ##  $ Replication: int  1 2 3 1 2 3 1 2 3 1 ...
    ##  $ Cq         : num  19.7 19.6 20 17.6 17.5 ...

``` r
summary(rep1) # Quick statistical summary
```

    ##   Treatment             Dpi             Replication       Cq       
    ##  Length:90          Length:90          Min.   :1    Min.   :15.98  
    ##  Class :character   Class :character   1st Qu.:1    1st Qu.:20.60  
    ##  Mode  :character   Mode  :character   Median :2    Median :21.72  
    ##                                        Mean   :2    Mean   :23.33  
    ##                                        3rd Qu.:3    3rd Qu.:26.77  
    ##                                        Max.   :3    Max.   :33.18

\###Make the factor variable

``` r
rep1$Treatment<-as.factor(rep1$Treatment)
rep1$Dpi<-as.factor(rep1$Dpi)
```

\##Calculate the ViralLoad using the regression equation derived by
(Feng , J.L et al., 2006 ) and log transform for normality.

[Feng, J.L et al.,
2006](https://academic.oup.com/abbs/article/38/10/669/217)

``` r
rep1<-rep1 %>% 
mutate (ViralLoad=((10^((Cq-49.153)/-3.93)))) %>% ##added a new column and named ViralLoad
mutate(logViralLoad=log10(ViralLoad))  ##added a new column and named logViralLoad
```

## Statistical analysis

``` r
model<-lm(logViralLoad~Treatment*Dpi,data=rep1) #linear model of interaction
car::Anova(model) #anova for linear model
```

    ## Anova Table (Type II tests)
    ## 
    ## Response: logViralLoad
    ##                Sum Sq Df  F value    Pr(>F)    
    ## Treatment      13.409  4  2050.38 < 2.2e-16 ***
    ## Dpi           104.420  5 12773.16 < 2.2e-16 ***
    ## Treatment:Dpi  14.435 20   441.45 < 2.2e-16 ***
    ## Residuals       0.098 60                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
lsmeans <- emmeans(model, ~Treatment|Dpi) # estimate lsmeans of Dpi within Treatment
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE,  Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC
```

    ## $emmeans
    ## Dpi = IL-1:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                7.48 0.0233 60     7.43     7.53  a    
    ##  B. subtilis            7.36 0.0233 60     7.31     7.41   b   
    ##  B. amyloliquifaciens   7.27 0.0233 60     7.22     7.31   b   
    ##  P. fluorescens         7.15 0.0233 60     7.11     7.20    c  
    ##  S. marcescens          6.97 0.0233 60     6.93     7.02     d 
    ## 
    ## Dpi = IL-10:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                8.43 0.0233 60     8.38     8.47  a    
    ##  S. marcescens          8.22 0.0233 60     8.17     8.26   b   
    ##  B. amyloliquifaciens   7.21 0.0233 60     7.16     7.26    c  
    ##  B. subtilis            7.00 0.0233 60     6.95     7.04     d 
    ##  P. fluorescens         6.66 0.0233 60     6.61     6.71      e
    ## 
    ## Dpi = IL-4:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                8.04 0.0233 60     8.00     8.09  a    
    ##  S. marcescens          7.08 0.0233 60     7.03     7.12   b   
    ##  B. amyloliquifaciens   6.87 0.0233 60     6.83     6.92    c  
    ##  B. subtilis            6.79 0.0233 60     6.74     6.84    c  
    ##  P. fluorescens         6.64 0.0233 60     6.59     6.68     d 
    ## 
    ## Dpi = IL-7:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                8.37 0.0233 60     8.33     8.42  a    
    ##  S. marcescens          7.61 0.0233 60     7.56     7.66   b   
    ##  P. fluorescens         7.15 0.0233 60     7.10     7.20    c  
    ##  B. amyloliquifaciens   7.08 0.0233 60     7.04     7.13    cd 
    ##  B. subtilis            7.01 0.0233 60     6.96     7.06     d 
    ## 
    ## Dpi = SL-10:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                5.70 0.0233 60     5.65     5.74  a    
    ##  P. fluorescens         5.25 0.0233 60     5.20     5.29   b   
    ##  B. subtilis            4.85 0.0233 60     4.80     4.90    c  
    ##  S. marcescens          4.30 0.0233 60     4.25     4.34     d 
    ##  B. amyloliquifaciens   4.15 0.0233 60     4.10     4.19      e
    ## 
    ## Dpi = SL-7:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  B. amyloliquifaciens   6.39 0.0233 60     6.35     6.44  a    
    ##  Control                5.87 0.0233 60     5.82     5.91   b   
    ##  S. marcescens          4.80 0.0233 60     4.75     4.84    c  
    ##  P. fluorescens         4.71 0.0233 60     4.66     4.76    cd 
    ##  B. subtilis            4.69 0.0233 60     4.64     4.73     d 
    ## 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 5 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## Dpi = IL-1:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  P. fluorescens - S. marcescens          0.1790 0.033 60   5.421  <.0001
    ##  B. amyloliquifaciens - S. marcescens    0.2952 0.033 60   8.940  <.0001
    ##  B. amyloliquifaciens - P. fluorescens   0.1162 0.033 60   3.520  0.0072
    ##  B. subtilis - S. marcescens             0.3859 0.033 60  11.689  <.0001
    ##  B. subtilis - P. fluorescens            0.2070 0.033 60   6.269  <.0001
    ##  B. subtilis - B. amyloliquifaciens      0.0908 0.033 60   2.749  0.0585
    ##  Control - S. marcescens                 0.5081 0.033 60  15.389  <.0001
    ##  Control - P. fluorescens                0.3291 0.033 60   9.968  <.0001
    ##  Control - B. amyloliquifaciens          0.2129 0.033 60   6.448  <.0001
    ##  Control - B. subtilis                   0.1221 0.033 60   3.699  0.0042
    ## 
    ## Dpi = IL-10:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  B. subtilis - P. fluorescens            0.3350 0.033 60  10.148  <.0001
    ##  B. amyloliquifaciens - P. fluorescens   0.5479 0.033 60  16.596  <.0001
    ##  B. amyloliquifaciens - B. subtilis      0.2129 0.033 60   6.448  <.0001
    ##  S. marcescens - P. fluorescens          1.5564 0.033 60  47.142  <.0001
    ##  S. marcescens - B. subtilis             1.2214 0.033 60  36.994  <.0001
    ##  S. marcescens - B. amyloliquifaciens    1.0085 0.033 60  30.546  <.0001
    ##  Control - P. fluorescens                1.7659 0.033 60  53.488  <.0001
    ##  Control - B. subtilis                   1.4309 0.033 60  43.340  <.0001
    ##  Control - B. amyloliquifaciens          1.2180 0.033 60  36.892  <.0001
    ##  Control - S. marcescens                 0.2095 0.033 60   6.346  <.0001
    ## 
    ## Dpi = IL-4:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  B. subtilis - P. fluorescens            0.1535 0.033 60   4.650  0.0002
    ##  B. amyloliquifaciens - P. fluorescens   0.2349 0.033 60   7.116  <.0001
    ##  B. amyloliquifaciens - B. subtilis      0.0814 0.033 60   2.466  0.1122
    ##  S. marcescens - P. fluorescens          0.4377 0.033 60  13.256  <.0001
    ##  S. marcescens - B. subtilis             0.2841 0.033 60   8.606  <.0001
    ##  S. marcescens - B. amyloliquifaciens    0.2027 0.033 60   6.140  <.0001
    ##  Control - P. fluorescens                1.4071 0.033 60  42.621  <.0001
    ##  Control - B. subtilis                   1.2536 0.033 60  37.971  <.0001
    ##  Control - B. amyloliquifaciens          1.1722 0.033 60  35.504  <.0001
    ##  Control - S. marcescens                 0.9695 0.033 60  29.364  <.0001
    ## 
    ## Dpi = IL-7:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  B. amyloliquifaciens - B. subtilis      0.0729 0.033 60   2.209  0.1903
    ##  P. fluorescens - B. subtilis            0.1399 0.033 60   4.239  0.0007
    ##  P. fluorescens - B. amyloliquifaciens   0.0670 0.033 60   2.030  0.2647
    ##  S. marcescens - B. subtilis             0.5988 0.033 60  18.138  <.0001
    ##  S. marcescens - B. amyloliquifaciens    0.5259 0.033 60  15.928  <.0001
    ##  S. marcescens - P. fluorescens          0.4589 0.033 60  13.899  <.0001
    ##  Control - B. subtilis                   1.3622 0.033 60  41.259  <.0001
    ##  Control - B. amyloliquifaciens          1.2892 0.033 60  39.050  <.0001
    ##  Control - P. fluorescens                1.2222 0.033 60  37.020  <.0001
    ##  Control - S. marcescens                 0.7634 0.033 60  23.122  <.0001
    ## 
    ## Dpi = SL-10:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  S. marcescens - B. amyloliquifaciens    0.1493 0.033 60   4.522  0.0003
    ##  B. subtilis - B. amyloliquifaciens      0.7040 0.033 60  21.323  <.0001
    ##  B. subtilis - S. marcescens             0.5547 0.033 60  16.802  <.0001
    ##  P. fluorescens - B. amyloliquifaciens   1.1009 0.033 60  33.346  <.0001
    ##  P. fluorescens - S. marcescens          0.9517 0.033 60  28.825  <.0001
    ##  P. fluorescens - B. subtilis            0.3969 0.033 60  12.023  <.0001
    ##  Control - B. amyloliquifaciens          1.5496 0.033 60  46.937  <.0001
    ##  Control - S. marcescens                 1.4003 0.033 60  42.415  <.0001
    ##  Control - B. subtilis                   0.8456 0.033 60  25.614  <.0001
    ##  Control - P. fluorescens                0.4487 0.033 60  13.590  <.0001
    ## 
    ## Dpi = SL-7:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  P. fluorescens - B. subtilis            0.0221 0.033 60   0.668  0.9625
    ##  S. marcescens - B. subtilis             0.1103 0.033 60   3.340  0.0121
    ##  S. marcescens - P. fluorescens          0.0882 0.033 60   2.672  0.0703
    ##  Control - B. subtilis                   1.1798 0.033 60  35.736  <.0001
    ##  Control - P. fluorescens                1.1578 0.033 60  35.068  <.0001
    ##  Control - S. marcescens                 1.0696 0.033 60  32.396  <.0001
    ##  B. amyloliquifaciens - B. subtilis      1.7048 0.033 60  51.638  <.0001
    ##  B. amyloliquifaciens - P. fluorescens   1.6828 0.033 60  50.970  <.0001
    ##  B. amyloliquifaciens - S. marcescens    1.5946 0.033 60  48.298  <.0001
    ##  B. amyloliquifaciens - Control          0.5250 0.033 60  15.902  <.0001
    ## 
    ## P value adjustment: tukey method for comparing a family of 5 estimates

\##Extract significance letters for plotting

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(Results_lsmeansEC$emmeans$Treatment, 
                               Results_lsmeansEC$emmeans$Dpi,
                               str_trim(Results_lsmeansEC$emmeans$.group))
colnames(sig.diff.letters) <- c("Treatment", 
                                "Dpi",
                                "Letters")
```

\##Data summarization

``` r
PVY1 <- rep1 %>%
  group_by(Treatment, Dpi) %>%   
  dplyr::summarize(
    avg.virus= mean(logViralLoad, na.rm=TRUE),
    n=n(),
    se = sd(logViralLoad)/sqrt(n)) %>%    # calculates the mean and se 
  left_join(sig.diff.letters) # join the significant letters with the PVY1 data frame.
```

    ## `summarise()` has grouped output by 'Treatment'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Treatment, Dpi)`

\##Classify the sample collection days based on leaf

``` r
PVY1$Leaf_sample <- ifelse(
  grepl("^IL", PVY1$Dpi), "Inoculated",  # Matches Dpi values starting with IL
  ifelse(grepl("^SL", PVY1$Dpi), "Systemic", NA)
)
```

## Visualization: Inoculated sample

``` r
inoculated<-filter(PVY1, Dpi=="IL-1" | Dpi=="IL-4" | Dpi=="IL-7" | Dpi=="IL-10")
inoculated$Treatment<-factor(inoculated$Treatment, levels=c("Control","P. fluorescens", "S. marcescens", "B. amyloliquifaciens", "B. subtilis"))
```

``` r
cbbPalette <- c( "#56B4E9","#F0E442","#CC79A7", "#E69F00", "#009E73")  ##Load colour blind friendly cbbPalette
PVY_IL_1<-ggplot(inoculated, aes(x =factor(Dpi,levels=c("IL-1", "IL-4", "IL-7", "IL-10")), y =avg.virus , fill = Treatment)) +
  stat_summary(fun=mean,geom="bar", position = position_dodge(width = 0.7), width = 0.6) +
   geom_errorbar(aes(ymin = avg.virus - se, 
                    ymax = avg.virus + se),
                width = 0.4,
                position = position_dodge(width = 0.7)) +
  geom_text(data = inoculated, aes(label = Letters, y = avg.virus),position = position_dodge(width = 0.7), vjust = -0.5)+
  theme_classic() +  #add the extracted significance letter to the bar
  labs(title= "" ,
       x = "",
       y = "Log Copy Number") +  # labels of x and y axis
  scale_fill_manual(values=cbbPalette,
                    labels=c("Control",expression(italic("P. fluorescens")), 
                             expression(italic("S. marcescens")),expression(italic("B. amyloliquifaciens")),expression(italic("B. subtilis"))))+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)), breaks = seq(1, 9, by = 1)  # This removes padding at the bottom of the plot
  )+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),  # Adjusts the size of the legend keys
        legend.title = element_blank(),  # Adjusts the legend title font size
        legend.text = element_text(size = 8),
               axis.text.x = element_text(color = "black", face="bold"),   # X-axis text color
        axis.text.y = element_text(color = "black", face="bold"),   # Y-axis text color
        axis.title.x = element_text(color = "black"),  # X-axis label color
        axis.title.y = element_text(color = "black", face="bold", size=11),  # Y-axis label color
        axis.line = element_line(color = "black")      # Axis lines color
  )

 PVY_IL_1 
```

![](Nb_PGPR+PVY_files/figure-gfm/PVY%20on%20IL_1st%20trial_NB-1.png)<!-- -->

## Visualization: Systemic sample

``` r
systemic<-filter(PVY1, Dpi=="SL-7" | Dpi=="SL-10")
systemic$Treatment<-factor(systemic$Treatment, levels=c("Control","P. fluorescens", "S. marcescens", "B. amyloliquifaciens", "B. subtilis"))
```

``` r
cbbPalette <- c( "#56B4E9","#F0E442","#CC79A7", "#E69F00", "#009E73")
PVY_SL_1<-ggplot(systemic, aes(x =factor(Dpi,levels=c("SL-7","SL-10")), y =avg.virus , fill = Treatment)) +
  stat_summary(fun=mean,geom="bar", position = position_dodge(width = 0.7), width = 0.6) +
   geom_errorbar(aes(ymin = avg.virus - se, 
                    ymax = avg.virus + se),
                width = 0.4,
                position = position_dodge(width = 0.7)) +
  geom_text(data = systemic, aes(label = Letters, y = avg.virus),position = position_dodge(width = 0.7), vjust = -0.5)+
  theme_classic() +
  labs(title= "" ,
       x = "",
       y = "Log Copy Number") +
  scale_fill_manual(values=cbbPalette,
                    labels=c("Control",expression(italic("P. fluorescens")), 
                             expression(italic("S. marcescens")),expression(italic("B. amyloliquifaciens")),expression(italic("B. subtilis"))))+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)), breaks = seq(1, 9, by = 1)  # This removes padding at the bottom of the plot
  )+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),  # Adjusts the size of the legend keys
        legend.title = element_blank(),  # Adjusts the legend title font size
        legend.text = element_text(size = 8),
                axis.text.x = element_text(color = "black", face="bold"),   # X-axis text color
        axis.text.y = element_text(color = "black", face="bold"),   # Y-axis text color
        axis.title.x = element_text(color = "black"),  # X-axis label color
        axis.title.y = element_text(color = "black", face="bold", size=11),  # Y-axis label color
        axis.line = element_line(color = "black")      # Axis lines color
  )

 PVY_SL_1 
```

![](Nb_PGPR+PVY_files/figure-gfm/PVY%20on%20SL_1st%20trial_NB-1.png)<!-- -->

## Combined figure

``` r
PVY.1st<-ggarrange(PVY_IL_1, PVY_SL_1,    #combine the  figures 
                   nrow=1,
                   ncol=2,
                   common.legend=T,    # Gives the common legend for both figures
                   widths=c(1.75,1))   # Makes the width of 1st fig 1.75 times more than 2nd fig
PVY.1st
```

![](Nb_PGPR+PVY_files/figure-gfm/Combine%20fig%201st%20trial-1.png)<!-- -->

# **Second trial**

\##Load and verify the data

``` r
rep2 <- read.csv("Data/Nb_PGPR+PVY_2nd_Reproducibility.csv", na.strings="na")
str(rep2)   # Check the structure of the data
```

    ## 'data.frame':    90 obs. of  4 variables:
    ##  $ Treatment  : chr  "Control" "Control" "Control" "Control" ...
    ##  $ Dpi        : chr  "IL-1" "IL-1" "IL-1" "IL-4" ...
    ##  $ Replication: int  1 2 3 1 2 3 1 2 3 1 ...
    ##  $ Cq         : num  20.9 21 21.2 24 24 ...

``` r
summary(rep2) # Quick statistical summary
```

    ##   Treatment             Dpi             Replication       Cq       
    ##  Length:90          Length:90          Min.   :1    Min.   :18.92  
    ##  Class :character   Class :character   1st Qu.:1    1st Qu.:22.48  
    ##  Mode  :character   Mode  :character   Median :2    Median :23.14  
    ##                                        Mean   :2    Mean   :24.94  
    ##                                        3rd Qu.:3    3rd Qu.:28.22  
    ##                                        Max.   :3    Max.   :32.46

\###Make the factor variable

``` r
rep2$Treatment<-as.factor(rep2$Treatment)
rep2$Dpi<-as.factor(rep2$Dpi)
```

\##Calculate the ViralLoad using the regression equation derived by
(Feng , J.L et al., 2006 ) and log transform for normality.

``` r
rep2<-rep2 %>% 
mutate (ViralLoad=((10^((Cq-49.153)/-3.93)))) %>% 
mutate(logViralLoad=log10(ViralLoad))
```

## Statistical analysis

``` r
model<-lm(logViralLoad~Treatment*Dpi,data=rep2)
car::Anova(model)
```

    ## Anova Table (Type II tests)
    ## 
    ## Response: logViralLoad
    ##               Sum Sq Df  F value    Pr(>F)    
    ## Treatment      2.953  4   649.19 < 2.2e-16 ***
    ## Dpi           71.942  5 12651.37 < 2.2e-16 ***
    ## Treatment:Dpi  6.385 20   280.72 < 2.2e-16 ***
    ## Residuals      0.068 60                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
lsmeans <- emmeans(model, ~Treatment|Dpi) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE,  Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC
```

    ## $emmeans
    ## Dpi = IL-1:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control               7.154 0.0195 60    7.115    7.193  a    
    ##  B. amyloliquifaciens  6.917 0.0195 60    6.878    6.956   b   
    ##  S. marcescens         6.840 0.0195 60    6.802    6.879   bc  
    ##  B. subtilis           6.773 0.0195 60    6.735    6.812    c  
    ##  P. fluorescens        6.622 0.0195 60    6.584    6.661     d 
    ## 
    ## Dpi = IL-10:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control               6.997 0.0195 60    6.958    7.036  a    
    ##  B. subtilis           6.783 0.0195 60    6.744    6.822   b   
    ##  B. amyloliquifaciens  6.760 0.0195 60    6.721    6.799   b   
    ##  P. fluorescens        6.710 0.0195 60    6.671    6.749   b   
    ##  S. marcescens         6.613 0.0195 60    6.574    6.652    c  
    ## 
    ## Dpi = IL-4:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  B. subtilis           6.887 0.0195 60    6.848    6.926  a    
    ##  B. amyloliquifaciens  6.762 0.0195 60    6.723    6.801   b   
    ##  S. marcescens         6.718 0.0195 60    6.679    6.757   b   
    ##  P. fluorescens        6.486 0.0195 60    6.447    6.525    c  
    ##  Control               6.400 0.0195 60    6.361    6.439     d 
    ## 
    ## Dpi = IL-7:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control               7.675 0.0195 60    7.636    7.714  a    
    ##  P. fluorescens        7.009 0.0195 60    6.970    7.048   b   
    ##  B. subtilis           6.777 0.0195 60    6.738    6.816    c  
    ##  B. amyloliquifaciens  6.557 0.0195 60    6.518    6.596     d 
    ##  S. marcescens         6.350 0.0195 60    6.311    6.389      e
    ## 
    ## Dpi = SL-10:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control               5.453 0.0195 60    5.414    5.492  a    
    ##  B. amyloliquifaciens  5.226 0.0195 60    5.187    5.265   b   
    ##  B. subtilis           4.671 0.0195 60    4.632    4.710    c  
    ##  S. marcescens         4.436 0.0195 60    4.397    4.475     d 
    ##  P. fluorescens        4.369 0.0195 60    4.330    4.408     d 
    ## 
    ## Dpi = SL-7:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  P. fluorescens        5.528 0.0195 60    5.489    5.567  a    
    ##  Control               5.323 0.0195 60    5.284    5.362   b   
    ##  S. marcescens         4.787 0.0195 60    4.748    4.826    c  
    ##  B. subtilis           4.655 0.0195 60    4.616    4.694     d 
    ##  B. amyloliquifaciens  4.565 0.0195 60    4.526    4.604      e
    ## 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 5 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## Dpi = IL-1:
    ##  contrast                              estimate     SE df t.ratio p.value
    ##  B. subtilis - P. fluorescens            0.1510 0.0275 60   5.483  <.0001
    ##  S. marcescens - P. fluorescens          0.2180 0.0275 60   7.916  <.0001
    ##  S. marcescens - B. subtilis             0.0670 0.0275 60   2.433  0.1205
    ##  B. amyloliquifaciens - P. fluorescens   0.2943 0.0275 60  10.689  <.0001
    ##  B. amyloliquifaciens - B. subtilis      0.1433 0.0275 60   5.206  <.0001
    ##  B. amyloliquifaciens - S. marcescens    0.0763 0.0275 60   2.772  0.0552
    ##  Control - P. fluorescens                0.5318 0.0275 60  19.314  <.0001
    ##  Control - B. subtilis                   0.3808 0.0275 60  13.831  <.0001
    ##  Control - S. marcescens                 0.3138 0.0275 60  11.397  <.0001
    ##  Control - B. amyloliquifaciens          0.2375 0.0275 60   8.625  <.0001
    ## 
    ## Dpi = IL-10:
    ##  contrast                              estimate     SE df t.ratio p.value
    ##  P. fluorescens - S. marcescens          0.0967 0.0275 60   3.512  0.0073
    ##  B. amyloliquifaciens - S. marcescens    0.1467 0.0275 60   5.329  <.0001
    ##  B. amyloliquifaciens - P. fluorescens   0.0500 0.0275 60   1.817  0.3734
    ##  B. subtilis - S. marcescens             0.1696 0.0275 60   6.161  <.0001
    ##  B. subtilis - P. fluorescens            0.0729 0.0275 60   2.649  0.0742
    ##  B. subtilis - B. amyloliquifaciens      0.0229 0.0275 60   0.832  0.9197
    ##  Control - S. marcescens                 0.3842 0.0275 60  13.954  <.0001
    ##  Control - P. fluorescens                0.2875 0.0275 60  10.442  <.0001
    ##  Control - B. amyloliquifaciens          0.2375 0.0275 60   8.625  <.0001
    ##  Control - B. subtilis                   0.2146 0.0275 60   7.793  <.0001
    ## 
    ## Dpi = IL-4:
    ##  contrast                              estimate     SE df t.ratio p.value
    ##  P. fluorescens - Control                0.0857 0.0275 60   3.111  0.0230
    ##  S. marcescens - Control                 0.3181 0.0275 60  11.551  <.0001
    ##  S. marcescens - P. fluorescens          0.2324 0.0275 60   8.440  <.0001
    ##  B. amyloliquifaciens - Control          0.3613 0.0275 60  13.122  <.0001
    ##  B. amyloliquifaciens - P. fluorescens   0.2757 0.0275 60  10.011  <.0001
    ##  B. amyloliquifaciens - S. marcescens    0.0433 0.0275 60   1.571  0.5215
    ##  B. subtilis - Control                   0.4869 0.0275 60  17.681  <.0001
    ##  B. subtilis - P. fluorescens            0.4012 0.0275 60  14.570  <.0001
    ##  B. subtilis - S. marcescens             0.1688 0.0275 60   6.130  <.0001
    ##  B. subtilis - B. amyloliquifaciens      0.1255 0.0275 60   4.559  0.0002
    ## 
    ## Dpi = IL-7:
    ##  contrast                              estimate     SE df t.ratio p.value
    ##  B. amyloliquifaciens - S. marcescens    0.2070 0.0275 60   7.516  <.0001
    ##  B. subtilis - S. marcescens             0.4266 0.0275 60  15.494  <.0001
    ##  B. subtilis - B. amyloliquifaciens      0.2197 0.0275 60   7.978  <.0001
    ##  P. fluorescens - S. marcescens          0.6590 0.0275 60  23.934  <.0001
    ##  P. fluorescens - B. amyloliquifaciens   0.4521 0.0275 60  16.418  <.0001
    ##  P. fluorescens - B. subtilis            0.2324 0.0275 60   8.440  <.0001
    ##  Control - S. marcescens                 1.3249 0.0275 60  48.114  <.0001
    ##  Control - B. amyloliquifaciens          1.1179 0.0275 60  40.598  <.0001
    ##  Control - B. subtilis                   0.8982 0.0275 60  32.620  <.0001
    ##  Control - P. fluorescens                0.6658 0.0275 60  24.180  <.0001
    ## 
    ## Dpi = SL-10:
    ##  contrast                              estimate     SE df t.ratio p.value
    ##  S. marcescens - P. fluorescens          0.0670 0.0275 60   2.433  0.1205
    ##  B. subtilis - P. fluorescens            0.3020 0.0275 60  10.966  <.0001
    ##  B. subtilis - S. marcescens             0.2349 0.0275 60   8.532  <.0001
    ##  B. amyloliquifaciens - P. fluorescens   0.8575 0.0275 60  31.142  <.0001
    ##  B. amyloliquifaciens - S. marcescens    0.7905 0.0275 60  28.708  <.0001
    ##  B. amyloliquifaciens - B. subtilis      0.5556 0.0275 60  20.176  <.0001
    ##  Control - P. fluorescens                1.0840 0.0275 60  39.366  <.0001
    ##  Control - S. marcescens                 1.0170 0.0275 60  36.933  <.0001
    ##  Control - B. subtilis                   0.7820 0.0275 60  28.400  <.0001
    ##  Control - B. amyloliquifaciens          0.2265 0.0275 60   8.224  <.0001
    ## 
    ## Dpi = SL-7:
    ##  contrast                              estimate     SE df t.ratio p.value
    ##  B. subtilis - B. amyloliquifaciens      0.0899 0.0275 60   3.265  0.0150
    ##  S. marcescens - B. amyloliquifaciens    0.2222 0.0275 60   8.070  <.0001
    ##  S. marcescens - B. subtilis             0.1323 0.0275 60   4.805  0.0001
    ##  Control - B. amyloliquifaciens          0.7583 0.0275 60  27.538  <.0001
    ##  Control - B. subtilis                   0.6684 0.0275 60  24.273  <.0001
    ##  Control - S. marcescens                 0.5360 0.0275 60  19.468  <.0001
    ##  P. fluorescens - B. amyloliquifaciens   0.9635 0.0275 60  34.992  <.0001
    ##  P. fluorescens - B. subtilis            0.8736 0.0275 60  31.727  <.0001
    ##  P. fluorescens - S. marcescens          0.7413 0.0275 60  26.922  <.0001
    ##  P. fluorescens - Control                0.2053 0.0275 60   7.454  <.0001
    ## 
    ## P value adjustment: tukey method for comparing a family of 5 estimates

\##Extract significance letters for plotting

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(Results_lsmeansEC$emmeans$Treatment, 
                               Results_lsmeansEC$emmeans$Dpi,
                               str_trim(Results_lsmeansEC$emmeans$.group))
colnames(sig.diff.letters) <- c("Treatment", 
                                "Dpi",
                                "Letters")
```

\##Data summarization

``` r
PVY2 <- rep2 %>%
  group_by(Treatment, Dpi) %>%
  dplyr::summarize(
    avg.virus= mean(logViralLoad, na.rm=TRUE),
    n=n(),
    se = sd(logViralLoad)/sqrt(n)) %>%
  left_join(sig.diff.letters)
```

    ## `summarise()` has grouped output by 'Treatment'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Treatment, Dpi)`

\##Classify the sample collection days based on leaf

``` r
PVY2$Leaf_sample <- ifelse(
  grepl("^IL", PVY2$Dpi), "Inoculated",  # Matches Dpi values starting with IL
  ifelse(grepl("^SL", PVY2$Dpi), "Systemic", NA)
)
```

## Visualization: Inoculated sample

``` r
inoculated<-filter(PVY2, Dpi=="IL-1" | Dpi=="IL-4" | Dpi=="IL-7" | Dpi=="IL-10")
inoculated$Treatment<-factor(inoculated$Treatment, levels=c("Control","P. fluorescens", "S. marcescens", "B. amyloliquifaciens", "B. subtilis"))
```

``` r
cbbPalette <- c( "#56B4E9","#F0E442","#CC79A7", "#E69F00", "#009E73")
PVY_IL_2<-ggplot(inoculated, aes(x =factor(Dpi,levels=c("IL-1", "IL-4", "IL-7", "IL-10")), y =avg.virus , fill = Treatment)) +
  stat_summary(fun=mean,geom="bar", position = position_dodge(width = 0.7), width = 0.6) +
   geom_errorbar(aes(ymin = avg.virus - se, 
                    ymax = avg.virus + se),
                width = 0.4,
                position = position_dodge(width = 0.7)) +
  geom_text(data = inoculated, aes(label = Letters, y = avg.virus),position = position_dodge(width = 0.7), vjust = -0.5)+
  theme_classic() +
  labs(title= "" ,
       x = "",
       y = "Log Copy Number") +
  scale_fill_manual(values=cbbPalette,
                    labels=c("Control",expression(italic("P. fluorescens")), 
                             expression(italic("S. marcescens")),expression(italic("B. amyloliquifaciens")),expression(italic("B. subtilis"))))+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)), breaks = seq(1, 9, by = 1)  # This removes padding at the bottom of the plot
  )+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),  # Adjusts the size of the legend keys
        legend.title = element_blank(),  # Adjusts the legend title font size
        legend.text = element_text(size = 8),
               axis.text.x = element_text(color = "black", face="bold"),   # X-axis text color and font
        axis.text.y = element_text(color = "black", face="bold"),   # Y-axis text color and font
        axis.title.x = element_text(color = "black"),  # X-axis label color
        axis.title.y = element_text(color = "black", face="bold", size=11),  # Y-axis label color, size and font
        axis.line = element_line(color = "black")      # Axis lines color
  )

 PVY_IL_2 
```

![](Nb_PGPR+PVY_files/figure-gfm/PVY%20on%20IL_2nd%20trial_NB-1.png)<!-- -->

## Visualization: Systemic sample

``` r
systemic<-filter(PVY2, Dpi=="SL-7" | Dpi=="SL-10")
systemic$Treatment<-factor(systemic$Treatment, levels=c("Control","P. fluorescens", "S. marcescens", "B. amyloliquifaciens", "B. subtilis"))
```

``` r
cbbPalette <- c( "#56B4E9","#F0E442","#CC79A7", "#E69F00", "#009E73")
PVY_SL_2<-ggplot(systemic, aes(x =factor(Dpi,levels=c("SL-7","SL-10")), y =avg.virus , fill = Treatment)) +
  stat_summary(fun=mean,geom="bar", position = position_dodge(width = 0.7), width = 0.6) +
   geom_errorbar(aes(ymin = avg.virus - se, 
                    ymax = avg.virus + se),
                width = 0.4,
                position = position_dodge(width = 0.7)) +
  geom_text(data = systemic, aes(label = Letters, y = avg.virus),position = position_dodge(width = 0.7), vjust = -0.5)+
  theme_classic() +
  labs(title= "" ,
       x = "",
       y = "Log Copy Number") +
  scale_fill_manual(values=cbbPalette,
                    labels=c("Control",expression(italic("P. fluorescens")), 
                             expression(italic("S. marcescens")),expression(italic("B. amyloliquifaciens")),expression(italic("B. subtilis"))))+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)), breaks = seq(1, 9, by = 1)  # This removes padding at the bottom of the plot
  )+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),  # Adjusts the size of the legend keys
        legend.title = element_blank(),  # Adjusts the legend title font size
        legend.text = element_text(size = 8),
                axis.text.x = element_text(color = "black", face="bold"),   # X-axis text color
        axis.text.y = element_text(color = "black", face="bold"),   # Y-axis text color
        axis.title.x = element_text(color = "black"),  # X-axis label color
        axis.title.y = element_text(color = "black", face="bold", size=11),  # Y-axis label color
        axis.line = element_line(color = "black")      # Axis lines color
  )

 PVY_SL_2 
```

![](Nb_PGPR+PVY_files/figure-gfm/PVY%20on%20SL_2ndtrial_NB-1.png)<!-- -->

## Combined figure

``` r
PVY.2nd<-ggarrange(PVY_IL_2, PVY_SL_2,    #combine figures
                   nrow=1,
                   ncol=2,
                   common.legend=T,
                   widths=c(1.75,1))   #make  width of 1st fig 1.75 times more than 2nd fig
PVY.2nd
```

![](Nb_PGPR+PVY_files/figure-gfm/Combine%20fig%202nd%20trial-1.png)<!-- -->

# **Third trial**

``` r
## Load data Rep 3
rep3 <- read.csv("Data/Nb_PGPR+PVY_3rd_Reproducibility.csv", na.strings="na")
####Verify the data
head(rep3)  # View the first few rows
```

    ##   Treatment  Dpi Replication    Cq
    ## 1   Control IL-1           1 22.09
    ## 2   Control IL-1           2 21.93
    ## 3   Control IL-1           3 21.78
    ## 4   Control IL-4           1 22.41
    ## 5   Control IL-4           2 22.17
    ## 6   Control IL-4           3 22.06

``` r
str(rep3)   # Check the structure of the data
```

    ## 'data.frame':    90 obs. of  4 variables:
    ##  $ Treatment  : chr  "Control" "Control" "Control" "Control" ...
    ##  $ Dpi        : chr  "IL-1" "IL-1" "IL-1" "IL-4" ...
    ##  $ Replication: int  1 2 3 1 2 3 1 2 3 1 ...
    ##  $ Cq         : num  22.1 21.9 21.8 22.4 22.2 ...

``` r
summary(rep3) # Quick statistical summary
```

    ##   Treatment             Dpi             Replication       Cq       
    ##  Length:90          Length:90          Min.   :1    Min.   :21.18  
    ##  Class :character   Class :character   1st Qu.:1    1st Qu.:22.43  
    ##  Mode  :character   Mode  :character   Median :2    Median :23.43  
    ##                                        Mean   :2    Mean   :25.25  
    ##                                        3rd Qu.:3    3rd Qu.:29.61  
    ##                                        Max.   :3    Max.   :32.35

\###Make the factor variable

``` r
rep3$Treatment<-as.factor(rep3$Treatment)
rep3$Dpi<-as.factor(rep3$Dpi)
```

\##Calculate the ViralLoad using the regression equation derived by
(Feng , J.L et al., 2006 ) and log transform for normality.

``` r
rep3<-rep3 %>% 
mutate (ViralLoad=((10^((Cq-49.153)/-3.93)))) %>% 
mutate(logViralLoad=log10(ViralLoad))
```

## Statistical analysis

``` r
model<-lm(logViralLoad~Treatment*Dpi,data=rep3)
car::Anova(model)
```

    ## Anova Table (Type II tests)
    ## 
    ## Response: logViralLoad
    ##               Sum Sq Df  F value    Pr(>F)    
    ## Treatment      2.230  4   386.37 < 2.2e-16 ***
    ## Dpi           74.915  5 10384.94 < 2.2e-16 ***
    ## Treatment:Dpi  6.571 20   227.73 < 2.2e-16 ***
    ## Residuals      0.087 60                       
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
lsmeans <- emmeans(model, ~Treatment|Dpi) # estimate lsmeans of variety within siteXyear
Results_lsmeansEC <- multcomp::cld(lsmeans, alpha = 0.05, reversed = TRUE, details = TRUE,  Letters = letters) # contrast with Tukey ajustment
Results_lsmeansEC
```

    ## $emmeans
    ## Dpi = IL-1:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                6.93 0.0219 60     6.88     6.97  a    
    ##  B. subtilis            6.82 0.0219 60     6.77     6.86   b   
    ##  S. marcescens          6.71 0.0219 60     6.67     6.76    c  
    ##  B. amyloliquifaciens   6.66 0.0219 60     6.61     6.70    cd 
    ##  P. fluorescens         6.58 0.0219 60     6.54     6.62     d 
    ## 
    ## Dpi = IL-10:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                7.11 0.0219 60     7.07     7.16  a    
    ##  B. amyloliquifaciens   6.74 0.0219 60     6.70     6.79   b   
    ##  S. marcescens          6.71 0.0219 60     6.66     6.75   b   
    ##  B. subtilis            6.44 0.0219 60     6.40     6.48    c  
    ##  P. fluorescens         6.36 0.0219 60     6.32     6.41    c  
    ## 
    ## Dpi = IL-4:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  B. subtilis            7.06 0.0219 60     7.01     7.10  a    
    ##  B. amyloliquifaciens   7.04 0.0219 60     6.99     7.08  a    
    ##  P. fluorescens         6.90 0.0219 60     6.86     6.95   b   
    ##  Control                6.85 0.0219 60     6.81     6.90   b   
    ##  S. marcescens          6.67 0.0219 60     6.63     6.71    c  
    ## 
    ## Dpi = IL-7:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  B. subtilis            7.04 0.0219 60     6.99     7.08  a    
    ##  Control                6.59 0.0219 60     6.54     6.63   b   
    ##  P. fluorescens         6.51 0.0219 60     6.46     6.55   bc  
    ##  S. marcescens          6.47 0.0219 60     6.43     6.52    c  
    ##  B. amyloliquifaciens   6.26 0.0219 60     6.21     6.30     d 
    ## 
    ## Dpi = SL-10:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  Control                5.90 0.0219 60     5.86     5.94  a    
    ##  S. marcescens          4.58 0.0219 60     4.54     4.62   b   
    ##  B. subtilis            4.41 0.0219 60     4.36     4.45    c  
    ##  B. amyloliquifaciens   4.37 0.0219 60     4.33     4.41    c  
    ##  P. fluorescens         4.34 0.0219 60     4.30     4.39    c  
    ## 
    ## Dpi = SL-7:
    ##  Treatment            emmean     SE df lower.CL upper.CL .group
    ##  P. fluorescens         5.23 0.0219 60     5.19     5.27  a    
    ##  B. subtilis            4.94 0.0219 60     4.90     4.99   b   
    ##  S. marcescens          4.94 0.0219 60     4.89     4.98   b   
    ##  Control                4.84 0.0219 60     4.80     4.88    c  
    ##  B. amyloliquifaciens   4.46 0.0219 60     4.41     4.50     d 
    ## 
    ## Confidence level used: 0.95 
    ## P value adjustment: tukey method for comparing a family of 5 estimates 
    ## significance level used: alpha = 0.05 
    ## NOTE: If two or more means share the same grouping symbol,
    ##       then we cannot show them to be different.
    ##       But we also did not show them to be the same. 
    ## 
    ## $comparisons
    ## Dpi = IL-1:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  B. amyloliquifaciens - P. fluorescens  0.07718 0.031 60   2.489  0.1069
    ##  S. marcescens - P. fluorescens         0.13232 0.031 60   4.266  0.0007
    ##  S. marcescens - B. amyloliquifaciens   0.05513 0.031 60   1.778  0.3959
    ##  B. subtilis - P. fluorescens           0.23494 0.031 60   7.576  <.0001
    ##  B. subtilis - B. amyloliquifaciens     0.15776 0.031 60   5.087  <.0001
    ##  B. subtilis - S. marcescens            0.10263 0.031 60   3.309  0.0132
    ##  Control - P. fluorescens               0.34521 0.031 60  11.131  <.0001
    ##  Control - B. amyloliquifaciens         0.26802 0.031 60   8.642  <.0001
    ##  Control - S. marcescens                0.21289 0.031 60   6.864  <.0001
    ##  Control - B. subtilis                  0.11026 0.031 60   3.555  0.0064
    ## 
    ## Dpi = IL-10:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  B. subtilis - P. fluorescens           0.07803 0.031 60   2.516  0.1006
    ##  S. marcescens - P. fluorescens         0.34266 0.031 60  11.049  <.0001
    ##  S. marcescens - B. subtilis            0.26463 0.031 60   8.533  <.0001
    ##  B. amyloliquifaciens - P. fluorescens  0.38083 0.031 60  12.280  <.0001
    ##  B. amyloliquifaciens - B. subtilis     0.30280 0.031 60   9.763  <.0001
    ##  B. amyloliquifaciens - S. marcescens   0.03817 0.031 60   1.231  0.7336
    ##  Control - P. fluorescens               0.74894 0.031 60  24.149  <.0001
    ##  Control - B. subtilis                  0.67091 0.031 60  21.633  <.0001
    ##  Control - S. marcescens                0.40628 0.031 60  13.100  <.0001
    ##  Control - B. amyloliquifaciens         0.36811 0.031 60  11.869  <.0001
    ## 
    ## Dpi = IL-4:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  Control - S. marcescens                0.18575 0.031 60   5.989  <.0001
    ##  P. fluorescens - S. marcescens         0.23494 0.031 60   7.576  <.0001
    ##  P. fluorescens - Control               0.04919 0.031 60   1.586  0.5119
    ##  B. amyloliquifaciens - S. marcescens   0.36726 0.031 60  11.842  <.0001
    ##  B. amyloliquifaciens - Control         0.18151 0.031 60   5.853  <.0001
    ##  B. amyloliquifaciens - P. fluorescens  0.13232 0.031 60   4.266  0.0007
    ##  B. subtilis - S. marcescens            0.38677 0.031 60  12.471  <.0001
    ##  B. subtilis - Control                  0.20102 0.031 60   6.482  <.0001
    ##  B. subtilis - P. fluorescens           0.15182 0.031 60   4.895  0.0001
    ##  B. subtilis - B. amyloliquifaciens     0.01951 0.031 60   0.629  0.9698
    ## 
    ## Dpi = IL-7:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  S. marcescens - B. amyloliquifaciens   0.21459 0.031 60   6.919  <.0001
    ##  P. fluorescens - B. amyloliquifaciens  0.24936 0.031 60   8.040  <.0001
    ##  P. fluorescens - S. marcescens         0.03478 0.031 60   1.121  0.7947
    ##  Control - B. amyloliquifaciens         0.32994 0.031 60  10.639  <.0001
    ##  Control - S. marcescens                0.11535 0.031 60   3.719  0.0039
    ##  Control - P. fluorescens               0.08058 0.031 60   2.598  0.0835
    ##  B. subtilis - B. amyloliquifaciens     0.78032 0.031 60  25.161  <.0001
    ##  B. subtilis - S. marcescens            0.56573 0.031 60  18.241  <.0001
    ##  B. subtilis - P. fluorescens           0.53096 0.031 60  17.120  <.0001
    ##  B. subtilis - Control                  0.45038 0.031 60  14.522  <.0001
    ## 
    ## Dpi = SL-10:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  B. amyloliquifaciens - P. fluorescens  0.02714 0.031 60   0.875  0.9049
    ##  B. subtilis - P. fluorescens           0.06361 0.031 60   2.051  0.2549
    ##  B. subtilis - B. amyloliquifaciens     0.03647 0.031 60   1.176  0.7649
    ##  S. marcescens - P. fluorescens         0.23749 0.031 60   7.658  <.0001
    ##  S. marcescens - B. amyloliquifaciens   0.21035 0.031 60   6.782  <.0001
    ##  S. marcescens - B. subtilis            0.17388 0.031 60   5.606  <.0001
    ##  Control - P. fluorescens               1.55725 0.031 60  50.212  <.0001
    ##  Control - B. amyloliquifaciens         1.53011 0.031 60  49.337  <.0001
    ##  Control - B. subtilis                  1.49364 0.031 60  48.161  <.0001
    ##  Control - S. marcescens                1.31976 0.031 60  42.554  <.0001
    ## 
    ## Dpi = SL-7:
    ##  contrast                              estimate    SE df t.ratio p.value
    ##  Control - B. amyloliquifaciens         0.38168 0.031 60  12.307  <.0001
    ##  S. marcescens - B. amyloliquifaciens   0.47668 0.031 60  15.370  <.0001
    ##  S. marcescens - Control                0.09500 0.031 60   3.063  0.0262
    ##  B. subtilis - B. amyloliquifaciens     0.48261 0.031 60  15.561  <.0001
    ##  B. subtilis - Control                  0.10093 0.031 60   3.254  0.0155
    ##  B. subtilis - S. marcescens            0.00594 0.031 60   0.191  0.9997
    ##  P. fluorescens - B. amyloliquifaciens  0.77099 0.031 60  24.860  <.0001
    ##  P. fluorescens - Control               0.38931 0.031 60  12.553  <.0001
    ##  P. fluorescens - S. marcescens         0.29432 0.031 60   9.490  <.0001
    ##  P. fluorescens - B. subtilis           0.28838 0.031 60   9.299  <.0001
    ## 
    ## P value adjustment: tukey method for comparing a family of 5 estimates

\##Extract significance letters for plotting

``` r
# Extracting the letters for the bars
sig.diff.letters <- data.frame(Results_lsmeansEC$emmeans$Treatment, 
                               Results_lsmeansEC$emmeans$Dpi,
                               str_trim(Results_lsmeansEC$emmeans$.group))
colnames(sig.diff.letters) <- c("Treatment", 
                                "Dpi",
                                "Letters")
```

\##Data summarization

``` r
PVY3 <- rep3 %>%
  group_by(Treatment, Dpi) %>%
  dplyr::summarize(
    avg.virus= mean(logViralLoad, na.rm=TRUE),
    n=n(),
    se = sd(logViralLoad)/sqrt(n)) %>%
  left_join(sig.diff.letters)
```

    ## `summarise()` has grouped output by 'Treatment'. You can override using the
    ## `.groups` argument.
    ## Joining with `by = join_by(Treatment, Dpi)`

\##Classify the sample collection days based on leaf

``` r
PVY3$Leaf_sample <- ifelse(
  grepl("^IL", PVY3$Dpi), "Inoculated",  # Matches Dpi values starting with IL
  ifelse(grepl("^SL", PVY3$Dpi), "Systemic", NA)
)
```

## Visualization: Inoculated sample

``` r
inoculated<-filter(PVY3, Dpi=="IL-1" | Dpi=="IL-4" | Dpi=="IL-7" | Dpi=="IL-10")
inoculated$Treatment<-factor(inoculated$Treatment, levels=c("Control","P. fluorescens", "S. marcescens", "B. amyloliquifaciens", "B. subtilis"))
```

``` r
cbbPalette <- c( "#56B4E9","#F0E442","#CC79A7", "#E69F00", "#009E73")
PVY_IL_3<-ggplot(inoculated, aes(x =factor(Dpi,levels=c("IL-1", "IL-4", "IL-7", "IL-10")), y =avg.virus , fill = Treatment)) +
  stat_summary(fun=mean,geom="bar", position = position_dodge(width = 0.7), width = 0.6) +
   geom_errorbar(aes(ymin = avg.virus - se, 
                    ymax = avg.virus + se),
                width = 0.4,
                position = position_dodge(width = 0.7)) +
  geom_text(data = inoculated, aes(label = Letters, y = avg.virus),position = position_dodge(width = 0.7), vjust = -0.5)+
  theme_classic() +
  labs(title= "" ,
       x = "",
       y = "Log Copy Number") +
  scale_fill_manual(values=cbbPalette,
                    labels=c("Control",expression(italic("P. fluorescens")), 
                             expression(italic("S. marcescens")),expression(italic("B. amyloliquifaciens")),expression(italic("B. subtilis"))))+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)), breaks = seq(1, 9, by = 1)  # This removes padding at the bottom of the plot
  )+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),  # Adjusts the size of the legend keys
        legend.title = element_blank(),  # Adjusts the legend title font size
        legend.text = element_text(size = 8),
               axis.text.x = element_text(color = "black", face="bold"),   # X-axis text color
        axis.text.y = element_text(color = "black", face="bold"),   # Y-axis text color
        axis.title.x = element_text(color = "black"),  # X-axis label color
        axis.title.y = element_text(color = "black", face="bold", size=11),  # Y-axis label color
        axis.line = element_line(color = "black")      # Axis lines color
  )

 PVY_IL_3
```

![](Nb_PGPR+PVY_files/figure-gfm/PVY%20on%20IL_3rd%20trial_NB-1.png)<!-- -->

## Visualization: Systemic sample

``` r
systemic<-filter(PVY3, Dpi=="SL-7" | Dpi=="SL-10")
systemic$Treatment<-factor(systemic$Treatment, levels=c("Control","P. fluorescens", "S. marcescens", "B. amyloliquifaciens", "B. subtilis"))
```

``` r
cbbPalette <- c( "#56B4E9","#F0E442","#CC79A7", "#E69F00", "#009E73")
PVY_SL_3<-ggplot(systemic, aes(x =factor(Dpi,levels=c("SL-7","SL-10")), y =avg.virus , fill = Treatment)) +
  stat_summary(fun=mean,geom="bar", position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = avg.virus - se, 
                    ymax = avg.virus + se),
                width = 0.4,
                position = position_dodge(width = 0.7)) +
  geom_text(data = systemic, aes(label = Letters, y = avg.virus+(3*se)),position = position_dodge(width = 0.7), vjust = -0.5)+
  theme_classic() +
  labs(title= "" ,
       x = "",
       y = "Log Copy Number") +
  scale_fill_manual(values=cbbPalette,
                    labels=c("Control",expression(italic("P. fluorescens")), 
                             expression(italic("S. marcescens")),expression(italic("B. amyloliquifaciens")),expression(italic("B. subtilis"))))+
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.1)), breaks = seq(1, 9, by = 1)  # This removes padding at the bottom of the plot
  )+
  theme(legend.position = "right",
        legend.key.size = unit(0.5, "cm"),  # Adjusts the size of the legend keys
        legend.title = element_blank(),  # Adjusts the legend title font size
        legend.text = element_text(size = 8),
                axis.text.x = element_text(color = "black", face="bold"),   # X-axis text color
        axis.text.y = element_text(color = "black", face="bold"),   # Y-axis text color
        axis.title.x = element_text(color = "black"),  # X-axis label color
        axis.title.y = element_text(color = "black", face="bold", size=11),  # Y-axis label color
        axis.line = element_line(color = "black")      # Axis lines color
  )

 PVY_SL_3
```

![](Nb_PGPR+PVY_files/figure-gfm/PVY%20on%20SL_3rdtrial_NB-1.png)<!-- -->

## Combined figure

``` r
PVY.3rd<-ggarrange(PVY_IL_3, PVY_SL_3,        #combine figure  
                   nrow=1,
                   ncol=2,
                   common.legend=T,
                   widths=c(1.75,1))       #Makes 1st fig 1.75 times wider than the 2nd fig. 
PVY.3rd
```

![](Nb_PGPR+PVY_files/figure-gfm/Combine%20fig%203rd%20trial-1.png)<!-- -->
