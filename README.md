# Metabolic-Elasticity-Score-Project

Metabolic Elasticity Score (MElaS) quantifies each metabolite's ability to adapt or recover after a disturbance such as an infection over time. We modified a formula derived from Zhou et al. that quantified gene elasticity at a transcript level. Here, we made a new method that provides a score for each metabolite using relative abundance from LC/MS. We then implemented this new tool to study the effects of Chagas disease over time at a metabolic level. To compute elasticity for each metabolite, natural log fold changes and p-values between two consecutive endpoints are required. MElaS_Compute_FC.Rmd will produce these FC and p-values and MElaS.Rmd will compute MElaS values. This code has an example of 3, 12 and 24 weeks because we want to calculate MElaS between 3 and 12 weeks and MElaS values between 12 and 24 weeks.

## MElaS_Compute_FC.Rmd:
- Reads DisTol_mergeddata.csv
    - feature table that includes metadata and feature table
    - Raw MS data processed through mzmine
    - TIC normalized data
- Median peak area values for each metabolite of all infected and uninfected samples at each 3, 12 and 24 weeks for each organ are computed
    - A value of 5.06 x 10⁻⁷, the minimum value in the dataset, is added to each median value
- Fold change between median infected over median uninfected at each endpoint is then calculated
    - an initial endpoint of 0 weeks is required so FC at 0 weeks is set to a value of 1
- Natural log fold changes between consecutive times are then determined for each metabolite
- P-values using Mann-Whitney U test to compare the differences of infected individuals between 3 versus 12 weeks and 12 versus 24 weeks for each metabolite are computed

These values will be required to calculate MElaS
lnFC_3v0
lnFC_12v0
lnFC_24v3
lnFC_24v12
adj.P.Val_24v3
adj.P.Val_24v12


## Compute_MElaS.Rmd
- R notebook that uses the code from (Zhou, et.al (2022)) with modifications
- Calculate Elastic Score between 3 and 12 weeks and between 12 and 24 weeks
- produce a csv file called MElaS_Values.csv containing MElaS values, FC at 0w, 3w, 12w, 24w lnFC at 3v0, 12v3 and 24v12 for all organs and conditions for analysis
  
### Paper and GitHub referenced in this code:
Zhou, Q., Yu, L., Cook, J.R., Qiang, L. & Sun, L. 2023, Cell Metabolism, 35, 1661-167.e6, doi: 10.1016/j.cmet.2023.08.001.
https://github.com/zhouqz/GElaS

## MElaSAnalysis.Rmd:
- Reads in MElaS_Values.csv
- Produce temporal patterns of each metabolite
- Chemical annotation of all metabolites and most and least elastic metabolites
- Impact of inoculum  and T. cruzi strain on MElaS scores by Benjamini–Hochberg FDR corrected Wilcoxon signed-rank test
    - Boxplots of log-scaled MElaS values by inoculum and strain
- Impact of inoculation concentration and T. cruzi strain on the proportion of elastic metabolites (number of metabolites that are elastic) by Benjamini–Hochberg FDR corrected Fisher’s Exact test
- Differences in MElaS values by sampling site by Benjamini–Hochberg FDR corrected Dunn’s post hoc multiple comparison test
    - Boxplots of MElaS by organ