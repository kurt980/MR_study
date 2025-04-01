# MR Study
Collecting code for analysis in the ABCD study

# Citations
If you consider using code from this Github repo please cite:

# Data

## ABCD Study
Data comes from the ABCD Study (Volkow et al., 2018): (https://abcdstudy.org/). The ABCD Study is the largest longitudinal study of the developmental trajectories of the brain and cognitive, social, and emotional growth of children in the United States. 

Among the 11,099 ABCD participants aged 9-10, we include the 6,381 (3,172 female) samples that have neuroimaging, genetic, and behavioral data collected at baseline in our study.

ABCD study collects data from multiple sites. The Data used for the 6,381 samples in the study contains 22 sites.
Neuroimaging data contains numeric SCFC values for 379 brain ROIs for 6,381 samples acquired from MRI; 
Genetic data contains 1.8 million SNPs per individual, stored as PLINK binary files (ped, map, bed, bim, fam); 
Behavioral data contains ADHD scores for 6,381 samples measured using CBCL standard, and binary ADHD diagnosis (0,1) for 6,381 samples.

Demographic data contains values for sex, race, ethnicity and age.

# Code

**association_scripts** code for running correlation, linear regression and other association analysis.
**gwas_scripts** code for running GWAS on SFC.
**heritability_scripts** code for running GREML for heritability estimates.
**onesample_mr_scripts** code for running one-sample MR (2SLS and 2SRI)
**twoample_mr_scripts** code for running two-sample MR (IVW, MR-CUE, MR-Egger)
**genetic_analysis_scripts**