# TS1015_sPLS-DA

#These data and scripts are supplement to Say, T. E., & Degnan, S. M. (2019). Interdependent photo- and chemosensory systems regulate larval settlement in a marine sponge. bioRxiv, 519512. doi:10.1101/519512

#To identify the genes that regulate larval settlement, based on light-dependent and age-related changes in gene expression, we used the multivariate sparse partial least squares discriminant analysis (sPLS-DA) (Lê Cao, Boitard, & Besse, 2011), implemented in the mixOmics package (Lê Cao, 2016) in R version 3.3.1 (R Core Team, 2016); see additional supplementary data on GitHub at https://github.com/tahshasay/TS1015_sPLS-DA. 

#This folder contains the files and scripts used to perform sPLS-DA on the CEL-Seq2 data.

Script:
'file name:  R_Script_sPLS-DA_TS1015.sh' - R commands used to perform sPLS-DA

This script is divided into different parts, outlined below:
Part i)	pls-da
Part ii)	sPLS-DA tuning
Part iii)	final sPLS-DA analysis
Part iv)	sPLS-DA plot
Part v)	heatmap


OUTPUT_Files:
"TS1015_Design.txt" - Experimental design file which includes the light regime and time (hours post emergence) for each sample. Light regimes include include natural day-night cycle (nat), constant dark (drk) and constant light (lgt). Time (hours post emergence) include 1 hpe (t1) and 9 hpe (t9). 

"TS1015_group.txt" - Experimental design file which includes the light regime and time (hours post emergence) for each sample. Light regimes include include natural day-night cycle (nat), constant dark (drk) and constant light (lgt). Time (hours post emergence) include 1 hpe (t1) and 9 hpe (t9).

"TS1015 plot of the classification error rate 2017.07.20.pdf" - plot of the error classification rate.

'TS1015_perf.plsda_error.rate_2017.07.20.txt' - list of error classification rates.

'TS1015_ptune.splsda.TS1015$choice.keepX_2017.07.20.txt' - optimal number of variables to select on the two components.

References:
Lê Cao, K.-A., Boitard, S., & Besse, P. (2011). Sparse PLS discriminant analysis: biologically relevant feature selection and graphical displays for multiclass problems. BMC Bioinformatics, 12(1), 253. doi:10.1186/1471-2105-12-253

Lê Cao, K.-A., Rohart F., Gautier B., Bartolo F., Gonz’lez I., Déjean S.,. (2016). mixOmics: omics data integration project. R package version 6.1.1. Retrieved from https://CRAN.R-project.org/package=mixOmics

R Core Team. (2016). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. Retrieved from https://www.R-project.org/

