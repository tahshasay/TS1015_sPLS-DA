# TS1015_sPLS-DA

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
