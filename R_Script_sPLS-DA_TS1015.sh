
####################################################################################
#These data and scripts are supplement to Say, T. E., & Degnan, S. M. (2019). 
#Interdependent photo- and chemosensory systems regulate larval settlement in a marine sponge. bioRxiv, 519512. doi:10.1101/519512

#To identify the genes that regulate larval settlement, based on light-dependent and age-related changes in gene expression, 
#we used the multivariate sparse partial least squares discriminant analysis (sPLS-DA) (Lê Cao, Boitard, & Besse, 2011), 
#implemented in the mixOmics package (Lê Cao, 2016) in R version 3.3.1 (R Core Team, 2016); see our customised script below.
#
#References: 
#Lê Cao, K.-A., Boitard, S., & Besse, P. (2011). Sparse PLS discriminant analysis: biologically relevant feature selection and graphical displays for multiclass problems. BMC Bioinformatics, 12(1), 253. doi:10.1186/1471-2105-12-253
#
#Lê Cao, K.-A., Rohart F., Gautier B., Bartolo F., Gonz’lez I., Déjean S.,. (2016). mixOmics: omics data integration project. R package version 6.1.1. Retrieved from https://CRAN.R-project.org/package=mixOmics
#
#R Core Team. (2016). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. Retrieved from https://www.R-project.org/
#
#
####################################################################################


#cd /opsin/u/tahsha/TSay_R_wd/splsda_wd # to change

# load required v of R 
# R v 3.3.1 has been used for TS1015 and SSCellTypes

#---------------
# mixOmics: v3.3.1
#---------------

#module load R/3.3.1-foss-2016a

#library('mixOmics', lib.loc="/opsin/u/tahsha/R/x86_64-redhat-linux-gnu-library/3.3") 
# above no longer works - new v R installed on Piwi and mixOmics is not compatible

#when the systems were updated, the default R on Piwi and Ago moved to 3.4.1. Using the module command you can load another v.
# So I asked nick to install v 3.3.1 on opsin on 19 Sept, 2017

library('mixOmics')

# call a nifty self contained script with: 
# nohup Rscript[command] <my R script> &

# The nohup command runs a command in the background and the ampersand is unix for return the prompt. This will redirect all standard output to a file nohup.out file placed in the directory you called it. 

# check my running processes
#top -u tahsha

# check for error messages
#less nohup.txt
# + should see  R  program running if it is working

# check library location
#.Library
#[1] "/usr/lib64/R/library"



####################################################################################
# R script for 'splsda' analysis using the mixomics package 
# Code by 
# http://mixomics.org/case-studies/spls-da-srbct-2/ 
# 
# and mixomics workshop material that attended 13 and 14 Aug, 2015. 
# Presented by Kim-Anh Le Cao
# File:
# /Users/tahshasay/Documents/PCFB_practical_computing_for_biologists/R workshop course material/MixOmics_Kim-Anh Le Cao_2015/mixOmics workshop TRI_updatedbyCOURSE_2015.08.17/R scripts/cs_multiblock.R
#
# Adapted by TSay for TS1015 in 2017
#
####################################################################################


setwd("/opsin/u/tahsha/TSay_R_wd/splsda_wd")
getwd()

# move files to piwi
# read in sample file
samples=read.table("TS1015_Design.txt", header=TRUE,row.names=1)
# or
# group file created above:
group=read.table("TS1015_group.txt", header=TRUE,row.names=1)


### to change
TS1015=read.table("/opsin/u/tahsha/TSay_R_wd/TS1015_transformed_data/TS1015_DESeq2_vsdMat_Blind_2016.11.21.txt", header=TRUE,row.names=1)
#header=TRUE indicates that the first line contains column names and row.names=1 means that the first column should be used as row names.
head(TS1015) # print the first 7 rows of the file
tail(TS1015)
dim(TS1015)  # print the dimensions of the tabl


# ----------------------------------------------------------
# load required packages
#library('mixOmics', lib.loc="/opsin/u/tahsha/R/x86_64-redhat-linux-gnu-library/3.3")
library('mixOmics')
# ----------------------------------------------------------

 #includes ggplot
#library("dplyr")tra

# samples <- samples[-28,] ##remove a dodgy sample
X <- TS1015
X <- t(X) #transpose the data


#Y <- samples$treatment
Y <- group$x
summary(Y)


#----#---- Part a) of every splsda script - setting up wd etc
#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----#----

# One practical way to choose the number of components is to run a PLS-DA model first with a large number of sPLS-DA components (e.g. ncomp = 10 ) using repeated cross-validation (here folds = 5 and nrepeat = 10 ), then use the function perf() which evaluates the performance of the model per component. This step will allow to choose the distance that minimises the classification error rate and the number of optimal components. Our experience has shown that usually ncomp = K-1 where K is the number of classes is often optimal, but this is highly dependent on the data.

#this chunk takes ~ 5 min to run

############
set.seed(32) # for reproducibility of the outputs of this code that performs random cross-validation sampling. To be removed in pro

############
TS1015.plsda.perf <- plsda(X, Y, ncomp = 10)

# to speed up computation in this example we can choose 5 folds repeated 10 times:
# however, 10 would probably be best so I will run 10 on piwi (faster than my computer). 
############
perf.plsda <- perf(TS1015.plsda.perf, validation = 'Mfold', folds = 6, progressBar = TRUE, nrepeat = 50)

# decided to re-do using 5 folds as per warning messages below. Consider a number of folds lower than the minimum in table(Y): 7

#############
# IMPORTANT: MUST CHOOSE WHICH CLASSIFICATION ERROR RATE TO USE!!! 
# want error rate to go down as no. comp increases
# centroids best for unbalanced data - but I will use max.dist
# should run perf multiple times and error rates should be averaged. 
#############

# The argument nrepeat indicates the number of cross-validation performed, and the performance are averaged across those repeats. Ideally you should choose the folds value so that the learning set is large enough, and so that the test set includes 􀁤 5 samples. Also consider increasing nrepeat when folds is small. Alternatively use leave-one-out cross validation validation = ’loo’ and nrepeat is not needed.


head(perf.plsda$error.rate)

write.table(perf.plsda$error.rate, file="TS1015_perf.plsda_error.rate_2017.07.20.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')

pdf("TS1015 plot of the classification error rate 2017.07.20.pdf", onefile=FALSE)
plot(perf.plsda, overlay = 'measure', sd=TRUE)
dev.off()

# Above is the plot of the classification error rate averaged across the 5 folds and the 10 repeated CV for all prediction distances (maximum, centroid and Mahalanobis). BER stands for balanced error rate, which accounts for unbalanced number of samples per group, , which is the case in this example. We can choose ncomp = 3 or 4 (depending on the standard deviation error bars).



# IMPORTANT: MUST CHOOSE WHICH CLASSIFICATION ERROR RATE TO USE!!! 



# We next fit a smaller PLS-DA model with ncomp = 3 and visualise the samples projected onto these components.
TS1015.plsda <- plsda(X, Y, ncomp = 3)


pdf("TS1015 PCA ncomp3 2017.07.20.pdf", onefile=FALSE)
plotIndiv(TS1015.plsda, comp = c(1,2),
group = group$x, ind.names = TRUE,
ellipse = TRUE, legend = TRUE, title = 'TS1015 ncomp3, PLSDA comp 1 - 2 2017.07.20')
dev.off()

# We observe a clear separation of the " treatement groups" compared to an unsupervised PCA sample plot. This is to be expected since the PLS-DA model includes the ~group information for each individual sample. Confidence ellipses for each class are plotted to highlight the strength of the discrimination (confidence level set to 0.95 per default). However, we observe many of the 2308 genes in X are noisy or uninformative to characterize the different classes. 

# NEXT
# The sPLS-DA analysis will help refine the sample clusters and select a small subset of variables relevant to discriminate each class.


#-----------------------------------------------------
#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink_2017.07.20.txt")
sessionInfo()
sink()


# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo_2017.07.20.txt")


####################################################################################
# bioinformatic analysis of Cel-seq2 data using sPLS-DA
# Code by Tahsha E. Say
#
# Part ii) 
#
#
####################################################################################

setwd("/opsin/u/tahsha/TSay_R_wd/splsda_wd")
getwd()

# move files to piwi
# read in sample file
samples=read.table("TS1015_Design.txt", header=TRUE,row.names=1)
# or
# group file created above:
group=read.table("TS1015_group.txt", header=TRUE,row.names=1)


### to change
TS1015=read.table("/opsin/u/tahsha/TSay_R_wd/TS1015_transformed_data/TS1015_DESeq2_vsdMat_Blind_2016.11.21.txt", header=TRUE,row.names=1)
#header=TRUE indicates that the first line contains column names and row.names=1 means that the first column should be used as row names.
head(TS1015) # print the first 7 rows of the file
tail(TS1015)
dim(TS1015)  # print the dimensions of the tabl


# load required packages
library('mixOmics', lib.loc="/opsin/u/tahsha/R/x86_64-redhat-linux-gnu-library/3.3")
 #includes ggplot
#library("dplyr")tra

# samples <- samples[-28,] ##remove a dodgy sample
X <- TS1015
X <- t(X) #transpose the data


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ to change
#Y <- samples$treatment
Y <- group$x
summary(Y)

###################
# sPLS-DA analysis
###################

# The PLS-DA model is built on all ncol(X) genes in , many of which may be uninformative to characterise the different classes. The sPLS-DA analysis will identify a small subset of genes that best discriminate the classes.

# Tuning sPLS-DA

# sPLS-DA selects the most predictive/discriminative features in the data that help classifying the samples.
# sPLS-DA is a special case of sparse PLS, where the l1 penalization is solely applied on the loading vector associated to the X data set, see [3] The parameters to choose by the user here is the number of components or dimensions ncomp and the number of variables to choose in the X data set


# keepX variables: the number of variables to select on each dimension (PC) 

# keepX. Using the function tune.splsda(), the tuning step outputs the number of variables to select, while the actual selected variables are output in the final model run on the entire data set. The tuning step will retain the first c(1:10, seq(20, 100, 10)) keepX values that gives the lowest average error rate and set that value for the next iteration. The tuning is being performed one component at the time inside the function and the optimal number of variables to select is automatically retrieved for each component. We set ncomp = 5 (ncomp + 1) and we used 10-fold cross validation (folds = 5 repeated 10 times) and specify a prediction distance (dist) to predict class membership across all CV runs. keepX is the number of variables to select on each dimesnsion.


#this chunk takes ~ 6 min to run
set.seed(32) # for reproducibility of the outputs of this code that performs random cross-validation sampling. To be removed in pro
# grid of possible keepX values that will be tested for comp 1 and comp 2
# T: think this corresponds to the grid along the X axis on plot!!!
list.keepX <- c(1:10, seq(20, 200, 10))
# [1]   1   2   3   4   5   6   7   8   9  10  20  30  40  50  60  70  80  90 100



# to speed up computation in this example we choose 5 folds repeated 10 times:
tune.splsda.TS1015 <- tune.splsda(X, Y, ncomp = 2, validation = 'Mfold', folds = 6,
progressBar = TRUE, dist = 'max.dist',
test.keepX = list.keepX, nrepeat = 100) #nrepeat 50-100
# Note: For a thorough tuning step, the following code should be repeated 50-100 times (as above).


head(tune.splsda.TS1015$error.rate)
# This output globally shows that 4 components are sufficient to achieve the lowest classification error rate in the sparse model:

# We display the mean classification error rate on each component. Each component is conditional on the previous components calculated with the optimal number of selected variables. The plot shows the optimal' number of variables to select as well as the number of components ncomp. The diamond on the plot indicates the best keepX value to achieve the lowest error rate per component.


#####################
# Splsda: Tuning step
#####################

tune.splsda.TS1015$choice.keepX

write.table(tune.splsda.TS1015$choice.keepX, file="TS1015_ptune.splsda.TS1015$choice.keepX_2017.07.20.txt", quote=FALSE, sep="\t", col.names = NA, row.names = TRUE, eol = '\n')

pdf("TS1015 tune.splsda 2017.07.20.pdf", onefile=FALSE)
plot(tune.splsda.TS1015, optimal = TRUE, sd = TRUE)
dev.off()


#---------- 
pca.bar <- pca(X, ncomp = 3, center = TRUE, scale = TRUE)

pdf("TS1015 PCA barplot 2017.07.20.pdf", onefile=FALSE)
plot(pca.bar)
dev.off()



# These results will depend on how fine your tuning grid list.keepX is, as well as the values chosen for folds and nrepeat. Therefore the performance of the final model, as well as the stability of the selected variables across the different folds should be examined, see below.

# The graphic above shows that the error rate decreases when more components are included in sPLS-DA. To obtain a more reliable estimation of the error rate, the number of repeats should be increased (between 50 to 100). This type of graph helps choosing not only the 'optimal' number of variables to select confirm the number of components ncomp . Indeed, when a sufficient number of components have been added, the error rate will stop decreasing. 

#The addition of the fourth component is probably not necessary here as the third and fourth component seem to leave to similar error rates.

#-----------------------------------------------------
#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink.txt")
sessionInfo()
sink()


# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


# ---------------------------------------------------
# Save results for use in subsequent parts/ scripts
save(splsda.TS1015, file = "Environment_TS1015_splsda_b_2017.06.19.RData")
# --------------------------------------------------- date will change automatically



####################################################################################
# bioinformatic analysis of Cel-seq2 data using sPLS-DA
# Code by Tahsha E. Say
#
# Part iii) 
#
#
####################################################################################

setwd("/opsin/u/tahsha/TSay_R_wd/splsda_wd")
getwd()

samples=read.table("TS1015_Design.txt", header=TRUE,row.names=1)
# or
# group file created above:
group=read.table("TS1015_group.txt", header=TRUE,row.names=1)


### to change
TS1015=read.table("/opsin/u/tahsha/TSay_R_wd/TS1015_transformed_data/TS1015_DESeq2_vsdMat_Blind_2016.11.21.txt", header=TRUE,row.names=1)
#header=TRUE indicates that the first line contains column names and row.names=1 means that the first column should be used as row names.
head(TS1015) # print the first 7 rows of the file
tail(TS1015)
dim(TS1015)  # print the dimensions of the tabl


# load required packages
library('mixOmics', lib.loc="/opsin/u/tahsha/R/x86_64-redhat-linux-gnu-library/3.3")
 #includes ggplot
#library("dplyr")tra

# samples <- samples[-28,] ##remove a dodgy sample
X <- TS1015
X <- t(X) #transpose the data


#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ to change
#Y <- samples$treatment
Y <- group$x
summary(Y)


#####################
# Final sPLS-DA model
#####################

# optimal number of variables to select on 3 comps:

# to change:
select.keepX = c(180,20) #from tuning step above
splsda.TS1015 <- splsda(X, Y, ncomp = 2, keepX = select.keepX)


pdf("TS1015 final sPLS-DA model -group 2017.07.25.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = c(1,2),
group = group$x, ind.names = TRUE,
ellipse = TRUE, legend = TRUE,
title = 'TS1015 final sPLS-DA comp 1 - 2')
dev.off()


pdf("TS1015 final sPLS-DA model -treatment 2017.07.25.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = c(1,2),
group = samples$treatment, ind.names = TRUE,
ellipse = TRUE, legend = TRUE,
title = 'TS1015 final sPLS-DA comp 1 - 2')
dev.off()


pdf("TS1015 final sPLS-DA model -time 2017.07.25.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = c(1,2),
group = samples$time, ind.names = TRUE,
ellipse = TRUE, legend = TRUE,
title = 'TS1015 final sPLS-DA comp 1 - 2')
dev.off()


pdf("carmels_TS1015_splsda_time_2016.06.12_2.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = 1:2,  ind.names = TRUE, group=samples$treatment,
add.legend=T, main = 'TS0616 XY-subspace', plot.ellipse= TRUE) # include sample names
dev.off()

comp = 1
loadingsX1 = abs(splsda.TS1015$loadings$X[, comp])
names(sort(loadingsX1, decreasing = T)[1:50])
PC1_transcripts <- which(splsda.TS1015$loadings$X[, comp] != 0)
write.csv( as.data.frame(PC1_transcripts), file="TS1015_PC1_transcripts_2017.07.25.csv")

comp = 2
PC2_transcripts <- which(splsda.TS1015$loadings$X[, comp] != 0)
PC2_transcripts_sorted <- names(sort(loadingsX1, decreasing = T)[1:50]) 
write.csv( as.data.frame(PC2_transcripts), file="TS1015_PC2_transcripts_2017.07.25.csv")


#-----------------------------------------------------
#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink.txt")
sessionInfo()
sink()


# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


# ---------------------------------------------------
# Save results for use in subsequent parts/ scripts
save(splsda.TS1015, file = "Environment_TS1015_splsda_c_2017.06.19.RData")
# --------------------------------------------------- 
# ----------------------------------------------------------------------------------------


####################################################################################
# bioinformatic analysis of Cel-seq2 data using sPLS-DA
# Code by Tahsha E. Say
#
# Part iv) component analysis
#
#
####################################################################################

setwd("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_q02.00_splsda_2017.06.14_tuning/splsda_2017.07.26_q03.01iii/Pretty_PCA") # to change

load(file = "/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_R_Scripts_wd_03.00/TS1015_q02.00_splsda_2017.06.14_tuning/splsda_2017.07.26_q03.01iii/Environment_TS1015_splsda_c_2017.06.19.RData") # to change

# help on editing plots
# http://mixomics.org/graphics/sample-plots/plotindiv/

# move files to piwi
# read in sample file
samples=read.table("../TS1015_Design.txt", header=TRUE,row.names=1)
# or
# group file created above:
group=read.table("../TS1015_group.txt", header=TRUE,row.names=1)

#deepskyblue
col.list <- c("deepskyblue3", "orange", "springgreen2", "snow4", "maroon2", "darkorchid2") ## to change


pdf("TS1015 final sPLS-DA model -group 2017.08.01_nolabels_pretty.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = c(1,2),
group = group$x, ind.names = FALSE,
ellipse = TRUE, legend = TRUE,
pch = 20, # dots
title = 'TS1015 final sPLS-DA',
legend.position = "right",
col.per.group = col.list,
)
dev.off()


# load required packages
library('mixOmics')

pdf("TS1015 final sPLS-DA model -group 2017.08.01_nolabels_style.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = c(1,2),
group = group$x, ind.names = FALSE,
ellipse = TRUE, legend = TRUE,
style="graphics",
pch = 20, # dots
legend.position = "right",
title = 'TS1015 final sPLS-DA comp 1 - 2'
)
dev.off()

pdf("TS1015 final sPLS-DA model -group 2017.08.01_nolabels_lattice.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = c(1,2),
group = group$x, ind.names = FALSE,
ellipse = TRUE, legend = TRUE,
pch = 20, # dots
title = 'TS1015 final sPLS-DA',
legend.position = "right",
style = 'lattice'
)
dev.off()


pdf("TS1015 final sPLS-DA model -time 2017.07.25.pdf", onefile=FALSE)
plotIndiv(splsda.TS1015, comp = c(1,2),
group = samples$time, ind.names = FALSE,
ellipse = TRUE, legend = TRUE,
pch = 20, # dots
title = 'TS1015 final sPLS-DA comp 1 - 2',
style="graphics"
)
dev.off()



########################################################
# customise the PCA plot using the ggplot function # # # 
########################################################

pdf("TS1015_2016.10.31_PCA_plot_rldt1_grid.pdf")
print(plotPCA(rldt1, intgroup=c("treatment")))		
# intgroup are the interesting groups for labeling the samples
# they tell the function to use them to choose colours

### remember to end the line with an operator so R knows that somepthing is coming (you are spreading the code over multiple lines)
  
data <- plotPCA(splsda.rld, intgroup=c("treatment"), returnData=TRUE)
percentVar <- round(100 * attr(data, "percentVar"))
ggplot(data, aes(PC1, PC2, color=treatment)) + 
geom_point(size=3) +
xlab(paste0("PC1: ",percentVar[1],"% variance")) +
ylab(paste0("PC2: ",percentVar[2],"% variance")) +
scale_colour_manual(values = col.list) +
#opts(panel.background = theme_rect(fill='white', colour='black')) + 
#opts(panel.grid.major = none, panel.grid.minor)
theme(axis.line = element_line(colour = "black"),
  panel.grid.major.x = element_blank(),
  panel.grid.major.y = element_blank(),
  panel.grid.minor.x = element_blank(),
  panel.grid.minor.y = element_blank(),
  panel.background = element_blank()
)

dev.off() 
          
          

####################################################################################
# bioinformatic analysis of Cel-seq2 data using sPLS-DA
# Code by Tahsha E. Say
#
# Part v) heatmap of sPLS-DA genes
#
#
####################################################################################


setwd("/opsin/u/tahsha/TSay_R_wd/f08.00_heatmap_splsda")
getwd()


# add library dependencies
library("DESeq2") # start R session and load the DESeq2_2016.11.28 package
library("ggplot2")


#######################
# to change
#######################
vsd=read.table("TS1015_PC1_transcripts_2017.07.25_vsd_Blind.txt", header=TRUE,row.names=1);

#vsdNTBL=read.table("/Users/tahshasay/Documents/Binary_Data_from_R/TS1015_splsda_q02.01iii_PC1_100_DEG_LRT_transformed_data/TS1015_splsda_q02.01iii_PC1_100_DEG_LRT_DESeq2_group_vsdMat_NotBlind_2016.10.24_Rv3.3.1.txt", header=TRUE,row.names=1)


vsdNTBL=read.table("TS1015_PC1_transcripts_2017.07.25_vsd_notBlind.txt", header=TRUE,row.names=1)


##############################
# re-order TRANSFORMED DATA used to create heatmaps to make ordering intuitive.
# best to have "nat" treatment in the middle (close to either light or dark dep on time)
# Here I use R to Reorder matrix columns by matching colnames to list of string
# Ref: http://stackoverflow.com/questions/25446714/r-reorder-matrix-columns-by-matching-colnames-to-list-of-string

col.order <- c("r0t1lgt1_15","r0t1lgt2_16","r0t1lgt3_17","r0t1lgt4_18","r0t1lgt5_19","r0t1lgt6_20","r0t1lgt7_21","r0t1nat1_01","r0t1nat2_02","r0t1nat3_03","r0t1nat4_04","r0t1nat5_05","r0t1nat6_06","r0t1nat7_07","r0t1drk1_08","r0t1drk2_09","r0t1drk3_10","r0t1drk4_11","r0t1drk5_12","r0t1drk6_13","r0t1drk7_14","r0t9lgt1_36","r0t9lgt2_37","r0t9lgt3_38","r0t9lgt4_39","r0t9lgt6_41","r0t9lgt7_42","r0t9lgt5_40","r0t9nat1_22","r0t9nat2_23","r0t9nat3_24","r0t9nat4_25","r0t9nat5_43","r0t9nat6_26","r0t9nat7_27","r0t9nat8_28","r0t9drk1_29","r0t9drk2_30","r0t9drk3_31","r0t9drk4_32","r0t9drk5_33","r0t9drk6_34","r0t9drk7_35")

vsd <- vsd[,col.order]
#rld <- rld[,col.order]

##########################################################################################
#
# Choose heatmap colours !!!!
#
##########################################################################################

library(ggplot2)
library(gplots) 
library(pheatmap)
library(RColorBrewer)
library("colorRamps") # for blue2red 


### remember to scale heat map by row (especially if you sorted by total expression level).  This will make the trend (down over different genes for each larva less obvious).  But will make your Left to Right trend MORE obvious.

##############################
# Diverging colours
##############################

#hmcol <- colorRampPalette(c("blue",'white','firebrick2'))(n=1000)
#hmcol = blue2red(400) # deceptive becuase lightt blue and yellow are sim exp level.
#hmcol <- colorRampPalette(rev(brewer.pal(9, "Spectral")))(1000)
# default used for heatmpa3.  hmcol <- colorRampPalette(c("blue",'white','red'))(n=1000)



##############################
# Sequential
##############################
#hmcol <- colorRampPalette(brewer.pal(9, "Blues")) (2000)

#---------------
library(viridis)
#hmcol <- rev(viridis(256))
hmcol <- rev(magma(256))
#hmcol <- rev(plasma(256))
#inferno
#---------------


#hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(100)
#hmcol <- colorRampPalette(brewer.pal(9, "Set3"))(100) ## to change
#hmcol <- colorRampPalette(brewer.pal(9, "BuGn"))(100) ## to change

#colnames(vsdMatBLT) <- with(colData(dds), paste(time, treatment, sep+" : "))
#cellheight = 6, works when pdf A4 pages but not otherwise



#----------------------------------------------------------------------------------------


########################################################
# Display row and color annotations
########################################################

# specify colours for annotating heatmap
#ann_col_colour = list(condition = c("orange", "green", "black"), time = c("snow2", "snow4"))

# read in sample file
Design=read.table("TS1015_Design.txt", header=TRUE,row.names=1)


annot_col = Design

# col <- c("black","orange", "springgreen4") ## to change

ann_colours = list(
    treatment = c(lgt = "orange", nat ="springgreen4", drk = "black"),
    time = c(t1 = "snow2", t9 = "snow4")
)


########################################################
# vsd    
########################################################

# add this before onefile
# width = 8, height = 10,  for big size

# if you want to see the row names
# cellheight = 5, show_rownames=TRUE,

pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_vsd_scalerow_clusteringboth.pdf",   onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(vsd, col=hmcol, 
border_color=NA, 
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=TRUE, 
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_vsd_scalerow_clusteringcols.pdf",   onefile=FALSE)
pheatmap(vsd, col=hmcol,  
border_color=NA, 
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=TRUE,
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 

# cluster rows
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_vsd_scalerow_clusteringrows.pdf",   onefile=FALSE)
pheatmap(vsd, col=hmcol,
border_color=NA,   
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=TRUE, 
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()


##########################################################################################
#
# Heatmap DEG. Visualising differential expression using heat maps
# # Selene said it is best to use vsd for heat maps
#
##########################################################################################

########################################################
# vsd     not blind
########################################################

pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_vsdnotBlind_scalerow_clusteringboth.pdf",   onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(vsdNTBL, col=hmcol,
border_color=NA,  
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=TRUE,
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_vsdnoBlind_scalerow_clusteringcols.pdf",   onefile=FALSE)
pheatmap(vsdNTBL, col=hmcol,  
border_color=NA, 
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=TRUE,
fontsize_row = 3,
fontsize_col = 6,  
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 

# cluster rows
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_vsdnoBlind_scalerow_clusteringrows.pdf",   onefile=FALSE)
pheatmap(vsdNTBL, col=hmcol,
border_color=NA,   
cluster_cols=FALSE, 
cluster_rows=TRUE, 
fontsize = 6, # for legends etc
show_rownames=TRUE,
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()


########################################################
# vsd    
########################################################

# add this before onefile
# width = 8, height = 10,  for big size

# if you want to see the row names
# cellheight = 5, show_rownames=FALSE,

pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_red_vsd_scalerow_clusteringboth.pdf",   onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(vsd, col=hmcol, 
border_color=NA, 
cellheight = 1,
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=FALSE,
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_red_vsd_scalerow_clusteringcols.pdf",   onefile=FALSE)
pheatmap(vsd, col=hmcol,  
border_color=NA, 
cellheight = 1,
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=FALSE,
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 

# cluster rows
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_red_vsd_scalerow_clusteringrows.pdf",   onefile=FALSE)
pheatmap(vsd, col=hmcol,
border_color=NA, 
cellheight = 1,  
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE, 
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()


##########################################################################################
#
# Heatmap DEG. Visualising differential expression using heat maps
# # Selene said it is best to use vsd for heat maps
#
##########################################################################################


########################################################
# vsd     not blind
########################################################

pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_red_vsdnotBlind_scalerow_clusteringboth.pdf",   onefile=FALSE)
#colnames(vsd) <- paste( rld$time, 1:23, sep="-" )
pheatmap(vsdNTBL, col=hmcol,
border_color=NA, 
cellheight = 1,  
cluster_cols=TRUE, 
cluster_rows=TRUE, 
show_rownames=FALSE,
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()

#cluster cols
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_red_vsdnoBlind_scalerow_clusteringcols.pdf",   onefile=FALSE)
pheatmap(vsdNTBL, col=hmcol,  
border_color=NA, 
cellheight = 1, 
cluster_cols=TRUE, 
cluster_rows=FALSE, 
show_rownames=FALSE,
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6,  
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off() 

# cluster rows
pdf("TS1015_splsda_q02.01iii_PC1_tune190_2017.08.03_red_vsdnoBlind_scalerow_clusteringrows.pdf",   onefile=FALSE)
pheatmap(vsdNTBL, col=hmcol,
border_color=NA,
cellheight = 1,    
cluster_cols=FALSE, 
cluster_rows=TRUE, 
show_rownames=FALSE,
fontsize = 6, # for legends etc
fontsize_row = 3,
fontsize_col = 6, 
annotation = annot_col,
annotation_colors = ann_colours,
#annotation_row = annot_row, 
scale="row")
dev.off()


########################################################
# print information needed for my excel sheet      # # # 
########################################################


# Capture the screen output into a character vector and use writeLines.
writeLines(capture.output(sessionInfo()), "sessionInfo.txt")


#‘sink’ diverts R output to a connection.
sink("sessionInfo_sink.txt")
sessionInfo()
sink()








