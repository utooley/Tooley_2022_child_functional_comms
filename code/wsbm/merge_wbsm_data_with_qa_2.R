library(dplyr)
library(psych)
library(mgcv)
library(stringi)
library(stringr)
library(summarytools)
library(lm.beta)
library(ggplot2)
require(reshape2)
# SETUP -------------------------------------------------------------------
#Cluster mounted locally on personal computer
netdatadir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/test_sample/"
trainingsampledir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/imageData/wsbm/training_sample/"
sublistdir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/parcellation/"
qadir="~/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_36p_gsr_multruns/"
demo_dir="~/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/code/analysis-nda17/data/"
analysis_dir="~/Desktop/cluster/jux/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/"

## TEMP
netdatadir="~/Downloads/"
sublistdir="~/Downloads/"
qadir="/Users/utooley/Downloads/"
demo_dir="~/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjData/"
analysis_dir="~/Downloads/"

# Read in files -----------------------------------------------------------

file1<-read.csv(paste0(netdatadir, "n64_test_sample_wsbm_yeo_comparisons_logevidence_k7.csv"))
demoData <- readRDS(paste0(demo_dir, "nda17.Rds"))
qaData <- read.csv(paste0(qadir,"XCP_QAVARS.csv"))
subjlist <- read.csv(paste0(sublistdir, "n64_test_sample.txt"))
wsbm_assign_variance <- read.csv(paste0(trainingsampledir, "variance_in_wsbm_assign_by_yeonet.csv"))

# Data Cleaning -----------------------------------------------------------

file1<-dplyr::rename(file1, ID=Var1)
demoData <- select(demoData,subjectid:race.ethnicity)
subjlist$ID<-subjlist
qaData<-dplyr::rename(qaData, ID=id0)
file1 <- dplyr::rename(file1, logevidence=Var6)

#remove filenames from subjlist
#make ID a character vector
subjlist$ID <- as.character(subjlist$ID)
file1$ID <- as.character(file1$ID)
qaData$ID <- as.character(qaData$ID)
demoData$ID <- as.character(demoData$subjectid)

#add 'sub' prefix to the demo data IDs so it matches
#demoData$ID <- paste0("sub-",demoData$ID)
demoData$ID <- gsub("_", "", demoData$ID)

#filter out runs other than run 1 from the QA data
qaData <- filter(qaData,id1=="run-01")
#filter out those who have no net stats data
file1 <- filter(file1,system_segreg_wsbm!=0 & !is.na(system_segreg_wsbm))

#join the files
master<-right_join(subjlist,file1, by="ID")
master <- right_join(qaData,master,by="ID")
master <- right_join(demoData,master,by="ID")

#make age in months
master$age <- master$age/12
#add participation coefficients together
master$sub_partcoef_wsbm <- master$sub_partcoef_pos_wsbm + master$sub_partcoef_neg_wsbm
master$sub_partcoef_yeo <- master$sub_partcoef_pos_yeo + master$sub_partcoef_neg_yeo

# Descriptives on WSBM Data ----------------------------------------------------
view(dfSummary(master))
describe(master$deviation_edge_weights_yeo)

l2 <- lm(sub_log_evidence ~ age+sex+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(avgweight ~ age+sex+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(mean_within_sys_wsbm ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(mean_between_sys_wsbm ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(system_segreg_wsbm ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_pos_wsbm ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_neg_wsbm ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

# Descriptives on Yeo Data ------------------------------------------
l2 <- lm(mean_within_sys_yeo ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(mean_between_sys_yeo ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(system_segreg_yeo ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_pos_yeo ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(sub_partcoef_neg_yeo ~ age+sex+avgweight+relMeanRMSMotion, data=master)
summary(l2)
lm.beta(l2)

# T.tests between Yeo and WSBM ------------------------------------------
t.test(master$deviation_edge_weights_wsbm, master$deviation_edge_weights_yeo,paired=TRUE)
t.test(master$mean_within_sys_yeo, master$mean_within_sys_wsbm,paired=TRUE)

t.test(master$sub_partcoef_pos_yeo, master$sub_partcoef_pos_wsbm,paired=TRUE)
t.test(master$sub_partcoef_neg_yeo, master$sub_partcoef_neg_wsbm,paired=TRUE) 
t.test(master$sub_partcoef_yeo, master$sub_partcoef_wsbm,paired=TRUE) 
#wsbm is higher

t.test(master$mean_within_sys_yeo, master$mean_within_sys_wsbm,paired=TRUE)
t.test(master$mean_between_sys_yeo, master$mean_between_sys_wsbm,paired=TRUE) 
t.test(master$system_segreg_wsbm, master$system_segreg_yeo,paired=TRUE) 
#wsbm is lower

# Make graphs -------------------------------------------------------------
## BAR GRAPH OF NODAL VARIANCE BY YEO SYSTEM
data=data.frame(t(wsbm_assign_variance))
colnames(data) <- "avgvariance"
#RGB colors of Yeo brain
yeo_in_order_cols <- c(rgb(120, 18, 134, maxColorValue = 255), rgb(70, 130, 180, maxColorValue = 255), rgb(0, 118, 14, maxColorValue = 255), rgb(196, 58, 250, maxColorValue = 255), 
                       rgb(220,248,164, maxColorValue = 255),rgb(230, 148, 34, maxColorValue = 255),  rgb(205, 62, 78, maxColorValue = 255))
Fig_variance <- ggplot(data=data, aes(1:7,avgvariance,)) +geom_bar(fill=yeo_in_order_cols, col="black",stat="identity", size=0)

## T-TESTS BETWEEN AND WITHIN SYSTEM CONNECTIVITY
df <- melt(data.frame(master$mean_within_sys_yeo, master$mean_within_sys_wsbm))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center", color="darkblue", fill="darkblue")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                        geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_bw()+scale_color_gradient()
df <- melt(data.frame(master$mean_between_sys_yeo, master$mean_between_sys_wsbm))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center", color="lightblue", fill="lightblue")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                        geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_bw()+scale_color_gradient()
df <- melt(data.frame(master$system_segreg_wsbm, master$system_segreg_yeo))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center", color="blue", fill="blue")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                          geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_bw()+scale_color_gradient()
## T-TESTS BETWEEN WSBM AND YEO EDGE DEVIATION NO DIFF
df <- melt(data.frame(master$deviation_edge_weights_wsbm,master$deviation_edge_weights_yeo))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center", color="darkblue", fill="lightblue")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                                         geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_bw()+scale_color_gradient()
## T-TESTS BETWEEN PARTICIPATION COEF
df <- melt(data.frame(master$sub_partcoef_yeo, master$sub_partcoef_wsbm))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center", color="darkblue", fill="darkblue")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                     geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_bw()+scale_color_gradient()

df <- melt(data.frame(master$sub_partcoef_neg_yeo, master$sub_partcoef_neg_wsbm))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center", color="lightgreen", fill="lightgreen")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
            
                                                                                                                                                        
                                                                                                                                                                                                                                                                                                    geom = "crossbar", width = 0.8, col=rgb(10, 58, 150, maxColorValue = 255))+theme_bw()+scale_color_gradient()
# MoveMe Function ---------------------------------------------------------


moveMe <- function(data, tomove, where = "last", ba = NULL) {
  temp <- setdiff(names(data), tomove)
  x <- switch(
    where,
    first = data[c(tomove, temp)],
    last = data[c(temp, tomove)],
    before = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)-1))]
    },
    after = {
      if (is.null(ba)) stop("must specify ba column")
      if (length(ba) > 1) stop("ba must be a single character string")
      data[append(temp, values = tomove, after = (match(ba, temp)))]
    })
  x
}