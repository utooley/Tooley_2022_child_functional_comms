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
netdatadir="~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/"
sublistdir="~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/subjLists/"
qadir="/Users/utooley/Desktop/cluster/jux/mackey_group/public_data/ABCD/bids_fast_track_site14/derivatives/xcpEngine_despike_onerun/"
demo_dir="~/Documents/projects/in_progress/spatial_topography_ABCD/data/subjData/"
analysis_dir="~/Desktop/cluster/picsl/mackey_group/Ursula/projects/in_progress/spatial_topography_parcellations_ABCD/data/"

# Read in files -----------------------------------------------------------

file1<-read.csv(paste0(netdatadir, "n27_one_run_only_net_meas_Schaefer400.csv"))
demoData <- readRDS(paste0(demo_dir, "nda17.Rds"))
qaData <- read.csv(paste0(qadir,"n27_quality_one_run_meanRMSonly.csv"))
qaData <- read.csv("~/Downloads/n27_quality_one_run_meanRMSonly.csv")
subjlist <- read.csv(paste0(sublistdir, "n27_cohort_file_one_run_only_21019.csv"))
#file2<-read.table(paste0(qadir, "group_bold.tsv"), sep = '\t', header = TRUE)

strengths<-read.csv(paste0(netdatadir, "n27_one_run_only_strengths.csv"))
strengths_null1<-read.csv(paste0(netdatadir, "n27_one_run_only_strengths_null1.csv"))
strengths_null2<-read.csv(paste0(netdatadir, "n27_one_run_only_strengths_null2.csv"))

# Data Cleaning -----------------------------------------------------------

file1<-dplyr::rename(file1, ID=Var1)
demoData <- select(demoData,subjectid:race.ethnicity)
subjlist<-dplyr::rename(subjlist, ID=id0)
qaData<-dplyr::rename(qaData, ID=id0)

#remove filenames from subjlist
#make ID a character vector
subjlist <- as.character(subjlist)
file1$ID <- as.character(file1$ID)
qaData$ID <- as.character(qaData$ID)
demoData$ID <- as.character(demoData$subjectid)

#add 'sub' prefix to the demo data IDs so it matches
demoData$ID <- paste0("sub-",demoData$ID)
demoData$ID <- gsub("_", "", demoData$ID)

#join the files
master<-right_join(subjlist,file1, by="ID")
master <- right_join(qaData,master,by="ID")
master <- right_join(demoData,master,by="ID")

#make age in months
master$age <- master$age/12

# Descriptives on Data ----------------------------------------------------
view(dfSummary(master))
describe(master$avgclustco_both)
describe(master$cpathleng)

l2 <- lm(avgclustco_both ~ age+sex+race.eth+avgweight+relMeanRMS, data=master)
summary(l2)
lm.beta(l2)

l2 <- lm(cpathleng ~ age+sex+race.eth+avgweight+relMeanRMS, data=master)
summary(l2)

l2 <- lm(cpathleng_thresh ~ age+sex+race.eth+avgweight+relMeanRMS, data=master)
summary(l2)
lm.beta(l2)

# Null model comparisons on data ------------------------------------------
t.test(master$avgclustco_both, master$avgclustco_both_null1,paired=TRUE)
t.test(master$avgclustco_both, master$avgclustco_both_null2,paired=TRUE)

t.test(master$cpathleng, master$cpathleng_null1,paired=TRUE)
t.test(master$cpathleng, master$cpathleng_null2,paired=TRUE)

t.test(master$cpathleng_thresh, master$cpathleng_thresh_null1,paired=TRUE)
t.test(master$cpathleng_thresh, master$cpathleng_thresh_null2,paired=TRUE)


# Make graphs -------------------------------------------------------------
## CLUST CO NULL MODELS
df <- melt(data.frame(master$avgclustco_both,master$avgclustco_both_null1, master$avgclustco_both_null2))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                                            geom = "crossbar", width = 0.8, col=rgb(101, 58, 150, maxColorValue = 255))+theme_bw()
## CPATHLENG NULL MODELS
df <- melt(data.frame(master$cpathleng_thresh,master$cpathleng_thresh_null1, master$cpathleng_thresh_null2))
colnames(df) <- c("which", "value")
print(df)
#dotplot
ggplot(df, aes(which, value))+ geom_dotplot(binaxis="y", mapping = aes(which, value),stackdir="center")+stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                                                                                                                     geom = "crossbar", width = 0.8, col=rgb(101, 58, 150, maxColorValue = 255))+theme_bw()

# Strength distribution ---------------------------------------------------
strengths_pos <- select(strengths, contains("pos_"))
strengths_neg <- select(strengths, contains("neg_"))
hist(unlist(strengths), breaks = 50, col="lightblue", ylim = c(0,1500), xlim=c(0,200))
hist(unlist(strengths_neg))
hist(unlist(strengths))


strengths_pos <- select(strengths_null1, contains("pos_"))
strengths_neg <- select(strengths_null1, contains("neg_"))
hist(unlist(strengths_pos))
hist(unlist(strengths_neg))
hist(unlist(strengths_null1), breaks=50, col="lightblue", ylim = c(0,1500), xlim=c(0,200))


strengths_pos <- select(strengths_null2, contains("pos_"))
strengths_neg <- select(strengths_null2, contains("neg_"))
hist(unlist(strengths))
hist(unlist(strengths_neg))
hist(unlist(strengths_null2), breaks=50, col="lightblue",  ylim = c(0,1500), xlim=c(0,200))

# Write out Data ----------------------------------------------------------

#write the network data file back into the output folder
write.csv(master,"~/Downloads/n47_within_between_Yeo_Schaefer400_with_qa.csv")
write.csv(master,paste0(netdatadir,"n47_within_between_Yeo_Schaefer400_with_QA.csv"))

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