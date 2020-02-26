###########
#### SETUP ###
###########
library(dplyr)
library(summarytools)

### DIRECTORIES ####
#Mackey Lab Computer
imageData <- read.delim("/Users/utooley/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/ImageAndDemoOnly/image03.txt")
demoData <- read.delim("/Users/utooley/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/ImageAndDemoOnly/pdem01.txt")
imageQCData <- read.delim("/Users/utooley/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/ImageAndDemoOnly/mriqc01.txt")

procdti4 <- read.delim("/Users/utooley/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/ProcDTIandMRI/pdti401.txt")
procdti <- read.delim("/Users/utooley/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/ProcDTIandMRI/pdti01.txt")
procfmri <- read.delim("/Users/utooley/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/ProcDTIandMRI/fmriresults01.txt")

outdir="~/Documents/projects/in_progress/spatial_topography_parcellations_ABCD/data/for_dylan/"
###########
#### DATA CLEANING ###
###########
procdti4 = procdti4[-1,]
procdti4 = droplevels(procdti4)
procdti4$dti_mean_motion <- as.numeric(as.character(procdti4$dti_mean_motion))

master <- arrange(procdti4, dti_mean_motion)
master$subjectkey <- as.character(master$subjectkey)
master$src_subject_id <- as.character(master$src_subject_id)

############
#### FIND LOW MOTION DTI SUBJECTS ###
############

#for prisma scanners
master <- master %>% filter(., scanner_type_pd=="Prisma")
summary(master$dti_mean_motion)
prisma_subjects <- master %>% filter(., dti_mean_motion < 0.75) %>% select(., subjectkey)

#pull their image addresses
imageData = imageData[-1,]
imageData = droplevels(imageData)
imageData$subjectkey <- as.character(imageData$subjectkey)

prisma <- imageData %>% filter(subjectkey %in% unlist(prisma_subjects))
diff_only <- prisma %>% filter(image_description=="ABCD-Diffusion-FM-AP"|image_description=="ABCD-Diffusion-FM-PA"|image_description=="ABCD-DTI")
write.csv(diff_only, paste0(outdir, "prisma_ABCD_subjects_lowest_motion.csv"))

#for prisma fit scanners
master <- arrange(procdti4, dti_mean_motion)
master$subjectkey <- as.character(master$subjectkey)
master$src_subject_id <- as.character(master$src_subject_id)
master <- master %>% filter(., scanner_type_pd=="Prisma_fit")
summary(master$dti_mean_motion)
prisma_fit_subjects <- master %>% filter(., dti_mean_motion < 0.75) %>% select(., subjectkey)

## pull image data and S3 links
prisma_fit <- imageData %>% filter(subjectkey %in% unlist(prisma_fit_subjects))
diff_only <- prisma_fit %>% filter(image_description=="ABCD-Diffusion-FM-AP"|image_description=="ABCD-Diffusion-FM-PA"|image_description=="ABCD-DTI")
write.csv(diff_only, paste0(outdir, "prisma_fit_ABCD_subjects_lowest_motion.csv"))

