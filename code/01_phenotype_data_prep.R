###############################
# Script to prepare the phenotype data set
# Generate Figure 2, Supplemenetal Figure 2
source("code/00_functions.R")
library_list<-c("tidyverse","labelled","stringr","plyr","dplyr","ggplot2","factoextra","survey","purrr","e1071","weights","corrplot")
lapply(library_list, require, character.only = TRUE)
# System parameters
options(scipen = 999, stringsAsFactors = F)
options(survey.lonely.psu="remove")

# read in phenotype data
pheno<-read.csv("data/metsleep_covariates_20220126.csv",header = T,sep=",")
# recode missing values in categorical variables
recode_var<-c("OCEA13","SLEA4","SLEA5","SLEA6","SLEA7","SLEA8","SLEA11")
pheno[,recode_var]<-sapply(pheno[,recode_var], function(x){ifelse(is.na(x),NA,ifelse(x%in%c(1:7),x,NA))})
# clean min_o2, avg_o2: if >100 set to missing, if <=0 set to missing
pheno[which(pheno$SLPA91>100),"SLPA91"]<-100
pheno[which(pheno$SLPA92<=0),"SLPA92"]<-NA

pheno$WEIGHT<-pheno$WEIGHT_FINAL_NORM_OVERALL
# Assess outliers in the phenotype data set
admin.var<-c("ID","PSU_ID","STRAT","WEIGHT_FINAL_NORM_OVERALL","WEIGHT")
pheno_outlier<-pheno[,!colnames(pheno)%in%admin.var]%>%
  mutate_if(is.numeric, remove_outliers, na.rm = TRUE)

# recode phenotype variables
pheno<-pheno%>% 
  dplyr::mutate(
    CENTER=factor(CENTER, levels = c("B","C","M","S"), labels=c("Bronx","Chicago","Miami","San Diego")),
    ALCOHOL_USE=factor(ALCOHOL_USE,levels=c(1,2,3),labels=c("never","former","current")),
    CIGARETTE_USE=factor(CIGARETTE_USE,levels=c(1,2,3),labels=c("never","former","current")),
    GENDER=factor(GENDER,levels = c("F","M"),labels = c("Female","Male")),
    OSA=ifelse(is.na(SLPA54)==0,ifelse(SLPA54==15,1,0),NA),
    OSA_status=factor(OSA,levels=c(0,1),labels = c("Non-OSA","OSA")),
    Per90=transformation(SLPA97),
    Per90_status=factor(Per90,levels=c(0,1),labels=c("Non-hypoxia","Hypoxia")),
    osa_ess=ifelse(is.na(SLPA54)|is.na(ESS),NA,ifelse(SLPA54>15&ESS>10,1,0)),
    snore=ifelse(SLEA13==4,1,ifelse(SLEA13%in%c(1,2,3),0,NA)), #always or almost always (6-7 nights a week)=1, less than 6-7 nights a weeek=0,dont know and others =NA
    pill=ifelse(!is.na(SLEA8),ifelse(SLEA8%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
    restless=ifelse(!is.na(SLEA11),ifelse(SLEA11%in%c("3","4",3,4),1,0),NA), # restless or very restless=1
    backtosleep=ifelse(!is.na(SLEA7),ifelse(SLEA7%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
    earlywake=ifelse(!is.na(SLEA6),ifelse(SLEA6%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
    freqwake=ifelse(!is.na(SLEA5),ifelse(SLEA5%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
    fallasleep=ifelse(!is.na(SLEA4),ifelse(SLEA4%in%c("4","5",4,5),1,0),NA), # 3 or more times a week=1
    restless_fct=factor(restless,levels=c(0,1),labels=c("No","Yes")),
    backtosleep_fct=factor(backtosleep,levels=c(0,1),labels=c("No","Yes")),
    earlywake_fct=factor(earlywake,levels=c(0,1),labels=c("No","Yes")),
    freqwake_fct=factor(freqwake,levels=c(0,1),labels=c("No","Yes")),
    fallasleep_fct=factor(fallasleep,levels=c(0,1),labels=c("No","Yes")),
    snore_fct=factor(snore,levels=c(0,1),labels=c("No","Yes")),
    essgt10_fct=factor(essgt10,levels=c(0,1),labels=c("No","Yes")),
    slpa36gt15_fct=factor(slpa36gt15,levels=c(0,1),labels=c("No","Yes")),
    slpa36gt5_fct=factor(slpa36gt5,levels=c(0,1),labels=c("No","Yes")),
    slpa54gt15_fct=factor(slpa54gt15,levels=c(0,1),labels=c("No","Yes")),
    slpa54gt5_fct=factor(slpa54gt5,levels=c(0,1),labels=c("No","Yes")),
    bkgrd1_clean=ifelse(BKGRD1_C7%in%c("","Q"),NA,as.character(BKGRD1_C7)),
    background=factor(bkgrd1_clean,levels=c(0:6),labels=c("Dominican","Central_American","Cuban","Mexican","Puerto_Rican", "South American","Multi")),
    # 0: Dominican
    # 1: Central American
    # 2: Cuban
    # 3: Mexican
    # 4: Puerto Rican
    # 5: South American
    # 6: More than one/Other heritage
    # Q: Unknown
    INCIDENT_DM_V1V2=ifelse(is.na(DIABETES3),NA, # if baseline is missing, then incident DM is missing
                            ifelse(is.na(DIABETES3_INDICATOR_V2),# if follow up is missing, then incident DM is missing
                                   NA,
                                   ifelse(DIABETES3%in%c(1,2)&DIABETES3_INDICATOR_V2==1,1,0))), # if baseline is 1 or 2, while follow up indicator is 1, then incident =1, otherwise incident =0
    DIABETES3_INDICATOR=ifelse(!is.na(DIABETES3),ifelse(DIABETES3==3,1,0),NA),
    DIABETES3_INDICATOR=factor(DIABETES3_INDICATOR,levels=c(0,1),labels=c("No","Yes")),
    DIABETES3_INDICATOR_V2=factor(DIABETES3_INDICATOR_V2,levels=c(0,1),labels=c("No","Yes")),
    OBSES=ifelse(!is.na(BMI),ifelse(BMI>=30,1,0),NA),
    OBSES_FCT=factor(OBSES,levels=c(0,1),labels=c("No","Yes")),
    HYPERTENSION=factor(HYPERTENSION2,levels=c(0,1),labels=c("No","Yes")), # overwrite hypertension with factor hypertension2
    HIGH_TOTAL_CHOL=factor(HIGH_TOTAL_CHOL,levels=c(0,1),labels=c("No","Yes")),
    MED_ANTIHYPERT=factor(MED_ANTIHYPERT,levels=c(0,1),labels=c("No","Yes")),
    MED_LLD=factor(MED_LLD,levels=c(0,1),labels=c("No","Yes")),
    INSULIN_FAST=ifelse(!is.na(INSULIN_FAST),ifelse(INSULIN_FAST%in%c("N","L",""),NA,as.numeric(INSULIN_FAST)),NA),
    LABA69=ifelse(!is.na(LABA69),ifelse(LABA69%in%c("N","L",""),NA,as.numeric(as.character(LABA69))),NA), # fix to as.numeric rounding up numbers
    LABA68=ifelse(!is.na(LABA68),ifelse(LABA68%in%c("N","L",""),NA,as.numeric(as.character(LABA68))),NA), # fix to as.numeric rounding up numbers
    LABA66=ifelse(!is.na(LABA66),ifelse(LABA66%in%c("N","L",""),NA,as.numeric(as.character(LABA66))),NA), # fix to as.numeric rounding up numbers
    MED_ANTIDIAB=factor(MED_ANTIDIAB,levels=c(0,1),labels=c("No","Yes")),
    HTN_NEW=ifelse(!is.na(HYPERTENSION2),
                   ifelse(!is.na(HYPERTENSION2_V2),
                          ifelse(HYPERTENSION2==0&HYPERTENSION2_V2==1,1,0),
                          0),
                   NA),
    incident_cvd=INCIDENT_CVD_V1V2,
    incident_dm=INCIDENT_DM_V1V2,
    incident_htn=HTN_NEW,
    baseline_htn=ifelse(!is.na(HYPERTENSION),ifelse(HYPERTENSION=="Yes",1,0),NA),
    baseline_dm=ifelse(!is.na(DIABETES3_INDICATOR),ifelse(DIABETES3_INDICATOR=="Yes",1,0),NA)
  )

# PCA for SDB variables
sdb_pca<-c("SLPA36","SLPA54","SLPA39","SLPA42","SLPA91","SLPA92","SLPA97","event_length_sec","hypoxicburden_harmonized")

sdb_active<-pheno[,c(sdb_pca,"ID","PSU_ID","STRAT","WEIGHT_FINAL_NORM_OVERALL")]%>%
  dplyr::rename(avgSpO2=SLPA92,
                minSpO2=SLPA91,
                avgEventLength=event_length_sec,
                REI3=SLPA54,
                REI0=SLPA36,
                HB=hypoxicburden_harmonized,
                Per90=SLPA97
  )
sdb_pca_rename<-c("avgSpO2","minSpO2","REI3","REI0","Per90","avgEventLength","HB")

# rank transform all traits
sdb_active[,sdb_pca_rename]<-apply(sdb_active[,sdb_pca_rename],2,rank_normalisation)
# survey weighted PCA
survey_design_sdb<-svydesign(id=~PSU_ID, strata=~STRAT,weights=~WEIGHT_FINAL_NORM_OVERALL,data=sdb_active)
res.pca <- svyprcomp(~avgSpO2+minSpO2+REI3+REI0+Per90+avgEventLength+HB, design=survey_design_sdb,scale=TRUE,scores=TRUE)

sdb_pca<-as.data.frame(res.pca$x)%>%
  cbind(.,res.pca$design$variables[-res.pca$naa,"ID"])

colnames(sdb_pca)[ncol(sdb_pca)]<-"ID"
sdb_pca<-sdb_pca[,c("ID","PC1","PC2","PC3")]
colnames(sdb_pca)[2:4]<-c("sdb_pc1","sdb_pc2","sdb_pc3")
sdb_pca[,2:4]<-sdb_pca[,2:4]%>%scale()
# flip the direction of SDB PC2
sdb_pca$sdb_pc2<-sdb_pca$sdb_pc2*(-1)

# merge with the phenotype dataset
pheno_pca<-merge(pheno,sdb_pca,by="ID",all.x=T)
# save the intermediate file
saveRDS(pheno_pca, file="data/phenotype_pca.rds",version=2) 

# Create Supplemental Figure S2 Scree Plot
fviz_eig(res.pca,addlabels = T)

# export the scree plot
jpeg("output/suppl_figs2.jpg")
print(fviz_eig(res.pca,addlabels = T))
dev.off()

# Generate Figure 2 Correlation Matrix
pheno_corr<-pheno_pca%>%
  dplyr::mutate(Age=AGE,
                REI0=SLPA36,
                REI3=SLPA54,
                HB=hypoxicburden_harmonized,
                Per90=SLPA97,
                avgSpO2=SLPA92,
                minSpO2=SLPA91,
                avgEventLength=event_length_sec,
                SDBPC1=sdb_pc1,
                SDBPC2=sdb_pc2
  )
survey.design.corr<-svydesign(id=~PSU_ID, strata=~STRAT,weights=~WEIGHT_FINAL_NORM_OVERALL,data=pheno_corr)

sdbpc_corr<-jtools::svycor(~Age+BMI+REI3+REI0+HB+avgEventLength+avgSpO2+minSpO2+Per90+SDBPC1+SDBPC2, design = survey.design.corr, digits = 4,na.rm=T, sig.stats = TRUE, bootn = 2000, mean1 = TRUE)
# to ensure the values are within [-1,1] range for plotting
sdbpc_corr$cors[sdbpc_corr$cors>1]<-1

figure2<-corrplot(sdbpc_corr$cors, method="number",type = "upper", order = "original", p.mat = sdbpc_corr$p.values,insig="pch",tl.col = "black", tl.srt = 90)

# export the correlation plot
jpeg("output/fig2.jpg", width = 800, height = 600)
print(corrplot(sdbpc_corr$cors, method="number",type = "upper", order = "original", p.mat = sdbpc_corr$p.values,insig="pch",tl.col = "black", tl.srt = 90))
# Sleep for a second to ensure the plot is written to the file
Sys.sleep(10)
dev.off()
