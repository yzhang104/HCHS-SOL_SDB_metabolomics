######################################
# Single metabolite analysis (SMA)
# Generate Table 1, Table 2, Supplemental Table S1
source("code/00_functions.R")
library_list<-c("tidyverse","labelled","stringr","plyr","dplyr","ggplot2","factoextra","survey","purrr","e1071","tableone")
lapply(library_list, require, character.only = TRUE)
# System parameters
options(scipen = 999, stringsAsFactors = F)
options(survey.lonely.psu="remove")
# read in phenotype and metabolomic data
pheno_pca<-readRDS("data/phenotype_pca.rds")
base::load("data/metabolomics_data.RData")
# Merge the phenotype and metabolomic data
# map  metabolites PARENT_SAMPLE_NAME ID with the SOL_ID (multiple PARENT_SAMPLE_NAME IDs mapped to the same SOL_ID)
metab_imp_cont_b1<-merge(id_list[,c("PARENT_SAMPLE_NAME","SOL_ID")],metab_imp_cont,by.y="ID",by.x="PARENT_SAMPLE_NAME",all.y=T)
multiple_SOL_ID_cont_b1<-metab_imp_cont_b1%>% #subset all obs with multiple obs for one SOL_ID
  group_by(SOL_ID) %>%
  filter(n()>1)
set.seed(100)
random_SOL_ID_cont_b1<-metab_imp_cont_b1%>% 
  group_by(SOL_ID)%>%
  sample_n(1) # randomly select one obs for each SOL_ID
pheno_var<-c("sdb_pc1","sdb_pc2","AGE","BMI","background","CENTER","GENDER","ID")
imp_only_data<-merge(pheno_pca[complete.cases(pheno_pca[,pheno_var]),],random_SOL_ID_cont_b1[,!colnames(random_SOL_ID_cont_b1)%in%c("PARENT_SAMPLE_NAME")],by.x="ID",by.y = "SOL_ID")
# add leading X to the chem_id col names
names(imp_only_data) <- gsub("^([0-9])", "X\\1", names(imp_only_data))
imp_only_id<-imp_only_data$ID

metab_imp_binary_b1<-merge(id_list[,c("PARENT_SAMPLE_NAME","SOL_ID")],metab_imp_binary,by.y="ID",by.x="PARENT_SAMPLE_NAME",all.y=T)
multiple_SOL_ID_binary_b1<-metab_imp_binary_b1%>% #subset all obs with multiple obs for one SOL_ID
  group_by(SOL_ID) %>%
  filter(n()>1)
set.seed(100)
random_SOL_ID_binary_b1<-metab_imp_binary_b1%>% 
  group_by(SOL_ID)%>%
  sample_n(1) # randomly select one obs for each SOL_ID
imp_bin_data<-merge(pheno_pca[complete.cases(pheno_pca[,pheno_var]),],random_SOL_ID_binary_b1[,!colnames(random_SOL_ID_binary_b1)%in%c("PARENT_SAMPLE_NAME")],by.x="ID",by.y = "SOL_ID")
# add leading X to the chem_id col names
names(imp_bin_data) <- gsub("^([0-9])", "X\\1", names(imp_bin_data))
imp_bin_id<-imp_bin_data$ID

# subset the final dataset by gender
imp_only_data_female<-imp_only_data[which(imp_only_data$GENDER=="Female"),]
imp_only_data_female_id<-imp_only_data_female$ID
imp_only_data_male<-imp_only_data[which(imp_only_data$GENDER=="Male"),]
imp_only_data_male_id<-imp_only_data_male$ID
imp_bin_data_female<-imp_bin_data[which(imp_bin_data$GENDER=="Female"),]
imp_bin_data_female_id<-imp_bin_data_female$ID
imp_bin_data_male<-imp_bin_data[which(imp_bin_data$GENDER=="Male"),]
imp_bin_data_male_id<-imp_bin_data_male$ID

# Batch 2
# map  metabolites PARENT_SAMPLE_NAME ID with the SOL_ID (multiple PARENT_SAMPLE_NAME IDs mapped to the same SOL_ID)
metab_imp_cont_b2<-merge(id_list[,c("PARENT_SAMPLE_NAME","SOL_ID")],metab_imp_cont_b2,by.y="ID",by.x="PARENT_SAMPLE_NAME",all.y=T)
multiple_SOL_ID_cont_b2<-metab_imp_cont_b2%>% #subset all obs with multiple obs for one SOL_ID
  group_by(SOL_ID) %>%
  filter(n()>1)
set.seed(100)
random_SOL_ID_cont_b2<-metab_imp_cont_b2%>% 
  group_by(SOL_ID)%>%
  sample_n(1) # randomly select one obs for each SOL_ID
imp_only_data_b2<-merge(pheno_pca[complete.cases(pheno_pca[,pheno_var]),],random_SOL_ID_cont_b2[,!colnames(random_SOL_ID_cont_b2)%in%c("PARENT_SAMPLE_NAME")],by.x="ID",by.y = "SOL_ID")
# add leading X to the chem_id col names
names(imp_only_data_b2) <- gsub("^([0-9])", "X\\1", names(imp_only_data_b2))
imp_only_id_b2<-imp_only_data_b2$ID

metab_imp_binary_b2<-merge(id_list[,c("PARENT_SAMPLE_NAME","SOL_ID")],metab_imp_binary_b2,by.y="ID",by.x="PARENT_SAMPLE_NAME",all.y=T)
multiple_SOL_ID_binary_b2<-metab_imp_binary_b2%>% #subset all obs with multiple obs for one SOL_ID
  group_by(SOL_ID) %>%
  filter(n()>1)
set.seed(100)
random_SOL_ID_binary_b2<-metab_imp_binary_b2%>% 
  group_by(SOL_ID)%>%
  sample_n(1) # randomly select one obs for each SOL_ID
imp_bin_data_b2<-merge(pheno_pca[complete.cases(pheno_pca[,pheno_var]),],random_SOL_ID_binary_b2[,!colnames(random_SOL_ID_binary_b2)%in%c("PARENT_SAMPLE_NAME")],by.x="ID",by.y = "SOL_ID")
# add leading X to the chem_id col names
names(imp_bin_data_b2) <- gsub("^([0-9])", "X\\1", names(imp_bin_data_b2))
imp_bin_id_b2<-imp_bin_data_b2$ID

# subset the final dataset by gender
imp_only_data_b2_female<-imp_only_data_b2[which(imp_only_data_b2$GENDER=="Female"),]
imp_only_data_b2_female_id<-imp_only_data_b2_female$ID
imp_only_data_b2_male<-imp_only_data_b2[which(imp_only_data_b2$GENDER=="Male"),]
imp_only_data_b2_male_id<-imp_only_data_b2_male$ID
imp_bin_data_b2_female<-imp_bin_data_b2[which(imp_bin_data_b2$GENDER=="Female"),]
imp_bin_data_b2_female_id<-imp_bin_data_b2_female$ID
imp_bin_data_b2_male<-imp_bin_data_b2[which(imp_bin_data_b2$GENDER=="Male"),]
imp_bin_data_b2_male_id<-imp_bin_data_b2_male$ID

pheno_tbl1<-pheno_pca%>%
  dplyr::mutate(
    REI3=SLPA54,
    REI0=SLPA36,
    Per90_status=factor(Per90,levels=c(0,1),labels=c("Non-hypoxia","Hypoxia")),
    baseline_dm=factor(baseline_dm,levels=c(0,1),labels=c("No","Yes")),
    baseline_htn=factor(baseline_htn,levels=c(0,1),labels=c("No","Yes")),
    incident_dm=factor(incident_dm,levels=c(0,1),labels=c("No","Yes")),
    incident_htn=factor(incident_htn,levels=c(0,1),labels=c("No","Yes")),
    fallasleep=factor(fallasleep,levels=c(0,1),labels=c("No","Yes")),
    freqwake=factor(freqwake,levels=c(0,1),labels=c("No","Yes")),
    earlywake=factor(earlywake,levels=c(0,1),labels=c("No","Yes")),
    backtosleep=factor(backtosleep,levels=c(0,1),labels=c("No","Yes")),
    pill=factor(pill,levels=c(0,1),labels=c("No","Yes")),
    restless=factor(restless,levels=c(0,1),labels=c("No","Yes")),
    snore=factor(snore,levels=c(0,1),labels=c("No","Yes")),
    sdbpc1_10=ntile(sdb_pc1,10),
    sdbpc2_10=ntile(sdb_pc2,10),
    sdbpc1_top10=ifelse(!is.na(sdbpc1_10),ifelse(sdbpc1_10==10,1,0),NA),
    sdbpc1_bottom10=ifelse(!is.na(sdbpc1_10),ifelse(sdbpc1_10==1,1,0),NA),
    sdbpc2_top10=ifelse(!is.na(sdbpc2_10),ifelse(sdbpc2_10==10,1,0),NA),
    sdbpc2_bottom10=ifelse(!is.na(sdbpc2_10),ifelse(sdbpc2_10==1,1,0),NA)
  )
# assign variable lables
tbl_list<-list(
  AGE="Age at baseline",
  GENDER="Gender",
  CENTER="Study center",
  background="Hispanic/Latino background",
  CIGARETTE_USE="Smoking status",
  ALCOHOL_USE="Alcohol drinking status",
  REI3="Respiratory Event Index (3% desat) (events/hr)",
  REI0="Respiratory Event Index (all desat) (events/hr)",
  OSA_status="OSA status",
  Per90_status="Hypoxia Status(>=5% sleep <90% saturation)",
  baseline_dm="Baseline Diabetes status (ADA)",
  baseline_htn="Baseline Hypertension status",
  incident_dm="Incident Diabetes status (ADA)",
  incident_htn="Incident Hypertension status",
  LABA70="Fasting glucose (mg/dL)",
  INSULIN_FAST="Fasting insulin (mU/L)",
  LABA69="LDL(mg/dL)",
  LABA66 ="TTotal cholesterol (mg/dL)",
  LABA68="HDL(mg/dL)",
  LABA67="Triglycerides (mg/dL)",
  SBPA5="Systolic Blood Pressure(mm Hg)",
  SBPA6="Diastolic Blood Pressure(mm Hg)",
  GPAQ_TOTAL_MET="Physical activity (MET-min/day)",
  AHEI2010="The Alternate Healthy Eating Index (2010)",
  ESS="Epworth Sleepiness Scale (ESS) total score",
  WHIIRS="Womens's Health Initiative Insomnia Rating Scale (WHIIRS) total score",
  fallasleep="Trouble falling asleep (3 or more times a week)",
  freqwake="Wake up several times at night (3 or more times a week)",
  earlywake="Wake up earlier than you plan (3 or more times a week)",
  backtosleep="Trouble getting back to sleep (3 or more times a week)",
  pill="Take sleeping pills",
  restless="Typical night’s sleep in past 4 weeks(restless or very restless)",
  snore="Self-reported snoring (6-7 nights a week)",
  essgt10="ESS>10",
  SLPDUR="Sleep duration",
  SLPA36="Respiratory Event Index (all desat)",
  SLPA91="Minimum SpO2%",
  SLPA92="Average SpO2%",
  SLPA97="Percentage sleep time with SpO2<90%",
  event_length_sec="Average length of each respiratory event (seconds)",
  SLPA111="Minimum resting heart rate during sleep",
  SLPA112="Maximum resting heart rate during sleep",
  SLPA113="Average resting heart rate during sleep",
  SLPA114="Standard deviation resting heart rate during sleep",
  hypoxicburden_harmonized="Hypoxic burden (%minute/hour)"
)

# create table one for b1 metabolomic data sets
pheno_b1<-pheno_tbl1[which(pheno_tbl1$ID%in%imp_only_data$ID),]
labelled::var_label(pheno_b1)<-tbl_list
pheno_b2<-pheno_tbl1[which(pheno_tbl1$ID%in%imp_only_data_b2$ID),]
labelled::var_label(pheno_b2)<-tbl_list
## export table 1
# variables
var<-c("AGE", "GENDER","background", "BMI", "ALCOHOL_USE", "CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010","OSA_status","REI3","REI0","event_length_sec","SLPA97", "Per90_status","hypoxicburden_harmonized"
       ,"SLPA91","SLPA92", "baseline_dm", "baseline_htn","incident_dm","incident_htn", "LABA67", "LABA68", "LABA69", "LABA70", "INSULIN_FAST", "HOMA_IR", "LABA66", "SBPA5", "SBPA6")
disturb_var<-c("WHIIRS","restless","pill","backtosleep","earlywake","freqwake","fallasleep","ESS","essgt10_fct","snore_fct","SLPDUR")
hr_var<-c("SLPA111","SLPA112","SLPA113","SLPA114")
# study design based on full data set
survey.design.total<-svydesign(id=~PSU_ID, strata=~STRAT,weights=~WEIGHT,data=pheno_tbl1)
survey.design.tbl1.pca<-subset(survey.design.total,!is.na(sdb_pc1))
survey.design.tbl1.b1<-subset(survey.design.tbl1.pca,ID%in%imp_only_data$ID&!is.na(BMI)&!is.na(CENTER)&!is.na(background))
survey.design.tbl1.b2<-subset(survey.design.tbl1.pca,ID%in%imp_only_data_b2$ID&!is.na(BMI)&!is.na(CENTER)&!is.na(background))
options(survey.lonely.psu="remove")  
 
# Export the baseline demographic characteristics tables (Table 1)
supp_tblone_weighted_b1<- print(svyCreateTableOne(vars = var, data = survey.design.tbl1.pca), missing=TRUE, varLabels = TRUE,digits =2,pDigits=2)
write.csv(supp_tblone_weighted_b1, file = "output/tables1/tables1_weighted_b1.csv",row.names = T)
supp_tblone_nonweighted_b1<- print(CreateTableOne(vars = var, data = subset(pheno_pca,!is.na(sdb_pc1))), missing=TRUE, varLabels = TRUE,digits =2,pDigits=2)
write.csv(supp_tblone_nonweighted_b1, file = "output/tables1/tables1_nonweighted_b1.csv",row.names = T)
supp_tblone_gender_weighted_b1<- print(svyCreateTableOne(vars = var, strata = "GENDER", data = survey.design.tbl1.pca), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(supp_tblone_gender_weighted_b1, file = "output/tables1/tables1_gender_weighted_b1.csv",row.names = T)
supp_tblone_gender_nonweighted_b1<- print(CreateTableOne(vars = var, strata = "GENDER", data = subset(pheno_pca,!is.na(sdb_pc1))), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(supp_tblone_gender_nonweighted_b1, file = "output/tables1/tables1_gender_nonweighted_b1.csv",row.names = T)
supp_tblone_osa_weighted_b1<- print(svyCreateTableOne(vars = var, strata = "OSA_status", data = survey.design.tbl1.pca), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(supp_tblone_osa_weighted_b1, file = "output/tables1/tables1_osa_weighted_b1.csv",row.names = T)
supp_tblone_osa_nonweighted_b1<- print(CreateTableOne(vars = var, strata = "OSA_status", data = subset(pheno_pca,!is.na(sdb_pc1))), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(supp_tblone_osa_nonweighted_b1, file = "output/tables1/tables1_osa_nonweighted_b1.csv",row.names = T)
# batch 1
tblone_weighted_b1<- print(svyCreateTableOne(vars = var, data = survey.design.tbl1.b1), missing=TRUE, varLabels = TRUE,digits =2,pDigits=2)
write.csv(tblone_weighted_b1, file = "output/table1/table1_weighted_b1.csv",row.names = T)
tblone_nonweighted_b1<- print(CreateTableOne(vars = var, data = pheno_b1), missing=TRUE, varLabels = TRUE,digits =2,pDigits=2)
write.csv(tblone_nonweighted_b1, file = "output/table1/table1_nonweighted_b1.csv",row.names = T)
tblone_gender_weighted_b1<- print(svyCreateTableOne(vars = var, strata = "GENDER", data = survey.design.tbl1.b1), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(tblone_gender_weighted_b1, file = "output/table1/table1_gender_weighted_b1.csv",row.names = T)
tblone_gender_nonweighted_b1<- print(CreateTableOne(vars = var, strata = "GENDER", data = pheno_b1), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(tblone_gender_nonweighted_b1, file = "output/table1/table1_gender_nonweighted_b1.csv",row.names = T)
# batch 2
tblone_weighted_b2<- print(svyCreateTableOne(vars = var, data = survey.design.tbl1.b2), missing=TRUE, varLabels = TRUE,digits =2,pDigits=2)
write.csv(tblone_weighted_b2, file = "output/table1/table1_weighted_b2.csv",row.names = T)
tblone_nonweighted_b2<- print(CreateTableOne(vars = var, data = pheno_b2), missing=TRUE, varLabels = TRUE,digits =2,pDigits=2)
write.csv(tblone_nonweighted_b2, file = "output/table1/table1_nonweighted_b2.csv",row.names = T)
tblone_gender_weighted_b2<- print(svyCreateTableOne(vars = var, strata = "GENDER", data = survey.design.tbl1.b2), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(tblone_gender_weighted_b2, file = "output/table1/table1_gender_weighted_b2.csv",row.names = T)
tblone_gender_nonweighted_b2<- print(CreateTableOne(vars = var, strata = "GENDER", data = pheno_b2), varLabels = TRUE,digits = 2,pDigits=2)
write.csv(tblone_gender_nonweighted_b2, file = "output/table1/table1_gender_nonweighted_b2.csv",row.names = T)

# Generate the baseline demographic characteristics tables for the top and bottom 10% SDB PCs (Table 2)
# SDB PC top/bottom 10% comparison
tblone_sdbpc1_top10_weighted<- print(svyCreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc1_top10", data = subset(survey.design.total,!is.na(sdb_pc1))), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc1_top10_weighted, file = "output/table2/table2_sdbpc1_top10_weighted.csv",row.names = T)
tblone_sdbpc1_top10_nonweighted<- print(CreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc1_top10", data = pheno_tbl1), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc1_top10_nonweighted, file = "output/table2/table2_sdbpc1_top10_nonweighted.csv",row.names = T)

tblone_sdbpc1_bottom10_weighted<- print(svyCreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc1_bottom10", data = subset(survey.design.total,!is.na(sdb_pc1))), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc1_bottom10_weighted, file = "output/table2/table2_sdbpc1_bottom10_weighted.csv",row.names = T)
tblone_sdbpc1_bottom10_nonweighted<- print(CreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc1_bottom10", data = pheno_tbl1), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc1_bottom10_nonweighted, file = "output/table2/table2_sdbpc1_bottom10_nonweighted.csv",row.names = T)

tblone_sdbpc2_top10_weighted<- print(svyCreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc2_top10", data = subset(survey.design.total,!is.na(sdb_pc2))), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc2_top10_weighted, file = "output/table2/table2_sdbpc2_top10_weighted.csv",row.names = T)
tblone_sdbpc2_top10_nonweighted<- print(CreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc2_top10", data = pheno_tbl1), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc2_top10_nonweighted, file = "output/table2/table2_sdbpc2_top10_nonweighted.csv",row.names = T)

tblone_sdbpc2_bottom10_weighted<- print(svyCreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc2_bottom10", data = subset(survey.design.total,!is.na(sdb_pc2))), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc1_bottom10_weighted, file = "output/table2/table2_sdbpc2_bottom10_weighted.csv",row.names = T)
tblone_sdbpc2_bottom10_nonweighted<- print(CreateTableOne(vars = c(var,disturb_var,hr_var), strata = "sdbpc2_bottom10", data = pheno_tbl1), missing=TRUE, varLabels = TRUE,digits = 3,pDigits=3)
write.csv(tblone_sdbpc2_bottom10_nonweighted, file = "output/table2/table2_sdbpc2_bottom10_nonweighted.csv",row.names = T)

merge_table1(tbl_nonweight=tblone_sdbpc1_top10_nonweighted,tbl_weight=tblone_sdbpc1_top10_weighted,out_path="output/table2/table2_sdbpc1_top10_merged.csv") 
merge_table1(tbl_nonweight=tblone_sdbpc2_top10_nonweighted,tbl_weight=tblone_sdbpc2_top10_weighted,out_path="output/table2/table2_sdbpc2_top10_merged.csv") 
merge_table1(tbl_nonweight=tblone_sdbpc1_bottom10_nonweighted,tbl_weight=tblone_sdbpc1_bottom10_weighted,out_path="output/table2/table2_sdbpc1_bottom10_merged.csv") 
merge_table1(tbl_nonweight=tblone_sdbpc2_bottom10_nonweighted,tbl_weight=tblone_sdbpc2_bottom10_weighted,out_path="output/table2/table2_sdbpc2_bottom10_merged.csv")

# export the processed metabolomics datasets
save(list=c("imp_only_data","imp_bin_data","imp_only_data_b2","imp_bin_data_b2"), file="data/merged_data.RData",version=2) 
