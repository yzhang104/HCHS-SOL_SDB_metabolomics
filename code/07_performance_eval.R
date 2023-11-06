#####################################
# Performance evaluation
# Generate Table 3, Figure 5, Supplemental Table s11, S12, S13, S14, Supplemental Figure S4
source("code/00_functions.R")
library_list<-c("tidyverse","labelled","stringr","plyr","dplyr","ggplot2","survey","purrr","RColorBrewer","ggrepel","gridExtra","ggpubr","ComplexHeatmap","circlize","pROC")
lapply(library_list, require, character.only = TRUE)
# System parameters
options(scipen = 999, stringsAsFactors = F)
options(survey.lonely.psu="remove")

pheno_pca<-readRDS("data/phenotype_pca.rds")
all_score_combined<-read.csv(file="data/sdb_mrs_scores.csv")
# read in time elapse between v1v2 and age calculated at visit 2
add_var<-read.csv("data/mcimet_covariates_20200716.csv",stringsAsFactors = F)%>%
  dplyr::select(ID,YRS_BTWN_V1V2,age_baseline_from_v2)
pheno_pca<-merge(pheno_pca,add_var,by="ID",all.x=T)%>%
  dplyr::mutate(AGE=age_baseline_from_v2,
                AGE_RAW=age_baseline_from_v2,
                time_y1y2=YRS_BTWN_V1V2)

# combine with phenotype data
pheno_all_score_combined<-merge(pheno_pca,all_score_combined,by="ID",all.x = T)%>%
  dplyr::mutate(HB=hypoxicburden_harmonized,
                REIgt15=slpa54gt15,
                REI3=SLPA54)
pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"sex_sdbpc1_mrs"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"female_sdbpc1_mrs"]
pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"sex_sdbpc2_mrs"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"female_sdbpc2_mrs"]
pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"sex_sdbpc1_mrs"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"male_sdbpc1_mrs"]
pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"sex_sdbpc2_mrs"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"male_sdbpc2_mrs"]
# log transform and z score HB
pheno_all_score_combined$mean_HB<-mean(log(pheno_all_score_combined$HB),na.rm =T)
pheno_all_score_combined$sd_HB<-sd(log(pheno_all_score_combined$HB),na.rm =T)
pheno_all_score_combined$log_HB<-(log(pheno_all_score_combined$HB)-pheno_all_score_combined$mean_HB)/pheno_all_score_combined$sd_HB

# standardize HB and REI
pheno_all_score_combined$BMI_RAW<-pheno_all_score_combined$BMI

pheno_all_score_combined[,c("REI3","HB","AGE","BMI")]<-pheno_all_score_combined[,c("REI3","HB","AGE","BMI")]%>%scale()
# create quartiles within combined batch
pheno_all_score_combined$sdbpc1_mrs_qt<-factor(ntile(pheno_all_score_combined$sdbpc1_mrs,4))
pheno_all_score_combined$sdbpc2_mrs_qt<-factor(ntile(pheno_all_score_combined$sdbpc2_mrs,4))
pheno_all_score_combined$osa_lasso_mrs_qt<-factor(ntile(pheno_all_score_combined$osa_lasso_mrs,4))
pheno_all_score_combined$female_sdbpc1_mrs_qt<-factor(ntile(pheno_all_score_combined$female_sdbpc1_mrs,4))
pheno_all_score_combined$female_sdbpc2_mrs_qt<-factor(ntile(pheno_all_score_combined$female_sdbpc2_mrs,4))
pheno_all_score_combined$male_sdbpc1_mrs_qt<-factor(ntile(pheno_all_score_combined$male_sdbpc1_mrs,4))
pheno_all_score_combined$male_sdbpc2_mrs_qt<-factor(ntile(pheno_all_score_combined$male_sdbpc2_mrs,4))

pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"sex_sdbpc1_mrs_qt"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"female_sdbpc1_mrs_qt"]
pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"sex_sdbpc2_mrs_qt"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Female"),"female_sdbpc2_mrs_qt"]
pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"sex_sdbpc1_mrs_qt"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"male_sdbpc1_mrs_qt"]
pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"sex_sdbpc2_mrs_qt"]<-pheno_all_score_combined[which(pheno_all_score_combined$GENDER=="Male"),"male_sdbpc2_mrs_qt"]

# association analysis: scores vs sdb pcs and incident outcomes within each batch
physio_biomarker_var<-c("REI3","log_HB","sdb_pc1","sdb_pc2")
all_biomarker_var<-c("REI3","log_HB","sdb_pc1","sdb_pc2","sdbpc1_mrs","sex_sdbpc1_mrs","sex_sdbpc2_mrs","sdbpc2_mrs","osa_lasso_mrs")
metab_biomarker_qt_var<-c("REI3","HB","sdb_pc1","sdb_pc2","sdbpc1_mrs","sex_sdbpc1_mrs","sex_sdbpc2_mrs","sdbpc2_mrs","osa_lasso_mrs")
all_biomarker_qt_var<-c("sdbpc1_mrs_qt","sdbpc2_mrs_qt","sex_sdbpc1_mrs_qt","sex_sdbpc2_mrs_qt","osa_lasso_mrs_qt")

covariate_both_list<-list(
  mdl1_both_covar=c("AGE","BMI","GENDER","CENTER","background"),
  mdl2_both_covar=c("AGE","BMI","GENDER","CENTER","background","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010")
)
covariate_gender_list<-list(
  mdl1_gender_covar=c("AGE","BMI","CENTER","background"),
  mdl2_gender_covar=c("AGE","BMI","CENTER","background","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010")
)

# Define the survey design
# survey design for SDB PCs as outcomes
sdbpc_survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_FINAL_NORM_OVERALL, data=subset(pheno_all_score_combined,!is.na(WEIGHT_FINAL_NORM_OVERALL)&!is.na(sdb_pc1)))
# survey design for incident outcomes
sdbpc_survey_design_v2=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_V2, data=subset(pheno_all_score_combined,!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(sdb_pc1)))
# survey design for incident htn outcome (drop obs with baseline htn)
sdbpc_survey_design_htn<-subset(sdbpc_survey_design_v2,baseline_htn==0&!is.na(incident_htn))
# survey design for incident dm outcome (drop obs with baseline dm)
sdbpc_survey_design_dm<-subset(sdbpc_survey_design_v2,baseline_dm==0&!is.na(incident_dm))
# survey design for incident dm outcome with baseline of normal glucose control (only include obs with non-diabetic, drop pre-diabetic obs at baseline)
sdbpc_survey_design_normvsdm<-subset(sdbpc_survey_design_v2,!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm))
# survey design for incident dm outcome with baseline of pre-dm (only include obs with pre-diabetic obs and drop obs with non-dm or missing at baseline)
sdbpc_survey_design_prevsdm<-subset(sdbpc_survey_design_v2,!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm))
# survey design for incident htn outcome among females
sdbpc_survey_design_htn_female<-subset(sdbpc_survey_design_v2,baseline_htn==0&!is.na(incident_htn)&GENDER=="Female")
# survey design for incident htn outcome among males
sdbpc_survey_design_htn_male<-subset(sdbpc_survey_design_v2,baseline_htn==0&!is.na(incident_htn)&GENDER=="Male")
# survey design for incident dm outcome among females
sdbpc_survey_design_dm_female<-subset(sdbpc_survey_design_v2,baseline_dm==0&!is.na(incident_dm)&GENDER=="Female")
# survey design for incident htn outcome among males
sdbpc_survey_design_dm_male<-subset(sdbpc_survey_design_v2,baseline_dm==0&!is.na(incident_dm)&GENDER=="Male")
# survey design for incident dm outcome among females with non-dm as control
sdbpc_survey_design_normvsdm_female<-subset(sdbpc_survey_design_v2,!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)&GENDER=="Female")
# survey design for incident dm outcome among females with pre-dm as control
sdbpc_survey_design_prevsdm_female<-subset(sdbpc_survey_design_v2,!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)&GENDER=="Female")
# survey design for incident htn outcome among males with non-dm as control
sdbpc_survey_design_normvsdm_male<-subset(sdbpc_survey_design_v2,!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)&GENDER=="Male")
# survey design for incident htn outcome among males with pre-dm as control
sdbpc_survey_design_prevsdm_male<-subset(sdbpc_survey_design_v2,!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)&GENDER=="Male")

# survey design for incident outcomes and MRS
metab_survey_design_v2=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT_NORM_OVERALL_V2, data=subset(pheno_all_score_combined,!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(sdbpc1_mrs)))
# survey design for incident htn outcome (drop obs with baseline htn)
metab_survey_design_htn<-subset(metab_survey_design_v2,baseline_htn==0&!is.na(incident_htn))
# survey design for incident dm outcome (drop obs with baseline dm)
metab_survey_design_dm<-subset(metab_survey_design_v2,baseline_dm==0&!is.na(incident_dm))
# survey design for incident dm outcome with baseline of normal glucose control (only include obs with non-diabetic, drop pre-diabetic obs at baseline)
metab_survey_design_normvsdm<-subset(metab_survey_design_v2,!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm))
# survey design for incident dm outcome with baseline of pre-dm (only include obs with pre-diabetic obs and drop obs with non-dm or missing at baseline)
metab_survey_design_prevsdm<-subset(metab_survey_design_v2,!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm))
# survey design for incident htn outcome among females
metab_survey_design_htn_female<-subset(metab_survey_design_v2,baseline_htn==0&!is.na(incident_htn)&GENDER=="Female")
# survey design for incident htn outcome among males
metab_survey_design_htn_male<-subset(metab_survey_design_v2,baseline_htn==0&!is.na(incident_htn)&GENDER=="Male")
# survey design for incident dm outcome among females
metab_survey_design_dm_female<-subset(metab_survey_design_v2,baseline_dm==0&!is.na(incident_dm)&GENDER=="Female")
# survey design for incident htn outcome among males
metab_survey_design_dm_male<-subset(metab_survey_design_v2,baseline_dm==0&!is.na(incident_dm)&GENDER=="Male")
# survey design for incident dm outcome among females with non-dm as control
metab_survey_design_normvsdm_female<-subset(metab_survey_design_v2,!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)&GENDER=="Female")
# survey design for incident dm outcome among females with pre-dm as control
metab_survey_design_prevsdm_female<-subset(metab_survey_design_v2,!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)&GENDER=="Female")
# survey design for incident htn outcome among males with non-dm as control
metab_survey_design_normvsdm_male<-subset(metab_survey_design_v2,!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)&GENDER=="Male")
# survey design for incident htn outcome among males with pre-dm as control
metab_survey_design_prevsdm_male<-subset(metab_survey_design_v2,!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)&GENDER=="Male")

# create index for the model list for all biomarkers loop in both gender combined
metab_both_mdl_index<-as.data.frame(matrix(1:(length(all_biomarker_var)*length(covariate_both_list)),ncol=length(all_biomarker_var),nrow=length(covariate_both_list)))
row.names(metab_both_mdl_index)<-names(covariate_both_list)
colnames(metab_both_mdl_index)<-all_biomarker_var

# create index for the model list for all physiological biomarkers loop in both gender combined
physio_both_mdl_index<-as.data.frame(matrix(1:(length(physio_biomarker_var)*length(covariate_both_list)),ncol=length(physio_biomarker_var),nrow=length(covariate_both_list)))
row.names(physio_both_mdl_index)<-names(covariate_both_list)
colnames(physio_both_mdl_index)<-physio_biomarker_var

# create index for the model list for all biomarkers loop in stratified gender
metab_gender_mdl_index<-as.data.frame(matrix(1:(length(all_biomarker_var)*length(covariate_gender_list)),ncol=length(all_biomarker_var),nrow=length(covariate_gender_list)))
row.names(metab_gender_mdl_index)<-names(covariate_gender_list)
colnames(metab_gender_mdl_index)<-all_biomarker_var

# create index for the model list for all biomarker(except for sdb pc1 and pc2) loop in both gender combined
metab_nosdbpc_both_mdl_index<-as.data.frame(matrix(1:(length(all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")])*length(covariate_both_list)),ncol=length(all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")]),nrow=length(covariate_both_list)))
row.names(metab_nosdbpc_both_mdl_index)<-names(covariate_both_list)
colnames(metab_nosdbpc_both_mdl_index)<-all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")]

# create index for the model list for all biomarker loop in both gender combined in quartiles
metab_both_qt_mdl_index<-as.data.frame(matrix(1:(length(all_biomarker_qt_var)*length(covariate_both_list)),ncol=length(all_biomarker_qt_var),nrow=length(covariate_both_list)))
row.names(metab_both_qt_mdl_index)<-names(covariate_both_list)
colnames(metab_both_qt_mdl_index)<-all_biomarker_qt_var

# create index for the model list for stratified gender in quartiles
metab_gender_qt_mdl_index<-as.data.frame(matrix(1:(length(all_biomarker_qt_var)*length(covariate_gender_list)),ncol=length(all_biomarker_qt_var),nrow=length(covariate_gender_list)))
row.names(metab_gender_qt_mdl_index)<-names(covariate_gender_list)
colnames(metab_gender_qt_mdl_index)<-all_biomarker_qt_var

# loop association for SDB PCs~MRS in batch 1
b1_sdbpc1_both_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_both_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&batch=="b1"),outcome_var="sdb_pc1")
b1_sdbpc1_female_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Female"&batch=="b1"),outcome_var="sdb_pc1")
b1_sdbpc1_male_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Male"&batch=="b1"),outcome_var="sdb_pc1")
b1_sdbpc2_both_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_both_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&batch=="b1"),outcome_var="sdb_pc2")
b1_sdbpc2_female_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Female"&batch=="b1"),outcome_var="sdb_pc2")
b1_sdbpc2_male_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Male"&batch=="b1"),outcome_var="sdb_pc2")

# extract output
b1_sdbpc1_both_output<-extract_mdl_sig_output(model_output_list=b1_sdbpc1_both_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc1",model_strata="both")
b1_sdbpc2_both_output<-extract_mdl_sig_output(model_output_list=b1_sdbpc2_both_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc2",model_strata="both")
b1_sdbpc1_female_output<-extract_mdl_sig_output(model_output_list=b1_sdbpc1_female_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc1",model_strata="female")
b1_sdbpc2_female_output<-extract_mdl_sig_output(model_output_list=b1_sdbpc2_female_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc2",model_strata="female")
b1_sdbpc1_male_output<-extract_mdl_sig_output(model_output_list=b1_sdbpc1_male_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc1",model_strata="male")
b1_sdbpc2_male_output<-extract_mdl_sig_output(model_output_list=b1_sdbpc2_male_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc2",model_strata="male")
b1_sdbpc_output_all<-rbind(b1_sdbpc1_both_output,b1_sdbpc2_both_output,b1_sdbpc1_female_output,b1_sdbpc2_female_output,b1_sdbpc1_male_output,b1_sdbpc2_male_output)%>%
  dplyr::select(-or,-lowerCI,-upperCI)%>%
  dplyr::mutate(batch="batch_1")

b1_sdbpc_output<-b1_sdbpc_output_all%>%
  dplyr::filter(!trait%in%c("HB","REI3","log_HB"))

##################
# loop association for SDB PCs~mrsx in batch 2
b2_sdbpc1_both_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_both_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&batch=="b2"),outcome_var="sdb_pc1")
b2_sdbpc1_female_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Female"&batch=="b2"),outcome_var="sdb_pc1")
b2_sdbpc1_male_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Male"&batch=="b2"),outcome_var="sdb_pc1")
b2_sdbpc2_both_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_both_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&batch=="b2"),outcome_var="sdb_pc2")
b2_sdbpc2_female_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Female"&batch=="b2"),outcome_var="sdb_pc2")
b2_sdbpc2_male_mdl_list<-loop_association_analysis(exposure_list=all_biomarker_var[!all_biomarker_var%in%c("sdb_pc1","sdb_pc2")],covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design,!is.na(sdbpc1_mrs)&GENDER=="Male"&batch=="b2"),outcome_var="sdb_pc2")
# extract output
b2_sdbpc1_both_output<-extract_mdl_sig_output(model_output_list=b2_sdbpc1_both_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc1",model_strata="both")
b2_sdbpc2_both_output<-extract_mdl_sig_output(model_output_list=b2_sdbpc2_both_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc2",model_strata="both")
b2_sdbpc1_female_output<-extract_mdl_sig_output(model_output_list=b2_sdbpc1_female_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc1",model_strata="female")
b2_sdbpc2_female_output<-extract_mdl_sig_output(model_output_list=b2_sdbpc2_female_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc2",model_strata="female")
b2_sdbpc1_male_output<-extract_mdl_sig_output(model_output_list=b2_sdbpc1_male_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc1",model_strata="male")
b2_sdbpc2_male_output<-extract_mdl_sig_output(model_output_list=b2_sdbpc2_male_mdl_list,model_index_matrix=metab_nosdbpc_both_mdl_index,model_outcome="sdb_pc2",model_strata="male")
b2_sdbpc_output_all<-rbind(b2_sdbpc1_both_output,b2_sdbpc2_both_output,b2_sdbpc1_female_output,b2_sdbpc2_female_output,b2_sdbpc1_male_output,b2_sdbpc2_male_output)%>%
  select(-or,-lowerCI,-upperCI)%>%
  dplyr::mutate(batch="batch_2")

b2_sdbpc_output<-b2_sdbpc_output_all%>%
  dplyr::filter(!trait%in%c("HB","REI3","log_HB"))

# loop association for incident outcomes~physiological biomarkers in combined batch (b1+b2) dataset
# incident htn
complete_inchtn_both_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_both_list,survey_design=sdbpc_survey_design_htn,outcome_var="incident_htn",offset_var="time_y1y2")
complete_inchtn_female_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_htn,GENDER=="Female"),outcome_var="incident_htn",offset_var="time_y1y2")
complete_inchtn_male_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_htn,GENDER=="Male"),outcome_var="incident_htn",offset_var="time_y1y2")
# combined control incident dm
complete_incdm_both_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_both_list,survey_design=sdbpc_survey_design_dm,outcome_var="incident_dm",offset_var="time_y1y2")
complete_incdm_female_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_dm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
complete_incdm_male_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_dm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")
# non-dm vs dm
complete_normvsdm_both_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_both_list,survey_design=sdbpc_survey_design_normvsdm,outcome_var="incident_dm",offset_var="time_y1y2")
complete_normvsdm_female_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_normvsdm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
complete_normvsdm_male_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_normvsdm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")
# pre-dm vs dm
complete_prevsdm_both_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_both_list,survey_design=sdbpc_survey_design_prevsdm,outcome_var="incident_dm",offset_var="time_y1y2")
complete_prevsdm_female_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_prevsdm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
complete_prevsdm_male_mdl_list<-loop_association_analysis_offset(exposure_list=physio_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(sdbpc_survey_design_prevsdm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")

# extract output
complete_inchtn_both_output<-extract_mdl_sig_output(model_output_list=complete_inchtn_both_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="incident_htn",model_strata="both")
complete_inchtn_female_output<-extract_mdl_sig_output(model_output_list=complete_inchtn_female_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="incident_htn",model_strata="female")
complete_inchtn_male_output<-extract_mdl_sig_output(model_output_list=complete_inchtn_male_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="incident_htn",model_strata="male")
complete_incdm_both_output<-extract_mdl_sig_output(model_output_list=complete_incdm_both_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="incident_dm",model_strata="both")
complete_incdm_female_output<-extract_mdl_sig_output(model_output_list=complete_incdm_female_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="incident_dm",model_strata="female")
complete_incdm_male_output<-extract_mdl_sig_output(model_output_list=complete_incdm_male_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="incident_dm",model_strata="male")
complete_normvsdm_both_output<-extract_mdl_sig_output(model_output_list=complete_normvsdm_both_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="norm_dm",model_strata="both")
complete_normvsdm_female_output<-extract_mdl_sig_output(model_output_list=complete_normvsdm_female_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="norm_dm",model_strata="female")
complete_normvsdm_male_output<-extract_mdl_sig_output(model_output_list=complete_normvsdm_male_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="norm_dm",model_strata="male")
complete_prevsdm_both_output<-extract_mdl_sig_output(model_output_list=complete_prevsdm_both_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="predm_dm",model_strata="both")
complete_prevsdm_female_output<-extract_mdl_sig_output(model_output_list=complete_prevsdm_female_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="predm_dm",model_strata="female")
complete_prevsdm_male_output<-extract_mdl_sig_output(model_output_list=complete_prevsdm_male_mdl_list,model_index_matrix=physio_both_mdl_index,model_outcome="predm_dm",model_strata="male")

complete_compare_dm_output<-rbind(complete_prevsdm_both_output,complete_prevsdm_female_output,complete_prevsdm_male_output,complete_normvsdm_both_output,complete_normvsdm_female_output,complete_normvsdm_male_output)%>%
  select(-beta,-lower95,-upper95)%>%
  dplyr::mutate(IRR=or)


complete_inc_output_all<-rbind(complete_inchtn_both_output,complete_inchtn_female_output,complete_inchtn_male_output,complete_incdm_both_output,complete_incdm_female_output,complete_incdm_male_output)%>%
  select(-beta,-lower95,-upper95)%>%
  dplyr::mutate(IRR=or)

complete_inc_output<-complete_inc_output_all

# loop association for incident outcomes~metabolomic biomarkers in combined batch
# incident htn
combined_inchtn_both_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_htn,outcome_var="incident_htn",offset_var="time_y1y2")
combined_inchtn_female_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_htn,GENDER=="Female"),outcome_var="incident_htn",offset_var="time_y1y2")
combined_inchtn_male_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_htn,GENDER=="Male"),outcome_var="incident_htn",offset_var="time_y1y2")
# incident dm
# combined control
combined_incdm_both_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_dm,outcome_var="incident_dm",offset_var="time_y1y2")
combined_incdm_female_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_dm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
combined_incdm_male_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_dm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")
# non-dm vs dm
combined_normvsdm_both_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_normvsdm,outcome_var="incident_dm",offset_var="time_y1y2")
combined_normvsdm_female_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_normvsdm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
combined_normvsdm_male_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_normvsdm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")
# pre-dm vs dm
combined_prevsdm_both_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_prevsdm,outcome_var="incident_dm",offset_var="time_y1y2")
combined_prevsdm_female_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_prevsdm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
combined_prevsdm_male_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_prevsdm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")

# extract output
combined_inchtn_both_output<-extract_mdl_sig_output(model_output_list=combined_inchtn_both_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="incident_htn",model_strata="both")
combined_inchtn_female_output<-extract_mdl_sig_output(model_output_list=combined_inchtn_female_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="incident_htn",model_strata="female")
combined_inchtn_male_output<-extract_mdl_sig_output(model_output_list=combined_inchtn_male_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="incident_htn",model_strata="male")
combined_incdm_both_output<-extract_mdl_sig_output(model_output_list=combined_incdm_both_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="incident_dm",model_strata="both")
combined_incdm_female_output<-extract_mdl_sig_output(model_output_list=combined_incdm_female_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="incident_dm",model_strata="female")
combined_incdm_male_output<-extract_mdl_sig_output(model_output_list=combined_incdm_male_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="incident_dm",model_strata="male")

combined_normvsdm_both_output<-extract_mdl_sig_output(model_output_list=combined_normvsdm_both_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="norm_dm",model_strata="both")
combined_normvsdm_female_output<-extract_mdl_sig_output(model_output_list=combined_normvsdm_female_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="norm_dm",model_strata="female")
combined_normvsdm_male_output<-extract_mdl_sig_output(model_output_list=combined_normvsdm_male_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="norm_dm",model_strata="male")
combined_prevsdm_both_output<-extract_mdl_sig_output(model_output_list=combined_prevsdm_both_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="predm_dm",model_strata="both")
combined_prevsdm_female_output<-extract_mdl_sig_output(model_output_list=combined_prevsdm_female_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="predm_dm",model_strata="female")
combined_prevsdm_male_output<-extract_mdl_sig_output(model_output_list=combined_prevsdm_male_mdl_list,model_index_matrix=metab_both_mdl_index,model_outcome="predm_dm",model_strata="male")

combined_compare_dm_output<-rbind(combined_prevsdm_both_output,combined_prevsdm_female_output,combined_prevsdm_male_output,combined_normvsdm_both_output,combined_normvsdm_female_output,combined_normvsdm_male_output)%>%
  select(-beta,-lower95,-upper95)%>%
  dplyr::mutate(IRR=or)

combined_inc_output_all<-rbind(combined_inchtn_both_output,combined_inchtn_female_output,combined_inchtn_male_output,combined_incdm_both_output,combined_incdm_female_output,combined_incdm_male_output)%>%
  select(-beta,-lower95,-upper95)%>%
  dplyr::mutate(IRR=or)
combined_inc_output<-combined_inc_output_all

# loop association for incident outcomes~metabolomic biomarker in quartiles in combined batch
combined_inchtn_both_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_htn,outcome_var="incident_htn",offset_var="time_y1y2")
combined_inchtn_female_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_htn,GENDER=="Female"),outcome_var="incident_htn",offset_var="time_y1y2")
combined_inchtn_male_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_htn,GENDER=="Male"),outcome_var="incident_htn",offset_var="time_y1y2")
combined_incdm_both_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_dm,outcome_var="incident_dm",offset_var="time_y1y2")
combined_incdm_female_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_dm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
combined_incdm_male_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_dm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")
# normal vs dm
combined_normvsdm_both_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_normvsdm,outcome_var="incident_dm",offset_var="time_y1y2")
combined_normvsdm_female_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_normvsdm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
combined_normvsdm_male_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_normvsdm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")
# pre-dm vs dm
combined_prevsdm_both_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_both_list,survey_design=metab_survey_design_prevsdm,outcome_var="incident_dm",offset_var="time_y1y2")
combined_prevsdm_female_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_prevsdm,GENDER=="Female"),outcome_var="incident_dm",offset_var="time_y1y2")
combined_prevsdm_male_qt_mdl_list<-loop_association_analysis_offset(exposure_list=all_biomarker_qt_var,covariate_list=covariate_gender_list,survey_design=subset(metab_survey_design_prevsdm,GENDER=="Male"),outcome_var="incident_dm",offset_var="time_y1y2")


# extract output value from quartiles model
inchtn_qt_both_output<-quartile_association_models_output_offset(model_list=combined_inchtn_both_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="both",
                outcome="inc_htn",
                baseline="normal")
incdm_qt_both_output<-quartile_association_models_output_offset(model_list=combined_incdm_both_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="both",
                outcome="inc_dm",
                baseline="norm+preDM")
inchtn_qt_female_output<-quartile_association_models_output_offset(model_list=combined_inchtn_female_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="female",
                outcome="inc_htn",
                baseline="normal")
incdm_qt_female_output<-quartile_association_models_output_offset(model_list=combined_incdm_female_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="female",
                outcome="inc_dm",
                baseline="norm+preDM")
inchtn_qt_male_output<-quartile_association_models_output_offset(model_list=combined_inchtn_male_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="male",
                outcome="inc_htn",
                baseline="normal")
incdm_qt_male_output<-quartile_association_models_output_offset(model_list=combined_incdm_male_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="male",
                outcome="inc_dm",
                baseline="norm+preDM")
normvsdm_qt_both_output<-quartile_association_models_output_offset(model_list=combined_normvsdm_both_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="both",
                outcome="inc_dm",
                baseline="norm")
normvsdm_qt_female_output<-quartile_association_models_output_offset(model_list=combined_normvsdm_female_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="female",
                outcome="inc_dm",
                baseline="norm")
normvsdm_qt_male_output<-quartile_association_models_output_offset(model_list=combined_normvsdm_male_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="male",
                outcome="inc_dm",
                baseline="norm")
prevsdm_qt_both_output<-quartile_association_models_output_offset(model_list=combined_prevsdm_both_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="both",
                outcome="inc_dm",
                baseline="preDM")
prevsdm_qt_female_output<-quartile_association_models_output_offset(model_list=combined_prevsdm_female_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="female",
                outcome="inc_dm",
                baseline="preDM")
prevsdm_qt_male_output<-quartile_association_models_output_offset(model_list=combined_prevsdm_male_qt_mdl_list,dataset="combined_batch")%>%
  dplyr::mutate(strata="male",
                outcome="inc_dm",
                baseline="preDM")
qt_output<-rbind(inchtn_qt_both_output,incdm_qt_both_output,inchtn_qt_female_output,incdm_qt_female_output,inchtn_qt_male_output,incdm_qt_male_output,normvsdm_qt_both_output,normvsdm_qt_female_output,normvsdm_qt_male_output,prevsdm_qt_both_output,prevsdm_qt_female_output,prevsdm_qt_male_output)

sdbpc_association<-rbind(b1_sdbpc_output_all,b2_sdbpc_output_all)
write.csv(sdbpc_association,file="output/table3.csv",row.names = F)
complete_inc_output_all<-complete_inc_output_all%>%
  dplyr::mutate(baseline=dplyr::case_when(
    outcome=="incident_htn" ~ "Normal",
    outcome=="incident_dm" ~ "Normal+Pre-diabetic"
  ))
complete_compare_dm_output<-complete_compare_dm_output%>%
  dplyr::mutate(baseline=dplyr::case_when(
    outcome=="predm_dm" ~ "Pre-diabetic",
    outcome=="norm_dm" ~ "Normal"
  ))
complete_compare_dm_output$outcome<-"incident_dm"
complete_output_all<-rbind(complete_inc_output_all,complete_compare_dm_output)
write.csv(complete_output_all,file="output/suppl_table_s11.csv",row.names = F)

combined_inc_output_all<-combined_inc_output_all%>%
  dplyr::mutate(baseline=dplyr::case_when(
    outcome=="incident_htn" ~ "Normal",
    outcome=="incident_dm" ~ "Normal+Pre-diabetic"
  ))
combined_compare_dm_output<-combined_compare_dm_output%>%
  dplyr::mutate(baseline=dplyr::case_when(
    outcome=="predm_dm" ~ "Pre-diabetic",
    outcome=="norm_dm" ~ "Normal"
  ))
combined_compare_dm_output$outcome<-"incident_dm"
combined_output_all<-rbind(combined_inc_output_all,combined_compare_dm_output)
write.csv(combined_inc_output_all,file="output/suppl_table_s12.csv",row.names = F)

write.csv(qt_output,file="output/suppl_table_s14.csv",row.names = F)

# Generate Figure 5
combined_inc_output_select<-combined_inc_output%>%
  dplyr::filter(!(trait=="sex_sdbpc1_mrs"&strata=="both"))%>%
  dplyr::filter(!(trait=="sex_sdbpc2_mrs"&strata=="both"))%>%
  dplyr::filter(trait%in%c("sdbpc1_mrs","sdbpc2_mrs","osa_lasso_mrs","log_HB","REI3"))
fig5<-ggplot(combined_inc_output_select, aes(trait, IRR, color = model, group = model)) +
  geom_errorbar(aes(min = lowerCI, max = upperCI), size = 1, 
                width = 0.25, position = position_dodge()) +
  geom_point(position = position_dodge(width = 0.5), color = "black") +
  facet_grid(row=vars(outcome),col=vars(strata))+
  geom_text(aes(label = p_label, color = model, group = model),position = position_dodge(0.9),vjust=0.5,hjust=-0.2)+
  geom_hline(yintercept=1, linetype="dashed",color = "red", size=0.5)+
  theme_bw()+ scale_color_brewer(palette="Dark2")+ scale_y_continuous(trans='log2')+ 
  theme(legend.position="top")+
  coord_flip()+
  labs(title = "Association between SDB phenotypes and incident cardiometabolic outcomes in the combined batch",
       caption = "*: p value is between 0.01 and 0.05; **: p value is between 0.001 and 0.01; ***: p value is below 0.001\nmodel_1:age,center,gender,latino background, BMI; model_2:model_1 covariates,cigarette use, alcohol use, physical activity, diet")
# export the association estimates plot
jpeg("output/fig5.jpg")
print(fig5)
dev.off()

###########
# Generate Supplemental Figure S4
qt_output_full<-qt_output
qt_output<-qt_output[which(qt_output$Exposure%in%c("osa_lasso_mrs_qt","sdbpc1_mrs_qt","sdbpc2_mrs_qt")),]
qt_output$outcome<-gsub("inc_","incident ",qt_output$outcome)
qt_output$outcome<-gsub("htn","hypertension",qt_output$outcome)
qt_output$outcome<-gsub("dm","diabetes",qt_output$outcome)
# qt_output$Exposure<-gsub("_MetInd_qt","_mrs",qt_output$Exposure)
qt_output$Exposure<-gsub("sdbpc1_mrs_qt","SDB PC1 MRS",qt_output$Exposure)
qt_output$Exposure<-gsub("sdbpc2_mrs_qt","SDB PC2 MRS",qt_output$Exposure)
qt_output$Exposure<-gsub("osa_lasso_mrs_qt","OSA LASSO MRS",qt_output$Exposure)
qt_output$IRR<-qt_output$OR
qt_output<-qt_output%>%
  dplyr::mutate(category=paste0(baseline,"-",outcome))
qt_plot<-ggplot(qt_output[which(qt_output$model%in%c("model_2")&qt_output$Exposure%in%c("SDB PC1 MRS","SDB PC2 MRS","OSA LASSO MRS")&qt_output$category%in%c("normal-incident hypertension","norm+preDM-incident diabetes")),], aes(Exposure, IRR, color = fct_rev(Quartile), group = fct_rev(Quartile))) +
  geom_errorbar(aes(min = lowerCI, max = upperCI), size = 0.75, 
                width = 0.25, position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8), color = "black")+
  geom_hline(yintercept=1, linetype="dashed",color = "red", size=0.5)+
  geom_text(aes(label = p_label, color = fct_rev(Quartile), group = fct_rev(Quartile)),position = position_dodge(0.7),hjust=-0.18,vjust=1.1,size=5)+
  # facet_wrap(.~dataset, nrow = 2) +
  facet_grid(cols=vars(strata),rows = vars(outcome))+
  theme_bw(base_size = 14)+ scale_color_brewer(palette="Dark2",limits = c("2nd", "3rd", "4th"),name = "Quartile") +
  scale_y_continuous(trans='log2')+
  coord_flip()
# export the quartile associations plot
jpeg("output/suppl_figs4.jpg", width = 800, height = 600)
print(qt_plot)
dev.off()

###########
# AUC Comparison
# Incident DM (normal+prediabetic)
# demographic predictors only model
combined_incdm_base_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+offset(log(time_y1y2)),design=metab_survey_design_dm,family=quasipoisson(link='log'))
# AHI model
combined_incdm_rei_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+REI3+offset(log(time_y1y2)),design=metab_survey_design_dm,family=quasipoisson(link='log'))
# HB model
combined_incdm_hb_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+log_HB+offset(log(time_y1y2)),design=metab_survey_design_dm,family=quasipoisson(link='log'))
# sdbpc1 met index model
combined_incdm_sdbpc1mrs_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+sdbpc1_mrs+offset(log(time_y1y2)),design=metab_survey_design_dm,family=quasipoisson(link='log'))
# sdbpc1 met index model
combined_incdm_sdbpc2mrs_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+sdbpc2_mrs+offset(log(time_y1y2)),design=metab_survey_design_dm,family=quasipoisson(link='log'))
# osa lasso met index model
combined_incdm_osalasso_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+osa_lasso_mrs+offset(log(time_y1y2)),design=metab_survey_design_dm,family=quasipoisson(link='log'))

# Incident DM (baseline: normal)
# demographic predictors only model
combined_normvsdm_base_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+offset(log(time_y1y2)),design=metab_survey_design_normvsdm,family=quasipoisson(link='log'))
# AHI model
combined_normvsdm_rei_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+REI3+offset(log(time_y1y2)),design=metab_survey_design_normvsdm,family=quasipoisson(link='log'))
# HB model
combined_normvsdm_hb_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+log_HB+offset(log(time_y1y2)),design=metab_survey_design_normvsdm,family=quasipoisson(link='log'))
# sdbpc1 met index model
combined_normvsdm_sdbpc1mrs_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+sdbpc1_mrs+offset(log(time_y1y2)),design=metab_survey_design_normvsdm,family=quasipoisson(link='log'))
# sdbpc2 met index model
combined_normvsdm_sdbpc2mrs_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+sdbpc2_mrs+offset(log(time_y1y2)),design=metab_survey_design_normvsdm,family=quasipoisson(link='log'))
# osa lasso met index model
combined_normvsdm_osalasso_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+osa_lasso_mrs+offset(log(time_y1y2)),design=metab_survey_design_normvsdm,family=quasipoisson(link='log'))

# Incident DM (baseline: prediabetic)
# demographic predictors only model
combined_prevsdm_base_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+offset(log(time_y1y2)),design=metab_survey_design_prevsdm,family=quasipoisson(link='log'))
# AHI model
combined_prevsdm_rei_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+REI3+offset(log(time_y1y2)),design=metab_survey_design_prevsdm,family=quasipoisson(link='log'))
# HB model
combined_prevsdm_hb_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+log_HB+offset(log(time_y1y2)),design=metab_survey_design_prevsdm,family=quasipoisson(link='log'))
# sdbpc1 met index model
combined_prevsdm_sdbpc1mrs_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+sdbpc1_mrs+offset(log(time_y1y2)),design=metab_survey_design_prevsdm,family=quasipoisson(link='log'))
# sdbpc1 met index model
combined_prevsdm_sdbpc2mrs_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+sdbpc2_mrs+offset(log(time_y1y2)),design=metab_survey_design_prevsdm,family=quasipoisson(link='log'))
# osa lasso met index model
combined_prevsdm_osalasso_mdl<-survey::svyglm(incident_dm~AGE+BMI+GENDER+CENTER+background+osa_lasso_mrs+offset(log(time_y1y2)),design=metab_survey_design_prevsdm,family=quasipoisson(link='log'))

# Incident HTN
# demographic predictors only model
combined_inchtn_base_mdl<-survey::svyglm(incident_htn~AGE+BMI+GENDER+CENTER+background,design=metab_survey_design_htn,family=quasipoisson(link='log'))
# AHI model
combined_inchtn_rei_mdl<-survey::svyglm(incident_htn~AGE+BMI+GENDER+CENTER+background+REI3,design=metab_survey_design_htn,family=quasipoisson(link='log'))
# HB model
combined_inchtn_hb_mdl<-survey::svyglm(incident_htn~AGE+BMI+GENDER+CENTER+background+log_HB,design=metab_survey_design_htn,family=quasipoisson(link='log'))
# sdbpc1 met index model
combined_inchtn_sdbpc1mrs_mdl<-survey::svyglm(incident_htn~AGE+BMI+GENDER+CENTER+background+sdbpc1_mrs,design=metab_survey_design_htn,family=quasipoisson(link='log'))
# sdbpc1 met index model
combined_inchtn_sdbpc2mrs_mdl<-survey::svyglm(incident_htn~AGE+BMI+GENDER+CENTER+background+sdbpc2_mrs,design=metab_survey_design_htn,family=quasipoisson(link='log'))
# osa lasso met index model
combined_inchtn_osalasso_mdl<-survey::svyglm(incident_htn~AGE+BMI+GENDER+CENTER+background+osa_lasso_mrs,design=metab_survey_design_htn,family=quasipoisson(link='log'))

model_order<-c("Demographics Only","Demographics+REI3","Demographics+HB","Demographics+OSA LASSO MRS","Demographics+SDB PC1 MRS","Demographics+SDB PC2 MRS")
# Incdident DM baseline (normal+prediabetic)
incdm_both_auc<-data.frame(model=model_order,lower95=rep(NA,6),auc=rep(NA,6),upper95=rep(NA,6))
# AUC
# base model
combined_incdm_base_pred<-predict(combined_incdm_base_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_dm!=1&!is.na(incident_dm)), type="response")
incdm_both_auc[1,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_dm!=1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_incdm_base_pred)[1:3]
# AHI model
combined_incdm_rei_pred<-predict(combined_incdm_rei_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_dm!=1&!is.na(incident_dm)), type="response")
incdm_both_auc[2,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_dm!=1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_incdm_rei_pred)[1:3]
# HB model
combined_incdm_hb_pred<-predict(combined_incdm_hb_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_dm!=1&!is.na(incident_dm)), type="response")
incdm_both_auc[3,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_dm!=1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_incdm_hb_pred)[1:3]
# OSA LASSO model
combined_incdm_osalasso_pred<-predict(combined_incdm_osalasso_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_dm!=1&!is.na(incident_dm)), type="response")
incdm_both_auc[4,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_dm!=1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_incdm_osalasso_pred)[1:3]
# SDB PC1 MRS model
combined_incdm_sdbpc1mrs_pred<-predict(combined_incdm_sdbpc1mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_dm!=1&!is.na(incident_dm)), type="response")
incdm_both_auc[5,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_dm!=1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_incdm_sdbpc1mrs_pred)[1:3]
# SDB PC2 MRS model
combined_incdm_sdbpc2mrs_pred<-predict(combined_incdm_sdbpc2mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_dm!=1&!is.na(incident_dm)), type="response")
incdm_both_auc[6,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_dm!=1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_incdm_sdbpc2mrs_pred)[1:3]

# Incident DM (baseline: normal)
normvsdm_both_auc<-data.frame(model=model_order,lower95=rep(NA,6),auc=rep(NA,6),upper95=rep(NA,6))

# AUC
# base model
combined_normvsdm_base_pred<-predict(combined_normvsdm_base_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)), type="response")
normvsdm_both_auc[1,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_normvsdm_base_pred)[1:3]
# AHI model
combined_normvsdm_rei_pred<-predict(combined_normvsdm_rei_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)), type="response")
normvsdm_both_auc[2,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_normvsdm_rei_pred)[1:3]
# HB model
combined_normvsdm_hb_pred<-predict(combined_normvsdm_hb_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)), type="response")
normvsdm_both_auc[3,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_normvsdm_hb_pred)[1:3]
# OSA LASSO model
combined_normvsdm_osalasso_pred<-predict(combined_normvsdm_osalasso_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)), type="response")
normvsdm_both_auc[4,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_normvsdm_osalasso_pred)[1:3]
# SDB PC1 MRS model
combined_normvsdm_sdbpc1mrs_pred<-predict(combined_normvsdm_sdbpc1mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)), type="response")
normvsdm_both_auc[5,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_normvsdm_sdbpc1mrs_pred)[1:3]
# SDB PC2 MRS model
combined_normvsdm_sdbpc2mrs_pred<-predict(combined_normvsdm_sdbpc2mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==1&!is.na(incident_dm)), type="response")
normvsdm_both_auc[6,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==1&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_normvsdm_sdbpc2mrs_pred)[1:3]

# Incident DM (baseline: prediabetic)
prevsdm_both_auc<-data.frame(model=model_order,lower95=rep(NA,6),auc=rep(NA,6),upper95=rep(NA,6))
# AUC
# base model
combined_prevsdm_base_pred<-predict(combined_prevsdm_base_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)), type="response")
prevsdm_both_auc[1,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==2&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_prevsdm_base_pred)[1:3]
# AHI model
combined_prevsdm_rei_pred<-predict(combined_prevsdm_rei_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)), type="response")
prevsdm_both_auc[2,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==2&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_prevsdm_rei_pred)[1:3]
# HB model
combined_prevsdm_hb_pred<-predict(combined_prevsdm_hb_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)), type="response")
prevsdm_both_auc[3,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==2&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_prevsdm_hb_pred)[1:3]
# OSA LASSO model
combined_prevsdm_osalasso_pred<-predict(combined_prevsdm_osalasso_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)), type="response")
prevsdm_both_auc[4,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==2&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_prevsdm_osalasso_pred)[1:3]
# SDB PC1 MRS model
combined_prevsdm_sdbpc1mrs_pred<-predict(combined_prevsdm_sdbpc1mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)), type="response")
prevsdm_both_auc[5,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==2&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_prevsdm_sdbpc1mrs_pred)[1:3]
# SDB PC2 MRS model
combined_prevsdm_sdbpc2mrs_pred<-predict(combined_prevsdm_sdbpc2mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&!is.na(DIABETES3)&DIABETES3==2&!is.na(incident_dm)), type="response")
prevsdm_both_auc[6,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&!is.na(pheno_all_score_combined$DIABETES3)&pheno_all_score_combined$DIABETES3==2&!is.na(pheno_all_score_combined$incident_dm)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_dm"],combined_prevsdm_sdbpc2mrs_pred)[1:3]

# Incident HTN
inchtn_both_auc<-data.frame(model=model_order,lower95=rep(NA,6),auc=rep(NA,6),upper95=rep(NA,6))
# AUC
# base model
combined_inchtn_base_pred<-predict(combined_inchtn_base_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_htn!=1&!is.na(incident_htn)), type="response")
inchtn_both_auc[1,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_htn!=1&!is.na(pheno_all_score_combined$incident_htn)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_htn"],combined_inchtn_base_pred)[1:3]
# AHI model
combined_inchtn_rei_pred<-predict(combined_inchtn_rei_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_htn!=1&!is.na(incident_htn)), type="response")
inchtn_both_auc[2,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_htn!=1&!is.na(pheno_all_score_combined$incident_htn)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_htn"],combined_inchtn_rei_pred)[1:3]
# HB model
combined_inchtn_hb_pred<-predict(combined_inchtn_hb_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_htn!=1&!is.na(incident_htn)), type="response")
inchtn_both_auc[3,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_htn!=1&!is.na(pheno_all_score_combined$incident_htn)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_htn"],combined_inchtn_hb_pred)[1:3]
# OSA LASSO model
combined_inchtn_osalasso_pred<-predict(combined_inchtn_osalasso_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_htn!=1&!is.na(incident_htn)), type="response")
inchtn_both_auc[4,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_htn!=1&!is.na(pheno_all_score_combined$incident_htn)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_htn"],combined_inchtn_osalasso_pred)[1:3]
# SDB PC1 MRS model
combined_inchtn_sdbpc1mrs_pred<-predict(combined_inchtn_sdbpc1mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_htn!=1&!is.na(incident_htn)), type="response")
inchtn_both_auc[5,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_htn!=1&!is.na(pheno_all_score_combined$incident_htn)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_htn"],combined_inchtn_sdbpc1mrs_pred)[1:3]
# SDB PC2 MRS model
combined_inchtn_sdbpc2mrs_pred<-predict(combined_inchtn_sdbpc2mrs_mdl,subset(pheno_all_score_combined,ID%in%all_score_combined$ID&!is.na(WEIGHT_NORM_OVERALL_V2)&baseline_htn!=1&!is.na(incident_htn)), type="response")
inchtn_both_auc[6,2:4]<-ci.auc(pheno_all_score_combined[which(pheno_all_score_combined$ID%in%all_score_combined$ID&pheno_all_score_combined$baseline_htn!=1&!is.na(pheno_all_score_combined$incident_htn)&!is.na(pheno_all_score_combined$WEIGHT_NORM_OVERALL_V2)),"incident_htn"],combined_inchtn_sdbpc2mrs_pred)[1:3]

inchtn_both_auc$outcome="Incident HTN"
inchtn_both_auc$baseline="Normal"
incdm_both_auc$outcome="Incident DM"
incdm_both_auc$baseline="Normal+PreDM"
normvsdm_both_auc$outcome="Incident DM"
normvsdm_both_auc$baseline="Normal"
prevsdm_both_auc$outcome="Incident DM"
prevsdm_both_auc$baseline="PreDM"

incident_both_auc<-rbind(inchtn_both_auc,incdm_both_auc,normvsdm_both_auc,prevsdm_both_auc)%>%
  dplyr::mutate(category=paste0(baseline,"-",outcome))
write.csv(incident_both_auc,file="output/suppl_tables13.csv",row.names = F)
