#################################
# Single Metabolite Association Analysis

# Generate Supplemental Table S2, S4, S5, S6, S7, Supplemental Figure S5

source("code/00_functions.R")
library_list<-c("tidyverse","labelled","stringr","dplyr","plyr","ggplot2","factoextra","survey","purrr","e1071","ggpubr")
lapply(library_list, require, character.only = TRUE)
# System parameters
options(scipen = 999, stringsAsFactors = F)
options(survey.lonely.psu="remove")
base::load("data/merged_data.RData")
base::load("data/metabolomics_data.RData")
pheno_pca<-readRDS("data/phenotype_pca.rds")
hchs_annot<-readxl::read_xlsx("data/2batch_combined_data_V1_only.xlsx",sheet="metabolites.info")
colnames(hchs_annot)<-tolower(colnames(hchs_annot))

# model 1: AGE, GENDER, CENTER, BACKGROUND BMI
# model 2: AGE, GENDER, CENTER, BACKGROUND BMI ALCOHOL_USE CIGARETTE_USE
# model 3: AGE GENDER CENTER BACKGROUND BMI ALCOHOL_USE CIGARETTE_USE 
# DIABETES2_INDICATOR HYPERTENSION LABA70 (fasting glucose) INSULIN_FAST
# HOMA_IR LABA68 (HDL) LABA69 (LDL) LABA66 (Total cholesterol) LABA67(Triglycerides)
# SBPA5(Systolic Blood Pressure) SBPA6(Diastolic Blood Pressure)

mdl1_covariates<-c("AGE", "GENDER", "CENTER", "background","BMI")
mdl2_covariates<-c(mdl1_covariates, "ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010")
mdl3_covariates<-c(mdl2_covariates,"DIABETES2_INDICATOR", "HYPERTENSION", "LABA70", "INSULIN_FAST", "HOMA_IR", "LABA68","LABA69", "LABA66", "LABA67", "SBPA5", "SBPA6")
admin_var<-c("PSU_ID", "STRAT", "WEIGHT","ID")
outcome_var<-c("sdb_pc1","sdb_pc2")

pheno_sma<-pheno_pca[complete.cases(pheno_pca[,mdl1_covariates]),c(mdl3_covariates,admin_var,outcome_var)]

# SMA in batch 1
imp_only_data_sma_b1<-imp_only_data%>%
  dplyr::select(ID, matches("^X"))%>%
  merge(.,pheno_sma,by="ID",all.y = T)
imp_bin_data_sma_b1<-imp_bin_data%>%
  dplyr::select(ID, matches("^X"))%>%
  merge(.,pheno_sma,by="ID",all.y = T)

# sdb pc1
# Continuous (imputed)
# Model 1
sdbpc1_imp_cont_mdl1<-svyreg_loop(data=imp_only_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_cont_mdl1_strat<-strat_svyreg_loop(data=imp_only_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc1_imp_cont_mdl2<-svyreg_loop(data=imp_only_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_cont_mdl2_strat<-strat_svyreg_loop(data=imp_only_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc1_imp_cont_mdl3<-svyreg_loop(data=imp_only_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_cont_mdl3_strat<-strat_svyreg_loop(data=imp_only_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# binary metabolites
# Model 1
sdbpc1_imp_bin_mdl1<-svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_bin_mdl1_strat<-strat_svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc1_imp_bin_mdl2<-svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_bin_mdl2_strat<-strat_svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc1_imp_bin_mdl3<-svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_bin_mdl3_strat<-strat_svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")

# sdb pc2
# Continuous (imputed)
# Model 1
sdbpc2_imp_cont_mdl1<-svyreg_loop(data=imp_only_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_cont_mdl1_strat<-strat_svyreg_loop(data=imp_only_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc2_imp_cont_mdl2<-svyreg_loop(data=imp_only_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_cont_mdl2_strat<-strat_svyreg_loop(data=imp_only_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc2_imp_cont_mdl3<-svyreg_loop(data=imp_only_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_cont_mdl3_strat<-strat_svyreg_loop(data=imp_only_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# binary metabolites
# Model 1
sdbpc2_imp_bin_mdl1<-svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_bin_mdl1_strat<-strat_svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc2_imp_bin_mdl2<-svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_bin_mdl2_strat<-strat_svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc2_imp_bin_mdl3<-svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_bin_mdl3_strat<-strat_svyreg_loop(data=imp_bin_data_sma_b1,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")

# Merge the results
sdbpc1_imp_mdl1<-rbind(sdbpc1_imp_cont_mdl1,sdbpc1_imp_bin_mdl1)%>%
  dplyr::mutate(model="model_1",
                strata="both",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl2<-rbind(sdbpc1_imp_cont_mdl2,sdbpc1_imp_bin_mdl2)%>%
  dplyr::mutate(model="model_2",
                strata="both",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl3<-rbind(sdbpc1_imp_cont_mdl3,sdbpc1_imp_bin_mdl3)%>%
  dplyr::mutate(model="model_3",
                strata="both",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl1_strat<-rbind(sdbpc1_imp_cont_mdl1_strat,sdbpc1_imp_bin_mdl1_strat)%>%
  dplyr::mutate(model="model_1",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl2_strat<-rbind(sdbpc1_imp_cont_mdl2_strat,sdbpc1_imp_bin_mdl2_strat)%>%
  dplyr::mutate(model="model_2",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl3_strat<-rbind(sdbpc1_imp_cont_mdl3_strat,sdbpc1_imp_bin_mdl3_strat)%>%
  dplyr::mutate(model="model_3",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl1<-rbind(sdbpc2_imp_cont_mdl1,sdbpc2_imp_bin_mdl1)%>%
  dplyr::mutate(model="model_1",
                strata="both",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl2<-rbind(sdbpc2_imp_cont_mdl2,sdbpc2_imp_bin_mdl2)%>%
  dplyr::mutate(model="model_2",
                strata="both",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl3<-rbind(sdbpc2_imp_cont_mdl3,sdbpc2_imp_bin_mdl3)%>%
  dplyr::mutate(model="model_3",
                strata="both",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl1_strat<-rbind(sdbpc2_imp_cont_mdl1_strat,sdbpc2_imp_bin_mdl1_strat)%>%
  dplyr::mutate(model="model_1",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl2_strat<-rbind(sdbpc2_imp_cont_mdl2_strat,sdbpc2_imp_bin_mdl2_strat)%>%
  dplyr::mutate(model="model_2",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl3_strat<-rbind(sdbpc2_imp_cont_mdl3_strat,sdbpc2_imp_bin_mdl3_strat)%>%
  dplyr::mutate(model="model_3",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
main_mdl_b1<-rbind(sdbpc1_imp_mdl1,sdbpc1_imp_mdl2,sdbpc1_imp_mdl3,sdbpc1_imp_mdl1_strat,sdbpc1_imp_mdl2_strat,sdbpc1_imp_mdl3_strat,sdbpc2_imp_mdl1,sdbpc2_imp_mdl2,sdbpc2_imp_mdl3,sdbpc2_imp_mdl1_strat,sdbpc2_imp_mdl2_strat,sdbpc2_imp_mdl3_strat)%>%
  dplyr::rename(chemid_x=metabolite)%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name")],by="chemid_x")

write.csv(main_mdl_b1,file="output/suppl_table_s4.csv",row.names = F)

# Interaction analysis
# metabolites that are significant in batch 1 model 1 stratified SMA
sig_strat_mdl1_sdbpc1_cont<-main_mdl_b1%>%
  dplyr::filter(!strata=="both"&model=="model_1"&p_val_fdr<0.05&trait=="sdb_pc1"&is_continuous=="Continuous")
sig_strat_mdl1_sdbpc1_binary<-main_mdl_b1%>%
  dplyr::filter(!strata=="both"&model=="model_1"&p_val_fdr<0.05&trait=="sdb_pc1"&is_continuous=="Dichotomized")
sig_strat_mdl1_sdbpc2_cont<-main_mdl_b1%>%
  dplyr::filter(!strata=="both"&model=="model_1"&p_val_fdr<0.05&trait=="sdb_pc2"&is_continuous=="Continuous")
sig_strat_mdl1_sdbpc2_binary<-main_mdl_b1%>%
  dplyr::filter(!strata=="both"&model=="model_1"&p_val_fdr<0.05&trait=="sdb_pc2"&is_continuous=="Dichotomized")

sdbpc1_imp_cont_mdl1_interaction_rep<-svyreg_loop_interaction(data=imp_only_data_sma_b1[,colnames(imp_only_data_sma_b1)%in%c("ID",sig_strat_mdl1_sdbpc1_cont$chemid_x,"AGE","GENDER","CENTER","background","BMI","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010","DIABETES2_INDICATOR","HYPERTENSION","LABA70","INSULIN_FAST","HOMA_IR","LABA68","LABA69" ,"LABA66","LABA67","SBPA5" , "SBPA6","PSU_ID","STRAT","WEIGHT","sdb_pc1","sdb_pc2")],covar=mdl1_covariates,end=(nrow(sig_strat_mdl1_sdbpc1_cont)+1),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,interaction = "GENDER")%>%
  dplyr::mutate(trait="sdb_pc1")

sdbpc1_imp_binary_mdl1_interaction_rep<-svyreg_loop_interaction(data=imp_bin_data_sma_b1[,colnames(imp_bin_data_sma_b1)%in%c("ID",sig_strat_mdl1_sdbpc1_binary$chemid_x,"AGE","GENDER","CENTER","background","BMI","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010","DIABETES2_INDICATOR","HYPERTENSION","LABA70","INSULIN_FAST","HOMA_IR","LABA68","LABA69" ,"LABA66","LABA67","SBPA5" , "SBPA6","PSU_ID","STRAT","WEIGHT","sdb_pc1","sdb_pc2")],covar=mdl1_covariates,end=(nrow(sig_strat_mdl1_sdbpc1_binary)+1),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,interaction = "GENDER")%>%
  dplyr::mutate(trait="sdb_pc1")

sdbpc2_imp_cont_mdl1_interaction_rep<-svyreg_loop_interaction(data=imp_only_data_sma_b1[,colnames(imp_only_data_sma_b1)%in%c("ID",sig_strat_mdl1_sdbpc2_cont$chemid_x,"AGE","GENDER","CENTER","background","BMI","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010","DIABETES2_INDICATOR","HYPERTENSION","LABA70","INSULIN_FAST","HOMA_IR","LABA68","LABA69" ,"LABA66","LABA67","SBPA5" , "SBPA6","PSU_ID","STRAT","WEIGHT","sdb_pc1","sdb_pc2")],covar=mdl1_covariates,end=(nrow(sig_strat_mdl1_sdbpc2_cont)+1),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,interaction = "GENDER")%>%
  dplyr::mutate(trait="sdb_pc2")

interaction_results_b1<-rbind(sdbpc1_imp_cont_mdl1_interaction_rep,sdbpc1_imp_binary_mdl1_interaction_rep,sdbpc2_imp_cont_mdl1_interaction_rep)%>%
  dplyr::mutate(cross_p_val_fdr=p.adjust(int_cross_p_val ,method="BH"),
                trait_p_val_fdr=p.adjust(trait_p_val ,method="BH"),
                cross_fdr_sig=ifelse(cross_p_val_fdr<0.05,"FDR<0.05","Not Sig"),
                dataset="batch_1")%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name","super_pathway","sub_pathway","hmdb")],by.x="metabolite",by.y="chemid_x")%>%
  dplyr::select(metabolite,n,trait_beta,trait_se,trait_p_val,trait_p_val_fdr,int_main_beta,int_main_se,int_main_p_val,int_cross_beta,int_cross_se,int_cross_p_val ,cross_term,sig_fdr,is_continuous , trait, cross_p_val_fdr, cross_fdr_sig,chemical_name,super_pathway,sub_pathway,hmdb)%>%
  dplyr::rename(chemid_x="metabolite",
                stratifier_beta="int_main_beta",
                stratifier_se="int_main_se",
                stratifier_p_val="int_main_p_val",
                cross_beta="int_cross_beta",
                cross_se="int_cross_se",
                cross_p_val="int_cross_p_val"
  )

# SMA in batch 2
imp_only_data_sma_b2<-imp_only_data_b2%>%
  dplyr::select(ID, matches("^X"))%>%
  merge(.,pheno_sma,by="ID",all.y = T)
imp_bin_data_sma_b2<-imp_bin_data_b2%>%
  dplyr::select(ID, matches("^X"))%>%
  merge(.,pheno_sma,by="ID",all.y = T)

# sdb pc1
# Continuous (imputed)
# Model 1
sdbpc1_imp_cont_mdl1_b2<-svyreg_loop(data=imp_only_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_cont_mdl1_strat_b2<-strat_svyreg_loop(data=imp_only_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc1_imp_cont_mdl2_b2<-svyreg_loop(data=imp_only_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_cont_mdl2_strat_b2<-strat_svyreg_loop(data=imp_only_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc1_imp_cont_mdl3_b2<-svyreg_loop(data=imp_only_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_cont_mdl3_strat_b2<-strat_svyreg_loop(data=imp_only_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# binary metabolites
# Model 1
sdbpc1_imp_bin_mdl1_b2<-svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_bin_mdl1_strat_b2<-strat_svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc1_imp_bin_mdl2_b2<-svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_bin_mdl2_strat_b2<-strat_svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc1_imp_bin_mdl3_b2<-svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F)
sdbpc1_imp_bin_mdl3_strat_b2<-strat_svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")

# sdb pc2
# Continuous 
# Model 1
sdbpc2_imp_cont_mdl1_b2<-svyreg_loop(data=imp_only_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_cont_mdl1_strat_b2<-strat_svyreg_loop(data=imp_only_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc2_imp_cont_mdl2_b2<-svyreg_loop(data=imp_only_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_cont_mdl2_strat_b2<-strat_svyreg_loop(data=imp_only_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc2_imp_cont_mdl3_b2<-svyreg_loop(data=imp_only_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_cont_mdl3_strat_b2<-strat_svyreg_loop(data=imp_only_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_only_data_sma_b2))),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# binary metabolites
# Model 1
sdbpc2_imp_bin_mdl1_b2<-svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_bin_mdl1_strat_b2<-strat_svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 2
sdbpc2_imp_bin_mdl2_b2<-svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_bin_mdl2_strat_b2<-strat_svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl2_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Model 3
sdbpc2_imp_bin_mdl3_b2<-svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F)
sdbpc2_imp_bin_mdl3_strat_b2<-strat_svyreg_loop(data=imp_bin_data_sma_b2,covar=mdl3_covariates,end=max(grep("^X", names(imp_bin_data_sma_b2))),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,stratifier = "GENDER")
# Merge the results
sdbpc1_imp_mdl1_b2<-rbind(sdbpc1_imp_cont_mdl1_b2,sdbpc1_imp_bin_mdl1_b2)%>%
  dplyr::mutate(model="model_1",
                strata="both",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl2_b2<-rbind(sdbpc1_imp_cont_mdl2_b2,sdbpc1_imp_bin_mdl2_b2)%>%
  dplyr::mutate(model="model_2",
                strata="both",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl3_b2<-rbind(sdbpc1_imp_cont_mdl3_b2,sdbpc1_imp_bin_mdl3_b2)%>%
  dplyr::mutate(model="model_3",
                strata="both",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl1_strat_b2<-rbind(sdbpc1_imp_cont_mdl1_strat_b2,sdbpc1_imp_bin_mdl1_strat_b2)%>%
  dplyr::mutate(model="model_1",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl2_strat_b2<-rbind(sdbpc1_imp_cont_mdl2_strat_b2,sdbpc1_imp_bin_mdl2_strat_b2)%>%
  dplyr::mutate(model="model_2",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc1_imp_mdl3_strat_b2<-rbind(sdbpc1_imp_cont_mdl3_strat_b2,sdbpc1_imp_bin_mdl3_strat_b2)%>%
  dplyr::mutate(model="model_3",
                trait="sdb_pc1",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl1_b2<-rbind(sdbpc2_imp_cont_mdl1_b2,sdbpc2_imp_bin_mdl1_b2)%>%
  dplyr::mutate(model="model_1",
                strata="both",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl2_b2<-rbind(sdbpc2_imp_cont_mdl2_b2,sdbpc2_imp_bin_mdl2_b2)%>%
  dplyr::mutate(model="model_2",
                strata="both",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl3_b2<-rbind(sdbpc2_imp_cont_mdl3_b2,sdbpc2_imp_bin_mdl3_b2)%>%
  dplyr::mutate(model="model_3",
                strata="both",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl1_strat_b2<-rbind(sdbpc2_imp_cont_mdl1_strat_b2,sdbpc2_imp_bin_mdl1_strat_b2)%>%
  dplyr::mutate(model="model_1",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl2_strat_b2<-rbind(sdbpc2_imp_cont_mdl2_strat_b2,sdbpc2_imp_bin_mdl2_strat_b2)%>%
  dplyr::mutate(model="model_2",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))
sdbpc2_imp_mdl3_strat_b2<-rbind(sdbpc2_imp_cont_mdl3_strat_b2,sdbpc2_imp_bin_mdl3_strat_b2)%>%
  dplyr::mutate(model="model_3",
                trait="sdb_pc2",
                metaboltie=as.character(metabolite))

main_mdl_b2<-rbind(sdbpc1_imp_mdl1_b2,sdbpc1_imp_mdl2_b2,sdbpc1_imp_mdl3_b2,sdbpc1_imp_mdl1_strat_b2,sdbpc1_imp_mdl2_strat_b2,sdbpc1_imp_mdl3_strat_b2,sdbpc2_imp_mdl1_b2,sdbpc2_imp_mdl2_b2,sdbpc2_imp_mdl3_b2,sdbpc2_imp_mdl1_strat_b2,sdbpc2_imp_mdl2_strat_b2,sdbpc2_imp_mdl3_strat_b2)%>%
  dplyr::rename(chemid_x=metabolite)%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name","hmdb")],by="chemid_x")

write.csv(main_mdl_b2, file = "output/suppl_table_s4_b2.csv",row.names = F)

# Replication in batch 2
main_mdl_b1<-main_mdl_b1%>%
  dplyr::mutate(lower95=beta-1.96*se,
                upper95=beta+1.96*se)
main_mdl_b2<-main_mdl_b2%>%
  dplyr::mutate(lower95=beta-1.96*se,
                upper95=beta+1.96*se)
both_metab_sdbpc1_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="both"&main_mdl_b1$p_val_fdr<0.05&main_mdl_b1$model=="model_1"&main_mdl_b1$trait=="sdb_pc1"),"chemid_x"]))
both_metab_sdbpc2_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="both"&main_mdl_b1$p_val_fdr<0.05&main_mdl_b1$model=="model_1"&main_mdl_b1$trait=="sdb_pc2"),"chemid_x"]))
female_metab_sdbpc1_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="Female"&main_mdl_b1$p_val_fdr<0.05&main_mdl_b1$model=="model_1"&main_mdl_b1$trait=="sdb_pc1"),"chemid_x"]))
female_metab_sdbpc2_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="Female"&main_mdl_b1$p_val_fdr<0.05&main_mdl_b1$model=="model_1"&main_mdl_b1$trait=="sdb_pc2"),"chemid_x"]))
male_metab_sdbpc1_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="Male"&main_mdl_b1$p_val_fdr<0.05&main_mdl_b1$model=="model_1"&main_mdl_b1$trait=="sdb_pc1"),"chemid_x"]))
male_metab_sdbpc2_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="Male"&main_mdl_b1$p_val_fdr<0.05&main_mdl_b1$model=="model_1"&main_mdl_b1$trait=="sdb_pc2"),"chemid_x"]))
gender_metab_sdbpc1_b1<-unique(c(female_metab_sdbpc1_b1,male_metab_sdbpc1_b1))
gender_metab_sdbpc2_b1<-unique(c(female_metab_sdbpc2_b1,male_metab_sdbpc2_b1))

batch1_sdbpc1_main<-main_mdl_b1[which(main_mdl_b1$chemid_x%in%both_metab_sdbpc1_b1&main_mdl_b1$strata=="both"&main_mdl_b1$trait=="sdb_pc1"),]
batch1_sdbpc2_main<-main_mdl_b1[which(main_mdl_b1$chemid_x%in%both_metab_sdbpc2_b1&main_mdl_b1$strata=="both"&main_mdl_b1$trait=="sdb_pc2"),]
batch1_sdbpc1_gender_main<-main_mdl_b1[which(main_mdl_b1$chemid_x%in%gender_metab_sdbpc1_b1&main_mdl_b1$trait=="sdb_pc1"),]
batch1_sdbpc2_gender_main<-main_mdl_b1[which(main_mdl_b1$chemid_x%in%gender_metab_sdbpc2_b1&main_mdl_b1$trait=="sdb_pc2"),]

batch2_sdbpc1_main<-main_mdl_b2[which(main_mdl_b2$chemid_x%in%both_metab_sdbpc1_b1&main_mdl_b2$strata=="both"&main_mdl_b2$trait=="sdb_pc1"),]
batch2_sdbpc2_main<-main_mdl_b2[which(main_mdl_b2$chemid_x%in%both_metab_sdbpc2_b1&main_mdl_b2$strata=="both"&main_mdl_b2$trait=="sdb_pc2"),]
batch2_sdbpc1_gender_main<-main_mdl_b2[which(main_mdl_b2$chemid_x%in%gender_metab_sdbpc1_b1&main_mdl_b2$trait=="sdb_pc1"),]
batch2_sdbpc2_gender_main<-main_mdl_b2[which(main_mdl_b2$chemid_x%in%gender_metab_sdbpc2_b1&main_mdl_b2$trait=="sdb_pc2"),]

batch1_sdbpc1_main$metabolite_mdl<-paste(batch1_sdbpc1_main$chemid_x,batch1_sdbpc1_main$model)
batch1_sdbpc2_main$metabolite_mdl<-paste(batch1_sdbpc2_main$chemid_x,batch1_sdbpc2_main$model)
batch2_sdbpc1_main$metabolite_mdl<-paste(batch2_sdbpc1_main$chemid_x,batch2_sdbpc1_main$model)
batch2_sdbpc2_main$metabolite_mdl<-paste(batch2_sdbpc2_main$chemid_x,batch2_sdbpc2_main$model)
batch2_sdbpc1_gender_main$metabolite_mdl<-paste(batch2_sdbpc1_gender_main$chemid_x,batch2_sdbpc1_gender_main$model,batch2_sdbpc1_gender_main$strata)
batch2_sdbpc2_gender_main$metabolite_mdl<-paste(batch2_sdbpc2_gender_main$chemid_x,batch2_sdbpc2_gender_main$model,batch2_sdbpc2_gender_main$strata)
batch1_sdbpc1_gender_main$metabolite_mdl<-paste(batch1_sdbpc1_gender_main$chemid_x,batch1_sdbpc1_gender_main$model,batch2_sdbpc1_gender_main$strata)
batch1_sdbpc2_gender_main$metabolite_mdl<-paste(batch1_sdbpc2_gender_main$chemid_x,batch1_sdbpc2_gender_main$model,batch2_sdbpc2_gender_main$strata)

# Interaction analysis
sdbpc1_imp_cont_mdl1_interaction_rep_b2<-svyreg_loop_interaction(data=imp_only_data_sma_b2[,colnames(imp_only_data_sma_b2)%in%c("ID",sig_strat_mdl1_sdbpc1_cont$chemid_x,"AGE","GENDER","CENTER","background","BMI","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010","DIABETES2_INDICATOR","HYPERTENSION","LABA70","INSULIN_FAST","HOMA_IR","LABA68","LABA69" ,"LABA66","LABA67","SBPA5" , "SBPA6","PSU_ID","STRAT","WEIGHT","sdb_pc1","sdb_pc2")],covar=mdl1_covariates,end=(nrow(sig_strat_mdl1_sdbpc1_cont)+1),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,interaction = "GENDER")%>%
  dplyr::mutate(trait="sdb_pc1")

sdbpc1_imp_binary_mdl1_interaction_rep_b2<-svyreg_loop_interaction(data=imp_only_data_sma_b2[,colnames(imp_only_data_sma_b2)%in%c("ID",sig_strat_mdl1_sdbpc1_binary$chemid_x,"AGE","GENDER","CENTER","background","BMI","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010","DIABETES2_INDICATOR","HYPERTENSION","LABA70","INSULIN_FAST","HOMA_IR","LABA68","LABA69" ,"LABA66","LABA67","SBPA5" , "SBPA6","PSU_ID","STRAT","WEIGHT","sdb_pc1","sdb_pc2")],covar=mdl1_covariates,end=(nrow(sig_strat_mdl1_sdbpc1_binary)+1),metab_is_cont=F,metab_is_complete=F,trait_original="sdb_pc1",trait_for_model="original",trait_as_predictor=F,interaction = "GENDER")%>%
  dplyr::mutate(trait="sdb_pc1")

sdbpc2_imp_cont_mdl1_interaction_rep_b2<-svyreg_loop_interaction(data=imp_only_data_sma_b2[,colnames(imp_only_data_sma_b2)%in%c("ID",sig_strat_mdl1_sdbpc2_cont$chemid_x,"AGE","GENDER","CENTER","background","BMI","ALCOHOL_USE","CIGARETTE_USE","GPAQ_TOTAL_MET","AHEI2010","DIABETES2_INDICATOR","HYPERTENSION","LABA70","INSULIN_FAST","HOMA_IR","LABA68","LABA69" ,"LABA66","LABA67","SBPA5" , "SBPA6","PSU_ID","STRAT","WEIGHT","sdb_pc1","sdb_pc2")],covar=mdl1_covariates,end=(nrow(sig_strat_mdl1_sdbpc2_cont)+1),metab_is_cont=T,metab_is_complete=F,trait_original="sdb_pc2",trait_for_model="original",trait_as_predictor=F,interaction = "GENDER")%>%
  dplyr::mutate(trait="sdb_pc2")

interaction_results_b2<-rbind(sdbpc1_imp_cont_mdl1_interaction_rep_b2,sdbpc1_imp_binary_mdl1_interaction_rep_b2,sdbpc2_imp_cont_mdl1_interaction_rep_b2)%>%
  dplyr::mutate(cross_p_val_fdr=p.adjust(int_cross_p_val ,method="BH"),
                trait_p_val_fdr=p.adjust(trait_p_val ,method="BH"),
                cross_fdr_sig=ifelse(cross_p_val_fdr<=0.05,"FDR<0.05","Not Sig"),
                dataset="batch_2")%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name","super_pathway","sub_pathway","hmdb")],by.x="metabolite",by.y="chemid_x")%>%
  dplyr::select(metabolite,n,trait_beta,trait_se,trait_p_val,trait_p_val_fdr,int_main_beta,int_main_se,int_main_p_val,int_cross_beta,int_cross_se,int_cross_p_val ,cross_term,sig_fdr,is_continuous , trait, cross_p_val_fdr, cross_fdr_sig,chemical_name,super_pathway,sub_pathway,hmdb)%>%
  dplyr::rename(chemid_x="metabolite",
                stratifier_beta="int_main_beta",
                stratifier_se="int_main_se",
                stratifier_p_val="int_main_p_val",
                cross_beta="int_cross_beta",
                cross_se="int_cross_se",
                cross_p_val="int_cross_p_val"
  )

interaction_results<-rbind(interaction_results_b1,interaction_results_b2)
write.csv(interaction_results, file = "output/suppl_table_s6.csv",row.names = F)

################################
# One-sided p:for beta of the same direction, adjusted p=unadjusted p/2, if beta have opposite signs, adjusted p = 1-unadjusted p/2
batch2_sdbpc1_main<-merge(batch2_sdbpc1_main,batch1_sdbpc1_main[,c("metabolite_mdl","trait","strata","beta","p_val_fdr")],by=c("metabolite_mdl","trait","strata"),suffixes = c("",".b1"))%>%
  dplyr::mutate(beta_ind=beta*beta.b1,
                onesided_p=ifelse(beta_ind>0,p_val/2,1-p_val/2))%>%
  dplyr::select(-beta.b1,-beta_ind)

batch2_sdbpc2_main<-merge(batch2_sdbpc2_main,batch1_sdbpc2_main[,c("metabolite_mdl","trait","strata","beta","p_val_fdr")],by=c("metabolite_mdl","trait","strata"),suffixes = c("",".b1"))%>%
  dplyr::mutate(beta_ind=beta*beta.b1,
                onesided_p=ifelse(beta_ind>0,p_val/2,1-p_val/2))%>%
  dplyr::select(-beta.b1,-beta_ind)

batch2_sdbpc1_gender_main<-merge(batch2_sdbpc1_gender_main,batch1_sdbpc1_gender_main[,c("metabolite_mdl","trait","strata","beta","p_val_fdr")],by=c("metabolite_mdl","trait","strata"),suffixes = c("",".b1"))%>%
  dplyr::mutate(beta_ind=beta*beta.b1,
                onesided_p=ifelse(beta_ind>0,p_val/2,1-p_val/2))%>%
  dplyr::select(-beta.b1,-beta_ind)
batch2_sdbpc2_gender_main<-merge(batch2_sdbpc2_gender_main,batch1_sdbpc2_gender_main[,c("metabolite_mdl","trait","strata","beta","p_val_fdr")],by=c("metabolite_mdl","trait","strata"),suffixes = c("",".b1"))%>%
  dplyr::mutate(beta_ind=beta*beta.b1,
                onesided_p=ifelse(beta_ind>0,p_val/2,1-p_val/2))%>%
  dplyr::select(-beta.b1,-beta_ind)

validated_batch2_main<-rbind(batch2_sdbpc1_main,batch2_sdbpc2_main)%>%
  dplyr::mutate(fdr_one_sided_p=p.adjust(onesided_p,method="BH"))

validated_batch2_gender<-rbind(batch2_sdbpc1_gender_main,batch2_sdbpc2_gender_main)%>%
  dplyr::mutate(fdr_one_sided_p=p.adjust(onesided_p,method="BH"))

batch2_sdbpc1_main<-merge(batch2_sdbpc1_main,validated_batch2_main[,c("metabolite_mdl","fdr_one_sided_p")],by="metabolite_mdl",all.x=T)
batch2_sdbpc2_main<-merge(batch2_sdbpc2_main,validated_batch2_main[,c("metabolite_mdl","fdr_one_sided_p")],by="metabolite_mdl",all.x=T)
batch2_sdbpc1_gender_main<-merge(batch2_sdbpc1_gender_main,validated_batch2_gender[,c("metabolite_mdl","fdr_one_sided_p","strata")],by=c("metabolite_mdl","strata"),all.x=T)
batch2_sdbpc2_gender_main<-merge(batch2_sdbpc2_gender_main,validated_batch2_gender[,c("metabolite_mdl","fdr_one_sided_p","strata")],by="metabolite_mdl",all.x=T)


batch1_batch2_sdbpc1_sma<-merge(batch1_sdbpc1_main[,c("chemid_x","beta","se","lower95","upper95","n","p_val","p_val_fdr","chemical_name","model","metabolite_mdl")],batch2_sdbpc1_main[,c("beta","se","lower95","upper95","n","p_val","onesided_p","fdr_one_sided_p","metabolite_mdl")],all.x=TRUE,by="metabolite_mdl",suffixes = c("_b1","_b2"))%>%
  merge(.,hchs_annot[,c("sub_pathway","super_pathway","chemical_name","hmdb")],by="chemical_name",all.x=T)
batch1_batch2_sdbpc2_sma<-merge(batch1_sdbpc2_main[,c("chemid_x","beta","se","lower95","upper95","n","p_val","p_val_fdr","chemical_name","model","metabolite_mdl")],batch2_sdbpc2_main[,c("beta","se","lower95","upper95","n","p_val","onesided_p","fdr_one_sided_p","metabolite_mdl")],all.y=TRUE,by="metabolite_mdl",suffixes = c("_b1","_b2"))%>%
  merge(.,hchs_annot[,c("sub_pathway","super_pathway","chemical_name","hmdb")],by="chemical_name",all.x=T)

batch1_batch2_sdbpc1_gender_sma<-merge(batch1_sdbpc1_gender_main[,c("chemid_x","beta","se","lower95","upper95","n","p_val","p_val_fdr","chemical_name","model","metabolite_mdl","strata")],batch2_sdbpc1_gender_main[,c("beta","se","lower95","upper95","n","p_val","onesided_p","fdr_one_sided_p","metabolite_mdl")],all.y=TRUE,by="metabolite_mdl",suffixes = c("_b1","_b2"))%>%
  merge(.,hchs_annot[,c("sub_pathway","super_pathway","chemical_name","hmdb")],by="chemical_name",all.x=T)%>%
  dplyr::mutate(trait="sdb_pc1")

batch1_batch2_sdbpc2_gender_sma<-merge(batch1_sdbpc2_gender_main[,c("chemid_x","beta","se","lower95","upper95","n","p_val","p_val_fdr","chemical_name","model","metabolite_mdl","strata")],batch2_sdbpc2_gender_main[,c("beta","se","lower95","upper95","n","p_val","onesided_p","fdr_one_sided_p","metabolite_mdl")],all.y=TRUE,by="metabolite_mdl",suffixes = c("_b1","_b2"))%>%
  merge(.,hchs_annot[,c("sub_pathway","super_pathway","chemical_name","hmdb")],by="chemical_name",all.x=T)%>%
  dplyr::mutate(trait="sdb_pc2")

batch1_batch2_gender_sma<-rbind(batch1_batch2_sdbpc1_gender_sma,batch1_batch2_sdbpc2_gender_sma)
write.csv(batch1_batch2_gender_sma, file = "output/suppl_table_s5.csv",row.names = F)

# export the batch 1 and batch 2 validated sma results
save(list=c("batch1_batch2_sdbpc1_sma","batch1_batch2_sdbpc2_sma"), file="data/b1b2_sma.RData",version=2) 


# reshape the output results
wide_batch1_batch2_sdbpc1_sma<-batch1_batch2_sdbpc1_sma%>%
  tidyr::separate(model, into = c("model", "measurement"), sep = "_")%>%
  tidyr::pivot_wider(
    id_cols = c(chemical_name, super_pathway, sub_pathway, hmdb, model),
    names_from = measurement,
    values_from = c(beta_b1, lower95_b1, upper95_b1, p_val_b1, p_val_fdr,onesided_p, fdr_one_sided_p)
  )%>%
  mutate(
    beta_lower_upper_b1_1 = paste(round(beta_b1_1,digits=2), "[", round(lower95_b1_1,digits=2), ", ", round(upper95_b1_1,digits=2), "]", sep = ""),
    beta_lower_upper_b1_2 = paste(round(beta_b1_2,digits=2), "[", round(lower95_b1_2,digits=2), ", ", round(upper95_b1_2,digits=2), "]", sep = ""),
    beta_lower_upper_b1_3 = paste(round(beta_b1_3,digits=2), "[", round(lower95_b1_3,digits=2), ", ", round(upper95_b1_3,digits=2), "]", sep = ""),
    trait="sdbpc1"
  ) %>%
  select(-model,-beta_b1_1, -lower95_b1_1, -upper95_b1_1,-beta_b1_2, -lower95_b1_2, -upper95_b1_2,-beta_b1_3, -lower95_b1_3, -upper95_b1_3)%>%
  select(trait,chemical_name,super_pathway,sub_pathway,hmdb,beta_lower_upper_b1_1,p_val_b1_1,p_val_fdr_1,onesided_p_1,fdr_one_sided_p_1,
         beta_lower_upper_b1_2,p_val_b1_2,p_val_fdr_2,onesided_p_2,fdr_one_sided_p_2,beta_lower_upper_b1_3,p_val_b1_3,p_val_fdr_3,onesided_p_3,fdr_one_sided_p_3)
wide_batch1_batch2_sdbpc2_sma<-batch1_batch2_sdbpc2_sma%>%
  tidyr::separate(model, into = c("model", "measurement"), sep = "_")%>%
  tidyr::pivot_wider(
    id_cols = c(chemical_name, super_pathway, sub_pathway, hmdb, model),
    names_from = measurement,
    values_from = c(beta_b1, lower95_b1, upper95_b1, p_val_b1, p_val_fdr,onesided_p, fdr_one_sided_p)
  )%>%
  mutate(
    beta_lower_upper_b1_1 = paste(round(beta_b1_1,digits=2), "[", round(lower95_b1_1,digits=2), ", ", round(upper95_b1_1,digits=2), "]", sep = ""),
    beta_lower_upper_b1_2 = paste(round(beta_b1_2,digits=2), "[", round(lower95_b1_2,digits=2), ", ", round(upper95_b1_2,digits=2), "]", sep = ""),
    beta_lower_upper_b1_3 = paste(round(beta_b1_3,digits=2), "[", round(lower95_b1_3,digits=2), ", ", round(upper95_b1_3,digits=2), "]", sep = ""),
    trait="sdbpc2"
  ) %>%
  select(-model,-beta_b1_1, -lower95_b1_1, -upper95_b1_1,-beta_b1_2, -lower95_b1_2, -upper95_b1_2,-beta_b1_3, -lower95_b1_3, -upper95_b1_3)%>%
  select(trait,chemical_name,super_pathway,sub_pathway,hmdb,beta_lower_upper_b1_1,p_val_b1_1,p_val_fdr_1,onesided_p_1,fdr_one_sided_p_1,beta_lower_upper_b1_2,
         p_val_b1_2,p_val_fdr_2,onesided_p_2,fdr_one_sided_p_2,beta_lower_upper_b1_3,p_val_b1_3,p_val_fdr_3,onesided_p_3,fdr_one_sided_p_3)

wide_batch1_batch2_sma<-rbind(wide_batch1_batch2_sdbpc1_sma,wide_batch1_batch2_sdbpc2_sma)

write.csv(wide_batch1_batch2_sma, file = "output/suppl_table_s2.csv",row.names = F)


# Secondary analysis: SMA against the original SDB phenotypes
both_metab_b1_mdl1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="both"&main_mdl_b1$p_val_fdr<0.05&main_mdl_b1$model=="model_1"),"chemid_x"]))
outcome_7_var<-c("hypoxicburden_harmonized","SLPA36","SLPA54","SLPA91","SLPA92","SLPA97","event_length_sec")

imp_only_data_sdb_sma_b1<-imp_only_data%>%
  dplyr::select(ID, matches("^X"))%>%
  merge(.,pheno_pca[complete.cases(pheno_pca[,mdl1_covariates]),c(mdl3_covariates,admin_var,outcome_7_var,outcome_var)],by="ID",all.y = T)%>%
  dplyr::rename(HB="hypoxicburden_harmonized",
                avgSpO2="SLPA92",
                minSpO2="SLPA91",
                REI0="SLPA36",
                Per90="SLPA97",
                REI3="SLPA54",
                avgEventLength="event_length_sec")
imp_bin_data_sdb_sma_b1<-imp_bin_data%>%
  dplyr::select(ID, matches("^X"))%>%
  merge(.,pheno_pca[complete.cases(pheno_pca[,mdl1_covariates]),c(mdl3_covariates,admin_var,outcome_7_var,outcome_var)],by="ID",all.y = T)%>%
  dplyr::rename(HB="hypoxicburden_harmonized",
                avgSpO2="SLPA92",
                minSpO2="SLPA91",
                REI0="SLPA36",
                Per90="SLPA97",
                REI3="SLPA54",
                avgEventLength="event_length_sec")
outcome_rename_7_var<-c("HB","avgSpO2","minSpO2","REI0","Per90","REI3","avgEventLength","sdb_pc1","sdb_pc2")

sma_7_cont_results<-list()
sma_7_bin_results<-list()
for (i in seq_along(outcome_rename_7_var)){
  sma_7_cont_results[[i]]<-svyreg_loop(data=imp_only_data_sdb_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_only_data_sdb_sma_b1))),metab_is_cont=T,metab_is_complete=F,trait_original=outcome_rename_7_var[i],trait_for_model="original",trait_as_predictor=F)
  sma_7_bin_results[[i]]<-svyreg_loop(data=imp_bin_data_sdb_sma_b1,covar=mdl1_covariates,end=max(grep("^X", names(imp_bin_data_sdb_sma_b1))),metab_is_cont=F,metab_is_complete=F,trait_original=outcome_rename_7_var[i],trait_for_model="original",trait_as_predictor=F)
}
sma_7_results_append<-dplyr::bind_rows(sma_7_cont_results, .id = "outcome")%>%
  rbind(.,dplyr::bind_rows(sma_7_bin_results, .id = "outcome"))%>%
  dplyr::mutate(beta_p=ifelse(beta>=0,-log10(p_val),-log10(p_val)*(-1)))%>%
  merge(.,hchs_b2_annot,by.x="metabolite",by.y="chemid_x",all.x=T)

sma_7_results_append$outcome<-outcome_rename_7_var[as.numeric(sma_7_results_append$outcome)]
sma_7_results_append<-sma_7_results_append[which(sma_7_results_append$metabolite%in%both_metab_b1_mdl1),]

# export the results
write.csv(sma_7_results_append,file="output/suppl_fig_s6_data.csv",row.names = F)

#######################
# Supplemental Figure S5
progestin_sig_cont<-main_mdl_b1%>%
  merge(.,hchs_b2_annot[,c("chemid_x","super_pathway","sub_pathway")])%>%
  dplyr::filter(model=="model_1"&p_val_fdr<0.05&sub_pathway%in%c("Pregnenolone Steroids","Progestin Steroids")&is_continuous=="Continuous")
progestin_sig_cont_list<-sub("X","",unique(progestin_sig_cont$chemid_x))
progestin_sig_bin<-main_mdl_b1%>%
  merge(.,hchs_b2_annot[,c("chemid_x","super_pathway","sub_pathway")])%>%
  dplyr::filter(model=="model_1"&p_val_fdr<0.05&sub_pathway%in%c("Pregnenolone Steroids","Progestin Steroids")&is_continuous=="Dichotomized")
progestin_sig_bin_list<-sub("X","",unique(progestin_sig_bin$chemid_x))

prog_cont_raw_b1<-metab_imp_cont_raw[,colnames(metab_imp_cont_raw)%in%c("ID",progestin_sig_cont_list)]%>%
  merge(.,id_list[,c("PARENT_SAMPLE_NAME","SOL_ID")],by.x="ID",by.y="PARENT_SAMPLE_NAME",all.x=T)%>%
  dplyr::mutate(dataset="batch_1")
prog_cont_raw_b2<-metab_imp_cont_b2_raw[,colnames(metab_imp_cont_b2_raw)%in%c("ID",progestin_sig_cont_list)]%>%
  merge(.,id_list[,c("PARENT_SAMPLE_NAME","SOL_ID")],by.x="ID",by.y="PARENT_SAMPLE_NAME",all.x=T)%>%
  dplyr::mutate(dataset="batch_2")

prog_cont_b1<-imp_only_data_sma_b1[which(imp_only_data_sma_b1$ID%in%imp_only_data$ID),colnames(imp_only_data_sma_b1)%in%c("ID",unique(progestin_sig_cont$chemid_x),"AGE","GENDER")]%>%
  dplyr::mutate(dataset="batch_1",
                gender="GENDER")
prog_bin_b1<-imp_bin_data_sma_b1[which(imp_bin_data_sma_b1$ID%in%imp_bin_data$ID),colnames(imp_bin_data_sma_b1)%in%c("ID",unique(progestin_sig_bin$chemid_x),"AGE","GENDER")]%>%
  dplyr::mutate(dataset="batch_1",
                gender="GENDER")
prog_cont_b2<-imp_only_data_sma_b2[which(imp_only_data_sma_b2$ID%in%imp_only_data_b2$ID),colnames(imp_only_data_sma_b2)%in%c("ID",unique(progestin_sig_cont$chemid_x),"AGE","GENDER")]%>%
  dplyr::mutate(dataset="batch_2",
                gender="GENDER")
prog_bin_b2<-imp_bin_data_sma_b2[which(imp_bin_data_sma_b2$ID%in%imp_bin_data_b2$ID),colnames(imp_bin_data_sma_b2)%in%c("ID",unique(progestin_sig_bin$chemid_x),"AGE","GENDER")]%>%
  dplyr::mutate(dataset="batch_2",
                gender="GENDER")

prog_cont_raw<-rbind(prog_cont_raw_b1,prog_cont_raw_b2)%>%
  merge(.,pheno_pca[,c("ID","AGE","GENDER")],by.x="SOL_ID",by.y="ID",all.x=T)%>%
  dplyr::mutate(
    ID=SOL_ID,
    age_group_5=dplyr::case_when(
      AGE<45 ~ 1,# age group: 1: <40 2: [40-50), 3: [50-55), 4:[55,60) 4:>=60
      AGE>=45&AGE<50 ~ 2,
      AGE>=50&AGE<55 ~ 3,
      AGE>=55&AGE<60 ~ 4,
      AGE>=60 ~ 5,
      TRUE ~ NA_real_
    ),
    age_group=ifelse(!is.na(age_group_5),ifelse(age_group_5<=3,"<50",">=50"),NA),
    age_group_5=factor(age_group_5,levels=c(1:5),labels=c("<45","45-50","50-55","55-60",">=60")),
    gender=GENDER)%>%
  dplyr::select(-SOL_ID)
colnames(prog_cont_raw)[colnames(prog_cont_raw)%in%progestin_sig_cont_list]<-paste0("X",colnames(prog_cont_raw)[colnames(prog_cont_raw)%in%progestin_sig_cont_list])
prog_cont_raw<-prog_cont_raw%>%
  tidyr::gather(.,key="chemid_x",value="concentration",2:7)%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name")],all.x=T)
# dplyr::mutate(chemical_name_level=factor(chemical_name))
prog_cont<-rbind(prog_cont_b1,prog_cont_b2)%>%
  dplyr::mutate(age_group_5=dplyr::case_when(
    AGE<45 ~ 1,# age group: 1: <40 2: [40-50), 3: [50-55), 4:[55,60) 4:>=60
    AGE>=45&AGE<50 ~ 2,
    AGE>=50&AGE<55 ~ 3,
    AGE>=55&AGE<60 ~ 4,
    AGE>=60 ~ 5,
    TRUE ~ NA_real_
  ),
  age_group=ifelse(!is.na(age_group_5),ifelse(age_group_5<=3,"<50",">=50"),NA),
  age_group_5=factor(age_group_5,levels=c(1:5),labels=c("<45","45-50","50-55","55-60",">=60")),
  gender=GENDER) 
prog_cont<-prog_cont[which(prog_cont$ID%in%unique(prog_cont_raw$ID)),]%>%
  tidyr::gather(.,key="chemid_x",value="concentration",2:7)%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name")],all.x=T)
prog_bin<-rbind(prog_bin_b1,prog_bin_b2)%>%
  dplyr::mutate(
    age_group_5=dplyr::case_when(
      AGE<45 ~ 1,# age group: 1: <40 2: [40-50), 3: [50-55), 4:[55,60) 4:>=60
      AGE>=45&AGE<50 ~ 2,
      AGE>=50&AGE<55 ~ 3,
      AGE>=55&AGE<60 ~ 4,
      AGE>=60 ~ 5,
      TRUE ~ NA_real_
    ),
    age_group=ifelse(!is.na(age_group_5),ifelse(age_group_5<=3,"<50",">=50"),NA),
    age_group_5=factor(age_group_5,levels=c(1:5),labels=c("<45","45-50","50-55","55-60",">=60")))
# Calculate percentage of 1 versos 0 in each batch, age group, by gender
prog_bin<-prog_bin[which(prog_bin$ID%in%unique(prog_cont_raw$ID)),]%>%
  tidyr::gather(.,key="chemid_x",value="concentration",2:3)%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name")],all.x=T)%>%
  dplyr::filter(!is.na(concentration))
# 5 age groups
prog_bin_age_5_numer<-prog_bin%>%dplyr::count(dataset,age_group_5,GENDER,chemid_x,concentration)%>%
  dplyr::mutate(index=paste0(dataset,age_group_5,GENDER,chemid_x))
prog_bin_age_5_denom<-prog_bin%>%
  group_by(dataset,age_group_5,GENDER,chemid_x)%>%
  dplyr::summarise(subtotal_n=n())%>%
  dplyr::mutate(index=paste0(dataset,age_group_5,GENDER,chemid_x))%>%
  ungroup()
prog_bin_age_5<-merge(prog_bin_age_5_numer,prog_bin_age_5_denom[,c("index","subtotal_n")],by="index",all.x=T)%>%
  dplyr::mutate(prop=n/subtotal_n,
                gender=GENDER)%>%
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name")],all.x=T)%>%
  dplyr::mutate(sd=sqrt(prop*(1-prop)/subtotal_n))

progesterone_raw_plot<-ggplot(prog_cont_raw, aes(x = age_group_5,y=concentration, fill = gender)) +
  geom_boxplot()  +theme_bw()+ xlab("Age")+ ylab("Concentration (Raw)")+scale_fill_brewer(palette="Set1")+theme(legend.position = "None")+
  facet_grid(col=vars(dataset),row=vars(chemical_name),scales="free_y")
progesterone_rank_plot<-ggplot(prog_cont, aes(x = age_group_5,y=concentration, fill = gender)) +
  geom_boxplot()  +theme_bw()+ xlab("Age")+ ylab("Concentration (Rank-normalized)")+scale_fill_brewer(palette="Set1")+theme(legend.position = "None")+
  facet_grid(col=vars(dataset),row=vars(chemical_name),scales="free_y")
progesterone_bin_plot<-ggplot(data = prog_bin_age_5[which(prog_bin_age_5$concentration==1),], aes(x = age_group_5, y=prop, group = gender)) + 
  geom_bar(
    aes(color = gender, fill = gender),
    stat = "identity", position = position_dodge(0.8),
    width = 0.7
  )+ geom_errorbar(aes(x=age_group_5,ymin=prop-1.96*sd,ymax=prop+1.96*sd, group = gender),width=0.4,position=position_dodge(.8)) +theme_bw()+ xlab("Age")+ ylab("Percentage of non-missing metabolite value")+scale_fill_brewer(palette="Set1")+theme(legend.position = "None")+
  facet_grid(col=vars(dataset),row=vars(chemical_name))
long_col_plot<-ggpubr::ggarrange(progesterone_raw_plot, progesterone_rank_plot, ncol = 2, common.legend = TRUE, legend="bottom")
short_col_plot<-ggpubr::ggarrange(progesterone_bin_plot,NULL, ncol = 1, heights = c(1,2.5))

# export the concentration plot
jpeg("output/suppl_figs5_a.jpg", width = 400, height = 800)
print(long_col_plot)
dev.off()
jpeg("output/suppl_figs5_b.jpg", width = 400, height = 400)
print(progesterone_bin_plot)
dev.off()