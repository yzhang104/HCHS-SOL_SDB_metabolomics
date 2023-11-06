###################################
# LASSO Regression and construction of SDB PC MRSs
# Generate supplemental Table S8, S9, Figure S3(need to add)
source("code/00_functions.R")
library_list<-c("tidyverse","stringr","plyr","dplyr","survey","purrr","glmnet","xlsx","gridExtra")
lapply(library_list, require, character.only = TRUE)
# System parameters
options(scipen = 999, stringsAsFactors = F)
options(survey.lonely.psu="remove")
base::load("data/merged_data.RData")

covar_included<-c("AGE","GENDER","CENTER","background","BMI","sdb_pc1","sdb_pc2","ID")
metab_included<-grep("^X", colnames(imp_only_data), value = TRUE)
imp_only_data<-imp_only_data[,colnames(imp_only_data)%in%c(covar_included,metab_included)]

imp_only_data_female<-imp_only_data[which(imp_only_data$GENDER=="Female"),]
imp_only_data_male<-imp_only_data[which(imp_only_data$GENDER=="Male"),]

# Prepare datasets for LASSO
# both genders
sdbpc1_cont_y<-imp_only_data$sdb_pc1
sdbpc2_cont_y<-imp_only_data$sdb_pc2
# convert categorical variables to dummy variables
cont_x_factors <- model.matrix(imp_only_data$sdb_pc1 ~ imp_only_data$GENDER+imp_only_data$background+imp_only_data$CENTER)
colnames(cont_x_factors)[2:ncol(cont_x_factors)]<-stringr::str_split(colnames(cont_x_factors)[2:ncol(cont_x_factors)], "\\$",2,simplify=T)[,2]
cont_x<- as.matrix(data.frame(imp_only_data[,which(!colnames(imp_only_data)%in%c("sdb_pc1","sdb_pc2","CENTER","background","GENDER","ID"))],cont_x_factors[,2:ncol(cont_x_factors),drop=F]))
# female 
imp_only_data_female_predictor<-as.data.frame(cont_x)%>%
  dplyr::filter(GENDERMale==0)%>%
  dplyr::select(-GENDERMale)
sdbpc1_cont_y_female<-imp_only_data_female$sdb_pc1
sdbpc2_cont_y_female<-imp_only_data_female$sdb_pc2
cont_x_female<- as.matrix(data.frame(imp_only_data_female_predictor[,which(!colnames(imp_only_data_female_predictor)%in%c("sdb_pc1","sdb_pc2","CENTER","background","GENDER","ID"))]))

# male
imp_only_data_male_predictor<-as.data.frame(cont_x)%>%
  dplyr::filter(GENDERMale==1)%>%
  dplyr::select(-GENDERMale)
sdbpc1_cont_y_male<-imp_only_data_male$sdb_pc1
sdbpc2_cont_y_male<-imp_only_data_male$sdb_pc2
cont_x_male<- as.matrix(data.frame(imp_only_data_male_predictor[,which(!colnames(imp_only_data_male_predictor)%in%c("sdb_pc1","sdb_pc2","CENTER","background","GENDER","ID"))]))

# LASSO regression
covar_lasso<-c("AGE","BMI","GENDERMale","backgroundCentral_American","backgroundCuban","backgroundMexican","backgroundPuerto_Rican","backgroundSouth.American","backgroundMulti","CENTERChicago","CENTERMiami","CENTERSan.Diego")
covar_lasso_gender<-c("AGE","BMI","backgroundCentral_American","backgroundCuban","backgroundMexican","backgroundPuerto_Rican","backgroundSouth.American","backgroundMulti","CENTERChicago","CENTERMiami","CENTERSan.Diego")

# both gender
# sdb pc1 in only continuous data
sdbpc1_cont_lasso<-lasso_function(x=cont_x,y=imp_only_data$sdb_pc1,binary_outcome=FALSE,var_not_penalized=covar_lasso,standardized=T)
# sdb pc2 in only continuous data
sdbpc2_cont_lasso<-lasso_function(x=cont_x,y=imp_only_data$sdb_pc2,binary_outcome=FALSE,var_not_penalized=covar_lasso,standardized=T)
# female
# sdb pc1 in only continuous data
sdbpc1_cont_lasso_female<-lasso_function(x=cont_x_female,y=imp_only_data_female$sdb_pc1,binary_outcome=FALSE,var_not_penalized=covar_lasso_gender,standardized=T)
# sdb pc2 in only continuous data
sdbpc2_cont_lasso_female<-lasso_function(x=cont_x_female,y=imp_only_data_female$sdb_pc2,binary_outcome=FALSE,var_not_penalized=covar_lasso_gender,standardized=T)
# male
# sdb pc1 in only continuous data
sdbpc1_cont_lasso_male<-lasso_function(x=cont_x_male,y=imp_only_data_male$sdb_pc1,binary_outcome=FALSE,var_not_penalized=covar_lasso_gender,standardized=T)
# sdb pc2 in only continuous data
sdbpc2_cont_lasso_male<-lasso_function(x=cont_x_male,y=imp_only_data_male$sdb_pc2,binary_outcome=FALSE,var_not_penalized=covar_lasso_gender,standardized=T)

# export coefficients from LASSO
##########
# export coefficients for all metabolites (including the non-selected ones)
sdbpc1_both<-as.data.frame(as.matrix(sdbpc1_cont_lasso$coef_lasso_all))%>%
  dplyr::mutate(metabolite=row.names(.))%>%
  dplyr::rename(sdbpc1_both_coeff="s1")
sdbpc2_both<-as.data.frame(as.matrix(sdbpc2_cont_lasso$coef_lasso_all))%>%
  dplyr::mutate(metabolite=row.names(.))%>%
  dplyr::rename(sdbpc2_both_coeff="s1")
sdbpc1_female<-as.data.frame(as.matrix(sdbpc1_cont_lasso_female$coef_lasso_all))%>%
  dplyr::mutate(metabolite=row.names(.))%>%
  dplyr::rename(sdbpc1_female_coeff="s1")
sdbpc2_female<-as.data.frame(as.matrix(sdbpc2_cont_lasso_female$coef_lasso_all))%>%
  dplyr::mutate(metabolite=row.names(.))%>%
  dplyr::rename(sdbpc2_female_coeff="s1")
sdbpc1_male<-as.data.frame(as.matrix(sdbpc1_cont_lasso_male$coef_lasso_all))%>%
  dplyr::mutate(metabolite=row.names(.))%>%
  dplyr::rename(sdbpc1_male_coeff="s1")
sdbpc2_male<-as.data.frame(as.matrix(sdbpc2_cont_lasso_male$coef_lasso_all))%>%
  dplyr::mutate(metabolite=row.names(.))%>%
  dplyr::rename(sdbpc2_male_coeff="s1")
sdb_lasso_coeff<-merge(sdbpc1_both,sdbpc2_both,by="metabolite",all.x=T,all.y = T)%>%
  merge(.,sdbpc1_female,by="metabolite",all.x=T,all.y = T)%>%
  merge(.,sdbpc2_female,by="metabolite",all.x=T,all.y = T)%>%
  merge(.,sdbpc1_male,by="metabolite",all.x=T,all.y = T)%>%
  merge(.,sdbpc2_male,by="metabolite",all.x=T,all.y = T)
sdb_lasso_coeff<-sdb_lasso_coeff%>%
  filter(!metabolite%in%c(covar_lasso,"(Intercept)"))%>% # exclude the covariates
  merge(.,hchs_b2_annot[,c("chemid_x","chemical_name")],by.x="metabolite",by.y="chemid_x",all.x=T)
write.csv(sdb_lasso_coeff,file = "output/suppl_table_s8.csv",row.names = F)

sdbpc1_cont_coef<-export_lasso_coeff(lasso_function_output=sdbpc1_cont_lasso,
                                     min_col_number=3,
                                     key_col="chemid_x",
                                     drop_covariates=covar_lasso,
                                     metab_reference_table=hchs_b2_annot[,c("chemid_x","chemical_name")])
sdbpc2_cont_coef<-export_lasso_coeff(lasso_function_output=sdbpc2_cont_lasso,
                                     min_col_number=3,
                                     key_col="chemid_x",
                                     drop_covariates=covar_lasso,
                                     metab_reference_table=hchs_b2_annot[,c("chemid_x","chemical_name")])
sdbpc1_cont_coef_female<-export_lasso_coeff(lasso_function_output=sdbpc1_cont_lasso_female,
                                            min_col_number=3,
                                            key_col="chemid_x",
                                            drop_covariates=covar_lasso,
                                            metab_reference_table=hchs_b2_annot[,c("chemid_x","chemical_name")])
sdbpc2_cont_coef_female<-export_lasso_coeff(lasso_function_output=sdbpc2_cont_lasso_female,
                                            min_col_number=3,
                                            key_col="chemid_x",
                                            drop_covariates=covar_lasso,
                                            metab_reference_table=hchs_b2_annot[,c("chemid_x","chemical_name")])
sdbpc1_cont_coef_male<-export_lasso_coeff(lasso_function_output=sdbpc1_cont_lasso_male,
                                          min_col_number=3,
                                          key_col="chemid_x",
                                          drop_covariates=covar_lasso,
                                          metab_reference_table=hchs_b2_annot[,c("chemid_x","chemical_name")])
sdbpc2_cont_coef_male<-export_lasso_coeff(lasso_function_output=sdbpc2_cont_lasso_male,
                                          min_col_number=3,
                                          key_col="chemid_x",
                                          drop_covariates=covar_lasso,
                                          metab_reference_table=hchs_b2_annot[,c("chemid_x","chemical_name")])

# Construct MRSs for each individuals in batch 1
sdbpc1_cont_score_b1<-data.frame(ID=imp_only_data$ID,score=summary_score(data=imp_only_data,coeff = sdbpc1_cont_coef))
sdbpc1_cont_score_male_b1<-data.frame(ID=imp_only_data_male$ID,score=summary_score(data=imp_only_data_male,coeff = sdbpc1_cont_coef_male))
sdbpc1_cont_score_female_b1<-data.frame(ID=imp_only_data_female$ID,score=summary_score(data=imp_only_data_female,coeff = sdbpc1_cont_coef_female))
sdbpc1_cont_score_bygender_b1<-plyr::rbind.fill(sdbpc1_cont_score_female_b1, sdbpc1_cont_score_male_b1)%>%
  merge(.,sdbpc1_cont_score_b1,by="ID",suffixes=c("bygender","bothgender"))
sdbpc1_cont_score_female_b1<-sdbpc1_cont_score_female_b1%>%
  dplyr::rename(female_sdbpc1_mrs="summary")
sdbpc1_cont_score_male_b1<-sdbpc1_cont_score_male_b1%>%
  dplyr::rename(male_sdbpc1_mrs="summary")

sdbpc2_cont_score_b1<-data.frame(ID=imp_only_data$ID,score=summary_score(data=imp_only_data,coeff = sdbpc2_cont_coef))
sdbpc2_cont_score_female_b1<-data.frame(ID=imp_only_data_female$ID,score=summary_score(data=imp_only_data_female,coeff = sdbpc2_cont_coef_female))
sdbpc2_cont_score_male_b1<-data.frame(ID=imp_only_data_male$ID,score=summary_score(data=imp_only_data_male,coeff = sdbpc2_cont_coef_male))

sdbpc2_cont_score_bygender_b1<-plyr::rbind.fill(sdbpc2_cont_score_female_b1, sdbpc2_cont_score_male_b1)%>%
  merge(.,sdbpc2_cont_score_b1,by="ID",suffixes=c("bygender","bothgender"))
sdbpc2_cont_score_female_b1<-sdbpc2_cont_score_female_b1%>%
  dplyr::rename(female_sdbpc2_mrs="summary")
sdbpc2_cont_score_male_b1<-sdbpc2_cont_score_male_b1%>%
  dplyr::rename(male_sdbpc2_mrs="summary")
# derive OSA MRS in batch 1
# read in the annotation file
annot<-read.csv("data/SOL_metabolomics_info_10202017.csv",header = TRUE, sep = ",",stringsAsFactors = FALSE)
# clean up the annotation file
colnames(annot)<-tolower(colnames(annot))
colnames(annot)[9]<-"metabolite"
annot[which(annot$metabolite=="dimethylarginineSDMA+ADMA_std"),"metabolite"]<-"dimethylarginineSDMA.ADMA_std"
annot[which(annot$metabolite=="methylglucopyranosidealpha+beta_std"),"metabolite"]<-"methylglucopyranosidealpha.beta_std"
annot<-annot[which((!stringr::str_detect(annot$metabolite,"^X"))),]
# trim the "_STD" in the annot file METABOLITE_STD column and matched metabolite list
annot$metabolite<-substr(annot$metabolite,1,nchar(annot$metabolite)-4)
osa_lasso_both_coeff<-xlsx::read.xlsx(file = "data/osamrs_supplemental_tables.xlsx",sheetName="LASSO_coeff_osa_both")%>%
  merge(.,hchs_annot[,c("chem_id","chemical_name")],by.x="biochemical",by.y="chemical_name",all.x=T)%>%
  dplyr::filter(!is.na(chem_id))%>%  
  dplyr::mutate(metabolite=paste0("X",chem_id))
osa_lasso_score_b1<-data.frame(ID=imp_only_data$ID, score=summary_score(data=imp_only_data, coeff = osa_lasso_both_coeff[,c("metabolite","coeff")]))
all_score_b1<-data.frame(ID=imp_only_data$ID,sdbpc1_mrs=sdbpc1_cont_score_b1$summary,sdbpc2_mrs=sdbpc2_cont_score_b1$summary,osa_lasso_mrs=osa_lasso_score_b1$summary)%>%
  merge(.,sdbpc1_cont_score_female_b1,by="ID",all.x=T)%>%
  merge(.,sdbpc1_cont_score_male_b1,by="ID",all.x=T)%>%
  merge(.,sdbpc2_cont_score_female_b1,by="ID",all.x=T)%>%
  merge(.,sdbpc2_cont_score_male_b1,by="ID",all.x=T)%>%
  dplyr::mutate(batch="b1")
# standardize scores within the batch
sdbpc_summary_zscore_b1<-all_score_b1[,c("sdbpc1_mrs","sdbpc2_mrs","osa_lasso_mrs","female_sdbpc1_mrs","male_sdbpc1_mrs","female_sdbpc2_mrs","male_sdbpc2_mrs")]%>%scale()
sdbpc_summary_zscore_b1<-data.frame(var=names(attributes(sdbpc_summary_zscore_b1)$'scaled:center'),mean=as.numeric(attributes(sdbpc_summary_zscore_b1)$'scaled:center'),
                                    sd=as.numeric(attributes(sdbpc_summary_zscore_b1)$'scaled:scale'),dataset="batch_1")
all_score_b1[,c("sdbpc1_mrs","sdbpc2_mrs","osa_lasso_mrs","female_sdbpc1_mrs","male_sdbpc1_mrs","female_sdbpc2_mrs","male_sdbpc2_mrs")]<-all_score_b1[,c("sdbpc1_mrs","sdbpc2_mrs","osa_lasso_mrs","female_sdbpc1_mrs","male_sdbpc1_mrs","female_sdbpc2_mrs","male_sdbpc2_mrs")]%>%scale()

# Construct MRSs for each individuals in batch 2
sdbpc1_cont_score_b2 <- map2_dfr(.x=imp_only_data_b2[sdbpc1_cont_coef$chemid_x], .y=sdbpc1_cont_coef$coeff, ~(.x * .y)) %>% 
  dplyr::mutate(sdbpc1_mrs = rowSums(., na.rm = TRUE)) 
sdbpc2_cont_score_b2 <- map2_dfr(.x=imp_only_data_b2[sdbpc2_cont_coef$chemid_x], .y=sdbpc2_cont_coef$coeff, ~(.x * .y)) %>% 
  dplyr::mutate(sdbpc2_mrs = rowSums(., na.rm = TRUE)) 
sdbpc1_cont_score_female_b2 <- map2_dfr(.x=imp_only_data_b2[sdbpc1_cont_coef_female$chemid_x], .y=sdbpc1_cont_coef_female$coeff, ~(.x * .y)) %>% 
  dplyr::mutate(female_sdbpc1_mrs = rowSums(., na.rm = TRUE)) 
sdbpc2_cont_score_female_b2 <- map2_dfr(.x=imp_only_data_b2[sdbpc2_cont_coef_female$chemid_x], .y=sdbpc2_cont_coef_female$coeff, ~(.x * .y)) %>% 
  dplyr::mutate(female_sdbpc2_mrs = rowSums(., na.rm = TRUE)) 
sdbpc1_cont_score_male_b2 <- map2_dfr(.x=imp_only_data_b2[sdbpc1_cont_coef_male$chemid_x], .y=sdbpc1_cont_coef_male$coeff, ~(.x * .y)) %>% 
  dplyr::mutate(male_sdbpc1_mrs = rowSums(., na.rm = TRUE)) 
sdbpc2_cont_score_male_b2 <- map2_dfr(.x=imp_only_data_b2[sdbpc2_cont_coef_male$chemid_x], .y=sdbpc2_cont_coef_male$coeff, ~(.x * .y)) %>% 
  dplyr::mutate(male_sdbpc2_mrs = rowSums(., na.rm = TRUE)) 
osa_lasso_mrs_b2 <- map2_dfr(.x=imp_only_data_b2[osa_lasso_both_coeff$metabolite], .y=osa_lasso_both_coeff$coeff, ~(.x * .y)) %>% 
  dplyr::mutate(osa_lasso_mrs = rowSums(., na.rm = TRUE)) 
all_score_b2<-bind_cols(sdbpc1_cont_score_b2$sdbpc1_mrs,sdbpc2_cont_score_b2$sdbpc2_mrs,sdbpc1_cont_score_female_b2$female_sdbpc1_mrs,sdbpc2_cont_score_female_b2$female_sdbpc2_mrs,sdbpc1_cont_score_male_b2$male_sdbpc1_mrs,sdbpc2_cont_score_male_b2$male_sdbpc2_mrs,osa_lasso_mrs_b2$osa_lasso_mrs)
colnames(all_score_b2)<-c("sdbpc1_mrs","sdbpc2_mrs","female_sdbpc1_mrs","female_sdbpc2_mrs","male_sdbpc1_mrs","male_sdbpc2_mrs","osa_lasso_mrs")
all_score_b2$ID<-imp_only_data_b2$ID
sdbpc_summary_zscore_b2<-all_score_b2[,c("sdbpc1_mrs","sdbpc2_mrs","osa_lasso_mrs","female_sdbpc1_mrs","male_sdbpc1_mrs","female_sdbpc2_mrs","male_sdbpc2_mrs")]%>%scale()
sdbpc_summary_zscore_b2<-data.frame(var=names(attributes(sdbpc_summary_zscore_b2)$'scaled:center'),mean=as.numeric(attributes(sdbpc_summary_zscore_b2)$'scaled:center'),
                                    sd=as.numeric(attributes(sdbpc_summary_zscore_b2)$'scaled:scale'),dataset="batch_2")
all_score_b2[,c("sdbpc1_mrs","sdbpc2_mrs","osa_lasso_mrs","female_sdbpc1_mrs","male_sdbpc1_mrs","female_sdbpc2_mrs","male_sdbpc2_mrs")]<-all_score_b2[,c("sdbpc1_mrs","sdbpc2_mrs","osa_lasso_mrs","female_sdbpc1_mrs","male_sdbpc1_mrs","female_sdbpc2_mrs","male_sdbpc2_mrs")]%>%scale()
all_score_b2$batch="b2"
# export the z score parameters for the scores
sdbpc_summary_zscore<-rbind(sdbpc_summary_zscore_b1,sdbpc_summary_zscore_b2)
write.csv(sdbpc_summary_zscore,file = "output/suppl_table_s9.csv",row.names = F)

# combine scores from batch 1 and 2
all_score_combined<-rbind(all_score_b1,all_score_b2)
write.csv(all_score_combined,file="data/sdb_mrs_scores.csv",row.names = F)


# Generate Supplemental Figure S3
# female specific SDB PC MRS
female_sdbpc1_cont_annot<-hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%sdbpc1_cont_coef_female$chemid_x),]%>%
  dplyr::mutate(strata="female_sdbpc1")
female_sdbpc2_cont_annot<-hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%sdbpc2_cont_coef_female$chemid_x),]%>%
  dplyr::mutate(strata="female_sdbpc2")

female_anno_metabolite_table<-rbind(female_sdbpc1_cont_annot,female_sdbpc2_cont_annot)

female_pathway_temp<-data.table::setDT(female_anno_metabolite_table)[,list(count = .N), by = .(strata,super_pathway)][,list(super_pathway = super_pathway, count = count,
                                                                                                                            percent_fmt = paste0(formatC(count*100/sum(count), digits = 2), "%"),
                                                                                                                            percent_metabolite_candidates = count/sum(count)
), by = strata]
female_pfilter_pathway_plot<-ggplot(data=female_pathway_temp, aes(x=strata, y= percent_metabolite_candidates, fill=super_pathway)) +   
  geom_bar(stat = "identity", width=0.7) +
  theme(legend.position = "None") +
  geom_text(aes(label = ifelse(percent_metabolite_candidates>0.01,percent_fmt,"")),position = position_stack(vjust = 0.5),size=3)+
  ggtitle("Female Only")+labs(x="Selected metabolites", y="% selected metabolite")

# male specific SDB PC MRS
male_sdbpc1_cont_annot<-hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%sdbpc1_cont_coef_male$chemid_x),]%>%
  dplyr::mutate(strata="male_sdbpc1")
male_sdbpc2_cont_annot<-hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%sdbpc2_cont_coef_male$chemid_x),]%>%
  dplyr::mutate(strata="male_sdbpc2")
male_anno_metabolite_table<-rbind(male_sdbpc1_cont_annot,male_sdbpc2_cont_annot)

male_pathway_temp<-data.table::setDT(male_anno_metabolite_table)[,list(count = .N), by = .(strata,super_pathway)][,list(super_pathway = super_pathway, count = count,
                                                                                                                        percent_fmt = paste0(formatC(count*100/sum(count), digits = 2), "%"),
                                                                                                                        percent_metabolite_candidates = count/sum(count)
), by = strata]
male_pfilter_pathway_plot<-ggplot(data=male_pathway_temp, aes(x=strata, y= percent_metabolite_candidates, fill=super_pathway)) +   
  geom_bar(stat = "identity", width=0.7) +
  theme(legend.position = "None") +
  geom_text(aes(label = ifelse(percent_metabolite_candidates>0.01,percent_fmt,"")),position = position_stack(vjust = 0.5),size=3)+
  ggtitle("Male Only")+labs(x="selected metabolites", y="% selected metabolite")



# SDB MRS metabolite + annotation
sdbpc1_cont_annot<-hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%sdbpc1_cont_coef$chemid_x),]%>%
  dplyr::mutate(strata="sdbpc1")

# AHI_cont metabolite +annotation
sdbpc2_cont_annot<-hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%sdbpc2_cont_coef$chemid_x),]%>%
  dplyr::mutate(strata="sdbpc2")

anno_metabolite_table<-rbind(sdbpc1_cont_annot,sdbpc2_cont_annot)

pathway_temp<-data.table::setDT(anno_metabolite_table)[,list(count = .N), by = .(strata,super_pathway)][,list(super_pathway = super_pathway, count = count,
                                                                                                              percent_fmt = paste0(formatC(count*100/sum(count), digits = 2), "%"),
                                                                                                              percent_metabolite_candidates = count/sum(count)
), by = strata]
pfilter_pathway_plot<-ggplot(data=pathway_temp, aes(x=strata, y= percent_metabolite_candidates, fill=super_pathway)) +   
  geom_bar(stat = "identity", width=0.7) +
  geom_text(aes(label = ifelse(percent_metabolite_candidates>0.01,percent_fmt,"")),position = position_stack(vjust = 0.5),size=3)+
  ggtitle("Sex-combined")+labs(x="selected metabolites", y="% selected metabolite")


ttl_pathway_temp<-hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%c(colnames(imp_only_data),colnames(imp_bin_data))),c("super_pathway","chemid_x")]%>%
  group_by(super_pathway)%>%
  dplyr::summarise(count=n())%>%
  dplyr::mutate(percent_fmt=paste0(formatC(100*count/sum(count),digits=2),'%'),
                percent_metabolite_candidates = count/sum(count),
                strata="Both")
ttl_pathway_plot<-ggplot(ttl_pathway_temp, aes(fill=super_pathway, y=percent_metabolite_candidates,x=strata)) +
  theme(legend.position = "None") +
  geom_bar(stat = "identity", width=0.7) +
  geom_text(aes(label = ifelse(percent_metabolite_candidates>0.01,percent_fmt,"")),position = position_stack(vjust = 0.5),size=3)+
  ggtitle("Total Mapped")+labs(x="All metabolites mapped", y="% mapped metabolites")

grid.arrange(ttl_pathway_plot,female_pfilter_pathway_plot,male_pfilter_pathway_plot,pfilter_pathway_plot,ncol=4,widths=c(1,1,1,2))

# export the scree plot
jpeg("output/suppl_figs3.jpg", width = 1100, height = 600)
print(grid.arrange(ttl_pathway_plot,female_pfilter_pathway_plot,male_pfilter_pathway_plot,pfilter_pathway_plot,ncol=4,widths=c(1,1,1,2)))
dev.off()