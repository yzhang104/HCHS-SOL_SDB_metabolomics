##########################
# Generate Figure 3, Figure 4, Figure S6 SMA results
# supplemental table s3
source("code/00_functions.R")
library_list<-c("tidyverse","labelled","stringr","plyr","dplyr","ggplot2","survey","purrr","RColorBrewer","ggrepel","gridExtra","ggpubr","ComplexHeatmap","circlize")
lapply(library_list, require, character.only = TRUE)
# System parameters
options(scipen = 999, stringsAsFactors = F)
options(survey.lonely.psu="remove")

main_mdl_b1<-read.csv("output/suppl_table_s4.csv")
base::load("data/metabolomics_data.RData")
base::load("data/b1b2_sma.RData")
base::load("data/merged_data.RData")
main_mdl_b1<-merge(main_mdl_b1,hchs_b2_annot[,c("chemid_x","super_pathway","sub_pathway","hmdb")],by="chemid_x",all.x=T)%>%
  dplyr::mutate(beta_p=ifelse(beta>=0,-log10(p_val),-log10(p_val)*(-1)))


model_list<-c("model_1","model_2","model_3")
outcome_list<-c("sdb_pc1","sdb_pc2")
strata_list<-c("Female","Male")
trait_list<-c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous","sdb_pc1","sdb_pc2")
# select the list of metabolites to be included in heatmaps in all three models
both_metab_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="both"&main_mdl_b1$p_val_fdr<0.05),"chemid_x"]))
female_metab_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="Female"&main_mdl_b1$p_val_fdr<0.05),"chemid_x"]))
male_metab_b1<-as.character(unique(main_mdl_b1[which(main_mdl_b1$strata=="Male"&main_mdl_b1$p_val_fdr<0.05),"chemid_x"]))
# merge the female and male metabolite candidate list
total_metab_b1<-union(union(both_metab_b1,female_metab_b1),male_metab_b1) 

# Generate Figure 3 Heatmap based on batch 1 SMA results

# Reformat the main result table to use with complexheatmap
both_beta_ind_b1<-data.frame()
both_p_ind_b1<-data.frame()
both_pfilter_ind_b1<-data.frame()
female_beta_ind_b1<-data.frame()
female_p_ind_b1<-data.frame()
female_pfilter_ind_b1<-data.frame()
male_beta_ind_b1<-data.frame()
male_p_ind_b1<-data.frame()
male_pfilter_ind_b1<-data.frame()
for (j in seq_along(model_list)) {
  # print(j)
  for (i in seq_along(outcome_list)){
    # print(i)
    # both
    both_temp_data_b1<-main_mdl_b1[which(main_mdl_b1$strata=="both"&main_mdl_b1$trait==outcome_list[i]&main_mdl_b1$model==model_list[j]),c("chemid_x","chemical_name","beta","beta_p","p_val_fdr","super_pathway","sub_pathway","is_continuous")]
    both_beta_data_b1<-both_temp_data_b1[,-c(4:5)]
    both_p_data_b1<-both_temp_data_b1[,-c(3,5)]
    both_pfilter_data_b1<-both_temp_data_b1[,-c(3:4)]
    colnames(both_beta_data_b1)[3]<-paste0(outcome_list[i])
    colnames(both_p_data_b1)[3]<-paste0(outcome_list[i])
    colnames(both_pfilter_data_b1)[3]<-paste0(outcome_list[i])
    # female
    female_temp_data_b1<-main_mdl_b1[which(main_mdl_b1$strata=="female"&main_mdl_b1$trait==outcome_list[i]&main_mdl_b1$model==model_list[j]),c("chemid_x","chemical_name","beta","beta_p","p_val_fdr","super_pathway","sub_pathway","is_continuous")]
    female_beta_data_b1<-female_temp_data_b1[,-c(4:5)]
    female_p_data_b1<-female_temp_data_b1[,-c(3,5)]
    female_pfilter_data_b1<-female_temp_data_b1[,-c(3:4)]
    colnames(female_beta_data_b1)[3]<-paste0(outcome_list[i])
    colnames(female_p_data_b1)[3]<-paste0(outcome_list[i])
    colnames(female_pfilter_data_b1)[3]<-paste0(outcome_list[i])
    # male
    male_temp_data_b1<-main_mdl_b1[which(main_mdl_b1$strata=="male"&main_mdl_b1$trait==outcome_list[i]&main_mdl_b1$model==model_list[j]),c("chemid_x","chemical_name","beta","beta_p","p_val_fdr","super_pathway","sub_pathway","is_continuous")]
    male_beta_data_b1<-male_temp_data_b1[,-c(4:5)]
    male_p_data_b1<-male_temp_data_b1[,-c(3,5)]
    male_pfilter_data_b1<-male_temp_data_b1[,-c(3:4)]
    colnames(male_beta_data_b1)[3]<-paste0(outcome_list[i])
    colnames(male_p_data_b1)[3]<-paste0(outcome_list[i])
    colnames(male_pfilter_data_b1)[3]<-paste0(outcome_list[i])
    if (i==1) {
      both_beta_temp_b1<-both_beta_data_b1
      both_p_temp_b1<-both_p_data_b1
      both_pfilter_temp_b1<-both_pfilter_data_b1
      female_beta_temp_b1<-female_beta_data_b1
      female_p_temp_b1<-female_p_data_b1
      female_pfilter_temp_b1<-female_pfilter_data_b1
      male_beta_temp_b1<-male_beta_data_b1
      male_p_temp_b1<-male_p_data_b1
      male_pfilter_temp_b1<-male_pfilter_data_b1
    } else {
      both_beta_temp_b1<-merge(both_beta_temp_b1,both_beta_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      both_p_temp_b1<-merge(both_p_temp_b1,both_p_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      both_pfilter_temp_b1<-merge(both_pfilter_temp_b1,both_pfilter_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      female_beta_temp_b1<-merge(female_beta_temp_b1,female_beta_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      female_p_temp_b1<-merge(female_p_temp_b1,female_p_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      female_pfilter_temp_b1<-merge(female_pfilter_temp_b1,female_pfilter_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      male_beta_temp_b1<-merge(male_beta_temp_b1,male_beta_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      male_p_temp_b1<-merge(male_p_temp_b1,male_p_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
      male_pfilter_temp_b1<-merge(male_pfilter_temp_b1,male_pfilter_data_b1,by=c("chemid_x","chemical_name","super_pathway","sub_pathway","is_continuous"))
    } 
  }
  if (nrow(both_beta_temp_b1)>0){
    both_beta_temp_b1$model<-model_list[j]
  }
  if (nrow(female_beta_temp_b1)>0){
    female_beta_temp_b1$model<-model_list[j]
  }
  if (nrow(male_beta_temp_b1)>0){
    male_beta_temp_b1$model<-model_list[j]
  }
  if (nrow(both_p_temp_b1)>0){
    both_p_temp_b1$model<-model_list[j]
  }
  if (nrow(female_p_temp_b1)>0){
    female_p_temp_b1$model<-model_list[j]
  }
  if (nrow(male_p_temp_b1)>0){
    male_p_temp_b1$model<-model_list[j]
  }
  if (nrow(both_pfilter_temp_b1)>0){
    both_pfilter_temp_b1$model<-model_list[j]
  }
  if (nrow(female_pfilter_temp_b1)>0){
    female_pfilter_temp_b1$model<-model_list[j]
  }
  if (nrow(male_pfilter_temp_b1)>0){
    male_pfilter_temp_b1$model<-model_list[j]
  }
  both_beta_ind_b1<-rbind(both_beta_ind_b1,both_beta_temp_b1)
  both_p_ind_b1<-rbind(both_p_ind_b1,both_p_temp_b1)
  both_pfilter_ind_b1<-rbind(both_pfilter_ind_b1,both_pfilter_temp_b1)
  female_beta_ind_b1<-rbind(female_beta_ind_b1,female_beta_temp_b1)
  female_p_ind_b1<-rbind(female_p_ind_b1,female_p_temp_b1)
  female_pfilter_ind_b1<-rbind(female_pfilter_ind_b1,female_pfilter_temp_b1)
  male_beta_ind_b1<-rbind(male_beta_ind_b1,male_beta_temp_b1)
  male_p_ind_b1<-rbind(male_p_ind_b1,male_p_temp_b1)
  male_pfilter_ind_b1<-rbind(male_pfilter_ind_b1,male_pfilter_temp_b1)
  # print(j)
}
p_both_md1_b1<-both_p_ind_b1[which(both_p_ind_b1$chemid_x%in%both_metab_b1&both_p_ind_b1$model=="model_1"),trait_list]%>%
  data.table::setorder(.,super_pathway,sub_pathway)%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")

p_female_md1_b1<-female_p_ind_b1[which(female_p_ind_b1$chemid_x%in%female_metab_b1&female_p_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
p_male_md1_b1<-male_p_ind_b1[which(male_p_ind_b1$chemid_x%in%male_metab_b1&male_p_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_both_md1_b1<-both_beta_ind_b1[which(both_beta_ind_b1$chemid_x%in%both_metab_b1&both_beta_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_female_md1_b1<-female_beta_ind_b1[which(female_beta_ind_b1$chemid_x%in%female_metab_b1&female_beta_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_male_md1_b1<-male_beta_ind_b1[which(male_beta_ind_b1$chemid_x%in%male_metab_b1&male_beta_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_both_md1_b1<-both_pfilter_ind_b1[which(both_pfilter_ind_b1$chemid_x%in%both_metab_b1&both_pfilter_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_female_md1_b1<-female_pfilter_ind_b1[which(female_pfilter_ind_b1$chemid_x%in%female_metab_b1&female_pfilter_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_male_md1_b1<-male_pfilter_ind_b1[which(male_pfilter_ind_b1$chemid_x%in%male_metab_b1&male_pfilter_ind_b1$model=="model_1"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)

# Model 2
p_both_md2_b1<-both_p_ind_b1[which(both_p_ind_b1$chemid_x%in%both_metab_b1&both_p_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
p_female_md2_b1<-female_p_ind_b1[which(female_p_ind_b1$chemid_x%in%female_metab_b1&female_p_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
p_male_md2_b1<-male_p_ind_b1[which(male_p_ind_b1$chemid_x%in%male_metab_b1&male_p_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_both_md2_b1<-both_beta_ind_b1[which(both_beta_ind_b1$chemid_x%in%both_metab_b1&both_beta_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_female_md2_b1<-female_beta_ind_b1[which(female_beta_ind_b1$chemid_x%in%female_metab_b1&female_beta_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_male_md2_b1<-male_beta_ind_b1[which(male_beta_ind_b1$chemid_x%in%male_metab_b1&male_beta_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_both_md2_b1<-both_pfilter_ind_b1[which(both_pfilter_ind_b1$chemid_x%in%both_metab_b1&both_pfilter_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_female_md2_b1<-female_pfilter_ind_b1[which(female_pfilter_ind_b1$chemid_x%in%female_metab_b1&female_pfilter_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_male_md2_b1<-male_pfilter_ind_b1[which(male_pfilter_ind_b1$chemid_x%in%male_metab_b1&male_pfilter_ind_b1$model=="model_2"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)

# Model 3
p_both_md3_b1<-both_p_ind_b1[which(both_p_ind_b1$chemid_x%in%both_metab_b1&both_p_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
p_female_md3_b1<-female_p_ind_b1[which(female_p_ind_b1$chemid_x%in%female_metab_b1&female_p_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
p_male_md3_b1<-male_p_ind_b1[which(male_p_ind_b1$chemid_x%in%male_metab_b1&male_p_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_both_md3_b1<-both_beta_ind_b1[which(both_beta_ind_b1$chemid_x%in%both_metab_b1&both_beta_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_female_md3_b1<-female_beta_ind_b1[which(female_beta_ind_b1$chemid_x%in%female_metab_b1&female_beta_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
beta_male_md3_b1<-male_beta_ind_b1[which(male_beta_ind_b1$chemid_x%in%male_metab_b1&male_beta_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_both_md3_b1<-both_pfilter_ind_b1[which(both_pfilter_ind_b1$chemid_x%in%both_metab_b1&both_pfilter_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_female_md3_b1<-female_pfilter_ind_b1[which(female_pfilter_ind_b1$chemid_x%in%female_metab_b1&female_pfilter_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
pfilter_male_md3_b1<-male_pfilter_ind_b1[which(male_pfilter_ind_b1$chemid_x%in%male_metab_b1&male_pfilter_ind_b1$model=="model_3"),trait_list]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous),
                sdb_pc1=ifelse(sdb_pc1<0.05,1,0),
                sdb_pc2=ifelse(sdb_pc2<0.05,1,0)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)

# set up grouping and order for the heatmap
#################################
# super pathway order
sp_color_both_b1 <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(levels(p_both_md1_b1$super_pathway)))
names(sp_color_both_b1) <- levels(p_both_md1_b1$super_pathway)
# sub pathway order
sub_pathway_order<-p_both_md1_b1[,c("sub_pathway","super_pathway")]%>%
  data.table::setorder(.,super_pathway,sub_pathway)

sub_pathway_orderd<-sub_pathway_order%>%
  data.table::setorder(.,super_pathway,sub_pathway)
sp_color_both_sub_b1 <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(nrow(sub_pathway_order))
names(sp_color_both_sub_b1) <- sub_pathway_orderd$sub_pathway
# dataset order
# list of continuous metabolites and binary metabolites
binary_metab_b1<-data.frame(metabolite=grep("^X", names(imp_bin_data), value = TRUE),is_continuous=rep("Dichotomized",length(grep("^X", names(imp_bin_data), value = TRUE))))
impute_metab_b1<-data.frame(metabolite=grep("^X", names(imp_only_data), value = TRUE),is_continuous=rep("Continuous",length(grep("^X", names(imp_only_data), value = TRUE))))
group_metab_b1<-rbind(binary_metab_b1,impute_metab_b1)
sp_color_dataset_b1<-colorRampPalette(RColorBrewer::brewer.pal(2,"Accent"))(length(levels(factor(group_metab_b1$is_continuous))))
names(sp_color_dataset_b1) <- levels(factor(group_metab_b1$is_continuous))

# Reorder sub_pathway according to the order they appear in the table (first sort by super_pathway)
p_both_md1_b1$sub_pathway_ordered<-factor(p_both_md1_b1$sub_pathway,as.character(unique(p_both_md1_b1$sub_pathway)))

# combine the super pathway and sub pathway annotation
row_ha_total_b1 <- HeatmapAnnotation(df = data.frame(super_pathway=p_both_md1_b1$super_pathway,
                                                     sub_pathway=p_both_md1_b1$sub_pathway_ordered,
                                                     continuous_or_dichotomized=p_both_md1_b1$metabolite_model
),
which="row", col = list(
  super_pathway = sp_color_both_b1,
  sub_pathway=sp_color_both_sub_b1,
  continuous_or_dichotomized=sp_color_dataset_b1
),
show_annotation_name = FALSE)

# legend
p_col_fun_b1 = colorRamp2(c(-10,-5, 0,5, 10), c("darkblue","blue", "white", "red", "darkred"))
beta_col_fun_b1=colorRamp2(c(-0.4,-0.2, 0,0.2, 0.4), c("darkblue","blue", "white","red", "darkred"))

# heatmap
ht_p_total_both_md1_b1 <- Heatmap(
  as.matrix(p_both_md1_b1[, c("sdb_pc1", "sdb_pc2")]),
  name = "-log10(p)",
  row_title = NULL,
  column_title = "Model 1",
  row_split = p_both_md1_b1$super_pathway,
  row_gap = unit(1, "mm"),
  left_annotation = row_ha_total_b1,
  border = TRUE,
  cluster_rows = F,
  cluster_columns = F,
  show_column_dend = FALSE,
  column_names_rot = 45,
  row_names_side = "left",
  row_dend_side = "right",
  column_names_gp = gpar(fontsize = c(10)),
  row_names_max_width = unit(12, "cm"), # give more room to the biochemical name
  row_names_gp = gpar(fontsize = 11),
  column_names_centered = T,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (as.matrix(pfilter_both_md1_b1[, c("sdb_pc1", "sdb_pc2")])[i, j] > 0)
      grid.text(
        ifelse(as.matrix(pfilter_both_md1_b1[, c("sdb_pc1", "sdb_pc2")])[i, j] > 0, "*", "") ,
        x,
        y,
        gp = gpar(fontsize = 20),
        just = c("centre", "centre")
      )
  },
  col = p_col_fun_b1
)

ht_p_total_both_md2_b1 <- Heatmap(
  as.matrix(p_both_md2_b1[, c("sdb_pc1", "sdb_pc2")]),
  name = "-log10(p)",
  row_title = NULL,
  column_title = "Model 2",
  row_split = p_both_md2_b1$super_pathway,
  row_gap = unit(1, "mm"),
  left_annotation = row_ha_total_b1,
  border = TRUE,
  cluster_rows = F,
  cluster_columns = F,
  show_column_dend = FALSE,
  column_names_rot = 45,
  row_names_side = "left",
  row_dend_side = "right",
  column_names_gp = gpar(fontsize = c(10)),
  row_names_max_width = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 11),
  column_names_centered = T,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (as.matrix(pfilter_both_md2_b1[, c("sdb_pc1", "sdb_pc2")])[i, j] > 0)
      grid.text(
        ifelse(as.matrix(pfilter_both_md2_b1[, c("sdb_pc1", "sdb_pc2")])[i, j] > 0, "*", "") ,
        x,
        y,
        gp = gpar(fontsize = 20),
        just = c("centre", "centre")
      )
  },
  col = p_col_fun_b1
)

ht_p_total_both_md3_b1 <- Heatmap(
  as.matrix(p_both_md3_b1[, c("sdb_pc1", "sdb_pc2")]),
  name = "-log10(p)",
  row_title = NULL,
  column_title = "Model 3",
  row_split = p_both_md3_b1$super_pathway,
  row_gap = unit(1, "mm"),
  left_annotation = row_ha_total_b1,
  border = TRUE,
  cluster_rows = F,
  cluster_columns = F,
  show_column_dend = FALSE,
  column_names_rot = 45,
  row_names_side = "left",
  row_dend_side = "right",
  column_names_gp = gpar(fontsize = c(10)),
  row_names_max_width = unit(12, "cm"),
  row_names_gp = gpar(fontsize = 11),
  # row_names_centered=T,
  column_names_centered = T,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (as.matrix(pfilter_both_md3_b1[, c("sdb_pc1", "sdb_pc2")])[i, j] > 0)
      grid.text(
        ifelse(as.matrix(pfilter_both_md3_b1[, c("sdb_pc1", "sdb_pc2")])[i, j] > 0, "*", "") ,
        x,
        y,
        gp = gpar(fontsize = 20),
        just = c("centre", "centre")
      )
  },
  col = p_col_fun_b1
)

total_p_ht_b1<-ht_p_total_both_md1_b1+ht_p_total_both_md2_b1+ht_p_total_both_md3_b1
draw(total_p_ht_b1, 
     column_title="Heatmap of unadjusted p values in both genders (Batch 1)", 
     row_title_side="right",
     row_title="Metabolites", 
     ht_gap=unit(0.5,"cm"),
     heatmap_legend_side = c("right"),
     merge_legends=T)

# export the scree plot
jpeg("output/fig3.jpg", width = 1100, height = 600)
print(draw(total_p_ht_b1, 
           column_title="Heatmap of unadjusted p values in both genders (Batch 1)", 
           row_title_side="right",
           row_title="Metabolites", 
           ht_gap=unit(0.5,"cm"),
           heatmap_legend_side = c("right"),
           merge_legends=T))
dev.off()

# export the biochemical information for the metabolites from the heatmap
write.csv(hchs_b2_annot[which(hchs_b2_annot$chemid_x%in%total_metab_b1),],file="output/suppl_table_s3.csv",row.names = F)

#####################
# Generate Figure 4 
validated_metabolite_pc1<-as.character(unique(batch1_batch2_sdbpc1_sma[which(batch1_batch2_sdbpc1_sma$fdr_one_sided_p<0.05&batch1_batch2_sdbpc1_sma$p_val_fdr<0.05),"chemid_x"]))
validated_metabolite_pc2<-as.character(unique(batch1_batch2_sdbpc2_sma[which(batch1_batch2_sdbpc2_sma$fdr_one_sided_p<0.05&batch1_batch2_sdbpc2_sma$p_val_fdr<0.05),"chemid_x"]))

batch1_batch2_sdbpc1_sma<-batch1_batch2_sdbpc1_sma%>%
  dplyr::mutate(p_label_b2=dplyr::case_when(
    fdr_one_sided_p > 0.05 ~ "",
    fdr_one_sided_p > 0.01 ~ "*",
    fdr_one_sided_p > 0.001 ~ "**",
    !is.na(fdr_one_sided_p) ~ "***",
    TRUE ~ NA_character_
  ),
  p_label_b1=dplyr::case_when(
    p_val_fdr > 0.05 ~ "",
    p_val_fdr > 0.01 ~ "*",
    p_val_fdr > 0.001 ~ "**",
    !is.na(p_val_fdr) ~ "***",
    TRUE ~ NA_character_
  ),
  trait="sdbpc1"
  )

batch1_batch2_sdbpc2_sma<-batch1_batch2_sdbpc2_sma%>%
  dplyr::mutate(p_label_b2=dplyr::case_when(
    fdr_one_sided_p > 0.05 ~ "",
    fdr_one_sided_p > 0.01 ~ "*",
    fdr_one_sided_p > 0.001 ~ "**",
    !is.na(fdr_one_sided_p) ~ "***",
    TRUE ~ NA_character_
  ),
  p_label_b1=dplyr::case_when(
    p_val_fdr > 0.05 ~ "",
    p_val_fdr > 0.01 ~ "*",
    p_val_fdr > 0.001 ~ "**",
    !is.na(p_val_fdr) ~ "***",
    TRUE ~ NA_character_
  ),
  trait="sdbpc2"
  )

batch1_batch2_sma<-rbind(batch1_batch2_sdbpc1_sma[,c("chemical_name","chemid_x","beta_b1","beta_b2","lower95_b1","upper95_b1","lower95_b2","upper95_b2","model","trait","p_val_fdr","fdr_one_sided_p","p_label_b1","p_label_b2")],batch1_batch2_sdbpc2_sma[,c("chemical_name","chemid_x","beta_b1","beta_b2","lower95_b1","upper95_b1","lower95_b2","upper95_b2","model","trait","p_val_fdr","fdr_one_sided_p","p_label_b1","p_label_b2")])%>%
  tidyr::gather(key="beta_batch",value="beta",c("beta_b1","beta_b2"))%>%
  tidyr::gather(key="lower95_batch",value="lower95",c("lower95_b1","lower95_b2"))%>%
  tidyr::gather(key="upper95_batch",value="upper95",c("upper95_b1","upper95_b2"))%>%
  tidyr::gather(key="p_label_batch",value="p_label",c("p_label_b1","p_label_b2"))%>%
  dplyr::mutate(beta_batch=sub(".*_", "", beta_batch),
                lower95_batch=sub(".*_", "", lower95_batch),
                upper95_batch=sub(".*_", "", upper95_batch),
                p_label_batch=sub(".*_", "", p_label_batch),
                model_batch=paste0(model,".",beta_batch))%>%
  filter(beta_batch==lower95_batch&lower95_batch==upper95_batch&upper95_batch==p_label_batch)

validated_sdbpc1_plot<-ggplot(batch1_batch2_sma[which(batch1_batch2_sma$chemid_x%in%validated_metabolite_pc1&batch1_batch2_sma$trait=="sdbpc1"),], aes(chemical_name, beta, color = model_batch, group = model_batch)) +
  geom_errorbar(aes(min = lower95, max = upper95), size = 0.75, 
                width = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8), color = "black",size=0.8) +
  geom_text(aes(label = p_label, color = model_batch, group = model_batch),position = position_dodge(width=0.8),hjust=-1.5,vjust=0.2)+
  geom_hline(yintercept=0, linetype="dashed",color = "red", size=0.5)+
  theme_bw()+ scale_color_brewer(palette="Paired")+
  coord_flip()+ ggtitle("SDB PC1")  +xlab("")+
  theme(plot.title = element_text(hjust = 1),axis.text.y = element_text(size = 11))

validated_sdbpc2_plot<-ggplot(batch1_batch2_sma[which(batch1_batch2_sma$chemid_x%in%validated_metabolite_pc2&batch1_batch2_sma$trait=="sdbpc2"),], aes(chemical_name, beta, color = model_batch, group = model_batch)) +
  geom_errorbar(aes(min = lower95, max = upper95), size = 0.75, 
                width = 0.5, position = position_dodge(width = 0.8)) +
  geom_point(position = position_dodge(width = 0.8), color = "black",size=0.8) +
  geom_text(aes(label = p_label, color = model_batch, group = model_batch),position = position_dodge(width=0.8),hjust=-1.5,vjust=0.2)+
  geom_hline(yintercept=0, linetype="dashed",color = "red", size=0.5)+
  theme_bw()+ scale_color_brewer(palette="Paired")+
  coord_flip()+ ggtitle("SDB PC2")  +xlab("")+
  theme(plot.title = element_text(hjust = 1),axis.text.y = element_text(size = 11))

validated_bar_plot <- ggarrange(validated_sdbpc1_plot, validated_sdbpc2_plot,
                                heights = c( 4,3),
                                ncol = 1,  common.legend = TRUE)
annotate_figure(validated_bar_plot,
                top = text_grob("Associations between SDB PCs and metabolites\nvalidated in HCHS/SOL batch 2", color = "black", face = "bold", size = 14))

# export the scree plot
jpeg("output/fig4.jpg", width = 600, height = 600)
print(validated_bar_plot)
dev.off()

#######################
# Supplemental Figure S6
sma_7_results_append<-read.csv(file="output/suppl_fig_s6_data.csv")
outcome_rename_7_var<-c("HB","avgSpO2","minSpO2","REI0","Per90","REI3","avgEventLength","sdb_pc1","sdb_pc2")
sma7_metab_b1<-as.character(unique(sma_7_results_append[which(sma_7_results_append$p_val_fdr<0.05),"metabolite"]))

# convert the result table to complexheatmap format
for (i in seq_along(outcome_rename_7_var)){
  both_temp_data_b1<-sma_7_results_append[which(sma_7_results_append$outcome==outcome_rename_7_var[i]),c("metabolite","chemical_name","beta","beta_p","p_val","p_val_fdr","super_pathway","sub_pathway","is_continuous")]
  both_beta_data_b1<-both_temp_data_b1[,-c(4:6)]
  both_p_data_b1<-both_temp_data_b1[,-c(3,5,6)]
  both_pfilter_data_b1<-both_temp_data_b1[,-c(3:5)]
  both_unjp_data_b1<-both_temp_data_b1[,-c(3,4,6)]
  colnames(both_beta_data_b1)[3]<-paste0(outcome_rename_7_var[i])
  colnames(both_p_data_b1)[3]<-paste0(outcome_rename_7_var[i])
  colnames(both_pfilter_data_b1)[3]<-paste0(outcome_rename_7_var[i])
  colnames(both_unjp_data_b1)[3]<-paste0(outcome_rename_7_var[i])
  if (i==1) {
    both_beta_temp_b1<-both_beta_data_b1
    both_p_temp_b1<-both_p_data_b1
    both_pfilter_temp_b1<-both_pfilter_data_b1
    both_unjp_temp_b1<-both_unjp_data_b1
  } else {
    both_beta_temp_b1<-merge(both_beta_temp_b1,both_beta_data_b1,by=c("metabolite","chemical_name","super_pathway","sub_pathway","is_continuous"))
    both_p_temp_b1<-merge(both_p_temp_b1,both_p_data_b1,by=c("metabolite","chemical_name","super_pathway","sub_pathway","is_continuous"))
    both_pfilter_temp_b1<-merge(both_pfilter_temp_b1,both_pfilter_data_b1,by=c("metabolite","chemical_name","super_pathway","sub_pathway","is_continuous"))
    both_unjp_temp_b1<-merge(both_unjp_temp_b1,both_unjp_data_b1,by=c("metabolite","chemical_name","super_pathway","sub_pathway","is_continuous"))
  }
}

p_both_md1_b1<-both_p_temp_b1[which(both_p_temp_b1$metabolite%in%sma7_metab_b1),c(outcome_rename_7_var,"super_pathway","sub_pathway","is_continuous","chemical_name")]%>%
  data.table::setorder(.,super_pathway,sub_pathway)%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")
beta_both_md1_b1<-both_beta_temp_b1[which(both_beta_temp_b1$metabolite%in%sma7_metab_b1),c(outcome_rename_7_var,"super_pathway","sub_pathway","is_continuous","chemical_name")]%>%
  data.table::setorder(.,super_pathway,sub_pathway)%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")
unjp_both_md1_b1<-both_unjp_temp_b1[which(both_unjp_temp_b1$metabolite%in%sma7_metab_b1),c(outcome_rename_7_var,"super_pathway","sub_pathway","is_continuous","chemical_name")]%>%
  data.table::setorder(.,super_pathway,sub_pathway)%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                metabolite_model=factor(is_continuous))%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")
pfilter_both_md1_b1<-both_pfilter_temp_b1[which(both_pfilter_temp_b1$metabolite%in%sma7_metab_b1),c(outcome_rename_7_var,"super_pathway","sub_pathway","is_continuous","chemical_name")]%>%
  dplyr::mutate(super_pathway=factor(super_pathway),
                sub_pathway=factor(sub_pathway),
                metabolite_model=factor(is_continuous)
  )%>%
  tibble::remove_rownames()%>%
  tibble::column_to_rownames(., var = "chemical_name")%>%
  data.table::setorder(.,super_pathway,sub_pathway)
unjp_both_md1_b1_bin<-ifelse(unjp_both_md1_b1[,1:7] < 0.01, 1, 0) # flag when unadjusted p<0.01
pfilter_both_md1_b1[,8:9]<-ifelse(pfilter_both_md1_b1[,8:9] < 0.05, 1, 0)
pfilter_both_md1_b1<-cbind(unjp_both_md1_b1_bin,pfilter_both_md1_b1[,8:ncol(pfilter_both_md1_b1)])
# Heatmap annotation
# super pathway order
sp_color_both_b1 <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(length(levels(p_both_md1_b1$super_pathway)))
names(sp_color_both_b1) <- levels(p_both_md1_b1$super_pathway)
# sub pathway order
sub_pathway_order<-p_both_md1_b1[,c("sub_pathway","super_pathway")]%>%
  data.table::setorder(.,super_pathway,sub_pathway)
sub_pathway_orderd<-sub_pathway_order%>%
  data.table::setorder(.,super_pathway,sub_pathway)
sp_color_both_sub_b1 <- colorRampPalette(RColorBrewer::brewer.pal(12,"Paired"))(nrow(sub_pathway_order))
names(sp_color_both_sub_b1) <- sub_pathway_orderd$sub_pathway
sp_color_dataset_b1<-colorRampPalette(RColorBrewer::brewer.pal(2,"Accent"))(length(levels(factor(group_metab_b1$is_continuous))))
names(sp_color_dataset_b1) <- levels(factor(group_metab_b1$is_continuous))

# Reorder sub_pathway according to the order they appear in the table (first sort by super_pathway)
p_both_md1_b1$sub_pathway_ordered<-factor(p_both_md1_b1$sub_pathway,as.character(unique(p_both_md1_b1$sub_pathway)))

# combine the super pathway and sub pathway annotation
row_ha_total_b1 <- HeatmapAnnotation(df = data.frame(super_pathway=p_both_md1_b1$super_pathway,
                                                     sub_pathway=p_both_md1_b1$sub_pathway_ordered,
                                                     continuous_or_dichotomized=p_both_md1_b1$metabolite_model
),
which="row", col = list(
  super_pathway = sp_color_both_b1,
  sub_pathway=sp_color_both_sub_b1,
  continuous_or_dichotomized=sp_color_dataset_b1
),
show_annotation_name = FALSE)

# legend
p_col_fun_b1 = colorRamp2(c(-5,-3, 0,3, 6), c("darkblue","blue", "white", "red", "darkred"))
beta_col_fun_b1=colorRamp2(c(-5,-1, 0,1, 4), c("darkblue","blue", "white","red", "darkred"))
# complex heatmap
suppl_fig6<-Heatmap(
  as.matrix(p_both_md1_b1[, 1:9]),
  name = "-log10(p)",
  row_title = NULL,
  column_title = "Model 1",
  row_split = p_both_md1_b1$super_pathway,
  row_gap = unit(1, "mm"),
  left_annotation = row_ha_total_b1,
  border = TRUE,
  cluster_rows = F,
  cluster_columns = F,
  show_column_dend = FALSE,
  column_names_rot = 0,
  row_names_side = "left",
  row_dend_side = "right",
  column_names_gp = gpar(fontsize = c(10)),
  row_names_max_width = unit(12, "cm"), # give more room to the biochemical name
  row_names_gp = gpar(fontsize = 11),
  # row_names_centered=T,
  column_names_centered = T,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (as.matrix(pfilter_both_md1_b1[, 1:9])[i, j] > 0)
      grid.text(
        ifelse(as.matrix(pfilter_both_md1_b1[, 1:9])[i, j] > 0, "*", "") ,
        x,
        y,
        gp = gpar(fontsize = 20),
        just = c("centre", "centre")
      )
  },
  col = p_col_fun_b1
)

# export the scree plot
jpeg("output/suppl_figs6.jpg", width = 1100, height = 600)
print(suppl_fig6)
dev.off()