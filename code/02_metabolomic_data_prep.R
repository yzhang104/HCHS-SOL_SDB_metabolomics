###############################
# Prepare the metabolomic datasets
source("code/00_functions.R")
# Read in batch 1 metabolomic data 
hchs_b1<-readxl::read_xlsx("data/2batch_combined_data_V1_only.xlsx",sheet="data")
hchs_annot<-readxl::read_xlsx("data/2batch_combined_data_V1_only.xlsx",sheet="metabolites.info")
colnames(hchs_annot)<-tolower(colnames(hchs_annot))
hchs_id<-readxl::read_xlsx("data/2batch_combined_data_V1_only.xlsx",sheet="sample.info")
hchs_b1<-hchs_b1%>%
  dplyr::rename(ID="PARENT_SAMPLE_NAME")

# Read in batch 2 metabolomic data
# Match batch 2 ID with batch 1 ID
hchs_b2v1<-readxl::read_xlsx("data/batch2_data.xlsx",sheet="batch2_batchnormalized_v1") # comparable to the original batch 1
hchs_b2v2<-readxl::read_xlsx("data/batch2_data.xlsx",sheet="batch2_batchnormalized_v2") # comparable to the original batch 1
hchs_b2<-hchs_b2v1 #only use batch 2 visit 1
id_list<-readxl::read_xlsx("data/IDlist_7204.xlsx")
# group PARENT_SAMPLE_NAME based on SOL_ID, if B01 and B02 are both present in each SOL_ID then drop the PARENT_SAMPLE_NAME
hchs_id_batch1<-id_list[which(id_list$BATCH=="B01"),]
hchs_id_batch2<-id_list[which(id_list$BATCH=="B02"),]
shared_id<-intersect(hchs_id_batch1$SOL_ID,hchs_id_batch2$SOL_ID) 
# SOL_ID:20698407 showed up in both batch 1 and batch 2 metabolomic datasets but not documented in the id_list file
# add this SOL_ID to the shared_id list
shared_id<-c(shared_id,"20698407")
hchs_id_batch2_only<-id_list[which((!id_list$SOL_ID%in%shared_id)&id_list$BATCH=="B02"&id_list$VISIT=="V01"),] 
hchs_b2_only<-hchs_b2[which(hchs_b2$PARENT_SAMPLE_NAME%in%hchs_id_batch2_only$PARENT_SAMPLE_NAME&(!hchs_b2$PARENT_SAMPLE_NAME%in%hchs_b1$ID)),]%>%
  dplyr::rename(ID="PARENT_SAMPLE_NAME") # drop samples when the same sample id was in both batch 1 and 2, and the same person was in both 1 and 2
hchs_b2_annot<-readxl::read_xlsx("data/batch2_data.xlsx",sheet="metabolites.info")%>%
  dplyr::mutate(chemid_x=paste0("X",CHEM_ID),
                CHEM_ID=as.character(CHEM_ID))
colnames(hchs_b2_annot)<-tolower(colnames(hchs_b2_annot))

########################
# Metabolomic data pre-processing for batch 1
# 1. Assess the missingness in the whole sample: 
# 1.1 Drop samples with more than 25% missing
hchs_b1 <- hchs_b1[rowSums(is.na(hchs_b1)) < (ncol(hchs_b1)-2)*0.25, ]
# 1.2 Drop metabolites with more than 75% missing
hchs_b1 <- hchs_b1[,colSums(is.na(hchs_b1)) < nrow(hchs_b1)*0.75]
# 2. Imputation and dichotomization
# 2.1 Subset metabolites into metab_cont (below metab_na_threshold_imp:25% missing) and metab_bin [metab_na_threshold_imp:25%,metab_na_threshold_drop:75%]
metab_imp_cont<-hchs_b1[,colSums(is.na(hchs_b1)) < nrow(hchs_b1)*0.25] # continuous dataset for imputation
metab_imp_binary<-hchs_b1[,colSums(is.na(hchs_b1)) >= nrow(hchs_b1)*0.25] # dichotomized dataset 
# 2.2 Turn metab_bin into binary data set (missing is 0 while nonmissing value is 1)
metab_imp_binary<-replace(data.frame(metab_imp_binary, stringsAsFactors = FALSE,check.names = F),!is.na(metab_imp_binary), 1)
metab_imp_binary<-replace(data.frame(metab_imp_binary, stringsAsFactors = FALSE,check.names = F),is.na(metab_imp_binary),0)
metab_imp_binary<-cbind(hchs_b1$ID,metab_imp_binary)%>%
  dplyr::rename(ID="hchs_b1$ID")
# 2.3 Impute the continuous data set with the min values
metab_imp_cont_ID<-metab_imp_cont
metab_imp_cont<-impute_func(data=metab_imp_cont[,!colnames(metab_imp_cont)%in%c("ID","LAB_ID")],method="min")%>%
  cbind(metab_imp_cont_ID$ID,.)%>%
  dplyr::rename(ID="metab_imp_cont_ID$ID")
# 3. drop unknown metabolites
metab_imp_cont<-metab_imp_cont[,!str_detect(colnames(metab_imp_cont),"^9999")]
metab_imp_binary<-metab_imp_binary[,!str_detect(colnames(metab_imp_binary),"^9999")]
# 4. drop xenobiotic metabolites
metab_imp_cont<-metab_imp_cont[,!colnames(metab_imp_cont)%in%hchs_b2_annot[which(hchs_b2_annot$super_pathway=="Xenobiotics"),]$chem_id]
metab_imp_binary<-metab_imp_binary[,!colnames(metab_imp_binary)%in%hchs_b2_annot[which(hchs_b2_annot$super_pathway=="Xenobiotics"),]$chem_id]
metab_list_b1_cont<-colnames(metab_imp_cont)[-1]
metab_list_b1_binary<-colnames(metab_imp_binary)[-1]
metab_imp_cont_raw<-metab_imp_cont

# 5. Rank normalization
# 5.1 Rank-normalize continuous metabolites
metab_imp_cont<-apply(metab_imp_cont[,-1],2,rank_normalisation)
metab_imp_cont<- cbind.data.frame(metab_imp_cont_ID$ID,metab_imp_cont)
colnames(metab_imp_cont)[1]<-"ID"

# Metabolomic data pre-processing for batch 2
# 1. Assess the missingness in the whole sample: 
# 1.1 Drop samples with more than 25% missing
hchs_b2_only <- hchs_b2_only[rowSums(is.na(hchs_b2_only)) < (ncol(hchs_b2_only)-2)*0.25, ]
# 1.2 Drop metabolites with more than 75% missing
hchs_b2_only <- hchs_b2_only[,colSums(is.na(hchs_b2_only)) < nrow(hchs_b2_only)*0.75]
 # 2. Imputation and dichotomization
# 2.1 Subset metabolites into metab_cont (below metab_na_threshold_imp:25% missing) and metab_bin [metab_na_threshold_imp:25%,metab_na_threshold_drop:75%]
metab_imp_cont_b2<-hchs_b2_only[,colSums(is.na(hchs_b2_only)) < nrow(hchs_b2_only)*0.25] # continuous dataset for imputation
metab_imp_binary_b2<-hchs_b2_only[,colSums(is.na(hchs_b2_only)) >= nrow(hchs_b2_only)*0.25] # dichotomized dataset 
# 2.2 Turn metab_bin into binary data set (missing is 0 while nonmissing value is 1)
metab_imp_binary_b2<-replace(data.frame(metab_imp_binary_b2, stringsAsFactors = FALSE,check.names = F),!is.na(metab_imp_binary_b2), 1)
metab_imp_binary_b2<-replace(data.frame(metab_imp_binary_b2, stringsAsFactors = FALSE,check.names = F),is.na(metab_imp_binary_b2),0)
metab_imp_binary_b2<-cbind(hchs_b2_only$ID,metab_imp_binary_b2)%>%
  dplyr::rename(ID="hchs_b2_only$ID")
# 2.3 Impute the continuous data set with the min values
metab_imp_cont_b2_ID<-metab_imp_cont_b2
metab_imp_cont_b2<-impute_func(data=metab_imp_cont_b2[,-1],method="min")%>%
  cbind(metab_imp_cont_b2_ID$ID,.)%>%
  dplyr::rename(ID="metab_imp_cont_b2_ID$ID")
# 3. drop unknown metabolites
metab_imp_cont_b2<-metab_imp_cont_b2[,!str_detect(colnames(metab_imp_cont_b2),"^9999")]
metab_imp_binary_b2<-metab_imp_binary_b2[,!str_detect(colnames(metab_imp_binary_b2),"^9999")]
# 4. drop xenobiotic metabolites
metab_imp_cont_b2<-metab_imp_cont_b2[,!colnames(metab_imp_cont_b2)%in%hchs_b2_annot[which(hchs_b2_annot$super_pathway=="Xenobiotics"),"chem_id"]]
metab_imp_binary_b2<-metab_imp_binary_b2[,!colnames(metab_imp_binary_b2)%in%hchs_b2_annot[which(hchs_b2_annot$super_pathway=="Xenobiotics"),"chem_id"]]
metab_list_b2_cont<-colnames(metab_imp_cont_b2)[-1]
metab_list_b2_binary<-colnames(metab_imp_binary_b2)[-1]
metab_imp_cont_b2_raw<-metab_imp_cont_b2
# 5. Rank normalization
# 5.1 Rank-normalize continuous metabolites
metab_imp_cont_b2<-apply(metab_imp_cont_b2[,-1],2,rank_normalisation)
metab_imp_cont_b2<- cbind.data.frame(metab_imp_cont_b2_ID$ID,metab_imp_cont_b2)
colnames(metab_imp_cont_b2)[1]<-"ID"

# 6. Drop unmapped metabolites between batch1 and batch2 
# Batch 1
# continuous metabolites
metab_imp_cont<-metab_imp_cont[,colnames(metab_imp_cont)%in%c("ID",metab_list_b2_cont)]
# binary metabolites
metab_imp_binary<-metab_imp_binary[,colnames(metab_imp_binary)%in%c("ID",metab_list_b2_binary)]
# Batch 2
# continuous metabolites
metab_imp_cont_b2<-metab_imp_cont_b2[,colnames(metab_imp_cont_b2)%in%c("ID",metab_list_b1_cont)]
# binary metabolites
metab_imp_binary_b2<-metab_imp_binary_b2[,colnames(metab_imp_binary_b2)%in%c("ID",metab_list_b1_binary)]


# export the processed metabolomics datasets
save(list=c("metab_imp_cont","metab_imp_binary","metab_imp_cont_b2","metab_imp_binary_b2","hchs_b2_annot","id_list","metab_imp_cont_raw","metab_imp_cont_b2_raw"), file="data/metabolomics_data.RData",version=2) 
