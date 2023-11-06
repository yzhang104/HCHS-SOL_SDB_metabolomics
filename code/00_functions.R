# function to check installed package, load packages
instant_pkgs<- function(pkgs) { 
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss,repos = "http://cran.us.r-project.org", method = "curl")
  }
  
  if (length(pkgs_miss) == 0) {
    message("\n ...Packages were already installed!\n")
  }
  
  # install packages not already loaded:
  pkgs_miss <- pkgs[which(!pkgs %in% installed.packages()[, 1])]
  if (length(pkgs_miss) > 0) {
    install.packages(pkgs_miss,repos = "http://cran.us.r-project.org", method = "curl")
  }
  
  # load packages not already loaded:
  attached <- search()
  attached_pkgs <- attached[grepl("package", attached)]
  need_to_attach <- pkgs[which(!pkgs %in% gsub("package:", "", attached_pkgs))]
  
  if (length(need_to_attach) > 0) {
    for (i in 1:length(need_to_attach)) require(need_to_attach[i], character.only = TRUE)
  }
  
  if (length(need_to_attach) == 0) {
    message("\n ...Packages were already loaded!\n")
  }
}

# remove outliers
remove_outliers <- function(x, na.rm = TRUE, ...) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
  H <- 3 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  y
}

# Function to transform oxygen saturation below 90%
transformation <- function(x){
  val <- rep(0, length(x))
  val[which(x >=5)] <- 1
  val[which(is.na(x))] <- NA
  return(val)
}

# Function to rank normalize var
rank_normalisation<-function(x){
  qnorm((rank(x,na.last="keep")-0.5)/length(x))
}

# Function to impute missing values
impute_func<-function(data,method){
  # mean: replace with mean
  # median: replace with median
  # min: replace with min
  # half: replace with half of the min value
  if (method=="mean"){
    imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),mean(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
  } else if (method=="median") {
    imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),median(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
  }  else if (method=="min") {
    imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),min(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
  } else if (method=="half") {
    imp_data<-data.frame(sapply(data,function(x) ifelse(is.na(x),0.5*min(x, na.rm = TRUE),x)),check.names = F,stringsAsFactors = F)
  } else stop('Method not supported')
  # else if (method=="knn") {
  #   imp_data<-VIM::kNN(data, variable = colnames(data), k = 5, numFun = median)
  # }
  return(imp_data)
}

# function to merge two table 1 (weighted and nonweighted) and to write into a csv
merge_table1<-function(tbl_nonweight,tbl_weight,out_path){
  # options(stringsAsFactors = FALSE)
  tbl_weight<-as.data.frame(tbl_weight, stringsAsFactors = FALSE)
  tbl_nonweight<-as.data.frame(tbl_nonweight, stringsAsFactors = FALSE)
  new_tbl<-tbl_weight # use the weighted table as base (including the p, test, missing)
  rownames(new_tbl)<-rownames(tbl_nonweight)
  if (ncol(tbl_weight)==2){
    end_col<-1
  } else {
    end_col<-ncol(tbl_weight)-2
  }
  for (i in 1:end_col){ # loop over each strata
    for (j in 1:nrow(tbl_weight)){ # loop over each row
      if (stringr::str_detect(rownames(tbl_weight)[j],"mean..SD")){ # if "mean  (SD)" is detected in the rowname
        new_tbl[j,i]<-new_tbl[j,i] # then don't modify the weighted table
      } else if (stringr::str_detect(new_tbl[j,i],"^\\s+$")) { # if the cell is empty (only contains white spaces)
        new_tbl[j,i]<-new_tbl[j,i]
      } else { 
        if (str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,2]==""){ # if no percentage is calculated (for the "n" row)
          new_tbl[j,i]<-tbl_nonweight[j,i] # then paste only the nonweighted count number
        } else { # otherwise stitch the count from the nonweighted and the percentage from the weighted together
          new_tbl[j,i]<-paste0(str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,1],"(",str_split_fixed(tbl_nonweight[j,i], '[(]', 2)[,2])
        }
        
      }
    }
  }
  write.csv(new_tbl, file = out_path)
}

# Function to loop regression model and export model summary using survey weight
svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID","SOL_ID","id")]
  out_nvar=ncol(dat)
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  out_nobs=rep(NA,out_nvar)
  if (trait_for_model=="binary") {
    trait<-"TraitBinary"
  } else if (trait_for_model=="original"){
    trait<-trait_original
  } else { trait<-"Trait" }
  
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                       ifelse(trait_for_model=="binary",FALSE,TRUE))
  
  
  for(i in 1:(ncol(dat))) {
    # met_df<-cbind(data[,i+2],data[,c("ID",covar)])
    if (!"ID"%in%colnames(data)){
      data<-data%>%
        dplyr::rename("ID"="SOL_ID")
    }
    met_df<-cbind(data[,colnames(dat)[i]],data[,c("ID",covar)])
    met_df_id<-met_df[complete.cases(met_df),"ID"]
    require(survey)
    
    if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
      survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
    } else {
      survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
    }
    
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    }
    
    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    out_se[i] = as.numeric(se[2])
    out_pvalue[i] = as.numeric(pval[2])
    out_nobs[i]=as.numeric(nobs(model))
  }
  regress_output<-data.frame(metabolite=colnames(dat),
                             beta=out_beta,
                             se=out_se,
                             n=out_nobs,
                             p_val=out_pvalue,
                             p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             p_val_fdr=p.adjust(out_pvalue,method="BH")
  )%>%
    dplyr::mutate(
      pval_Bonf_neglog=-log10(p_val_Bonf),
      pval_fdr_neglog=-log10(p_val_fdr),
      sig_Bonf=ifelse(p_val_Bonf<0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(p_val_fdr<0.05,"FDR<0.05","Not Sig"),
      is_continuous=metab_is_cont
    )
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
                  sig_fdr_named=ifelse(p_val_fdr_named<0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
                  sig_fdr_named=NA)
  regress_output<-rbind(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Continuous","Dichotomized"))
  }
  
  return(regress_output)
}

# Function to loop regression model and export model summary using survey weight
strat_svyreg_loop<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,stratifier){
  if (trait_for_model=="binary") {
    trait<-"TraitBinary"
  }  else if (trait_for_model=="original"){
    trait<-trait_original
  } else { trait<-"Trait" }
  
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("SOL_ID","ID","LAB_ID")]
  
  covar<-covar[!covar%in%stratifier]
  group_by_value<-unique(data[,stratifier])
  out_nvar<-ncol(dat)
  out_beta<-rep(NA,out_nvar)
  out_se<-rep(NA,out_nvar)
  out_pvalue<-rep(NA,out_nvar)
  out_nobs<-rep(NA,out_nvar)
  group<-rep(NA,out_nvar)
  regress_output<-data.frame()
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont<-ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                        ifelse(trait_for_model=="binary",FALSE,TRUE))
  
  
  for (j in seq_along(group_by_value)){
    newdata<-data[which(data[,stratifier]==as.character(group_by_value[j])),]
    # survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
    if (!"ID"%in%colnames(data)){
      data<-data%>%
        dplyr::rename("ID"="SOL_ID")
      newdata<-newdata%>%
        dplyr::rename("ID"="SOL_ID")
    }
    dat<-newdata[,1:end]
    dat<-dat[,!colnames(dat)%in%c("ID","LAB_ID")]
    
    if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
      survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
    } else {
      survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
    }
    
    for(i in 1:(ncol(dat))) {
      # print(i)
      met_df<-cbind(newdata[,colnames(newdata)%in%colnames(dat)[i]],newdata[,c("ID",covar)])
      met_df_id<-met_df[complete.cases(met_df),"ID"]
      if (trait_as_predictor==T){
        if (metab_is_cont==T){
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
        } else {
          model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      }else {
        if (trait_is_cont==T){
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse= "+"))),design=subset(survey_design,ID%in%met_df_id))
        }else{
          model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar,collapse = "+"))),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
        }
      }
      
      Vcov <- vcov(model, useScale = FALSE)
      beta<- coef(model)
      se<- sqrt(diag(vcov(model, useScale = FALSE)))
      zval<- beta / se
      pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
      out_beta[i]=as.numeric(beta[2])
      out_se[i] = as.numeric(se[2])
      out_pvalue[i] = as.numeric(pval[2])
      out_nobs[i]=as.numeric(nobs(model))
      group[i]=as.character(group_by_value[j])
      # print(i)
    }
    strata_output<-data.frame(metabolite=colnames(dat),
                              beta=out_beta,
                              se=out_se,
                              n=out_nobs,
                              p_val=out_pvalue,
                              p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                              p_val_fdr=p.adjust(out_pvalue,method="BH"),
                              strata=group)%>%
      dplyr::mutate(
        pval_Bonf_neglog=-log10(p_val_Bonf),
        pval_fdr_neglog=-log10(p_val_fdr),
        sig_Bonf=ifelse(p_val_Bonf<0.05,"Bonf<0.05","Not Sig"),
        sig_fdr=ifelse(p_val_fdr<0.05,"FDR<0.05","Not Sig"),
        is_continuous=metab_is_cont
      )
    regress_output<-rbind(regress_output,strata_output)
  }
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::group_by(strata)%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(p_val,method="BH"),
                  sig_fdr_named=ifelse(p_val_fdr_named<0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
                  sig_fdr_named=NA)
  regress_output<-plyr::rbind.fill(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Continuous","Dichotomized"))
  }
  
  return(regress_output)
}

# Function to loop the interaction regression model and export model summary using survey weight
svyreg_loop_interaction<-function(data,covar,end,metab_is_cont,metab_is_complete,trait_original,trait_binary,trait_for_model,trait_as_predictor,interaction){
  dat<-data[,1:end]
  dat<-dat[,!colnames(dat)%in%c("LAB_ID","ID")]
  if (is.vector(dat)){
    dat<-data.frame(dat)
    colnames(dat)[1]<-colnames(data)[2]
  }
  out_nvar=ncol(dat) 
  out_beta=rep(NA,out_nvar)
  out_se=rep(NA,out_nvar)
  out_pvalue=rep(NA,out_nvar)
  out_nobs=rep(NA,out_nvar)
  int_beta=rep(NA,out_nvar)
  int_se=rep(NA,out_nvar)
  int_pvalue=rep(NA,out_nvar)
  cross_beta=rep(NA,out_nvar)
  cross_se=rep(NA,out_nvar)
  cross_pvalue=rep(NA,out_nvar)
  
  # drop the interaction term from the covariate list
  covar2<-covar[!covar%in%interaction]
  
  if (trait_for_model=="binary") {
    trait<-"TraitBinary"
  }  else if (trait_for_model=="original"){
    trait<-trait_original
  } else { trait<-"Trait" }
  
  # if the original trait variable is binary or the user pick the binary version to use in the model
  trait_is_cont=ifelse(length(unique(data[!is.na(data[,trait]),trait]))==2,FALSE,
                       ifelse(trait_for_model=="binary",FALSE,TRUE))
  
  if (all(grep("[A-Z]", names(data[,-1])) == 0)) { # if column names are all in lower case
    survey_design=svydesign(id=~psu_id, strata=~strat, weights=~weight, data=data)
  } else {
    survey_design=svydesign(id=~PSU_ID, strata=~STRAT, weights=~WEIGHT, data=data)
  }
  
  for(i in 1:(ncol(dat))) {
    met_df<-cbind(data[,colnames(dat)[i]],data[,c("ID",covar)])
    met_df_id<-met_df[complete.cases(met_df),"ID"]
    if (trait_as_predictor==T){
      if (metab_is_cont==T){
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2,collapse= "+"),"+",trait,"*",interaction)),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(colnames(dat)[i],"~",trait,"+",paste(covar2,collapse = "+"),"+",trait,"*",interaction)),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    } else {
      if (trait_is_cont==T){
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2,collapse= "+"),"+",colnames(dat)[i],"*",interaction)),design=subset(survey_design,ID%in%met_df_id))
      } else {
        model<- svyglm(as.formula(paste0(trait,"~",colnames(dat)[i],"+",paste(covar2,collapse = "+"),"+",colnames(dat)[i],"*",interaction)),design=subset(survey_design,ID%in%met_df_id),family=quasipoisson(link='log'))
      }
    }
    
    Vcov <- vcov(model, useScale = FALSE)
    beta<- coef(model)
    se<- sqrt(diag(vcov(model, useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[2])
    out_se[i] = round(as.numeric(se[2]),digits = 3)
    out_pvalue[i] = as.numeric(pval[2])
    out_nobs[i]=as.numeric(nobs(model))
    int_beta[i]=as.numeric(beta[length(beta)-1])
    int_se[i] = round(as.numeric(se[length(se)-1]),digits = 3)
    int_pvalue[i] = round(as.numeric(pval[length(pval)-1]),digits=9)
    cross_beta[i]=as.numeric(beta[length(beta)])
    cross_se[i] = round(as.numeric(se[length(se)]),digits = 3)
    cross_pvalue[i] = as.numeric(pval[length(pval)])
  }
  regress_output<-data.frame(metabolite=colnames(dat),
                             trait_beta=out_beta,
                             trait_se=out_se,
                             n=out_nobs,
                             trait_p_val=out_pvalue,
                             trait_p_val_Bonf=p.adjust(out_pvalue,method="bonferroni"),
                             trait_p_val_fdr=p.adjust(out_pvalue,method="BH"),
                             int_main_beta=int_beta,
                             int_main_se=int_se,
                             int_main_p_val=int_pvalue,
                             int_cross_beta=cross_beta,
                             int_cross_se=cross_se,
                             int_cross_p_val=cross_pvalue,
                             interaction_term=names(beta)[length(beta)-1],
                             cross_term=names(beta)[length(beta)]
                             
  )%>%
    dplyr::mutate(
      pval_Bonf_neglog=-log10(trait_p_val_Bonf),
      pval_fdr_neglog=-log10(trait_p_val_fdr),
      sig_Bonf=ifelse(trait_p_val_Bonf<0.05,"Bonf<0.05","Not Sig"),
      sig_fdr=ifelse(trait_p_val_fdr<0.05,"FDR<0.05","Not Sig"),
      is_continuous=metab_is_cont
    )
  identified_data<-regress_output%>%
    dplyr::filter(!stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=p.adjust(trait_p_val,method="BH"),
                  sig_fdr_named=ifelse(p_val_fdr_named<0.05,"FDR<0.05","Not Sig"))
  unidentified_data<-regress_output%>%
    dplyr::filter(stringr::str_detect(metabolite,"^X"))%>%
    dplyr::mutate(p_val_fdr_named=NA,
                  sig_fdr_named=NA)
  regress_output<-rbind(identified_data,unidentified_data)
  if (metab_is_complete==T){
    regress_output$is_continuous<-"Complete-cases"
  } else {
    regress_output$is_continuous<-factor(regress_output$is_continuous,levels=c(TRUE,FALSE),labels=c("Continuous","Dichotomized"))
  }
  
  return(regress_output)
}

# LASSO Regression
lasso_function<-function(x,y,binary_outcome,var_not_penalized,standardized){
  # Note alpha=1 for lasso only and can blend with ridge penalty down to
  # alpha=0 ridge only.
  penalize_vector<-data.frame(factor=rep(1,ncol(x)),predictor=colnames(x))
  penalize_vector$penalize<-ifelse(penalize_vector$predictor%in%var_not_penalized,0,1)
  
  if (binary_outcome==TRUE) {
    glmmod<-glmnet(x,y,alpha=1,family="binomial",penalty.factor = penalize_vector$penalize,standardize=standardized)
    lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
    set.seed(123)
    glmmod_cv <- cv.glmnet(x, y, alpha = 1, type.measure="auc",lambda = lambdas_to_try, nfolds = 10,family="binomial",penalty.factor = penalize_vector$penalize,standardize=standardized)
    # Plot cross-validation results
    plot(glmmod_cv)
    plot_cv<-recordPlot()
    best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
    mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
    se1_lambda<-glmmod_cv$lambda.1se
    plot(glmmod, xvar="lambda")
    abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
    abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
    plot_coef_lambda<-recordPlot()
    plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
      dplyr::mutate(row=row.names(.))%>%
      dplyr::rename(value="s1")%>%
      dplyr::filter(row != "(Intercept)"&value!=0) %>%
      ggplot(aes(value, reorder(row, value), color = value > 0)) +
      geom_point(show.legend = FALSE) +
      ggtitle("Influential variables") +
      xlab("Coefficient") +
      ylab(NULL) 
    coef<-coef(glmmod_cv, s = "lambda.min")
    coef_lasso_all<-coef
    coef_lasso<-coef[which(coef[,1]!=0),1]
    if (length(coef_lasso)==(length(var_not_penalized)+1)){
      set.seed(1234)
      glmmod_cv <- cv.glmnet(x, y, alpha = 1, type.measure="auc",lambda = lambdas_to_try, nfolds = 10,family="binomial",penalty.factor = penalize_vector$penalize,standardize=standardized)
      # Plot cross-validation results
      plot(glmmod_cv)
      plot_cv<-recordPlot()
      best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
      mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
      se1_lambda<-glmmod_cv$lambda.1se
      plot(glmmod, xvar="lambda")
      abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
      abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
      plot_coef_lambda<-recordPlot()
      plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
        dplyr::mutate(row=row.names(.))%>%
        dplyr::rename(value="s1")%>%
        dplyr::filter(row != "(Intercept)"&value!=0) %>%
        ggplot(aes(value, reorder(row, value), color = value > 0)) +
        geom_point(show.legend = FALSE) +
        ggtitle("Influential variables") +
        xlab("Coefficient") +
        ylab(NULL) 
      coef<-coef(glmmod_cv, s = "lambda.min")
      coef_lasso_all<-coef
      coef_lasso<-coef[which(coef[,1]!=0),1]
    }
   output<-list(plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, 
                 best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,coef_lasso=coef_lasso,coef_lasso_all=coef_lasso_all,plot_influ=plot_influ)
  } else {
    glmmod<-glmnet(x,y,alpha=1,penalty.factor = penalize_vector$penalize,standardize=standardized)
    lambdas_to_try <- 10^seq(-3, 5, length.out = 100)
    set.seed(123)
    glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
    # Plot cross-validation results
    plot(glmmod_cv)
    plot_cv<-recordPlot()
    best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
    mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
    se1_lambda<-glmmod_cv$lambda.1se
    plot(glmmod, xvar="lambda")
    abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
    abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
    plot_coef_lambda<-recordPlot()
    plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
      dplyr::mutate(row=row.names(.))%>%
      dplyr::rename(value="s1")%>%
      dplyr::filter(row != "(Intercept)"&value!=0) %>%
      ggplot(aes(value, reorder(row, value), color = value > 0)) +
      geom_point(show.legend = FALSE) +
      ggtitle("Influential variables") +
      xlab("Coefficient") +
      ylab(NULL)+ theme(axis.text.y = element_text(size = 7))
    coef<-coef(glmmod_cv, s = "lambda.min")
    coef_lasso<-coef[which(coef[,1]!=0),1]
    coef_lasso_all<-coef
    if (length(coef_lasso)==(length(var_not_penalized)+1)){
      set.seed(1026)
      glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
      # Plot cross-validation results
      plot(glmmod_cv)
      plot_cv<-recordPlot()
      best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
      mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
      se1_lambda<-glmmod_cv$lambda.1se
      plot(glmmod, xvar="lambda")
      abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
      abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
      plot_coef_lambda<-recordPlot()
      plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
        dplyr::mutate(row=row.names(.))%>%
        dplyr::rename(value="s1")%>%
        # broom::tidy() %>%
        dplyr::filter(row != "(Intercept)"&value!=0) %>%
        ggplot(aes(value, reorder(row, value), color = value > 0)) +
        geom_point(show.legend = FALSE) +
        ggtitle("Influential variables") +
        xlab("Coefficient") +
        ylab(NULL) 
      coef<-coef(glmmod_cv, s = "lambda.min")
      coef_lasso_all<-coef
      coef_lasso<-coef[which(coef[,1]!=0),1]
      if (length(coef_lasso)==(length(var_not_penalized)+1)){
        set.seed(10000)
        glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
        # Plot cross-validation results
        plot(glmmod_cv)
        plot_cv<-recordPlot()
        best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
        mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
        se1_lambda<-glmmod_cv$lambda.1se
        plot(glmmod, xvar="lambda")
        abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
        abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
        plot_coef_lambda<-recordPlot()
        plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
          dplyr::mutate(row=row.names(.))%>%
          dplyr::rename(value="s1")%>%
          # broom::tidy() %>%
          dplyr::filter(row != "(Intercept)"&value!=0) %>%
          ggplot(aes(value, reorder(row, value), color = value > 0)) +
          geom_point(show.legend = FALSE) +
          ggtitle("Influential variables") +
          xlab("Coefficient") +
          ylab(NULL) 
        coef<-coef(glmmod_cv, s = "lambda.min")
        coef_lasso_all<-coef
        coef_lasso<-coef[which(coef[,1]!=0),1]
        if (length(coef_lasso)==(length(var_not_penalized)+1)){
          set.seed(100)
          glmmod_cv <- cv.glmnet(x, y, alpha = 1, lambda = lambdas_to_try,nfolds = 10,penalty.factor = penalize_vector$penalize,standardize=standardized)
          # Plot cross-validation results
          plot(glmmod_cv)
          plot_cv<-recordPlot()
          best_lambda<-glmmod_cv$lambda.min # find the lambda that minimizes the CV error
          mean_lambda<-glmmod_cv$cvm #cvm:mean cross-validated error
          se1_lambda<-glmmod_cv$lambda.1se
          plot(glmmod, xvar="lambda")
          abline(v = log(glmmod_cv$lambda.min), col = "red", lty = "dashed")
          abline(v = log(glmmod_cv$lambda.1se), col = "blue", lty = "dashed")
          plot_coef_lambda<-recordPlot()
          plot_influ<-data.frame(as.matrix(coef(glmmod_cv, s = "lambda.min"))) %>%
            dplyr::mutate(row=row.names(.))%>%
            dplyr::rename(value="s1")%>%
            # broom::tidy() %>%
            dplyr::filter(row != "(Intercept)"&value!=0) %>%
            ggplot(aes(value, reorder(row, value), color = value > 0)) +
            geom_point(show.legend = FALSE) +
            ggtitle("Influential variables") +
            xlab("Coefficient") +
            ylab(NULL) 
          coef<-coef(glmmod_cv, s = "lambda.min")
          coef_lasso_all<-coef
          coef_lasso<-coef[which(coef[,1]!=0),1]
        }
      }
    }
    output<-list(plot_coef_lambda=plot_coef_lambda,plot_cv=plot_cv, best_lambda=best_lambda,mean_lambda=mean_lambda,se1_lambda=se1_lambda,
                 coef_lasso=coef_lasso,coef_lasso_all=coef_lasso_all,plot_influ=plot_influ)
    
  }
  return(output)
}

# function to calculate a summary score based on metabolite concentration and coefficient table
summary_score<-function(data, coeff){
  final<-data.frame()
  for (i in 1:nrow(coeff)){
    # print(i)
    component<-coeff[i,1]
    coefficient<-coeff[i,2]
    temp<-data.frame(metab=as.numeric(data[,colnames(data)%in%component]),coef=rep(coefficient,nrow(data)))
    temp$product<-temp$metab*temp$coef
    
    if (i==1){
      final<-data.frame(temp[,3])
    } else {
      final<-data.frame(cbind(final,temp[,3]))
    }
    colnames(final)[i]<-component
    # print(i)
  }
  if (ncol(final)==1) {
    score<-final
    colnames(score)[1]<-"summary"
  } else {score<-data.frame(summary=rowSums(final,na.rm = T))}
  return(score)
}
# function to export the coefficient and map to chemical names
export_lasso_coeff<-function(lasso_function_output,min_col_number,key_col,drop_covariates,metab_reference_table){
  if (length(lasso_function_output$coef_lasso)>3){
    lasso_coeff<-data.frame(lasso_function_output$coef_lasso)
    lasso_coeff$feature_id<-rownames(lasso_coeff)
    colnames(lasso_coeff)<-c("coeff",key_col)
    lasso_coeff<-lasso_coeff[which(!lasso_coeff[,key_col]%in%c("(Intercept)",drop_covariates)),] 
    lasso_coeff<-merge(lasso_coeff,metab_reference_table,by=key_col)
    return(lasso_coeff)
  }
}

# Loop over association analysis gender stratified model using the gender-combined MRS for incident htn
loop_association_analysis<-function(exposure_list,covariate_list,survey_design,outcome_var,dataset){
  output_mdl_list<-list()
  number=1
  model<-list()
  for (j in seq_along(exposure_list)){
    # print(paste("j=",j))
    for (k in seq_along(covariate_list)){
      # print(paste("k=",k))
      if (!is.na(covariate_list[k])){
        if (!is.na(exposure_list[j])){
          if (!is.null(survey_design)){
            if (length(unique(survey_design$variables[,outcome_var]))==2){
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),design = survey_design, family=quasipoisson(link='log'))
            } else {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),design = survey_design)
            }
          } else {
            if (length(unique(dataset[,outcome_var]))==2) {
              model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset, family=binomial(link='logit'))
            } else {
              model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset)
            }
          }
        } else {
          if (!is.null(survey_design)){
            if (length(unique(survey_design$variables[,outcome_var]))==2){
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),design = survey_design, family=quasipoisson(link='log'))
            } else {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),design = survey_design)
            }
          } else {
            if (length(unique(dataset[,outcome_var]))==2) {
              model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset, family=binomial(link='logit'))
            } else {
              model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset)
            }
          }
        }
      } else {
        if (!is.na(exposure_list[j])){
          if (!is.null(survey_design)){
            if (length(unique(survey_design$variables[,outcome_var]))==2) {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j])),design = survey_design, family=quasipoisson(link='log'))
            } else {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j])),design = survey_design)
            }
          } else {
            if (length(unique(dataset[,outcome_var]))==2){
              model<-glm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset, family=binomial(link='logit'))
            } else (
              model<-lm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset)
              
            )
          }
        } 
      }
      # print(paste("k=",k))
      output_mdl_list[[number]]<-model
      number=number+1
    }
    # print(paste("j=",j))
  }
  return(output_mdl_list)
}

# Loop over association analysis gender stratified model using the gender-combined MRS for incident htn
loop_association_analysis_offset<-function(exposure_list,covariate_list,survey_design,outcome_var,dataset,offset_var){
  output_mdl_list<-list()
  number=1
  model<-list()
  for (j in seq_along(exposure_list)){
    # print(paste("j=",j))
    for (k in seq_along(covariate_list)){
      # print(paste("k=",k))
      if (!is.na(covariate_list[k])){
        if (!is.na(exposure_list[j])){
          if (!is.null(survey_design)){
            if (length(unique(survey_design$variables[,outcome_var]))==2){
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design, family=quasipoisson(link='log'))
            } else {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design)
            }
          } else {
            if (length(unique(dataset[,outcome_var]))==2) {
              model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset, family=binomial(link='logit'))
            } else {
              model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+",exposure_list[j])),data = dataset)
            }
          }
        } else {
          if (!is.null(survey_design)){
            if (length(unique(survey_design$variables[,outcome_var]))==2){
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+offset(log(",offset_var,"))")),design = survey_design, family=quasipoisson(link='log'))
            } else {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"),"+offset(log(",offset_var,"))")),design = survey_design)
            }
          } else {
            if (length(unique(dataset[,outcome_var]))==2) {
              model<-glm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset, family=binomial(link='logit'))
            } else {
              model<-lm(as.formula(paste0(outcome_var,"~",paste0(covariate_list[[k]],collapse="+"))),data = dataset)
            }
          }
        }
      } else {
        if (!is.na(exposure_list[j])){
          if (!is.null(survey_design)){
            if (length(unique(survey_design$variables[,outcome_var]))==2) {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design, family=quasipoisson(link='log'))
            } else {
              model<-survey::svyglm(as.formula(paste0(outcome_var,"~",exposure_list[j],"+offset(log(",offset_var,"))")),design = survey_design)
            }
          } else {
            if (length(unique(dataset[,outcome_var]))==2){
              model<-glm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset, family=binomial(link='logit'))
            } else (
              model<-lm(as.formula(paste0(outcome_var,"~",exposure_list[j])),data = dataset)
              
            )
          }
        } 
      }
      # print(paste("k=",k))
      output_mdl_list[[number]]<-model
      number=number+1
    }
    # print(paste("j=",j))
  }
  return(output_mdl_list)
}

# extract model summary in the format of a data frame after loop_association_analysis
extract_mdl_sig_output<-function(model_output_list,model_index_matrix,model_outcome,model_strata) {
  out_beta<-rep(NA,length(model_output_list))
  out_se<-rep(NA,length(model_output_list))
  out_pvalue<-rep(NA,length(model_output_list))
  out_nobs<-rep(NA,length(model_output_list))
  out_index<-rep(NA,length(model_output_list))
  out_trait<-rep(NA,length(model_output_list))
  out_model<-rep(NA,length(model_output_list))
  for (i in seq_along(model_output_list)){
    # print(i)
    Vcov <- vcov(model_output_list[[i]], useScale = FALSE)
    beta<- coef(model_output_list[[i]])
    se<- sqrt(diag(vcov(model_output_list[[i]], useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    out_beta[i]=as.numeric(beta[length(beta)])
    out_se[i] = round(as.numeric(se[length(se)]),digits = 3)
    out_pvalue[i] = as.numeric(pval[length(pval)])
    out_nobs[i]=as.numeric(nobs(model_output_list[[i]]))
    out_index[i]=paste(i)
    out_trait[i]=colnames(model_index_matrix)[ceiling(i/nrow(model_index_matrix))]
    if (i%%nrow(model_index_matrix!=0)){
      out_model[i]=rownames(model_index_matrix)[i%%nrow(model_index_matrix)]
    } else {
      out_model[i]=rownames(model_index_matrix)[nrow(model_index_matrix)]
    }
    # print(i)
  }
  
  extracted_output<-data.frame(trait=out_trait,
                               model=out_model,
                               index=out_index,
                               beta=out_beta,
                               or=exp(out_beta),
                               se=out_se,
                               lower95=out_beta-1.96*out_se,
                               upper95=out_beta+1.96*out_se,
                               lowerCI=exp(out_beta-1.96*out_se),
                               upperCI=exp(out_beta+1.96*out_se),
                               n=out_nobs,
                               p_val=out_pvalue,
                               outcome=model_outcome,
                               strata=model_strata
  )%>%
    dplyr::mutate(p_label=dplyr::case_when(
      p_val > 0.05 ~ "",
      p_val > 0.01 ~ "*",
      p_val > 0.001 ~ "**",
      !is.na(p_val) ~ "***",
      TRUE ~ NA_character_
    ))
}

# extract model output of quartile regression of multiple exposures and models
quartile_association_models_output_offset<-function(model_list,dataset){
  qt2_beta<-rep(NA,length(model_list))
  qt2_se<-rep(NA,length(model_list))
  qt2_pvalue<-rep(NA,length(model_list))
  qt3_beta<-rep(NA,length(model_list))
  qt3_se<-rep(NA,length(model_list))
  qt3_pvalue<-rep(NA,length(model_list))
  qt4_beta<-rep(NA,length(model_list))
  qt4_se<-rep(NA,length(model_list))
  qt4_pvalue<-rep(NA,length(model_list))
  exposure<-rep(NA,length(model_list))
  covariates<-rep(NA,length(model_list))
  number<-1
  for (i in 1:length(model_list)){
    Vcov <- vcov(model_list[[i]], useScale = FALSE)
    beta<- coef(model_list[[i]])
    se<- sqrt(diag(vcov(model_list[[i]], useScale = FALSE)))
    zval<- beta / se
    pval<- 2 * pnorm(abs(zval), lower.tail = FALSE)
    qt2_beta[number]=as.numeric(beta[length(beta)-2])
    qt2_se[number] = as.numeric(se[length(beta)-2])
    qt2_pvalue[number] = as.numeric(pval[length(beta)-2])
    qt3_beta[number]=as.numeric(beta[length(beta)-1])
    qt3_se[number] = as.numeric(se[length(beta)-1])
    qt3_pvalue[number] = as.numeric(pval[length(beta)-1])
    qt4_beta[number]=as.numeric(beta[length(beta)])
    qt4_se[number] = as.numeric(se[length(beta)])
    qt4_pvalue[number] = as.numeric(pval[length(beta)])
    exposure[number]=substr(names(beta)[length(names(beta))],1,nchar(names(beta)[length(names(beta))])-1)
    predictors=paste(model_list[[i]]$formula[3])
    # covariates[number]=stringr::str_extract(predictors,".*(?=[+])")
    pattern <- "^(.*?)\\+(.*?\\+[^+]+)$"
    covariates[number]=sub(pattern, "\\1", predictors)
    if (covariates[number] == predictors) {
      covariates[number] <- "None"
    }
    number<-number+1
  }
  qt2_output<-data.frame(beta=qt2_beta,se=qt2_se,pval=qt2_pvalue,Quartile=rep("2nd",length(model_list)),Exposure=exposure,Covariate=covariates)
  qt3_output<-data.frame(beta=qt3_beta,se=qt3_se,pval=qt3_pvalue,Quartile=rep("3rd",length(model_list)),Exposure=exposure,Covariate=covariates)        
  qt4_output<-data.frame(beta=qt4_beta,se=qt4_se,pval=qt4_pvalue,Quartile=rep("4th",length(model_list)),Exposure=exposure,Covariate=covariates)        
  qt_output<-rbind(qt2_output,qt3_output,qt4_output)%>%
    dplyr::mutate(OR=exp(beta),
                  lowerCI=exp(beta-1.96*se),
                  upperCI=exp(beta+1.96*se),
                  p_label=dplyr::case_when(
                    pval > 0.05 ~ "",
                    pval > 0.01 ~ "*",
                    pval > 0.001 ~ "**",
                    !is.na(pval) ~ "***",
                    TRUE ~ NA_character_
                  ),
                  # model=rep(c("Model 1","Model 2","Model 3"),length(model_list)),
                  dataset=paste(dataset))
  model_names<-paste0("model_",seq_along(unique(qt_output$Covariate)))
  qt_output<-qt_output%>%
    dplyr::mutate(model=rep(model_names,(nrow(qt_output)/length(model_names))))
}   
