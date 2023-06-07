##### Representative Script #####

# ----------------------------------
# ----------------------------------
# The following was intended to give an overview of the steps included in the processing and analysis of the
# data across the project.

## Note the data wrangling/transformation that takes place between the following steps
## Not all wrangling steps were included to maintain clarity of the key processing and analysis steps 
# ----------------------------------
# ----------------------------------

# Libraries -----
library(tidyverse)
library(mosaic)
library(ggplot2)
library(ggpubr)
library(caret)
library(tidymodels)
library(ranger)
library(neuralnet)
library(mboost)
library(qvalue)
library(limma)
library(parallel)
library(doParallel)
library(e1071)
library(purrr)
library(kableExtra)
library(factoextra)
library(xgboost)
library(naivebayes)
library(kknn)
library(naivebayes)
library(patchwork)
library(discrim)
library(klaR)
library(dplyr)
library(reshape2)
library(ROCR) 
library(pROC)
library(grid)
library(gridExtra)
library(MLeval)
library(randomForest)
library(data.table)
library(ggrepel)
library(DALEXtra)
library(tibble)
library(explore)

### Functions and theme -----
# theme
theme_set(theme_minimal())

# parallel processing -----
doParallel::registerDoParallel(cores = 4)

#Helpful functions ------
`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)
back_2_front <- function(df){
  df <- df[,c(ncol(df),1:(ncol(df)-1))]
  return(df)
}
flextable_only <- function(table){
  table_flex <- flextable::flextable(table)%>%
    flextable::fontsize(size = 8)%>%
    flextable::autofit()%>%
    flextable:: height( height = .5)%>%
    flextable::theme_vanilla() %>% 
    flextable::bold(j=1)
}

# find out RAM usage
objects <- data.frame('object' = ls()) %>% 
  dplyr::mutate(size_unit = object %>%
                  sapply(. %>% get() %>%
                           object.size %>%
                           format(., unit = 'auto')),
                size = as.numeric(sapply(strsplit(size_unit, split = ' '), 
                                         FUN = function(x) x[1])),
                unit = factor(sapply(strsplit(size_unit, split = ' '), 
                                     FUN = function(x) x[2]), 
                              levels = c('Gb', 'Mb', 'Kb', 'bytes'))) %>% 
  dplyr::arrange(unit, dplyr::desc(size)) %>% 
  dplyr::select(-size_unit)


# limma's differential abundance functions
# differential abundance
limma_fun <- function(matrix_AB, no., var1, var2){
  Group <- factor(colnames(matrix_AB), levels = c(`var1`, `var2`))
  design <- model.matrix (~Group)
  colnames(design) <- c('var1', 'var1vsvar2')
  eset <- matrix_AB
  fit <- lmFit(eset, design)
  fit <- eBayes(fit)
  toptable <- topTable(fit, coef = 'var1vsvar2', adjust = 'BH', number = no.)
  toptable <- as.data.frame(toptable)
  toptable$Peak_ID_Ion <- rownames(toptable)
  toptable <- toptable[,c(ncol(toptable),1:(ncol(toptable)-1))]
  #toptable <- inner_join(peak_metadata_sim, toptable, by='Peak_ID_Ion')
  toptable$Sig <- 0
  toptable$Sig <- ifelse(toptable$adj.P.Val <0.05, '< 0.05', '> 0.05') 
  toptable$Sig_Names <-0
  toptable$Sig_Names <- ifelse(toptable$Sig =='< 0.05' ,toptable$Peak_ID_Ion, '')
  return(toptable)
}

### Set up example dataset ----
# example dataframe
omic_df <- data.frame(ID = paste0('Gene_', 1:200),
                      hc_1= sample(0:500, size=200, replace=T),
                      hc_2= sample(0:600, size=200, replace=T),
                      hc_3= sample(0:700, size=200, replace=T),
                      sepsis_1= sample(500:3300, size=200, replace=T), # I've made these obviously different to show an 'ideal' exmaple
                      sepsis_2= sample(500:2000, size=200, replace=T),
                      sepsis_3= sample(500:2200, size=200, replace=T),
                      gout_1 = sample(1000:2000, size=200, replace=T),
                      gout_2 = sample(1000:2800, size=200, replace=T),
                      gout_3 = sample(1000:3000, size=200, replace=T),
                      AS_1 = sample(0:2000, size=200, replace=T),
                      AS_2 = sample(0:2500, size=200, replace=T),
                      AS_3 = sample(0:1500, size=200, replace=T))

# increase number of 0s as example
omic_df[2:ncol(omic_df)][omic_df[2:ncol(omic_df)] < 100] <- 0
rownames(omic_df) <-  omic_df$ID
omic_df <- omic_df[,-1]
omic_df <- as.data.frame(t(omic_df))
omic_df$Samples <- rownames(omic_df)
omic_df <- back_2_front(omic_df)

## Processing data ##
# Removing features with high proportion (80%) of samples with NAs #
number_samples <- length(omic_df$Samples)
omic_df <- as.data.frame(t(omic_df[-1]))
omic_df$Prop_Zero <- (rowSums(is.na(omic_df))/number_samples)*100
omic_df <- back_2_front(omic_df)

omic_df <- omic_df %>%
  filter(omic_df$Prop_Zero < 20)

omic_df <- omic_df[,-1]
omic_df <- as.matrix(omic_df)

## Imputing for missing values (missing not at random) using half minimum imputaion ##
# imputation function: Take the missing values for each feature 
# where these are expected to be missing due to levels below quantification #
# replacezero <- function(x) "[<-"(x, !x|is.na(x), 
#                                  min(x[x>0], 
#                                      na.rm=TRUE)/runif(1,min=2))
# 
# 
# omic_df <- as.data.frame(apply(omic_df, 1, replacezero)) # apply the function to the omic dataframe

## Remove batch effects ## ----
# batch <- samples_df$Batch
# omic_df <- removeBatchEffect(omic_df, batch=batch, covariates=NULL)
# 
# ## Normalise using cyclic loess function where high variation observed across features in dataset
# ## Metabolomic data ##
# omic_df <- normalizeCyclicLoess(omic_df)

### Transcriptomic data ### ----
## Normalisation of transcriptomic data ##
# BSData = readBeadSummaryData(dataFile = dataFile,
#                              qcFile = qcFile, controlID = "ProbeID",
#                              skip = 0, qc.skip = 0, qc.columns = list(exprs = "AVG_Signal",
#                                                                       Detection = "Detection Pval"))
# BSData.quantile = normaliseIllumina(BSData,
#                                     method = "quantile", 
#                                     transform = "log2")
# omic_df = data.frame(exprs(BSData.quantile),
#                  type = "Quantile",
#                  allpositive = allpositive)
# 
# ## Aggregate probes using limma's avereps function
# omic_df <-avereps(omic_df, ID=omic_df$Probe_ID)

### Exploratory Data Analysis ### -----
## Principal Component Analysis ## 
omic_df <- omic_df %>% 
  filter(rownames(omic_df) %notlike% 'Prop') 

pca_omic <- omic_df # ensure numeric only
scaled_intensities <- scale(pca_omic)
scaled_intensities[do.call(cbind, lapply(scaled_intensities, is.nan))] <- 0
scaled_intensities<- as.data.frame(scaled_intensities)
pca_data <- prcomp(scaled_intensities)
pca_coord <- data.frame(pca_data$x)
var_explained <- pca_data$sdev^2/sum(pca_data$sdev^2)
var_explained[1:5] # observe amount of variance for first 5 PCs

# generate scree plot
scree <- fviz_eig(pca_data) 

# generate PCA plot
pca_coord %>% 
  mutate(Sample = as.factor(gsub("(.*)_.*","\\1",  rownames(pca_coord)))) %>% 
  ggplot() + 
  geom_point(size=2, alpha=0.7, 
             aes(x=PC1,y=PC2, colour= Sample, fill= Sample))+
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       colour='Sample Type',
       fill= 'Sample Type')+
  geom_hline(yintercept = 0,
             colour='navy',
             linetype='dashed')+
  geom_vline(xintercept = 0,
             colour='navy',
             linetype='dashed')+
  theme_minimal()+
  theme(axis.text= element_text(size=14),
        axis.title= element_text(size=16, face='bold'),
    legend.key.size = unit(01, 'cm'), 
        legend.key.height = unit(1, 'cm'),
        legend.key.width = unit(1, 'cm'), 
        legend.title = element_text(size=16, face='bold'),
        legend.text = element_text(size=14)) 

## Differential Analysis ##
# visualise the number of samples in each condition of interest

omic_df$Condition <- as.factor(gsub("(.*)_.*","\\1",  rownames(omic_df)))
omic_df%>% 
  ggplot(aes(x=Condition, fill=Condition)) +
  geom_bar() 

# differential analysis function using the limma package
omic_df <- omic_df %>% 
  filter(rownames(omic_df) %notlike% 'Prop') 
number_of_features <- length(prelimma_omic[1])

limma_omic_DA <- limma_fun(prelimma_omic, number_of_features, 'Condition_1', 'Condition_2')

# volcano plot
limma_omic_DA%>%
  ggplot(aes(x=logFC, y=-log10(P.Value))) +
  geom_point (size=3,alpha=0.7,
              aes(colour=Sig, 
                  group=Sig)) +
  theme_minimal() +
  labs (x='LogFC',
        y='-Log p-value',
        colour='Unadjusted \np-value')+
  geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Names_2),
                  box.padding =1,
                  size=5,
                  max.overlaps = Inf,
                  position = position_jitter(seed = 1),
                  arrow = arrow(length = unit(0.0035, "npc"))) +
  geom_vline(xintercept = c(-0.2,0.2), linetype='dashed')+
  geom_hline(yintercept = -log10(0.05), linetype='dashed')+
  geom_text(aes(x=-1.5, y= -log10(0.05)), 
            label='unadjusted p-value = 0.05',
            vjust=2.1, colour='dark grey')+
  
  theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))+
  scale_color_brewer(palette = "Set1",direction=-1)

## Correlation analysis ##
ints_nested <- omics_df %>%
  group_by(Feature) %>%
  nest()

ints_unnested <- omics_df %>%
  unnest(cols=c())
identical(omics_df, ints_unnested)
ints_lm <- ints_nested %>%
  mutate(model = purrr::map(data, ~lm(formula = Feature_Level~Disease_Activity, data = .x)))
model_coef_nested <- ints_lm %>%
  mutate(coef=map(model, ~tidy(.x)))
model_coef <- model_coef_nested %>%
  unnest(coef)
model_perf_nested <- ints_lm %>%
  mutate(fit=map(model, ~glance(.x)))
model_perf <- model_perf_nested%>%
  unnest(fit)
best_fit <- model_perf %>%
  top_n(n=4, wt=r.squared)
bestest_fit <- with(model_perf,model_perf[order(-r.squared),])
best_augmented <- bestest_fit %>% 
  mutate(augmented = map(model, ~augment(.x))) %>% 
  unnest(augmented) 
best_augmented$adj_p <- p.adjust(best_augmented$p.value, method='BH')

# Plot the features that correlated with disease activity measure with adjusted p-value < 0.05
# if no features shown then use unadjusted p-value < 0.05 with caveat of non significance.
best_augmented %>%
  filter(adj_p < 0.05) %>% 
  ggplot(aes(x = Disease_Activity, y=Feature_Level)) +
  geom_point(size=1, alpha=0.7) + 
  stat_cor(vjust=1, hjust=0,
           size=5)+
  geom_smooth(method='lm',
              colour='red')+
  facet_wrap(~Feature, scales = "free_y")+
  theme(
    strip.text.x= element_text(#face = "bold",
      size=16),
    title=element_text(size=16),
    axis.title = element_text(size=14,face='bold'),
    axis.text = element_text(size=12))


## Supervised machine learning ##
# Generating an omic profile associated with the classification of samples across conditions of interest

# For holdout cross validation, split omic_df into training and testing subsets
# example below uses 70:30 train:test split #
set.seed(42)
index <- createDataPartition(omic_df$Condition, p = 0.7, list = FALSE)
train_data <- omic_df[index, ]
test_data  <- omic_df[-index, ]

# Feature selection using recursive feature elimination #
# stringent function
ft_sel <- function(train_df){
  for (i in 1:10){ # run the RFE process 10 times  
    profile_function <- function(Profile_number){
      options(warn=-1)
      subsets <- c(1:10) 
      set.seed(42)
      # 10-fold cross validation repeated 10 times for each run
      ctrl <- rfeControl(functions = rfFuncs,
                         method = "repeatedcv",
                         number = 10, 
                         repeats=10,
                         verbose = FALSE)
      
      profile <- rfe(x=train_df[,-1], y=train_df$Condition,
                     sizes = subsets,
                     rfeControl = ctrl)
      
      profile <- profile$variables
      profile <- profile %>%
        arrange(-Overall) %>%
        distinct(var, .keep_all=TRUE) %>%
        filter(Overall > 1) # filter features with relative importance > 1
      
      profile_df <- profile %>%
        group_by(var) %>%
        nest() %>%
        mutate(Profile=1) 
    }
    profile <- profile_function(i)
    assign(paste0("profile_test", i), profile)
  }
  df <- do.call(rbind, mget(ls(pattern = "profile_test")))
  
  df_2 <- df %>%
    group_by(var) %>%
    filter(n() >= 5) %>% # from each of the runs, select features that are included in at least half of the runs
    unnest()
  df_2 <- df_2 %>%
    aggregate(.~var, mean, na.rm=TRUE) %>%
    mutate(Round_overall=round(Overall,1)) %>%
    arrange(-Overall)
}

ft_sel_results <- ft_sel(train_data)

ft_sel_results %>%
  ggplot(aes(x=Overall, y=reorder(var, Overall)))+
  geom_col(aes(fill=Overall))+
  geom_text(aes(label = round(Overall, 2)),
            vjust=0.8,
            hjust=1.2,
            color="white",
            size=6) +
  #scale_color_colorblind()+
  scale_fill_continuous(low='light green', high='navy')+
  theme_minimal()+
  theme(legend.position = 'none',
        axis.title.y = element_text(size = 20),
        axis.title.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        title=element_text(size=22)
  )+
  labs(x='Relative Importance',
       y='Feature')

# Keep only selected features
train_data_2 <- train_data %>% 
  dplyr::select(Condition, ft_sel_results$var) %>%
  as.data.frame()

# Algorithm Selection
# Function for algorithm tuning and selection
multi_mods_test <- function(){
  control <- trainControl(method="repeatedcv",
                          number=10,
                          repeats=10,
                          summaryFunction=twoClassSummary,
                          savePredictions = TRUE,
                          classProbs = TRUE,
                          verboseIter = TRUE,
                          search = 'random')
  # train the SVM model
  set.seed(42)
  modelSvm <- caret::train(Condition~., data=train_data_2,
                           method='svmRadial',
                           metric='ROC',
                           tuneLength=10,
                           trControl= control)
  
  # train the LogM model
  set.seed(42)
  modelglm <- caret::train(Condition~., data=train_data_2,
                           method='glm',
                           metric='ROC',
                           tuneLength=10,
                           trControl= control)
  # train the LogM model
  set.seed(42)
  modelglm_boost <- caret::train(Condition~., data=train_data_2,
                                 method='glmboost',
                                 metric='ROC',
                                 tuneLength=10,
                                 trControl= control)
  # train the RF model
  set.seed(42)
  model_rf <- caret::train(Condition~., data=train_data_2,
                           method='ranger',
                           metric='ROC',
                           tuneLength=10,
                           trControl= control)
  # train the XGBoost model
  set.seed(42)
  model_xgb <- caret::train(Condition~., data=train_data_2,
                            method='xgbTree',
                            metric='ROC',
                            tuneLength=10,
                            trControl= control)
  # train the KNN model
  set.seed(42)
  model_knn <- caret::train(Condition~., data=train_data_2,
                            method='kknn',
                            metric='ROC',
                            tuneLength=10,
                            trControl= control)
  # train the naive Bayes model
  set.seed(42)
  model_naivebayes <- caret::train(Condition~., data=train_data_2,
                                   method='naive_bayes',
                                   metric='ROC',
                                   tuneLength=10,
                                   trControl= control)
  
  
  # collect resamples
  results <- resamples(list(SVM = modelSvm, RF= model_rf, LRM=modelglm,GLMB = modelglm_boost,
                            XGBoost =model_xgb, KNN=model_knn,
                            Naive.Bayes= model_naivebayes)) #LRM= modelLRM))
  
  
  # summarize the distributions
  summary(results)
  # boxplots of results
  plott <- bwplot(results)
  comp_roc <- evalm(list(modelSvm, 
                         model_rf, 
                         modelglm, 
                         modelglm_boost,
                         model_xgb,  
                         model_knn,
                         model_naivebayes),
                    gnames=c('KNN', 'NB', 'RF','LRM','GLMB', 'SVM', 'XGB'))
  
  ml_eval_output <- as.data.frame(comp_roc$stdres)
  ml_eval_output$Measure <- rownames(ml_eval_output)
  ml_eval_output <- back_2_front(ml_eval_output)
  ml_eval_output <- flextable_only(ml_eval_output)
  listt <- list(comp_roc, ml_eval_output)
  return(listt)
}

# Compare algorithms
alg_sel <- multi_mods_test()
met_table <- alg_sel[[2]]
met_roc <- alg_sel[[1]]$roc
met_prg <- alg_sel[[1]]$prg


# Select final model and increase number of repeats in cross validation to 100
# Model generation
set.seed(42)
fit_control <- trainControl(method="repeatedcv",
                            number=10,
                            repeats=100,
                            summaryFunction=twoClassSummary,
                            savePredictions = TRUE,
                            classProbs = TRUE,
                            verboseIter = TRUE,
                            search = 'random')

set.seed(42)
omic_model <- caret::train(Response~., data=train_data_2,
                           method='xgbTree',
                           metric='ROC',
                           tuneLength=15,
                           trControl= fit_control)

ggplot(omic_model) # check the trained model's initial performance

# Evaluate model's performance in the training subset using the ML_eval package
train_results <- evalm(omic_model)
train_roc <- train_results$roc
train_prg <- train_results$prg

# Evaluate final tuned model in testing subset
# testing subset
# function
model_performance <- function(model, test_df){
  predictions <- predict(model, test_df)
  confusionMatrix(predictions, test_df$Condition)
  con_matr <- confusionMatrix(predictions, test_df$Condition)
  con_stats <- con_matr$overall
  
  pr <- prediction(as.numeric(predictions), as.numeric(test_df$Condition))
  prf <- performance(pr, measure = "tpr", x.measure = "fpr")
  auc <- performance(pr, measure = "auc")
  auc_val <- auc@y.values[[1]]
  result.predicted.prob <- predict(model, test_df, type="prob") # Prediction
  result.roc <- roc(test_df$Good, result.predicted.prob$Good) # Apply the positive class e.g. a 'Good' response
  list_pred <- list(model, con_stats, result.roc, con_matr)
  return(list_pred)
}

test_performance <- model_performance(omic_model, test_data)
test_roc <- test_performance[[3]]
test_roc$auc

# Feature Interpretation using the DALEX and DALEXtra packages
model_explainer <- DALEX::explain(model = omic_model,  
                              data = omic_df[, -1],
                              y = omic_df$Condition, 
                              type='classification')

## attempt accumulated local effects plots to explain interaction of features in the plot
pdp_model <- model_profile(explainer = model_explainer,
                       type = "partial",
                       variables = names(omic_df)[2:ncol(omic_df)])
ld_model <- model_profile(explainer = model_explainer,
                       type       = "conditional",
                       variables = names(omic_df)[2:ncol(omic_df)])
ale_model <- model_profile(explainer = model_explainer, 
                        type='accumulated',
                        variables = names(omic_df)[2:ncol(omic_df)])

shapleys_model <- predict_parts(explainer = model_explainer, 
                             type='shap',
                             new_observation =  omic_df, 
                             B=25) # the B=25 is the number of random orderings of the explanatory variables

pdp_model$agr_profiles$`_label_` = "partial dependence"
ld_model$agr_profiles$`_label_` = "local dependence"
ale_model$agr_profiles$`_label_` = "accumulated local"

partials <- plot(pdp_model, ld_model, ale_model)

breakdown <- model_explainer %>% 
  predict_parts(new_observation = omic_df) %>% 
  plot()

SHAP <- plot(shapleys_rf)

# Explain the variance
# include only the selected features
omic_df_cut <- omic_df %>% 
  dplyr::select(c(Sample, ft_sel_results$var))

# Select patient factors of interest
patient_factors <- patient_metadata%>%
  dplyr::select(Sample,Age, Sex, aCCP.Status, RhF.Status, Smoking.Status )

form <- ~ Age + (1|Sex) + (1|aCCP.Status) + (1|RhF.Status) + (1|Smoking.Status) 

varPart <- fitExtractVarPartModel(omic_df_cut, form, patient_factors)

vp <- sortCols(varPart)
explain_each <- plotPercentBars( vp[] ) +
  theme(text=element_text(size=16))
var_plot <- plotVarPart( vp )+
  theme(text=element_text(size=16))


# Model comparison
# generate the above objects for each model of interest
test_roc_model_1$auc
test_roc_model_2$auc
test_roc_model_3$auc

roc_comp <- rbind.data.frame(test_roc_model_1,
                             test_roc_model_2, 
                             test_roc_model_3)

# Plot the ROC curves for each model as a single figure
roc_comparison <- roc_comp%>%
  ggplot(aes(x = One_Minus_Spec, y = Sensitivity, group_by=Model)) +
  geom_path(aes(colour=Model), size=1)+
  geom_abline(linetype='solid', colour='grey') +
  coord_equal()+
  theme_minimal()+
  labs(x='False positive rate',
       y='True positive',
       colour='')+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_rect(colour = "black", fill=NA, size=1),
        legend.text = element_text(size=16),
        axis.title = element_text(size = 16),
        axis.text = element_text(size = 16))

# Statistical testing for difference between models' performance
#delong

roc_comparison_1 <- roc.test(test_roc_model_1, test_roc_model_2)
roc_comparison_2 <- roc.test(test_roc_model_1, test_roc_model_3)
roc_comparison_3 <- roc.test(test_roc_model_2, test_roc_model_3)


roc_table_1 <- (cbind.data.frame(roc_comparison_1$roc2$auc,
                                 roc_comparison_1$roc1$auc,
                                 round(roc_comparison_1$p.value,3)))
names(roc_table_1) <- c('AUC-ROC 1', 'AUC-ROC 2', 'p-value')
roc_table_2 <- (cbind.data.frame(roc_comparison_2$roc2$auc,
                                 roc_comparison_2$roc1$auc,
                                 round(roc_comparison_2$p.value,3)))
names(roc_table_2) <- c('AUC-ROC 1', 'AUC-ROC 2', 'p-value')
roc_table_3 <- (cbind.data.frame(roc_comparison_3$roc2$auc,
                                 roc_comparison_3$roc1$auc,
                                 round(roc_comparison_3$p.value,3)))
names(roc_table_3) <- c('AUC-ROC 1', 'AUC-ROC 2', 'p-value')
roc_table_complete <- rbind.data.frame(roc_table_1, roc_table_2, roc_table_3)
rownames(roc_table_complete) <-  c('Model 1 vs Model 2',
                                   'Model 1 vs Model 3',
                                   'Model 2 vs Model 3')

plottt <- roc_table_complete %>%
  kable() %>%
  kable_styling(bootstrap_options=c("striped", "hover", "responsive"))

ggplot() +
  theme_void()+
  annotate(geom='table',
           x=1, y=1, 
           label=list(roc_table_complete))

# ----------------------------------
# ----------------------------------
# ----------------------------------
