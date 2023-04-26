rm(list = ls())

# Packages nedeed
list.of.packages <- c("RhpcBLASctl", "dplyr", "BGLR", "rBayesianOptimization", "devtools")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if(!("SKM" %in% installed.packages()[,"Package"])){
  devtools::install_github("brandon-mosqueda/SKM")
}

library(dplyr)
library(SKM)
library(RhpcBLASctl)
blas_set_num_threads(2)

# Cambiar la dirección de los datos 
#Dir <- "C:/Users/Osval/Documents/Berna/Practicas profesionales/Hyperparameter tuning/To Suecia"
Dir <- "/home/omontesinos/OAML/PLS"
setwd(Dir)

load("Pheno.Rdata", verbose = TRUE) #Pheno Data
load("Markers.Rdata", verbose = TRUE)
Pheno <- Y
Markers <- X

dim(Pheno)
# Number of rows with any missing value:
sum(apply(is.na(Pheno), 1, sum) != 0)
apply(is.na(Pheno), 2, sum)
str(Pheno)
head(Pheno)

# NA's
dim(Markers)
sum(apply(is.na(Markers), 2, sum) == 0)
rownames(Markers)
head(Markers[,1:10])
sum(apply(Markers, 2, typeof) == "double")

pos_na <- (apply(is.na(Pheno), 1, sum) != 0)
Pheno <- Pheno[!pos_na,]

# Markers
dim(Markers)
cbind(rownames(Markers), unique(Pheno$GID))

# Data preparation of Env, G & GE
Line <- model.matrix(~0 + GID, data = Pheno)
pos <- gsub("GID", "",colnames(Line))
Markers <- Markers[pos,] #each Line's column correspond to each Geno's row
Geno = tcrossprod(as.matrix(Markers))/dim(Markers)[2]
dim(Geno); sum(rownames(Geno) != pos)
Env <- model.matrix(~0 + site, data = Pheno)
Geno <- t(chol(Geno)) 
LineG <- Line %*% Geno
LinexGenoxEnv <- model.matrix(~ 0 + LineG:Env)

# Predictor Variable
X <- cbind(Env, LineG, LinexGenoxEnv)

Traits <- colnames(Pheno)[-c(1,2)]

for(j in 1:length(Traits)){
  # Identify the j-th Trait
  y <- Pheno[, Traits[j]]
  
  # Note that y is a continuous numeric vector
  class(y)
  typeof(y)
  
  # 5-Fold CV Partition
  set.seed(2022)
  folds <- cv_kfold(records_number = nrow(X), k = 5)
  
  # A data frame that will contain the variables:
  ## (Number) Fold, Line, Env, (testing values) Observed and Predicted (values)
  Predictions <- data.frame()
  
  # Model training and predictions of the ith partition
  for (i in seq_along(folds)) {
    cat("\n")
    cat("\t\t\t\t*** Fold:", i, " ***\n")
    fold <- folds[[i]]
    
    #Identify the training and testing sets
    X_training <- X[fold$training, ]
    X_testing <- X[fold$testing, ]
    y_training <- y[fold$training]
    y_testing <- y[fold$testing]
    
    # Model training
    model <- partial_least_squares(
      x = X_training,
      y = y_training
    )
    
    #Prediction of the testing set
    ncomp_line <- model$optimal_components_num
    predictions <- predict(model, X_testing, components_num = ncomp_line)
    #plot(y_testing, predictions$predicted)
    #plot(y_testing, predictions2$predicted)
    #cbind(predictions$predicted, predictions2$predicted)
    numeric_summary(y_testing, predictions$predicted)
    
    # Predictions for the Fold
    FoldPredictions <- data.frame(
      Fold = i,
      Line = Pheno$GID[fold$testing],
      Env = Pheno$site[fold$testing],
      Observed = y_testing,
      Predicted = predictions$predicted
    )
    
    Predictions <- rbind(Predictions, FoldPredictions)
    
    cat("Optimum Components number: ", ncomp_line)
  }
  
  NewDir <- paste("5-Fold CV", Traits[j], sep = "/")
  mkdir(NewDir)
  
  head(Predictions)
  unique(Predictions$Fold)
  
  # Summaries
  summaries <- gs_summaries(Predictions) 
  
  # Elements of summaries
  names(summaries)
  
  # Prediction of j-th Trait:
  data.table::fwrite(
    Predictions,
    file = paste(NewDir,"Predictions.csv", sep = "/"),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  
  # Summaries: j-th Trait
  ## by Enviroment
  data.table::fwrite(
    summaries$env,
    file = paste(NewDir,"Sum_Env.csv", sep = "/"),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  ## by line
  data.table::fwrite(
    summaries$line,
    file = paste(NewDir,"Sum_Line.csv", sep = "/"),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  ## by fold
  data.table::fwrite(
    summaries$fold,
    file = paste(NewDir,"Sum_Fold.csv", sep = "/"),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
}

############################### Multitrait ####################################
# Predictor Variable
X <- cbind(Env, LineG, LinexGenoxEnv)
y <- Pheno[, Traits]

# 5-Fold CV Partition
set.seed(2022)
folds <- cv_kfold(records_number = nrow(X), k = 5)

# A data frame that will contain the variables:
## (Number) Fold, Line, Env, (testing values) Observed and Predicted (values)
Predictions <- list()
for(i in names(y)){
  Predictions[[i]] <- data.frame()
}

# Model training and predictions of the ith partition
for (i in seq_along(folds)) {
  cat("\n")
  cat("\t\t\t\t*** Fold:", i, " ***\n")
  fold <- folds[[i]]
  
  #Identify the training and testing sets
  X_training <- X[fold$training, ]
  X_testing <- X[fold$testing, ]
  y_training <- y[fold$training,]
  y_testing <- y[fold$testing,]
  
  # Model training
  model <- partial_least_squares(
    x = X_training,
    y = y_training
  )
  
  #Prediction of the testing set
  ncomp_line <- model$optimal_components_num
  X_testing2 <- scale(X_testing, center=colMeans(X_testing)-colMeans(X_training), scale=F)
  predictions <- predict(model, X_testing2, components_num = ncomp_line, format = "data.frame")
  #plot(y_testing[, c("BLUE__starch")], predictions$BLUE__starch$predicted)
  #plot(y_testing[, c("BLUE__starch")], predictions2$BLUE__starch$predicted)
  #cbind(y_testing[, c("BLUE__starch")], predictions$BLUE__starch$predicted, predictions2$BLUE__starch$predicted)
  #numeric_summary(y_testing[, c("BLUE__40mm")], predictions$BLUE__40mm$predicted)
  #numeric_summary(y_testing[, c("BLUE__40mm")], predictions2$BLUE__40mm$predicted)
  
  # Predictions for the Fold: Bayesian Optimization Tuning
  for(j in colnames(y)){
    FoldPredictionsVi <- paste("FoldPredictions", j, sep = "")
    
    assign(FoldPredictionsVi,
           data.frame(
             Fold = i,
             Line = Pheno$GID[fold$testing],
             Env = Pheno$site[fold$testing],
             Observed = y_testing[, j],
             Predicted = predictions[, j]))
    Predictions[[j]] <- rbind(Predictions[[j]], get(FoldPredictionsVi))
  }
  cat("Optimum Components number: ", ncomp_line)
}

NewDir <- "5-Fold CV/Multitrait"
mkdir(NewDir)

# Generated files
for(k in colnames(y)){
  mkdir(paste(NewDir, k, sep = "/"))
  Predictions_Vi <- paste("Predictions", k, sep = "_")
  print(Predictions_Vi)
  
  # Prediction of Variable i
  data.table::fwrite(
    Predictions[[k]],
    file = paste(NewDir,"/",k,"/",Predictions_Vi,".csv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  
  # Summaries: No Tuning
  sum <- paste("Summaries", k, sep = "_")
  assign(sum, gs_summaries(Predictions[[k]]))
  ## by Enviroment
  data.table::fwrite(
    get(sum)$env,
    file = paste(NewDir,"/",k,"/",sum,"_Env.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  ## by line
  data.table::fwrite(
    get(sum)$line,
    file = paste(NewDir,"/",k,"/",sum,"_Line.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
  ## by fold
  data.table::fwrite(
    get(sum)$fold,
    file = paste(NewDir,"/",k,"/",sum,"_Fold.csv", sep = ""),
    row.names = FALSE,
    quote = FALSE,
    na = "NA"
  )
}
