# Copyright (C) 2022 The University of Edinburgh 
# Author Stuart Aitken MRC IGC s.aitken@ed.ac.uk
# All Rights Reserved.
# Funded by the Medical Research Council
# https://www.ed.ac.uk/mrc-human-genetics-unit

## Code for the IMPROVE HPO term classifier
## Requires files generated by the script IMPROVE_resource.r
## generates Fig 3 when run on DDD data; runs on demo resource data;

library(hash);
library(entropy);
library(e1071);
library(naivebayes);
library(GenSA);
library(gplots);
library(ROCR);

## demo data
source(file="IMPROVE_functions.r");
load(file="informativePhenotypicTerms.RData"); ## provided, see also IMPROVE_resource.r
load(file="hpou_demo.RData"); ## generated using IMPROVE_resource.r from database of HPO annotations
load(file="pheno_demo.RData");
proband_to_diagnosis <- data.frame(pheno);
proband_to_diagnosis$GENE <- factor(proband_to_diagnosis$GENE);

## data structures used by improve: diagns_to_probands and geneTotals
allgenes <- sort(unique(proband_to_diagnosis[,'GENE']));
length(allgenes); ## 2 for demo data; 77 for DDD data
diagns_to_probands <- vector('list',length=length(allgenes));
names(diagns_to_probands) <- allgenes;
for(i in 1:length(allgenes)) {
    diagns_to_probands[[i]] <- list(as.character(unique(proband_to_diagnosis[proband_to_diagnosis[,'GENE']==allgenes[i],'ID'])));
}
geneTotals <- table(proband_to_diagnosis$GENE);

## run model comparison and ROC analysis of improve_HPO() on proband*HPO matrix
## and improve_nBayes() on continuous and nominal data
## growth and development must be normalised continuous data
## gender data is nominal, all data input must have a column 'GENE' giving the diagnosis
resultsList <- vector('list',length=8);
names(resultsList) <- c('HPO_test_on_training','HPO_leave_one_out',
                        'Growth_test_on_training','Growth_leave_one_out',
                        'Development_test_on_training','Development_leave_one_out',
                        'Gender_test_on_training','Gender_leave_one_out');

resultsList[[1]] <- run_improve(proband_to_diagnosis,
                                geneTotals,
                                improve_HPO,  ## method to generate HPO gene models and compute likelihood
                                filename='HPO_testOnTraining_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=FALSE); ## test on training

resultsList[[2]] <- run_improve(proband_to_diagnosis,
                                geneTotals,
                                improve_HPO,  ## method to generate HPO gene models and compute likelihood
                                filename='HPO_LOO_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=TRUE); ## cross-validation

resultsList[[3]] <- run_improve(data_growth,
                                geneTotals,
                                improve_nBayes,  ## method to generate gene models from continuous data and compute likelihood
                                filename='Growth_testOnTraining_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=FALSE); ## test on training

resultsList[[4]] <- run_improve(data_growth,
                                geneTotals,
                                improve_nBayes,  ## method to generate gene models from continuous data and compute likelihood
                                filename='Growth_LOO_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=TRUE); ## cross-validation

resultsList[[5]] <- run_improve(data_development,
                                geneTotals,
                                improve_nBayes,  ## method to generate gene models from continuous data and compute likelihood
                                filename='Development_testOnTraining_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=FALSE); ## test on training

resultsList[[6]] <- run_improve(data_development,
                                geneTotals,
                                improve_nBayes,  ## method to generate gene models from continuous data and compute likelihood
                                filename='Development_LOO_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=TRUE); ## cross-validation

resultsList[[7]] <- run_improve(data_gender,
                                geneTotals,
                                improve_nBayes_nominal,  ## method to generate gene models from nominal data and compute likelihood
                                filename='Gender_testOnTraining_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=FALSE); ## test on training

resultsList[[8]] <- run_improve(data_gender,
                                geneTotals,
                                improve_nBayes_nominal,  ## method to generate gene models from nominal data and compute likelihood
                                filename='Gender_LOO_ROC', ## file name for ROC curve pdf
                                PLOT=TRUE,LOO=TRUE); ## cross-validation


## summarise performance [F1 and AUC presented in Figure 3]
performanceList <- vector('list',length=length(resultsList));
names(performanceList) <- names(resultsList);
for(i in 1:length(resultsList)) { 
    performanceList[[i]] <- mkQuantPerformance(resultsList[[i]],
                                               geneTotalsp=geneTotals,
                                               dataf=proband_to_diagnosis);
}
performanceList; ## demo HPO data to allow code to run only
'performanceList
$test_on_training
  TP FP prec    recall   F1 AUC
X  6  4  0.6 1.0000000 0.75   1
Y  2  0  1.0 0.3333333 0.50   1

$leave_one_out
  TP FP prec    recall  F1 AUC
X  4  6  0.4 0.6666667 0.5   0  ## ROC not meaningful here
Y  0  2  0.0 0.0000000 NaN   0'

## DDD study data
head(performanceList[[1]])
'        TP FP       prec    recall        F1       AUC
ADNP    10 45 0.18181818 0.2702703 0.2173913 0.8571782
AHDC1    5 58 0.07936508 0.2941176 0.1250000 0.9122798
ANKRD11 38 65 0.36893204 0.4691358 0.4130435 0.8855348
ARID1B  31 46 0.40259740 0.4366197 0.4189189 0.8804218';

save(performanceList,file='performanceList_IMPROVE_DD.RDATA');

## optimally combine test-on-training classifiers
priors_equalising <- log(table(data_growth$GENE)[rownames(geneTotals)]/dim(data_growth)[1])-log(1-table(data_growth$GENE)[rownames(geneTotals)]/dim(data_growth)[1])

## make object with log likelihoods input to optimisation function
mkOpt <- function(i,resultsObj=resultsList) {
    limitFn <- function(x,abslim=7) { ## limit extremes of logl to ~ 1000 exp(7) = 1096.633
        return(apply(array(x),1, function(y) {
            if(y<= -abslim) { return(-abslim);
            } else if(y>= abslim) { return(abslim);
            } else { return(y); }}));
    }    
    gnm <- rownames(geneTotals)[i];
    allLL <- cbind(limitFn(resultsObj[[7]]$quant[[gnm]] - priors_equalising[gnm]),
                   limitFn(resultsObj[[3]]$quant[[gnm]] - priors_equalising[gnm]),
                   limitFn(resultsObj[[5]]$quant[[gnm]] - priors_equalising[gnm]),
                   limitFn(resultsObj[[1]]$quant[[gnm]] - priors_equalising[gnm]));
    
    colnames(allLL) <- c('Gender','Growth','Development','HPO');
    rownames(allLL) <- rownames(data_growth);
    posIds <- rownames(data_growth)[as.character(data_growth$GENE)==gnm];
    negIds <- setdiff(rownames(data_growth),posIds);
    if((noPos<-length(posIds))+length(negIds)!=dim(data_growth)[1]) { cat('error!!'); }
    if(sum(posIds==rownames(data_growth)[resultsObj[[1]]$test[[gnm]]])!=noPos) { cat('error!!!'); }
    return(list(ll=allLL,posIds=posIds,negIds=negIds,noPos=noPos,test=resultsObj[[1]]$test[[gnm]],
                prior=priors_equalising[gnm],i=i));
}
## function to optimise [oo is set as global variable]
optLLFn <- function(wts,RETURN_VALUES=FALSE) {
    allLLwtd <- apply(sweep(oo$ll,2,wts[-1],'*'),1,sum) + oo$prior + wts[1];
    TP <- sum(allLLwtd[oo$posIds]>0);
    FP <- sum(allLLwtd[oo$negIds]>0);
    ALLP <- sum(allLLwtd>0)
    prec <- TP/ALLP;
    rec <- TP/oo$noPos;
    f1 <- 2*(prec * rec)/(prec + rec);
    if(RETURN_VALUES) { return(list(TP=TP,FP=FP,prec=prec,rec=rec,f1=f1,auc=mkAUC(oo,allLLwtd))); }
    return(1-f1);
}

## optimise test-on-training by simulated annealing in GenSA
## setup GenSA run, call mkOpt [note: optLLFn accesses oo as global variable]
tol <- 1e-13;
global.min <- 0;                        ## target for minimisation
lower <- c(-7, 0.1,  0.1, 0.1, 0.1);    ## lower and upper limits
upper <- c( 7,  2,   10,  10,  10);
params <- array(dim=c(length(geneTotals),5));
rownames(params) <- rownames(geneTotals);
colnames(params) <- c('pseudoPrior','wtGender','wtGrowth','wtDev','wtHPO');
optPerformance <- array(dim=c(length(geneTotals),6));
rownames(optPerformance) <- rownames(geneTotals);
colnames(optPerformance) <- c('TP','FP','prec','recall','F1','AUC');
for(i in 1:length(geneTotals)) {
    cat(paste('\n',rownames(optPerformance)[i]));
    set.seed(2021) # 
    oo <- NULL;
    oo <- mkOpt(i,resultsObj=resultsList); ## object accessed by optLLFn << test-on-training
    optmn <- GenSA(lower = lower, upper = upper, fn = optLLFn,
                   par = c(0,1,1,1,5), ##start-more reproducible results when changing seed
                   control=list(maxit=1500,threshold.stop=global.min+tol,verbose=TRUE));
    params[i,] <- optmn$par;
    optPerformance[i,] <- unlist(optLLFn(optmn$par,RETURN_VALUES=T));
    cat('\n'); print(optPerformance[i,]);
}

optimisedTTR <- list(params=params,performance=optPerformance);


save(optimisedTTR,file='optimisedTTR_IMPROVE_DD.RDATA');
