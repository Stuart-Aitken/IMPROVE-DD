# Copyright (C) 2022 The University of Edinburgh 
# Author Stuart Aitken MRC HGU IGC s.aitken@ed.ac.uk
# All Rights Reserved.
# Funded by the Medical Research Council
# https://www.ed.ac.uk/mrc-human-genetics-unit

## functions used in the resource generation script and IMPROVE classifier 

## convert proband HPO annotations "HP:0000494; HP:0001999; HP:0100543"
## to HPO / ontoCAT format "HP_0000494","HP_0001999","HP_0100543"
mkUnderscore <- function(l) {
  apply(array(unlist(strsplit(l,'\\|'))),1,
        xFn <- function(x) { paste(unlist(strsplit(x,'\\:')),collapse='_'); });
}
## expand term id to all parents [include original term by default]
expandToParents <- function(id,inc=TRUE) {  
    if(inc) {
        return(union(id,applyListFn(getAllTermParentsById(OBOHPO,id),getAccession)));
    } else {
        return(union(NULL,applyListFn(getAllTermParentsById(OBOHPO,id),getAccession)))
    }
}
## ontoCAT functions return lists, define a sapply like function 
applyListFn <- function(l,fn) {  
  len <- length(l);
  if(len>0) {
    a <- array(dim=len);
    for(i in 1:len) {
      a[i] <- fn(l[[i]]);
    }
    return(a);
  } else return(NULL);
}

getChildTermCounts <- function(term,hpoep=hpoe) {
    ins <- intersect(sapply(getTermChildrenById(OBOHPO,term),getAccession),colnames(hpoep));
    if(length(ins)>1) {
        tc <- apply(hpoep[,ins],2,sum)
        tca <- array(tc);
        rownames(tca) <- names(tc);
        return(tca);
    } else if(length(ins)==1) {
        tca <- array(sum(hpoep[,ins]));
        rownames(tca) <- ins;
        return(tca);
    }
  return(NULL);
}

getDuplicatesAndParents <- function(informativeTerms) {
    duplicatesIndx <- 0;
    duplicates <- list();
    parentsIndx <- 0;
    parents <- list();
    for(i in 1:length(informativeTerms)) {
        ti <- informativeTerms[[i]];
        cat(paste('\n',names(informativeTerms)[i]));
        if(!is.null(ti)) {
            for(j in 1:length(informativeTerms)) {
                tj <- informativeTerms[[j]];
                if(!is.null(tj)) {
                    for(k in 1:length(ti)) {
                        for(l in 1:length(tj)) {
                            if(rownames(ti)[k] %in% expandToParents(rownames(tj)[l])) {
                                if(i!=j && rownames(ti)[k]==rownames(tj)[l]) { ## top level differs but duplicated child
                                    upperT <- sort(c(names(informativeTerms)[i],names(informativeTerms)[j]));
                                    duplicatesIndx <- duplicatesIndx+1;
                                    duplicates[[duplicatesIndx]] <- list(upper1=upperT[1],
                                                                         upper2=upperT[2],
                                                                         child=rownames(ti)[k]);
                                } else { ## rownames(tj)[l] under [j] has parent rownames(ti)[k] under [i]
                                    parentsIndx <- parentsIndx+1;
                                    parents[[parentsIndx]] <- list(upper1=names(informativeTerms)[i],
                                                                   parent1=rownames(ti)[k],
                                                                   upper2=names(informativeTerms)[j],
                                                                   child2=rownames(tj)[l]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    return(list(parents=parents,duplicates=duplicates));
}

## IMPOVE classifier functions
## HPO gene models
geneModelFn <- function(selgenes,hpop,d2p) { 
    prsi <- array(dim=c(length(selgenes),dim(hpop)[2]));
    rownames(prsi) <- selgenes;
    colnames(prsi) <- colnames(hpop);
    for(i in 1:length(selgenes)) {  ## note all of diagns_to_probands
        prsi[i,] <- apply(hpop[unlist(d2p[selgenes[i]]),],2,sum);
        prsi[i,] <- (prsi[i,] + 1) /(sum(prsi[i,]) + 100);  ### normalise by no. annotations + pseudo count
    }
    return(prsi);
}

## improve_HPO: generate HPO gene models for gene and "other" in data$GENE
## from proband*HPO matrix hpou which must be compiled before running improve
## [see IMPROVE_resouce.r]
## return log likelihoods for each proband and model
## either testing on training data or by generating model omitting each individual
## in turn and testing that individual (LOO=F or T)
## informativePhenotypicTerms accessed
improve_HPO <- function(data,priors,LOO,dtest=NULL) {
    if(sum(priors>0)>0) { cat('\nerror priors should be logged'); return(NULL); }
    gns <- c(setdiff(unique(as.character(data$GENE)),"other"),"other");
    rownames(priors) <- gns;
    trmns <- c(names(informativePhenotypicTerms), unlist(sapply(informativePhenotypicTerms,names)));
    trmns <- trmns[!is.null(trmns)];
    uqtopinform <- unique(trmns);
    probandIDs <- data$ID;
    ## subset of hpoe setting usage of toplevel terms (names(informativePhenotypicTerms)) to have
    ## usage 0 where one or more specific child term is used (toplevel term is ** X-other ** cf X)
    ## retaining completeness of usage
    hpou_subset <- hpou[probandIDs,uqtopinform]*1;
    for(i in 1:length(informativePhenotypicTerms)) {
        if(!is.null(informativePhenotypicTerms[[i]])) {
            hpou_subset[apply(hpou_subset[,names(informativePhenotypicTerms[[i]])],1,sum)>0, ##probands with any term in list
                        names(informativePhenotypicTerms)[i]] <- 0; ## set top level term usage to 0
        }
    }
    pbndsi <- intersect(unlist(diagns_to_probands[gns[1]]),probandIDs);
    pbndso <- setdiff(probandIDs,pbndsi);
    diagns_to_probands_subset <- list(list(pbndsi),list(pbndso));
    names(diagns_to_probands_subset) <- c(gns[1],'other');
    prary <- array(dim=c(dim(hpou_subset)[1],2));
    colnames(prary) <- gns;
    rownames(prary) <- rownames(hpou_subset);
    if(!LOO) {
        prsinz <- geneModelFn(gns, hpou_subset, diagns_to_probands_subset); ##
        no_terms <- array(dim=c(dim(hpou_subset)[1],dim(prsinz)[1]));
        rownames(no_terms) <- rownames(hpou_subset);
        for(i in 1:dim(hpou_subset)[1]) { ## all probands
            for(j in 1:dim(prsinz)[1]) {  ## all genes
                pij <- hpou_subset[i,] * prsinz[j,];   ## probabilities of terms that actually occur 
                pij <- pij[pij!=0];                    ## (ignoring absent terms with probs)
                no_terms[i,j] <- length(pij);
                ## product of probs [scaled to 10 terms] and prior (Mitchell)
                prary[i,j] <- (sum(log(pij)) * (10/no_terms[i,j])) + priors[j];
            }
        }
        if(!is.null(dtest)) {
            cat('\n ** returning results on test data only **\n');
            hpou_subset_test <- hpou[rownames(dtest),uqtopinform]*1; 
            for(i in 1:length(informativePhenotypicTerms)) {
                if(!is.null(informativePhenotypicTerms[[i]])) {
                    hpou_subset_test[apply(hpou_subset_test[,names(informativePhenotypicTerms[[i]])],1,sum)>0, ##probands with any term in list
                                     names(informativePhenotypicTerms)[i]] <- 0; ## set top level term usage to 0
                }
            }
            prarytest <- array(dim=c(dim(hpou_subset_test)[1],2));
            colnames(prarytest) <- gns;
            rownames(prarytest) <- rownames(hpou_subset_test);               
            no_terms <- array(dim=c(dim(hpou_subset_test)[1],dim(prsinz)[1]));
            rownames(no_terms) <- rownames(hpou_subset_test);
            for(i in 1:dim(hpou_subset_test)[1]) { ## all probands
                for(j in 1:dim(prsinz)[1]) {  ## all genes
                    pij <- hpou_subset_test[i,] * prsinz[j,];   ## probabilities of terms that actually occur 
                    pij <- pij[pij!=0];                          ## (ignoring absent terms with probs)
                    no_terms[i,j] <- length(pij);
                    prarytest[i,j] <- (sum(log(pij)) * (10/no_terms[i,j])) + priors[j];
                    ## product of probs and prior (Mitchell)
                }
            }
            return(prarytest);
        }
    } else { ## leave one out cross-validation
        for(i in 1:dim(hpou_subset)[1]) { ## all probands
            ## removing probandIDs[i] from diag 2 proband is sufficient to derive bootstrap probs
            diagns_to_probandsi <- vector('list',length=length(diagns_to_probands_subset));
            names(diagns_to_probandsi) <- names(diagns_to_probands_subset);
            for(j in 1:length(diagns_to_probandsi)) {
                diagns_to_probandsi[[j]] <- list(setdiff(unlist(diagns_to_probands_subset[[j]]),probandIDs[i]));
            }
            
            prsinz <- geneModelFn(gns, hpou_subset, diagns_to_probandsi);
            no_terms <- array(dim=dim(prsinz)[1]);
            for(j in 1:dim(prsinz)[1]) {  ## all genes
                pij <- hpou_subset[i,] * prsinz[j,];   ## probabilities of terms that actually occur 
                pij <- pij[pij!=0];                    ## (ignoring absent terms with probs)
                no_terms[j] <- length(pij);
                ## product of probs [scaled to 10 terms] and prior (Mitchell)
                prary[i,j] <- (sum(log(pij)) * (10/no_terms[j])) + priors[j]; 
            }
        }
    }
    return(prary);
}

## improve_nBayes_nominal: generate naive Bayes gene models for nominal data
## calls improve_nBayes with alternative parameters
improve_nBayes_nominal <- function(data,priors,LOO,dtest=NULL) {
    improve_nBayes(data,priors,LOO,dtest=dtest,
                   treatAsContinuous=FALSE,PLOT_DENSITY=F);
}
## improve_nBayes: generate naive Bayes gene models for gene and "other" in data$GENE
## from continuous data (treatAsContinuous=TRUE) or nominal (treatAsContinuous=FALSE)
## return log likelihoods for each proband and model
## either testing on training data or by generating model omitting each individual
## in turn and testing that individual (LOO=F or T)
improve_nBayes <- function(data,priors,LOO,dtest=NULL,
                           treatAsContinuous=TRUE,laplace=0.5,PLOT_DENSITY=F) {
    if(!'GENE'%in%colnames(data)) { cat('error - GENE is not a column'); return(NULL); }
    allcols <- colnames(data); 
    colsTest <- setdiff(allcols,"GENE");
    if(treatAsContinuous) {
        if(!LOO) {            
            nbayes <- naive_bayes(GENE ~ ., data = data[,allcols],
                                  usekernel = T,bw='nrd0',adjust=1.5, kernel = "gaussian");
            if(PLOT_DENSITY) {
                plot(nbayes, prob = "conditional");
                legend('topleft',legend=paste(unique(get_cond_dist(nbayes)),collapse=';'));
            }
            if(!is.na(priors)[1]) {
                ## check priors in input match those calculated by naive_bayes
                if(sum(log(nbayes$prior)==priors)!=length(priors)) { cat('\npriors differ'); return(NULL); } else { cat('.'); }
            } else { cat('prior is NA'); }            
            if(!is.null(dtest)) {
                cat('\n ** returning results on test data only **\n');
                prary     <- predict(nbayes, dtest[,colsTest], type = "prob");
            } else {
                prary     <- predict(nbayes, data[,colsTest], type = "prob");
            }
        } else {
            gns <- rownames(priors);
            prary <- array(dim=c(dim(data)[1],length(gns)));
            rownames(prary) <- rownames(data);
            colnames(prary) <- gns;
            for(i in 1:dim(data)[1]) { ## prary is pbnd * gene array (all genes)
                nbayes <- naive_bayes(GENE ~ ., data = data[-i,allcols],
                                      usekernel = T,bw='nrd0',adjust=1.5, kernel = "gaussian");
                if(PLOT_DENSITY && data[i,]$GENE!="other") {
                    plot(nbayes, prob = "conditional");
                    legend('topleft',legend=paste(unique(get_cond_dist(nbayes)),collapse=';'));
                }
                prary[i,] <- predict(nbayes, data[i,colsTest], type = "prob");
            }
        }
    } else {
        data1 <- data.frame(data[,allcols],stringsAsFactors=FALSE);
        colnames(data1) <- allcols;
        if(!LOO) {
            if(!is.null(dtest)) {
                cat('\n ** returning results on test data only **\n');
                data1test <- data.frame(dtest[,colsTest],stringsAsFactors=TRUE);
                colnames(data1test) <- colsTest;
            } else {
                data1test <- data.frame(data[,colsTest],stringsAsFactors=TRUE);
                colnames(data1test) <- colsTest;
            }
            nbayes <- naive_bayes(GENE ~ ., data = data1,laplace=laplace); ##
            prary  <- predict(nbayes, data1test, type = "prob");
            if(!is.na(priors)[1]) { if(sum(log(nbayes$prior)==priors)!=length(priors)) { cat('\npriors differ'); print(priors); } else { cat('.'); } } else { cat('prior is NA'); }
        } else {
            gns <- rownames(priors);
            prary <- array(dim=c(dim(data)[1],length(gns)));
            rownames(prary) <- rownames(data);
            colnames(prary) <- gns;
            data1test <- data.frame(data[,colsTest],stringsAsFactors=TRUE);
            colnames(data1test) <- colsTest;
            for(i in 1:dim(data)[1]) { ## prary is pbnd * gene array (all genes)
                nbayes    <- naive_bayes(GENE ~ ., data = data1[-i,],laplace=laplace); ##
                prary[i,] <- predict(nbayes, data1test, type = "prob")[i,];
            }
        }
    }
    return(log(prary[,rownames(priors)]));
}


## run improve() (pryFn) on all genes testing models for gene vs not-gene
## test on training or leave-one-out (LOO=F or T)
## performs AUC analysis from classifier margin
## returns only log likelihoods of test data if dtest is not NULL after training on data
## hpou is used by improve() [not passed as an argument]
run_improve <- function(data,gnsCounts,pryFn,priors=NA,
                        filename='test',PLOT=TRUE,LOO=FALSE,dtest=NULL) {
    gns <- rownames(gnsCounts);
    resultsQuant <- vector('list',length=length(gns));
    names(resultsQuant) <- gns;
    resultsTest <- vector('list',length=length(gns));
    names(resultsTest) <- gns;
    resultsPrary <- vector('list',length=length(gns));
    names(resultsPrary) <- gns;
   
    priors_used <- array(dim=c(length(gns),2));
    colnames(priors_used) <- c('gene','other');
    rownames(priors_used) <- gns;
    for(i in 1:length(gns)) { ## for all genes in gnsCounts
        di <- data;  ## re-label as gene gns[i] and "other"
        levels(di$GENE) <- c(levels(di$GENE),"other");
        di[as.character(di$GENE)!=gns[i],]$GENE <- "other";
        di$GENE <- droplevels(di$GENE); ## prior as fractions
        if(is.na(priors)[1]) {
            priorsi <- table(di$GENE);
            priorsia <- array(log(priorsi/sum(priorsi)));
            rownames(priorsia) <- rownames(priorsi);
            priorsia <- priorsia[c(gns[i],'other')];
        } else {  ## or prior as provided
            priorsia <- priors;
        }
        priors_used[gns[i],1] <- priorsia[1];
        priors_used[gns[i],2] <- priorsia[2];
        prary <- pryFn(di,priorsia,LOO,dtest); ## call pryFn()
        if(is.null(dtest)) {
            prdiff <- array(prary[,gns[i]]-prary[,'other']); ## margin
            rownames(prdiff) <- di$ID;
            decn <- array(dim=length(prdiff),FALSE); ## decision
            rownames(decn) <- di$ID;
            decn[prdiff>0] <- TRUE;
            hasgene <- array(dim=length(prdiff),FALSE); ## positive cases
            rownames(hasgene) <- di$ID;
            hasgene[as.character(di$GENE)==gns[i]] <- TRUE;
            resultsQuant[[i]] <- prdiff;
            resultsTest[[i]]  <- hasgene;
            resultsPrary[[i]] <- prary;
        } else {
            resultsPrary[[i]] <- prary;
        }
    }
    if(is.null(dtest)) {
        pred <- prediction(sapply(resultsQuant, function(x) { list(exp(x-max(x)));}),
                           sapply(resultsTest,function(x) { list(x*2-1);}));
    
        if(PLOT) {
            perf <- performance(pred, "tpr", "fpr");
            pdf(file=paste(filename,'.pdf',sep=''),width=8,height=8);
            plot(perf,lty=1,col="grey40",colorize=F,lwd= 0.5, main= "ROC for all genes");
            perf <- performance(pred, "auc");
            AUC <- unlist(perf@y.values);
            
            par(mfrow=c(2,2));
            for(i in 1:length(gns)) {
                pred <- prediction(list(exp(resultsQuant[[i]]-max(resultsQuant[[i]]))),
                                   list(resultsTest[[i]]*2-1));
                perf <- performance(pred, "tpr", "fpr");
                plot(perf,colorize=TRUE,main=paste(gns[i],'AUC',signif(AUC[i],2)));
            }
            dev.off();
        } else {
            perf <- performance(pred, "auc");
            AUC <- unlist(perf@y.values);
        }
        AUC <- array(AUC);
        rownames(AUC) <- gns;
    } else {
        AUC <- NA;
    }
    return(list(praryList=resultsPrary,AUC=AUC,
                priors=priors_used,
                quant=resultsQuant,test=resultsTest));
}

mkAUC <- function(oo,x) {
    quant <- exp(x-max(x))
    pred <- prediction(list(quant),list(oo$test*2-1));
    perf <- performance(pred, "auc");
    return(unlist(perf@y.values));
}

statsFromResultsObjFn <- function(obj,i,prior=0,geneTotalsp=geneTotals,dataf=proband_to_diagnosis) {
    gnm <- rownames(geneTotalsp)[i];
    allLL <- obj$quant[[gnm]] + prior;
    rownames(allLL) <- array(dataf$ID);
    posIds <- dataf$ID[as.character(dataf$GENE)==gnm];
    negIds <- setdiff(dataf$ID,posIds);
    if((noPos<-length(posIds))+length(negIds)!=dim(dataf)[1]) { cat('error!!'); }
    if(sum(posIds==dataf$ID[obj$test[[gnm]]])!=noPos) { cat('error!!!'); }
    TP <- sum(allLL[posIds]>0);
    FP <- sum(allLL[negIds]>0);
    ALLP <- sum(allLL>0)
    prec <- TP/ALLP;
    rec <- TP/noPos;
    f1 <- 2*(prec * rec)/(prec + rec);
    ooo <- NULL;
    ooo$test <- obj$test[[i]];
    return(list(TP=TP,FP=FP,prec=prec,rec=rec,f1=f1,auc=mkAUC(ooo,allLL))); 
}

mkQuantPerformance <-function(statsGvsO,prior=0,geneTotalsp=geneTotals,dataf=proband_to_diagnosis) {
    performancep <- array(dim=c(length(geneTotalsp),6));
    rownames(performancep) <- rownames(geneTotalsp);
    colnames(performancep) <- c('TP','FP','prec','recall','F1','AUC');
    for(i in 1:length(geneTotalsp)) {
        performancep[i,] <- unlist(statsFromResultsObjFn(statsGvsO,i,prior=prior,geneTotalsp=geneTotalsp,dataf=dataf));
    }
    if(sum(performancep[,6]==statsGvsO$AUC)!=length(geneTotalsp)) { cat('error'); }
    return(performancep);
}
