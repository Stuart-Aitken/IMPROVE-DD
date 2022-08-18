# Copyright (C) 2022 The University of Edinburgh 
# Author Stuart Aitken MRC HGU IGC s.aitken@ed.ac.uk
# All Rights Reserved.
# Funded by the Medical Research Council
# https://www.ed.ac.uk/mrc-human-genetics-unit

## Code to create informative phenotypic terms (IPTs) and the proband*HPO array hpou
## from a database of HPO terms. These are used by the IMPROVE HPO classifier.
## Generates resources for TF IDF IC (Fig 2) and classification (Fig 3) when run on DDD data.
## Generates demo resource data to allow TF IDF IC and classifier code to be run.
## These steps are run once for the database of HPO terms at hand.
## IPTs generated from the DDD study data are provided in: informativePhenotypicTerms.RData

## 

library(hash);
library(rJava);
library(ontoCAT);


source(file="IMPROVE_functions.r");

## Depends on Human phenotype ontology https://hpo.jax.org/
## Version 02/08/2021 available here:
## https://bioportal.bioontology.org/ontologies/HP
## should be saved as: hp_2021_08_02.obo
OBOHPO <- getOntology("hp_2021_08_02.obo"); 


## For demonstration, define HPO annotations (pheno$child_hpo) for 12 probands
## For a new database, create an equivalent data structure of HPO annotations, individual ids and diagnoses. Note that diagnoses are not used here but in the classification step.
pheno <- NULL;
pheno$child_hpo <- array(c("HP:0000494; HP:0000176; HP:0100851; HP:0000750; HP:0000286; HP:0001263; HP:0001999; HP:0100543",
                     "HP:0002558; HP:0000158; HP:0010535; HP:0000403", 
                     "HP:0009843; HP:0008734; HP:0007018",
                     "HP:0000684; HP:0005709; HP:0000954; HP:0000430; HP:0002007; HP:0001374; HP:0009467; HP:0000347",
                     "HP:0002019; HP:0001328; HP:0002342; HP:0000750; HP:0001548; HP:0011407; HP:0004233",
                     "HP:0010529; HP:0000494; HP:0011343",
                     "HP_0000707; HP_0001939; HP_0000818; HP_0000769; HP_0001507; HP_0001608; HP_0002715; HP_0001626; HP_0000119; HP_0002664; HP_0000598; HP_0025354; HP_0000478; HP_0001197; HP_0000152; HP_0033127; HP_0001871; HP_0001574; HP_0040064; HP_0002086; HP_0025142; HP_0025031; HP_0012447",
                     
                     "HP_0011146; HP_0020219; HP_0033259; HP_0002197; HP_0100543; HP_0001328; HP_0000736; HP_0002360; HP_0100851; HP_0000752; HP_0004305; HP_0011443; HP_0002493; HP_0001288; HP_0011282; HP_0002118; HP_0001256; HP_0002342; HP_0010864; HP_0001270; HP_0000717; HP_0007370; HP_0001344; HP_0011342; HP_0011343; HP_0011344; HP_0010993; HP_0004323; HP_0000098; HP_0004322; HP_0001511; HP_0002597",
                     
                     "HP_0011025; HP_0025015; HP_0001627; HP_0000079; HP_0000078; HP_0031704; HP_0000370; HP_0000377; HP_0000357; HP_0000539; HP_0000508; HP_0000504; HP_0000496; HP_0000553; HP_0004328; HP_0008056; HP_0004329; HP_0000490; HP_0000316; HP_0000464; HP_0001965",
                     
                     "HP_0009116; HP_0011821; HP_0000309; HP_0001999; HP_0000306; HP_0000290; HP_0000235; HP_0002648; HP_0040194; HP_0010938; HP_0005288; HP_0031815; HP_0000422; HP_0000463; HP_0003196; HP_0005484; HP_0000534; HP_0000202; HP_0000174; HP_0000164; HP_0011337; HP_0000286; HP_0012471",
                     
                     "HP_0000178; HP_0000233; HP_0000177; HP_0200007; HP_0200006; HP_0003549; HP_0011805; HP_0009115; HP_0001367; HP_0100261; HP_0009122; HP_0000925; HP_0000765; HP_0011314; HP_0001388; HP_0001371; HP_0001382; HP_0001252; HP_0001276; HP_0030084; HP_0011122; HP_0001597; HP_0011356; HP_0011354; HP_0011355; HP_0001000; HP_0000499; HP_0011362; HP_0010720; HP_0100037; HP_0009815; HP_0009810",
                     
                     "HP_0040070; HP_0006496; HP_0040069; HP_0100491; HP_0045060; HP_0005927; HP_0100871; HP_0009484; HP_0005922; HP_0005656; HP_0001763; HP_0001780; HP_0006494; HP_0001159; HP_0011927; HP_0100807; HP_0005918; HP_0004207; HP_0001211; HP_0001172; HP_0004097; HP_0011024; HP_0004298; HP_0025033; HP_0012719; HP_0011458"));
pheno$ID <-  paste('id',seq(from=1,to=length(pheno$child_hpo)),sep=''); ## make ids
pheno$GENE <- c('X','X','X','X','X','X','Y','Y','Y','Y','Y','Y');       ## diagnoses used in classification
rownames(pheno$child_hpo) <- pheno$ID;
save(pheno,file="pheno_demo_data.RData");

## 1. create the arrays of HPO data (individuals * HPO terms)
allHPO <- NULL;  # all HPO terms used in data
for(i in 1:length(pheno$child_hpo)) {
  allHPO <- union(allHPO,
                  mkUnderscore(unlist(strsplit(pheno$child_hpo[i],'; '))));
}
length(allHPO)==181
allHPO <- allHPO[!is.na(allHPO) & allHPO!="NA"]; # eliminate "NA"
length(allHPO)==181
length(unique(allHPO))==181

## hpo : proband * HPO array for direct proband-term annotations
hpo <- array(dim=c(length(allHPO),length(pheno$child_hpo)),F); # array of HPO terms used * probands
rownames(hpo) <- allHPO;  # rownames HPO terms
colnames(hpo) <- pheno$ID;
noOfSwaps <- 0;
for(i in 1:length(pheno$child_hpo)) {  # for each proband (column)
  hpoi <- mkUnderscore(unlist(strsplit(pheno$child_hpo[i],'; '))); # get HPO terms
  hpoi <- hpoi[!is.na(hpoi) & hpoi!="NA"]; # eliminate "NA"
  if(length(hpoi)>0) {
      hpo[hpoi,i] <- T;  # set rows (HPO terms) in column (proband) i to true for [mapped] terms
  }
}

hpo <- t(hpo);
mean((apply(hpo,2,sum))); ## 1.049724 uses of HPO terms in demo data [22.6 in DDD data]
mean((apply(hpo,1,sum))); ## 15.83333  # HPO terms / case in demo data [7.0 in DDD data]
dim(hpo) == c(12,181);

## hpoe : proband * HPO array for propagated proband-term annotations
## annotations are propagated to parent terms
## {T/F} indicates term (column) holds of proband (row)

allIDS <- getAllTermIds(OBOHPO)
allHPOExpanded <- NULL;
for(i in 1:length(allHPO)) {
  if(allHPO[i] %in% allIDS) {
    allHPOExpanded <- union(allHPOExpanded,expandToParents(allHPO[i]));
  } else { cat(paste('not in HPO',allHPO[i],';\n')); }
}

sum(is.na(allHPOExpanded))==0
sum(is.na(allHPOExpanded=="NA"))==0

## array of HPO terms used EXPANDED * probands
unk <- 0;
hpoe <- array(dim=c(length(allHPOExpanded),length(pheno$child_hpo)),F);
rownames(hpoe) <- allHPOExpanded;
colnames(hpoe) <- pheno$ID;
for(i in 1:length(pheno$child_hpo)) {
    hpoei <- mkUnderscore(unlist(strsplit(pheno$child_hpo[i],'; ')));
    if(is.na(hpoei) || hpoei=="NA") { cat('error'); return(NULL); }
    if(length(hpoei)>0) {
        hpoeiex <- NULL;
        for(j in 1:length(hpoei)) {
            if(hpoei[j] %in% allIDS) {
                hpoeiex <- union(hpoeiex,expandToParents(hpoei[j]));
            } else { unk <- unk+1; cat(paste('not in HPO',hpoei[j],';\n')); }
        }
    }
    hpoe[hpoeiex,i] <- T;
}
i==12;
unk==0;
hpoe <- t(hpoe);
dim(hpoe)
##[1]   12 315
## hpoe is further modified after creating informativePhenotypicTerms, see hpou below
## hpoe is used in TF-IDF-IC analysis
save(hpoe,file="hpoe_demo.RData");

## 2. Create the informative phenotypic terms list:
## a list of entries for top-level terms (names of the list informativeSpecificTerms and of sub lists)
## where top-level term is expanded, entry is a list of counts, names are informative child terms
## This step requires the full hpoe matrix from DDD data to generate informativePhenotypicTerms.RData
## Note: This should be skipped when using the demo data as the result is provided:
## informativePhenotypicTerms.RData

## toplevel terms below root
toplevel <- sapply(getTermChildrenById(OBOHPO,'HP_0000118'),getAccession);
toplevel <- intersect(toplevel,colnames(hpoe));

## no. annotations of child terms
toplevelcounts <- getChildTermCounts('HP_0000118');
sum(rownames(toplevelcounts)==toplevel)==length(toplevel);

## expand top-level terms with >= upper annotation threshold; retain terms with >= retain threshold
upper_annotation_threshold <- 1500; ## approx 10%; 
retain_threshold <- 250;            ## approx 2%;
informativeSpecificTerms <- vector('list',length=length(toplevelcounts));
names(informativeSpecificTerms) <- rownames(toplevelcounts);
for(i in 1:length(toplevelcounts)) {
    if(toplevelcounts[i]>=upper_annotation_threshold) {  ## expand
        tci <- getChildTermCounts(rownames(toplevelcounts)[i]);
        informativeSpecificTerms[[i]] <- array(tci[tci>=retain_threshold]);
        rownames(informativeSpecificTerms[[i]]) <- names(tci[tci>=retain_threshold]);
    }    
}
## while terms with frequency greater than upper_annotation_threshold remain expand those terms
while(max(sapply(informativeSpecificTerms,function(x) { if(!is.null(x)) { return(max(unlist(x))); }; 0; })) > upper_annotation_threshold) { 
    for(i in 1:length(informativeSpecificTerms)) {
        if(!is.null(informativeSpecificTerms[[i]])) {
            tci <- informativeSpecificTerms[[i]];
            tcidelete <- NULL;
            tciadd <- NULL;
            for(j in 1:length(tci)) {
                if(tci[j]>=upper_annotation_threshold) {  ## expand
                    tcj <- getChildTermCounts(rownames(tci)[j]);
                    tcjsel <- tcj[tcj>=retain_threshold];
                    if(length(tcjsel)>0) {
                        tcidelete <- c(tcidelete,rownames(tci)[j]);
                        tciadd <- c(tciadd,tcjsel);
                    }
                    tcikeep <- !(rownames(tci)%in%tcidelete);
                    informativeSpecificTerms[[i]] <- array(c(tci[tcikeep],tciadd));
                    rownames(informativeSpecificTerms[[i]]) <- names(c(tci[tcikeep],tciadd));
                    informativeSpecificTerms[[i]] <- unique(informativeSpecificTerms[[i]]);
                }
            }
        }    
    }
}
## terms are reduced to a set in use
uqtopinform <-  unique(c(names(informativeSpecificTerms), array(unlist(sapply(informativeSpecificTerms,names)))));
length(uqtopinform)==157; ## from DDD HPO data

## run checkIPTList - find terms with multiple paths from top-level parent
iST <- checkIPTList(informativeSpecificTerms);
## flatten iST$multi to array columns: upper-term no. 1 , upper-term no. 2 , common child term 
iST_multis <- t(sapply(iST$multi,function(x){ c((x[[1]]),(x[[2]]),x[[3]]) }));
colnames(iST_multis) <- c('upper1','upper2','child');
dim(iST_multis)[1]==length(iST$multi);
dim(iST_multis) == c(42,  3)
iST_multis <- unique(iST_multis); 
dim(iST_multis) == c(21,  3);
rev(sort(table(iST_multis[,3])))

## We list the follow steps used to generate
## the informativePhenotypicTerms data structure derived from DDD data.
## This step can be automated, see B. below line 256.

## 2A. A single top-level parent was selected by inspection for terms with multiple paths to top-level.
iST_multis <- iST_multis[order(iST_multis[,2]),]
txt <- apply(iST_multis,1,function(x) {
    paste('\n',getLabel(getTermById(OBOHPO,x[[1]])),'\t',
          getLabel(getTermById(OBOHPO,x[[2]])),'\t',
          getLabel(getTermById(OBOHPO,x[[3]])),
          '\n',x[[1]],'\t',
          x[[2]],'\t',
          x[[3]],
          sep='');
})
cat(txt);
## top-level, child-to-be-deleted [manually curated from txt0 txt and HPO, multiple parents identified below]
delete_multi <- matrix(c('HP_0033127','HP_0005484', ##HP_0005484 has 3 top-level parents, delete 2
                              'HP_0000707','HP_0005484',                             
                              'HP_0001574','HP_0000534', 
                              'HP_0033127','HP_0009116', 
                              'HP_0033127','HP_0011821',
                              'HP_0033127','HP_0000235', 
                              'HP_0033127','HP_0002648', 
                              'HP_0033127','HP_0040194',                               
                              'HP_0033127','HP_0040069',                             
                              'HP_0033127','HP_0040070',                              
                              'HP_0033127','HP_0045060',                              
                              'HP_0033127','HP_0001159',                              
                              'HP_0033127','HP_0011927',                              
                              'HP_0033127','HP_0001780',                            
                              'HP_0033127','HP_0100807',                                
                              'HP_0033127','HP_0005918',                               
                              'HP_0033127','HP_0004207',                               
                              'HP_0033127','HP_0001211',
                              'HP_0033127','HP_0001172',
                              'HP_0033127','HP_0004097'),nrow=20,byrow=T);

## apply deletion in delete_multi to informativeSpecificTerms
## gives informativeSpecificNonMPTerms
informativeSpecificNonMPTerms <- vector('list',length=length(informativeSpecificTerms));
names(informativeSpecificNonMPTerms) <- names(informativeSpecificTerms);
for(i in 1:length(informativeSpecificTerms)) {
    if(!is.null(informativeSpecificTerms[[i]])) {
        tmsi <- informativeSpecificTerms[[i]];
        delni <- delete_multi[delete_multi[,1]==names(informativeSpecificTerms)[i],];
        if(length(delni)>0) {
            if(length(delni)==2) {
                informativeSpecificNonMPTerms[[i]] <- tmsi[setdiff(rownames(tmsi),delni[2])];
            } else {
                informativeSpecificNonMPTerms[[i]] <- tmsi[setdiff(rownames(tmsi),delni[,2])];
            }
        } else {
            informativeSpecificNonMPTerms[[i]] <- tmsi;
        }
    }
}

## save as 'informativePhenotypicTerms'
informativePhenotypicTerms <- informativeSpecificNonMPTerms;
save(informativePhenotypicTerms,file="informativePhenotypicTerms.RData")
## end of A.

## 2B. automatic processing of multi-path terms
iST_multis <- iST_multis[order(iST_multis[,3]),]
uqtms <- unique(iST_multis[,3]);

delete_multi <- NULL;
for(i in 1:length(uqtms)) {
    uqpi <- unique(c(iST_multis[iST_multis[,3]==uqtms[i],1:2]))[-1]
    if((li<-length(uqpi))>0) {
        nrws <- array(dim=c(li,2));
        nrws[,2] <- uqtms[i];
        nrws[,1] <- uqpi;
        delete_multi <- rbind(delete_multi,nrws);
    } 
}

## apply deletion in delete_multi to informativeSpecificTerms
## gives informativeSpecificNonMPTerms
informativeSpecificNonMPTerms <- vector('list',length=length(informativeSpecificTerms));
names(informativeSpecificNonMPTerms) <- names(informativeSpecificTerms);
for(i in 1:length(informativeSpecificTerms)) {
    if(!is.null(informativeSpecificTerms[[i]])) {
        tmsi <- informativeSpecificTerms[[i]];
        delni <- delete_multi[delete_multi[,1]==names(informativeSpecificTerms)[i],];
        if(length(delni)>0) {
            if(length(delni)==2) {
                informativeSpecificNonMPTerms[[i]] <- tmsi[setdiff(rownames(tmsi),delni[2])];
            } else {
                informativeSpecificNonMPTerms[[i]] <- tmsi[setdiff(rownames(tmsi),delni[,2])];
            }
        } else {
            informativeSpecificNonMPTerms[[i]] <- tmsi;
        }
    }
}


## save as 'informativePhenotypicTerms' - gives the same set of IPTs as A.
informativePhenotypicTerms <- informativeSpecificNonMPTerms;
save(informativePhenotypicTerms,file="informativePhenotypicTerms.RData")


## 3. Create hpou: a modified hpoe where annotations to IPTs are edited out
## for individuals with annotations to an IPT that is a child of that IPT
## use the informativePhenotypicTerms provided, or that generated from a
## different set of HPO annotations after completing 2A. or 2B. above.
load("informativePhenotypicTerms.RData")

## First, find any parent-child links within informativePhenotypicTerms
iSNDT <- checkIPTList(informativePhenotypicTerms);
## tidy data structure
iSNDT_parens <- t(sapply(iSNDT$parents,function(x){ c((x[[1]]),(x[[2]]),x[[3]],x[[4]]) }));
colnames(iSNDT_parens) <- c('upper1','parent1','upper2','child2');
dim(iSNDT_parens)[1]==length(iSNDT$parents);
dim(iSNDT_parens) == c(164,4); ## 164 for informativePhenotypicTerms
iSNDT_parens <- unique(iSNDT_parens);
dim(iSNDT_parens) == c(164,4); ##
sum(del<-apply(iSNDT_parens,1,function(x) { x[1]==x[3] && x[2]==x[4];}))==135;
iSNDT_parens <- iSNDT_parens[!del,];
dim(iSNDT_parens) == c(29,4); ## 

iSNDT_child_parent_correlated <- iSNDT_parens[,c(4,2)]
colnames(iSNDT_child_parent_correlated) <- c('child','parent');

## make the edits
hpou <- hpoe;
for(i in 1:dim(iSNDT_child_parent_correlated)[1]) {
    cat(paste('\n',i,sum(hpou[,iSNDT_child_parent_correlated[i,2]]),sum(hpou[,iSNDT_child_parent_correlated[i,1]]))); ## in rows where child is TRUE set parent FALSE
    hpou[hpou[,iSNDT_child_parent_correlated[i,1]],iSNDT_child_parent_correlated[i,2]] <- FALSE;
    cat(paste('\t',i,sum(hpou[,iSNDT_child_parent_correlated[i,2]]),sum(hpou[,iSNDT_child_parent_correlated[i,1]])));
}

## hpou is used in the IMPROVE classifier
save(hpou,file="hpou_demo.RData");

## done
