# Copyright (C) 2022 The University of Edinburgh 
# Author Stuart Aitken MRC HGU IGC s.aitken@ed.ac.uk
# All Rights Reserved.
# Funded by the Medical Research Council
# https://www.ed.ac.uk/mrc-human-genetics-unit

## term frequency; inverse document frequency and information content
## following the metrics defined by Groza et al. 2015 AJHG 97 111-124
## generates Fig 2 when run on DDD data; runs on demo resource data


library(gplots);

## code will run on demo data [results not that meaningful]
load(file="informativePhenotypicTerms.RData"); ## provided, see also IMPROVE_resource.r
load(file="hpoe_demo.RData"); ## generated using IMPROVE_resource.r from database of HPO annotations
load(file="pheno_demo.RData");
proband_to_diagnosis <- data.frame(pheno);
proband_to_diagnosis$GENE <- factor(proband_to_diagnosis$GENE);


## label included in pdf file names
label <- '_demo_data';

## term frequency
TF <- array(dim=c((noDiagns<-length(diagns<-sort(unique(as.character(proband_to_diagnosis$GENE))))),dim(hpoe)[2]));
colnames(TF) <- colnames(hpoe);
rownames(TF) <- diagns;
for(i in 1:noDiagns) {
    pbndsDi <- proband_to_diagnosis[proband_to_diagnosis$GENE==diagns[i],]$ID;
    TF[i,] <- apply(hpoe[pbndsDi,],2,sum);
}
dim(TF)
##[1]   77 5153
cat(paste('summary TF',dim(TF)[1],dim(TF)[2],'\n',sep=' '));
print(summary(apply(TF,2,sum))); ##

## IDF = log(Total_No_Diseases/(No diseases where term is mentioned in abstract))
IDF <- log(noDiagns/apply(TF,2,function(x){ sum(x>0); }));
TFIDF <- sweep(TF,2,IDF,'*');

cat(paste('summary IDF',length(IDF),'\n',sep=' '));
print(summary(IDF));

## IC = -log(freq)
IC <- -log((x<-apply(hpoe,2,sum))/max(x));
TFIDF_IC <- sweep(TFIDF,2,IC,'*');

TFIDF_plot <- TFIDF;
TFIDF_plot <- TFIDF_plot[,apply(TFIDF,2,function(x){sum(is.na(x));})<noDiagns];
cat(paste('summary c(TFIDF_plot)',dim(TFIDF_plot)[1],dim(TFIDF_plot)[2],'\n',sep=' '));
print(summary(c(TFIDF_plot)));

prod(dim(TFIDF_plot))-sum(c(TFIDF_plot)==0,na.rm=T);

termns <- c(names(informativeSpecificNonduplicateTerms), unlist(sapply(informativeSpecificNonduplicateTerms,names)));
termns <- termns[!is.null(termns)];
termsNonDupl <- unique(termns);

topN <- 10;
topN3 <- length(termsNonDupl);

TFIDF_IC_plot <- TFIDF_IC;
TFIDF_IC_plot <- TFIDF_IC_plot[,apply(TFIDF_IC,2,function(x){sum(is.na(x));})<noDiagns];
cat('dim plottable TFIDF_IC');
print(dim(TFIDF_IC_plot));
##[1]   77 2541
cat(paste('summary c(TFIDF_IC_plot)',dim(TFIDF_IC_plot)[1],dim(TFIDF_IC_plot)[2],'\n',sep=' '));
print(summary(c(TFIDF_IC_plot)));

## heatmap of rank of IPTs
topRank <- 800;
TFIDF_IC_rank <- array(dim=dim(TFIDF_IC_plot),0);
rownames(TFIDF_IC_rank) <- rownames(TFIDF_IC_plot);
for(i in 1:noDiagns) {
    oi <- colnames(TFIDF_IC_plot)[order(TFIDF_IC_plot[i,],decreasing=T)];
    TFIDF_IC_rank[i,oi%in%termsNonDupl] <- 1;
}

scols <- NULL;
noBlocks <- ceiling(dim(TFIDF_IC_rank)[2]/100); ## col colours blocks of 100
for(i in 1:noBlocks) { scols <- c(scols,array(dim=100,rainbow(noBlocks)[i])); }
scols <- scols[1:dim(TFIDF_IC_rank)[2]];

scolstopRank <- NULL;
for(i in 1:(topRank/100)) { scolstopRank <- c(scolstopRank,array(dim=100,rainbow((topRank/100))[i])); }


pdf(file=paste('TFIDF_IC_rank',label,'.pdf',sep=''),width=12,height=8);
heatmap.2(TFIDF_IC_rank,trace='none',cexRow=0.4,col=c('ivory','black'),ColSideColors=scols,
          cexCol=0.1,Colv=NA,Rowv=NA,labCol=array(dim=dim(TFIDF_IC_rank)[1],''));
dev.off();

odr1 <- apply(TFIDF_IC_rank,1,function(x){ median(seq(from=1,to=length(x))[x==1]);});
pdf(file=paste('TFIDF_IC_rank_sorted',label,'.pdf',sep=''),width=12,height=8);
heatmap.2(TFIDF_IC_rank[order(odr1,decreasing=F),],trace='none',cexRow=0.4,col=c('ivory','black'),ColSideColors=scols,
          cexCol=0.1,Colv=NA,Rowv=NA,labCol=array(dim=dim(TFIDF_IC_rank)[1],''));
dev.off();

pdf(file=paste('TFIDF_IC_rank_top',label,'.pdf',sep=''),width=12,height=8);
heatmap.2(TFIDF_IC_rank[,1:topRank],trace='none',cexRow=0.4,col=c('ivory','black'),ColSideColors=scolstopRank,
          cexCol=0.1,Colv=NA,Rowv=NA,labCol=array(dim=dim(TFIDF_IC_rank)[1],''));
dev.off();

pdf(file=paste('TFIDF_IC_rank_sorted_top',label,'.pdf',sep=''),width=12,height=8);
heatmap.2(TFIDF_IC_rank[order(odr1,decreasing=F),1:topRank],trace='none',cexRow=0.4,col=c('ivory','black'),ColSideColors=scolstopRank,
              cexCol=0.1,Colv=NA,Rowv=NA,labCol=array(dim=dim(TFIDF_IC_rank)[1],''));
dev.off();

pdf(file=paste('TFIDF_IC_rank_sum',label,'.pdf',sep=''),width=12,height=8);
par(mfrow=c(2,1));
plot(apply(TFIDF_IC_rank,2,sum),xlab='',main='',ylab='',type='h')
plot(apply(TFIDF_IC_rank,2,sum),xlab='',main='',ylab='',type='h',col=scols);
dev.off();

pdf(file=paste('TFIDF_IC_rank_sum_top',label,'.pdf',sep=''),width=12,height=8);
par(mfrow=c(2,1));
plot(apply(TFIDF_IC_rank[,1:topRank],2,sum),xlab='',main='',ylab='',type='h');
plot(apply(TFIDF_IC_rank[,1:topRank],2,sum),xlab='',main='',ylab='',type='h',col=scolstopRank)
dev.off();

## barplot of no. probands on ranking
pdf(file=paste('TFIDF_IC_rank_sum_top_barplot',label,'.pdf',sep=''),width=12,height=3);
barplot(rev(table(proband_to_diagnosis$GENE)[order(odr1,decreasing=F)]),col='black',las=2,cex.names=0.2)
dev.off()

topN_per_gene <- array(dim=c(dim(TFIDF_IC)[1],topN));
rownames(topN_per_gene) <- rownames(TFIDF_IC);
for(i in 1:dim(TFIDF_IC)[1]) {
    topN_per_gene[i,] <- colnames(TFIDF_IC)[order(TFIDF_IC[i,],decreasing=T)][1:topN];
}

pdf(file=paste('IC_hist_of_terms',label,'.pdf',sep=''),width=4,height=8);
par(mfrow=c(2,1))
hist(apply(topN_per_gene,1,function(x){ mean(IC[x])}),main=paste('Mean IC of top',topN,'terms\nper gene'),xlab='Information content',col=rainbow((topRank/100))[1]);
hist(IC[termsNonDupl],main='IC of informative terms',xlab='Information content',col=rainbow((topRank/100))[2]);
dev.off();

stats <- array(dim=c(noDiagns,2));
rownames(stats) <- rownames(topN_per_gene);
colnames(stats) <- c('topN','termsNonDupl');
for(i in 1:noDiagns) {
    stats[i,1] <- mean(apply(hpoe[proband_to_diagnosis[proband_to_diagnosis$GENE==rownames(topN_per_gene)[i],]$ID,topN_per_gene[i,]]*1,1,sum));
    stats[i,2] <- mean(apply(hpoe[proband_to_diagnosis[proband_to_diagnosis$GENE==rownames(topN_per_gene)[i],]$ID,termsNonDupl]*1,1,sum));
}


pdf(file=paste('hist_terms_per_case',label,'.pdf',sep=''),width=8,height=8);
par(mfrow=c(2,2));
hist(stats[,1],main=paste('No. HPO terms/case\ntop',topN,'terms per gene'),
     xlab='Mean terms/case per gene',col=rainbow((topRank/100))[1]);
hist(stats[,2],main=paste('No. HPO terms/case\ninformative terms'),
     xlab='Mean terms/case per gene',col=rainbow((topRank/100))[2]);
dev.off();

    
