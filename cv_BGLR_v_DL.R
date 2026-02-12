####################################################################################################
################--------Filter phenotypic data and calculate adjusted means(BLUPs)-----#############
####################################################################################################
library(lmerTest)
library(tidyr)
library(stringr)
################################import data######################################
pheno_2014=read.csv("g2f_2014_hybrid_data_clean.csv",h=T,na.strings='NA',sep=",",check.names=F)
pheno_2015=read.csv("g2f_2015_hybrid_data_clean.csv",h=T,na.strings='NA',sep=",",check.names=F)
pheno_2016=read.csv("g2f_2016_hybrid_data_clean.csv",h=T,na.strings='NA',sep=",",check.names=F)
pheno_2017=read.csv("g2f_2017_hybrid_data_clean.csv",h=T,na.strings='NA',sep=",",check.names=F)
################################recombine########################################
column_id=c(1,2,4,5,6,7,20,22:34)
pheno_master=rbind(pheno_2014[,column_id],pheno_2015[,column_id],pheno_2016[,column_id],pheno_2017[,column_id])
################################check names######################################
pedigree_list=grep("check",pheno_master$Source,ignore.case = T)
check_variety=as.character(unique(pheno_master[pedigree_list,4]))
conv=paste("check",seq(1,length(check_variety)))
pheno_master$Pedigree=as.character(pheno_master$Pedigree)
pheno_master[pedigree_list,4]=conv[match(pheno_master[pedigree_list,4],check_variety)]
################################data strcuture####################################
pheno_master$`﻿Year`=as.factor(pheno_master$`﻿Year`)
pheno_master$`Field-Location`=as.factor(pheno_master$`Field-Location`)
pheno_master$Block=as.factor(pheno_master$Block)
pheno_master$Pedigree=as.factor(pheno_master$Pedigree)
pheno_master$DTA=as.numeric(as.Date(as.character(pheno_master$`Anthesis [date]`), format="%m/%d/%Y")-as.Date(as.character(pheno_master$`Date Planted`), format="%m/%d/%Y"))
pheno_master$env=as.factor(paste(pheno_master$`﻿Year`,pheno_master$`Field-Location`, sep="_"))
pheno_master$Pedigree=toupper(pheno_master$Pedigree)
################################table missing data from each env####################
# na_rate=table(is.na(pheno_master$`Grain Yield [bu/A]`),pheno_master$env)
# freq=as.data.frame(prop.table(na_rate, 2))
# freq=freq[freq$Var1=="TRUE",]
# ################################match lines among env#############################
# pheno_master_list=as.data.frame(str_split_fixed(pheno_master$Pedigree, "/",2))
# pheno_masteringeno=pheno_master[(pheno_master_list$V1%in%line_filtered)&(pheno_master_list$V1%in%line_filtered),]
# table(pheno_masteringeno$env)
# testdata <- split(pheno_masteringeno[,c(4,22)],pheno_masteringeno[,c(4,22)][2])
# intersect(unlist(testdata$`2014_IAH3`[1]),unlist(testdata$`2014_ILH1`[1]))
# test1=(pheno_master[pheno_master$env=="2014_IAH1a",4])
# test2=(pheno_master[pheno_master$env=="2014_NYH2",4])
# test3=(pheno_master[pheno_master$env=="2014_IAH3",4])
# test4=(pheno_master[pheno_master$env=="2014_ILH1",4])
# # commonlist=Reduce(intersect, list(test1,test2,test3,test4))
###########################select IA, NY, and IL env as test data (6)###############
pheno_selected=droplevels(pheno_master[pheno_master$env%in%c("2014_IAH1a","2014_NYH2","2014_IAH3","2014_ILH1"),])
colnames(pheno_selected)[20]="yield"
################calcuate BLUPs###################
lm =lmer(yield~(1|Pedigree)+(1|env)+(1|env:Pedigree)+(1|Block),data=pheno_selected)
lmblup = ranef(lm)
temp=cbind(rownames(lmblup$`env:Pedigree`),lmblup$`env:Pedigree`)

####################################################################################################
############--------Filter genotypic data and create genotypic data for hybrids--------#############
####################################################################################################
library(NAM)
library(adegenet)
library(ggplot2)
##################################import genotypic data#########################
G = Import_data(file.choose(),type = "HapMap")
gen = G$gen
chr = G$chr
#######################impute genotypic data (forwards markov model)#############################
G_imputed=markov(gen,chr)
#####################################filter SNPs##############################
maize_G=as.data.frame(G_imputed)
G_info=read.csv(file.choose(), sep="\t",stringsAsFactors=FALSE,check.names = F)[-(2:11)]
line_tailered=toupper(gsub(":.*","",colnames(G_info[,-1])))
entry_index=which(!duplicated(gsub(":.*","",line_tailered)))
G_filtered=maize_G[entry_index,]
line_filtered=line_tailered[entry_index]
line_filtered=toupper(line_filtered)
rownames(G_filtered)=line_filtered
colnames(G_filtered)=G_info[,1]
#####################################match genomic and phenotypic data##############################
pheno_BLUP=temp[-(grep("CHECK",temp[,1])),][,1:2]
colnames(pheno_BLUP)=c("line","yield")
rownames(pheno_BLUP)=NULL
pheno_BLUP=as.data.frame(cbind(str_split_fixed(pheno_BLUP$line, ":",2),pheno_BLUP[,2]))
colnames(pheno_BLUP)=c("env","line","yield")
datagFULL=pheno_BLUP[order(pheno_BLUP$line),]
datagFULL=datagFULL[order(datagFULL$env),]
datagFULL$yield=as.numeric(as.character(datagFULL$yield))
dataEnv1=datagFULL[datagFULL$env=="2014_IAH1a",]
dataEnv2=datagFULL[datagFULL$env=="2014_NYH2",]
dataEnv3=datagFULL[datagFULL$env=="2014_IAH3",]
dataEnv4=datagFULL[datagFULL$env=="2014_ILH1",]
testa=(datagFULL[datagFULL$env=="2014_IAH1a",2])
testb=(datagFULL[datagFULL$env=="2014_NYH2",2])
testc=(datagFULL[datagFULL$env=="2014_IAH3",2])
testd=(datagFULL[datagFULL$env=="2014_ILH1",2])
corelist=Reduce(intersect, list(testa,testb,testc,testd))
pheno_BLUP_1=pheno_BLUP[pheno_BLUP$line%in%corelist,]
pheno_list_1=as.data.frame(str_split_fixed(unique(pheno_BLUP_1$line), "/",2))
pheno_BLUP_2=pheno_BLUP_1[((as.character(pheno_list_1$V1))%in%line_filtered)&((as.character(pheno_list_1$V2))%in%line_filtered),]
pheno_list_2=as.data.frame(str_split_fixed(unique(pheno_BLUP_2$line), "/",2))
geno_1=G_filtered[line_filtered%in%(c(as.character(pheno_list_2$V1),as.character(pheno_list_2$V2))),]
#####################################create new SNP data for hybrids##############################
geno_2=matrix(data=NA,nrow=nrow(pheno_list_2),ncol=ncol(geno_1))
for (i in 1:nrow(pheno_list_2)){
  print(i)  
  geno_temp=rbind(unlist(geno_1[rownames(geno_1)==pheno_list_2[i,1],]),unlist(geno_1[rownames(geno_1)==pheno_list_2[i,2],]))
  geno_2[i,]=colSums(geno_temp)/2
}
rownames(geno_2)=droplevels(unique(pheno_BLUP_2$line))
markerready=as.data.frame(geno_2)
colnames(markerready)=G_info[,1]
#####################################composit phenotypic data#####################################
dataEnv1=dataEnv1[match(rownames(markerready),dataEnv1$line),]
rownames(dataEnv1)=NULL
dataEnv1$line=rownames(markerready)
dataEnv2$env=c("2014_IAH1a")
dataEnv2=dataEnv2[match(rownames(markerready),dataEnv2$line),]
rownames(dataEnv2)=NULL
dataEnv2$line=rownames(markerready)
dataEnv2$env=c("2014_NYH2")
dataEnv3=dataEnv3[match(rownames(markerready),dataEnv3$line),]
rownames(dataEnv3)=NULL
dataEnv3$line=rownames(markerready)
dataEnv3$env=c("2014_IAH3")
dataEnv4=dataEnv4[match(rownames(markerready),dataEnv4$line),]
rownames(dataEnv4)=NULL
dataEnv4$line=rownames(markerready)
dataEnv4$env=c("2014_ILH1")
newdatagFULL=rbind(dataEnv1,dataEnv2,dataEnv3,dataEnv4)
##############################stratify population##############################
S1 <- nrow(markerready)
PopSNP <- vector(mode="list",S1)
for(i in 1:S1){
  listdata=as.integer(markerready[i,])
  PopSNP[[i]] <- new("SNPbin", listdata)
  print(i)
}
Popgen=new("genlight",PopSNP,pop=rownames(markerready))
grp=find.clusters(Popgen,stat=c("BIC"))
newdatagFULL$cluster=rep(as.numeric(grp$grp),4)
#####################################PCA plot#####################################################
pca_matrix = scale(markerready)
pca_matrix_noNA = apply(pca_matrix, 2, function (x) {ifelse(is.na(x), 0, x)})
eig = prcomp(pca_matrix_noNA, center = T)
eigenvalues = (eig$sdev)^2
pc1_var = round((eigenvalues[1]/sum(eigenvalues)) * 100, 1)
pc2_var = round((eigenvalues[2]/sum(eigenvalues)) * 100, 1)
pc3_var = round((eigenvalues[3]/sum(eigenvalues)) * 100, 1)
scores = eig$x
PCV1=scores[,1]
PCV2=scores[,2]
PCV3=scores[,3]
GSpcadata=cbind(PCV1,PCV2,PCV3)
GSpcadata=as.data.frame(cbind(GSpcadata,grp$grp))
GSpcadata$cluster=as.factor(GSpcadata$V4)
group.col=c("Blue","Red","orange","darkslategrey","lightseagreen","darkgoldenrod3","green3","cyan","burlywood4","darkred","orangered","purple","ivory4","plum4","black","springgreen2","hotpink4","lightskyblue4","darkorange","tan3","orchid4","sienna4")
png(file="PCA_pop_1vs2.png",width=10,height=8,units="in",res=500)
q=ggplot(GSpcadata,aes(x=PCV1, y=PCV2,color=cluster))+
  geom_point(size=5,stroke = 0.5, aes(fill =cluster))+theme_bw()+xlab(paste("PC1(",pc1_var,"%)",sep=""))+ylab(paste("PC2(",pc2_var,"%)",sep=""))+
  theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 1),
        axis.title=element_text(size=16,face="bold"))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  theme(legend.justification=c(1,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())
q
dev.off()
png(file="PCA_pop_2vs3.png",width=10,height=8,units="in",res=500)
q=ggplot(GSpcadata,aes(x=PCV2, y=PCV3,color=cluster))+
  geom_point(size=5,stroke = 0.5, aes(fill =cluster))+theme_bw()+xlab(paste("PC2(",pc2_var,"%)",sep=""))+ylab(paste("PC3(",pc3_var,"%)",sep=""))+
  theme(axis.text=element_text(size=10),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  theme(axis.line = element_line(color="black", size = 1),
        axis.title=element_text(size=16,face="bold"))+
  theme(panel.border= element_blank(),panel.background = element_blank())+
  theme(legend.justification=c(1,1),legend.title=element_text(size=15),
        legend.text=element_text(size=18),legend.box.just = "left")+
  theme(legend.key = element_blank())
q 
dev.off()

##################################################################################################
############################------Bayeisan model using BGLR package------#########################
##################################################################################################
###############################  setup iteration number########################
nIter<-15000; burnIn<-5000
##############################modeling with stratified different snp number##########
library(BGLR)
setwd("work space for batch data")
X<-as.matrix(markerready)
y<-newdatagFULL$yield
############# choose snp number included and start############
nofsnps=80556
nrep=5
corrsnew_env1<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env2<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env3<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env4<-matrix(data=NA,nrow=nrep,ncol=2)
predy5nrep<-matrix(data=NA,nrow=length(y),ncol=nrep*5)
# number of reps/iterations (5 in the example(z in 1:5))
for (z in 1:nrep){
  # create progress bar
  pb <- txtProgressBar(min = 1, max =nrep, style = 3)
  print(paste("sample snp",z,sep=" "))
  Xnew=X[,sample(ncol(X),nofsnps,replace=TRUE)]
  Xnew<-scale(Xnew,center=F,scale=T)
  G<-tcrossprod(Xnew)/ncol(Xnew)
  G0=kronecker(matrix(nrow=nlevels(newdatagFULL$env),ncol=nlevels(newdatagFULL$env),1),G)
  ###G+G*Env
  ETA<-list(list(~factor(newdatagFULL$env)-1,model='FIXED'),
                 list(K=G0, model='RKHS'))
  for(i in 1:4){
    tmp <- rep(0,4) ; tmp[i] <- 1
    G1 <- kronecker(diag(tmp),G)
    ETA[[(i+2)]] <- list(K=G1, model='RKHS')
  }
  #################actual 5-fold cross validation#####################
  corrs5<-matrix(data=NA,nrow=5,ncol=nlevels(newdatagFULL$env))
  predy5=matrix(data=NA,nrow=length(y),ncol=5)
  for (l in 1:5) {
    print(paste("sample lines",l,sep=" "))
    n<-length(rownames(markerready))
    group<-vector(mode="numeric", length=n)
        for (m in 1:10){
      tempphenodata=newdatagFULL[as.numeric(newdatagFULL$env)==1,]
      nincluster=length(which(tempphenodata$cluster==m))
      floorno=floor(nincluster/5)
      samplepool=c(rep(1:5,each=floorno),sample(1:5,(nincluster-floorno*5),replace=F))
      listemp=sample(samplepool,nincluster,replace=F)
      group[which(tempphenodata$cluster==m)]=listemp
    }
    newdatagFULL$group=rep(group,nlevels(newdatagFULL$env))
    ypred_cv<-matrix(data=NA,nrow=nrow(newdatagFULL),ncol=1)
    
    for (k in 1:5)
    { 
      # Reading Phenotypic data
      ycv<-y
            for (j in 1:nrow(newdatagFULL)) {
        if(newdatagFULL$group[j] == k) { ycv[j]<-NA } 
      }
            gsBayesBcv<-BGLR(y=ycv,ETA=ETA,nIter=nIter,burnIn=burnIn,verbose=F)
      predGScv<-gsBayesBcv$yHat
            for (j in 1:nrow(newdatagFULL)) {
        if(newdatagFULL$group[j] == k) { ypred_cv[j]<-predGScv[j] } 
      }
      print(k)
    }
    (corr_cv_env1<-cor(y[which(newdatagFULL$env=="2014_IAH1a")],
                       ypred_cv[which(newdatagFULL$env=="2014_IAH1a")],method='pearson',use="complete.obs"))
    (corr_cv_env2<-cor(y[which(newdatagFULL$env=="2014_NYH2")],
                       ypred_cv[which(newdatagFULL$env=="2014_NYH2")],method='pearson',use="complete.obs"))
    (corr_cv_env3<-cor(y[which(newdatagFULL$env=="2014_IAH3")],
                       ypred_cv[which(newdatagFULL$env=="2014_IAH3")],method='pearson',use="complete.obs"))
    (corr_cv_env4<-cor(y[which(newdatagFULL$env=="2014_ILH1")],
                       ypred_cv[which(newdatagFULL$env=="2014_ILH1")],method='pearson',use="complete.obs"))
    corrs5[l,]<-cbind(corr_cv_env1,corr_cv_env2,corr_cv_env3,corr_cv_env4)
    predy5[,l]=ypred_cv
  }
  corrsnew_env1[z,1]=mean(corrs5[,1])
  corrsnew_env1[z,2]=sd(corrs5[,1])
  corrsnew_env2[z,1]=mean(corrs5[,2])
  corrsnew_env2[z,2]=sd(corrs5[,2])
  corrsnew_env3[z,1]=mean(corrs5[,3])
  corrsnew_env3[z,2]=sd(corrs5[,3])
  corrsnew_env4[z,1]=mean(corrs5[,4])
  corrsnew_env4[z,2]=sd(corrs5[,4])
  predy5nrep[,(z*5-4):(z*5)]=predy5
  setTxtProgressBar(pb, z)
}
########mean accuracies for each environment########################
meanEnv1=mean(corrsnew_env1[,1])
meanEnv2=mean(corrsnew_env2[,1])
meanEnv3=mean(corrsnew_env3[,1])
meanEnv4=mean(corrsnew_env4[,1])
myPA=t(rbind(levels(pheno_BLUP_2$env),c(meanEnv1,meanEnv2,meanEnv3,meanEnv4)))
write.csv(cbind(newdatagFULL$line,as.character(newdatagFULL$env),newdatagFULL$yield,predy5nrep),row.names = F,"predy_GBLUP_80556_GY.csv")


##################################################################################################
#################################-----mixed model using rrblup package----########################
##################################################################################################
library(rrBLUP)
Pheno=newdatagFULL
y =Pheno$yield
G_ready=as.matrix(markerready)
#########choose your snp number and start#############################
nofsnps=80556
nrep=5
corrsnew_env1<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env2<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env3<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env4<-matrix(data=NA,nrow=nrep,ncol=2)
###########sample SNPS##############
predy5nrep<-matrix(data=NA,nrow=length(y),ncol=nrep*5)
for (z in 1:nrep){
  pb <- txtProgressBar(min = 1, max =nrep, style = 3)
  print(paste("sample snp",z,sep=" "))
  G_readynew=G_ready[,sample(ncol(G_ready),nofsnps,replace=TRUE)]
  corrs5<-matrix(data=NA,nrow=5,ncol=4)
  predy5=matrix(data=NA,nrow=length(y),ncol=5)
  for (l in 1:5) {
    print(paste("sample lines",l,sep=" "))
    group<-vector(mode="numeric",length=length(rownames(markerready)))
    for (m in 1:10){
      tempphenodata=newdatagFULL[as.numeric(newdatagFULL$env)==1,]
      nincluster=length(which(tempphenodata$cluster==m))
      floorno=floor(nincluster/5)
      samplepool=c(rep(1:5,each=floorno),sample(1:5,(nincluster-floorno*5),replace=F))
      listemp=sample(samplepool,nincluster,replace=F)
      group[which(tempphenodata$cluster==m)]=listemp
    }
    newdatagFULL$group=rep(group,nlevels(newdatagFULL$env))
    n=nrow(tempphenodata)
    ypred_cv<-matrix(data=NA,nrow=nrow(newdatagFULL),ncol=1)
    for (o in c(levels(newdatagFULL$env))){
      print(paste("env",o,sep=" "))
      cor_cv=c()
      ypred_cv_env<-matrix(data=NA,nrow=n,ncol=1)
      pheno_env=newdatagFULL[newdatagFULL$env==o,]
      for (k in 1:5)
      { 
        Post_trn_pheno=which(pheno_env$group!=k)
        trait=pheno_env[Post_trn_pheno,colnames(newdatagFULL)=="yield"]
        m_train=G_readynew[Post_trn_pheno,]
        m_valid=G_readynew[-Post_trn_pheno,]
        #rr-BLUP
        trait_M=mixed.solve(trait, Z=m_train, K=NULL,SE=F,return.Hinv=F)
        traitu=trait_M$u
        e=as.matrix(traitu)
        pred_trait_valid=m_valid%*%e
        pred_trait=as.vector(pred_trait_valid[,1])+as.vector(trait_M$beta)
        trait_valid=trait=pheno_env[-Post_trn_pheno,colnames(pheno_env)=="yield"]
        trait_accuracy=cor(trait_valid,pred_trait_valid,use="complete")
        ypred_cv_env[which(pheno_env$group==k),1]=pred_trait
        cor_cv[k]=trait_accuracy
        print(k)
      }
      ypred_cv[which(newdatagFULL$env==o),1]=ypred_cv_env
    }
    (corr_cv_env1<-cor(y[which(newdatagFULL$env=="2014_IAH1a")],
                       ypred_cv[which(newdatagFULL$env=="2014_IAH1a")],method='pearson',use="complete.obs"))
    (corr_cv_env2<-cor(y[which(newdatagFULL$env=="2014_NYH2")],
                       ypred_cv[which(newdatagFULL$env=="2014_NYH2")],method='pearson',use="complete.obs"))
    (corr_cv_env3<-cor(y[which(newdatagFULL$env=="2014_IAH3")],
                       ypred_cv[which(newdatagFULL$env=="2014_IAH3")],method='pearson',use="complete.obs"))
    (corr_cv_env4<-cor(y[which(newdatagFULL$env=="2014_ILH1")],
                       ypred_cv[which(newdatagFULL$env=="2014_ILH1")],method='pearson',use="complete.obs"))
    corrs5[l,]<-cbind(corr_cv_env1,corr_cv_env2,corr_cv_env3,corr_cv_env4)
    predy5[,l]=ypred_cv
  }
  corrsnew_env1[z,1]=mean(corrs5[,1])
  corrsnew_env1[z,2]=sd(corrs5[,1])
  corrsnew_env2[z,1]=mean(corrs5[,2])
  corrsnew_env2[z,2]=sd(corrs5[,2])
  corrsnew_env3[z,1]=mean(corrs5[,3])
  corrsnew_env3[z,2]=sd(corrs5[,3])
  corrsnew_env4[z,1]=mean(corrs5[,4])
  corrsnew_env4[z,2]=sd(corrs5[,4])
  predy5nrep[,(z*5-4):(z*5)]=predy5
  setTxtProgressBar(pb, z)
}
########mean accuracies for each environment########################
meanEnv1=mean(corrsnew_env1[,1])
meanEnv2=mean(corrsnew_env2[,1])
meanEnv3=mean(corrsnew_env3[,1])
meanEnv4=mean(corrsnew_env4[,1])
myPA=t(rbind(levels(pheno_BLUP_2$env),c(meanEnv1,meanEnv2,meanEnv3,meanEnv4)))
write.csv(cbind(newdatagFULL$line,as.character(newdatagFULL$env),newdatagFULL$yield,predy5nrep),row.names = F,"predy_rrblup_80556_GY.csv")


##################################################################################################
#####################-------Machine learning model using rrblup package-------####################
##################################################################################################
library(tensorflow)
library(keras)
#########Genomic data and phenotipic data#############################
Pheno=newdatagFULL
y =Pheno$yield
G_ready=as.matrix(markerready)
#########choose your snp number and start#############################
nofsnps=80556
nrep=5
corrsnew_env1<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env2<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env3<-matrix(data=NA,nrow=nrep,ncol=2)
corrsnew_env4<-matrix(data=NA,nrow=nrep,ncol=2)
###########sample SNPS##############
predy5nrep<-matrix(data=NA,nrow=length(y),ncol=nrep*5)
for (z in 1:nrep){
  pb <- txtProgressBar(min = 1, max =nrep, style = 3)
  print(paste("sample snp",z,sep=" "))
  G_readynew=G_ready[,sample(ncol(G_ready),nofsnps,replace=TRUE)]
  G_readynew<-scale(G_readynew,center=F,scale=T)
  G<-tcrossprod(G_readynew)/ncol(G_readynew)
  ########Creating the desing matrices ########################
  LG=t(chol(G))
  Z1G=model.matrix(~0+as.factor(Pheno$line))
  ZE=model.matrix(~0+as.factor(Pheno$env))
  Z1G=Z1G%*%LG
  Z2GE=model.matrix(~0+as.factor(Pheno$line):as.factor(Pheno$env))
  G2=kronecker(diag(4),data.matrix(G))
  LG2=t(chol(G2))
  Z2GE=Z2GE%*%LG2
  ###Defining the number of epoch and units#####################
  units_M=50
  epochs_M=20
  X = cbind(ZE, Z1G, Z2GE)
  #############sample lines and cross validation########################
  n=dim(X)[1]
  corrs5<-matrix(data=NA,nrow=5,ncol=nlevels(newdatagFULL$env))
  predy5=matrix(data=NA,nrow=length(y),ncol=5)
  for (l in 1:5) {
    print(paste("sample lines",l,sep=" "))
    group<-vector(mode="numeric",length=length(rownames(markerready)))
    for (m in 1:10){
      tempphenodata=newdatagFULL[as.numeric(newdatagFULL$env)==1,]
      nincluster=length(which(tempphenodata$cluster==m))
      floorno=floor(nincluster/5)
      samplepool=c(rep(1:5,each=floorno),sample(1:5,(nincluster-floorno*5),replace=F))
      listemp=sample(samplepool,nincluster,replace=F)
      group[which(tempphenodata$cluster==m)]=listemp
    }
    newdatagFULL$group=rep(group,nlevels(newdatagFULL$env))
    ypred_cv<-matrix(data=NA,nrow=nrow(newdatagFULL),ncol=1)
    
    cor_cv=c()
    ypred_cv=matrix(data=NA,nrow=nrow(newdatagFULL),ncol=1)
    for (k in 1:5)
    { 
      Post_trn=which(newdatagFULL$group!=k) ##80%training###
      X_tr = X[Post_trn,]
      X_ts = X[-Post_trn,]
      y_tr = scale(y[Post_trn])
      Mean_trn=mean(y[Post_trn])
      SD_trn=sd(y[Post_trn])
      y_ts = (y[-Post_trn]- Mean_trn)/SD_trn
      #########Model fitting in Keras################################
      #install_ensorflow()
      model= keras_model_sequential()
      #########Layers specification ################################
      model %>%
        layer_dense(
          units =units_M,
          activation = "relu",
          input_shape = c(dim(X_tr)[2])) %>%
        layer_dropout(rate = 0.3) %>% ###Input Layer
        layer_dense(units = units_M, activation = "relu") %>%
        layer_dropout(rate = 0.3) %>% ###Hidden layer 1
        layer_dense(units = units_M, activation =  "relu") %>%
        layer_dropout(rate = 0.3) %>% ####Hidden layer 2
        layer_dense(units = 1) ####Output layer
      model %>% compile(
        loss = "mean_squared_error",
        optimizer = optimizer_adam(),
        metrics = c("mean_squared_error"))
      history <- model %>% fit(
        X_tr, y_tr, epochs = epochs_M, batch_size = 30,
        verbose = FALSE)
      #######Evaluating the performance of the model###################
      pf = model %>% evaluate(x = X_ts, y = y_ts, verbose = 0)
      y_p = model %>% predict(X_ts)
      y_p=y_p*SD_trn+ Mean_trn
      y_ts=y_ts
      y_ts=y_ts*SD_trn+ Mean_trn
      ###############Observed and predicted values of the testing set#
      Y_all_tst = data.frame(cbind(y_ts, y_p))
      cor_cv[k]=cor(Y_all_tst[,1],Y_all_tst[,2])
      ypred_cv[which(newdatagFULL$group==k),]=y_p
      print(k)
    }
    
    (corr_cv_env1<-cor(y[which(newdatagFULL$env=="2014_IAH1a")],
                       ypred_cv[which(newdatagFULL$env=="2014_IAH1a")],method='pearson',use="complete.obs"))
    (corr_cv_env2<-cor(y[which(newdatagFULL$env=="2014_NYH2")],
                       ypred_cv[which(newdatagFULL$env=="2014_NYH2")],method='pearson',use="complete.obs"))
    (corr_cv_env3<-cor(y[which(newdatagFULL$env=="2014_IAH3")],
                       ypred_cv[which(newdatagFULL$env=="2014_IAH3")],method='pearson',use="complete.obs"))
    (corr_cv_env4<-cor(y[which(newdatagFULL$env=="2014_ILH1")],
                       ypred_cv[which(newdatagFULL$env=="2014_ILH1")],method='pearson',use="complete.obs"))
    
    corrs5[l,]<-cbind(corr_cv_env1,corr_cv_env2,corr_cv_env3,corr_cv_env4)
    predy5[,l]=ypred_cv
  }
  corrsnew_env1[z,1]=mean(corrs5[,1])
  corrsnew_env1[z,2]=sd(corrs5[,1])
  
  corrsnew_env2[z,1]=mean(corrs5[,2])
  corrsnew_env2[z,2]=sd(corrs5[,2])
  
  corrsnew_env3[z,1]=mean(corrs5[,3])
  corrsnew_env3[z,2]=sd(corrs5[,3])
  
  corrsnew_env4[z,1]=mean(corrs5[,4])
  corrsnew_env4[z,2]=sd(corrs5[,4])
  
  predy5nrep[,(z*5-4):(z*5)]=predy5
  setTxtProgressBar(pb, z)
  
}
########mean accuracies for each environment########################
meanEnv1=mean(corrsnew_env1[,1])
meanEnv2=mean(corrsnew_env2[,1])
meanEnv3=mean(corrsnew_env3[,1])
meanEnv4=mean(corrsnew_env4[,1])
myPA=t(rbind(levels(pheno_BLUP_2$env),c(meanEnv1,meanEnv2,meanEnv3,meanEnv4)))
write.csv(cbind(newdatagFULL$line,as.character(newdatagFULL$env),newdatagFULL$yield,predy5nrep),row.names = F,"predy_DL_80556_GY.csv")
