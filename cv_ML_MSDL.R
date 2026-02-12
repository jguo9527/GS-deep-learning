##################################################################################################
#####################-------Machine learning model using tensorflow package-------################
##################################################################################################
library(tensorflow)
library(keras)
library(beepr)
k_clear_session()
# Disable GPU                                                                   
Sys.setenv("CUDA_VISIBLE_DEVICES" = -1)   
library(mice)
options(scipen=999)
#########Genomic data and phenotipic data#############################
newdatagFULL=read.csv(file.choose(),check.names = F)
markerready=read.table(file.choose(),check.names=F)
Pheno=newdatagFULL
str(Pheno)
pheno_impute=mice(as.data.frame(Pheno[,-c(1,2,9)]), m=5, maxit = 50, method = 'pmm', seed = 500)
pheno_imputed=complete(pheno_impute,5)
y_imputed=pheno_imputed[,c(19,7,10,8)]
#y=Pheno[,c(22,13,10,11)]
G_ready=as.matrix(markerready)
#########choose your snp number and start#############################
nofsnps=ncol(markerready)
###########sample SNPS##############
G_readynew=G_ready[,sample(ncol(G_ready),nofsnps,replace=TRUE)]
G_readynew<-scale(G_readynew,center=F,scale=T)
G<-tcrossprod(G_readynew)/ncol(G_readynew)
corrsnew_env=matrix(data=NA,nrow=length(unique(newdatagFULL$Env)),ncol=2)
########Creating the desing matrices ########################
LG=t(chol(G))
Z1G=model.matrix(~0+as.factor(Pheno$ENTRY))
ZE=model.matrix(~0+as.factor(Pheno$Env))
Z1G=Z1G%*%LG
Z2GE=model.matrix(~0+as.factor(Pheno$ENTRY):as.factor(Pheno$Env))
G2=kronecker(diag(4),data.matrix(G))
LG2=t(chol(G2))
Z2GE=Z2GE%*%LG2
###Defining the number of epoch and units#####################
nt=ncol(y_imputed)
ne=length(unique(newdatagFULL$Env))
X = cbind(ZE, Z1G, Z2GE)
#####################################training################################################
nrep=25
stratified=T
# my_units_M=c(3,5,8,seq(from=10,to=70,by=5))
# my_epochs_M=c(3,5,8,seq(from=10,to=80,by=10))

my_epochs_M=41
my_units_M=9

my_units_level=nlevels(as.factor(my_units_M))
my_epochs_level=nlevels(as.factor(my_epochs_M))
epochs_curve=as.data.frame(matrix(data=NA,nrow=0,ncol=7))
myPA=as.data.frame(matrix(data=NA,nrow=ne*nt*2,ncol=1))
#pb <- txtProgressBar(min = 1, max =nlevels(as.factor(my_epochs_M)), style = 3)
for (t in 1:my_epochs_level){
  epochs_M=my_epochs_M[t]
  for (s in 1:nlevels(as.factor(my_units_M))){
    units_curve=as.data.frame(matrix(data=NA,nrow=0,ncol=7))
    units_M=my_units_M[s]
    predy<-as.data.frame(matrix(data=NA,nrow=nrow(X),ncol=nrep*nt))
    for (h in 1:ncol(y_imputed)){
      y=y_imputed[,h]
      rep<-matrix(data=NA,nrow=nrep,ncol=ne)
      predyrep=matrix(data=NA,nrow=length(y),ncol=25)
      for (l in 1:nrep) {
        n<-length(rownames(markerready))
        if (stratified==T){
          group<-vector(mode="numeric",length=length(rownames(markerready)))
          for (m in 1:10){
            tempphenodata=newdatagFULL[1:length(unique(newdatagFULL$ENTRY)),]
            nincluster=length(which(tempphenodata$cluster==m))
            floorno=floor(nincluster/5)
            samplepool=c(rep(1:5,each=floorno),sample(1:5,(nincluster-floorno*5),replace=F))
            listemp=sample(samplepool,nincluster,replace=F)
            group[which(tempphenodata$cluster==m)]=listemp
          }
        }else{
          group=sample(1:5,n,replace = T)
        }
        newdatagFULL$group=rep(group,ne)
        ypred_cv<-matrix(data=NA,nrow=nrow(newdatagFULL),ncol=1)
        cor_cv=c()
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
            layer_dropout(rate = 0.25) %>%
            layer_dense(units = units_M, activation ="relu") %>%
            layer_dropout(rate = 0.25) %>%
            layer_dense(units = units_M, activation ="relu") %>%
            layer_dropout(rate = 0.25) %>%
            layer_dense(units = units_M, activation ="relu") %>%
            layer_dropout(rate = 0.25) %>%
            layer_dense(units = 1)
          model %>% compile(
            loss = "mean_squared_error",
            optimizer = optimizer_adam(),
            metrics = c("mean_squared_error"))
          history <- model %>% fit(
            X_tr, y_tr, epochs = epochs_M, batch_size = 40,
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
        }
        (corr_cv_env3<-cor(y[which(newdatagFULL$Env=="2016_1")],
                           ypred_cv[which(newdatagFULL$Env=="2016_1")],method='pearson',use="complete.obs"))
        (corr_cv_env1<-cor(y[which(newdatagFULL$Env=="2016_2")],
                           ypred_cv[which(newdatagFULL$Env=="2016_2")],method='pearson',use="complete.obs"))
        (corr_cv_env4<-cor(y[which(newdatagFULL$Env=="2017_1")],
                           ypred_cv[which(newdatagFULL$Env=="2017_1")],method='pearson',use="complete.obs"))
        (corr_cv_env2<-cor(y[which(newdatagFULL$Env=="2017_2")],
                           ypred_cv[which(newdatagFULL$Env=="2017_2")],method='pearson',use="complete.obs"))
        rep[l,]<-cbind(corr_cv_env1,corr_cv_env2,corr_cv_env3,corr_cv_env4)
        predyrep[,l]=ypred_cv
      }
      myPA[(ne*h-3):(ne*h),]=c(mean(rep[,1]),mean(rep[,2]),mean(rep[,3]),mean(rep[,4]))
      myPA[(ne*h+13):(ne*h+16),]=c(sd(rep[,1]),sd(rep[,2]),sd(rep[,3]),sd(rep[,4]))
      
      predy[,(25*h-24):(25*h)]=predyrep
      
      rep_curve=rowMeans(rep)
      units_curve[1:nrep,h+3]=rep_curve
    }
    units_curve[,3]=seq(from=1,to=nrep,by=1)
    units_curve[,2]=my_units_M[s]
    units_curve[,1]=my_epochs_M[t]
    epochs_curve=rbind(epochs_curve,units_curve)
    k_clear_session()
    print(paste("neuron",my_units_M[s],"epochs",my_epochs_M[t],sep=" "))
  }
  #setTxtProgressBar(pb, t)
}
beep("fanfare")


write.csv(epochs_curve, row.names = F,"curve_SMDL_stratfied_1.csv")

write.csv(cbind(newdatagFULL$ENTRY,as.character(newdatagFULL$Env),y,predy),
          row.names = F,"predy_SMDL_stratfied(41,9).csv")

