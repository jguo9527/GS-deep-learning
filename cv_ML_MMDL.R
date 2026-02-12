##################################################################################################
#####################-------Machine learning model using tensorflow package-------####################
##################################################################################################
#install_tensorflow(version = "gpu")
library(tensorflow)
library(beepr)
library(keras)
k_clear_session()
library(mice)
options(scipen=999)
#########Genomic data and phenotipic data#############################
newdatagFULL=read.csv(file.choose(),row.names = 1)
markerready=read.table(file.choose())
Pheno=newdatagFULL
str(Pheno)
pheno_impute=mice(as.data.frame(Pheno[,-c(1,2,9)]), m=5, maxit = 50, method = 'pmm', seed = 500)
pheno_imputed=complete(pheno_impute,5)
y=pheno_imputed[,c(19,7,10,8)]
y_imputed=y
#y=Pheno[,c(22,10,13,11)]
nofsnps=27957
G_ready=as.matrix(markerready)
G_readynew=G_ready[,sample(ncol(G_ready),nofsnps,replace=TRUE)]
G_readynew<-scale(G_readynew,center=F,scale=T)
G<-tcrossprod(G_readynew)/ncol(G_readynew)
########Creating the desing matrices ########################
LG=t(chol(G))
Z1G=model.matrix(~0+as.factor(Pheno$ENTRY))
ZE=model.matrix(~0+as.factor(Pheno$Env))
Z1G=Z1G%*%LG
Z2GE=model.matrix(~0+as.factor(Pheno$ENTRY):as.factor(Pheno$Env))
G2=kronecker(diag(4),data.matrix(G))
LG2=t(chol(G2))
Z2GE=Z2GE%*%LG2
X = cbind(ZE, Z1G, Z2GE)
n=dim(X)[1]
nt=dim(y)[2]
nenv=length(unique(newdatagFULL$Env))
###Defining the number of epoch and units####################
stratified=T
nrep=25
#my_units_M=c(seq(from=10,to=100,by=10))
#my_epochs_M=c(seq(from=20,to=100,by=20))

my_epochs_M=51
my_units_M=32

my_units_level=nlevels(as.factor(my_units_M))
my_epochs_level=nlevels(as.factor(my_epochs_M))
epochs_curve=as.data.frame(matrix(data=NA,nrow=0,ncol=7))
myPA=matrix(data=NA,nrow=nenv*nt*2,ncol=1)
#pb <- txtProgressBar(min = 1, max =nlevels(as.factor(my_epochs_M)), style = 3)
for (t in 1:my_epochs_level){
  epochs_M=my_epochs_M[t]
  for (s in 1:nlevels(as.factor(my_units_M))){
    units_curve=as.data.frame(matrix(data=NA,nrow=0,ncol=7))
    units_M=my_units_M[s]
    corrs<-as.data.frame(matrix(data=NA,nrow=nenv*nt,ncol=nrep))
    predy<-as.data.frame(matrix(data=NA,nrow=nrow(y),ncol=nrep*nt))
    for (l in 1:nrep){
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
        group=sample(1:5,length(rownames(markerready)),replace = T)
      }
       newdatagFULL$group=rep(group,length(unique(newdatagFULL$Env)))
        ypred_cv=matrix(data=NA,nrow=nrow(newdatagFULL),ncol=nt)
        ###############################cross-validation############################    
        for (k in 1:5)
        { 
          Post_trn=which(newdatagFULL$group!=k) ##80%training###
          X_tr = X[Post_trn,]
          X_ts = X[-Post_trn,]
          y_tr = scale(y[Post_trn,])
          Mean_trn=apply(y[Post_trn,],2,mean,na.rm=TRUE)
          SD_trn=apply(y[Post_trn,],2,sd,na.rm=TRUE)
          
          y_ts=matrix(NA,ncol=nt,nrow=dim(X_ts)[1])
          
          for(u in 1:nt){
            y_ts[,u] = (y[-Post_trn,u]- Mean_trn[u])/SD_trn[u]
          }
          #########Model fitting in Keras################################
          #install_tensorflow()
          model= keras_model_sequential()
          #########Layers specification ################################
          # add covariates (independent variables)
          input <- layer_input(shape=dim(X_tr)[2],name="covars")
          # add hidden layers
          base_model <- input %>%
            layer_dense(units = units_M, activation ="relu") %>%
            layer_dropout(rate = 0.25) %>%
            layer_dense(units = units_M, activation ="relu") %>%
            layer_dropout(rate = 0.25) %>%
            layer_dense(units = units_M, activation ="relu") %>%
            layer_dropout(rate = 0.25) %>%
            layer_dense(units = units_M, activation ="relu") %>%
            layer_dropout(rate = 0.25)
          # add output 1
          yhat1 <- base_model %>%
            layer_dense(units = 1, name="yhat1")
          # add output 2
          yhat2 <- base_model %>%
            layer_dense(units = 1, name="yhat2")
          # add output 3
          yhat3 <- base_model %>%
            layer_dense(units = 1, name="yhat3")
          # add output 4
          yhat4 <- base_model %>%
            layer_dense(units = 1, name="yhat4")
          # build multi-output model
          model <- keras_model(input,list(yhat1,yhat2,yhat3,yhat4)) %>%
            compile(optimizer = "rmsprop",
                    loss="mse",
                    metrics="mae",
                    loss_weights=c(0.25,0.25,0.25,0.25))
          # fit model
          model_fit <- model %>%
            fit(x=X_tr,
                y=list(y_tr[,1],y_tr[,2],y_tr[,3],y_tr[,4]),
                epochs=epochs_M,
                batch_size = 50,
                verbose=0)
          # predict values for test set
          Yhat <- predict(model, X_ts) %>%
            data.frame() %>%
            setNames(colnames(y_tr))
          predB=Yhat
          y_p=predB
          for (v in 1:nt){
            y_p[,v]=y_p[,v]*SD_trn[v]+ Mean_trn[v]
            y_ts[,v]=y_ts[,v]*SD_trn[v]+ Mean_trn[v]
          }
          ###############Observed and predicted values of the testing set#
          ypred_cv[which(newdatagFULL$group==k),]=as.matrix(y_p)
        }
        Y_all_tst=cbind(y,ypred_cv)
        corr_env1<-cor(Y_all_tst[which(Pheno$Env=="2016_1"),],method='pearson',use="complete.obs")
        corr_cv_env1=diag(corr_env1[(nt+1):(2*nt),1:nt])
        corr_env2<-cor(Y_all_tst[which(Pheno$Env=="2016_2"),],method='pearson',use="complete.obs")
        corr_cv_env2=diag(corr_env2[(nt+1):(2*nt),1:nt])
        corr_env3<-cor(Y_all_tst[which(Pheno$Env=="2017_1"),],method='pearson',use="complete.obs")
        corr_cv_env3=diag(corr_env3[(nt+1):(2*nt),1:nt])
        corr_env4<-cor(Y_all_tst[which(Pheno$Env=="2017_2"),],method='pearson',use="complete.obs")
        corr_cv_env4=diag(corr_env4[(nt+1):(2*nt),1:nt])
        
         units_curve[l,4:(3+nt)]=colMeans(as.data.frame(rbind(corr_cv_env2,corr_cv_env4,corr_cv_env1,corr_cv_env3)))
         corrs[,l]<-c(stack(as.data.frame(rbind(corr_cv_env2,corr_cv_env4,corr_cv_env1,corr_cv_env3)))[1])
        predy[,(nt*l-3):(nt*l)]=as.matrix(ypred_cv)
      #corrsnew[1:16,l]=rowMeans(rep)
      #corrsnew[17:32,l]=c(apply(rep,1, sd))
    }
    units_curve[,3]=seq(from=1,to=nrep,by=1)
    units_curve[,2]=my_units_M[s]
    units_curve[,1]=my_epochs_M[t]
    epochs_curve=rbind(epochs_curve,units_curve)
    k_clear_session()
    print(paste("neuron",my_units_M[s],"epochs",my_epochs_M[t],sep=" "))
  }
   # setTxtProgressBar(pb, t)  
}
myPA[1:(nenv*nt),1]=c(apply(corrs,1, mean))
myPA[((nenv*nt)+1):(nenv*nt*2),1]=c(apply(corrs,1, sd))
beep("fanfare")


# write.csv(epochs_curve, row.names = F,"curve_MMDL_stratfied.csv")
# write.csv(myPA, row.names = F,"PA_MMDL_stratfied.csv")
write.csv(cbind(newdatagFULL$ENTRY,as.character(newdatagFULL$Env),y,predy),
          row.names = F,"predy_MMDL_stratfied(51,32).csv")
