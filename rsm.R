library("rsm")
ML_data=read.csv(file.choose())
colnames(ML_data)=c("epochs","neuron","rep","gy","HI","SF","TGW")
table(ML_data$neuron)
table(ML_data$epochs)
#str(ML_data)
ML_RSM_gy <- rsm(gy ~ SO(epochs, neuron), 
                 data = ML_data[(ML_data$epochs<50)&(ML_data$neuron<55),])
ML_RSM_HI <- rsm(HI ~ SO(epochs, neuron), 
                 data = ML_data[(ML_data$epochs<60)&(ML_data$neuron<55),])
ML_RSM_SF <- rsm(SF ~ SO(epochs, neuron), 
                 data = ML_data[ML_data$epochs<60,])
ML_RSM_TGW <- rsm(TGW ~ SO(epochs, neuron), 
                  data = ML_data[ML_data$epochs<60,])
summary(ML_RSM_gy)
summary(ML_RSM_HI)
summary(ML_RSM_SF)
summary(ML_RSM_TGW)
#contour(ML_RSM, ~ epochs+neuron, image = TRUE,at = summary(ML_RSM)$canonical$'xs')


png(file="SM_stratified.png",width=12.57703333,height=8,units="in",res=300)
par(mfrow = c(2,2),
    oma = c(0,0,0,0) + 0.1,
    mar = c(2,1,0,0) + 0.1)
persp(ML_RSM_gy, ~ epochs+neuron, image = TRUE,
      at = c(summary(ML_RSM_gy)$canonical$'xs'),
      theta=30,col.lab=33,contour="colors",cex.lab=1.3,cex.axis=1.2)
persp(ML_RSM_HI, ~ epochs+neuron, image = TRUE,
      at = c(summary(ML_RSM_HI)$canonical$'xs'),
      theta=30,col.lab=33,contour="colors",cex.lab=1.3,cex.axis=1.2)
persp(ML_RSM_SF, ~ epochs+neuron, image = TRUE,
      at = c(summary(ML_RSM_SF)$canonical$'xs'),
      theta=30,col.lab=33,contour="colors",cex.lab=1.3,cex.axis=1.2)
persp(ML_RSM_TGW, ~ epochs+neuron, image = TRUE,
      at = c(summary(ML_RSM_TGW)$canonical$'xs'),
      theta=30,col.lab=33,contour="colors",cex.lab=1.3,cex.axis=1.2)
dev.off()



steepest(ML_RSM)
canonical.path(ML_RSM, dist = seq(-10, 10, by = 0.5))


ag1 <- aggregate(. ~ epochs+neuron, ML_data[,-3], function(x) c(mean = mean(x), sd = sd(x)))
ag1
ag2 <- aggregate(. ~ epochs, ML_data[,-3], function(x) c(mean = mean(x), sd = sd(x)))
ag3 <- aggregate(. ~ neuron, ML_data[,-3], function(x) c(mean = mean(x), sd = sd(x)))
