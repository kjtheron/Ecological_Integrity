# Raster and shapefile processing
library(rgdal)
library(rgeos)
library(raster)
library(RStoolbox)
library(maptools)
library(spatialEco)
library(spdep)
library(glcm)

# Machine learning
library(caret)
library(MLmetrics)

# Parallel processing
library(doParallel)
library(snow)

# Data manipulation and plotting
library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(corrplot)
library(ggpubr)

# Multivariate analysis
library(vegan)
library(SSDM)
library(gdm)
library(ResistanceGA)
library(leastcostpath)
library(gdistance)
library(betapart)
library(eulerr)

# Spatial autocorrelation
library(ape)
library(ade4)

# Stats
library(lme4)
library(glmmADMB)
library(MuMIn)
library(robustHD)
library(ggcorrplot)
library(car)
library(vcd)

# Set working directory
setwd("~/")

############### Process Raster data ###############
# Load Sentinel data
b2<-raster("Raster_Data/Imagery/Sentinel/T36JTN_20200129T074049_B02_10m.jp2")
b3<-raster("Raster_Data/Imagery/Sentinel/T36JTN_20200129T074049_B03_10m.jp2")
b4<-raster("Raster_Data/Imagery/Sentinel/T36JTN_20200129T074049_B04_10m.jp2")
b8<-raster("Raster_Data/Imagery/Sentinel/T36JTN_20200129T074049_B08_10m.jp2")
Sent<-stack(b2,b3,b4,b8)
rm(b2,b3,b4,b8)

# Reduce extent
Estates<-readOGR("Shapefiles/Study_Extent.shp")
Sent<-crop(Sent,Estates)
Sent<-mask(Sent,Estates)
writeRaster(Sent,filename="Products/Raster_Products/Sent.tif",format="GTiff",overwrite=TRUE)

# Load and process SUDEM
SUDEM<-raster("Raster_Data/DEM/SUDEM5_10m_projected_mask.tif")
SUDEM<-projectRaster(SUDEM,Sent)
SUDEM<-crop(SUDEM,Estates)
SUDEM<-mask(SUDEM,Estates)
writeRaster(SUDEM,filename="Products/Raster_Products/SUDEM.tif",format="GTiff",overwrite=TRUE)

# Load SA_NLC
SA_NLC<-raster("Raster_Data/SA_NLC_2018/SA_NLC_2018_projected_mask_10m.tif")
SA_NLC<-projectRaster(SA_NLC,Sent)
SA_NLC<-crop(SA_NLC,Estates)
SA_NLC<-mask(SA_NLC,Estates)
writeRaster(SA_NLC,filename="Products/Raster_Products/SA_NLC.tif",format="GTiff",overwrite=TRUE)

# Load and process Drainage Raster
Drain<-raster("Raster_Data/DEM/SUDEM_Drainage_10m.tif")
Drain<-projectRaster(Drain,Sent)
Drain<-crop(Drain,Estates)
Drain<-mask(Drain,Estates)
writeRaster(Drain,filename="Products/Raster_Products/Drain.tif",format="GTiff",overwrite=TRUE)

# Load Max_NBR
Max_NBR<-raster("Raster_Data/Fire/Max_NBR_2017_2020.tif")
Max_NBR<-resample(Max_NBR,Sent)
Max_NBR<-crop(Max_NBR,Estates)
Max_NBR<-mask(Max_NBR,Estates)
writeRaster(Max_NBR,filename="Products/Raster_Products/Max_NBR.tif",format="GTiff",overwrite=TRUE)
rm(Sent,Estates,SUDEM,SA_NLC,Drain,Max_NBR)

############### Calculate raster variables ###############
Sent_SR<-brick("Products/Raster_Products/Sent.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
SA_NLC<-raster("Products/Raster_Products/SA_NLC.tif")
Drain<-raster("Products/Raster_Products/Drain.tif")

# Calculate NDVI
names(Sent_SR)<-c("Blue","Green","Red","NIR")
Sent_NDVI<-(Sent_SR[["NIR"]]-Sent_SR[["Red"]])/(Sent_SR[["NIR"]]+Sent_SR[["Red"]])
writeRaster(Sent_NDVI,filename="Products/Raster_Products/Sent_NDVI.tif",format="GTiff",overwrite=TRUE)

# PCA
PCA<-rasterPCA(Sent_SR,nComp=1)
PCA<-PCA$map
writeRaster(PCA,filename="Products/Raster_Products/PCA.tif",format="GTiff",overwrite=TRUE)

# Calculate slope and aspect
derivedVars<-terrain(SUDEM,opt=c("slope","aspect"),neighbors=8,unit="degrees")
SUDEM_Slope<-derivedVars[["slope"]]
SUDEM_Aspect<-derivedVars[["aspect"]]
writeRaster(SUDEM_Slope,filename="Products/Raster_Products/SUDEM_Slope.tif",format="GTiff",overwrite=TRUE)
writeRaster(SUDEM_Aspect,filename="Products/Raster_Products/SUDEM_Aspect.tif",format="GTiff",overwrite=TRUE)
rm(SUDEM,derivedVars)

# Calculate distance matrix
Drain[Drain>1]<-1
Estate_1<-readOGR("Shapefiles/Study_Estate1_Extent.shp")
Drain_Estate1<-crop(Drain,Estate_1)
Drain_Estate1<-mask(Drain_Estate1,Estate_1)
Drain_Estate1<-raster::distance(Drain_Estate1)
Drain_Estate1<-crop(Drain_Estate1,Sent_SR[[1]])
Estate_2<-readOGR("Shapefiles/Study_Estate2_Extent.shp")
Drain_Estate2<-crop(Drain,Estate_2)
Drain_Estate2<-mask(Drain_Estate2,Estate_2)
Drain_Estate2<-raster::distance(Drain_Estate2)
Drain_Estate2<-crop(Drain_Estate2,Sent_SR[[1]])
Estate_3<-readOGR("Shapefiles/Study_Estate3_Extent.shp")
Drain_Estate3<-crop(Drain,Estate_3)
Drain_Estate3<-mask(Drain_Estate3,Estate_3)
Drain_Estate3<-raster::distance(Drain_Estate3)
Drain_Estate3<-crop(Drain_Estate3,Sent_SR[[1]])
Estate_4<-readOGR("Shapefiles/Study_Estate4_Extent.shp")
Drain_Estate4<-crop(Drain,Estate_4)
Drain_Estate4<-mask(Drain_Estate4,Estate_4)
Drain_Estate4<-raster::distance(Drain_Estate4)
Drain_Estate4<-crop(Drain_Estate4,Sent_SR[[1]])
Drain_Dist<-merge(Drain_Estate1,Drain_Estate2,Drain_Estate3,Drain_Estate4)
Drain_Dist<-crop(Drain_Dist,Sent_SR[[1]])
Drain_Dist<-mask(Drain_Dist,Sent_SR[[1]])
writeRaster(Drain_Dist,filename="Products/Raster_Products/Drain_Dist.tif",format="GTiff",overwrite=TRUE)
rm(Estate_1,Estate_2,Estate_3,Estate_4,Drain_Estate1,Drain_Estate2,Drain_Estate3,Drain_Estate4,Drain)

# Create Stack for classification 
Raster_data<-stack(Sent_SR,Sent_NDVI,PCA,SA_NLC,SUDEM_Slope)
writeRaster(Raster_data,filename="Products/Raster_Products/Raster_data.tif",format="GTiff",overwrite=TRUE)
rm(Raster_data,Sent_SR,Sent_NDVI,PCA,SA_NLC,SUDEM_Slope,SUDEM_Aspect,Drain_Dist)

############### Classify satellite imagery ###############
#https://valentinitnelav.github.io/satellite-image-classification-r/

# Load processed data
Raster_data<-stack("Products/Raster_Products/Raster_data.tif")
names(Raster_data)<-c("Blue","Green","Red","NIR","NDVI","PCA","Land_Use","Slope")

# Add training data
poly<-readOGR("Shapefiles/Training_Data.shp",stringsAsFactors=FALSE)
setDT(poly@data)

# Manipulating training data
poly_utm<-sp::spTransform(poly,CRSobj=Raster_data[[1]]@crs)
template_rst<-raster(extent(Raster_data[[1]]),resolution=10,crs=projection(Raster_data[[1]]))
poly_utm_rst<-raster::rasterize(poly_utm,template_rst,field='ID')
poly_dt<-as.data.table(rasterToPoints(poly_utm_rst))
setnames(poly_dt,old="layer",new="ID_cls")
points<-SpatialPointsDataFrame(coords=poly_dt[,.(x,y)],data=poly_dt,proj4string=poly_utm_rst@crs)
rm(poly_utm,template_rst,poly_utm_rst)

# Extract raster values
df<-Raster_data %>%
  raster::extract(y=points) %>%  
  as.data.frame %>% 
  mutate(ID_cls=points@data$ID_cls) %>%
  left_join(y=unique(poly@data),by=c("ID_cls"="ID")) %>% 
  mutate(ID_cls=NULL) %>% 
  mutate(Class=factor(Class))
setDT(df)
rm(points,poly,poly_dt)
write.csv(df,"Products/Excel_Sheets_Products/Extracted_Raster_Data_Evaluation.csv",row.names=FALSE)

# Inspect class balance
plot(df$Class,main="Class frequencies")

# Split dataset
set.seed(321)
df_Split<-createDataPartition(df$Class,p=0.7,list=FALSE)
df_train<-df[df_Split,]
plot(df_train$Class,main="Train Class frequencies")
df_test<-df[-df_Split,]
plot(df_test$Class,main="Test Class frequencies")
rm(df_Split,df)

# Classification CART     
system.time(model_rpart<-caret::train(Class~Blue+Green+Red+NIR+NDVI+PCA+Land_Use+Slope, 
                                      method="rpart",data=na.omit(df_train),tuneLength=10,metric="Kappa",
                                      trControl=trainControl(summaryFunction=multiClassSummary)))
# Saving model
saveRDS(model_rpart,file="Products/Models/Model_rpart.rds")
# Reading model
model_rpart<-readRDS("Products/Models/Model_rpart.rds")
# Confusion matrix
predict_rpart<-predict(model_rpart,na.omit(df_test))
confusionMatrix(data=predict_rpart,na.omit(df_test$Class),mode="prec_recall")
# Variable importance
plot(varImp(model_rpart),main="rpart variable importance")

# Predict and save
rpart_map<-raster::predict(object=Raster_data,model=model_rpart,type="raw")
rpart_map<-focal(rpart_map,matrix(1,3,3),fun=modal) 
writeRaster(rpart_map,filename="Products/Raster_Products/Map_rpart.tif",format="GTiff",overwrite=TRUE)
rm(predict_rpart,rpart_map)

# Classification kknn
set.seed(321)
cl <- makeCluster(3/4 * detectCores())
registerDoParallel(cl)
system.time(model_kknn<-caret::train(Class~Blue+Green+Red+NIR+NDVI+PCA+Land_Use+Slope, 
                                     method="kknn",data=na.omit(df_train),tuneLength=10,metric="Kappa",
                                     trControl=trainControl(summaryFunction=multiClassSummary)))
stopCluster(cl); remove(cl)
registerDoSEQ()
# Saving model
saveRDS(model_kknn,file="Products/Models/Model_kknn.rds")
# Reading model
model_kknn<-readRDS("Products/Models/Model_kknn.rds")
# Confusion matrix
predict_kknn<-predict(model_kknn,na.omit(df_test))
confusionMatrix(data=predict_kknn,na.omit(df_test$Class),mode="prec_recall")
# Variable importance
plot(varImp(model_kknn),main="kknn variable importance")
# Predict and save
kknn_map<-raster::predict(object=Raster_data,model=model_kknn,type="raw")
kknn_map<-focal(kknn_map,matrix(1,3,3),fun=modal) 
writeRaster(kknn_map,filename="Products/Raster_Products/Map_kknn.tif",format="GTiff",overwrite=TRUE)
rm(predict_kknn,kknn_map)

# Classification ranger
system.time(model_ranger<-caret::train(Class~Blue+Green+Red+NIR+NDVI+PCA+Land_Use+Slope, 
                                       method="ranger",importance="impurity",data=na.omit(df_train),tuneLength=10,metric="Kappa",
                                       trControl=trainControl(summaryFunction=multiClassSummary,classProbs=TRUE)))
# Saving model
saveRDS(model_ranger,file="Products/Models/Model_ranger.rds")
# Reading model
model_ranger<-readRDS("Products/Models/Model_ranger.rds")
# Confusion matrix
predict_ranger<-predict(model_ranger,na.omit(df_test))
confusionMatrix(data=predict_ranger,na.omit(df_test$Class),mode="prec_recall")
# Variable importance
plot(varImp(model_ranger),main="ranger variable importance")
# Predict and save
ranger_map<-raster::predict(object=Raster_data,model=model_ranger,type="raw")
ranger_map<-focal(ranger_map,matrix(1,3,3),fun=modal) 
writeRaster(ranger_map,filename="Products/Raster_Products/Map_ranger.tif",format="GTiff",overwrite=TRUE)
rm(predict_ranger,ranger_map)

#Compare model performance 7x6
model_list<-list(rpart=model_rpart,kknn=model_kknn,ranger=model_ranger)
resamples<-caret::resamples(model_list)
bwplot(resamples)
rm(model_list,model_kknn,model_ranger,model_rpart,df_train,df_test,resamples,Raster_data)

############### SSDM ###############
# https://cran.r-project.org/web/packages/SSDM/SSDM.pdf
# https://cran.r-project.org/web/packages/SSDM/vignettes/SSDM.html
# https://doi.org/10.1111/2041-210X.12841

# Load and calculating response variables
Data<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
groups<-c("Low","Intermediate","Intermediate","Low","Low","High","Low","Low","Intermediate","Intermediate","Low","Low",
          "Intermediate","Low","Low","Low","Low","Intermediate","Low","Low","High","Intermediate","Low","Intermediate",
          "Intermediate","High","Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","Intermediate",
          "Low","High","Intermediate","Low","High","High","Low","Intermediate","Low","Intermediate","High","Intermediate",
          "Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","High","High","Intermediate","High",
          "Low","Intermediate","Low","Intermediate","Low")
DataT<-t(Data)
DataT<-as.data.frame(cbind(DataT,"Groups"=groups))
Low<-subset(DataT,Groups=="Low")
Low$Groups=NULL
Low<-as.data.frame(t(Low))
indx<-sapply(Low,is.factor)
Low[indx]<-lapply(Low[indx],function(x) as.numeric(as.character(x)))
Intermediate<-subset(DataT,Groups=="Intermediate")
Intermediate$Groups=NULL
Intermediate<-as.data.frame(t(Intermediate))
indx<-sapply(Intermediate,is.factor)
Intermediate[indx]<-lapply(Intermediate[indx],function(x) as.numeric(as.character(x)))
High<-subset(DataT,Groups=="High")
High$Groups=NULL
High<-as.data.frame(t(High))
indx<-sapply(High,is.factor)
High[indx]<-lapply(High[indx],function(x) as.numeric(as.character(x)))
rm(DataT,indx,groups)

# Merge GPS data with multivariate abundance data
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Plantation<-GPS$Plantation
GPS<-SpatialPointsDataFrame(GPS[,2:3],GPS,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
GPS<-spTransform(GPS,CRS="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Proj_Cor<-as.data.frame(GPS@coords)

# Removing species with 3 or less occurrences
Data<-Data[,colSums(Data!=0)>3]
Low<-Low[,colSums(Low!=0)>3]
Intermediate<-Intermediate[,colSums(Intermediate!=0)>3]
High<-High[,colSums(High!=0)>3]
Data<-as.data.frame(cbind(Plantation,Proj_Cor,Data))
Low<-as.data.frame(cbind(Plantation,Proj_Cor,Low))
Intermediate<-as.data.frame(cbind(Plantation,Proj_Cor,Intermediate))
High<-as.data.frame(cbind(Plantation,Proj_Cor,High))
rm(GPS,Proj_Cor,Plantation)

# Transform multivariate abundance data to occurrence matrix
Data<-Data %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
Low<-Low %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
Intermediate<-Intermediate %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
High<-High %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name) %>%
  subset(Presence==1) %>%
  dplyr::select(-Presence)
write.csv(Data,"Products/Excel_Sheets_Products/Occurrence_Data.csv",row.names=FALSE)
write.csv(Low,"Products/Excel_Sheets_Products/Occurrence_Low_Data.csv",row.names=FALSE)
write.csv(Intermediate,"Products/Excel_Sheets_Products/Occurrence_Intermediate_Data.csv",row.names=FALSE)
write.csv(High,"Products/Excel_Sheets_Products/Occurrence_High_Data.csv",row.names=FALSE)
rm(Data,Low,Intermediate,High)

#Create design variable (corridor width)
Ranger_map<-raster("Products/Raster_Products/Map_ranger.tif")
Ranger_map[Ranger_map==2]<-NA
Ranger_map[Ranger_map==1]<-1
Ranger_map[Ranger_map==3]<-1
Ranger_map[Ranger_map==4]<-1
Ranger_map[Ranger_map==5]<-1
Ranger_map[Ranger_map==6]<-1
Ranger_map[Ranger_map==7]<-1
Estate_1<-readOGR("Shapefiles/Study_Estate1_Extent.shp")
Ranger_map_Estate1<-crop(Ranger_map,Estate_1)
Ranger_map_Estate1<-mask(Ranger_map_Estate1,Estate_1)
Ranger_map_Estate1<-raster::distance(Ranger_map_Estate1)
Ranger_map_Estate1<-crop(Ranger_map_Estate1,Estate_1)
Ranger_map_Estate1<-mask(Ranger_map_Estate1,Estate_1)
Estate_2<-readOGR("Shapefiles/Study_Estate2_Extent.shp")
Ranger_map_Estate2<-crop(Ranger_map,Estate_2)
Ranger_map_Estate2<-mask(Ranger_map_Estate2,Estate_2)
Ranger_map_Estate2<-raster::distance(Ranger_map_Estate2)
Ranger_map_Estate2<-crop(Ranger_map_Estate2,Estate_2)
Ranger_map_Estate2<-mask(Ranger_map_Estate2,Estate_2)
Estate_3<-readOGR("Shapefiles/Study_Estate3_Extent.shp")
Ranger_map_Estate3<-crop(Ranger_map,Estate_3)
Ranger_map_Estate3<-mask(Ranger_map_Estate3,Estate_3)
Ranger_map_Estate3<-raster::distance(Ranger_map_Estate3)
Ranger_map_Estate3<-crop(Ranger_map_Estate3,Estate_3)
Ranger_map_Estate3<-mask(Ranger_map_Estate3,Estate_3)
Estate_4<-readOGR("Shapefiles/Study_Estate4_Extent.shp")
Ranger_map_Estate4<-crop(Ranger_map,Estate_4)
Ranger_map_Estate4<-mask(Ranger_map_Estate4,Estate_4)
Ranger_map_Estate4<-raster::distance(Ranger_map_Estate4)
Ranger_map_Estate4<-crop(Ranger_map_Estate4,Estate_4)
Ranger_map_Estate4<-mask(Ranger_map_Estate4,Estate_4)
Cor_Width<-merge(Ranger_map_Estate1,Ranger_map_Estate2,Ranger_map_Estate3,Ranger_map_Estate4)
writeRaster(Cor_Width,filename="Products/Raster_Products/Cor_Width.tif",format="GTiff",overwrite=TRUE)
rm(Estate_1,Estate_2,Estate_3,Estate_4,Ranger_map_Estate1,Ranger_map_Estate2,Ranger_map_Estate3,Ranger_map_Estate4,Ranger_map,Cor_Width)

# Load environmental raster data
Env<-load_var(path="~/.../Products/Raster_Products/",
              files=c("Sent_NDVI.tif","Map_ranger.tif","Drain_Dist.tif","SUDEM_Aspect.tif","Cor_Width.tif","Max_NBR.tif"),
              format=".tif",Norm=TRUE,categorical=c("Map_ranger.tif"))

# Load occurrence data
Occ<-load_occ(path="~/.../Products/Excel_Sheets_Products",
              Env=Env,file="Occurrence_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")

# Perform stacked species distribution models
SSDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                      Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",tmp=TRUE,rep=10,metric="Kappa",
                      save=TRUE,name="SSDM",path="~/.../Products/",
                      cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))

# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/SSDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Drainage","Aspect","Dist Edge","NBR")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(variable.importance)

# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/SSDM/Stack/Tables/StackEval.csv",row.names = 1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(evaluation)

# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/SSDM/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#5x4
rm(correlation,Occ)

# Load occurrence data Low
Occ_Low<-load_occ(path="~/.../Products/Excel_Sheets_Products",
                  Env=Env,file="Occurrence_Low_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")

# Perform stacked species distribution models
SSDM_Low<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ_Low,Env=Env,cores=4,method="pSSDM",
                          Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",tmp=TRUE,rep=10,metric="Kappa",
                          save=TRUE,name="SSDM_Low",path="~/.../Products/",
                          cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))

# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/SSDM_Low/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Drainage","Aspect","Dist Edge","NBR")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(variable.importance)

# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/SSDM_Low/Stack/Tables/StackEval.csv",row.names = 1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(evaluation)

# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/SSDM_Low/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#5x4
rm(correlation,Occ_Low,SSDM_Low)

# Load occurrence data Intermediate
Occ_Intermediate<-load_occ(path="~/.../Products/Excel_Sheets_Products",
                           Env=Env,file="Occurrence_Intermediate_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")

# Perform stacked species distribution models
SSDM_Intermediate<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ_Intermediate,Env=Env,cores=4,method="pSSDM",
                                   Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",tmp=TRUE,rep=10,metric="Kappa",
                                   save=TRUE,name="SSDM_Intermediate",path="~/.../Products/",
                                   cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))

# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/SSDM_Intermediate/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Drainage","Aspect","Dist Edge","NBR")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(variable.importance)

# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/SSDM_Intermediate/Stack/Tables/StackEval.csv",row.names = 1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(evaluation)

# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/SSDM_Intermediate/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#5x4
rm(correlation,Occ_Intermediate,SSDM_Intermediate)

# Load occurrence data High
Occ_High<-load_occ(path="~/.../Products/Excel_Sheets_Products",
                   Env=Env,file="Occurrence_High_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")


# Perform stacked species distribution models
SSDM_High<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ_High,Env=Env,cores=4,method="pSSDM",
                           Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",tmp=TRUE,rep=10,metric="Kappa",
                           save=TRUE,name="SSDM_High",path="~/.../Products/",
                           cv="holdout",cv.param=c(0.7,1),ensemble.metric=c("AUC"),ensemble.thresh=c(0.7))

# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/SSDM_High/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Land Use","Drainage","Aspect","Dist Edge","NBR")
variable.importance<-dplyr::add_rownames(variable.importance,var="Variables")
ggplot(variable.importance) +
  geom_bar(aes(x=Variables,y=Mean),stat="identity",alpha=0.7) +
  ylab("Variable relative contribution (%)\n") +
  xlab("\nVariables") +
  geom_errorbar(aes(x=Variables,ymin=Mean-SD,ymax=Mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(variable.importance)

# Evaluate model
evaluation<-as.data.frame(t(read.csv("Products/SSDM_High/Stack/Tables/StackEval.csv",row.names = 1)))
row.names(evaluation)<-c("Rich Error","Prediction","Kappa","Specificity","Sensitivity","Jaccard")
evaluation<-dplyr::add_rownames(evaluation,var="Variables")
ggplot(evaluation) +
  geom_bar(aes(x=Variables,y=mean),stat="identity",alpha=0.7) +
  ylab("Evaluation metric contribution\n") +
  xlab("\nMetrics") +
  geom_errorbar(aes(x=Variables,ymin=mean-SD,ymax=mean+SD),alpha=0.9,width=0.2) +
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#7x5
rm(evaluation)

# Algorithm correlation
correlation<-as.data.frame(read.csv("Products/SSDM_High/Stack/Tables/AlgoCorr.csv"))
correlation<-pivot_longer(data=correlation, 
                          cols=-c(1), 
                          names_to="Variable", 
                          values_to="Correlation")
ggplot(data=correlation,mapping=aes(x=X,y=Variable,fill=Correlation)) +
  geom_tile(aes(fill=Correlation)) +
  geom_text(aes(label = round(Correlation, 3))) +
  scale_fill_gradient(low = "white", high = "red") +
  xlab(label="Variable")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
#5x4
rm(correlation,Occ_High,Env,SSDM_High)

############### Prepare data for Circuitscape ###############
# Load diversity map and select highest 95% of suitable pixel values
Diversity_map<-raster("Products/SSDM/Stack/Rasters/Diversity.tif")
Focal_Areas<-Diversity_map[[1]]>=quantile(Diversity_map,.95)
Focal_Areas[Focal_Areas==0]<-NA
Diversity_map_Low<-raster("Products/SSDM_Low/Stack/Rasters/Diversity.tif")
Focal_Areas_Low<-Diversity_map_Low[[1]]>=quantile(Diversity_map_Low,.95)
Focal_Areas_Low[Focal_Areas_Low==0]<-NA
Diversity_map_Intermediate<-raster("Products/SSDM_Intermediate/Stack/Rasters/Diversity.tif")
Focal_Areas_Intermediate<-Diversity_map_Intermediate[[1]]>=quantile(Diversity_map_Intermediate,.95)
Focal_Areas_Intermediate[Focal_Areas_Intermediate==0]<-NA
Diversity_map_High<-raster("Products/SSDM_High/Stack/Rasters/Diversity.tif")
Focal_Areas_High<-Diversity_map_High[[1]]>=quantile(Diversity_map_High,.95)
Focal_Areas_High[Focal_Areas_High==0]<-NA

# Convert to vector
Vector_Focal_Areas<-rasterToPolygons(Focal_Areas)
Vector_Focal_Areas$Obj_ID<-1:nrow(Vector_Focal_Areas)
Vector_Focal_Areas_Low<-rasterToPolygons(Focal_Areas_Low)
Vector_Focal_Areas_Low$Obj_ID<-1:nrow(Vector_Focal_Areas_Low)
Vector_Focal_Areas_Intermediate<-rasterToPolygons(Focal_Areas_Intermediate)
Vector_Focal_Areas_Intermediate$Obj_ID<-1:nrow(Vector_Focal_Areas_Intermediate)
Vector_Focal_Areas_High<-rasterToPolygons(Focal_Areas_High)
Vector_Focal_Areas_High$Obj_ID<-1:nrow(Vector_Focal_Areas_High)

# Dissolve vector
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas,dissolve=TRUE)
Vector_Focal_Areas_dis<-disaggregate(Vector_Focal_Areas_agg)
Vector_Focal_Areas_buf<-buffer(Vector_Focal_Areas_dis,width=0.001,dissolve=TRUE)
Vector_Focal_Areas<-disaggregate(Vector_Focal_Areas_buf)
rm(Vector_Focal_Areas_agg,Vector_Focal_Areas_dis,Vector_Focal_Areas_buf)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_Low,dissolve=TRUE)
Vector_Focal_Areas_dis<-disaggregate(Vector_Focal_Areas_agg)
Vector_Focal_Areas_buf<-buffer(Vector_Focal_Areas_dis,width=0.001,dissolve=TRUE)
Vector_Focal_Areas_Low<-disaggregate(Vector_Focal_Areas_buf)
rm(Vector_Focal_Areas_agg,Vector_Focal_Areas_dis,Vector_Focal_Areas_buf)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_Intermediate,dissolve=TRUE)
Vector_Focal_Areas_dis<-disaggregate(Vector_Focal_Areas_agg)
Vector_Focal_Areas_buf<-buffer(Vector_Focal_Areas_dis,width=0.001,dissolve=TRUE)
Vector_Focal_Areas_Intermediate<-disaggregate(Vector_Focal_Areas_buf)
rm(Vector_Focal_Areas_agg,Vector_Focal_Areas_dis,Vector_Focal_Areas_buf)
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas_High,dissolve=TRUE)
Vector_Focal_Areas_dis<-disaggregate(Vector_Focal_Areas_agg)
Vector_Focal_Areas_buf<-buffer(Vector_Focal_Areas_dis,width=0.001,dissolve=TRUE)
Vector_Focal_Areas_High<-disaggregate(Vector_Focal_Areas_buf)
rm(Vector_Focal_Areas_agg,Vector_Focal_Areas_dis,Vector_Focal_Areas_buf)

# Select only large focal areas
Vector_Focal_Areas_features<-as(Vector_Focal_Areas,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
Vector_Focal_Areas_features_100<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=10000) #_100_pixel region
Vector_Focal_Areas_features_50<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=5000) #_50_pixel region
Vector_Focal_Areas_features_Low<-as(Vector_Focal_Areas_Low,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features_Low$Obj_ID<-1:nrow(Vector_Focal_Areas_features_Low)
Vector_Focal_Areas_features_Low$Area_sqm<-area(Vector_Focal_Areas_features_Low)
Vector_Focal_Areas_features_Low<-st_as_sf(Vector_Focal_Areas_features_Low)
Vector_Focal_Areas_features_Low_100<-Vector_Focal_Areas_features_Low %>%
  filter(Area_sqm>=10000) #_100_pixel region
Vector_Focal_Areas_features_Low_50<-Vector_Focal_Areas_features_Low %>%
  filter(Area_sqm>=5000) #_50_pixel region
Vector_Focal_Areas_features_Intermediate<-as(Vector_Focal_Areas_Intermediate,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features_Intermediate$Obj_ID<-1:nrow(Vector_Focal_Areas_features_Intermediate)
Vector_Focal_Areas_features_Intermediate$Area_sqm<-area(Vector_Focal_Areas_features_Intermediate)
Vector_Focal_Areas_features_Intermediate<-st_as_sf(Vector_Focal_Areas_features_Intermediate)
Vector_Focal_Areas_features_Intermediate_100<-Vector_Focal_Areas_features_Intermediate %>%
  filter(Area_sqm>=10000) #_100_pixel region
Vector_Focal_Areas_features_Intermediate_50<-Vector_Focal_Areas_features_Intermediate %>%
  filter(Area_sqm>=5000) #_50_pixel region
Vector_Focal_Areas_features_High<-as(Vector_Focal_Areas_High,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features_High$Obj_ID<-1:nrow(Vector_Focal_Areas_features_High)
Vector_Focal_Areas_features_High$Area_sqm<-area(Vector_Focal_Areas_features_High)
Vector_Focal_Areas_features_High<-st_as_sf(Vector_Focal_Areas_features_High)
Vector_Focal_Areas_features_High_100<-Vector_Focal_Areas_features_High %>%
  filter(Area_sqm>=10000) #_100_pixel region
Vector_Focal_Areas_features_High_50<-Vector_Focal_Areas_features_High %>%
  filter(Area_sqm>=5000) #_50_pixel region

# Save to aci file for Circuitscape
write.dir<-"~/.../Circuitscape/"
rst_template<-raster(crs=projection(Diversity_map),ext=extent(Diversity_map),resolution=20)
Focal_Areas_100<-rasterize(Vector_Focal_Areas_features_100,rst_template)
Focal_Areas_50<-rasterize(Vector_Focal_Areas_features_50,rst_template)
Estate_1<-readOGR("Shapefiles/Study_Estate1_Extent.shp")
Focal_Areas_Estate_1<-crop(Focal_Areas_100,Estate_1)
Focal_Areas_Estate_1<-mask(Focal_Areas_Estate_1,Estate_1)
writeRaster(Focal_Areas_Estate_1,filename=paste0(write.dir,"Focal_Areas_Estate_1.asc"),overwrite=TRUE)
Estate_2<-readOGR("Shapefiles/Study_Estate2_Extent.shp")
Focal_Areas_Estate_2<-crop(Focal_Areas_100,Estate_2)
Focal_Areas_Estate_2<-mask(Focal_Areas_Estate_2,Estate_2)
writeRaster(Focal_Areas_Estate_2,filename=paste0(write.dir,"Focal_Areas_Estate_2.asc"),overwrite=TRUE)
Estate_3<-readOGR("Shapefiles/Study_Estate3_Extent.shp")
Focal_Areas_Estate_3<-crop(Focal_Areas_100,Estate_3)
Focal_Areas_Estate_3<-mask(Focal_Areas_Estate_3,Estate_3)
writeRaster(Focal_Areas_Estate_3,filename=paste0(write.dir,"Focal_Areas_Estate_3.asc"),overwrite=TRUE)
Estate_4<-readOGR("Shapefiles/Study_Estate4_Extent.shp")
Focal_Areas_Estate_4<-crop(Focal_Areas_50,Estate_4)
Focal_Areas_Estate_4<-mask(Focal_Areas_Estate_4,Estate_4)
writeRaster(Focal_Areas_Estate_4,filename=paste0(write.dir,"Focal_Areas_Estate_4.asc"),overwrite=TRUE)
Focal_Areas_Low_100<-rasterize(Vector_Focal_Areas_features_Low_100,rst_template)
Focal_Areas_Low_50<-rasterize(Vector_Focal_Areas_features_Low_50,rst_template)
Focal_Areas_Estate_1_Low<-crop(Focal_Areas_Low_100,Estate_1)
Focal_Areas_Estate_1_Low<-mask(Focal_Areas_Estate_1_Low,Estate_1)
writeRaster(Focal_Areas_Estate_1_Low,filename=paste0(write.dir,"Focal_Areas_Estate_1_Low.asc"),overwrite=TRUE)
Focal_Areas_Estate_2_Low<-crop(Focal_Areas_Low_100,Estate_2)
Focal_Areas_Estate_2_Low<-mask(Focal_Areas_Estate_2_Low,Estate_2)
writeRaster(Focal_Areas_Estate_2_Low,filename=paste0(write.dir,"Focal_Areas_Estate_2_Low.asc"),overwrite=TRUE)
Focal_Areas_Estate_3_Low<-crop(Focal_Areas_Low_100,Estate_3)
Focal_Areas_Estate_3_Low<-mask(Focal_Areas_Estate_3_Low,Estate_3)
writeRaster(Focal_Areas_Estate_3_Low,filename=paste0(write.dir,"Focal_Areas_Estate_3_Low.asc"),overwrite=TRUE)
Focal_Areas_Estate_4_Low<-crop(Focal_Areas_Low_50,Estate_4)
Focal_Areas_Estate_4_Low<-mask(Focal_Areas_Estate_4_Low,Estate_4)
writeRaster(Focal_Areas_Estate_4_Low,filename=paste0(write.dir,"Focal_Areas_Estate_4_Low.asc"),overwrite=TRUE)
Focal_Areas_Intermediate_100<-rasterize(Vector_Focal_Areas_features_Intermediate_100,rst_template)
Focal_Areas_Intermediate_50<-rasterize(Vector_Focal_Areas_features_Intermediate_50,rst_template)
Focal_Areas_Estate_1_Intermediate<-crop(Focal_Areas_Intermediate_100,Estate_1)
Focal_Areas_Estate_1_Intermediate<-mask(Focal_Areas_Estate_1_Intermediate,Estate_1)
writeRaster(Focal_Areas_Estate_1_Intermediate,filename=paste0(write.dir,"Focal_Areas_Estate_1_Intermediate.asc"),overwrite=TRUE)
Focal_Areas_Estate_2_Intermediate<-crop(Focal_Areas_Intermediate_100,Estate_2)
Focal_Areas_Estate_2_Intermediate<-mask(Focal_Areas_Estate_2_Intermediate,Estate_2)
writeRaster(Focal_Areas_Estate_2_Intermediate,filename=paste0(write.dir,"Focal_Areas_Estate_2_Intermediate.asc"),overwrite=TRUE)
Focal_Areas_Estate_3_Intermediate<-crop(Focal_Areas_Intermediate_100,Estate_3)
Focal_Areas_Estate_3_Intermediate<-mask(Focal_Areas_Estate_3_Intermediate,Estate_3)
writeRaster(Focal_Areas_Estate_3_Intermediate,filename=paste0(write.dir,"Focal_Areas_Estate_3_Intermediate.asc"),overwrite=TRUE)
Focal_Areas_Estate_4_Intermediate<-crop(Focal_Areas_Intermediate_50,Estate_4)
Focal_Areas_Estate_4_Intermediate<-mask(Focal_Areas_Estate_4_Intermediate,Estate_4)
writeRaster(Focal_Areas_Estate_4_Intermediate,filename=paste0(write.dir,"Focal_Areas_Estate_4_Intermediate.asc"),overwrite=TRUE)
Focal_Areas_High_100<-rasterize(Vector_Focal_Areas_features_High_100,rst_template)
Focal_Areas_High_50<-rasterize(Vector_Focal_Areas_features_High_50,rst_template)
Focal_Areas_Estate_1_High<-crop(Focal_Areas_High_100,Estate_1)
Focal_Areas_Estate_1_High<-mask(Focal_Areas_Estate_1_High,Estate_1)
writeRaster(Focal_Areas_Estate_1_High,filename=paste0(write.dir,"Focal_Areas_Estate_1_High.asc"),overwrite=TRUE)
Focal_Areas_Estate_2_High<-crop(Focal_Areas_High_100,Estate_2)
Focal_Areas_Estate_2_High<-mask(Focal_Areas_Estate_2_High,Estate_2)
writeRaster(Focal_Areas_Estate_2_High,filename=paste0(write.dir,"Focal_Areas_Estate_2_High.asc"),overwrite=TRUE)
Focal_Areas_Estate_3_High<-crop(Focal_Areas_High_100,Estate_3)
Focal_Areas_Estate_3_High<-mask(Focal_Areas_Estate_3_High,Estate_3)
writeRaster(Focal_Areas_Estate_3_High,filename=paste0(write.dir,"Focal_Areas_Estate_3_High.asc"),overwrite=TRUE)
Focal_Areas_Estate_4_High<-crop(Focal_Areas_High_50,Estate_4)
Focal_Areas_Estate_4_High<-mask(Focal_Areas_Estate_4_High,Estate_4)
writeRaster(Focal_Areas_Estate_4_High,filename=paste0(write.dir,"Focal_Areas_Estate_4_High.asc"),overwrite=TRUE)
rm(Focal_Areas_Estate_4,Focal_Areas_Estate_4_High,Focal_Areas_Estate_4_Intermediate,Focal_Areas_Estate_4_Low,Focal_Areas_High,
   Focal_Areas_High_100,Focal_Areas_50,Focal_Areas_Intermediate,Focal_Areas_Intermediate_100,Focal_Areas_Intermediate_50,
   Focal_Areas_Low,Focal_Areas_Low_100,Focal_Areas_Low_50,Vector_Focal_Areas,Vector_Focal_Areas_features,Focal_Areas,
   Focal_Areas_100,Focal_Areas_Estate_1,Focal_Areas_Estate_1_High,Focal_Areas_Estate_1_Intermediate,Focal_Areas_Estate_1_Low,
   Focal_Areas_Estate_2,Focal_Areas_Estate_2_High,Focal_Areas_Estate_2_Intermediate,Focal_Areas_Estate_2_Low,Focal_Areas_Estate_3,
   Focal_Areas_Estate_3_High,Focal_Areas_Estate_3_Intermediate,Focal_Areas_Estate_3_Low,Focal_Areas_High_50,Vector_Focal_Areas_features_100,
   Vector_Focal_Areas_features_50,Vector_Focal_Areas_features_High,Vector_Focal_Areas_features_High_100,Vector_Focal_Areas_features_High_50,
   Vector_Focal_Areas_features_Intermediate,Vector_Focal_Areas_features_Intermediate_100,Vector_Focal_Areas_features_Intermediate_50,
   Vector_Focal_Areas_features_Low,Vector_Focal_Areas_features_Low_100,Vector_Focal_Areas_features_Low_50,Vector_Focal_Areas,
   Vector_Focal_Areas_High,Vector_Focal_Areas_Intermediate,Vector_Focal_Areas_Low)

# Save Diversity map to aci file for Circuitscape
Diversity_map<-reclassify(Diversity_map,cbind(-3.4e+38, -9999))
Diversity_map<-resample(Diversity_map,rst_template)
Diversity_map_Estate_1<-crop(Diversity_map,Estate_1)
Diversity_map_Estate_1<-mask(Diversity_map_Estate_1,Estate_1)
writeRaster(Diversity_map_Estate_1,filename =paste0(write.dir,"Diversity_Estate_1.asc"),overwrite=TRUE)
Diversity_map_Estate_2<-crop(Diversity_map,Estate_2)
Diversity_map_Estate_2<-mask(Diversity_map_Estate_2,Estate_2)
writeRaster(Diversity_map_Estate_2,filename =paste0(write.dir,"Diversity_Estate_2.asc"),overwrite=TRUE)
Diversity_map_Estate_3<-crop(Diversity_map,Estate_3)
Diversity_map_Estate_3<-mask(Diversity_map_Estate_3,Estate_3)
writeRaster(Diversity_map_Estate_3,filename =paste0(write.dir,"Diversity_Estate_3.asc"),overwrite=TRUE)
Diversity_map_Estate_4<-crop(Diversity_map,Estate_4)
Diversity_map_Estate_4<-mask(Diversity_map_Estate_4,Estate_4)
writeRaster(Diversity_map_Estate_4,filename =paste0(write.dir,"Diversity_Estate_4.asc"),overwrite=TRUE)
Diversity_map_Low<-reclassify(Diversity_map_Low,cbind(-3.4e+38, -9999))
Diversity_map_Low<-resample(Diversity_map_Low,rst_template)
Diversity_map_Estate_1_Low<-crop(Diversity_map_Low,Estate_1)
Diversity_map_Estate_1_Low<-mask(Diversity_map_Estate_1_Low,Estate_1)
writeRaster(Diversity_map_Estate_1_Low,filename =paste0(write.dir,"Diversity_Estate_1_Low.asc"),overwrite=TRUE)
Diversity_map_Estate_2_Low<-crop(Diversity_map_Low,Estate_2)
Diversity_map_Estate_2_Low<-mask(Diversity_map_Estate_2_Low,Estate_2)
writeRaster(Diversity_map_Estate_2_Low,filename =paste0(write.dir,"Diversity_Estate_2_Low.asc"),overwrite=TRUE)
Diversity_map_Estate_3_Low<-crop(Diversity_map_Low,Estate_3)
Diversity_map_Estate_3_Low<-mask(Diversity_map_Estate_3_Low,Estate_3)
writeRaster(Diversity_map_Estate_3_Low,filename =paste0(write.dir,"Diversity_Estate_3_Low.asc"),overwrite=TRUE)
Diversity_map_Estate_4_Low<-crop(Diversity_map_Low,Estate_4)
Diversity_map_Estate_4_Low<-mask(Diversity_map_Estate_4_Low,Estate_4)
writeRaster(Diversity_map_Estate_4_Low,filename =paste0(write.dir,"Diversity_Estate_4_Low.asc"),overwrite=TRUE)
Diversity_map_Intermediate<-reclassify(Diversity_map_Intermediate,cbind(-3.4e+38, -9999))
Diversity_map_Intermediate<-resample(Diversity_map_Intermediate,rst_template)
Diversity_map_Estate_1_Intermediate<-crop(Diversity_map_Intermediate,Estate_1)
Diversity_map_Estate_1_Intermediate<-mask(Diversity_map_Estate_1_Intermediate,Estate_1)
writeRaster(Diversity_map_Estate_1_Intermediate,filename =paste0(write.dir,"Diversity_Estate_1_Intermediate.asc"),overwrite=TRUE)
Diversity_map_Estate_2_Intermediate<-crop(Diversity_map_Intermediate,Estate_2)
Diversity_map_Estate_2_Intermediate<-mask(Diversity_map_Estate_2_Intermediate,Estate_2)
writeRaster(Diversity_map_Estate_2_Intermediate,filename =paste0(write.dir,"Diversity_Estate_2_Intermediate.asc"),overwrite=TRUE)
Diversity_map_Estate_3_Intermediate<-crop(Diversity_map_Intermediate,Estate_3)
Diversity_map_Estate_3_Intermediate<-mask(Diversity_map_Estate_3_Intermediate,Estate_3)
writeRaster(Diversity_map_Estate_3_Intermediate,filename =paste0(write.dir,"Diversity_Estate_3_Intermediate.asc"),overwrite=TRUE)
Diversity_map_Estate_4_Intermediate<-crop(Diversity_map_Intermediate,Estate_4)
Diversity_map_Estate_4_Intermediate<-mask(Diversity_map_Estate_4_Intermediate,Estate_4)
writeRaster(Diversity_map_Estate_4_Intermediate,filename =paste0(write.dir,"Diversity_Estate_4_Intermediate.asc"),overwrite=TRUE)
Diversity_map_High<-reclassify(Diversity_map_High,cbind(-3.4e+38, -9999))
Diversity_map_High<-resample(Diversity_map_High,rst_template)
Diversity_map_Estate_1_High<-crop(Diversity_map_High,Estate_1)
Diversity_map_Estate_1_High<-mask(Diversity_map_Estate_1_High,Estate_1)
writeRaster(Diversity_map_Estate_1_High,filename =paste0(write.dir,"Diversity_Estate_1_High.asc"),overwrite=TRUE)
Diversity_map_Estate_2_High<-crop(Diversity_map_High,Estate_2)
Diversity_map_Estate_2_High<-mask(Diversity_map_Estate_2_High,Estate_2)
writeRaster(Diversity_map_Estate_2_High,filename =paste0(write.dir,"Diversity_Estate_2_High.asc"),overwrite=TRUE)
Diversity_map_Estate_3_High<-crop(Diversity_map_High,Estate_3)
Diversity_map_Estate_3_High<-mask(Diversity_map_Estate_3_High,Estate_3)
writeRaster(Diversity_map_Estate_3_High,filename =paste0(write.dir,"Diversity_Estate_3_High.asc"),overwrite=TRUE)
Diversity_map_Estate_4_High<-crop(Diversity_map_High,Estate_4)
Diversity_map_Estate_4_High<-mask(Diversity_map_Estate_4_High,Estate_4)
writeRaster(Diversity_map_Estate_4_High,filename =paste0(write.dir,"Diversity_Estate_4_High.asc"),overwrite=TRUE)
rm(Diversity_map,Diversity_map_Estate_1,Diversity_map_Estate_1_High,Diversity_map_Estate_1_Intermediate,Diversity_map_Estate_1_Low,
   Diversity_map_Estate_2,Diversity_map_Estate_2_High,Diversity_map_Estate_2_Intermediate,Diversity_map_Estate_2_Low,Diversity_map_Estate_3,
   Diversity_map_Estate_3_High,Diversity_map_Estate_3_Intermediate,Diversity_map_Estate_3_Low,Diversity_map_Estate_4,
   Diversity_map_Estate_4_High,Diversity_map_Estate_4_Intermediate,Diversity_map_Estate_4_Low,Diversity_map_High,
   Diversity_map_Intermediate,Diversity_map_Low)

# Load and transform NDVI
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI.tif")
Sent_NDVI<-resample(Sent_NDVI,rst_template)
Plot.trans(PARM =c(3,100),Resistance=Sent_NDVI,transformation="Reverse Ricker")
NDVI_trans<-Resistance.tran("Reverse Ricker",shape=3,max=100,scale=NULL,r=Sent_NDVI)
NDVI_trans_Estate_1<-crop(NDVI_trans,Estate_1)
NDVI_trans_Estate_1<-mask(NDVI_trans_Estate_1,Estate_1)
writeRaster(NDVI_trans_Estate_1,filename =paste0(write.dir,"NDVI_trans_Estate_1.asc"),overwrite=TRUE)
NDVI_trans_Estate_2<-crop(NDVI_trans,Estate_2)
NDVI_trans_Estate_2<-mask(NDVI_trans_Estate_2,Estate_2)
writeRaster(NDVI_trans_Estate_2,filename =paste0(write.dir,"NDVI_trans_Estate_2.asc"),overwrite=TRUE)
NDVI_trans_Estate_3<-crop(NDVI_trans,Estate_3)
NDVI_trans_Estate_3<-mask(NDVI_trans_Estate_3,Estate_3)
writeRaster(NDVI_trans_Estate_3,filename =paste0(write.dir,"NDVI_trans_Estate_3.asc"),overwrite=TRUE)
NDVI_trans_Estate_4<-crop(NDVI_trans,Estate_4)
NDVI_trans_Estate_4<-mask(NDVI_trans_Estate_4,Estate_4)
writeRaster(NDVI_trans_Estate_4,filename =paste0(write.dir,"NDVI_trans_Estate_4.asc"),overwrite=TRUE)

# Load and reclassify classification
Map_ranger<-raster("Products/Raster_Products/Map_ranger.tif")
Map_ranger<-resample(Map_ranger,rst_template)
'''
ID  Land use      Conductivity value (0-100)
1 = Agriculture = 40
2 = Grassland   = 100
3 = Ground      = 70
4 = Plantation  = 10
5 = Shrubland   = 30
6 = Water       = 0
7 = Woodland    = 20
'''
Map_ranger[Map_ranger==1]<-40
Map_ranger[Map_ranger==2]<-100
Map_ranger[Map_ranger==3]<-70
Map_ranger[Map_ranger==4]<-10
Map_ranger[Map_ranger==5]<-30
Map_ranger[Map_ranger==6]<-0
Map_ranger[Map_ranger==7]<-20
Map_ranger_Estate_1<-crop(Map_ranger,Estate_1)
Map_ranger_Estate_1<-mask(Map_ranger_Estate_1,Estate_1)
writeRaster(Map_ranger_Estate_1,filename =paste0(write.dir,"Map_ranger_Estate_1.asc"),overwrite=TRUE)
Map_ranger_Estate_2<-crop(Map_ranger,Estate_2)
Map_ranger_Estate_2<-mask(Map_ranger_Estate_2,Estate_2)
writeRaster(Map_ranger_Estate_2,filename =paste0(write.dir,"Map_ranger_Estate_2.asc"),overwrite=TRUE)
Map_ranger_Estate_3<-crop(Map_ranger,Estate_3)
Map_ranger_Estate_3<-mask(Map_ranger_Estate_3,Estate_3)
writeRaster(Map_ranger_Estate_3,filename =paste0(write.dir,"Map_ranger_Estate_3.asc"),overwrite=TRUE)
Map_ranger_Estate_4<-crop(Map_ranger,Estate_4)
Map_ranger_Estate_4<-mask(Map_ranger_Estate_4,Estate_4)
writeRaster(Map_ranger_Estate_4,filename =paste0(write.dir,"Map_ranger_Estate_4.asc"),overwrite=TRUE)
rm(Estate_1,Estate_2,Estate_3,Estate_4,Map_ranger,Map_ranger_Estate_1,Map_ranger_Estate_2,Map_ranger_Estate_3,Map_ranger_Estate_4,
   NDVI_trans,NDVI_trans_Estate_1,NDVI_trans_Estate_2,NDVI_trans_Estate_3,NDVI_trans_Estate_4,rst_template,Sent_NDVI,write.dir)

############### Calculate proportions of movement habitat over land use types ###############
# Load land use data
Map_ranger<-raster("Products/Raster_Products/Map_ranger.tif")

# Load, merge and average a connectivity maps
Cur_Map_Estate1_SSDM<-raster("Circuitscape/Estate1/Estate1_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_SSDM<-raster("Circuitscape/Estate2/Estate2_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_SSDM<-raster("Circuitscape/Estate3/Estate3_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_SSDM<-raster("Circuitscape/Estate4/Estate4_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_SSDM<-merge(Cur_Map_Estate1_SSDM,Cur_Map_Estate2_SSDM,Cur_Map_Estate3_SSDM,Cur_Map_Estate4_SSDM)
Cur_Map_SSDM<-resample(Cur_Map_SSDM,Map_ranger)
rm(Cur_Map_Estate1_SSDM,Cur_Map_Estate2_SSDM,Cur_Map_Estate3_SSDM,Cur_Map_Estate4_SSDM)
Cur_Map_Estate1_NDVI<-raster("Circuitscape/Estate1_NDVI/Estate1_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_NDVI<-raster("Circuitscape/Estate2_NDVI/Estate2_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_NDVI<-raster("Circuitscape/Estate3_NDVI/Estate3_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_NDVI<-raster("Circuitscape/Estate4_NDVI/Estate4_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_NDVI<-merge(Cur_Map_Estate1_NDVI,Cur_Map_Estate2_NDVI,Cur_Map_Estate3_NDVI,Cur_Map_Estate4_NDVI)
Cur_Map_NDVI<-resample(Cur_Map_NDVI,Map_ranger)
rm(Cur_Map_Estate1_NDVI,Cur_Map_Estate2_NDVI,Cur_Map_Estate3_NDVI,Cur_Map_Estate4_NDVI)
Cur_Map_Estate1_Ranger<-raster("Circuitscape/Estate1_Ranger/Estate1_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_Ranger<-raster("Circuitscape/Estate2_Ranger/Estate2_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_Ranger<-raster("Circuitscape/Estate3_Ranger/Estate3_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_Ranger<-raster("Circuitscape/Estate4_Ranger/Estate4_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Ranger<-merge(Cur_Map_Estate1_Ranger,Cur_Map_Estate2_Ranger,Cur_Map_Estate3_Ranger,Cur_Map_Estate4_Ranger)
Cur_Map_Ranger<-resample(Cur_Map_Ranger,Map_ranger)
rm(Cur_Map_Estate1_Ranger,Cur_Map_Estate2_Ranger,Cur_Map_Estate3_Ranger,Cur_Map_Estate4_Ranger)
Cur_Map<-stack(Cur_Map_SSDM,Cur_Map_NDVI,Cur_Map_Ranger)
ncores<-4
beginCluster(ncores)
Cur_Map_Mean<-clusterR(Cur_Map,mean,args=list(na.rm=TRUE))
endCluster()
Cur_Map_Mean[Cur_Map_Mean<0]<-NA
writeRaster(Cur_Map_Mean,filename="Products/Raster_Products/Cur_Map.tif",format="GTiff",overwrite=TRUE)
rm(Cur_Map,Cur_Map_SSDM,Cur_Map_NDVI,Cur_Map_Ranger,ncores)
Cur_Map_Estate1_SSDM_Low<-raster("Circuitscape/Estate1_Low/Estate1_Low_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_SSDM_Low<-raster("Circuitscape/Estate2_Low/Estate2_Low_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_SSDM_Low<-raster("Circuitscape/Estate3_Low/Estate3_Low_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_SSDM_Low<-raster("Circuitscape/Estate4_Low/Estate4_Low_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_SSDM_Low<-merge(Cur_Map_Estate1_SSDM_Low,Cur_Map_Estate2_SSDM_Low,Cur_Map_Estate3_SSDM_Low,Cur_Map_Estate4_SSDM_Low)
Cur_Map_SSDM_Low<-resample(Cur_Map_SSDM_Low,Map_ranger)
rm(Cur_Map_Estate1_SSDM_Low,Cur_Map_Estate2_SSDM_Low,Cur_Map_Estate3_SSDM_Low,Cur_Map_Estate4_SSDM_Low)
Cur_Map_Estate1_NDVI_Low<-raster("Circuitscape/Estate1_NDVI_Low/Estate1_Low_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_NDVI_Low<-raster("Circuitscape/Estate2_NDVI_Low/Estate2_Low_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_NDVI_Low<-raster("Circuitscape/Estate3_NDVI_Low/Estate3_Low_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_NDVI_Low<-raster("Circuitscape/Estate4_NDVI_Low/Estate4_Low_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_NDVI_Low<-merge(Cur_Map_Estate1_NDVI_Low,Cur_Map_Estate2_NDVI_Low,Cur_Map_Estate3_NDVI_Low,Cur_Map_Estate4_NDVI_Low)
Cur_Map_NDVI_Low<-resample(Cur_Map_NDVI_Low,Map_ranger)
rm(Cur_Map_Estate1_NDVI_Low,Cur_Map_Estate2_NDVI_Low,Cur_Map_Estate3_NDVI_Low,Cur_Map_Estate4_NDVI_Low)
Cur_Map_Estate1_Ranger_Low<-raster("Circuitscape/Estate1_Ranger_Low/Estate1_Low_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_Ranger_Low<-raster("Circuitscape/Estate2_Ranger_Low/Estate2_Low_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_Ranger_Low<-raster("Circuitscape/Estate3_Ranger_Low/Estate3_Low_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_Ranger_Low<-raster("Circuitscape/Estate4_Ranger_Low/Estate4_Low_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Ranger_Low<-merge(Cur_Map_Estate1_Ranger_Low,Cur_Map_Estate2_Ranger_Low,Cur_Map_Estate3_Ranger_Low,Cur_Map_Estate4_Ranger_Low)
Cur_Map_Ranger_Low<-resample(Cur_Map_Ranger_Low,Map_ranger)
rm(Cur_Map_Estate1_Ranger_Low,Cur_Map_Estate2_Ranger_Low,Cur_Map_Estate3_Ranger_Low,Cur_Map_Estate4_Ranger_Low)
Cur_Map_Low<-stack(Cur_Map_SSDM_Low,Cur_Map_NDVI_Low,Cur_Map_Ranger_Low)
ncores<-4
beginCluster(ncores)
Cur_Map_Mean_Low<-clusterR(Cur_Map_Low,mean,args=list(na.rm=TRUE))
endCluster()
Cur_Map_Mean_Low[Cur_Map_Mean_Low<0]<-NA
writeRaster(Cur_Map_Mean_Low,filename="Products/Raster_Products/Cur_Map_Low.tif",format="GTiff",overwrite=TRUE)
rm(Cur_Map_Low,Cur_Map_SSDM_Low,Cur_Map_NDVI_Low,Cur_Map_Ranger_Low,ncores)
Cur_Map_Estate1_SSDM_Intermediate<-raster("Circuitscape/Estate1_Intermediate/Estate1_Intermediate_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_SSDM_Intermediate<-raster("Circuitscape/Estate2_Intermediate/Estate2_Intermediate_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_SSDM_Intermediate<-raster("Circuitscape/Estate3_Intermediate/Estate3_Intermediate_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_SSDM_Intermediate<-raster("Circuitscape/Estate4_Intermediate/Estate4_Intermediate_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_SSDM_Intermediate<-merge(Cur_Map_Estate1_SSDM_Intermediate,Cur_Map_Estate2_SSDM_Intermediate,Cur_Map_Estate3_SSDM_Intermediate,Cur_Map_Estate4_SSDM_Intermediate)
Cur_Map_SSDM_Intermediate<-resample(Cur_Map_SSDM_Intermediate,Map_ranger)
rm(Cur_Map_Estate1_SSDM_Intermediate,Cur_Map_Estate2_SSDM_Intermediate,Cur_Map_Estate3_SSDM_Intermediate,Cur_Map_Estate4_SSDM_Intermediate)
Cur_Map_Estate1_NDVI_Intermediate<-raster("Circuitscape/Estate1_NDVI_Intermediate/Estate1_Intermediate_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_NDVI_Intermediate<-raster("Circuitscape/Estate2_NDVI_Intermediate/Estate2_Intermediate_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_NDVI_Intermediate<-raster("Circuitscape/Estate3_NDVI_Intermediate/Estate3_Intermediate_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_NDVI_Intermediate<-raster("Circuitscape/Estate4_NDVI_Intermediate/Estate4_Intermediate_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_NDVI_Intermediate<-merge(Cur_Map_Estate1_NDVI_Intermediate,Cur_Map_Estate2_NDVI_Intermediate,Cur_Map_Estate3_NDVI_Intermediate,Cur_Map_Estate4_NDVI_Intermediate)
Cur_Map_NDVI_Intermediate<-resample(Cur_Map_NDVI_Intermediate,Map_ranger)
rm(Cur_Map_Estate1_NDVI_Intermediate,Cur_Map_Estate2_NDVI_Intermediate,Cur_Map_Estate3_NDVI_Intermediate,Cur_Map_Estate4_NDVI_Intermediate)
Cur_Map_Estate1_Ranger_Intermediate<-raster("Circuitscape/Estate1_Ranger_Intermediate/Estate1_Intermediate_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_Ranger_Intermediate<-raster("Circuitscape/Estate2_Ranger_Intermediate/Estate2_Intermediate_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_Ranger_Intermediate<-raster("Circuitscape/Estate3_Ranger_Intermediate/Estate3_Intermediate_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_Ranger_Intermediate<-raster("Circuitscape/Estate4_Ranger_Intermediate/Estate4_Intermediate_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Ranger_Intermediate<-merge(Cur_Map_Estate1_Ranger_Intermediate,Cur_Map_Estate2_Ranger_Intermediate,Cur_Map_Estate3_Ranger_Intermediate,Cur_Map_Estate4_Ranger_Intermediate)
Cur_Map_Ranger_Intermediate<-resample(Cur_Map_Ranger_Intermediate,Map_ranger)
rm(Cur_Map_Estate1_Ranger_Intermediate,Cur_Map_Estate2_Ranger_Intermediate,Cur_Map_Estate3_Ranger_Intermediate,Cur_Map_Estate4_Ranger_Intermediate)
Cur_Map_Intermediate<-stack(Cur_Map_SSDM_Intermediate,Cur_Map_NDVI_Intermediate,Cur_Map_Ranger_Intermediate)
ncores<-4
beginCluster(ncores)
Cur_Map_Mean_Intermediate<-clusterR(Cur_Map_Intermediate,mean,args=list(na.rm=TRUE))
endCluster()
Cur_Map_Mean_Intermediate[Cur_Map_Mean_Intermediate<0]<-NA
writeRaster(Cur_Map_Mean_Intermediate,filename="Products/Raster_Products/Cur_Map_Intermediate.tif",format="GTiff",overwrite=TRUE)
rm(Cur_Map_Intermediate,Cur_Map_SSDM_Intermediate,Cur_Map_NDVI_Intermediate,Cur_Map_Ranger_Intermediate,ncores)
Cur_Map_Estate1_SSDM_High<-raster("Circuitscape/Estate1_High/Estate1_High_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_SSDM_High<-raster("Circuitscape/Estate2_High/Estate2_High_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_SSDM_High<-raster("Circuitscape/Estate3_High/Estate3_High_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_SSDM_High<-raster("Circuitscape/Estate4_High/Estate4_High_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_SSDM_High<-merge(Cur_Map_Estate1_SSDM_High,Cur_Map_Estate2_SSDM_High,Cur_Map_Estate3_SSDM_High,Cur_Map_Estate4_SSDM_High)
Cur_Map_SSDM_High<-resample(Cur_Map_SSDM_High,Map_ranger)
rm(Cur_Map_Estate1_SSDM_High,Cur_Map_Estate2_SSDM_High,Cur_Map_Estate3_SSDM_High,Cur_Map_Estate4_SSDM_High)
Cur_Map_Estate1_NDVI_High<-raster("Circuitscape/Estate1_NDVI_High/Estate1_High_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_NDVI_High<-raster("Circuitscape/Estate2_NDVI_High/Estate2_High_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_NDVI_High<-raster("Circuitscape/Estate3_NDVI_High/Estate3_High_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_NDVI_High<-raster("Circuitscape/Estate4_NDVI_High/Estate4_High_NDVI_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_NDVI_High<-merge(Cur_Map_Estate1_NDVI_High,Cur_Map_Estate2_NDVI_High,Cur_Map_Estate3_NDVI_High,Cur_Map_Estate4_NDVI_High)
Cur_Map_NDVI_High<-resample(Cur_Map_NDVI_High,Map_ranger)
rm(Cur_Map_Estate1_NDVI_High,Cur_Map_Estate2_NDVI_High,Cur_Map_Estate3_NDVI_High,Cur_Map_Estate4_NDVI_High)
Cur_Map_Estate1_Ranger_High<-raster("Circuitscape/Estate1_Ranger_High/Estate1_High_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate2_Ranger_High<-raster("Circuitscape/Estate2_Ranger_High/Estate2_High_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate3_Ranger_High<-raster("Circuitscape/Estate3_Ranger_High/Estate3_High_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Estate4_Ranger_High<-raster("Circuitscape/Estate4_Ranger_High/Estate4_High_Ranger_Pairwise_Cor4_Coductive_cg_amg_cum_curmap.asc")
Cur_Map_Ranger_High<-merge(Cur_Map_Estate1_Ranger_High,Cur_Map_Estate2_Ranger_High,Cur_Map_Estate3_Ranger_High,Cur_Map_Estate4_Ranger_High)
Cur_Map_Ranger_High<-resample(Cur_Map_Ranger_High,Map_ranger)
rm(Cur_Map_Estate1_Ranger_High,Cur_Map_Estate2_Ranger_High,Cur_Map_Estate3_Ranger_High,Cur_Map_Estate4_Ranger_High)
Cur_Map_High<-stack(Cur_Map_SSDM_High,Cur_Map_NDVI_High,Cur_Map_Ranger_High)
ncores<-4
beginCluster(ncores)
Cur_Map_Mean_High<-clusterR(Cur_Map_High,mean,args=list(na.rm=TRUE))
endCluster()
Cur_Map_Mean_High[Cur_Map_Mean_High<0]<-NA
writeRaster(Cur_Map_Mean_High,filename="Products/Raster_Products/Cur_Map_High.tif",format="GTiff",overwrite=TRUE)
rm(Cur_Map_High,Cur_Map_SSDM_High,Cur_Map_NDVI_High,Cur_Map_Ranger_High,ncores)

# Select top 85% of suitable pixel values
Cur_Map_Mean<-Cur_Map_Mean[[1]]>=quantile(Cur_Map_Mean,.85)
Cur_Map_Mean[Cur_Map_Mean==0]<-NA
Cur_Map_Mean_Low<-Cur_Map_Mean_Low[[1]]>=quantile(Cur_Map_Mean_Low,.85)
Cur_Map_Mean_Low[Cur_Map_Mean_Low==0]<-NA
Cur_Map_Mean_Intermediate<-Cur_Map_Mean_Intermediate[[1]]>=quantile(Cur_Map_Mean_Intermediate,.85)
Cur_Map_Mean_Intermediate[Cur_Map_Mean_Intermediate==0]<-NA
Cur_Map_Mean_High<-Cur_Map_Mean_High[[1]]>=quantile(Cur_Map_Mean_High,.85)
Cur_Map_Mean_High[Cur_Map_Mean_High==0]<-NA

# Clip land use to focal areas
Map_ranger<-crop(Map_ranger,Cur_Map_Mean)
Map_ranger<-mask(Map_ranger,Cur_Map_Mean)
Map_ranger_Low<-crop(Map_ranger,Cur_Map_Mean_Low)
Map_ranger_Low<-mask(Map_ranger_Low,Cur_Map_Mean_Low)
Map_ranger_Intermediate<-crop(Map_ranger,Cur_Map_Mean_Intermediate)
Map_ranger_Intermediate<-mask(Map_ranger_Intermediate,Cur_Map_Mean_Intermediate)
Map_ranger_High<-crop(Map_ranger,Cur_Map_Mean_High)
Map_ranger_High<-mask(Map_ranger_High,Cur_Map_Mean_High)

# Calculate area per land use
Areas<-as.data.frame(tapply(area(Map_ranger),Map_ranger[],sum)) %>%
  rename(sqM=`tapply(area(Map_ranger), Map_ranger[], sum)`) %>%
  mutate(Land_use=c("Agriculture","Grassland","Ground","Plantation","Shrubland","Water","Woodland"))
write.csv(Areas,"tmp/Areas.csv",row.names=FALSE)
rm(Map_ranger,Cur_Map_Mean,Areas)
Areas_Low<-as.data.frame(tapply(area(Map_ranger_Low),Map_ranger_Low[],sum)) %>%
  rename(sqM=`tapply(area(Map_ranger_Low), Map_ranger_Low[], sum)`) %>%
  mutate(Land_use=c("Agriculture","Grassland","Ground","Plantation","Shrubland","Water","Woodland"))
write.csv(Areas_Low,"tmp/Areas_Low.csv",row.names=FALSE)
rm(Map_ranger_Low,Cur_Map_Mean_Low,Areas_Low)
Areas_Intermediate<-as.data.frame(tapply(area(Map_ranger_Intermediate),Map_ranger_Intermediate[],sum)) %>%
  rename(sqM=`tapply(area(Map_ranger_Intermediate), Map_ranger_Intermediate[], `) %>%
  mutate(Land_use=c("Agriculture","Grassland","Ground","Plantation","Shrubland","Water","Woodland"))
write.csv(Areas_Intermediate,"tmp/Areas_Intermediate.csv",row.names=FALSE)
rm(Map_ranger_Intermediate,Cur_Map_Mean_Intermediate,Areas_Intermediate)
Areas_High<-as.data.frame(tapply(area(Map_ranger_High),Map_ranger_High[],sum)) %>%
  rename(sqM=`tapply(area(Map_ranger_High), Map_ranger_High[], sum)`) %>%
  mutate(Land_use=c("Agriculture","Grassland","Ground","Plantation","Shrubland","Water","Woodland"))
write.csv(Areas_High,"tmp/Areas_High.csv",row.names=FALSE)
rm(Map_ranger_High,Cur_Map_Mean_High,Areas_High)

# Plot area per land use 7.x5
Areas<-read.csv("tmp/Areas.csv",stringsAsFactors=FALSE)
Areas$sqHa<-Areas$sqM/10000
ggplot(Areas,aes(y=sqHa,x=Land_use)) + 
  ylab("Area in hectares") +
  xlab("Land use") +
  ggtitle("Importance of land use types for functional landscape connectivity") +
  geom_bar(stat="identity") +
  theme(plot.title=element_text(hjust = 0.5,size=14),
        axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        legend.title=element_text(size=13,colour="black"),
        legend.text=element_text(size=13,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(Areas)
Areas_Low<-read.csv("tmp/Areas_Low.csv",stringsAsFactors=FALSE)
Areas_Low$sqHa<-Areas_Low$sqM/10000
ggplot(Areas_Low,aes(y=sqHa,x=Land_use)) + 
  ylab("Area in hectares") +
  xlab("Land use") +
  ggtitle("Importance of land use types to improve functional connectivity
          for low conservation value grasshoppers") +
  geom_bar(stat="identity") +
  theme(plot.title=element_text(hjust = 0.5,size=14),
        axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        legend.title=element_text(size=13,colour="black"),
        legend.text=element_text(size=13,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(Areas_Low)
Areas_Intermediate<-read.csv("tmp/Areas_Intermediate.csv",stringsAsFactors=FALSE)
Areas_Intermediate$sqHa<-Areas_Intermediate$sqM/10000
ggplot(Areas_Intermediate,aes(y=sqHa,x=Land_use)) + 
  ylab("Area in hectares") +
  xlab("Land use") +
  ggtitle("Importance of land use types to improve functional connectivity
          for intermediate conservation value grasshoppers") +
  geom_bar(stat="identity") +
  theme(plot.title=element_text(hjust = 0.5,size=14),
        axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        legend.title=element_text(size=13,colour="black"),
        legend.text=element_text(size=13,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(Areas_Intermediate)
Areas_High<-read.csv("tmp/Areas_High.csv",stringsAsFactors=FALSE)
Areas_High$sqHa<-Areas_High$sqM/10000
ggplot(Areas_High,aes(y=sqHa,x=Land_use)) + 
  ylab("Area in hectares") +
  xlab("Land use") +
  ggtitle("Importance of land use types to improve functional connectivity
          for high conservation value grasshoppers") +
  geom_bar(stat="identity") +
  theme(plot.title=element_text(hjust = 0.5,size=14),
        axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        legend.title=element_text(size=13,colour="black"),
        legend.text=element_text(size=13,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())
rm(Areas_High)

############### Least cost paths ###############
# Create cost surfaces
ras_temp_Estate1<-raster(nrows=686,ncols=647,xmn=200040,xmx=212980,ymn=6706980,ymx=6720700,
                         res=10,crs=("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
ras_temp_Estate2<-raster(nrows=527,ncols=415,xmn=213040,xmx=221340,ymn=6719960,ymx=6730500,
                         res=10,crs=("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
ras_temp_Estate3<-raster(nrows=491,ncols=656,xmn=230140,xmx=243260,ymn=6756480,ymx=6766300,
                         res=10,crs=("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
ras_temp_Estate4<-raster(nrows=652,ncols=596,xmn=250160,xmx=262080,ymn=6752400,ymx=6765440,
                         res=10,crs=("+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs"))
Diversity_map_Estate_1<-raster("Circuitscape/Diversity_Estate_1.asc")
Diversity_map_Estate_1<-resample(Diversity_map_Estate_1,ras_temp_Estate1)
Diversity_map_Estate_1_cs<-gdistance::transition(Diversity_map_Estate_1,transitionFunction=mean,directions=8)
Diversity_map_Estate_2<-raster("Circuitscape/Diversity_Estate_2.asc")
Diversity_map_Estate_2<-resample(Diversity_map_Estate_2,ras_temp_Estate2)
Diversity_map_Estate_2_cs<-gdistance::transition(Diversity_map_Estate_2,transitionFunction=mean,directions=8)
Diversity_map_Estate_3<-raster("Circuitscape/Diversity_Estate_3.asc")
Diversity_map_Estate_3<-resample(Diversity_map_Estate_3,ras_temp_Estate3)
Diversity_map_Estate_3_cs<-gdistance::transition(Diversity_map_Estate_3,transitionFunction=mean,directions=8)
Diversity_map_Estate_4<-raster("Circuitscape/Diversity_Estate_4.asc")
Diversity_map_Estate_4<-resample(Diversity_map_Estate_4,ras_temp_Estate4)
Diversity_map_Estate_4_cs<-gdistance::transition(Diversity_map_Estate_4,transitionFunction=mean,directions=8)
Diversity_map_Estate_1_High<-raster("Circuitscape/Diversity_Estate_1_High.asc")
Diversity_map_Estate_1_High<-resample(Diversity_map_Estate_1_High,ras_temp_Estate1)
Diversity_map_Estate_1_High_cs<-gdistance::transition(Diversity_map_Estate_1_High,transitionFunction=mean,directions=8)
Diversity_map_Estate_2_High<-raster("Circuitscape/Diversity_Estate_2_High.asc")
Diversity_map_Estate_2_High<-resample(Diversity_map_Estate_2_High,ras_temp_Estate2)
Diversity_map_Estate_2_High_cs<-gdistance::transition(Diversity_map_Estate_2_High,transitionFunction=mean,directions=8)
Diversity_map_Estate_3_High<-raster("Circuitscape/Diversity_Estate_3_High.asc")
Diversity_map_Estate_3_High<-resample(Diversity_map_Estate_3_High,ras_temp_Estate3)
Diversity_map_Estate_3_High_cs<-gdistance::transition(Diversity_map_Estate_3_High,transitionFunction=mean,directions=8)
Diversity_map_Estate_4_High<-raster("Circuitscape/Diversity_Estate_4_High.asc")
Diversity_map_Estate_4_High<-resample(Diversity_map_Estate_4_High,ras_temp_Estate4)
Diversity_map_Estate_4_High_cs<-gdistance::transition(Diversity_map_Estate_4_High,transitionFunction=mean,directions=8)
Diversity_map_Estate_1_Intermediate<-raster("Circuitscape/Diversity_Estate_1_Intermediate.asc")
Diversity_map_Estate_1_Intermediate<-resample(Diversity_map_Estate_1_Intermediate,ras_temp_Estate1)
Diversity_map_Estate_1_Intermediate_cs<-gdistance::transition(Diversity_map_Estate_1_Intermediate,transitionFunction=mean,directions=8)
Diversity_map_Estate_2_Intermediate<-raster("Circuitscape/Diversity_Estate_2_Intermediate.asc")
Diversity_map_Estate_2_Intermediate<-resample(Diversity_map_Estate_2_Intermediate,ras_temp_Estate2)
Diversity_map_Estate_2_Intermediate_cs<-gdistance::transition(Diversity_map_Estate_2_Intermediate,transitionFunction=mean,directions=8)
Diversity_map_Estate_3_Intermediate<-raster("Circuitscape/Diversity_Estate_3_Intermediate.asc")
Diversity_map_Estate_3_Intermediate<-resample(Diversity_map_Estate_3_Intermediate,ras_temp_Estate3)
Diversity_map_Estate_3_Intermediate_cs<-gdistance::transition(Diversity_map_Estate_3_Intermediate,transitionFunction=mean,directions=8)
Diversity_map_Estate_4_Intermediate<-raster("Circuitscape/Diversity_Estate_4_Intermediate.asc")
Diversity_map_Estate_4_Intermediate<-resample(Diversity_map_Estate_4_Intermediate,ras_temp_Estate4)
Diversity_map_Estate_4_Intermediate_cs<-gdistance::transition(Diversity_map_Estate_4_Intermediate,transitionFunction=mean,directions=8)
Diversity_map_Estate_1_Low<-raster("Circuitscape/Diversity_Estate_1_Low.asc")
Diversity_map_Estate_1_Low<-resample(Diversity_map_Estate_1_Low,ras_temp_Estate1)
Diversity_map_Estate_1_Low_cs<-gdistance::transition(Diversity_map_Estate_1_Low,transitionFunction=mean,directions=8)
Diversity_map_Estate_2_Low<-raster("Circuitscape/Diversity_Estate_2_Low.asc")
Diversity_map_Estate_2_Low<-resample(Diversity_map_Estate_2_Low,ras_temp_Estate2)
Diversity_map_Estate_2_Low_cs<-gdistance::transition(Diversity_map_Estate_2_Low,transitionFunction=mean,directions=8)
Diversity_map_Estate_3_Low<-raster("Circuitscape/Diversity_Estate_3_Low.asc")
Diversity_map_Estate_3_Low<-resample(Diversity_map_Estate_3_Low,ras_temp_Estate3)
Diversity_map_Estate_3_Low_cs<-gdistance::transition(Diversity_map_Estate_3_Low,transitionFunction=mean,directions=8)
Diversity_map_Estate_4_Low<-raster("Circuitscape/Diversity_Estate_4_Low.asc")
Diversity_map_Estate_4_Low<-resample(Diversity_map_Estate_4_Low,ras_temp_Estate4)
Diversity_map_Estate_4_Low_cs<-gdistance::transition(Diversity_map_Estate_4_Low,transitionFunction=mean,directions=8)

# Convert focal habitat to points
Focal_Areas_Estate_1<-raster("Circuitscape/Focal_Areas_Estate_1.asc")
Focal_Areas_Estate_1<-rasterToPolygons(Focal_Areas_Estate_1,dissolve=TRUE) 
Focal_Areas_Estate_1<-gCentroid(Focal_Areas_Estate_1,byid=TRUE)
Focal_Areas_Estate_2<-raster("Circuitscape/Focal_Areas_Estate_2.asc")
Focal_Areas_Estate_2<-rasterToPolygons(Focal_Areas_Estate_2,dissolve=TRUE) 
Focal_Areas_Estate_2<-gCentroid(Focal_Areas_Estate_2,byid=TRUE)
Focal_Areas_Estate_3<-raster("Circuitscape/Focal_Areas_Estate_3.asc")
Focal_Areas_Estate_3<-rasterToPolygons(Focal_Areas_Estate_3,dissolve=TRUE) 
Focal_Areas_Estate_3<-gCentroid(Focal_Areas_Estate_3,byid=TRUE)
Focal_Areas_Estate_4<-raster("Circuitscape/Focal_Areas_Estate_4.asc")
Focal_Areas_Estate_4<-rasterToPolygons(Focal_Areas_Estate_4,dissolve=TRUE) 
Focal_Areas_Estate_4<-gCentroid(Focal_Areas_Estate_4,byid=TRUE)
Focal_Areas_Estate_1_High<-raster("Circuitscape/Focal_Areas_Estate_1_High.asc")
Focal_Areas_Estate_1_High<-rasterToPolygons(Focal_Areas_Estate_1_High,dissolve=TRUE) 
Focal_Areas_Estate_1_High<-gCentroid(Focal_Areas_Estate_1_High,byid=TRUE)
Focal_Areas_Estate_2_High<-raster("Circuitscape/Focal_Areas_Estate_2_High.asc")
Focal_Areas_Estate_2_High<-rasterToPolygons(Focal_Areas_Estate_2_High,dissolve=TRUE) 
Focal_Areas_Estate_2_High<-gCentroid(Focal_Areas_Estate_2_High,byid=TRUE)
Focal_Areas_Estate_3_High<-raster("Circuitscape/Focal_Areas_Estate_3_High.asc")
Focal_Areas_Estate_3_High<-rasterToPolygons(Focal_Areas_Estate_3_High,dissolve=TRUE) 
Focal_Areas_Estate_3_High<-gCentroid(Focal_Areas_Estate_3_High,byid=TRUE)
Focal_Areas_Estate_4_High<-raster("Circuitscape/Focal_Areas_Estate_4_High.asc")
Focal_Areas_Estate_4_High<-rasterToPolygons(Focal_Areas_Estate_4_High,dissolve=TRUE) 
Focal_Areas_Estate_4_High<-gCentroid(Focal_Areas_Estate_4_High,byid=TRUE)
Focal_Areas_Estate_1_Intermediate<-raster("Circuitscape/Focal_Areas_Estate_1_Intermediate.asc")
Focal_Areas_Estate_1_Intermediate<-rasterToPolygons(Focal_Areas_Estate_1_Intermediate,dissolve=TRUE) 
Focal_Areas_Estate_1_Intermediate<-gCentroid(Focal_Areas_Estate_1_Intermediate,byid=TRUE)
Focal_Areas_Estate_2_Intermediate<-raster("Circuitscape/Focal_Areas_Estate_2_Intermediate.asc")
Focal_Areas_Estate_2_Intermediate<-rasterToPolygons(Focal_Areas_Estate_2_Intermediate,dissolve=TRUE) 
Focal_Areas_Estate_2_Intermediate<-gCentroid(Focal_Areas_Estate_2_Intermediate,byid=TRUE)
Focal_Areas_Estate_3_Intermediate<-raster("Circuitscape/Focal_Areas_Estate_3_Intermediate.asc")
Focal_Areas_Estate_3_Intermediate<-rasterToPolygons(Focal_Areas_Estate_3_Intermediate,dissolve=TRUE) 
Focal_Areas_Estate_3_Intermediate<-gCentroid(Focal_Areas_Estate_3_Intermediate,byid=TRUE)
Focal_Areas_Estate_4_Intermediate<-raster("Circuitscape/Focal_Areas_Estate_4_Intermediate.asc")
Focal_Areas_Estate_4_Intermediate<-rasterToPolygons(Focal_Areas_Estate_4_Intermediate,dissolve=TRUE) 
Focal_Areas_Estate_4_Intermediate<-gCentroid(Focal_Areas_Estate_4_Intermediate,byid=TRUE)
Focal_Areas_Estate_1_Low<-raster("Circuitscape/Focal_Areas_Estate_1_Low.asc")
Focal_Areas_Estate_1_Low<-rasterToPolygons(Focal_Areas_Estate_1_Low,dissolve=TRUE) 
Focal_Areas_Estate_1_Low<-gCentroid(Focal_Areas_Estate_1_Low,byid=TRUE)
Focal_Areas_Estate_2_Low<-raster("Circuitscape/Focal_Areas_Estate_2_Low.asc")
Focal_Areas_Estate_2_Low<-rasterToPolygons(Focal_Areas_Estate_2_Low,dissolve=TRUE) 
Focal_Areas_Estate_2_Low<-gCentroid(Focal_Areas_Estate_2_Low,byid=TRUE)
Focal_Areas_Estate_3_Low<-raster("Circuitscape/Focal_Areas_Estate_3_Low.asc")
Focal_Areas_Estate_3_Low<-rasterToPolygons(Focal_Areas_Estate_3_Low,dissolve=TRUE) 
Focal_Areas_Estate_3_Low<-gCentroid(Focal_Areas_Estate_3_Low,byid=TRUE)
Focal_Areas_Estate_4_Low<-raster("Circuitscape/Focal_Areas_Estate_4_Low.asc")
Focal_Areas_Estate_4_Low<-rasterToPolygons(Focal_Areas_Estate_4_Low,dissolve=TRUE) 
Focal_Areas_Estate_4_Low<-gCentroid(Focal_Areas_Estate_4_Low,byid=TRUE)

# Calculate least cost paths
Overall_Estate_1<-create_FETE_lcps(Diversity_map_Estate_1_cs,Focal_Areas_Estate_1,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_1_Density<-create_lcp_density(Overall_Estate_1,Diversity_map_Estate_1,rescale=FALSE)
Overall_Estate_1_Density[Overall_Estate_1_Density==0]<-NA
Overall_Estate_2<-create_FETE_lcps(Diversity_map_Estate_2_cs,Focal_Areas_Estate_2,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_2_Density<-create_lcp_density(Overall_Estate_2,Diversity_map_Estate_2,rescale=FALSE)
Overall_Estate_2_Density[Overall_Estate_2_Density==0]<-NA
Overall_Estate_3<-create_FETE_lcps(Diversity_map_Estate_3_cs,Focal_Areas_Estate_3,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_3_Density<-create_lcp_density(Overall_Estate_3,Diversity_map_Estate_3,rescale=FALSE)
Overall_Estate_3_Density[Overall_Estate_3_Density==0]<-NA
Overall_Estate_4<-create_FETE_lcps(Diversity_map_Estate_4_cs,Focal_Areas_Estate_4,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_4_Density<-create_lcp_density(Overall_Estate_4,Diversity_map_Estate_4,rescale=FALSE)
Overall_Estate_4_Density[Overall_Estate_4_Density==0]<-NA
Overall_Estate_1_High<-create_FETE_lcps(Diversity_map_Estate_1_High_cs,Focal_Areas_Estate_1_High,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_1_High_Density<-create_lcp_density(Overall_Estate_1_High,Diversity_map_Estate_1_High,rescale=FALSE)
Overall_Estate_1_High_Density[Overall_Estate_1_High_Density==0]<-NA
Overall_Estate_2_High<-create_FETE_lcps(Diversity_map_Estate_2_High_cs,Focal_Areas_Estate_2_High,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_2_High_Density<-create_lcp_density(Overall_Estate_2_High,Diversity_map_Estate_2_High,rescale=FALSE)
Overall_Estate_2_High_Density[Overall_Estate_2_High_Density==0]<-NA
Overall_Estate_3_High<-create_FETE_lcps(Diversity_map_Estate_3_High_cs,Focal_Areas_Estate_3_High,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_3_High_Density<-create_lcp_density(Overall_Estate_3_High,Diversity_map_Estate_3_High,rescale=FALSE)
Overall_Estate_3_High_Density[Overall_Estate_3_High_Density==0]<-NA
Overall_Estate_4_High<-create_FETE_lcps(Diversity_map_Estate_4_High_cs,Focal_Areas_Estate_4_High,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_4_High_Density<-create_lcp_density(Overall_Estate_4_High,Diversity_map_Estate_4_High,rescale=FALSE)
Overall_Estate_4_High_Density[Overall_Estate_4_High_Density==0]<-NA
Overall_Estate_1_Intermediate<-create_FETE_lcps(Diversity_map_Estate_1_Intermediate_cs,Focal_Areas_Estate_1_Intermediate,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_1_Intermediate_Density<-create_lcp_density(Overall_Estate_1_Intermediate,Diversity_map_Estate_1_Intermediate,rescale=FALSE)
Overall_Estate_1_Intermediate_Density[Overall_Estate_1_Intermediate_Density==0]<-NA
Overall_Estate_2_Intermediate<-create_FETE_lcps(Diversity_map_Estate_2_Intermediate_cs,Focal_Areas_Estate_2_Intermediate,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_2_Intermediate_Density<-create_lcp_density(Overall_Estate_2_Intermediate,Diversity_map_Estate_2_Intermediate,rescale=FALSE)
Overall_Estate_2_Intermediate_Density[Overall_Estate_2_Intermediate_Density==0]<-NA
Overall_Estate_3_Intermediate<-create_FETE_lcps(Diversity_map_Estate_3_Intermediate_cs,Focal_Areas_Estate_3_Intermediate,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_3_Intermediate_Density<-create_lcp_density(Overall_Estate_3_Intermediate,Diversity_map_Estate_3_Intermediate,rescale=FALSE)
Overall_Estate_3_Intermediate_Density[Overall_Estate_3_Intermediate_Density==0]<-NA
Overall_Estate_4_Intermediate<-create_FETE_lcps(Diversity_map_Estate_4_Intermediate_cs,Focal_Areas_Estate_4_Intermediate,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_4_Intermediate_Density<-create_lcp_density(Overall_Estate_4_Intermediate,Diversity_map_Estate_4_Intermediate,rescale=FALSE)
Overall_Estate_4_Intermediate_Density[Overall_Estate_4_Intermediate_Density==0]<-NA
Overall_Estate_1_Low<-create_FETE_lcps(Diversity_map_Estate_1_Low_cs,Focal_Areas_Estate_1_Low,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_1_Low_Density<-create_lcp_density(Overall_Estate_1_Low,Diversity_map_Estate_1_Low,rescale=FALSE)
Overall_Estate_1_Low_Density[Overall_Estate_1_Low_Density==0]<-NA
Overall_Estate_2_Low<-create_FETE_lcps(Diversity_map_Estate_2_Low_cs,Focal_Areas_Estate_2_Low,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_2_Low_Density<-create_lcp_density(Overall_Estate_2_Low,Diversity_map_Estate_2_Low,rescale=FALSE)
Overall_Estate_2_Low_Density[Overall_Estate_2_Low_Density==0]<-NA
Overall_Estate_3_Low<-create_FETE_lcps(Diversity_map_Estate_3_Low_cs,Focal_Areas_Estate_3_Low,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_3_Low_Density<-create_lcp_density(Overall_Estate_3_Low,Diversity_map_Estate_3_Low,rescale=FALSE)
Overall_Estate_3_Low_Density[Overall_Estate_3_Low_Density==0]<-NA
Overall_Estate_4_Low<-create_FETE_lcps(Diversity_map_Estate_4_Low_cs,Focal_Areas_Estate_4_Low,cost_distance=FALSE,parallel=TRUE)
Overall_Estate_4_Low_Density<-create_lcp_density(Overall_Estate_4_Low,Diversity_map_Estate_4_Low,rescale=FALSE)
Overall_Estate_4_Low_Density[Overall_Estate_4_Low_Density==0]<-NA
writeRaster(Overall_Estate_1_Density,filename="Products/Raster_Products/Overall_Estate_1_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_2_Density,filename="Products/Raster_Products/Overall_Estate_2_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_3_Density,filename="Products/Raster_Products/Overall_Estate_3_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_4_Density,filename="Products/Raster_Products/Overall_Estate_4_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_1_High_Density,filename="Products/Raster_Products/Overall_Estate_1_High_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_2_High_Density,filename="Products/Raster_Products/Overall_Estate_2_High_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_3_High_Density,filename="Products/Raster_Products/Overall_Estate_3_High_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_4_High_Density,filename="Products/Raster_Products/Overall_Estate_4_High_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_1_Intermediate_Density,filename="Products/Raster_Products/Overall_Estate_1_Intermediate_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_2_Intermediate_Density,filename="Products/Raster_Products/Overall_Estate_2_Intermediate_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_3_Intermediate_Density,filename="Products/Raster_Products/Overall_Estate_3_Intermediate_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_4_Intermediate_Density,filename="Products/Raster_Products/Overall_Estate_4_Intermediate_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_1_Low_Density,filename="Products/Raster_Products/Overall_Estate_1_Low_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_2_Low_Density,filename="Products/Raster_Products/Overall_Estate_2_Low_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_3_Low_Density,filename="Products/Raster_Products/Overall_Estate_3_Low_Density.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_4_Low_Density,filename="Products/Raster_Products/Overall_Estate_4_Low_Density.tif",format="GTiff",overwrite=TRUE)

# Prune networks
Overall_Estate_1_Density<-raster("Products/Raster_Products/Overall_Estate_1_Density.tif")
Overall_Estate_1_Density[Overall_Estate_1_Density<50]<-NA
Overall_Estate_2_Density<-raster("Products/Raster_Products/Overall_Estate_2_Density.tif")
Overall_Estate_2_Density[Overall_Estate_2_Density<15]<-NA
Overall_Estate_3_Density<-raster("Products/Raster_Products/Overall_Estate_3_Density.tif")
Overall_Estate_3_Density[Overall_Estate_3_Density<50]<-NA
Overall_Estate_4_Density<-raster("Products/Raster_Products/Overall_Estate_4_Density.tif")
Overall_Estate_4_Density[Overall_Estate_4_Density<1]<-NA
Overall_Estate_1_High_Density<-raster("Products/Raster_Products/Overall_Estate_1_High_Density.tif")
Overall_Estate_1_High_Density[Overall_Estate_1_High_Density<30]<-NA
Overall_Estate_2_High_Density<-raster("Products/Raster_Products/Overall_Estate_2_High_Density.tif")
Overall_Estate_2_High_Density[Overall_Estate_2_High_Density<10]<-NA
Overall_Estate_3_High_Density<-raster("Products/Raster_Products/Overall_Estate_3_High_Density.tif")
Overall_Estate_3_High_Density[Overall_Estate_3_High_Density<30]<-NA
Overall_Estate_4_High_Density<-raster("Products/Raster_Products/Overall_Estate_4_High_Density.tif")
Overall_Estate_4_High_Density[Overall_Estate_4_High_Density<1]<-NA
Overall_Estate_1_Intermediate_Density<-raster("Products/Raster_Products/Overall_Estate_1_Intermediate_Density.tif")
Overall_Estate_1_Intermediate_Density[Overall_Estate_1_Intermediate_Density<45]<-NA
Overall_Estate_2_Intermediate_Density<-raster("Products/Raster_Products/Overall_Estate_2_Intermediate_Density.tif")
Overall_Estate_2_Intermediate_Density[Overall_Estate_2_Intermediate_Density<10]<-NA
Overall_Estate_3_Intermediate_Density<-raster("Products/Raster_Products/Overall_Estate_3_Intermediate_Density.tif")
Overall_Estate_3_Intermediate_Density[Overall_Estate_3_Intermediate_Density<45]<-NA
Overall_Estate_4_Intermediate_Density<-raster("Products/Raster_Products/Overall_Estate_4_Intermediate_Density.tif")
Overall_Estate_4_Intermediate_Density[Overall_Estate_4_Intermediate_Density<3]<-NA
Overall_Estate_1_Low_Density<-raster("Products/Raster_Products/Overall_Estate_1_Low_Density.tif")
Overall_Estate_1_Low_Density[Overall_Estate_1_Low_Density<45]<-NA
Overall_Estate_2_Low_Density<-raster("Products/Raster_Products/Overall_Estate_2_Low_Density.tif")
Overall_Estate_2_Low_Density[Overall_Estate_2_Low_Density<10]<-NA
Overall_Estate_3_Low_Density<-raster("Products/Raster_Products/Overall_Estate_3_Low_Density.tif")
Overall_Estate_3_Low_Density[Overall_Estate_3_Low_Density<45]<-NA
Overall_Estate_4_Low_Density<-raster("Products/Raster_Products/Overall_Estate_4_Low_Density.tif")
Overall_Estate_4_Low_Density[Overall_Estate_4_Low_Density<1]<-NA

# Merge
Overall_Density<-raster::merge(Overall_Estate_1_Density,Overall_Estate_2_Density,Overall_Estate_3_Density,Overall_Estate_4_Density)
writeRaster(Overall_Density,filename="Products/Raster_Products/Overall_Density.tif",format="GTiff",overwrite=TRUE)
Overall_Density_High<-raster::merge(Overall_Estate_1_High_Density,Overall_Estate_2_High_Density,Overall_Estate_3_High_Density,Overall_Estate_4_High_Density)
writeRaster(Overall_Density_High,filename="Products/Raster_Products/Overall_Density_High.tif",format="GTiff",overwrite=TRUE)
Overall_Density_Intermediate<-raster::merge(Overall_Estate_1_Intermediate_Density,Overall_Estate_2_Intermediate_Density,Overall_Estate_3_Intermediate_Density,Overall_Estate_4_Intermediate_Density)
writeRaster(Overall_Density_Intermediate,filename="Products/Raster_Products/Overall_Density_Intermediate.tif",format="GTiff",overwrite=TRUE)
Overall_Density_Low<-raster::merge(Overall_Estate_1_Low_Density,Overall_Estate_2_Low_Density,Overall_Estate_3_Low_Density,Overall_Estate_4_Low_Density)
writeRaster(Overall_Density_Low,filename="Products/Raster_Products/Overall_Density_Low.tif",format="GTiff",overwrite=TRUE)
rm(Diversity_map_Estate_1,Diversity_map_Estate_1_cs,Focal_Areas_Estate_1,Overall_Estate_1,Overall_Estate_1_Density,
   Diversity_map_Estate_2,Diversity_map_Estate_2_cs,Focal_Areas_Estate_2,Overall_Estate_2,Overall_Estate_2_Density,
   Diversity_map_Estate_3,Diversity_map_Estate_3_cs,Focal_Areas_Estate_3,Overall_Estate_3,Overall_Estate_3_Density,
   Diversity_map_Estate_4,Diversity_map_Estate_4_cs,Focal_Areas_Estate_4,Overall_Estate_4,Overall_Estate_4_Density,
   Diversity_map_Estate_1_High,Diversity_map_Estate_1_High_cs,Focal_Areas_Estate_1_High,Overall_Estate_1_High,Overall_Estate_1_High_Density,
   Diversity_map_Estate_2_High,Diversity_map_Estate_2_High_cs,Focal_Areas_Estate_2_High,Overall_Estate_2_High,Overall_Estate_2_High_Density,
   Diversity_map_Estate_3_High,Diversity_map_Estate_3_High_cs,Focal_Areas_Estate_3_High,Overall_Estate_3_High,Overall_Estate_3_High_Density,
   Diversity_map_Estate_4_High,Diversity_map_Estate_4_High_cs,Focal_Areas_Estate_4_High,Overall_Estate_4_High,Overall_Estate_4_High_Density,
   Diversity_map_Estate_1_Intermediate,Diversity_map_Estate_1_Intermediate_cs,Focal_Areas_Estate_1_Intermediate,Overall_Estate_1_Intermediate,Overall_Estate_1_Intermediate_Density,
   Diversity_map_Estate_2_Intermediate,Diversity_map_Estate_2_Intermediate_cs,Focal_Areas_Estate_2_Intermediate,Overall_Estate_2_Intermediate,Overall_Estate_2_Intermediate_Density,
   Diversity_map_Estate_3_Intermediate,Diversity_map_Estate_3_Intermediate_cs,Focal_Areas_Estate_3_Intermediate,Overall_Estate_3_Intermediate,Overall_Estate_3_Intermediate_Density,
   Diversity_map_Estate_4_Intermediate,Diversity_map_Estate_4_Intermediate_cs,Focal_Areas_Estate_4_Intermediate,Overall_Estate_4_Intermediate,Overall_Estate_4_Intermediate_Density,
   Diversity_map_Estate_1_Low,Diversity_map_Estate_1_Low_cs,Focal_Areas_Estate_1_Low,Overall_Estate_1_Low,Overall_Estate_1_Low_Density,
   Diversity_map_Estate_2_Low,Diversity_map_Estate_2_Low_cs,Focal_Areas_Estate_2_Low,Overall_Estate_2_Low,Overall_Estate_2_Low_Density,
   Diversity_map_Estate_3_Low,Diversity_map_Estate_3_Low_cs,Focal_Areas_Estate_3_Low,Overall_Estate_3_Low,Overall_Estate_3_Low_Density,
   Diversity_map_Estate_4_Low,Diversity_map_Estate_4_Low_cs,Focal_Areas_Estate_4_Low,Overall_Estate_4_Low,Overall_Estate_4_Low_Density,
   Overall_Density,Overall_Density_High,Overall_Density_Intermediate,Overall_Density_Low,ras_temp_Estate1,ras_temp_Estate2,ras_temp_Estate3,ras_temp_Estate4)

############### Identify restoration areas ###############
# Load least cost networks
Overall_Density<-raster("Products/Raster_Products/Overall_Density.tif")
Overall_Density_High<-raster("Products/Raster_Products/Overall_Density_High.tif")
Overall_Density_Intermediate<-raster("Products/Raster_Products/Overall_Density_Intermediate.tif")
Overall_Density_Low<-raster("Products/Raster_Products/Overall_Density_Low.tif")

# Load current maps and select lowest 85%
Extent<-shapefile("Shapefiles/Plantations_Clipped.shp",stringsAsFactors=FALSE)
Cur_Map<-raster("Products/Raster_Products/Cur_Map.tif")
Cur_Map_Low<-raster("Products/Raster_Products/Cur_Map_Low.tif")
Cur_Map_Intermediate<-raster("Products/Raster_Products/Cur_Map_Intermediate.tif")
Cur_Map_High<-raster("Products/Raster_Products/Cur_Map_High.tif")
Cur_Map<-crop(Cur_Map,Extent)
Cur_Map<-mask(Cur_Map,Extent)
Cur_Map_Low<-crop(Cur_Map_Low,Extent)
Cur_Map_Low<-mask(Cur_Map_Low,Extent)
Cur_Map_Intermediate<-crop(Cur_Map_Intermediate,Extent)
Cur_Map_Intermediate<-mask(Cur_Map_Intermediate,Extent)
Cur_Map_High<-crop(Cur_Map_High,Extent)
Cur_Map_High<-mask(Cur_Map_High,Extent)
Cur_Map<-Cur_Map[[1]]<=quantile(Cur_Map,.85)
Cur_Map[Cur_Map==0]<-NA
Cur_Map_Low<-Cur_Map_Low[[1]]<=quantile(Cur_Map_Low,.85)
Cur_Map_Low[Cur_Map_Low==0]<-NA
Cur_Map_Intermediate<-Cur_Map_Intermediate[[1]]<=quantile(Cur_Map_Intermediate,.85)
Cur_Map_Intermediate[Cur_Map_Intermediate==0]<-NA
Cur_Map_High<-Cur_Map_High[[1]]<=quantile(Cur_Map_High,.85)
Cur_Map_High[Cur_Map_High==0]<-NA

# Select network over low current areas
Overall_Density<-crop(Overall_Density,Cur_Map)
Overall_Density<-mask(Overall_Density,Cur_Map)
Overall_Density_High<-crop(Overall_Density_High,Cur_Map_High)
Overall_Density_High<-mask(Overall_Density_High,Cur_Map_High)
Overall_Density_Intermediate<-crop(Overall_Density_Intermediate,Cur_Map_Intermediate)
Overall_Density_Intermediate<-mask(Overall_Density_Intermediate,Cur_Map_Intermediate)
Overall_Density_Low<-crop(Overall_Density_Low,Cur_Map_Low)
Overall_Density_Low<-mask(Overall_Density_Low,Cur_Map_Low)
rm(Cur_Map_High,Cur_Map_Intermediate,Cur_Map_Low,Cur_Map)

# Clip paths over plantations
Map_ranger<-raster("Products/Raster_Products/Map_ranger.tif")
Map_ranger<-crop(Map_ranger,Extent)
Map_ranger<-mask(Map_ranger,Extent)
Map_ranger[Map_ranger==4]<-NA
Overall_Density<-crop(Overall_Density,Map_ranger)
Overall_Density<-mask(Overall_Density,Map_ranger)
Overall_Density_High<-crop(Overall_Density_High,Map_ranger)
Overall_Density_High<-mask(Overall_Density_High,Map_ranger)
Overall_Density_Intermediate<-crop(Overall_Density_Intermediate,Map_ranger)
Overall_Density_Intermediate<-mask(Overall_Density_Intermediate,Map_ranger)
Overall_Density_Low<-crop(Overall_Density_Low,Map_ranger)
Overall_Density_Low<-mask(Overall_Density_Low,Map_ranger)
rm(Extent,Map_ranger)

# Subset to include highly used corridors only
Overall_Density[Overall_Density<522]<-NA
Overall_Density_High[Overall_Density_High<233]<-NA
Overall_Density_Intermediate[Overall_Density_Intermediate<654]<-NA
Overall_Density_Low[Overall_Density_Low<620]<-NA
writeRaster(Overall_Density,filename="Products/Raster_Products/Overall_Density_Restore.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Density_High,filename="Products/Raster_Products/Overall_Density_High_Restore.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Density_Intermediate,filename="Products/Raster_Products/Overall_Density_Intermediate_Restore.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Density_Low,filename="Products/Raster_Products/Overall_Density_Low_Restore.tif",format="GTiff",overwrite=TRUE)

# Convert to points for visualization
Overall_Density_Point<-rasterToPoints(Overall_Density,spatial=TRUE)
writeOGR(Overall_Density_Point,"Products/Shapefiles/","Overall_Density_Point_Restore",driver="ESRI Shapefile")
Overall_Density_High_Point<-rasterToPoints(Overall_Density_High,spatial=TRUE)
writeOGR(Overall_Density_High_Point,"Products/Shapefiles/","Overall_Density_High_Point_Restore",driver="ESRI Shapefile")
Overall_Density_Intermediate_Point<-rasterToPoints(Overall_Density_Intermediate,spatial=TRUE)
writeOGR(Overall_Density_Intermediate_Point,"Products/Shapefiles/","Overall_Density_Intermediate_Point_Restore",driver="ESRI Shapefile")
Overall_Density_Low_Point<-rasterToPoints(Overall_Density_Low,spatial=TRUE)
writeOGR(Overall_Density_Low_Point,"Products/Shapefiles/","Overall_Density_Low_Point_Restore",driver="ESRI Shapefile")
Overall_Density<-raster("Products/Raster_Products/Overall_Density.tif")
Overall_Density_High<-raster("Products/Raster_Products/Overall_Density_High.tif")
Overall_Density_Intermediate<-raster("Products/Raster_Products/Overall_Density_Intermediate.tif")
Overall_Density_Low<-raster("Products/Raster_Products/Overall_Density_Low.tif")
Overall_Density_Point<-rasterToPoints(Overall_Density,spatial=TRUE)
writeOGR(Overall_Density_Point,"Products/Shapefiles/","Overall_Density_Point",driver="ESRI Shapefile")
Overall_Density_High_Point<-rasterToPoints(Overall_Density_High,spatial=TRUE)
writeOGR(Overall_Density_High_Point,"Products/Shapefiles/","Overall_Density_High_Point",driver="ESRI Shapefile")
Overall_Density_Intermediate_Point<-rasterToPoints(Overall_Density_Intermediate,spatial=TRUE)
writeOGR(Overall_Density_Intermediate_Point,"Products/Shapefiles/","Overall_Density_Intermediate_Point",driver="ESRI Shapefile")
Overall_Density_Low_Point<-rasterToPoints(Overall_Density_Low,spatial=TRUE)
writeOGR(Overall_Density_Low_Point,"Products/Shapefiles/","Overall_Density_Low_Point",driver="ESRI Shapefile")
rm(Overall_Density_Point,Overall_Density_High_Point,Overall_Density_Intermediate_Point,Overall_Density_Low_Point,
   Overall_Density,Overall_Density_High,Overall_Density_Intermediate,Overall_Density_Low)

############### Turnover vs nestedness ###############
# Load biological data and subset per plantation farm
Bio<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
Bio<-Bio[,colSums(Bio!=0)>3]
Bio[Bio!=0]<-1
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Plantation<-GPS$Plantation
Bio<-as.data.frame(cbind("Plantation"=Plantation,Bio))
Estate1<-subset(Bio,Plantation=="Estate1")
Estate1$Plantation=NULL
Estate2<-subset(Bio,Plantation=="Estate2")
Estate2$Plantation=NULL
Estate3<-subset(Bio,Plantation=="Estate3")
Estate3$Plantation=NULL
Estate4<-subset(Bio,Plantation=="Estate4")
Estate4$Plantation=NULL
Bio$Plantation=NULL
rm(GPS,Plantation)

# Calculate turnover and nestedness
Estate1_core<-beta.multi(betapart.core(Estate1),index.family="jac")
Estate1_core<-as.data.frame(cbind("Beta diversity"=Estate1_core$beta.JAC,
                                  "Nestedness"=Estate1_core$beta.JNE,
                                  "Turnover"=Estate1_core$beta.JTU,
                                  "Estate"="Estate1"))
Estate2_core<-beta.multi(betapart.core(Estate2),index.family="jac")
Estate2_core<-as.data.frame(cbind("Beta diversity"=Estate2_core$beta.JAC,
                                  "Nestedness"=Estate2_core$beta.JNE,
                                  "Turnover"=Estate2_core$beta.JTU,
                                  "Estate"="Estate2"))
Estate3_core<-beta.multi(betapart.core(Estate3),index.family="jac")
Estate3_core<-as.data.frame(cbind("Beta diversity"=Estate3_core$beta.JAC,
                                  "Nestedness"=Estate3_core$beta.JNE,
                                  "Turnover"=Estate3_core$beta.JTU,
                                  "Estate"="Estate3"))
Estate4_core<-beta.multi(betapart.core(Estate4),index.family="jac")
Estate4_core<-as.data.frame(cbind("Beta diversity"=Estate4_core$beta.JAC,
                                  "Nestedness"=Estate4_core$beta.JNE,
                                  "Turnover"=Estate4_core$beta.JTU,
                                  "Estate"="Estate4"))
Bio_core<-beta.multi(betapart.core(Bio),index.family="jac")
Bio_core<-as.data.frame(cbind("Beta diversity"=Bio_core$beta.JAC,
                              "Nestedness"=Bio_core$beta.JNE,
                              "Turnover"=Bio_core$beta.JTU,
                              "Estate"="All"))
betaDiv<-as.data.frame(rbind(Estate1_core,Estate2_core,Estate3_core,Estate4_core,Bio_core))
write.csv(betaDiv,"Products/Excel_Sheets_Products/betaDiv_Data.csv",row.names=FALSE)
rm(Areas,betaDiv,Bio,Bio_core,Estate1,Estate1_core,Estate2,Estate2_core,Estate3,Estate3_core,Estate4,Estate4_core)

############### Drivers of turnover gdm ###############
# Load data environmental data
Cur_Map<-raster("Products/Raster_Products/Cur_Map.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist.tif")
Map_ranger<-raster("Products/Raster_Products/Map_ranger.tif")
Max_NBR<-stack("Products/Raster_Products/Max_NBR.tif")
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
Env<-stack(Sent_NDVI,Aspect,Max_NBR,Cur_Map,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
rm(Sent_NDVI,Aspect,Max_NBR,Cur_Map,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
names(Env)
names(Env)<-c("NDVI","Aspect","Fire_Hist","Func_Conn","Land_Use","Drain_Dist","Cor_Width","Elv")

# Spatial autocorrelation
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Bio<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
Bio<-Bio[,colSums(Bio!=0)>3]
provi<-deldir::deldir(GPS[2:3])
provi.neig<-neig(edges=as.matrix(provi$delsgs[,5:6]))
maf.listw<-spdep::nb2listw(neig2nb(provi.neig))
multispati.rtest((dudi.pca(Bio,scannf=FALSE)),maf.listw) #Spatial autocorrelation
rm(provi,provi.neig,maf.listw)

# Load biological data
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Site_ID<-GPS$Site_ID
GPS<-SpatialPointsDataFrame(GPS[,2:3],GPS,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
GPS<-spTransform(GPS,CRS="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Proj_Cor<-as.data.frame(GPS@coords)
Bio<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
Bio<-Bio[,colSums(Bio!=0)>3]
Bio[Bio!=0]<-1
Bio<-as.data.frame(cbind(Site_ID,Proj_Cor,Bio))
rm(GPS,Proj_Cor,Site_ID)

# Prepare data for gdm
Data<-formatsitepair(bioData=Bio,bioFormat=1,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# Variable importance
model_1_Imp<-gdm.varImp(spTable=Data,geo=TRUE,nPerm=100,fullModelOnly=FALSE)
barplot(sort(model_1_Imp[[2]][,1],decreasing=T))
model_1_Imp[[1]]
model_1_Imp[[2]]
model_1_Imp[[3]]
model_1_Imp[[4]]

# Select only significant variables
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Env<-stack(Cor_Width,SUDEM)
rm(Cor_Width,SUDEM)
names(Env)
names(Env)<-c("Cor_Width","Elv")
Data<-formatsitepair(bioData=Bio,bioFormat=1,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# GDM
model_1<-gdm(data=Data,geo=TRUE)
summary.gdm(model_1)
plot(model_1,plot.layout=c(3,2))
model_1.splineDat<-isplineExtract(model_1)
par(mfrow=c(2,2))
plot(model_1.splineDat$x[,"Cor_Width"],model_1.splineDat$y[,"Cor_Width"],ylim=c(0,0.7),
     lwd=2,type="l",xlab="Edge distance",ylab="Ecological distance")
plot(model_1.splineDat$x[,"Elv"],model_1.splineDat$y[,"Elv"],ylim=c(0,0.7),
     lwd=2,type="l",xlab="Elevation",ylab="Ecological distance")
plot(model_1.splineDat$x[,"Geographic"],model_1.splineDat$y[,"Geographic"],ylim=c(0,0.7),
     lwd=2,type="l",xlab="Geographic distance",ylab="Ecological distance")
plot(x=model_1[["ecological"]],y=model_1[["observed"]],ylim=c(0,1),
     type="n",xlab="Pred. ecological distance",ylab="Obs. compositional dissimilarity")
points(x=model_1[["ecological"]],y=model_1[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_1[["ecological"]]),to=max(model_1[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)

# Variable importance
model_1_Imp<-gdm.varImp(spTable=Data,geo=TRUE,nPerm=100,fullModelOnly=FALSE)
barplot(sort(model_1_Imp[[2]][,1],decreasing=T))
model_1_Imp[[1]]
model_1_Imp[[2]]
model_1_Imp[[3]]
model_1_Imp[[4]]
rm(overlayX,overlayY,model_1_Imp,Env,Data)

# Create variable combinations
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
# Create a dummy raster variable
geo_dummy<-raster(extent(SUDEM),resolution=10,crs=projection(SUDEM))
geo_dummy<-setValues(geo_dummy,1)
geo_dummy<-crop(geo_dummy,SUDEM)
geo_dummy<-mask(geo_dummy,SUDEM)
Env<-stack(Cor_Width,SUDEM)
names(Env)<-c("Cor_Width","Elv")
All<-formatsitepair(bioData=Bio,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                    siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Cor_Width<-formatsitepair(bioData=Bio,bioFormat=2,predData=Cor_Width,sppColumn="Species",dist="jaccard",
                          siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
SUDEM<-formatsitepair(bioData=Bio,bioFormat=2,predData=SUDEM,sppColumn="Species",dist="jaccard",
                      siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
geo_dummy<-formatsitepair(bioData=Bio,bioFormat=2,predData=geo_dummy,sppColumn="Species",dist="jaccard",
                          siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
All<-na.omit(All)
Cor_Width<-na.omit(Cor_Width)
SUDEM<-na.omit(SUDEM)
geo_dummy<-na.omit(geo_dummy)

# Calculate explained deviance
cor_width<-gdm(data=Cor_Width)
cor_width<-round(cor_width[["explained"]],digits=2)
elevation<-gdm(data=SUDEM)
elevation<-round(elevation[["explained"]],digits=2)
geo_dist<-gdm(data=geo_dummy,geo=TRUE)
geo_dist<-round(geo_dist[["explained"]],digits=2)
cor_width_geo_dist<-gdm(data=Cor_Width,geo=TRUE)
cor_width_geo_dist<-round(cor_width_geo_dist[["explained"]],digits=2)
elevation_geo_dist<-gdm(data=SUDEM,geo=TRUE)
elevation_geo_dist<-round(elevation_geo_dist[["explained"]],digits=2)
cor_width_elevation<-gdm(data=All)
cor_width_elevation<-round(cor_width_elevation[["explained"]],digits=2)
cor_width_elevation_geo_dist<-gdm(data=All,geo=TRUE)
cor_width_elevation_geo_dist<-round(cor_width_elevation_geo_dist[["explained"]],digits=2)

# Produce euler plot
fit1<-euler(c("A"=cor_width,"B"=elevation,"C"=geo_dist,
              "A&B"=cor_width_elevation,"A&C"=cor_width_geo_dist,"B&C"=elevation_geo_dist,
              "A&B&C"=cor_width_elevation_geo_dist),shape="ellipse")
plot(fit1,quantities=TRUE,
     labels=c("Edge dist","Elv","Geo dist"),
     fill=c("Grey40","Grey70","Grey100"))
rm(fit1,cor_width_elevation_geo_dist,cor_width_elevation,elevation_geo_dist,cor_width_geo_dist,geo_dist,
   elevation,cor_width,geo_dummy,SUDEM,Cor_Width,Env,Bio,All)


# Load and group biological data
Data<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
groups<-c("Low","Intermediate","Intermediate","Low","Low","High","Low","Low","Intermediate","Intermediate","Low","Low",
          "Intermediate","Low","Low","Low","Low","Intermediate","Low","Low","High","Intermediate","Low","Intermediate",
          "Intermediate","High","Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","Intermediate",
          "Low","High","Intermediate","Low","High","High","Low","Intermediate","Low","Intermediate","High","Intermediate",
          "Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","High","High","Intermediate","High",
          "Low","Intermediate","Low","Intermediate","Low")
DataT<-t(Data)
DataT<-as.data.frame(cbind(DataT,"Groups"=groups))
Low<-subset(DataT,Groups=="Low")
Low$Groups=NULL
Low<-as.data.frame(t(Low))
indx<-sapply(Low,is.factor)
Low[indx]<-lapply(Low[indx],function(x) as.numeric(as.character(x)))
Intermediate<-subset(DataT,Groups=="Intermediate")
Intermediate$Groups=NULL
Intermediate<-as.data.frame(t(Intermediate))
indx<-sapply(Intermediate,is.factor)
Intermediate[indx]<-lapply(Intermediate[indx],function(x) as.numeric(as.character(x)))
High<-subset(DataT,Groups=="High")
High$Groups=NULL
High<-as.data.frame(t(High))
indx<-sapply(High,is.factor)
High[indx]<-lapply(High[indx],function(x) as.numeric(as.character(x)))
rm(DataT,indx,groups,Data)

# Spatial autocorrelation
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
provi<-deldir::deldir(GPS[2:3])
provi.neig<-neig(edges=as.matrix(provi$delsgs[,5:6]))
maf.listw<-spdep::nb2listw(neig2nb(provi.neig))
Low<-Low[,colSums(Low!=0)>3]
multispati.rtest((dudi.pca(Low,scannf=FALSE)),maf.listw) #Spatial autocorrelation
Intermediate<-Intermediate[,colSums(Intermediate!=0)>3]
multispati.rtest((dudi.pca(Intermediate,scannf=FALSE)),maf.listw) #Spatial autocorrelation
High<-High[,colSums(High!=0)>3]
multispati.rtest((dudi.pca(High,scannf=FALSE)),maf.listw) #Spatial autocorrelation
rm(provi,provi.neig,maf.listw)

# Merge GPS with biological data
Site_ID<-GPS$Site_ID
GPS<-SpatialPointsDataFrame(GPS[,2:3],GPS,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
GPS<-spTransform(GPS,CRS="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Proj_Cor<-as.data.frame(GPS@coords)
Low<-as.data.frame(cbind(Site_ID,Proj_Cor,Low))
Intermediate<-as.data.frame(cbind(Site_ID,Proj_Cor,Intermediate))
High<-as.data.frame(cbind(Site_ID,Proj_Cor,High))
Low[Low!=0]<-1
Intermediate[Intermediate!=0]<-1
High[High!=0]<-1
rm(GPS,Proj_Cor,Site_ID)

# Load data environmental data for LOW conservation value species
Cur_Map_Low<-raster("Products/Raster_Products/Cur_Map_Low.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist.tif")
Map_ranger<-raster("Products/Raster_Products/Map_ranger.tif")
Max_NBR<-stack("Products/Raster_Products/Max_NBR.tif")
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
Env<-stack(Sent_NDVI,Aspect,Max_NBR,Cur_Map_Low,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
rm(Sent_NDVI,Aspect,Max_NBR,Cur_Map_Low,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
names(Env)
names(Env)<-c("NDVI","Aspect","Fire_Hist","Func_Conn","Land_Use","Drain_Dist","Cor_Width","Elv")

# Prepare data for gdm LOW
Data<-formatsitepair(bioData=Low,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# Variable importance 11.5x5
model_1_Imp<-gdm.varImp(spTable=Data,geo=TRUE,nPerm=100,fullModelOnly=FALSE)
barplot(sort(model_1_Imp[[2]][,1],decreasing=T))
model_1_Imp[[1]]
model_1_Imp[[2]]
model_1_Imp[[3]]
model_1_Imp[[4]]

# Select only significant variables
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Env<-stack(Cor_Width,SUDEM)
rm(Cor_Width,SUDEM)
names(Env)
names(Env)<-c("Cor_Width","Elv")
Data<-formatsitepair(bioData=Low,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# GDM 7x7
model_2<-gdm(data=Data,geo=TRUE)
summary.gdm(model_2)
plot(model_2,plot.layout=c(3,2))
model_2.splineDat<-isplineExtract(model_2)
par(mfrow=c(2,2))
plot(model_2.splineDat$x[,"Cor_Width"],model_2.splineDat$y[,"Cor_Width"],ylim=c(0,0.8),
     lwd=2,type="l",xlab="Edge distance",ylab="Ecological distance")
plot(model_2.splineDat$x[,"Elv"],model_2.splineDat$y[,"Elv"],ylim=c(0,0.7),
     lwd=2,type="l",xlab="Elevation",ylab="Ecological distance")
plot(model_2.splineDat$x[,"Geographic"],model_2.splineDat$y[,"Geographic"],ylim=c(0,0.8),
     lwd=2,type="l",xlab="Geographic distance",ylab="Ecological distance")
plot(x=model_2[["ecological"]],y=model_2[["observed"]],ylim=c(0,1),
     type="n",xlab="Pred. ecological distance",ylab="Obs. compositional dissimilarity")
points(x=model_2[["ecological"]],y=model_2[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_2[["ecological"]]),to=max(model_2[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)

# Variable importance
model_1_Imp<-gdm.varImp(spTable=Data,geo=TRUE,nPerm=100,fullModelOnly=FALSE)
barplot(sort(model_1_Imp[[2]][,1],decreasing=T))
model_1_Imp[[1]]
model_1_Imp[[2]]
model_1_Imp[[3]]
model_1_Imp[[4]]
rm(overlayX,overlayY,model_1_Imp,Env,Data)

# Create variable combinations
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
# Create a dummy raster variable
geo_dummy<-raster(extent(SUDEM),resolution=10,crs=projection(SUDEM))
geo_dummy<-setValues(geo_dummy,1)
geo_dummy<-crop(geo_dummy,SUDEM)
geo_dummy<-mask(geo_dummy,SUDEM)
Env<-stack(Cor_Width,SUDEM)
names(Env)<-c("Cor_Width","Elv")
All<-formatsitepair(bioData=Low,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                    siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Cor_Width<-formatsitepair(bioData=Low,bioFormat=2,predData=Cor_Width,sppColumn="Species",dist="jaccard",
                          siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
SUDEM<-formatsitepair(bioData=Low,bioFormat=2,predData=SUDEM,sppColumn="Species",dist="jaccard",
                      siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
geo_dummy<-formatsitepair(bioData=Low,bioFormat=2,predData=geo_dummy,sppColumn="Species",dist="jaccard",
                          siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
All<-na.omit(All)
Cor_Width<-na.omit(Cor_Width)
SUDEM<-na.omit(SUDEM)
geo_dummy<-na.omit(geo_dummy)

# Calculate explained deviance
cor_width<-gdm(data=Cor_Width)
cor_width<-round(cor_width[["explained"]],digits=2)
elevation<-gdm(data=SUDEM)
elevation<-round(elevation[["explained"]],digits=2)
geo_dist<-gdm(data=geo_dummy,geo=TRUE)
geo_dist<-round(geo_dist[["explained"]],digits=2)
cor_width_geo_dist<-gdm(data=Cor_Width,geo=TRUE)
cor_width_geo_dist<-round(cor_width_geo_dist[["explained"]],digits=2)
elevation_geo_dist<-gdm(data=SUDEM,geo=TRUE)
elevation_geo_dist<-round(elevation_geo_dist[["explained"]],digits=2)
cor_width_elevation<-gdm(data=All)
cor_width_elevation<-round(cor_width_elevation[["explained"]],digits=2)
cor_width_elevation_geo_dist<-gdm(data=All,geo=TRUE)
cor_width_elevation_geo_dist<-round(cor_width_elevation_geo_dist[["explained"]],digits=2)

# Produce euler plot 5.5x5.5
fit1<-euler(c("A"=cor_width,"B"=elevation,"C"=geo_dist,
              "A&B"=cor_width_elevation,"A&C"=cor_width_geo_dist,"B&C"=elevation_geo_dist,
              "A&B&C"=cor_width_elevation_geo_dist),shape="ellipse")
plot(fit1,quantities=TRUE,
     labels=c("Dist edge","Elv","Geo dist"),
     fill=c("Grey40","Grey70","Grey100"))
rm(fit1,cor_width_elevation_geo_dist,cor_width_elevation,elevation_geo_dist,cor_width_geo_dist,geo_dist,
   elevation,cor_width,geo_dummy,SUDEM,Cor_Width,Env,Low,All)

# Load data environmental data for INTERMEDIATE conservation value species
Cur_Map_Intermediate<-raster("Products/Raster_Products/Cur_Map_Intermediate.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist.tif")
Map_ranger<-raster("Products/Raster_Products/Map_ranger.tif")
Max_NBR<-stack("Products/Raster_Products/Max_NBR.tif")
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
Env<-stack(Sent_NDVI,Aspect,Max_NBR,Cur_Map_Intermediate,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
rm(Sent_NDVI,Aspect,Max_NBR,Cur_Map_Intermediate,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
names(Env)
names(Env)<-c("NDVI","Aspect","Fire_Hist","Func_Conn","Land_Use","Drain_Dist","Cor_Width","Elv")

# Prepare data for gdm Intermediate
Data<-formatsitepair(bioData=Intermediate,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# Variable importance 11.5x5
model_1_Imp<-gdm.varImp(spTable=Data,geo=TRUE,nPerm=100,fullModelOnly=FALSE)
barplot(sort(model_1_Imp[[2]][,1],decreasing=T))
model_1_Imp[[1]]
model_1_Imp[[2]]
model_1_Imp[[3]]
model_1_Imp[[4]]

# Select only significant variables
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Env<-SUDEM
rm(SUDEM)
Data<-formatsitepair(bioData=Intermediate,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# GDM 7x3.5
model_3<-gdm(data=Data,geo=TRUE)
summary.gdm(model_3)
plot(model_3,plot.layout=c(2,2))
model_3.splineDat<-isplineExtract(model_3)
par(mfrow=c(1,2))
plot(model_3.splineDat$x[,"extract(predData, bioData$cellName)"],model_3.splineDat$y[,"extract(predData, bioData$cellName)"],ylim=c(0,0.4),
     lwd=2,type="l",xlab="Elv",ylab="Ecological distance")
plot(x=model_3[["ecological"]],y=model_3[["observed"]],ylim=c(0,1),
     type="n",xlab="Pred. ecological distance",ylab="Obs. compositional dissimilarity")
points(x=model_3[["ecological"]],y=model_3[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_3[["ecological"]]),to=max(model_3[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)
rm(overlayX,overlayY,model_1_Imp,Env,Data)

# Create variable combinations
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
# Create a dummy raster variable
geo_dummy<-raster(extent(SUDEM),resolution=10,crs=projection(SUDEM))
geo_dummy<-setValues(geo_dummy,1)
geo_dummy<-crop(geo_dummy,SUDEM)
geo_dummy<-mask(geo_dummy,SUDEM)
Env<-SUDEM
names(Env)<-c("Elv")
All<-formatsitepair(bioData=Intermediate,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                    siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
SUDEM<-formatsitepair(bioData=Intermediate,bioFormat=2,predData=SUDEM,sppColumn="Species",dist="jaccard",
                      siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
geo_dummy<-formatsitepair(bioData=Intermediate,bioFormat=2,predData=geo_dummy,sppColumn="Species",dist="jaccard",
                          siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
All<-na.omit(All)
SUDEM<-na.omit(SUDEM)
geo_dummy<-na.omit(geo_dummy)

# Calculate explained deviance
elevation<-gdm(data=SUDEM)
elevation<-round(elevation[["explained"]],digits=2)
geo_dist<-gdm(data=geo_dummy,geo=TRUE)
geo_dist<-round(geo_dist[["explained"]],digits=2)
elevation_geo_dist<-gdm(data=SUDEM,geo=TRUE)
elevation_geo_dist<-round(elevation_geo_dist[["explained"]],digits=2)

# Produce euler plot 5.5x5.5
fit1<-euler(c("A"=elevation,"B"=geo_dist,
              "A&B"=elevation_geo_dist),shape="ellipse")
plot(fit1,quantities=TRUE,
     labels=c("Elv","Geo dist"),
     fill=c("Grey40","Grey70","Grey100"))
rm(fit1,elevation_geo_dist,geo_dist,elevation,geo_dummy,SUDEM,Cor_Width,Env,Intermediate,All)

# Load data environmental data for HIGH conservation value species
Cur_Map_High<-raster("Products/Raster_Products/Cur_Map_High.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist.tif")
Map_ranger<-raster("Products/Raster_Products/Map_ranger.tif")
Max_NBR<-stack("Products/Raster_Products/Max_NBR.tif")
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
Env<-stack(Sent_NDVI,Aspect,Max_NBR,Cur_Map_High,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
rm(Sent_NDVI,Aspect,Max_NBR,Cur_Map_High,Map_ranger,Drain_Dist,Cor_Width,SUDEM)
names(Env)
names(Env)<-c("NDVI","Aspect","Fire_Hist","Func_Conn","Land_Use","Drain_Dist","Cor_Width","Elv")

# Prepare data for gdm High
Data<-formatsitepair(bioData=High,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# Variable importance 11.5x5
model_1_Imp<-gdm.varImp(spTable=Data,geo=TRUE,nPerm=100,fullModelOnly=FALSE)
barplot(sort(model_1_Imp[[2]][,1],decreasing=T))
model_1_Imp[[1]]
model_1_Imp[[2]]
model_1_Imp[[3]]
model_1_Imp[[4]]

# Select only significant variables
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
Env<-stack(Aspect,SUDEM)
rm(Aspect,SUDEM)
names(Env)
names(Env)<-c("Aspect","Elv")
Data<-formatsitepair(bioData=High,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# GDM 7x7
model_4<-gdm(data=Data,geo=TRUE)
summary.gdm(model_4)
plot(model_4,plot.layout=c(2,2))
model_4.splineDat<-isplineExtract(model_4)
par(mfrow=c(2,2))
plot(model_4.splineDat$x[,"Aspect"],model_4.splineDat$y[,"Aspect"],ylim=c(0,1.8),
     lwd=2,type="l",xlab="Aspect",ylab="Ecological distance")
plot(model_4.splineDat$x[,"Elv"],model_4.splineDat$y[,"Elv"],ylim=c(0,1.8),
     lwd=2,type="l",xlab="Elevation",ylab="Ecological distance")
plot(x=model_4[["ecological"]],y=model_4[["observed"]],ylim=c(0,1),
     type="n",xlab="Pred. ecological distance",ylab="Obs. compositional dissimilarity")
points(x=model_4[["ecological"]],y=model_4[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_4[["ecological"]]),to=max(model_4[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)

# Variable importance
model_1_Imp<-gdm.varImp(spTable=Data,geo=TRUE,nPerm=100,fullModelOnly=FALSE)
barplot(sort(model_1_Imp[[2]][,1],decreasing=T))
model_1_Imp[[1]]
model_1_Imp[[2]]
model_1_Imp[[3]]
model_1_Imp[[4]]
rm(overlayX,overlayY,model_1_Imp,Env,Data)

# Create variable combinations
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
# Create a dummy raster variable
geo_dummy<-raster(extent(SUDEM),resolution=10,crs=projection(SUDEM))
geo_dummy<-setValues(geo_dummy,1)
geo_dummy<-crop(geo_dummy,SUDEM)
geo_dummy<-mask(geo_dummy,SUDEM)
Env<-stack(Aspect,SUDEM)
names(Env)<-c("Aspect","Elv")
All<-formatsitepair(bioData=High,bioFormat=2,predData=Env,sppColumn="Species",dist="jaccard",
                    siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Aspect<-formatsitepair(bioData=High,bioFormat=2,predData=Aspect,sppColumn="Species",dist="jaccard",
                       siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
SUDEM<-formatsitepair(bioData=High,bioFormat=2,predData=SUDEM,sppColumn="Species",dist="jaccard",
                      siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
geo_dummy<-formatsitepair(bioData=High,bioFormat=2,predData=geo_dummy,sppColumn="Species",dist="jaccard",
                          siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
All<-na.omit(All)
Aspect<-na.omit(Aspect)
SUDEM<-na.omit(SUDEM)
geo_dummy<-na.omit(geo_dummy)

# Calculate explained deviance
aspect<-gdm(data=Aspect)
aspect<-round(aspect[["explained"]],digits=2)
elevation<-gdm(data=SUDEM)
elevation<-round(elevation[["explained"]],digits=2)
geo_dist<-gdm(data=geo_dummy,geo=TRUE)
geo_dist<-round(geo_dist[["explained"]],digits=2)
aspect_geo_dist<-gdm(data=Aspect,geo=TRUE)
aspect_geo_dist<-round(aspect_geo_dist[["explained"]],digits=2)
elevation_geo_dist<-gdm(data=SUDEM,geo=TRUE)
elevation_geo_dist<-round(elevation_geo_dist[["explained"]],digits=2)
aspect_elevation<-gdm(data=All)
aspect_elevation<-round(aspect_elevation[["explained"]],digits=2)
aspect_elevation_geo_dist<-gdm(data=All,geo=TRUE)
aspect_elevation_geo_dist<-round(aspect_elevation_geo_dist[["explained"]],digits=2)

# Produce euler plot 5.5x5.5
fit1<-euler(c("A"=aspect,"B"=elevation,"C"=geo_dist,
              "A&B"=aspect_elevation,"A&C"=aspect_geo_dist,"B&C"=elevation_geo_dist,
              "A&B&C"=aspect_elevation_geo_dist),shape="ellipse")
plot(fit1,quantities=TRUE,
     labels=c("Aspect","Elv","Geo dist"),
     fill=c("Grey40","Grey70","Grey100"))
rm(fit1,aspect_elevation_geo_dist,aspect_elevation,elevation_geo_dist,aspect_geo_dist,geo_dist,
   elevation,aspect,geo_dummy,SUDEM,Aspect,Env,High,All)


# Plot GDM Variable importance together
par(mar=c(3.3,3.3,1,1))
par(oma=c(0,0,0,0))
par(mfrow=c(7,2))
#OVERALL
plot(model_1.splineDat$x[,"Cor_Width"],model_1.splineDat$y[,"Cor_Width"],ylim=c(0,1),
     lwd=2,type="l",xlab="",ylab="")
mtext("a.1)",adj=0,line=0,cex=0.7)
mtext("Edge distance",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(model_1.splineDat$x[,"Elv"],model_1.splineDat$y[,"Elv"],ylim=c(0,1),
     lwd=2,type="l",xlab="",ylab="")
mtext("a.2)",adj=0,line=0,cex=0.7)
mtext("Elevation",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(model_1.splineDat$x[,"Geographic"],model_1.splineDat$y[,"Geographic"],ylim=c(0,1),
     lwd=2,type="l",xlab="",ylab="")
mtext("a.3)",adj=0,line=0,cex=0.7)
mtext("Geographic distance",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(x=model_1[["ecological"]],y=model_1[["observed"]],ylim=c(0,1),
     type="n",xlab="",ylab="")
mtext("a.4)",adj=0,line=0,cex=0.7)
mtext("Pred ecol dist",side=1,line=2,cex=0.7)
mtext("Obs comp diss",side=2,line=2,cex=0.7)
points(x=model_1[["ecological"]],y=model_1[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_1[["ecological"]]),to=max(model_1[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)
#LOW
plot(model_2.splineDat$x[,"Cor_Width"],model_2.splineDat$y[,"Cor_Width"],ylim=c(0,1),
     lwd=2,type="l",xlab="",ylab="")
mtext("b.1)",adj=0,line=0,cex=0.7)
mtext("Edge distance",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(model_2.splineDat$x[,"Elv"],model_2.splineDat$y[,"Elv"],ylim=c(0,1),
     lwd=2,type="l",xlab="",ylab="")
mtext("b.2)",adj=0,line=0,cex=0.7)
mtext("Elevation",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(model_2.splineDat$x[,"Geographic"],model_2.splineDat$y[,"Geographic"],ylim=c(0,1),
     lwd=2,type="l",xlab="",ylab="")
mtext("b.3)",adj=0,line=0,cex=0.7)
mtext("Geographic distance",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(x=model_2[["ecological"]],y=model_2[["observed"]],ylim=c(0,1),
     type="n",xlab="",ylab="")
mtext("b.4)",adj=0,line=0,cex=0.7)
mtext("Pred ecol dist",side=1,line=2,cex=0.7)
mtext("Obs comp diss",side=2,line=2,cex=0.7)
points(x=model_2[["ecological"]],y=model_2[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_2[["ecological"]]),to=max(model_2[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)
#INTERMEDIATE
plot(model_3.splineDat$x[,"extract(predData, bioData$cellName)"],model_3.splineDat$y[,"extract(predData, bioData$cellName)"],ylim=c(0,1),
     lwd=2,type="l",xlab="",ylab="")
mtext("c.1)",adj=0,line=0,cex=0.7)
mtext("Elevation",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(x=model_3[["ecological"]],y=model_3[["observed"]],ylim=c(0,1),
     type="n",xlab="",ylab="")
mtext("c.2)",adj=0,line=0,cex=0.7)
mtext("Pred ecol dist",side=1,line=2,cex=0.7)
mtext("Obs comp diss",side=2,line=2,cex=0.7)
points(x=model_3[["ecological"]],y=model_3[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_3[["ecological"]]),to=max(model_3[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)
#HIGH
plot(model_4.splineDat$x[,"Aspect"],model_4.splineDat$y[,"Aspect"],ylim=c(0,2),
     lwd=2,type="l",xlab="",ylab="")
mtext("d.1)",adj=0,line=0,cex=0.7)
mtext("Aspect",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(model_4.splineDat$x[,"Elv"],model_4.splineDat$y[,"Elv"],ylim=c(0,2),
     lwd=2,type="l",xlab="",ylab="")
mtext("d.2)",adj=0,line=0,cex=0.7)
mtext("Elevation",side=1,line=2,cex=0.7)
mtext("Ecol dist",side=2,line=2,cex=0.7)
plot(x=model_4[["ecological"]],y=model_4[["observed"]],ylim=c(0,1),
     type="n",xlab="",ylab="")
mtext("d.3)",adj=0,line=0,cex=0.7)
mtext("Pred ecol dist",side=1,line=2,cex=0.7)
mtext("Obs comp diss",side=2,line=2,cex=0.7)
points(x=model_4[["ecological"]],y=model_4[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_4[["ecological"]]),to=max(model_4[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)
par(mar=c(5.1,4.1,4.1,2.1))
rm(overlayX,overlayY,model_1,model_1.splineDat,model_2,model_2.splineDat,model_3,model_3.splineDat,model_4,model_4.splineDat)

############### Drivers of nestedness ###############
#Load raster data
Max_NBR<-raster("Products/Raster_Products/Max_NBR.tif")
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist.tif")
SUDEM_Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
Func_Con<-raster("Products/Raster_Products/Cur_Map.tif")
Func_Con_Low<-raster("Products/Raster_Products/Cur_Map_Low.tif")
Func_Con_Intermediate<-raster("Products/Raster_Products/Cur_Map_Intermediate.tif")
Func_Con_High<-raster("Products/Raster_Products/Cur_Map_High.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")

#Extract values over GPS
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Plantation<-GPS$Plantation
GPS<-SpatialPointsDataFrame(GPS[,2:3],GPS,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
GPS<-spTransform(GPS,CRS="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Max_NBR<-raster::extract(Max_NBR,GPS,fun=mean,buffer=50)
Sent_NDVI<-raster::extract(Sent_NDVI,GPS,fun=mean,buffer=50)
Drain_Dist<-raster::extract(Drain_Dist,GPS,fun=mean,buffer=50)
SUDEM_Aspect<-raster::extract(SUDEM_Aspect,GPS,fun=mean,buffer=50)
Cor_Width<-raster::extract(Cor_Width,GPS,fun=mean,buffer=50)
Func_Con<-raster::extract(Func_Con,GPS,fun=mean,buffer=50)
Func_Con_Low<-raster::extract(Func_Con_Low,GPS,fun=mean,buffer=50)
Func_Con_Intermediate<-raster::extract(Func_Con_Intermediate,GPS,fun=mean,buffer=50)
Func_Con_High<-raster::extract(Func_Con_High,GPS,fun=mean,buffer=50)
SUDEM<-raster::extract(SUDEM,GPS,fun=mean,buffer=50)
df<-as.data.frame(cbind("Fire"=Max_NBR,"NDVI"=Sent_NDVI,"Drain"=Drain_Dist,"Aspect"=SUDEM_Aspect,
                        "Edge_Dist"=Cor_Width,"Func_Con"=Func_Con,"Func_Con_Low"=Func_Con_Low,
                        "Func_Con_Intermediate"=Func_Con_Intermediate,"Func_Con_High"=Func_Con_High,"Elv"=SUDEM))
rm(Max_NBR,Sent_NDVI,Drain_Dist,SUDEM_Aspect,Cor_Width,Func_Con,Func_Con_Low,Func_Con_Intermediate,Func_Con_High,SUDEM)

# Load and calculating response variables
Data<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
groups<-c("Low","Intermediate","Intermediate","Low","Low","High","Low","Low","Intermediate","Intermediate","Low","Low",
          "Intermediate","Low","Low","Low","Low","Intermediate","Low","Low","High","Intermediate","Low","Intermediate",
          "Intermediate","High","Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","Intermediate",
          "Low","High","Intermediate","Low","High","High","Low","Intermediate","Low","Intermediate","High","Intermediate",
          "Intermediate","Intermediate","Intermediate","Intermediate","Intermediate","High","High","Intermediate","High",
          "Low","Intermediate","Low","Intermediate","Low")
DataT<-t(Data)
DataT<-as.data.frame(cbind(DataT,"Groups"=groups))
Low<-subset(DataT,Groups=="Low")
Low$Groups=NULL
Low<-as.data.frame(t(Low))
indx<-sapply(Low,is.factor)
Low[indx]<-lapply(Low[indx],function(x) as.numeric(as.character(x)))
Low_Richness<-rowSums(Low!=0)
Low_exShan<-exp(diversity(Low,index="shannon"))
Intermediate<-subset(DataT,Groups=="Intermediate")
Intermediate$Groups=NULL
Intermediate<-as.data.frame(t(Intermediate))
indx<-sapply(Intermediate,is.factor)
Intermediate[indx]<-lapply(Intermediate[indx],function(x) as.numeric(as.character(x)))
Intermediate_Richness<-rowSums(Intermediate!=0)
Intermediate_exShan<-exp(diversity(Intermediate,index="shannon"))
High<-subset(DataT,Groups=="High")
High$Groups=NULL
High<-as.data.frame(t(High))
indx<-sapply(High,is.factor)
High[indx]<-lapply(High[indx],function(x) as.numeric(as.character(x)))
High_Richness<-rowSums(High!=0)
High_exShan<-exp(diversity(High,index="shannon"))
Richness<-rowSums(Data!=0)
exShan<-exp(diversity(Data,index="shannon"))
df<-as.data.frame(cbind(Plantation,"Richness"=Richness,"Low_Richness"=Low_Richness,"Intermediate_Richness"=Intermediate_Richness,"High_Richness"=High_Richness,
                        "exShannon"=exShan,"Low_exShan"=Low_exShan,"Intermediate_exShan"=Intermediate_exShan,"High_exShan"=High_exShan,df))
rm(DataT,Richness,exShan,Low_Richness,Low_exShan,Intermediate_Richness,
   Intermediate_exShan,High_Richness,High_exShan,indx,groups,High,Low,Intermediate)

#Plot overall curve
Curve<-specaccum(Data,method="rarefaction")
plot(Curve,main="Species rarefaction",ylab="Rarefaction")
#Plot per plantation
hopperPlant<-as.data.frame(cbind(Data,"Plantation"=Plantation))
#Subset
Estate1<-subset(hopperPlant,Plantation=="Estate1")
Estate1$Plantation=NULL
Estate2<-subset(hopperPlant,Plantation=="Estate2")
Estate2$Plantation=NULL
Estate3<-subset(hopperPlant,Plantation=="Estate3")
Estate3$Plantation=NULL
Estate4<-subset(hopperPlant,Plantation=="Estate4")
Estate4$Plantation=NULL
#Plotting
Estate1Curve<-specaccum(Estate1,method="rarefaction")
Estate2Curve<-specaccum(Estate2,method="rarefaction")
Estate3Curve<-specaccum(Estate3,method="rarefaction")
Estate4Curve<-specaccum(Estate4,method="rarefaction")
plot(Estate1Curve,main="Species rarefaction",ylab="Rarefaction")
plot(Estate2Curve,add=TRUE,col='red')
plot(Estate3Curve,add=TRUE,col='blue')
plot(Estate4Curve,add=TRUE,col='green')
legend(x="bottomright",legend=c("Estate 1","Estate 2","Estate 3","Estate 4"),fill=c("black","red","blue","green"))
#Clean environment
rm(Curve,Plantation,hopperPlant,Estate1,Estate2,Estate3,Estate4,Estate1Curve,Estate2Curve,Estate3Curve,Estate4Curve)

# Normality testing
shapiro.test(df$Richness) 
hist(df$Richness)     # Fine
shapiro.test(df$exShan)
hist(df$exShan)       # Fine
shapiro.test(df$Low_Richness)
hist(df$Low_Richness) # Fine
shapiro.test(df$Low_exShan)
hist(df$Low_exShan)   # Fine
shapiro.test(df$Intermediate_Richness)
hist(df$Intermediate_Richness) # Fine
shapiro.test(df$Intermediate_exShan)
hist(df$Intermediate_exShan)   # Fine
shapiro.test(df$High_Richness)
hist(df$High_Richness)
summary(goodfit(df$High_Richness,type="poisson",method="MinChisq"))
shapiro.test(df$High_exShan)
hist(df$High_exShan)
summary(goodfit(df$High_exShan,type="poisson",method="MinChisq"))

# Spatial autocorrelation
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Data.dist.inv<-1/(as.matrix(dist(cbind(GPS$Long_X,GPS$Lat_Y))))
Data.dist.inv[is.infinite(Data.dist.inv)]<-0
Moran.I(df$Richness,Data.dist.inv)
Moran.I(df$exShannon,Data.dist.inv)
Moran.I(df$Low_Richness,Data.dist.inv)
Moran.I(df$Low_exShan,Data.dist.inv)
Moran.I(df$Intermediate_Richness,Data.dist.inv)
Moran.I(df$Intermediate_exShan,Data.dist.inv)
Moran.I(df$High_Richness,Data.dist.inv)
Moran.I(df$High_exShan,Data.dist.inv)
rm(GPS,Data.dist.inv)

#Check for correlation between features
features<-df[,10:19]
ggcorrplot(cor(features,method="spearman"),
           type="lower",lab=TRUE,lab_size=3,tl.cex= 10)
df$NDVI=NULL
df$Func_Con=NULL
df$Func_Con_Low=NULL
df$Func_Con_Intermediate=NULL
df$Func_Con_High=NULL
rm(features)

# Standerdise variables
df_stand<-standardize(df[,2:14])
df_stand<-as.data.frame(cbind(df_stand,"Plantation"=df$Plantation))

# Model Richness
Model_Rich<-lmer(Richness~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),data=df_stand)
vif(Model_Rich)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.359000  1.029662  1.090805  1.226544  1.486527 
'''
options(na.action="na.fail")
Model_Rich_Dredge<-dredge(Model_Rich,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Rich_Subset<-subset(Model_Rich_Dredge,delta<2)
Model_Rich_Ave<-model.avg(Model_Rich_Subset)
importance(Model_Rich_Ave)
confint(Model_Rich_Ave)
summary(Model_Rich_Ave)
rm(Model_Rich_Dredge,Model_Rich_Subset,Model_Rich,Model_Rich_Ave)
# Plot 5.5x5.5
ggplot(data=df,aes(x=Edge_Dist,y=Richness))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera species richness\n",x="\nDistance to edge (m)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())

# Model exShannon
Model_exShan<-lmer(exShannon~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),data=df_stand)
vif(Model_exShan)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.532133  1.024122  1.093195  1.309799  1.727296 
'''
options(na.action="na.fail")
Model_exShan_Dredge<-dredge(Model_exShan,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_exShan_Subset<-subset(Model_exShan_Dredge,delta<2)
Model_exShan_Top<-get.models(Model_exShan_Subset,1)[[1]]
confint(Model_exShan_Top)
summary(Model_exShan_Top)
rm(Model_exShan_Dredge,Model_exShan_Subset,Model_exShan,Model_exShan_Top)
# Plot
ggplot(data=df,aes(x=Edge_Dist,y=exShannon))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Caelifera exponent Shannon diversity\n",x="\nDistance to edge (m)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())

# Model LOW Richness
Model_Low_Rich<-lmer(Low_Richness~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),data=df_stand)
vif(Model_Low_Rich)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.363694  1.029256  1.090731  1.228733  1.493079 
'''
options(na.action="na.fail")
Model_Low_Rich_Dredge<-dredge(Model_Low_Rich,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Low_Rich_Subset<-subset(Model_Low_Rich_Dredge,delta<2)
Model_Low_Rich_Ave<-model.avg(Model_Low_Rich_Subset)
importance(Model_Low_Rich_Ave)
confint(Model_Low_Rich_Ave)
summary(Model_Low_Rich_Ave)
rm(Model_Low_Rich_Dredge,Model_Low_Rich_Subset,Model_Low_Rich,Model_Low_Rich_Ave)
# Plot
ggplot(data=df,aes(x=Edge_Dist,y=Low_Richness))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Low conservation value species richness\n",x="\nDistance to edge (m)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())

# Model LOW exShannon
Model_Low_exShannon<-lmer(Low_exShan~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),data=df_stand)
vif(Model_Low_exShannon)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.479321  1.024680  1.091709  1.283887  1.655060 
'''
options(na.action="na.fail")
Model_Low_exShannon_Dredge<-dredge(Model_Low_exShannon,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Low_exShannon_Subset<-subset(Model_Low_exShannon_Dredge,delta<2)
Model_Low_exShannon_Top<-get.models(Model_Low_exShannon_Subset,1)[[1]]
confint(Model_Low_exShannon_Top)
summary(Model_Low_exShannon_Top)
rm(Model_Low_exShannon_Dredge,Model_Low_exShannon_Subset,Model_Low_exShannon,Model_Low_exShannon_Top)
# Plot
ggplot(data=df,aes(x=Edge_Dist,y=Low_exShan))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Low conservation value Shannon diversity\n",x="\nDistance to edge (m)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())

# Model INTERMEDIATE Richness
Model_Intermediate_Rich<-lmer(Intermediate_Richness~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),data=df_stand)
vif(Model_Intermediate_Rich)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.767695  1.024012  1.103625  1.429411  2.028319
'''
options(na.action="na.fail")
Model_Intermediate_Rich_Dredge<-dredge(Model_Intermediate_Rich,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Intermediate_Rich_Subset<-subset(Model_Intermediate_Rich_Dredge,delta<2)
Model_Intermediate_Rich_Ave<-model.avg(Model_Intermediate_Rich_Subset)
importance(Model_Intermediate_Rich_Ave)
confint(Model_Intermediate_Rich_Ave)
summary(Model_Intermediate_Rich_Ave)
rm(Model_Intermediate_Rich_Dredge,Model_Intermediate_Rich_Subset,Model_Intermediate_Rich,Model_Intermediate_Rich_Ave)

# Model INTERMEDIATE exShannon
Model_Intermediate_exShannon<-lmer(Intermediate_exShan~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),data=df_stand)
vif(Model_Intermediate_exShannon)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.835046  1.024149  1.107434  1.464527  2.108105 
'''
options(na.action="na.fail")
Model_Intermediate_exShannon_Dredge<-dredge(Model_Intermediate_exShannon,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_Intermediate_exShannon_Subset<-subset(Model_Intermediate_exShannon_Dredge,delta<2)
Model_Intermediate_exShannon_Ave<-model.avg(Model_Intermediate_exShannon_Subset)
importance(Model_Intermediate_exShannon_Ave)
confint(Model_Intermediate_exShannon_Ave)
summary(Model_Intermediate_exShannon_Ave)
rm(Model_Intermediate_exShannon_Dredge,Model_Intermediate_exShannon_Subset,Model_Intermediate_exShannon,Model_Intermediate_exShannon_Ave)
# Plot
ggplot(data=df,aes(x=Elv,y=Intermediate_exShan))+
  geom_point()+
  geom_smooth(method="lm",colour="black")+
  labs(y="Intermediate conservation value Shannon diversity\n",x="\nElevation (m)")+
  theme(axis.title=element_text(size=14),
        axis.text=element_text(size=13,colour="black"),
        axis.ticks=element_line(size=0.8,colour="black"),
        axis.line=element_line(size=0.8,colour="black"),
        panel.grid=element_blank(),
        panel.background=element_blank())

# Model HIGH Richness
Model_High_Rich<-glmmadmb(High_Richness~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),
                          data=df,zeroInflation=TRUE,family="poisson")
vif(Model_High_Rich)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.016142  1.008139  1.005499  1.029271  1.019168 
'''
options(na.action="na.fail")
Model_High_Rich_Dredge<-dredge(Model_High_Rich,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_High_Rich_Subset<-subset(Model_High_Rich_Dredge,delta<2)
Model_High_Rich_Ave<-model.avg(Model_High_Rich_Subset)
importance(Model_High_Rich_Ave)
confint(Model_High_Rich_Ave)
summary(Model_High_Rich_Ave)
rm(Model_High_Rich_Dredge,Model_High_Rich_Subset,Model_High_Rich,Model_High_Rich_Ave)

# Model HIGH exShannon
Model_High_exShannon<-glmmadmb(High_exShan~Fire+Drain+Aspect+Edge_Dist+Elv+(1|Plantation),
                          data=df,zeroInflation=TRUE,family="gamma")
vif(Model_High_exShannon)
'''
     Fire     Drain    Aspect Edge_Dist       Elv 
 1.000000  1.000396  1.000949  1.000373  1.000550 
'''
options(na.action="na.fail")
Model_High_exShannon_Dredge<-dredge(Model_High_exShannon,evaluate=TRUE,rank=AICc)
options(na.action="na.omit")
Model_High_exShannon_Subset<-subset(Model_High_exShannon_Dredge,delta<2)

Model_High_exShannon_Ave<-model.avg(Model_High_exShannon_Subset)
importance(Model_High_exShannon_Ave)
confint(Model_High_exShannon_Ave)
summary(Model_High_exShannon_Ave)
rm(Model_High_exShannon_Dredge,Model_High_exShannon_Subset,Model_High_exShannon,Model_High_exShannon_Ave)
