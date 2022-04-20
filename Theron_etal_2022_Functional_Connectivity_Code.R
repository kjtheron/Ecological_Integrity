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
library(lattice)
library(ggcorrplot)

# Multivariate analysis
library(vegan)
library(mvabund)
library(SSDM)
library(gdm)
library(ResistanceGA)
library(leastcostpath)
library(gdistance)
library(robustHD)

# Spatial autocorrelation
library(ape)
library(ade4)

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

# Load and calculate response variable
Data<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")

# Merge GPS data with multivariate abundance data
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Plantation<-GPS$Plantation
GPS<-SpatialPointsDataFrame(GPS[,2:3],GPS,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
GPS<-spTransform(GPS,CRS="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Proj_Cor<-as.data.frame(GPS@coords)

# Removing species with 3 or less occurrences
Data<-Data[,colSums(Data!=0)>3]
Data<-as.data.frame(cbind(Plantation,Proj_Cor,Data))
rm(GPS,Proj_Cor,Plantation)

# Transform multivariate abundance data to occurrence matrix
Data<-Data %>%
  dplyr::select(-Plantation) %>%
  pivot_longer(cols=-(1:2),values_to="Presence") %>%
  mutate(Presence=pmin(Presence,1)) %>%
  arrange(name) %>%
  rename(Species=name)
write.csv(Data,"Products/Excel_Sheets_Products/Occurrence_Data.csv",row.names=FALSE)
rm(Data)

#Create design variable (corridor width)
Ranger_map<-raster("Products/Raster_Products/Map_ranger.tif")
Ranger_map[Ranger_map==1]<-1
Ranger_map[Ranger_map==2]<-NA
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

# Restrict SSDM to grassland corridors
Ranger_map<-raster("Products/Raster_Products/Map_ranger.tif")
Ranger_map[Ranger_map==1]<-NA
Ranger_map[Ranger_map==4]<-NA
Ranger_map[Ranger_map==6]<-NA
Ranger_map[Ranger_map==7]<-NA
NDVI<-raster("Products/Raster_Products/Sent_NDVI.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist.tif")
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width.tif")
Max_NBR<-raster("Products/Raster_Products/Max_NBR.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM.tif")
NDVI<-crop(NDVI,Ranger_map)
NDVI<-mask(NDVI,Ranger_map)
Drain_Dist<-crop(Drain_Dist,Ranger_map)
Drain_Dist<-mask(Drain_Dist,Ranger_map)
Aspect<-crop(Aspect,Ranger_map)
Aspect<-mask(Aspect,Ranger_map)
Cor_Width<-crop(Cor_Width,Ranger_map)
Cor_Width<-mask(Cor_Width,Ranger_map)
Max_NBR<-crop(Max_NBR,Ranger_map)
Max_NBR<-mask(Max_NBR,Ranger_map)
SUDEM<-crop(SUDEM,Ranger_map)
SUDEM<-mask(SUDEM,Ranger_map)
Plantation_Mask<-raster("Products/Raster_Products/Plantation_Mask.tif")
Ranger_map<-mask(Ranger_map,Plantation_Mask,inverse=TRUE)
NDVI<-mask(NDVI,Plantation_Mask,inverse=TRUE)
Drain_Dist<-mask(Drain_Dist,Plantation_Mask,inverse=TRUE)
Aspect<-mask(Aspect,Plantation_Mask,inverse=TRUE)
Cor_Width<-mask(Cor_Width,Plantation_Mask,inverse=TRUE)
Max_NBR<-mask(Max_NBR,Plantation_Mask,inverse=TRUE)
SUDEM<-mask(SUDEM,Plantation_Mask,inverse=TRUE)
writeRaster(Ranger_map,filename="Products/Raster_Products/Map_ranger_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(NDVI,filename="Products/Raster_Products/Sent_NDVI_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Drain_Dist,filename="Products/Raster_Products/Drain_Dist_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Aspect,filename="Products/Raster_Products/SUDEM_Aspect_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Cor_Width,filename="Products/Raster_Products/Cor_Width_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Max_NBR,filename="Products/Raster_Products/Max_NBR_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(SUDEM,filename="Products/Raster_Products/SUDEM_Grass.tif",format="GTiff",overwrite=TRUE)
rm(Ranger_map,NDVI,Drain_Dist,Aspect,Cor_Width,Max_NBR,SUDEM,Plantation_Mask)

# Load environmental raster data
Env<-load_var(path="~/Products/Raster_Products/",
              files=c("Sent_NDVI_Grass.tif","Drain_Dist_Grass.tif","SUDEM_Aspect_Grass.tif","Cor_Width_Grass.tif","Max_NBR_Grass.tif","SUDEM_Grass.tif"),
              format=".tif",Norm=TRUE)

# Load occurrence data
Occ<-load_occ(path="~/Products/Excel_Sheets_Products",
              Env=Env,file="Occurrence_Data.csv",Xcol="Long_X",Ycol="Lat_Y",Spcol="Species",sep=",")

# Perform stacked species distribution models
SSDM<-stack_modelling(algorithms=c("RF","SVM"),Occurrences=Occ,Env=Env,cores=4,method="pSSDM",
                      Xcol="Long_X",Ycol="Lat_Y",Pcol="Presence",Spcol="Species",tmp=TRUE,rep=10,
                      save=TRUE,name="SSDM",path="~/Products/")

# Variable importance 
variable.importance<-as.data.frame(t(read.csv("Products/SSDM/Stack/Tables/VarImp.csv",row.names=1)))
row.names(variable.importance)<-c("NDVI","Drainage","Aspect","Dist Edge","NBR","Elv")
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
rm(correlation,Occ,Env,SSDM)

############### Prepare data for Circuitscape ###############
# Load diversity map and select highest 95% of suitable pixel values
Diversity_map<-raster("Products/SSDM/Stack/Rasters/Diversity.tif")
Focal_Areas<-Diversity_map[[1]]>=quantile(Diversity_map,.95)
Focal_Areas[Focal_Areas==0]<-NA

# Convert to vector
Vector_Focal_Areas<-rasterToPolygons(Focal_Areas)
Vector_Focal_Areas$Obj_ID<-1:nrow(Vector_Focal_Areas)

# Dissolve vector
Vector_Focal_Areas_agg<-aggregate(Vector_Focal_Areas,dissolve=TRUE)
Vector_Focal_Areas_dis<-disaggregate(Vector_Focal_Areas_agg)
Vector_Focal_Areas_buf<-buffer(Vector_Focal_Areas_dis,width=0.001,dissolve=TRUE)
Vector_Focal_Areas<-disaggregate(Vector_Focal_Areas_buf)
rm(Vector_Focal_Areas_agg,Vector_Focal_Areas_dis,Vector_Focal_Areas_buf)

# Select only large focal areas to sustain grasshopper populations
# https://doi.org/10.1007/s10531-005-2356-1
Vector_Focal_Areas_features<-as(Vector_Focal_Areas,"SpatialPolygonsDataFrame")
Vector_Focal_Areas_features$Obj_ID<-1:nrow(Vector_Focal_Areas_features)
Vector_Focal_Areas_features$Area_sqm<-raster::area(Vector_Focal_Areas_features)
Vector_Focal_Areas_features<-st_as_sf(Vector_Focal_Areas_features)
Vector_Focal_Areas_features_100<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=10000) #_100_pixel region
Vector_Focal_Areas_features_50<-Vector_Focal_Areas_features %>%
  filter(Area_sqm>=5000) #_50_pixel region
Vector_Focal_Areas_features_100<-as(Vector_Focal_Areas_features_100,"Spatial")
Vector_Focal_Areas_features_50<-as(Vector_Focal_Areas_features_50,"Spatial")
rm(Vector_Focal_Areas_features,Focal_Areas,Vector_Focal_Areas)

# Save to aci file for Circuitscape
write.dir<-"~/Circuitscape/"
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
rm(Focal_Areas_Estate_1,Focal_Areas_Estate_2,Focal_Areas_Estate_3,Focal_Areas_Estate_4,
   Focal_Areas_100,Focal_Areas_50,Vector_Focal_Areas_features_100,Vector_Focal_Areas_features_50)

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
rm(Diversity_map,Diversity_map_Estate_1,Diversity_map_Estate_2,Diversity_map_Estate_3,Diversity_map_Estate_4)

# Load and transform NDVI
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI_Grass.tif")
Sent_NDVI<-resample(Sent_NDVI,rst_template)
#Plot.trans(PARM =c(3,100),Resistance=Sent_NDVI,transformation="Reverse Ricker")
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
rm(Sent_NDVI,NDVI_trans,NDVI_trans_Estate_1,NDVI_trans_Estate_2,NDVI_trans_Estate_3,NDVI_trans_Estate_4)

# Load and reclassify classification
Map_ranger<-raster("Products/Raster_Products/Map_ranger_Grass.tif")
Map_ranger<-resample(Map_ranger,rst_template,method="ngb")
#ID  Land use      Conductivity value (0-100)
#2 = Grassland   = 100
#3 = Ground      = 70
#5 = Shrubland   = 30
Map_ranger[Map_ranger==2]<-100
Map_ranger[Map_ranger==3]<-70
Map_ranger[Map_ranger==5]<-30
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
rm(Estate_1,Estate_2,Estate_3,Estate_4,Map_ranger,Map_ranger_Estate_1,Map_ranger_Estate_2,
   Map_ranger_Estate_3,Map_ranger_Estate_4,rst_template,write.dir)

############### Average current maps ###############
# Load land use data
Map_ranger<-raster("Products/Raster_Products/Map_ranger_Grass.tif")

# Load and merge connectivity maps
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

# Average current maps
Cur_Map<-stack(Cur_Map_SSDM,Cur_Map_NDVI,Cur_Map_Ranger)
Cur_Map_Mean<-calc(Cur_Map,fun=mean,na.rm=TRUE)
Cur_Map_Mean[Cur_Map_Mean<0]<-NA
Cur_Map_Mean[Cur_Map_Mean==0]<-NA
writeRaster(Cur_Map_Mean,filename="Products/Raster_Products/Cur_Map_Grass.tif",format="GTiff",overwrite=TRUE)
rm(Cur_Map,Cur_Map_SSDM,Cur_Map_NDVI,Cur_Map_Ranger)

# Load and merge focal areas
Focal_Areas_Estate_1<-raster("Circuitscape/Focal_Areas_Estate_1.asc")
Focal_Areas_Estate_2<-raster("Circuitscape/Focal_Areas_Estate_2.asc")
Focal_Areas_Estate_3<-raster("Circuitscape/Focal_Areas_Estate_3.asc")
Focal_Areas_Estate_4<-raster("Circuitscape/Focal_Areas_Estate_4.asc")
Focal_Areas<-merge(Focal_Areas_Estate_1,Focal_Areas_Estate_2,Focal_Areas_Estate_3,Focal_Areas_Estate_4)
Focal_Areas<-resample(Focal_Areas,Map_ranger)
rm(Focal_Areas_Estate_1,Focal_Areas_Estate_2,Focal_Areas_Estate_3,Focal_Areas_Estate_4,Map_ranger)

# Mask out focal areas from maps
Cur_Map_Mean<-crop(Cur_Map_Mean,Focal_Areas)
Cur_Map_Mean<-mask(Cur_Map_Mean,Focal_Areas,inverse=TRUE)

# Save focal areas as shapefile
Focal_Areas<-rasterToPolygons(Focal_Areas,dissolve=TRUE) 
writeOGR(Focal_Areas,"Products/Shapefiles/","Focal_Areas_Grass",driver="ESRI Shapefile")
rm(Focal_Areas)

# Clip per plantation for visualization
Estate_1<-readOGR("Shapefiles/Study_Estate1_Extent.shp")
Cur_Map_Mean_Estate_1<-crop(Cur_Map_Mean,Estate_1)
Cur_Map_Mean_Estate_1<-mask(Cur_Map_Mean_Estate_1,Estate_1)
Estate_2<-readOGR("Shapefiles/Study_Estate2_Extent.shp")
Cur_Map_Mean_Estate_2<-crop(Cur_Map_Mean,Estate_2)
Cur_Map_Mean_Estate_2<-mask(Cur_Map_Mean_Estate_2,Estate_2)
Estate_3<-readOGR("Shapefiles/Study_Estate3_Extent.shp")
Cur_Map_Mean_Estate_3<-crop(Cur_Map_Mean,Estate_3)
Cur_Map_Mean_Estate_3<-mask(Cur_Map_Mean_Estate_3,Estate_3)
Estate_4<-readOGR("Shapefiles/Study_Estate4_Extent.shp")
Cur_Map_Mean_Estate_4<-crop(Cur_Map_Mean,Estate_4)
Cur_Map_Mean_Estate_4<-mask(Cur_Map_Mean_Estate_4,Estate_4)
writeRaster(Cur_Map_Mean_Estate_1,filename="Products/Raster_Products/Cur_Map_Mean_Estate_1_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Cur_Map_Mean_Estate_2,filename="Products/Raster_Products/Cur_Map_Mean_Estate_2_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Cur_Map_Mean_Estate_3,filename="Products/Raster_Products/Cur_Map_Mean_Estate_3_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Cur_Map_Mean_Estate_4,filename="Products/Raster_Products/Cur_Map_Mean_Estate_4_Grass.tif",format="GTiff",overwrite=TRUE)
rm(Cur_Map_Mean,Estate_1,Estate_2,Estate_3,Estate_4,Cur_Map_Mean_Estate_1,Cur_Map_Mean_Estate_2,Cur_Map_Mean_Estate_3,Cur_Map_Mean_Estate_4)

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
rm(ras_temp_Estate1,ras_temp_Estate2,ras_temp_Estate3,ras_temp_Estate4)

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
writeRaster(Overall_Estate_1_Density,filename="Products/Raster_Products/Overall_Estate_1_Density_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_2_Density,filename="Products/Raster_Products/Overall_Estate_2_Density_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_3_Density,filename="Products/Raster_Products/Overall_Estate_3_Density_Grass.tif",format="GTiff",overwrite=TRUE)
writeRaster(Overall_Estate_4_Density,filename="Products/Raster_Products/Overall_Estate_4_Density_Grass.tif",format="GTiff",overwrite=TRUE)
rm(Diversity_map_Estate_1_cs,Diversity_map_Estate_2_cs,Diversity_map_Estate_3_cs,Diversity_map_Estate_4_cs,
   Focal_Areas_Estate_1,Focal_Areas_Estate_2,Focal_Areas_Estate_3,Focal_Areas_Estate_4,Overall_Estate_1,
   Overall_Estate_2,Overall_Estate_3,Overall_Estate_4,Diversity_map_Estate_1,Diversity_map_Estate_2,
   Diversity_map_Estate_3,Diversity_map_Estate_4)

# Prune networks
Overall_Estate_1_Density<-raster("Products/Raster_Products/Overall_Estate_1_Density_Grass.tif")
Overall_Estate_1_Density[Overall_Estate_1_Density<2]<-NA
Overall_Estate_2_Density<-raster("Products/Raster_Products/Overall_Estate_2_Density_Grass.tif")
Overall_Estate_2_Density[Overall_Estate_2_Density<2]<-NA
Overall_Estate_3_Density<-raster("Products/Raster_Products/Overall_Estate_3_Density_Grass.tif")
Overall_Estate_3_Density[Overall_Estate_3_Density<2]<-NA
Overall_Estate_4_Density<-raster("Products/Raster_Products/Overall_Estate_4_Density_Grass.tif")
Overall_Estate_4_Density[Overall_Estate_4_Density<2]<-NA

# Merge
Overall_Density<-raster::merge(Overall_Estate_1_Density,Overall_Estate_2_Density,Overall_Estate_3_Density,Overall_Estate_4_Density)
writeRaster(Overall_Density,filename="Products/Raster_Products/Overall_Density_Grass.tif",format="GTiff",overwrite=TRUE)
rm(Overall_Density,Overall_Estate_1_Density,Overall_Estate_2_Density,Overall_Estate_3_Density,Overall_Estate_4_Density)

############### Identify restoration areas ###############
# Load least cost networks
Overall_Density<-raster("Products/Raster_Products/Overall_Density_Grass.tif")

# Load current maps and clip to estates
Extent<-shapefile("Shapefiles/Plantations_Clipped.shp",stringsAsFactors=FALSE)
Cur_Map<-raster("Products/Raster_Products/Cur_Map_Grass.tif")
Cur_Map<-crop(Cur_Map,Extent)
Cur_Map<-mask(Cur_Map,Extent)
Estate_1<-readOGR("Shapefiles/Study_Estate1_Extent.shp")
Estate_2<-readOGR("Shapefiles/Study_Estate2_Extent.shp")
Estate_3<-readOGR("Shapefiles/Study_Estate3_Extent.shp")
Estate_4<-readOGR("Shapefiles/Study_Estate4_Extent.shp")
Estate_1_Cur_Map<-crop(Cur_Map,Estate_1)
Estate_1_Cur_Map<-mask(Estate_1_Cur_Map,Estate_1)
Estate_2_Cur_Map<-crop(Cur_Map,Estate_2)
Estate_2_Cur_Map<-mask(Estate_2_Cur_Map,Estate_2)
Estate_3_Cur_Map<-crop(Cur_Map,Estate_3)
Estate_3_Cur_Map<-mask(Estate_3_Cur_Map,Estate_3)
Estate_4_Cur_Map<-crop(Cur_Map,Estate_4)
Estate_4_Cur_Map<-mask(Estate_4_Cur_Map,Estate_4)

# Select lowest 85%
Estate_1_Cur_Map<-Estate_1_Cur_Map[[1]]<=quantile(Estate_1_Cur_Map,.85)
Estate_1_Cur_Map[Estate_1_Cur_Map==0]<-NA
Estate_2_Cur_Map<-Estate_2_Cur_Map[[1]]<=quantile(Estate_2_Cur_Map,.85)
Estate_2_Cur_Map[Estate_2_Cur_Map==0]<-NA
Estate_3_Cur_Map<-Estate_3_Cur_Map[[1]]<=quantile(Estate_3_Cur_Map,.85)
Estate_3_Cur_Map[Estate_3_Cur_Map==0]<-NA
Estate_4_Cur_Map<-Estate_4_Cur_Map[[1]]<=quantile(Estate_4_Cur_Map,.85)
Estate_4_Cur_Map[Estate_4_Cur_Map==0]<-NA
rm(Estate_1,Estate_2,Estate_3,Estate_4,Cur_Map,Extent)

# Select the network over low current areas
Estate_1_Overall_Density<-crop(Overall_Density,Estate_1_Cur_Map)
Estate_1_Overall_Density<-mask(Estate_1_Overall_Density,Estate_1_Cur_Map)
Estate_2_Overall_Density<-crop(Overall_Density,Estate_2_Cur_Map)
Estate_2_Overall_Density<-mask(Estate_2_Overall_Density,Estate_2_Cur_Map)
Estate_3_Overall_Density<-crop(Overall_Density,Estate_3_Cur_Map)
Estate_3_Overall_Density<-mask(Estate_3_Overall_Density,Estate_3_Cur_Map)
Estate_4_Overall_Density<-crop(Overall_Density,Estate_4_Cur_Map)
Estate_4_Overall_Density<-mask(Estate_4_Overall_Density,Estate_4_Cur_Map)
rm(Overall_Density,Estate_1_Cur_Map,Estate_2_Cur_Map,Estate_3_Cur_Map,Estate_4_Cur_Map)

# Subset to include highly used corridors only
Estate_1_Overall_Density<-Estate_1_Overall_Density[[1]]>=quantile(Estate_1_Overall_Density,.95)
Estate_1_Overall_Density[Estate_1_Overall_Density==0]<-NA
Estate_2_Overall_Density<-Estate_2_Overall_Density[[1]]>=quantile(Estate_2_Overall_Density,.95)
Estate_2_Overall_Density[Estate_2_Overall_Density==0]<-NA
Estate_3_Overall_Density<-Estate_3_Overall_Density[[1]]>=quantile(Estate_3_Overall_Density,.95)
Estate_3_Overall_Density[Estate_3_Overall_Density==0]<-NA
Estate_4_Overall_Density<-Estate_4_Overall_Density[[1]]>=quantile(Estate_4_Overall_Density,.95)
Estate_4_Overall_Density[Estate_4_Overall_Density==0]<-NA

# Save
Overall_Density<-merge(Estate_1_Overall_Density,Estate_2_Overall_Density,Estate_3_Overall_Density,Estate_4_Overall_Density)
writeRaster(Overall_Density,filename="Products/Raster_Products/Overall_Density_Restore_Grass.tif",format="GTiff",overwrite=TRUE)
rm(Estate_1_Overall_Density,Estate_2_Overall_Density,Estate_3_Overall_Density,Estate_4_Overall_Density)

# Convert to points for visualization
Overall_Density_Point<-rasterToPoints(Overall_Density,spatial=TRUE)
writeOGR(Overall_Density_Point,"Products/Shapefiles/","Overall_Density_Point_Restore_Grass",driver="ESRI Shapefile")
Overall_Density<-raster("Products/Raster_Products/Overall_Density_Grass.tif")
Overall_Density_Point<-rasterToPoints(Overall_Density,spatial=TRUE)
writeOGR(Overall_Density_Point,"Products/Shapefiles/","Overall_Density_Point_Grass",driver="ESRI Shapefile")
rm(Overall_Density_Point,Overall_Density)

############### Drivers of turnover gdm ###############
# Load data environmental data
Cur_Map<-raster("Products/Raster_Products/Cur_Map_Grass.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width_Grass.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist_Grass.tif")
Max_NBR<-stack("Products/Raster_Products/Max_NBR_Grass.tif")
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI_Grass.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM_Grass.tif")
Aspect<-raster("Products/Raster_Products/SUDEM_Aspect_Grass.tif")
Env<-stack(Sent_NDVI,Aspect,Max_NBR,Cur_Map,Drain_Dist,Cor_Width,SUDEM)
rm(Sent_NDVI,Aspect,Max_NBR,Cur_Map,Drain_Dist,Cor_Width,SUDEM)
names(Env)
names(Env)<-c("NDVI","Aspect","Fire_Hist","Func_Conn","Drain_Dist","Cor_Width","Elv")

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
Cor_Width<-raster("Products/Raster_Products/Cor_Width_Grass.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM_Grass.tif")
Env<-stack(Cor_Width,SUDEM)
rm(Cor_Width,SUDEM)
names(Env)
names(Env)<-c("Cor_Width","Elv")
Data<-formatsitepair(bioData=Bio,bioFormat=1,predData=Env,sppColumn="Species",dist="jaccard",
                     siteColumn="Site_ID",XColumn="Long_X",YColumn="Lat_Y")
Data<-na.omit(Data)

# GDM
model_1<-gdm(data=Data,geo=TRUE)
summary(model_1)
plot(model_1,plot.layout=c(3,2))
model_1.splineDat<-isplineExtract(model_1)
par(mfrow=c(4,1))
par(mar=c(4,4,2,2))
par(oma=c(1,1,1,1))
plot(model_1.splineDat$x[,"Cor_Width"],model_1.splineDat$y[,"Cor_Width"],ylim=c(0,0.7),
     lwd=2,type="l",xlab="Edge distance",ylab="Dissimilarity")
plot(model_1.splineDat$x[,"Elv"],model_1.splineDat$y[,"Elv"],ylim=c(0,0.7),
     lwd=2,type="l",xlab="Elevation",ylab="Dissimilarity")
plot(model_1.splineDat$x[,"Geographic"],model_1.splineDat$y[,"Geographic"],ylim=c(0,0.7),
     lwd=2,type="l",xlab="Geographic distance",ylab="Dissimilarity")
plot(x=model_1[["ecological"]],y=model_1[["observed"]],ylim=c(0,1),
     type="n",xlab="Pred. dissimilarity",ylab="Obs. dissimilarity")
points(x=model_1[["ecological"]],y=model_1[["observed"]],cex=0.5,lwd=0.5,col="black")
overlayX<-seq(from=min(model_1[["ecological"]]),to=max(model_1[["ecological"]]),length=200)
overlayY<-1-exp(-overlayX)
lines(overlayX,overlayY,lwd=2)
par(mar=c(5.1,4.1,4.1,2.1))
rm(model_1_Imp,overlayX,overlayY,model_1,model_1.splineDat,Env,Data,Bio)

############### Community analysis manyglm ###############
#Load raster data
Max_NBR<-raster("Products/Raster_Products/Max_NBR_Grass.tif")
Sent_NDVI<-raster("Products/Raster_Products/Sent_NDVI_Grass.tif")
Drain_Dist<-raster("Products/Raster_Products/Drain_Dist_Grass.tif")
SUDEM_Aspect<-raster("Products/Raster_Products/SUDEM_Aspect_Grass.tif")
Cor_Width<-raster("Products/Raster_Products/Cor_Width_Grass.tif")
Func_Con<-raster("Products/Raster_Products/Cur_Map_Grass.tif")
SUDEM<-raster("Products/Raster_Products/SUDEM_Grass.tif")

#Extract values over GPS
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
GPS<-SpatialPointsDataFrame(GPS[,2:3],GPS,proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
GPS<-spTransform(GPS,CRS="+proj=utm +zone=36 +south +datum=WGS84 +units=m +no_defs")
Max_NBR<-raster::extract(Max_NBR,GPS,fun=mean,buffer=50)
Sent_NDVI<-raster::extract(Sent_NDVI,GPS,fun=mean,buffer=50)
Drain_Dist<-raster::extract(Drain_Dist,GPS,fun=mean,buffer=50)
SUDEM_Aspect<-raster::extract(SUDEM_Aspect,GPS,fun=mean,buffer=50)
Cor_Width<-raster::extract(Cor_Width,GPS,fun=mean,buffer=50)
Func_Con<-raster::extract(Func_Con,GPS,fun=mean,buffer=50)
SUDEM<-raster::extract(SUDEM,GPS,fun=mean,buffer=50)
Env<-as.data.frame(cbind("Fire"=Max_NBR,"NDVI"=Sent_NDVI,"Drain"=Drain_Dist,"Aspect"=SUDEM_Aspect,
                        "Edge_Dist"=Cor_Width,"Func_Con"=Func_Con,"Elv"=SUDEM))
Env[is.na(Env)]<-0
rm(Max_NBR,Sent_NDVI,Drain_Dist,SUDEM_Aspect,Cor_Width,Func_Con,SUDEM,GPS)

#Check for correlation between features
ggcorrplot(cor(Env,method="spearman"),
           type="lower",lab=TRUE,lab_size=3,tl.cex= 10)
Env$Fire=NULL

# Standerdise variables and add random effect
Env<-standardize(Env)
GPS<-read.csv("Excel_Sheets/Site_GPS_coordinates.csv")
Env<-as.data.frame(cbind(Env,"Plantation"=GPS$Plantation))
rm(GPS)

# Load biological data
Bio<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
Bio<-Bio[,colSums(Bio!=0)>3]

# Load trait data
Traits<-read.csv("Excel_Sheets/Caelifera_assemblage_data.csv")
Traits<-Traits[,colSums(Traits!=0)>3]
groups<-c("Low","Intermediate","Low","Low","High","Low","Low","Intermediate","Intermediate","Low","Low",
          "Intermediate","Low","Low","Low","Low","Intermediate","Low","Low","Intermediate","Low","Intermediate",
          "Intermediate","Intermediate","Intermediate","Low","High","Intermediate","Low","High","High","Low",
          "Intermediate","Low","Intermediate","High","Intermediate","Intermediate","Intermediate","Intermediate",
          "Intermediate","Intermediate","High","High","Intermediate","High","Low","Intermediate","Low","Intermediate","Low")
Traits<-t(Traits)
Traits<-as.data.frame(cbind(Traits,"Groups"=groups))
Traits<-as.data.frame(Traits$Groups)
colnames(Traits)<-"Groups"
rm(groups)

# Create dummy variables
Traits<-data.frame(predict(dummyVars("~Groups",data=Traits),newdata=Traits))
colnames(Traits)<-c("High","Intermediate","Low")

# Define model
mvSpecies<-mvabund(Bio)
model<-manyglm(mvSpecies~.,family="negative binomial",data=Env)

# Inspect model
plot(model)
resP<-residuals(model)
qqnorm(resP)
qqline(resP,col="red")
rm(resP,mvSpecies)

# Model statistics
anova(model,resamp="monte-carlo",cor.type="shrink",test="score")
'''
Time elapsed: 0 hr 3 min 11 sec
Analysis of Variance Table

Model: mvSpecies ~ .

Multivariate test:
            Res.Df Df.diff  score Pr(>score)    
(Intercept)     50                              
NDVI            49       1  31.52      0.230    
Drain           48       1  30.93      0.303    
Aspect          47       1  23.92      0.817    
Edge_Dist       46       1  41.26      0.042 *  
Func_Con        45       1  58.81      0.003 ** 
Elv             44       1  65.38      0.001 ***
Plantation      41       3 127.13      0.009 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Arguments:
 Test statistics calculated assuming correlated response via ridge regularization 
 P-value calculated using 999 iterations via monte-carlo resampling.
'''
rm(model)

# 4th corner
corner<-traitglm(Bio,Env,Traits,family="negative.binomial",method="manyglm")
plot(corner)
qqnorm(residuals(corner)); abline(c(0,1),col="red")
anova(corner,resamp="monte-carlo",test="score")
'''
Multivariate test:
                          Res.Df Df.diff score Pr(>score)  
Main effects only           1735                           
env:trait (fourth corner)   1717      18 33.77       0.04 *
'''

# Plot 4th corner
fourth.corner<-as.data.frame(corner$fourth.corner)
fourth.corner<-fourth.corner %>%
  dplyr::select(c("Edge_Dist","Func_Con","Elv"))
rownames(fourth.corner)<-c("High","Intermediate","Low")
a = max(abs(fourth.corner))
colort = colorRampPalette(c("blue","white","red"))
plot.4th = levelplot(t(as.matrix(fourth.corner)),xlab="Environmental Variables",
ylab ="Conservation groupings",col.regions=colort(100),at=seq(-a,a,length=100),
scales = list(x=list(rot=45)))
print(plot.4th)
rm(a,colort,plot.4th,corner,fourth.corner,Bio,Env,Traits)
