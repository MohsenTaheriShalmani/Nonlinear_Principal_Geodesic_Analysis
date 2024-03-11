library(shapes)
library(rgl)
library(ks)
library(matrixcalc) 
library(mvShapiroTest)

#####################################################################################################
#####################################################################################################
#set working directory to file location

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#####################################################################################################
#####################################################################################################
# functions 

# function to calculate Euclidean distance
euclideanDistance<-function(x,y){
  return(sqrt((x[1]-y[1])^2+(x[2]-y[2])^2+(x[3]-y[3])^2))
}

# functon to project a point on a line (to improve up and down spokes' lengths)
projectPointOnLine<-function(a,b,p){
  ap = p-a
  ab = b-a
  projection = a + sum(ap*ab)/sum(ab*ab) * ab
  return(projection)
}

# function to calculate Mahalanobis distance (Hotelling Metric)
MahalanobisDistance<-function(k,nx,ny,X,Y){
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  return(T2)
}

# Hotelling T2 test
HotellingT2<-function(X,Y){
  if(dim(X)[2]!=dim(Y)[2]){
    cat("Dimention Error!\n")
    break
  }
  k<-dim(X)[2]
  nx<-dim(X)[1]
  ny<-dim(Y)[1]
  Sx<-cov(X)
  Sy<-cov(Y)
  meanX<-colMeans(X)
  meanY<-colMeans(Y)
  n<-nx+ny-1
  S<-((nx-1)*Sx+(ny-1)*Sy)/(nx+ny-2) #S=pooled covariance matrix
  T2<-t(meanX-meanY)%*%solve(S*(1/nx+1/ny))%*%(meanX-meanY)
  F_value<-((n-k)/(k*(n-1)))*T2
  df1<-k
  df2<-n-k
  p_value<-1-pf(F_value,df1,df2)
  return(p_value)
}

# library(BisRNA)
# The fisher.method function is also available in BisRNA library
fisher.method<-function (pvalues)
{
  df <- 2 * length(pvalues)
  global_pValue<-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = FALSE)
  # global_pValue<-1-pchisq(-2 * sum(log(pvalues), na.rm = TRUE), df, lower.tail = TRUE)
  return(global_pValue)
}
#####################################################################################################
#####################################################################################################
# read SPHARM-PDM data 

#choose two groups to compare from {"CG", "PD", "PDD1", "PDD2", "PDD3", "PDD4", "PDD_1_2_3_4"} e.g. "CG" and "PD"

G1<-"CG"
G2<-"PD"

#read meshes
meshPoints_G1 <- read.csv(file = paste("../files/meshPoints_",G1,".csv",sep = ""),check.names = FALSE, header=TRUE, sep=",")
meshPoints_G2 <- read.csv(file = paste("../files/meshPoints_",G2,".csv",sep = ""),check.names = FALSE, header=TRUE, sep=",")
#View(meshPoints_G1)
#View(meshPoints_G2)

# connections of triangular mesh
PolygonsCsv <- read.csv("../files/PolygonsCsv.csv")
polyMatrix<-cbind((PolygonsCsv$point1+1),(PolygonsCsv$point2+1),(PolygonsCsv$point3+1))


#number of samples
nSamplesG1<-dim(meshPoints_G1)[2]
nSamplesG1
nSamplesG2<-dim(meshPoints_G2)[2]
nSamplesG2

#####################################################################################################
#####################################################################################################
# extract all SPHARM-PDM points

numberOfPoints<-length(meshPoints_G1[,1])/3

spharmPDM_G1<-array(NA, dim = c(numberOfPoints,3,nSamplesG1))
k<-1
for (i in names(meshPoints_G1)){
  spharmPDM_G1[,,k]<-matrix(meshPoints_G1[[i]], ncol = 3, byrow = TRUE)
  k<-k+1
}

spharmPDM_G2<-array(NA, dim = c(numberOfPoints,3,nSamplesG2))
k<-1
for (i in names(meshPoints_G2)){
  spharmPDM_G2[,,k]<-matrix(meshPoints_G2[[i]], ncol = 3, byrow = TRUE)
  k<-k+1
}

#####################################################################################################
#####################################################################################################
# t-test for SPHARM-PDM centroid size

allSpharm<-list()
allSpharm$group<-factor(c(rep(G1,nSamplesG1),rep(G2,nSamplesG2)))
allSpharm$x<-abind(spharmPDM_G1,spharmPDM_G2)
plot(allSpharm$group,centroid.size(allSpharm$x), main="Centroid size of SPHARM-PDM", col = c("blue","red"),xlab="",ylab="")
testResult<-t.test(centroid.size(spharmPDM_G1),centroid.size(spharmPDM_G2)) 
testResult
mean(centroid.size(spharmPDM_G1))
mean(centroid.size(spharmPDM_G2))
sd(centroid.size(spharmPDM_G1))
sd(centroid.size(spharmPDM_G2))

#####################################################################################################
#####################################################################################################
# compare SPHARM-PDM means visually shape analysis

# find and compare mean shapes
# align SPHARM-PDMs it takes a minute!
procG1<-procGPA(spharmPDM_G1,scale = T) #scale=T for shape analysis
procG2<-procGPA(spharmPDM_G2,scale = T)
procMeans<-procOPA(procG1$mshape,procG2$mshape,scale = F)

#plot
vertsG1 <- rbind(t(as.matrix(procMeans$Ahat)),1)
trglsG1 <- as.matrix(t(polyMatrix))
tmeshG1 <- tmesh3d(vertsG1, trglsG1)
vertsG2 <- rbind(t(as.matrix(procMeans$Bhat)),1)
trglsG2 <- as.matrix(t(polyMatrix))
tmeshG2 <- tmesh3d(vertsG2, trglsG2)
wire3d(tmeshG1, col="blue")         # G1 mean as a wire mesh
# shade3d(tmeshG1, col="blue")
# wire3d(tmeshG2, col="red")        # G2 mean as a mesh
shade3d(tmeshG2, col="red",alpha=0.8)

#####################################################################################################
#####################################################################################################
# compare SPHARM-PDM means visually size-and-shape analysis

# find and compare mean shapes
# align SPHARM-PDMs it takes a minute!
procG1<-procGPA(spharmPDM_G1,scale = F) #scale=F for size-and-shape analysis
procG2<-procGPA(spharmPDM_G2,scale = F)
procMeans<-procOPA(procG1$mshape,procG2$mshape,scale = F)

#plot
vertsG1 <- rbind(t(as.matrix(procMeans$Ahat)),1)
trglsG1 <- as.matrix(t(polyMatrix))
tmeshG1 <- tmesh3d(vertsG1, trglsG1)
vertsG2 <- rbind(t(as.matrix(procMeans$Bhat)),1)
trglsG2 <- as.matrix(t(polyMatrix))
tmeshG2 <- tmesh3d(vertsG2, trglsG2)
wire3d(tmeshG1, col="blue")         # G1 mean as a wire mesh
# shade3d(tmeshG1, col="blue")
# wire3d(tmeshG2, col="red")        # G2 mean as a mesh
shade3d(tmeshG2, col="red",alpha=0.8)

#####################################################################################################
#####################################################################################################
# Hotelling T2 Test for SPHARM-PDM shape analysis

# align SPHARM-PDMs it takes a minute!
all_spharmPDM<-abind(spharmPDM_G1,spharmPDM_G2)
procSPHARM_PDM<-procGPA(all_spharmPDM, scale = T) #scale=T for shape analysis
alignedSpharmPDM_G1<-procSPHARM_PDM$rotated[,,1:nSamplesG1]
alignedSpharmPDM_G2<-procSPHARM_PDM$rotated[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]

# Hotelling test
pValues_Hotelling_SPHARM_PDM<-c()
for (i in 1:numberOfPoints) {
  pValues_Hotelling_SPHARM_PDM<-c(pValues_Hotelling_SPHARM_PDM,
                                  HotellingT2(t(alignedSpharmPDM_G1[i,,]),t(alignedSpharmPDM_G2[i,,])))
}

alpha<-0.05
significantpValues_HotellingSpharmPDM<-which(pValues_Hotelling_SPHARM_PDM<=alpha)
significantpValues_HotellingSpharmPDM

FDR<-0.05
pValues_Hotelling_SPHARM_PDM_BH<-p.adjust(pValues_Hotelling_SPHARM_PDM,method = "BH")
significantpValues_HotellingSpharmPDM_BH<-which(pValues_Hotelling_SPHARM_PDM_BH<=FDR)
significantpValues_HotellingSpharmPDM_BH 

#plot significant without correction
open3d()
plot3d(procSPHARM_PDM$mshape,type="s", radius = 0.1,col = "blue",expand = 10, box=FALSE,add = TRUE)
plot3d(procSPHARM_PDM$mshape[significantpValues_HotellingSpharmPDM,],type="s",
       radius = 0.3 ,col = "red",expand = 10,box=FALSE,add = TRUE)
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="white",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#plot significant by FDR correction
open3d()
plot3d(procSPHARM_PDM$mshape,type="s", radius = 0.1,col = "blue",expand = 10, box=FALSE,add = TRUE)
plot3d(procSPHARM_PDM$mshape[significantpValues_HotellingSpharmPDM_BH,],type="s",
       radius = 0.3 ,col = "red",expand = 10,box=FALSE,add = TRUE)
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="white",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#####################################################################################################
#####################################################################################################
# Global hypothesis test for parametric approach for shape 

# Bonferroni test for spharm
alpha<-0.05
alpha/length(pValues_Hotelling_SPHARM_PDM)
min(pValues_Hotelling_SPHARM_PDM)

# Fisher test for SPHARM-PDM 
fisher.method(pValues_Hotelling_SPHARM_PDM)
fisher.method(pValues_Hotelling_SPHARM_PDM_BH)

#plot SPHARM-PDM p-values and adjusted p-values
plot(1:length(pValues_Hotelling_SPHARM_PDM),sort(p.adjust(pValues_Hotelling_SPHARM_PDM,method = "bonferroni")),col="orange", ylim = c(0,1),xlab = "",ylab = "" )
par(new=TRUE)
plot(1:length(pValues_Hotelling_SPHARM_PDM),sort(pValues_Hotelling_SPHARM_PDM),col="red", ylim = c(0,1),xlab = "",ylab = "" )
par(new=TRUE)
plot(1:length(pValues_Hotelling_SPHARM_PDM_BH),sort(pValues_Hotelling_SPHARM_PDM_BH),col="blue", ylim = c(0,1),xlab = "",ylab = "")
title(main = "SPHARM-PDM p-values and adjusted p-values", sub = "",
      xlab = "", ylab = "",
      cex.main = 1,   font.main= 4, col.main= "black",
      cex.sub = 0.75, font.sub = 3, col.sub = "green",
      col.lab ="darkblue"
)
legend("bottomright",c("p-values","Benjamini-Hochberg","Bonferroni"),fill=c("red","blue","orange"),cex = 0.6,text.width = 150)



#####################################################################################################
#####################################################################################################
# Hotelling T2 Test for SPHARM-PDM shape size-and-shape analysis

# align SPHARM-PDMs it takes a minute!
all_spharmPDM<-abind(spharmPDM_G1,spharmPDM_G2)
procSPHARM_PDM<-procGPA(all_spharmPDM, scale = F) #scale=F for shape analysis
alignedSpharmPDM_G1<-procSPHARM_PDM$rotated[,,1:nSamplesG1]
alignedSpharmPDM_G2<-procSPHARM_PDM$rotated[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]

# Hotelling test
pValues_Hotelling_SPHARM_PDM<-c()
for (i in 1:numberOfPoints) {
  pValues_Hotelling_SPHARM_PDM<-c(pValues_Hotelling_SPHARM_PDM,
                                  HotellingT2(t(alignedSpharmPDM_G1[i,,]),t(alignedSpharmPDM_G2[i,,])))
}

alpha<-0.05
significantpValues_HotellingSpharmPDM<-which(pValues_Hotelling_SPHARM_PDM<=alpha)
significantpValues_HotellingSpharmPDM

FDR<-0.05
pValues_Hotelling_SPHARM_PDM_BH<-p.adjust(pValues_Hotelling_SPHARM_PDM,method = "BH")
significantpValues_HotellingSpharmPDM_BH<-which(pValues_Hotelling_SPHARM_PDM_BH<=FDR)
significantpValues_HotellingSpharmPDM_BH 

#plot significant without correction
open3d()
plot3d(procSPHARM_PDM$mshape,type="s", radius = 0.1,col = "blue",expand = 10, box=FALSE,add = TRUE)
plot3d(procSPHARM_PDM$mshape[significantpValues_HotellingSpharmPDM,],type="s",
       radius = 0.3 ,col = "red",expand = 10,box=FALSE,add = TRUE)
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="white",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#plot significant by FDR correction
open3d()
plot3d(procSPHARM_PDM$mshape,type="s", radius = 0.1,col = "blue",expand = 10, box=FALSE,add = TRUE)
plot3d(procSPHARM_PDM$mshape[significantpValues_HotellingSpharmPDM_BH,],type="s",
       radius = 0.3 ,col = "red",expand = 10,box=FALSE,add = TRUE)
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="white",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)


#####################################################################################################
#####################################################################################################
# Global hypothesis test for parametric approach size-and-shape

# Bonferroni test for spharm
alpha<-0.05
alpha/length(pValues_Hotelling_SPHARM_PDM_BH)
min(pValues_Hotelling_SPHARM_PDM_BH)

# Fisher test for SPHARM-PDM 
fisher.method(pValues_Hotelling_SPHARM_PDM_BH)

#plot SPHARM-PDM p-values and adjusted p-values
plot(1:length(pValues_Hotelling_SPHARM_PDM),sort(p.adjust(pValues_Hotelling_SPHARM_PDM,method = "bonferroni")),col="orange", ylim = c(0,1),xlab = "",ylab = "" )
par(new=TRUE)
plot(1:length(pValues_Hotelling_SPHARM_PDM),sort(pValues_Hotelling_SPHARM_PDM),col="red", ylim = c(0,1),xlab = "",ylab = "" )
par(new=TRUE)
plot(1:length(pValues_Hotelling_SPHARM_PDM_BH),sort(pValues_Hotelling_SPHARM_PDM_BH),col="blue", ylim = c(0,1),xlab = "",ylab = "")
title(main = "SPHARM-PDM p-values and adjusted p-values", sub = "",
      xlab = "", ylab = "",
      cex.main = 1,   font.main= 4, col.main= "black",
      cex.sub = 0.75, font.sub = 3, col.sub = "green",
      col.lab ="darkblue"
)
legend("bottomright",c("p-values","Benjamini-Hochberg","Bonferroni"),fill=c("red","blue","orange"),cex = 0.6,text.width = 150)


#####################################################################################################
#####################################################################################################
########### Hypothesis test independent from the alignment by distance matrix #######################
#####################################################################################################
#####################################################################################################
# functions

# function to adjust upper triangle portion of the test matrix
adjustByUpperTriangle <- function(pvalMatrix, method) {
  A<-pvalMatrix
  d<-p.adjust(A[upper.tri(A)],method = method)
  B<-array(0, dim = dim(A))
  k<-1
  k2<-1
  for (j in 2:dim(A)[1]) {
    for (i in 1:k2) {
      B[i,j]<-d[k]
      k<-k+1
    }
    k2<-k2+1
  }
  return(B+t(B)+diag(dim(B)[1]))
}

# function to calculate euclidean distance tensor
euclideanDistanceTensor <- function(landmarkGroup) {
  q<-dim(landmarkGroup)[1]
  m<-dim(landmarkGroup)[2]
  n<-dim(landmarkGroup)[3]
  tensor<-array(NA, dim=c(q,q,n))
  for (k in 1:n) {
    for (i in 1:q) {
      for (j in 1:q) {
        tensor[i,j,k]<-euclideanDistance(landmarkGroup[i,,k],landmarkGroup[j,,k])
      }
    }
    cat("Sample",k,"from",n,"samples is done.\n")
    
  }
  return(tensor)
}

# function to calculate geodesic distance tensor
geodesicDistanceTensor <- function(directionGroup) {
  q<-dim(directionGroup)[1]
  m<-dim(directionGroup)[2]
  n<-dim(directionGroup)[3]
  tensor<-array(NA, dim=c(q,q,n))
  for (k in 1:n) {
    for (i in 1:q) {
      for (j in 1:q) {   
        tensor[i,j,k]<-acos(pmin(pmax(
          as.numeric(directionGroup[i,,k]%*%directionGroup[j,,k]) ,-1.0), 1.0)) #pmin pmax is to make sure value is between [-1,1]
      }
    }
    cat("Sample",k,"from",n,"samples is done.\n")
    
  }
  return(tensor)
}


# function to calculate p-value matrix of two distance tensors based on T-test
calculatePvalMatrix <- function(tensor1, tensor2) {
  q<-dim(tensor1)[1]
  pvalMatrix<-array(NA, dim=c(q,q))
  for (i in 1:q) {
    for (j in 1:q) {
      if(i!=j){
        pvalMatrix[i,j]<-t.test(tensor1[i,j,],tensor2[i,j,])$p.value
      }else{
        pvalMatrix[i,j]<-1
      }
    }
    cat("Element",i,"from",q,"elements is done.\n")
  } 
  return(pvalMatrix)
}

#####################################################################################################
#####################################################################################################
# Hypothesis test independent from alignment for SPHARM-PDM by removing scale

# align SPHARM-PDMs it takes a minute!
all_spharmPDM<-abind(spharmPDM_G1,spharmPDM_G2)
procSPHARM_PDM<-procGPA(all_spharmPDM, scale = T) #scale=T for shape analysis
alignedSpharmPDM_G1<-procSPHARM_PDM$rotated[,,1:nSamplesG1]
alignedSpharmPDM_G2<-procSPHARM_PDM$rotated[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]


# comment out lines below to calculate SPHARM-PDM distance tensors and p-value matrix!!!

# # takes some time !!!
# spharmTensorScaleTG1<-euclideanDistanceTensor(alignedSpharmPDM_G1)
# spharmTensorScaleTG2<-euclideanDistanceTensor(alignedSpharmPDM_G2)
# pvalMatrix_spharmScaleT<-calculatePvalMatrix(spharmTensorScaleTG1, spharmTensorScaleTG2)
# save(pvalMatrix_spharmScaleT,file = "../files/pvalMatrix_spharmScaleT.Rdata")

load("../files/pvalMatrix_spharmScaleT.Rdata")

# pvalue matrix and adjustment
pvalMatrix_spharm<-pvalMatrix_spharmScaleT
pvalMatrix_spharm_BH<-adjustByUpperTriangle(pvalMatrix_spharmScaleT,method = "BH")


# calculate kde of mirrored p-values
h<-0.05
pvalMatrix<-pvalMatrix_spharm
# pvalMatrix<-pvalMatrix_spharm_BH  #comment out to use adjusted p-value matrix
a<-pvalMatrix
diag(a)<-NA
kdeValues<-c()
for (i in 1:dim(pvalMatrix)[1]) {
  tempData<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
  fhat<-kde(x = tempData,xmin =-1, xmax = 1, positive=F, h=h)
  kdeValues<-c(kdeValues,predict(fhat, x=c(0)))
}
kdeValues

# max possible KDE 
maxPossible<-1/(h*sqrt(2*pi))
maxPossible

# Plot the point with maximum kde value
library(fields) #To plot color bar
i<-which.max(kdeValues)
nLevel<-10
mirroredPval<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
fhat2<-kde(x = mirroredPval,xmin =-1, xmax = 1, positive=F,h = h)
plot(fhat2$eval.points,fhat2$estimate,type = "l", col=3,main=paste("Landmark", i, sep = " "),
     xlab="Mirrored p-values",ylab="KDE",ylim=c(0,maxPossible))
colfunc<-colorRampPalette(c("white","yellow","red","darkred"))
image.plot(legend.only=T, zlim=c(1,nLevel), col=colfunc(nLevel),legend.lab = "Level of significance")

nLevel<-10
shape4plot<-procSPHARM_PDM$mshape
open3d()
for (i in 1:numberOfPoints) {
  temp<-kdeValues[i]*nLevel/maxPossible
  myCol<-colfunc(nLevel)[ceiling(temp)]
  plot3d(shape4plot[c(i,i),],type="s", size=1.2,col = myCol,expand = 10,box=FALSE,add = TRUE)
}
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="black",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)


# Plot the point with maximum kde value
i<-which.max(kdeValues)
mirroredPval<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
fhat2<-kde(x = mirroredPval,xmin =-1, xmax = 1, positive=F,h = h)
plot(fhat2$eval.points,fhat2$estimate,type = "l", col=3,main=paste("Landmark", i, sep = " "),
     xlab="Mirrored p-values",ylab="KDE",ylim=c(0,maxPossible))
colfunc<-colorRampPalette(c("black","yellow"))
image.plot(legend.only=T, zlim=c(1,10), col=colfunc(2),legend.lab = "Level of significance")

#plot with threshold maxPossible/2
shape4plot<-procSPHARM_PDM$mshape
threshold<-maxPossible/2
open3d()
for (i in 1:numberOfPoints) {
  if(kdeValues[i]>threshold){
    myCol<-"yellow"
  }else{
    myCol<-"black"
  }
  plot3d(shape4plot[c(i,i),],type="s", size=1.2,col = myCol,expand = 10,box=FALSE,add = TRUE)
}
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="black",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)


#####################################################################################################
#####################################################################################################
# Hypothesis test independent from alignment for SPHARM-PDM without removing scale

# comment out lines below to calculate SPHARM-PDM distance tensors and p-value matrix!!!

# # takes some time !!!
# spharmTensorScaleFG1<-euclideanDistanceTensor(spharmPDM_G1)
# spharmTensorScaleFG2<-euclideanDistanceTensor(spharmPDM_G2)
# pvalMatrix_spharmScaleF<-calculatePvalMatrix(spharmTensorScaleFG1, spharmTensorScaleFG2)


load("../files/pvalMatrix_spharmScaleF.Rdata")


pvalMatrix_spharm<-pvalMatrix_spharmScaleF
pvalMatrix_spharm_BH<-adjustByUpperTriangle(pvalMatrix_spharm,method = "BH")
pvalMatrix_spharm_BH

# calculate kde of mirrored p-values
h<-0.05
pvalMatrix<-pvalMatrix_spharm
# pvalMatrix<-pvalMatrix_spharm_BH  #comment out to use adjusted p-value matrix
a<-pvalMatrix
diag(a)<-NA
kdeValues<-c()
for (i in 1:dim(pvalMatrix)[1]) {
  tempData<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
  fhat<-kde(x = tempData,xmin =-1, xmax = 1, positive=F, h=h)
  kdeValues<-c(kdeValues,predict(fhat, x=c(0)))
}
kdeValues

# max possible KDE 
maxPossible<-1/(h*sqrt(2*pi))
maxPossible

# Plot the point with maximum kde value
library(fields) #To plot color bar
i<-which.max(kdeValues)
nLevel<-10
mirroredPval<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
fhat2<-kde(x = mirroredPval,xmin =-1, xmax = 1, positive=F,h = h)
plot(fhat2$eval.points,fhat2$estimate,type = "l", col=3,main=paste("Landmark", i, sep = " "),
     xlab="Mirrored p-values",ylab="KDE",ylim=c(0,maxPossible))
colfunc<-colorRampPalette(c("white","yellow","red","darkred"))
image.plot(legend.only=T, zlim=c(1,nLevel), col=colfunc(nLevel),legend.lab = "Level of significance")

nLevel<-10
shape4plot<-procSPHARM_PDM$mshape
open3d()
for (i in 1:numberOfPoints) {
  temp<-kdeValues[i]*nLevel/maxPossible
  myCol<-colfunc(nLevel)[ceiling(temp)]
  plot3d(shape4plot[c(i,i),],type="s", size=1.2,col = myCol,expand = 10,box=FALSE,add = TRUE)
}
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="black",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)


# Plot the point with maximum kde value
i<-which.max(kdeValues)
mirroredPval<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
fhat2<-kde(x = mirroredPval,xmin =-1, xmax = 1, positive=F,h = h)
plot(fhat2$eval.points,fhat2$estimate,type = "l", col=3,main=paste("Landmark", i, sep = " "),
     xlab="Mirrored p-values",ylab="KDE",ylim=c(0,maxPossible))
colfunc<-colorRampPalette(c("black","yellow"))
image.plot(legend.only=T, zlim=c(1,10), col=colfunc(2),legend.lab = "Level of significance")

#plot with threshold maxPossible/2
shape4plot<-procSPHARM_PDM$mshape
threshold<-maxPossible/2
open3d()
for (i in 1:numberOfPoints) {
  if(kdeValues[i]>threshold){
    myCol<-"yellow"
  }else{
    myCol<-"black"
  }
  plot3d(shape4plot[c(i,i),],type="s", size=1.2,col = myCol,expand = 10,box=FALSE,add = TRUE)
}
verts <- rbind(t(as.matrix(procSPHARM_PDM$mshape)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
shade3d(tmesh, col="black",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)




#####################################################################################################
#####################################################################################################
# Generalized Shapiro-Wilk test
# Testing multivariate normality on aligned G1 and G2 SPHARM-PDMs

# align SPHARM-PDMs it takes a minute!
all_spharmPDM<-abind(spharmPDM_G1,spharmPDM_G2)
procSPHARM_PDM<-procGPA(all_spharmPDM, scale = F)
alignedSpharmPDM_G1<-procSPHARM_PDM$rotated[,,1:nSamplesG1]
alignedSpharmPDM_G2<-procSPHARM_PDM$rotated[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]

testOfNormalitySpharmG1<-c()
for (i in 1:numberOfPoints) {
  testOfNormalitySpharmG1<-c(testOfNormalitySpharmG1,
                             mvShapiro.Test(rbind(t(alignedSpharmPDM_G1[i,,])))$p.value)
}
testOfNormalitySpharmG2<-c()
for (i in 1:numberOfPoints) {
  testOfNormalitySpharmG2<-c(testOfNormalitySpharmG2,
                             mvShapiro.Test(rbind(t(alignedSpharmPDM_G2[i,,])))$p.value)
}

# percentage of points that failed in Shapiro-Wilk multivariate normality test
percentageOfNormalPointsG1<-sum(testOfNormalitySpharmG1<=0.05)/numberOfPoints
percentageOfNormalPointsG1  
percentageOfNormalPointsG2<-sum(testOfNormalitySpharmG2<=0.05)/numberOfPoints
percentageOfNormalPointsG2  

#####################################################################################################
#####################################################################################################
# Permutation as non-parametric approach for SPHARM-PDM size-and-shape
# The result is available in appendixA

# comment out below lines to recalculate the permutation. It takes some time!!!

# align SPHARM-PDMs it takes a minute!
all_spharmPDM<-abind(spharmPDM_G1,spharmPDM_G2)
procSPHARM_PDM<-procGPA(all_spharmPDM, scale = F) #scale=F for shape analysis
alignedSpharmPDM_G1<-procSPHARM_PDM$rotated[,,1:nSamplesG1]
alignedSpharmPDM_G2<-procSPHARM_PDM$rotated[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]

T0<-array(NA,dim = c(1,numberOfPoints))  

for (i in 1:numberOfPoints) {
  T0[1,i]<-MahalanobisDistance(k=3,nx=nSamplesG1,ny=nSamplesG2,
                               t(alignedSpharmPDM_G1[i,,]),t(alignedSpharmPDM_G2[i,,]))
}

# takes time!!!
pooledPositions<-abind(alignedSpharmPDM_G1,alignedSpharmPDM_G2)
nPerm<-10000
TdistanceMatrix<-array(NA,dim = c(numberOfPoints,nPerm))
for (j in 1:nPerm) {
  permNumbers<-sample(1:(nSamplesG1+nSamplesG2))
  groupA_elementNumbers<-permNumbers[c(1:nSamplesG1)]
  groupB_elementNumbers<-permNumbers[c((nSamplesG1+1):(nSamplesG1+nSamplesG2))]
  
  for (i in 1:numberOfPoints) {
    TdistanceMatrix[i,j]<-MahalanobisDistance(k=3,nx=nSamplesG1,ny=nSamplesG2,t(pooledPositions[i,,groupA_elementNumbers]),
                                              t(pooledPositions[i,,groupB_elementNumbers]))
  }
  
  # for (i in 1:numberOfPoints) {
  #   TdistanceMatrixSpharm1[i,j]<-euclideanDistance(colMeans(t(pooledPositions[i,,groupA_elementNumbers]))
  #                                           ,colMeans(t(pooledPositions[i,,groupB_elementNumbers])))
  # }
  
  cat("Permutation number",j,"is done.\n")
  
}

N<-dim(TdistanceMatrix)[2]
pValuesPermutation<-array(NA,dim = c(1,numberOfPoints))
for(i in 1:numberOfPoints){
  pValuesPermutation[1,i]<-(0.5+sum(TdistanceMatrix[i,]>=T0[1,i]))/(N+1)
}


pValuesPermutation_BH<-p.adjust(pValuesPermutation,method = "BH")

significantPvalues<-which(pValuesPermutation<=0.05)
significantPvalues
significantPvalues_BH<-which(pValuesPermutation_BH<=0.05)
significantPvalues_BH
significantPvalues_BH2<-which(pValuesPermutation_BH<=0.01)
significantPvalues_BH2

#plot without correction
shape4plot<-procSPHARM_PDM$mshape
open3d()
plot3d(shape4plot,type="s",  radius = 0.4,col = "#90a0ff",expand = 10,box=FALSE,add = TRUE)
plot3d(shape4plot[significantPvalues,],type="s", radius = 0.45,col = "red",expand = 10,box=FALSE,add = TRUE)
vertsG1 <- rbind(t(as.matrix(shape4plot)),1)
trglsG1 <- as.matrix(t(polyMatrix))
tmeshG1 <- tmesh3d(vertsG1, trglsG1)
shade3d(tmeshG1, col="#90a0ff",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#plot with correction FDR=0.05
open3d()
plot3d(shape4plot,type="s",  radius = 0.4,col = "#90a0ff",expand = 10,box=FALSE,add = TRUE)
plot3d(shape4plot[significantPvalues_BH,],type="s", radius = 0.45,col = "red",expand = 10,box=FALSE,add = TRUE)
vertsG1 <- rbind(t(as.matrix(shape4plot)),1)
trglsG1 <- as.matrix(t(polyMatrix))
tmeshG1 <- tmesh3d(vertsG1, trglsG1)
shade3d(tmeshG1, col="#90a0ff",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#plot with correction FDR=0.01
open3d()
plot3d(shape4plot,type="s",  radius = 0.4,col = "#90a0ff",expand = 10,box=FALSE,add = TRUE)
plot3d(shape4plot[significantPvalues_BH2,],type="s", radius = 0.45,col = "red",expand = 10,box=FALSE,add = TRUE)
vertsG1 <- rbind(t(as.matrix(shape4plot)),1)
trglsG1 <- as.matrix(t(polyMatrix))
tmeshG1 <- tmesh3d(vertsG1, trglsG1)
shade3d(tmeshG1, col="#90a0ff",alpha=1)
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)
