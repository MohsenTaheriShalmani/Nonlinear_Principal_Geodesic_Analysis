library(shapes)
library(rgl)
library(RiemBase)

#####################################################################################################
#####################################################################################################
#set working directory to file location

setwd(dirname(rstudioapi::getSourceEditorContext()$path))

#####################################################################################################
#####################################################################################################
# functions 

# calculate unit normal vector of triangle mesh
unitNormalVec <- function(point1,point2,point3) {
  a<-point2-point1
  b<-point3-point1
  
  normalVec<-c((a[2]*b[3]-a[3]*b[2]),-(a[1]*b[3]-a[3]*b[1]),(a[1]*b[2]-a[2]*b[1]))
  triangleArea<-sqrt(sum(normalVec^2))/2
  unitNormal<-normalVec/sqrt(sum(normalVec^2))
  
  my_list <- list("unitNormal" = unitNormal, "trianleArea" = triangleArea)
  
  return(my_list)
  
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

#####################################################################################################
#####################################################################################################
# read SPHARM-PDM and s-rep data 

#choose two groups from {"CG", "PD", "PDD1", "PDD2", "PDD3", "PDD4"} to compare (e.g. "CG" and "PD")
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
nSamplesG2<-dim(meshPoints_G2)[2]

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
# rigid Procrustes alignment

#samples are already aligned by SPHARM-PDM toolbox so we can skip this part
# takes a minute !!!
procSpharm<-procGPA(abind(spharmPDM_G1,spharmPDM_G2),scale = F)
spharmPDM_G1<-procSpharm$rotated[,,1:nSamplesG1]
spharmPDM_G2<-procSpharm$rotated[,,(nSamplesG1+1):(nSamplesG1+nSamplesG2)]


#####################################################################################################
#####################################################################################################
# find unit normal vectors and triangle size of meshes 

unitNormalMesh_G1<-array(NA,dim = c(dim(polyMatrix)[1],3,nSamplesG1))
triangleSizeMesh_G1<-array(NA,dim = c(dim(polyMatrix)[1],nSamplesG1))
for (k in 1:nSamplesG1) {
  for (i in 1:dim(PolygonsCsv)[1]) {
    p1<-polyMatrix[i,1]
    p2<-polyMatrix[i,2]
    p3<-polyMatrix[i,3]
    point1<-spharmPDM_G1[p1,,k]
    point2<-spharmPDM_G1[p2,,k]
    point3<-spharmPDM_G1[p3,,k]
    tempNormalInfo<-unitNormalVec(point1,point2,point3)
    unitNormalMesh_G1[i,,k]<-tempNormalInfo$unitNormal
    triangleSizeMesh_G1[i,k]<-tempNormalInfo$trianleArea
  }
}

unitNormalMesh_G2<-array(NA,dim = c(dim(polyMatrix)[1],3,nSamplesG2))
triangleSizeMesh_G2<-array(NA,dim = c(dim(polyMatrix)[1],nSamplesG2))
for (k in 1:nSamplesG2) {
  for (i in 1:dim(PolygonsCsv)[1]) {
    p1<-polyMatrix[i,1]
    p2<-polyMatrix[i,2]
    p3<-polyMatrix[i,3]
    point1<-spharmPDM_G2[p1,,k]
    point2<-spharmPDM_G2[p2,,k]
    point3<-spharmPDM_G2[p3,,k]
    tempNormalInfo<-unitNormalVec(point1,point2,point3)
    unitNormalMesh_G2[i,,k]<-tempNormalInfo$unitNormal
    triangleSizeMesh_G2[i,k]<-tempNormalInfo$trianleArea
  }
}


#plot i-th correponding vector of two groups (e.g. vector no. 127)
i<-127 #change this to see other vectors
spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
plot3d(t(unitNormalMesh_G1[i,,]),type="p",col = "blue",expand = 10,box=TRUE,add = TRUE)
plot3d(t(unitNormalMesh_G2[i,,]),type="p",col = "red",expand = 10,box=TRUE,add = TRUE)

# plot mesh and normal vectors of a sample
sampleNo<-1  #change this to see normal vectors of other samples 
verts <- rbind(t(as.matrix(spharmPDM_G1[,,sampleNo])),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
wire3d(tmesh, col="blue")  #wire mesh
# shade3d(tmesh, col="blue")  #surface mech
for (i in 1:dim(PolygonsCsv)[1]) {
  p1<-polyMatrix[i,1]
  p2<-polyMatrix[i,2]
  p3<-polyMatrix[i,3]
  point1<-spharmPDM_G1[p1,,sampleNo]
  point2<-spharmPDM_G1[p2,,sampleNo]
  point3<-spharmPDM_G1[p3,,sampleNo]
  triangleMean<-(point1+point2+point3)/3
  triangleVecTip<-triangleMean+unitNormalVec(point1,point2,point3)$unitNormal
  
  plot3d(rbind(triangleMean,triangleVecTip),type="l",col = "red",expand = 10,box=FALSE,add = TRUE)
}


#####################################################################################################
#####################################################################################################
# test for unit normal vectors

pValues_Hotelling_SPHARM_Directions<-c()
for(i in 1:dim(polyMatrix)[1]){
  
  allDirTemp<-t(cbind(unitNormalMesh_G1[i,,],unitNormalMesh_G2[i,,]))
  data1 <- list()
  for (j in 1:dim(allDirTemp)[1]){
    data1[[j]] <-allDirTemp[j,]
  }
  data2 <- riemfactory(data1, name="sphere")
  ### Compute Fre'chet Mean
  out1<- rbase.mean(data2)
  mu_g<-as.vector(out1$x)
  
  R <- rotMat(mu_g,c(0,0,1))
  shiftedG1i<-R%*%unitNormalMesh_G1[i,,]
  shiftedG2i<-R%*%unitNormalMesh_G2[i,,]
  
  theta1<-acos(shiftedG1i[3,])
  logPointG1ix<-shiftedG1i[1,]*theta1/sin(theta1)
  logPointG1iy<-shiftedG1i[2,]*theta1/sin(theta1)
  logPointG1iz<-rep(1,nSamplesG1)
  logG1i<-cbind(logPointG1ix,logPointG1iy,logPointG1iz)
  
  theta2<-acos(shiftedG2i[3,])
  logPointG2ix<-shiftedG2i[1,]*theta2/sin(theta2)
  logPointG2iy<-shiftedG2i[2,]*theta2/sin(theta2)
  logPointG2iz<-rep(1,nSamplesG2)
  logG2i<-cbind(logPointG2ix,logPointG2iy,logPointG2iz)
  
  pValues_Hotelling_SPHARM_Directions<-c(pValues_Hotelling_SPHARM_Directions,HotellingT2(logG1i[,1:2],logG2i[,1:2]))
  
  cat(i,"vectors of",dim(polyMatrix)[1],"is done.\n")
}


alpha<-0.05
FDR<-0.15

sum(pValues_Hotelling_SPHARM_Directions<=alpha)
pValues_Hotelling_SPHARM_Directions_BH<-p.adjust(pValues_Hotelling_SPHARM_Directions, method = "BH")
sum(pValues_Hotelling_SPHARM_Directions_BH<=FDR)

sigNormals<-which(pValues_Hotelling_SPHARM_Directions<=alpha)
sigNormals
sigNormals_BH<-which(pValues_Hotelling_SPHARM_Directions_BH<=FDR)
sigNormals_BH


#plot significant normals without correction
shape4plot<-procSpharm$mshape
verts <- rbind(t(as.matrix(shape4plot)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
open3d()
shade3d(tmesh, col="blue",alpha=0.5)
verts2 <- rbind(t(as.matrix(shape4plot)),1)
trgls2 <- as.matrix(t(polyMatrix[sigNormals,]))
tmesh2 <- tmesh3d(verts2, trgls2)
shade3d(tmesh2, col="orange",alpha=1)
for (i in sigNormals) {
  p1<-polyMatrix[i,1]
  p2<-polyMatrix[i,2]
  p3<-polyMatrix[i,3]
  point1<-shape4plot[p1,]
  point2<-shape4plot[p2,]
  point3<-shape4plot[p3,]
  triangleMean<-(point1+point2+point3)/3
  triangleVecTip<-triangleMean+unitNormalVec(point1,point2,point3)$unitNormal
  plot3d(rbind(triangleMean,triangleVecTip),type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
}
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#plot significant normals with FDR correction
shape4plot<-procSpharm$mshape
verts <- rbind(t(as.matrix(shape4plot)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
open3d()
shade3d(tmesh, col="blue",alpha=0.5)
verts2 <- rbind(t(as.matrix(shape4plot)),1)
trgls2 <- as.matrix(t(polyMatrix[sigNormals_BH,]))
tmesh2 <- tmesh3d(verts2, trgls2)
shade3d(tmesh2, col="orange",alpha=1)
for (i in sigNormals_BH) {
  p1<-polyMatrix[i,1]
  p2<-polyMatrix[i,2]
  p3<-polyMatrix[i,3]
  point1<-shape4plot[p1,]
  point2<-shape4plot[p2,]
  point3<-shape4plot[p3,]
  triangleMean<-(point1+point2+point3)/3
  triangleVecTip<-triangleMean+unitNormalVec(point1,point2,point3)$unitNormal
  plot3d(rbind(triangleMean,triangleVecTip),type="l",lwd = 3,col = "red",expand = 10,box=FALSE,add = TRUE)
}
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#####################################################################################################
#####################################################################################################
# test for triangle size

pValues_Ttest_meshTriangleSizes<-c()
for(i in 1:dim(polyMatrix)[1]){
  
  pValues_Ttest_meshTriangleSizes<-c(pValues_Ttest_meshTriangleSizes,
                                     t.test(triangleSizeMesh_G1[i,],triangleSizeMesh_G2[i,])$p.value)
}

alpha<-0.05
FDR<-0.15

sum(pValues_Ttest_meshTriangleSizes<=alpha)
pValues_Ttest_meshTriangleSizes_BH<-p.adjust(pValues_Ttest_meshTriangleSizes, method = "BH")
sum(pValues_Ttest_meshTriangleSizes_BH<=FDR)

sigSize<-which(pValues_Ttest_meshTriangleSizes<=alpha)
sigSize
sigSize_BH<-which(pValues_Ttest_meshTriangleSizes_BH<=FDR)
sigSize_BH

#plot significant normals without correction
verts <- rbind(t(as.matrix(shape4plot)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
open3d()
shade3d(tmesh, col="blue",alpha=0.03)
verts2 <- rbind(t(as.matrix(shape4plot)),1)
trgls2 <- as.matrix(t(polyMatrix[sigSize,]))
tmesh2 <- tmesh3d(verts2, trgls2)
shade3d(tmesh2, col="red")
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)

#plot significant normals with FDR correction
verts <- rbind(t(as.matrix(shape4plot)),1)
trgls <- as.matrix(t(polyMatrix))
tmesh <- tmesh3d(verts, trgls)
open3d()
shade3d(tmesh, col="blue",alpha=0.03)
verts2 <- rbind(t(as.matrix(shape4plot)),1)
trgls2 <- as.matrix(t(polyMatrix[sigSize_BH,]))
tmesh2 <- tmesh3d(verts2, trgls2)
shade3d(tmesh2, col="red")
decorate3d(xlab = "x", ylab = "y", zlab = "z",
           xlim = c(-10,15),ylim =c(-20,20),zlim = c(-5,5),
           box = TRUE, axes = TRUE, main = NULL, sub = NULL,
           top = TRUE, aspect = FALSE, expand = 1.1)



