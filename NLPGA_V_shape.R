library(shapes)
library(RiemBase)
library(rgl)


#function to plot PNS circle on the sphere
drawCircleOnUnitSphere <- function(normalVec, radius, drawsphere) {
  if(sum(normalVec^2)<0.999 | sum(normalVec^2)>1.001){
    cat("Error: the normalVec is not normal!\n")
    break
  }
  if(sum(abs(normalVec)==c(1,0,0))==3 | sum(abs(normalVec)==c(0,1,0))==3 | sum(abs(normalVec)==c(0,0,1))==3){
    cat("Error: The orth axis is exact!\n")
    break
  }
  if(drawsphere==TRUE){
    spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
  }
  r<-radius
  v<-normalVec
  v3x<-v[1]
  v3y<-v[2]
  v3z<-v[3]
  # Calculate v1.
  a <- v3z
  b <- 0
  c <- -v3x
  v1x <- a/sqrt(a^2+b^2+c^2)
  v1y <- b/sqrt(a^2+b^2+c^2)
  v1z <- c/sqrt(a^2+b^2+c^2)
  # Calculate v2 as cross product of v3 and v1.
  # Since v1y is 0, it could be removed from the following calculations. Keeping it for consistency.
  v2x <- v3y * v1z - v3z * v1y
  v2y <- v3z * v1x - v3x * v1z
  v2z <- v3x * v1y - v3y * v1x
  # For each circle point.
  z<-sqrt(1^2-r^2)
  C<-v*z
  cx<-C[1]
  cy<-C[2]
  cz<-C[3]
  n<-300
  theta<-seq(0, 2*pi, len=n)
  px <- cx + r * (v1x * cos(theta) + v2x * sin(theta))
  py <- cy + r * (v1y * cos(theta) + v2y * sin(theta))
  pz <- cz + r * (v1z * cos(theta) + v2z * sin(theta))
  lines3d(px,py,pz)
}


#Generate points on Sphere S^2 
#v-shape
ndata = 1000  
theta = seq(from=-0.45,to=0.45,length.out=ndata)*1.3*pi 
tmpx = 3*cos(theta) + rnorm(ndata,sd=0.1)
tmpy = sin(theta) + rnorm(ndata,sd=0.1)
points<-c()
for (i in 1:ndata){
  tgt <- c(tmpx[i],0.5*tmpy[i],1)
  points<-rbind(points,tgt/sqrt(sum(tgt^2)))
}

#plot data on sphere
plot3d(points,type="p",expand=10, add=TRUE)
spheres3d(x=0,y=0,z=0,radius = 1,col="lightblue",alpha=0.1)

#Find Fre'chet mean
data <- list()
for (i in 1:dim(points)[1]) {
  data[[i]]<-points[i,]
}
data2 <- riemfactory(data, name="sphere")
out1 <- rbase.mean(data2)
FrechetMean<-t(out1$x) 
#Make sure Fre'chet mean is on the sphere
if(sum(FrechetMean^2)!=1){
  FrechetMean<-FrechetMean/sqrt(sum(FrechetMean^2))
}

#PNS and PNG
pnsSmall<-pns(t(points), sphere.type = "small")
pnsGreat<-pns(t(points), sphere.type = "great")
orthAxisGreat<-pnsGreat$PNS$orthaxis[[1]]
orthAxisSmall<-pnsSmall$PNS$orthaxis[[1]]
pnsMeanSmall<-pnsSmall$PNS$mean
pnsMeanGreat<-pnsGreat$PNS$mean

#plot PNS, PNG and Fre'chet mean
spheres3d(x=0,y=0,z=0,radius = 1,col="lightblue",alpha=0.1)
plot3d(points,type="p",expand=10, add=TRUE)
plot3d(rbind(FrechetMean,FrechetMean),type="s",col = "red",radius = 0.02,expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(pnsSmall$PNS$mean,pnsSmall$PNS$mean),type="s",col = "yellow",radius = 0.02,expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(pnsGreat$PNS$mean,pnsGreat$PNS$mean),type="s",col = "orange",radius = 0.02,expand = 10,box=TRUE,add = TRUE)
drawCircleOnUnitSphere(normalVec = orthAxisGreat,radius = pnsGreat$PNS$radii[2],drawsphere=TRUE)
drawCircleOnUnitSphere(normalVec = orthAxisSmall,radius = pnsSmall$PNS$radii[2],drawsphere=TRUE)

#map points to the tangent plane at the north pole
u<-t(points)
mu_g<-FrechetMean
R <- rotMat(mu_g,c(0,0,1))
shifted_u<-R%*%u
nSamples<-dim(u)[2]
theta<-acos(shifted_u[3,])
logPointx<-shifted_u[1,]*theta/sin(theta)
logPointy<-shifted_u[2,]*theta/sin(theta)
logPointz<-rep(1,nSamples)
logPoints<-cbind(logPointx,logPointy,logPointz)

#PGA i.e. PCA on tangent space 
log2D<-logPoints[,1:2]
dataCentered<-log2D-colMeans(log2D)  #centered data by PCA mean
pcaDataCentered<-prcomp(dataCentered)  #PCA
rotationMatrix<-pcaDataCentered$rotation  #matrix that its columns contain the eigenvectors
pcaPoints<-pcaDataCentered$x

#plot centered data on tangent space
plot(pcaPoints,col="green", xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2), xlab = "",ylab = "")

#curve fitting by rotating the PCA coordinate 
nPoints<-100
x1<-seq(from = -pi/2, to = pi/2,length.out = nPoints)
newDATA<-data.frame(x=x1)
polyDegree<-3  #choose the degree of the polynomial
gap<-10  #skip angle
allFirstEigenProp<-c()
for (k in seq(0,360,gap)) {
  tempTheta<-k*2*pi/360
  cat("Rotation angle theta:",tempTheta,"\n")
  planeRotationMatrix<-matrix(c(cos(tempTheta),-sin(tempTheta),sin(tempTheta),cos(tempTheta)),ncol = 2,byrow = T)
  pcaPointsTemp<-pcaPoints%*%planeRotationMatrix
  par(new=T)
  plot(pcaPointsTemp,col="green", xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2), xlab = "",ylab = "")
  
  x<-pcaPointsTemp[,1]
  y<-pcaPointsTemp[,2]
  fit <- lm(y ~ poly(x, degree = polyDegree,raw = TRUE),
            data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
  
  # ploynomial function
  f_fit <- function(x) {
    result<-0
    for (i in polyDegree:0) {
      result<-result+fit$coefficients[i+1]*x^i
    }
    return(result)
  }
  
  #distance function to calculate the distance between two points on the curve
  distanceFuncfit <- function(x) { 
    temp<-0
    for (i in polyDegree:1) {
      temp<-temp+i*fit$coefficients[i+1]*x^(i-1)
    }
    return(sqrt(1+temp^2))
  }
  
  #find minimum 
  curveY<-predict(fit, newdata = newDATA)
  par(new=T)
  plot(x1,curveY,type = "o", xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2), xlab = "",ylab = "")
  
  sumDistances<-rep(0, length(x1))
  for (i in 1:length(x1)) {
    for (j in 1:length(x)) {
      sumDistances[i]<-sumDistances[i]+abs(integrate(distanceFuncfit,x1[i],x[j])$value)
    }
  }
  indexMeanPoint<-which.min(sumDistances)
  glmMean<-c(x1[indexMeanPoint],f_fit(x1[indexMeanPoint]))
  
  par(new=T)
  plot(glmMean[1],glmMean[2] , pch=19,col="red", cex=2,
       xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2), xlab = "",ylab = "")
  
  glmResiduals1<-c()
  for (i in 1:length(x)) {
    glmResiduals1<-c(glmResiduals1,integrate(distanceFuncfit,glmMean[1],x[i])$value)
  }
  glmResiduals2<-fit$residuals
  glmResMatrix<-rbind(glmResiduals1,glmResiduals2)
  
  s<-prcomp(t(glmResMatrix),scale. = F)
  pScores<-s$x        #The coordinates of the observations on the principal components
  dim(s$x)
  eigenVectors<-s$rotation
  eigenValues<-(s$sdev)^2
  eigen_prop<- eigenValues/sum(eigenValues)*100
  
  cat("Percentage of eigenmodes contribution:",eigen_prop[1],"%",eigen_prop[2],"%.\n")
  allFirstEigenProp<-c(allFirstEigenProp,eigen_prop[1]) #we will rank performance by the first eigenmode contribution
}


#choose the best rotation and redo the procedure
max(allFirstEigenProp)
k<-which.max(allFirstEigenProp)
tempTheta<-seq(0,360,gap)[k]*2*pi/360   
planeRotationMatrix<-matrix(c(cos(tempTheta),-sin(tempTheta),sin(tempTheta),cos(tempTheta)),ncol = 2,byrow = T)
pcaPointsTemp<-pcaPoints%*%planeRotationMatrix
#find the equation of the best fitted polynomial
x<-pcaPointsTemp[,1]
y<-pcaPointsTemp[,2]
fit <- lm(y ~ poly(x, degree = polyDegree,raw = TRUE), data=as.data.frame(cbind(x,y))) #NB!!! raw=F calculate orthogonal polynomial
#function of the polynomial
f_fit <- function(x) { 
  result<-0
  for (i in polyDegree:0) {
    result<-result+fit$coefficients[i+1]*x^i
  }
  return(result)
}
#distance function to calculate the distance between two points on the curve
distanceFuncfit <- function(x) { 
  temp<-0
  for (i in polyDegree:1) {
    temp<-temp+i*fit$coefficients[i+1]*x^(i-1)
  }
  return(sqrt(1+temp^2))
}
#find mean
sumDistances<-rep(0, length(x1))
for (i in 1:length(x1)) {
  for (j in 1:length(x)) {
    sumDistances[i]<-sumDistances[i]+abs(integrate(distanceFuncfit,x1[i],x[j])$value)
  }
}
indexMeanPoint<-which.min(sumDistances)
glmMean<-c(x1[indexMeanPoint],f_fit(x1[indexMeanPoint]))

#plot the best fitted curve
curveY<-predict(fit, newdata = newDATA)
plot(pcaPointsTemp[,1],pcaPointsTemp[,2],col="green", 
     xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2),main = "NLPCA on tangent space",
     xlab = "x",ylab = "f(x)")
curveY<-predict(fit, newdata = newDATA)
par(new=T)
plot(x1,curveY,type = "o", xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2), xlab = "",ylab = "") #the same as predict

projectedData<-cbind(x,f_fit(x))
par(new=T)
plot(projectedData[,1],projectedData[,2] , pch=20,col="blue",
     xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2), xlab = "",ylab = "")
par(new=T)
plot(glmMean[1],glmMean[2] , pch=19,col="red", cex=2,
     xlim = c(-pi/2,pi/2),ylim = c(-pi/2,pi/2), xlab = "",ylab = "")
legend("bottomright", 
       legend = c("NLPCA mean", "Projected data", "Data", "Fitted curve"), 
       col = c("red", "blue","green","black"), 
       pch = c(19,20,1,1),
       lty = c(NA, NA, NA, 2),
       pt.cex = c(1.5,1.5,1,1), 
       cex = 0.9, 
       text.col = "black", 
       horiz = F)


#NLPCA residuals
glmResiduals1<-c()
for (i in 1:length(x)) {
  glmResiduals1<-c(glmResiduals1,integrate(distanceFuncfit,glmMean[1],x[i])$value)
}
glmResiduals2<-fit$residuals
glmResMatrix<-rbind(glmResiduals1,glmResiduals2)

#map NLPCA mean from tangent space to the sphere 
rotatedBackGlmMean<-(glmMean%*%solve(planeRotationMatrix))%*%solve(rotationMatrix)+colMeans(log2D)
v1<-rotatedBackGlmMean[1]
v2<-rotatedBackGlmMean[2]
sizeV<-sqrt(v1^2+v2^2)
expMapGlmMean<-c(v1*sin(sizeV)/sizeV, v2*sin(sizeV)/sizeV, cos(sizeV))
rotatedBackExpMapGlmMean<-expMapGlmMean%*%R

#map data from tangent space to the sphere 
rotatedBackProjectedData<-(projectedData%*%solve(planeRotationMatrix))%*%solve(rotationMatrix)+colMeans(log2D)
w1<-rotatedBackProjectedData[,1]
w2<-rotatedBackProjectedData[,2]
sizeW<-sqrt(w1^2+w2^2)
expMapProjectedData<-rbind(w1*sin(sizeW)/sizeW, w2*sin(sizeW)/sizeW, cos(sizeW))
rotatedBackExpMapprojectedData<-t(expMapProjectedData)%*%R

#plot PNS, PNG, Fre'chet and NLPCA mean
open3d()
spheres3d(x = 0, y = 0, z = 0, radius = 1,col = "lightblue", alpha=0.1)
plot3d(t(u),type="p",expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(FrechetMean,FrechetMean),type="s",col = "red",radius = 0.02,expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(pnsSmall$PNS$mean,pnsSmall$PNS$mean),type="s",col = "yellow",radius = 0.02,expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(pnsGreat$PNS$mean,pnsGreat$PNS$mean),type="s",col = "orange",radius = 0.02,expand = 10,box=TRUE,add = TRUE)
plot3d(rbind(rotatedBackExpMapGlmMean,rotatedBackExpMapGlmMean),type="s",col ="#1E90FF",radius = 0.02,expand = 10,box=TRUE,add = TRUE)
plot3d(rotatedBackExpMapprojectedData,type="p",col = "blue",expand = 10,box=TRUE,add = TRUE)
drawCircleOnUnitSphere(normalVec = orthAxisGreat,radius = pnsGreat$PNS$radii[2],drawsphere=TRUE)
drawCircleOnUnitSphere(normalVec = orthAxisSmall,radius = pnsSmall$PNS$radii[2],drawsphere=TRUE)

#compare performance
#PNS small circle
s_CG1<-prcomp(t(pnsSmall$resmat),scale. = F)
pScores_CG1<-s_CG1$x        #The coordinates of the observations on the principal components
eigenVectors_CG1<-s_CG1$rotation
eigenValues_CG1<-(s_CG1$sdev)^2
eigen_prop_CG1<- eigenValues_CG1/sum(eigenValues_CG1)*100
xx1=barplot(eigen_prop_CG1[1:3],xlab = 'No. of Eigenmodes=20', ylab = 'Percentage contribution', col = "blue",
            main= "Plot of Eigenmode contributions for CPNG control group ",border = par("fg"),ylim =c(0,110))
y_CG1<-round(eigen_prop_CG1[1:3], digits = 2)
text(x=xx1,y=y_CG1,labels=as.character(y_CG1),pos = 3 , col = "blue", cex = 0.9)
cum_prop_CG1 <- cumsum(eigenValues_CG1)/sum(eigenValues_CG1)*100
lines(cum_prop_CG1)

#PNG
s_CG2<-prcomp(t(pnsGreat$resmat),scale. = F)
pScores_CG2<-s_CG2$x        #The coordinates of the observations on the principal components
dim(s_CG2$x)
eigenVectors_CG2<-s_CG2$rotation
eigenValues_CG2<-(s_CG2$sdev)^2
eigen_prop_CG2<- eigenValues_CG2/sum(eigenValues_CG2)*100
xx2=barplot(eigen_prop_CG2[1:3],xlab = 'No. of Eigenmodes=20', ylab = 'Percentage contribution', col = "blue",
            main= "Plot of Eigenmode contributions for CPNG control group ",border = par("fg"),ylim =c(0,110))
y_CG2<-round(eigen_prop_CG2[1:3], digits = 2)
text(x=xx2,y=y_CG2,labels=as.character(y_CG2),pos = 3 , col = "blue", cex = 0.9)
cum_prop_CG2 <- cumsum(eigenValues_CG2)/sum(eigenValues_CG2)*100
lines(cum_prop_CG2)

#NLPCA
s_CG3<-prcomp(t(glmResMatrix),scale. = F)
pScores_CG3<-s_CG3$x        #The coordinates of the observations on the principal components
dim(s_CG3$x)
eigenVectors_CG3<-s_CG3$rotation
eigenValues_CG3<-(s_CG3$sdev)^2
eigen_prop_CG3<- eigenValues_CG3/sum(eigenValues_CG3)*100
xx3=barplot(eigen_prop_CG3[1:3],xlab = 'No. of Eigenmodes=20', ylab = 'Percentage contribution', col = "blue",
            main= "Plot of Eigenmode contributions for CPNG control group ",border = par("fg"),ylim =c(0,110))
y_CG3<-round(eigen_prop_CG3[1:3], digits = 2)
text(x=xx3,y=y_CG3,labels=as.character(y_CG3),pos = 3 , col = "blue", cex = 0.9)
cum_prop_CG3 <- cumsum(eigenValues_CG3)/sum(eigenValues_CG3)*100
lines(cum_prop_CG3)

#plot the performance comparison
library(dplyr)
library(ggplot2)
df2 <- data.frame(Method=rep(c("Small cirlce", "NLPCA","Great circle"), each=3),
                  Eigenmodes=rep(c(1,2,3),3),
                  Percentage_Contribution=c(cum_prop_CG1,100,cum_prop_CG3,100,cum_prop_CG2,100))
head(df2)

ggplot(df2, aes(x=Eigenmodes, y=Percentage_Contribution, group=Method)) +
  geom_line(aes(color=Method))+
  geom_point(aes(color=Method))+
  ggtitle("Performance comparison") + 
  theme(plot.title = element_text( size = 12,hjust = 0.5)) +
  xlab("Eigenmodes") + ylab("Percentage contribution")








