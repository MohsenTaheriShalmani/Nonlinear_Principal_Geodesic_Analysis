library(shapes)
library(rgl)

#####################################################################################################
#####################################################################################################
# functions

# Euclidean Distance
euclideanDistance_2D<-function(x,y){
  return(sqrt((x[1]-y[1])^2+(x[2]-y[2])^2))
}

# Angle between 2 vectors in 2D between [0,2pi]
angleBetween3points_2D <- function(vertex,point2,point3) {
  
  v1<-point2-vertex
  v2<-point3-vertex
  
  dot <- v1[1]*v2[1] + v1[2]*v2[2]      # dot product
  det <- v1[1]*v2[2] - v1[2]*v2[1]      # determinant
  angle <- atan2(-det, -dot)+pi         #atan2(sin, cos)
  
  return(angle)
  
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

# function to adjust upper triangle portion of the p-value matrix
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

#####################################################################################################
#####################################################################################################
# Hypothesis test HotellingT2 without removing scale (size-and-shape) 

group1<-panf.dat
group2<-panm.dat

#plot both groups 
plotshapes(abind(group1,group2))

# without removing scale for size-and-shape analysis
pooledGroup<-abind(group1,group2)
procPooled<-procGPA(pooledGroup, scale = F)
group1<-procPooled$rotated[,,1:dim(group1)[3]]
group2<-procPooled$rotated[,,(dim(group1)[3]+1):(dim(group1)[3]+dim(group2)[3])]

#plot aligned groups
par(new=F)
xlim <- ylim <- c(-120,120)
plot(1, type="n", main = "Overlaid size-and-shapes" ,
     xlim = xlim, ylim = ylim,xlab = "x", ylab = "y")
par(new=T)
for (i in 1:dim(group1)[3]) {
  plot(group1[,,i], pch=20, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.4) ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y",cex=2)
  par(new=T)
}
for (i in 1:dim(group2)[3]) {
  plot(group2[,,i], pch=20, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.4),
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y",cex=2)
  par(new=T)
}
legend("bottomright", legend=c("panf", "panm"),pch = c(19,19),
       col=c("blue","red"), cex=0.6,pt.cex = 1.5, horiz=F)

# HotellingT2 test
pvalues_Hotelling<-c()
for (i in 1:dim(group1)[1]) {
  pvalues_Hotelling<-c(pvalues_Hotelling,HotellingT2(t(group1[i,,]),t(group2[i,,])))
}
pvalues_Hotelling

# significant landmarks
alpha<-0.05
FDR<-0.15
significant<-which(pvalues_Hotelling<=alpha)
significant
significant_BH<-which(p.adjust(pvalues_Hotelling,method = "BH")<=FDR)
significant_BH

#plot significants
shape4plot<-procPooled$mshape
par(new=F)
plot(shape4plot, pch=3, main = "Significant Landmarks",
     xlim = xlim, ylim = ylim,
     xlab = "x", ylab = "y",cex=1)
par(new=T)
plot(shape4plot[significant_BH,], pch=10,
     xlim = xlim, ylim = ylim,
     xlab = "x", ylab = "y",cex=2)
legend("bottomright", legend=c("Non-significant", "Significant"), pch = c(3,10),
        cex=0.55,pt.cex = 1, horiz=T)

#####################################################################################################
#####################################################################################################
# Hypothesis test with HotellingT^2 and removing scale

group1<-panf.dat
group2<-panm.dat

# Remove scale for shape analysis
pooledGroup<-abind(group1,group2)
procPooled<-procGPA(pooledGroup, scale = T)
group1<-procPooled$rotated[,,1:dim(group1)[3]]
group2<-procPooled$rotated[,,(dim(group1)[3]+1):(dim(group1)[3]+dim(group2)[3])]

#plot aligned groups
par(new=F)
xlim <- ylim <- c(-120,120)
plot(1, type="n", main = "Overlaid shapes" ,
     xlim = xlim, ylim = ylim,xlab = "x", ylab = "y")
par(new=T)
for (i in 1:dim(group1)[3]) {
  plot(group1[,,i], pch=20, col=rgb(red = 0, green = 0, blue = 1, alpha = 0.4) ,
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y",cex=2)
  par(new=T)
}
for (i in 1:dim(group2)[3]) {
  plot(group2[,,i], pch=20, col=rgb(red = 1, green = 0, blue = 0, alpha = 0.4),
       xlim = xlim, ylim = ylim,
       xlab = "x", ylab = "y",cex=2)
  par(new=T)
}
legend("bottomright", legend=c("panf", "panm"),pch = c(19,19),
       col=c("blue","red"), cex=0.6,pt.cex = 1.5, horiz=F)

# HotellingT2 test
pvalues_Hotelling<-c()
for (i in 1:dim(group1)[1]) {
  pvalues_Hotelling<-c(pvalues_Hotelling,HotellingT2(t(group1[i,,]),t(group2[i,,])))
}
pvalues_Hotelling

#significant landmarks
alpha<-0.05
FDR<-0.15
significant<-which(pvalues_Hotelling<=alpha)
significant
significant_BH<-which(p.adjust(pvalues_Hotelling,method = "BH")<=FDR)
significant_BH

#plot significants
shape4plot<-procPooled$mshape
par(new=F)
plot(shape4plot, pch=3, main = "Significant Landmarks",
     xlim = xlim, ylim = ylim,
     xlab = "x", ylab = "y",cex=1)
par(new=T)
plot(shape4plot[significant_BH,], pch=10,
     xlim = xlim, ylim = ylim,
     xlab = "x", ylab = "y",cex=2)
legend("bottomright", legend=c("Non-significant", "Significant"), pch = c(3,10),
       cex=0.55,pt.cex = 1, horiz=T)

#####################################################################################################
#####################################################################################################
# Hypothesis test with distance, without removing scale

group1<-panf.dat
group2<-panm.dat

# calculate 3d tensors of distances
q<-dim(group1)[1]
m<-dim(group1)[2]
n1<-dim(group1)[3]
tensor1<-array(NA, dim=c(q,q,n1))
for (k in 1:n1) {
  for (i in 1:q) {
    for (j in 1:q) {
      tensor1[i,j,k]<-euclideanDistance_2D(group1[i,,k],group1[j,,k])
    }
  }
}
n2<-dim(group2)[3]
tensor2<-array(NA, dim=c(q,q,n2))
for (k in 1:n2) {
  for (i in 1:q) {
    for (j in 1:q) {
      tensor2[i,j,k]<-euclideanDistance_2D(group2[i,,k],group2[j,,k])
    }
  }
}

# p-value matrix
pvalMatrix<-array(NA, dim=c(q,q))
for (i in 1:q) {
  for (j in 1:q) {
    if(i!=j){
      pvalMatrix[i,j]<-t.test(tensor1[i,j,],tensor2[i,j,])$p.value
    }else{
      pvalMatrix[i,j]<-1
    }
  }
}

pvalMatrix
#adjusted p-value matrix
pvalMatrix_BH<-adjustByUpperTriangle(pvalMatrix,method = "BH")
pvalMatrix_BH

# calculate kde of mirrored p-values
library(ks)
h<-0.05
a<-pvalMatrix
# a<-pvalMatrix_BH #comment out to use adjusted p-value matrix
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

#plot landmarks' KDE
par(mfrow=c(2,4))
h=0.05
for (i in 1:dim(group1)[1]) {
  mirroredPval<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
  fhat2<-kde(x = mirroredPval,xmin =-1, xmax = 1, positive=F,h = h)
  plot(fhat2$eval.points,fhat2$estimate,type = "l", col=3,main=paste("Landmark", i, sep = " "),
       xlab="Mirrored p-values",ylab="KDE",ylim=c(0,maxPossible))
}
par(mfrow=c(1,1))

# plot significant landmark by the size of kde
library(fields) #To plot color bar
colfunc<-colorRampPalette(c("red","yellow"))
nLevel<-10
shape4plot<-procGPA(abind(group1,group2))$mshape
for (i in 1:dim(group1)[1]) {
  temp<-kdeValues[i]*nLevel/maxPossible
  myCol<-colfunc(nLevel)[ceiling(temp)]
  plot(shape4plot[i,1],shape4plot[i,2], pch=19, 
       col=myCol,
       xlim = xlim, ylim = ylim,
       xlab = "", ylab = "",cex=2)
  par(new=T)
}
image.plot(legend.only=T, zlim=c(1,nLevel), col=colfunc(nLevel),legend.lab = "Level of significance")


#####################################################################################################
#####################################################################################################
# Hypothesis test with distance and removing scale

group1<-panf.dat
group2<-panm.dat

# Remove scale for shape analysis
pooledGroup<-abind(group1,group2)
procPooled<-procGPA(pooledGroup, scale = T)
group1<-procPooled$rotated[,,1:dim(group1)[3]]
group2<-procPooled$rotated[,,(dim(group1)[3]+1):(dim(group1)[3]+dim(group2)[3])]

# calculate 3d tensors of distances
q<-dim(group1)[1]
m<-dim(group1)[2]
n1<-dim(group1)[3]
tensor1<-array(NA, dim=c(q,q,n1))
for (k in 1:n1) {
  for (i in 1:q) {
    for (j in 1:q) {
      tensor1[i,j,k]<-euclideanDistance_2D(group1[i,,k],group1[j,,k])
    }
  }
}
n2<-dim(group2)[3]
tensor2<-array(NA, dim=c(q,q,n2))
for (k in 1:n2) {
  for (i in 1:q) {
    for (j in 1:q) {
      tensor2[i,j,k]<-euclideanDistance_2D(group2[i,,k],group2[j,,k])
    }
  }
}

# p-value matrix
pvalMatrix<-array(NA, dim=c(q,q))
for (i in 1:q) {
  for (j in 1:q) {
    if(i!=j){
      pvalMatrix[i,j]<-t.test(tensor1[i,j,],tensor2[i,j,])$p.value
    }else{
      pvalMatrix[i,j]<-1
    }
  }
}

pvalMatrix
#adjusted p-value matrix
pvalMatrix_BH<-adjustByUpperTriangle(pvalMatrix,method = "BH")
pvalMatrix_BH

# calculate kde of mirrored p-values
library(ks)
h<-0.05
a<-pvalMatrix
# a<-pvalMatrix_BH #comment out to use adjusted p-value matrix
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

#plot landmarks' KDE
par(mfrow=c(2,4))
h=0.05
for (i in 1:dim(group1)[1]) {
  mirroredPval<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
  fhat2<-kde(x = mirroredPval,xmin =-1, xmax = 1, positive=F,h = h)
  plot(fhat2$eval.points,fhat2$estimate,type = "l", col=3,main=paste("Landmark", i, sep = " "),
       xlab="Mirrored p-values",ylab="KDE",ylim=c(0,maxPossible))
}
par(mfrow=c(1,1))

# plot significant landmark by the size of kde
library(fields) #To plot color bar
colfunc<-colorRampPalette(c("red","yellow"))
nLevel<-10
shape4plot<-procGPA(abind(group1,group2))$mshape
for (i in 1:dim(group1)[1]) {
  temp<-kdeValues[i]*nLevel/maxPossible
  myCol<-colfunc(nLevel)[ceiling(temp)]
  plot(shape4plot[i,1],shape4plot[i,2], pch=19, 
       col=myCol,
       xlim = xlim, ylim = ylim,
       xlab = "", ylab = "",cex=2)
  par(new=T)
}
image.plot(legend.only=T, zlim=c(1,nLevel), col=colfunc(nLevel),legend.lab = "Level of significance")

#####################################################################################################
#####################################################################################################
# Hypothesis test with angle

group1<-panf.dat
group2<-panm.dat

#angles of a group of shapes
anglesMatrixG1<-array(NA,dim=c(dim(group1)[3],3*choose(dim(group1)[1],3)))
for (t in 1:dim(group1)[3]) {
  tempShape<-group1[,,t]
  points<-c(1:dim(tempShape)[1])
  h<-1
  for (i in points) {
    tips<-points[!points %in% i]
    constantTips<-c()
    for (j in tips) {
      constantTip<-j
      constantTips<-c(constantTips,constantTip)
      remainTips<-tips[!tips %in% constantTips]
      for (k in remainTips) {
        tempAngle<-angleBetween3points_2D(vertex = tempShape[i,],
                                          point2 = tempShape[constantTip,],
                                          point3 = tempShape[k,])
        anglesMatrixG1[t,h]<-tempAngle
        h<-h+1
      }
      constantTips<-c(constantTips,constantTip)
    }
  }
}

anglesMatrixG2<-array(NA,dim=c(dim(group2)[3],3*choose(dim(group2)[1],3)))
for (t in 1:dim(group2)[3]) {
  tempShape<-group2[,,t]
  points<-c(1:dim(tempShape)[1])
  h<-1
  for (i in points) {
    tips<-points[!points %in% i]
    constantTips<-c()
    for (j in tips) {
      constantTip<-j
      constantTips<-c(constantTips,constantTip)
      remainTips<-tips[!tips %in% constantTips]
      for (k in remainTips) {
        tempAngle<-angleBetween3points_2D(vertex = tempShape[i,],
                                          point2 = tempShape[constantTip,],
                                          point3 = tempShape[k,])
        anglesMatrixG2[t,h]<-tempAngle
        h<-h+1
      }
      constantTips<-c(constantTips,constantTip)
    }
  }
}

library(RiemBase)
pvalueVector<-rep(NA,dim(anglesMatrixG1)[2])
for (i in 1:dim(anglesMatrixG1)[2]) {
  
  angle2directionG1<-cbind(cos(anglesMatrixG1[,i]),sin(anglesMatrixG1[,i]))
  angle2directionG2<-cbind(cos(anglesMatrixG2[,i]),sin(anglesMatrixG2[,i]))
  allDirTemp<-rbind(angle2directionG1,angle2directionG2)
  data1 <- list()
  for (j in 1:dim(allDirTemp)[1]){
    data1[[j]] <-allDirTemp[j,]
  }
  data2 <- riemfactory(data1, name="sphere")
  ### Compute Fre'chet Mean
  out1<- rbase.mean(data2)
  mu_g<-as.vector(out1$x)
  mu_g
  
  R <- rotMat(mu_g,c(0,1))
  shifted1<-t(R%*%t(angle2directionG1))
  shifted2<-t(R%*%t(angle2directionG2))
  
  theta1<-acos(shifted1[,2])
  logPointG1x<-shifted1[,1]*theta1/sin(theta1)
  logPointG1y<-rep(1,dim(shifted1)[1])
  theta2<-acos(shifted2[,2])
  logPointG2x<-shifted2[,1]*theta2/sin(theta2)
  logPointG2y<-rep(1,dim(shifted2)[1])
  
  pvalueVector[i]<-t.test(logPointG1x,logPointG2x)$p.value
  
}


pvalMatrix<-matrix(pvalueVector,nrow =dim(group1)[1] ,byrow = T)
pvalMatrix_BH<-matrix(p.adjust(pvalueVector,method = "BH"),nrow =dim(group1)[1] ,byrow = T)

# calculate kde of mirrored p-values
library(ks)
h<-0.05
a<-pvalMatrix
# a<-pvalMatrix_BH #comment out to use adjusted p-value matrix
diag(a)<-NA
kdeValues<-c()
kdeValuesScaled<-c()
for (i in 1:dim(pvalMatrix)[1]) {
  tempData<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
  fhat<-kde(x = tempData,xmin =-1, xmax = 1, positive=F, h=h)
  kdeValues<-c(kdeValues,predict(fhat, x=c(0)))
}
kdeValues

# max possible KDE 
maxPossible<-1/(h*sqrt(2*pi))
maxPossible

#plot landmarks' KDE
par(mfrow=c(2,4))
h=0.05
for (i in 1:dim(group1)[1]) {
  mirroredPval<-c(a[i,][which(!is.na(a[i,]))],-(a[i,][which(!is.na(a[i,]))]))
  fhat2<-kde(x = mirroredPval,xmin =-1, xmax = 1, positive=F,h = h)
  plot(fhat2$eval.points,fhat2$estimate,type = "l", col=3,main=paste("Landmark", i, sep = " "),
       xlab="Mirrored p-values",ylab="KDE",ylim=c(0,maxPossible))
}
par(mfrow=c(1,1))

# plot significant landmark by the size of kde
library(fields) #To plot color bar
colfunc<-colorRampPalette(c("red","yellow"))
nLevel<-10
shape4plot<-procGPA(abind(group1,group2))$mshape
for (i in 1:dim(group1)[1]) {
  temp<-kdeValues[i]*nLevel/maxPossible
  myCol<-colfunc(nLevel)[ceiling(temp)]
  plot(shape4plot[i,1],shape4plot[i,2], pch=19, 
       col=myCol,
       xlim = xlim, ylim = ylim,
       xlab = "", ylab = "",cex=2)
  par(new=T)
}
image.plot(legend.only=T, zlim=c(1,nLevel), col=colfunc(nLevel),legend.lab = "Level of significance")
