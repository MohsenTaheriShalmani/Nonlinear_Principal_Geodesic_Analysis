library(shapes)
library(rgl)

#####################################################################################################
#####################################################################################################
#functions

# Elimination
elimination <- function(G1,G2) {
  if(dim(G1)[1]!=dim(G2)[1]){
    cat("Error: Number of landmarks in two groups are not equal!!!/n")
    break
  }
  
  k<-dim(G1)[1]
  x<-1:dim(G1)[1]
  Index<-1
  sortedLandmarks<-c()
  distances<-c()
  while(length(x)>2){
    minProcDistance<-Inf
    proc1<-procGPA(G1[x,,], scale = F)
    proc2<-procGPA(G2[x,,], scale = F)
    temp2<-procdist(proc1$mshape,proc2$mshape)
    distances<-c(distances,temp2)
    for (i in x) {
      temp<-x[!x %in% i]
      print(temp)
      proc1<-procGPA(G1[temp,,], scale = F)
      proc2<-procGPA(G2[temp,,], scale = F)
      temp3<-procdist(proc1$mshape,proc2$mshape)
      if(temp3<minProcDistance){
        minProcDistance<-temp3
        Index<-i
      }
    }
    cat("Min proc distance is:",minProcDistance,"\n")
    sortedLandmarks<-c(sortedLandmarks,Index)
    x<-x[!x %in% sortedLandmarks]
  }
  
  sortedLandmarks<-c(sortedLandmarks,x)
  distances<-c(distances,c(0,0))
  pointsAndDis<-cbind(sortedLandmarks,distances)  #landmarks and distances
  
  return(pointsAndDis)
  
}

#####################################################################################################
#####################################################################################################
# Example 1. 7-landmarks female gorilla skull vs. modified female gorilla skull 

group1<-gorf.dat[c(1:2,4:8),,]
changeMatrix<-t(matrix(c(50,20,0,0,0,0,0,0,10,0,0,0,0,0),nrow = 2))
group2<-array(NA, dim = dim(group1))
for (i in 1:dim(group1)[3]) {
  group2[,,i]<-group1[,,i]+changeMatrix
}

#elimination
pointsAndDis<-elimination(G1 = group1,G2 = group2)
pointsAndDis

# threshold<-mean(pointsAndDis[,2])
threshold<-0.03
suspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]>=threshold),1]
suspiciousPoints
unsuspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]<threshold),1]
unsuspiciousPoints

#plot suspiciousPoints (red color indicate significat points)
shape4plot<-procGPA(abind(group1,group2))$mshape
par(new=F)
plot(shape4plot, pch=3, col="blue" , main = "gorf vs. Modified gorf",
     xlim = c(-170,170), ylim = c(-170,170),
     xlab = "x", ylab = "y")
par(new=T)
plot(shape4plot[suspiciousPoints,], pch=16, col="red" ,
     xlim = c(-170,170), ylim = c(-170,170),
     xlab = "x", ylab = "y")
legend("bottomright", legend=c("unsuspicious", "suspicious"),pch = c(3,16),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)

k<-dim(group1)[1]
weights<-c(1:k)
weights[unsuspiciousPoints]<-0.0001
weights[-unsuspiciousPoints]<-1

covW<-weights*diag(k)

pooledGroup<-abind(group1,group2)
procPooledW<-procWGPA(pooledGroup,fixcovmatrix =covW ,scale = T)

groupW1<-procPooledW$rotated[,,1:dim(group1)[3]]
groupW2<-procPooledW$rotated[,,(dim(group1)[3]+1):(dim(group1)[3]+dim(group2)[3])]

par(new=F)
for (i in 1:dim(groupW1)[3]) {
  plot(groupW1[,,i], pch=3, col="blue" ,
       xlim = c(-170,170), ylim = c(-170,180),
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(groupW2)[3]) {
  plot(groupW2[,,i], pch=4, col="red" ,
       xlim = c(-170,170), ylim = c(-170,180),
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(3,4),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)

#####################################################################################################
#####################################################################################################
# Example 2. female gorilla skull vs male gorilla skull

group1<-gorf.dat
group2<-gorm.dat

#elimination
pointsAndDis<-elimination(G1 = group1,G2 = group2)
pointsAndDis

threshold<-mean(pointsAndDis[,2])
threshold<-0.03
suspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]>=threshold),1]
suspiciousPoints
unsuspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]<threshold),1]
unsuspiciousPoints

#plot suspiciousPoints (red color indicate significat points)
shape4plot<-procGPA(abind(group1,group2))$mshape
par(new=F)
plot(shape4plot, pch=3, col="blue" , main = "gorf vs. gorm",
     xlim = c(-170,170), ylim = c(-170,170),
     xlab = "x", ylab = "y")
par(new=T)
plot(shape4plot[suspiciousPoints,], pch=16, col="red" ,
     xlim = c(-170,170), ylim = c(-170,170),
     xlab = "x", ylab = "y")
legend("bottomright", legend=c("unsuspicious", "suspicious"),pch = c(3,16),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)

k<-dim(group1)[1]
weights<-c(1:k)
weights[unsuspiciousPoints]<-0.0001
weights[-unsuspiciousPoints]<-1

covW<-weights*diag(k)

pooledGroup<-abind(group1,group2)
procPooledW<-procWGPA(pooledGroup,fixcovmatrix =covW ,scale = T)

groupW1<-procPooledW$rotated[,,1:dim(group1)[3]]
groupW2<-procPooledW$rotated[,,(dim(group1)[3]+1):(dim(group1)[3]+dim(group2)[3])]

par(new=F)
for (i in 1:dim(groupW1)[3]) {
  plot(groupW1[,,i], pch=3, col="blue" ,
       xlim = c(-170,170), ylim = c(-190,170),
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(groupW2)[3]) {
  plot(groupW2[,,i], pch=4, col="red" ,
       xlim = c(-170,170), ylim = c(-190,170),
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(3,4),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)


#####################################################################################################
#####################################################################################################
# Example 3. female chimpanzee skull vs male chimpanzee skull

group1<-panf.dat
group2<-panm.dat

#elimination
pointsAndDis<-elimination(G1 = group1,G2 = group2)
pointsAndDis

threshold<-mean(pointsAndDis[,2])
threshold<-0.03
suspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]>=threshold),1]
suspiciousPoints
unsuspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]<threshold),1]
unsuspiciousPoints

#plot suspiciousPoints (red color indicate significat points)
shape4plot<-procGPA(abind(group1,group2))$mshape
par(new=F)
plot(shape4plot, pch=3, col="blue" , main = "panf vs. panm",
     xlim = c(-170,170), ylim = c(-170,170),
     xlab = "x", ylab = "y")
par(new=T)
plot(shape4plot[suspiciousPoints,], pch=16, col="red" ,
     xlim = c(-170,170), ylim = c(-170,170),
     xlab = "x", ylab = "y")
legend("bottomright", legend=c("unsuspicious", "suspicious"),pch = c(3,16),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)

k<-dim(group1)[1]
weights<-c(1:k)
weights[unsuspiciousPoints]<-0.0001
weights[-unsuspiciousPoints]<-1

covW<-weights*diag(k)

pooledGroup<-abind(group1,group2)
procPooledW<-procWGPA(pooledGroup,fixcovmatrix =covW ,scale = T)

groupW1<-procPooledW$rotated[,,1:dim(group1)[3]]
groupW2<-procPooledW$rotated[,,(dim(group1)[3]+1):(dim(group1)[3]+dim(group2)[3])]

par(new=F)
for (i in 1:dim(groupW1)[3]) {
  plot(groupW1[,,i], pch=3, col="blue" ,
       xlim = c(-170,170), ylim = c(-190,170),
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(groupW2)[3]) {
  plot(groupW2[,,i], pch=4, col="red" ,
       xlim = c(-170,170), ylim = c(-190,170),
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(3,4),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)



#####################################################################################################
#####################################################################################################
# Example 4. schizophrenia vs control
data(schizophrenia)
schizophrenia$group
group1<-schizophrenia$x[,,1:14]
group2<-schizophrenia$x[,,15:28]

#elimination
pointsAndDis<-elimination(G1 = group1,G2 = group2)
pointsAndDis

threshold<-mean(pointsAndDis[,2])
threshold<-0.03
suspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]>=threshold),1]
suspiciousPoints
unsuspiciousPoints<-pointsAndDis[which(pointsAndDis[,2]<threshold),1]
unsuspiciousPoints

#plot suspiciousPoints (red color indicate significat points)
shape4plot<-procGPA(abind(group1,group2))$mshape
par(new=F)
plot(shape4plot, pch=3, col="blue" ,main = "Schizophrenia vs. Control",
     xlim = c(-1,1), ylim = c(-1,1),
     xlab = "x", ylab = "y")
par(new=T)
plot(shape4plot[suspiciousPoints,], pch=16, col="red" ,
     xlim = c(-1,1), ylim = c(-1,1),
     xlab = "x", ylab = "y")
legend("bottomright", legend=c("unsuspicious", "suspicious"),pch = c(3,16),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)

k<-dim(group1)[1]
weights<-c(1:k)
weights[unsuspiciousPoints]<-0.0001
weights[-unsuspiciousPoints]<-1

covW<-weights*diag(k)

pooledGroup<-abind(group1,group2)
procPooledW<-procWGPA(pooledGroup,fixcovmatrix =covW ,scale = T)

groupW1<-procPooledW$rotated[,,1:dim(group1)[3]]
groupW2<-procPooledW$rotated[,,(dim(group1)[3]+1):(dim(group1)[3]+dim(group2)[3])]

par(new=F)
for (i in 1:dim(groupW1)[3]) {
  plot(groupW1[,,i], pch=3, col="blue" ,
       xlim = c(-1,1), ylim = c(-1,1),
       xlab = "x", ylab = "y")
  par(new=T)
}
for (i in 1:dim(groupW2)[3]) {
  plot(groupW2[,,i], pch=4, col="red" ,
       xlim = c(-1,1), ylim = c(-1,1),
       xlab = "x", ylab = "y")
  par(new=T)
}
legend("bottomright", legend=c("Group A", "Group B"),pch = c(3,4),
       col=c("blue","red"), cex=0.6,pt.cex = 1, horiz=T)

