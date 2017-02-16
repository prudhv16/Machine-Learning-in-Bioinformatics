library(caret);library(randomForest);library(e1071);library(kernlab);library(doMC);library(foreach);library(RColorBrewer)
require(car);require(Hmisc);require(corrplot);require(neuralnet);require(plyr);require(flux)
setwd("/Users/Jing/Desktop/I529_Final")


tr <- read.csv('trainingData-releaseV2.csv')
dim(tr)

summary(tr[,c(42:56)])
png("test1.png", width=9, height=7, units="in", res=300)
par(mfrow=c(1,2))
boxplot(tr[,c(42:56)],notch=T,col=c(42:56),main='Boxplots of training data')
grid(10, 10, lwd = 2)
hist(tr[,42],col=rgb(1,0,0,0.5),breaks=20,prob=T,xlab='',ylab='',main='')
lines(density(tr[,42]),lwd=3,col=rgb(1,0,0,0.7),lty=5)
hist(tr[,43],col=rgb(0,0.5,0,0.5),add=T,breaks=20,prob=T)
lines(density(tr[,43]),lwd=3,col=rgb(0,0.5,0,0.7),lty=5)
hist(tr[,44],col=rgb(0,0,1,0.5),add=T,breaks=20,prob=T)
lines(density(tr[,44]),lwd=3,col=rgb(0,0,0.7,0.5),lty=5)
title('Histograms for training')
dev.off()

png("Scatter.png", width=9, height=7, units="in", res=300)
scatterplotMatrix(tr[,42:45])
dev.off()

png("heat.png", width=9, height=7, units="in", res=300)
M<-rcorr(as.matrix(tr[,42:56]), type="pearson")$r
corrplot(M, method = "circle")
dev.off()


dim(tr[,-1:-41])
tr1 = data.frame(t(tr))
dim(tr1)
M=as.matrix(dist(tr1[-1:-41,]))
#M=as.matrix(dist(tr[,-1:-41]))
dim(M)

# Determine number of clusters
wss <- (nrow(tr[,-1:-41])-1)*sum(apply(tr[,-1:-41],2,var))
for (i in 2:20) {
  wss[i] <- sum(kmeans(tr[,-1:-41],centers=i,nstart=25, iter.max=1000)$withinss)
}
png("kmeans1.png", width=9, height=7, units="in", res=300)
plot(1:20, wss, type="b",lwd=3, xlab="Number of Clusters",ylab="Within groups sum of squares")
grid(10, 10, lwd = 2)
title('k-means clustering for features')
dev.off()

centers=6
cl<-kmeans(M, centers, iter.max = 100, nstart = 5)
png("kmeans2.png", width=9, height=7, units="in", res=300)
plot(M[,c(2,5)], col = cl$cluster,cex=1.5,pch=cl$cluster)
title('Clustering analysis')
text(M[,c(2,5)]+0.5, labels=colnames(M),col = cl$cluster)
grid(10, 10, lwd = 2)
dev.off()

# Cluster sizes
sort(table(cl$clust))
clust <- names(sort(table(cl$clust)))

table(tr[,14])
# First cluster
row.names(M[cl$clust==clust[1],])
row1<-row.names(M[cl$clust==clust[1],])
# Second Cluster
row.names(M[cl$clust==clust[2],])
row2<-row.names(M[cl$clust==clust[2],])
# Third Cluster
row.names(M[cl$clust==clust[3],])
row3<-row.names(M[cl$clust==clust[3],])
# Fourth Cluster
row.names(M[cl$clust==clust[4],])
row4<-row.names(M[cl$clust==clust[4],])
# Fifth Cluster
row.names(M[cl$clust==clust[5],])
row5<-row.names(M[cl$clust==clust[5],])
row.names(M[cl$clust==clust[6],])

#random Forest
#registerDoMC(cores=6)
#y1 <- as.factor(paste('X.', tr[,14], sep = ''))

rf <- randomForest(tr[,-1:-41], y1, ntree=2500)
rf.mod = grow(rf, 50)
Importance <- importance(rf.mod)
Importance <- as.matrix(Importance)
Importance <- as.matrix(sort(Importance[,1], decreasing=TRUE))
plot(Importance)

png("RandomForest.png", width=9, height=7, units="in", res=300)
op = par(mfrow=c(1,2))
plot(Importance, main="Importance of features")
points(13,0.6, type='h', lwd=2)
grid(10, 10, lwd = 2)
dev.off()

#SVM
set.seed(3)
index <- sample(1:nrow(tr),round(0.75*nrow(tr)))
rownames(Importance)[c(1:8)]
eight<-rownames(Importance)[c(1:8)]
dat.log <- tr[,eight]
colnames(dat.log)

train <- dat.log[index,]
test <- dat.log[-index,]
y.train<-tr[index,14]
y.test<-tr[-index,14]
length(y.test)
length(y.train)
table(y.train)
table(y.test)

datain.train<-cbind(train,y.train)
n <- names(datain.train)
f <- as.formula(paste("y.train ~", paste(n[!n %in% "y.train"], collapse = " + ")))
svm.fit=svm(f,data=tr,kernel ="radial",cost=1)
summary(svm.fit)
tune.out=tune(svm,f,data=datain.train, kernel ="radial"
              ,ranges = list( cost=c(0.01,0.1,1,10,100)))
bestmod=tune.out$best.model
svm.probs = predict(bestmod,test)
table(test[,14],svm.probs)


svm.pred=rep(0,48)
svm.pred[svm.probs=="RESISTANT"] = "RESISTANT"
svm.pred[svm.probs=="CR"] = "CR"
table(svm.pred,y.test)
mean(svm.pred==y.test)





