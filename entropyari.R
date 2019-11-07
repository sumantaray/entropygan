library('ClusterR')
library('EntropyEstimation')
library('MASS')
library(foreach)
library(doParallel)
library('Rfast')
library(markovchain) 
library(infotheo)
library('caret')
library('ggplot2')
library('RMThreshold')
library(GLDEX)
library(fossil)

set.seed(476576)
covar<-matrix(0,50,50)
feanumber=50
for(i in 1:50)
{
  for(j in 1:50)
  {
    covar[i,j]<-0.5^{abs(i-j)}
  }
}

nf=50
ariminrenyi=matrix(0,nrow=100,ncol=1)
arirenyi=matrix(0,nrow=100,ncol=12)
aritsallis=matrix(0,nrow=100,ncol=12)
arishanon=matrix(0,nrow=100,ncol=1)
cluscen=3

jointprob_3 <-function(u) {
  ncol=length(unique(u[,1]))
  nrow=length(unique(u[,2]))*length(unique(u[,3]))
  probtable=matrix(0,nrow,ncol)
  a=unique(u[,1])
  for(i in 1:length(a))
  {
    index=which(u[,1]==a[i])
    pmat=as.matrix(u[index,2:3])
    b=unique(pmat[,1])
    count=1
    for(j in 1: length(b))
    {
      index1=which(pmat[,1]==b[j])
      pmat1=pmat[index1,2]
      s1=aggregate(data.frame(count = pmat1), list(value = pmat1), length)
      lens1=length(s1$count)
      probtable[(count:(count+lens1-1)),i]=(s1$count)/nrow(data)
      count=count+lens1
    }
  }
  probtable
}

jointprob_2 <-function(u) {
  ncol=length(unique(u[,1]))
  nrow=length(unique(u[,2]))
  probtable=matrix(0,nrow,ncol)
  a=unique(u[,1])
  for(i in 1:length(a))
  {
    index=which(u[,1]==a[i])
    pmat1=as.matrix(u[index,2])
    s1=aggregate(data.frame(count = pmat1), list(value = pmat1), length)
    lens1=length(s1$count)
    probtable[1:lens1,i]=(s1$count)/nrow(data)
  }  
  probtable
}
mu_3classnonoverlap <- read.csv("mu_3classnonoverlap.csv", header=FALSE)

ptm <- proc.time()
for(x in 1:100)
{
mu=as.matrix(mu_3classnonoverlap)
gen.mix <- function(n, k, mu, sig) {
  library(MASS)
  
  d <- length(mu[1,])  # number of dimensions
  result <- matrix(rep(NA,n*d), ncol=d)
  colnames(result) <- paste0("X",1:d)
  
  for(i in 1:n) {
    result[i,] <- mvrnorm(1, mu = mu[k[i],], Sigma=sig[,,k[i]])
  }
  
  result
}
sigs <- array(rep(NA,feanumber*feanumber*3), c(feanumber,feanumber,3))  # 3D matrix
sigs[,,1] <- covar
sigs[,,2]<-covar
sigs[,,3]<-covar
#sigs[,,4]<-covar
#sigs[,,5]<-covar
pi<-c(0.6,0.2,0.2)
classesnormal <- sample(1:3, 500, replace=TRUE, prob=pi)
mydata <- gen.mix(500, classesnormal, mu, sigs)
gmm = GMM(mydata, 3, dist_mode = "eucl_dist", seed_mode = "random_subset", km_iter = 10,em_iter = 10, verbose = F)          
siggmm <- array(rep(NA,feanumber*feanumber*3), c(feanumber,feanumber,3))  # 3D matrix
siggmm[,,1] <- diag(gmm$covariance_matrices[1,])
siggmm[,,2]<-diag(gmm$covariance_matrices[2,])
siggmm[,,3]<-diag(gmm$covariance_matrices[3,])
#siggmm[,,4]<-diag(gmm$covariance_matrices[4,])
#siggmm[,,5]<-diag(gmm$covariance_matrices[5,])
classes <- sample(1:3, 500, replace=TRUE, prob=(gmm$weights))
gmmdata <- gen.mix(500, classes, gmm$centroids, siggmm)
data<-matrix(0,nrow=500,ncol=300)
data[,1:feanumber]<-gmmdata
for(i in (feanumber+1):300)
{rn<-matrix(0,nrow=500,ncol=1)
rn[,1]<-runif(500)
data[,i]<-add.Gaussian.noise(as.matrix(rn), mean = 0, stddev = 1, symm = FALSE)
}

mun=as.matrix(3*mu_3classnonoverlap)
pi<-c(0.6,0.2,0.2)
classesn <- sample(1:3, 500, replace=TRUE, prob=pi)
newdata<- gen.mix(500, classesn, mun, sigs)
class1=as.matrix(classes)
indexorg=which(class1==1)
indexnew=which(classesn==2)
data[indexorg[1:50],1:50]=newdata[indexorg[1:50],1:50]
indexorg=which(class1==2)
indexnew=which(classesn==3)
data[indexorg[1:25],1:50]=newdata[indexorg[1:25],1:50]
indexorg=which(class1==3)
indexnew=which(classesn==1)
data[indexorg[1:25],1:50]=newdata[indexorg[1:25],1:50]




data<-discretize(data)
classc<-as.matrix(classes)
n <- nrow(data)
col<-ncol(data)
count=ncol(data)





# Renyi min entropy
fea<- matrix(0, nrow=1,ncol =nf)
cl <- makeCluster(7)
registerDoParallel(cl)
parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
{
  u<-as.matrix(cbind(classc,data[,j]))
  tmat=jointprob_2(u)
  transvecfea=-(log(sum(rowMaxs(tmat,value=TRUE))))
}

mimc1<-as.matrix(parl)
fea[1,1]<-which(mimc1==min(mimc1[which(mimc1>0)]))[1]
idx=which(mimc1==min(mimc1[which(mimc1>0)]))[1]

for(i in 1:(nf-1))
{  
  parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
  {
    u2<-as.matrix(cbind(classc,data[,idx],data[,j]))
    tmat=jointprob_3(u2)
    feat=-(log(sum(rowMaxs(tmat,value=TRUE))))
  }
  feat<-as.matrix(parl)
  feat[fea[1:i]]<-Inf
  idx=which(feat==min(feat[which(feat>0)]))[1]
  fea[1,(i+1)]=idx
}
stopCluster(cl)
clus<-kmeans(data[,fea],centers=cluscen)
ariminrenyi[x,1]<-adj.rand.index(classc, clus$cluster)



# Renyi entropy with q
cl <- makeCluster(7)
registerDoParallel(cl)
for(m in 1:12)
{qindex=c(0.1,0.3,0.5,0.7,1.3,1.5,1.7,2,2.5,3,5,7)
q=qindex[m]
fea<- matrix(0, nrow=1,ncol =nf)
parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
{
  u<-as.matrix(cbind(classc,data[,j]))
  tmat=jointprob_2(u)
  tmat=tmat^q
  transvecfea=(q/(1-q))*(log(sum((rowsums(tmat))^(1/q))))
}

mimc1<-as.matrix(parl)
fea[1,1]<-which(mimc1==min(mimc1[which(mimc1>0)]))[1]
idx=which(mimc1==min(mimc1[which(mimc1>0)]))[1]

for(i in 1:(nf-1))
{  
  parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
  {
    u2<-as.matrix(cbind(classc,data[,idx],data[,j]))
    tmat=jointprob_3(u2)
    tmat=tmat^q
    transvecfea=(q/(1-q))*(log(sum((rowsums(tmat))^(1/q))))
  }
  feat<-as.matrix(parl)
  feat[fea[1:i]]<-Inf
  idx=which(feat==min(feat[which(feat>0)]))[1]
  fea[1,(i+1)]=idx
}
clus<-kmeans(data[,fea],centers=cluscen)
arirenyi[x,m]<-adj.rand.index(classc, clus$cluster)
}
stopCluster(cl)

# Tsallis Entropy
cl <- makeCluster(7)
registerDoParallel(cl)
for(m in 1:12)
{qindex=c(0.1,0.3,0.5,0.7,1.3,1.5,1.7,2,2.5,3,5,7)
q=qindex[m]
fea<- matrix(0, nrow=1,ncol =nf)
parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
{
  u<-as.matrix(cbind(classc,data[,j]))
  tmat=jointprob_2(u)
  tmatq=tmat^q
  nem=sum(tmatq)
  s1=aggregate(data.frame(count = data[,j]), list(value = data[,j]), length)
  lens1=length(s1$count)
  probnew=(s1$count)/nrow(data)
  denom=sum(probnew^q)
  trans=(1/(q-1))*(1-(nem/denom))
}

mimc1<-as.matrix(parl)
fea[1,1]<-which(mimc1==min(mimc1[which(mimc1>0)]))[1]
idx=which(mimc1==min(mimc1[which(mimc1>0)]))[1]

for(i in 1:(nf-1))
{  
  parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast")) %dopar%
  {
    u2<-as.matrix(cbind(classc,data[,idx],data[,j]))
    tmat=jointprob_3(u2)
    tmatq=tmat^q
    nem=sum(tmatq)
    u3<-as.matrix(cbind(data[,idx],data[,j]))
    probnew=jointprob_2(u3)
    denom=sum(probnew^q)
    trans=(1/(q-1))*(1-(nem/denom))
  }
  feat<-as.matrix(parl)
  feat[fea[1:i]]<-Inf
  #idx=which.min(feat)
  idx=which(feat==min(feat[which(feat>0)]))[1]
  fea[1,(i+1)]=idx
}
clus<-kmeans(data[,fea],centers=cluscen)
aritsallis[x,m]<-adj.rand.index(classc, clus$cluster)
}
stopCluster(cl)


# Shanon Entropy
cl <- makeCluster(7)
registerDoParallel(cl)
fea<- matrix(0, nrow=1,ncol =nf)
parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast","GLDEX")) %dopar%
{
  u<-as.matrix(cbind(classc,data[,j]))
  tmat=jointprob_2(u)
  tmat=as.vector(tmat)
  tmat=fun.zero.omit(tmat)
  ent1=-sum(tmat*log2(tmat))
  s1=aggregate(data.frame(count = data[,j]), list(value = data[,j]), length)
  lens1=length(s1$count)
  probnew=(s1$count)/nrow(data)
  ent2=-sum(probnew*log2(probnew))
  trans=ent1-ent2
}

mimc1<-as.matrix(parl)
fea[1,1]<-which(mimc1==min(mimc1[which(mimc1>0)]))[1]
idx=which(mimc1==min(mimc1[which(mimc1>0)]))[1]

for(i in 1:(nf-1))
{  
  parl<-foreach(j=1:count, .combine=c,.packages= c("Rfast","GLDEX")) %dopar%
  {
    u2<-as.matrix(cbind(classc,data[,idx],data[,j]))
    tmat=jointprob_3(u2)
    tmat=as.vector(tmat)
    tmat=fun.zero.omit(tmat)
    ent1=-sum(tmat*log2(tmat))
    u3<-as.matrix(cbind(data[,idx],data[,j]))
    probnew=jointprob_2(u3)
    probnew=as.vector(probnew)
    probnew=fun.zero.omit(probnew)
    ent2=-sum(probnew*log2(probnew))
    trans=ent1-ent2
  }
  feat<-as.matrix(parl)
  feat[fea[1:i]]<-Inf
  #idx=which.min(feat)
  idx=which(feat==min(feat[which(feat>0)]))[1]
  fea[1,(i+1)]=idx
}
clus<-kmeans(data[,fea],centers=cluscen)
arishanon[x,1]<-adj.rand.index(classc, clus$cluster)
stopCluster(cl)
print(x)
}
registerDoSEQ()
proc.time() - ptm
write.table(arirenyi,file="/home/snehalika/Desktop/entropyfeature/arirenyi.csv",sep=",",row.names = FALSE,col.names = FALSE)
write.table(ariminrenyi,file="/home/snehalika/Desktop/entropyfeature/ariminrenyi.csv",sep=",",row.names = FALSE,col.names = FALSE)
write.table(aritsallis,file="/home/snehalika/Desktop/entropyfeature/aritsallis.csv",sep=",",row.names = FALSE,col.names = FALSE)
write.table(arishanon,file="/home/snehalika/Desktop/entropyfeature/arishanon.csv",sep=",",row.names = FALSE,col.names = FALSE)



