## test for the ratio at which the relationship is undetectable
# gaussian and rnaseq

microb = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Arghavan/MetaSample2/samples.txt")
microb = as.matrix(microb)
microb <- microb[,-1]
microb_rows = dim(microb)[1]/2  # divide by two to choose only the first class
microb = microb[1:microb_rows,] 
sd <- apply(microb,2,sd)
microb <- microb[,which(sd!=0)]  # choose only OTUs with non-zero standard deviation
## The top three are 417,541,42
microb <- microb[,c(417,541)] # two OTUs with min numbers of zeros

gaussian = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/params_sb_high_1_ma_trn.txt",sep = ",")
gaussian <- as.matrix(gaussian)
gaussian_rows <- dim(gaussian)[1]/2    # divide by two to choose only the first class
gaussian <- gaussian[1:gaussian_rows,c(4,7,10)]

rnaseq = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/params_sb_high_1_trn.txt",sep = ",")
rnaseq = as.matrix(rnaseq)
rnaseq_rows = dim(rnaseq)[1]/2   # divide by two to choose only the first class
rnaseq = rnaseq[1:rnaseq_rows,c(4,7,10)]  

y <- microb[,1]
x1 <- gaussian[,1]
x2 <- gaussian[,2]
x3 <- gaussian[,3]
data = data.frame(y,x1,x2,x3)
#x4 <- gaussian[,4]
#x5 <- gaussian[,5]
#x6 <- gaussian[,6]
#x7 <- gaussian[,7]
#x8 <- gaussian[,8]
#x9 <- gaussian[,9]
#x10 <- gaussian[,10]

#data = data.frame(y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)

# For small cases
microb1<-microb2<-microb3<-microb4<-microb5<-microb6<-microb7<-microb8<-microb9<-microb10<-microb
microb1[,1] <-predict(lm(y~x1,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb2[,1] <-predict(lm(y~x1+x2,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb3[,1] <-predict(lm(y~x1+x2+x3,data = data, subset = y>0),newdata = as.data.frame(gaussian))

microb4[,1] <-predict(lm(y~x1+x2+x3+x4,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb5[,1] <-predict(lm(y~x1+x2+x3+x4+x5,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb6[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb7[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb8[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7+x8,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb9[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb10[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data = data, subset = y>0),newdata = as.data.frame(gaussian))


####  for high number of n, use for-loop instead
cor_g<-numeric(dim(gaussian)[2])
cor_rna<-numeric(dim(rnaseq)[2])
  for (n in 1:100){
    new_microb <- microb
    m <- lm(y~gaussian[,1:n],subset = y>0)
    new_microb[,1] <- predict(m, newdata=as.data.frame(gaussian))
    cor_g[n] <- myCCA(new_microb, gaussian)
    cor_rna[n] <- myCCA(new_microb,rnaseq)
  }


myCCA <- function(microb_new,gaussian) {
library(PMA)
set.seed(1105)
ccaPerm = CCA.permute(x = microb_new, z = gaussian,
                      typex = "standard", typez = "standard", 
                      nperms = 30, niter = 5, standardize = T,trace = F)
penXtemp = ccaPerm$bestpenaltyx
penZtemp = ccaPerm$bestpenaltyz
ccaRslt = CCA(x = microb_new, z = gaussian,
              typex = "standard", typez = "standard",
              penaltyx = penXtemp, penaltyz = penZtemp,
              K = 2, niter = 5, standardize = T)
#sum(ccaRslt$u != 0)
#sum(ccaRslt$v != 0)

ccaScoreU = microb_new %*% ccaRslt$u
ccaScoreV = gaussian %*% ccaRslt$v
ccaScores = cbind(ccaScoreU, ccaScoreV)
colnames(ccaScores) = c("U1", "U2", "V1", "V2")
ccaScores = as.data.frame(ccaScores)
return(cor(ccaScores$U1,ccaScores$V1))
}





