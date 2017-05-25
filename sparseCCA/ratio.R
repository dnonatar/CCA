## test for the ratio at which the relationship is undetectable
# gaussian

microb = read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Arghavan/MetaSample2/samples.txt")
microb = as.matrix(microb)
microb <- microb[,-1]
microb_rows = dim(microb)[1]/2   # divide by two to choose only the first class
microb = microb[1:microb_rows,] 
sd <- apply(microb,2,sd)
microb <- microb[,which(sd!=0)]  # choose only OTUs with non-zero standard deviation
microb <- microb[,c(417,541,42)] # two OTUs with min numbers of zeros

gaussian <- read.table("/home/ratanond/Desktop/Masters_Project/Synthetic/Eunji/rnaseq_cls/rnaseq_cls/out/myparams_ma_trn.txt",sep = ",")
gaussian <- as.matrix(gaussian)
gaussian_rows <- dim(gaussian)[1]/2  
gaussian <- gaussian[1:gaussian_rows,1:10] # try with 10 genes

y <- microb[,1:3]
x <- gaussian
x1 <- gaussian[,1]
x2 <- gaussian[,2]
x3 <- gaussian[,3]
x4 <- gaussian[,4]
x5 <- gaussian[,5]
x6 <- gaussian[,6]
x7 <- gaussian[,7]
x8 <- gaussian[,8]
x9 <- gaussian[,9]
x10 <- gaussian[,10]

data = data.frame(y,x1,x2,x3,x4,x5,x6,x7,x8,x9,x10)

microb1<-microb2<-microb3<-microb4<-microb5<-microb6<-microb7<-microb8<-microb9<-microb10<-microb
microb1[,1] <-predict(lm(y~x1,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb2[,1] <-predict(lm(y~x1+x2,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb3[,1] <-predict(lm(y~x1+x2+x3,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb4[,1] <-predict(lm(y~x1+x2+x3+x4,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb5[,1:2] <-predict(lm(y~x1+x2+x3+x4+x5,data = data, subset = y[,1]>0&y[,2]>0),newdata = as.data.frame(gaussian))
microb6[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb7[,1:2] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7,data = data, subset = y[,1]>0&y[,2]>0),newdata = as.data.frame(gaussian))
microb8[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7+x8,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb9[,1] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data = data, subset = y>0),newdata = as.data.frame(gaussian))
microb10[,1:2] <-predict(lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9+x10,data = data, subset = y[,1]>0&y[,2]>0),newdata = as.data.frame(gaussian))


for(k in 1:k) {
  for (n in 1:n){
    m <- lm(y[])
    
    
  }
}


library(PMA)

myCCA <- function(microb_new,gaussian) {
set.seed(1105)
ccaPerm = CCA.permute(x = microb_new, z = gaussian,
                      typex = "standard", typez = "standard", 
                      nperms = 30, niter = 5, standardize = T)
penXtemp = ccaPerm$bestpenaltyx
penZtemp = ccaPerm$bestpenaltyz
ccaRslt = CCA(x = microb_new, z = gaussian,
              typex = "standard", typez = "standard",
              penaltyx = penXtemp, penaltyz = penZtemp,
              K = 2, niter = 5, standardize = T)
sum(ccaRslt$u != 0)
sum(ccaRslt$v != 0)

ccaScoreU = microb_new %*% ccaRslt$u
ccaScoreV = gaussian %*% ccaRslt$v
ccaScores = cbind(ccaScoreU, ccaScoreV)
colnames(ccaScores) = c("U1", "U2", "V1", "V2")
ccaScores = as.data.frame(ccaScores)
return(cor(ccaScores$U1,ccaScores$V1))
}



