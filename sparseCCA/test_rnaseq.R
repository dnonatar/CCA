library(PMA)
set.seed(1105)

ccaPerm = CCA.permute(x = microb3, z = rnaseq,
                      typex = "standard", typez = "standard", 
                      nperms = 30, niter = 5, standardize = T, trace = F)
penXtemp = ccaPerm$bestpenaltyx
penZtemp = ccaPerm$bestpenaltyz
ccaRslt = CCA(x = microb3, z = rnaseq,
              typex = "standard", typez = "standard",
              penaltyx = penXtemp, penaltyz = penZtemp,
              K = 2, niter = 5, standardize = T)
sum(ccaRslt$u != 0)
sum(ccaRslt$v != 0)

ccaScoreU = microb3 %*% ccaRslt$u
ccaScoreV = rnaseq %*% ccaRslt$v
ccaScores = cbind(ccaScoreU, ccaScoreV)
colnames(ccaScores) = c("U1", "U2", "V1", "V2")
ccaScores = as.data.frame(ccaScores)
#number of each type should be flexible
ccaScores$type = c(rep("class 1", 17), rep("class 2", 17))

U1_cor <- function(g) {return(cor(ccaScores$U1,g))}
V1_cor <- function(g) {return(cor(ccaScores$V1,g))}
rnaseq_cor <- apply(rnaseq,2,U1_cor)
microb3_cor <- apply(microb3,2,V1_cor)

rnaseq_cor
microb3_cor