library(PMA)
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
#number of each type should be flexible
ccaScores$type = c(rep("class 1", 17), rep("class 2", 17))

