if(!file.exists("/home/ratanond/prototype.db")){
  stop("Database does not exist!")
}

db = dbConnect(SQLite(), dbname = "/home/ratanond/prototype.db")
writeLines("Starting Sparse Canonical Correlation Analysis")

while(TRUE) {
res <- dbGetQuery(db, "SELECT * from ccajobsss where (status=='Incomplete')")

if(nrow(res) > 0){

row <- dbGetQuery(db, "SELECT name, project, rand, id from ccajobsss where (status=='Incomplete') limit 1")
    #args = c(row[1,1],row[1,2])
    randname = paste(row[1,3],row[1,4], sep="")
    writeLines(paste("Processing job",randname,sep=": "))


result <- tryCatch({

#######################################
######## For Microbiota OTU ###########
#######################################

microBiomeSampFileNames = read.table("/home/ratanond/Desktop/Thesis/Prototype/app/microbiota/sample_filenames.csv", 
                                     stringsAsFactors=FALSE, sep=',',
                                     comment.char="%", quote='"',
                                     header=TRUE) 

# create a list of length 12. Each slot corresponds to the data in each file in metabolic_profiles folder.
microBiomeSampList = NULL
for(f in 1:length(microBiomeSampFileNames$filename)){
  microBiomeSampList[[f]] = 
    read.table(paste("/home/ratanond/Desktop/Thesis/Prototype/app/microbiota/",
                     microBiomeSampFileNames$filename[f], sep=""), 
               stringsAsFactors=TRUE, sep='\t',
               comment.char="%", quote='"',
               header=TRUE) 
}


microBiomeStack = data.frame("id" = microBiomeSampFileNames$sample[1],
                             microBiomeSampList[[1]])
for(i in 2:length(microBiomeSampList)){
  microBiomeStack = rbind(microBiomeStack, 
                          data.frame("id" = microBiomeSampFileNames$sample[i],
                                     microBiomeSampList[[i]]))
}  

#### In the end, we get a dataframe containing the data we will be using.####

# give new column names
names(microBiomeStack) = c("id", "subHier1", "subHier2", 
                           "subHier3", "hits")
sampNames = unique(microBiomeStack$id)

seedLev1Names = unique(microBiomeStack$subHier1)

#cluster hits by id and subHier1 (?)
seedLev1Counts = aggregate(hits ~ id + subHier1,
                           data = microBiomeStack,
                           FUN = sum)

# need packege reshape2 in order to use melt and dcast
library(reshape2)
seedLev1Counts = dcast(seedLev1Counts, subHier1 ~ id,
                       fill = 0, value.var = "hits") ## value.var needs a character
rownames(seedLev1Counts) = seedLev1Counts[, "subHier1"]
seedLev1Counts = seedLev1Counts[, -1]
# in the end we get a dataframe containing the counts of each gene for each subHeir1

## Or do directly:
# seedLev1Counts = dcast(formula = subHier1~id, value.var = "hits", 
#                        data = microBiomeStack, fun = sum)

## For system level2:

#why paste0 step ??
microBiomeStack$subHier2 = paste0(microBiomeStack$subHier1, "_",
                                  microBiomeStack$subHier2)
seedLev2Counts = dcast(data = microBiomeStack, subHier2 ~ id,
                       value.var = "hits", fun = sum)
rownames(seedLev2Counts) = seedLev2Counts[,"subHier2"]
seedLev2Counts = seedLev2Counts[,-1]
# in the end we get a dataframe containing the counts of each gene for each subHeir1_subHeir2

## For system level3:
seedLev3Counts = dcast(microBiomeStack, subHier3 ~ id,
                       value.var = "hits", fun = sum)
rownames(seedLev3Counts) = seedLev3Counts[,"subHier3"]
seedLev3Counts = seedLev3Counts[,-1]
# in the end we get a dataframe containing the counts of each gene for each subHeir3

## Do count OTUs normalization follow the way below
## https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html
library(metagenomeSeq)
set.seed(1105)
seedLev2Mr = newMRexperiment(seedLev2Counts)
seedLev2P = cumNormStat(seedLev2Mr)
seedLev2Mr = cumNorm(seedLev2Mr, p = seedLev2P)
seedLev2Norm = MRcounts(seedLev2Mr, norm = T, log = T)

test = seedLev2Counts
test[test == 0] = NA
test = na.omit(test)
seedLev2TestMr = newMRexperiment(test)
seedLev2TestP = cumNormStat(seedLev2TestMr)
seedLev2TestMr = cumNorm(seedLev2TestMr, p = seedLev2TestP)
seedLev2TestNorm = MRcounts(seedLev2TestMr, norm = T, log = T)

library(ggplot2)
temp = melt(seedLev2Counts)
temp$variable = factor(temp$variable,
                       levels = sort(as.character(temp$variable)))
seedLev2Box = ggplot(data = temp, aes(x = variable, y = value)) +
  geom_boxplot(aes(fill = variable)) +
  scale_y_log10(breaks = c(1,10,100,1000), name = "Raw Counts") + 
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        plot.margin = unit(c(1,1,1,1), "cm"))
seedLev2Box

temp = melt(seedLev2TestNorm)
temp$variable = factor(temp$Var2,
                       levels = sort(as.character(temp$Var2)))
seedLev2TestBox = ggplot(data = temp,
                         aes(x = Var2, y = value)) +
  geom_boxplot(aes(fill = Var2)) + ylab("Normalized") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 13),
        plot.margin = unit(c(1,1,1,1), "cm"))
seedLev2TestBox


###############################
####### For Microarray ########
###############################

iddoSampPheno = read.table("/home/ratanond/Desktop/Thesis/Prototype/app/microarray/file_names.csv", stringsAsFactors=FALSE,
                           sep=',', comment.char="#", header=TRUE)
iddoSampPheno[1,3] = 4 ## Fix a typo

iddoSampPheno$name = paste0(iddoSampPheno$group, iddoSampPheno$id)

## Import data: making the data set accessible
iddoSampRawVals = NULL
iddoSampBkgd = NULL
iddoSampFlags = NULL

for(i in iddoSampPheno$filename){
  iddoSample =
    read.table(paste("/home/ratanond/Desktop/Thesis/Prototype/app/microarray/", i, sep=''),
               stringsAsFactors=FALSE, sep='\t', comment.char="#", quote='"',
               skip=7, nrows=54359, header=TRUE, fill=TRUE)
  ## One sample
  dim(iddoSample)
  
  ## Only 2:4 should vary between samples, 5 is close, the rest are fixed identifiers
  iddoSample = iddoSample[,c(1:4,12,16,21,27,32,24,25,33,34,35)]
  names(iddoSample)=c("ID_REF", "Raw_intensity", "VALUE", "Quality_flag",
                      "Bkgd_mean", "Probe_name",  "Annotation_NCBI_Acc",
                      "Probe_type","Annotation_Pub_Probe_Targets",
                      "Annotation_OGS", "Annotation_UniGene",
                      "Annotation_Molecular_Function",
                      "Annotation_Biological_Process",
                      "Annotation_Cellular_Component")
  
  iddoSampRawVals = cbind(iddoSampRawVals, iddoSample$Raw_intensity)
  iddoSampBkgd = cbind(iddoSampBkgd, iddoSample$Bkgd_mean)
  iddoSampFlags = cbind(iddoSampFlags, iddoSample$Quality_flag)
}
iddoAnno = iddoSample$Annotation_OGS


sampAnno = cbind(idRef = iddoSample$ID_REF,
                 probName = iddoSample$Probe_name, iddoAnno)
rownames(sampAnno) = sampAnno[,1]

## The damir list does not include the paired data
pairedSamples = which(iddoSampPheno$name %in% colnames(seedLev2Norm))
pairedPheno = iddoSampPheno[pairedSamples,]


rownames(pairedPheno) = pairedPheno[,"name"]
# change BMS group to BF
pairedPheno[which(pairedPheno[,"group"] == "BMS"),"group"] ="BF"
pairedExpr = iddoSampRawVals[,pairedSamples]
colnames(pairedExpr) = rownames(pairedPheno)
rownames(pairedExpr) = iddoSample$ID_REF

####in the end we get the matrix showing intensity of each gene (values from iddoSample$Raw_intensity) #####

## Remove the rows which are less than 0,
# pairedExpr[which(pairedExpr==-9999)] = NA ## Automatically vectorized
pairedExpr[which(pairedExpr <= 0)] = NA
remaFeat = which(rowSums(is.na(pairedExpr)) == 0)
pairedExpr = pairedExpr[remaFeat, ]
sampAnno = sampAnno[remaFeat,]

## Create ExpressionSet object:
pairedPheno = AnnotatedDataFrame(as.data.frame(pairedPheno))
pairedFeat = AnnotatedDataFrame(as.data.frame(sampAnno))
# pairedExprBkCor = oligo::backgroundCorrect(pairedExpr, method = "rma")
eset = new("ExpressionSet", exprs = log(pairedExpr, base = 2),
           phenoData = pairedPheno,
           featureData = pairedFeat)


library(oligo)

oligo::MAplot(eset, group = as.factor(pData(eset)$group),
              trans = identity, pairs = T) ## "main" doesn't work here

oligo::hist(eset, trans = identity,
            xlab = expression(paste(log[2], " of Raw Intensity")),
            ylab = "",
            main = "Density Estimation of Raw Intensity for Each Sample",
            cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.8)

# oligo::boxplot(eset, trans = identity, cex.axis = 1, las = 2,
#                main = "")
# title(ylab = expression(paste(log[2], " of Raw Intensity")), line = 2)
## Or use # oldpar = par(); par(mgp = c(2,1,0))

# exprs returns a matrix of expression values
temp = melt(exprs(eset))
temp$Var2 = factor(temp$Var2, sort(as.character(unique(temp$Var2))))

geneExpRawBox = ggplot(data = temp,
                       aes(x = Var2, y = value)) +
  stat_boxplot(geom = 'errorbar', coef = 5, width = 0.3) +
  geom_boxplot(aes(fill = Var2), outlier.size = NA, coef = 0) +  
  ylab(expression(paste(log[2], " of Raw Intensity"))) + 
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))
print(geneExpRawBox)

### Normalization ###
## To do bkgd correction and normalization:
library(affyPLM)
set.seed(1105) # in case, seed is needed
esetNorm = normalize.ExpressionSet.quantiles(eset = eset)

pairedExprSumm = as.data.frame(exprs(esetNorm))
pairedExprSumm = cbind(geneNames = pairedFeat@data[,3], pairedExprSumm)
pairedExprSumm = pairedExprSumm[base::which(pairedExprSumm[,1] != "NULL"),]
# cluster based on gene names (dot means include all variables)
pairedExprSumm = aggregate(.~geneNames, data = pairedExprSumm,
                           FUN = "median") # It is robust to "quantiles" method
rownames(pairedExprSumm) = pairedExprSumm[,1]
pairedExprSumm = pairedExprSumm[,-1]
pairedExprSumm = pairedExprSumm[,colnames(seedLev1Counts)]

head(pairedExprSumm)

# outputMark = lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
#                     detach, character.only=TRUE, unload = T)
## If not detach all other packages, the graph will be wrong if run from beginning.

library(oligo)
oligo::MAplot(esetNorm, group = as.factor(pData(esetNorm)$group),
              trans = identity, pairs = T)

## Within knitr, we cannot detach all the packages any more.
## Thus use ggplot to draw boxplot myself, density plot form ggplot shows no advantages.
oligo::hist(esetNorm, trans = identity,
            xlab = "Normalized Intensity",
            ylab = "",
            main = "Density Estimation of Normalized Intensity for Each Sample", 
            cex.lab = 1.8, cex.axis = 1.5, cex.main = 1.8)

# oligo::boxplot(esetNorm, trans = identity, cex.axis = 1, las = 2,
#                main = "" )
# title(ylab = "Normalized Intensity", line = 2)
## Or use # oldpar = par(); par(mgp = c(2,1,0))

temp = melt(exprs(esetNorm))
temp$Var2 = factor(temp$Var2, sort(as.character(unique(temp$Var2))))

geneExpNormBox = ggplot(data = temp,
                        aes(x = Var2, y = value)) +
  stat_boxplot(geom = 'errorbar', coef = 5, width = 0.3) +
  geom_boxplot(aes(fill = Var2), outlier.size = NA, coef = 0) + 
  ylab("Normalized Values") +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.title.y = element_text(size = 18),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 14))
geneExpNormBox


### Immunology ###

wtf = read.table("/home/ratanond/Desktop/Thesis/Prototype/app/immunology/PANTHER_BP.txt", sep = '\t', header = TRUE,
                 stringsAsFactors = FALSE)
# this is the group we're going to work with

# oldpar = par()
# par(mar = c(16,3,3,1))
# # par(mar = c(6,3,3,1))
# panthRank = sort(wtf$Count, index.return = TRUE,
#                  decreasing = TRUE)$ix
# barplot(wtf$Count[panthRank[1:20]], space = 0,
#         main = "Probes assigned to Panther MF's")
# axis(1,1:20-.5, wtf$Term[panthRank[1:20]], las = 2)
# par(oldpar)


immunStuff = which(wtf$Term == "BP00148:Immunity and defense")
immunGenes1List = wtf$Genes[immunStuff]  

## Upgrade to R.3.3.2, no need for the following. 
# lapply(paste('package:', names(sessionInfo()$otherPkgs), sep=""),
#        detach, character.only=TRUE, unload = TRUE) 
# ## If not remove packages, "nchar()" will result in: 
# ## "could not find symbol "keepNA" in environment of the generic function".
# ## It is due to the conflicts of packages and current R version.
# ## If want to permanently solve it, needs to upgrade R itself.

breakList = NULL
for(i in 1:wtf$Count[immunStuff]){
  deliniation = regexpr(", ", immunGenes1List)
  breakList = c(breakList, substr(immunGenes1List, 1, deliniation-1))
  immunGenes1List = substr(immunGenes1List, deliniation + 2,
                           nchar(immunGenes1List))
}
breakList = c(breakList, immunGenes1List) # don't leave the last one behind
#THIS THING ABOVE IS PANTERH IMMUNITY LIST

immunGenes1 = breakList


## ROBB JENIFER IMMUNITY LIST
immunGenes2 = read.table("/home/ratanond/Desktop/Thesis/Prototype/app/immunology/robbs_imunology_genes.csv", sep=',',
                         header=TRUE, skip=0, comment.char = "#",
                         stringsAsFactors=FALSE)

## There's some genes on Robb's list that aren't actually on the codelink.
## So for those, we're just not bothering to look for them obviously
immunGenes2 = immunGenes2$Official.Gene.Name[immunGenes2$codelink == 1]


## LAST ROBB JENIFER IMMUNITY LIST (SECOND LIST) 
## here, we're adding the third list from Robb/Jennifer
immunGenes3 = read.table("/home/ratanond/Desktop/Thesis/Prototype/app/immunology/robbs_imunology_genes_TWO_from_jennifer.csv",
                         sep=',', header=T, skip=0, comment.char = "#",
                         stringsAsFactors=FALSE)

names(immunGenes3)
immunGenes3 = immunGenes3$Official.Gene.Name

## Combine the names of all immunology genes:
immunGeneNames = unique(c(immunGenes1, immunGenes2, immunGenes3))
## Extract the immunology gene list from the observed gene names
immunGeneExpr = pairedExprSumm[which(rownames(pairedExprSumm)
                                     %in% immunGeneNames),]
# length(immunGeneExpr) ## 811


### Start to do sparse CCA for Immunology and SEEDlevel2 ###

library(PMA)
set.seed(1105)
ccaPermImmun = CCA.permute(x = t(seedLev2TestNorm), z = t(immunGeneExpr),
                           typex = "standard", typez = "standard", 
                           nperms = 30, niter = 5)
penXtemp = ccaPermImmun$bestpenaltyx
penZtemp = ccaPermImmun$bestpenaltyz
ccaRsltImmun = CCA(x = t(seedLev2TestNorm), z = t(immunGeneExpr),
                   typex = "standard", typez = "standard",
                   penaltyx = penXtemp, penaltyz = penZtemp,
                   K = 2, niter = 5)
sum(ccaRsltImmun$u != 0)
sum(ccaRsltImmun$v != 0)

ccaImmunScoreU = t(seedLev2TestNorm) %*% ccaRsltImmun$u
ccaImmunScoreV = t(immunGeneExpr) %*% ccaRsltImmun$v
ccaImmunScores = cbind(ccaImmunScoreU, ccaImmunScoreV)
colnames(ccaImmunScores) = c("U1", "U2", "V1", "V2")
ccaImmunScores = as.data.frame(ccaImmunScores)
ccaImmunScores$type = c(rep("BF", 6), rep("FF", 6))


library(ggplot2)
myCCAPlot = function(x = U1, y = U2, col = V1, shape = type, data = ccaImmunScores,
                     xyName = "Seedlevel2", coloName = "Immunology",
                     textVjust = -1.0, elliLev = 0.6, ...){
  jitterPara = list(...)
  if(!"height" %in% names(jitterPara)){
    jitterPara = c(jitterPara, height = 0.01) 
  } else if(!"width" %in% names(jitterPara)){
    jitterPara = c(jitterPara, width = 0.01)
  }
  x = deparse(substitute(x))
  y = deparse(substitute(y))
  col = deparse(substitute(col))
  shape = deparse(substitute(shape))
  myPlot1 = ggplot(data, aes(x = data[,x], y = data[,y],
                             col = data[,col], shape = data[,shape])) +
    geom_point(size = 4) +
    scale_color_continuous(name = paste0("First Component \nScores of ",
                                         coloName),
                           low = "blue", high = "red") +
    geom_text(aes(label = rownames(data)),
              col = "black", size = 5, vjust = textVjust,
              position = do.call("position_jitter", args = jitterPara)) +
    ## The position_jitter will make the values within a group 
    ## a litter bit separate.
    ## On ther other hand, position_dodge will separate the values between groups.
    scale_x_continuous(paste0("First Component Scores of ",
                              xyName)) +
    scale_y_continuous(paste0("Second Component Scores of ",
                              xyName)) +
    labs(title = paste0("Sparse CCA Scores for ", xyName, " as Base")) +
    theme(legend.title = element_text(size = 12),
          plot.title = element_text(size = 16, vjust = 2.0, face = "bold"),
          legend.text = element_text(size = 10)) +
    stat_ellipse(aes(fill = data[,shape]), level = elliLev, alpha = 0.2,
                 geom = "polygon", linetype = 2) +
    scale_fill_discrete(name = "Feeding Type",
                        labels = c("Breastfeeding", "Formula Feeding")) +
    scale_shape_discrete(name = "Feeding Type",
                         labels = c("Breastfeeding", "Formula Feeding"))
  myPlot1
}

## Write output files (1 csv + 2 png)
pdf(file.path("/home/ratanond/Desktop/Thesis/Prototype/app/ccaresults",paste(randname,"SeedLevel2_Scores.pdf", sep="")),bg="transparent")
myCCAPlot()
dev.off()
pdf(file.path("/home/ratanond/Desktop/Thesis/Prototype/app/ccaresults",paste(randname,"Immunology_Scores.pdf", sep="")),bg="transparent")
myCCAPlot(V1, V2, U1, xyName = "Immunology", coloName = "Seedlevel2")
dev.off()

outfile = file.path("/home/ratanond/Desktop/Thesis/Prototype/app/ccaresults",  ## directory is on your side
     	          paste(randname,"-", "ccaImmunScores", ".csv",sep=""))
write.csv(x = ccaImmunScores, file = outfile)

dbSendQuery(conn = db, sprintf("update ccajobsss set status='Complete' where rand='%s'",row[1,3]))
}, error = function(err) {
	dbSendQuery(conn = db, sprintf("update ccajobsss set status='Errored' where rand='%s'",row[1,3]))

})

## close bracket for if(res)
}
Sys.sleep(1)
## close bracket for while(TRUE) loop
}

dbDisconnect(db)
