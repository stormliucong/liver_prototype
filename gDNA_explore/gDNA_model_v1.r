library(DESeq2) ### 1.12.4
library(glmnet) ### 2.0-5
#source("https://bioconductor.org/biocLite.R")
#biocLite("DESeq2")
#################################################
# STEP 1 get INPUT
#################################################
# column 1: SEQ_ID.
# column 2: 血液编号
# column 3: 提取DNA类型.
# column 4: 收样医院
# column 5: Age: <35,35-55,55
# column 6: Diagnois.
# column 7: featureCount_file.
# column 8-?: Other covariates. TBD.
# seperate by "\t"
# no head.
cat("STEP1: get input.\n")
input <- read.table("training_sample.txt",sep = "\t")
#input <- input[input$V4 != 'NULL',]
seq_id <- paste("SEQ",input$V1,sep = "")
# train <- sample(1:dim(input)[1],size = 60)
train <- sample(1:dim(input)[1])


#################################################
# STEP 2 go through files to build big feature count matrix.
#################################################
if(!exists("big_mat")){
  file_name_list <- as.character(input$V7)
  N <- length(file_name_list)
  # get gene_id from file_1.
  # use genebody feature count.
  file_name_list <- paste(file_name_list,"/",seq_id,".genebody",sep = "")
  file_1 <- file_name_list[1]
  gene_id <- as.character(read.table(file_1,head=TRUE,sep = "\t",skip = 1)[,1])
  
  count <- NULL
  for(i in 1:length(seq_id)){
    file <- file_name_list[i]
    count_now <- read.table(file,head=TRUE,sep = "\t",skip = 1)[,7]
    count <- cbind(count,count_now)
  }
  big_mat <- data.frame(gene_id,count)
  colnames(big_mat) <- c("gene_name",as.character(seq_id))
  #save(big_mat,file = "big_mat.RData")
}else{
  #load("big_mat.RData")
}
cat("STEP2: load big_mat.\n")


#################################################
# STEP 3 prepare condition matrix for Deseq2.
#################################################

# countData is the input format for Deseq2.
countData <- big_mat[,-1]
# trainning sample.
countData <- countData[,train]
rownames(countData) <- big_mat$gene_name
countData <- countData[which(rowMeans(countData) > 0),]

# read condition data.
colData <- as.matrix(input[train,c(3)])
# group <- cut(as.numeric(colData$V4), 
#              breaks = c(-Inf, 35, 55, Inf), 
#              labels = c("35-", "35-55", "55+"), 
#              right = FALSE)
rownames(colData) <- seq_id[train]
colnames(colData) <- c("condition")

cat("STEP3: get colData and countData for DESeq2.\n")

#################################################
# STEP 4 Run Deseq2.
#################################################
dds <- DESeqDataSetFromMatrix(countData=countData,
                              colData=colData,
                              design=~condition)
dds <- DESeq(dds)
res <- results(dds)
resOrdered <- res[order(res$padj),]
#resSig <- subset(resOrdered, padj < 0.05)
resSig <- resOrdered[1:100,]
# record useful parameters.
gene.list <- rownames(resSig)
geoMeans <- apply(countData,1,function(x) exp(mean(log(x+0.1))))
dispFit <- mcols(dds)$dispFit

#rld <- rlog(dds, blind=FALSE,dispFit =dispFit) 
# rld is very slow. It might takes several hours.
# we could afford to take 30~40 subset of samples.
rld <- rlog(dds[,sample(1:dim(dds)[2],20)], blind=FALSE) 
# vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
betaPriorVar <- attr(rld,"betaPriorVar")
intercept <- mcols(rld)$rlogIntercept
rld <- rlog(dds, blind=F, betaPriorVar=betaPriorVar, intercept=intercept, fitType="parametric") 

cat("STEP4: DESeq2 finished.\n")

#################################################
# STEP 5 Run glmnet
#################################################
rld <- assay(rld)[match(gene.list, names(geoMeans)),]
cvfit <- cv.glmnet(x = t(rld), y = colData[,1],family = "binomial", type.measure = "class")
pred <- function(x){
  counts <- read.delim(x, skip=1, row.names=1)
  counts <- counts[names(geoMeans),6]
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(condition="unknown"), design=~1) 
  dds <- estimateSizeFactors(dds, geoMeans=geoMeans) 
  mcols(dds)$dispFit <- dispFit
  rld <- rlog(dds, blind=F, betaPriorVar=betaPriorVar, intercept=intercept, fitType="parametric") 
  rld <- assay(rld) [match(gene.list, names(geoMeans)),]
  return( predict(cvfit, newx=t(rld), s="lambda.1se", type="response") )
}
#################################################
# STEP 6 Save the Model
#################################################
save(betaPriorVar,
     cvfit,
     dispFit,
     gene.list,
     geoMeans,
     intercept,dds,
     pred=pred,file = "Model_gDNA_v1.rda")


