library(DESeq2) ### 1.12.4
library(glmnet) ### 2.0-5
count <- function(x){
  counts <- read.delim(x, skip=1, row.names=1)
  counts <- counts[gene_list,6]
  return(counts)
}

sizeCount <- function(x){
  counts <- read.delim(x, skip=1, row.names=1)
  counts <- as.data.frame(counts[names(geoMeans),6])
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(condition="unknown"), design=~1) 
  dds <- estimateSizeFactors(dds, geoMeans=geoMeans) # where the problem happens.
  #mcols(dds)$dispFit <- dispFit
  #rld <- rlog(dds, blind=F, betaPriorVar=betaPriorVar, intercept=intercept, fitType="parametric") 
  #rld <- assay(rld) [match(gene_list, names(geoMeans)),]
  #return(rld)
  return(as.numeric(dds$sizeFactor))
}

rldCount <- function(x){
  counts <- read.delim(x, skip=1, row.names=1)
  counts <- as.data.frame(counts[names(geoMeans),6])
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(condition="unknown"), design=~1) 
  dds <- estimateSizeFactors(dds, geoMeans=geoMeans) # where the problem happens.
  mcols(dds)$dispFit <- dispFit
  rld <- rlog(dds, blind=F, betaPriorVar=betaPriorVar, intercept=intercept, fitType="parametric") 
  rld <- assay(rld) [match(gene_list, names(geoMeans)),]
  return(rld)
}

pred <- function(x){
  counts <- read.delim(x, skip=1, row.names=1)
  counts <- counts[names(geoMeans),6]
  dds <- DESeqDataSetFromMatrix(countData=counts, colData=data.frame(condition="unknown"), design=~1) 
  dds <- estimateSizeFactors(dds, geoMeans=geoMeans) 
  mcols(dds)$dispFit <- dispFit
  rld <- rlog(dds, blind=F, betaPriorVar=betaPriorVar, intercept=intercept, fitType="parametric") 
  rld <- assay(rld) [match(gene_list, names(geoMeans)),]
  return( predict(cvfit, newx=t(rld), s="lambda.1se", type="response") )
}

load("/database/scripts/get_score/liver_prototype/Model_liver_v1.rda") # load the model 
b <- as.matrix(coef(cvfit,s="lambda.1se"))
gene_list <- names(b[which(b!=0),][-1])
input <- read.table("prediction_file_list.txt",sep = "\t")
# column 1: SEQ_ID.
# column 2: 血液编号.
# column 3: DNA_type
# column 4: Diagnois.
# column 5: featureCount_file.
# seperate by "\t"
# no head.
seq_id <- paste("SEQ",input$V1,sep = "")
file_name_list <- as.character(input$V7)
file_name_list <- paste(file_name_list,"/",seq_id,".genebody",sep = "")
colData_test <- input
raw_out <- NULL
rld_out <- NULL
pred_out <- NULL
size_out <- NULL
for(i in 1:length(file_name_list))
{
  size_out <- rbind(size_out,sizeCount(file_name_list[i]))
  raw_out <- rbind(raw_out,count(file_name_list[i]))
  rld_out <- rbind(rld_out,rldCount(file_name_list[i]))
  pred_out <- rbind(pred_out,pred(file_name_list[i]))
  cat(i,"\n")
}
colnames(raw_out) <- gene_list
colnames(rld_out) <- gene_list
colnames(pred_out) <- "score"
colnames(size_out) <- "factorSize"

result <- cbind(colData_test,pred_out,raw_out,rld_out,size_out)
write.table(result,file="gDNA_explore.txt",sep="\t")