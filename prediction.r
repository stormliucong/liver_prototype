if(!exists("train")){
  library(DESeq2) ### 1.12.4
  library(glmnet) ### 2.0-5
  load("Model_liver_v1.rda") # load the model 
  input <- read.table("prediction_file_list.txt",sep = "\t")
  # column 1: SEQ_ID.
  # column 2: 血液编号.
  # column 3: 收样医院
  # column 4: Age: <35,35-55,55
  # column 5: Diagnois.
  # column 6: featureCount_file.
  # column 6-?: Other covariates. TBD.
  # seperate by "\t"
  # no head.
  seq_id <- paste("SEQ",input$V1,sep = "")
  file_name_list <- as.character(input$V6)
  file_name_list <- paste(file_name_list,"/",seq_id,".genebody",sep = "")
  colData_test <- input
}else{
  file_name_list <- as.character(input$V6)[-train]
  file_name_list <- paste(file_name_list,"/",seq_id[-train],".genebody",sep = "")
  colData_test <- input[-train,]
}

out <- rep(NA, length(file_name_list))
for(i in 1:length(file_name_list))
{
  out[i] <- pred( file_name_list[i])
  cat(i,"\n")
}
score <- out
result <- cbind(colData_test,score)
write.table(result,"pred_result.txt",sep = "\t",row.names = F)
