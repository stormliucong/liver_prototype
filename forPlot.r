# compare median value.
require(ggplot2)
require(DESeq2) ### 1.12.4
require(glmnet) ### 2.0-5
load("Model_liver_v1.rda")
res <- results(dds)
d <- plotCounts(dds,gene = which(res$padj < 0.00001),intgroup = "condition",replaced = T)
