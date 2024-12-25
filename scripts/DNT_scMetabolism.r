library(scMetabolism)
library(ggplot2)
library(rsvd)


countexp.Seurat <- readRDS("Pancancer_T/rds/T11.rds")

scRNA_metabolism <- sc.metabolism.Seurat(obj = countexp.Seurat, method = "VISION", imputation = F, ncores = 24, metabolism.type = "KEGG")

DimPlot.metabolism(obj = countexp.Seurat, 
                   pathway = rownames(scRNA_metabolism@assays[["METABOLISM"]][["score"]])[1:20], 
                   dimention.reduction.type = "umap", 
                   dimention.reduction.run = F, size = 1,
                   phenotype = "T11_type", norm = "y")


