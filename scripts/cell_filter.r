# library(tidyverse)
library(data.table)
library(Seurat)
library(getopt)

spec <- matrix(c(
    "name", "n", 1, "character",
    "help", "h", 0, "logical"
), byrow = TRUE, ncol = 4)
opt <- getopt(spec)

name <- opt$name

# l <- fread("level2.sample_list", header = FALSE)
l <- fread(name, header = FALSE)
print(name)

for (i in 1:nrow(l)){
    print(l[i])

    data <- readRDS(paste0("/mnt/xu/scRNAseq/SingleCell_PublicData/sample_rds_level2/",l[i],".RDS"))

    row1 = data@assays$RNA[1,]
    col1 = data@assays$RNA[,1]
    Non1 <- row1[!as.integer(row1) ==  row1]
    Non2 <- col1[!as.integer(col1) ==  col1]
    
    error_flag <- FALSE
    # 检查是否是原始counts
    if (length(Non1) > 0) {
        error_flag <- TRUE
        warning(paste("Raw Counts not found:", paste(Non1, collapse = ", ")))
    }
    if (length(Non2) > 0) {
        error_flag <- TRUE
        warning(paste("Raw Counts not found:", paste(Non2, collapse = ", ")))
    }
    
    ### filter
    tryCatch({
        Cell2 <- WhichCells(data, expression = (CD3D > 0 | CD3E > 0 | CD3G > 0))
        
        if (length(Cell2) > 0) {
            #取出符合表达的细胞
            fil2 <- subset(data, cells = Cell2)
            #保存
            saveRDS(fil2, file = paste0("/work/haoq/Pancancer_T/rds/level2_T_rds/",l[i],".rds"))
            print("File: Saved.")
        } else {
                print("File: no cell. Skipping.")
                error_flag <- TRUE
        }

    }, error = function(e) {
        print(paste("error:", conditionMessage(e), ". Skipping."))
        error_flag <- TRUE
    })

    if (error_flag) {
        next
    }
}
