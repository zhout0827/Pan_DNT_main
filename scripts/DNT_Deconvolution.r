library(CARD)
library(MuSiC)
library(Seurat)
library(patchwork)
library(tidyverse)
library(progress)

####CARD
scRNA <- schard::h5ad2seurat("/work/haoq/Pancancer_T/h5ad/ALL_Epi_DNT_T_M_counts.h5ad")
sc_counts <- as.matrix(scRNA[['RNA']]@counts)
sc_meta <- scRNA@meta.data %>% dplyr::select(batch,orig.ident,cell_type)

N <- 100
pb <- progress_bar$new(
 format = '[:bar] :percent eta: :eta',
 total = nrow(folders), clear = FALSE, width = 80
)
####Read10X
stRNAlist <- list()
CARD_obj_list <- list()
for (i in 1:nrow(folders)) {
  group <- folders$group[i]
  type <- folders$type[i]
  # 构建预期的文件夹路径
  expected_folder_path <- file.path(".", paste0(group,folders$dir[i]))
  # 判断是否含有lowres文件
  if (type == "A") {
  # 读取数据
    stRNA = Load10X_Spatial(expected_folder_path,slice = folders$sample[i])
  } else {
    img_dir <- file.path(expected_folder_path,"spatial")
    img = Read10X_Image(img_dir,image.name = "tissue_hires_image.png")
    stRNA = Load10X_Spatial(expected_folder_path, image = img,slice = folders$sample[i])
  }
    
    stRNA <- SCTransform(stRNA, assay = "Spatial", verbose = FALSE)
    stRNA <- RunPCA(stRNA, assay = "SCT", verbose = FALSE) 
    pc.num=1:20
    
    stRNA <- FindNeighbors(stRNA, reduction = "pca", dims = pc.num)
    stRNA <- FindClusters(stRNA, verbose = FALSE)
   
    stRNA <- RunUMAP(stRNA, reduction = "pca", dims =  pc.num)
   
    coords <- GetTissueCoordinates(stRNA,cols=c("row","col"),scale=NULL)
    colnames(coords)=c('x','y')
    coords <- coords[,1:2]

    sp_counts <- GetAssayData(stRNA,layer = "counts")

    if(!ncol(sp_counts) == nrow(coords)){
      sp_counts <- sp_counts[,intersect(rownames(coords),colnames(sp_counts))]
    }

    CARD_obj = createCARDObject(sc_count = sc_counts,
                                sc_meta = sc_meta,
                                spatial_count = sp_counts,
                                spatial_location = coords,
                                ct.varname = "cell_type",
                                ct.select = unique(sc_meta$cell_type), 
                                sample.varname = "orig.ident")

    CARD_obj = CARD_deconvolution(CARD_object = CARD_obj)
    CARD_obj_list[[i]] <- CARD_obj
    stRNAlist[[i]] <- stRNA

pb$tick()
}

saveRDS(CARD_obj_list,"CARD_obj_list.Rds")
saveRDS(stRNAlist,"stRNAlist.Rds")

stRNA[['CARD']] <- CreateAssayObject(counts = t(composition))
DefaultAssay(stRNA) <- "CARD"
SpatialFeaturePlot(stRNA, 
                   features = c("T1-RGS1+CREM+","T11-ITGA1+","T14-HLA-DR+"), 
                   image.alpha = 0,
                   pt.size.factor = 1.6,
                   combine = TRUE
                   )


####10X Visium HD
stRNA <- Load10X_Spatial("/work/zhout/pancancer//DNT/spatial/VisiumHD",bin.size = c(8,16))
DefaultAssay(hd_data) <- "Spatial.008um"

###RCTD
scRNA <- schard::h5ad2seurat("/work/haoq/Pancancer_T/h5ad/ALL_Epi_DNT_T_M_counts.h5ad")
sc_counts <- as.matrix(scRNA[['RNA']]@counts)
sc_nUMI = colSums(sc_counts)
cellType=data.frame(barcode=colnames(scRNA), celltype=scRNA$cell_type)
names(cellType) = c('barcode', 'cell_type')
cell_types = cellType$cell_type; names(cell_types) <- cellType$barcode 
cell_types = as.factor(cell_types)
reference = Reference(sc_counts, cell_types, sc_nUMI)

sp_counts<-Read10X_h5("pancancer/DNT/spatial/VisiumHD/binned_outputs/square_008um/filtered_feature_bc_matrix.h5")
coords<-read_parquet("pancancer/DNT/spatial/VisiumHD/binned_outputs/square_008um/spatial/tissue_positions.parquet",as_data_frame = TRUE)
coords <- as.data.frame(coords)
rownames(coords)<-coords$barcode
coords<-coords[colnames(sp_counts),]
coords<-coords[,3:4]
sp_nUMI <- colSums(sp_counts)

puck <- SpatialRNA(coords, sp_counts, sp_nUMI)
barcodes <- colnames(puck@counts)

myRCTD <- create.RCTD(puck,reference,max_cores = 48)

myRCTD <- run.RCTD(myRCTD,doublet_mode = 'doublet')