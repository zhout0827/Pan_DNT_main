#LoadSeuratlibrary
library(Seurat)

input_file<-commandArgs(trailingOnly=TRUE)[1]

if(!file.exists(input_file)){
stop("Filedoesnotexist:",input_file)
}

load(file=input_file)

seurat_objs<-Filter(Negate(is.null),Filter(function(x)inherits(get(x),"Seurat"),ls()))

if(length(seurat_objs)==0){
stop("NoSeuratobjectfoundintheloadedenvironment.")
}

print(seurat_objs)

seurat_obj_name<-seurat_objs[1]

output_dir<-"output"
if(!dir.exists(output_dir)){
dir.create(output_dir,recursive=TRUE)
}

seurat_obj<-get(seurat_obj_name)
write.csv(t(as.matrix(seurat_obj[['RNA']]@counts)),file=paste0(output_dir,"/count.csv"))
