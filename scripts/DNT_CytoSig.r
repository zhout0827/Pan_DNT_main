
T11 <- readRDS("/work/haoq/Pancancer_T/rds/T11.rds")
dir_for_data <- 
DropletUtils::write10xCounts(
    x = T11@assays$RNA@counts,
    path = file.path(dir_for_data, "01.cellranger_output"),
    version = "3"
)
dir.create(file.path(dir_for_data, "02.output"), recursive = TRUE)
system("python3 CytoSig_run.py -i /work/zhout/pancancer/DNT/CytoSig/01.cellranger_output/ -o /work/zhout/pancancer/DNT/CytoSig/02.output -c 0 -z 1", wait = FALSE)

####Grm2 VS Grm1
cellinfo <- T11@meta.data
targets<-data.table(FileName=rownames(cellinfo),Target=cellinfo$T11_type)

lev<-unique(targets$Target)
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f)
colnames(design) <- lev 

eset=getAUC(anaAUC)
eset=eset[,targets$FileName]
eset<-t(scale(t(eset)))


cont.wt <- makeContrasts("Grm2-Grm1",levels=design) 
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)

tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
write.table(tT,file="T11_cytokine_2vs1_logFC.csv",sep="\t",quote=F)

####Tumor VS Nontumor
targets<-data.table(FileName=rownames(cellinfo),Target=cellinfo$State)

lev<-unique(targets$Target)
f <- factor(targets$Target, levels=lev) 
design <- model.matrix(~0+f) 
colnames(design) <- lev 

eset <- t(zscore_all) 
eset <- eset[,targets$FileName]
eset <- t(scale(t(eset))) 

cont.wt <- makeContrasts("Tumor-Nontumor",levels=design)
fit <- lmFit(eset, design)
fit2 <- contrasts.fit(fit, cont.wt) 
fit2 <- eBayes(fit2) 
tT=topTable(fit2, adjust="BH",sort.by="logFC",n=Inf)

tT = subset(tT, select=c("adj.P.Val","P.Value","logFC"))
colnames(tT)=c("FDR","P.Value","logFC")
write.table(tT,file="T11_cytokine_TvsN_logFC.csv",sep="\t",quote=F)