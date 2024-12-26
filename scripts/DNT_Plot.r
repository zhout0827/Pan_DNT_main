library(readxl)
library(ggplot2)
library(patchwork)
library(dplyr)
library(reshape2)
library(grid)
library(cowplot)
library(gridExtra)
library(data.table)
library(limma)
library(viridis)
library(RColorBrewer)
library(tidyverse)
library(patchwork)
library(networkD3)
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(Seurat)
library(ggpubr)
library(ggrepel)
library(wesanderson)

# mycolor <- wes_palette("FantasticFox1", 8, type ="continuous")
# mycolor <- sample(mycolor, length(mycolor))

#####TCR sanky
df <- fread("DNT_TCR_sanky_input.csv")

p<-ggplot(df, aes(x = x,node= node, next_x= next_x, next_node= next_node,
                   fill = factor(x),
                   label= node)
                   ) +
    geom_sankey(flow.alpha = 0.5,
                # flow.fill = 'grey',
                flow.color = 'grey80', #条带描边色
                # node.fill =mycolor, #节点填充色
                smooth= 8,
                width= 0.08) +
    geom_sankey_text(size = 3.2,color= "black")+
    theme_void()+
#     scale_fill_manual(values = c("#2FA26E","#ECAC03","#A1A377","#BD2923"))+
scale_fill_manual(values = c("#DB8E2D","#DFD356","#A5382A","#63AAC6"))+
    # coord_flip() +
    theme(legend.position = 'none')
p
# ggsave("test.pdf",p,width = 6,height = 20)
ggsave("TCR_sanky.pdf",p,width = 6,height = 5)


#####ROE
meta <- fread("script/DNT_metadata6_new_cancertype.csv")

meta$S %>% table
cluster.order <- c("T13_TRDV1+","T12_TRDV2+","T14_HLA-DR+","T6_CXCL1+CXCL8+","T9_CX3CR1+","T11_ITGA1+","T10_GZMK+","T7_PTGDS+","T4_CXCL13+PDCD1+","T5_FOXP3+IL2RA+","T1_RGS1+CREM+","T3_CD40LG+CCR6+","T2_LEF1+CCR7+","T8_IFNG+")


res.chisq <- chisq.test(table(meta[["DNT_C"]], meta[["Type"]]))
R.oe <- (res.chisq$observed) / (res.chisq$expected)
R.oe_Type <- R.oe[cluster.order, c("Tumor","Nontumor")]

write.table(R.oe_Type, file = "R_oe_Type.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

col_fun2 = colorRamp2(c(0.3,0.5,1,1.5,1.8), 
                        c("#FFFDE0", "#F4E2A5", "#E7B57A","#D57D57","#B7493A"))

pdf('ROE_heatmap.pdf', height=4, width=5*0.6)

Heatmap(R.oe,
        col = col_fun2,
        cluster_rows = FALSE,
        row_order = c(5,4,6,11,3,14, 2,9,10,12,1,8,7,13),
        # cluster_columns = FALSE,
        # row_km = 2,
        # column_km = 4,
        cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(sprintf("%.2f",R.oe[i,j]), x, y, gp = gpar(fontsize = 6,col="black"))
        }
)

dev.off()

####rader

scale_0_1 <- function(x) {
  (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

df <- rbind(fread("Pancancer_T/score/abT_AUCell_Tstate_APC.tsv"),
fread("Pancancer_T/score/gdT_AUCell_Tstate_APC.tsv"))

df1 <- df%>% 
            group_by(DNT_C) %>%
            summarise(
                    Quiescence =    mean(Quiescence),
                    Regulating =    mean(Regulating),
                    Proliferation = mean(Proliferation),
                    Cytotoxicity =  mean(Cytotoxicity),
                    Progenitor_exhaustion = mean(Progenitor_exhaustion),
                    Terminal_exhaustion =   mean(Terminal_exhaustion),
                    Helper =        mean(Helper),
                    Senescence = mean(Senescence),
                    APC =        mean(APC)
                        ) %>% ungroup()
df1

df3 <- df1 %>% melt

df3$DNT_C <- factor(df3$DNT_C, levels = c("T1_RGS1+CREM+","T2_LEF1+CCR7+","T3_CD40LG+CCR6+","T4_CXCL13+PDCD1+","T5_FOXP3+IL2RA+","T6_CXCL1+CXCL8+","T7_PTGDS+","T8_IFNG+","T9_CX3CR1+","T10_GZMK+","T11_ITGA1+","T12_TRDV2+","T13_TRDV1+","T14_HLA-DR+"))

coord_radar <- function (theta = "x", start = 0, direction = 1) 
{  theta <- match.arg(theta, c("x", "y"))
r <- if (theta == "x") 
  "y"
else "x"
ggproto("CoordRadar", CoordPolar, theta = theta, r = r, start = start, 
        direction = sign(direction),
        is_linear = function(coord) TRUE)}

myAngle<- 360- 360 * nrow(df3) /nrow(df3)  

p<-ggplot(data=df3,aes(x=DNT_C, y=value,group=variable,fill=variable)) + 
  geom_polygon(aes(color=variable),alpha=0.2)+
  geom_point(aes(x=DNT_C, y=value,group=variable,fill=variable,color=variable),size=2.5,shape=21)+
  coord_radar()+
  theme_bw() +
facet_grid(vars(variable),scales = "free_y")+
  theme(axis.text.x=element_text(size = 6,colour="black",angle = myAngle),
        axis.title=element_blank(),
        # panel.grid.major = element_line(color="grey80"),
        axis.line = element_blank(),
        panel.border = element_blank(),
        # axis.ticks.x =  element_line(color="black"),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()
        )+
  guides(fill="none",color="none")
p

ggsave("DNT_rader.pdf", p, width = 8, height = 20)


#####T11 SCENIC
# load data
tT <- fread("T11_TF_17vs1_logFC.csv")
df <- tT[grepl("+",rownames(tT),fixed = TRUE),]
df$TF <- rownames(df)
df$abs_logFC <- abs(df$logFC)
df <- filter(df,abs_logFC > 0.1)
df <- df[order(df$logFC,decreasing = TRUE),]
df$TF <- factor(df$TF, levels=rev(df$TF))
select_tf <- c("CTCF(+)","ATF3(+)","CREM(+)","BCL11B(+)","")
df <- df[select_tf,]

p <- ggplot(df,aes(x = TF,y = logFC))+ 
            geom_segment(aes(x=TF,xend=TF,y=0,yend=logFC,color= logFC),linewidth = 1.5)+
            geom_col(aes(fill = logFC), width = 0.1)+
            geom_point(aes(size = log10pvalue, 
            color = logFC))+
            scale_size_continuous(range = c(2, 7))+
            scale_color_continuous_c4a_div('army_rose', mid = 0, reverse = T)+
            scale_fill_continuous_c4a_div('army_rose', mid = 0, reverse = T)+
            scale_y_continuous(breaks = seq(-0.5,0.5,by = 0.25),
            labels = seq(-0.5,0.5,by = 0.25))+
            ylab('')+
            theme_classic()+
            theme(
                  axis.text = element_text(size = 12), #轴标签大小调整
                  axis.title.x = element_text(size = 13), #x轴标题大小调整
                  legend.title = element_text(size = 13), #图例标题大小调整
                  legend.text = element_text(size = 12) #图例标签大小调整
                  )+
            coord_flip()

ggsave(p,filename = "SCENIC_TF_Grm2_VS_Grm1.pdf",width = 8,height = 8)

####vocanol
df_plot <- fread("T14_vocanol_input.csv")
p<-ggplot(
    data=df_plot, aes(x=avg_log2FC_N, y=avg_log2FC_A, color = avg_log2FC_N)) + 
    geom_point(aes(x=avg_log2FC_N, y=avg_log2FC_A, color = avg_log2FC_N, size = -log10(p_val_adj_N)))+
    theme_bw()+
    xlab("Log2FC [Tumor vs Normal]")+
    ylab("Log2FC [Tumor vs Adjacent]")+
    scale_color_gradientn(
        colors = c("#0072B5","#e0d1d1c0","#BC3C29"),  # 定义颜色
        values = rescale(c(-1.4, 0.3, 2)),  # 定义数值对应的颜色节点，
        limits = c(min(df_plot$avg_log2FC_N), max(df_plot$avg_log2FC_N))  # 设置连续变量的范围
    )+
    scale_size(range = c(2,5))+
    theme( 
        panel.border = element_rect(linewidth = 2,fill = NA),
        plot.margin = unit(c(1, 0.5, 0.5, 1), "cm"),
        axis.ticks = element_line(linewidth=1, color = "black"),
        axis.ticks.length = unit(0.2, "cm"),
        axis.text.x=element_text(face="plain",size=12,angle=0,color="Black",vjust=1),
        axis.text.y=element_text(face="plain",size=12,angle=0,color="Black"),
        axis.title.y = element_text(size = 14, face="plain", color = "Black"),
        axis.title.x = element_text(size = 14, face="plain", color = "Black"),
        panel.grid= element_blank()
    )+
    geom_vline(aes(xintercept=0.25),color = 'grey', linetype="dashed")+
    geom_vline(aes(xintercept=-0.25),color = 'grey', linetype="dashed")+
    geom_hline(aes(yintercept=0.25),color = 'grey', linetype="dashed")+
    geom_hline(aes(yintercept=-0.25),color = 'grey', linetype="dashed")+
    geom_text_repel(
        data= df_plot %>% filter(avg_log2FC_N > 1.5 & avg_log2FC_A > 0),
        aes(label = rowname))+
    geom_text_repel(
        data= df_plot %>% filter(avg_log2FC_N < -2 & avg_log2FC_A < 0),
        aes(label = rowname))
p
ggsave("T14_vocanol.pdf",p, width = 7, height = 5)


####boxplot

rds1 <- readRDS("T11.Rds")

p<-ggplot(rds1@meta.data %>% filter(subcluster_type %in% c("Grm1","Grm17")), aes(x = subcluster_type, y = Hypoxia.HIF_regulated))+
        geom_boxplot(aes(fill = subcluster_type), outlier.size = 0,outlier.shape = NA,
        position = position_dodge(0.8), lwd = 0.7, width = 0.6)+
            scale_fill_manual(values = c('#d33f6a',
                                  '#E17B84',
                                  '#8dd593',
                                  '#f0b98d',
                                  '#ef9708',
                                  '#b5bbe3',
                                  '#4a6fe3',
                                '#023fa5',
                                '#c6dec7',
                                '#9cded6',
                                '#BEF197',
                                '#BFC1D3',
                                '#d6bcc0',
                                '#30E0C0'))+
        theme_bw()+
        ylab("Scores")+
        xlab(" ")+
        theme( 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(linewidth=1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x=element_text(face="plain",size=11,angle=90,color="Black",vjust=0.5,hjust=1),
        axis.text.y=element_text(face="plain",size=11,angle=0,color="Black"),
        axis.title.y = element_text(size = 13, face="plain", color = "Black"),
        axis.title.x = element_text(size = 13, face="plain", color = "Black"),
        panel.grid= element_blank()
        )+
        guides(fill="none")p<-ggplot(rds1@meta.data %>% filter(subcluster_type %in% c("Trm1","Trm17")), aes(x = subcluster_type, y = Hypoxia.HIF_regulated))+
        geom_boxplot(aes(fill = subcluster_type), outlier.size = 0,outlier.shape = NA,
        position = position_dodge(0.8), lwd = 0.7, width = 0.6)+
            scale_fill_manual(values = c('#d33f6a',
                                  '#E17B84',
                                  '#8dd593',
                                  '#f0b98d',
                                  '#ef9708',
                                  '#b5bbe3',
                                  '#4a6fe3',
                                '#023fa5',
                                '#c6dec7',
                                '#9cded6',
                                '#BEF197',
                                '#BFC1D3',
                                '#d6bcc0',
                                '#30E0C0'))+
        theme_bw()+
        ylab("Scores")+
        xlab(" ")+
        theme( 
        panel.border = element_blank(),
        axis.line = element_line(color = "black", size = 1),
        axis.ticks = element_line(linewidth=1, color = "black"),
        axis.ticks.length = unit(0.15, "cm"),
        axis.text.x=element_text(face="plain",size=11,angle=90,color="Black",vjust=0.5,hjust=1),
        axis.text.y=element_text(face="plain",size=11,angle=0,color="Black"),
        axis.title.y = element_text(size = 13, face="plain", color = "Black"),
        axis.title.x = element_text(size = 13, face="plain", color = "Black"),
        panel.grid= element_blank()
        )+
        guides(fill="none")

ggsave("T11_hypoxia_boxplot.pdf",p, width = 7, height = 5)

####survival
library(survival)
library(survminer)

df <- 
fread("/ICB_cohort/CIBERSORTx_Adjusted.txt") %>%
  left_join(fread("ICB_cohort/BLCA_Sanjeev_2018_metadata.csv"), by = c("Mixture"="V1")
) %>% 
  mutate(
    Res = case_when(Response %in% c("CR","PR") ~ "Response",
                    Response %in% c("PD","SD") ~ "Non_Response"
                      )) %>%
    filter(!is.na(Res))

res.cut <- surv_cutpoint(df ,
                        time = "os",
                        event = "censOS",
                        variables = "T14_HLA-DR+")
fit <- survfit(Surv(os, censOS) ~ T14, data = res.cat)
ggsurv <- ggsurvplot(
    fit, 
    data =  res.cat,
    size = 1,
    palette = c("#A6B0B1", "#2E9FDF"),
    pval = TRUE, 
    risk.table = TRUE, 
    xlab = "Time in Days", 
    ylab = "Overall Survival",
    legend.title = "BLCA ICB cohort",
    legend.labs = c("High","Low"),
    risk.table.title = "Overall Survival",
    ggtheme = theme_bw()
    )

p <- ggsurv$plot +
 guides(
   color = guide_legend(
     override.aes = list(size=1),
     keywidth= unit(0.5,'cm'),
     keyheight= unit(0.5,'cm')
   )
 )+
 theme( 
    panel.border = element_blank(),
    axis.line = element_line(color = "black", size = 1),
    plot.margin = unit(c(1, 0.5, 0.5, 1)/2, "cm"),
    axis.ticks = element_line(linewidth=1, color = "black"),
    axis.ticks.length = unit(0.15, "cm"),
    axis.text.x=element_text(face="plain",size=12,angle=0,color="Black",vjust=1),
    axis.text.y=element_text(face="plain",size=12,angle=0,color="Black"),
    axis.title.y = element_text(size = 15, face="plain", color = "Black"),
    axis.title.x = element_text(size = 15, face="plain", color = "Black"),
    panel.grid= element_blank()
 )

ggsave(filename = 'T14_ICB_surv.pdf',p,width = 4,height = 5)