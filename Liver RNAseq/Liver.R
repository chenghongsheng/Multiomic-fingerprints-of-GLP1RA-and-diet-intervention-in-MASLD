setwd("your-working-directory")

library(DESeq2)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(biomaRt)
library(pheatmap)
library(ViSEAGO)
library(EnhancedVolcano)
library(stringr)
library(ggfortify)
library(ggforce)
library(GSVA)
library(rstatix)

##### load data ######
#metadata
metadata <- read.delim("./metadata.txt", header=T, row.names = 1)

metadata_liver<-subset(metadata,metadata$Tissue %in% c('Liver'))
#metadata_liver <- metadata_liver[,c("Group","Tissue")]
rownames(metadata_liver) <- sub("-", ".", rownames(metadata_liver))
metadata_liver$Group <- factor(metadata_liver$Group,
                               levels=c('Control','Control diet + Semaglutide',
                                        'LIDPAD','LIDPAD + Semaglutide','LIDPAD + diet reversion'))

#metadata_liver <- metadata_liver[-c(8,22), ]

#count
count <- read.table("./Fcount.txt", row.names=1, header = T)
count<-count[,c(6:53)]

rownames(metadata) <- sub("-", ".", rownames(metadata))
count<-count[,rownames(metadata)]

count_liver<-count[,rownames(metadata_liver)]

#sanity check
all(colnames(count_liver) %in% rownames(metadata_liver))
all(rownames(metadata) == colnames(count))

#### set color ######


color4liver<-c('black','grey50','red','#EBB8DD','#5C62D6')
pie(rep(1,5),col=color4liver)

##### DESeq_liver-first trial #####
#conducts normalization
dds_liver <- DESeqDataSetFromMatrix(countData = count_liver,
                                    colData = metadata_liver,
                                    design = ~ Group)
dds_liver<-DESeq(dds_liver)
vsd_liver<-varianceStabilizingTransformation(dds_liver, blind=TRUE)

pcadat_liver<-plotPCA(vsd_liver,intgroup="Group",ntop=3000,returnData=T)

all(rownames(pcadat_liver)==rownames(metadata_liver)) #sanity check
pcadat_liver$Group<-metadata_liver$Group

percentVar.vsd_liver<-round(100*attr(pcadat_liver,"percentVar"))

pcadat_liver$name<-str_split_i(pcadat_liver$name,'_',1)
pcadat_liver$name<-str_split_i(pcadat_liver$name,'\\.',2)


pca_liver<-ggplot(pcadat_liver, aes(PC1, PC2, fill=Group)) +
  geom_point(size=5,pch=21,stroke=0.5)+ 
  geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd_liver[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd_liver[2],"% variance")) + 
  #stat_ellipse()+
   #geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = color4liver)+
  scale_color_manual(values = color4liver)+
  coord_fixed()


pca_liver


###### DEGs Liver ####
Hsa.dataset<-useDataset('mmusculus_gene_ensembl',mart=useMart("ensembl"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name',"gene_biotype"), 
               filters='ensembl_gene_id',
               values=rownames(count),mart=Hsa.dataset)
genesymbols <- tapply(Genemap$external_gene_name, 
                      Genemap$ensembl_gene_id, paste, collapse="; ")

#1 control vs LIDPAD
res.liver.ctrlvLIDPAD <- as.data.frame(results(dds_liver, contrast=c('Group','Control','LIDPAD'),alpha=0.05)) 
res.liver.ctrlvLIDPAD<-res.liver.ctrlvLIDPAD[order(res.liver.ctrlvLIDPAD$padj),]
res.liver.ctrlvLIDPAD$symbol<-genesymbols[rownames(res.liver.ctrlvLIDPAD)]
sig.liver.ctrlvLIDPAD<-res.liver.ctrlvLIDPAD[which(res.liver.ctrlvLIDPAD$padj<0.05 &
                                       (abs(res.liver.ctrlvLIDPAD$log2FoldChange)>1)),]

write.csv(sig.liver.ctrlvLIDPAD,"sig.liver.ctrlvLIDPAD.csv")
write.csv(res.liver.ctrlvLIDPAD,"res.liver.ctrlvLIDPAD.csv")

#2 control vs ctrl_semaglutide (CSG)
res.liver.ctrlvCSG <- as.data.frame(results(dds_liver, contrast=c('Group','Control','Control diet + Semaglutide'),alpha=0.05)) 
res.liver.ctrlvCSG<-res.liver.ctrlvCSG[order(res.liver.ctrlvCSG$padj),]
res.liver.ctrlvCSG$symbol<-genesymbols[rownames(res.liver.ctrlvCSG)]
sig.liver.ctrlvCSG<-res.liver.ctrlvCSG[which(res.liver.ctrlvCSG$padj<0.05 &
                                                     (abs(res.liver.ctrlvCSG$log2FoldChange)>1)),]

write.csv(sig.liver.ctrlvCSG,"sig.liver.ctrlvCSG.csv")
write.csv(res.liver.ctrlvCSG,"res.liver.ctrlvCSG.csv")

#3 LIDPAD vs LP_semaglutide (LPSG)
res.liver.LIDPADvLPSG <- as.data.frame(results(dds_liver, contrast=c('Group','LIDPAD','LIDPAD + Semaglutide'),alpha=0.05)) 
res.liver.LIDPADvLPSG<-res.liver.LIDPADvLPSG[order(res.liver.LIDPADvLPSG$padj),]
res.liver.LIDPADvLPSG$symbol<-genesymbols[rownames(res.liver.LIDPADvLPSG)]
sig.liver.LIDPADvLPSG<-res.liver.LIDPADvLPSG[which(res.liver.LIDPADvLPSG$padj<0.05 &
                                                   (abs(res.liver.LIDPADvLPSG$log2FoldChange)>1)),]

write.csv(sig.liver.LIDPADvLPSG,"sig.liver.LIDPADvLPSG.csv")
write.csv(res.liver.LIDPADvLPSG,"res.liver.LIDPADvLPSG.csv")

#4 LIDPAD vs LP_dietreversion (LPDR)
res.liver.LIDPADvLPDR <- as.data.frame(results(dds_liver, contrast=c('Group','LIDPAD','LIDPAD + diet reversion'),alpha=0.05)) 
res.liver.LIDPADvLPDR<-res.liver.LIDPADvLPDR[order(res.liver.LIDPADvLPDR$padj),]
res.liver.LIDPADvLPDR$symbol<-genesymbols[rownames(res.liver.LIDPADvLPDR)]
sig.liver.LIDPADvLPDR<-res.liver.LIDPADvLPDR[which(res.liver.LIDPADvLPDR$padj<0.05 &
                                                     (abs(res.liver.LIDPADvLPDR$log2FoldChange)>1)),]

write.csv(sig.liver.LIDPADvLPDR,"sig.liver.LIDPADvLPDR.csv")
write.csv(res.liver.LIDPADvLPDR,"res.liver.LIDPADvLPDR.csv")

#5 ctrl_semaglutide (CSG) vs LP_semaglutide (LPSG)
res.liver.CSGvLPSG <- as.data.frame(results(dds_liver, contrast=c('Group','Control diet + Semaglutide','LIDPAD + Semaglutide'),alpha=0.05)) 
res.liver.CSGvLPSG<-res.liver.CSGvLPSG[order(res.liver.CSGvLPSG$padj),]
res.liver.CSGvLPSG$symbol<-genesymbols[rownames(res.liver.CSGvLPSG)]
sig.liver.CSGvLPSG<-res.liver.CSGvLPSG[which(res.liver.CSGvLPSG$padj<0.05 &
                                                     (abs(res.liver.CSGvLPSG$log2FoldChange)>1)),]

write.csv(sig.liver.CSGvLPSG,"sig.liver.CSGvLPSG.csv")
write.csv(res.liver.CSGvLPSG,"res.liver.CSGvLPSG.csv")

#6 LP_semaglutide (LPSG) vs LP_dietreversion (LPDR)
res.liver.LPSGvLPDR <- as.data.frame(results(dds_liver, contrast=c('Group','LIDPAD + Semaglutide','LIDPAD + diet reversion'),alpha=0.05)) 
res.liver.LPSGvLPDR<-res.liver.LPSGvLPDR[order(res.liver.LPSGvLPDR$padj),]
res.liver.LPSGvLPDR$symbol<-genesymbols[rownames(res.liver.LPSGvLPDR)]
sig.liver.LPSGvLPDR<-res.liver.LPSGvLPDR[which(res.liver.LPSGvLPDR$padj<0.05 &
                                               (abs(res.liver.LPSGvLPDR$log2FoldChange)>1)),]

write.csv(sig.liver.LPSGvLPDR,"sig.liver.LPSGvLPDR.csv")
write.csv(res.liver.LPSGvLPDR,"res.liver.LPSGvLPDR.csv")



#plot heatmap
liver.DEG.all<-c(rownames(sig.liver.ctrlvLIDPAD),rownames(sig.liver.ctrlvCSG),
                 rownames(sig.liver.LIDPADvLPSG),rownames(sig.liver.LIDPADvLPDR),
                 rownames(sig.liver.CSGvLPSG),
                 rownames(sig.liver.LPSGvLPDR))
liver.DEG.all<-liver.DEG.all[!duplicated(liver.DEG.all)]

hm_mat.liver<-assay(vsd_liver)
hm_mat.liver<-subset(hm_mat.liver,rownames(hm_mat.liver) %in% liver.DEG.all)

coldat_hm.liver<-metadata_liver[order(metadata_liver$Group),,drop=F]

hm_mat.liver<-hm_mat.liver[,rownames(coldat_hm.liver)]
hm_mat.liver<-as.data.frame(hm_mat.liver)

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
break_hm = seq(-2, 2,length.out=100)

my_color_annotation.liver<-list(Group= c('Control'='black','LIDPAD'='red',
                                         'Control diet + Semaglutide'='grey50',
                                         'LIDPAD + Semaglutide'='#EBB8DD',
                                         'LIDPAD + diet reversion'='#5C62D6'))

coldat_hm.liver<-as.data.frame(dplyr::select(coldat_hm.liver,-Tissue))

pheatmap(hm_mat.liver,scale="row",border_color = NA,color = hm_color,
                     show_rownames = F,show_colnames = T,
                     cluster_rows = T,cluster_cols =F,
                     annotation_col = coldat_hm.liver,
                     breaks = break_hm,
                     annotation_colors = my_color_annotation.liver,
                     clustering_distance_rows = "correlation",
                     clustering_distance_cols = "correlation",
                     angle_col = 45,cutree_rows = 5, gaps_col = 10,
                     fontsize_col = 5) 

##### Remove outliers DEG ####
metadata_liver.sub<-metadata_liver 
metadata_liver.sub$SN<-rownames(metadata_liver.sub)

metadata_liver.sub$SN<-str_split_i(metadata_liver.sub$SN,'_',1)
metadata_liver.sub$SN<-str_split_i(metadata_liver.sub$SN,'\\.',2)


out<-c(1541,1551,
       #1542,
       1543,1555, 1556,1557,1559)

metadata_liver.sub<-subset(metadata_liver.sub,!metadata_liver.sub$SN %in% out)
count_liver.sub<-count_liver[,rownames(metadata_liver.sub)]

all(rownames(metadata_liver.sub)==colnames(count_liver.sub))

# deseq2 analysis
dds_liver.sub <- DESeqDataSetFromMatrix(countData = count_liver.sub,
                                    colData = metadata_liver.sub,
                                    design = ~ Group)
dds_liver.sub<-DESeq(dds_liver.sub)
vsd_liver.sub<-varianceStabilizingTransformation(dds_liver.sub, blind=TRUE)

pcadat_liver.sub<-plotPCA(vsd_liver.sub,intgroup="Group",returnData=T)

all(rownames(pcadat_liver.sub)==rownames(metadata_liver.sub)) #sanity check
pcadat_liver.sub$Group<-metadata_liver.sub$Group

percentVar.vsd_liver.sub<-round(100*attr(pcadat_liver.sub,"percentVar"))

pcadat_liver.sub$name<-str_split_i(pcadat_liver.sub$name,'_',1)
pcadat_liver.sub$name<-str_split_i(pcadat_liver.sub$name,'\\.',2)

pcadat_liver.sub<-subset(pcadat_liver.sub,!pcadat_liver.sub$name %in% c('1550','1548'))


pca_liver.sub<-ggplot(pcadat_liver.sub, aes(PC1, PC2, fill=Group)) +
  geom_point(size=5,pch=21,stroke=0.5)+ 
  #geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd_liver.sub[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd_liver.sub[2],"% variance")) + 
  #tat_ellipse()+
  geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = color4liver)+
  scale_color_manual(values = color4liver)+
  coord_fixed()+
  scale_x_continuous(limits = c(-18,28))+
  scale_y_continuous(limits = c(-25,18))

pdf('pca_liver.sub.pdf')
pca_liver.sub
dev.off()

###### DEGs Liver #### #need to rerun
#1  LIDPAD vs Control
res.liver.LIDPADvCtrl <- as.data.frame(results(dds_liver.sub, contrast=c('Group','LIDPAD','Control'),alpha=0.05)) 
res.liver.LIDPADvCtrl<-res.liver.LIDPADvCtrl[order(res.liver.LIDPADvCtrl$padj),]
res.liver.LIDPADvCtrl$symbol<-genesymbols[rownames(res.liver.LIDPADvCtrl)]
sig.liver.LIDPADvCtrl<-res.liver.LIDPADvCtrl[which(res.liver.LIDPADvCtrl$padj<0.05 &
                                                     (abs(res.liver.LIDPADvCtrl$log2FoldChange)>1)),]

write.csv(sig.liver.LIDPADvCtrl,"sig.liver.LIDPADvCtrl.csv")
write.csv(res.liver.LIDPADvCtrl,"res.liver.LIDPADvCtrl.csv")


#2  ctrl_semaglutide (CSG) vs control
res.liver.CSGvctrl <- as.data.frame(results(dds_liver.sub, contrast=c('Group','Control diet + Semaglutide','Control'),alpha=0.05)) 
res.liver.CSGvctrl<-res.liver.CSGvctrl[order(res.liver.CSGvctrl$padj),]
res.liver.CSGvctrl$symbol<-genesymbols[rownames(res.liver.CSGvctrl)]
sig.liver.CSGvctrl<-res.liver.CSGvctrl[which(res.liver.CSGvctrl$padj<0.05 &
                                               (abs(res.liver.CSGvctrl$log2FoldChange)>1)),]

write.csv(sig.liver.CSGvctrl,"sig.liver.CSGvctrl.csv")
write.csv(res.liver.CSGvctrl,"res.liver.CSGvctrl.csv")

#3 LP_semaglutide (LPSG) vs LIDPAD 
res.liver.LPSGvLIDPAD <- as.data.frame(results(dds_liver.sub, contrast=c('Group','LIDPAD + Semaglutide','LIDPAD'),alpha=0.05)) 
res.liver.LPSGvLIDPAD<-res.liver.LPSGvLIDPAD[order(res.liver.LPSGvLIDPAD$padj),]
res.liver.LPSGvLIDPAD$symbol<-genesymbols[rownames(res.liver.LPSGvLIDPAD)]
sig.liver.LPSGvLIDPAD<-res.liver.LPSGvLIDPAD[which(res.liver.LPSGvLIDPAD$padj<0.05 &
                                                     (abs(res.liver.LPSGvLIDPAD$log2FoldChange)>1)),]

write.csv(sig.liver.LPSGvLIDPAD,"sig.liver.LPSGvLIDPAD.csv")
write.csv(res.liver.LPSGvLIDPAD,"res.liver.LPSGvLIDPAD.csv")

#4 LP_dietreversion (LPDR) vs LIDPAD 
res.liver.LPDRvLIDPAD <- as.data.frame(results(dds_liver.sub, contrast=c('Group','LIDPAD + diet reversion','LIDPAD'),alpha=0.05)) 
res.liver.LPDRvLIDPAD<-res.liver.LPDRvLIDPAD[order(res.liver.LPDRvLIDPAD$padj),]
res.liver.LPDRvLIDPAD$symbol<-genesymbols[rownames(res.liver.LPDRvLIDPAD)]
sig.liver.LPDRvLIDPAD<-res.liver.LPDRvLIDPAD[which(res.liver.LPDRvLIDPAD$padj<0.05 &
                                                     (abs(res.liver.LPDRvLIDPAD$log2FoldChange)>1)),]

write.csv(sig.liver.LPDRvLIDPAD,"sig.liver.LPDRvLIDPAD.csv")
write.csv(res.liver.LPDRvLIDPAD,"res.liver.LPDRvLIDPAD.csv")


#5  LP_semaglutide (LPSG) vs ctrl_semaglutide (CSG) 
res.liver.LPSGvCSG <- as.data.frame(results(dds_liver.sub, contrast=c('Group','LIDPAD + Semaglutide','Control diet + Semaglutide'),alpha=0.05)) 
res.liver.LPSGvCSG<-res.liver.LPSGvCSG[order(res.liver.LPSGvCSG$padj),]
res.liver.LPSGvCSG$symbol<-genesymbols[rownames(res.liver.LPSGvCSG)]
sig.liver.LPSGvCSG<-res.liver.LPSGvCSG[which(res.liver.LPSGvCSG$padj<0.05 &
                                               (abs(res.liver.LPSGvCSG$log2FoldChange)>1)),]

write.csv(sig.liver.LPSGvCSG,"sig.liver.LPSGvCSG.csv")
write.csv(res.liver.LPSGvCSG,"res.liver.LPSGvCSG.csv")

#6 LP_semaglutide (LPSG) vs LP_dietreversion (LPDR)
res.liver.LPSGvLPDR <- as.data.frame(results(dds_liver.sub, contrast=c('Group','LIDPAD + Semaglutide','LIDPAD + diet reversion'),alpha=0.05)) 
res.liver.LPSGvLPDR<-res.liver.LPSGvLPDR[order(res.liver.LPSGvLPDR$padj),]
res.liver.LPSGvLPDR$symbol<-genesymbols[rownames(res.liver.LPSGvLPDR)]
sig.liver.LPSGvLPDR<-res.liver.LPSGvLPDR[which(res.liver.LPSGvLPDR$padj<0.05 &
                                                 (abs(res.liver.LPSGvLPDR$log2FoldChange)>1)),]

write.csv(sig.liver.LPSGvLPDR,"sig.liver.LPSGvLPDR.csv")
write.csv(res.liver.LPSGvLPDR,"res.liver.LPSGvLPDR.csv")

#plot heatmap
liver.DEG.all<-c(rownames(sig.liver.LIDPADvCtrl),rownames(sig.liver.CSGvctrl),
                 rownames(sig.liver.LPSGvLIDPAD),rownames(sig.liver.LPDRvLIDPAD),
                 rownames(sig.liver.LPSGvCSG),
                 rownames(sig.liver.LPSGvLPDR))
liver.DEG.all<-liver.DEG.all[!duplicated(liver.DEG.all)]

hm_mat.liver<-assay(vsd_liver.sub)
hm_mat.liver<-subset(hm_mat.liver,rownames(hm_mat.liver) %in% liver.DEG.all)

coldat_hm.liver<-metadata_liver.sub[order(metadata_liver.sub$Group),,drop=F]
coldat_hm.liver<-coldat_hm.liver[-c(6,8),'Group',drop=F]

hm_mat.liver<-hm_mat.liver[,rownames(coldat_hm.liver)]
hm_mat.liver<-as.data.frame(hm_mat.liver)


hm.liver<-pheatmap(hm_mat.liver,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm.liver,
         breaks = break_hm,
         annotation_colors = my_color_annotation.liver,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45,
         cutree_rows = 6,
         gaps_col = 6) 


hc<-hclust(as.dist(1-cor(t(hm_mat.liver),method = "pearson")),method = "complete")
all(hc$order==hm.liver$tree_row$order) #all TRUE
hm_cluster <- as.data.frame(cutree(tree = hc, k = 6)) #k must be the same as the number of cutree rows
table(hm_cluster)


order.row<-hm.liver$tree_row$order
hm_mat.liver <- hm_mat.liver[order.row,]

hm_cluster<-hm_cluster[rownames(hm_mat.liver),,drop=F]
colnames(hm_cluster)<-c('cluster')
all(rownames(hm_mat.liver)==rownames(hm_cluster))
hm_cluster$cluster<-factor(hm_cluster$cluster)
rowdat_hm.liver<-hm_cluster

#check heatmap
my_color_annotation.liver<-list(Group= c('Control'='black','LIDPAD'='red',
                                         'Control diet + Semaglutide'='grey50',
                                         'LIDPAD + Semaglutide'='#EBB8DD',
                                         'LIDPAD + diet reversion'='#5C62D6'),
                                cluster=c('1'='#E41A1C','2'='#377EB8','3'='#4DAF4A',
                                          '4'='#984EA3','5'="#FF7F00",'6'='#FFFF33'
                                          ))

pheatmap(hm_mat.liver,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm.liver,
         annotation_row = rowdat_hm.liver,
         breaks = break_hm,
         annotation_colors = my_color_annotation.liver,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45, gaps_col = 6) 

rowdat_hm.liver<-hm_cluster %>% 
  mutate(cluster=case_when(cluster=='5'~'1',
                           cluster=='4'~'2',
                           cluster=='1'~'3',
                           cluster=='6'~'4',
                           cluster=='2'~'5',
                           cluster=='3'~'6'
                           ))

rowdat_hm.liver<-rowdat_hm.liver %>%
  arrange(cluster)

hm_mat.liver<-hm_mat.liver[rownames(rowdat_hm.liver),]

hm.liver_reclus <-pheatmap(hm_mat.liver,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = F,
         cluster_rows = F,cluster_cols =F,
         annotation_col = coldat_hm.liver,
         annotation_row = rowdat_hm.liver,
         breaks = break_hm,
         annotation_colors = my_color_annotation.liver,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45, gaps_col = 6) 

dev.new()
pdf("hm.plot.liver.pdf")
hm.liver_reclus
dev.off()

##### liver functional analysis #####
clus1_liver<-subset(rowdat_hm.liver,rowdat_hm.liver$cluster=='1')
clus2_liver<-subset(rowdat_hm.liver,rowdat_hm.liver$cluster=='2')
clus3_liver<-subset(rowdat_hm.liver,rowdat_hm.liver$cluster=='3')
clus4_liver<-subset(rowdat_hm.liver,rowdat_hm.liver$cluster=='4')
clus5_liver<-subset(rowdat_hm.liver,rowdat_hm.liver$cluster=='5')
clus6_liver<-subset(rowdat_hm.liver,rowdat_hm.liver$cluster=='6')


expressed_genes_liver<-as.data.frame(count_liver.sub)
expressed_genes_liver[expressed_genes_liver=="0"]<-NA
expressed_genes_liver<-expressed_genes_liver[complete.cases(expressed_genes_liver),]
expressed_genes_liver<-rownames(expressed_genes_liver)
write.table(expressed_genes_liver,'expressed_genes_liver.txt')

background<-scan("expressed_genes_liver.txt",
                 quiet=TRUE,
                 what="")

Ensembl<-ViSEAGO::Ensembl2GO()
ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "mmusculus_gene_ensembl",
  Ensembl)

#Clus1 gene (liver) - biological process used for ontology
write.table(rownames(clus1_liver),'clus1_liver.txt')

clus1_liver<-scan(
  "clus1_liver.txt",
  quiet=TRUE,
  what="")

liver_clus1<-ViSEAGO::create_topGOdata(
  geneSel=clus1_liver,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus1<-topGO::runTest(
  liver_clus1,
  algorithm ="classic",
  statistic = "fisher")

liver_clus1_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("liver_clus1","classic_clus1")))

liver_clus1_res<-as.data.frame(liver_clus1_sResults@data)
write.table(liver_clus1_res,"liver_clus1_res.txt")

#Clus2 gene (liver) - biological process used for ontology
write.table(rownames(clus2_liver),'clus2_liver.txt')

clus2_liver<-scan(
  "clus2_liver.txt",
  quiet=TRUE,
  what="")

liver_clus2<-ViSEAGO::create_topGOdata(
  geneSel=clus2_liver,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus2<-topGO::runTest(
  liver_clus2,
  algorithm ="classic",
  statistic = "fisher")

liver_clus2_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("liver_clus2","classic_clus2")))

liver_clus2_res<-as.data.frame(liver_clus2_sResults@data)
write.table(liver_clus2_res,"liver_clus2_res.txt")

#Clus3 gene (liver) - biological process used for ontology
write.table(rownames(clus3_liver),'clus3_liver.txt')

clus3_liver<-scan(
  "clus3_liver.txt",
  quiet=TRUE,
  what="")

liver_clus3<-ViSEAGO::create_topGOdata(
  geneSel=clus3_liver,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus3<-topGO::runTest(
  liver_clus3,
  algorithm ="classic",
  statistic = "fisher")

liver_clus3_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("liver_clus3","classic_clus3")))

liver_clus3_res<-as.data.frame(liver_clus3_sResults@data)
write.table(liver_clus3_res,"liver_clus3_res.txt")

#Clus4 gene (liver) - biological process used for ontology
write.table(rownames(clus4_liver),'clus4_liver.txt')

clus4_liver<-scan(
  "clus4_liver.txt",
  quiet=TRUE,
  what="")

liver_clus4<-ViSEAGO::create_topGOdata(
  geneSel=clus4_liver,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus4<-topGO::runTest(
  liver_clus4,
  algorithm ="classic",
  statistic = "fisher")

liver_clus4_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("liver_clus4","classic_clus4")))

liver_clus4_res<-as.data.frame(liver_clus4_sResults@data)
write.table(liver_clus4_res,"liver_clus4_res.txt")

#Clus5 gene (liver) - biological process used for ontology
write.table(rownames(clus5_liver),'clus5_liver.txt')

clus5_liver<-scan(
  "clus5_liver.txt",
  quiet=TRUE,
  what="")

liver_clus5<-ViSEAGO::create_topGOdata(
  geneSel=clus5_liver,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus5<-topGO::runTest(
  liver_clus5,
  algorithm ="classic",
  statistic = "fisher")

liver_clus5_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("liver_clus5","classic_clus5")))

liver_clus5_res<-as.data.frame(liver_clus5_sResults@data)
write.table(liver_clus5_res,"liver_clus5_res.txt")


#clus6 gene (liver) - biological process used for ontology
write.table(rownames(clus6_liver),'clus6_liver.txt')

clus6_liver<-scan(
  "clus6_liver.txt",
  quiet=TRUE,
  what="")

liver_clus6<-ViSEAGO::create_topGOdata(
  geneSel=clus6_liver,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus6<-topGO::runTest(
  liver_clus6,
  algorithm ="classic",
  statistic = "fisher")

liver_clus6_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("liver_clus6","classic_clus6")))

liver_clus6_res<-as.data.frame(liver_clus6_sResults@data)
write.table(liver_clus6_res,"liver_clus6_res.txt")


# below are not used #
#plot top10 GO terms
#clus1
liver_clus1_res$condition.genes_frequency <- str_replace(liver_clus1_res$condition.genes_frequency, "\\(\\d+/\\d+\\)", "") #removes fraction behind
liver_clus1_res$condition.genes_frequency <- str_replace(liver_clus1_res$condition.genes_frequency, "%", "") #removes fraction behind

liver_clus1_res <- liver_clus1_res[order(as.numeric(as.character(liver_clus1_res$condition.genes_frequency)),decreasing = T),  ]
top_clus1 <- head(liver_clus1_res, 10)

#clus2
liver_clus2_res$condition.genes_frequency <- str_replace(liver_clus2_res$condition.genes_frequency, "\\(\\d+/\\d+\\)", "") #removes fraction behind
liver_clus2_res$condition.genes_frequency <- str_replace(liver_clus2_res$condition.genes_frequency, "%", "") #removes fraction behind

liver_clus2_res <- liver_clus2_res[order(as.numeric(as.character(liver_clus2_res$condition.genes_frequency)),decreasing = T),  ]
top_clus2 <- head(liver_clus2_res, 10)

#clus3
liver_clus3_res$condition.genes_frequency <- str_replace(liver_clus3_res$condition.genes_frequency, "\\(\\d+/\\d+\\)", "") #removes fraction behind
liver_clus3_res$condition.genes_frequency <- str_replace(liver_clus3_res$condition.genes_frequency, "%", "") #removes fraction behind

liver_clus3_res <- liver_clus3_res[order(as.numeric(as.character(liver_clus3_res$condition.genes_frequency)),decreasing = T),  ]
top_clus3 <- head(liver_clus3_res, 10)

#clus4
liver_clus4_res$condition.genes_frequency <- str_replace(liver_clus4_res$condition.genes_frequency, "\\(\\d+/\\d+\\)", "") #removes fraction behind
liver_clus4_res$condition.genes_frequency <- str_replace(liver_clus4_res$condition.genes_frequency, "%", "") #removes fraction behind

liver_clus4_res <- liver_clus4_res[order(as.numeric(as.character(liver_clus4_res$condition.genes_frequency)),decreasing = T),  ]
top_clus4 <- head(liver_clus4_res, 10)

#clus5
liver_clus5_res$condition.genes_frequency <- str_replace(liver_clus5_res$condition.genes_frequency, "\\(\\d+/\\d+\\)", "") #removes fraction behind
liver_clus5_res$condition.genes_frequency <- str_replace(liver_clus5_res$condition.genes_frequency, "%", "") #removes fraction behind

liver_clus5_res <- liver_clus5_res[order(as.numeric(as.character(liver_clus5_res$condition.genes_frequency)),decreasing = T),  ]
top_clus5 <- head(liver_clus5_res, 10)


trial <- ggplot(top_clus5) + 
  geom_bar(aes(x = reorder(term,condition.genes_frequency), y = as.numeric(condition.genes_frequency), fill = (`condition.-log10_pvalue`)), stat = "identity",size=2) +
  coord_cartesian() +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(size=rel(1.15)),
        axis.title = element_text(size=rel(1.15))) +
  xlab("GO Terms") +
  ylab("Gene Number") +
  ggtitle("Top 10 GO Terms (Cluster 5)") +
  theme(plot.title = element_text(hjust=0.5, 
                                  face = "bold")) +
  scale_fill_gradientn(name = "-log10(pvalue)",colours = c("#132B43","#56B1F7")) +
  theme(legend.title = element_text(size=rel(1.15),
                                    hjust=0.5, 
                                    face="bold"))

tiff('liver_clus5_GOBP.tiff',width=3500,height=1000,units='px',res=300,compression='lzw')
trial
dev.off()


##### WCGNA data #####

nmf.df<-as.data.frame(assay(vsd_liver.sub))
nmf.df$symbol<-genesymbols[rownames(nmf.df)]
nmf.df<-nmf.df %>%
  subset(!duplicated(symbol)) %>%
  subset(!symbol=="")

rownames(nmf.df)<-nmf.df$symbol
nmf.df$symbol<-NULL


nmf.df$var<-rowVars(as.matrix(nmf.df))
nmf.df <- nmf.df %>%
  arrange(desc(var)) %>%
  slice_head(n=2000) %>%
  dplyr::select(-c(var))

metadata_liver.sub<-metadata_liver.sub %>%
  mutate(newID=str_replace_all(MouseID,'_',''),
         newID=str_sub(newID, end = -6))

all(colnames(nmf.df)==rownames(metadata_liver.sub))

colnames(nmf.df)<-metadata_liver.sub$newID

write.csv(nmf.df,'liver.RNAseq.top2000.csv')
