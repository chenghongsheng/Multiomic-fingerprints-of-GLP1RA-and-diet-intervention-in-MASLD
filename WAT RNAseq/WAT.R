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


##### load data ######
#metadata
metadata <- read.delim("./metadata.txt", header=T, row.names = 1)

metadata_eWAT<-subset(metadata,metadata$Tissue %in% c('eWAT'))
metadata_eWAT <- metadata_eWAT[,c("Group","Tissue")]
rownames(metadata_eWAT) <- sub("-", ".", rownames(metadata_eWAT))
#metadata_eWAT <- metadata_eWAT[-c(12), ]
metadata_eWAT$Group<-factor(metadata_eWAT$Group,
                            levels=c('Control diet','Control diet + Semaglutide',
                                     'LIDPAD','LIDPAD + Semaglutide','LIDPAD + diet reversion'))

#count
count <- read.table("./Fcount.txt", row.names=1, header = T)
count<-count[,c(6:53)]

rownames(metadata) <- sub("-", ".", rownames(metadata))
count<-count[,rownames(metadata)]

count_eWAT<-count[,rownames(metadata_eWAT)]

#sanity check

all(rownames(metadata_eWAT) == colnames(count_eWAT))

#### set color ######

color4eWAT<-c('black','grey50','red','#EBB8DD','#5C62D6')
pie(rep(1,5),col=color4eWAT)

##### DESeq_eWAT #####
#conducts normalization
dds_eWAT <- DESeqDataSetFromMatrix(countData = count_eWAT,
                                   colData = metadata_eWAT,
                                   design = ~ Group)
dds_eWAT<-DESeq(dds_eWAT)
vsd_eWAT<-varianceStabilizingTransformation(dds_eWAT, blind=TRUE)
pcadat_eWAT<-plotPCA(vsd_eWAT,intgroup="Group",returnData=T)
all(rownames(pcadat_eWAT)==rownames(metadata_eWAT)) #sanity check
pcadat_eWAT$Group<-metadata_eWAT$Group

#pcadat_eWAT<-subset(pcadat_eWAT,!pcadat_eWAT$name %in% out2)


pcadat_eWAT$name<-str_split_i(pcadat_eWAT$name,'_',1)
pcadat_eWAT$name<-str_split_i(pcadat_eWAT$name,'\\.',2)

percentVar.vsd_eWAT<-round(100*attr(pcadat_eWAT,"percentVar"))

ggplot(pcadat_eWAT, aes(PC1, PC2, fill=Group)) +
  geom_point(size=5,pch=21,stroke=0.5)+ 
  #geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd_eWAT[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd_eWAT[2],"% variance")) +
  geom_text(aes(label=name))+
  #stat_ellipse()+
  # geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = color4eWAT)+
  scale_color_manual(values = color4eWAT)+
  coord_fixed()


##### remove outliers #####
out<-c(1567,
  1573,1570,1574,
       1576,1577,1580,
       1584,1585)

metadata_eWAT$name<-rownames(metadata_eWAT)
metadata_eWAT$name<-str_split_i(metadata_eWAT$name,'_',1)
metadata_eWAT$name<-str_split_i(metadata_eWAT$name,'\\.',2)

metadata_eWAT.sub<-subset(metadata_eWAT,!metadata_eWAT$name %in% out)
count_eWAT.sub <- count_eWAT[,rownames(metadata_eWAT.sub)]

all(rownames(metadata_eWAT.sub)==colnames(count_eWAT.sub))

dds_eWAT.sub <- DESeqDataSetFromMatrix(countData = count_eWAT.sub,
                                   colData = metadata_eWAT.sub,
                                   design = ~ Group)
dds_eWAT.sub<-DESeq(dds_eWAT.sub)
vsd_eWAT.sub<-varianceStabilizingTransformation(dds_eWAT.sub, blind=TRUE)
pcadat_eWAT.sub<-plotPCA(vsd_eWAT.sub,intgroup="Group",returnData=T)
all(rownames(pcadat_eWAT.sub)==rownames(metadata_eWAT.sub)) #sanity check
pcadat_eWAT.sub$Group<-metadata_eWAT.sub$Group
pcadat_eWAT.sub<-pcadat_eWAT.sub[-c(4),]

percentVar.vsd_eWAT.sub<-round(100*attr(pcadat_eWAT.sub,"percentVar"))

pca_eWAT.sub<-ggplot(pcadat_eWAT.sub, aes(PC1, PC2, fill=Group)) +
  geom_point(size=5,pch=21,stroke=0.5)+ 
  #geom_text(aes(label=name))+
  xlab(paste0("PC1: ",percentVar.vsd_eWAT.sub[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar.vsd_eWAT.sub[2],"% variance")) +
  #geom_text(aes(label=name))+
  geom_mark_ellipse(aes(fill = group,color = group))+
  #stat_ellipse()+
  # geom_mark_ellipse(aes(fill = group,color = group))+
  theme(aspect.ratio = 1)+
  scale_fill_manual(values = color4eWAT)+
  scale_color_manual(values = color4eWAT)+
  coord_fixed()+
  scale_x_continuous(limits = c(-18,20))+
  scale_y_continuous(limits = c(-23,13))

pdf('pca_eWAT.sub.pdf')
pca_eWAT.sub
dev.off()

##### DEGs eWAT #####
Mmu.dataset<-useDataset('mmusculus_gene_ensembl',mart=useMart("ensembl"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name',"gene_biotype"), 
               filters='ensembl_gene_id',
               values=rownames(count),mart=Hsa.dataset)
genesymbols <- tapply(Genemap$external_gene_name, 
                      Genemap$ensembl_gene_id, paste, collapse="; ")

#1  LIDPAD vs control
res.eWAT.LIDPADvctrl <- as.data.frame(results(dds_eWAT.sub, contrast=c('Group','LIDPAD','Control diet'),alpha=0.05)) 
res.eWAT.LIDPADvctrl<-res.eWAT.LIDPADvctrl[order(res.eWAT.LIDPADvctrl$padj),]
res.eWAT.LIDPADvctrl$symbol<-genesymbols[rownames(res.eWAT.LIDPADvctrl)]
sig.eWAT.LIDPADvctrl<-res.eWAT.LIDPADvctrl[which(res.eWAT.LIDPADvctrl$padj<0.05 &
                                                   (abs(res.eWAT.LIDPADvctrl$log2FoldChange)>1)),]

write.csv(sig.eWAT.LIDPADvctrl,"sig.eWAT.LIDPADvctrl.csv")
write.csv(res.eWAT.LIDPADvctrl,"res.eWAT.LIDPADvctrl.csv")

#2 ctrl_semaglutide (CSG) vs control 
res.eWAT.CSGvctrl <- as.data.frame(results(dds_eWAT.sub, contrast=c('Group','Control diet + Semaglutide','Control diet'),alpha=0.05)) 
res.eWAT.CSGvctrl<-res.eWAT.CSGvctrl[order(res.eWAT.CSGvctrl$padj),]
res.eWAT.CSGvctrl$symbol<-genesymbols[rownames(res.eWAT.CSGvctrl)]
sig.eWAT.CSGvctrl<-res.eWAT.CSGvctrl[which(res.eWAT.CSGvctrl$padj<0.05 &
                                             (abs(res.eWAT.CSGvctrl$log2FoldChange)>1)),]

write.csv(sig.eWAT.CSGvctrl,"sig.eWAT.CSGvctrl.csv")
write.csv(res.eWAT.CSGvctrl,"res.eWAT.CSGvctrl.csv")


#3 LP_semaglutide (LPSG) vs LIDPAD
res.eWAT.LPSGvLIDPAD <- as.data.frame(results(dds_eWAT.sub, contrast=c('Group','LIDPAD + Semaglutide','LIDPAD'),alpha=0.05)) 
res.eWAT.LPSGvLIDPAD<-res.eWAT.LPSGvLIDPAD[order(res.eWAT.LPSGvLIDPAD$padj),]
res.eWAT.LPSGvLIDPAD$symbol<-genesymbols[rownames(res.eWAT.LPSGvLIDPAD)]
sig.eWAT.LPSGvLIDPAD<-res.eWAT.LPSGvLIDPAD[which(res.eWAT.LPSGvLIDPAD$padj<0.05 &
                                                   (abs(res.eWAT.LPSGvLIDPAD$log2FoldChange)>1)),]

write.csv(sig.eWAT.LPSGvLIDPAD,"sig.eWAT.LPSGvLIDPAD.csv")
write.csv(res.eWAT.LPSGvLIDPAD,"res.eWAT.LPSGvLIDPAD.csv")

#4  LP_dietreversion (LPDR) vs LIDPAD 
res.eWAT.LPDRvLIDPAD <- as.data.frame(results(dds_eWAT.sub, contrast=c('Group','LIDPAD + diet reversion','LIDPAD'),alpha=0.05)) 
res.eWAT.LPDRvLIDPAD<-res.eWAT.LPDRvLIDPAD[order(res.eWAT.LPDRvLIDPAD$padj),]
res.eWAT.LPDRvLIDPAD$symbol<-genesymbols[rownames(res.eWAT.LPDRvLIDPAD)]
sig.eWAT.LPDRvLIDPAD<-res.eWAT.LPDRvLIDPAD[which(res.eWAT.LPDRvLIDPAD$padj<0.05 &
                                                   (abs(res.eWAT.LPDRvLIDPAD$log2FoldChange)>1)),]

write.csv(sig.eWAT.LPDRvLIDPAD,"sig.eWAT.LPDRvLIDPAD.csv")
write.csv(res.eWAT.LPDRvLIDPAD,"res.eWAT.LPDRvLIDPAD.csv")

#5  LP_semaglutide (LPSG) vs ctrl_semaglutide (CSG)
res.eWAT.LPSGvCSG <- as.data.frame(results(dds_eWAT.sub, contrast=c('Group','LIDPAD + Semaglutide','Control diet + Semaglutide'),alpha=0.05)) 
res.eWAT.LPSGvCSG<-res.eWAT.LPSGvCSG[order(res.eWAT.LPSGvCSG$padj),]
res.eWAT.LPSGvCSG$symbol<-genesymbols[rownames(res.eWAT.LPSGvCSG)]
sig.eWAT.LPSGvCSG<-res.eWAT.LPSGvCSG[which(res.eWAT.LPSGvCSG$padj<0.05 &
                                             (abs(res.eWAT.LPSGvCSG$log2FoldChange)>1)),]

write.csv(sig.eWAT.LPSGvCSG,"sig.eWAT.LPSGvCSG.csv")
write.csv(res.eWAT.LPSGvCSG,"res.eWAT.LPSGvCSG.csv")


#6 LP_semaglutide (LPSG) vs LP_dietreversion (LPDR)
res.eWAT.LPSGvLPDR <- as.data.frame(results(dds_eWAT.sub, contrast=c('Group','LIDPAD + Semaglutide','LIDPAD + diet reversion'),alpha=0.05)) 
res.eWAT.LPSGvLPDR<-res.eWAT.LPSGvLPDR[order(res.eWAT.LPSGvLPDR$padj),]
res.eWAT.LPSGvLPDR$symbol<-genesymbols[rownames(res.eWAT.LPSGvLPDR)]
sig.eWAT.LPSGvLPDR<-res.eWAT.LPSGvLPDR[which(res.eWAT.LPSGvLPDR$padj<0.05 &
                                               (abs(res.eWAT.LPSGvLPDR$log2FoldChange)>1)),]

write.csv(sig.eWAT.LPSGvLPDR,"sig.eWAT.LPSGvLPDR.csv")
write.csv(res.eWAT.LPSGvLPDR,"res.eWAT.LPSGvLPDR.csv")


#plot heatmap
eWAT.DEG.all<-c(rownames(sig.eWAT.CSGvctrl),rownames(sig.eWAT.LIDPADvctrl),
                rownames(sig.eWAT.LPDRvLIDPAD),rownames(sig.eWAT.LPSGvCSG),
                rownames(sig.eWAT.LPSGvLIDPAD),rownames(sig.eWAT.LPSGvLPDR))
eWAT.DEG.all<-eWAT.DEG.all[!duplicated(eWAT.DEG.all)]

hm_mat.eWAT<-assay(vsd_eWAT.sub)
hm_mat.eWAT<-subset(hm_mat.eWAT,rownames(hm_mat.eWAT) %in% eWAT.DEG.all)

coldat_hm.eWAT<-metadata_eWAT.sub[order(metadata_eWAT.sub$Group),c(1),drop=F]
coldat_hm.eWAT<-coldat_hm.eWAT[-3,,drop=F]

hm_mat.eWAT<-hm_mat.eWAT[,rownames(coldat_hm.eWAT)]
hm_mat.eWAT<-as.data.frame(hm_mat.eWAT)

hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
break_hm = seq(-2, 2,length.out=100)

my_color_annotation.eWAT<-list(Group= c('Control diet'='black','LIDPAD'='red',
                                         'Control diet + Semaglutide'='grey50',
                                         'LIDPAD + Semaglutide'='#EBB8DD',
                                         'LIDPAD + diet reversion'='#5C62D6'))

hm.eWAT<-pheatmap(hm_mat.eWAT,scale="row",border_color = NA,color = hm_color,
                       show_rownames = F,show_colnames = T,
                       cluster_rows = T,cluster_cols =F,
                       annotation_col = coldat_hm.eWAT,
                       breaks = break_hm,
                       annotation_colors = my_color_annotation.eWAT,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       angle_col = 45,gaps_col = 6) 



hc<-hclust(as.dist(1-cor(t(hm_mat.eWAT),method = "pearson")),method = "complete")
all(hc$order==hm.eWAT$tree_row$order) #all TRUE
hm_cluster <- as.data.frame(cutree(tree = hc, k = 2))
table(hm_cluster)
 
  
order.row<-hm.eWAT$tree_row$order
hm_mat.eWAT <- hm_mat.eWAT[order.row,]

hm_cluster<-hm_cluster[rownames(hm_mat.eWAT),,drop=F]
colnames(hm_cluster)<-c('cluster')
all(rownames(hm_mat.eWAT)==rownames(hm_cluster))

pheatmap(hm_mat.eWAT,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = T,
         cluster_rows = T,cluster_cols =F,
         annotation_col = coldat_hm.eWAT,
         annotation_row = hm_cluster,
         breaks = break_hm,
         annotation_colors = my_color_annotation.eWAT,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45,gaps_col = 6)


rowdat_hm.eWAT<-hm_cluster
rowdat_hm.eWAT$cluster<-factor(rowdat_hm.eWAT$cluster)

my_color_annotation.eWAT<-list(Group= c('Control diet'='black','LIDPAD'='red',
                                        'Control diet + Semaglutide'='grey50',
                                        'LIDPAD + Semaglutide'='#EBB8DD',
                                        'LIDPAD + diet reversion'='#5C62D6'),
                               cluster=c('1'='#E41A1C','2'='#377EB8')
                               )

hm.plot.eWAT<-pheatmap(hm_mat.eWAT,scale="row",border_color = NA,color = hm_color,
                       show_rownames = F,show_colnames = T,
                       cluster_rows = T,cluster_cols =F,
                       annotation_col = coldat_hm.eWAT,
                       annotation_row = rowdat_hm.eWAT,
                       breaks = break_hm,
                       annotation_colors = my_color_annotation.eWAT,
                       clustering_distance_rows = "correlation",
                       clustering_distance_cols = "correlation",
                       angle_col = 45,gaps_col = 6)

dev.new()
pdf('hm.plot.eWAT.pdf')
hm.plot.eWAT
dev.off()


######plot heatmap add sample size######
out2<-c(1567,
       1573,1570,1574,
       1576,1580,
       1584,1585)

hm_mat2.eWAT<-assay(vsd_eWAT)
hm_mat2.eWAT<-subset(hm_mat2.eWAT,rownames(hm_mat2.eWAT) %in% eWAT.DEG.all)

coldat2_hm.eWAT<-metadata_eWAT[order(metadata_eWAT$Group),,drop=F]
coldat2_hm.eWAT<-subset(coldat2_hm.eWAT,!coldat2_hm.eWAT$name %in% out2)
coldat2_hm.eWAT<-coldat2_hm.eWAT[-3,1,drop=F]

hm_mat2.eWAT<-hm_mat2.eWAT[,rownames(coldat2_hm.eWAT)]
hm_mat2.eWAT<-as.data.frame(hm_mat2.eWAT)
hm_mat2.eWAT <- hm_mat2.eWAT[order.row,]

hm.plot.eWAT<-pheatmap(hm_mat2.eWAT,scale="row",border_color = NA,color = hm_color,
         show_rownames = F,show_colnames = T,
         cluster_rows = F,cluster_cols =F,
         annotation_col = coldat2_hm.eWAT,
         annotation_row = rowdat_hm.eWAT,
         breaks = break_hm,
         annotation_colors = my_color_annotation.eWAT,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45,gaps_col = 6)

dev.new()
pdf('hm.plot.eWAT.pdf')
hm.plot.eWAT
dev.off()



##### eWAT functional analysis #####
clus1_eWAT<-subset(rowdat_hm.eWAT,rowdat_hm.eWAT$cluster=='1')
clus2_eWAT<-subset(rowdat_hm.eWAT,rowdat_hm.eWAT$cluster=='2')

expressed_genes_eWAT<-as.data.frame(count_eWAT.sub)
expressed_genes_eWAT[expressed_genes_eWAT=="0"]<-NA
expressed_genes_eWAT<-expressed_genes_eWAT[complete.cases(expressed_genes_eWAT),]
expressed_genes_eWAT<-rownames(expressed_genes_eWAT)
write.table(expressed_genes_eWAT,'expressed_genes_eWAT.txt')

background<-scan("expressed_genes_eWAT.txt",
                 quiet=TRUE,
                 what="")

Ensembl<-ViSEAGO::Ensembl2GO()
ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "mmusculus_gene_ensembl",
  Ensembl)

#Clus1 gene (eWAT) - biological process used for ontology
write.table(rownames(clus1_eWAT),'clus1_eWAT.txt')

clus1_eWAT<-scan(
  "clus1_eWAT.txt",
  quiet=TRUE,
  what="")

eWAT_clus1<<-ViSEAGO::create_topGOdata(
  geneSel=clus1_eWAT,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus1_eWAT<-topGO::runTest(
  eWAT_clus1,
  algorithm ="classic",
  statistic = "fisher")

eWAT_clus1_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("eWAT_clus1","classic_clus1_eWAT")))

eWAT_clus1_res<-as.data.frame(eWAT_clus1_sResults@data)
write.table(eWAT_clus1_res,"eWAT_clus1_res.txt")

#Clus2 gene (eWAT) - biological process used for ontology
write.table(rownames(clus2_eWAT),'clus2_eWAT.txt')

clus2_eWAT<-scan(
  "clus2_eWAT.txt",
  quiet=TRUE,
  what="")

eWAT_clus2<<-ViSEAGO::create_topGOdata(
  geneSel=clus2_eWAT,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_clus2_eWAT<-topGO::runTest(
  eWAT_clus2,
  algorithm ="classic",
  statistic = "fisher")

eWAT_clus2_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("eWAT_clus2","classic_clus2_eWAT")))

eWAT_clus2_res<-as.data.frame(eWAT_clus2_sResults@data)
write.table(eWAT_clus2_res,"eWAT_clus2_res.txt")


##### WCGNA data #####

nmf.df<-as.data.frame(assay(vsd_eWAT.sub))
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

metadata_eWAT.sub2 <- read.delim("./metadata.txt")

metadata_eWAT.sub2 <- metadata_eWAT.sub2 %>%
    mutate(SampleID=str_replace(SampleID,'-','.'))%>%
  subset(SampleID %in% colnames(nmf.df))


metadata_eWAT.sub2<-metadata_eWAT.sub2 %>%
  mutate(newID=str_replace_all(MouseID,'_',''),
         newID=str_sub(newID, end = -4))

all(colnames(nmf.df)==metadata_eWAT.sub2$SampleID)

colnames(nmf.df)<-metadata_eWAT.sub2$newID

write.csv(nmf.df,'eWAT.RNAseq.top2000.csv')

