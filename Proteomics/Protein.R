setwd("your-working-directory")

library(DEqMS)
library(dplyr)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(limma)
library(ggforce)
library(matrixStats)
library(biomaRt)
library(ViSEAGO)

############## load data ######################################

protein_intensity1 <- read.delim("./protein_intensity1.txt")
protein_intensity2 <- read.delim("./protein_intensity2.txt")
metadata<- read.delim("./metadata.txt")

############## set color ######################################

color4liver<-c('black','grey80','red','#EBB8DD','#5C62D6')
pie(rep(1,5),col=color4liver)

############## merge data ######################################
intensity <- full_join(protein_intensity1,protein_intensity2,by="Protein.Ids")

colnames(intensity)

#create abundance matrix
abundance<-intensity[,c(6:53,58:63)]
abundance$ID<-intensity$Protein.Ids
rownames(abundance)<-abundance$ID
abundance$ID<-NULL
abundance[is.na(abundance)]<-0

#create rowdata
rowdat<-intensity[,c(1:5)]
colnames(rowdat)<-c("Protein.Group", "Protein.Ids", "Protein.Names","Genes","First.Protein.Description")
rownames(rowdat)<-rowdat$Protein.Ids
rowdat$Protein.Ids<-NULL

all(rownames(rowdat)==rownames(abundance))

#create coldata
coldat <- metadata
coldat$Treatment <-factor(coldat$Treatment,levels=c('Ctrl', 'CS', 'LP', 'LPS', 'Rev'))
coldat$ID<-str_replace(coldat$ID,"//+",'//.')

all(coldat$ID==colnames(abundance)) #TRUE

coldat$unique.sampleID<-paste0(coldat$sampleID,'_',coldat$Batch)
coldat2<-coldat %>%
  


#combine replicates 
abundance<-as.data.frame(t(abundance))
abundance$ID<-rownames(abundance)

abundance<-left_join(abundance,coldat,by='ID')

abundance$unique.sampleID<-paste0(abundance$sampleID,'_',abundance$Batch)

abundance.merged<-abundance %>%
  select(-c(ID, sampleID, Diet, Rep,Batch)) %>%
  tidyr::gather('proteinID','intensity',1:3980) %>%
  group_by(unique.sampleID,proteinID) %>%
  summarise(sum.int=sum(intensity)) %>%
  tidyr::spread(unique.sampleID,sum.int) %>%
  as.data.frame(.)

rownames(abundance.merged)<-abundance.merged$proteinID
abundance.merged$proteinID<-NULL

all(rownames(abundance.merged)==rownames(rowdat)) #sanity check

abundance.merged<-abundance.merged[rownames(rowdat),]


coldat2<-coldat[!duplicated(coldat$unique.sampleID),] %>%
  select(-c(ID,sampleID,Diet,Rep)) %>%
  arrange(Treatment)

rownames(coldat2)<-coldat2$unique.sampleID

all(colnames(abundance.merged)==rownames(coldat2)) #sanity check

abundance.merged<-abundance.merged[,rownames(coldat2)]


############## remove low expression genes ######################################
abundance.merged$count<-rowSums(abundance.merged==0)
abundance.merged<-subset(abundance.merged,abundance.merged$count<14)
abundance.merged$count<-NULL

#protein.matrix <- log2(as.matrix(abundance)) #log transformation

#boxplot(protein.matrix,las=2) 
#weird results for sema 4, sema 8 and sema 15 
#suspect haemoglobin contamination, perform corection using haemoglobin expression

############## correction ######################################

hemo<-rownames(rowdat[grep(rowdat$First.Protein.Description,pattern='Hemoglobin'),])
hemo.exp<-subset(abundance.merged,rownames(abundance.merged) %in% hemo)
hemo.exp<-as.data.frame(t(hemo.exp))
hemo.exp$mean<-rowMeans(as.matrix(hemo.exp))
all(rownames(hemo.exp)==coldat2$unique.sampleID)

mm<-model.matrix(~Treatment,coldat2)

abundance.merged.bhc<-removeBatchEffect(abundance.merged,batch=coldat2$Batch,covariates = hemo.exp$mean,design=mm)
abundance.merged.bhc<-abundance.merged.bhc+abs(min(abundance.merged.bhc))+1

protein.matrix.bhc <- log2(as.matrix(abundance.merged.bhc)) #log transformation

boxplot(protein.matrix.bhc,las=2)

############## PCA ######################################
#PCA
pca<-prcomp(t(protein.matrix.bhc))
pca_df<-as.data.frame(pca$x)

pca_df$Treatment<-coldat2$Treatment
pca_df$batch<-coldat2$Batch
pca_df$ID<-rownames(pca_df)

pca.proportionvariances <- round((pca$sdev^2) / (sum(pca$sdev^2))*100,2)
pca.proportionvariances <- paste(colnames(pca_df),"(",paste(as.character(pca.proportionvariances),"%",")", sep=""))


ggplot(pca_df, aes(PC1, PC2, fill=Treatment)) +
  geom_point(size=5,pch=21,stroke=1.5) +
  geom_text(aes(label=ID))+
  xlab(paste0(pca.proportionvariances[1],"% variance")) +
  ylab(paste0(pca.proportionvariances[2],"% variance")) + 
  scale_fill_manual(values = color4liver)+
  #geom_mark_ellipse(aes(fill = Treatment))+
  #xlim(c(-35,35))+
  #ylim(c(-25,20))+
  theme(aspect.ratio = 1)+
  coord_fixed() #Sema1_1 is an outlier

#remove other outliers
out<-c('Sema1_1','Sema23_1','Sema14_1')

coldat2.sel<-subset(coldat2, !rownames(coldat2) %in% out)
protein.matrix.bhc.sel<-protein.matrix.bhc[,rownames(coldat2.sel)]

pca<-prcomp(t(protein.matrix.bhc.sel))
pca_df<-as.data.frame(pca$x)

pca_df$Treatment<-coldat2.sel$Treatment
pca_df$batch<-coldat2.sel$Batch
pca_df$ID<-rownames(pca_df)

pca.proportionvariances <- round((pca$sdev^2) / (sum(pca$sdev^2))*100,2)
pca.proportionvariances <- paste(colnames(pca_df),"(",paste(as.character(pca.proportionvariances),"%",")", sep=""))


pca.proteomic<-ggplot(pca_df, aes(PC1, PC2, fill=Treatment)) +
  geom_point(size=5,pch=21,stroke=1.5) +
  geom_text(aes(label=ID))+
  xlab(paste0(pca.proportionvariances[1]," variance")) +
  ylab(paste0(pca.proportionvariances[2]," variance")) + 
  scale_fill_manual(values = color4liver)+
  geom_mark_ellipse(aes(fill = Treatment))+
  xlim(c(-2,2))+
  ylim(c(-1.25,0.5))+
  theme(aspect.ratio = 1)+
  coord_fixed() 

pdf('pca.proteomic.pdf')
pca.proteomic
dev.off()


################differential analysis########################
hemo.exp.sub<-subset(hemo.exp,rownames(hemo.exp) %in% rownames(coldat2.sel))
all(rownames(hemo.exp.sub)==rownames(coldat2.sel))
coldat2.sel$hemo<-hemo.exp.sub$mean

abundance.merged.bhc.sel<-as.data.frame(abundance.merged.bhc[,rownames(coldat2.sel)])
all(colnames(abundance.merged.bhc.sel)==rownames(coldat2.sel))

pep.count.table = data.frame(count = rowMins(as.matrix(abundance.merged.bhc.sel)),
                              row.names = rownames(abundance.merged.bhc.sel))

design = model.matrix(~0+Treatment,coldat2.sel)
con<-c('TreatmentCS-TreatmentCtrl',
       'TreatmentLP-TreatmentCtrl',
       'TreatmentLPS-TreatmentLP',
       'TreatmentLPS-TreatmentCS',
       'TreatmentRev-TreatmentLP',
       'TreatmentLPS-TreatmentRev') 

contrast <- makeContrasts(contrasts = con,levels=design)
fit1 = lmFit(abundance.merged.bhc.sel,design = design)
fit2 = contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)
fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
fit4 = spectraCounteBayes(fit3)

# Visualize the fit curve
VarianceBoxplot(fit4, n=30,
                xlab="peptide count")
VarianceScatterplot(fit4)

#if you are not sure which coef_col refers to the specific contrast,type
head(fit4$coefficients)

###### CS vs Ctrl ####
res.CSvCtrl = outputResult(fit4,coef_col = 1)
res.CSvCtrl$symbol<-rowdat[rownames(res.CSvCtrl),]$Genes
sig.CSvCtrl<-res.CSvCtrl[which(res.CSvCtrl$sca.adj.pval<0.1 & abs(res.CSvCtrl$logFC)>1),]
write.csv(res.CSvCtrl,file='res.CSvCtrl.csv')
write.csv(sig.CSvCtrl,file='sig.CSvCtrl.csv')

# LP vs Ctrl ####
res.LPvCtrl = outputResult(fit4,coef_col = 2)
res.LPvCtrl$symbol<-rowdat[rownames(res.LPvCtrl),]$Genes
sig.LPvCtrl<-res.LPvCtrl[which(res.LPvCtrl$sca.adj.pval<0.1 & abs(res.LPvCtrl$logFC)>1),]
write.csv(res.LPvCtrl,file='res.LPvCtrl.csv')
write.csv(sig.LPvCtrl,file='sig.LPvCtrl.csv')


# LPS vs LP ####
res.LPSvLP = outputResult(fit4,coef_col = 3)
res.LPSvLP$symbol<-rowdat[rownames(res.LPSvLP),]$Genes
sig.LPSvLP<-res.LPSvLP[which(res.LPSvLP$sca.adj.pval<0.1 & abs(res.LPSvLP$logFC)>1),]
write.csv(res.LPSvLP,file='res.LPSvLP.csv')
write.csv(sig.LPSvLP,file='sig.LPSvLP.csv')

# LPS vs CS ####
res.LPSvCS = outputResult(fit4,coef_col = 4)
res.LPSvCS$symbol<-rowdat[rownames(res.LPSvCS),]$Genes
sig.LPSvCS<-res.LPSvCS[which(res.LPSvCS$sca.adj.pval<0.1 & abs(res.LPSvCS$logFC)>1),]
write.csv(res.LPSvCS,file='res.LPSvCS.csv')
write.csv(sig.LPSvCS,file='sig.LPSvCS.csv')

# Rev vs LP ####
res.RevvLP = outputResult(fit4,coef_col = 5)
res.RevvLP$symbol<-rowdat[rownames(res.RevvLP),]$Genes
sig.RevvLP<-res.RevvLP[which(res.RevvLP$sca.adj.pval<0.1 & abs(res.RevvLP$logFC)>1),]
write.csv(res.RevvLP,file='res.RevvLP.csv')
write.csv(sig.RevvLP,file='sig.RevvLP.csv')

# LPS vs Rev ####
res.LPSvRev = outputResult(fit4,coef_col = 6)
res.LPSvRev$symbol<-rowdat[rownames(res.LPSvRev),]$Genes
sig.LPSvRev<-res.LPSvRev[which(res.LPSvRev$sca.adj.pval<0.1 & abs(res.LPSvRev$logFC)>1),]
write.csv(res.LPSvRev,file='res.LPSvRev.csv')
write.csv(sig.LPSvRev,file='sig.LPSvRev.csv')


######## volcano plot ######

library(EnhancedVolcano)

volp_df<-subset(res.RepvKO,!rownames(res.RepvKO) %in% rownames(sig.RepvKO))
volp_df<-as.data.frame(rbind(volp_df,sig.RepvKO))

Volcano_RepvKO<-EnhancedVolcano(volp_df, #res from DESEQ2
                                lab = NA,x="logFC",y="sca.P.Value",
                                title = "Rep vs KO", pCutoff = 0.05, FCcutoff = 1,
                                pointSize = 3.0,labSize = 6.0,
                                col=c('grey50', 'grey50', 'grey50', 'red3'),
                                legendPosition = "none",
                                ylab = bquote(~-Log[10]~italic(pval)),
                                colAlpha = 0.5)

tiff('Volcano_RepvKO.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
Volcano_RepvKO
dev.off()


volp_df<-subset(res.RepvWT,!rownames(res.RepvWT) %in% rownames(sig.RepvWT))
volp_df<-as.data.frame(rbind(volp_df,sig.RepvWT))

Volcano_RepvWT<-EnhancedVolcano(volp_df, #res from DESEQ2
                                lab = NA,x="logFC",y="sca.P.Value",
                                title = "Rep vs WT", pCutoff = 0.05, FCcutoff = 1,
                                pointSize = 3.0,labSize = 6.0,
                                col=c('grey50', 'grey50', 'grey50', 'red3'),
                                legendPosition = "none",
                                ylab = bquote(~-Log[10]~italic(pval)),
                                colAlpha = 0.5)

tiff('Volcano_RepvWT.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
Volcano_RepvWT
dev.off()


volp_df<-subset(res.ARKOvWT,!rownames(res.ARKOvWT) %in% rownames(sig.ARKOvWT))
volp_df<-as.data.frame(rbind(volp_df,sig.ARKOvWT))

Volcano_ARKOvWT<-EnhancedVolcano(volp_df, #res from DESEQ2
                                 lab = NA,x="logFC",y="sca.P.Value",
                                 title = "ARKO vs WT", pCutoff = 0.05, FCcutoff = 1,
                                 pointSize = 3.0,labSize = 6.0,
                                 col=c('grey50', 'grey50', 'grey50', 'red3'),
                                 legendPosition = "none",
                                 ylab = bquote(~-Log[10]~italic(pval)),
                                 colAlpha = 0.5)

tiff('Volcano_ARKOvWT.tiff',width=2000,height=2000,units='px',res=300,compression='lzw')
Volcano_ARKOvWT
dev.off()



######### heatmap ###########

coldat_hm<-coldat2.sel[,1,drop=F]

DEP<-c(rownames(sig.CSvCtrl),
       rownames(sig.LPSvCS),
       rownames(sig.LPSvLP),
       rownames(sig.LPSvRev),
       rownames(sig.LPvCtrl),
       rownames(sig.RevvLP))

DEP<-DEP[!duplicated(DEP)]

hm_mat<-abundance.merged.bhc.sel[DEP,rownames(coldat_hm)]


#heatmap construction
break_hm = seq(-3, 3,length.out=100)
hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

my_color_annotation<-list(Treatment= c('Ctrl'='black',
                                         'CS'='grey50',
                                         'LP'='red',
                                         'LPS'='#EBB8DD',
                                         'Rev'='#5C62D6'))


pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
                    show_rownames = T,show_colnames = T,
                    cluster_rows = T,cluster_cols = F,
                    annotation_col = coldat_hm,
                    breaks = break_hm,
                    annotation_colors = my_color_annotation,
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "correlation",
                    angle_col = 45,cutree_rows = 5)


#get heatmap clusters (4 clusters)
hc<-hclust(as.dist(1-cor(t(hm_mat),method = "pearson")),method = "complete")
all(hc$order==hm_allDEP$tree_row$order) #all TRUE
hm_cluster <- as.data.frame(cutree(tree = hc, k = 5))
table(hm_cluster)
#1   2   3   4   5 
#36  17  70 172  68  

plot(x = hc, labels =  row.names(hc), cex = 0.5)
rect.hclust(tree = hc, k =5, which = 1:5, border = 1:5, cluster = hm_cluster$`cutree(tree = hc, k = 5)`)

#heatmap reconstruction
rowdat_hm<-hm_cluster
rm(hm_cluster)
colnames(rowdat_hm)<-c('cluster')
rowdat_hm$cluster<-factor(rowdat_hm$cluster)

rowdat_hm<-rowdat_hm %>% 
  mutate(cluster=case_when(cluster=='5'~'1',
                           cluster=='1'~'2',
                           cluster=='4'~'3',
                           cluster=='2'~'4',
                           cluster=='3'~'5'
  ))


brewer.pal(5,"Set2")
pie(rep(1,5),col =brewer.pal(5,"Set2") )

my_color_annotation<-list(Treatment = c('Ctrl'='black',
                                       'CS'='grey80',
                                       'LP'='red',
                                       'LPS'='#EBB8DD',
                                       'Rev'='#5C62D6'),
                          cluster=c('1'="#FC8D62",'2'="#66C2A5",
                                         '3'="#E78AC3",'4'="#8DA0CB",
                                         '5'='#A6D854'))

hm_allDEP<-pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
                    show_rownames = F,show_colnames = F,
                    cluster_rows = T,cluster_cols = F,
                    annotation_col = coldat_hm,
                    annotation_row = rowdat_hm,
                    breaks = break_hm,
                    annotation_colors = my_color_annotation,
                    clustering_distance_rows = "correlation",
                    clustering_distance_cols = "correlation",
                    angle_col = 45,gaps_col = 9)

dev.new()
pdf("hm_allDEP.pdf")
hm_allDEP
dev.off()

##############Enrichment analysis############

#obtain genes from each cluster for functional enrichment
Mmu.dataset<-useDataset('mmusculus_gene_ensembl',mart=useMart("ensembl"))
Genemap<-getBM(attributes = c('ensembl_gene_id','external_gene_name',
                              'uniprotswissprot','uniprot_gn_id'), 
               filters='uniprot_gn_id',
               values=rownames(rowdat),mart=Mmu.dataset)
genesymbols <- tapply(Genemap$ensembl_gene_id, 
                      Genemap$uniprot_gn_id, paste, collapse="; ")


genesymbols2 <- tapply(Genemap$external_gene_name, 
                      Genemap$uniprot_gn_id, paste, collapse="; ")
#clus1
Clus1_orange<-subset(rowdat_hm,rowdat_hm$cluster=="1")
Clus1_orange$ensembl<-genesymbols[rownames(Clus1_orange)]
Clus1_orange<-unlist(str_split(Clus1_orange$ensembl,pattern=';'))
Clus1_orange<-na.omit(Clus1_orange[!duplicated(Clus1_orange)])

#clus2
Clus2_green<-subset(rowdat_hm,rowdat_hm$cluster=="2")
Clus2_green$ensembl<-genesymbols[rownames(Clus2_green)]
Clus2_green<-unlist(str_split(Clus2_green$ensembl,pattern=';'))
Clus2_green<-na.omit(Clus2_green[!duplicated(Clus2_green)])

#clus3
Clus3_pink<-subset(rowdat_hm,rowdat_hm$cluster=="3")
Clus3_pink$ensembl<-genesymbols[rownames(Clus3_pink)]
Clus3_pink<-unlist(str_split(Clus3_pink$ensembl,pattern=';'))
Clus3_pink<-na.omit(Clus3_pink[!duplicated(Clus3_pink)])

#clsu4
Clus4_blue<-subset(rowdat_hm,rowdat_hm$cluster=="4")
Clus4_blue$ensembl<-genesymbols[rownames(Clus4_blue)]
Clus4_blue<-unlist(str_split(Clus4_blue$ensembl,pattern=';'))
Clus4_blue<-na.omit(Clus4_blue[!duplicated(Clus4_blue)])

#clus5
Clus5_lightgreen<-subset(rowdat_hm,rowdat_hm$cluster=="5")
Clus5_lightgreen$ensembl<-genesymbols[rownames(Clus5_lightgreen)]
Clus5_lightgreen<-unlist(str_split(Clus5_lightgreen$ensembl,pattern=';'))
Clus5_lightgreen<-na.omit(Clus5_lightgreen[!duplicated(Clus5_lightgreen)])



#background gene
expressed_genes<-rownames(abundance.merged.bhc.sel)
expressed_genes<-genesymbols[expressed_genes]
expressed_genes<-unlist(str_split(expressed_genes,pattern=';'))
expressed_genes<-na.omit(expressed_genes[!duplicated(expressed_genes)])

write.table(expressed_genes,'expressed_genes.txt')

background<-scan("expressed_genes.txt",
                 quiet=TRUE,
                 what="")

Ensembl<-ViSEAGO::Ensembl2GO()
ViSEAGO::available_organisms(Ensembl)

myGENE2GO<-ViSEAGO::annotate(
  "mmusculus_gene_ensembl",
  Ensembl)


#Clus1_orange
write.table(Clus1_orange,'sig_Clus1_orange.txt')

sig_Clus1_orange<-scan(
  "sig_Clus1_orange.txt",
  quiet=TRUE,
  what="")

BP_Clus1_orange<-ViSEAGO::create_topGOdata(
  geneSel=sig_Clus1_orange,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_Clus1_orange<-topGO::runTest(
  BP_Clus1_orange,
  algorithm ="classic",
  statistic = "fisher")

BP_Clus1_orange_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_Clus1_orange","classic_Clus1_orange")))

BP_Clus1_orange_res<-as.data.frame(BP_Clus1_orange_sResults@data)
write.table(BP_Clus1_orange_res,"BP_Clus1_orange_res.txt")



#Clus2_green
write.table(Clus2_green,'sig_Clus2_green.txt')

sig_Clus2_green<-scan(
  "sig_Clus2_green.txt",
  quiet=TRUE,
  what="")

BP_Clus2_green<-ViSEAGO::create_topGOdata(
  geneSel=sig_Clus2_green,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_Clus2_green<-topGO::runTest(
  BP_Clus2_green,
  algorithm ="classic",
  statistic = "fisher")

BP_Clus2_green_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_Clus2_green","classic_Clus2_green")))

BP_Clus2_green_res<-as.data.frame(BP_Clus2_green_sResults@data)
write.table(BP_Clus2_green_res,"BP_Clus2_green_res.txt")


#Clus3_pink
write.table(Clus3_pink,'sig_Clus3_pink.txt')

sig_Clus3_pink<-scan(
  "sig_Clus3_pink.txt",
  quiet=TRUE,
  what="")

BP_Clus3_pink<-ViSEAGO::create_topGOdata(
  geneSel=sig_Clus3_pink,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_Clus3_pink<-topGO::runTest(
  BP_Clus3_pink,
  algorithm ="classic",
  statistic = "fisher")

BP_Clus3_pink_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_Clus3_pink","classic_Clus3_pink")))

BP_Clus3_pink_res<-as.data.frame(BP_Clus3_pink_sResults@data)
write.table(BP_Clus3_pink_res,"BP_Clus3_pink_res.txt")


#Clus4_blue
write.table(Clus4_blue,'sig_Clus4_blue.txt')

sig_Clus4_blue<-scan(
  "sig_Clus4_blue.txt",
  quiet=TRUE,
  what="")

BP_Clus4_blue<-ViSEAGO::create_topGOdata(
  geneSel=sig_Clus4_blue,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_Clus4_blue<-topGO::runTest(
  BP_Clus4_blue,
  algorithm ="classic",
  statistic = "fisher")

BP_Clus4_blue_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_Clus4_blue","classic_Clus4_blue")))

BP_Clus4_blue_res<-as.data.frame(BP_Clus4_blue_sResults@data)
write.table(BP_Clus4_blue_res,"BP_Clus4_blue_res.txt")


#Clus5_lightgreen
write.table(Clus5_lightgreen,'sig_Clus5_lightgreen.txt')

sig_Clus5_lightgreen<-scan(
  "sig_Clus5_lightgreen.txt",
  quiet=TRUE,
  what="")

BP_Clus5_lightgreen<-ViSEAGO::create_topGOdata(
  geneSel=sig_Clus5_lightgreen,
  allGenes=background,
  gene2GO=myGENE2GO, 
  ont="BP",
  nodeSize=5)

classic_Clus5_lightgreen<-topGO::runTest(
  BP_Clus5_lightgreen,
  algorithm ="classic",
  statistic = "fisher")

BP_Clus5_lightgreen_sResults<-ViSEAGO::merge_enrich_terms(
  Input=list(
    condition=c("BP_Clus5_lightgreen","classic_Clus5_lightgreen")))

BP_Clus5_lightgreen_res<-as.data.frame(BP_Clus5_lightgreen_sResults@data)
write.table(BP_Clus5_lightgreen_res,"BP_Clus5_lightgreen_res.txt")
