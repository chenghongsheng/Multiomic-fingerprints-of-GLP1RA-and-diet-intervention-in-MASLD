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
library(MetaboAnalystR)
library(lipidr)
library(rstatix)

list.files()

### load files #####

coldat <- read.delim("./coldat.txt", row.names=1)
coldat$Treatment<-factor(coldat$Treatment,levels=c('Blank','QC','Ctrl', 'CS','LP', 'LPS',  'Rev'))
coldat<-coldat[order(coldat$Treatment),,drop=F]

Metabolomics <- read.delim("./Metabolomics.txt")                       
Metabolomics$No.<-NULL

rowdat<-Metabolomics[,1:5]
rownames(rowdat)<-rowdat$Name
rowdat$Name<-NULL


met_df<-Metabolomics
rownames(met_df)<-met_df$Name  
met_df$Name  <- NULL
met_df<-met_df[,-c(1:4)]
met_df<-met_df[,rownames(coldat)]

### set color #####

color4liver<-c('black','grey50','red','#EBB8DD','#5C62D6')
pie(rep(1,5),col=color4liver)


color<-c('black','grey80','red','#EBB8DD','#5C62D6')
pie(rep(1,5),col=color)


### normalized AUC #####
met_df<-as.data.frame(met_df-met_df$Blank_03)
met_df$Blank_03<-NULL

coldat.sub<-coldat[colnames(met_df),,drop=F]
coldat.sub$Treatment<-factor(coldat.sub$Treatment,levels=c('QC','Ctrl', 'CS','LP', 'LPS',  'Rev'))

met_df[met_df<0]<-0
met_df$count<-rowSums(met_df[,c(7:30)]==0)
met_df<-subset(met_df,met_df$count<6)
met_df$count<-NULL

rowdat.sub<-subset(rowdat,rownames(rowdat) %in% rownames(met_df))

met.matrix<-met_df+1

met.matrix <- log2(as.matrix(met.matrix)) #log transformation
boxplot(met.matrix,las=2) 


### PCA #####

pca<-prcomp(t(met.matrix))
pca_df<-as.data.frame(pca$x)

all(rownames(pca_df)==rownames(coldat.sub))

pca_df$Treatment<-coldat.sub$Treatment
pca_df$name<-rownames(pca_df)

pca.proportionvariances <- round((pca$sdev^2) / (sum(pca$sdev^2))*100,2)
pca.proportionvariances <- paste(colnames(pca_df),"(",paste(as.character(pca.proportionvariances),"%",")", sep=""))


ggplot(pca_df, aes(PC1, PC2, fill=Treatment)) +
  geom_point(size=5,pch=21,stroke=1.5) +
  geom_text(aes(label=name))+
  xlab(paste0(pca.proportionvariances[1]," variance")) +
  ylab(paste0(pca.proportionvariances[2]," variance")) + 
  #scale_fill_manual(values = color)+
  geom_mark_ellipse(aes(fill = Treatment))+
  #xlim(c(-35,35))+
  #ylim(c(-25,20))+
  theme(aspect.ratio = 1)+
  coord_fixed() #Sema4, Sema_11, Sema_5, Sema_15, Sema22 are potential outliers # remove QC


# remove QC & potential outliers
out<- c('QC_1','QC_2','QC_3','QC_4','QC_5','QC_6','Sema_4', 'Sema_11', 'Sema_5', 'Sema_15', 'Sema_22')

met.matrix<-as.data.frame(met.matrix ) %>%
  dplyr::select(-c(out))


coldat.sub<-subset(coldat.sub,rownames(coldat.sub) %in% colnames(met.matrix))
coldat.sub$Treatment<-factor(coldat.sub$Treatment,levels=c('Ctrl', 'CS','LP', 'LPS',  'Rev'))


pca<-prcomp(t(met.matrix))
pca_df<-as.data.frame(pca$x)

all(rownames(pca_df)==rownames(coldat.sub))

pca_df$Treatment<-coldat.sub$Treatment
pca_df$name<-rownames(pca_df)

pca.proportionvariances <- round((pca$sdev^2) / (sum(pca$sdev^2))*100,2)
pca.proportionvariances <- paste(colnames(pca_df),"(",paste(as.character(pca.proportionvariances),"%",")", sep=""))


pca<-ggplot(pca_df, aes(PC1, PC2, fill=Treatment)) +
  geom_point(size=5,pch=21,stroke=1.5) +
  #geom_text(aes(label=name))+
  xlab(paste0(pca.proportionvariances[1]," variance")) +
  ylab(paste0(pca.proportionvariances[2]," variance")) + 
  scale_fill_manual(values = color)+
  geom_mark_ellipse(aes(fill = Treatment))+
  xlim(c(-35,30))+
  ylim(c(-25,30))+
  theme(aspect.ratio = 1)+
  coord_fixed()

dev.new()
pdf('pca_metabolomics.pdf')
pca
dev.off()


### DeqMS #####

met_df.sub<-met_df[,rownames(coldat.sub)]
met_df.sub<-met_df.sub+1


pep.count.table = data.frame(count = rowMins(as.matrix(met_df.sub)),
                             row.names = rownames(met_df.sub))

design = model.matrix(~0+Treatment,coldat.sub)
con<-c('TreatmentCS-TreatmentCtrl',
       'TreatmentLP-TreatmentCtrl',
       'TreatmentLPS-TreatmentLP',
       'TreatmentLPS-TreatmentCS',
       'TreatmentRev-TreatmentLP',
       'TreatmentLPS-TreatmentRev') 

contrast <- makeContrasts(contrasts = con,levels=design)
fit1 = lmFit(met_df.sub,design = design)
fit2 = contrasts.fit(fit1,contrasts = contrast)
fit3 <- eBayes(fit2)
fit3$count = pep.count.table[rownames(fit3$coefficients),"count"]
fit4 <- spectraCounteBayes(fit3)

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

### heatmap #####

coldat_hm<-coldat.sub

DEM<-c(rownames(sig.CSvCtrl),
       rownames(sig.LPSvCS),
       rownames(sig.LPSvLP),
       rownames(sig.LPSvRev),
       rownames(sig.LPvCtrl),
       rownames(sig.RevvLP))

DEM<-DEM[!duplicated(DEM)]

hm_mat<-met_df.sub[DEM,rownames(coldat_hm)]

#heatmap construction
break_hm = seq(-2, 2,length.out=100)
hm_color<- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

my_color_annotation<-list(Treatment= c('Ctrl'='black',
                                       'CS'='grey50',
                                       'LP'='red',
                                       'LPS'='#EBB8DD',
                                       'Rev'='#5C62D6'))



hm_allDEM<-pheatmap(hm_mat,scale="row",border_color = NA,color = hm_color,
         show_rownames = T,show_colnames = T,
         cluster_rows = T,cluster_cols = F,
         annotation_col = coldat_hm,
         breaks = break_hm,
         annotation_colors = my_color_annotation,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         angle_col = 45,gaps_col = 8)

dev.new()
pdf("hm_allDEM.pdf")
hm_allDEM
dev.off()

### metaboanalystR ##### #use https://www.metaboanalyst.ca/ModuleView.xhtml 
metabolite_name <- read.delim("C:/Users/hscheng/OneDrive - Nanyang Technological University/Project/Semaglutide/Metabolomics/Metabolomic/metabolite_name2.txt",row.names = 1)
all(metabolite_name$Name==Metabolomics$Name)

# rev
res.rev.metab<-res.RevvLP
colnames(res.rev.metab)[7]<-'Name'
df<-subset(metabolite_name,metabolite_name$Name %in% rownames(res.rev.metab))

res.rev.metab<-res.rev.metab %>%
  left_join(.,df,by='Name') %>%
  arrange(sca.adj.pval)

write.csv(res.rev.metab,'res.rev.metab.csv')


sig.rev.metab <- subset(res.rev.metab,abs(res.rev.metab$logFC)>1 & res.rev.metab$sca.adj.pval<0.1)
sig.rev.metab$count <- rep(1,nrow(sig.rev.metab))

rev.class.pie<-ggplot(sig.rev.metab,aes(x='',y=count,fill=Class))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  scale_fill_manual(values=c('red','orange','purple','grey90'))+
  theme_void()

dev.new()
pdf('rev.class.pie.pdf')
rev.class.pie
dev.off()



rev.subclass.pie<-ggplot(subset(sig.rev.metab,sig.rev.metab$Class=='Lipids'),
                                aes(x='',y=count,fill=Subclass_2))+
  geom_bar(stat="identity", width=1) +
  coord_polar("y", start=0)+
  #scale_fill_manual(values=c('red','orange','purple','grey90'))+
  theme_void()

dev.new()
pdf('rev.subclass.pie.pdf')
rev.subclass.pie
dev.off()


#selected metabolites
sel<-subset(sig.rev.metab,sig.rev.metab$Class=='Lipids')$Name
#[1] "1-Oleoyl-2-hydroxy-sn-glycero-3-PE"           "(-)-Menthylacetate"                          
#[3] "Brassylic acid"                               "3-(Tetradecanoyloxy)hexadecanoic acid"       
#[5] "9-OAHSA"                                      "1-linoleoyl-sn-glycero-3-phosphoethanolamine"
#[7] "Dodecanedioicacid"                            "9-PAHPA"       

sel.met_df.sub<-subset(met_df.sub,rownames(met_df.sub) %in% sel)
sel.met_df.sub<-as.data.frame(t(sel.met_df.sub))

all(rownames(coldat.sub)==rownames(sel.met_df.sub))

sel.met_df.sub$Treatment <-coldat.sub$Treatment


  

sel.met_df.sub %>%
  pairwise_t_test(`Brassylic acid`~Treatment)
#.y.            group1 group2    n1    n2            p p.signif       p.adj p.adj.signif
#* <chr>          <chr>  <chr>  <int> <int>        <dbl> <chr>          <dbl> <chr>       
#  1 Brassylic acid Ctrl   CS         3     5 0.842        ns       1           ns          
#  2 Brassylic acid Ctrl   LP         3     3 0.000000758  ****     0.00000531  ****        
#  3 Brassylic acid CS     LP         5     3 0.000000153  ****     0.00000138  ****        
#  4 Brassylic acid Ctrl   LPS        3     4 0.000000245  ****     0.00000196  ****        
#  5 Brassylic acid CS     LPS        5     4 0.0000000381 ****     0.000000381 ****        
#  6 Brassylic acid LP     LPS        3     4 0.806        ns       1           ns          
#  7 Brassylic acid Ctrl   Rev        3     4 0.0816       ns       0.245       ns          
#  8 Brassylic acid CS     Rev        5     4 0.0335       *        0.134       ns          
#  9 Brassylic acid LP     Rev        3     4 0.00000519   ****     0.000026    ****        
# 10 Brassylic acid LPS    Rev        4     4 0.00000146   ****     0.00000873  ****


Brassylicacid.plot<-ggplot(sel.met_df.sub,aes(x=Treatment,y=`Brassylic acid`,fill=Treatment))+
  stat_summary(fun.data  ='mean_se',geom='errorbar',width=0.5)+
  stat_summary(fun  ='mean',geom='crossbar')+
  geom_dotplot(stackdir = 'center',binaxis = 'y')+
  scale_fill_manual(values=color4liver)+
  theme(legend.position="none")

dev.new()
pdf('Brassylic acid.plot.pdf')
Brassylicacid.plot
dev.off()


sel.met_df.sub %>%
  pairwise_t_test(Dodecanedioicacid~Treatment)
# A tibble: 10 x 9
#.y.               group1 group2    n1    n2          p p.signif     p.adj p.adj.signif
#* <chr>             <chr>  <chr>  <int> <int>      <dbl> <chr>        <dbl> <chr>       
#  1 Dodecanedioicacid Ctrl   CS         3     5 0.992      ns       0.992     ns          
#2 Dodecanedioicacid Ctrl   LP         3     3 0.0000125  ****     0.000113  ***         
#  3 Dodecanedioicacid CS     LP         5     3 0.00000359 ****     0.0000359 ****        
#  4 Dodecanedioicacid Ctrl   LPS        3     4 0.000159   ***      0.00112   **          
#  5 Dodecanedioicacid CS     LPS        5     4 0.000044   ****     0.000352  ***         
#  6 Dodecanedioicacid LP     LPS        3     4 0.076      ns       0.152     ns          
#7 Dodecanedioicacid Ctrl   Rev        3     4 0.00849    **       0.0339    *           
#  8 Dodecanedioicacid CS     Rev        5     4 0.00358    **       0.0179    *           
#  9 Dodecanedioicacid LP     Rev        3     4 0.00141    **       0.00848   **          
#  10 Dodecanedioicacid LPS    Rev        4     4 0.0441     *        0.132     ns  

Dodecanedioicacid.plot<-ggplot(sel.met_df.sub,aes(x=Treatment,y=Dodecanedioicacid,fill=Treatment))+
  stat_summary(fun.data  ='mean_se',geom='errorbar',width=0.5)+
  stat_summary(fun  ='mean',geom='crossbar')+
  geom_dotplot(stackdir = 'center',binaxis = 'y')+
  scale_fill_manual(values=color4liver)+
  theme(legend.position="none")

dev.new()
pdf('Dodecanedioicacid.plot.pdf')
Dodecanedioicacid.plot
dev.off()


#1-Oleoyl-2-hydroxy-sn-glycero-3-PE
sel.met_df.sub %>%
  pairwise_t_test(`1-Oleoyl-2-hydroxy-sn-glycero-3-PE`~Treatment)
#A tibble: 10 x 9
#.y.                                group1 group2    n1    n2           p p.signif      p.adj p.adj.signif
#* <chr>                              <chr>  <chr>  <int> <int>       <dbl> <chr>         <dbl> <chr>       
#  1 1-Oleoyl-2-hydroxy-sn-glycero-3-PE Ctrl   CS         3     5 0.604       ns       1          ns          
#2 1-Oleoyl-2-hydroxy-sn-glycero-3-PE Ctrl   LP         3     3 0.00000239  ****     0.0000191  ****        
#  3 1-Oleoyl-2-hydroxy-sn-glycero-3-PE CS     LP         5     3 0.000000316 ****     0.00000316 ****        
#  4 1-Oleoyl-2-hydroxy-sn-glycero-3-PE Ctrl   LPS        3     4 0.00772     **       0.0386     *           
#  5 1-Oleoyl-2-hydroxy-sn-glycero-3-PE CS     LPS        5     4 0.00105     **       0.00631    **          
#  6 1-Oleoyl-2-hydroxy-sn-glycero-3-PE LP     LPS        3     4 0.00018     ***      0.00126    **          
#  7 1-Oleoyl-2-hydroxy-sn-glycero-3-PE Ctrl   Rev        3     4 0.674       ns       1          ns          
#8 1-Oleoyl-2-hydroxy-sn-glycero-3-PE CS     Rev        5     4 0.304       ns       0.913      ns          
#9 1-Oleoyl-2-hydroxy-sn-glycero-3-PE LP     Rev        3     4 0.00000207  ****     0.0000186  ****        
#  10 1-Oleoyl-2-hydroxy-sn-glycero-3-PE LPS    Rev        4     4 0.0118      *        0.0473     *   


LysoPE_18.1.plot<-ggplot(sel.met_df.sub,aes(x=Treatment,y=`1-Oleoyl-2-hydroxy-sn-glycero-3-PE`,fill=Treatment))+
  stat_summary(fun.data  ='mean_se',geom='errorbar',width=0.5)+
  stat_summary(fun  ='mean',geom='crossbar')+
  geom_dotplot(stackdir = 'center',binaxis = 'y')+
  scale_fill_manual(values=color4liver)+
  theme(legend.position="none")

dev.new()
pdf('LysoPE_18.1.plot.pdf')
LysoPE_18.1.plot
dev.off()



sel.met_df.sub %>%
  pairwise_t_test(`1-linoleoyl-sn-glycero-3-phosphoethanolamine`~Treatment)
#A tibble: 10 x 9
#.y.                                          group1 group2    n1    n2         p p.signif    p.adj p.adj.signif
#* <chr>                                        <chr>  <chr>  <int> <int>     <dbl> <chr>       <dbl> <chr>       
#  1 1-linoleoyl-sn-glycero-3-phosphoethanolamine Ctrl   CS         3     5 0.0493    *        0.197    ns          
#  2 1-linoleoyl-sn-glycero-3-phosphoethanolamine Ctrl   LP         3     3 0.00144   **       0.0115   *           
#  3 1-linoleoyl-sn-glycero-3-phosphoethanolamine CS     LP         5     3 0.0000124 ****     0.000124 ***         
#  4 1-linoleoyl-sn-glycero-3-phosphoethanolamine Ctrl   LPS        3     4 0.535     ns       1        ns          
#  5 1-linoleoyl-sn-glycero-3-phosphoethanolamine CS     LPS        5     4 0.00835   **       0.0501   ns          
#  6 1-linoleoyl-sn-glycero-3-phosphoethanolamine LP     LPS        3     4 0.00295   **       0.0206   *           
#  7 1-linoleoyl-sn-glycero-3-phosphoethanolamine Ctrl   Rev        3     4 0.87      ns       1        ns          
#  8 1-linoleoyl-sn-glycero-3-phosphoethanolamine CS     Rev        5     4 0.0238    *        0.119    ns          
#  9 1-linoleoyl-sn-glycero-3-phosphoethanolamine LP     Rev        3     4 0.00117   **       0.0105   *           
#  10 1-linoleoyl-sn-glycero-3-phosphoethanolamine LPS    Rev        4     4 0.621     ns       1        ns    


LysoPE_18.2.plot<-ggplot(sel.met_df.sub,aes(x=Treatment,y=`1-linoleoyl-sn-glycero-3-phosphoethanolamine`,fill=Treatment))+
  stat_summary(fun.data  ='mean_se',geom='errorbar',width=0.5)+
  stat_summary(fun  ='mean',geom='crossbar')+
  geom_dotplot(stackdir = 'center',binaxis = 'y')+
  scale_fill_manual(values=color4liver)+
  theme(legend.position="none")

dev.new()
pdf('LysoPE_18.2.plot.pdf')
LysoPE_18.2.plot
dev.off()



sel.met_df.sub %>%
  pairwise_t_test(`9-OAHSA`~Treatment)
#A tibble: 10 x 9
#.y.     group1 group2    n1    n2           p p.signif      p.adj p.adj.signif
#* <chr>   <chr>  <chr>  <int> <int>       <dbl> <chr>         <dbl> <chr>       
#  1 9-OAHSA Ctrl   CS         3     5 0.987       ns       1          ns          
#  2 9-OAHSA Ctrl   LP         3     3 0.00000276  ****     0.0000193  ****        
#  3 9-OAHSA CS     LP         5     3 0.000000737 ****     0.00000664 ****        
#  4 9-OAHSA Ctrl   LPS        3     4 0.00000238  ****     0.0000191  ****        
#  5 9-OAHSA CS     LPS        5     4 0.000000509 ****     0.00000509 ****        
#  6 9-OAHSA LP     LPS        3     4 0.679       ns       1          ns          
#  7 9-OAHSA Ctrl   Rev        3     4 0.00102     **       0.0051     **          
#  8 9-OAHSA CS     Rev        5     4 0.000328    ***      0.00197    **          
#  9 9-OAHSA LP     Rev        3     4 0.00155     **       0.00619    **          
# 10 9-OAHSA LPS    Rev        4     4 0.00205     **       0.00619    **    



OAHSA.plot<-ggplot(sel.met_df.sub,aes(x=Treatment,y=`9-OAHSA`,fill=Treatment))+
  stat_summary(fun.data  ='mean_se',geom='errorbar',width=0.5)+
  stat_summary(fun  ='mean',geom='crossbar')+
  geom_dotplot(stackdir = 'center',binaxis = 'y')+
  scale_fill_manual(values=color4liver)+
  theme(legend.position="none")

dev.new()
pdf('9-OAHSA.plot.pdf')
OAHSA.plot
dev.off()


sel.met_df.sub %>%
  pairwise_t_test(`9-PAHPA`~Treatment)

# A tibble: 10 x 9
#.y.     group1 group2    n1    n2         p p.signif    p.adj p.adj.signif
# <chr>   <chr>  <chr>  <int> <int>     <dbl> <chr>       <dbl> <chr>       
#  1 9-PAHPA Ctrl   CS         3     5 0.982     ns       1        ns          
#  2 9-PAHPA Ctrl   LP         3     3 0.0599    ns       0.24     ns          
#  3 9-PAHPA CS     LP         5     3 0.0366    *        0.22     ns          
#  4 9-PAHPA Ctrl   LPS        3     4 0.0675    ns       0.24     ns          
#  5 9-PAHPA CS     LPS        5     4 0.0387    *        0.22     ns          
#  6 9-PAHPA LP     LPS        3     4 0.839     ns       1        ns          
#  7 9-PAHPA Ctrl   Rev        3     4 0.0000544 ****     0.000489 ***         
#  8 9-PAHPA CS     Rev        5     4 0.0000136 ****     0.000136 ***         
#  9 9-PAHPA LP     Rev        3     4 0.00342   **       0.0239   *           
#  10 9-PAHPA LPS    Rev        4     4 0.00126   **       0.0101   *

PAHPA.plot<-ggplot(sel.met_df.sub,aes(x=Treatment,y=`9-PAHPA`,fill=Treatment))+
  stat_summary(fun.data  ='mean_se',geom='errorbar',width=0.5)+
  stat_summary(fun  ='mean',geom='crossbar')+
  geom_dotplot(stackdir = 'center',binaxis = 'y')+
  scale_fill_manual(values=color4liver)+
  theme(legend.position="none")

dev.new()
pdf('9-PAHPA.plot.pdf')
PAHPA.plot
dev.off()




# sema
res.sema.metab<-res.LPSvLP
colnames(res.sema.metab)[7]<-'Name'
df<-subset(metabolite_name,metabolite_name$Name %in% rownames(res.sema.metab))

res.sema.metab<-res.sema.metab %>%
  left_join(.,df,by='Name') %>%
  arrange(sca.adj.pval)

write.csv(res.sema.metab,'res.sema.metab.csv')

sig.sema.metab <- subset(res.sema.metab,abs(res.sema.metab$logFC)>1 & res.sema.metab$sca.adj.pval<0.1)


### omu ##### 
omu_df<- met_df.sub
omu_df$Name<-rownames(omu_df)
df<-subset(metabolite_name,metabolite_name$Name %in% rownames(omu_df))

omu_df<-omu_df %>%
  left_join(.,df,by='Name') 


omu_df<-assign_hierarchy(count_data = omu_df, keep_unknowns = TRUE, identifier = "KEGG")

write.csv(omu_df,'metabolite_name2.txt')

# WGCNA ####


met.matrix$var<-rowVars(as.matrix(met.matrix))

summary(met.matrix$var)

met.matrix <- met.matrix %>%
  arrange(desc(var)) %>%
  #slice_head(n=2000) %>%
  dplyr::select(-c(var))

colnames(met.matrix)<-str_replace(colnames(met.matrix),'_','')

write.csv(met.matrix,'met.matrix.csv')

