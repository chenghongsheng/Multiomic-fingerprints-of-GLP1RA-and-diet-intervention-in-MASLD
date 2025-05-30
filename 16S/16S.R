setwd("your-working-directory")

library(phyloseq)
library(ggplot2)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(vegan)
library(pairwiseAdonis)
library(ALDEx2)
library(ggpubr)
library(rstatix)
library(emmeans)
library(stringr)
library(pheatmap)


################# set unique ID ########
giveUniqueID<-function(x){
  x<-x %>%
    mutate(ID = case_when(Genus=='na'~paste0('Family_',Family),
                          TRUE ~ Genus))
  x<-x %>%
    mutate(ID = case_when(ID=='Family_na'~paste0('Order_',Order),
                          TRUE ~ ID))
  x<-x %>%
    mutate(ID = case_when(ID=='Order_na'~paste0('Class_',Class),
                          TRUE ~ ID))
  
  x<-x %>%
    mutate(ID = case_when(ID=='Class_na'~paste0('Phylum_',Phylum),
                          TRUE ~ ID))
  
  x<-x %>%
    mutate(ID = case_when(ID=='Phylum_na'~paste0('Kingdom_',Kingdom),
                          TRUE ~ ID))
}

#######load data############
#create metadata
metadata <- read.delim("./metadata.txt", row.names=1)
metadata<-sample_data(metadata)

load("../sema.mice_F240R230.RData") #Rdata from DADA2 based on raw sequencing data
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, metadata,dna)

taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

ps_percent <- transform_sample_counts(ps, function(x) 100*x/sum(x))

saveRDS(ps_percent,file="./ps_percent.RDS")
saveRDS(ps,file="./ps.RDS")

ps.genus<-tax_glom(ps,taxrank = 'Genus',NArm = F)

ps.genus_percent <- transform_sample_counts(ps.genus, function(x) 100*x/sum(x))

saveRDS(ps.genus_percent,file="../ps.genus_percent.RDS")
saveRDS(ps.genus,file="../ps.genus.RDS")

bar_Genus <- plot_bar(ps.genus_percent, fill = "Genus")
relabundance<-as.data.frame(bar_Genus$data)
relabundance<-subset(relabundance,relabundance$Kingdom=="Bacteria" & !relabundance$Abundance==0)

relabundance<-relabundance %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(percent=sum(Abundance))

relabundance<-spread(relabundance,'Sample','percent')

relabundance[,c(1:6)][is.na(relabundance[,c(1:6)])]<-'na'

relabundance<-giveUniqueID(relabundance)
write.csv(relabundance,'relabundance.csv')


coldat<-as.data.frame(bar_Genus$data) %>%
  select(Sample, Week, Sema_ID, Treatment,Group)
coldat<-coldat[!duplicated(coldat$Sample),]


coldat$Group<-factor(coldat$Group,levels=c('Ctrl_8','Ctrl_12',
                                           'CS_8' ,'CS_12',
                                           'LP_8','LP_12', 
                                           'LPS_8','LPS_12', 
                                           'Rev_8','Rev_12' ))

coldat$Week<-factor(coldat$Week)
coldat<-coldat[order(coldat$Group,coldat$Week),]


relabundance<-gather(relabundance,"Sample",'percent',7:54)
relabundance<-full_join(relabundance,coldat,by='Sample')

n <- length(levels(factor(relabundance$ID)))
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(888)
col<-sample(col_vector, n,replace = T)
pie(rep(1,n), col=col)

relabundance<-arrange(relabundance,rev(Group),rev(Week),desc(percent))

df.mean <-relabundance %>%
  select(ID,percent) %>%
  group_by(ID)%>%
  mutate(percent=replace_na(percent,0)) %>%
  summarise(mean=mean(percent)) %>%
  arrange(desc(mean))

relabundance$ID<-factor(relabundance$ID,levels =df.mean$ID[!duplicated(df.mean$ID)])


Genus_bar_FMT_used <-ggplot(relabundance,aes(fill=ID,y=percent,x=Sample))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',show.legend=T,width = 0.9)+
  scale_fill_manual(values=col)+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "None",
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  xlab("")+
  ylab("Relative Abundance (%)")

dev.new()
pdf('Genus_bar_FMT_used.pdf',width = 20)
Genus_bar_FMT_used
dev.off()

#Top20
relabundance_long<-relabundance %>%
  select(Genus,Family,Order,Class,Phylum,Kingdom,Sample,percent,ID) %>%
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  spread('Sample','percent')

colnames(relabundance_long)

relabundance_long[,c(8:55)][is.na(relabundance_long[,c(8:55)])]<-0
relabundance_long$sum<-rowSums(as.matrix(relabundance_long[,c(8:55)]))
relabundance_long<-relabundance_long %>%
  arrange(desc(sum))

Top10 = relabundance_long$ID[1:10]

relabundance.top10<-subset(relabundance,relabundance$ID  %in% Top10)

Genus_bar_FMT.top10 <-ggplot(relabundance.top10,aes(fill=ID,y=percent,x=Sample))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',show.legend=T)+
  scale_fill_manual(values=col)+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "None",
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  xlab("")+
  ylab("Relative Abundance (%)")

dev.new()
pdf('Genus_bar_FMT.top10.pdf',width = 15)
Genus_bar_FMT.top10
dev.off()

### set color #####

color<-c(rep('black',2),rep('grey80',2),
         rep('red',2),rep('#EBB8DD',2),
         rep('#5C62D6',2))
pie(rep(1,10),col=color)


#######alpha diversity############
#alpha-diversity

alpha<-plot_richness(ps) #do not use Chao1 index from DADA2
diversity_index<-as.data.frame(alpha$data)

diversity_index$Group <-factor(diversity_index$Group,levels=c('Ctrl_8','Ctrl_12',
                                                     'CS_8' ,'CS_12',
                                                     'LP_8','LP_12', 
                                                     'LPS_8','LPS_12', 
                                                     'Rev_8','Rev_12' ))

diversity_index<-diversity_index %>%
  select(-c(se)) %>%
  spread(variable,value)

write.csv(diversity_index,file="diversity_index.csv")

Chao1 <- ggplot(diversity_index,aes(x=Group,y=Chao1,color=Group))+
  geom_point(size=3)+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  stat_summary(fun = mean,geom = "crossbar", width = 0.5,color='black')+
  scale_y_log10()+
  scale_color_manual(values=color)+
  ggtitle("bar is mean")+
  xlab("")+
  ylab("Chao1")

dev.new()
pdf('Chao1.pdf',width=10)
Chao1
dev.off()


Shannon <- ggplot(diversity_index,aes(x=Group,y=Shannon,color=Group))+
  geom_point(size=3)+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  stat_summary(fun = mean,geom = "crossbar", width = 0.5,color='black')+
  scale_y_log10()+
  scale_color_manual(values=color)+
  ggtitle("bar is mean")+
  xlab("")+
  ylab("Shannon")

dev.new()
pdf('Shannon.pdf',width=10)
Shannon
dev.off()



InvSimpson <- ggplot(diversity_index,aes(x=Group,y=InvSimpson,color=Group))+
  geom_point(size=3)+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        axis.line=element_line(linewidth =1,color="black"))+
  facet_grid(.~Group,scales = "free_x",space = "free_x")+
  stat_summary(fun = mean,geom = "crossbar", width = 0.5,color='black')+
  scale_y_log10()+
  scale_color_manual(values=color)+
  ggtitle("bar is mean")+
  xlab("")+
  ylab("InvSimpson")

dev.new()
pdf('InvSimpson.pdf',width=10)
InvSimpson
dev.off()


#ANOVA 

diversity.nodonor<-subset(diversity,!diversity$treatment=='donor')
diversity.nodonor$newID<-paste0(str_split_i(diversity.nodonor$ID,'-',2),
                                "_",
                                str_split_i(diversity.nodonor$ID,'-',3))
diversity.nodonor$treatment<-factor(diversity.nodonor$treatment)
diversity.nodonor$week<-factor(diversity.nodonor$week)
diversity.nodonor$diet<-factor(diversity.nodonor$diet)

diversity_index$Week<-factor(diversity_index$Week)
diversity_index$Treatment<-factor(diversity_index$Treatment,levels=c('Ctrl','CS','LP','LPS','Rev'))


Shannon_aov<-anova_test(diversity_index,dv=Shannon, wid=Sema_ID,within=Week,between = Treatment)
Shannon_aov #non-sig
#       Effect DFn DFd     F     p p<.05   ges
#1      Treatment   4  19 1.672 0.198       0.189
#2           Week   1  19 3.836 0.065       0.064
#3 Treatment:Week   4  19 2.558 0.072       0.154


Chao1_aov<-anova_test(diversity_index,dv=Chao1, wid=Sema_ID,within=Week,between = Treatment)
Chao1_aov
# Effect DFn DFd     F     p p<.05   ges
#1      Treatment   4  19 1.250 0.324       0.118
#2           Week   1  19 1.107 0.306       0.028
#3 Treatment:Week   4  19 1.386 0.276       0.126

InvSimpson_aov<-anova_test(diversity_index,dv=InvSimpson, wid=Sema_ID,within=Week,between = Treatment)
InvSimpson_aov #non-sig
# Effect DFn DFd     F     p p<.05   ges
#1      Treatment   4  19 1.321 0.298       0.138
#2           Week   1  19 1.314 0.266       0.028
#3 Treatment:Week   4  19 2.759 0.058       0.197


########### beta diversity ################
# Plot PCoA of Sample using bray-curtis dissimilarity 
# need to collapse to genus level 

ps.genus_percent@sam_data$Week<-factor(ps.genus_percent@sam_data$Week)
ps.genus_percent@sam_data$Treatment<-factor(ps.genus_percent@sam_data$Treatment,
                                            levels=c('Ctrl','CS','LP','LPS','Rev'))

ps.genus.ord <-ordinate(ps.genus_percent, "PCoA", "bray")

PCoA <- plot_ordination(ps.genus_percent, ps.genus.ord, type="samples", color="Week", title="Treatment effect")


PCoA  <- PCoA  + 
  geom_point(aes(color=Week),size=5)+
  #scale_color_manual(values=color)+
  theme(text = element_text(size=15),axis.text.x = element_text(angle=45,hjust=1))+
  facet_wrap(~Treatment,1)+
  xlab("PCoA 1 (29%)")+
  ylab("PCoA 2 (20.8%)")

dev.new()
pdf('PCoA.beta.diversity.pdf',width=10)
PCoA 
dev.off()


#Permanova using BC dissimilarity
metadata <- as(sample_data(ps.genus_percent), "data.frame")


metadata$Week<-factor(as.numeric(metadata$Week))

metadata$Sema_ID<-factor(metadata$Sema_ID)

dist<-distance(ps.genus_percent, method="bray")

permanova_ps <- adonis2(dist~ Treatment*Week , 
                        strata=metadata$Sema_ID,
                        data=metadata, 
                        permutations=9999)
permanova_ps

#Number of permutations: 9999

#adonis2(formula = dist ~ Treatment * Week, data = metadata, permutations = 9999, strata = metadata$Sema_ID)
#                Df SumOfSqs      R2      F Pr(>F)  
#Treatment       4   2.6212 0.24289 3.5203 0.0113 *
#  Week          1   0.2337 0.02165 1.2554 0.0697 .
#Treatment:Week  4   0.8629 0.07996 1.1589 0.0217 *
#  Residual     38   7.0736 0.65549                
#Total          47  10.7914 1.00000      

pairwise.adonis(dist, metadata$Treatment:metadata$Week,p.adjust.m ="BH")
#                    pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
#1     control:0 vs LIDPAD:0  1 0.10864676 1.1433537 0.08084035   0.360 0.415384615    
#18    control:4 vs LIDPAD:4  1 0.22565544 2.0816521 0.13802547   0.007 0.012600000   .
#31    control:8 vs LIDPAD:8  1 0.21886220 2.5380117 0.16334212   0.017 0.027321429   .
#40  control:12 vs LIDPAD:12  1 0.28598754 3.2892517 0.20192773   0.003 0.006136364   *
#45  control:16 vs LIDPAD:16  1 0.27419331 2.9450190 0.18469837   0.001 0.003750000   *

pairwise.adonis2(dist ~ Group,data=metadata,p.adjust.m ="BH" )

#$Ctrl_8_vs_Ctrl_12
#Df SumOfSqs     R2     F Pr(>F)
#Group     1  0.22519 0.1183 0.805  0.668
#Residual  6  1.67841 0.8817             
#Total     7  1.90361 1.0000   

#$LP_8_vs_LP_12
#Df SumOfSqs      R2      F Pr(>F)
#Group     1  0.08827 0.06524 0.4188   0.91
#Residual  6  1.26477 0.93476              
#Total     7  1.35304 1.00000  

#$CS_8_vs_CS_12
#Df SumOfSqs      R2      F Pr(>F)
#Group     1  0.17087 0.08906 0.9777  0.381
#Residual 10  1.74759 0.91094              
#Total    11  1.91846 1.00000  

#$LPS_8_vs_LPS_12
#Df SumOfSqs      R2      F Pr(>F)
#Group     1  0.24628 0.16821 1.6178  0.193
#Residual  8  1.21788 0.83179              
#Total     9  1.46416 1.00000 

# $Rev_8_vs_Rev_12
# Df SumOfSqs      R2     F Pr(>F)  
# Group     1  0.36596 0.23904 2.513  0.014 *
#  Residual  8  1.16500 0.76096               
# Total     9  1.53096 1.00000               



################# differential analysis ########


bar_Genus <- plot_bar(ps.genus, fill = "Genus")

rownames(coldat)<-coldat$Sample

count<-as.data.frame(bar_Genus$data)
count<-subset(count,count$Kingdom=="Bacteria" & !count$Abundance==0)

count<-count %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(count=sum(Abundance))

count<-spread(count,'Sample','count')

count[,c(1:6)][is.na(count[,c(1:6)])]<-'na'

count<-giveUniqueID(count)

count<-as.data.frame(count[,c(coldat$Sample,'ID')])
count[is.na(count)]<-0
rownames(count)<-count$ID
count$ID<-NULL


#interpret aldex2 results#
# rab.all is median log2 relative abundance
# diff.btw is log2 fold difference
# diff.win is log2 variance
# effect is the ratio of the two
# we.ep - Expected P value of Welch’s t test
# we.eBH - Expected Benjamini-Hochberg corrected P value of Welch’s t test
# wi.ep - Expected P value of Wilcoxon rank test
# wi.eBH - Expected Benjamini-Hochberg corrected P value of Wilcoxon test


# ctrl_12 vs ctrl_8 ####
coldata_C12vC8 <- subset(coldat, coldat$Group %in% c("Ctrl_8","Ctrl_12"))
coldata_C12vC8$Group<-factor(coldata_C12vC8$Group,levels=c('Ctrl_8','Ctrl_12'))
count_C12vC8<- count[,rownames(coldata_C12vC8)]
aldex_C12vC8 <- aldex(count_C12vC8, coldata_C12vC8$Group, mc.samples=1000, test="t", effect=TRUE,
                           include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_C12vC8<-subset(aldex_C12vC8,abs(aldex_C12vC8$effect)>0.8 & 
                          abs(aldex_C12vC8$diff.btw)>1)

write.csv(aldex_C12vC8,'aldex_C12vC8.csv' )
write.csv(sig_C12vC8,'sig_C12vC8.csv' )


# CS_12 vs CS_8 ####
coldata_CS12vCS8 <- subset(coldat, coldat$Group %in% c("CS_8","CS_12"))
coldata_CS12vCS8$Group<-factor(coldata_CS12vCS8$Group,levels=c('CS_8','CS_12'))
count_CS12vCS8<- count[,rownames(coldata_CS12vCS8)]

aldex_CS12vCS8 <- aldex(count_CS12vCS8, coldata_CS12vCS8$Group, mc.samples=1000, test="t", effect=TRUE,
                      include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_CS12vCS8<-subset(aldex_CS12vCS8,abs(aldex_CS12vCS8$effect)>0.8 & 
                     abs(aldex_CS12vCS8$diff.btw)>1)

write.csv(aldex_CS12vCS8,'aldex_CS12vCS8.csv' )
write.csv(sig_CS12vCS8,'sig_CS12vCS8.csv' )


# LP_12 vs LP_8 ####
coldata_LP12vLP8 <- subset(coldat, coldat$Group %in% c("LP_8","LP_12"))
coldata_LP12vLP8$Group<-factor(coldata_LP12vLP8$Group,levels=c('LP_8','LP_12'))
count_LP12vLP8<- count[,rownames(coldata_LP12vLP8)]

aldex_LP12vLP8 <- aldex(count_LP12vLP8, coldata_LP12vLP8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LP12vLP8<-subset(aldex_LP12vLP8,abs(aldex_LP12vLP8$effect)>0.8 & 
                       abs(aldex_LP12vLP8$diff.btw)>1)

write.csv(aldex_LP12vLP8,'aldex_LP12vLP8.csv' )
write.csv(sig_LP12vLP8,'sig_LP12vLP8.csv' )


# LPS_12 vs LPS_8 ####
coldata_LPS12vLPS8 <- subset(coldat, coldat$Group %in% c("LPS_8","LPS_12"))
coldata_LPS12vLPS8$Group<-factor(coldata_LPS12vLPS8$Group,levels=c('LPS_8','LPS_12'))
count_LPS12vLPS8<- count[,rownames(coldata_LPS12vLPS8)]

aldex_LPS12vLPS8 <- aldex(count_LPS12vLPS8, coldata_LPS12vLPS8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LPS12vLPS8<-subset(aldex_LPS12vLPS8,abs(aldex_LPS12vLPS8$effect)>0.8 & 
                       abs(aldex_LPS12vLPS8$diff.btw)>1)

write.csv(aldex_LPS12vLPS8,'aldex_LPS12vLPS8.csv' )
write.csv(sig_LPS12vLPS8,'sig_LPS12vLPS8.csv' )


# Rev_12 vs Rev_8 ####
coldata_Rev12vRev8 <- subset(coldat, coldat$Group %in% c("Rev_8","Rev_12"))
coldata_Rev12vRev8$Group<-factor(coldata_Rev12vRev8$Group,levels=c('Rev_8','Rev_12'))
count_Rev12vRev8<- count[,rownames(coldata_Rev12vRev8)]

aldex_Rev12vRev8 <- aldex(count_Rev12vRev8, coldata_Rev12vRev8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_Rev12vRev8<-subset(aldex_Rev12vRev8,abs(aldex_Rev12vRev8$effect)>0.8 & 
                       abs(aldex_Rev12vRev8$diff.btw)>1)

write.csv(aldex_Rev12vRev8,'aldex_Rev12vRev8.csv' )
write.csv(sig_Rev12vRev8,'sig_Rev12vRev8.csv' )


# CS_8 vs Ctrl_8 ####
coldata_CS8vC8 <- subset(coldat, coldat$Group %in% c("Ctrl_8","CS_8"))
coldata_CS8vC8$Group<-factor(coldata_CS8vC8$Group,levels=c('Ctrl_8','CS_8'))
count_CS8vC8<- count[,rownames(coldata_CS8vC8)]

aldex_CS8vC8 <- aldex(count_CS8vC8, coldata_CS8vC8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_CS8vC8<-subset(aldex_CS8vC8,abs(aldex_CS8vC8$effect)>0.8 & 
                       abs(aldex_CS8vC8$diff.btw)>1)

write.csv(aldex_CS8vC8,'aldex_CS8vC8.csv' )
write.csv(sig_CS8vC8,'sig_CS8vC8.csv' )


# LP_8 vs Ctrl_8 ####
coldata_LP8vC8 <- subset(coldat, coldat$Group %in% c("LP_8","Ctrl_8"))
coldata_LP8vC8$Group<-factor(coldata_LP8vC8$Group,levels=c('Ctrl_8','LP_8'))
count_LP8vC8<- count[,rownames(coldata_LP8vC8)]

aldex_LP8vC8 <- aldex(count_LP8vC8, coldata_LP8vC8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LP8vC8<-subset(aldex_LP8vC8,abs(aldex_LP8vC8$effect)>0.8 & 
                       abs(aldex_LP8vC8$diff.btw)>1)

write.csv(aldex_LP8vC8,'aldex_LP8vC8.csv' )
write.csv(sig_LP8vC8,'sig_LP8vC8.csv' )


# LPS_8 vs CS_8 ####
coldata_LPS8vCS8 <- subset(coldat, coldat$Group %in% c("CS_8","LPS_8"))
coldata_LPS8vCS8$Group<-factor(coldata_LPS8vCS8$Group,levels=c('CS_8','LPS_8'))
count_LPS8vCS8<- count[,rownames(coldata_LPS8vCS8)]

aldex_LPS8vCS8 <- aldex(count_LPS8vCS8, coldata_LPS8vCS8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LPS8vCS8<-subset(aldex_LPS8vCS8,abs(aldex_LPS8vCS8$effect)>0.8 & 
                       abs(aldex_LPS8vCS8$diff.btw)>1)

write.csv(aldex_LPS8vCS8,'aldex_LPS8vCS8.csv' )
write.csv(sig_LPS8vCS8,'sig_LPS8vCS8.csv' )



# LPS_8 vs LP_8 ####
coldata_LPS8vLP8 <- subset(coldat, coldat$Group %in% c("LP_8","LPS_12"))
coldata_LPS8vLP8$Group<-factor(coldata_LPS8vLP8$Group,levels=c('LP_8','LPS_12'))
count_LPS8vLP8<- count[,rownames(coldata_LPS8vLP8)]

aldex_LPS8vLP8 <- aldex(count_LPS8vLP8, coldata_LPS8vLP8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LPS8vLP8<-subset(aldex_LPS8vLP8,abs(aldex_LPS8vLP8$effect)>0.8 & 
                       abs(aldex_LPS8vLP8$diff.btw)>1)

write.csv(aldex_LPS8vLP8,'aldex_LPS8vLP8.csv' )
write.csv(sig_LPS8vLP8,'sig_LPS8vLP8.csv' )



# Rev_8 vs LP_8 ####
coldata_Rev8vLP8 <- subset(coldat, coldat$Group %in% c("LP_8","Rev_8"))
coldata_Rev8vLP8$Group<-factor(coldata_Rev8vLP8$Group,levels=c('LP_8','Rev_8'))
count_Rev8vLP8<- count[,rownames(coldata_Rev8vLP8)]

aldex_Rev8vLP8 <- aldex(count_Rev8vLP8, coldata_Rev8vLP8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_Rev8vLP8<-subset(aldex_Rev8vLP8,abs(aldex_Rev8vLP8$effect)>0.8 & 
                       abs(aldex_Rev8vLP8$diff.btw)>1)

write.csv(aldex_Rev8vLP8,'aldex_Rev8vLP8.csv' )
write.csv(sig_Rev8vLP8,'sig_Rev8vLP8.csv' )


# Rev_8 vs LPS_8 ####
coldata_Rev8vLPS8 <- subset(coldat, coldat$Group %in% c("LPS_8","Rev_8"))
coldata_Rev8vLPS8$Group<-factor(coldata_Rev8vLPS8$Group,levels=c('LPS_8','Rev_8'))
count_Rev8vLPS8<- count[,rownames(coldata_Rev8vLPS8)]

aldex_Rev8vLPS8 <- aldex(count_Rev8vLPS8, coldata_Rev8vLPS8$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_Rev8vLPS8<-subset(aldex_Rev8vLPS8,abs(aldex_Rev8vLPS8$effect)>0.8 & 
                       abs(aldex_Rev8vLPS8$diff.btw)>1)

write.csv(aldex_Rev8vLPS8,'aldex_Rev8vLPS8.csv' )
write.csv(sig_Rev8vLPS8,'sig_Rev8vLPS8.csv' )


# CS_12 vs Ctrl_12 ####
coldata_CS12vC12 <- subset(coldat, coldat$Group %in% c("Ctrl_12","CS_12"))
coldata_CS12vC12$Group<-factor(coldata_CS12vC12$Group,levels=c('Ctrl_12','CS_12'))
count_CS12vC12<- count[,rownames(coldata_CS12vC12)]

aldex_CS12vC12 <- aldex(count_CS12vC12, coldata_CS12vC12$Group, mc.samples=1000, test="t", effect=TRUE,
                      include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_CS12vC12<-subset(aldex_CS12vC12,abs(aldex_CS12vC12$effect)>0.8 & 
                     abs(aldex_CS12vC12$diff.btw)>1)

write.csv(aldex_CS12vC12,'aldex_CS12vC12.csv' )
write.csv(sig_CS12vC12,'sig_CS12vC12.csv' )


# LP_12 vs Ctrl_12 ####
coldata_LP12vC12 <- subset(coldat, coldat$Group %in% c("LP_12","Ctrl_12"))
coldata_LP12vC12$Group<-factor(coldata_LP12vC12$Group,levels=c('Ctrl_12','LP_12'))
count_LP12vC12<- count[,rownames(coldata_LP12vC12)]

aldex_LP12vC12 <- aldex(count_LP12vC12, coldata_LP12vC12$Group, mc.samples=1000, test="t", effect=TRUE,
                      include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LP12vC12<-subset(aldex_LP12vC12,abs(aldex_LP12vC12$effect)>0.8 & 
                     abs(aldex_LP12vC12$diff.btw)>1)

write.csv(aldex_LP12vC12,'aldex_LP12vC12.csv' )
write.csv(sig_LP12vC12,'sig_LP12vC12.csv' )


# LPS_12 vs CS_12 ####
coldata_LPS12vCS12 <- subset(coldat, coldat$Group %in% c("CS_12","LPS_12"))
coldata_LPS12vCS12$Group<-factor(coldata_LPS12vCS12$Group,levels=c('CS_12','LPS_12'))
count_LPS12vCS12<- count[,rownames(coldata_LPS12vCS12)]

aldex_LPS12vCS12 <- aldex(count_LPS12vCS12, coldata_LPS12vCS12$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LPS12vCS12<-subset(aldex_LPS12vCS12,abs(aldex_LPS12vCS12$effect)>0.8 & 
                       abs(aldex_LPS12vCS12$diff.btw)>1)

write.csv(aldex_LPS12vCS12,'aldex_LPS12vCS12.csv' )
write.csv(sig_LPS12vCS12,'sig_LPS12vCS12.csv' )



# LPS_12 vs LP_12 ####
coldata_LPS12vLP12 <- subset(coldat, coldat$Group %in% c("LP_12","LPS_12"))
coldata_LPS12vLP12$Group<-factor(coldata_LPS12vLP12$Group,levels=c('LP_12','LPS_12'))
count_LPS12vLP12<- count[,rownames(coldata_LPS12vLP12)]

aldex_LPS12vLP12 <- aldex(count_LPS12vLP12, coldata_LPS12vLP12$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_LPS12vLP12<-subset(aldex_LPS12vLP12,abs(aldex_LPS12vLP12$effect)>0.8 & 
                       abs(aldex_LPS12vLP12$diff.btw)>1)

write.csv(aldex_LPS12vLP12,'aldex_LPS12vLP12.csv' )
write.csv(sig_LPS12vLP12,'sig_LPS12vLP12.csv' )



# Rev_12 vs LP_12 ####
coldata_Rev12vLP12 <- subset(coldat, coldat$Group %in% c("LP_12","Rev_12"))
coldata_Rev12vLP12$Group<-factor(coldata_Rev12vLP12$Group,levels=c('LP_12','Rev_12'))
count_Rev12vLP12<- count[,rownames(coldata_Rev12vLP12)]

aldex_Rev12vLP12 <- aldex(count_Rev12vLP12, coldata_Rev12vLP12$Group, mc.samples=1000, test="t", effect=TRUE,
                        include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_Rev12vLP12<-subset(aldex_Rev12vLP12,abs(aldex_Rev12vLP12$effect)>0.8 & 
                       abs(aldex_Rev12vLP12$diff.btw)>1)

write.csv(aldex_Rev12vLP12,'aldex_Rev12vLP12.csv' )
write.csv(sig_Rev12vLP12,'sig_Rev12vLP12.csv' )


# Rev_12 vs LPS_12 ####
coldata_Rev12vLPS12 <- subset(coldat, coldat$Group %in% c("LPS_12","Rev_12"))
coldata_Rev12vLPS12$Group<-factor(coldata_Rev12vLPS12$Group,levels=c('LPS_12','Rev_12'))
count_Rev12vLPS12<- count[,rownames(coldata_Rev12vLPS12)]

aldex_Rev12vLPS12 <- aldex(count_Rev12vLPS12, coldata_Rev12vLPS12$Group, mc.samples=1000, test="t", effect=TRUE,
                         include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)
sig_Rev12vLPS12<-subset(aldex_Rev12vLPS12,abs(aldex_Rev12vLPS12$effect)>0.8 & 
                        abs(aldex_Rev12vLPS12$diff.btw)>1)

write.csv(aldex_Rev12vLPS12,'aldex_Rev12vLPS12.csv' )
write.csv(sig_Rev12vLPS12,'sig_Rev12vLPS12.csv' )


################# heatmap using pheatmap ##################
library(pheatmap)
library(RColorBrewer)

Dgenus <- c(rownames(sig_C12vC8),
            rownames(sig_CS12vC12),
            rownames(sig_CS12vCS8),
            rownames(sig_CS8vC8),
            rownames(sig_LP12vC12),
            rownames(sig_LP12vLP8),
            rownames(sig_LP8vC8),
            rownames(sig_LPS12vCS12),
            rownames(sig_LPS12vLP12),
            rownames(sig_LPS12vLPS8),
            rownames(sig_LPS8vCS8),
            rownames(sig_LPS8vLP8),
            rownames(sig_Rev12vLP12),
            rownames(sig_Rev12vLPS12),
            rownames(sig_Rev12vRev8),
            rownames(sig_Rev8vLP8),
            rownames(sig_Rev8vLPS8)
            )

Dgenus <- Dgenus[!duplicated(Dgenus)]

bar_Genus <- plot_bar(ps.genus_percent, fill = "Genus")

relabundance<-as.data.frame(bar_Genus$data)
relabundance<-subset(relabundance,relabundance$Kingdom=="Bacteria" & !relabundance$Abundance==0)

relabundance<-relabundance %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(percent=sum(Abundance))

relabundance<-spread(relabundance,'Sample','percent')

relabundance[,c(1:6)][is.na(relabundance[,c(1:6)])]<-'na'

relabundance<-giveUniqueID(relabundance)

relabundance<-as.data.frame(relabundance[,c(coldat$Sample,'ID')])
relabundance[is.na(relabundance)]<-0
rownames(relabundance)<-relabundance$ID
relabundance$ID<-NULL


rel.Dgenus <-subset(relabundance,rownames(relabundance) %in% Dgenus)

coldat.hm<-coldat[,'Group',drop=F]
coldat.hm<-coldat.hm[order(coldat.hm$Group),,drop=F]
rel.Dgenus<-rel.Dgenus[,rownames(coldat.hm)]

hm_color<- colorRampPalette(rev(brewer.pal(11,'RdBu')))(100)

breaks=seq(-2,2,length.out=100)


Dgenus_hm<-pheatmap(rel.Dgenus,scale="row",border_color = NA,color = hm_color,
                    show_rownames = T,show_colnames = F,
                    cluster_rows = T,cluster_cols = F,
                    #annotation_row = rowdat,
                    breaks = breaks,
                    annotation_col = coldat.hm,
                    #annotation_colors = my_color_annotation,
                    clustering_distance_rows = "correlation",
                    cellwidth = 10,cellheight = 10,
                    angle_col = 45,gaps_col = c(4,8,14,20,24,28,33,38,43))

dev.new()
pdf('Dgenus_hm.pdf',width = 12)
Dgenus_hm
dev.off()


# abundance bar-reversion ####
bar_Genus <- plot_bar(ps.genus_percent, fill = "Genus")
relabundance<-as.data.frame(bar_Genus$data)
relabundance<-subset(relabundance,relabundance$Kingdom=="Bacteria" & !relabundance$Abundance==0)

relabundance<-relabundance %>% #collapse genus
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Sample) %>%
  summarise(percent=sum(Abundance))

relabundance<-spread(relabundance,'Sample','percent')

relabundance[,c(1:6)][is.na(relabundance[,c(1:6)])]<-'na'

relabundance<-giveUniqueID(relabundance)
relabundance[,7:54][is.na(relabundance[,7:54])]<-0
relabundance$mean<-rowMeans(as.matrix(relabundance[,7:54]))
relabundance$count<-rowSums(!as.matrix(relabundance[,7:54])=='0')

relabundance <- relabundance %>%
  mutate(ID2=case_when(mean>0.1 & count>5 ~ID,
                       T~'others'))


coldat<-as.data.frame(bar_Genus$data) %>%
  select(Sample, Week, Sema_ID, Treatment,Group)
coldat<-coldat[!duplicated(coldat$Sample),]


coldat$Group<-factor(coldat$Group,levels=c('Ctrl_8','Ctrl_12',
                                           'CS_8' ,'CS_12',
                                           'LP_8','LP_12', 
                                           'LPS_8','LPS_12', 
                                           'Rev_8','Rev_12' ))

coldat$Week<-factor(coldat$Week)
coldat<-coldat[order(coldat$Group,coldat$Week),]



relabundance<-gather(relabundance,"Sample",'percent',7:54)
relabundance<-full_join(relabundance,coldat,by='Sample')
relabundance.rev<-subset(relabundance,
                         relabundance$Group %in% c('Rev_8','Rev_12') & !relabundance$percent==0)


relabundance.rev<-relabundance.rev %>%
  group_by(Genus,Family,Order,Class,Phylum,Kingdom,Group,ID,ID2) %>%
  summarise(percent=mean(percent))

relabundance.rev<-arrange(relabundance.rev,desc(percent))

n <- length(levels(factor(relabundance.rev$ID2)))-1
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
set.seed(888)
col<-c(sample(col_vector, n,replace = T),'grey50')
pie(rep(1,n+1), col=col)

relabundance.rev$ID2<-factor(relabundance.rev$ID2,levels =relabundance.rev$ID2[!duplicated(relabundance.rev$ID2)])
relabundance.rev$ID2<-factor(relabundance.rev$ID2,levels=c(levels(relabundance.rev$ID2)[-14], 'others'))
levels(relabundance.rev$ID2)



Genus.bar.rev<-ggplot(relabundance.rev,aes(fill=ID2,y=percent,x=Group))+
  geom_bar(position=position_fill(reverse = T), stat = 'identity',
           show.legend=T)+
  scale_fill_manual(values=col)+
  guides(fill=guide_legend(reverse=T))+
  theme(text = element_text(size=10),axis.text.x = element_text(angle=45,hjust=1),
        panel.background = element_blank(),
        panel.border = element_blank(),
        #legend.position = "None",
        axis.line=element_line(linewidth =1,color="black"))+
  #facet_grid(.~Group,scales = "free_x",space = "free_x")+
  xlab("")+
  ylab("Relative Abundance (%)")


dev.new()
pdf('Genus.bar.rev.pdf',width = 12)
Genus.bar.rev
dev.off()

