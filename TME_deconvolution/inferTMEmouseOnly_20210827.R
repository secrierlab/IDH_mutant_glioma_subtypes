########################
##### This script defines signatures and calculates the TME composition in mouse samples.

library(pheatmap)
library(reshape)
library(GSVA)
library(ggpubr)
library(gdata)
library(erer)

load("immuneListALL.BRAIN.update01042021.RData")

##### Read lineage signatures:
load("lineageList.plusPaoloOPC.20210719.RData")

# Read E2F targets:
e2f <- read.table("E2Ftargets_MSigDB.txt", header=FALSE, stringsAsFactors = FALSE)$V1

lineageList.noe2f <- lineageList
for (i in 1:length(lineageList.noe2f)) {
  lineageList.noe2f[[i]] <- setdiff(lineageList.noe2f[[i]],e2f)
}

lineageList <- lineageList.noe2f

##################################
######### Mouse data
####################################

### First, read in TPM counts:
expr <- read.csv("../rnaseq_latest/Adi_IDH_r123.tpm.csv")

### Also sample info:
info <- read.xls("../20210203_Adi_RNAseq_Master_modi.xlsx")
toremove <- info[which((info$Remove == "Yes")&(!grepl("Adi_9",info$File_name))&(!grepl("Adi_10",info$File_name))),]

load("mouse.mat.expr.log.keepAllC9.RData")
load("mouse.mat.expr.log.20210301.RData")


##############################
#### Lineage signatures
##############################

### First, filter out genes that are not highly enough expressed across samples:
keep <- info[which((info$Remove == "Yes")&(!grepl("Adi_9",info$File_name))&(!grepl("Adi_10",info$File_name))&(info$Condition %in% c("C1","C2","C1A","C1B","C4"))),
                   c("File_name","Condition")]
lineagemarkers <- unlist(lineageList)
mat.expre.log.keep <- mat.expr.log[which(rownames(mat.expr.log) %in% lineagemarkers),]
                                   #keep$File_name]

# Check expression of every marker and remove the ones very lowly expressed:
markerstorem <- NULL
for (i in 1:nrow(mat.expre.log.keep)) {
  if (max(mat.expre.log.keep[i,]<2)) {
    markerstorem <- c(markerstorem,rownames(mat.expre.log.keep)[i])
  }
}

# Also remove markers with lowest variance:
# vars <- sort(apply(mat.expre.log.keep, 1, function(x) var(x)))
# q<-quantile(vars)
# # Remove markers with variance < 1.25 (25% quantile):
# markerstorem2 <- names(vars[vars<q[1]])

#allmarkerstorem <- unique(c(markerstorem,markerstorem2))

### Keeping only markers that passed the cut-off:
# for (i in 1:(length(lineageList)-3)) {
#   lineageList[[i]] <- setdiff(lineageList[[i]],markerstorem)
# }
save(lineageList, file="lineageListFiltered.20210827.RData")

write.list(lineageList, "lineageSignatures.noE2Ftargets.csv")


library(fgsea)

df.fgsea.mouse <-NULL
for (s in colnames(mat.expr.log)) {
  m <- mat.expr.log[,s]
  names(m) <- rownames(mat.expr.log)
  fgseaRes <- fgsea(pathways = lineageList, 
                    stats    = m,
                    minSize  = 2,
                    maxSize  = 300)
  df.fgsea <- data.frame(fgseaRes)
  df.fgsea$Sample <- s
  df.fgsea.mouse <- rbind(df.fgsea.mouse, df.fgsea)
}
# merge with C1A/B info:
df.fgsea.mouse.plusCondition <- merge(df.fgsea.mouse,
                                      info[,c("File_name","Condition")],
                                      by.x="Sample", by.y="File_name",
                                      all.x=FALSE, all.y=FALSE)

df.fgsea.mouse.selected <- df.fgsea.mouse.plusCondition[which(df.fgsea.mouse.plusCondition$Condition
                                                            %in% c("C1","C1A","C1B","C2","C4","C9")),]
df.fgsea.mouse.selected$star <- sapply(df.fgsea.mouse.selected$padj,
                                       function(x) ifelse(x<0.0001,"***",ifelse(x<0.001,"**",
                                                                                ifelse(x<0.05,"*",""))))
df.fgsea.mouse.selected$SampleID <- apply(df.fgsea.mouse.selected[,c("Sample","Condition")],
                                          1, function(x)
                                          paste(c(strsplit(x[1],"_")[[1]][1:2],x[2]),collapse="_"))
                                          

library(ggplot2)
library(viridis)
library(hrbrthemes)
library(reshape2)
library(forcats)

df.fg <- df.fgsea.mouse.selected[order(df.fgsea.mouse.selected$Condition),]
df.fg$SampleID <- factor(df.fg$SampleID,
                         levels=unique(df.fg$SampleID))

pdf("plots.lineage27082021///lineage_signatures.fGSEA.pdf",w=12,h=3)
ggplot(df.fg, aes(SampleID, pathway,fill= NES)) + 
  geom_tile()+
  geom_text(aes(label=star), color="black", size=3) +
  scale_fill_distiller(palette = "RdBu")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  # theme(axis.title.x=element_blank(),
  #       axis.text.x=element_blank(),
  #       axis.ticks.x=element_blank())
dev.off()

## Also GSVA scores:
nada.mouse <- gsva(as.matrix(mat.expr.log),lineageList, 
                         method="gsva", kcdf="Gaussian",
                         mx.diff=TRUE, abs.ranking=FALSE)
df.scores.nada <- data.frame(t(nada.mouse))
df.scores.nada$Sample <-rownames(df.scores.nada)
df.scores.nada.gsva <- merge(df.scores.nada, info[,c("File_name","Condition")],
                             by.x="Sample", by.y="File_name",
                             all.x=FALSE, all.y=FALSE)
df.scores.nada.gsva <- df.scores.nada.gsva[which(!is.na(df.scores.nada.gsva$Condition)),]
df.gsva.nada <- melt(df.scores.nada.gsva)
df.gsva.nada.selected <- df.gsva.nada[which(df.gsva.nada$Condition 
                                            %in% c("C1","C1A","C1B","C2","C4","C9")),]
df.gsva.nada.selected$SampleID <- apply(df.gsva.nada.selected[,c("Sample","Condition")],
                                          1, function(x)
                                            paste(c(strsplit(x[1],"_")[[1]][1:2],x[2]),collapse="_"))

df.scores.nada.gsva <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition
                                                 %in% c("C1","C1A","C1B","C2","C4","C9")),]
df.scores.nada.gsva$SampleID <- apply(df.scores.nada.gsva[,c("Sample","Condition")],
                                        1, function(x)
                                          paste(c(strsplit(x[1],"_")[[1]][1:2],x[2]),collapse="_"))

# matc1 <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C1"),-c(1,6,7)]
# matc1a <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C1A"),-c(1,6,7)]
# matc1b <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C1B"),-c(1,6,7)]
# matc2 <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C2"),-c(1,6,7)]
# matc4 <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C4"),-c(1,6,7)]
# rownames(matc1) <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C1"),]$SampleID
# rownames(matc1a) <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C1A"),]$SampleID
# rownames(matc1b) <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C1B"),]$SampleID
# rownames(matc2) <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C2"),]$SampleID
# rownames(matc4) <- df.scores.nada.gsva[which(df.scores.nada.gsva$Condition == "C4"),]$SampleID
# pc1 <- pheatmap(matc1)
# pc1a <- pheatmap(matc1a)
# pc1b <- pheatmap(matc1b)
# pc2 <- pheatmap(matc2)
# pc4 <- pheatmap(matc4)

df.gs <- df.gsva.nada.selected[order(df.gsva.nada.selected$Condition),]
df.gs$SampleID <- factor(df.gs$SampleID,
                         levels=unique(df.gs$SampleID))


pdf("plots.lineage27082021//lineage_signatures.GSVA.pdf",w=11,h=3)
ggplot(df.gs, aes(SampleID, variable,fill= value)) + 
  geom_tile()+
  #geom_text(aes(label=star), color="black", size=3) +
  scale_fill_distiller(palette = "RdBu")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
# theme(axis.title.x=element_blank(),
#       axis.text.x=element_blank(),
#       axis.ticks.x=element_blank())
dev.off()

## Now compare using boxplots:

df.fg$ConditionBin <- sapply(df.fg$Condition, function(x)
                             ifelse (x %in% c("C1A","C1B"),"C1",x))
df.gs$ConditionBin <- sapply(df.gs$Condition, function(x)
  ifelse (x %in% c("C1A","C1B"),"C1",x))
my_comparisons <- list(c("C1","C2"),
                       c("C2","C4"),
                       c("C1","C4"),
                       c("C1","C9"))
pdf("plots.lineage27082021//violinplots_lineage_sigs_NES.pdf",w=10,h=7)
ggviolin(df.fg, x = "ConditionBin", 
         y = "NES",
         palette = c("grey","#1b9e77", "#2d6353", "#d95f02","#7570b3"),
         fill = "ConditionBin",alpha = 0.7,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~pathway,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Normalised enrichment score")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()

pdf("plots.lineage27082021//violinplots_lineage_sigs_GSVA.pdf",w=10,h=7)
ggviolin(df.gs, x = "ConditionBin", 
         y = "value",
         palette = c("grey","#1b9e77", "#2d6353", "#d95f02","#7570b3"),
         fill = "ConditionBin",alpha = 0.7,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=2)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("GSVA score")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()



#### Balloon plots:
library(ggpubr)

# Define color palette
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")

library(dplyr)
df.select <- df.fg[,c("SampleID","ConditionBin","pathway","NES")]
df.select.summarised <- data.frame(df.select%>%
                                     group_by(ConditionBin,pathway)%>% 
                                     summarise(Mean=mean(NES), Median=median(NES)))

pdf("plots.lineage27082021//ballonplot.lineage.NES.mean.pdf",w=6,h=3)
ggballoonplot(df.select.summarised,
              x = "pathway", y = "ConditionBin",
              fill = "Mean") +
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()

pdf("plots.lineage27082021//ballonplot.lineage.NES.median.pdf",w=6,h=3)
ggballoonplot(df.select.summarised,
              x = "pathway", y = "ConditionBin",
              fill = "Median") +
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()

df.select <- df.gs[,c("SampleID","ConditionBin","variable","value")]
df.select.summarised <- data.frame(df.select%>%
                                     group_by(ConditionBin,variable)%>% 
                                     summarise(Mean=mean(value), Median=median(value)))

pdf("plots.lineage27082021//ballonplot.lineage.GSVA.mean.pdf",w=6,h=3)
ggballoonplot(df.select.summarised,
              x = "variable", y = "ConditionBin",
              fill = "Mean") +
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()

pdf("plots.lineage27082021//ballonplot.lineage.GSVA.median.pdf",w=6,h=3)
ggballoonplot(df.select.summarised,
              x = "variable", y = "ConditionBin",
              fill = "Median") +
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()


#################

# Ternary plot:
library(ggtern)
df.nada.melt <- melt(df.scores.nada.gsva)
pdf("plots.20210211/ternaryPlot.2.pdf")
ggtern(data = df.scores.nada.gsva, 
       aes(x = Astrocytes, y = Oligodendrocytes, z = OPC)) +
  tern_limits(T=1, L=1, R=1)+
  geom_point(aes(fill = Condition),
             size = 4, 
             shape = 21, 
             color = "black") + 
  labs(fill = "Condition") + 
  #theme_tern_rgbw() + 
  theme(legend.position      = c(0,1), 
        legend.justification = c(0, 1))
dev.off()


#####################################
#####################################
# Now score immune infiltrates by GSVA:
load("mouse.mat.expr.log.noSamplesRemoved.RData")
imm.scores.mouse <- gsva(as.matrix(mat.expr.log),immuneList.all, 
                   method="gsva", kcdf="Gaussian",
                   mx.diff=TRUE, abs.ranking=FALSE)
df.scores.mouse <- data.frame(t(imm.scores.mouse))
df.scores.mouse$Sample <-rownames(df.scores.mouse)
df.scores.mouse <- merge(df.scores.mouse, info[,c("File_name",
                                                  "Sample_name",
                                                  "Condition",
                                                  "Grp","Tumour_origin")],
                         all.x=TRUE, all.y=FALSE, by.x="Sample",
                         by.y="File_name")
#save(df.scores.mouse, file="df.scores.gsva.MOUSE.plusC9.RData")
#write.csv(df.scores.mouse, "TMEscores.mouse.ssGSEA.20210301.csv",row.names=FALSE)

#load("df.scores.gsva.MOUSE.20210521.RData")

# keep only C1,C2,C4,C9:
df.scores.mouse <- df.scores.mouse[which(df.scores.mouse$Condition %in% c("C1A","C1B","C2","C4","C9")),]

################
# Read samples of interest:
samps <- rev(read.table("samplesOrderTMEheatmap.txt", header = FALSE, stringsAsFactors = FALSE)$V1)
df.scores.mouse.keep <- df.scores.mouse[which(df.scores.mouse$Sample_name %in% samps),]

selectedprogrammes <- c("CD8_Tcells",
                        "CD8_Effector.Memory","CD8_Chemokine.IFNg","CD8_Interferon",
                        "CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                        "CD4_Treg","T_regs",
                        "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                        "Neutrophils","Exhaustion")

df.melt <- melt(df.scores.mouse.keep,id.vars = c("Sample","Sample_name","Condition","Grp","Tumour_origin"))
df.select <- df.melt[which(df.melt$variable %in% selectedprogrammes),c("Sample_name","variable","value")]
df.select$Compartment <- sapply(df.select$variable,
                                           function(x) ifelse (x %in% c("Macrophages_M1","Macrophages_M2",
                                                                        "Microglia_Consensus","MDM_Consensus","Neutrophils"),"MyeloidCells",
                                                               ifelse(x %in% c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                                                               "CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
                                                                               "CD8_Tcells",
                                                                               "Memory_Tcells","T_regs","CD4_Treg"),"Lymphocytes","Exhaustion")))

df.compartment<- df.select#[which((df.select.summarised$Condition %in% c("C1A","C1B"))&
#         (df.select.summarised$Compartment !="Other")),] 

df.compartment$variable <- factor(df.compartment$variable,
                                  levels=c("CD8_Tcells",
                                           "CD8_Effector.Memory","CD8_Chemokine.IFNg","CD8_Interferon",
                                           "CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                           "CD4_Treg","T_regs",
                                           "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                                           "Neutrophils","Exhaustion"))
df.select.TMEsummary <- data.frame(df.compartment%>%
                                     group_by(Sample_name,Compartment)%>% 
                                     summarise(Mean=mean(value)))
df.mouse.tmesum <- cast(Sample_name~Compartment,value='Mean',
                       fun.aggregate='mean',
                       data = df.select.TMEsummary)
df.mouse.tmesum$ImmuneScore <- apply(df.mouse.tmesum[,c("Lymphocytes","MyeloidCells")],1,
                                  function(x) mean(x[1],x[2]))
write.csv(df.mouse.tmesum, file="mouse.TMEscores.csv",row.names=FALSE)
df.toplot <- df.mouse.tmesum[,-1]
rownames(df.toplot) <- df.mouse.tmesum$Sample_name

##############################
### Generate heat map for mouse:
pdf("TMEscores.mouse.pdf",onefile=FALSE,w=12,h=5)
pheatmap(t(df.toplot[samps,]),cluster_rows = FALSE,
         show_colnames = TRUE)
dev.off()

df.melt <- melt(df.scores.mouse[which(df.scores.mouse$Condition
                                      %in% c("C1A","C1B","C2","C9")),])
df.melt[which(df.melt$Condition %in% c("C1A","C1B")),]$Condition = "C1"
df.melt$Condition <- factor(df.melt$Condition,
                            levels=c("C1","C2","C9"))
my_comparisons <- list(
                       c("C1","C2"),
                       c("C1","C9"),
                       c("C2","C9"))

df.melt.selected <- df.melt[which(df.melt$variable %in%
                c("CD8_Tcells",
                  "Memory_Tcells",
                  "T_cell_activation",
                  "APC_activation",
                  "Exhaustion",
                  "T_regs",
                  "MDM_Consensus",
                  "Bcells","Macrophages_M1",
                  "MDM_Paolo","Microglia_Paolo")),]
df.melt.selected$variable <- factor(df.melt.selected$variable,
                                    levels=c("CD8_Tcells",
                                             "CD4_Tcells",
                                             "Memory_Tcells",
                                             "T_cell_activation",
                                             "APC_activation",
                                             "Exhaustion",
                                             "T_regs",
                                             "MDM_Consensus",
                                             "Bcells","Macrophages_M1",
                                             "MDM_Paolo","Microglia_Paolo"
                                    ))

df.melt.selected$Condition <- factor(df.melt.selected$Condition,
                            levels=c("C1","C2","C9"))


pdf("plots.C9/TMEcompared.MOUSE.GSVA.allCells.C1C2C9.pdf",w=20,h=10)
ggviolin(df.melt, x = "Condition", 
         y = "value",
         palette = c("#1b9e77",  "#d95f02","grey"),
         fill = "Condition",alpha = 0.7,
         add = "boxplot", add.params = list(fill = "white"))+
  facet_wrap(~variable,scales = "free",nrow=5)+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
  xlab("")+
  ylab("Score")+
  theme(axis.text.x=element_blank(),
        axis.ticks.x = element_blank())
dev.off()


#############################################
#### Create balloon plot of TME in mouse:

# selectedprogrammes.old <- c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
#                         "CD8_Tcells",
#                         "Memory_Tcells","CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
#                         "T_regs","CD4_Treg",
#                         "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
#                         "Neutrophils")

selectedprogrammes <- c("CD8_Tcells",
                        "CD8_Effector.Memory",
                        "CD4_Tcells","CD4_Effector.Memory",
                            "T_regs",
                            "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                            "Neutrophils","Exhaustion")

library(ggpubr)

# Define color palette
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")

library(dplyr)
df.melt <- melt(df.scores.mouse[which(df.scores.mouse$Condition
                                      %in% c("C1A","C1B","C2","C9")),])
df.select <- df.melt[which(df.melt$variable %in% selectedprogrammes),]
df.select.summarised <- data.frame(df.select%>%
  group_by(Condition,variable)%>% 
  summarise(Mean=mean(value), Median=median(value)))
#

df.select.summarised$Compartment <- sapply(df.select.summarised$variable,
                                           function(x) ifelse (x %in% c("Macrophages_M1","Macrophages_M2",
                                                                        "Microglia_Consensus","MDM_Consensus","Neutrophils"),"MyeloidCells",
                                           ifelse(x %in% c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                                           "CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
                                                           "CD8_Tcells",
                                                           "Memory_Tcells","T_regs","CD4_Treg"),"Lymphocytes","Exhaustion")))

df.compartment<- df.select.summarised[which((df.select.summarised$Condition %in% c("C1A","C1B","C9"))&
                             (df.select.summarised$Compartment !="Other")),] 

df.compartment$variable <- factor(df.compartment$variable,
                                        levels=c("CD8_Tcells","CD8_Effector.Memory","CD4_Tcells","CD4_Effector.Memory",
                                                 "T_regs",
                                                 "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                                                 "Neutrophils","Exhaustion"))

pdf("plots.C9/ballonplot.mouseTME.byCompartment.GSVA.plusC9.pdf",w=6,h=3)
ggballoonplot(df.compartment,
              x = "variable", y = "Condition",
              fill = "Mean") +
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()

####
## Only select matched cases:

c9info <- read.xls("20210116_Adi_RNAseq_CNA_Master-4-onlyprimary-C9.xlsx")[,c("File_name","Tumour.origin.","TumourID","Condition")]

df.toplot$Sample <- rownames(df.toplot)
df.toplot.matched <- merge(df.toplot,
                         c9info,
                         by.x=c("Sample"),
                         by.y=c("File_name"),
                         all.x=FALSE,all.y=FALSE)
df.toplot.matched <- merge(df.toplot.matched,
                           unique(df.melt[,c("Sample","Condition")]),
                           by.x="Sample",by.y="Sample")
annot <- df.toplot.matched[,c("Condition.y","TumourID")]
rownames(annot) <- df.toplot.matched$Sample
df.toplot.matched.pheat <- df.toplot.matched[,2:50]
rownames(df.toplot.matched.pheat) <- df.toplot.matched$Sample

pdf("plots.C9/heatmap.onlyMatchedCases.pdf",w=8,h=10)
 pheatmap(t(df.toplot.matched.pheat),
           annotation=annot, 
           #annotation_colors = ann_colors,
           show_colnames = TRUE, show_rownames = TRUE)
dev.off()

pdf("plots.C9/heatmap.onlyMatchedCases.selectedProgrammes.pdf",w=8,h=6)
pheatmap(t(df.toplot.matched.pheat[,selectedprogrammes]),
         annotation=annot, 
         #annotation_colors = ann_colors,
         show_colnames = TRUE, show_rownames = TRUE)
dev.off()

## Ordered:

df.toplot.ordered <- df.toplot.matched[with(df.toplot.matched,order(TumourID, Condition.y)),]
annot <- df.toplot.ordered[,c("Condition.y","TumourID")]
rownames(annot) <- df.toplot.ordered$Sample
df.toplot.ordered.pheat <- df.toplot.ordered[,2:50]
rownames(df.toplot.ordered.pheat) <- df.toplot.ordered$Sample

pdf("plots.C9/heatmap.onlyMatchedCases.ordered.pdf",w=8,h=10)
pheatmap(t(df.toplot.ordered.pheat),
         annotation=annot, cluster_cols = FALSE,
         #annotation_colors = ann_colors,
         show_colnames = TRUE, show_rownames = TRUE)
dev.off()

pdf("plots.C9/heatmap.onlyMatchedCases.selectedProgrammes.ordered.pdf",w=8,h=6)
pheatmap(t(df.toplot.ordered.pheat[,selectedprogrammes]),
         annotation=annot, cluster_cols = FALSE,
         #annotation_colors = ann_colors,
         show_colnames = TRUE, show_rownames = TRUE)
dev.off()


