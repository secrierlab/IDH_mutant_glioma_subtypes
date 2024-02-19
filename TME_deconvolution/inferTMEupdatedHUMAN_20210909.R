########################
##### This script defines signatures and calculates the TME composition in mouse samples.

library(GSVA)
library(gdata)
library(pheatmap)
library(erer)

load("immuneListALL.BRAIN.update01042021.RData")
load("df.scores.gsva.HUMAN.TCGA.11052021.RData")
load("tcga.gbm.mut.RData")
load("tcga.lgg.mut.RData")

mat.expr.gbm <- tcga.gbm[,-c(1:6,ncol(tcga.gbm))]
rownames(mat.expr.gbm) <- tcga.gbm$Sample
mat.expr.lgg <- tcga.lgg[,-c(1:6,ncol(tcga.gbm))]
rownames(mat.expr.lgg) <- tcga.lgg$Sample
mat.expr.tcga <- rbind(mat.expr.gbm,mat.expr.lgg)

mat.expr.tcga.log <- log2(mat.expr.tcga+1)

mat.exhaustion <- mat.expr.tcga.log[,which(colnames(mat.expr.tcga.log) %in% immuneList.all$Exhaustion)]
df.exh <- data.frame(mat.exhaustion)
df.exh$Sample <- rownames(df.exh)

#Merge with CD8/CD4Tcells:
df.mergedexh <- merge(df.exh, df.scores.human.tcga[,c("Sample","CD8_Tcells",
                                                 "CD8_Effector.Memory",
                                                 "CD4_Tcells","CD4_Effector.Memory")],
                      by.x="Sample",by.y="Sample",
                      all.x=FALSE, all.y=FALSE
)
df.m <- df.mergedexh[,-1]

library(corrplot)
testRes = cor.mtest(df.m, conf.level = 0.95)
M = cor(df.m)
pdf("plots.manuscript.20210909/correlations.exhaustionCD4CD8reduced.human.pdf",w=15,h=15)
corrplot.mixed(M, order = 'AOE')
dev.off()


load("TMEinfiltration_TCGA_GBMplusLGG.RData")

### Load C1A/B annotation:
tcgaannot <- read.table("c1abhuman/OS.Comp8_C1A_vs_C1B_0.05_1.gsva_MCP_0_LGG_TRIPLE_MUT.coxph_thrcox.0.001.up_dw.txt",
                        sep="\t",header=TRUE)
tcgaannot$Category <- sapply(tcgaannot$status_C1AC1B,
                             function(x) ifelse(x=="high","C1A","C1B"))
df.scores.human.tcga <- merge(df.scores.human.tcga,
                              tcgaannot[,c("Samples","Category")],
                              by.x="Case",by.y="Samples",
                              all.x=FALSE,all.y=FALSE)

### Generate heat map for human:
selectedprogrammes <- setdiff(colnames(df.scores.human.tcga)[3:51],
                              c("MDM_Joyce","Microglia_Joyce","Monocytes_Becher",
                                "Macrophage_Venteicher","Microglia_Venteicher",
                                "Microglia_Becher","Metabolic.escape","MDM.Joyce.mouse.from.STM",
                                "CD4_Cytotoxicity_Chemokine","CD8_ND","HighNK","CD8_Chemokine.IFNg",
                                "CD8_Cytotoxicity.NK","MDM_Becher","MDM1_Becher",
                                "MDM2_Becher","MDM3_Becher","MDM4_Becher"))
df.toplot <- df.scores.human.tcga[,which(colnames(df.scores.human.tcga) %in% selectedprogrammes)]
rownames(df.toplot) <- df.scores.human.tcga$Sample
annot <- df.scores.human.tcga[,rev(c(
  "Category","Exhaustion"))]
rownames(annot) <- rownames(df.toplot)
annot <- annot[,rev(c("Category","Exhaustion"))]
# ann_colors_orig <- list(Type=c(LGG="#ECBEB4",GBM="#7EB09B"),
#                         m_TP53=c("1"="black","0"="lightgrey"),
#                         m_ATRX=c("1"="black","0"="lightgrey"),
#                         `m_IDH1.R132H`=c("1"="black","0"="lightgrey"),
#                         `X1p.19q.codeletion`=c("codel"="black","non-codel"="lightgrey"))
# 

annot$Category <- factor(annot$Category,
                             levels=c("C1B","C1A"))
  
pdf("plots.20210702/heatmap.HUMAN.C1AC1B.pdf",onefile=FALSE,w=10,h=6)
pall <- pheatmap(t(df.toplot),
                 annotation_col=annot, 
                 cutree_cols = 2,
                 #annotation_colors = ann_colors_orig,
                 show_colnames = FALSE)
dev.off()


groups <- cutree(pall$tree_col,k = 2)
table(groups)
df.groups <- data.frame(groups)
df.groups$Sample <- rownames(df.groups)
colnames(df.groups) <- c("Group","Sample")

df.toplot$Sample <- rownames(df.toplot)
df.tme <- merge(df.toplot, df.groups,
                by.x="Sample",by.y="Sample",
                all.x=FALSE, all.y=FALSE)

annot$Sample <- rownames(annot)
df.annot <- merge(annot, df.groups,
                by.x="Sample",by.y="Sample",
                all.x=FALSE, all.y=FALSE)
tb <- table(df.annot[,c("Category","Group")])
fisher.test(tb)
# not signif

#################################################
### Generate baloon plot for human infiltrates:
selectedprogrammes_old <- c("CD8_Tcells",
                        "CD8_Effector.Memory",
                        "CD4_Tcells","CD4_Effector.Memory",
                        "T_regs",
                        "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                        "Neutrophils","Exhaustion")
selectedprogrammes <- c("CD8_Tcells",
                            "CD8_Effector.Memory","CD8_Chemokine.IFNg","CD8_Interferon",
                            "CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                            "CD4_Treg","T_regs",
                            "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                            "Neutrophils","Exhaustion")

library(ggpubr)
df.tme <- merge(tme.TCGA, tcgaannot[,c("Samples","Category")],
                by.x="Case", by.y="Samples",
                all.x=FALSE, all.y=FALSE)


# Define color palette
my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")

library(dplyr)
library(reshape)
# selectedprogrammes <- c("CD8_Tcells","CD8_Effector.Memory","CD4_Tcells","CD4_Effector.Memory",
#                         "T_regs",
#                         "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
#                         "Neutrophils","Exhaustion")
df.melt <- melt(df.tme,id.vars = c("Case","Sample","Category"))
df.select <- df.melt[which(df.melt$variable %in% selectedprogrammes),c("Case","Category","variable","value")]
df.select.summarised <- data.frame(df.select%>%
                                     group_by(Category,variable)%>% 
                                     summarise(Mean=mean(value), Median=median(value)))
#

df.select.summarised$Compartment <- sapply(df.select.summarised$variable,
                                           function(x) ifelse (x %in% c("Macrophages_M1","Macrophages_M2",
                                                                        "Microglia_Consensus","MDM_Consensus","Neutrophils"),"MyeloidCells",
                                                               ifelse(x %in% c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                                                               "CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
                                                                               "CD8_Tcells",
                                                                               "Memory_Tcells","T_regs","CD4_Treg"),"Lymphocytes","Exhaustion")))

df.compartment<- df.select.summarised#[which((df.select.summarised$Condition %in% c("C1A","C1B"))&
                                     #         (df.select.summarised$Compartment !="Other")),] 

df.compartment$variable <- factor(df.compartment$variable,
                                  levels=c("CD8_Tcells",
                                           "CD8_Effector.Memory","CD8_Chemokine.IFNg","CD8_Interferon",
                                           "CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                           "CD4_Treg","T_regs",
                                           "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                                           "Neutrophils","Exhaustion"))

pdf("figs_manuscript/ballonplot.humanTME.byC1ABstatus_fromSurvivalAnalysis.GSVA_extendedTMElist.pdf",w=6,h=3)
ggballoonplot(df.compartment,
              x = "variable", y = "Category",
              fill = "Mean") +
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()

## Compare as boxplots:
# selectedprogrammes <- setdiff(colnames(df.scores.human.tcga)[3:51],
#                               c("MDM_Joyce","Microglia_Joyce","Monocytes_Becher",
#                                 "Macrophage_Venteicher","Microglia_Venteicher",
#                                 "Microglia_Becher","Metabolic.escape","MDM.Joyce.mouse.from.STM",
#                                 "CD4_Cytotoxicity_Chemokine","CD8_ND","HighNK","CD8_Chemokine.IFNg",
#                                 "CD8_Cytotoxicity.NK","MDM_Becher","MDM1_Becher",
#                                 "MDM2_Becher","MDM3_Becher","MDM4_Becher"))
df.melt <- melt(df.tme,id.vars = c("Case","Sample","Category"))
df.select <- df.melt[which(df.melt$variable %in% selectedprogrammes),c("Case","Category","variable","value")]

pdf("figs_manuscript/boxplots.humanTME.compared_byC1ABstatus_updatedTMElist.pdf",w=10,h=16)
my_comparisons <- list( c("C1A","C1B" ))
ggboxplot(df.select, x = "Category", y = "value",
          color = "Category", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")+ # Add pairwise comparisons p-value
  facet_wrap(~variable,scale="free",nrow=4)
dev.off()

### Define immune score, lymphoid and myeloid:

df.select$Compartment <- sapply(df.select$variable,
                                function(x) ifelse (x %in% c("Macrophages_M1","Macrophages_M2",
                                                             "Microglia_Consensus","MDM_Consensus","Neutrophils"),"MyeloidCells",
                                                    ifelse(x %in% c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                                                    "CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
                                                                    "CD8_Tcells",
                                                                    "Memory_Tcells","T_regs","CD4_Treg"),"Lymphocytes","Exhaustion")))
df.select.TMEsummary <- data.frame(df.select%>%
                                     group_by(Case,Compartment)%>% 
                                     summarise(Mean=mean(value)))
df.tcga.tmesum <- cast(Case~Compartment,value='Mean',
                       fun.aggregate='mean',
                       data = df.select.TMEsummary)

# read guidantonio table:
tb.tcga.guid <- read.xls("C1a-C1b-new-patient-annotations-20210715.xlsx")
selectedprogrammes <- c("CD8_Tcells",
                        "CD8_Effector.Memory","CD8_Chemokine.IFNg","CD8_Interferon",
                        "CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                        "CD4_Treg","T_regs",
                        "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                        "Neutrophils","Exhaustion")
df.melt <- melt(tme.TCGA,id.vars = c("Case","Sample"))
df.select <- df.melt[which(df.melt$variable %in% selectedprogrammes),c("Case","variable","value")]
df.select$Compartment <- sapply(df.select$variable,
                                function(x) ifelse (x %in% c("Macrophages_M1","Macrophages_M2",
                                                             "Microglia_Consensus","MDM_Consensus","Neutrophils"),"MyeloidCells",
                                                    ifelse(x %in% c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                                                    "CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
                                                                    "CD8_Tcells",
                                                                    "Memory_Tcells","T_regs","CD4_Treg"),"Lymphocytes","Exhaustion")))
df.select.TMEsummary <- data.frame(df.select%>%
                                     group_by(Case,Compartment)%>% 
                                     summarise(Mean=mean(value)))
df.tcga.tmesum <- cast(Case~Compartment,value='Mean',
                       fun.aggregate='mean',
                       data = df.select.TMEsummary)


tb.tcga.guid <- merge(tb.tcga.guid, 
                      df.tcga.tmesum,
                      by.x="Samples",by.y="Case",
                      all.x=FALSE,all.y=FALSE)
tb.tcga.guid$ImmuneScore <- apply(tb.tcga.guid[,c("Lymphocytes","MyeloidCells")],1,
                                   function(x) mean(x[1],x[2]))
tb.tcga.guid <- merge(tb.tcga.guid, 
                      tme.TCGA[,-2],
                      by.x="Samples",by.y="Case",
                      all.x=FALSE,all.y=FALSE)
write.csv(tb.tcga.guid, row.names = FALSE,
          file="C1a-C1b-new-patient-annotations-20210726.plusTME.csv")

##############
### Check other clinical variables:
library(TCGAbiolinks)
query <- GDCquery(project = "TCGA-LGG", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab")
GDCdownload(query)
clinical.BCRtab.all <- GDCprepare(query)
names(clinical.BCRtab.all)

query <- GDCquery(project = "TCGA-LGG", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  file.type = "clinical_patient_lgg")
GDCdownload(query)
clinical <- data.frame(GDCprepare(query)$clinical_patient_lgg)[-c(1:3),]

df.tme$Case <- sapply(df.tme$Sample,
                      function(x) paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
clin.merged <- merge(clinical,
                     df.tme[,c("Case","Category")],
                     by.x="bcr_patient_barcode", by.y="Case",
                     all.x=FALSE, all.y=FALSE)
fisher.test(table(clin.merged[,c("Category","gender")]))

fisher.test(table(clin.merged[,c("Category","tumor_grade")]))
#p-value = 0.009659
#odds ratio 
#2.414761
tumor_grade
# Category G2 G3
# C1A 51 22
# C1B 42 44
### 2.4 enrichment of C1A G2 tumours compared to G3 - might explain better prognosis?


query <- GDCquery(project = "TCGA-LGG", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  file.type = "clinical_follow_up_v1.0_lgg")
GDCdownload(query)
clinical <- data.frame(GDCprepare(query)$clinical_follow_up_v1.0_lgg)[-c(1:3),]

clin.merged <- merge(clinical,
                     df.tme[,c("Case","Category")],
                     by.x="bcr_patient_barcode", by.y="Case",
                     all.x=FALSE, all.y=FALSE)
tb1 <- table(clin.merged[,c("Category","treatment_outcome_at_tcga_followup")])
fisher.test(matrix(c(18,46,26,62),nrow=2))

tb2 <- table(clin.merged[,c("Category","treatment_outcome_first_course")])
fisher.test(matrix(c(22,51,33,75),nrow=2))

tb3 <- table(clin.merged[,c("Category","radiation_treatment_adjuvant")])
fisher.test(tb3[,c(4,5)])
# p-value = 0.0001526
# alternative hypothesis: true odds ratio is not equal to 1
# 95 percent confidence interval:
#   1.753410 7.218937
# sample estimates:
#   odds ratio 
# 3.520689 
# radiation_treatment_adjuvant
# Category NO YES
# C1A 34  42
# C1B 21  92
### 3.5 fold enrichment of adjuvant radiation treatment in C1B

query <- GDCquery(project = "TCGA-LGG", 
                  data.category = "Clinical",
                  data.type = "Clinical Supplement", 
                  data.format = "BCR Biotab",
                  file.type = "clinical_drug_lgg")
GDCdownload(query)
clinical <- data.frame(GDCprepare(query)$clinical_drug_lgg)[-c(1:3),]

clin.merged <- merge(clinical,
                     df.tme[,c("Case","Category")],
                     by.x="bcr_patient_barcode", by.y="Case",
                     all.x=FALSE, all.y=FALSE)


### Aneuploidy from https://www.sciencedirect.com/science/article/pii/S1535610818301119#sec5

library(gdata)
aneu <- read.xls("TCGA_aneuploidy.xlsx")
aneu$Case <- sapply(aneu$Sample, function(x)
  paste(strsplit(x,"-")[[1]][1:3],collapse="-"))
aneu.merged <- merge(aneu,
                     df.tme[,c("Case","Category")],
                     by.x="Case", by.y="Case",
                     all.x=FALSE, all.y=FALSE)
tb1 <- table(aneu.merged[,c("Category","Genome_doublings")])
fisher.test(tb1[,1:2])
# no difference in WGD

wilcox.test(aneu.merged[which(aneu.merged$Category=="C1A"),]$AneuploidyScore.AS,
aneu.merged[which(aneu.merged$Category=="C1B"),]$AneuploidyScore.AS)
# no difference

pdf("figs_manuscript/aneuploidy.humanC1AB.pdf")
my_comparisons <- list( c("C1A","C1B" ))
ggboxplot(aneu.merged, x = "Category", y = "AneuploidyScore.AS.",
          color = "Category", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
dev.off()

pdf("figs_manuscript/aneuploidy_del.humanC1AB.pdf")
my_comparisons <- list( c("C1A","C1B" ))
ggboxplot(aneu.merged, x = "Category", y = "AS_del",
          color = "Category", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
dev.off()

pdf("figs_manuscript/aneuploidy_amp.humanC1AB.pdf")
my_comparisons <- list( c("C1A","C1B" ))
ggboxplot(aneu.merged, x = "Category", y = "AS_amp",
          color = "Category", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
dev.off()

pdf("figs_manuscript/nonsilentmutations.humanC1AB.pdf")
my_comparisons <- list( c("C1A","C1B" ))
ggboxplot(aneu.merged, x = "Category", y = "Non.silentMutationsperMb",
          color = "Category", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
dev.off()

pdf("figs_manuscript/purity.humanC1AB.pdf")
aneu.merged$Purity <- as.numeric(aneu.merged$Purity)
my_comparisons <- list( c("C1A","C1B" ))
ggboxplot(aneu.merged, x = "Category", y = "Purity",
          color = "Category", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
dev.off()

pdf("figs_manuscript/leukocytes.humanC1AB.pdf")
aneu.merged$Leuk <- as.numeric(aneu.merged$Leuk)
my_comparisons <- list( c("C1A","C1B" ))
ggboxplot(aneu.merged, x = "Category", y = "Leuk",
          color = "Category", palette =c("#00AFBB", "#E7B800"),
          add = "jitter")+
  stat_compare_means(comparisons = my_comparisons,label="p.signif")
dev.off()
