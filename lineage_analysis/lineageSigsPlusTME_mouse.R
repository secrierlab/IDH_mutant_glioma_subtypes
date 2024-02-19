########################
##### This script defines signatures and calculates the TME composition in mouse samples.

library(pheatmap)
library(reshape)
library(GSVA)
library(ggpubr)
library(gdata)
library(erer)
library(ComplexHeatmap)


##### Read lineage signatures:
load("lineagesigs_suva1.RData")
load("lineagesigs_suva2.RData")
load("lineagesigs_verhaark.RData")

### The Venteicher sigs:
lineage <- read.xls("lineageSigs_2200221.xlsx",
                    header=TRUE, stringsAsFactors=FALSE)

lineageList <- NULL
lineageList$Astro.Venteicher <- setdiff(unique(trimws(lineage$Astro_Venteicher)),"")
lineageList$Oligo.Venteicher <- setdiff(unique(trimws(lineage$Oligo_Venteicher)),"")
lineageList$Stemness.Venteicher <- setdiff(unique(trimws(lineage$Stemness_Venteicher)),"")

lineage.IDH <- list(setdiff(unique(trimws(lineage$IDH.A_accountedByGenetics)),""),
                    setdiff(unique(trimws(lineage$IDH.A_.notaccounted.ByGenetics)),""),
                    setdiff(unique(trimws(lineage$IDH.O_accountedByGenetics)),""),
                    setdiff(unique(trimws(lineage$IDH.O_notaccountedByGenetics)),""))
names(lineage.IDH) <- c("IDH.A_AccGenetics",
                        "IDH.A_.NotAccGenetics",
                        "IDH.O_AccGenetics",
                        "IDH.O_NotAccGenetics")
                    

### The Parada sigs:
lineage <- read.xls("luisparada/paradaSigs.xlsx",
                    header=TRUE, stringsAsFactors=FALSE)

lineage.parada <- NULL
lineage.parada$Astro.Parada <- toupper(setdiff(unique(trimws(lineage$Astro)),""))
lineage.parada$OLS.Parada <- toupper(setdiff(unique(trimws(lineage$OLS)),""))
lineage.parada$Neurons.Parada <- toupper(setdiff(unique(trimws(lineage$Neurons)),""))


### Nada sigs:

lin.nada <- read.xls("nada/celltypes1.xlsx")
lineage.nada <- NULL
lineage.nada$Astro.Nada <- toupper(setdiff(unique(trimws(lin.nada$Astrocytes)),""))
lineage.nada$Oligo.Nada <- toupper(setdiff(unique(trimws(lin.nada$Oligodendrocytes)),""))


### Build the proneural/neural/mes classif:
# sigs.classif <- list()
# sigs.classif$Proneural <- c("DLL3","NKX2-2","SOX2","ERBB3","OLIG2")
# sigs.classif$Neural <- c("FBXO3","GABRB2","SNCG","MBP")
# sigs.classif$Classical <- c("FGFR3","PDGFA","EGFR","AKT2","NES")
# sigs.classif$Mesenchymal <- c("CASP1","CASP4","CASP5","CASP8",
#                               "ILR4","CHI3L1","TRADD","TLR2","TLR4","RELB")
# save(sigs.classif, file="sigs.classif.gbm.RData")

### Verhaak sigs:
lin.verhaak <- read.xls("verhaak_cancercell2017/verhaakSigs.xlsx")
lineage.verhaak <- NULL
lineage.verhaak$Classical <- setdiff(unique(trimws(lin.verhaak$Classical)),"")
lineage.verhaak$Classical <- setdiff(unique(trimws(lin.verhaak$Proneural)),"")
lineage.verhaak$Mesenchymal <- setdiff(unique(trimws(lin.verhaak$Mesenchymal)),"")



lists <- list(lineagesigs_verhaark,lineageList,lineage.parada,
              lineagesigs_suva2, lineagesigs_suva1,
           lineage.nada,lineage.IDH)
names(lists) <- c("Verhaak","Venteicher","Parada","Suva2","Suva1","Nada","IDH_genetics")

##################################
######### Mouse data
####################################

### First, read in expression data:
#load("mouse.mat.expr.log.keepAllC9.RData")

expr <- read.csv("rnaseq.mouse.latest/Adi_IDH_r1234.tpm.csv")
x1 <- gsub("\\.", "-", colnames(expr))
x2 <- gsub("X2", "2", x1)
colnames(expr) <- x2


### Also sample info:
info <- read.xls("latest/20220315_Adi_RNAseq_Master.xlsx")
info <- info[which(info$Remove != "Yes"),]

mat.expr <- expr[,which(colnames(expr) %in% info$File_name)]
mat.expr <- mat.expr[ !duplicated(expr$Gene_name), ]
rownames(mat.expr) <- toupper(expr$Gene_name[!duplicated(expr$Gene_name)])
mat.expr.log <- log2(mat.expr+1)

#####################
### Produce all plots:

methods = c("gsva","ssgsea")

i<-0

for (lin.mouse in list(lineagesigs_verhaark,lineageList,lineage.parada,
                       lineagesigs_suva2, lineagesigs_suva1,
                       lineage.nada,lineage.IDH)) {
  
  i<-i+1
  print(i)
  #names(lin.mouse) <- names(lists[i])
  
  print("Genes missing:")
  print(setdiff(unlist(lin.mouse),rownames(mat.expr)))
  
  setwd("plots.lineagesMouse.May2022/")
  setwd(names(lists[i]))
  
  wid1=2
  hei1=5
  wid2=5
  hei2=4
  if (names(lists[i]) %in% c("Suva1","IDH_genetics")) {
    wid1=3
    hei1=5
    wid2=10
    hei2=3
  }
  
  for (m in methods) {
    lineage.mouse <- gsva(as.matrix(mat.expr.log),lin.mouse, 
                             method=m, kcdf="Gaussian",
                             mx.diff=TRUE, abs.ranking=FALSE)
    df.scores.lineage <- data.frame(t(lineage.mouse))
    df.scores.lineage$Sample <- rownames(df.scores.lineage)
    df.scores.lineage <- merge(df.scores.lineage, info[,c("File_name","Condition","Grp")],
                                 by.x="Sample", by.y="File_name",
                                 all.x=FALSE, all.y=FALSE)
    df.scores.lineage.keep <- df.scores.lineage[which(df.scores.lineage$Condition %in% c("C1A","C1B","C2","C4")),]
    
    
    df.scores.lineage.keep$Grp <-factor(df.scores.lineage.keep$Grp,
                                             levels=c("C1","C2","C4"))

    pdf(paste0(names(lists[i]),".",m,".lineageHeatmap.mouse.C1C2C4.pdf"),w=wid1,h=hei1)
    draw(Heatmap(as.matrix(df.scores.lineage.keep[,c(2:(ncol(df.scores.lineage.keep)-2))]), 
                 split=df.scores.lineage.keep$Grp, show_row_names=FALSE,
                 cluster_columns =FALSE), heatmap_legend_side="left")
    dev.off()
    
    df.scores.lineage.keep2 <- df.scores.lineage.keep[which(df.scores.lineage.keep$Grp=="C1"),]
    pdf(paste0(names(lists[i]),".",m,".lineageHeatmap.mouse.C1AC1B.pdf"),w=wid1,h=hei1)
    draw(Heatmap(as.matrix(df.scores.lineage.keep2[,c(2:(ncol(df.scores.lineage.keep2)-2))]), 
                 split=df.scores.lineage.keep2$Condition, show_row_names=FALSE,
                 cluster_columns =FALSE), heatmap_legend_side="left")
    dev.off()
    
    df.gs <- melt(df.scores.lineage.keep)
    my_comparisons <- list(c("C1","C2"),
                        c("C2","C4"),
                        c("C1","C4"))
    pdf(paste0(names(lists[i]),".",m,".lineageViolinplot.mouse.C1C2C4.pdf"),w=wid2,h=hei2)
    print(ggviolin(df.gs, x = "Grp", 
             y = "value",
             palette = c("#1b9e77", "#d95f02","#7570b3"),
             fill = "Grp",alpha = 0.7,
             add = "boxplot", add.params = list(fill = "white"))+
      facet_wrap(~variable,scales = "free",nrow=1)+
      stat_compare_means(comparisons = my_comparisons, label = "p.signif")+
      xlab("")+
      ylab(paste0(m,"score"))+
      theme(axis.text.x=element_blank(),
            axis.ticks.x = element_blank()))
    dev.off()
    
    df.gs2 <- melt(df.scores.lineage.keep2)
    my_comparisons2 <- list(c("C1A","C1B"))
    pdf(paste0(names(lists[i]),".",m,".lineageViolinplot.mouse.C1AC1B.pdf"),w=wid2,h=hei2)
    print(ggviolin(df.gs2, x = "Condition", 
             y = "value",
             fill = "Condition",alpha = 0.7,
             add = "boxplot", add.params = list(fill = "white"))+
      facet_wrap(~variable,scales = "free",nrow=1)+
      stat_compare_means(comparisons = my_comparisons2, label = "p.signif")+
      xlab("")+
      ylab(paste0(m,"score"))+
      theme(axis.text.x=element_blank(),
            axis.ticks.x = element_blank()))
    dev.off()
    
  }
  
  setwd("../../.")

}

### Finally, score TME:
load("immuneListALL.BRAIN.update01042021.RData")
immuneList.all$CD8_Treg <- c(immuneList.all$T_regs,"CD8A","CD8B")
imm.scores.mouse <- gsva(as.matrix(mat.expr.log),immuneList.all, 
                         method="gsva", kcdf="Gaussian",
                         mx.diff=TRUE, abs.ranking=FALSE)
df.scores.mouse <- data.frame(t(imm.scores.mouse))
df.scores.mouse$Sample <-rownames(df.scores.mouse)
df.scores.mouse <- merge(df.scores.mouse, info[,c("File_name","Condition","Grp")],
                         by.x="Sample", by.y="File_name",
                         all.x=FALSE, all.y=FALSE)
df.scores.mouse.keep <- df.scores.mouse[which(df.scores.mouse$Condition %in% c("C1A","C1B")),]


selectedprogrammes <- c("CD8_Tcells",
                        "CD8_Effector.Memory","CD8_Treg",
                        "CD8_Interferon",
                        "CD4_Tcells","CD4_Effector.Memory","CD4_Treg",
                        "CD4_Interferon",
                        "Microglia_Consensus","MDM_Consensus","Macrophages_M1","Macrophages_M2",
                        "Neutrophils","Exhaustion")
df.melt <- melt(unique(df.scores.mouse.keep[,c("Sample",selectedprogrammes,"Condition")]),
                id.vars = c("Sample","Condition"))

library(dplyr)
df.select.toplot <- as.data.frame(df.melt%>%
                                    group_by(Condition,variable)%>% 
                                    summarise(Enrichment=mean(value)))

df.select.toplot<-data.frame(df.select.toplot)
df.select.toplot$Condition <- as.character(df.select.toplot$Condition)

df.select.toplot$Condition <- factor(df.select.toplot$Condition,
                                     levels=c("C1A",
                                              "C1B"))

my_cols <- c("#0D0887FF", "#6A00A8FF", "#B12A90FF",
             "#E16462FF", "#FCA636FF", "#F0F921FF")
pdf("baloonplots.mouse/TME.baloonplot.mouse.C1AC1B.pdf",w=6,h=3)
ggballoonplot(df.select.toplot,
              x = "variable", y = "Condition",
              fill = "Enrichment")+
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()


## Also broad programmes:
df.melt$Compartment <- sapply(df.melt$variable,
                                           function(x) ifelse (x %in% c("Macrophages_M1","Macrophages_M2",
                                                                        "Microglia_Consensus","MDM_Consensus","Neutrophils"),"MyeloidCells",
                                                               ifelse(x %in% c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                                                               "CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
                                                                               "CD8_Tcells","CD8_Treg",
                                                                               "Memory_Tcells","T_regs","CD4_Treg"),"Lymphocytes","Exhaustion")))


df.select.toplot <- as.data.frame(df.melt%>%
                                    group_by(Condition,Compartment)%>% 
                                    summarise(Enrichment=mean(value)))

df.select.toplot<-data.frame(df.select.toplot)
df.select.toplot$Compartment <- as.character(df.select.toplot$Compartment)

df.select.toplot$Compartment <- factor(df.select.toplot$Compartment,
                                     levels=c("Lymphocytes","MyeloidCells",
                                              "Exhaustion"))

pdf("baloonplots.mouse/TME.baloonplot.mouse.C1AC1B_broadProgrammes.pdf",w=4,h=3)
ggballoonplot(df.select.toplot,
              x = "Compartment", y = "Condition",
              fill = "Enrichment")+
  scale_fill_gradientn(colors = my_cols)+
  guides(size = FALSE)
dev.off()


# Print out the data:
df.scores.mouse$Lymphocytes <- apply(df.scores.mouse[,c("CD4_Tcells","CD4_Effector.Memory","CD4_Cytotoxicity_Chemokine","CD4_Interferon",
                                                                  "CD8_Interferon","CD8_Effector.Memory","CD8_Chemokine.IFNg",
                                                                  "CD8_Tcells","CD8_Treg",
                                                                  "Memory_Tcells","T_regs","CD4_Treg")], 1, function(x) mean(x))
df.scores.mouse$MyeloidCells <- apply(df.scores.mouse[,c("Macrophages_M1","Macrophages_M2",
                                                                   "Microglia_Consensus","MDM_Consensus","Neutrophils")], 1, function(x) mean(x))

df.scores.mouse$ImmuneScore <- (df.scores.mouse$Lymphocytes+df.scores.mouse$MyeloidCells)/2
write.csv(df.scores.mouse, file="mouse.TME.09082022.csv", quote=FALSE, row.names=FALSE)
