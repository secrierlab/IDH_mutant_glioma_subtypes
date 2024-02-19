### Check whether PRC2, EZH2 and SUZ12 targets have corresponding changes in gene expression.

library(reshape)
library(gdata)

convertMouseGeneList2 <- function(x){
  
  require("biomaRt")
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl",version="102")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl",version="102")
  
  humanx = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  
  return(humanx)
  
}

# Read DE signatures for C1A/B:
sigs <- read.csv("../DESeq2_annotMar2022_v2_summary/Adi.ALL.log2FC.csv")
sigs$Adi.Comp1_C1_vs_C2.gz

# Read methylation DMRs:
dmrs <- read.delim("../project_adi-main/methylation_scripts/annotations_charts_ml_methylation/annotation_dmrs_adi_pval0.05_lfc1.txt")
# Select only promoters:
dmrs <- dmrs[which(dmrs$annot.type == "mm10_genes_promoters"),]

mapped <- convertMouseGeneList2(unique(dmrs$annot.symbol))
save(mapped, file="mapped.dmrs.mousetohuman.RData")

dmrs <- merge(dmrs, mapped,
              by.x="annot.symbol", by.y="MGI.symbol",
              all.x=FALSE, all.y=FALSE)

### Compile targets from different lists:
alltargets.up <- read.xls("Adi.Comp1_C1_vs_C2.TF.byType.up.xlsx")
alltargets.down <- read.xls("Adi.Comp1_C1_vs_C2.TF.byType.down.xlsx")
suz12.targets.up <- alltargets.up[which(grepl("SUZ12",toupper(alltargets.up$Term))),]
suz12.targets.down <- alltargets.down[which(grepl("SUZ12",toupper(alltargets.down$Term))),]
ezh2.targets.up <- alltargets.up[which(grepl("EZH2",toupper(alltargets.up$Term))),]
ezh2.targets.down <- alltargets.down[which(grepl("EZH2",toupper(alltargets.down$Term))),]

ezh2.up <-  unique(unlist(sapply(ezh2.targets.up$Genes, function(x) setdiff(unique(strsplit(x,";")[[1]]),""))))
suz12.up <- unique(unlist(sapply(suz12.targets.up$Genes, function(x) setdiff(unique(strsplit(x,";")[[1]]),""))))
ezh2.down <-  unique(unlist(sapply(ezh2.targets.down$Genes, function(x) setdiff(unique(strsplit(x,";")[[1]]),""))))
suz12.down <- unique(unlist(sapply(suz12.targets.down$Genes, function(x) setdiff(unique(strsplit(x,";")[[1]]),""))))
ezh2 <- c(ezh2.up,ezh2.down)
suz12 <- c(suz12.up, suz12.down)
prc2targets <- read.delim("PRC_EZH_SUZ12_targets/BENPORATH_PRC2_TARGETS.v2022.1.Hs.tsv")
prc2 <- setdiff(unique(strsplit(prc2targets$BENPORATH_PRC2_TARGETS[19],",")[[1]]),"")

dmr.prc2 <- dmrs[which(dmrs$HGNC.symbol %in% prc2),]
dmr.ezh2 <- dmrs[which((toupper(dmrs$HGNC.symbol) %in% ezh2)|(toupper(dmrs$annot.symbol) %in% ezh2)),]
dmr.suz12 <- dmrs[which((toupper(dmrs$HGNC.symbol) %in% suz12)|(toupper(dmrs$annot.symbol) %in% suz12)),]

## Among these up/downregulated targets, how many are overexpressed/underexpressed?

overexp <- sigs[which(sigs$Adi.Comp1_C1_vs_C2.gz>0),]$Gene
underexp <- sigs[which(sigs$Adi.Comp1_C1_vs_C2.gz<0),]$Gene

length(intersect(overexp, dmr.ezh2$annot.symbol))
length(intersect(underexp, dmr.ezh2$annot.symbol))

length(intersect(overexp, dmr.suz12$annot.symbol))
length(intersect(underexp, dmr.suz12$annot.symbol))

## Merge sigs table with dmr annot:
dmr.ezh2$TF <- "EZH2"
dmr.suz12$TF <- "SUZ12"
dmr.prc2$TF <- "PRC2"

dmrs <- rbind(dmr.ezh2,dmr.suz12,dmr.prc2)

sigs.plusdmr <- merge(sigs, unique(dmrs[,c("annot.symbol","log2FoldChange",
                                    "DynamicallyDetected","status","HGNC.symbol","TF")]),
                      by.x="Gene", by.y="annot.symbol",
                      all.x=FALSE, all.y=TRUE
                      )

library(pheatmap)
library(ggplot2)

pdf("EZH_SUZ12_targets.updated.dec2023/hypermethylatedTargets.DEexpression.pdf",w=20,h=3)
ggplot(sigs.plusdmr[which(sigs.plusdmr$status=="DMRs-up"),], 
       aes(x=Gene, y=TF, fill=Adi.Comp1_C1_vs_C2.gz))+
  geom_tile()+
  scale_fill_gradient2(low="#1D3557", high="#E63946", mid="white",
                       midpoint=0,na.value="black")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
  #theme_classic()
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/hypomethylatedTargets.DEexpression.pdf",w=10,h=3)
ggplot(sigs.plusdmr[which(sigs.plusdmr$status=="DMRs-down"),], 
       aes(x=Gene, y=TF, fill=Adi.Comp1_C1_vs_C2.gz))+
  geom_tile()+
  scale_fill_gradient2(low="#1D3557", high="#E63946", mid="white",
                       midpoint=0,na.value="black")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle = 45, hjust=1))
dev.off()


sigs.plusdmr$Expression <- sapply(sigs.plusdmr$Adi.Comp1_C1_vs_C2.gz,
                                  function(x) ifelse(is.na(x),"UNCHANGED",ifelse(x>0,"UP",ifelse(x<0,"DOWN","UNCHANGED"))))

tb.ezh2 <- table(sigs.plusdmr[which(sigs.plusdmr$TF=="EZH2"),c("status","Expression")])
tb.suz12 <- table(sigs.plusdmr[which(sigs.plusdmr$TF=="SUZ12"),c("status","Expression")])
tb.prc2<- table(sigs.plusdmr[which(sigs.plusdmr$TF=="PRC2"),c("status","Expression")])

write.csv(sigs.plusdmr[which(sigs.plusdmr$TF=="EZH2"),c("Gene","status","Expression")],
          file="EZH2targets.expressionMethylation.csv")
write.csv(sigs.plusdmr[which(sigs.plusdmr$TF=="PRC2"),c("Gene","status","Expression")],
          file="PRC2targets.expressionMethylation.csv")
write.csv(sigs.plusdmr[which(sigs.plusdmr$TF=="SUZ12"),c("Gene","status","Expression")],
          file="SUZ12targets.expressionMethylation.csv")


round(prop.table(tb.ezh2),2)
round(prop.table(tb.suz12),2)

require(vcd)

pdf("EZH_SUZ12_targets.updated.dec2023/ezh2targets.DMRsandExpression.pdf")
mosaic(tb.ezh2, shade=T, legend=T)
assoc(tb.ezh2, shade=T, legend=T)
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/suz12targets.DMRsandExpression.pdf")
mosaic(tb.suz12, shade=T, legend=T)
assoc(tb.suz12, shade=T, legend=T)
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/prc2targets.DMRsandExpression.pdf")
mosaic(tb.prc2, shade=T, legend=T)
assoc(tb.prc2, shade=T, legend=T)
dev.off()

# no significance by Fisher's exact test

pdf("EZH_SUZ12_targets.updated.dec2023/ezh2targets.DMRsandExpression.spineplot.pdf")
spineplot(tb.ezh2)
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/prc2targets.DMRsandExpression.spineplot.pdf")
spineplot(tb.prc2)
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/suz12targets.DMRsandExpression.spineplot.pdf")
spineplot(tb.suz12)
dev.off()


pdf("EZH_SUZ12_targets.updated.dec2023/ezh2targets.DMRsandExpression.contingencyPlot.pdf",w=4,h=3)
ggplot(melt(tb.ezh2),aes(x = status, y = Expression)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(label = value), color = "black", fontface = "bold", size = 6)+
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  theme_bw()
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/ezh2targets.DMRsandExpression.contingencyPlot.percentage.pdf",w=4,h=3)
ggplot(melt(round(prop.table(tb.ezh2),2)),aes(x = status, y = Expression)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(label = value), color = "black", fontface = "bold", size = 6)+
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  theme_bw()
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/suz12targets.DMRsandExpression.contingencyPlot.pdf",w=4,h=3)
ggplot(melt(tb.suz12),aes(x = status, y = Expression)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(label = value), color = "black", fontface = "bold", size = 6)+
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  theme_bw()
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/suz12targets.DMRsandExpression.contingencyPlot.percentage.pdf",w=4,h=3)
ggplot(melt(round(prop.table(tb.suz12),2)),aes(x = status, y = Expression)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(label = value), color = "black", fontface = "bold", size = 6)+
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  theme_bw()
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/prc2targets.DMRsandExpression.contingencyPlot.pdf",w=4,h=3)
ggplot(melt(tb.prc2),aes(x = status, y = Expression)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(label = value), color = "black", fontface = "bold", size = 6)+
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  theme_bw()
dev.off()

pdf("EZH_SUZ12_targets.updated.dec2023/prc2targets.DMRsandExpression.contingencyPlot.percentage.pdf",w=4,h=3)
ggplot(melt(round(prop.table(tb.prc2),2)),aes(x = status, y = Expression)) +
  geom_tile(aes(fill = value))+
  geom_text(aes(label = value), color = "black", fontface = "bold", size = 6)+
  scale_fill_gradient2(low="#1D3557", high="#E63946",
                       midpoint=0,na.value="black")+
  theme_bw()
dev.off()

#### Output top 20 DMR up, expression down:
sigs.plusdmr.dmrupexpdown <- sigs.plusdmr[which((sigs.plusdmr$status=="DMRs-up")&
                                                  (sigs.plusdmr$Expression=="DOWN")),]
ezh.out <- sigs.plusdmr.dmrupexpdown[which(sigs.plusdmr.dmrupexpdown$TF=="EZH2"),]
write.csv(ezh.out[rev(order(ezh.out$log2FoldChange)),][1:20,],
          file="EZH2.targets.DMRupEXPRdown.top20.csv", row.names=FALSE)
prc.out <- sigs.plusdmr.dmrupexpdown[which(sigs.plusdmr.dmrupexpdown$TF=="PRC2"),]
write.csv(prc.out[rev(order(prc.out$log2FoldChange)),][1:20,],
          file="PRC2.targets.DMRupEXPRdown.top20.csv", row.names=FALSE)
SUZ12.out <- sigs.plusdmr.dmrupexpdown[which(sigs.plusdmr.dmrupexpdown$TF=="SUZ12"),]
write.csv(SUZ12.out[rev(order(SUZ12.out$log2FoldChange)),][1:20,],
          file="SUZ12.targets.DMRupEXPRdown.top20.csv", row.names=FALSE)

#### Plot heat map of only DMRs that are up when expr is down:
pdf("EZH_SUZ12_targets.updated.dec2023/DMRupEXPRdown_Targets.pdf",w=25,h=3)
ggplot(sigs.plusdmr.dmrupexpdown[which(sigs.plusdmr.dmrupexpdown$status=="DMRs-up"),], 
       aes(x=Gene, y=TF, fill=Adi.Comp1_C1_vs_C2.gz))+
  geom_tile()+
  scale_fill_gradient2(low="#AD343E", high="#1D3557", mid="white",
                       midpoint=0,na.value="black")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle =90, hjust=1))
#theme_classic()
dev.off()


# DMR up EXP down targets of EZH2 and SUZ12:
length(unique(sigs.plusdmr.dmrupexpdown[which(sigs.plusdmr.dmrupexpdown$TF %in% c("SUZ12","EZH2","PRC2")),]$Gene))
### only 212 for EZH2 and SUZ12 alone, otherwise 243

## Now select only TSs or DDR regulators:
ts <- read.delim("Census_allThu Dec 21 09_48_56 2023.tsv")
ts <- toupper(unique(ts[which(grepl("TSG",ts$Role.in.Cancer)),]$Gene.Symbol))
ddr <- read.xls("~/Desktop/UCL/supervisedProjects/DanielJacobson/Rotation1-master/DDRpathways.xlsx")
ddr.genes <- toupper(setdiff(unique(ddr$Gene.ID),""))
mousegenes <- convertMouseGeneList2(unique(c(ts,ddr.genes)))
extra <- c("Atr", "Fbxo41", "Lmx1a", "Glipr2", "Plk1")


# Keep only TFs and DDR regulators and plot heat map of only DMRs that are up when expr is down:
sigs.plusdmr.dmrupexpdown.keep <- sigs.plusdmr.dmrupexpdown[which(toupper(sigs.plusdmr.dmrupexpdown$Gene) %in% toupper(c(mousegenes$MGI.symbol,extra))),]
pdf("EZH_SUZ12_targets.updated.dec2023/EZH2.SUZ12.PRC2_DMRupEXPRdown_Targets_TS_DDR.plusextra.pdf",w=8,h=2)
ggplot(sigs.plusdmr.dmrupexpdown.keep[which(sigs.plusdmr.dmrupexpdown.keep$status=="DMRs-up"),], 
       aes(x=Gene, y=TF, fill=Adi.Comp1_C1_vs_C2.gz))+
  geom_tile()+
  scale_fill_gradient2(low="#AD343E", high="#1D3557", mid="white",
                       midpoint=0,na.value="black")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle =90, hjust=1))
#theme_classic()
dev.off()

# The same but excluding PRC2 targets:
sigs.plusdmr.dmrupexpdown.keep2 <- sigs.plusdmr.dmrupexpdown.keep[which(sigs.plusdmr.dmrupexpdown.keep$TF %in% c("SUZ12","EZH2")),]

pdf("EZH_SUZ12_targets.updated.dec2023/EZH2.SUZ12_DMRupEXPRdown_Targets_TS_DDR.plusextra.pdf",w=8,h=2)
ggplot(sigs.plusdmr.dmrupexpdown.keep2[which(sigs.plusdmr.dmrupexpdown.keep2$status=="DMRs-up"),], 
       aes(x=Gene, y=TF, fill=Adi.Comp1_C1_vs_C2.gz))+
  geom_tile()+
  scale_fill_gradient2(low="#AD343E", high="#1D3557", mid="white",
                       midpoint=0,na.value="black")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle =90, hjust=1))
#theme_classic()
dev.off()

# Keep only TFs and DDR regulators and plot heat map of only DMRs that are up when expr is down, without the extra genes:
sigs.plusdmr.dmrupexpdown.keep <- sigs.plusdmr.dmrupexpdown[which(toupper(sigs.plusdmr.dmrupexpdown$Gene) %in% toupper(c(mousegenes$MGI.symbol))),]
pdf("EZH_SUZ12_targets.updated.dec2023/EZH2.SUZ12.PRC2_DMRupEXPRdown_Targets_TS_DDR.pdf",w=8,h=2)
ggplot(sigs.plusdmr.dmrupexpdown.keep[which(sigs.plusdmr.dmrupexpdown.keep$status=="DMRs-up"),], 
       aes(x=Gene, y=TF, fill=Adi.Comp1_C1_vs_C2.gz))+
  geom_tile()+
  scale_fill_gradient2(low="#AD343E", high="#1D3557", mid="white",
                       midpoint=0,na.value="black")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle =90, hjust=1))
#theme_classic()
dev.off()

# The same but excluding PRC2 targets:
sigs.plusdmr.dmrupexpdown.keep2 <- sigs.plusdmr.dmrupexpdown.keep[which(sigs.plusdmr.dmrupexpdown.keep$TF %in% c("SUZ12","EZH2")),]

pdf("EZH_SUZ12_targets.updated.dec2023/EZH2.SUZ12_DMRupEXPRdown_Targets_TS_DDR.pdf",w=8,h=2)
ggplot(sigs.plusdmr.dmrupexpdown.keep2[which(sigs.plusdmr.dmrupexpdown.keep2$status=="DMRs-up"),], 
       aes(x=Gene, y=TF, fill=Adi.Comp1_C1_vs_C2.gz))+
  geom_tile()+
  scale_fill_gradient2(low="#AD343E", high="#1D3557", mid="white",
                       midpoint=0,na.value="black")+
  xlab("")+
  ylab("")+
  theme(axis.text.x = element_text(angle =90, hjust=1))
#theme_classic()
dev.off()
