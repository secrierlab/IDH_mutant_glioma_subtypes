library(biomaRt)
library(sva)
library(edgeR)
library(pheatmap)
library(sva)
library(glmnet)
library(ComplexHeatmap)
library(data.table)
library(circlize)
library(stringr)
library(stringi)
library(gdata)


source("convert_mouse_to_human.R")


# Upload the RNA-seq values from Human
load('../../tme/tcga.lgg.mut.RData')
load('../../tme/tcga.gbm.mut.RData')

input_ge<-rbind(tcga.gbm,tcga.lgg)

col_to_select<-c(colnames(input_ge[,c(1:6)]),'Type')

#load annotation data
load("data/annotation_gbm_lgg.RData")
#load annotation 
load("data/CGGA_ann.RData")
#load microarray for metanalysis
load("data/database_microarray_brain.RData")

### Read in DEGs:
genes <- read.csv("../../DESeq2_annotMar2022_v2_summary/Adi.DE.signature.csv")
c1c2 <- unique(setdiff(genes$Adi.Comp1_C1_vs_C2_DEG,""))

ortholog_table_ensembl<-unique(convertMouseGeneList2(c1c2))
hs_genes_sig<-unique(ortholog_table_ensembl[,2])
save(hs_genes_sig, file="data/hs_genes_sig.RData")
  	
#
  	# Create matrix h.sapiens --NOT-- processed with combat
	#
	col_to_select<-c(colnames(input_ge[,c(1:6)]),'Type')
  matrix_to_plot_hs_nocombat<-t(input_ge[,-c(which(colnames(input_ge)%in%col_to_select))])
	colnames(matrix_to_plot_hs_nocombat)<-unlist(lapply(X=strsplit(input_ge$Sample,split='\\-'),FUN=function(X)paste(X[1:3],collapse='-')))

	matrix_to_plot_hs_nocombat2<-matrix_to_plot_hs_nocombat[which(rownames(matrix_to_plot_hs_nocombat)%in%hs_genes_sig),]
	
	source("support_functions_heatmaps.R")
        source("support_functions_MutBurden_GE_METH.R")
        source("create_matrix_meth.R")
        source("support_function_heatmap_mouse.R")

	###############################################
	# Plot human results
	################################################ 
	

	annotation_all<-ann3[,c("Samples","type","IDH1_R132H","TP53","ATRX","PDGFRA","IDH_codel_subtype","chr7_gain_chr_10loss","histology","clusters")]

	## Next, load TCGA C1AB annot and add here:
	c1abannot <- read.xls("../../tme/TCGA_GBM_LGG_GeneSig_median_cBio_wCNV_r1234.IDHmut.annot.plusTME.xlsx")
	
	annotation_all <- merge(annotation_all, 
	                        c1abannot[,c("patient","Mar2022.Comp8_C1A_vs_C1B_DN","C1AB")], 
	                        by.x="Samples", by.y="patient",
	                        all.x=TRUE, all.y=FALSE)
	
	# No Combat: all cancers
	output_hs<-'heatmapTCGA.nocombat.plusC1AB.pdf'
  	output_hs_var<-'heatmapTCGA.nocombat.var.plusC1AB.pdf'

	plot_heatmap_hs(matrix_to_plot_hs_nocombat2,annotation_all,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=3, png_output=F)	

	onlylgg <- annotation_all[which(annotation_all$type == "LGG"),]$Samples
	matrix_to_plot_hs_nocombat2.lgg<-matrix_to_plot_hs_nocombat2[,which(colnames(matrix_to_plot_hs_nocombat2) %in% onlylgg)]
	annotation_all.lgg <- annotation_all[which(annotation_all$Samples %in% colnames(matrix_to_plot_hs_nocombat2.lgg)),]
	output_hs<-'heatmapTCGA.nocombat.plusC1AB.onlyLGG.pdf'
	plot_heatmap_hs(matrix_to_plot_hs_nocombat2.lgg,annotation_all.lgg,output_hs=output_hs,output_hs2=output_hs_var,column_title="DEGs mm in GBM and LGG,no combat",width=10,split=TRUE,km=TRUE,number_clusters_row=2, png_output=F)	
	
	
		library(ConsensusClusterPlus)
	results = ConsensusClusterPlus(matrix_to_plot_hs_nocombat2,
				       maxK=4,
				       reps=50,
				       pItem=0.8,
				       pFeature=1,
				       title="test",
				       clusterAlg="km",
				       seed=1262118388.71279,plot="pdf")

	classesdf<-data.frame(ID=names(results[[3]]$consensusClass),
			      classes=results[[3]]$consensusClass)
	
	library(glmnet)
	cvfit = cv.glmnet(t(matrix_to_plot_hs_nocombat2),classesdf$classes,type.measure="class",alpha=1,family="multinomial")
	test<-coef(cvfit, s = "lambda.min") 

	resLasso<-function(x){
	idx<-which(x!=0)
	cg<-names(x[idx,])[-1]
	return(cg)
	}

	signatureGenes<-unique(unlist(lapply(test,resLasso)))

	output_hs<-'method.fc.combatsig.3cl.pdf'
	output_hs_var<-'method.fc.nocombat.var.test.pdf'
		       
  plot_heatmap_hs(matrix_to_plot_hs_nocombat2[rownames(matrix_to_plot_hs_nocombat2)%in%signatureGenes,],annotation_all,output_hs=output_hs,output_hs2=output_hs_var,png_output=F,column_title="DEGs mm in GBM and LGG,no combat",width=10,km=F,split=T,how_split=classesdf$classes,number_clusters_row=3)

  plot_heatmap_hs(matrix_to_plot_hs_nocombat2[rownames(matrix_to_plot_hs_nocombat2)%in%signatureGenes,],
                  annotation_all,output_hs=output_hs,output_hs2=output_hs_var,
                  column_title="DEGs mm in GBM and LGG,combat",width=10,split=TRUE,
                  km=T,number_clusters_row=3, png_output=F)	
  

### chose the 3 cluster solution as the best
