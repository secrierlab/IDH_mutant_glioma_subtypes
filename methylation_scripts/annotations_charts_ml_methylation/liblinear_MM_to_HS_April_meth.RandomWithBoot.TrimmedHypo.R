library("biomaRt")
library("data.table")
library(cowplot)
library(patchwork)
require(caret)
require(LiblineaR)
source("libLinear_mouse_hsApril.R")
library(plyr)
library(dplyr)
library(ggplot2)
library(factoextra)
library(ggpubr)

collapseMethValues<-function(ml_results,method=c("75","MLS","mean")){
  
  if(method=="75"){
    
    ml_temp<-ml_results[[2]]
    ml_mat1<-melt(ml_temp[,c(1:9)])
    colnames(ml_mat1)[3]<-"prob"
    
    ml_mat2<-unique(melt(ml_temp[,-c(2:9)]))
    custom_fun_quantile<-function(X){quantile(X)[4]}
    ml_mat2_agg<-aggregate(.~ID, ml_mat2[,-2], custom_fun_quantile)
    colnames(ml_mat2_agg)[2]<-"meth_levels"
    
    ml_final<-merge(ml_mat1,ml_mat2_agg)
    
  }
  
  if(method=="MLS"){
    
    ml_temp<-ml_results[[2]]
    ml_mat1<-melt(ml_temp[,c(1:9)])
    colnames(ml_mat1)[3]<-"prob"
    
    ml_mat2<-unique(melt(ml_temp[,-c(2:9)]))
    ml_mat2_agg<-aggregate(.~ID, ml_mat2[,-2], sum)
    colnames(ml_mat2_agg)[2]<-"meth_levels"
    
    ml_final<-merge(ml_mat1,ml_mat2_agg)
    
  }
  
  if(method=="mean"){
    
    ml_temp<-ml_results[[2]]
    ml_mat1<-melt(ml_temp[,c(1:9)])
    colnames(ml_mat1)[3]<-"prob"
    
    ml_mat2<-unique(melt(ml_temp[,-c(2:9)]))
    custom_fun_quantile<-function(X){sum}
    ml_mat2_agg<-aggregate(.~ID, ml_mat2[,-2], mean)
    colnames(ml_mat2_agg)[2]<-"meth_levels"
    
    ml_final<-merge(ml_mat1,ml_mat2_agg)
    
  }
  
  return(ml_final)
}

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
load("Mouse_beta_matrix_adi_pval0.05_lfc1.RData")

source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/clean_methylation_data.R")

#
# Read the DMRs regions
#
dmrMouse<-read.delim(file="/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/annotation_dmrs_adi_pval0.05_lfc1.txt",stringsAsFactors = F)
dmrMouse_rid<-unique(dmrMouse[,c(1,2,3,19,17)])
dmrMouse_stats<-dmrMouse[,c(19,13,14,17,18)]


#
# Convert mouse genes in human genes
#
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl",ssl.verifypeer = FALSE)
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",ssl.verifypeer = FALSE)
# httr::set_config(httr::config(ssl_verifypeer = FALSE))

# human = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
# mouse = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")

dmrs_up_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-up",19])
hs_dmrs_up_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_up_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

dmrs_down_mm<-unique(dmrMouse[dmrMouse$status=="DMRs-down",19])
hs_dmrs_down_from_mm = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = dmrs_down_mm , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

#
# Prepare human matrices
#

#look the script: Methylation_barplot_mouse_human.R for pre-processing methylation data
setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")
# see script: Methylation_barplot_mouse_human.random.R to see how I have created the RData
load("GBM_Mouse_beta_matrix_adi_pval0.05_lfc1.RANDOMGENES.RData")
load("LGG_Mouse_beta_matrix_adi_pval0.05_lfc1.RANDOMGENES.RData")

#
# Download the regions of the genes for each genes: is a variation of the script without the randomization, to make fast the analysis
#
#httr::set_config(httr::config(ssl_verifypeer = FALSE))         

ensembl = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
tssAll <- getBM(attributes = c("transcription_start_site", "chromosome_name","transcript_start", "transcript_end", "strand",  "ensembl_gene_id", "ensembl_transcript_id", "external_gene_name"),filters = "external_gene_name", values = unique(c(hs_meth_in_mouse_gbm[,4],hs_meth_in_mouse_lgg[,4])),mart = ensembl)

genes_dmr_mouse<-c(hs_dmrs_up_from_mm[,2],hs_dmrs_down_from_mm[,2])
ngenes_dmr_mouse<-length(genes_dmr_mouse)

LGG_ALL_STATS<-data.frame()
GBM_LGG_ALL_STATS<-data.frame()

for(ngenes_dmr_mouse in c(length(genes_dmr_mouse),round(length(genes_dmr_mouse)/seq(2,60,by=2)))){

for(i in 1:100){
	
	print(i)

	genes_to_use_all<-sample(setdiff(x=c(hs_meth_in_mouse_lgg$genes.Gene_Symbol,hs_meth_in_mouse_gbm$genes.Gene_Symbol),y=genes_dmr_mouse),ngenes_dmr_mouse)
	
	tss_to_use<-tssAll[which(tssAll$external_gene_name%in%genes_to_use_all),]
	hs_meth_in_mouse_gbm_filt<-clean_methylation_data_for_rand(mat_meth=hs_meth_in_mouse_gbm,tss=tss_to_use)

	hs_meth_in_mouse_lgg_filt<-clean_methylation_data_for_rand(mat_meth=hs_meth_in_mouse_lgg,tss=tss_to_use)

	setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human")
	load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/annotation_gbm_lgg_July2021.RData")

	ann_tot<-ann3

	#  Step 4: Combine all methylation
	hs_meth_in_mouse_gbm_filt2<-hs_meth_in_mouse_gbm_filt[which(hs_meth_in_mouse_gbm_filt[,4]%in%hs_meth_in_mouse_lgg_filt[,4]),]
  	hs_meth_in_mouse_lgg_filt2<-hs_meth_in_mouse_lgg_filt[which(hs_meth_in_mouse_lgg_filt[,4]%in%hs_meth_in_mouse_gbm_filt[,4]),]


	combined_meth<-combineMethMatricesRandomAnalysis(lgg_matrix=hs_meth_in_mouse_lgg_filt2,gbm_matrix=hs_meth_in_mouse_gbm_filt2)
	colnames(hs_meth_in_mouse_lgg_filt)<-gsub(unlist(lapply(strsplit(colnames(hs_meth_in_mouse_lgg_filt),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")
	colnames(hs_meth_in_mouse_gbm_filt)<-gsub(unlist(lapply(strsplit(colnames(hs_meth_in_mouse_gbm_filt),split="\\."),FUN=function(x){paste(x[1:3],collapse="-")})),pattern="\\.",replacement="-")
 	colnames(combined_meth)<-unlist(lapply(strsplit(colnames(combined_meth),split="\\-"),FUN=function(x){paste(x[1:3],collapse="-")}))

	#
	# Approach similar to the paper between dogs and human
	#


	setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/")

	#Measure the average Beta values for the genes in mus musculus
	meth_data_mouse<-setDT(res_meth3)
	meth_data_mouse2 <- data.frame(meth_data_mouse[, lapply(.SD, mean), by = gene_symbol])

	#httr::set_config(httr::config(ssl_verifypeer = FALSE))

	# human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
	#human = useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
	# mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
	#mouse = useEnsembl(biomart = "ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")

	#genes_mm_map_hs = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = meth_data_mouse2[,1] , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)

	#Measure the average Beta values for the genes in homo sapiens
	colnames(combined_meth)[1]<-"gene_symbol"
	meth_data_hs<-setDT(combined_meth)
	meth_data_hs2 <- data.frame(meth_data_hs[, lapply(.SD, mean), by = gene_symbol])
	colnames(meth_data_hs2)<-gsub(colnames(meth_data_hs2),pattern="\\.",replacement="-")

	#common_genes_hs_mm<-intersect(meth_data_hs2[,1],genes_mm_map_hs[,2])

	#mouse_genes_to_consider<-genes_mm_map_hs[,1]
	#human_genes_to_consider<-genes_mm_map_hs[,2]

	#meth_data_mouse3<-merge(genes_mm_map_hs,meth_data_mouse2[meth_data_mouse2[,1]%in%mouse_genes_to_consider,],by.x="MGI.symbol",by.y="gene_symbol")
	#meth_data_hs3<-meth_data_hs2[meth_data_hs2[,1]%in%human_genes_to_consider,]

	# Now merge the two data.frame because is the only way to get the same number of columns
	# Warning: i will not use more the mouse gene symbol but only the human gene symbol to be consistent

	#merged_matrix<-merge(meth_data_mouse3,meth_data_hs3,by.x="HGNC.symbol",by.y="gene_symbol")

	human_tab<-meth_data_hs2

	input_for_ml<-t(human_tab[,-1])
	colnames(input_for_ml)<-human_tab[,1]
	input_for_ml2<-data.frame(ID=rownames(input_for_ml),input_for_ml)
	input_for_ml3<-merge(ann_tot,input_for_ml2,by.x="Case",by.y="ID",all.x=F,all.y=F)
	colnames(input_for_ml3)<-gsub(colnames(input_for_ml3),pattern="\\.",replacement="-")

	mouse_tab<-meth_data_mouse2 #remove columns with gene-symbol and control EC

	input_mm<-t(mouse_tab[,-c(1,2)]) #remove columns with gene-symbol and control EC
	colnames(input_mm)<-mouse_tab[,1]
	input_mm2<-data.frame(ID=c("mut","mut","mut","wt","wt","wt","wt"),input_mm)
	colnames(input_mm2)[2:ncol(input_mm2)]<-gsub(colnames(input_mm2)[2:ncol(input_mm2)],pattern="\\.",replacement="-")
	colnames(input_mm2)[1]<-"status_mouse"
	rownames(input_mm2)<-colnames(mouse_tab)[-c(1,2)]

	list_ann<-c("DNA_methylation_clusters","IDH_TOT")

	#
	# Update the annotation with the most recent version
	#

	source("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September/libLinear_mouse_hsApril.RANDOMGENES.R")
	#load("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/compare_mouse_signatures_with_human/annotation_gbm_lgg_July2021.RData")
	#colnames(ann3)[1]<-"Samples"

	#input_ml4<-input_for_ml3[,-c(2:19)]
	#input_ml5<-merge(ann3,input_ml4,by.x="Samples",by.y="Case")
	ann3$IDH_TOT<-rep("0",nrow(ann3))
	ann3[which(ann3$IDH1==1 | ann3$IDH2==1),"IDH_TOT"]<-1
	input_ml5<-merge(ann3[,-1],input_for_ml3[,-c(1,3:19)],by.x="SAMPLE_ID",by.y="SAMPLE_ID")

	input_ml5_new<-input_ml5[sample(nrow(input_ml5),replace=T),]

	temp_meth<-input_ml5_new[,c(20:ncol(input_ml5_new))]
	#temp_meth[temp_meth>0.5]<-
	
	#custom_trim_function<-function(X2){
	#	X2[which(X2<=0.5)]<-0
	#	return(X2)
	#}
	
	#temp_meth<-apply(temp_meth,2,custom_trim_function)
	temp_meth[temp_meth<=0.5]<-0
	temp_meth[is.na(temp_meth)]<-0

	input_ml5_new2<-data.frame(input_ml5_new[,c(1:19)],temp_meth)

	lgg_ml_res<-libLinear_mouse_hsAprilRANDOM(input_ml5_new2,input_mm2,list_ann=list_ann,type_tumours="LGG",prefix_output="LGG",suffix_output="DMRs_hyper_hypo")
	gbm_lgg_ml_res<-libLinear_mouse_hsAprilRANDOM(input_ml5_new2,input_mm2,list_ann=list_ann,type_tumours=c("GBM","LGG"),prefix_output="GBM_LGG",suffix_output="DMRs_hyper_hypo")

	lgg_acc<-data.frame(ngenes=ngenes_dmr_mouse,niter=as.character(i),t(data.frame(t(data.frame(lgg_ml_res[[3]]))[1,])))
	gbm_lgg_acc<-data.frame(ngenes=ngenes_dmr_mouse,niter=as.character(i),t(data.frame(t(data.frame(gbm_lgg_ml_res[[3]]))[1,])))

	if(i==1){
	  
	input_ml5_new_temp<-input_ml5_new2[input_ml5_new2$type%in%"LGG",]
	  
	input_pca<-input_ml5_new_temp[,19:ncol(input_ml5_new_temp)]
	input_pca[is.na(input_pca)]<-0
	res.pca <- prcomp(input_pca[,-1], scale = F)
	
	pdf(paste(paste("PCA_BOOT_RANDOMGENES_TRIMMED_iter",i,ngenes_dmr_mouse,sep="."),".pdf",sep=""))
	p<-fviz_pca_ind(res.pca,
	                col.ind = input_pca[,1], # color by groups
	                palette = c("#00AFBB",  "#FC4E07"),
	                legend.title = "Groups",
	                label="none",
	                pointsize = 2.5,
	                title=paste("LGG",i,ngenes_dmr_mouse))
	print(p)
	
	input_boxplot<-melt(input_pca)
	
	p2<-ggboxplot(input_boxplot, "IDH_TOT", "value",
	          fill = "IDH_TOT", palette = c("#E7B800","#00AFBB"))+ stat_compare_means(method = "wilcox.test")
	print(p2)

	p3<-ggboxplot(input_boxplot[which(input_boxplot[,3]!=0),], "IDH_TOT", "value",
		                        fill = "IDH_TOT", palette = c("#E7B800","#00AFBB"))+ stat_compare_means(method = "wilcox.test")
        print(p3)

	dev.off()
	
	}
	
	LGG_ALL_STATS<-rbind(LGG_ALL_STATS,lgg_acc)
	GBM_LGG_ALL_STATS<-rbind(GBM_LGG_ALL_STATS,gbm_lgg_acc)

setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")

}

}

LGG_ALL_STATS[,1]<-as.character(LGG_ALL_STATS[,1])
GBM_LGG_ALL_STATS[,1]<-as.character(GBM_LGG_ALL_STATS[,1])

lgg_melt<-melt(LGG_ALL_STATS)
lgg_melt$ngenes<-factor(lgg_melt$ngenes,levels=rev(c(length(genes_dmr_mouse),round(length(genes_dmr_mouse)/seq(2,60,by=2)))))

gbm_lgg_melt<-melt(GBM_LGG_ALL_STATS)
gbm_lgg_melt$ngenes<-factor(gbm_lgg_melt$ngenes,levels=rev(c(length(genes_dmr_mouse),round(length(genes_dmr_mouse)/seq(2,60,by=2)))))

save(list=c("lgg_melt","gbm_lgg_melt"),file="liblinear_randomgenes_methylationWithBootTrimmed.RData")

pdf("boxplot_stats_liblinear_randomgenesWithBootTrimmed.pdf")

# Stat Summary contains the accuracy obtained using the original datasets (see liblinear_MM_to_HS_April_meth.R)

p<-ggplot(lgg_melt, aes(x=variable, y=value, fill=variable)) + 
	geom_boxplot()+geom_hline(yintercept=c(0.79,0.98), linetype="dashed",color = "red", size=1)+
	stat_summary(geom = 'text', label = c(0.79,0.98), fun.y = max, vjust = -1)+ylab("accuracy")+ggtitle("Accuracy 10000 iterations LGG")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

p2<-ggplot(gbm_lgg_melt, aes(x=variable, y=value, fill=variable)) +
        geom_boxplot()+geom_hline(yintercept=c(0.79,0.97), linetype="dashed",color = "red", size=1)+
	stat_summary(geom = 'text', label = c(0.79,0.97), fun.y = max, vjust = -1)+ylab("accuracy")+ggtitle("Accuracy 10000 iterations GBM+LGG")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

print(plot_grid(p,p2))

lgg_melt$ngenes<-as.factor(lgg_melt$ngenes)

prand<-ggplot(lgg_melt, aes(x=variable, y=value, fill=ngenes)) +
        geom_boxplot()+ylab("accuracy")+ggtitle("Accuracy 100 iterations LGG")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

prand2<-ggplot(gbm_lgg_melt, aes(x=variable, y=value, fill=ngenes)) +
        geom_boxplot()+ylab("accuracy")+ggtitle("Accuracy 100 iterations GBM+LGG")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

# ph<-ggplot(lgg_melt, aes(niter, ngenes, fill= value)) +
#       geom_tile()+scale_fill_gradient2(low = "white", mid = "orange2", high = "blue2", midpoint = .25)+ggtitle("Accuracy 10000 iterations LGG")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+facet_wrap(~variable,nrow=2)

# ph2<-ggplot(gbm_lgg_melt, aes(niter, ngenes, fill= value)) +
#  geom_tile()+scale_fill_gradient2(low = "white", mid = "orange2", high = "blue2", midpoint = .25)+ggtitle("Accuracy 10000 iterations GBM+ LGG")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+theme_bw()+facet_wrap(~variable,nrow=2)

print(prand)
print(prand2)

dev.off()





# 
# setwd("/mnt/data/lab/gmt_data/data_brain/gmt_analysis/meth_analysis/analysis_September")
# 
# list_res_ml_to_use<-c("lgg_ml_res","gbm_lgg_ml_res")
# names(list_res_ml_to_use)<-c("LGG","GBM_LGG")
# 
# for(lrmt in 1:length(list_res_ml_to_use)){
#  	 
#     current_ml_res<-names(list_res_ml_to_use)[lrmt]
#     input_ml_to_pca<-get(list_res_ml_to_use[[lrmt]])[[2]][,c(1:9)]
# 
#     brain_melt_ml<-melt(input_ml_to_pca)
#     brain_melt_ml[,2]<-gsub(brain_melt_ml[,2],pattern="variable_",replacement="")
#     colnames(brain_melt_ml)[3]<-"prob"
#     brain_melt_ml[,1]<-gsub(brain_melt_ml[,1],pattern="\\.",replacement="-")
#     
#     input_for_pca<-merge(brain_melt_ml,input_ml5,by.x="ID",by.y="Samples")
#     ann_pca<-input_for_pca[,c(1:21)]
#     data_pca<-input_for_pca[,-c(1:21)]
#     
#     data_pca[is.na(data_pca)]<-0
#     res.pca <- prcomp(data_pca,  scale = TRUE)
#     res.pca2<- data.frame(res.pca$x[,c(1:2)])
#     res.pca2$samples<-input_for_pca[,1]
#     
#     res.pca3<-merge(res.pca2,ann_pca,by.x="samples",by.y="ID")
#     
#     pdf(paste(current_ml_res,"PCA_mm_to_hs_prob_meth.pdf",sep="_"),width=16)
#     
#     res.pca3$DNA_methylation_clusters<-gsub(res.pca3$DNA_methylation_clusters,pattern="no_ann",replacement="LGr_Unkown")
#     
#     p <- ggplot(data = res.pca3, aes(x = PC1, y = PC2,colour=prob)) + geom_point() + scale_colour_gradient2(low = "gainsboro", mid = "orange2", high = "red2", midpoint = 0.5) + geom_hline(yintercept=0, linetype="dashed", color = "black") + geom_vline(xintercept = 0, linetype="dashed", color = "black")+ggtitle(paste("n.samples:",length(unique(res.pca3[,1]))))
#     
#     p2 <- ggplot(data = res.pca3, aes(x = PC1, y = PC2,colour=IDH1)) + geom_point(aes(colours=IDH1))+scale_color_manual(breaks = c("0", "1"),values=c("green4", "gold2")) + geom_hline(yintercept=0, linetype="dashed", color = "black") + geom_vline(xintercept = 0, linetype="dashed", color = "black")+ggtitle(paste("n.samples:",length(unique(res.pca3[,1]))))
#     
#     p3 <- ggplot(data = res.pca3, aes(x = PC1, y = PC2,colour=prob)) + geom_point() + scale_colour_gradient2(low = "gainsboro", mid = "orange2", high = "red2", midpoint = 0.5) + geom_hline(yintercept=0, linetype="dashed", color = "black") + geom_vline(xintercept = 0, linetype="dashed", color = "black")+ggtitle(paste("n.samples:",length(unique(res.pca3[,1]))))
#     
#     ngroups<-length(table(res.pca3$DNA_methylation_clusters))
#     p4 <- ggplot(data = res.pca3, aes(x = PC1, y = PC2,colour=DNA_methylation_clusters)) + geom_point(aes(colours=DNA_methylation_clusters))+scale_shape_manual(values = c(15,16,17,18,19,0,1))+ scale_color_manual(breaks = c("LGm1", "LGm2", "LGm3","LGm4","LGm5","LGm6","LGr_Unkown"),values=c("brown2", "darkorchid3", "dodgerblue4","olivedrab3","olivedrab3","olivedrab3","black")) + geom_hline(yintercept=0, linetype="dashed", color = "black") + geom_vline(xintercept = 0, linetype="dashed", color = "black")+ggtitle(paste("n.samples:",length(unique(res.pca3[,1]))))
#     
#     p5 <- ggplot(data = res.pca3, aes(x = PC1, y = PC2,colour=prob)) + geom_point() + scale_colour_gradient2(low = "gainsboro", mid = "orange2", high = "red2", midpoint = 0.5) + geom_hline(yintercept=0, linetype="dashed", color = "black") + geom_vline(xintercept = 0, linetype="dashed", color = "black")+ggtitle(paste("n.samples:",length(unique(res.pca3[,1]))))
#     
#     p6 <- ggplot(data = res.pca3, aes(x = PC1, y = PC2,colour=type)) + geom_point(aes(colours=type))+scale_color_manual(breaks = c("LGG", "GBM"),values=c("skyblue1", "red2")) + geom_hline(yintercept=0, linetype="dashed", color = "black") + geom_vline(xintercept = 0, linetype="dashed", color = "black")+ggtitle(paste("n.samples:",length(unique(res.pca3[,1]))))+ggtitle(paste("n.samples:",length(unique(res.pca3[,1]))))
#     
#     
#     print(p+p2)
#     print(p3+p4)
#     print(p5+p6)
#     
#     dev.off()
# 
# }
# 
# library(ggpubr)
# 
# pdf("LGG_boxplot_mm_to_hs_prob_meth.pdf")
# 
# lgg_75<-collapseMethValues(lgg_ml_res,method="75")
# lgg_75$variable<-gsub(lgg_75$variable,pattern="variable_",replacement="")
# lgg_75$variable<-factor(lgg_75$variable,levels=c("X0","X1","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6"))
# 
# lgg_MLS<-collapseMethValues(lgg_ml_res,method="MLS")
# lgg_MLS$variable<-gsub(lgg_MLS$variable,pattern="variable_",replacement="")
# lgg_MLS$variable<-factor(lgg_MLS$variable,levels=c("X0","X1","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6"))
# 
# lgg_mean<-collapseMethValues(lgg_ml_res,method="mean")
# lgg_mean$variable<-gsub(lgg_mean$variable,pattern="variable_",replacement="")
# lgg_mean$variable<-factor(lgg_mean$variable,levels=c("X0","X1","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6"))
# 
# list_lgg<-vector(mode="list",3)
# 
# list_lgg[[1]]<-lgg_75
# list_lgg[[2]]<-lgg_MLS
# list_lgg[[3]]<-lgg_mean
# names(list_lgg)<-c("75quantile","MLS","mean")
# 
# nsamples<-paste("N.patients",paste(lgg_ml_res[[3]],collapse="-"),sep=": ")
# 
# for(i in 1:length(list_lgg)){
# 
# current_lgg<-list_lgg[[i]]
# current_measure<-names(list_lgg)[i]
#   
# current_lgg<-current_lgg[-which(current_lgg[,2]%in%"X0"),]
# 
# p3<-ggplot(current_lgg, aes(x=variable, y=prob))+geom_boxplot()+ggtitle(nsamples)+ylab("Probability")+
#   geom_jitter(shape=16, position=position_jitter(0.2),aes(color = meth_levels))+scale_colour_gradient(low="blue",high="red")+theme_bw() 
# 
# print(p3)
# 
# cont_tables<-split(lgg_ml_res[[5]],lgg_ml_res[[5]]$comparison)
# 
# cont_table_DNAmeth<-cont_tables[[1]]
# cont_table_IDH<-cont_tables[[2]]
# 
# pct1<-ggplot(cont_table_IDH, aes(predictions, real)) +
#   geom_tile(aes(fill = value)) + 
#   geom_text(aes(label = round(value, 1))) +
#   scale_fill_gradient(low = "white", high = "red")+ggtitle(paste("IDH model accuracy:","\n",round(lgg_ml_res[[4]]["IDH1",1],3)))
# 
# pct2<-ggplot(cont_table_DNAmeth, aes(predictions, real)) +
#   geom_tile(aes(fill = value)) + 
#   geom_text(aes(label = round(value, 1))) +
#   scale_fill_gradient(low = "white", high = "red")+ggtitle(paste("DNA methylation accuracy:","\n",round(lgg_ml_res[[4]]["DNA_methylation_clusters",1],3)))
# 
# print(plot_grid(pct1,pct2))
# 
# psm<-ggplot(current_lgg, aes(x=meth_levels, y=prob))+geom_point()+geom_smooth(method = "loess", se = TRUE)+facet_wrap(~variable)+theme_bw()+ylab("Probability")+xlab(paste("Methylation levels","(",current_measure,")"))
# print(psm)
# 
# psm2<-ggplot(current_lgg, aes(x=meth_levels, y=prob))+geom_point()+geom_smooth(method = "loess", se = TRUE)+theme_bw()+ylab("Probability")+xlab(paste("Methylation levels","(",current_measure,")"))
# print(psm2)
# 
# hist_meth<-ggplot(current_lgg, aes(x=meth_levels)) + geom_histogram(bins=100)+xlab(paste("Methylation levels","(",current_measure,")"))
# print(hist_meth)
# 
# }
# 
# dev.off()
# 
# 
# # 
# # Plot results LGG 
# # 
# 
# pdf("GBM_LGG_boxplot_mm_to_hs_prob_meth.pdf")
# 
# gbm_lgg_75<-collapseMethValues(gbm_lgg_ml_res,method="75")
# gbm_lgg_75$variable<-gsub(gbm_lgg_75$variable,pattern="variable_",replacement="")
# gbm_lgg_75$variable<-factor(gbm_lgg_75$variable,levels=c("X0","X1","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6"))
# 
# gbm_lgg_MLS<-collapseMethValues(gbm_lgg_ml_res,method="MLS")
# gbm_lgg_MLS$variable<-gsub(gbm_lgg_MLS$variable,pattern="variable_",replacement="")
# gbm_lgg_MLS$variable<-factor(gbm_lgg_MLS$variable,levels=c("X0","X1","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6"))
# 
# gbm_lgg_mean<-collapseMethValues(gbm_lgg_ml_res,method="mean")
# gbm_lgg_mean$variable<-gsub(gbm_lgg_mean$variable,pattern="variable_",replacement="")
# gbm_lgg_mean$variable<-factor(gbm_lgg_mean$variable,levels=c("X0","X1","LGm1","LGm2","LGm3","LGm4","LGm5","LGm6"))
# 
# list_gbm_lgg<-vector(mode="list",3)
# 
# list_gbm_lgg[[1]]<-gbm_lgg_75
# list_gbm_lgg[[2]]<-gbm_lgg_MLS
# list_gbm_lgg[[3]]<-gbm_lgg_mean
# names(list_gbm_lgg)<-c("75quantile","MLS","mean")
# 
# nsamples<-paste("N.patients",paste(gbm_lgg_ml_res[[3]],collapse="-"),sep=": ")
# 
# for(i in 1:length(list_gbm_lgg)){
#   
#   current_gbm_lgg<-list_gbm_lgg[[i]]
#   current_measure<-names(list_gbm_lgg)[i]
#   
#   current_gbm_lgg<-current_gbm_lgg[-which(current_gbm_lgg[,2]%in%"X0"),]
#   
#   p3<-ggplot(current_gbm_lgg, aes(x=variable, y=prob))+geom_boxplot()+ggtitle(nsamples)+ylab("Probability")+
#     geom_jitter(shape=16, position=position_jitter(0.2),aes(color = meth_levels))+scale_colour_gradient(low="blue",high="red")+theme_bw() 
#   
#   print(p3)
#   
#   cont_tables<-split(gbm_lgg_ml_res[[5]],gbm_lgg_ml_res[[5]]$comparison)
#   
#   cont_table_DNAmeth<-cont_tables[[1]]
#   cont_table_IDH<-cont_tables[[2]]
#   
#   pct1<-ggplot(cont_table_IDH, aes(predictions, real)) +
#     geom_tile(aes(fill = value)) + 
#     geom_text(aes(label = round(value, 1))) +
#     scale_fill_gradient(low = "white", high = "red")+ggtitle(paste("IDH model accuracy:","\n",round(gbm_lgg_ml_res[[4]]["IDH1",1],3)))
#   
#   pct2<-ggplot(cont_table_DNAmeth, aes(predictions, real)) +
#     geom_tile(aes(fill = value)) + 
#     geom_text(aes(label = round(value, 1))) +
#     scale_fill_gradient(low = "white", high = "red")+ggtitle(paste("DNA methylation accuracy:","\n",round(gbm_lgg_ml_res[[4]]["DNA_methylation_clusters",1],3)))
#   
#   print(plot_grid(pct1,pct2))
#   
#   psm<-ggplot(current_gbm_lgg, aes(x=meth_levels, y=prob))+geom_point()+geom_smooth(method = "loess", se = TRUE)+facet_wrap(~variable)+theme_bw()+ylab("Probability")+xlab(paste("Methylation levels","(",current_measure,")"))
#   print(psm)
#   
#   psm2<-ggplot(current_gbm_lgg, aes(x=meth_levels, y=prob))+geom_point()+geom_smooth(method = "loess", se = TRUE)+theme_bw()+ylab("Probability")+xlab(paste("Methylation levels","(",current_measure,")"))
#   print(psm2)
#   
#   hist_meth<-ggplot(current_gbm_lgg, aes(x=meth_levels)) + geom_histogram(bins=100)+xlab(paste("Methylation levels","(",current_measure,")"))
#   print(hist_meth)
#   
# }
# 
# dev.off()
# 
# 
