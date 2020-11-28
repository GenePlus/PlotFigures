library(ggplot2)
library(ggpubr)
library(ggsci)
library(pheatmap)
library(ComplexHeatmap)
library(Biostrings)
library(rtracklayer)
library(ComplexHeatmap)
library(reshape2)


#-----------------------
#     Figure 1a
#-----------------------
PCA=read.table("Data/Data_Fig2a_PCA.eigenvec",header=T)
pdf("Data/Fig2a_PCA.pdf",height=3.5,width=4.5)
ggplot(PCA, aes(x=PCA1,y=PCA2)) +
  geom_point(aes(color=Group1),size=0.7)+scale_color_d3()+
  labs(x="PC1",y="PC2")+theme_bw()+theme(legend.title = element_blank())
dev.off()

#-----------------------
#     Figure 1b
#-----------------------
raw_fst=read.table("Data/Data_Fig2b_fst.txt",header=T)
fst=as.matrix(raw_fst/1000)
pheatmap(fst,cluster_rows = F,cluster_cols = F, display_numbers = T,fontsize = 10,
	filename = "Data/Fig2b_fst.pdf", width = 4.5, height = 3.5)

#-----------------------
#     Figure 4a
#-----------------------
Clin<-read.csv("Data/Data_Fig4a_Clin.csv",header=T,stringsAsFactors =F,sep=",") 
Mut<-read.table("Data/Data_Fig4a_MAF.txt",header=T,stringsAsFactors =F,sep="\t",comment.char="") 
	
#--set TMB on top
qctmb<-as.data.frame(t(table(Mut[,c(1,4)])))#row=pat id; col=variant type
qctmb<-dcast(qctmb,Sample~Function,fill="0",drop=F)
rownames(qctmb)<-qctmb[,1]
qctmb<-qctmb[,-1]
TMBp<-na.omit(qctmb)
TMBp=as.matrix(TMBp)
storage.mode(TMBp)="numeric"
TMBp=as.data.frame(TMBp)

#--set the variant type
maf.dat1= Mut[,c("Gene", "Sample","Function")]
ma<-maf.dat1$Variant_Classification
maf.dat.mat1 = dcast(data = maf.dat1, Gene ~ Sample,
                                 fun.aggregate=function(x)paste(x,collapse=";"),
                                 value.var ="Function",fill="",drop = FALSE)
rownames(maf.dat.mat1 ) =maf.dat.mat1 [,1]
maf.dat.mat1 = maf.dat.mat1 [,-1]
mat<-maf.dat.mat1
mat[is.na(mat)]=""
gene_order<-apply(mat,1,function(x)as.character(length(which(x==""))))
gene_order<-gene_order[order(as.numeric(gene_order))]
mat<-mat[names(gene_order),]

#--set the grid
alter_fun = function(x, y, w, h, v) {
  n=sum(v)
  h=h*0.9
  grid.rect(x, y, w-unit(0.40,"mm"),h-unit(0.40,"mm"), gp = gpar(fill ="#CCCCCCDF",col=NA))
  if(v["CNV"])  grid.rect(x, y - h*0.5 + 0.95:n/n*h, w*0.85, 1/n*h,gp = gpar(fill = nocol[names(which(v))], col = NA),just = "top")
  if(v["frameshift"]) grid.rect(x, y - h*0.5 + 0.95:n/n*h, w*0.85, 1/n*h,gp = gpar(fill = nocol[names(which(v))], col = NA), just = "top")
  if(v["indel"])  grid.rect(x, y - h*0.5 + 0.95:n/n*h, w*0.85, 1/n*h,gp = gpar(fill = nocol[names(which(v))], col = NA), just = "top")
  if(v["missense"])  grid.rect(x, y - h*0.5 + 0.95:n/n*h, w*0.85, 1/n*h,gp = gpar(fill = nocol[names(which(v))], col = NA), just = "top")
  if(v["nonsense"])  grid.rect(x, y - h*0.5 + 0.95:n/n*h, w*0.85, 1/n*h,gp = gpar(fill = nocol[names(which(v))], col = NA), just = "top")
  if(v["splice"])  grid.rect(x, y - h*0.5 + 0.95:n/n*h, w*0.85, 1/n*h,gp = gpar(fill = nocol[names(which(v))], col = NA), just = "top")
  if(v["SV"])  grid.rect(x, y - h*0.5 + 0.95:n/n*h, w*0.85, 1/n*h,gp = gpar(fill = nocol[names(which(v))], col = NA), just = "top")
}


#--set the annotation
nocol<-c("CNV"="#9912A896","frameshift"="#E3A452D4","indel"="#E6738DCA","missense"="dodgerblue3","nonsense"="#C9BD14BE",
"others"="grey","splice"="indianred1","SV"="#228B22")
hat<-HeatmapAnnotation(TMB=anno_barplot(TMBp,gp=gpar(fill=nocol,
                                                    col="white"),border=F,bar_width=0.9),
						df=Clin[,c("Gender","Stage","pathology")],
                       col=list("Gender"=c("male"="steelblue1","female"="pink1"),
                                "Stage"=c("I_II"="#4EEE94","III"="#00CD00","IV"="#228B22" ,"NA"="grey"),
                                "pathology"=c("LUAD"="#D62FA79B","LUSC"= "#92C3E0DA","NSCLC-NOS"="#EEC591","other"="grey")),
                       gap=unit(0.55,"mm"),gp = gpar(col = "white"),
                       annotation_height = unit(c(2.0,0.3,0.3,0.3), c("cm","cm","cm","cm")),
                       show_annotation_name=T, 
                       annotation_name_side="left",
                       height = unit(4.1, "cm"))

colnol=c("CNV","frameshift","indel","missense","nonsense","splice","others","SV")


genelist=c("EGFR","TP53","LRP1B","KRAS","CDKN2A",	"CDH23",	"MET",
           "GRIN2A",	"AR",	"NTRK3",	"ATM",	"FAT1",	"NF1",	"APC",	"SMAD4",
           "JAK2"	"ALK",	"CDKN2B",	"CREBBP",	"EPHA5",	"FAT2",	"MLL3",	"NOTCH4",
           "PIK3CA",	"RB1",	"ROS1",	"TET2",	"NSD1",	"ATRX",	"MLL",	"CTNNB1",
           "EPHA3",	"ERBB2",	"FGFR3",	"CUL3",	"RUNX1",	"BRCA2",	"BRAF",	"RET" )
		   
mat1=mat[geneList,]
#mat1<-mat[1:30,] # visualize the top 30 rows
#--plot the figure
ht<-oncoPrint(mat1,get_type =function(x) strsplit(x, ";")[[1]],
              alter_fun = alter_fun,
              col=nocol,
              row_order = NULL, column_order = NULL,
			  show_column_names = TRUE,
			  row_names_gp = gpar(fontsize = 7),column_names_gp = gpar(fontsize = 7),
              show_pct = TRUE, right_annotation = rowAnnotation(rbar = anno_oncoprint_barplot()),
              remove_empty_columns=F,
              top_annotation=hat,
              heatmap_legend_param = list(title = "Alternations",at = colnol,  labels = colnol))
pdf("Data/Fig4a_heatmap.pdf",height=6,width=12)			  
draw(ht,heatmap_legend_side ="left",annotation_legend_side="left")
dev.off()
