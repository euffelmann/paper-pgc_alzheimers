# =============================================================
#
#  magma_gs.sh
#  Perform MAGMA gene set and gene property analyses
#  Table of Contents:
#  1) Gene set analysis using MSIGDB C7+C2+C5 + SYNGO
#  2) Gene set clustering
#  3) FUMA cell type analysis
#  4) Microglia state analysis
#  5) AD case-control cell type analysis

#
# =============================================================


### 1) Gene set analysis using MSIGDB C7+C2+C5 + SYNGO
#set up
#download magma_v1.10 https://cncr.nl/research/magma/
#download MSigDB v2024.1.Hs.symbols C2, C5, C7  https://www.gsea-msigdb.org/gsea/index.jsp
#convert the gene symbols to ENSG using biomaRt to link ENSG to symbol
#remove gene-sets with fewer than 10 genes
#download SynGO data https://www.syngoportal.org/data/SynGO_bulk_download_release_20231201.zip
#remove gene-sets with fewer than 10 genes
#combine MSigDB and SynGO genesets

#run gene set analysis using MAGMA results from FLAMES
./magma --gene-results magma.genes.raw --set-annot  MSigDB2024_C257_ENSG_syngo_20231201_ontologies.gmt --out MSigDB_SynGO

#select FDR <0.05 gene sets in R
cat > FDR.R << END
library(data.table)
a<-fread("MSigDB_SynGO.gsa.out")
a$FDR<-p.adjust(a$P, method="BH")
a<-a[a$FDR<0.05,]
a<-a[order(a$P),]
fwrite(a, file="FDRsig.txt", col.names=T, row.names=F, sep=" ", na=NA, quote=F)
END
Rscript FDR.R
awk '{print $8}' FDRsig.txt | tail -n +2 > fdrgs.txt



#### 2) Gene set clustering

#run pairwise conditioning
cat fdrgs.txt | while read i
do 
./magma --gene-results magma.genes.raw --set-annot MSigDB2024_C257_ENSG_syngo_20231201_ontologies.gmt --out FDR_${i} --model analyse=file, fdrgs.txt condition-hide=${i}
done

#summarise into one file
echo FULL_NAME ConditionedOn TYPE NGENES BETA BETA_STD SE P > condFDR.txt
for i in FDR_*gsa.out
do
gs=$(echo $i | sed 's/FDR_//g' | sed 's/.gsa.out//g')
echo $gs
tail -n +6 $i | grep -v "$gs" | awk -v gse="$gs" '{print $8,gse,$2,$3,$4,$5,$6,$7}' >> condFDR.txt
done


#cluster based on proportion of pairwise explained significance
cat > cluster.R << END
library(data.table)
library(igraph)
library(ape)
library(svglite)

#read in conditional and marginal P-value
a<-fread("condFDR.txt")
#read in marginal P
b<-fread("MSigDB_SynGO.gsa.out")
b<-b[,c("FULL_NAME","P")]
colnames(b)[2]<-"marginal"
a<-merge(a,b,by="FULL_NAME",all.x=T)
#amount explained by conditioned on gs
a$ce<- -log10(a$marginal) - -log10(a$P)
#get proportion of marginal explained by conditioned on gs
a$prop<- a$ce/-log10(a$marginal)

#node size is marginal P
nodeinfo<-a[,c("FULL_NAME","marginal")]
#set edge width as proportion of significance explained by each gs
a<-a[,c("FULL_NAME","ConditionedOn","prop")]
colnames(a)<-c("to","from","prop")

#get node size as -log10(P) marginal
colnames(nodeinfo)<-c("id","set_size")
nodeinfo<-nodeinfo[!duplicated(nodeinfo),]
nodeinfo$set_size<- -log10(nodeinfo$set_size)

#only use edges where proportion of significance is >0.3
net <- graph_from_data_frame(d=a, vertices=nodeinfo, directed=F) 
net.sp <- delete_edges(net, E(net)[prop<0.3])
net.sp <- simplify(net.sp, remove.multiple = F, remove.loops = T) 
ceb <- cluster_edge_betweenness(net.sp) 

#plot clustering
svglite("ClusteredGS30.svg")
plot(ceb, net.sp, edge.width=E(net.sp)$prop, edge.label=round(E(net.sp)$prop,3), edge.label.color="black", edge.label.cex=.2, edge.arrow.size=.05, vertex.label.cex=.1, vertex.label.color="black",vertex.label.cex=.2, rescale=T, vertex.size=V(net.sp)$set_size) 
dev.off()
#save groups
g<-as.data.table(cbind(ceb$names,ceb$membership))
fwrite(g, file="cebgroups03.txt", col.names=F, row.names=F, na=NA, quote=F, sep=" ")
END
Rscript cluster.R


## 3) FUMA cell type analysis
#https://fuma.ctglab.nl/
#upload to FUMA with the following options:
#datasets
#59_Siletti_Hippocampus.HiT.CA4-DGC_Human_2022_level2
#60_Siletti_Hippocampus.HiH.HiT.Sub_Human_2022_level2
#61_Siletti_Hippocampus.HiH.CA1_Human_2022_level2
#62_Siletti_Hippocampus.HiH.CA1-3_Human_2022_level2
#63_Siletti_Hippocampus.HiH.CA1-CA3_Human_2022_level2
#64_Siletti_Hippocampus.HiH.DG-CA4_Human_2022_level2
#65_Siletti_Hippocampus.HiB-RostralCA1-CA3_Human_2022_level2
#66_Siletti_Hippocampus.HiH.CA2-3_Human_2022_level2
#67_Siletti_Hippocampus.HiB.RostralCA1-2_Human_2022_level2
#68_Siletti_Hippocampus.HiB.RostralDG-CA4_Human_2022_level2
#69_Siletti_Hippocampus.HiB.RostralCA3_Human_2022_level2
#Allen_Human_MTG_level2
#GSE168408_Human_Prefrontal_Cortex_level2_Adult
#PsychENCODE_Adult
#bonferroni correction
#step 2 and step 3




### 4) Microglia state analysis
#Download Sun 2023 microglia states from supplementary table 2
#https://www.sciencedirect.com/science/article/pii/S0092867423009716
#convert gene symbols to ENSG using using biomaRt to link ENSG to symbol

#MG6 - Stress signature
#MG0 - homeostatic
#MG8 - Inflammatory II
#MG10 - Inflammatory III
#MG7 - Glycolytic
#MG2 - Inflammatory I
#MG4 - Lipid processing
#MG11 - Antiviral
#MG1 - Neuronal surveillance
#MG5 - Phagocytic
#MG12 - Cycling
#MG3 - Ribosome biogenesis

# Run MAGMA analysis
./magma --gene-results magma.genes.raw --set-annot Sun2023_MicrogliaState.gmt --out SunMicroState





### 5) AD case-control cell type analysis

#Download AD differential expression values from Supplementary file 3 of Nakatsuka 2025
#https://www.biorxiv.org/content/10.1101/2024.10.15.618577v2
#convert gene symbols to ENSG using using biomaRt to link ENSG to symbol
#one file per differential direction, each column is a cell type, each row is a gene, cells= -log10(P) of differential expression for that gene in that cell type
#use -log10(P) as gene property for analysis
#calculate average differential expression for each gene (rowMean) for conditional analysis


#run MAGMA
for j in {UP,DOWN}
do
./magma --gene-results magma.genes.raw --gene-covar ${j}_PVal_ENSG.txt --model direction=greater --out ${j}_PVal
./magma --gene-results magma.genes.raw --gene-covar ${j}_PVal_ENSG.txt --model condition-hide=Average direction=greater --out ${j}_PVal_condition
done



