library(data.table)
library(ggplot2)
library(patchwork)

#Read in supplementary data Nc
a<-fread("SuppDataNc.txt")
#remove trailing genes and replace with ...
a$Gene<-gsub("/.*","/..",a$Gene)

#novel
n<-a[a$Novel==1,]
n<-n[,-c("chromosome","start","end","size_bp","Novel")]
n[,3:22] <- lapply(n[,3:22], as.numeric)
nl = melt(n, id.vars = c("Gene","locus"),measure.vars = colnames(n)[-c(1,2)])
nl$p<--log10(nl$value)
nl<-nl[order(nl$locus),]
nl$id<-paste("Locus",nl$locus,nl$Gene,sep=" ")
nl$id<-factor(nl$id,levels=unique(nl$id))
nl$value<-signif(nl$value, digits=3)
nl$sig<-"none"
nl$sig[nl$value<1e-5]<-"sug"
nl$sig[nl$value<5e-8]<-"sig"
nl$value<-as.character(nl$value)
nl$value[is.na(nl$value)]<-"NA"

ggheatmap <- ggplot(nl, aes(variable, id, fill = sig))+
  scale_y_discrete(limits=rev) +
  scale_x_discrete(position = "top") +
  geom_tile(color = "black")+
  scale_fill_manual(values=c("white","#FF9900","#FF3333"),
                    breaks=c("none","sug","sig"),
                    labels=c(">1e-5","Suggestive (<1e-5)","Significant (<5e-8)"))+
  theme_minimal()+ # minimal theme
  theme(axis.text.x=element_text(vjust = 0, size = 7, hjust = 0.1, angle=45), axis.title.x=element_blank()) +
  theme(axis.text.y = element_text(vjust = 0.4, size = 7, hjust = 1.05), axis.title.y=element_blank())+
  geom_text(aes(variable, id, label = value), color = "black", size = 1.8) +
  guides(fill="none") +
  theme(plot.background = element_rect(fill = 'white', colour = 'white'),axis.line = element_line(colour = "black"))

ggheatmap
ggsave("novelloci.png",ggheatmap,scale=1.3,height=200,width = 160, dpi = 300, units = "mm", device='png')

#known
n<-a[a$Novel==0,]
n<-n[,-c("chromosome","start","end","size_bp","Novel")]
n[,3:22] <- lapply(n[,3:22], as.numeric)
nl = melt(n, id.vars = c("Gene","locus"),measure.vars = colnames(n)[-c(1,2)])
nl$p<--log10(nl$value)
nl<-nl[order(nl$locus),]
nl$id<-paste("Locus",nl$locus,nl$Gene,sep=" ")
nl$id<-factor(nl$id,levels=unique(nl$id))
nl$value<-signif(nl$value, digits=3)
nl$sig<-"none"
nl$sig[nl$value<1e-5]<-"sug"
nl$sig[nl$value<5e-8]<-"sig"

ggheatmap <- ggplot(nl, aes(variable, id, fill = sig))+
  scale_y_discrete(limits=rev) +
  scale_x_discrete(position = "top") +
  geom_tile(color = "black")+
  scale_fill_manual(values=c("white","#FF9900","#FF3333"),
                    breaks=c("none","sug","sig"),
                    labels=c(">1e-5","Suggestive (<1e-5)","Significant (<5e-8)"))+
  theme_minimal()+ # minimal theme
  theme(axis.text.x=element_text(vjust = 0, size = 7, hjust = 0.1, angle=45), axis.title.x=element_blank()) +
  theme(axis.text.y = element_text(vjust = 0.4, size = 7, hjust = 1.05), axis.title.y=element_blank())+
  geom_text(aes(variable, id, label = value), color = "black", size = 1.8) +
  guides(fill="none") +
  theme(plot.background = element_rect(fill = 'white', colour = 'white'),axis.line = element_line(colour = "black"))
ggheatmap
ggsave("knownloci.png",ggheatmap,scale=1.3,height=200, width = 160, dpi = 300, units = "mm", device='png')
