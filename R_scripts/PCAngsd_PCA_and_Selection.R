############################################################################################
#PCAngsd PCA
############################################################################################


#load data and do eigendecomposition
angsd_pca <- read.table("~/angsd_output.cov", quote="\"", comment.char="")
m<-as.matrix(angsd_pca)
e<-eigen(m)

#plot variances explained
variances<-vector()
for (i in 1:10){
  variances[i]<-round(e$values[i]/sum(e$values)*100,digits=2)}
barplot(variances)

#Check if PC1 is correlated with depth (needs dataframe "macro" with morphometric data from plsr analysis)
cor.test(e$vectors[,2],macro$Depth,method="spearman")


#plot according to depth but 90m with other symbol (needs dataframe "macro" with morphometric data from plsr analysis)
pdf("~/genomic_pc1_against_pc2_90_diff.pdf")
pch_vec<-rep(19,length(macro[which(macro$Sequenced==1),37]))
pch_vec[which(macro[which(macro$Sequenced==1),37]=="#16642AFF")]<-8
plot(cex=2,e$vectors[,c(1,2)],ylab=paste0("PC2 (",round(e$values[2]/sum(e$values)*100,digits=2),"%)"),xlab=paste0("PC1 (",round(e$values[1]/sum(e$values)*100,digits=2),"%)"),col=macro[which(macro$Sequenced==1),37],pch=pch_vec,lwd=3)
legend("bottomleft",y.intersp=1.2,pt.cex=2,pt.lwd=3,pch=c(rep(19,5),8),bty="n",legend=c("4m","12m","20m","40m","60m","90m"),col=c("#BBDAADFF","#97C38BFF","#74AC69FF","#468E3DFF","#247524FF","#16642AFF"))
dev.off()

############################################################################################
#PCAngsd Selection scan
############################################################################################

#Now load selection scan data
selection <- read.delim("~/selectionscan_output.txt", header=FALSE)
colnames(selection)<-c("Chr","Pos","S")

#transform chisquare statistic to pvalues and log10 of pvalues
pvalues<-pchisq(selection$S, df=1,lower.tail = F)
log<-(-log10(pvalues))
selection<-cbind(selection,pvalues,log)

#Make a qqplot similar to the PCAngsd selection scan paper to evaluate appropriatedness of the selection scan (looks good)
png("~/Dropbox/Depth_Transect/genomics/Pcangsd/qqnormown.png")
#Calculate expectations
exp.pvalues<-(rank(selection$pvalues, ties.method="first")+.5)/(length(selection$pvalues)+1)
#Make plot
plot(-log10(exp.pvalues), -log10(selection$pvalues),xlim=c(0,6),ylim=c(0,10),pch=19,col="darkgreen",xlab="Expected -log10(p)",ylab="Observed -log10(p)")
abline(0,1,col="red",lwd=2)
dev.off()


#fdr correction
fdr<-((p.adjust(selection$pvalues, method ="fdr", n = length(selection$Chr))))
fdr_log<-(-log10(p.adjust(selection$pvalues, method ="fdr", n = length(selection$Chr))))


#prepare data for creating a manhattan plot with function manhattan_custom() from manhattan_custom_plot.R
manhattan_data<-cbind(manhattan_data,fdr,fdr_log)
selection<-cbind(selection,fdr,fdr_log)

#create vector with names
names<-vector()
for(i in 1:length(which(manhattan_data$fdr<0.05))){
  names[i]<-paste0("snp",i)
}

#get indices of signifcatn snps
vec<-rep(0,length(manhattan_data$CHR))
indices<-which(manhattan_data$fdr<0.05)
for(i in 1:length(indices)){vec[indices[i]]<-names[i]}
manhattan_data<-cbind(manhattan_data,vec)
manhattan_data$vec<-as.character(manhattan_data$vec)

#make the plot
pdf("/Users/david/dropbox/Depth_Transect/genomics/selscan/plots/fastpca_fdr_manhattan_durchsichtig.pdf",width=20,height=6)
manhattan_custom(manhattan_data,p="fdr",snp="vec",highlight = names,chrlabs=as.character(chr_lab),suggestiveline=F,genomewideline=F,col=c(add.alpha("#BBDAADFF",0.8),add.alpha("#16642AFF",0.8)))
abline(h=(-log10(0.05)),col="black",lty=2)
dev.off()


#make a bedfile with significant snps
sign<-selection[which(selection$fdr<=0.05),]
bed<-cbind(sign$Chr,sign$Pos-1,sign$Pos)
write.table(bed,file="~/selection_fdr.bed",sep="\t",row.names = F, col.names = F,quote=F)
