#load libraries
library(readr)
library(vioplot)
#load data
Depth_Transect_Data <- read_delim("~/Depth_Transect_Morphology_R.csv",  ";", escape_double = FALSE, trim_ws = TRUE)
#subset to only the caught macrophthalmus (verified with genomic data)
only_macro<-Data_SC[which(Data_SC$`Species Field`=="Coregonus macrophthalmus"),]
#size correction of the morphology data
for(i in c(1:6,8:23)){
  only_macro[,2+i]<-(lm(log10(only_macro[,2+i])~log10(only_macro$SLmean))$resid)}
#Remove standard length
macro<-only_macro[,-9]
#extract relevant columns with morphometric traits
morpho<-as.matrix(macro[,c(3:24)])

#perform plsr
partial_lsr<-plsr(macro$Depth~morpho,data=as.data.frame(morpho),validation="LOO")

#extract values of components
values<-partial_lsr$scores

#automatic component selection
cval<-crossval(partial_lsr, segments = 10)

#estimated MSEP as functions of the number of components (Figure 2)
pdf("~/MSEP.pdf",width=12)
plot(MSEP(cval), legendpos="topright")
selectNcomp(cval, method = "randomization", plot = TRUE,xlim=c(0,22),main="randomization")
selectNcomp(cval, method = "onesigma", plot = TRUE,xlim=c(0,22),main="onesigma")
dev.off()


#components 1 and 2
pdf("~/comps_1:2.pdf")
vs<-explvar(partial_lsr)
plot(values[,1],values[,2],cex=2,col=macro$color,pch=19,xlab=paste0("Component 1 (",round(vs[1],digits=1),"%)"),ylab=paste0("Component 2 (",round(vs[2],digits=1),"%)"),lwd=3)
plot(values[,3],values[,4],cex=2,col=macro$color,pch=19,xlab=paste0("Component 3 (",round(vs[3],digits=1),"%)"),ylab=paste0("Component 4 (",round(vs[4],digits=1),"%)"),lwd=3)

#highlighting the 90m fish with other symbol
pch_vec<-rep(19,length(macro$color))
pch_vec[which(macro$color=="#16642AFF")]<-8
plot(values[,1],values[,2],cex=2,col=macro$color,pch=pch_vec,xlab=paste0("Component 1 (",round(vs[1],digits=1),"%)"),ylab=paste0("Component 2 (",round(vs[2],digits=1),"%)"),lwd=3)

dev.off()

#loading plot
pdf("~/loadings.pdf",width=12)
par(mar=c(7, 4.1, 4.1, 2.1))
fake<-rep(" ",22)
plot(partial_lsr, "loadings", comps = 1:2, legendpos = "topleft",labels=fake,xlab=" ",ylab="Loading")
names<-colnames(morpho)
names<-gsub("mean","",names)
axis(at=1:22,side=1,labels=names,cex=0.5,las=2)
dev.off()



#components against depth
pdf("~/comp-vs-depth.pdf")
cor<-cor.test(values[,1],-macro$Depth,method="spearman")
plot(values[,1],-macro$Depth,col=macro$color,pch=19,xlab=paste0("Component 1 (",round(vs[1],digits=1),"%)"),ylab="Depth",cex=2,lwd=3)

cor<-cor.test(values[,2],-macro$Depth,method="spearman")
plot(values[,2],-macro$Depth,col=macro$color,pch=19,xlab="Depth",ylab="Comp 2",main=paste0("rho=",round(cor$estimate,digits=2)," , p=",round(cor$p.value,digits=6)),cex=2)

dev.off()
