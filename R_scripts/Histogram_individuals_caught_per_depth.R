#create a vector with number per individuals caught per 12 hours
hist_d<-rev(c(15,19,31,18,11,8.5,0))
#labels (corresponding to the nets
labels<-rev(c("4","12","20","40","60","90","120"))
#set the spaces between the bars (roughly equaling distance between nets)
spaces<-rev(c(8,8,20,20,30,30,4))

#create pdf with plot
pdf("~/depth_histogram.pdf")
barplot(hist_d,col=rev(c("#BBDAADFF","#97C38BFF","#74AC69FF","#468E3DFF","#247524FF","#16642AFF","black")),xlab="Number of individuals",ylab="Depth [m]",names=labels,width=1,horiz=T,space=spaces/30,xlim=c(0,32))
box()
dev.off()
