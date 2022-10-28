#Define colors for species
library(RColorBrewer)
colors<-brewer.pal(12, name="Paired")
greys<-brewer.pal(9, name="Greys")
GF = colors[4]
SF = colors[8]
BF = colors[2]
K = greys[3]


#load data (fastpca result fdr)
freq_fastpca<- as.data.frame(read.csv("~/frequency_table.txt", sep=""))

#collapse close positions that are within 5000000bp from each other
new<-freq_fastpca[1,]
ld=5000000
for(i in 2:(length(freq_fastpca$chromo))){
  
  if(freq_fastpca[i,]$chromo==freq_fastpca[i-1,]$chromo){
    if( (freq_fastpca[i,]$position  <= (freq_fastpca[i-1,]$position+ld))){
      print(paste("Position",freq_fastpca[i,]$chromo,freq_fastpca[i,]$position,"has been deleted"))
    }else{new<-rbind(new,freq_fastpca[i,])}
  }else{new<-rbind(new,freq_fastpca[i,])}
  
}
freq_fastpca<-new


#add color column
freq_fastpca$col<-"grey"
freq_fastpca$lwd<-1


#highlight potentially introgressed positions
freq_fastpca[which(freq_fastpca$kilch>0.05 & freq_fastpca$gangfisch<0.01 & freq_fastpca$sandfelchen<0.01 & freq_fastpca$blaufelchen<0.01 & freq_fastpca$shallow<freq_fastpca$deep ),]$col<-"#16642AFF"
freq_fastpca[which(freq_fastpca$kilch<0.95 & freq_fastpca$gangfisch>0.99 & freq_fastpca$sandfelchen>0.99 & freq_fastpca$blaufelchen>0.99 & freq_fastpca$shallow>freq_fastpca$deep ),]$col<-"#16642AFF"


#make lines of interesting snps wider
freq_fastpca[which(freq_fastpca$col!="grey"),]$lwd<-2

#count how many kilch introgressed snps
kilch_intro<-length(which(freq_fastpca$col=="#16642AFF"))
kilch_intro/total


#make a plot with frequencies across all populations (figure 3)
pdf("~/Dropbox/Depth_Transect/genomics/frequencies/fastpca_fdr_frequencies.pdf",width=14,height=6)

#make emtpy plot
label1<-c("4m","12m","20m","40m","60m","90m","Kilch_pre","Gang_pre","Sand_pre","Blau_pre")
plot(0,xlim=c(1,10),pch=NA,ylim=c(0,1),axes=F,xlab=" ",ylab="Allele frequency")
axis(side=1,labels = label1,at=c(1:10),las=2)
axis(side=2)
box()


for(i in 1:length(freq_fastpca[,1])){
  #complile frequencies
  freq1<-as.vector(c(freq_fastpca[i,7],freq_fastpca[i,8],freq_fastpca[i,9],freq_fastpca[i,10],freq_fastpca[i,11],freq_fastpca[i,12],freq_fastpca[i,]$kilch,freq_fastpca[i,]$gangfisch,freq_fastpca[i,]$sandfelchen,freq_fastpca[i,]$blaufelchen))
  #draw lines
  points(freq1[1:6],type="l",col=freq_fastpca[i,]$col,lwd=freq_fastpca[i,]$lwd)
  points(x=c(6:10),y=freq1[6:10],type="l",col=freq_fastpca[i,]$col,lty=2,lwd=freq_fastpca[i,]$lwd)
  #add points and symbols
  points(freq1[1:6],pch=20,col=c("#BBDAADFF","#97C38BFF","#74AC69FF","#468E3DFF","#247524FF","#16642AFF"))
  points(x=c(7:10),y=freq1[7:10],pch=3,col=c(K,GF,SF,BF),cex=1.5,lwd=2)}


dev.off()




#########################################
#permutation test
#########################################
#load genome-wide allele frequencies for each population
genomewide_frequencies <- read.csv("~/genomewide_frequencies.txt", sep="")

#create vectors
permutations_darkgreen<-vector()

#set number of permutations
iterations=10000


########################################################################################################
for (k in 1:iterations){
  
  #number of total positions under selection
  n_observed=length(freq_fastpca$chromo)
  #create vector for results
  result<-vector()
  
  #sample different positions until with have sampled as many as are under selection 
  while(length(result)<n_observed){
    #sample a snps
    i<-sample(1:length(genomewide_frequencies$chr),1)
    int<-genomewide_frequencies[i,]
    
    #add color column
    int$col<-"grey"
    int$lwd<-1
    
    #sampled snps consistent with kilch introgression?
    if(  length(which(int$kilch>0.05 & int$gangfisch<0.01 & int$sandfelchen<0.01 & int$blaufelchen<0.01 & shallow<deep )) > 0  ){
      int[which(int$kilch>0.05 & int$gangfisch<0.01 & int$sandfelchen<0.01 & int$blaufelchen<0.01 & shallow<deep ),]$col<-"darkgreen"}
    if( length(which(int$kilch<0.95 & int$gangfisch>0.99 & int$sandfelchen>0.99 & int$blaufelchen>0.99 & shallow>deep)) > 0){
      int[which(int$kilch<0.95 & int$gangfisch>0.99 & int$sandfelchen>0.99 & int$blaufelchen>0.99 & shallow>deep ),]$col<-"darkgreen"}
  
    #store the result
    result[length(result)+1]<-int$col
  }
  #store the proportion of selected sites with kilch introgression of each permutation
  permutations_darkgreen[k]<-length(which(result=="darkgreen"))/length(result)
}
########################################################################################################

#calculate p value
#observed
observed_darkgreen_fastpca<-length(which(freq_fastpca$col=="#16642AFF"))/length(freq_fastpca$col)
#compare to expected and get p
length(which(permutations_darkgreen>=observed_darkgreen_fastpca))/length(permutations_darkgreen)
