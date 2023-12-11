#genome summary plots
#09 March 2021

library(vioplot)


###########################################################################################
###########################################################################################
#shows how many scaffolds there are on each chr
#also how many bases assigned to each chr are placed on the chr vs unplaced

meta<-read.csv("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.meta.txt", header=T);str(meta)
par(mfrow=c(2,1), cex=1)

meta$Chr <- sapply(strsplit(meta$Scaffold, "_" ), FUN="[", 1);head(meta)
meta$Chr[grep("or", meta$Chr)] = "ChrUnknown"
meta$Chr[meta$Chr=="ChrUnknown"] = "Unk."
plot_order = c(1,12,15:21,2:11,13:14,23,24,22,25)
barplot(table(meta$Chr)[plot_order], cex.names=1, main="Scaffolds per chromosome", xlab="", ylab="Number of scaffolds", las=2)
text(x=(0:24)*1.2+.7, y=15, table(meta$Chr)[plot_order], cex=.75)
clist = split(meta, meta$Chr);clist

sumData <- data.frame("Chr"=names(clist), "Len"=NA, "Longest"=0);sumData
idx=0
for (i in names(clist)){
	idx = idx + 1
	sumData$Len[sumData$Chr==i] = sum(clist[[idx]]$total_chr_len)	
	sumData$Longest[sumData$Chr==i] = max(clist[[idx]]$scaf_len)
} 
sumData
sumData$percent = round(sumData$Longest/sumData$Len*100,1)
i
rownames(sumData) = sumData$Chr;sumData

barplot(sumData$Len[plot_order]/1000000, col="white", names=sumData$Chr[plot_order], xlab="", ylab="Megabases", main="Chromosome length", las=2)
barplot(sumData$Longest[plot_order]/1000000, col="gray", add=T, axes=F)
#text(paste0(sumData$percent[plot_order], "%"), x=(0:24)*1.2+.7, y=20, cex=.7, srt=90)
#legend(x=12.2, y=206, legend=c("The longest scaffold on each chromosome is shaded grey."), cex=1, box.col="white")


?text()

#Chr13bedfile = meta[meta$Chr==13,c(1,4,5)]
#Chr13bedfile$start = 0;Chr13bedfile
#write.table(Chr13bedfile, file="/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v4.chr13.bed", sep=",", row.name=F, col.names=F, quote=F)


###########################################################################################
###########################################################################################
#shows the genetic map compared to the physical map (code similar to DTmap.rqtl2.r, but this uses a the output file of that: Meriones_chromonome.geneticmap.csv)

allMap = read.table("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.geneticMap.csv", header=T, sep=",");str(allMap)

allMap$Chrfactor=as.factor(allMap$Chr);str(allMap)
allMap$Chrnumeric = sub("Chr", "", allMap$Chr);head(allMap)
allMap$Chrnumeric = ifelse(as.character(allMap$Chrnumeric)=="X", 22, as.numeric(as.character(allMap$Chrnumeric)));cat("there will be warning here - nothing to worry about - the ifelse replaces all the NAs with 'X'.\n");head(allMap)
allMap$ChrLen = NA
for(i in 1:nrow(allMap)){
	chr = allMap$Chr[i]
	allMap$ChrLen[i] = sumData$Len[sumData$Chr==chr]
}
allMap$ChrLenMb = allMap$ChrLen/1000000
allMap$ChrPosMb = allMap$ChrPos/1000000 #for the most part this is true...
head(allMap)

str(allMap)
par(mfrow=c(1,1))
plot(x=allMap$Chrnumeric, y=allMap$cM, pch="_", bty="L", axes=FALSE, ylab="cM or Mb", xlab="Linkage Groups and Chromosomes", main="Comparison of genetic and physcial maps")
axis(side=1, labels=c(1:21,"X"), at=c(1:22)+0.25, cex.axis=1);#the LG names that correspond to chromosome names.
axis(side=2, labels=c(0,50,100,150,200,250,300,350), at=c(0,50,100,150,200,250,300,350))
points(allMap$Chrnumeric+.5, y=allMap$ChrPosMb, pch="_")
segments(x0=allMap$Chrnumeric+0.05, y0=allMap$cM, x1=allMap$Chrnumeric+.45, y1=allMap$ChrPosMb, col="grey")
for (LG in 1:22){
	#the linkage group bars:
	x0 = LG
	x1 = LG
	y0 = 0
	y1 = max(allMap$cM[allMap$Chrnumeric==LG])
	segments(x0=x0, x1=x1, y0=y0, y1=y1)

	#the chromosome outline:
	x0 = LG+0.45
	x1 = LG+0.55
	y0 = 0
	y1 = max(allMap$ChrLenMb[allMap$Chrnumeric==LG])
	rect(xleft=x0, xright=x1, ybottom=y0, ytop=y1)
}
legend(x=17, y=350, legend=c("Linkage groups to the left, chromosomes to the right", "There are 2492 markers spread across 22 linkage groups"))

par(mfrow=c(5,5), cex=.9)
for (i in unique(allMap$Chr)[c(1,12,15:21,2:11,13,14,22)]){
	plot(x=allMap$ChrPosMb[allMap$Chr==i], y=allMap$cM[allMap$Chr==i], type="p", col=as.factor(allMap$Scaf[allMap$Chr==i]), main=i, ylab="cM", xlab="Mb", pch=19, lwd=2, cex=1.1)
	points(x=allMap$ChrPosMb[allMap$Chr==i], y=allMap$cM[allMap$Chr==i], type="l", lwd=2)
	
}

###########################################################################################
###########################################################################################
#shows the distribution of R, GC, and Gene density across a chr. Also now Entropy and Linguistic complexity  
read_in_R_GC_GeneDens<-function(FileHandle){
	df = read.table(FileHandle, header=T, skip=0, sep=",");dim(df)
	colnames(df) = c("Scaf", "start", "stop","ChrLen", "num_Ns", "GC", "R", "GeneDens");head(df)
	df$Chr = sapply(strsplit(df$Scaf, "_" ), FUN="[", 1)
	df$Chr[grep("or", df$Chr)] = "Unk."
	df$Chr[grep("Unk", df$Chr)] = "Unk."
	df$Scaf = sub("^.*ed_", "", df$Scaf);head(df)
	df$Chrnumeric = sub("Chr", "", sub("_.*$", "", df$Chr));head(df)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="MT", "24", df$Chrnumeric)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="Y", "23", df$Chrnumeric)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="Unk.", "25", df$Chrnumeric)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="X", "22", df$Chrnumeric)
	df$Chrnumeric = as.numeric(df$Chrnumeric)
	df$Chrfactor = as.factor(df$Chrnumeric)
	df$log10R = log10(abs(df$R))
	df$log10R[2622260] = 0 #this keeps the Y chr in the R plot. doesn't actually plot anything for ChrY though, just holds the place.
	df$log10R[2701069] = 0 #this keeps the MT chr in the R plot. doesn't actually plot anything for ChrY though, just holds the place.
	df$GeneDens[2701069] = 0 #this keeps the MT chr in the R plot. doesn't actually plot anything for ChrY though, just holds the place.
	df$log10R[2702315] = 0 #this keeps the Unk chr in the R plot. doesn't actually plot anything for ChrY though, just holds the place.

	return(df)
}




#dfEntropy10k = read.table("/Volumes/MulleySeqs_Working/Genome_comparison/Nessie/Meriones_chromonome_v6.1.nessieOut_E_l10000_s10000_asDF.csv", header=F, skip=1, sep=",");dim(dfEntropy10k)



read_in_Nessie<-function(FileHandle, type="Entropy"){
	df = read.table(FileHandle, header=F, skip=1, sep=",");dim(df)
	colnames(df) = c("Chr", "Scaf", "start", type);head(df)#, "A", "C", "G", "T")
	rows = nrow(df);rows
	df = rbind.data.frame(df, c("ChrMT", "ChrMT", 0, 0))
	df$start = as.numeric(df$start)
	df[,4] = as.numeric(df[,4])
	df$Chr[grep("or", df$Chr)] = "Unk."
	df$Chr[grep("Unk", df$Chr)] = "Unk."
	df$Scaf = sub("^.*ed_", "", df$Scaf);head(df)
	df$Chrnumeric = sub("Chr", "", sub("_.*$", "", df$Chr));head(df)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="Y", "23", df$Chrnumeric)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="MT", "24", df$Chrnumeric)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="Unk.", "25", df$Chrnumeric)
	df$Chrnumeric = ifelse(as.character(df$Chrnumeric)=="X", "22", df$Chrnumeric)
	df$Chrnumeric = as.numeric(df$Chrnumeric)
	df$Chrfactor = as.factor(df$Chrnumeric)
	return(df)
}



dfR_GC_GeneDens = read_in_R_GC_GeneDens("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.R_GC_Genedens.csv")
#dfR_GC_GeneDens2 = dfR_GC_GeneDens[!is.na(dfR_GC_GeneDens$log10R) & ! (dfR_GC_GeneDens$log10R==-Inf),] #just for the Recombination rate - apparently the vioplot doesn't like the NAs in there.

dfEntropy = read_in_Nessie("/Volumes/MulleySeqs_Working/Genome_comparison/Nessie/Meriones_chromonome_v6.1.nessieOut_E_l10000_s1000_asDF.csv", "Entropy")
dfLinComp = read_in_Nessie("/Volumes/MulleySeqs_Working/Genome_comparison/Nessie/Meriones_chromonome_v6.1.nessieOut_L_l10000_s1000_asDF.csv", "Complexity")


if(FALSE){
	
	#compare all chromosomes
	par(mfrow=c(5, 1))
	vioplot(dfR_GC_GeneDens$GC~dfR_GC_GeneDens$Chrnumeric, xlab="Chromosome", main="GC content", ylab="GC content", ylim=c(0,1), colMed=c(rep("white", 4), "red", rep("white",7), "red", rep("white", 100)), names=c(1:21, "X", "Y", "MT", "Unk."))
	vioplot(dfR_GC_GeneDens2$log10R~dfR_GC_GeneDens2$Chrnumeric, xlab="Chromosome", main="Recombination rate", ylab="cM/Mb log10", colMed=c(rep("white", 4), "red", rep("white",7), "red", rep("white", 100)), names=c(1:21, "X", "Y", "MT", "Unk."))
	vioplot(dfR_GC_GeneDens$GeneDens~dfR_GC_GeneDens$Chrnumeric, xlab="Chromosome", main="Gene density", ylab="Genes/Mb", colMed=c(rep("white", 4), "red", rep("white",7), "red", rep("white", 100)), names=c(1:21, "X", "Y", "MT", "Unk."))
	vioplot(dfEntropy$Entropy~dfEntropy$Chrnumeric, xlab="Chromosome", ylab="", main="Entropy", colMed=c(rep("white", 4), "red", rep("white",7), "red", rep("white", 100)), names=c(1:21, "X", "Y", "MT", "Unk."), ylim=c(.6,1))
	
	#vioplot(-log10(0.99999999-dfEntropy$Entropy)~dfEntropy$Chrnumeric, xlab="Chromosome", ylab="", main="Scaled Entropy", colMed=c(rep("white", 4), "red", rep("white",7), "red", rep("white", 100)), names=c(1:21, "X", "Y"))
	
	
	vioplot(dfLinComp$Complexity~dfLinComp$Chrnumeric, xlab="Chromosome", ylab="", main="Linguistic Complexity", colMed=c(rep("white", 4), "red", rep("white",7), "red", rep("white", 100)), names=c(1:21, "X", "Y", "MT", "Unk."))
	
	
	
	
	
	
	
	
	
	
	
	par(mfrow==c(1,1))
	vioplot(-log10(0.999999999-dfLinComp$Complexity)~dfLinComp$Chrnumeric, xlab="Chromosome", ylab="", main="Linguistic Complexity", colMed=c(rep("white", 4), "red", rep("white",7), "red", rep("white", 100)), names=c(1:21, "X", "Y", "MT", "Unk."))
	
	
	if(FALSE){ # because there are so many window - everything looks significant
	mod1 = aov(dfEntropy$Entropy~dfEntropy$Chrfactor2)
	summary(mod1)
	anova(mod1)
	TukeyHSD(mod1)
	levels(dfEntropy$Chrfactor2)
	}
	

}

############
#now plot chr-by-chr instead of all together.
dfR_GC_GeneDens = dfR_GC_GeneDens[order(dfR_GC_GeneDens$Chr, dfR_GC_GeneDens$Scaf, dfR_GC_GeneDens$start),]
dfRs<-split(dfR_GC_GeneDens, f=dfR_GC_GeneDens$Chr)
dfEntropy = dfEntropy[order(dfEntropy$Chr, dfEntropy$Scaf, dfEntropy$start),]
dfEs<-split(dfEntropy, f=dfEntropy$Chr)
dfLinComp = dfLinComp[order(dfLinComp$Chr, dfLinComp$Scaf, dfLinComp$start),]
dfLs<-split(dfLinComp, f=dfLinComp$Chr)


#end up with 
for(i in 1:length(dfRs)){
	dfRs[[i]] = dfRs[[i]][order(dfRs[[i]]$ChrLen, decreasing=TRUE),]
	Scafs = unique(dfRs[[i]]$Scaf);Scafs
	ScafCol_factors = cbind.data.frame("ScafCol"=as.factor(1:length(Scafs)), "Scaf"=Scafs);ScafCol_factors
	dfRs[[i]] = merge(x=dfRs[[i]], y=ScafCol_factors, by="Scaf", all.x=TRUE, all.y=FALSE)
	dfRs[[i]] = dfRs[[i]][order(dfRs[[i]]$ChrLen, decreasing=TRUE),]

}

ChrLens = unique(dfR_GC_GeneDens[c(1,4)]);ChrLens

for(i in 1:length(dfEs)){
	dfEs[[i]] = merge(x=dfEs[[i]], y=ChrLens, by="Scaf", all.x=TRUE, all.y=FALSE)
	dfEs[[i]] = dfEs[[i]][order(dfEs[[i]]$ChrLen, decreasing=TRUE),]
	Scafs = unique(dfEs[[i]]$Scaf)
	ScafCol_factors = cbind.data.frame("ScafCol"=as.factor(1:length(Scafs)), "Scaf"=Scafs)
	dfEs[[i]] = merge(x=dfEs[[i]], y=ScafCol_factors, by="Scaf", all.x=TRUE, all.y=FALSE)
	dfEs[[i]] = dfEs[[i]][order(dfEs[[i]]$ChrLen, decreasing=TRUE),]

}

for(i in 1:length(dfLs)){
	dfLs[[i]] = merge(x=dfLs[[i]], y=ChrLens, by="Scaf", all.x=TRUE, all.y=FALSE)
	dfLs[[i]] = dfLs[[i]][order(dfLs[[i]]$ChrLen, decreasing=TRUE),]
	Scafs = unique(dfLs[[i]]$Scaf)
	ScafCol_factors = cbind.data.frame("ScafCol"=as.factor(1:length(Scafs)), "Scaf"=Scafs)
	dfLs[[i]] = merge(x=dfLs[[i]], y=ScafCol_factors, by="Scaf", all.x=TRUE, all.y=FALSE)
	dfLs[[i]] = dfLs[[i]][order(dfLs[[i]]$ChrLen, decreasing=TRUE),]

}


palette(c("white", 'lightgrey'))#, 'lightsalmon', 'lightblue', 'lightyellow', 'palegreen', "red1", "peru"))


cbind.data.frame(1:25, names(dfRs))

#		i		Chr
#1     1        Chr1
#12   12        Chr2
#15   15        Chr3
#16   16        Chr4
#17   17        Chr5
#18   18        Chr6
#19   19        Chr7
#20   20        Chr8
#21   21        Chr9
#2     2       Chr10
#3     3       Chr11
#4     4       Chr12
#5     5       Chr13
#6     6       Chr14
#7     7       Chr15
#8     8       Chr16
#9     9       Chr17
#10   10       Chr18
#11   11       Chr19
#13   13       Chr20
#14   14       Chr21
#23   23        ChrX
#24   24        ChrY
#22   22       ChrMT
#25   25        Unk.



cbind(names(dfEs), names(dfRs), names(dfLs))
i=5
for(i in 1:length(names(dfRs))){


	Chr =dfRs[[i]]$Chr[1];Chr	
	cat(Chr, "\n")
	start = 1                  #for start of chr:   1
	end = max(nrow(dfEs[[i]]));end  #for end of the chr: max(nrow(L[[i]]))


	#pdf(file=paste0("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/summary_figures/", dfRs[[i]]$Chr[1], ".GC_R_GeneDens_Entropy_LinComp.pdf"), width = 7, height = 12)#for pdf:w=7, h=12, for png, w=700, h=1200
	
	png(file=paste0("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/summary_figures/", dfRs[[i]]$Chr[1], ".GC_R_GeneDens_Entropy_LinComp.png"), width = 700, height = 1200)#for pdf:w=7, h=12, for png, w=700, h=1200
	par(mfrow=c(4,1), cex=1)
	
	
	plot(y=dfRs[[i]]$GC*100,x=1:nrow(dfRs[[i]]), ylim=c(0,100), pch=19, cex=.4, col=rgb(0,0,0,.1), main="GC content", cex.main=1.5, xlab="Position along chromosome", ylab="GC percent", bty="L", cex.lab = 1, axes=F, type="n", xlim = c(start, end));
	segments(x0=1:nrow(dfRs[[i]]), y0=0, y1=100, col=as.factor(dfRs[[i]]$ScafCol))
	points(y=dfRs[[i]]$GC*100,x=1:nrow(dfRs[[i]]), type="p", pch=19, cex=.4, col=rgb(0,0,0,.1));
	axis(side=1, at=seq(0, nrow(dfRs[[i]]), by=10000), labels=seq(0, nrow(dfRs[[i]]), by=10000)*1000)
	axis(side=2, at=c(0:10)*10, labels=c(0:10)*10)
	
	mtext(sub("_.*$", "", dfRs[[i]]$Chr[1]), side=3, line=2, cex=2, adj=0)
	
	
	if(FALSE){
		plot(y=dfRs[[i]]$R,x=1:nrow(dfRs[[i]]), ylim=c(0,10), pch=19, cex=.4, col=rgb(0,0,0,.1), cex.main=2, xlab="Position along chromosome", ylab="cm/Mb", main="Recombination Rate", bty="L", cex.lab = 1, axes=F, type="n", xlim = c(start, end), cex.main=1.5);
		segments(x0=1:nrow(dfRs[[i]]), y0=0, y1=10, col=as.factor(dfRs[[i]]$ScafCol))
		points(y=dfRs[[i]]$R,x=1:nrow(dfRs[[i]]), type="p", pch=19, cex=.4, col=rgb(0,0,0,.1));
		axis(side=1, at=seq(0, nrow(dfRs[[i]]), by=10000), labels=seq(0, nrow(dfRs[[i]]), by=10000)*1000)
		axis(side=2, at=c(0:10), labels=c(0:10))
	}
	
	plot(y=dfRs[[i]]$GeneDens,x=1:nrow(dfRs[[i]]), ylim=c(0,1), pch=19, cex=.4, col=rgb(0,0,0,.1), cex.main=1.5, xlab="Position along chromosome", ylab="Genic bases/Mb", main="Gene Density", bty="L", cex.lab = 1, axes=F, type="n", xlim = c(start, end));
	segments(x0=1:nrow(dfRs[[i]]), y0=0, y1=1, col=as.factor(dfRs[[i]]$ScafCol))
	points(y=dfRs[[i]]$GeneDens,x=1:nrow(dfRs[[i]]),  type="p", pch=19, cex=.4, col=rgb(0,0,0,.1));
	axis(side=1, at=seq(0, nrow(dfRs[[i]]), by=10000), labels=seq(0, nrow(dfRs[[i]]), by=10000)*1000)
	axis(side=2, at=c(0:10)/10, labels=c(0:10)/10)

	
	if(FALSE){#the scaled entropy
		plot(y=-log10(1-dfEs[[i]]$Entropy), x=1:nrow(dfEs[[i]]), pch=19, cex=.4, ylim=c(0,6), col=rgb(0,0,0,.1), cex.main=1.5, xlab="Position along chromosome", ylab="", bty="L", cex.lab = 1.5, axes=F, type="n", xlim = c(start, end), main="Entropy");
		segments(x0=1:nrow(dfEs[[i]]), y0=0, y1=5, col=as.factor(dfEs[[i]]$ScafCol))
		points(y=-log10(0.9999999-dfEs[[i]]$Entropy),x=1:nrow(dfEs[[i]]),  type="p", pch=19, cex=.4, col=rgb(0,0,0,.1));
		axis(side=1, at=seq(0, nrow(dfEs[[i]]), by=10000), labels=seq(0, nrow(dfEs[[i]]), by=10000)*1000)
		axis(side=2, at=c(1,2,3,4,5), labels=c(0.9,0.99,0.999,0.9999, 0.99999), las=1)
	}
	
	plot(y=dfEs[[i]]$Entropy, x=1:nrow(dfEs[[i]]), pch=19, cex=.4, ylim=c(0.9,1), col=rgb(0,0,0,.1), cex.main=1.5, xlab="Position along chromosome", ylab="Entropy", bty="L", cex.lab = 1, axes=F, type="n", xlim = c(start, end), main="Entropy");
	segments(x0=1:nrow(dfEs[[i]]), y0=0, y1=5, col=as.factor(dfEs[[i]]$ScafCol))
	points(y=dfEs[[i]]$Entropy,x=1:nrow(dfEs[[i]]),  type="p", pch=19, cex=.4, col=rgb(0,0,0,.1)); 
	axis(side=1, at=seq(0, nrow(dfEs[[i]]), by=10000), labels=seq(0, nrow(dfEs[[i]]), by=10000)*1000)
	axis(side=2, at=c(.90, .92, .94, .96, .98,1), labels=c(.90, .92, .94, .96, .98,1), las=1)

	
	plot(y=dfLs[[i]]$Complexity,x=1:nrow(dfLs[[i]]), pch=19, cex=.4, ylim=c(0,1), col=rgb(0,0,0,.1), cex.main=1.5, xlab="Position along chromosome", ylab="Linguistic complexity", main="Lingiustic Complexity", bty="L", cex.lab = 1, axes=F, type="n", xlim = c(start, end));
	segments(x0=1:nrow(dfLs[[i]]), y0=0, y1=5, col=as.factor(dfLs[[i]]$ScafCol))
	points(y=dfLs[[i]]$Complexity,x=1:nrow(dfLs[[i]]),  type="p", pch=19, cex=.4, col=rgb(0,0,0,.1));
	axis(side=1, at=seq(0, nrow(dfLs[[i]]), by=10000), labels=seq(0, nrow(dfLs[[i]]), by=10000)*1000)
	axis(side=2, at=c(0:10)/10, labels=c(0:10)/10)

	dev.off()
	
	
	
}



unique(dfLs[[i]][,c(1,7)])





png(file=paste0("~/Dropbox/Post_Doc/Manuscripts/Gerbil_genome/figures/Entropy_Comparison.png"), width = 1000, height = 1200)#for pdf:w=7, h=12, for png, w=700, h=1200
par(mfrow=c(6,1), cex=1)
for(i in c(1,17,5,14,23,24)){ 
	Chr =dfEs[[i]]$Chr[1];Chr	
	cat(Chr, "\n")
	start = 1                  #for start of chr:   1
	end = max(nrow(dfEs[[1]]));end  #for end of the chr: max(nrow(L[[i]]))


	plot(y=dfEs[[i]]$Entropy, x=1:nrow(dfEs[[i]]), pch=19, cex=.4, ylim=c(0.9,1), col=rgb(0,0,0,.1), cex.main=1.5, xlab="Position (Mbp)", ylab="Entropy", bty="L", cex.lab = 1, axes=F, type="n", xlim = c(start, end), main=paste0(dfEs[[i]]$Chr[1]));
	segments(x0=1:nrow(dfEs[[i]]), y0=0, y1=5, col=as.factor(dfEs[[i]]$ScafCol))
	points(y=dfEs[[i]]$Entropy,x=1:nrow(dfEs[[i]]),  type="p", pch=19, cex=.4, col=rgb(0,0,0,.1)); 
	axis(side=1, at=seq(0, nrow(dfEs[[i]]), by=10000), labels=seq(0, nrow(dfEs[[i]]), by=10000)*1000/1e6)
	axis(side=2, at=c(.86, .88, .90, .92, .94, .96, .98,1), labels=c(.86, .88, .90, .92, .94, .96, .98, 1), las=1)
}
dev.off()


i=14







if(FALSE){

###########################################################################################
###########################################################################################
#This is the translation table between the HiFiAsm scaffolds, the intermediate DT scaffolds, the final DT scaffolds, and the chromonome scaffolds:
#it may be worth actualy calculating out the bases where each starts and stops all in the final Chr coordinates, but since they are chained together and that isn't a trivial calculation, I'm going to remove those bases here. This will tell you which Chr a given scaf from a given assembly is on, but not where it starts and stops. I'll leave that complexity for another time.
#as a hint if I ever come back to that: HiFiAsm scaffolds ("HiFiAsmScaf") were assembled into Dovetail's intermediate scaffolds ("IntermediateScaf") which were assembled into Dovetail final scaffolds ("DTscaf") from which the first bit of the names was extracted to get the scaffold ids ("Scaf"). These were then assembled into Chrs based on the genetic map and pools etc. So the coordinates will chain together, but to get a HiFiAsm scaf from the final "Chr" genome, you'll have to back calculate where the DT final Scaf is on that Chr, where the DTintermediateScaf is on that DTscaf, and where the HiFiAsmScaf is on that DTintermediateScaf... 


hictable<-read.table("/Volumes/MulleySeqs_Working/RAW/Dovetail/OmniC/hic.table.txt", header=F)
names(hictable)<-c("DTintermediateScaf", "HiFiAsmScaf", "HifiStart", "HifiEnd", "HiFiOri", "interStart", "interEnd")
hictable<-hictable[,c("DTintermediateScaf", "HiFiAsmScaf")];head(hictable)
agp<-read.table("/Volumes/MulleySeqs_Working/RAW/Dovetail/OmniC/final_assembly.agp", header=F, skip=2)
names(agp)<-c("DTfinalScaf", "DTstart", "DTend", "unknown", "letter", "DTintermediateScaf", "IscafStart", "IscafEnd", "IscafOri")
agp<-agp[,c("DTfinalScaf", "DTintermediateScaf")];head(agp)
agp<-agp[agp$DTintermediateScaf!=100,]#this is the record of where the runs of N's are that DT added when they joined scaffolds
t2<-merge(agp, hictable, by.x="DTintermediateScaf", by.y="DTintermediateScaf", all.x=T, all.y=T);head(t2)
t2$Scaffold = sub("__[0-9]+_contigs__length_[0-9]+","",t2$DTfinalScaf, perl=T);head(t2)
t3 = merge(x=t2, y=meta[,c("Scaffold", "Chr", "pool", "placed")], by.x="Scaffold", by.y= "Scaffold", all.x=T, all.y=T);t3 = t3[,c(5,1,3,2,4,6,7)];t3=t3[order(t3$HiFiAsmScaf),];t3=t3[order(as.numeric(sub("Scaffold_", "", t3$Scaffold))),];t3=t3[order(t3$pool),];t3=t3[order(as.numeric(t3$Chr)),];head(t3)#there are a lot of sorts here, coudl probably have put them all in one order() call, but whatever, this gets me what I want. 
write.table(t3, file="/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Scaffold_links_through_the_genome_versions.csv", sep=",", row.name=F, col.names=F, quote=F)

t3$ChrOrder = 1
t3$ChrOrder = ifelse((as.numeric(t3$Chr) %in% c(1:21)), as.numeric(t3$Chr), 24)
t3$ChrOrder = ifelse(t3$Chr=="X", 22, t3$ChrOrder)
t3$ChrOrder = ifelse(t3$Chr=="Y", 23, t3$ChrOrder)
t3

###########################################################################################
###########################################################################################
par(mfrow=c(2,1))

#This shows the number of reads mapped to each chr:
depth=read.table("/Volumes/MulleySeqs_Working/RAW/Dovetail/OmniC/DTG-OmniC-157_R1_001.fastq.gz.bam.idxstats.txt", header=F)
names(depth) = c("HiFiAsmScaf", "HiFiScafLength", "ReadsMapped", "ReadsUnmapped")
t4 = merge(x=depth, y=t3, by.x="HiFiAsmScaf", by.y="HiFiAsmScaf", all.x=T, all.y=T)
str(t4)
t5 = aggregate(cbind(t4$ReadsMapped, t4$HiFiScafLength), by=list(t4$ChrOrder), FUN=sum);colnames(t5) =c("Chr", "ReadsMapped", "Length");str(t5)
t5$coverage = t5$ReadsMapped / t5$Length;str(t5)
barplot(t5$coverage, names=c(1:21, "X", "Y", "Unk."), main="Average reads per base (~'coverage') for each Chr", ylab=  "reads mapped / chr length")





#This is plots out the mapping quality of reads from the omniC bam file for each Chr
readdepth <-read.table("/Volumes/MulleySeqs_Working/RAW/Dovetail/OmniC/Read_ChrPosQual_s1.01.txt", header=F)
names(readdepth)<-c("HiFiAsmScaf", "Pos", "MapQual")
readdepth<-merge(x= readdepth, y=t3, by="HiFiAsmScaf");str(readdepth)

vioplot(readdepth$MapQual~readdepth$ChrOrder, xlab="Chromosome", names=c(1:21, "X", "Y", "Unk."), ylab="Mapping quality", main="OmniC read-mapping quality on various chromosomes")




err <- read.csv("/Volumes/MulleySeqs_Working/RAW/Dovetail/OmniC/mHiC/bedtools_error_flags.txt", header=T);str(err)
err<-err[order(err$flag),]
err$type = as.factor(err$type)
err$flag = as.factor(err$flag)
str(err)


}