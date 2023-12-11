#written 15-06-2022 by Tom Brekke
#plots out the number of genes per chromosome as a function of the chrlength. 
library(vioplot)


chrLens = read.table("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.chrLengths.csv", header=T, sep=",");names(chrLens)<-c("Scaf", "Len");chrLens$Chr <- sapply(strsplit(chrLens$Scaf, "_"), FUN="[", 1);chrLens$Chr[grep("or", chrLens$Chr)] = "ChrUnknown"; head(chrLens)
#this bit sums the lengths across various scaffolds that are all on teh same chr and returns the length of the total scaf.
chrLensList = split(chrLens, chrLens$Chr);head(chrLensList)
CL = data.frame("Chr"=rep(NA, 25), "Len" = rep(NA, 25), "ChrNumeric" = c(1,10,11,12,13,14,15,16,17,18,19,2,20,21,3,4,5,6,7,8,9,24,25,22,23), "NumGenes" = rep(NA, 25), "NumUniqueGenes"=rep(NA, 25));CL
for (i in 1:length(chrLensList)){
	ChrDF = chrLensList[[i]]
	chr = ChrDF$Chr[1]
	len = sum(ChrDF$Len)
	CL$Chr[i] = chr
	CL$Len[i] = len
}
CL = CL[order(CL$ChrNumeric),]
CL




#now get the gene counts. 
#this "gff" file started as a gff, and I extraced the chr and the name of the gene.
#then count number of genes per chr and number of unique genes per chr
gff = read.table("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Merv6.1.cols_1_9.gff", header=F, sep=",");names(gff)<-c("Scaf", "Gene");gff$Chr <- sapply(strsplit(gff$Scaf, "_" ), FUN="[", 1); gff$Chr[grep("or", gff$Chr)] = "ChrUnknown";    head(gff)





gffList = split(gff, gff$Chr);
for(gffDF in gffList){
	Chr = gffDF$Chr[1];Chr
	NumGenes = nrow(gffDF);
	NumUniqueGenes = length(unique(toupper(gffDF$Gene)))
	CL$NumGenes[CL$Chr==Chr] = NumGenes
	CL$NumUniqueGenes[CL$Chr==Chr] = NumUniqueGenes

}

CL$LenMbp = CL$Len / 1000000
CL$DuplicationRate = 1- (CL$NumUniqueGenes/CL$NumGenes)
CL

newx = seq(min(CL$LenMbp), max(CL$LenMbp), length.out=100);newx

dev.new()


par(mfrow=c(1,3), cex=1.1)
plot(x=CL$LenMbp, y=CL$NumGenes, type="n", bty="L", xlim=c(-10,230),ylim=c(0,4300), xlab="Chromsome Length (Mbp)", ylab="Number of Genes annotated", main = "Number of Genes");
NG.lm = lm(NumGenes~LenMbp, data=CL);NG.lm;
#for the confidence intervals:
preds = predict(NG.lm, newdata=data.frame(LenMbp=newx), interval="confidence");preds
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
if(TRUE){#points or "CHRX" labels. True for point, false for chrnames
	points( x=CL$LenMbp, y=CL$NumGenes, pch=19, cex=.7)
	text("Chr13", x=94, y=3150, cex=.8)
	text("ChrX", x=174, y=4300, cex=.8)
	text("ChrY", x=80, y=1600, cex=.8)
}else{
	text(CL$Chr, x=CL$LenMbp, y=CL$NumGenes)
}
abline(NG.lm)
#mtext("A", side=3, line=2, cex=2, at=-70)

plot(x=CL$LenMbp, y=CL$NumUniqueGenes, type="n", bty="L", xlim=c(-10,230),ylim=c(0,2000), xlab="Chromsome Length (Mbp)", ylab="Number of Unique Genes annotated", main = "Number of Unique Genes");
NUG.lm = lm(NumUniqueGenes~LenMbp, data=CL)
preds = predict(NUG.lm, newdata=data.frame(LenMbp=newx), interval="confidence");preds
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
if(TRUE){
	points( x=CL$LenMbp, y=CL$NumUniqueGenes, pch=19, cex=.7)
	text("Chr13", x=125, y=30, cex=.8)
	text("ChrX", x=195, y=600, cex=.8)
	text("ChrY", x=72, y=28, cex=.8)
}else{	
	text(CL$Chr, x=CL$LenMbp, y=CL$NumUniqueGenes);
}
abline(NUG.lm)
#mtext("B", side=3, line=2, cex=2, at=-70)



#these are the numbers of GC genes on each chromosome - pulled from the R script called "get_outlier_genes_from_Roddys_data.R" hard-coded here to save time.
CL$GCGeneCount = c(40, 29, 21, 19, 24, 22, 19, 8, 28, 8, 17, 13, 29, 10, 6, 8, 7, 12, 13, 6, 12, 33+5+5, 7+9, 0, 0);CL
CL$GCdens = CL$GCGeneCount / (CL$Len/1e6) 
CL$color = rep("black", nrow(CL));CL$color[CL$Chr=="Chr13"]=  "red"



CL
CLNoNa= CL[1:(nrow(CL)),];CLNoNa

plot(x=CLNoNa$Len/1e6, y=CLNoNa$GCGeneCount, pch=19, cex=1, bty="L", xlim=c(0,200), main="GC-rich genes", xlab="Chromosome length (Mbp)", ylab="Number of GC-rich genes");
CLNoNa$LenMbp = CLNoNa$Len/1e6
mod.lm = lm(GCGeneCount ~ LenMbp, data=CLNoNa)
newx = seq(min(CLNoNa$LenMbp, na.rm=T), max(CLNoNa$LenMbp, na.rm=T), length.out=100);newx
preds = predict(mod.lm, newdata = data.frame(LenMbp = newx), interval="confidence");preds
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
points( x=CLNoNa$Len/1e6, y=CLNoNa$GCGeneCount, pch=19, cex=1)
text(x=100, y=30.5, "Chr13", cex=.8)
text(x=175, y=43, "ChrX", cex=.8)
text(x=80, y=18, "ChrY", cex=.8)
abline(mod.lm)



if(FALSE){
	hist (CL$GCdens, breaks=seq(0,.29,by=.02), main="GC-rich genes per Mbp", xlab="GC-rich gene density");
	text(x=.27, y=.8, "Chr13")
	text(x=.23, y=.8, "Chr9\nChrX")
	text(x=.21, y=.8, "Chr1\nChr21")
	text(x=.19, y=.8, "Chr19\nChrY")
	text(x=.17, y=.8, "Chr5")
	text(x=.15, y=.8, "Chr2\nChr6\nChr7\nChr11\nChr18")
	text(x=.13, y=.8, "Chr4\nChr12")
	text(x=.11, y=.8, "Chr3\nChr14\nChr16")
	text(x=.09, y=.8, "Chr20")
	text(x=.07, y=.8, "Chr8\nChr10\nChr15\nChr17")
}






plot(x=CL$LenMbp, y= CL$DuplicationRate, type="n", bty="L", xlim=c(-10,230),ylim=c(0,1), xlab="Chromsome Length (Mbp)", ylab="Duplication rate: 1 - (Number of Unique Genes / Number of Genes)", main = "Gene duplication rate")
DR.lm = lm(DuplicationRate~LenMbp, data=CL);
preds = predict(DR.lm, newdata=data.frame(LenMbp=newx), interval="confidence");preds
polygon(c(rev(newx), newx), c(rev(preds[ ,3]), preds[ ,2]), col = 'grey80', border = NA)
lines(newx, preds[ ,3], lty = 'dashed', col = 'red')
lines(newx, preds[ ,2], lty = 'dashed', col = 'red')
text(CL$Chr, x=CL$LenMbp, y=CL$DuplicationRate)
abline(DR.lm)

#extract the residual for each Chr from the regression line:
CL = cbind.data.frame(CL, "NumGenesResid"=NG.lm$residuals);CL
CL = cbind.data.frame(CL, "NumUniqueGenesResid"=NUG.lm$residuals);CL
CL = cbind.data.frame(CL, "DuplicationRateResid"=DR.lm$residuals);CL



#########################################################################################################
#analysis of Toby's TE stuff
par(mfrow=c(1,1))
#making a plot of observed TE bases per 100kb
TE<-read.table("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/TE_Analysis_v5.2_from_Hayward_Baril/chromosomeTEdistWithChrMeans.dataset.txt", header=T);names(TE)<-c("Scaf", "window", "obsDen", "expDen", "winno", "xstart", "xend", "mean", "stdev");TE$Chr = sapply(strsplit(TE$Scaf, "_" ), FUN="[", 1);head(TE);
TE$Chr[TE$Chr=="ChrUnknown"] = "Unk."
TE$Chr[TE$Chr=="Chr5|6"] = "Unk."
TE$Chr[TE$Chr=="Chr8|10"] = "Unk."
TE$Chr[TE$Chr=="Chr12|17"] = "Unk."
TE$Chr[TE$Chr=="Chr15|18"] = "Unk."
TE$Chr[TE$Chr=="Chr19|20"] = "Unk."

ChrDF = data.frame("Chr"=unique(TE$Chr), "ChrNumeric"=c(1,10,11,12,24,13,14,15,16,17,18,19,2,20,21,3,4,5,6,7,8,9,22,23));ChrDF[order(ChrDF$ChrNumeric),]
TE = merge(x=TE, y=ChrDF, by="Chr", all.x=T, all.y=F)
vioplot(TE$obsDen~TE$ChrNumeric, names=ChrDF$Chr[order(ChrDF$ChrNumeric)], xlab="", ylab="", main="Repetitive DNA density", las=2) #ylab="Observed repetitive bases per 100kb"
TE$Chr = as.factor(TE$Chr)
TE.lm = lm(TE$obsDen~TE$Chr, data=TE)
TE.lm
anova(aov(TE.lm))
TukeyHSD(aov(TE.lm))

for( Chr in ChrDF$Chr){ 
cat(Chr, mean(TE$obsDen[TE$Chr==Chr], na.rm=T), "\n")
}
