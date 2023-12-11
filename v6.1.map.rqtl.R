#Written by TDB on 04-02-2021
#sorting out the genetic map with the dovetail genome
#Apparently rqtl2 doesn't do the map-making diagnostic bits, so I'll need to do r/qtl1. 

wd = "/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/stacks/rqtl/"
setwd(wd)
library(qtl)

mergeCHR<-function(object, chr1, chr1dir = "+", chr2, chr2dir="+", start=500){
	#put them in order as they would be joined, 
	if(!(chr1 %in% names(object$geno)) | !(chr2 %in% names(object$geno))){
		cat("\n\n!!one of the chrs doesn't exist - try again!!\n\n")
		return(object)
	}
	i=start
	if (chr1dir == "-"){
		object<-flip.order(object, chr1)
	}
	if (chr2dir =="-"){
		object<-flip.order(object, chr2)
	}
	markers<-markernames(object, chr2)
	for (marker in markers){
		i=i+1
		object<-movemarker(object, marker, chr1, i)		
	}	
	return(object)
}


Ripple<-function(object, chr, window=4, method="countxo", error.prob=0.0001, breakAt=10){
	rip<-ripple(object, chr=chr, window=window, method=method, error.prob=error.prob)
	a=1
	while(as.data.frame(summary(rip))$obligXO[1]!=as.data.frame(summary(rip))$obligXO[2]){
		cat("switching order", a,"\n")
		object<-switch.order(object, chr=chr, rip[2,])
		rip<-ripple(object, chr=chr, window=window, method=method,error.prob=error.prob)
		a=a+1
		if(a>breakAt){cat("It doesn't seem to be converging...\n");break}
	}
	return(object)
}


#read in the cross and look at the marker data:

MerMap<-read.cross(format="csvr", file="Meriones_chromonome_v6.1.geneticMap.rqtl.csvr", alleles=c("A", "B"), genotypes=c("A", "H", "B"), na.strings=c("UN","NA"), dir=wd, estimate.map=FALSE, convertXdata=FALSE, crosstype="f2"); summary(MerMap) #don't use the X option for make_rqtl_input_file.py script.

names(MerMap$geno)

xchr=c("X")
xmarkers = markernames(MerMap, chr=xchr);length(xmarkers)

nind(MerMap)
plotMissing(MerMap)

nt.byind <-ntyped(MerMap, "ind")
barplot(nt.byind[order(nt.byind)])
indstodrop<-names(nt.byind[nt.byind < 1500]);indstodrop #1500 for 90%missing, 2000 for 50%missing
indstokeep<-names(nt.byind)[!(names(nt.byind)%in%(indstodrop))];length(indstokeep)
MerMap<-subset(MerMap, ind=indstokeep);summary(MerMap)

nt.bymar <- ntyped(MerMap, "mar")
barplot(nt.bymar[order(nt.bymar)])

marsleft=length(markernames(MerMap));marsleft;par(mfrow=c(1,2), las=1);plot(ntyped(MerMap, "ind"), ylab="No. typed markers", main="No. genotypes by individual", pch=19, cex=.4, ylim=c(0,marsleft)) ;plot(ntyped(MerMap, "mar"), ylab="No. typed individuals",  main="No. genotypes by marker", ylim=c(0,150), pch=19, cex=.2); summary(MerMap);summaryMap(MerMap)



gt<-geno.table(MerMap)
gt[gt$P.value <0.05/totmar(MerMap),] #could drop ones that don't follow hwe, but that would throw out the X and I have the genome - so I don't need to throw out transmission ratio distorters



#the number of duplicate markers, no need to get rid of them here - they are locations in the genome
print(dup<-findDupMarkers(MerMap, exact.only=F)); todrop<-unlist(dup);length(todrop)
#MerMap<-drop.markers(MerMap, todrop)
summary(MerMap)


MerMap <- est.rf(MerMap)
#rf <- pull.rf(MerMap);lod <- pull.rf(MerMap, what="lod");dev.new();plot(as.numeric(rf), as.numeric(lod), xlab="Recombination fraction", ylab="LOD score", pch=20, cex=.4);dev.new()
todrop<-c("Chr14_11250691");MerMap<-drop.markers(MerMap, todrop) # potentially switched alleles, but switching doesn't fix anything. Just drop it.


#don't need to form linkage groups - There are already in chromosomes.
#now go by pool and esitmate rf:
plotRF(MerMap, chr=14)


#newmap.haldane<-est.map(MerMap, verbose=True, map.function="haldane")
#newmap.kosambi<-est.map(MerMap, verbose=TRUE, map.function="kosambi") #seems better than haldane - at least LGs are shorter. Still super long though.
#newmap.morgan<-est.map(MerMap, verbose=TRUE, map.function="morgan") # seems to take ages
newmap.cf<-est.map(MerMap, verbose=TRUE, map.function="c-f")
#plotMap(newmap.kosambi, newmap.haldane)
#plotMap(newmap.cf, newmap.haldane)

#plotMap(newmap.cf)
MerMap<-replace.map(MerMap, map=newmap.cf)


g<-pull.geno(MerMap)
gfreq <- apply(g, 1, function(a) table(factor(a, levels=1:3)))
gfreq <- t(t(gfreq)/colSums(gfreq))
par(mfrow=c(1,3), las=1)
for (i in 1:3){
	plot(gfreq[i,], ylab="Genotype frequency", main=c("AA", "AB", "BB")[i], ylim=c(0,1))
}


MerMap_backup<-MerMap
MerMap<-MerMap_backup;summaryMap(MerMap)


##########################################################
#then re-name chromsoomes:
chr_names = c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12", "Chr13", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19", "Chr20", "Chr21", "ChrX")
names(MerMap$geno) = chr_names
chrnames(MerMap)


write.cross(MerMap, format="csvr", filestem="Meriones_chromonome_v6.1.geneticmap_with_genotypes")
save.image(file="v6.1map.rqtl.1.RData")
load("v6.1map.rqtl.1.RData");summary(MerMap)


######################################################################################
######################################################################################
######################################################################################
######################################################################################
######################################################################################
#do import the map - the old seems to have some issue with male XBY genotypes. better to just import it correctly  - the DTmap.csv has correct genotypes.
library(qtl)
wd = "/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/stacks/rqtl/"
setwd(wd)
MerMap<-read.cross("csvr", dir=wd, file="Meriones_chromonome_v6.1.geneticmap_with_genotypes.csv",  na.strings=c("-", "NA"), alleles=c("A", "B"), genotypes=c("AA", "AB", "BB"), convertXdata=T);summary(MerMap)


######################################################################################
#now make the agp file that makes marey maps easy to plot. 
#this is a bit convoluted, but should get there in the end. 





pmap <-read.csv("Meriones_chromonome_v6.1.geneticMap.pmap.csv", header=F);names(pmap)=c("Marker", "Chr", "ChrPos");head(pmap)
gmap <- pull.map(MerMap, as.table=T); gmap$Marker = rownames(gmap); head(gmap)

allMap=merge(x=gmap, y=pmap, by="Marker");str(allMap);names(allMap) = c("Marker", "chr", "cM", "Chr", "ChrPos");head(allMap)
GC = read.csv("../../../Meriones_chromonome_v6.1.GC.txt", header=T, sep="	");head(GC)
allMap = merge(x=allMap, y=GC, by.x="Chr", by.y="name", all.x=T, all.y=T);names(allMap) = c("Chr", "Marker", "chr", "cM", "ChrPos", "ChrLen", "ChrNumN", "GC"); head(allMap)
#genome_links = read.csv("../../../../Meriones_Final/Scaffold_links_through_the_genome_versions_v5.csv", header=T, sep=",");genome_links<-na.omit(genome_links[,c(2,4,5)]);names(genome_links) = c("Chr", "Pool", "Scaffold");genome_links=unique.data.frame(genome_links);head(genome_links) #this has all the scaffolds from each genome version and which other scaffolds they become. It does not have the base positions of those - for that you need to find each individual agp file. The columns we care about are the pool, chromosomes, and annotation_jasmine-ban1818-mb-hirise-31o9t_12-07-2020__final_assembly.fasta

#allMap<-merge(x=allMap, y=genome_links, by=c("Chr", "Scaffold"), all.x=T, all.y=F)

allMap$ChrLenMb = allMap$ChrLen/1000000
allMap$ChrPosMb = allMap$ChrPos/1000000
head(allMap)

allMap$Chrfactor=as.factor(allMap$chr);str(allMap)
allMap$Chrnumeric = ifelse(as.character(allMap$chr)=="X", 22, as.numeric(as.character(allMap$Chrfactor)));str(allMap)
table(allMap$Chr, allMap$Chrnumeric)
head(allMap)

allMap=na.omit(allMap)
allMap = allMap[order(allMap$Chr, allMap$cM),];(head(allMap))


#################################################

#here's a quick look at the marey map of 7 - there seems to be something going on with one end of it - an inversion? perhaps that Dovetail made?
Chr7<-pull.map(MerMap, "7", T);head(Chr7)
Chr7$marker=rownames(Chr7);
Chr7$Scaf = sub("_[0-9]*$", "", Chr7$marker)
Chr7$Mb = as.numeric(sub("^.*_", "", Chr7$marker))
Chr7$cumMb = Chr7$Mb
Chr7$chr=7
Chr7$ScafDir = "+"
Chr7$ChrLen = 122584303
head(Chr7)


plot(x=Chr7$cumMb, y=Chr7$pos, col=as.factor(Chr7$Scaf), type="l")


head(allMap)






#plot the pmap vs the gmap
par(mfrow=c(1,1))
plot(x=allMap$Chrnumeric, y=allMap$cM, ylim=c(0,350), pch="_", bty="L", axes=FALSE, ylab="cM or Mb", xlab="Linkage Groups and Chromosomes", main="v6.1")
axis(side=1, labels=c(1:21,"X"), at=c(1:22)+0.25, cex.axis=1);#the LG names that correspond to chromosome names.
axis(side=2, labels=c(0,50,100,150,200,250,300,350), at=c(0,50,100,150,200,250,300,350))
points(allMap$Chrnumeric+.5, y=allMap$ChrPosMb, pch="_")
segments(x0=allMap$Chrnumeric+0.05, y0=allMap$cM, x1=allMap$Chrnumeric+.45, y1=allMap$ChrPosMb)
for (LG in 1:22){
	#the linkage group bars:
	x0 = LG
	x1 = LG
	y0 = 0
	y1 = max(allMap$cM[allMap$Chrnumeric==LG], na.rm=T)
	segments(x0=x0, x1=x1, y0=y0, y1=y1)

	#the chromosome outline:
	x0 = LG+0.45
	x1 = LG+0.55
	y0 = 0
	y1 = max(allMap$ChrLenMb[allMap$Chrnumeric==LG], na.rm=T) 
	#get_max_for_chr(LG, allMap)# this isn't the top - the top is based on the length of the czome...
	rect(xleft=x0, xright=x1, ybottom=y0, ytop=y1)
}
legend(x=8, y=350, legend=c("Linkage groups to the left, chromosomes to the right", "2492 markers spread across 22 linkage groups"))
nrow(allMap)



#################################################
#here is a plot of cM vs Mb:
par(mfrow=c(4,6))
for (Chr in unique(allMap$Chr)){
	plot(y=allMap$cM[allMap$Chr == Chr], x=allMap$ChrPos[allMap$Chr== Chr], main= Chr, ylab="centiMorgans", xlab="Megabases", type="b", ylim=c(0,max(allMap$cM[allMap$Chr == Chr], na.rm=T)), xlim=c(0,allMap$ChrLen[allMap$Chr== Chr][1]))
}
#dev.new()


summaryMap(MerMap)
#par(mfrow=c(1,1))
#plotGeno(MerMap, chr="1", ind=c(1:130),cutoff=10, include.xo=T, min.sep=2)






#################################################
#Here is the output files, a brief genetic map, and a complex one with all the metadata.

#then this needs to be done with the map too: 
genetic_map = allMap[c("Marker", "Chr", "ChrPos", "cM", "ChrLen")];genetic_map; #head(allMap)
genetic_map = genetic_map[order(genetic_map$Chr, genetic_map$ChrPos),];head(genetic_map)
genetic_map = genetic_map[c("Marker", "Chr", "ChrPos", "cM")];head(genetic_map)
write.table(x=genetic_map, file="Meriones_chromonome_v6.1.geneticmap.csv", sep=",", row.names=F, quote=F)
