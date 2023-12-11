#extracting the 300 outlier genes that Roddy discusses in his paper

setwd("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/GC_outlier_genes/")




###################################################################################################
###################################################################################################
###################################################################################################


#recreating Roddy's image to identify outlier genes:
df<-read.csv("dataset_7_normalised_rates_per_orthogroup.csv")

#recreate figure 3b:
plot(x=df$dS_SW_meriones_unguiculatus, y=df$dS_WS_meriones_unguiculatus, xlim=c(0,6), ylim=c(0,6))

#outliers are 2.5x the average:
ds_WS_threshold = 2.5 * mean(df$dS_WS_meriones_unguiculatus);ds_WS_threshold
ds_SW_threshold = 2.5 * mean(df$dS_SW_meriones_unguiculatus);ds_SW_threshold
#both are 2.5 - that must have been what the standardization did....
abline(h=ds_WS_threshold, v=ds_SW_threshold)

points(x=df$dS_SW_meriones_unguiculatus[df$dS_SW_meriones_unguiculatus>2.5], y=df$dS_WS_meriones_unguiculatus[df$dS_SW_meriones_unguiculatus>2.5], xlim=c(0,6), ylim=c(0,6), col="blue")

points(x=df$dS_SW_meriones_unguiculatus[df$dS_WS_meriones_unguiculatus>2.5], y=df$dS_WS_meriones_unguiculatus[df$dS_WS_meriones_unguiculatus>2.5], xlim=c(0,6), ylim=c(0,6), col="red")


points(x=df$dS_SW_meriones_unguiculatus[df$dS_WS_meriones_unguiculatus>2.5 & df$dS_SW_meriones_unguiculatus>2.5], y=df$dS_WS_meriones_unguiculatus[df$dS_WS_meriones_unguiculatus>2.5 & df$dS_SW_meriones_unguiculatus>2.5], xlim=c(0,6), ylim=c(0,6), col="purple")

gerbil_outliers = df$name[df$dS_WS_meriones_unguiculatus>2.5];length(gerbil_outliers);gerbil_outliers #387 of these

gerbil_outlier_genes = tolower(gsub('.{4}$', '', gerbil_outliers));gerbil_outlier_genes;length(gerbil_outliers)


write.table(x=paste0("Note=Similar to ", gerbil_outlier_genes, ":"), file = "Gene_names.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

#then use case insensitive grep on that new list against the current version gff file as such:
#grep "\tgene\t" ../Meriones_chromonome_v6.1.gff | grep  -i -F -f Gene_names.txt >  ./RoddysGenes.v6.1.gff 


#PROBLEM in version 5.1: of the 387 identified outliers, only 294 are annotated in our gerbil genome.. SOLVED: annotated the genes by hand and added them into a file incorporated below
#use a mummer alignment between mouse and gerbil to convert the mouse locations to gerbil coordinates. 
#for later versions, I should be ale to re-coordinated the gff file...? maybe??

#PROBLEM in version 6.1: the grep of 387 outliers returns 344 genes (some are duplicates though: 324 unique)... so still need to annotate by hand it looks like. 
#missing genes are these: 
#cut -f 9 RoddysGenes.v6.1.gff | cut -f 1 -d ":" | cut -f 2 -d";" > genes_found.txt
#grep -iv -F -f genes_found.txt Gene_names.txt | cut -f 3 -d " " | cut -f 1 -d ":" > genes_missing.txt
#Take that to biomart to get the chr,start,stop of those genes from mouse and use the mouse-merionse mummer alignment to identify the approximate locations of them in gerbil.

###################################################################################################
###################################################################################################
###################################################################################################

extract_gene_name <- function(df_line){
	#print(df_line)
	string = df_line[9]
	#search for 'Note=Similar to XXXXXX:' and extract the XXXXX bit
	#split on ; and keep #2, then split on : and keep #1:
	note = strsplit(x=string, split=';')[[1]][2]
	shortnote = strsplit(x=note, split=":")[[1]][1];shortnote
	gene = strsplit(x=shortnote, split=" ")[[1]][3];gene
	return(tolower(gene))
}
#use as such: gff$gene = apply(X = gff, MARGIN = 1, FUN = extract_gene_name)

#Does the genome plot: 
gffComplete<-read.table("../Meriones_chromonome_v6.1.gff", header=F, sep="	");names(gffComplete)<-c("Chr", "source", "type", "start", "stop", "dot", "dir", "dot2", "notes");gffComplete=gffComplete[gffComplete$type=="gene",]; gffComplete = gffComplete[order(gffComplete$Chr, gffComplete$start),];head(gffComplete);dim(gffComplete)

gff<-read.table("RoddysGenes.v6.1.gff", header=F, sep="	", quote="");names(gff)<-c("Scaf", "source", "type", "start", "stop", "dot", "dir", "dot2", "notes");gff = gff[gff$type=="gene",];gff$Chr = unlist(lapply(strsplit(gff$Scaf, split="_"), `[[`, 1));gff$gene = apply(X = gff, MARGIN = 1, FUN = extract_gene_name);str(gff)




genes_missing = gerbil_outlier_genes[!(gerbil_outlier_genes%in%gff$gene) ];genes_missing
genes_missing_annotated_by_eye = read.table("Missing_genes_annotated_by_eye_v6.1.txt", header=T);names(genes_missing_annotated_by_eye) <- c("musChr", "musStart", "musStop", "gene", "id", "Chr", "start");genes_missing_annotated_by_eye$gene <- tolower(genes_missing_annotated_by_eye$gene);genes_missing_annotated_by_eye<-genes_missing_annotated_by_eye[,c("Chr", "start", "gene")];
head(genes_missing_annotated_by_eye);


dim(gffComplete);head(gffComplete)
gffComplete<-merge(x=gffComplete, y=genes_missing_annotated_by_eye, by=c("Chr", "start"), all.x=T, all.y=T); 
gffComplete<-gffComplete[order(gffComplete$Chr, gffComplete$start),];
dim(gffComplete);head(gffComplete)


dim(gff)
gff<-merge(x=gff, y=genes_missing_annotated_by_eye, by=c("Chr", "start", "gene"), all.x=T, all.y=T);
head(gff, n=25)
dim(gff)


table(gff$gene %in% gerbil_outlier_genes) #needs to be 100%true
gff$gene[!(gff$gene %in% gerbil_outlier_genes)]#these are the problem ones: should be an empty list


table(gerbil_outlier_genes%in%gff$gene) #ideal would be 100%true, but as of now: 3 False, 384 True - these couldn't be annotated by hand for various reasons. mnot many, so move on.
gerbil_outlier_genes[!(gerbil_outlier_genes%in%gff$gene)]

chrlen<-read.table("Meriones_chromonome_v6.1.chrLengths_CHRSONLY.csv", header=T, sep=",");chrlen$ChrN = 1:24; str(chrlen)
 cent_locs<-read.csv("../cent_locations_rough_v6.1.bed", header=F);cent_locs<-cent_locs[,-c(4,5)];names(cent_locs)=c("Chr", "centStart", "centStop");head(cent_locs)
 chrlen<-merge(x=chrlen, y=cent_locs, all.x=T, all.y=T, by="Chr");names(chrlen)<-c("Chr", "Len", "ChrN", "centStart", "centStop");chrlen
alldata = merge(x=gff, y=chrlen, by=c('Chr'), all.x=T, all.y=F)

ChrN_df<-unique(cbind.data.frame(alldata$Chr, alldata$ChrN));names(ChrN_df)<-c("Chr", "ChrN");ChrN_df<-rbind.data.frame(ChrN_df, c("ChrY", 23), c("Chr20", 20));ChrN_df$ChrN<-as.numeric(ChrN_df$ChrN);ChrN_df<-ChrN_df[order(ChrN_df$ChrN),];str(ChrN_df);ChrN_df

gff<-merge(x=gff, y=ChrN_df, by="Chr", all.x=T, all.y=F)
gff_complete<-merge(x=gffComplete, y=ChrN_df, by="Chr", all.x=T, all.y=F)

teloSeqs=c("CTAACC","TAACCC","AACCCT","ACCCTA","CCCTAA","CCTAAC","GGTTAG","GTTAGG","TTAGGG","TAGGGT","AGGGTT","GGGTTA")
telomeres = read.csv("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/TandemRepeatFinder/Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000.simplified.csv");dim(telomeres);telomeres=telomeres[telomeres$consensus %in% teloSeqs,];dim(telomeres);telomeres<-merge(x=telomeres, y=ChrN_df, by="Chr", all.x=T, all.y=F);telomeres<-telomeres[order(telomeres$Chr, telomeres$Start),];head(telomeres)
hist(telomeres$numCopies, breaks=2000, xlim=c(0,100), main="Number of copies in telomere repeat arrays")
abline(v=70, col="red")
telomeres=telomeres[telomeres$numCopies>70,]

outliers = as.data.frame(table(alldata$Chr));outliers;dim(outliers)
chrlen<-merge(x=chrlen, y=outliers, by.x='Chr', by.y="Var1",all=T)
chrlen<-chrlen[order(chrlen$ChrN),];
chrlen$centMean = rowMeans(chrlen[c("centStart", "centStop")]);chrlen
r=1.2
s=1

table(alldata$Chr)#number of outlier genes on each chr
table(gffComplete$Chr)#number of total genes on each chr.

barplot(chrlen$Len/1000000, names=chrlen$Chr, col="white", main="Meriones genome with high-GC genes marked", ylab="position in genome (Mbp)", ylim=c(0,220));
#genes
segments(y0=alldata$start/1000000, x0=(r*alldata$ChrN)-s, y1=alldata$start/1000000, x1=(r*(alldata$ChrN+.83)-s), col="red")
#telomeres
segments(y0=telomeres$Start/1000000, x0=(r* telomeres$ChrN)-s, y1= telomeres$Start/1000000, x1=(r*(telomeres$ChrN+.83)-s), col=rgb(0,1,0,1))
#centromeres
rect(ybottom=chrlen$centStart/1000000, xleft=(r*chrlen$ChrN)-s, ytop=chrlen$centStop/1000000, xright=(r*(chrlen$ChrN+.83)-s), col=rgb(0,0,1,.1), border=rgb(0,0,1,1))

text(chrlen$Freq, x=(r*chrlen$ChrN)-s/2, y=210)
legend(x=11, y=190, legend=c("GC outlier genes", "Centromeres", "Telomere repeats (>70x copies)"), fill=c("red", rgb(0,0,1,.3), rgb(0,1,0,.3)))
#dev.new()








###################################################################################################
###################################################################################################
###################################################################################################
#calculates the genome stats.

#stats 
#read in gff file, extract genes, randomly draw the same number of genes and compare how far apart they are. (by chromosome? or can I do the whole gneome?)




#compare_distances_all <- function(gff){#takes in a gff file and returns a matrix of the all-by-all distances between each gene in it. 
	# mat = matrix(data=NA, nrow=nrow(gff), ncol=nrow(gff))
	
	# for(i in 1:nrow(gff)){
		# #extract the focal gene.
 		# for(j in i:nrow(gff)){#compare it to every other gene
			# if(gff$Chr[i] == gff$Chr[j]){ #only record if genes are on teh same chr
				# if(gff$start[j] != gff$start[i]){	#leave NA if they are the same gene (i.e. i=j)
					# mat[i,j] = abs(gff$start[i] - gff$start[j])
				# }
			# }
		# } 
	# }
	# return(mean(mat, na.rm=T)/1000000)
	# #return(mat)
# }



#this one is very slow! the next is far quicker.
compare_distances_nearestNeighbor <- function(gff){#takes in a gff file and returns a matrix of the all-by-all distances between each gene in it. 
	mat = matrix(data=NA, nrow=nrow(gff), ncol=nrow(gff))
	
	for(i in 1:nrow(gff)){
		#extract the focal gene.
		for(j in 1:nrow(gff)){#compare it to every other gene
			if(gff$Chr[j] == gff$Chr[i]){ #only record if genes are on the same chr
				#cat(i, j, fchr, fstart, tchr, tstart, "\n")
				if(gff$start[i] != gff$start[j]){	#leave NA if they are the same gene (i.e. start[i]=start[j])
					mat[i,j] = abs(gff$start[j] - gff$start[i])
				}
			}
		} 
	}
	mat = mat[, !apply(is.na(mat), 2, all)] #the matrix probably has some columns of all NAs, which means the next call would generate a bunch of Infs . drop them out here instead. The all-NA cols are b/c of chrs with a single gene on them
	mins = apply(mat,2,min, na.rm=T)
	return(mean(mins, na.rm=T)/1000000)
}


compare_distances_nearestNeighbor_quicker <- function(gff){#takes in a gff file, orders it, and calculates teh distance between each gene, i, and it's neighbours, i-1 and i+1. then chooses which of those is smaller. returns teh average of the distances. 
	#make an empty vector: list of distances same length as genes
	vec = rep(NA, nrow(gff))
	gff = gff[order(gff$Chr, gff$start),]#gff should be ordered, but just in case it isn't
	#for each gene, compare dist to gene [i-1] and to gene [i+1] - need to sort both edge cases specially: first and last entires, if they are on different chrs, return NA
	#edgecase 1, only compare gene i and gene i+1:
	vec[1] = ifelse(gff$Chr[1]==gff$Chr[2], gff$start[2] - gff$start[1], NA) #as long as they are on the same chr, take the difference, otherwise, return NA
	#mid cases, compere gene i to gene i-1 and to gene i+1 and take the minimum:
	for (i in 2:(nrow(gff)-1)){
		dist1 = ifelse(gff$Chr[i]==gff$Chr[i-1], gff$start[i] - gff$start[i-1], NA)
		dist2 = ifelse(gff$Chr[i]==gff$Chr[i+1], gff$start[i+1] - gff$start[i], NA)
		vec[i] = min(dist1, dist2)
	}
	#edgecase2, only compare gene i and i-1:
	vec[nrow(gff)] = ifelse(gff$Chr[nrow(gff)]==gff$Chr[(nrow(gff)-1)], gff$start[nrow(gff)] - gff$start[(nrow(gff)-1)], NA) #as long as they are on the same chr, take teh distance between gene i and i-1 where i is the final entry
	
	return(mean(vec, na.rm=T)/1000000)
}

#are they the same?
compare_distances_nearestNeighbor(gff)
compare_distances_nearestNeighbor_quicker(gff)
#yes. carry on.

reps = 1000000
dist = rep(NA, reps)
n = nrow(gff)
x = 1:nrow(gffComplete)
observed = compare_distances_nearestNeighbor_quicker(gff)

for(i in 1:reps){#takes my desektop about an hour to run 1e6 simulations with *_quicker, takes two days with the original version.
	if(i %%100 == 0){
		cat(i/reps*100, "%\r")
	}
	sIDX = sample(x=x, size=n)
	gffSample = gffComplete[sIDX,]
	dist[i] = compare_distances_nearestNeighbor_quicker(gffSample)
}
dist_nearestNeighbor = na.omit(dist)


plot_distNearestNeighbor_perms = function(){
	observed = compare_distances_nearestNeighbor_quicker(gff)
	hist(dist_nearestNeighbor, breaks=100, xlim=c(1.5,3.7), main=paste0("Average distance to nearest neighbor\n", reps, " permutations"), xlab="Average distance to the nearest neighor (Mbp)")
	abline(v=observed, col="red", lwd=2)

	pval = sum(dist_nearestNeighbor < observed) / length(dist_nearestNeighbor); pval#the proportion of measurments that are less than the observed out of the entire distances list gives us a p-value
	#or how many sigmas away from the mean also works:
	m = mean(dist_nearestNeighbor);m
	sigma = sd(dist_nearestNeighbor);sigma
	sds = abs(observed - m)/sigma
	text(x=2.2, y=45000*1, paste0("p = ", round(pval, 2)), cex=.7)
	text(x=2.2, y=40000*1, paste0("mean = ", round(m, 2)), cex=.7)
	text(x=2.2, y=35000*1, paste0("observed = ", round(observed, 2)), cex=.7)
	text(x=2.2, y=30000*1, paste0(round(sds, digits=2), "sd from the mean"), cex=.75)
}
plot_distNearestNeighbor_perms()


par(cex=1)


save.image(file="get_outlier_genes_from_Roddys_data.RData")
load(file="get_outlier_genes_from_Roddys_data.RData")






###################################################################################################
###################################################################################################
###################################################################################################
#compare clusters

if(FALSE){
compare_clusters <- function(gff, extracts){#takes in a gff file and returns a matrix of the all-by-all distances between each gene in it. 
	#the gff needs to be the entire gff.
	#extracts is an ordered vector of the genes to extract.

	#look for adjacencies
	count = 0
	run=1
	runs = NULL
	flag=TRUE #this will discount runs of many adjacent genes and count them as one cluster. i.e. only the first pair will add to teh count.
	
	for (i in 1:(length(extracts)-1)){
		if((extracts[i+1] == extracts[i]+1) & (gff$Chr[i+1] == gff$Chr[i])){ #makes sure that adjacent genes are also on the same chr, so that the final gene on chr1 and first gene on chr2 don't count as adjacent
			run = run + 1
			if(flag){
				count = count + 1
				flag = FALSE
			}
		}else{
			if(! flag){
				runs=c(runs, run)
			}#counts the number of genes in the run
			flag=TRUE
			run=1
		}
		
	}
	if(! flag){#have to check that the final entry wasn't a cluster
		runs=c(runs, run)
	}
	return(list(count=count, runs=runs))
}


reps = 1e6
dist = rep(NA, reps);
n = nrow(gff);n
x = 1:nrow(gffComplete);

observed_extracts =  (1:nrow(gffComplete))[gffComplete$notes %in% gff$notes]
observed = compare_clusters(gffComplete, observed_extracts)


nClusts = NULL
runLengths = NULL
avRunLengths = NULL
for(i in 1:reps){
	if(i %%100 == 0){
		cat("\t", i/reps*100, "%\r")
		#if(i %%10000 == 0){
		#	hist(nClusts, main="Histogram of number of clusters of 2 or more genes", xlim=c(-1,30), breaks=-1:30)
		#abline(v=observed[[1]])
		#}
	}	
	extracts = sample(x=x, size=n); 
	extracts = extracts[order(extracts)];
	out = compare_clusters(gffComplete, extracts);out
	nClusts = c(nClusts, out[[1]])
	avRunLengths = c(avRunLengths, ifelse(out[[2]], mean(out[[2]]), NA))
	runLengths = c(runLengths, out[[2]])
}

hist(nClusts, main=c("Histogram of number of clusters of 2 or more genes", paste0(reps, " permutations")), xlim=c(-1,50), breaks=-1:30, xlab="number of clusters of 2+ genes")
abline(v=observed[[1]])

pval = sum(nClusts > observed[[1]]) / length(nClusts); pval#the proportion of measurments that are less than the observed out of the entire distances list gives us a p-value
#or how many sigmas away from the mean also works:
m = mean(nClusts);m
sigma = sd(nClusts);sigma
sds = abs(observed[[1]] - m)/sigma
text(x=20, y=150000, paste0("p = ", pval))
text(x=20, y=140000, paste0("Mean = ", round(m, 2)))
text(x=20, y=130000, paste0("observed = ", observed[[1]]))
text(x=20, y=120000, paste0(round(sds, 1), "sd from the mean"))



#hist(runLengths, main="Histogram of size of clusters of 2 or more genes", xlim=c(0,4))
#abline(v=observed[[2]])

#pval = sum(runLengths > observed[[2]]) / length(runLengths); pval#the proportion of measurments that are less than the observed out of the entire distances list gives us a p-value
#or how many sigmas away from the mean also works:
#m = mean(runLengths);m
#sigma = sd(runLengths);sigma
#abs(observed[[2]] - m)/sigma


hist(avRunLengths, main=c("Histogram of average size of clusters of 2 or more genes", paste0(reps, " permutations")), xlim=c(0,4))
abline(v=mean(observed[[2]]))

pval = sum(avRunLengths > mean(observed[[2]])) / length(avRunLengths); pval#the proportion of measurments that are less than the observed out of the entire distances list gives us a p-value
#or how many sigmas away from the mean also works:
m = mean(avRunLengths);m
sigma = sd(avRunLengths);sigma
sds = abs(mean(observed[[2]]) - m)/sigma
text(x=.9, y=2000000, paste0("p = ", round(pval, 5)))
text(x=.9, y=1850000, paste0("Mean = ", round(m, 2)))
text(x=.9, y=1700000, paste0("observed = ", round(mean(observed[[2]]), 2)))
text(x=.9, y=1550000, paste0(round(sds, 1), "sd from the mean"))
}


save.image(file="get_outlier_genes_from_Roddys_data2.RData")



###################################################################################################
###################################################################################################
###################################################################################################
#calculate location of clusters in relation to telomeres. 

load("get_outlier_genes_from_Roddys_data2.RData")
ls()
calc_dist_to_cent<-function(gff_entry){
	#takes in a gene as a line in a gff file and then calculates the distance to the centromere along the arm as a percentage of the distance between the telo and the centromere.
	
	#get base-pairs locations
	Chr = gff_entry[1];Chr
	if(! Chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19", "Chr20", "Chr21", "Chr22", "ChrX")){return(NA)}#drop out things on Chr13, the Y, and unplaced scaffolds. i.e anything that is on a chr that doesn't have a centromere
	Start = as.numeric(gff_entry[2]);Start
	Cent_loc = as.numeric(chrlen[chrlen$Chr==Chr, 7]);Cent_loc
	ChrLen = as.numeric(chrlen[chrlen$Chr==Chr,2]);ChrLen
	
	#get length of chr arm (must decide what side of centromere it's on)
	Arm_len = as.numeric(ifelse(Start>Cent_loc, ChrLen-Cent_loc, Cent_loc));Arm_len
	
	#calc length in bps from cent
	dist_from_cent = abs(Cent_loc - Start);dist_from_cent

	#calc length in percent from cent
	percent_from_cent = dist_from_cent / Arm_len * 100; percent_from_cent
	
	#return percent
	return(percent_from_cent)
	
}

gffComplete$percent_to_centromere = apply(X=gffComplete, MARGIN=1, FUN=calc_dist_to_cent)
gff$percent_to_centromere=apply(X=gff, MARGIN=1, FUN=calc_dist_to_cent)

plot_distAlongCzomeArm = function(){
	hist(x=gffComplete$percent_to_centromere, breaks=100, main=c("Gene location"), xlab="Chromosome arm\nCentromere (0) to Telomere (100)", freq=F, col=rgb(0,0,1,.2), ylim=c(0,.11))
	hist(gff$percent_to_centromere, breaks=100, freq=F, add=T, col=rgb(1,0,0,.2))
	legend(x=10, y=.1, legend=c(paste0("All genes, n=",nrow(gffComplete)), paste0("GC outliers, n=", nrow(gff))), fill=c(rgb(0,0,1,.2), rgb(1,0,0,.2)), cex=.7)
	#is the observed distribution different than the genome-wide dist? yes, resoundingly so.
	pval = t.test(x=gffComplete$percent_to_centromere, y=gff$percent_to_centromere)$p.value;pval
	text(y=0.08, x=70, paste0("p = ", round(pval, 5)), cex=.7)
	wilcox.test(x=gffComplete$percent_to_centromere, y=gff$percent_to_centromere)
}
plot_distAlongCzomeArm()


reps = 1e6
dist = rep(NA, reps)
n = nrow(gff)
x = 1:nrow(gffComplete)
observed = mean(gff$percent_to_centromere, na.rm=T);observed


for(i in 1:reps){
	if(i %%1000 == 0){
		cat(i/reps*100, "%\r")
	}
	sIDX = sample(x=x, size=n);sIDX
	dist[i] = mean(gffComplete$percent_to_centromere[sIDX], na.rm=T);dist[i]
}
dist_along_CzomeArm = na.omit(dist)


plot_distAlongCzomeArm_perms = function(){
	observed = mean(gff$percent_to_centromere, na.rm=T);observed
	hist(dist_along_CzomeArm, breaks=20, main=c("Average gene location", paste0(reps, " permutations")), xlab=c("Chromosome arm", "Centromere (0) to Telomere (100)"), xlim=c(0,100))
	abline(v=observed, col="red", lwd=2)
	pval = sum(dist_along_CzomeArm > observed) / length(dist_along_CzomeArm);pval#the proportion of measurments that are less than the observed out of the entire distances list gives us a p-value
	m = mean(dist_along_CzomeArm);m
	sigma = sd(dist_along_CzomeArm);sigma
	sds = abs(observed - m)/sigma
	text(x=20, y=23000*8, paste0("p = ", round(pval, 5)), cex=.7)
	text(x=20, y=19500*8, paste0("Mean = ", round(m, 2)), cex=.7)
	text(x=20, y=16000*8, paste0("observed = ", round(observed, 2)), cex=.7)
	text(x=20, y=12500*8, paste0(round(sds, 1), "sd from the mean"), cex=.7)
}
plot_distAlongCzomeArm_perms()






calc_dist_to_nearest_telo<-function(gff_entry){
	#extract gene location in bp
	Chr = gff_entry[1];Chr
	if(! Chr %in% c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5", "Chr6", "Chr7", "Chr8", "Chr9", "Chr10", "Chr11", "Chr12", "Chr14", "Chr15", "Chr16", "Chr17", "Chr18", "Chr19", "Chr20", "Chr21", "Chr22", "ChrX")){return(NA)}#drop out things on Chr13, the Y, and unplaced scaffolds. i.e anything that is on a chr that doesn't have a centromere
	Start = as.numeric(gff_entry[2]);Start
	
	
	#find closest telomere
	telos = telomeres[telomeres[,1]==Chr,];telos
	old_dist = 100000000000000
	for(i in 1:nrow(telos)){
		dist = abs(telos[i,2]- Start);dist
		if(dist < old_dist){#as long as you're getting closer to a tel0, keep going, once it's getting further away, then bail out.
			old_dist=dist;dist
		}else{
			break
		}
	}
	return(old_dist)
}
calc_dist_to_nearest_telo(gff[1,])







gffComplete$dist_to_nearest_telo = apply(X=gffComplete, MARGIN=1, FUN=calc_dist_to_nearest_telo)
gff$dist_to_nearest_telo=apply(X=gff, MARGIN=1, FUN=calc_dist_to_nearest_telo)



par(mfrow=c(2,1))

plot_distNearestTelo = function(log=FALSE){
	if(log){
		hist(x=log10(gffComplete$dist_to_nearest_telo), breaks=seq(2.5, 8.0, by=.05), main=c("Distance to nearest telomere repeat"), xlab="Distance to nearest telomere repeat (log10 bp)", freq=F, col=rgb(0,0,1,.2), ylim=c(0,2));hist(log10(gff$dist_to_nearest_telo), breaks=seq(2.8, 8.8, by=.05), freq=F, add=T, col=rgb(1,0,0,.2))
		pval = t.test(x=log10(gffComplete$dist_to_nearest_telo), y=log10(gff$dist_to_nearest_telo))$p.value;pval
		legend(x=3, y=2, legend=c(paste0("All genes, n=",nrow(gffComplete)), paste0("GC outliers, n=", nrow(gff))), fill=c(rgb(0,0,1,.2), rgb(1,0,0,.2)), cex=.7)
		text(x=6.5, y=1.5, paste0("p = ", round(pval, 6)), cex=.7)
	}else{
		hist(x=gffComplete$dist_to_nearest_telo/1000000, breaks=100, main=c("Distance to nearest telomere repeat"), xlab="Distance to nearest telomere repeat (Mbp)", freq=F, col=rgb(0,0,1,.2), xlim=c(0,1e2), ylim=c(0,.2))
		hist(gff$dist_to_nearest_telo/1000000, breaks=100, freq=F, add=T, col=rgb(1,0,0,.2))
		legend(x=10, y=.2, legend=c(paste0("All genes, n=",nrow(gffComplete)), paste0("GC outliers, n=", nrow(gff))), fill=c(rgb(0,0,1,.2), rgb(1,0,0,.2)), cex=.7)
		#is the observed distribution different than the genome-wide dist? yes, resoundingly so.
		pval = t.test(x=gffComplete$dist_to_nearest_telo/1000000, y=gff$dist_to_nearest_telo/1000000)$p.value;pval #the log-transformmed is porbably better -  more normal-shaped data.
		text(x=85, y=0.15, paste0("p = ", round(pval, 6)), cex=.7)
		#wilcox.test(x=gffComplete$dist_to_nearest_telo/1000000, y=gff$dist_to_nearest_telo/1000000)
	}
}
plot_distNearestTelo()
plot_distNearestTelo(log=T)



reps = 1e6
dist = rep(NA, reps)
n = nrow(gff)
x = 1:nrow(gffComplete)
observed = mean(gff$dist_to_nearest_telo, na.rm=T);observed


for(i in 1:reps){
	if(i %%1000 == 0){
		cat(i/reps*100, "%\r")
	}
	sIDX = sample(x=x, size=n);sIDX
	dist[i] = mean(gffComplete$dist_to_nearest_telo[sIDX], na.rm=T);dist[i]
}
dist_nearestTelo = na.omit(dist)


plot_distNearestTelo_perms = function(){
	observed = mean(gff$dist_to_nearest_telo, na.rm=T);observed
	hist(dist_nearestTelo/1000000, breaks=50, main=c("Average distance to nearest telomere repeat", paste0(reps, " permutations")), xlab="Distance to the closest telomere repeat (Mbp)", xlim=c(10,20))
	abline(v=observed/1e6, col="red", lwd=2)
	pval = sum(dist_nearestTelo < observed) / length(dist_nearestTelo);pval#the proportion of measurments that are less than the observed out of the entire distances list gives us a p-value
	m = mean(dist_nearestTelo);m
	sigma = sd(dist_nearestTelo);sigma
	sds = abs(observed - m)/sigma
	text(x=18, y=27000*2, paste0("p = ", round(pval, 5)), cex=.7)
	text(x=18, y=24000*2, paste0("Mean = ", round(m/1e6, 2)), cex=.7)
	text(x=18, y=21000*2, paste0("observed = ", round(observed/1e6, 2)), cex=.7)
	text(x=18, y=18000*2, paste0(round(sds, 1), "sd from the mean"), cex=.7)
}
plot_distAlongCzomeArm_perms()
plot_distNearestTelo_perms()


save.image(file="get_outlier_genes_from_Roddys_data3.RData")
###################################################################################################
###################################################################################################
###################################################################################################

load("get_outlier_genes_from_Roddys_data3.RData")
#calculate Recombintation

thin_map<-function(map, bp_window=1000){
#keeps only one marker in each bp_window. I.e: only one of two markers closer than bp_window will be kept. 
	map<-map[order(map$ChrN, map$ChrPos),];head(map)
	map$drop=0;head(map)
	old_Chr = "NA"
	old_pos = -2000
	for (i in 1:nrow(map)){
		Chr = map$Chr[i];Chr
		if( Chr == old_Chr){
			pos = map$ChrPos[i];pos
			if((pos - old_pos) < bp_window){
				map$drop[i] = 1
			}else{
				old_pos = pos;old_pos
			}	
		}else{
			old_pos = -2000
		}		
		old_Chr = Chr;old_Chr
	}
	return(map[map$drop==0,])
}

map<-read.csv("../Meriones_chromonome_v6.1.geneticmap.csv");map$R = NA; map$Mbp = map$ChrPos/1000000;map<-merge(x=map, y=ChrN_df, by="Chr", all.x=T, all.y=F);map<-map[order(map$ChrN, map$ChrPos),];table(map$Chr);map<-thin_map(map, 1000);table(map$Chr);map_list<-split(map, f=map$Chr);#map_list


window_main = 8
for (i in 1:length(map_list)){
	Chr_df = map_list[[i]];Chr_df
	#slide a window across the df and regress the cM by the Mbp
	window = min(window_main, nrow(Chr_df)) #in case there are fewer than window_main markers on the chr
	
	#plot(x=Chr_df$ChrPos, y=Chr_df$cM, main=paste(map_list[[i]][1,1]))

	for(j in 1:nrow(Chr_df)){
		small_DF = Chr_df[j:(j+window-1),];small_DF
		small_DF = small_DF[rowSums(is.na(small_DF)) != ncol(small_DF),];small_DF #drops rows that have NAs in them
		if(nrow(small_DF)<3){
			#midPos = NA
			R=NA
		}else{
			#midPos = mean(small_DF$ChrPos, na.rm=T);
			mod = lm(small_DF$cM ~ small_DF$Mbp);mod
			R = mod$coefficients[2];R #gets the slope of the model, which is the Recombintaion rate
			x0 = small_DF$ChrPos[1];x0
			x0b = small_DF$ChrPos[2];x0b#this is the start of the next marker, used to get the non-overlapping window size
			x1 = small_DF$ChrPos[nrow(small_DF)];x1
			y0 = (x0 * R) + mod$coefficients[1];y0
			y1 = (x1 * R) + mod$coefficients[1];y1
			#segments(x0=x0, x1=x1, y1=y1, y0=y0)
		}
		#cat(R, "\n")
		map_list[[i]]$R[j] = R
		map_list[[i]]$OwindowStart[j] = x0 
		map_list[[i]]$OwindowEnd[j] = x1
		map_list[[i]]$NOwindowLen[j] = x0b - x0 #the distance between one marker and the next
		#map_list[[i]]$midPos[j] = midPos
	}
	for(k in 1:(nrow(map_list[[i]])-1)){#takes the average of the overlapping windows for each inter-marker space.
		map_list[[i]]$finalR[k] = mean(map_list[[i]]$R[map_list[[i]]$OwindowStart <= map_list[[i]]$ChrPos[k] & map_list[[i]]$OwindowEnd >= map_list[[i]]$ChrPos[k+1]], na.rm=T)
	}
}
map = do.call(rbind, map_list);map = map[order(map$ChrN, map$ChrPos),];tail(map)
map$finalR = ifelse(is.na(map$R), NA, map$finalR) #keeps windows with NA from having a finalR.

head(map) #OwindowSt* is for overlapping window: includes 8 markers. 
		  #NOwindow* is for non-overlapping window, the region between Pos[i] and Pos[i+1]
genomeR = weighted.mean(map$finalR, map$NOwindowLen, na.rm=T);genomeR
chrlen$aveR = NA

for(i in 1:23){
	Chr = chrlen$Chr[i]
	chrlen$aveR[i]=weighted.mean(map$finalR[map$Chr==Chr], map$NOwindowLen[map$Chr==Chr])
}
chrlen

rawHotspots = na.omit(map[map$finalR>5*genomeR,]);rawHotspots#every marker that is above the genome average - need to find each hotspot and get the range and a midpoint
#some chrs have multiple hotspots: always more that 1MB apart. so get range of each set of hotspot markers that are within 10MB of eachother

hotspots = data.frame("Chr"=rep(NA, 1000), "ChrN"=rep(NA, 1000), "Start"=rep(NA, 1000), "Stop"=rep(NA, 1000));hotspots
#for each marker in rawHotspots, if it is an outlier, grab it and check the next, if it's an outlier grab it and check the next etc. when one isn't, write out the hotspot and reset everything

Start = NULL
Stop = NULL
ChrN = NULL
prevChr="Chr1"
prevFlag = FALSE
count = 0
for (i in 1:nrow(map)){
	Chr = map$Chr[i];Chr
	cat(Chr, "\n")
	if (Chr != prevChr){#reset everything for each chromsome
		Start = NULL
		Stop = NULL
		ChrN = NULL
	}
	if (! is.na(map$finalR[i]) & map$finalR[i] > 5*genomeR){#if a marker exceedes the genomeR, grab it and add it into the hotspot
		Start = min(Start, map$ChrPos[i]);Start
		Stop = max(Stop, map$ChrPos[i] + map$NOwindowLen[i]);Stop
		ChrN = map$ChrN[i];ChrN
		flag = TRUE;flag
	}else{
		flag = FALSE
	}
	if( !(flag) & prevFlag){#meaning we just left a hotspot - the current marker is not a hotspot.
		count = count + 1;count
		hotspots$Chr[count] = Chr
		hotspots$ChrN[count] = ChrN
		hotspots$Start[count] = Start	
		hotspots$Stop[count] = Stop
		Start = NULL
		Stop = NULL
		ChrN = NULL
	}
	
	prevFlag = flag;prevFlag
	prevChr = Chr;prevChr
}

i=i+1

hotspots=na.omit(hotspots);
hotspots


par(mfrow=c(5,5))
for(Chr in chrlen$Chr[1:22]){ 
#plot(x=map$Mbp[map$Chr==Chr], y=map$R[map$Chr==Chr], xlim=c(0,chrlen$Len[chrlen$Chr==Chr]/1e6), ylim=c(0,150), type="p", pch=19, cex=.5, main=Chr, xlab="Mbp", ylab="R (cM/Mbp)");
plot(x=map$Mbp[map$Chr==Chr], y=map$cM[map$Chr==Chr], xlim=c(0,chrlen$Len[chrlen$Chr==Chr]/1e6), type="l", main=Chr, ylab="cM", xlab="Mbp")
h = hotspots[hotspots$Chr==Chr,]
if(nrow(h)>0){
rect(xleft=h$Start/1e6, xright=h$Stop/1e6, ytop=500, ybottom=0, col=rgb(1,0,0, .4), border=NA)}
}



#now the df called map:
head(map) # has a column called "finalR" which is the R between the ChrPos and the ChrPos of the next row (which is not windowEnd! - windowEnd is overlapping window endpoint - 7 markers down the line...)

hist(map$finalR)#this is just th edifferent finalR's - it does not weight them by their length, so not to be used for finding the 95% CI. 
#Chr,Pos,R for every pos in the genome

#get weighted average R: for genome and then for each chr


plot_genome_summary = function(){
barplot(chrlen$Len[1:24]/1000000, names=chrlen$Chr[1:24], col="white", main="Meriones genome", ylab="Position on chromosome (Mbp)", ylim=c(0,220), cex.axis = .8, cex.names=.7, las=2);


#genes
segments(y0=alldata$start/1000000, x0=(r*alldata$ChrN)-s, y1=alldata$start/1000000, x1=(r*(alldata$ChrN+.83)-s), col=rgb(1, .3, .4, 1))
#telomeres
segments(y0=telomeres$Start/1000000, x0=(r* telomeres$ChrN)-s, y1=telomeres$Start/1000000, x1=(r*(telomeres$ChrN+.83)-s), col=rgb(.1,.8,.4,1))
#centromeres
rect(ybottom=chrlen$centStart/1000000, xleft=(r*chrlen$ChrN)-s, ytop=chrlen$centStop/1000000, xright=(r*(chrlen$ChrN+.83)-s), col=rgb(0,0,1,.1), border=rgb(0,0,1,1))


#recombintaion hotspots
xleft = (r*(hotspots$ChrN+.83)-s)
xright = (r* hotspots$ChrN)-s
ybottom = hotspots$Start/1000000
ytop = hotspots$Stop/1000000
rect(xleft=xleft, xright = xright, ybottom= ybottom, ytop = ytop, col=rgb(0,.5,0.5,1))

#recombination rate
points(x=((r*map$ChrN)-s)+(map$R*1/max(map$R, na.rm=T)), y=map$Mbp, pch=19, cex=.5, col="black")#as.factor(map$ChrN))

text(chrlen$Freq, x=(r*chrlen$ChrN)-s/2, y=210)
legend(x=11, y=190, legend=c("GC outlier genes", "Centromeres", "Telomere repeats (>70x copies)", "Recombination hotspots"), fill=c(rgb(1, .3, .4, 1), rgb(0,0, 1,.1), rgb(.1,.8,.4,1), rgb(0, .5, 0.5, 1)))

}
plot_genome_summary()



####need to identify R hotspots and then do the permutation test for gc genes near R hotspots.

hotspots#has the list of hotspots.
#now do the permutations.

calc_dist_to_nearest_hotspot<-function(gff_entry){
	#extract gene location in bp
	Chr = gff_entry[1];Chr
	if(! Chr %in% unique(hotspots$Chr)){return(NA)}#drop out things on that don't have hotspots or genetic markers
	Start = as.numeric(gff_entry[2]);Start
	
	
	#find closest telomere
	hs = hotspots[hotspots$Chr==Chr,];hs
	old_dist = 100000000000000
	for(i in 1:nrow(hs)){
		if(Start<hs$Stop[i] & Start>hs$Start[i]){#it's within the hotspot
			return(1)
		} 
		distStart = abs(hs$Start[i]- Start)
		distStop = abs(hs$Stop[i] - Start)
		dist = min(distStart, distStop);dist
		if(dist < old_dist){#as long as you're getting closer to a hs, keep going, once it's getting further away, then bail out.
			old_dist=dist;dist
		}else{
			break
		}
	}
	return(old_dist)
}


gffComplete$dist_to_nearest_Rhs = apply(X=gffComplete, MARGIN=1, FUN=calc_dist_to_nearest_hotspot);head(gffComplete)
gff$dist_to_nearest_Rhs=apply(X=gff, MARGIN=1, FUN=calc_dist_to_nearest_hotspot);head(gff)




plot_DistRHotspot = function(log=FALSE){
	if(log){
		hist(x=log10(gff$dist_to_nearest_Rhs), breaks=seq(0, 8.5, by=.1), main=c("Distance to nearest recombination hotspot"), xlab="Distance to nearest recombination hotspot (log10 bp)", freq=F, col=rgb(0,0,1,.2));hist(log10(gffComplete$dist_to_nearest_Rhs), breaks=seq(0, 8.5, by=.1), freq=F, add=T, col=rgb(1,0,0,.2))
		pval = t.test(x=log10(gffComplete$dist_to_nearest_Rhs), y=log10(gff$dist_to_nearest_Rhs))$p.value;pval
		legend(x=.01, y=1, legend=c(paste0("All genes, n=",nrow(gffComplete)), paste0("GC outliers, n=", nrow(gff))), fill=c(rgb(0,0,1,.2), rgb(1,0,0,.2)), cex=.7)
		text(x=5, y=.7, paste0("p = ", round(pval, 6)), cex=.7)
	}else{
		hist(x=gffComplete$dist_to_nearest_Rhs/1e6, breaks=100, main="Distance to nearest recombination hotspot", xlab="Distance to nearest recombination hotspot (Mbp)", freq=F, col=rgb(0,0,1,.2), ylim = c(0,.1))
		hist(gff$dist_to_nearest_Rhs/1e6, breaks=100, freq=F, add=T, col=rgb(1,0,0,.2))
		legend(x=15, y=.1, legend=c(paste0("All genes, n=",nrow(gffComplete)), paste0("GC outliers, n=", nrow(gff))), fill=c(rgb(0,0,1,.2), rgb(1,0,0,.2)), cex=.7)
		#is the observed distribution different than the genome-wide dist? yes, resoundingly so.
		pval = t.test(x=gffComplete$dist_to_nearest_Rhs/1e6, y=gff$dist_to_nearest_Rhs/1e6)$p.value;pval
		text(x=100, y=0.08, paste0("p = ", round(pval, 6)), cex=.7)
		#wilcox.test(x=gffComplete$dist_to_nearest_Rhs/1000000, y=gff$dist_to_nearest_Rhs/1000000)		
	}
}




reps = 1e6
dist = rep(NA, reps)
n = nrow(gff)
x = 1:nrow(gffComplete)
observed = mean(gff$dist_to_nearest_Rhs, na.rm=T);observed


for(i in 1:reps){
	if(i %%1000 == 0){
		cat(i/reps*100, "%\r")
	}
	sIDX = sample(x=x, size=n);sIDX
	dist[i] = mean(gffComplete$dist_to_nearest_Rhs[sIDX], na.rm=T);dist[i]
}
dist_nearestRhs = na.omit(dist)

plot_DistRHotspot_perms = function(){
	par(cex=1)
	observed = mean(gff$dist_to_nearest_Rhs, na.rm=T);observed
	hist(dist_nearestRhs/1000000, breaks=100, main=c(paste0("Average distance to nearest \nrecombination hotspot\n", reps, " permutations")), xlab="Distance to nearest recombination hotspot (Mbp)", xlim=c(18,38))
	abline(v=observed/1e6, col="red", lwd=2)
	pval = sum(dist_nearestRhs < observed) / length(dist_nearestRhs);pval#the proportion of measurments that are less than the observed out of the entire distances list gives us a p-value
	m = mean(dist_nearestRhs);m
	sigma = sd(dist_nearestRhs);sigma
	sds = abs(observed - m)/sigma
	
	text(x=35, y=40000, paste0("p = ", round(pval, 5)), cex=.7)
	text(x=35, y=35000, paste0("Mean = ", round(m/1e6, 2)), cex=.7)
	text(x=35, y=30000, paste0("observed = ", round(observed/1e6, 2)), cex=.7)
	text(x=35, y=25000, paste0(round(sds, 1), "sd from the mean"), cex=.7)
	par(cex=1)
}
plot_distNearestTelo_perms()

plot_DistRHotspot_perms()



plot_RecomHS_dist = function(){
	RHSdat = as.data.frame(table(hotspots$Chr));dim(RHSdat)
	RHSdat$Var1 = as.character(RHSdat$Var1)
	RHSdat = rbind.data.frame(RHSdat, c("Chr2", 0));RHSdat
	RHSdat = rbind.data.frame(RHSdat, c("Chr18", 0));RHSdat
	RHSdat = rbind.data.frame(RHSdat, c("Chr21", 0));RHSdat
	RHSdat = rbind.data.frame(RHSdat, c("ChrX", 0));RHSdat
	RHSdat$Freq = as.numeric(RHSdat$Freq);str(RHSdat)
	
	
	hist(RHSdat$Freq, breaks=c(-1:10), axes=F, xlab="Number of recombination hotspots", ylab="Number of chromsomes", ylim=c(0,8), main="Recombination hotspots per chromsome")	
	axis(side=1, at = 0:10-.5, labels=0:10)
	axis(side=2, at = seq(0,8, by=2), labels=seq(0,8, by=2))
	text(x=7, y=7.5, paste0("N = 22"))
	text(x=7, y=6, paste0("Mean = ", round(mean(RHSdat$Freq), 2)))
	text(x=7, y=4.5, paste0("StDev = ", round(sd(RHSdat$Freq), 2)))
	

}

plot_RecomHS_dist()


save.image(file="get_outlier_genes_from_Roddys_data4.RData")






setwd("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/GC_outlier_genes/")
load("get_outlier_genes_from_Roddys_data4.RData")









########################################################################################################################
########################################################################################################################
########################################################################################################################
#here is the main figure  - all pulled from above. 
par(mfrow=c(4,2), cex=1, cex.main=1.01)

########################################################################################################################
#Panel 0: Genome overveiw with GC genes, Recombination rates, intersitial telomeres, and centromeres annotated.
#plot_genome_summary() #needs a full-length plot, one little square is too little.

########################################################################################################################
#Panel 1: Are GC genes clustered? Yes, binomial dist with observed
#plot(x=0, y=1, type="n", bty="n", axes=FALSE, xlab="", ylab="", main="")
plot_distNearestNeighbor_perms()
mtext("A", side=3, line=2, at= 1, cex=2)

########################################################################################################################
#Panel 2: Where are the recombination hotspots?
plot_RecomHS_dist()
mtext("B", side=3, line=2, at= -4, cex=2)
########################################################################################################################
#Panel 4: Are GC genes near recombination hotspots? Yes
plot_DistRHotspot(log=TRUE)
mtext("C", side=3, line=2, at= -2, cex=2)
plot_DistRHotspot_perms()
mtext("D", side=3, line=2, at= 12, cex=2)

########################################################################################################################
#Panel 2: Are GC genes near centromeres? Yes 
plot_distAlongCzomeArm()
mtext("E", side=3, line=2, at= -24, cex=2)

plot_distAlongCzomeArm_perms()
mtext("F", side=3, line=2, at= -28, cex=2)

########################################################################################################################
#Panel 3: Are GC genes near telomere sites (interstital or otherwise)? Yes
plot_distNearestTelo(log=TRUE)
mtext("G", side=3, line=2, adj=-.22, cex=2)
plot_distNearestTelo_perms()
mtext("H", side=3, line=2, at= 7, cex=2)
















