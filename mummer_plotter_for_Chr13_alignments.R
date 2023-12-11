#written 20 May 2022 by TDB
#plots out coords file of the Chr13 mummer alignments. 
library(bedr)

#nucmer --maxmatch --nosimplify -t 2 -p Chr13ScafsvChr13.maxmatch.nosimplify ../Chr13_repeat_structure/Meriones_chromonome_v6.1.AllChr13Scaffolds.fasta ../v6.1_seqs/Chr13.v6.1.fasta

#delta-filter -l 20000 Chr13ScafsvChr13.maxmatch.nosimplify.delta > Chr13ScafsvChr13.maxmatch.nosimplify.delta.l20000.filter
#delta-filter -l 10000 Chr13ScafsvChr13.maxmatch.nosimplify.delta > Chr13ScafsvChr13.maxmatch.nosimplify.delta.l10000.filter
#delta-filter -l 1000 Chr13ScafsvChr13.maxmatch.nosimplify.delta > Chr13ScafsvChr13.maxmatch.nosimplify.delta.l1000.filter

#show-coords -Tl Chr13ScafsvChr13.maxmatch.nosimplify.delta.l20000.filter > Chr13ScafsvChr13.maxmatch.nosimplify.delta.l20000.filter.Tl.coords
#show-coords -Tl Chr13ScafsvChr13.maxmatch.nosimplify.delta.l10000.filter > Chr13ScafsvChr13.maxmatch.nosimplify.delta.l10000.filter.Tl.coords
#show-coords -Tl Chr13ScafsvChr13.maxmatch.nosimplify.delta.l1000.filter > Chr13ScafsvChr13.maxmatch.nosimplify.delta.l1000.filter.Tl.coords



#for S in ../v6.1_seqs/Chr13*fasta
#do
#	Q=${S/.v6.1.fasta/}; 
#	Q=${Q/..\/v6.1_seqs\//}; 
#	run_cmd.sh "nucmer --maxmatch --nosimplify -t 2 -p Chr13Scafsv$Q.maxmatch.nosimplify ../Chr13_repeat_structure/Meriones_chromonome_v6.1.AllChr13Scaffolds.fasta $S"  Chr13Scafsv$Q.maxmatch.nosimplify.delta;
#	run_cmd.sh "delta-filter -l 1000 Chr13Scafsv$Q.maxmatch.nosimplify.delta > Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter" Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter; 
#	run_cmd.sh "show-coords -Tl Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter > Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter.Tl.coords" Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter.Tl.coords;
#done
#run_cmd.sh "cat Chr13Scafsv*.maxmatch.nosimplify.delta.l1000.filter.Tl.coords > Chr13ScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT" Chr13ScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT



get_chr_length = function(ChrID, coord=coords){
	if( ! is.na(coord$RLen[coord$Ref==ChrID][1])){
		return(coord$RLen[coord$Ref==ChrID][1])
	}
	return(coord$QLen[coord$Query==ChrID][1])
}

plot_mummer <-function(coords, REF, QUE ,IDYthreshold=0, LENthreshold = 1000, xlimMin = 0, xlimMax=get_chr_length(REF, coord=coords), ylimMin=0, ylimMax=get_chr_length(QUE, coord=coords)){
	REFLEN = xlimMax - xlimMin #get_chr_length(REF, coord=coords)
	QUELEN = ylimMax - ylimMin #get_chr_length(QUE, coord=coords)
	#maxplotLim = max(REFLEN, QUELEN)
	
	coords = coords[coords$Ref==REF & coords$Query==QUE & coords$IDY >= IDYthreshold & coords$QAliLen >= LENthreshold, ]
	coords = coords[coords$RStart>(xlimMin*0.8) & coords$QStart>(ylimMin*0.8) & coords$RStop<(xlimMax*1.2) & coords$QStop<(ylimMax*1.2) ,]
	
	if (REF ==QUE){
		#mainTitle = paste("Self-alignment of", REF)
		mainTitle = paste(REF)
	}else{
		mainTitle = paste("Alignments between", REF, "and", QUE)
	}
	xlab = c(paste0("Position in Mbp"))#, paste0("total length = ", round(get_chr_length(REF, coord=coords)/1e6,1), "Mbp"))
	ylab = c(paste0("Position in Mbp"))#, paste0("total length = ", round(get_chr_length(QUE, coord=coords)/1e6,1), "Mbp"))
	plot(x=0, y=0, xlim=c(xlimMin/1e6, xlimMax/1e6), ylim=c(ylimMin/1e6, ylimMax/1e6), type="n", bty="L", xlab=xlab, ylab=ylab, axes=F)
	title(main=mainTitle)
	#title(main=c("","","","", "Grey box denotes the length of the scaffolds in the X and Y directions", paste0("Minimum alignment identity: ", IDYthreshold, "%"), paste0("Minimum alignment length: ", LENthreshold, "bp")), cex.main=.8)
	axis(side=1, at = (c(0/1e6, round(REFLEN/2e6, 1), round(2*REFLEN/2e6, 1)) + xlimMin/1e6))
	axis(side=2, at = (c(0/1e6, round(QUELEN/2e6, 1), round(2*QUELEN/2e6, 1)) + ylimMin/1e6))
	#rect(xleft = 0, ybottom = 0 , xright = get_chr_length(REF, coord=coords), ytop = get_chr_length(QUE, coord=coords), col=rgb(0,0,0,0.05), lwd=0)
	x0 = coords$RStart/1e6
	y0 = coords$QStart/1e6
	x1 = coords$RStop/1e6
	y1 = coords$QStop/1e6
	segments(x0 = x0, y0 = y0, x1 = x1, y1 = y1)
}


read_in_coords<-function(fileHandle, drop=FALSE){
	coords<-read.table(fileHandle, header=F, skip=4, sep="\t")
	colnames(coords) = c("QStart", "QStop", "RStart", "RStop", "QAliLen", "RAliLen", "IDY", "QLen" ,"RLen", "Query", "Ref");
	#need to get rid of 1:1 alignments when calculating duplication lelves, but not when plotting self-alingmnets - those aren't repetitve, they're just the self-align. But this was moved to the get_merged_beds bit, so drop should really always be true ehre
	if(drop){coords = coords[!(coords$QStart == coords$RStart & coords$QStop == coords$RStop & coords$Query == coords$Ref),]}
	return(coords)
}

get_merged_bed_of_repeats<-function(coords){
	coords = coords[!(coords$QStart == coords$RStart & coords$QStop == coords$RStop & coords$Query == coords$Ref),]
	bed<-cbind.data.frame("chr"=coords$Query, "start"=coords$QStart, "end"=coords$QStop)
	bed.sorted<-bedr.sort.region(bed, check.chr=F)
	bed.sorted.merged<-bedr.merge.region(bed.sorted, check.chr=F)
	return(bed.sorted.merged)
}

get_redundant_bits_of_chr<-function(bed.merged, ChrLen, minLenThreshold=0){
	bed.merged$Len = abs(bed.merged$end - bed.merged$start);head(bed.merged)
	repeatLen = sum(bed.merged$Len[bed.merged$Len>minLenThreshold]);repeatLen
	return(repeatLen / ChrLen *100)
}



##################################################################
###########using the "canonical" repeat to fish out percent of the genome that is that repeat
#here is the 170kb fragment, self-aligned:
coords170self<-read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Chr13_repeat_structure/Chr13_170kb_repeat_v_SELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
plot_mummer(coords170self, "Chr13:31508-204060", "Chr13:31508-204060", 0, 1000)

#here is the rest of Chr13 - will be useful to have (see below for teh bash code to prepare this file.):
coords13 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr13ScafsvChr13various/Chr13ScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT", drop=FALSE)
bed13.merged<-get_merged_bed_of_repeats(coords13)
head(bed13.merged)

#the scaffold lengths will also be useful.
LenDF = read.table("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.chrLengths.csv", sep=",", header=T);head(LenDF)

#here's the entire length of chr13:
Chr13Len = sum(c(2679916, 17095644, 16671041, 9569972, 7902876, 6143218, 5705833, 4935230, 3650033, 2883794, 2347156, 1703228, 1573360, 1369987, 1209380, 925606, 914264, 908068, 853260, 835279, 775897, 701984, 701463, 669220, 575929, 533563, 532724, 496881, 472376, 470050, 457961, 434242, 392047, 384531, 371784, 332107, 314591, 310064, 294655, 270248, 260653, 259644, 239478, 237422, 237028, 206944, 205817, 203597, 203151, 193665, 188782, 188145, 177885, 175620, 174130, 172853, 170959, 167363, 166017, 164385, 161338, 150568, 149329, 144469, 127471, 120360, 117036, 114300, 112281, 105674, 100321, 99644, 98203, 94581, 93246, 90409, 90357, 90277, 89034, 87185, 86295, 85060, 84320, 84268, 84227, 81987, 80657, 80184, 78313, 76513, 74806, 71733, 70814, 68359, 67811, 66952, 60814, 60274, 59817, 54774, 54696, 52703, 50854, 47742, 45263, 43915, 43673, 42036, 38234, 31296, 31105, 29840, 29710, 28696, 28403, 24856, 24727, 23870, 23855, 23044, 22393));Chr13Len



#the 170kb fragment aligned to the chr13scaffolds:
	coords170<-read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr13ScafsvChr13_170kb_repeat.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
	bed170.merged<-get_merged_bed_of_repeats(coords170)

#39% of chr13 is comprised of these 170kb repeats:
	get_redundant_bits_of_chr(bed170.merged, Chr13Len) 

#but 93% of the chromosome is comprised of repeated bits of elsewhere in the chromosome:
	get_redundant_bits_of_chr(bed13.merged, Chr13Len) 

#so I must be missing something - some repeat that is also at high frequency. 


#here is another way to visualise it: 
	par(mfrow=c(1,2))
	plot_mummer(coords170, "Chr13:31508-204060", "Chr13_unordered_Scaffold_72", 0, 1000)
#the 170kb is from bases 83,000-100,000, and only occurs once in the S72
	coordsChr13S72<-read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr13ScafsvChr13various/Chr13ScafsvChr13_unordered_Scaffold_72.maxmatch.nosimplify.delta.l1000.filter.TL.coords", drop=FALSE)
	plot_mummer(coordsChr13S72, "Chr13_unordered_Scaffold_72", "Chr13_unordered_Scaffold_72", 0, 1000)
#but S72 is very self-repetitive: 76% is elsewhere on S72:
	coordsChr13S72only <- coordsChr13S72[coordsChr13S72$Query == "Chr13_unordered_Scaffold_72",];dim(coordsChr13S72only)
	bedChr13S72only.merged<-get_merged_bed_of_repeats(coordsChr13S72only);
	get_redundant_bits_of_chr(bedChr13S72only.merged, LenDF$Pos[LenDF$Chr=="Chr13_unordered_Scaffold_72"])
#and 99% of it is redundant with something else on chr13:
	get_redundant_bits_of_chr(bed13.merged[bed13.merged$chr=="Chr13_unordered_Scaffold_72",], LenDF$Pos[LenDF$Chr=="Chr13_unordered_Scaffold_72"])
#so I must be missing something, and I think I can get it by looking closesly at the repeat in S72: 
	plot_mummer(coordsChr13S72, "Chr13_unordered_Scaffold_72", "Chr13_unordered_Scaffold_72", 0, 1000);abline(v=163000);abline(v=187900);

#starts at 163000, ends at 187900
189700-163000 #=26700

#extract that and do the mummer align:
#extract_fasta_region.py -f ../Meriones_chromonome_v6.1.fasta -chr Chr13_unordered_Scaffold_72 -s 163000 -e 1879000 -o Chr13_26kb_repeat.fasta
#nucmer --maxmatch --nosimplify -t 2 -p Chr13ScafsvChr13_26kb_repeat.maxmatch.nosimplify. Meriones_chromonome_v6.1.AllChr13Scaffolds.fasta Chr13_26kb_repeat.fasta 
#delta-filter -l 1000 Chr13ScafsvChr13_26kb_repeat.maxmatch.nosimplify.delta > Chr13ScafsvChr13_26kb_repeat.maxmatch.nosimplify.delta.l1000.filter
#show-coords -Tl Chr13ScafsvChr13_26kb_repeat.maxmatch.nosimplify.delta.l1000.filter > Chr13ScafsvChr13_26kb_repeat.maxmatch.nosimplify.delta.l1000.filter.Tl.coords


coords26 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Chr13_repeat_structure/Chr13ScafsvChr13_26kb_repeat.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
plot_mummer(coords26, "Chr13_unordered_Scaffold_72:163000-187900", "Chr13_unordered_Scaffold_72", 0, 1000)
head(coords26);dim(coords26)
bed26.merged = get_merged_bed_of_repeats(coords26)
get_redundant_bits_of_chr(bed26.merged, Chr13Len)
get_redundant_bits_of_chr(bed26.merged[bed26.merged$chr=="Chr13_unordered_Scaffold_72",],LenDF$Pos[LenDF$Chr=="Chr13_unordered_Scaffold_72"])
#so 40% of S72 is this 26k repeat, but only 3% of all of Chr13. So that's not a major explanation. Try again





#try  Chr13_unordered_Scaffold_108 - it is 99% repetitive, but only 16% is the 170kb and 0% is the 26kb
	plot_mummer(coords170, "Chr13:31508-204060", "Chr13_unordered_Scaffold_108", 0, 1000) #not much is the 170kb

	coordsChr13S108<-read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr13ScafsvChr13various/Chr13ScafsvChr13_unordered_Scaffold_108.maxmatch.nosimplify.delta.l1000.filter.TL.coords", drop=FALSE)
	plot_mummer(coordsChr13S108, "Chr13_unordered_Scaffold_108", "Chr13_unordered_Scaffold_108", 0, 1000)
#but S108 is NOT very self-repetitive: 5% is elsewhere on S108:
	coordsChr13S108only <- coordsChr13S108[coordsChr13S108$Query == "Chr13_unordered_Scaffold_108",];dim(coordsChr13S108only)
	bedChr13S108only.merged<-get_merged_bed_of_repeats(coordsChr13S108only);
	get_redundant_bits_of_chr(bedChr13S108only.merged, LenDF$Pos[LenDF$Chr=="Chr13_unordered_Scaffold_108"])
#and 99% of it is redundant with something else on chr13:
	get_redundant_bits_of_chr(bed13.merged[bed13.merged$chr=="Chr13_unordered_Scaffold_108",], LenDF$Pos[LenDF$Chr=="Chr13_unordered_Scaffold_108"])
#13% is the amount of S108 that appears on chr13 (incudling itself: S108 is up 0.1% of all of Chr13): 
	coordsChr13S108asRef <- coordsChr13S108[coordsChr13S108$Ref == "Chr13_unordered_Scaffold_108",];dim(coordsChr13S108asRef)
	bedChr13S108asRef.merged<-get_merged_bed_of_repeats(coordsChr13S108asRef);
	get_redundant_bits_of_chr(bedChr13S108asRef.merged, Chr13Len)



#so I must be missing something, and I think I can get it by looking closesly at the repeat in S108: 
	plot_mummer(coords13, "Chr13_unordered_Scaffold_108", "Chr13_unordered_Scaffold_28", 0, 1000, ylimMax = 4167760, ylimMin = 4167760-2000000);abline(v=163000);abline(v=187900);abline(v=c(29000,43500))


	plot_mummer(coords13, "Chr13_unordered_Scaffold_108", "Chr13_unordered_Scaffold_28", 0, 1000);#abline(h= c(4167760,4167760-2000000));abline(v=c(31000,40000))




#not worth continuing



#by scaffold now instead of the whole chr:
bed_list = split(bed170.merged, bed170.merged$chr)
redundant_by_scaf = data.frame("Chr"=rep(NA, length(bed_list)), "Len"=rep(NA,  length(bed_list)), "Repetitive"=rep(NA,  length(bed_list)),"R170k"=rep(NA,  length(bed_list)))
i=0
for (bit in bed_list){
	i=i+1
	redundant_by_scaf$Chr[i] = bit$chr[1]
	redundant_by_scaf$Len[i] = LenDF$Pos[LenDF$Chr==bit$chr[1]]
	redundant_by_scaf$Repetitive[i] = get_redundant_bits_of_chr(bed13.merged[bed13.merged$chr==bit$chr[1],], LenDF$Pos[LenDF$Chr==bit$chr[1]])	
	redundant_by_scaf$R170k[i] = get_redundant_bits_of_chr(bit, LenDF$Pos[LenDF$Chr==bit$chr[1]])
}
redundant_by_scaf = redundant_by_scaf[order(redundant_by_scaf$R170k),]
redundant_by_scaf

plot(x=log10(redundant_by_scaf$Len), y=redundant_by_scaf$R170k)
plot(x=redundant_by_scaf$Repetitive, y=redundant_by_scaf$R170k)#nearly everything is 100% repeated, but it ranges from 0-100% are copies of the 170kb bit.





dim(coordsChr13S72)



##################################################################
###########How much of Chr13 is made up of repeats of anywhere else on Chr13? Any other repeat is sufficient: i.e. what is 1-(unique bits of Chr13)?

#here's the bash code to get all this ready:
#for S in ../v6.1_seqs/Chr13_unordered_Scaffold_*fasta;
# do Q=${S/.v6.1.fasta/};  
#Q=${Q/..\/v6.1_seqs\//};  
#run_cmd.sh "nucmer --maxmatch --nosimplify -t 2 -p Chr13Scafsv$Q.maxmatch.nosimplify ../Chr13_repeat_structure/Meriones_chromonome_v6.1.AllChr13Scaffolds.fasta $S"  Chr13Scafsv$Q.maxmatch.nosimplify.delta;
# run_cmd.sh "delta-filter -l 1000 Chr13Scafsv$Q.maxmatch.nosimplify.delta > Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter" Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter;
#  run_cmd.sh "show-coords -Tl Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter > Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter.Tl.coords" Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter.Tl.coords; 
#done; 

#for F in Chr13ScafsvChr13*.maxmatch.nosimplify.delta.l1000.filter.Tl.coords
#do
#run_cmd.sh "awk 'NR >4 { print }' $F >> Chr13ScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT" 1
#done

#see below for all of this


##################################################################
###########So how does that compare to other chromosomes: 
###########How much of Chr1 is made up of repeats of anywhere else on Chr1, etc.? 




#coords1 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr1.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords2 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr2.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords3 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr3.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords4 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr4.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords5 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr5.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords6 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr6.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords7 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr7.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords8 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr8.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords9 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr9.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coords10 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr10.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coords11 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr11.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coords12 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr12.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coords13 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr13ScafsvChr13various/Chr13ScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT", drop=FALSE)
#coords14 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr14.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coords15 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr15.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coords16 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr16.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords17 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr17.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords18 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr18.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords19 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr19.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coords20 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr20.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coords21 = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/Chr21.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
#coordsX = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/ChrX.vSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords", drop=FALSE)
coordsY = read_in_coords("/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/mummer/ChrYScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT", drop=FALSE)


bed1.merged<-get_merged_bed_of_repeats(coords1)
bed2.merged<-get_merged_bed_of_repeats(coords2)
bed3.merged<-get_merged_bed_of_repeats(coords3)
bed4.merged<-get_merged_bed_of_repeats(coords4)
bed5.merged<-get_merged_bed_of_repeats(coords5)
bed6.merged<-get_merged_bed_of_repeats(coords6)
bed7.merged<-get_merged_bed_of_repeats(coords7)
bed8.merged<-get_merged_bed_of_repeats(coords8)
bed9.merged<-get_merged_bed_of_repeats(coords9)
bed10.merged<-get_merged_bed_of_repeats(coords10)
bed11.merged<-get_merged_bed_of_repeats(coords11)
bed12.merged<-get_merged_bed_of_repeats(coords12)
bed13.merged<-get_merged_bed_of_repeats(coords13)
bed14.merged<-get_merged_bed_of_repeats(coords14)
bed15.merged<-get_merged_bed_of_repeats(coords15)
bed16.merged<-get_merged_bed_of_repeats(coords16)
bed17.merged<-get_merged_bed_of_repeats(coords17)
bed18.merged<-get_merged_bed_of_repeats(coords18)
bed19.merged<-get_merged_bed_of_repeats(coords19)
bed20.merged<-get_merged_bed_of_repeats(coords20)
bed21.merged<-get_merged_bed_of_repeats(coords21)
bedX.merged<-get_merged_bed_of_repeats(coordsX)
bedY.merged<-get_merged_bed_of_repeats(coordsY)


get_redundant_bits_of_chr(bed1.merged, coords1$RLen[1])
get_redundant_bits_of_chr(bed2.merged, coords2$RLen[1])
get_redundant_bits_of_chr(bed3.merged, coords3$RLen[1])
get_redundant_bits_of_chr(bed4.merged, coords4$RLen[1])
get_redundant_bits_of_chr(bed5.merged, coords5$RLen[1])
get_redundant_bits_of_chr(bed6.merged, coords6$RLen[1])
get_redundant_bits_of_chr(bed7.merged, coords7$RLen[1])
get_redundant_bits_of_chr(bed8.merged, coords8$RLen[1])
get_redundant_bits_of_chr(bed9.merged, coords9$RLen[1])
get_redundant_bits_of_chr(bed10.merged, coords10$RLen[1], 1000)
get_redundant_bits_of_chr(bed11.merged, coords11$RLen[1])
get_redundant_bits_of_chr(bed12.merged, coords12$RLen[1])
get_redundant_bits_of_chr(bed13.merged, 107817972, 1000)
get_redundant_bits_of_chr(bed14.merged, coords14$RLen[1])
get_redundant_bits_of_chr(bed15.merged, coords15$RLen[1])
get_redundant_bits_of_chr(bed16.merged, coords16$RLen[1])
get_redundant_bits_of_chr(bed17.merged, coords17$RLen[1])
get_redundant_bits_of_chr(bed18.merged, coords18$RLen[1])
get_redundant_bits_of_chr(bed19.merged, coords19$RLen[1])
get_redundant_bits_of_chr(bed20.merged, coords20$RLen[1])
get_redundant_bits_of_chr(bed21.merged, coords21$RLen[1])
get_redundant_bits_of_chr(bedX.merged, coordsX$RLen[1])
get_redundant_bits_of_chr(bedY.merged, 4450620+8639947+11283957+27601515+23978+35123045)




par(mfrow=c(3,3))
par(cex=.8)
plot_mummer(coords10, "Chr10", "Chr10", 0, 1000, xlimMax = 38e6, ylimMax = 38e6, xlimMin=36e6, ylimMin=36e6)
plot_mummer(coords16, "Chr16", "Chr16", 0, 1000, xlimMax = 12e6, ylimMax = 12e6, xlimMin=10e6, ylimMin=10e6)
plot_mummer(coords21, "Chr21", "Chr21", 0, 1000, xlimMax = 18e6, ylimMax = 18e6, xlimMin=16e6, ylimMin=16e6)
#plot_mummer(coords10, "Chr10", "Chr10", 0, 1000, xlimMax = 32.5e6, ylimMax = 32.5e6, xlimMin=30e6, ylimMin=30e6)
#plot_mummer(coords10, "Chr10", "Chr10", 0, 1000, xlimMax = 102.5e6, ylimMax = 102.5e6, xlimMin=100e6, ylimMin=100e6)
#plot_mummer(coords11, "Chr11", "Chr11", 0, 1000, xlimMax = 2e6, ylimMax = 2e6, xlimMin=0e6, ylimMin=0e6)
#plot_mummer(coords11, "Chr11", "Chr11", 0, 1000, xlimMax = 32.5e6, ylimMax = 32.5e6, xlimMin=30e6, ylimMin=30e6)
#plot_mummer(coords11, "Chr11", "Chr11", 0, 1000, xlimMax = 102.5e6, ylimMax = 102.5e6, xlimMin=100e6, ylimMin=100e6)
#plot_mummer(coords12, "Chr12", "Chr12", 0, 1000, xlimMax = 52.5e6, ylimMax = 52.5e6, xlimMin=50e6, ylimMin=50e6)
plot_mummer(coords13, "Chr13", "Chr13", 0, 1000, xlimMax = 2e6, ylimMax = 2e6)
plot_mummer(coords13, "Chr13_unordered_Scaffold_28", "Chr13_unordered_Scaffold_28", 0, 1000, xlimMax = 2e6, ylimMax = 2e6)
plot_mummer(coords13, "Chr13_unordered_Scaffold_33", "Chr13_unordered_Scaffold_33", 0, 1000, xlimMax = 7e6, ylimMax = 7e6, xlimMin=5e6, ylimMin=5e6)
#plot_mummer(coords13, "Chr13_unordered_Scaffold_37", "Chr13_unordered_Scaffold_37", 0, 1000, xlimMax = 2.5e6, ylimMax = 2.5e6)
#plot_mummer(coords15, "Chr15", "Chr15", 0, 1000, xlimMax = 27e6, ylimMax = 27e6, xlimMin=24e6, ylimMin=24e6)
plot_mummer(coordsY, "ChrY_unordered_Scaffold_23", "ChrY_unordered_Scaffold_23", 0, 1000, xlimMax = 2e6, ylimMax = 2e6)
plot_mummer(coordsY, "ChrY_unordered_Scaffold_25", "ChrY_unordered_Scaffold_25", 0, 1000, xlimMax = 2e6, ylimMax = 2e6)
plot_mummer(coordsY, "ChrY_unordered_Scaffold_32", "ChrY_unordered_Scaffold_32", 0, 1000, xlimMax = 2e6, ylimMax = 2e6)

#hist(coords10$QAliLen)
