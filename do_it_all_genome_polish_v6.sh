#!/bin/bash

#this takes verson 5.2 and makes 6 then 6.1 (6.1 is 6.0 with chr5 flipped the right way round)
#then runs all subsequent analyses



#once Meriones_chromonome_v6.0.assembly.csv, Meriones_chromonome_v6.0.breaks.csv, and Meriones_chromonome_v6.0.scaffolds_to_drop.csv have been created (by hand), see the README for a list of the changes. 
#this bit assembles version 6.0:
#also add in contig_107 and contig_2709 from assembly.fasta, the ONT Flye assembly as these have the parahox region on them. Rename to Chr13_un



########################################################################################
#build genome
########################################################################################
run_cmd.sh "assemble_genome_and_recoordinate_gff.py -a /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.0.assembly.csv -f /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/old_versions/V5.2/Meriones_chromonome_v5.2.fasta  -g /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/old_versions/V5.2/Meriones_chromonome_v5.2.sorted.gff -b /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.0.breaks.csv -d /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.0.scaffolds_to_drop.csv -o /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final//old_versions/V6.0/Meriones_chromonome_v6.0 && cat /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/old_versions/V5.2/Mitochondiral_Genes/MT_CONSENSUS.fa >> /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final//old_versions/V6.0/Meriones_chromonome_v6.0.fasta && cat /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/old_versions/V5.2/Mitochondiral_Genes/MT_CONSENSUS.gff >> /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final//old_versions/V6.0/Meriones_chromonome_v6.0.gff && cat /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/old_versions/V5.2/ONT_denovo/All/assembly.contig_107_and_contig_2907.fasta >> /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final//old_versions/V6.0/Meriones_chromonome_v6.0.fasta" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final//old_versions/V6.0/Meriones_chromonome_v6.0.fasta
run_cmd.sh "assemble_genome_and_recoordinate_gff.py -a /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.assembly.csv -f /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/old_versions/V6.0/Meriones_chromonome_v6.0.fasta  -g /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/old_versions/V6.0/Meriones_chromonome_v6.0.gff   -o /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1 " /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta


run_cmd.sh "makeblastdb -in /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta -dbtype nucl" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta.nhr
run_cmd.sh "bwa index /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta"  /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta.bwt

mkdir /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/v6.1_seqs/
run_cmd.sh "grep \">\" Meriones_chromonome_v6.1.fasta | cut -f 1 -d \" \" | cut -f 2 -d \">\" > /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/v6.1_seqs/Meriones_chromonome_v6.1.chrNames.txt" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/v6.1_seqs/Meriones_chromonome_v6.1.chrNames.txt

#there are a couple things I need each record individually for:
for line in $(cat /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/v6.1_seqs/Meriones_chromonome_v6.1.chrNames.txt)
do
run_cmd.sh "extract_fasta_region.py -f /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta -chr \"$line\" -o \"/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/v6.1_seqs/$line.v6.1.fasta\""   "/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/v6.1_seqs/$line.v6.1.fasta"
done	



########################################################################################
#genetic map
########################################################################################
mkdir /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/
mkdir /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/stacks/
run_cmd.sh "Meriones_Stacks_7.sh" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/stacks/run_stacks2_complete.txt
run_cmd.sh "v6.1.map.rqtl.R" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/stacks/population_map.txt #just put a file in here to make sure it didn't try to run itself. 
echo "the previous command needs to be run by hand"
read var
cp /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/stacks/rqtl/Meriones_chromonome_v6.1.geneticmap.csv /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/




########################################################################################
#GC data etc
########################################################################################
#now start collecting genome-wide stats on variouse things. 
#run_cmd.sh "CalcGC_v1.py -f /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta -o /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.GC.txt" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.GC.txt
#run_cmd.sh "CalcGC_v1.py -f /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta -o /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.1kbwindow.GC.txt -w 1000" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.1kbwindow.GC.txt
#run_cmd.sh "CalcGC_v1.py -f /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta -o /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.10kbwindow.GC.txt -w 10000" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.10kbwindow.GC.txt
run_cmd.sh "extract_chr_lengths_simple.py /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta > /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.chrLengths.csv" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.chrLengths.csv
run_cmd.sh "Calc_R_GC_Gene_density.py -f /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta -g /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.gff -m /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.geneticmap.csv -o /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.R_GC_GeneDens.csv" /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.R_GC_GeneDens.csv



########################################################################################
#Tandem Repeat Finder
########################################################################################
mkdir /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/TandemRepeatFinder
cd /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/TandemRepeatFinder
ln -s ../v6.1_seqs/*fasta ./
counter="0"
mybatch="3"
for Chr in *fasta
do
	counter=$(expr $counter + 1)
	if [ $counter -gt $mybatch ]
	then
		wait
		counter="0"
	fi
	run_cmd.sh "/Users/bss81d/software/trf409.macosx \"$Chr\" 2 7 7 80 10 50 2000 -d -h -l 5" "$Chr.2.7.7.80.10.50.2000.dat" &
done
wait
cat *.2.7.7.80.10.50.2000.dat > Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000.dat
run_cmd.sh "TRFoutput_reformatter.py --dat Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000.dat --centromeres ../cent_locations_rough_v6.1.bed -o Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000 -s" Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000.simplified.csv 
run_cmd.sh "TRFoutput_reformatter.py --dat Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000.dat --centromeres ../cent_locations_rough_v6.1.bed -o Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000 " Meriones_chromonome_v6.1.fasta.2.7.7.80.10.50.2000.csv 


########################################################################################
#Nessie
########################################################################################
cd /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/
run_cmd.sh "Nessie_runner.sh" 1 





########################################################################################
#RepeatMasker
########################################################################################
#need to run repeatmasker - on SuperComputingWales, but this is the script:
if false
then
	run_cmd.sh "Meriones_chromonome_RepeatMasker.sh" 1



	cd /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/RepeatMasker/

	#once that finishes, scp all files back to local computer - repeat masker has been run on each fasta entry in isolation, so they need to be concatenated together and some summary data needs to be extracted.
	#cat the repeat-masked genome: #two entries had no repetitive bits found: ChrUnknown_unordered_Scaffold_150 and  Chr13_unordered_Scaffold_195 so those original seqs need to be added in as well.
	cat *.fasta.masked ../v6.1_seqs/ChrUnknown_unordered_Scaffold_150.v6.1.fasta ../v6.1_seqs/Chr13_unordered_Scaffold_195.v6.1.fasta | seqkit sort -N -n > ../Meriones_chromonome_v6.1.masked.fasta
	#then combine the summary tables/gffs etc 
	grep -v "#" *gff > ../Meriones_chromonome_v6.1.repeats.gff
	run_cmd.sh "RepeatMasker_tbl_summarizer.py -d ./ -o ../Meriones_chromonome_v6.1.masked.summary.tbl" ../Meriones_chromonome_v6.1.masked.summary.tbl

fi




########################################################################################
#mummer
########################################################################################
#on SCW, but here is the commands.

#run_cmd.sh "nucmer --mum -t 40 -p MER6.1vSELF /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta" MER6.1vSELF.delta
#needs to me maxmatch, not mum for the self-alignment to find repeats. errors out with maxmatch, try --batch 1 instead as per here: https://github.com/mummer4/mummer/releases
#needs to have --nosimplify in there to keep the repeats. 
#run_cmd.sh "nucmer --maxmatch --nosimplify -t 60 -p MER6.1vSELF.maxmatch.noSimplify /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta" MER6.1vSELF.maxmatch.noSimplify.delta


#--maxmatch --nosimplify should do, btut doesn't seem to work on the whole genome - maybe chr5 is messing it up?? 
#this ran for the chr13 scaffolds aligned to the chr13, might work for all chr13 scafs aligned to themselves:
run_cmd.sh "nucmer --maxmatch --nosimplify -t 2 -p Chr13ScafsvChr13.maxmatch.nosimplify ../Chr13_repeat_structure/Meriones_chromonome_v6.1.AllChr13Scaffolds.fasta ../v6.1_seqs/Chr13.v6.1.fasta" Chr13ScafsvChr13.maxmatch.nosimplify.delta
run_cmd.sh "nucmer --maxmatch --nosimplify -t 2 -p Chr13ScafsvChr13_unordered_Scaffold_33.maxmatch.nosimplify ../Chr13_repeat_structure/Meriones_chromonome_v6.1.AllChr13Scaffolds.fasta ../v6.1_seqs/Chr13_unordered_Scaffold_33.v6.1.fasta" Chr13ScafsvChr13_unordered_Scaffold_33.maxmatch.nosimplify.delta



for S in ../v6.1_seqs/Chr13*fasta
do
	Q=${S/.v6.1.fasta/}; 
	Q=${Q/..\/v6.1_seqs\//}; 
	run_cmd.sh "nucmer --maxmatch --nosimplify -t 2 -p Chr13Scafsv$Q.maxmatch.nosimplify ../Chr13_repeat_structure/Meriones_chromonome_v6.1.AllChr13Scaffolds.fasta $S"  Chr13Scafsv$Q.maxmatch.nosimplify.delta;
	run_cmd.sh "delta-filter -l 1000 Chr13Scafsv$Q.maxmatch.nosimplify.delta > Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter" Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter; 
	run_cmd.sh "show-coords -Tl Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter > Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter.Tl.coords" Chr13Scafsv$Q.maxmatch.nosimplify.delta.l1000.filter.Tl.coords;
done
run_cmd.sh "cat Chr13Scafsv*.maxmatch.nosimplify.delta.l1000.filter.Tl.coords > Chr13ScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT" Chr13ScafsvSELF.maxmatch.nosimplify.delta.l1000.filter.Tl.coords.CAT


run_cmd.sh "nucmer --mum -t 40 -p MER6.1vPachy /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta /Users/bss81d/Dropbox/Post_Doc/Genomes/Pachyuromys_duprasi/Pachyuromys_duprasi_HiC.fasta" MER6.1vPachy.delta
run_cmd.sh "nucmer --mum -t 40 -p MER6.1vPsam /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta /Users/bss81d/Dropbox/Post_Doc/Genomes/Psammomys_obesus/Thybert/mPsaObe1_REL_1811_v2.fa" MER6.1vPsam.delta
run_cmd.sh "nucmer --mum -t 40 -p MER6.1vMusc /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta /Users/bss81d/Dropbox/Post_Doc/Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa " MER6.1vMusc.delta
run_cmd.sh "nucmer --maxmatch --nosimplify -t 60 -p MER6.1vPsam.maxmatch.noSimplify /Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta /Users/bss81d/Dropbox/Post_Doc/Genomes/Psammomys_obesus/Thybert/mPsaObe1_REL_1811_v2.fa" MER6.1vPsam.maxmatch.noSimplify.delta


#run_cmd.sh "delta-filter -1 -l 500 MER6.1vSELF.delta > MER6.1vSELF.1l500.filter" MER6.1vSELF.1l500.filter
run_cmd.sh "delta-filter    -l 1000 MER6.1vSELF.maxmatch.delta > MER6.1vSELF.maxmatch.l1000.filter" MER6.1vSELF.maxmatch.l1000.filter
run_cmd.sh "delta-filter -1 -l 10000 MER6.1vPachy.delta > MER6.1vPachy.1l10000.filter" MER6.1vPachy.1l10000.filter
run_cmd.sh "delta-filter -1 -l 10000 MER6.1vPsam.delta > MER6.1vPsam.1l10000.filter" MER6.1vPsam.1l10000.filter
run_cmd.sh "delta-filter -1 -l 500 MER6.1vMusc.delta > MER6.1vMusc.1l500.filter" MER6.1vMusc.1l500.filter

run_cmd.sh "show-coords -T MER6.1vSELF.1l10000.filter > MER6.1vSELF.1l10000.filter.T.coords" MER6.1vSELF.1l10000.filter.T.coords
run_cmd.sh "show-coords -T MER6.1vPachy.1l10000.filter> MER6.1vPachy.1l10000.filter.T.coords" MER6.1vPachy.1l10000.filter.T.coords
run_cmd.sh "show-coords -T MER6.1vPsam.1l10000.filter > MER6.1vPsam.1l10000.filter.T.coords" MER6.1vPsam.1l10000.filter.T.coords
run_cmd.sh "show-coords -T MER6.1vMusc.1l500.filter > MER6.1vMusc.1l500.filter.T.coords" MER6.1vMusc.1l500.filter.T.coords



#can do mummer plots:
	#mummerplot MER6.1vPsam.1l10000.filter
#or circle comparisons with my processing script Genome_alignments_plot_v61.pde:
	#copy to processing directory, then run format_data.py to sort it out, i.e.:
	#./format_data.py -c MER6.1vPsam.1l10000.filter.T.coords -o MER6.1vPsam.1l10000.filter.T.coords.formatted.csv
	#./format_data.py -c MER6.1vPachy.1l10000.filter.T.coords -o MER6.1vPachy.1l10000.filter.T.coords.formatted.csv

	

fi




########################################################################################
#NTRprism analysis for centromere repeats
########################################################################################
mkdir NTRprism
cd ./NTRprism

	run_cmd.sh "extract_fasta_region.py -f ../Meriones_chromonome_v6.1.fasta -chrList ../cent_locations_rough_v6.1.bed -o Meriones_chromonome_v6.1.centromeres.fasta" ./Meriones_chromonome_v6.1.centromeres.fasta
	run_cmd.sh "perl ~/software/NTRprism_v0.1/NTRprism_ProcessFasta_v0.1.pl Meriones_chromonome_v6.1.centromeres.fasta Meriones_chromonome_v6.1.centromeres 1 5000 10 5" ./Meriones_chromonome_v6.1.centromeres.region_ChrY_unordered_Scaffold_35.1500000.2800000.span5000.k5.mincount10.bin1.txt
	run_cmd.sh "perl ~/software/NTRprism_v0.1/NTRprism_ProcessFasta_v0.1.pl Meriones_chromonome_v6.1.centromeres.fasta Meriones_chromonome_v6.1.centromeres 100 5000 10 5" ./Meriones_chromonome_v6.1.centromeres.region_ChrY_unordered_Scaffold_35.1500000.2800000.span5000.k5.mincount10.bin100.txt
	#run by hand to make plots for each chr:
	#~/scripts/NTRprism_PlotSpectrum


cd ../






########################################################################################
#Plots for manuscript
########################################################################################
#Number of contigs per chromosome
mkdir summary_figures
cd ./summary_figures
 #should probably get the codes in here - there are 4 r scripts that I used for most of the plotsls
 
 
 #should be run by hand:
if false
then
	Figure_5_GeneContent_by_ChrLength.R  #makes many of the plots for Figure 5
	get_outlier_genes_from_Roddys_data.R #does the GC outlier gene analysis. Also calculates recombination rate in a slididng window
	#Genome_summary_figures_v6.1.R  #many interesting plots as I built and refined the genome, but none in the final manuscript
	mummer_plotter_for_Chr13_alignments.R #plots out the mummer self-alignments
fi



cd ../













