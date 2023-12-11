#!/bin/bash


##7 is for the Meriones_chromonome_v6.1 genome.
#updated 16-12-2020
#and again on 04/08/2020 to deal with dovetail changing the scaffold names on me :(
#and again on 21/03/2022 for v6.1

#this is the script that runs the GBS pipeline for the Meriones. 
#written by TDB on 6 November 2017. 



	#original raw gzipped data should be in:	 /UDrive/CNS/ResearchData/BIOL-GROUP-Mulley/tom_brekke/raw/Meriones_GBS_LGC_1a/RE_processed/
	#original raw gzipped data should be in:	 /UDrive/CNS/ResearchData/BIOL-GROUP-Mulley/tom_brekke/raw/Meriones_GBS_LGC_1b/RE_processed/
	#original raw gzipped data should be in:	 /UDrive/CNS/ResearchData/BIOL-GROUP-Mulley/tom_brekke/raw/Meriones_GBS_LGC_2/AdapterClipped/  #didn't come RE_processed, not sure why... But it also looks like RE processing doesn't remove much, so it shouldn't be much of a problem. - any wonky reads won't form stacks and will be dropped later anyway.

		#some of the raw data needs to be manipulated - 
	#1b needs to be renamed - they included sample names in the titles. 
		#here is how I renamed them:
		#for files in folders:
			#for i in S*E*/*; do echo $i; j=${i/\/[[:digit:]]*-E//\E}; mv $i $j; done
		#for folders:
			#for i in S*E*; do echo $i; j=${i/_[[:digit:]]*-E/_E};  mv $i $j; done	
	#2 needs to be as well:
		#this renames the files
			#for i in Sample_[[:digit:]]_*/*; do echo $i; j=${i/\/[[:digit:]]_/\/}; echo $j; mv $i $j; done
			#for i in Sample_[[:digit:]][[:digit:]]_*/*; do echo $i; j=${i/\/[[:digit:]][[:digit:]]_/\/}; echo $j; mv $i $j;done
			#for i in Sample_[[:digit:]][[:digit:]]-*/*; do echo $i; j=${i/\/[[:digit:]][[:digit:]]-/\/}; echo $j; mv $i $j;done
			#for i in Sample_[[:digit:]][[:digit:]][[:digit:]]-*/*; do echo $i; j=${i/\/[[:digit:]][[:digit:]][[:digit:]]-/\/}; echo $j; mv $i $j;done
			
		#this renames the dirs
			#for i in Sample_[[:digit:]][[:digit:]][[:digit:]]-*; do echo $i; j=${i/Sample_[[:digit:]][[:digit:]][[:digit:]]-/Sample_};echo $j ;mv $i $j;done
			#for i in Sample_[[:digit:]][[:digit:]][[:digit:]]_*; do echo $i; j=${i/Sample_[[:digit:]][[:digit:]][[:digit:]]/Sample};echo $j ;mv $i $j;done
			#for i in Sample_[[:digit:]][[:digit:]]-*; do echo $i; j=${i/Sample_[[:digit:]][[:digit:]]-/Sample_};echo $j ;mv $i $j;done
			#for i in Sample_[[:digit:]][[:digit:]]_*; do echo $i; j=${i/Sample_[[:digit:]][[:digit:]]/Sample};echo $j ;mv $i $j;done
			#for i in Sample_[[:digit:]]_*; do echo $i; j=${i/Sample_[[:digit:]]/Sample};echo $j ;mv $i $j;done
		#then some files I missed - the '-' ones...:
			#for i in Sample_*/[[:digit:]][[:digit:]][[:digit:]]*; do echo $i; j=${i/\/[[:digit:]][[:digit:]][[:digit:]]_/\/}; echo $j; mv $i $j; done
		


################################################
###  SET THESE PRIOR TO RUNNING THE SCRIPT:  ###
################################################

#This is the location of the working directory:

homedir="/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/genetic_map/stacks"

#These are the location of the raw data - many of these were re-seqeunced. All combined in the SRA upload.:
	raw_data_dir_1a="/Volumes/MulleySeqs_Working/RAW/Meriones_GBS_Lucigen/Meriones_GBS_LGC_1a_30-Oct-2017/AdapterClipped"
	raw_data_dir_1b="/Volumes/MulleySeqs_Working/RAW/Meriones_GBS_Lucigen/Meriones_GBS_LGC_1b_08-Jan-2018/AdapterClipped"
	raw_data_dir_2="/Volumes/MulleySeqs_Working/RAW/Meriones_GBS_Lucigen/Meriones_GBS_LGC_2_03-Apr-2018/AdapterClipped"

#for the genome - used to link radtags to a genome contig.
GENOME_LOCATION="/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta"
run_cmd.sh "bwa index $GENOME_LOCATION" $GENOME_LOCATION.amb

################################################
###                                         ####
################################################



#align all fastqs, then merge bams together as needed. 
if [ ! -d $homedir ]
then
	mkdir $homedir
fi	
cd $homedir
if [ ! -d stacks_out ]
then 	
	mkdir raw
	mkdir raw/1a
	mkdir raw/1b
	mkdir raw/2
	mkdir stacks_out
	mkdir alignments
	mkdir alignments/1a
	mkdir alignments/1b
	mkdir alignments/2
	mkdir merged_alignments
	mkdir alignments/control_files/ #each is named something like "Sample_EEEE0-1F.sorted.RG.control.txt"  and in it is a list of the file names for that individual - all were sqeuenced twice, some three times, this just allows me to merge all the sequencing efforts from each individual for the alignment.
fi
cd alignments

#for each directory in the 3 sequencing runs: ($d)
for d in $raw_data_dir_1a $raw_data_dir_1b $raw_data_dir_2  #complete file paths ending in "1a", "1b", or "2"
do
	#echo "waiting 1"
	#read var
	run=$(echo $d | cut -d"_" -f 7) #set to "1a", "1b", or "2"
	#echo $run
	cd $d # $d is a full pathway
	#for each sample ($s) in the sequencing run directory (the sample is a directory as well):
	for s in Sample*
	do
		#echo "waiting 2"
		#read var
		cd $d/$s #$d is a full pathway so this works
		#for each set of reads in the sample directory
		for R1 in *R1_clipped.fastq.bz2
		do
			#echo "waiting 3"
			#read var
			cd $d/$s
			R2=${R1/R1_clipped.fastq.bz2/R2_clipped.fastq.bz2}
#			echo $R1 
#			echo $s.sam
#			echo $s.sorted.bam
			if [ ! -f $homedir/alignments/$run/$s.sorted.bam ]
			then
				#echo "waiting 4"
				#read var
				if [ ! -f $homedir/alignments/$run/$s.sorted.RG.bam ]
				then
					#echo "waiting 5"
					#read var
					run_cmd.sh "bunzip2 -kc $R1 > $homedir/raw/$run/${R1/.bz2/} " $homedir/raw/$run/${R1/.bz2/} &
					#echo "done"
					#read var
					run_cmd.sh "bunzip2 -kc $R2 > $homedir/raw/$run/${R2/.bz2/} " $homedir/raw/$run/${R2/.bz2/} 
					#echo "done"
					#read var
					wait
					run_cmd.sh "bwa mem -t 2 $GENOME_LOCATION $homedir/raw/$run/${R1/.bz2/} $homedir/raw/$run/${R2/.bz2/} > $homedir/alignments/$run/$s.sam" $homedir/alignments/$run/$s.sam 
					run_cmd.sh "rm $homedir/raw/$run/${R1/.bz2/}" 1 &
					run_cmd.sh "rm $homedir/raw/$run/${R2/.bz2/}" 1 &
					#echo "done"
					#read var 				
					cd $homedir/alignments/$run/
					run_cmd.sh "Convert_sam_to_sorted_bam.sh $s.sam" $homedir/alignments/$run/$s.sorted.bam &
				fi
			fi
			#echo "done"
			#read var
		done	
	done
done	


wait

cd $homedir/alignments/1a
ID="ID:1a"
if [ ! -f ../readGroupsAdded.flag ]
then
	for b in *sorted.bam
	do
		run_cmd.sh "samtools addreplacerg -O BAM -r $ID -o ${b/.bam/.RG.bam} $b " ${b/.bam/.RG.bam} && rm $b
	done	
	cd $homedir/alignments/1b
	ID="ID:1b"
	for b in *sorted.bam
	do
		run_cmd.sh "samtools addreplacerg -O BAM -r $ID -o ${b/.bam/.RG.bam} $b " ${b/.bam/.RG.bam} && rm $b
	done
	cd $homedir/alignments/2
	ID="ID:2"
	for b in *sorted.bam
	do
		run_cmd.sh "samtools addreplacerg -O BAM -r $ID -o ${b/.bam/.RG.bam} $b " ${b/.bam/.RG.bam} && rm $b
	done
	date >> $homedir/alignments/ReadGroupsAdded.flag
fi



#now change the sorted.RG.sorted.bam to sorted.bam for all three dirs and run the following:

#next merge the bam files that need be merged
cd $homedir/alignments/1a 


#this goes through and makes a control file for each sample that has each sequencing run in it, then merge all those files. .
for i in *.sorted.RG.bam #name of the sequencing reads
do
	echo 1a/$i > $homedir/alignments/control_files/${i/.bam/.control.txt}
	#check if the sample has been sequenced a second time in 1b
	if [ -f $homedir/alignments/1b/$i ] 
	then
		echo 1b/$i >> $homedir/alignments/control_files/${i/.bam/.control.txt}
	fi

	#folder 2 can have both a dash or a dot in its sequenceing:
	if [ -f $homedir/alignments/2/$i ] 
	then
		echo 2/$i >> $homedir/alignments/control_files/${i/.bam/.control.txt}
	fi
	doti=${i/-/.}
	if [ -f $homedir/alignments/2/$doti ]
	then
		echo 2/$doti >> $homedir/alignments/control_files/${i/.bam/.control.txt}
	fi
done



cd $homedir/alignments
flag=0
for cf in ./control_files/*
do
	let "flag=flag+1"
	#echo $cf
	outfh1=${cf/.\/control_files\/Sample_/}
	outfh=${outfh1/.control.txt/.merged.bam}
	#echo $outfh 
	run_cmd.sh "samtools merge -b $cf $homedir/merged_alignments/$outfh" $homedir/merged_alignments/$outfh &
	if test $flag -eq 3  #this sends three of the samtools merge calls out, then waits for them to finish. "3" can instead be set to the number of processors available.
	then
		wait
		flag=0
	fi
done
wait

##########################################################################################
#Next step: stacks2!
##########################################################################################
cd $homedir


#makes the population maps and sex.csv file:
mkdir rqtl
run_cmd.sh "printf  \"id,sex,pgm\nEEEE0-1F,f,na\nSSSS0-1M,m,na\nESES1-1F,f,0\nESES1-2M,m,0\nESES1-3M,m,0\nESES1-4M,m,0\nESES1-5M,m,0\nESES1-6M,m,0\nESES10-1M,m,0\nESES10-2F,f,0\nESES10-3M,m,0\nESES10-4F,f,0\nESES10-5F,f,0\nESES11-1F,f,0\nESES11-2M,m,0\nESES11-3,m,0\nESES11-4,f,0\nESES11-5,m,0\nESES11-6,m,0\nESES12-1F,f,0\nESES12-2,m,0\nESES12-3,m,0\nESES12-4,m,0\nESES12-5M,m,0\nESES12-6F,f,0\nESES13-1F,f,0\nESES13-2F,f,0\nESES13-3F,f,0\nESES13-4F,f,0\nESES13-6,m,0\nESES14-10F,f,0\nESES14-1F,f,0\nESES14-2M,m,0\nESES14-3M,m,0\nESES14-4F,f,0\nESES14-5M,m,0\nESES14-6M,m,0\nESES14-7F,f,0\nESES14-8M,m,0\nESES14-9F,f,0\nESES15-1F,f,0\nESES15-2F,f,0\nESES15-3F,f,0\nESES15-4M,m,0\nESES15-5M,m,0\nESES15-6M,m,0\nESES15-7M,m,0\nESES16-1F,f,0\nESES16-2M,m,0\nESES16-3F,f,0\nESES16-4F,f,0\nESES16-5M,m,0\nESES16-6M,m,0\nESES16-7M,m,0\nESES16-8F,f,0\nESES17-1M,m,0\nESES17-2F,f,0\nESES17-3M,m,0\nESES17-4F,f,0\nESES17-5F,f,0\nESES17-6F,f,0\nESES17-7F,f,0\nESES18-1F,f,0\nESES18-2F,f,0\nESES18-3F,f,0\nESES18-4F,m,0\nESES18-5M,m,0\nESES18-6M,f,0\nESES18-7M,m,0\nESES19-1,m,0\nESES19-2,m,0\nESES19-3,m,0\nESES19-4,f,0\nESES19-5,m,0\nESES19-6,f,0\nESES19-7,m,0\nESES19-8,f,0\nESES19-9,m,0\nESES2-1M,m,0\nESES2-2F,f,0\nESES2-3M,m,0\nESES2-4M,m,0\nESES2-5M,m,0\nESES20-1,f,0\nESES20-2,f,0\nESES20-3,m,0\nESES20-4,m,0\nESES20-5,m,0\nESES20-6,f,0\nESES20-7,f,0\nESES21-1,m,0\nESES21-2,m,0\nESES21-3,f,0\nESES21-4,f,0\nESES21-5,f,0\nESES21-6,m,0\nESES22-10,f,0\nESES22-1,m,0\nESES22-2,m,0\nESES22-3,f,0\nESES22-4,m,0\nESES22-5,f,0\nESES22-6,m,0\nESES22-7,m,0\nESES22-8,f,0\nESES22-9,f,0\nESES3-1M,m,0\nESES3-2F,m,0\nESES3-3F,f,0\nESES3-4M,m,0\nESES4-1F,f,0\nESES4-2F,f,0\nESES4-3M,m,0\nESES5-1M,m,0\nESES5-2M,m,0\nESES5-3M,m,0\nESES5-4M,m,0\nESES5-5F,f,0\nESES5b-1M,m,0\nESES5b-2M,m,0\nESES5b-3,m,0\nESES6-1M,m,0\nESES6-2M,m,0\nESES6-3M,m,0\nESES6-4F,f,0\nESES6-5F,f,0\nESES6-6M,m,0\nESES6-7F,f,0\nESES7-1F,f,0\nESES7-2F,f,0\nESES7-3F,f,0\nESES7-4F,f,0\nESES8-1F,f,0\nESES8-2F,f,0\nESES9-1F,f,0\nESES9-2M,m,0\nESES9-3F,f,0\nESES9-4M,f,0\nESES9-5M,m,0\" > ./rqtl/sex.csv " ./rqtl/sex.csv
run_cmd.sh "printf \"EEEE0-1F.sorted.RG.merged\tparent \nEESS2-1F.sorted.RG.merged\tprogeny\nEESS2-3F.sorted.RG.merged\tprogeny\nEESS2-4M.sorted.RG.merged\tprogeny\nEESS2-6M.sorted.RG.merged\tprogeny\nEESS3-1F.sorted.RG.merged\tprogeny\nEESS3-2F.sorted.RG.merged\tprogeny\nEESS3-3M.sorted.RG.merged\tprogeny\nEESS3-4M.sorted.RG.merged\tprogeny\nESES1-1F.sorted.RG.merged\tprogeny\nESES1-2M.sorted.RG.merged\tprogeny\nESES1-3M.sorted.RG.merged\tprogeny\nESES1-4M.sorted.RG.merged\tprogeny\nESES1-5M.sorted.RG.merged\tprogeny\nESES1-6M.sorted.RG.merged\tprogeny\nESES10-1M.sorted.RG.merged\tprogeny\nESES10-2F.sorted.RG.merged\tprogeny\nESES10-3M.sorted.RG.merged\tprogeny\nESES10-4F.sorted.RG.merged\tprogeny\nESES10-5F.sorted.RG.merged\tprogeny\nESES11-1F.sorted.RG.merged\tprogeny\nESES11-2M.sorted.RG.merged\tprogeny\nESES11-3.sorted.RG.merged\tprogeny\nESES11-4.sorted.RG.merged\tprogeny\nESES11-5.sorted.RG.merged\tprogeny\nESES11-6.sorted.RG.merged\tprogeny\nESES12-1F.sorted.RG.merged\tprogeny\nESES12-2.sorted.RG.merged\tprogeny\nESES12-3.sorted.RG.merged\tprogeny\nESES12-4.sorted.RG.merged\tprogeny\nESES12-5M.sorted.RG.merged\tprogeny\nESES12-6F.sorted.RG.merged\tprogeny\nESES13-1F.sorted.RG.merged\tprogeny\nESES13-2F.sorted.RG.merged\tprogeny\nESES13-3F.sorted.RG.merged\tprogeny\nESES13-4F.sorted.RG.merged\tprogeny\nESES13-6.sorted.RG.merged\tprogeny\nESES14-10F.sorted.RG.merged\tprogeny\nESES14-1F.sorted.RG.merged\tprogeny\nESES14-2M.sorted.RG.merged\tprogeny\nESES14-3M.sorted.RG.merged\tprogeny\nESES14-4F.sorted.RG.merged\tprogeny\nESES14-5M.sorted.RG.merged\tprogeny\nESES14-6M.sorted.RG.merged\tprogeny\nESES14-7F.sorted.RG.merged\tprogeny\nESES14-8M.sorted.RG.merged\tprogeny\nESES14-9F.sorted.RG.merged\tprogeny\nESES15-1F.sorted.RG.merged\tprogeny\nESES15-2F.sorted.RG.merged\tprogeny\nESES15-3F.sorted.RG.merged\tprogeny\nESES15-4M.sorted.RG.merged\tprogeny\nESES15-5M.sorted.RG.merged\tprogeny\nESES15-6M.sorted.RG.merged\tprogeny\nESES15-7M.sorted.RG.merged\tprogeny\nESES16-1F.sorted.RG.merged\tprogeny\nESES16-2M.sorted.RG.merged\tprogeny\nESES16-3F.sorted.RG.merged\tprogeny\nESES16-4F.sorted.RG.merged\tprogeny\nESES16-5M.sorted.RG.merged\tprogeny\nESES16-6M.sorted.RG.merged\tprogeny\nESES16-7M.sorted.RG.merged\tprogeny\nESES16-8F.sorted.RG.merged\tprogeny\nESES17-1M.sorted.RG.merged\tprogeny\nESES17-2F.sorted.RG.merged\tprogeny\nESES17-3M.sorted.RG.merged\tprogeny\nESES17-4F.sorted.RG.merged\tprogeny\nESES17-5F.sorted.RG.merged\tprogeny\nESES17-6F.sorted.RG.merged\tprogeny\nESES17-7F.sorted.RG.merged\tprogeny\nESES18-1F.sorted.RG.merged\tprogeny\nESES18-2F.sorted.RG.merged\tprogeny\nESES18-3F.sorted.RG.merged\tprogeny\nESES18-4F.sorted.RG.merged\tprogeny\nESES18-5M.sorted.RG.merged\tprogeny\nESES18-6M.sorted.RG.merged\tprogeny\nESES18-7M.sorted.RG.merged\tprogeny\nESES19-1.sorted.RG.merged\tprogeny\nESES19-2.sorted.RG.merged\tprogeny\nESES19-3.sorted.RG.merged\tprogeny\nESES19-4.sorted.RG.merged\tprogeny\nESES19-5.sorted.RG.merged\tprogeny\nESES19-6.sorted.RG.merged\tprogeny\nESES19-7.sorted.RG.merged\tprogeny\nESES19-8.sorted.RG.merged\tprogeny\nESES19-9.sorted.RG.merged\tprogeny\nESES2-1M.sorted.RG.merged\tprogeny\nESES2-2F.sorted.RG.merged\tprogeny\nESES2-3M.sorted.RG.merged\tprogeny\nESES2-4M.sorted.RG.merged\tprogeny\nESES2-5M.sorted.RG.merged\tprogeny\nESES20-1.sorted.RG.merged\tprogeny\nESES20-2.sorted.RG.merged\tprogeny\nESES20-3.sorted.RG.merged\tprogeny\nESES20-4.sorted.RG.merged\tprogeny\nESES20-5.sorted.RG.merged\tprogeny\nESES20-6.sorted.RG.merged\tprogeny\nESES20-7.sorted.RG.merged\tprogeny\nESES21-1.sorted.RG.merged\tprogeny\nESES21-2.sorted.RG.merged\tprogeny\nESES21-3.sorted.RG.merged\tprogeny\nESES21-4.sorted.RG.merged\tprogeny\nESES21-5.sorted.RG.merged\tprogeny\nESES21-6.sorted.RG.merged\tprogeny\nESES22-10.sorted.RG.merged\tprogeny\nESES22-1.sorted.RG.merged\tprogeny\nESES22-2.sorted.RG.merged\tprogeny\nESES22-3.sorted.RG.merged\tprogeny\nESES22-4.sorted.RG.merged\tprogeny\nESES22-5.sorted.RG.merged\tprogeny\nESES22-6.sorted.RG.merged\tprogeny\nESES22-7.sorted.RG.merged\tprogeny\nESES22-8.sorted.RG.merged\tprogeny\nESES22-9.sorted.RG.merged\tprogeny\nESES3-1M.sorted.RG.merged\tprogeny\nESES3-2F.sorted.RG.merged\tprogeny\nESES3-3F.sorted.RG.merged\tprogeny\nESES3-4M.sorted.RG.merged\tprogeny\nESES4-1F.sorted.RG.merged\tprogeny\nESES4-2F.sorted.RG.merged\tprogeny\nESES4-3M.sorted.RG.merged\tprogeny\nESES5-1M.sorted.RG.merged\tprogeny\nESES5-2M.sorted.RG.merged\tprogeny\nESES5-3M.sorted.RG.merged\tprogeny\nESES5-4M.sorted.RG.merged\tprogeny\nESES5-5F.sorted.RG.merged\tprogeny\nESES5b-1M.sorted.RG.merged\tprogeny\nESES5b-2M.sorted.RG.merged\tprogeny\nESES5b-3.sorted.RG.merged\tprogeny\nESES6-1M.sorted.RG.merged\tprogeny\nESES6-2M.sorted.RG.merged\tprogeny\nESES6-3M.sorted.RG.merged\tprogeny\nESES6-4F.sorted.RG.merged\tprogeny\nESES6-5F.sorted.RG.merged\tprogeny\nESES6-6M.sorted.RG.merged\tprogeny\nESES6-7F.sorted.RG.merged\tprogeny\nESES7-1F.sorted.RG.merged\tprogeny\nESES7-2F.sorted.RG.merged\tprogeny\nESES7-3F.sorted.RG.merged\tprogeny\nESES7-4F.sorted.RG.merged\tprogeny\nESES8-1F.sorted.RG.merged\tprogeny\nESES8-2F.sorted.RG.merged\tprogeny\nESES9-1F.sorted.RG.merged\tprogeny\nESES9-2M.sorted.RG.merged\tprogeny\nESES9-3F.sorted.RG.merged\tprogeny\nESES9-4M.sorted.RG.merged\tprogeny\nESES9-5M.sorted.RG.merged\tprogeny\nSSSS0-1M.sorted.RG.merged\tparent\" > population_map.txt" population_map.txt
run_cmd.sh "printf \"EEEE0-1F.sorted.RG.merged\t1\nESES1-1F.sorted.RG.merged\t1\nESES1-2M.sorted.RG.merged\t1\nESES1-3M.sorted.RG.merged\t1\nESES1-4M.sorted.RG.merged\t1\nESES1-5M.sorted.RG.merged\t1\nESES1-6M.sorted.RG.merged\t1\nESES10-1M.sorted.RG.merged\t1\nESES10-2F.sorted.RG.merged\t1\nESES10-3M.sorted.RG.merged\t1\nESES10-4F.sorted.RG.merged\t1\nESES10-5F.sorted.RG.merged\t1\nESES11-1F.sorted.RG.merged\t1\nESES11-2M.sorted.RG.merged\t1\nESES11-3.sorted.RG.merged\t1\nESES11-4.sorted.RG.merged\t1\nESES11-5.sorted.RG.merged\t1\nESES11-6.sorted.RG.merged\t1\nESES12-1F.sorted.RG.merged\t1\nESES12-2.sorted.RG.merged\t1\nESES12-3.sorted.RG.merged\t1\nESES12-4.sorted.RG.merged\t1\nESES12-5M.sorted.RG.merged\t1\nESES12-6F.sorted.RG.merged\t1\nESES13-1F.sorted.RG.merged\t1\nESES13-2F.sorted.RG.merged\t1\nESES13-3F.sorted.RG.merged\t1\nESES13-4F.sorted.RG.merged\t1\nESES13-6.sorted.RG.merged\t1\nESES14-10F.sorted.RG.merged\t1\nESES14-1F.sorted.RG.merged\t1\nESES14-2M.sorted.RG.merged\t1\nESES14-3M.sorted.RG.merged\t1\nESES14-4F.sorted.RG.merged\t1\nESES14-5M.sorted.RG.merged\t1\nESES14-6M.sorted.RG.merged\t1\nESES14-7F.sorted.RG.merged\t1\nESES14-8M.sorted.RG.merged\t1\nESES14-9F.sorted.RG.merged\t1\nESES15-1F.sorted.RG.merged\t1\nESES15-2F.sorted.RG.merged\t1\nESES15-3F.sorted.RG.merged\t1\nESES15-4M.sorted.RG.merged\t1\nESES15-5M.sorted.RG.merged\t1\nESES15-6M.sorted.RG.merged\t1\nESES15-7M.sorted.RG.merged\t1\nESES16-1F.sorted.RG.merged\t1\nESES16-2M.sorted.RG.merged\t1\nESES16-3F.sorted.RG.merged\t1\nESES16-4F.sorted.RG.merged\t1\nESES16-5M.sorted.RG.merged\t1\nESES16-6M.sorted.RG.merged\t1\nESES16-7M.sorted.RG.merged\t1\nESES16-8F.sorted.RG.merged\t1\nESES17-1M.sorted.RG.merged\t1\nESES17-2F.sorted.RG.merged\t1\nESES17-3M.sorted.RG.merged\t1\nESES17-4F.sorted.RG.merged\t1\nESES17-5F.sorted.RG.merged\t1\nESES17-6F.sorted.RG.merged\t1\nESES17-7F.sorted.RG.merged\t1\nESES18-1F.sorted.RG.merged\t1\nESES18-2F.sorted.RG.merged\t1\nESES18-3F.sorted.RG.merged\t1\nESES18-4F.sorted.RG.merged\t1\nESES18-5M.sorted.RG.merged\t1\nESES18-6M.sorted.RG.merged\t1\nESES18-7M.sorted.RG.merged\t1\nESES19-1.sorted.RG.merged\t1\nESES19-2.sorted.RG.merged\t1\nESES19-3.sorted.RG.merged\t1\nESES19-4.sorted.RG.merged\t1\nESES19-5.sorted.RG.merged\t1\nESES19-6.sorted.RG.merged\t1\nESES19-7.sorted.RG.merged\t1\nESES19-8.sorted.RG.merged\t1\nESES19-9.sorted.RG.merged\t1\nESES2-1M.sorted.RG.merged\t1\nESES2-2F.sorted.RG.merged\t1\nESES2-3M.sorted.RG.merged\t1\nESES2-4M.sorted.RG.merged\t1\nESES2-5M.sorted.RG.merged\t1\nESES20-1.sorted.RG.merged\t1\nESES20-2.sorted.RG.merged\t1\nESES20-3.sorted.RG.merged\t1\nESES20-4.sorted.RG.merged\t1\nESES20-5.sorted.RG.merged\t1\nESES20-6.sorted.RG.merged\t1\nESES20-7.sorted.RG.merged\t1\nESES21-1.sorted.RG.merged\t1\nESES21-2.sorted.RG.merged\t1\nESES21-3.sorted.RG.merged\t1\nESES21-4.sorted.RG.merged\t1\nESES21-5.sorted.RG.merged\t1\nESES21-6.sorted.RG.merged\t1\nESES22-10.sorted.RG.merged\t1\nESES22-1.sorted.RG.merged\t1\nESES22-2.sorted.RG.merged\t1\nESES22-3.sorted.RG.merged\t1\nESES22-4.sorted.RG.merged\t1\nESES22-5.sorted.RG.merged\t1\nESES22-6.sorted.RG.merged\t1\nESES22-7.sorted.RG.merged\t1\nESES22-8.sorted.RG.merged\t1\nESES22-9.sorted.RG.merged\t1\nESES3-1M.sorted.RG.merged\t1\nESES3-2F.sorted.RG.merged\t1\nESES3-3F.sorted.RG.merged\t1\nESES3-4M.sorted.RG.merged\t1\nESES4-1F.sorted.RG.merged\t1\nESES4-2F.sorted.RG.merged\t1\nESES4-3M.sorted.RG.merged\t1\nESES5-1M.sorted.RG.merged\t1\nESES5-2M.sorted.RG.merged\t1\nESES5-3M.sorted.RG.merged\t1\nESES5-4M.sorted.RG.merged\t1\nESES5-5F.sorted.RG.merged\t1\nESES5b-1M.sorted.RG.merged\t1\nESES5b-2M.sorted.RG.merged\t1\nESES5b-3.sorted.RG.merged\t1\nESES6-1M.sorted.RG.merged\t1\nESES6-2M.sorted.RG.merged\t1\nESES6-3M.sorted.RG.merged\t1\nESES6-4F.sorted.RG.merged\t1\nESES6-5F.sorted.RG.merged\t1\nESES6-6M.sorted.RG.merged\t1\nESES6-7F.sorted.RG.merged\t1\nESES7-1F.sorted.RG.merged\t1\nESES7-2F.sorted.RG.merged\t1\nESES7-3F.sorted.RG.merged\t1\nESES7-4F.sorted.RG.merged\t1\nESES8-1F.sorted.RG.merged\t1\nESES8-2F.sorted.RG.merged\t1\nESES9-1F.sorted.RG.merged\t1\nESES9-2M.sorted.RG.merged\t1\nESES9-3F.sorted.RG.merged\t1\nESES9-4M.sorted.RG.merged\t1\nESES9-5M.sorted.RG.merged\t1\nSSSS0-1M.sorted.RG.merged\t1\n\" > population_map_onePop.txt " population_map_onePop.txt


run_cmd.sh "ref_map.pl --samples $homedir/merged_alignments/ --popmap $homedir/population_map.txt -o $homedir/stacks_out -T 3" $homedir/stacks_out/gstacks.log #-s sorted.RG.merged doesn't seem to work. modified the population map instead. #DONT USE "remove PCR duplicates""! its GBS - all the instert sizes are the same - between two cutsites...
run_cmd.sh "/usr/local/bin/populations -P $homedir/stacks_out -O $homedir/rqtl/ -M $homedir/population_map_OnePop.txt --threads 3 --ordered-export --vcf -R 0.90 --merge-sites -e MslI --min-mac 30 --min-maf 0.1" $homedir/rqtl/populations.log
run_cmd.sh "make_rqtl_input_file.py -v $homedir/rqtl/populations.snps.vcf -o $homedir/rqtl/Meriones_chromonome_v6.1.geneticMap -t .sorted.RG.merged --parents EEEE0-1F,SSSS0-1M -p $homedir/rqtl/sex.csv 2>> $homedir/rqtl/make_rqtl_metadata.txt" 1



date >> "$homedir/run_stacks2_complete.txt"





