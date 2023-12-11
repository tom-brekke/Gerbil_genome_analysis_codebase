#!/bin/bash

#this was written by TDB on 05 July 2021
#is should go through a list of defined genomes and run the nessie program to analyze entropy for each
#maybe also linguistic complexity


#genomes here: 
#full path names
declare -a genomes=(
					"/Volumes/MulleySeqs_Working/Assemblies/Meriones_Final/Meriones_chromonome_v6.1.fasta"
					#"~/Dropbox/Post_Doc/Genomes/Psammomys_obesus/Second_sandrat/mPsaObe1.curated_primary.fa" #not to chromosome scale
				#	"~/Dropbox/Post_Doc/Genomes/Pachyuromys_duprasi/Pachyuromys_duprasi_HiC.fasta"
					#"~/Dropbox/Post_Doc/Genomes/Rhombomys_opimus/gerbilAsm.fasta" #not to chromosome scale
					#"~/Dropbox/Post_Doc/Genomes/Mus_musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa" #throwing a segfault and not getting the whole genome. not sure what that's about.
				#	"~/Dropbox/Post_Doc/Genomes/Rattus_norvegicus/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
				#	"/Users/bss81d/Dropbox/Post_Doc/Genomes/Peromyscus_maniculatus/GCA_003704035.1_HU_Pman_2.1_genomic.fasta" #deer mouse?
				#	"/Users/bss81d/Dropbox/Post_Doc/Genomes/Microtus_ochrogaster/GCF_000317375.1_MicOch1.0_genomic.fasta" #vole?
				#	"/Users/bss81d/Dropbox/Post_Doc/Genomes/Cricetulus_griseus/CriGri_chromosomes.fasta" #Chinese hamster?
				#	"/Users/bss81d/Dropbox/Post_Doc/Genomes/Psammomys_obesus/Thybert/mPsaObe1_REL_1811_v2.fa"
					)


outDir="/Volumes/MulleySeqs_Working/Genome_comparison/Nessie/"
cd $outDir
pwd
declare -i length=10000
declare -i step=1000




#here is the loop for the Nessie to run the entropy calculation:
for step in 1000 10000
do
	for genome in "${genomes[@]}"
	do
		genome_file_stem=${genome/sta/}
		genome_file_stem=${genome_file_stem/.fa/}
		genome_name=${genome_file_stem##*/} #drops substring from start of string to last occurence of substring: i.e. the path to the genome.
		#echo $genome_file_stem
		#echo $genome_name
		outfile=$genome_name.nessieOut_E_l"$length"_s"$step".txt
		date 
		run_cmd.sh "~/software/nessie-master/nessie -I $genome -O $outfile -E -l $length -s $step" $outfile
		run_cmd.sh "Nessie_to_df.py -N $outfile -o ${outfile/.txt/_asDF.csv}" ${outfile/.txt/_asDF.csv}
	done
done




#here is the loop for linguistic complexity - much longer! 

# 
for step in 1000 10000
do
	for genome in "${genomes[@]}"
 	do
 		genome_file_stem=${genome/sta/}
 		genome_file_stem=${genome_file_stem/.fa/}
 		genome_name=${genome_file_stem##*/}  #drops substring from start of string to last occurence of substring: i.e. the path to the genome.
		#echo $genome_file_stem
	 	#echo $genome_name
	 	outfile=$genome_name.nessieOut_L_l"$length"_s"$step".txt 
		date 
	 	run_cmd.sh "~/software/nessie-master/nessie -I $genome -O $outfile -L -l $length -s $step" $outfile
 		run_cmd.sh "Nessie_to_df.py -N $outfile -o ${outfile/.txt/_asDF.csv}" ${outfile/.txt/_asDF.csv}
 	done
done
# 
