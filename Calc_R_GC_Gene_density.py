#!/usr/bin/env Python3

##########################################################################################
#docstring
##########################################################################################
desc='''
Written by TDB on 12 February 2021
This is a more general script that get_R_GC_Gene_density_of_SFs.py. It should work for any genome with a GFF and map.
This script goes through a genome, gff file, and genetic map for a species and
calculates the gene density, recombination rate, and gc content for each window specified by -w. 

The map needs the following columns in this order. A header line is necessary, but the what each column is named doesn't matter:
SnpID, chr, base_pos, female_cM, male_cM, ave_cM
Only 'female_cM' will be used to calculate the recombintaion rate (so as to include X chrs)


This script pulls many bits from my CalcGC_v1.py script
'''


##########################################################################################
#imports
##########################################################################################
import argparse
import sys
from Bio import SeqIO
import re
from math import ceil
import pprint

##########################################################################################
#globals
##########################################################################################
gff_dict_gene = {} #this will become a nested dict with ECRid as the outer key and start:end as the inner dict as such: {ECRid: {start: end}}
#gff_dict_cds = {} #this will become a nested dict with ECRid as the outer key and start:end as the inner dict as such: {ECRid: {start: end}}
map_dict      = {} #this will become a nested dict as well, with chr as the outer key and start as the next inner key. 'End' and 'R' will be keys in the innermost dict as such: {chr: {start: {"end": end, "R":r}}}. THIS IS IN GENOME COORDINATES, NOT ECR COORDINATES!!!
chrlen_dict   = {} #has the lengths of each chromosome in it, needed to access the final base in each chr. 
##########################################################################################
#functions
##########################################################################################
def process_command_line():
	'''
	This just uses the argparse module to get the command line arguments -
	looking for a fasta, a GFF, a Map, a window size, and possibly an outfile
	'''
	parser = argparse.ArgumentParser(description=desc)                                            
	parser.add_argument("--fasta",                    "-f", type=str, required=True,  help="fasta file. A genome to chromosome scale. NOT repeat-masked. Will need to have an associated *.chrLengths.csv file from extract_chr_lengths_simple.py *.fasta > *.chrLengths.csv", metavar="*.fasta")
	parser.add_argument("--GFF",                      "-g", type=str, required=True,  help="Gff3 file with annotations of cds to calculate gene density")
	parser.add_argument("--Map",                      "-m", type=str, required=False,  help="Genetic map file to calculate recombination rate: SnpID, chr, base_pos, (female_)cM")
	parser.add_argument("--geneDensityWindowSize",    "-d", type=int, required=False,  help="This is the size of the window that gene density will be calculated for. Default = %(default)s", default = 1000000)
	parser.add_argument("--Window",                   "-w", type=int, required=False,  help="Used to run a sliding window of the specified size across the sequences. For recombination and GC content this is the window size, for gene density the window size is set with -d. This is also the step size for R, GC, and GeneDens. Default = %(default)s", default = 1000)
	parser.add_argument("--Out",                      "-o", type=str, required=False, help="Output file name. Omit for stdout")
	args = parser.parse_args()
	return(args)	

def read_in_gff(GFF_FH):
	'''
	This takes the GFF file and creates a dictionary of dictionaries called gff_dict 
	that has ECRid as the outer key and then start: end as the inner dictionary. 
	i.e.: {ECRid: {start: end}} 
	input:   the GFF file handle
	output:  updates the global gff_dict dictionary
	returns: 0 if everything worked, 1 if there was an error
	'''
	print("Parsing GFF file:							", GFF_FH, file = sys.stderr)
	try:
		with open(GFF_FH, 'r') as GFF:
			for line in GFF:
				if not line.startswith("#"):
					if re.search("\tgene\t", line):
						#print(line)
						fields =line.split("	")
						ECRid  = fields[0]
						start  = int(fields[3])
						end    = int(fields[4])
						gff_dict_gene[ECRid] = gff_dict_gene.setdefault(ECRid, {start:end}) #makes the entry if it doesn't exist with an empty dict as the value, if it does exist, leaves it alone
						gff_dict_gene[ECRid].update({start:end})
#					elif re.search("	CDS	", line):
#						fields =line.split("	")
#						ECRid  = fields[0]
#						start  = int(fields[3])
#						end    = int(fields[4])		
#						gff_dict_cds[ECRid] = gff_dict_cds.setdefault(ECRid, {start:end}) #makes the entry if it doesn't exist with an empty dict as the value, if it does exist, leaves it alone
#						gff_dict_cds[ECRid].update({start:end})
	#	'''#'''
		return(0)		
	except:
		return(1)
		
def read_genome_calcAll(FASTA_FH, windowSize, OUT, geneDensityWindowSize):
		'''
	This reads in the genome, breaks it into bits of size 'windowSize' and calculates the 
	GC content of the window, the gene density surrounding the window, and the 
	recombination rate in the window.
	inputs:  the fasta file handle, the window size, and the OUTput file handle
	output:  prints to the outfile (or stdout) all the data properly formatted
	returns: 0 if everything worked, 1 if something went wrong.		 
	'''
		print("Genome file:								", FASTA_FH, file=sys.stderr, end="\n")
		print("Output at: 								", OUT.name, file=sys.stderr)
		print("chr", "start", "stop", "ChrLen", "num_Ns", "GC", "R", "GeneDensity", sep=",", file = OUT)
	#try:
		with open(FASTA_FH, 'r') as FASTA:
			for record in SeqIO.parse(FASTA, 'fasta'):									#for each fasta entry 
				print("Parsing genome and calculating GC and gene density stats. Working on:	", record.name, "			 ", sep=" ", file=sys.stderr, end="\r")			
				fragLen      = len(record.seq)
				#recordname = record.name.lstrip("Chr") #was for the genetic map that had Chrs as"1" "2" etc. But the new map uses CHr1, Chr2 etc, so not needed. 
				#print(record.name, recordname, file=sys.stderr)
				for i in range( ceil( len( record.seq ) / windowSize ) ):			#this breaks it down to the window
					start        = (i * windowSize) + 1								#get start & stop
					stop         = min((i + 1) * windowSize, len(record.seq))						
					seq          = record.seq[start:stop + 1:1].upper() 				#get the sequence for GC calculation	
					GC           = calc_GC(seq) 										#get sequence GC
					GeneDensGene = calc_GeneDens(record.name, start, stop, fragLen, geneDensityWindowSize, gff_dict_gene)	#get gene density
					R            = recoord_R(record.name, start, stop)	#get recombination rate
					print(record.name,  start, stop, len(record.seq), seq.count("N"), GC, R, GeneDensGene, sep=",", file = OUT)#for each bit windowSize across the entry:		
		print("\n", file=sys.stderr)
		return(0)

	
def calc_GeneDens(CHR, start, stop, fragLen, size, dict):
	'''
	This takes in the location of the fragment that we're looking at and calculates the average gene density for it.
	Will count all coding bases 500,000 bases before the midpoint and 500,000 bases after the midpoint
	inputs: chr, start and stop
	returns: gene density in units of 'thousands of genic bases per million bases which is ~= to genes per million bases'
	'''
	differential     = size / 2 
	mid_base         = int((start + stop) / 2)
	start_base       = max(0, mid_base - differential) #the start of the chromosome or 0.5Mb before the window
	end_base         = min(mid_base + differential, fragLen) #either the end of the chr or 0.5Mb past the end of the fragment
	focalAreaLength  = end_base - start_base # to divide out the length and make sure the denominator is actually Megabases
	coding_bases     = get_coding_bases(CHR, start_base, end_base, dict)
	if coding_bases == "NA":
		return(coding_bases)
	else:
		return(coding_bases / focalAreaLength)	#this makes the units likelihood of a base being genic


def get_coding_bases(CHR, start_base, end_base, dict):
	'''
	This should go through the gff_dict and count how many coding bases 
	are within the start:end range
	'''
	bases = 0
	if CHR not in dict: # this means there are no coding bases in the ECR
		return("NA")
	else:
		for START in dict[CHR]:
			END = dict[CHR][START]
			if (START < start_base and END > end_base):
				bases = bases + ( end_base - start_base) #if the gene encompasses the whole ROI, add the length of the roi.
			elif ( START < start_base and END > start_base ): #if the gene starts before the roi and crosses the start base (if it is on the far side of the end base, the first if would have been triggered instead)
				bases = bases + ( END - start_base)
			elif ( START > start_base and END < end_base ): #if the fragment is completely within the ROI
				bases = bases + (END - START)
			elif ( START < end_base and END > end_base): #if the fragment starts in the ROI and leaves it. 
				bases = bases + ( end_base - START)			
	return(bases)			   

def recoord_R(chr, start, stop):
	'''
	this bit will take in the chr, start, and stop coordinates and return the recombination 
	rate for those bases.
	'''
	R_at_Start,breakpoint  = extract_R(chr, start) #this grabs the breakpoint which is a stop base of the earliest fragment
	R_at_End,_    = extract_R(chr, stop)
	if R_at_Start == "NA" or R_at_End =="NA":
		return("NA")
	elif R_at_Start == R_at_End: 
		return(round(R_at_Start, 6))
	else: #get the breakpoint and calculate the weightings for a weighted average
		weight_Start = breakpoint - start 
		weight_End   = stop - breakpoint
		R = (R_at_Start * weight_Start) + (R_at_End * weight_End) / (weight_End + weight_Start) #if the window size is big, it may cross over multiple blocks of recombination. That would be bad and lead to the length of the middle one being counted as the R of the last one. But unless the window is really big, it should be fine.  
		return(round(R, 6))


def extract_R(chr, pos):
	'''
	grabs the recombination rate from the dictionary given a chr and a position
	inputs: chr and pos
	outputs: none
	returns: R at that position
	'''
	if chr in map_dict:
		for start in map_dict[chr]:				
			if pos > start and pos <= map_dict[chr][start]["stop"]:
				return(map_dict[chr][start]["R"], map_dict[chr][start]["stop"])
		return(0, map_dict[chr][start]["stop"]) #so in this case, return 0 for and the mid-point base. It will likely be the R_at_End that triggers this and we should assume that the R there is 0. 
	else:
		return("NA", 0)
		
def calc_GC(str):
	'''
	This takes in a string of DNA and returns the GC percent as a decimal rounded to 4 places
	If the length of the string is 0, or if the string is all Ns it returns 'NA'
	'''
	if len(str)==0:
		return("NA")
	else:	
		GC = str.count("G") + str.count("C")
		AT = str.count("A") + str.count("T")
		N = len(str) - GC - AT
		if N == len(str):
			perGC = "NA"
		else:
			perGC = round(GC / ( GC + AT ), 4)
		return( perGC )
			
def read_in_map(MAP_FH):
	'''
	This bit reads in the genetic map file (which must be in genome coordinates), 
	NOT ECR coordinates and then calculates the recombination rate between every marker.
	These recombination rates can then be extracted for each ECR.
	inputs: genetic map filehandle
	outputs: updated dictionary "map_dict" with the following structure Chr:{start:{"stop":stop, "rate":R}}
	returns: 0 or 1 for success or failure. 
	'''
	print("Parsing the genetic map and calculating recombination rate. File:	", MAP_FH, file=sys.stderr)
	with open(MAP_FH, 'r') as MAP:
		prevchr = "Chr1"
		prevpos = -1
		prevcM = 0
		counter = 0
		for line in MAP:
			counter += 1
			if counter > 1:
				line=line.strip("\n")
				bits = line.split(",")  
				mname = bits[0]
				chr = bits[1] 
				pos = int(float(bits[2]))
				cM = float(bits[3]) 
				#print(mname, chr, pos, cM, prevpos, prevcM)
				R = calc_R(chr, prevchr, pos, prevpos, cM, prevcM)
				map_dict[chr] = map_dict.setdefault(chr, {0:{"stop":pos, "R":"NA"}})
				if prevchr == chr:
					map_dict[chr].update({prevpos:{"stop": pos, "R": R}})
				#print(chr, prevpos, pos, R)
				prevchr = chr
				prevpos = pos
				prevcM = cM
	#pprint.pprint(map_dict)
	return(0) #0 means success	
		
def calc_R(chr, prevchr, pos, prevpos, cM, prevcM):
	'''
	this bit will take in a bunch of data and return a recombination rate in cM/Mb
	inputs: chr the chr that we're on
			prevchr the previous chromosome, if this is not the same, then return NA
			pos and prevpos are positions to give the denominator
			cM and prevcM are the centimorgans for the numerator
	outputs: none
	returns: the recombination rate in cM/Mb		
	'''
	if chr == prevchr:
		R = (cM - prevcM) / (pos - prevpos) * 1000000
		return(round(R, 6))
	else:		
		return("NA")	
		
##########################################################################################
#main body
##########################################################################################
def main():
	args = process_command_line()
	#print(args)
	
	if args.Out:
		OUT = open(args.Out, "w")
		OUT_FH = args.Out
	else:
		OUT = sys.stdout 
		OUT_FH = "stdout"
	'''
	first read in the gff file and make it into a dictionary
	it will be used to look up the number of coding bases in a region
	'''
	if read_in_gff(args.GFF): #this will return 0 if it worked, if it returns anything else, that means there was an error
		print("\n	Error reading in the GFF file. Exiting now", file=sys.stderr)
		return(1)


	'''
	second read in the genetic map and calculate the recombination rate across the SFs
	store recombination rate in some meaningful way...?
	'''	
	if args.Map:		
		if read_in_map(args.Map):
			print("\n	Error reading in the map file. Exiting now", file=sys.stderr)	
			return(2)


	'''
	finally read in the genome and start scanning across it in the windows.
	calculate the gc content as I go and also the gene density, and finally the
	'''	
	if read_genome_calcAll(args.fasta, args.Window, OUT, args.geneDensityWindowSize):
		print("\n	Error reading in the genome or calculating something. Exiting now.", file=sys.stderr)
		return(4)	#returning !0 means something went wrong  
	
	
	
	if args.Out:	
		OUT.close()	
	return(0) #returning 0 means success of this subroutine and all that is in it.

if __name__ =='__main__':
	status = main() # will return 0 upon success
	sys.exit(status)



