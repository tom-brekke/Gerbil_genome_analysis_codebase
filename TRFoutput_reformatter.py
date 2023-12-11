#!/usr/bin/env Python3


import argparse
import sys
import re
import pprint

cent_dict = {}

def process_command_line():
	'''
	This just uses the argparse module to get the command line arguments -
	looking for a *.dat file and an outfile.
	'''
	parser = argparse.ArgumentParser(description="written by TDB\n20-09-2021\nTakes a *dat file from TandemRepeatFinder and reformats it into a proper data frame with the chr on every line. Also exports fastas for every 6, 36, 127, and 1747bp sequences.")                                              
	parser.add_argument("--dat", "-d", type=str, required=True, help="The output file from TandemRepeatFinder with flags d and h used.", metavar="*.dat")
	parser.add_argument("--centromeres", type=str, required=False, help="The centromere locations file. Default = %(default)s", default = "cent_locations_v5.2.csv", metavar="*.csv")
	parser.add_argument("--mincopies", type=int, required=False, help="The minimum number of copies that a repeat has to have to be output. There are many repeats that are not centromere repeats and are at low copy number (1-50) which won't align with the true centromere repeats, this is the easiest way to get rid of them. Default = %(default)s", default = 100)
	parser.add_argument("--Out", "-o", type=str, required=True, help="Output file stem. There are a couple output files: the reformatted dataframe, and a fasta of all the 6bp, 36bp, 127bp, and 1747bp sequences.")
	parser.add_argument("--consensus", "-c", default=False, action='store_const', const=True, required=False, help="Include to extract only consensus sequences of each length, not every sequence. Will save a huge amount of time on the next alignment step, but may miss variants that are common but not represented in the consensus sequences.")
	parser.add_argument("--simplify", "-s", default=False, action='store_const', const=True, required=False, help="TandemRepeatFinder will report series of repeats of period sizes, for instance if a true repeat is 20 copies of 6bp long, it will report the 20x6bp repeat, but also a 10x12bp repeat and perhaps a 5x24bp repeat. In other words - it will report the exact same sequence multiple times. Include this flag to attempt to choose the best (typically shortest period size) of these to output. Omitting this flag will output everything.")

	args = parser.parse_args()
	return(args)	


def main():
	args = process_command_line()
	#These are base location sof the centromeres. they are approximate, but that is sufficeient for my purposes now.
	#cent_dict_v5 = {"Chr1":(119800000,124000000),"Chr2":(101100000,105200000),"Chr3":(87500000,88900000),"Chr4":(73300000,78000000),"Chr5":(75800000,77250000),"Chr5_unordered_Scaffold_22":(0,30010000),"Chr6":(73900000,75600000),"Chr7":(70800000,72100000),"Chr8":(64500000,65800000),"Chr9":(79800000,81200000),"Chr10":(44700000,46200000),"Chr11":(51900000,53300000),"Chr12":(47500000,49000000),"Chr13_unordered_Scaffold_34":(2000000,6000000),"Chr14":(56100000,60000000),"Chr15":(50000000,51500000),"Chr16":(38300000,39800000),"Chr17":(0,3400000),"Chr18":(71800000,76500000),"Chr19":(67800000,72450000),"Chr20":(68000000,73400000),"Chr21":(56000000,58800000),"ChrX":(154000000,160100000),"ChrY_unordered_Scaffold_35":(1500000,2800000),}
	#cent_dict_v52 = {"Chr1":(75369399,79569399),"Chr2":(86002158,90102158),"Chr3":(87500000,88900000),"Chr4":(63826194,68526194),"Chr5":(75800000,77250000),"Chr5_unordered_Scaffold_22":(0,30010000),"Chr6":(61714327,63414327),"Chr7":(50484302,51784302),"Chr8":(57913762,59213762),"Chr9":(39429165,40829165),"Chr10":(44700000,46200000),"Chr11":(51900000,53300000),"Chr12":(47500000,49000000),"Chr13_unordered_Scaffold_34":(2000000,6000000),"Chr14":(31019591,34919591),"Chr15":(27421940,28921940),"Chr16":(38300000,39800000),"Chr17":(0,3400000),"Chr18":(865389,5565389),"Chr19":(0,4319935),"Chr20":(0,5204060),"Chr21":(0,2726396),"ChrX":(0,6106605),"ChrY_unordered_Scaffold_35":(1500000,2800000),}
	read_in_cent_dict(args.centromeres)

	
	if args.simplify:
		return(simple(args))
	else:
		return(complete(args))
		
		

def read_in_cent_dict(centFH):
	print("Reading in the centromere locations file:", centFH, sep=" ", file=sys.stderr)
	with open(centFH) as CENT:
		lineCounter = 0
		for line in CENT:
			lineCounter += 1
			if lineCounter < 1:
				continue
			bits = line.split("	")
			if bits[0] not in cent_dict:
				#cent_dict[bits[0]] = int(bits[1]), int(bits[2])
				cent_dict[bits[0]] = {}
			if int(bits[1]) not in cent_dict[bits[0]]:
				cent_dict[bits[0]][int(bits[1])] = int(bits[2])
	#pprint.pprint(cent_dict)				

def simple(args):
	print("Reading in the data file:", args.dat, "and processing it using the -s flag.", sep=" ", file=sys.stderr)

	#print(args)
	dict = {}
	
	#read dat into a dictionary
	with open(args.dat, "r") as DAT:
		oldStop = 0
		skipped = 0
		added = 0
		for line in DAT:			
			line=line.strip("\n")
			#print(line, file=sys.stderr)
			if line.startswith("Sequence"):
				Chr = line.split(" ")[1].rstrip("	EXTRACTED")
				if Chr =="Chr":
					Chr="ChrX" #because removing the EXTRACTED will drop off the X too..
				#print(Chr, file=sys.stderr)
				oldStop = 0
			elif re.match("^\d", line):	
				#prints new df
				#print(line, file=sys.stderr)
				Start,Stop,period,numCopies,consensusSize,percentMatch,percentIndel,alignmentScore,A,C,G,T,entropy,consensus,alignment = line.split(" ")
				Start = int(Start)
				Stop = int(Stop)
				if Start > oldStop:
					added = added + 1
					if Chr not in dict:
						dict[Chr] = {Start:line}
					else:
						dict[Chr][Start] = line	
				else:
					skipped = skipped + 1	
				oldStop = max(Stop, oldStop)
				#extract fasta of the repeats for alignment and DNA logo/consensus calling.

	#print dictionary
	if args.consensus:
		OUT     = open(args.Out+".consensus.simplified.csv", "w")
		OUT6    = open(args.Out+".6bp.consensus.simplified.fasta", "w")
		OUT36   = open(args.Out+".36bp.consensus.simplified.fasta", "w")
		OUT127  = open(args.Out+".127bp.consensus.simplified.fasta","w")
		OUT1747 = open(args.Out+".1747bp.consensus.simplified.fasta","w")
	else:
		OUT     = open(args.Out+".simplified.csv", "w")
		OUT6    = open(args.Out+".6bp.simplified.fasta", "w")
		OUT36   = open(args.Out+".36bp.simplified.fasta", "w")
		OUT127  = open(args.Out+".127bp.simplified.fasta","w")
		OUT1747 = open(args.Out+".1747bp.simplified.fasta","w")

	#pprint.pprint(dict.keys())
	print("Writing the data files", sep=" ", file=sys.stderr)
	print("Chr,Start,Stop,period,numCopies,consensusSize,percentMatch,percentIndel,alignmentScore,A,C,G,T,entropy,consensus,alignment", file=OUT)
	for Chr in sorted(dict):
		for start in sorted(dict[Chr]):		
			#print(Chr, start, file=sys.stderr)
			print(Chr, ",".join(dict[Chr][start].split(" ")), sep=",", file=OUT)
			extract_fasta(Chr, dict[Chr][start], OUT6, OUT36, OUT127, OUT1747, args.consensus, args.mincopies)
			
	OUT.close()	
	OUT6.close()	
	OUT36.close()	
	OUT127.close()	
	OUT1747.close()	
	
	
	return(0) #returning 0 means success of this subroutine and all that is in it.
	
	
	
def complete(args):	
	#print(args)
	print("Reading in the data file:", args.dat, "and processing it without the -s flag.", sep=" ", file=sys.stderr)

	if args.consensus:
		OUT     = open(args.Out+".consensus.csv", "w")
		OUT6    = open(args.Out+".6bp.consensus.6bp.fasta", "w")
		OUT36   = open(args.Out+".36bp.consensus.fasta", "w")
		OUT127  = open(args.Out+".127bp.consensus.fasta","w")
		OUT1747 = open(args.Out+".1747bp.consensus.fasta","w")
	else:
		OUT     = open(args.Out+".csv", "w")
		OUT6    = open(args.Out+".6bp.fasta", "w")
		OUT36   = open(args.Out+".36bp.fasta", "w")
		OUT127  = open(args.Out+".127bp.fasta","w")
		OUT1747 = open(args.Out+".1747bp.fasta","w")



	#read dat
	print("Chr,Start,Stop,period,numCopies,consensusSize,percentMatch,percentIndel,alignmentScore,A,C,G,T,entropy,consensus,alignment", file=OUT)
	with open(args.dat, "r") as DAT:
		for line in DAT:			
			line=line.strip("\n")
			#print(line, file=sys.stderr)
			if line.startswith("Sequence"):
				Chr = line.split(" ")[1].rstrip("	EXTRACTED")
				if Chr =="Chr":
					Chr="ChrX" #because removing the EXTRACTED will drop off the X too..
				#print(Chr, file=sys.stderr)
				#print(Chr, file=sys.stderr)
			elif re.match("^\d", line):	
				#prints new df
				print(Chr, ",".join(line.split(" ")), sep=",", file=OUT)
				#extract fasta of the repeats for alignment and DNA logo/consensus calling.
				extract_fasta(Chr, line, OUT6, OUT36, OUT127, OUT1747, args.consensus, args.mincopies)
				
				
			
	OUT.close()	
	OUT6.close()	
	OUT36.close()	
	OUT127.close()	
	OUT1747.close()	
	
	
	return(0) #returning 0 means success of this subroutine and all that is in it.

def extract_fasta(Chr, line, OUT6, OUT36, OUT127, OUT1747, print_consensus, mincopies):
	#example line:
	#3679 3724 20 2.5 17 80 16 56 82 17 0 0 0.67 AAAAAAAAACAAAACAA AAAACAAAACAAAACAAAAAACCAAAAACAACAACAAAAAAAAAAA
	#print(line, file=sys.stderr)
	start, stop, period, copies, consensusSize, percentMatches, percentIndels, alignmentScore, A, C, G, T, entropy, consensus, fullSequence  = line.split(" ")
	start = int(start)
	stop = int(stop)
	copies = int(float(copies))
	flag = 0 # 0 or 1 for whether it worked or didn't
	
	if Chr not in cent_dict: #bail out if there is no location for a centromere for the chromosome we're on. happens for unplaced scaffolds which need to be ignored.
		return(1) #bail out - no fastas to write
	#print(Chr, cent_dict[Chr], start, stop, copies)
	for startpos in cent_dict[Chr]:
		if ((start > startpos) & (stop < cent_dict[Chr][startpos]) & (copies > mincopies)):#only print repeats that are in the cent sequence. and at high copy number. 
			#print("yes")
			period = int(period)
			#print(period, file=sys.stderr)
	
			if period == 1747: #((period > 1700) & (test_remainder(period, 1747, 20))):
				#print(">", Chr,":", start, "-", stop, " period=", period, " copies=", copies, " perT=", T, "\n", consensus, sep="", end="\n", file=OUT1747)
				print_seqs(Chr, start, stop, period, consensusSize, copies, fullSequence, consensus, OUT1747, print_consensus)

		
			elif period == 127: #(period > 115) & (test_remainder(period, 127, 5)):
				#print(">", Chr,":", start, "-", stop, " period=", period, " copies=", copies, " perT=", T,  "\n", consensus, sep="", end="\n", file=OUT127)
				print_seqs(Chr, start, stop, period, consensusSize, copies, fullSequence, consensus, OUT127, print_consensus)


			elif period == 36 or period == 37:#(period > 28) & (test_remainder(period, 36, 4)):
				if (not re.search("TTAGGG", consensus) and not re.search("TCTCTCTC", consensus) and not re.search("CACTCACT", consensus)):		#because the telomere repeats can be classed as period 36 and are TTAGGG
					#print(">", Chr,":", start, "-", stop, " period=", period, " copies=", copies, " perT=", T,  "\n", consensus, sep="", end="\n", file=OUT36)
					print_seqs(Chr, start, stop, period, consensusSize, copies, fullSequence, consensus, OUT36, print_consensus)

				
				
			elif period == 6: #test_remainder(period, 6, 1):
				#print(">", Chr,":", start, "-", stop, " period=", period, " copies=", copies, " perT=", T,  "\n", consensus, sep="", end="\n", file=OUT6)
				print_seqs(Chr, start, stop, period, consensusSize, copies, fullSequence, consensus, OUT6, print_consensus)
		
			else:
			#	print("none", file=sys.stderr)
				flag = 1 # not a repeat I'm interested in.
		else:
			flag = 1
	
	return(flag)


def print_seqs(Chr, start, stop, period, consensusSize, copies, fullSequence, consensus, OUT, print_consensus):
	if print_consensus: #prints just the consensus - printing all the alignments would capture all the variation, but the alignment of them all takes quite a long time.
		T = calcT(consensus)
		if T < 0.30:
			consensus = revcom(consensus)
			T = calcT(consensus)
		print(">", Chr,":", start, "-", stop, " period=", period, " consensusSize=", consensusSize, " copies=", copies, " perT=", T,  "\n", consensus, sep="", end="\n", file=OUT)
	else:
		copies = int(float(copies))
		period = int(period)
		consensusSize = int(consensusSize)
		rep_len = consensusSize # change between consensusSize and period
		for i in range(0,copies):
			b = i * rep_len
			e = (i+1) * rep_len
			if e < len(fullSequence):
				seq = fullSequence[b:e:1]
				T = calcT(seq)
				if T< 0.30: #put everything in the forward direction. arbitraily defined as T% greater than 30%
					seq = revcom(seq)
					T = calcT(seq)
				print(">", Chr,":", start + b, "-", start + e, " period=", period, " consensusSize=", consensusSize, " copy=", i+1, "of", copies, " perT=", T,  "\n", seq, sep="", end="\n", file=OUT)



def test_remainder(testValue, truePeriod, Threshold):
	#print(testValue, truePeriod, Threshold, file=sys.stderr)
	Remain = testValue % truePeriod
	#print("	", Remain,file=sys.stderr)
	if Remain < Threshold:
		return(True)
	elif (truePeriod - Remain < Threshold):
		return(True)
	else:
		return(False)

def revcom(seq):
	d = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N"}
	seq = seq.upper()[::-1]
	return("".join([d[base] for base in seq]))
	



def calcT(seq):
	A = seq.upper().count("A")
	C = seq.upper().count("C")
	G = seq.upper().count("G")
	T = seq.upper().count("T")
	if A+C+G+T == 0: #runs of all N's  - shouldn't be (m)any, and this just kinda bypasses them.
		return(.25)
	return(round( (T)/(A+C+G+T), 2))
	
	
if __name__ =='__main__':
	status = main() # will return 0 upon success
	sys.exit(status)

