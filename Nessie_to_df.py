#!/usr/bin/env Python3

##########################################################################################
#docstring
##########################################################################################
desc='''
Written by TDB on 14 June 2021
Should rip through the output of Nessie and make the pseudo-fasta format into a proper data table. 

'''


##########################################################################################
#imports
##########################################################################################
import argparse
import sys
import pprint
import re

##########################################################################################
#functions
##########################################################################################
def process_command_line():
	'''
	This just uses the argparse module to get the command line arguments -
	looking for a fasta, a GFF, a Map, a window size, and possibly an outfile
	'''
	parser = argparse.ArgumentParser(description=desc)                                            
	parser.add_argument("--Nessie",                    "-N", type=str, required=True,  help="Output from the nessie command", metavar="*.txt")
	parser.add_argument("--Out",                      "-o", type=str, required=False, help="Output file name. Omit for stdout")
	args = parser.parse_args()
	return(args)	


def main():
	args = process_command_line()
	#print(args)
	
	if args.Out:
		OUT = open(args.Out, "w")
		OUT_FH = args.Out
	else:
		OUT = sys.stdout 
		OUT_FH = "stdout"
	
	if re.search("_E_", args.Nessie):
		type = "entropy,A,C,G,T"
	elif re.search("_L_", args.Nessie):
		type = "LinComp"
	else:
		print("the regex search didn't work. Stopping now.", file=sys.stderr)
		exit()	
	
	with open(args.Nessie, "r") as NESSIE:
		CHR = "not_yet_set"
		SHORTCHR = "not_yet_set"
		offset = "not_yet_set"
		print("Chr,Scaf,Start", type, file=OUT, sep=",")
		for line in NESSIE:
			line=line.strip("\n")
			if not (line.startswith("#")):
				if line.startswith(">"):
					CHR = line.lstrip(">")
					CHR = CHR.rstrip("\n")
					CHR = CHR.split(" ")[0]
					match = re.search("(Chr\w+?)[_,]", CHR)
					if match:
						SHORTCHR = match.group(1)
					else:
						SHORTCHR = CHR	
				elif line.startswith("@"): # for some reason they report the start in blocks, there will be an @0-1234004 and then data blocks at 1000, 2000, 3000 etc. once it gets to 12340004, there is another @1235004-2345005 and the next data will start again at 0. Seems a bit dumb, but easy to reset everything by the offset. 
					offset = int(line.strip("@").split("-")[0])
				else:
					bits =line.split("	")
					start = int(bits[0]) + offset 
					data = bits[1]
					print(SHORTCHR, CHR, start, data, sep=",", end="\n", file=OUT)	
	
	if args.Out:	
		OUT.close()	
	return(0) #returning 0 means success of this subroutine and all that is in it.


if __name__ =='__main__':
	status = main() # will return 0 upon success
	sys.exit(status)

