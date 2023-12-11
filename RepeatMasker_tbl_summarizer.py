#!/usr/bin/env Python3


import argparse
import sys
import re
import os

def process_command_line():
	'''
	This just uses the argparse module to get the command line arguments -
	looking for a fasta that has been aligned.
	'''
	parser = argparse.ArgumentParser(description="written by TDB\n07-10-2021\nTakes a *fasta file from clustalo and reformats it by wrapping the dna around the base")                                              
	parser.add_argument("--directory", "-d", type=str, required=True, help="The RepeatMasker directory which contains *.tbl files to reformat and summarize", metavar="directory location")
	parser.add_argument("--Out", "-o", type=str, required=False, help="Output file name. Omit for stdout.")
	args = parser.parse_args()
	return(args)	
	
	
def extract_file_data(file):
	with open(file, "r") as FH:
		for line in FH:
			line = line.strip("\n")
			if line.startswith("file name:"):	
				scaf               = line[11:].strip()[:-11:1]
				Chr                = scaf.split("_")[0]
				line               = FH.readline()
				totalLength        = int(FH.readline().strip("\n")[14:].lstrip(" ").split()[0])
				GC                 = float(FH.readline().strip("\n")[9:].lstrip(" ").split()[0])
				line               = FH.readline().strip("\n")[13:].lstrip(" ").split()
				maskedCount        = int(line[0])
			elif line.startswith("SINEs"):
				SineCount          = int(line.split()[1])
				SineLength         = int(line.split()[2])
			elif line.startswith("LINEs"):
				LineCount          = int(line.split()[1])
				LineLength         = int(line.split()[2])
			elif line.startswith("LTR elements"):
				LTRCount           = int(line.split()[2])
				LTRLength          = int(line.split()[3])
			elif line.startswith("DNA elements"):
				DNACount           = int(line.split()[2])
				DNALength          = int(line.split()[3].strip("bp"))
			elif line.startswith("Unclassified"):
				UnclassCount       = int(line.split()[1])
				UnclassLength      = int(line.split()[2])
			elif line.startswith("Total interspersed repeats"):
				InterspersedLength = int(line.split()[3])
			elif line.startswith("Small RNA"):
				sRNACount          = int(line.split()[2])
				sRNALength         = int(line.split()[3])
			elif line.startswith("Satellites"):
				SatelliteCount     = int(line.split()[1])
				SatelliteLength    = int(line.split()[2])
			elif line.startswith("Simple repeats"):
				SimpleCount        = int(line.split()[2])
				SimpleLength       = int(line.split()[3])
			elif line.startswith("Low complexity"):
				LowComplexCount    = int(line.split()[2])
				LowComplexLength   = int(line.split()[3])
	#Chr = scaf		
	if Chr not in repeat_dict:
		repeat_dict[Chr] = {"totalLength":totalLength, "GC":GC, "maskedCount":maskedCount, "SineCount":SineCount, "SineLength":SineLength, "LineCount":LineCount, "LineLength":LineLength, "LTRCount":LTRCount, "LTRLength":LTRLength, "DNACount":DNACount, "DNALength":DNALength, "UnclassCount":UnclassCount, "UnclassLength":UnclassLength, "InterspersedLength":InterspersedLength, "SatelliteCount":SatelliteCount, "SatelliteLength":SatelliteLength, "SimpleCount":SimpleCount, "SimpleLength":SimpleLength, "LowComplexCount":LowComplexCount, "LowComplexLength":LowComplexLength}
	else: #update the entries with the other scaffolds on the same chr.
		old_len                                = repeat_dict[Chr]["totalLength"]
		repeat_dict[Chr]["totalLength"]        = old_len + totalLength
		repeat_dict[Chr]["GC"]                 = (repeat_dict[Chr]["GC"] * old_len) + (GC * totalLength) / repeat_dict[Chr]["totalLength"] #the weighted average of the GC contents should be the new GC content.
		repeat_dict[Chr]["maskedCount"]        = repeat_dict[Chr]["maskedCount"] + maskedCount
		repeat_dict[Chr]["SineCount"]          = repeat_dict[Chr]["SineCount"] + SineCount
		repeat_dict[Chr]["SineLength"]         = repeat_dict[Chr]["SineLength"] + SineLength
		repeat_dict[Chr]["LineCount"]          = repeat_dict[Chr]["LineCount"] + LineCount
		repeat_dict[Chr]["LineLength"]         = repeat_dict[Chr]["LineLength"] + LineLength
		repeat_dict[Chr]["LTRCount"]           = repeat_dict[Chr]["LTRCount"] + LTRCount
		repeat_dict[Chr]["LTRLength"]          = repeat_dict[Chr]["LTRLength"] + LTRLength
		repeat_dict[Chr]["DNACount"]           = repeat_dict[Chr]["DNACount"] + DNACount
		repeat_dict[Chr]["DNALength"]          = repeat_dict[Chr]["DNALength"] + DNALength
		repeat_dict[Chr]["UnclassCount"]       = repeat_dict[Chr]["UnclassCount"] + UnclassCount
		repeat_dict[Chr]["UnclassLength"]      = repeat_dict[Chr]["UnclassLength"] + UnclassLength
		repeat_dict[Chr]["InterspersedLength"] = repeat_dict[Chr]["InterspersedLength"] + InterspersedLength
		repeat_dict[Chr]["SatelliteCount"]     = repeat_dict[Chr]["SatelliteCount"] + SatelliteCount
		repeat_dict[Chr]["SatelliteLength"]    = repeat_dict[Chr]["SatelliteLength"] + SatelliteLength
		repeat_dict[Chr]["SimpleCount"]        = repeat_dict[Chr]["SimpleCount"] + SimpleCount
		repeat_dict[Chr]["SimpleLength"]       = repeat_dict[Chr]["SimpleLength"] + SimpleLength
		repeat_dict[Chr]["LowComplexCount"]    = repeat_dict[Chr]["LowComplexCount"] + LowComplexCount
		repeat_dict[Chr]["LowComplexLength"]   = repeat_dict[Chr]["LowComplexLength"] + LowComplexLength




#get all tbl files in the directory
args = process_command_line()

if args.Out:
	OUTFH = open(args.Out, "w")
else:
	OUTFH = sys.stdout	



repeat_dict = {}
for file in os.listdir(args.directory):
	if file.endswith(".tbl"):
		extract_file_data(file)
		
		
print("Chr", "Length", "GC", "maskedCount", "SineCount", "SineLength", "LineCount", "LineLength", "LTRCount", "LTRLength", "DNACount", "DNALength", "UnclassCount", "UnclassLength", "InterspersedLength", "SatelliteCount", "SatelliteLength", "SimpleCount", "SimpleLength", "LowComplexCount", "LowComplexLength", sep=",", file=OUTFH)
for Chr in repeat_dict:
	print(Chr, repeat_dict[Chr]["totalLength"], repeat_dict[Chr]["GC"], repeat_dict[Chr]["maskedCount"], repeat_dict[Chr]["SineCount"], repeat_dict[Chr]["SineLength"], repeat_dict[Chr]["LineCount"], repeat_dict[Chr]["LineLength"], repeat_dict[Chr]["LTRCount"], repeat_dict[Chr]["LTRLength"], repeat_dict[Chr]["DNACount"], repeat_dict[Chr]["DNALength"], repeat_dict[Chr]["UnclassCount"], repeat_dict[Chr]["UnclassLength"], repeat_dict[Chr]["InterspersedLength"], repeat_dict[Chr]["SatelliteCount"], repeat_dict[Chr]["SatelliteLength"], repeat_dict[Chr]["SimpleCount"], repeat_dict[Chr]["SimpleLength"], repeat_dict[Chr]["LowComplexCount"], repeat_dict[Chr]["LowComplexLength"], sep=",", file=OUTFH)
		
		
		
		
		
		
		