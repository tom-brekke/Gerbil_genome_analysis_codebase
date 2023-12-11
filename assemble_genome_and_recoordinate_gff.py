#!/usr/bin/env Python3


import argparse
import sys
import pprint
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

parser = argparse.ArgumentParser(description="written by TDB on 10-02-2021. To assemble a genome from scaffolds and recoordinates the associated gff file appropriately.")
parser.add_argument("--Fasta", "-f", type=str, required=True, help="The fasta file for the original genome.", metavar="*.fasta")
parser.add_argument("--GFF", "-g", type=str, required=False, help="GFF file for the original genome. Optional.", metavar="*.gff")
#Map was sorted in DTMap.qtl, so this shouldn't be run: parser.add_argument("--Map", "-m", type=str, required=False, help="A genetic map file for the original genome. Columns must be: SnpID, chr, base_pos, (female_)cM", metavar="*.csv")
parser.add_argument("--Assembly", "-a", type=str, required=True, help="The *.assembly.csv for the genome assembly you are interested in. This is a csv with 7 columns where the first base is '1': Scaffold, scaf_start, scaf_end, Chr, chr_Start, chr_end, orientation. The assembly file needs in the 'Scaffold' column to have all of the genome scaffolds including any scaffolds that will appear as a result of all the breaks. It must include scaffolds that need to be dropped. The genome will be printed in the order that the chromosomes appear in this file.")
parser.add_argument("--Breaks", "-b", type=str, required=False, help="The *.breaks.csv for the genome assembly you are interested in. This is a csv with 6 columns: Scaffold, scaf_start, scaf_end, Chr, chr_Start, chr_end. Each scaffold should be there twice and will be broken into the output chromosomes.")
parser.add_argument("--scaffolds_to_Drop", "-d", type=str, required=False, help="A file of scaffold names to remove from the assembly and gff. One name per line. Probably used for the case wehre there were 5 scaffs that were all the mt genome and they needed to be removed.")
parser.add_argument("--Pools", "-p", type=str, required=False, help="The *pools_assignments.csv for the genome assembly you are interested in. This is the output from chromosome_pools_Dovetail.r and has the pool assignments for each scaffold in it - used here to name unplaced scaffolds.")
parser.add_argument("--scaffold_name_to_Trim", "-s", type=str, required=False, help="If scaffold names have an extra bit (i.e. '__1_contigs__length_1234567') that you don't want to use, type it here and it will be removed from the end of the name only. Regex allowed. Perhaps: __\d+_contigs__length_\d+")
parser.add_argument("--Out", "-o", type=str, required=False, help="Output filestem, will make a <stem>.fasta, a <stem>.gff, and a meta-data file: <stem>.meta.txt.")

args = parser.parse_args()


def parse_scaffolds_to_Drop(scaffolds_to_Drop):
	print("Parsing the scaffolds to drop file:	", scaffolds_to_Drop.strip("\n"), file=sys.stderr)
	with open (scaffolds_to_Drop, "r") as StD:
		for line in StD:
			line=line.strip("\n")
			if line not in scaffolds_to_drop_dict:
				scaffolds_to_drop_dict[line] = 0
			else:
				print(scaffolds_to_Drop, "has duplicates. This isn't really a problem, just something you should probably know.", file=sys.stderr)
				
def parse_assembly(Assembly_FH):
	placed = True
	print("Parsing the assembly file:	", Assembly_FH, file=sys.stderr)
	with open(Assembly_FH, "r") as AFH:
		counter=0
		for line in AFH:
			counter += 1
			if counter > 1:
				line=line.strip("\n")
				#pprint.pprint(line)
				scaffName, scafStart, scafEnd, Chr, start, end, dir = line.split(",")
				scafStart = int(scafStart)
				scafEnd = int(scafEnd)
				name=""
				if not start =="NA": #this means that the scaffold is placed in the chromosome
					start = int(start)
					end = int(end)
					placed = True
					name = "Chr"+ Chr.lstrip("Chr")
				else: #this means it's unplaced
					start = scafStart
					end = scafEnd
					dir = "+"
					placed = False
					name = "Chr"+ Chr + "_unordered_" + scaffName
				if name in scaffolds_to_drop_dict:
					placed = "dropped_from_assembly"
				if scaffName not in assembly_dict:
					assembly_dict[scaffName] = {"Chr":Chr, "scaffStart":scafStart, "scaffEnd":scafEnd, "scaffDir":dir, "ChrStart":start, "ChrEnd":end, "Seq":"", "Pool":"", "placed":placed, "Name":name}
					if name not in chr_dict:
						chr_dict[name] = {}
						if name not in scaffolds_to_drop_dict:
							chr_names.append(name)	
					chr_dict[name][start] = scaffName
				else:
					print("The scaffold", scaffName, "is in the genome twice! Check the assembly file:", args.Assembly, "! Exiting now", file=sys.stderr)
					exit()	
	#pprint.pprint(assembly_dict)
	
def parse_pools(Pools_FH): # this file has a line for every marker, but whatever, it also has what we need, just don't add to dicts for each line.
	print("Parsing the pools file:		", Pools_FH, file=sys.stderr)
	pool_dict = {"pool01":"X", "pool02":"2", "pool03":"1or3", "pool04":"5or6", "pool05":"4", "pool06":"9", "pool07":"8or10", "pool08":"7", "pool09":"11", "pool10":"13", "pool11":"12or17", "pool12":"14", "pool13":"Y", "pool14":"15or18", "pool15":"19or20", "pool16":"16", "pool17":"21"}
	with open(Pools_FH, "r") as PFH:
		counter=0
		for line in PFH:
			counter += 1
			if counter > 1:
				line=line.strip("\n")
				#pprint.pprint(line)
				Scaf,Len,dump,Ns,GC,nonN,pool01,pool02,pool03,pool04,pool05,pool06,pool07,pool08,pool09,pool10,pool11,pool12,pool13,pool14,pool15,pool16,pool17,Pool,Pools,color = line.split(",")#nearly all of these won't be used, but whatever
				#print(Scaf)
				if Scaf in assembly_dict:
					assembly_dict[Scaf]["Pool"]=Pool
				elif Scaf in breaks_dict:
					#these are scaffolds that were broken. I need to look up their new names in breaks dict and apply the proper pool to the new names in assembly_dict
					for outputScaf in breaks_dict[Scaf]:
						assembly_dict[outputScaf]["Pool"]=Pool
				else:
					print("There's something wrong: the pools file has scaffolds that are not present in either the assembly or breaks files. Exiting now.", file=sys.stderr)
					exit()
	#pprint.pprint(assembly_dict)				


							
def recoord_item(scaffName, geneStart, geneEnd, geneDir):
	Chr = assembly_dict[scaffName]["Chr"]
	#print(assembly_dict[scaffName]["scaffDir"], file=sys.stderr)
	newStart = "notCorrectYet"
	newEnd = "notCorrectYet"
	newDir = geneDir
	if  assembly_dict[scaffName]["scaffDir"] == "+" and geneDir == "+":
		newStart = geneStart + assembly_dict[scaffName]["ChrStart"]
		newEnd   = geneEnd   + assembly_dict[scaffName]["ChrStart"]
		newDir   = "+" #+
		#print(scaffName, "+/+", file=sys.stderr)
	elif  assembly_dict[scaffName]["scaffDir"] == "+" and geneDir == "-":
		newStart = geneStart + assembly_dict[scaffName]["ChrStart"]
		newEnd   = geneEnd   + assembly_dict[scaffName]["ChrStart"]
		newDir   = "-" #-
		#print(scaffName, "+/-", file=sys.stderr)
	elif  assembly_dict[scaffName]["scaffDir"] == "-" and geneDir == "+":
		newStart = assembly_dict[scaffName]["ChrEnd"] - geneEnd   + 1  #needs the +1 to avoid an off-by-one error but only for the -/+ and -/- cases. 
		newEnd   = assembly_dict[scaffName]["ChrEnd"] - geneStart + 1  #needs the +1 to avoid an off-by-one error but only for the -/+ and -/- cases. 
		#print(scaffName, "-/+", file=sys.stderr)
		newDir   = "-" #-
	elif  assembly_dict[scaffName]["scaffDir"] == "-" and geneDir == "-":
		#print(scaffName, "-/-", file=sys.stderr)
		newStart = assembly_dict[scaffName]["ChrEnd"] - geneEnd   + 1  #needs the +1 to avoid an off-by-one error but only for the -/+ and -/- cases. 
		newEnd   = assembly_dict[scaffName]["ChrEnd"] - geneStart + 1  #needs the +1 to avoid an off-by-one error but only for the -/+ and -/- cases. 
		newDir   = "+" #+
	return(Chr, newStart, newEnd, newDir)


def recoord_item_break(scaffNameOld, geneStartOld, geneEndOld):#this takes chrs that were broken and shifts the gene start appropriately. Then returns the new chr name, the new start, and the new stop
	#use breaks dict to find which new scaf the old ones go on. 
	for scaffNameNew in breaks_dict[scaffNameOld]:
		breakStart = breaks_dict[scaffNameOld][scaffNameNew][0]
		breakEnd   = breaks_dict[scaffNameOld][scaffNameNew][1]
		if geneStartOld > breakStart and geneStartOld < breakEnd:
			#this will be true when the gene start is indeed in the scaffNameNew we're working with. Now we need to test if the end base is too - otherwise we have split a gene with our break which would be bad.
			if geneEndOld > breakStart and geneEndOld < breakEnd:
				#this means that both the start and end of the gene are within one of the new breaks. Great! now convert the old start to a new start, convert the old end to a new end, and return those with the new scaffold name.
				newStart = geneStartOld - breakStart + 1 #needs the +1 to avoid an off-by-one error.
				newEnd = geneEndOld - breakStart + 1     #needs the +1 to avoid an off-by-one error. 
				#print(scaffNameOld, geneStartOld, geneEndOld, scaffNameNew, newStart, newEnd, file=sys.stderr)
				return(scaffNameNew, newStart, newEnd)
			else:
				print("One of the breaks splits up a gene! Not sure what to do here: which scaffold should that gene go on?? Should there really be a break there?? Going to panic and exit.", file=sys.stderr)
				exit()
		#scaffName, geneStart, geneEnd = recoord_item_break(scaffName, geneStart, geneEnd)




def recoordinate_gff(GFF_FH):
	GFF_OUT=open(args.Out+".gff",'w')
	print("Recoordinating gff:		", GFF_FH, file=sys.stderr)
	print("Output will go to:		", args.Out+".gff", file=sys.stderr)
	with open(GFF_FH, 'r') as GFH:
		counter = 0
		flag = 0
		for line in GFH:
			line=line.strip("\n")
			if line.startswith("##") and flag == 0:
				print(line, file = GFF_OUT)		
				print("#This GFF has been recoordinated by the script 'assemble_genome_and_recoordinate_gff.py' written by Thomas D Brekke", file = GFF_OUT)
				flag = 1		
			elif line.startswith("##") and flag == 1:
				a=1
			elif line.startswith("#"):	
				print(line, file = GFF_OUT)		
			else:
				scaffName,Gnomon,type,geneStart,geneEnd,dot1,geneDir,dot2,desc = line.split("	")
				if scaffName in scaffolds_to_drop_dict:
					continue
				if scaffName in assembly_dict:
					Chr, newGeneStart, newGeneEnd, newGeneDir = recoord_item(scaffName, int(geneStart), int(geneEnd), geneDir)
					line = "	".join([assembly_dict[scaffName]["Name"],Gnomon,type,str(newGeneStart),str(newGeneEnd),dot1,newGeneDir,dot2,desc])
					#print(scaffName,geneStart,geneEnd, geneDir, assembly_dict[scaffName]["scaffDir"], file=sys.stderr)	
				elif scaffName in breaks_dict:
					scaffName, geneStart, geneEnd = recoord_item_break(scaffName, int(geneStart), int(geneEnd))
					Chr, newGeneStart, newGeneEnd, newGeneDir = recoord_item(scaffName, int(geneStart), int(geneEnd), geneDir)
					line = "	".join([assembly_dict[scaffName]["Name"],Gnomon,type,str(newGeneStart),str(newGeneEnd),dot1,newGeneDir,dot2,desc])
				else:
					print("There's a problem with the gff recoordinate: the GFF file has scaffolds in neither the breaks nor the assembly file!", scaffName, "is the problem. Exiting now.", file=sys.stderr)
					exit()
				print(line, file=GFF_OUT)	
	GFF_OUT.close()

def recoordinate_map(MAP_FH): #this seems to not be used. 
	MAP_OUT=open(args.Out+".geneticmap.csv",'w')
	print("Recoordinating map:		", MAP_FH, file=sys.stderr)
	print("Output will go to:		", args.Out+".geneticmap.csv", file=sys.stderr)
	with open(MAP_FH, 'r') as MFH:
		counter = 0
		flag = 0
		for line in MFH:
			line=line.strip("\n")
			flag = flag + 1
			if flag == 1:
				print(line, file=MAP_OUT)
			elif flag > 1:
				Marker,Chr,ChrPos,cM = line.split(",")
				scaffName = "_".join(Marker.split("_")[0:2])
				#print(Marker, file=sys.stderr)
				#print(scaffName, file=sys.stderr)
				if scaffName in assembly_dict:
					Chr, newStart, newEnd, newDir = recoord_item(scaffName, int(ChrPos), int(ChrPos), "+") #
					#newMarkerName = "Chr"+Chr+"_"+newStart
					line = ",".join([Marker,Chr,str(newStart),cM])
					#print(scaffName,geneStart,geneEnd, geneDir, assembly_dict[scaffName]["scaffDir"], file=sys.stderr)	
				print(line, file=MAP_OUT)	
	MAP_OUT.close()



def parse_scaffolds(Fasta_FH):
	print("Parsing the fasta file:		", Fasta_FH, file=sys.stderr)
	with open(Fasta_FH, "r") as FFH:
		counter=0
		for record in SeqIO.parse(FFH, 'fasta'):
			counter+=1
			chr = record.name
			if args.scaffold_name_to_Trim:
					chr = trim(chr)
			print("	Chromosome: ", chr, "			", end="\r", file=sys.stderr)
			if chr in breaks_dict:
				#print(chr, "is in the breaks dict!")
				#now I need to take the entry and make new entries broken at the appropriate spot: one new entry for each value of the breaks_dict[chr]
				for outputChr in breaks_dict[chr]:
					start = int(breaks_dict[chr][outputChr][0])
					end = int(breaks_dict[chr][outputChr][1])
					#print(type(start), type(end))
					seqStr = str(record.seq[start:end])
					#print(seqStr[0:10], seqStr[-10:], file=sys.stderr) #this is a test the my indexing base-0 vs base-1 is done correctly.
					
					if assembly_dict[outputChr]["scaffDir"] =="-":
						assembly_dict[outputChr]["Seq"] = Seq(seqStr).reverse_complement()
					else:	
						assembly_dict[outputChr]["Seq"] = Seq(seqStr)
			else:
				if assembly_dict[chr]["Seq"] == "": #this is because the psammomys genome, there are two of each sequence name -one of the sequence, and the next of 1000 Ns which is the spacer for the next sequence. Need to skip that run of Ns and not replace it here. This may cause problems with genomes formatted in other ways.
					if assembly_dict[chr]["scaffDir"] =="-":
						assembly_dict[chr]["Seq"] = record.seq.reverse_complement()
					else:
						assembly_dict[chr]["Seq"] = record.seq
	print("", file=sys.stderr)
	
	
def trim(ch):
	pattern = "(.*)"+args.scaffold_name_to_Trim
	ch = re.fullmatch(pattern, ch)[1]		
	return(ch)		

def assemble_genome():
	FASTA_OUT=open(args.Out+".fasta",'w')
	print("Assembling the genome to:	", args.Out+".fasta", file=sys.stderr)
	for chr in chr_names:
		seq = ""
		for start in sorted(chr_dict[chr]): #sorted chr_dict makes it so that the order of the scaffolds is the order of the start positions, rather than the order of a dict (which is random)
			scaf = chr_dict[chr][start]
			seq = seq + "".join(["N" for i in range(100)]) + assembly_dict[scaf]["Seq"]
		seq=seq.lstrip("N")
		SR = SeqRecord(seq, id=assembly_dict[scaf]["Name"], description=assembly_dict[scaf]["Pool"])
		seqRecords.append(SR)	
	SeqIO.write(seqRecords, FASTA_OUT, "fasta")
	FASTA_OUT.close()
	
def write_metadata():
	META_OUT=open(args.Out+".meta.txt", "w")
	print("Writing metadata to:	", args.Out+".meta.txt", file=sys.stderr)
	print("Scaffold", "scaf_len", "Chr", "start", "stop", "orientation", "total_chr_len", "pool", "placed", sep=",", file=META_OUT)
	for scaf in assembly_dict:
		chrLen="NA"
		chr =  assembly_dict[scaf]['Chr']
		scaffEnd = assembly_dict[scaf]['scaffEnd']
		chrStart = assembly_dict[scaf]['ChrStart']
		chrEnd = assembly_dict[scaf]['ChrEnd']
		scaffDir = assembly_dict[scaf]['scaffDir']
		chrLen = calc_Chr_len(chr)
		pool = assembly_dict[scaf]['Pool']
		if pool=="":
			pool = "Unknown"
		placed = assembly_dict[scaf]["placed"]
		print(scaf, scaffEnd, chr, chrStart, chrEnd, scaffDir, chrLen, pool, placed, sep=",", file=META_OUT)  
	META_OUT.close()

def calc_Chr_len(Chr):
	scaf_set = get_all_scafs_on_chr(Chr)
	chr_len = 0
	for scaf in scaf_set:
		chr_len = chr_len + assembly_dict[scaf]['scaffEnd'] # scaffEnd is always the length of the scaffold
	chr_len = chr_len + (100 * (len(scaf_set)-1))		 
	return(chr_len)	

def get_all_scafs_on_chr(chr):
	scaf_set = set()
	#two places to look: in the Desc of assembly_dict and in the Chr or assembly dict
	for scaf in assembly_dict:
		if assembly_dict[scaf]["Chr"]==chr:
			scaf_set.add(scaf)
	return(scaf_set)

def parse_breaks(Breaks_FH):
	placed = True
	print("\nParsing the breaks file:	", Breaks_FH, file=sys.stderr)
	with open(Breaks_FH, "r") as BFH:
		counter=0
		for line in BFH:
			counter += 1
			if counter > 1:#skip the header line
				line=line.strip("\n")
				#pprint.pprint(line)
				#collect the scaffold to break ("inputScaf"), the locations from it to take ("inputStart" and "inputEnd") and the new scaffold name "outputScaf"
				inputScaf, inputStart, inputEnd, outputScaf, outputStart, outputEnd, = line.split(",")
				inputStart = int(inputStart)
				inputEnd = int(inputEnd)
				#print(line, file=sys.stderr)
				if inputScaf not in breaks_dict:
					breaks_dict[inputScaf] = {outputScaf:tuple([inputStart,inputEnd])}
				else:
					if outputScaf not in breaks_dict[inputScaf]:
						breaks_dict[inputScaf][outputScaf] = tuple([inputStart,inputEnd])
					else: 
						print("there is something wrong with the breaks file:", Breaks_FH, ". Multiple output scaffold names are the same - input scaffolds can exist multiply (and must at least be in there twice - or you will be dropping out parts of the genome...), but each output scaffold name must be unique. Exiting now.", sep="", file=sys.stderr)
						exit()
	#pprint.pprint(breaks_dict)
	

scaffolds_to_drop_dict = {} 
seqRecords = []	
chr_dict = {}
chr_names = [] # a list to store the names of the chrs in the proper order.
assembly_dict = {}
breaks_dict = {} #this has a key for every scaffold that needs be broken. the values are dictionaries with keys as the output chr names and the values as a start,end tuple where start,end are the coords from the input scaff. 
if args.scaffolds_to_Drop:
	parse_scaffolds_to_Drop(args.scaffolds_to_Drop)
if args.Breaks:
	parse_breaks(args.Breaks)
parse_assembly(args.Assembly)
if args.Pools:
	parse_pools(args.Pools)
parse_scaffolds(args.Fasta)
if args.GFF:
	recoordinate_gff(args.GFF)
#if args.Map:
#	recoordinate_map(args.Map)
#pprint.pprint(assembly_dict["Contig141263_ctg1"])
#pprint.pprint(assembly_dict["Contig141542_ctg1"])
#pprint.pprint(len(assembly_dict["Contig141542_ctg1"]["Seq"]))
#pprint.pprint(chr_dict)
assemble_genome()
write_metadata()


