#!/usr/bin/env Python3

import argparse
import re
import select
import sys
import pprint


'''
class TimeoutExpired(Exception):
    pass
def input_with_timeout(prompt, timeout):
    sys.stdout.write(prompt)
    sys.stdout.flush()
    ready, _, _ = select.select([sys.stdin], [],[], timeout)
    if ready:
        return sys.stdin.readline().rstrip('\n') # expect stdin to be line-buffered
    raise TimeoutExpired
'''
    


parser = argparse.ArgumentParser(description="written by TDB on 09-12-2014. Updated by TDB on 22-11-2017 for use with any VCF. And updated again by TDB on 04 Feb 2021 to phase things by the parents genotypes.")                                              
parser.add_argument("--VCF", "-v", type=str, required=True, help="VCF file with genotypes. Default: %(default)s", default="gstacks.vcf", metavar="*.vcf")
parser.add_argument("--Pheno", "-p", type=str, required=False, help="If phenotypes need be added use this. It needs to be a csv formatted file with phenotypes in columns, first column must be the names of the individuals EXACTLY the same as in the VCF OR -t must be used, see below. This is a good way to add in sex, pgm, and any QTL trait data for rqtl", metavar="*.csv")
parser.add_argument("--Out", "-o", type=str, required=True, help="Output file stem. A csvr file will be output with the genotypes for r/qtl called '<stem>.rqtl.csvr' and a physical map will be output called '<stem>.pmap.csv'. ")
parser.add_argument("--individual_name_to_Trim", "-t", type=str, required=False, help="If vcf individual names have an extra bit (i.e. 'sorted.bam') that you don't want to use, type it here and it will be removed. Flat string only - no Regex. Perhaps: .sorted.RG.merged", metavar="STR")
parser.add_argument("--scaffold_name_to_Trim", "-s", type=str, required=False, help="If vcf scaffold names have an extra bit (i.e. '__1_contigs__length_1234567') that you don't want to use, type it here and it will be removed from the end of the name only. Regex allowed. Perhaps: __\d+_contigs__length_\d+", metavar="STR")
parser.add_argument("--only_first_SNP", default=None, action='store_const', const=True, required=False, help="Include to only use the first SNP per radtag. DO NOT USE WITH GENOME-ALIGNED DATA. Only for de novo RADs.")
parser.add_argument("--parents", default=None, type=str, required=False, help="If this is for a genetic map, it may help to only output SNPs that are fixed oppositely in the parents. Include one mother and one father as a comma-separated list and the script will only export sites that are alternatively fixed in those two individuals. Also - if this is included, it will polarise the allele calls so that A is from parent0 and B is from parent1")
parser.add_argument("--xchr", "-x", default=None, type=str, required=False, help="If you know which scaffolds are X linked, put them here as a comma separated list. It will rename them all to \"X\" and unless you have it recorded elsewhere you will loose the prior scaffold names. This may make the import to rqtl easier or it might just confuse things. Use with caution. Example: Scaffold_15,Scaffold_23,Scaffold_27")
args = parser.parse_args()

vcf = args.VCF
Pheno = args.Pheno
Out = args.Out+".rqtl.csvr"
pmap = args.Out+".pmap.csv"
OUT=open(Out, 'w')
PMAP=open(pmap, 'w')
iTrim = args.individual_name_to_Trim
sTrim = args.scaffold_name_to_Trim
parents = args.parents.split(",")
#pprint.pprint(parents)
parent_indexes = [] #This will hold the indexes of the parents in the name vector, needs to be global as it's accessed in a variety of definitions. Set in format_Header()
names=[]
cM=[0.0]
xcMs=[0.0]
old_ch=[""]
if args.xchr:
	xchrList = args.xchr.split(",")
else:
	xchrList = []
print("\n\n###########################\nRunning make_rqtl2_input_file.py", file=sys.stderr)
print("physical map will be output to		", pmap, file=sys.stderr)
print("genotypes will be output to		", Out, file=sys.stderr)



def format_Header(str): #first make the header row for the table: do this from the id's of the vcf file: ID,,ind1,ind2,ind3... followed by a row of sexes: sex,,1,1,1,2,2,1,1...
	print("Formatting header", file=sys.stderr)
	fields=str.split()
	print("ID,,", sep=",", end="", file=OUT)
	cols=0
	
	#need to get the indexes of the parents:
	f2 = [i.replace(iTrim,"") for i in fields]
	#pprint.pprint(f2)
	parent_indexes.extend([f2.index(p) for p in parents])
	#pprint.pprint(parent_indexes)
	
	#first line: individual names and make a list of the individuals to use later
	for idx in range(9, len(fields)):
		if idx not in parent_indexes:
			name=fields[idx]
			cols+=1
			if iTrim: #this gets rid of extra shit at the ends of VCF names if it's included
			#	print(name, file=sys.stderr)
				name=name.replace(iTrim,"") 
				#	print(name, iTrim, "\n", file=sys.stderr)
			names.append(name)

			print(",", name, sep="", end="", file=OUT)
	print("\n", end="", file=OUT)
	return(cols)

	
		
	
def	format_Pheno(cols):
	#second line: individual numbers - just need something
	print("Formatting phenotypes", file=sys.stderr)

	print("IndNum,,", sep=",", end="", file=OUT)
	for i in range(0,cols):		
		print(",", i+1, sep="", end="", file=OUT)
	print("\n", end="", file=OUT)
	#print(names, file=sys.stderr)
	#now any more phenotype lines that there may be:
	if Pheno: #skip this bit if there is no phenotype file
		Pheno_dict={}
		#parse it into dictionaries for each phenotype: this dict has each ID as a key and then a list of phenos, the list of phenos will be in the same order/position as the phenoNames
		with open(Pheno, 'r') as PHENO:	
			for line in PHENO:
				line=line.strip("\n")
				if line.startswith("id") or line.startswith("ID"):
					phenoNames=line.split(",")[1:]
				else:	
					fields=line.split(",")				
					Pheno_dict[fields[0]]=fields[1:]
		#pprint.pprint(Pheno_dict)
		#pprint.pprint(phenoNames)
		#now go through and print each phenotype based for each individual.
		for i in range(0,len(phenoNames)):
			print(phenoNames[i], end=",,", file=OUT)
			for j in range(0,len(names)):
				print(",",Pheno_dict[names[j]][i], sep="",end="", file=OUT)
			print("\n", end="", file=OUT)
			
def extract_parent_genos_and_compare(fields):
	#find column index of the parents
		#this is done in the header definition and stored in the global list parent_indexes, [idx1, idx2]
	#extract the genos at those indexes
	#pprint.pprint(parent_indexes)
	geno0 = fields[parent_indexes[0]].split(":")[0]
	geno1 = fields[parent_indexes[1]].split(":")[0]
	good_genos = ["0/0", "1/1"]			
	#print(geno0, geno1, file=sys.stderr)

	if geno0 in good_genos and geno1 in good_genos and geno0 != geno1:
		return(True)
	else:
		return(False)



def format_Geno(str, print_info):
	if print_info: # just makes it so the information 'formatting genotypes' is only printed once.
		print_info = False
		print("Formatting genotypes", file=sys.stderr)
	fields=str.split("	")
	#print(fields, file=sys.stderr)
	#test if geno is alternatively fixed in parents:
	if extract_parent_genos_and_compare(fields):

		ch=fields[0]
		if sTrim:
			pattern = "(.*)"+sTrim
			#print(ch, sTrim, pattern)
			ch = re.fullmatch(pattern, ch)[1]
			#print(ch)
		if ch in ch_dict: #the ch_dict keeps count of how many snps there are on each chr
			ch_dict[ch]+=1
			if args.only_first_SNP:
				return(print_info) #this bit here forces the script to skip 2nd and 3rd... SNPs per each radtag - may be messing up later steps which try to find those steps. 
		else:
			ch_dict[ch]=1
		pos=fields[1]
		name=ch+"_"+pos # this is a marker name. 
		#write PMAP physical map - pos needs to be in Mbp
		#three columns: marker name, chr, pos

		print(name, ch, pos, sep=",", file=PMAP)
	
		#write geno for parents and F2s
		if ch != old_ch[0]: #these are dummy centiMorgans for rqtl. rqtl needs a column of them, and this will set the first to 1, then 2, etc. rqtl will estimate the actual cMs. 
			cM[0] = -1
			#print(cM, xcMs)
			if ch in xchrList:
				cM[0] = xcMs[0] 
		cM[0] = round(cM[0] + 1, 0)
		old_ch[0] = ch
		if ch in xchrList:
			xcMs[0] = cM[0]
			print(name, "X", cM[0], sep=",", end="", file=OUT)
		else:
			print(name, ch, cM[0], sep=",", end="", file=OUT)	#make sure the name is written in both the geno file and the founder geno file.		
		
	#	j=0
		p0geno = fields[parent_indexes[0]].split(":")[0]
		p1geno = fields[parent_indexes[1]].split(":")[0]
		
		for  idx in range(9, len(fields)):
			fi = fields[idx]
	#		j+=1   #not sure what this does...?
			if fi==".":
				f=fi
			else:		
				f=fi.split(":")[0]
		
			if f==p0geno:
				g="A"
			elif f=="0/1":	
				g="H"
			elif f==p1geno:
				g="B"
			elif f=="./.":
				g="NA"	
			elif f==".":
				g="NA"	
			else:
				g="UN"
			if idx not in parent_indexes:	 				#these are the F2s (and maybe F1s if I forgot to take them out earlier...)	Don't bother printing the parents - we know their genotypes: AA for p0 and BB for p1
				print(",", g, end="", sep="", file=OUT)
					
		print("\n", end="", sep="", file=OUT)				
		return(print_info)






ch_dict={}
print_info = True
with open(vcf, 'r') as VCF:
	for line in VCF:
		if not line.startswith("##"):
			if line.startswith("#"):
				cols=format_Header(line)		#first make the header row for the table: do this from the id's of the vcf file: ID,,ind1,ind2,ind3... followed by a row of sexes: sex,,1,1,1,2,2,1,1...
				format_Pheno(cols)			#second make the phenotype row for the table: do this from the Pheno dataset or lacking that, make a row of for each phenotype in the input file: pheno1,,x,y,z...	
			else:
				#print(line.split("	")[0:2], end=" ", file=sys.stderr)
				print_info = format_Geno(line, print_info)		#third convert the vcf into genotypesA,H,B and NA: markerName,pos(1),A,B,H,NA,A...



#some stats:
print(len(ch_dict), "scaffolds have", sum(ch_dict.values()), "markers on them in the following distribution:", file=sys.stderr)
for ch in ch_dict:
	print(ch, "has", ch_dict[ch], "markers on it.", sep=" ", file=sys.stderr)
	
	
OUT.close()
PMAP.close() 
 
