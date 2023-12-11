#!/usr/bin/env Python3

#written by TDB on 18-10-2017
#extracts chr lengths from fasta (1 line)
from Bio import SeqIO
import sys
import re
fasta=sys.argv[1]
FASTA=open(fasta, 'r')
print("Chr,Pos", file=sys.stdout)
for record in SeqIO.parse(FASTA, 'fasta'):
	print(record.name, len(record.seq),sep=",", file=sys.stdout)

FASTA.close()
