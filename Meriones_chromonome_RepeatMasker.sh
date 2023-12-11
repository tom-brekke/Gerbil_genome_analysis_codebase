#!/bin/bash --login

##job resources
#SBATCH --ntasks=40
#SBATCH --time=0-72:00
#SBATCH --mem-per-cpu=8000
#SBATCH --job-name=RepMask
#SBATCH --account=scw1055
#SBATCH --output=RM.out.%J
#SBATCH --error=RM.err.%J
#SBATCH --partition=compute

module purge
module load RepeatMasker/4.1.0

#genome="Meriones_chromonome_v5.Fasta"
outdir="/scratch/b.bss81d/Assemblies/Meriones_Chromonome/v5.1_fasta_entries/"

cd $outdir

for genome in `ls Chr*.fasta`
do
	run_cmd.sh "RepeatMasker -pa 40 -species rodentia -gff -dir .  $genome" $genome.masked
done
