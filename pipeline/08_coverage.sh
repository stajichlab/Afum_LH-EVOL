#!/usr/bin/bash
#SBATCH -p short -N 1 -n 2

module unload miniconda2
module load miniconda3
source activate mosdepth

if [ ! -f genome/chroms.bed ]; then
	awk 'BEGIN{OFS="\t"} {print $1,1,$2}' genome/FungiDB-39_AfumigatusAf293_Genome.*.fai > genome/chroms.bed
fi

for file in aln/*.cram; 
do 
	b=$(basename $file .cram)
	mosdepth -b genome/chroms.bed -t 2 -f genome/FungiDB-39_AfumigatusAf293_Genome.fasta depth/$b $file
done

