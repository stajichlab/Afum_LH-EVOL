#!/usr/bin/bash

#SBATCH --mem=24gb --ntasks 24 --nodes 1
#SBATCH --time=2:00:00 -p short
#SBATCH -J maketree --out logs/make_tree.log


CPU=2
if [ $SLURM_CPUS_ON_NODE ]; then
 CPU=$SLURM_CPUS_ON_NODE
fi

if [[ -f config.txt ]]; then
    source config.txt
else
    echo "Need a config.txt"
    exit
fi

if [[ -z $REFNAME ]]; then
    REFNAME=REF
fi
module load parallel
module load perl/5.20.2
module load bcftools/1.9
module load samtools/1.9
module load IQ-TREE
module load fasttree

print_fas() {
    printf ">%s\n%s\n" $1 $(bcftools query -e 'INFO/AF < 0.1' -s $1 -f '[%TGT]' $2)
}

export -f print_fas

for TYPE in SNP INDEL
do
    root=$FINALVCF/$PREFIX.selected.$TYPE
    FAS=$TREEDIR/$PREFIX.$TYPE.mfa

    if [ -f $root.vcf ]; then
	bgzip $root.vcf
	tabix $root.vcf.gz
    fi
    vcf=$root.vcf.gz
    # for the clade-specific runs don't include the Af293 ref seq
    #printf ">%s\n%s\n" $REFNAME $(bcftools query -e 'INFO/AF < 0.1' -f '%REF' ${vcf}) > $FAS
    # note that if you change back to use the Af293 ref then you needto change this next line from ">" to ">>" to
    # append
    parallel -j $CPU print_fas ::: $(bcftools query -l ${vcf}) ::: $vcf > $FAS

    perl -ip -e 'if(/^>/){s/[\(\)#]/_/g; s/_+/_/g } else {s/[\*\.]/-/g }' $FAS

    tab=$root.bcftools.tab

    if [ ! -f $tab ]; then
	bcftools query -H -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' ${vcf} > $tab 
	bcftools query -H -e 'INFO/AF < 0.1' -f '%CHROM\t%POS\t%INFO/AF\t%INFO/AC\t%INFO/AN\t%REF\t%ALT{0}[\t%TGT]\n' ${vcf} > $root.highfreq.bcftools.tab 
    fi

    if [[ $TYPE == "SNP" && ! -f $TREEDIR/$PREFIX.$TYPE.fasttree.tre ]]; then
	FastTreeMP -gtr -gamma -nt < $FAS > $TREEDIR/$PREFIX.$TYPE.fasttree.tre
    fi
    
done
