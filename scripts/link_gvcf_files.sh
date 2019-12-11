#!/usr/bin/bash
mkdir -p gvcf
cd gvcf
for strain in $(cat ../strains.txt); do ln -s /bigdata/stajichlab/shared/projects/Afum_popgenome/variantcall/Variants/$strain.* .; done
cd ..
