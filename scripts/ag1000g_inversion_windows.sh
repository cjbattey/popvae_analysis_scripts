#!/bin/bash
#run locator in windows across the genome

cd /home/cbattey2/popvae/

step=200000
window_size=200000
start=12000000
stop=50000000
#get chromosome length
header=`tabix -H /home/data_share/ag1000g/phase2/ag1000g.phase2.ar1.pass.biallelic.2L.vcf.gz | grep "##contig=<ID=2L,length="`
length=`echo $header | awk '{sub(/.*=/,"");sub(/>/,"");print}'`

#subset vcf by region and run locator
endwindow=$((start+window_size))
for startwindow in `seq $start $step $stop`
do
	echo "processing $startwindow to $endwindow"
	tabix -h /home/data_share/ag1000g/phase2/ag1000g.phase2.ar1.pass.biallelic.2L.vcf.gz \
	2L:$startwindow\-$endwindow > data/ag1000g/tmp.vcf

	python scripts/popvae.py \
	--infile data/ag1000g/tmp.vcf \
	--out out/ag1000g/inversion_windows/$startwindow\-$endwindow \
	--PCA \
	--seed 666 \
	--patience 50

	endwindow=$((endwindow+step))
	rm data/ag1000g/tmp.vcf
done
