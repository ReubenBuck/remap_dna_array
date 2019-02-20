#!/bin/bash

# other potential flags
# header lines to skip
# 


FILE=HillsFeline_250K_iSelect_scoreResults_designBumpv3.csv
OUTNAME=250k

INPATH=/home/buckley/Desktop/hills_remap/Cat_Array_SNP_Design
TMP=/home/buckley/Desktop/hills_remap/tmp_result
OUTPATH=/home/buckley/Desktop/hills_remap/out

DBIN=/home/buckley/Documents/ref_genomes/fasta/mammals/cat/Felis_catus_9.0/Felis_catus_9.0.fa

tail -n +2 $INPATH/$FILE | cut -f2 -d"," | cut -f1 -d"[" > $TMP/${OUTNAME}_front_seq.txt

tail -n +2 $INPATH/$FILE | cut -f2 -d"," | cut -f2 -d"]" > $TMP/${OUTNAME}_back_seq.txt

tail -n +2 $INPATH/$FILE | cut -f15 -d"," > $TMP/probe.name.txt



touch $TMP/query.fa
LNS=$(wc -l $TMP/probe.name.txt | cut -f1 -d" ")

for i in $(seq 1 $LNS) ; do

	NAME=$(sed "${i}q;d" $TMP/probe.name.txt)
	FRONT=$(sed "${i}q;d" $TMP/${OUTNAME}_front_seq.txt)
	BACK=$(sed "${i}q;d" $TMP/${OUTNAME}_back_seq.txt)
	
	echo ">${NAME}_front" >> $TMP/query.fa
	echo $FRONT >> $TMP/query.fa

	echo ">${NAME}_back" >> $TMP/query.fa
	echo $BACK >> $TMP/query.fa
done


# next is blast step

# build database
makeblastdb -dbtype nucl -in $DBIN -out $TMP/blast.db

blastn -query $TMP/q_test.fa \
	-db $TMP/blast.db \
	-out $OUTPATH/$OUTNAME.blast_out.txt \
	-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore btop" \
	-culling_limit 2 \
	-num_threads 16


