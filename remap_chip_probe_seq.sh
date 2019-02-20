#!/bin/bash

FILE=HillsFeline_340K_20024177_A1.csv
OUTNAME=340KHILLS

INPATH=/home/buckley/Desktop/hills_remap/probe_seq
TMP=/home/buckley/Desktop/hills_remap/tmp_result
OUTPATH=/home/buckley/Desktop/hills_remap/out

DBIN=/home/buckley/Documents/ref_genomes/fasta/mammals/cat/Felis_catus_9.0/Felis_catus_9.0.fa

Rscript /home/buckley/Desktop/hills_remap/scripts/csv_2_fasta.R $INPATH/$FILE $TMP/$OUTNAME

makeblastdb -dbtype nucl -in $DBIN -out $TMP/blast.db

(
blastn -query $TMP/${OUTNAME}_Aseq.fa \
	-db $TMP/blast.db \
	-out $OUTPATH/$OUTNAME.blast_A_out.txt \
	-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand btop" \
	-culling_limit 2 \
	-num_threads 16
	) &

(
blastn -query $TMP/${OUTNAME}_Bseq.fa \
	-db $TMP/blast.db \
	-out $OUTPATH/$OUTNAME.blast_B_out.txt \
	-outfmt "6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore sstrand btop" \
	-culling_limit 2 \
	-num_threads 16
	) &

wait

Rscript /home/buckley/Desktop/hills_remap/scripts/filter_mapped_probes.R $INPATH/$FILE $OUTPATH/$OUTNAME.blast_A_out.txt $OUTPATH/$OUTNAME.blast_B_out.txt $OUTPATH/$OUTNAME.map.tsv 