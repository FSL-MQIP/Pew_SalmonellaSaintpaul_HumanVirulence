#!/bin/bash

#The not-so fancy bash script for blast
for db in *.fas
do
makeblastdb -in $db -dbtype 'nucl'
for genome in *.fasta
do
blastn -query $genome -db $db -evalue 1e-10 -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -out "$genome"_"$db".txt
done
done