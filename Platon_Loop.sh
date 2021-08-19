#!/bin/bash

#The not-so fancy bash script for blast
for genome in *.fasta
do
platon --db ./db --mode accuracy --output "$genome"_output/ --threads 5 $genome
done