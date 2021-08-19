#!usr/bin/perl
# This is a script to go over several files with SRA IDs and then create assemblies using SKESA

# Get all files with SRA Ids from directory
my @txt = glob "*srr*.txt";

# Read each file at a time
foreach $filename (@txt){

# Open file using file handler
open INPUT, "<", $filename;

# Check each line at a time                         
while(defined($line=<INPUT>)){

chomp($line);

# Run fastq-dump
$output= "" . $line . "_skesa_contigs.fasta";

open FASTQ, "|/programs/skesa.centos6.9/skesa --sra_run $line --cores 5 --memory 32 --use_paired_ends --contigs_out $output";
close FASTQ;
}
close INPUT;
}
