# Bolonlab_barcode_ORF_assembly
Association of ORF mutations with barcodes from paired-end sequencing data

00_PE_splice.pl - takes PE fastq files as input, outputs file with 
               joined sequences ("read1,read2") minimal quality check (bases
               with PHRED < 10 noted as N)

01_bcconsen.pl - takes sorted bc file (bc,r2 format) and outputs consensus
               r2 sequence for each barcode as well as mutations relative
               to a parental sequence descriped in input file.

01b_bccheck.pl - takes 01 file and outputs number of bc for each possible point 
               mutation. uses same input description file as 01.pl

02_bcclean.pl - takes 01 net assembly file and reorder's based on ORF variants
               removes ORF variants from list that were unplanned/unintended

03_bclistident.pl - takes cleaned 02 bc file and outputs duplicate bcs as 
                    std out > reject file

03b_bchamcheck.pl - NOTE: only suggested for files with <40,000 bc
                  takes some time to run on larger files
                  takes cleaned02 bc file and calculates minimum hamming
                  distance between all barcodes. Stdout is list of all bc
                  and min hamming distances. Also makes reject_file with
                  with barcodes that differ by less than cutoff (input param)

04_bcprune.pl - takes cleaned02 bc file and reject03 file, removes barcodes
               in reject list and outputs cleaned04 bc file. Stdout is list
               of rejected barcodes and the position and ORF variant. Cutoff
               parameter is in input (min ham). The cleaned04 bc file is
               good to use for analyzing frequencies from SE reads of bcs.
