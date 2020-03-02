#!/usr/bin/perl
#
# read in paired end fastq files
# check that R1 and R2 are from same cluster
# output reads with abbreviated quality indicator encrypted in sequences
# updated to output without encoding (qscore>10 as is, <10 as N)
#

sub nt_qual_encode {
  my($qval, $letin, $letout);
  $letin = $_[0] ;
  $qval = $_[1] ;
  $letout = "z" ;
  if ($qval>9) {
    if ($letin eq "A") {$letout = "A"} ;
    if ($letin eq "C") {$letout = "C"} ;
    if ($letin eq "G") {$letout = "G"} ;
    if ($letin eq "T") {$letout = "T"} ;
    if ($letin eq "N") {$letout = "N"} ;
  } else {
    $letout = "N" ;
  }
  $letout ;
}

#sub nt_qual_encode {
#  my($qval, $letin, $letout);
#  $letin = $_[0] ;
#  $qval = $_[1] ;
#  $letout = "z" ;
#  if ($qval>29) {
#    $letout =  uc($letin) ;
#  } elsif ($qval>19) {
#    $letout = lc($letin) ;
#  } elsif ($qval>14) {
#    if ($letin eq "A") {$letout = "B"} ;
#    if ($letin eq "C") {$letout = "D"} ;
#    if ($letin eq "G") {$letout = "H"} ;
#    if ($letin eq "T") {$letout = "U"} ;
#  } elsif ($qval>9) {
#    if ($letin eq "A") {$letout = "b"} ;
#    if ($letin eq "C") {$letout = "d"} ;
#    if ($letin eq "G") {$letout = "h"} ;
#    if ($letin eq "T") {$letout = "u"} ;
#  } elsif ($qval>4) {
#    if ($letin eq "A") {$letout = "1"} ;
#    if ($letin eq "C") {$letout = "2"} ;
#    if ($letin eq "G") {$letout = "3"} ;
#    if ($letin eq "T") {$letout = "4"} ;
#  } else {
#    $letout = "n" ;
#  }
#  $letout ;
#}

if ($#ARGV != 5) {
  print "usage: script.pl READ1.fastq READ2.fastq seq_outfile sum_outfile read1_length read2_length\n";
  exit;
}

$R1_file = $ARGV[0];
$R2_file = $ARGV[1];
$outfile = $ARGV[2];
$sumfile = $ARGV[3];
$R1_len = $ARGV[4];
$R2_len = $ARGV[5];

print "R1 file: $R1_file\n" ;
print "R2 file: $R2_file\n" ;

open(R1F, $R1_file) ;
open(R2F, $R2_file) ;
open(SUMF, ">$sumfile") ;
open(OUTF, ">$outfile") ;

# zero array to store quality information for each read
for ($i=0; $i<$R1_len; $i++) {
  for ($j=0; $j<41; $j++) {
    $R1_qsum[$i][$j] = 0 ;
  }
}
for ($i=0; $i<$R2_len; $i++) {
  for ($j=0; $j<41; $j++) {
    $R2_qsum[$i][$j] = 0 ;
  }
}

$countlines = 0;
$check = 0 ;
# Read in fastq files in synchrony 
while ($line1 = <R1F>) {
  $line2 = <R2F> ;
  chomp($line1) ;
  chomp($line2) ;
  $char1 = substr($line1,0,1);
  $char2 = substr($line2,0,1);
  $countlines++ ;
  $check++ ;
  if ($check == 1) {
    @spline1 = split (/ /, $line1) ;
    $R1_label = $spline1[0] ;
    @spline2 = split (/ /, $line2) ;
    $R2_label = $spline2[0] ;
    if ($R1_label ne $R2_label) {
      print "Read1 and Read2 labels differ\n" ;
      print "  $R1_label\n" ;
      print "  $R2_label\n" ;
      exit ;
    }
    if ($char1 ne "@") {
      print "A error at line $countlines\n";
      exit ;
    }
  }
  if ($check == 2) {
    $R1_seq = $line1 ;
    $R2_seq = $line2 ;
    $R1_curlen = length($R1_seq) ;
    $R2_curlen = length($R2_seq) ;
    if ($R1_curlen != $R1_len) {
      print "Read 1 input length $R1_len does NOT match actual length $R1_curlen\n" ;
      exit;
    }
    if ($R2_curlen != $R2_len) {
      print "Read 2 input length $R2_len does NOT match actual length $R2_curlen\n" ;
      exit;
    }
  }
  if ($check == 3) {
    if ($char1 ne "+") {
      print "B error at line $countlines\n";
      exit ;
    }
  }
  if ($check == 4) {
    $check = 0 ;
    $R1_qual = $line1 ;
    $R2_qual = $line2 ;
    # encode quality for each read into sequence and output
    $R1_40 = substr($R1_seq, 0, 40) ;
    $R1_eseq = "" ;
    for ($i=0; $i<$R1_len; $i++) {
      $curnt = substr($R1_seq,$i,1) ;
      $let = substr($R1_qual,$i,1) ;
      $curqual = ord($let) - 33 ;
      $curent = &nt_qual_encode($curnt, $curqual) ;
      $R1_eseq = $R1_eseq . $curent ;
      $R1_qsum[$i][$curqual]++ ;
    }
    $R2_eseq = "" ;
    for ($i=0; $i<$R2_len; $i++) {
      $curnt = substr($R2_seq,$i,1) ;
      $let = substr($R2_qual,$i,1) ;
      $curqual = ord($let) - 33 ;
      $curent = &nt_qual_encode($curnt, $curqual) ;
      $R2_eseq = $R2_eseq . $curent ;
      $R2_qsum[$i][$curqual]++ ;
    }
    print OUTF "$R1_eseq\," ;
    print OUTF "$R2_eseq\n" ;
  }
}

print SUMF "READ1 Quality summary\n" ;
for ($j=0; $j<41; $j++) {
  print SUMF "quality $j\," ;
  for ($i=0; $i<$R1_len; $i++) {
    print SUMF "$R1_qsum[$i][$j]\," ;
  }
  print SUMF "\n" ;
}
print SUMF "\n" ;
print SUMF "READ2 Quality summary\n" ;
for ($j=0; $j<41; $j++) {
  print SUMF "quality $j\," ;
  for ($i=0; $i<$R2_len; $i++) {
    print SUMF "$R2_qsum[$i][$j]\," ;
  }
  print SUMF "\n" ;
}
print SUMF "\n" ;

close(R1F) ;
close(R2F) ;
close(SUMF) ;
close(OUTF) ;

