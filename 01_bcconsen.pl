#!/usr/bin/perl
#
# read in sorted sequence file (sequences only)
# count number of each unique barcodes in R1
# for each unique barcode output:
#   barcode sequence, N aa substitutions, first aa sub, first non wt codon
#
# updated to start considering r2 at end of primer
#

use POSIX ;

sub translate {
  my($i,$n_nuke,$n_aa,$aa_seq,$cur_codon,$nuke_seq,$cur_aa,$start) ;
  $nuke_seq = lc($_[0]) ;
  $n_nuke = length($nuke_seq) ;
  $n_aa = floor($n_nuke/3) ;
  $aa_seq = "" ;
  for ($i=0; $i<$n_aa; $i++) {
    $start = $i*3 ;
    $cur_codon = substr($nuke_seq,$start,3) ;
    $cur_aa = &codontoaa($cur_codon) ;
    $aa_seq = $aa_seq . $cur_aa ;
  }
  $aa_seq ;
}

sub codontoaa {
  my($q,$aa,$codon) ;
  $codon = lc($_[0]) ;
  $codon = substr($codon,0,3) ;
  $aa = "X" ;
  if ($codon eq "aaa") {$aa = "K"} ;
  if ($codon eq "aac") {$aa = "N"} ;
  if ($codon eq "aag") {$aa = "K"} ;
  if ($codon eq "aat") {$aa = "N"} ;
  if ($codon eq "aca") {$aa = "T"} ;
  if ($codon eq "acc") {$aa = "T"} ;
  if ($codon eq "acg") {$aa = "T"} ;
  if ($codon eq "act") {$aa = "T"} ;
  if ($codon eq "aga") {$aa = "R"} ;
  if ($codon eq "agc") {$aa = "S"} ;
  if ($codon eq "agg") {$aa = "R"} ;
  if ($codon eq "agt") {$aa = "S"} ;
  if ($codon eq "ata") {$aa = "I"} ;
  if ($codon eq "atc") {$aa = "I"} ;
  if ($codon eq "atg") {$aa = "M"} ;
  if ($codon eq "att") {$aa = "I"} ;
  if ($codon eq "caa") {$aa = "Q"} ;
  if ($codon eq "cac") {$aa = "H"} ;
  if ($codon eq "cag") {$aa = "Q"} ;
  if ($codon eq "cat") {$aa = "H"} ;
  if ($codon eq "cca") {$aa = "P"} ;
  if ($codon eq "ccc") {$aa = "P"} ;
  if ($codon eq "ccg") {$aa = "P"} ;
  if ($codon eq "cct") {$aa = "P"} ;
  if ($codon eq "cga") {$aa = "R"} ;
  if ($codon eq "cgc") {$aa = "R"} ;
  if ($codon eq "cgg") {$aa = "R"} ;
  if ($codon eq "cgt") {$aa = "R"} ;
  if ($codon eq "cta") {$aa = "L"} ;
  if ($codon eq "ctc") {$aa = "L"} ;
  if ($codon eq "ctg") {$aa = "L"} ;
  if ($codon eq "ctt") {$aa = "L"} ;
  if ($codon eq "gaa") {$aa = "E"} ;
  if ($codon eq "gac") {$aa = "D"} ;
  if ($codon eq "gag") {$aa = "E"} ;
  if ($codon eq "gat") {$aa = "D"} ;
  if ($codon eq "gca") {$aa = "A"} ;
  if ($codon eq "gcc") {$aa = "A"} ;
  if ($codon eq "gcg") {$aa = "A"} ;
  if ($codon eq "gct") {$aa = "A"} ;
  if ($codon eq "gga") {$aa = "G"} ;
  if ($codon eq "ggc") {$aa = "G"} ;
  if ($codon eq "ggg") {$aa = "G"} ;
  if ($codon eq "ggt") {$aa = "G"} ;
  if ($codon eq "gta") {$aa = "V"} ;
  if ($codon eq "gtc") {$aa = "V"} ;
  if ($codon eq "gtg") {$aa = "V"} ;
  if ($codon eq "gtt") {$aa = "V"} ;
  if ($codon eq "taa") {$aa = "*"} ;
  if ($codon eq "tac") {$aa = "Y"} ;
  if ($codon eq "tag") {$aa = "*"} ;
  if ($codon eq "tat") {$aa = "Y"} ;
  if ($codon eq "tca") {$aa = "S"} ;
  if ($codon eq "tcc") {$aa = "S"} ;
  if ($codon eq "tcg") {$aa = "S"} ;
  if ($codon eq "tct") {$aa = "S"} ;
  if ($codon eq "tga") {$aa = "*"} ;
  if ($codon eq "tgc") {$aa = "C"} ;
  if ($codon eq "tgg") {$aa = "W"} ;
  if ($codon eq "tgt") {$aa = "C"} ;
  if ($codon eq "tta") {$aa = "L"} ;
  if ($codon eq "ttc") {$aa = "F"} ;
  if ($codon eq "ttg") {$aa = "L"} ;
  if ($codon eq "ttt") {$aa = "F"} ;
  $aa ;
}


$len_index = 18 ;

if ($#ARGV !=2) {
  print "usage: script.pl input_file bc_ORF_outfile obs_cut\n" ;
  exit ;
}

$indesc = $ARGV[0] ;
$outfile = $ARGV[1] ;
$obs_cut = $ARGV[2] ;

open(IND, $indesc) ;

$line = <IND> ;
chomp($line) ;
@spline = split(/,/, $line) ;
print "test: $spline[0]\n" ;
if (index($spline[0],"sorted_seq_file") < 0 ) {
  print "Error input file expect sorted_seq_file at line 1\n" ;
  exit
}
$infile = $spline[1] ;

$line = <IND> ;
chomp($line) ;
@spline = split(/,/, $line) ;
if ($spline[0] ne "R2_anneal_primer_seq" ) {
  print "Error input file expect R2_anneal_primer_seq at line 2\n" ;
  exit
}
$primer_seq = uc($spline[1]) ;
$primer_len = length($primer_seq) ;

$line = <IND> ;
chomp($line) ;
@spline = split(/,/, $line) ;
if ($spline[0] ne "R2_length_to_analyze" ) {
  print "Error input file expect R2_length_to_analyze at line 3\n" ;
  exit
}
$r2_len = $spline[1] ;

$line = <IND> ;
chomp($line) ;
@spline = split(/,/, $line) ;
if ($spline[0] ne "rand_aastart-aaend" ) {
  print "Error input file expect rand_aastart-aaend at line 4\n" ;
  exit
}
$rand_aastart = $spline[1] ;
$rand_aaend = $spline[2] ;

$line = <IND> ;
chomp($line) ;
@spline = split(/,/, $line) ;
if ($spline[0] ne "nuke_start-wt_orf" ) {
  print "Error input file expect nuke_start-wt_orf at line 5\n" ;
  exit
}
$wt_start = $spline[1] ;
$wt_orf = uc($spline[2]) ;

#print "wt orf: $wt_orf\n" ;

$primer_start = index($wt_orf,$primer_seq) ;
$primer_end = $primer_start + $primer_len ;
$diff = $primer_end - $wt_start ;
$remainder = $diff % 3 ;
$primer_shift = 0 ;
if ($remainder !=0) {$primer_shift = 3 - $remainder} ;
$primer_end_aa = ($diff+4)/3 ;
print "primer end, next aa: $primer_end_aa\n" ;

$r2_end_len = $r2_len - $primer_len ;

print "primer end: $primer_end\n" ;
if ($primer_start < 0) {
  print "Error primer $primer_seq not in wt ORF\n" ;
  exit ;
}

$wt_r2 = substr($wt_orf,$primer_end,$r2_end_len) ;

print "wt r2: $wt_r2\n" ;

$wt_coding_start = $primer_end + $primer_shift ;
$r2_coding_start = $primer_len + $primer_shift ;
$r2_coding_len = $r2_end_len - $primer_shift ;
$wt_r2_coding = substr($wt_orf,$wt_coding_start,$r2_coding_len) ;
print "wt_r2_coding: $wt_r2_coding\n" ;
$aa_r2_coding = &translate($wt_r2_coding) ;
$n_aa = length($aa_r2_coding) ;
print "translation of r2: $aa_r2_coding\n" ;

close(IND) ;

open(INF, $infile) ;
open(OUTF, ">$outfile") ;

$old_bc = "ZZZ" ;
$count = 0 ;
$max = 0 ;
while ($line = <INF>) {
  chomp($line) ;
  @spline = split(/,/, $line) ;
  $cur_bc = substr($spline[0],0,18) ;
  $r2_tran_len =$r2_end_len-$primer_shift ;
  $cur_orf = substr($spline[1],$r2_coding_start,$r2_coding_len) ;
  if ($cur_bc eq $old_bc) {
    $orf[$count] = $cur_orf ;
    $count++;
  } else {
    if ($count>=$obs_cut) {
#   barcode sequence, N aa substitutions, first aa sub, first non wt codon
      print OUTF "$old_bc,$count," ;
      for ($i=0; $i<$r2_coding_len; $i++) {
        $nA[$i] = 0 ;
        $nC[$i] = 0 ;
        $nG[$i] = 0 ;
        $nT[$i] = 0 ;
        $nN[$i] = 0 ;
      }
      for ($i=0; $i<$count; $i++) {
        for ($j=0; $j<$r2_coding_len; $j++) {
          $nt = substr($orf[$i],$j,1) ;
          if ($nt eq "A") {$nA[$j]++} ;
          if ($nt eq "C") {$nC[$j]++} ;
          if ($nt eq "G") {$nG[$j]++} ;
          if ($nt eq "T") {$nT[$j]++} ;
          if ($nt eq "N") {$nN[$j]++} ;
        }
      }
# All reads for this index are input; analyze and output consensus
      $consensus = "" ;
      for ($j=0; $j<$r2_coding_len; $j++) {
        $nuke = "N" ;
        $tot = $nA[$j] + $nC[$j] + $nG[$j] + $nT[$j] ;
        $cut = floor($tot/2) ;
        if ($cut < 2) {$cut = 3} ;
        $nuke = "N" ;
        if ($nA[$j] >= $cut) {$nuke = "A"} ;
        if ($nC[$j] >= $cut) {$nuke = "C"} ;
        if ($nG[$j] >= $cut) {$nuke = "G"} ;
        if ($nT[$j] >= $cut) {$nuke = "T"} ;
        $consensus = $consensus . $nuke ;
      }
#      print "con: $consensus\n" ;
#      print "wt : $wt_r2_coding\n" ;
      $ndiffcodons = 0 ;
      $ndiffaas = 0 ;
      $aamut = "" ;
      $codon_mut = "" ;
      for ($i=0; $i<$n_aa; $i++) {
        $start = $i*3 ;
        $cur_codon = substr($consensus,$start,3) ;
        $wt_codon = substr($wt_r2_coding,$start,3) ;
        $cur_aa = &codontoaa($cur_codon) ;
        $wt_aa = &codontoaa($wt_codon) ;
        if ($wt_aa ne $cur_aa) {
          $ndiffaas++ ;
          $pos = $i + $primer_end_aa ;
          $temp = $wt_aa . $pos . $cur_aa ;
          $aamut = $aamut . $temp ;
        }
        if ($wt_codon ne $cur_codon) {
          $ndiffcodons++ ;
          $pos = $i + $primer_end_aa ;
          $temp = $wt_codon . $pos . $cur_codon ;
          $codon_mut = $codon_mut . $temp ;
        }
      }
      print OUTF "$ndiffaas\,$aamut\,$ndiffcodons\,$codon_mut\n"
    }
    $count = 0 ;
    $orf[$count] = $cur_orf ;
    $count++;
    $old_bc = $cur_bc ;
  }
}
     
close(INF);
close(OUTF);
