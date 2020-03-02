#!/usr/bin/perl
#
# inputs 01_bc_ORF_association file (ordered alphabetically by bc seq)
# outputs bc_ORF organized relative to mutations to parental sequence
# updated to organize by amino acid then codon

use POSIX ;
use List::Util qw(first) ;

sub aatoi {
  my($iaa,$aa) ;
  $aa = uc($_[0]) ;
  $iaa = 99 ;
  if ($aa eq "A") {$iaa = 0} ;
  if ($aa eq "C") {$iaa = 1} ;
  if ($aa eq "D") {$iaa = 2} ;
  if ($aa eq "E") {$iaa = 3} ;
  if ($aa eq "F") {$iaa = 4} ;
  if ($aa eq "G") {$iaa = 5} ;
  if ($aa eq "H") {$iaa = 6} ;
  if ($aa eq "I") {$iaa = 7} ;
  if ($aa eq "K") {$iaa = 8} ;
  if ($aa eq "L") {$iaa = 9} ;
  if ($aa eq "M") {$iaa = 10} ;
  if ($aa eq "N") {$iaa = 11} ;
  if ($aa eq "P") {$iaa = 12} ;
  if ($aa eq "Q") {$iaa = 13} ;
  if ($aa eq "R") {$iaa = 14} ;
  if ($aa eq "S") {$iaa = 15} ;
  if ($aa eq "T") {$iaa = 16} ;
  if ($aa eq "V") {$iaa = 17} ;
  if ($aa eq "W") {$iaa = 18} ;
  if ($aa eq "Y") {$iaa = 19} ;
  if ($aa eq "*") {$iaa = 20} ;
  $iaa ;
}

sub itoaa {
  my($i,$aa) ;
  $i = $_[0] ;
  $aa = "Z" ;
  if ($i == 0) {$aa = "A"} ;
  if ($i == 1) {$aa = "C"} ;
  if ($i == 2) {$aa = "D"} ;
  if ($i == 3) {$aa = "E"} ;
  if ($i == 4) {$aa = "F"} ;
  if ($i == 5) {$aa = "G"} ;
  if ($i == 6) {$aa = "H"} ;
  if ($i == 7) {$aa = "I"} ;
  if ($i == 8) {$aa = "K"} ;
  if ($i == 9) {$aa = "L"} ;
  if ($i == 10) {$aa = "M"} ;
  if ($i == 11) {$aa = "N"} ;
  if ($i == 12) {$aa = "P"} ;
  if ($i == 13) {$aa = "Q"} ;
  if ($i == 14) {$aa = "R"} ;
  if ($i == 15) {$aa = "S"} ;
  if ($i == 16) {$aa = "T"} ;
  if ($i == 17) {$aa = "V"} ;
  if ($i == 18) {$aa = "W"} ;
  if ($i == 19) {$aa = "Y"} ;
  if ($i == 20) {$aa = "*"} ;
  $aa ;
}

sub itocodon {
  my($i,$aa,$codon) ;
  $i = $_[0] ;
  if ($i == 1) {$codon = "aaa"}
  if ($i == 2) {$codon = "aac"}
  if ($i == 3) {$codon = "aag"}
  if ($i == 4) {$codon = "aat"}
  if ($i == 5) {$codon = "aca"}
  if ($i == 6) {$codon = "acc"}
  if ($i == 7) {$codon = "acg"}
  if ($i == 8) {$codon = "act"}
  if ($i == 9) {$codon = "aga"}
  if ($i == 10) {$codon = "agc"}
  if ($i == 11) {$codon = "agg"}
  if ($i == 12) {$codon = "agt"}
  if ($i == 13) {$codon = "ata"}
  if ($i == 14) {$codon = "atc"}
  if ($i == 15) {$codon = "atg"}
  if ($i == 16) {$codon = "att"}
  if ($i == 17) {$codon = "caa"}
  if ($i == 18) {$codon = "cac"}
  if ($i == 19) {$codon = "cag"}
  if ($i == 20) {$codon = "cat"}
  if ($i == 21) {$codon = "cca"}
  if ($i == 22) {$codon = "ccc"}
  if ($i == 23) {$codon = "ccg"}
  if ($i == 24) {$codon = "cct"}
  if ($i == 25) {$codon = "cga"}
  if ($i == 26) {$codon = "cgc"}
  if ($i == 27) {$codon = "cgg"}
  if ($i == 28) {$codon = "cgt"}
  if ($i == 29) {$codon = "cta"}
  if ($i == 30) {$codon = "ctc"}
  if ($i == 31) {$codon = "ctg"}
  if ($i == 32) {$codon = "ctt"}
  if ($i == 33) {$codon = "gaa"}
  if ($i == 34) {$codon = "gac"}
  if ($i == 35) {$codon = "gag"}
  if ($i == 36) {$codon = "gat"}
  if ($i == 37) {$codon = "gca"}
  if ($i == 38) {$codon = "gcc"}
  if ($i == 39) {$codon = "gcg"}
  if ($i == 40) {$codon = "gct"}
  if ($i == 41) {$codon = "gga"}
  if ($i == 42) {$codon = "ggc"}
  if ($i == 43) {$codon = "ggg"}
  if ($i == 44) {$codon = "ggt"}
  if ($i == 45) {$codon = "gta"}
  if ($i == 46) {$codon = "gtc"}
  if ($i == 47) {$codon = "gtg"}
  if ($i == 48) {$codon = "gtt"}
  if ($i == 49) {$codon = "taa"}
  if ($i == 50) {$codon = "tac"}
  if ($i == 51) {$codon = "tag"}
  if ($i == 52) {$codon = "tat"}
  if ($i == 53) {$codon = "tca"}
  if ($i == 54) {$codon = "tcc"}
  if ($i == 55) {$codon = "tcg"}
  if ($i == 56) {$codon = "tct"}
  if ($i == 57) {$codon = "tga"}
  if ($i == 58) {$codon = "tgc"}
  if ($i == 59) {$codon = "tgg"}
  if ($i == 60) {$codon = "tgt"}
  if ($i == 61) {$codon = "tta"}
  if ($i == 62) {$codon = "ttc"}
  if ($i == 63) {$codon = "ttg"}
  if ($i == 64) {$codon = "ttt"}
  $codon ;
}

sub codontoi {
  my($q,$aa,$codon) ;
  $codon = lc($_[0]) ;
  $codon = substr($codon,0,3) ;
  if ($codon eq "aaa") {$aa = 1} ;
  if ($codon eq "aac") {$aa = 2} ;
  if ($codon eq "aag") {$aa = 3} ;
  if ($codon eq "aat") {$aa = 4} ;
  if ($codon eq "aca") {$aa = 5} ;
  if ($codon eq "acc") {$aa = 6} ;
  if ($codon eq "acg") {$aa = 7} ;
  if ($codon eq "act") {$aa = 8} ;
  if ($codon eq "aga") {$aa = 9} ;
  if ($codon eq "agc") {$aa = 10} ;
  if ($codon eq "agg") {$aa = 11} ;
  if ($codon eq "agt") {$aa = 12} ;
  if ($codon eq "ata") {$aa = 13} ;
  if ($codon eq "atc") {$aa = 14} ;
  if ($codon eq "atg") {$aa = 15} ;
  if ($codon eq "att") {$aa = 16} ;
  if ($codon eq "caa") {$aa = 17} ;
  if ($codon eq "cac") {$aa = 18} ;
  if ($codon eq "cag") {$aa = 19} ;
  if ($codon eq "cat") {$aa = 20} ;
  if ($codon eq "cca") {$aa = 21} ;
  if ($codon eq "ccc") {$aa = 22} ;
  if ($codon eq "ccg") {$aa = 23} ;
  if ($codon eq "cct") {$aa = 24} ;
  if ($codon eq "cga") {$aa = 25} ;
  if ($codon eq "cgc") {$aa = 26} ;
  if ($codon eq "cgg") {$aa = 27} ;
  if ($codon eq "cgt") {$aa = 28} ;
  if ($codon eq "cta") {$aa = 29} ;
  if ($codon eq "ctc") {$aa = 30} ;
  if ($codon eq "ctg") {$aa = 31} ;
  if ($codon eq "ctt") {$aa = 32} ;
  if ($codon eq "gaa") {$aa = 33} ;
  if ($codon eq "gac") {$aa = 34} ;
  if ($codon eq "gag") {$aa = 35} ;
  if ($codon eq "gat") {$aa = 36} ;
  if ($codon eq "gca") {$aa = 37} ;
  if ($codon eq "gcc") {$aa = 38} ;
  if ($codon eq "gcg") {$aa = 39} ;
  if ($codon eq "gct") {$aa = 40} ;
  if ($codon eq "gga") {$aa = 41} ;
  if ($codon eq "ggc") {$aa = 42} ;
  if ($codon eq "ggg") {$aa = 43} ;
  if ($codon eq "ggt") {$aa = 44} ;
  if ($codon eq "gta") {$aa = 45} ;
  if ($codon eq "gtc") {$aa = 46} ;
  if ($codon eq "gtg") {$aa = 47} ;
  if ($codon eq "gtt") {$aa = 48} ;
  if ($codon eq "taa") {$aa = 49} ;
  if ($codon eq "tac") {$aa = 50} ;
  if ($codon eq "tag") {$aa = 51} ;
  if ($codon eq "tat") {$aa = 52} ;
  if ($codon eq "tca") {$aa = 53} ;
  if ($codon eq "tcc") {$aa = 54} ;
  if ($codon eq "tcg") {$aa = 55} ;
  if ($codon eq "tct") {$aa = 56} ;
  if ($codon eq "tga") {$aa = 57} ;
  if ($codon eq "tgc") {$aa = 58} ;
  if ($codon eq "tgg") {$aa = 59} ;
  if ($codon eq "tgt") {$aa = 60} ;
  if ($codon eq "tta") {$aa = 61} ;
  if ($codon eq "ttc") {$aa = 62} ;
  if ($codon eq "ttg") {$aa = 63} ;
  if ($codon eq "ttt") {$aa = 64} ;
  $aa ;
}

sub codontoaa {
  my($q,$aa,$codon) ;
  $codon = lc($_[0]) ;
  $codon = substr($codon,0,3) ;
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


if ($#ARGV !=2) {
  print "usage: script.pl description_file outputfile bc_obs_cutoff\n" ;
  exit ;
}
$dfile = $ARGV[0] ;
$outfile = $ARGV[1] ;
$cutoff = $ARGV[2] ;

open(DF, $dfile) ;
$line = <DF> ;
chomp($line) ;
@spline = split (/,/, $line) ;
$nf = $spline[1] ;
$npositions = 0 ;
for ($i=0; $i<$nf; $i++) {
  $line = <DF> ;
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $fname[$i]=$spline[0] ;

  $aa_start[$i] = $spline[1] ;
  $aa_end[$i] = $spline[2] ;
  @list_position = (@list_position,($aa_start[$i]..$aa_end[$i])) ;
  $wt_seq[$i] = $spline[3] ;
}

$npositions = @list_position ;
#zero arrays for storing bc info
for ($i=0; $i<=$npositions; $i++) {
  for ($j=0; $j<=64; $j++) {
    $bc_list[$i][$j] = "" ;
  }
}
$nflaggedbc = 0 ;

for my $i (0..$nf) {
  open(BCF, $fname[$i]) ;
  while ($line = <BCF>) {
    chomp($line) ;
    @spline = split (/,/, $line) ;
    $bc = $spline[0] ;
    $bc_reads = $spline[1] ;
    $naamut = $spline[2] ;
    $aa_desc = $spline[3] ;
    $ncodmut = $spline[4] ;
    $cod_desc = $spline[5] ;
    if ($bc_reads >= $cutoff) {
      if ($ncodmut == 0) {
# wildtype sequence store for output
        $bc_list[$npositions][0] = $bc_list[$npositions][0] . "," . $bc ;
      } elsif ($ncodmut == 1) {
# one codon differs from parental sequence, store for output
        $aa_position = $cod_desc ;
        $aa_position =~ s/[^0-9]//g ;
        if ($aa_position >= $aa_start[$i]) {
          if ($aa_position <= $aa_end[$i]) {
# mutation in expected range, store for output
            $l = length($cod_desc) ;
            $mutcod = substr($cod_desc,$l-3,3) ;
            if (index($mutcod, 'N') == -1) {
              if ($naamut == 1) {
                $jmut = &codontoi($mutcod) ;
                $imut = first { $list_position[$_] eq $aa_position } 0..$#list_position;
                $bc_list[$imut][$jmut] = $bc_list[$imut][$jmut] . "," . $bc ;
              } elsif ($naamut == 0) {
                $mutaa = &codontoaa($mutcod) ;
                $imut = first { $list_position[$_] eq $aa_position } 0..$#list_position;
                $jmut = &codontoi($mutcod) ;
                $bc_list[$imut][$jmut] = $bc_list[$imut][$jmut] . "," . $bc ;
              }
            }
          }
        }
      }
    } 
  }
  close(BCF) ;
}

open(OUTF, ">$outfile") ;
print OUTF "Position, aa, codon, bc_list\n" ;
print OUTF "WILDTYPE, WILDTYPE, WILDTYPE" ;
print OUTF "$bc_list[$npositions][0]\n" ;
for ($i=0; $i<$npositions; $i++) {
  $cur_pos = $list_position[$i] ;
  @aa_list = qw(* W Y F L M I V C A G P S T N Q H R K D E ) ;
  for $aa (@aa_list) {
    for ($j=1; $j<=64; $j++) {
      $cur_cod = &itocodon($j) ;
      $cur_cod = uc($cur_cod) ;
      $cur_aa = &codontoaa($cur_cod) ;
      if ($cur_aa eq $aa) {
        print OUTF "$cur_pos\,$cur_aa\,$cur_cod" ;
        print OUTF "$bc_list[$i][$j]\n" ;
      }
    }
  }
}
close(OUTF) ;
