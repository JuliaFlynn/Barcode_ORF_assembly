#!/usr/bin/perl
#

use List::Util qw[min max];

if ($#ARGV !=2) {
  print "usage: script.pl input_bc_file output_reject_file ham_cutoff\n" ;
  exit ;
}

$infile = $ARGV[0] ;
$outfile = $ARGV[1] ;
$cutoff = $ARGV[2] ;

$nbc=0 ;
open (INF, $infile) ;

$line = <INF> ;

@list_bc = () ;
@list_ham = () ;
@list_bc2 = () ;

while ($line = <INF>) {
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $cur_pos = $spline[0] ;
  $cur_cod = $spline[1] ;
  $cur_aa = $spline[2] ;
  $nspline = @spline ;
  for ($i=3; $i<$nspline; $i++) {
    push @list_bc, $spline[$i] ;
  }
}

close(INF) ;

$nlist_bc = @list_bc ;
for ($i=0; $i<$nlist_bc; $i++) {
  $list_ham[$i]=20 ;
  $list_bc2[$i]="A" ;
}

@list_remove = () ;

@list_int = (0..17) ;

for ($i=0; $i<$nlist_bc; $i++) {
  $cur_bc1 = $list_bc[$i] ;
  for ($count=0; $count<18; $count++) {
    $b1[$count] = substr($cur_bc1,$count,1) ;
  }
  for ($j=$i+1; $j<$nlist_bc; $j++) {
    $cur_bc2 = $list_bc[$j] ;
    $cur_ham = max($list_ham[$i],$list_ham[$j]) ;
    $test_ham = 0 ;
    for $k (@list_int) {
      $b2 = substr($cur_bc2,$k,1) ;
      if ($b1[$k] ne $b2) {$test_ham++} ;
      if ($test_ham>$cur_ham) {
        last ;
      }
    }
#    print "Atest $i $j $test_ham\n" ;
    if ($test_ham < $list_ham[$i]) {
#      print "$test_ham $cur_ham $cur_bc2 $cur_bc1\n" ;
      $list_bc2[$i] = $cur_bc2 ;
      $list_ham[$i] = $test_ham ;
    }
    if ($test_ham < $list_ham[$j]) {
#      print "test $i $test_ham $cur_ham $cur_bc2 $cur_bc1\n" ;
      $list_bc2[$j] = $cur_bc1 ;
      $list_ham[$j] = $test_ham ;
    }
  }
}

open(OUTF, ">$outfile") ;
for ($i=0; $i<$nlist_bc; $i++) {
  print "$list_bc[$i]\,$list_ham[$i]\,$list_bc2[$i]\n" ;
  if ($list_ham[$i] < $cutoff) {
    print OUTF "$list_bc[$i]\,$list_ham[$i]\,$list_bc2[$i]\n" ;
  }
}
