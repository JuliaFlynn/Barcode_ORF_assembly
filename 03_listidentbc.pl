#!/usr/bin/perl
#

use POSIX ;
use List::Util qw(first) ;

if ($#ARGV !=0) {
  print "usage: script.pl bc02file\n" ;
  exit ;
}

$infile = $ARGV[0] ;

@wt_list = () ;
@bc_list = () ;

open(INF, $infile) ;
$line = <INF> ;

while ($line = <INF>) {
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $n_curlist = @spline ;
  for ($j=3; $j<$n_curlist; $j++) {
    push @bc_list, $spline[$j] ;
  }
}

close(INF) ;

my %bc_list_counts;
for (@bc_list) {$bc_list_counts{$_}++ } ;
@dupes = grep { $bc_list_counts{$_} > 1 } keys %bc_list_counts;
$ndupes = @dupes ;
for ($i=0; $i<$ndupes; $i++) {
  print "$dupes[$i]\,0\n" ;
}
