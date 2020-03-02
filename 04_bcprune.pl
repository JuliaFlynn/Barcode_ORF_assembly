#!/usr/bin/perl
#
# input 02_bclist, and reject_list files
# remove barcodes in reject list
# output rejected barcodes that link to same ORF - may want to 
#   update reject list to remove most abundant of these

if ($#ARGV !=3) {
  print "usage: script.pl input_02_file reject_03_file output_04_file ham_cut\n" ;
  exit ;
}

$infile = $ARGV[0] ;
$rejectfile = $ARGV[1] ;
$outfile = $ARGV[2] ;
$cutoff = $ARGV[3] ;

@list_reject = () ;

open(RF, $rejectfile) ;
while ($line = <RF>) {
  chomp($line) ;
  @spline = split (/,/, $line) ;
  if ($spline[1]<=$cutoff) {
    push @list_reject, $spline[0] ;
# note no need to push $spline[2] - it will be $spline[1] another time
#    push @list_reject, $spline[2] ;
  }
}
close(RF) ;
$nreject = @list_reject ;

print "n reject $nreject\n" ;
#for ($i=0; $i<@list_reject; $i++) {
#  print "$i $list_reject[$i]\n" ;
#}

open(INF, $infile) ;
open(OUTF, ">$outfile") ;

$line = <INF> ;
chomp($line) ;
print OUTF "$line\n" ;

while ($line = <INF>) {
  chomp($line) ;
  @spline = split (/,/, $line) ;
  $cur_pos = $spline[0] ;
  $cur_cod = $spline[1] ;
  $cur_aa = $spline[2] ;
  $nspline = @spline ;
  @list_bc = () ;
  for ($i=3; $i<$nspline; $i++) {
    if ($spline[$i] ~~ @list_reject) {
      print "ejected $spline[$i]\, $cur_pos\, $cur_cod\n" ;
    } else {
      push @list_bc, $spline[$i] ;
    }
  }
  print OUTF "$cur_pos\,$cur_cod\,$cur_aa" ;
  for ($i=0; $i<@list_bc; $i++) {
    print OUTF ",$list_bc[$i]" ;
  }
  print OUTF "\n" ;
}

close(INF) ;
close(OUTF) ;

