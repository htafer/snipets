#!/usr/bin/perl

# PURPOSE: We need to enumrate the CDS out of EVM
# USAGE: enumerateExon.pl < bla.gff > corrected.gff
# AUTHOR: htafer@gmail.com


use strict;
use warnings;

main();

sub main{
  my $count=0;
  while(my $line=<STDIN>){
    if($line=~/Name/){
      $count = 0;
    }
    elsif($line=~/\tCDS\t/){
      $count++;
      my $pattern = "ID=cds".$count;
      $line=~s/ID=cds/$pattern/;
    }
    print $line;
  }

}
