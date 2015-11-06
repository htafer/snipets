#!/usr/bin/env perl                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
# PURPOSE: From the swissprot and ipr file, generate protein names that can be then fed to the testfeatures.pl script                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
#perl ./parseNameFromIPRSP.pl -i debFab.merged.tsv -s ~/share/database/swissProtFungi.fasta -b debFab.swissProt  > debFab.names                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
# USAGE: gff2cdsfa.pl -fasta t/VV10.fa -gff t/t1.gff                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
# AUTHOR: malcolm.cook@gmail.com                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              

use strict;
use warnings;
use Getopt::Std;

use Data::Dumper;

my $nameHash;
my %opts;
my $ACCID="bla";

sub init{
  use Getopt::Std;
  getopt("i:s:b",\%opts);
  usage() if(!(defined($opts{s})) && !(defined($opts{i})) && !(defined($opts{b})));
}

sub usage{
  print STDERR << "EOF";
Input:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
-i ipr file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
-s swissProt file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
-b swissprot blast result                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
usage:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
in /media/work/genomes/debFab/annotation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
perl ./parseNameFromIPRSP.pl -i debFab.merged.tsv -s ~/share/database/swissProtFungi.fasta -b debFab.swissProt  > debFab.names                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
EOF                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
  
  exit;
}
init();
main();

sub main{
  parseBlastFile();
  parseIPRFile();
  getName();
  #parse blastfile                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  #parse inteproscan                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
  #join files                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
  #print table                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
}

sub getName{
  for my $name (keys %{$nameHash}){
    if(defined($nameHash->{$name}->{BLAST})){
      print $name,"\t",$nameHash->{$name}->{BLAST},"\n";
    }
    elsif(defined($nameHash->{$name}->{TIGRFAM})){
      my $desc = $nameHash->{$name}->{TIGRFAM};
      $desc=~s/.+://;
      $desc=~s/region/domain/g;
      if($desc=~/(domain|repeat)/){
        $desc=~s/domain.*/domain-containing protein/;
        $desc=~s/repeat.*/repeat-containing protein/;
        print $name,"\t",$desc,"\n";
      }
    }
    elsif(defined($nameHash->{$name}->{Pfam})){
      my $desc = $nameHash->{$name}->{Pfam};
      $desc=~s/.+://;
      $desc=~s/region/domain/g;
      if($desc=~/(domain|repeat)/){
        $desc=~s/domain.*/domain-containing protein/;
        $desc=~s/repeat.*/repeat-containing protein/;
        print $name,"\t",$desc,"\n";
      }
    }
  }
}



sub  parseIPRFile{
  open(IPR,"<$opts{i}");
  while(my $line=<IPR>){
    next unless $line=~/(Pfam|TIGRFAM)/;
    my @F=split(/\t/,$line);
    if(!defined($nameHash->{$F[0]}->{$F[1]})){
      $nameHash->{$F[0]}->{$F[3]}=$F[5];
    }
  }
  close(FILE);
}


sub parseBlastFile{
  open(SWISSPROT,"<$opts{s}");
  #Parse swissprot                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  my $swissProtHash;
  while(my $line=<SWISSPROT>){
    next unless($line=~/>/);
    $line=~s/>//;
    $line=~s/ OS=.+//;
    $line=~/([^\s]+)\s(.+)/;
    $swissProtHash->{$1}=$2;
  }
  close(SWISSPROT);
  #Get annotation id <->swiss prot gene correspondance                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
  open(FILE,"<$opts{b}");
  #Get full gene name from the swissprot database;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
  while(my $line=<FILE>){
    next if($line=~/#/);
    my @F = split(/\t/,$line);
    if(!defined($nameHash->{$F[0]}) && abs(1-$F[1]/$F[5])<=0.1 && abs($F[3]-$F[2])/$F[1] >= 0.7 && abs($F[6]-$F[7])/$F[5] >=0.7 && $F[8]>=70){
      if(defined($swissProtHash->{$F[4]})){
        $nameHash->{$F[0]}->{BLAST}=$swissProtHash->{$F[4]};
      }
    }
  }
  close(FILE);
}



























