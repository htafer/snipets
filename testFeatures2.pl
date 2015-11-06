#!/usr/bin/perl 
# takes a gff / fasta and functional annotation file and turn it into a ncbi-compatible tbl file for genome submission                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
#FROM SWISSPROT 'cat res | grep -v "#" | perl -lane 'BEGIN{$name; }if(!defined($name->{$F[0]}) && abs(1-$F[1]/$F[5])<=0.1 && abs($F[3]-$F[4])/$F[1] >= 0.7 && abs($F[6]-$F[7])/$F[5] >=0.7 && $F[8]>=70){$name->{$F[0]}=1; my $id = `grep "$F[4]" swissProtFungi.fasta`;$id=~s/>[^ ]+//; $id=~s/OS=.+//; print "$F[0] $id" }' '                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               

#FROM interprocsan annotation                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
#'cat debFabLargeContigs.gff  | perl ../enumerateExon.pl | perl ../testFeatures2.pl -F debFabLargeContigs.fsa -g /dev/stdin -a AC631 -f ../debFab.merged.tsv > debFabLargeContigs.tbl && linux64.tbl2asn -t template.sbt -p . -a s -V v -j "[organism=Debaryomyces fabryi][strain=CBS789]" -i debFabLargeContigs.fsa'                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         

use strict;
use warnings;
use Getopt::Std;
use List::MoreUtils 'uniq';

use Bio::SeqIO;
use Bio::DB::SeqFeature::Store;
use Bio::SeqFeatureI;

use Data::Dumper;


my %opts;
my $ACCID="bla";

sub init{
  use Getopt::Std;
  getopt("F:g:a:f:n",\%opts);
  usage() if(!(defined($opts{F}) && defined($opts{g})));
}

sub usage{
  print STDERR << "EOF";
This program return the sequence and structure of fusion transcripts                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          
Input:                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        
-f functional annotation file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
-F fasta file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
-a accessionID                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
-g gff file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   
-n name file                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
Usage                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
cd /media/work/genomes/debFab/annotation/asm                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
Get name dictionary                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
-h         : this (help) message                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              
EOF                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
  exit;
}
init();
main();

sub main{
  #Set accession ID                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
  $ACCID=(defined($opts{a}) ? $opts{a} : "bla");
  my $funcHash = parseMerged();
  my $nameHash = parseNames();
  #Set feature DB                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
  my $db = Bio::DB::SeqFeature::Store->new(-fasta => $opts{F},
                                           -gff   => $opts{g},
                                           -adaptor=> 'memory') or die $!;
  my @seqIds = $db->seq_ids();
  my $count=0;
  for my $seqId (sort {$a cmp $b} @seqIds){
    print ">Features $seqId\n";
    my @ss = $db->features(-type =>  'gene', -seqid=> $seqId);
    #Sort the features based on the contig and the start position                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
    @ss = sort {$a->start() <=> $b->start()} @ss;
    #Now Print the element                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
    for my $gene (@ss) {
      $count++;
      my $gC = sprintf("%05d",$count);
      printTBL($gene,$gC,$funcHash,$nameHash);
    }
  }
}

sub parseNames{
  open(FILE,$opts{n}) || die "cant get $opts{n}";
  my $nameHash;
  while(my $line=<FILE>){
    chomp($line);
    my @data = split("\t",$line);
    $nameHash->{$data[0]} = $data[1];
  }
  return $nameHash;
}


sub parseMerged{
  open(FILE,$opts{f}) || die "cant get $opts{f}";
  my $funcHash;
  while(my $line=<FILE>){
    my @data = split("\t",$line);
    if($line=~/IPR/){
      push @{$funcHash->{$data[0]}->{InterPro}} , $data[11];
    }
    if($line=~/Pfam/){
      push @{$funcHash->{$data[0]}->{PFAM}} , $data[4];
    }
    if($line=~/PIRSF/){
      push @{$funcHash->{$data[0]}->{PIR}} , $data[4];
    }
    if($line=~/TIGRFAM/){
      push @{$funcHash->{$data[0]}->{TIGRFAM}} , $data[4];
    }
    if($line=~/(GO:\d+[^\s]*\s)/){
      my $gos=$1;
      $gos=~s/GO://g;
      my @dataGO=split(/\|/,$gos);
      push @{$funcHash->{$data[0]}->{GO}},@dataGO;
    }
  }
  return $funcHash;
}

sub printTBL{
  my $gene = shift;
  my $gC = shift;
  my $funcHash=shift;
  my $nameHash=shift;
  my $geneOutput;
  $geneOutput->{partial5}="";
  $geneOutput->{partial3}="";
  #Coding Gene                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  if(defined($gene->get_SeqFeatures('mRNA'))){
    for my $mRNA ($gene->get_SeqFeatures('mRNA')){
      #We need to get the phase in order to translate correctly                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               
      my @exons = $mRNA->get_SeqFeatures('CDS');
      @exons = sort {$a cmp $b} @exons;
      my $offset = (defined($exons[0]->phase) ? $exons[0]->phase()+1 : 1);
      my $seq = $mRNA->spliced_cds()->translate(-offset => $offset)->seq;

      $geneOutput->{partial5}="<" if(substr($seq,0,1) ne "M");
      $geneOutput->{partial3}=">" if(substr($seq,-1,1) ne "*");
    }
    #Now we can start printing the tbl format                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
    ####################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    #Print Gene        #                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    ####################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    print $geneOutput->{partial5},
      ($gene->strand<0?$gene->end():$gene->start()),
      "\t",$geneOutput->{partial3},
      ,($gene->strand>0?$gene->end():$gene->start()),
      "\tgene\n";
    print "\t\t\tlocus_tag\t$ACCID"."_".$gC,"\n";

    ####################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    #Print mRNA        #                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    ####################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
    for my $mRNA ($gene->get_SeqFeatures('mRNA')){
      printMRNA($mRNA,$gC,$nameHash,$funcHash,$geneOutput);
      ####################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      #Print CDS         #                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      ####################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    
      printCDS($mRNA,$nameHash,$funcHash,$geneOutput,$gC);
    }
  }
}


sub printMRNA{
  my $mRNA=shift;
  my $gC = shift;
  my $nameHash=shift;
  my $funcHash=shift;
  my $geneOutput=shift;
  my $mRNAID = $mRNA->load_id;
  #We need to extract the exon                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                
  my @exons = $mRNA->get_SeqFeatures('exon');
  @exons = sort {$a cmp $b} @exons;
  if(scalar(@exons) == 1) {
    print $geneOutput->{partial5},
      ($exons[0]->strand<0?$exons[0]->end():$exons[0]->start()),
        "\t",$geneOutput->{partial3},
          ($exons[0]->strand>0?$exons[0]->end():$exons[0]->start()),
            "\tmRNA\n";
    print "\t\t\tprotein_id\tgnl|ncbi|$ACCID"."_".$gC,"\n";
    print "\t\t\ttranscript_id\tgnl|ncbi|mRNA.$ACCID"."_".$gC,"\n";
    print "\t\t\tcodon_start\t",$exons[0]->phase()+1  ,"\n" if(defined($exons[0]->phase) && $exons[0]->phase!=0);
    #Print Product                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    my $product = (defined($nameHash->{$mRNAID}) ? $nameHash->{$mRNAID} : "hypothetical protein");
    print "\t\t\tproduct\t$product\n";
  }
  else{
    print $geneOutput->{partial5},($exons[0]->strand<0?$exons[0]->end():$exons[0]->start()),
      "\t",($exons[0]->strand>0?$exons[0]->end():$exons[0]->start()),
        "\tmRNA\n";
    for(my $i=1; $i<=$#exons-1;$i++){
      print "",($exons[$i]->strand<0?$exons[$i]->end():$exons[$i]->start()),"\t",($exons[$i]->strand>0?$exons[$i]->end():$exons[$i]->start()),"\n";
    }
    print "",($exons[-1]->strand<0?$exons[-1]->end():$exons[-1]->start()),"\t",$geneOutput->{partial3},($exons[-1]->strand>0?$exons[-1]->end():$exons[-1]->start()),"\n";
    print "\t\t\tprotein_id\tgnl|ncbi|$ACCID"."_".$gC,"\n";
    print "\t\t\ttranscript_id\tgnl|ncbi|mRNA.$ACCID"."_".$gC,"\n";
    print "\t\t\tcodon_start\t",$exons[0]->phase()+1 ,"\n" if(defined($exons[0]->phase) && $exons[0]->phase!=0);
    #Print Product                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    my $product = (defined($nameHash->{$mRNAID}) ? $nameHash->{$mRNAID} : "hypothetical protein");
    print "\t\t\tproduct\t$product\n";
  }
}



sub printCDS{
  my $mRNA = shift;
  my $nameHash = shift;
  my $funcHash = shift;
  my $geneOutput = shift;
  my $gC = shift;
  my $mRNAID=$mRNA->load_id;
  my @CDS = $mRNA->CDS;
  @CDS = sort {$a cmp $b} @CDS;
  if(scalar(@CDS) == 1){
    print $geneOutput->{partial5},
      ($CDS[0]->strand<0?$CDS[0]->end():$CDS[0]->start()),
        "\t",$geneOutput->{partial3},
          ($CDS[0]->strand>0?$CDS[0]->end():$CDS[0]->start()),
            "\tCDS\n";
    print "\t\t\tprotein_id\tgnl|ncbi|$ACCID"."_".$gC,"\n";
    print "\t\t\ttranscript_id\tgnl|ncbi|mRNA.$ACCID"."_".$gC,"\n";
    print "\t\t\tcodon_start\t",$CDS[0]->phase()+1,"\n" if($CDS[0]->phase);
    #Print Product                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            
    my $product = (defined($nameHash->{$mRNAID}) ? $nameHash->{$mRNAID} : "hypothetical protein");
    print "\t\t\tproduct\t$product\n";
    if(defined($funcHash->{$mRNAID})){
      for my $func(keys %{$funcHash->{$mRNAID}}){
        for my $element(uniq @{$funcHash->{$mRNAID}->{$func}}){
          print "\t\t\t";
          print "db_xref\t$func:$element\n";
        }
      }
    }
  }
  else{
    print $geneOutput->{partial5},($CDS[0]->strand<0?$CDS[0]->end():$CDS[0]->start()),
      "\t",($CDS[0]->strand>0?$CDS[0]->end():$CDS[0]->start()),
        "\tCDS\n";
    for(my $i=1; $i<=$#CDS-1;$i++){
      print "",($CDS[$i]->strand<0?$CDS[$i]->end():$CDS[$i]->start()),"\t",($CDS[$i]->strand>0?$CDS[$i]->end():$CDS[$i]->start()),"\n";
    }
    print "",($CDS[-1]->strand<0?$CDS[-1]->end():$CDS[-1]->start()),"\t",$geneOutput->{partial3},($CDS[-1]->strand>0?$CDS[-1]->end():$CDS[-1]->start()),"\n";
    print "\t\t\tprotein_id\tgnl|ncbi|$ACCID"."_".$gC,"\n";
    print "\t\t\ttranscript_id\tgnl|ncbi|mRNA.$ACCID"."_".$gC,"\n";
    print "\t\t\tcodon_start\t",$CDS[0]->phase()+1,"\n" if($CDS[0]->phase);
    my $product = (defined($nameHash->{$mRNAID}) ? $nameHash->{$mRNAID} : "hypothetical protein");
    print "\t\t\tproduct\t$product\n";
    if(defined($funcHash->{$mRNAID})){
      for my $func (keys %{$funcHash->{$mRNAID}}){
        for my $element (uniq @{$funcHash->{$mRNAID}->{$func}}){
          print "\t\t\t";
          print "db_xref\t$func:$element\n";
        }
      }
    }
  }
}

  







