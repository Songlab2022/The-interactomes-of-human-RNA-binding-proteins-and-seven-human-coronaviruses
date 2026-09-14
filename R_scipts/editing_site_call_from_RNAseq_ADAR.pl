#! /bin/env perl 

use strict;
use Getopt::Long; 
my ($in,$out,$help);
GetOptions("in=s"=>\$in,"out=s"=>\$out,"help"=>\$help);
if($in and $out){

my $referencegenome = "";
my $alu_editing_sites = "";

$in=substr($in,-1) eq '/'?$in:$in.'/';
opendir IN,$in;
$out=substr($out,-1) eq '/'?$out:$out.'/';
foreach(readdir IN){
  if(m/\.bam$/){
    my $inbam=$in.$_;
    my $key=substr($_,0,index($_,'.'));
    my $outk=$out.$key;
    my $outfk=$outk.'/'.$key;
    `mkdir $outk` if (!-e $outk);
    `samtools index $inbam`;
    `perl Query_Editing_Level.pl $alu_editing_sites $inbam $outfk.all.levels.txt`;
   }
}
 closedir IN;
}elsif($help){
   print "The pipline for RNA-seq editing sites calling!\n\n
   --in:the input path;\n
   --out:the result path;\n
   --help:the help!\n\n";
}else{
   print "Please read the help!\n\n";
}

