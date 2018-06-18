#!/usr/bin/perl -w
#author: Lin Li
#date: 2018-06-16
#functions: QTG_summarizer
#input files: VCF file
#output files: Statistical file
#format as follows:
#column 1-chromosome
#column 2-start position
#column 3-end position
#column 4-name
#modified at 20180423 to add new function that calculates all the statistics within window spanning similar number of markers 
#rather than similar genomic length

use strict;
use warnings;
use Carp qw< croak >;
use List::Util qw< sum >;
#use Statistics::Distributions qw< chisqrprob >;

#require lilinperl;

&main(); # the main interface that the program starts from

sub main(){
	#input parameters
	my ($chrname,$tstat,$outputfile)=@ARGV;
	unless($outputfile)
	{
		print "Syntax:\n";
		print "./QTG_summarizer.pl <ChrName file> <Statistic> <output file> \n";
		print "Note:\n";
		print "to do final QTG mining.\n";
		exit(0);
	}
	my @namearr=getDatafromFile($chrname);
	#my $finalbinresfile=$outputfile;
	#my $QTGfile=$outputfile."_target";
	open OUT,">$outputfile" or die "Cannot create the output file $!";
	#open SUMOUT,">$QTGfile" or die "Cannot create QTG target region file $!";
	#print SUMOUT "chr,position,$tstat,Mark\n";
	print OUT "chr,position,$tstat,Mark,Coord\n";
	
	my $i=0;
	my @tchrarr=();
	my @tposiarr=();
	my @tstatarr=();
	my @tmarkarr=();
	my @tcoordarr=();
	my @tmaxarr=();
	
	my $f=0;
	while($i<@namearr){
		my $tempfile="res_".trim($namearr[$i]).".txt_res";
		my @temparr=getDatafromFile($tempfile);
		my $j=1;
		my $c=11;
		if($tstat eq "SNPindex"){
			$c=12;
		}elsif($tstat eq "Pvalue"){
			$c=14;
		}elsif($tstat eq "ED4"){
			$c=21;
		}
		my $chrabspos=0;
		for (my $m=0;$m<$i;$m++){
			$chrabspos+=$tmaxarr[$m]#+1000000;
		}
		print $chrabspos,"\n";
		while($j<@temparr){
			#print $temparr[$j];
			my @tarr=split(",",trim($temparr[$j]));
			$tchrarr[$f]=$tarr[0];
			$tposiarr[$f]=$tarr[1];
			if($c==21){
				$tstatarr[$f]=$tarr[$c-10]**4;
			}else{
				$tstatarr[$f]=$tarr[$c];
			}
			if(!defined($tmaxarr[$i])){
				$tmaxarr[$i]=$tarr[1];
			}elsif($tmaxarr[$i]<$tarr[1]){
				$tmaxarr[$i]=$tarr[1];
			}
			$tcoordarr[$f]=$tarr[1]+$chrabspos;
			$tmarkarr[$f]=0;
			$f++;
			$j++;
		}
		$i++;
	}
	#my @sorted_indexes = sort { $tstatarr[$b] <=> $tstatarr[$a] } 0..$#tstatarr;
	#my $pstat95=$tstatarr[$sorted_indexes[int(@sorted_indexes*0.05-0.5)]];
	#if($pstat95>3*median(@tstatarr)){
		$f=0;
		while($f<@tstatarr){
			#if($tstatarr[$f]>=$pstat95){
			#	$tmarkarr[$f]=1;
			#}
			print OUT $tchrarr[$f],",",$tposiarr[$f],",",$tstatarr[$f],",",$tmarkarr[$f],",",$tcoordarr[$f],"\n";
			$f++;
		}
	#}
	close(OUT);
	#close(SUMOUT);
}

#-----------------------------------------------------------------------
#获取文本文件的数据
#
#
#-----------------------------------------------------------------------
sub getDatafromFile(){
  my($filename) = @_;
  # Initialize variables
  my @filedata = ();
  unless(open(GET_FILE_DATA, $filename)) {
      print STDERR "Cannot open file \"$filename\"\n\n";
      exit;
  }
  @filedata = <GET_FILE_DATA>;
  close GET_FILE_DATA;
  return @filedata;	
}

sub trim()
{
	my $string=shift;
	$string=~s/^\s+//;
	$string=~s/\s+$//;
	return $string;
}
