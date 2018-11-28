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
	my ($statfiles,$outputfile)=@ARGV;
	unless($outputfile)
	{
		print "Syntax:\n";
		print "./QTG_Miner.pl <Statistic files> <output file> \n";
		print "Note:\n";
		print "to do final QTG mining.\n";
		exit(0);
	}
	my @namearr=getDatafromFile($statfiles);
	#my $finalbinresfile=$outputfile;
	#my $QTGfile=$outputfile."_target";
	open OUT,">$outputfile" or die "Cannot create the output file $!";
	#open SUMOUT,">$QTGfile" or die "Cannot create QTG target region file $!";
	#print SUMOUT "chr,position,$tstat,Mark\n";
	print OUT "QTGID,Chr,Left,Peak,Right\n";
	my $qtgid=0;
	my $i=0;
	my $tmax=0;
	my $tmposi=0;
	while($i<@namearr){
		my $tempfile=trim($namearr[$i]);
		my @temparr=getDatafromFile($tempfile);
		my $j=1;
		while($j<@temparr){
			#print $temparr[$j];
			my @tarr=split(",",trim($temparr[$j]));
			if($tarr[3]==0){
				$j++;
			}else{
				$tmax=$tarr[5];
				$tmposi=$j;
				my ($leftcoord,$peakcoord,$rightcoord)=(0,0,0);
				my $QTGchr="";
				my $m=$j+1;
				while($m<@temparr){
					my @ttarr=split(",",trim($temparr[$m]));
					if($ttarr[3]==1){
						if($tmax<$ttarr[5]){
							$tmax=$ttarr[5];
							$tmposi=$m;
							$peakcoord=$ttarr[1];
							$QTGchr=$ttarr[0];
						}
					}else{
						last;
					}
					$m++;
				}
				$j=$m;
				#searching the upstream left boundary
				my $upi=$tmposi-1;
				my $leftposi=$tmposi;
				while($upi>=1){
					my @ttarr=split(",",trim($temparr[$upi]));
					my $iszero=0;
					my $is0posi=0;
					my $is0dist=0;
					if($ttarr[3]==1 && $tmax-$ttarr[5]>=$ttarr[7] && $iszero==0){
						$leftposi=$upi;
						$leftcoord=$ttarr[1];
						last;
					}elsif($ttarr[3]==0 && $tmax-$ttarr[5]>=$ttarr[7] && $iszero==0){
						$leftposi=$upi;
						$leftcoord=$ttarr[1];
						last;
					}elsif($ttarr[3]==0 && $tmax-$ttarr[5]<$ttarr[7] && $iszero==0){
						$iszero=1;
						$is0posi=$upi;
						$is0dist=$tmax-$ttarr[5];
						$leftcoord=$ttarr[1];
					}elsif($ttarr[3]==0 && $tmax-$ttarr[5]<$ttarr[7] && $iszero==1){
						if($is0dist<$tmax-$ttarr[5]){
							$is0dist=$tmax-$ttarr[5];
							$is0posi=$upi;
							$leftcoord=$ttarr[1];
						}
					}elsif($ttarr[3]==1 && $iszero==1){
						$leftposi=$is0posi;
						last;
					}else{
						$leftposi=$upi;
						$leftcoord=$ttarr[1];
					}
					$upi--;
				}
				#searching the downstream right boundary
				my $downi=$tmposi+1;
				my $rightposi=$tmposi;
				while($downi<@temparr){
					my @ttarr=split(",",trim($temparr[$downi]));
					my $iszero=0;
					my $is0posi=0;
					my $is0dist=0;
					if($ttarr[3]==1 && $tmax-$ttarr[5]>=$ttarr[7] && $iszero==0){
						$rightposi=$downi;
						$rightcoord=$ttarr[1];
						last;
					}elsif($ttarr[3]==0 && $tmax-$ttarr[5]>=$ttarr[7] && $iszero==0){
						$rightposi=$downi;
						$rightcoord=$ttarr[1];
						last;
					}elsif($ttarr[3]==0 && $tmax-$ttarr[5]<$ttarr[7] && $iszero==0){
						$iszero=1;
						$is0posi=$downi;
						$is0dist=$tmax-$ttarr[5];
						$rightcoord=$ttarr[1];
					}elsif($ttarr[3]==0 && $tmax-$ttarr[5]<$ttarr[7] && $iszero==1){
						if($is0dist<$tmax-$ttarr[5]){
							$is0dist=$tmax-$ttarr[5];
							$is0posi=$downi;
							$rightcoord=$ttarr[1];
						}
					}elsif($ttarr[3]==1 && $iszero==1){
						$rightposi=$is0posi;
						
						last;
					}else{
						$rightposi=$downi;
						$rightcoord=$ttarr[1];
					}
					$downi++;
				}
				$qtgid++;
				print OUT $qtgid,",",$QTGchr,",",$leftcoord,",",$peakcoord,",",$rightcoord,"\n"; 
			}
			
		}
		$i++;
	}
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
