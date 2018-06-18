#!/usr/bin/perl -w
#author: Lin Li
#date: 2017-04-04
#functions: QTG-Seq parser
#input files: VCF file
#output files: Statistical file
#format as follows:
#column 1-chromosome
#column 2-start position
#column 3-end position
#column 4-name


use strict;
use warnings;
use Carp qw< croak >;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;

#require lilinperl;

&main(); # the main interface that the program starts from

sub main(){
	#input parameters
	my ($vcf,$rdepth,$outputfile)=@ARGV;
	unless($outputfile)
	{
		print "Syntax:\n";
		print "./QTG_pretreatment.pl <VCF file> <Read Cov> <output file> \n";
		print "Note:\n";
		print "to do the pretreatment on the QTG-Seq data.\n";
		exit(0);
	}
	my @tarr=getDatafromFile($vcf);
	open OUT,">$outputfile" or die "Cannot create the output file $!";
	my $filteredSNPfile="temp_".time().".vcf";
	open SNPOUT,">$filteredSNPfile" || die "Cannot create temp filtered SNP file:$!";
	
	my $i=0;
	while($i<@tarr){
		if($tarr[$i] =~/^#/){
			#print OUT $tarr[$i];
			$i++;
		}else{
			my @trow=split("\t",trim($tarr[$i]));
			if(($trow[3] eq "A" || $trow[3] eq "T" || $trow[3] eq "C" || $trow[3] eq "G" ) && ($trow[4] eq "A" || $trow[4] eq "T" || $trow[4] eq "C" || $trow[4] eq "G")){
				my $SNPstr=$trow[0].",".$trow[1];
				my $flag=0;
				my $flag2=0;
				#A and a allele in High pool and Low pool
				my $HA=-1;
				my $Ha=-1;
				my $LA=-1;
				my $La=-1;
				if($trow[9]=~/(0|1)\/(0|1)\:\d+,\d+\:\d+\:\d+\:\d+,\d+,\d+/){
					my @garr1=split(/\:/,trim($trow[9]));
					if($garr1[2]<$rdepth){
						$SNPstr.=",NA,NA";
						$flag++;
					}else{
						#$SNPstr.=",".$garr1[1];
						my @Harr=split(/,/,$garr1[1]);
						$HA=$Harr[0];
						$Ha=$Harr[1];
						if($garr1[1]=~/^0,/){
							$flag2++;
						}
					}
				}else{
					$flag++;
				}
				if($trow[10]=~/(0|1)\/(0|1)\:\d+,\d+\:\d+\:\d+\:\d+,\d+,\d+/){
					my @garr2=split(/\:/,trim($trow[10]));
					if($garr2[2]<$rdepth){
						$SNPstr.=",NA,NA";
						$flag++;
					}else{
						#$SNPstr.=",".$garr2[1];
						my @Larr=split(/,/,$garr2[1]);
						$LA=$Larr[0];
						$La=$Larr[1];
						if($garr2[1]=~/^0,/){
							$flag2++;
						}						
					}
				}else{
					$flag++;
				}
				if($flag==0 && $flag2==0){
					my $p1=chi_squared_test(observed => [$HA,$Ha],expected => [($HA+$Ha)/2,($HA+$Ha)/2]);
					my $p2=chi_squared_test(observed => [$LA,$La],expected => [($LA+$La)/2,($LA+$La)/2]);				
					if(($p1>=0.05 && $p2>=0.05)){
						print SNPOUT $tarr[$i];
					}elsif( ($p1<0.05 && $p2<0.05 && (($HA/($HA+$Ha)>0.5 && $LA/($LA+$La)<0.5) || ($HA/($HA+$Ha)<0.5 && $LA/($LA+$La)>0.5)))){
						if($HA<$Ha){
							my @garr1=split(/\:/,trim($trow[9]));
							my @garr2=split(/\:/,trim($trow[10]));
							my $modifiedStr=$trow[0]."\t".$trow[1]."\t".$trow[2]."\t".$trow[3]."\t".$trow[4]."\t".$trow[5]."\t".$trow[6]."\t".$trow[7]."\t".$trow[8]."\t".$garr1[0].":".$Ha.",".$HA.":".$garr1[2].":".$garr1[3].":".$garr1[4]."\t".$garr2[0].":".$La.",".$LA.":".$garr2[2].":".$garr2[3].":".$garr2[4]."\n";
							print SNPOUT $modifiedStr;
						}else{
							print SNPOUT $tarr[$i];
						}
					}
				}
			}
			$i++;
		}
	}
	close(SNPOUT);
	@tarr=(); #release the old memory
	@tarr=getDatafromFile($filteredSNPfile);
	$i=0;
	my $upi=-1;
	my $downi=-1;
	while($i<@tarr){
		if($i==0){
			$downi=1;
			my ($tfreq1,$tfreq2)=getFreq($tarr[$i]);
			my ($downfreq1,$downfreq2)=getFreq($tarr[$downi]);
			if(abs($tfreq1-$downfreq1)<0.1 && abs($tfreq2-$downfreq2)<0.1){
				print OUT $tarr[$i];
				$upi=$i;
			}else{
				if(($tfreq1-0.5)*($tfreq2-0.5)<($downfreq1-0.5)*($downfreq2-0.5)){
					print OUT $tarr[$downi];
					$upi=$downi;
					$i=$downi;
				}else{
					print OUT $tarr[$i];
					$upi=$i;					
				}
			}
		}elsif($i==@tarr-1){
			my ($tfreq1,$tfreq2)=getFreq($tarr[$i]);
			my ($upfreq1,$upfreq2)=getFreq($tarr[$upi]);
			if(abs($tfreq1-$upfreq1)<0.1 && abs($tfreq2-$upfreq2)<0.1){
				print OUT $tarr[$i];
				$upi=$i;
			}
		}else{
			$downi=$i+1;
			my ($tfreq1,$tfreq2)=getFreq($tarr[$i]);
			my ($downfreq1,$downfreq2)=getFreq($tarr[$downi]);
			my ($upfreq1,$upfreq2)=getFreq($tarr[$upi]);
			if(abs($tfreq1-$downfreq1)<0.1 && abs($tfreq2-$downfreq2)<0.1 && abs($tfreq1-$upfreq1)<0.1 && abs($tfreq2-$upfreq2)<0.1){
				print OUT $tarr[$i];
				$upi=$i;
			}
		}
		$i++;
	}
	close(OUT);
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

sub getFreq(){
  my($tstr) = @_;
  # Initialize variables
  #my @freqarr = ();
  #print "tstr in getFreq sub:",$tstr,"\n";
  my @trow=split("\t",trim($tstr));
  my @garr1=split(/\:/,trim($trow[9]));
  my @garr2=split(/\:/,trim($trow[10]));
  my @Harr=split(/,/,$garr1[1]);
  my @Larr=split(/,/,$garr2[1]);
  my $tfreq1=$Harr[0]/($Harr[0]+$Harr[1]);
  my $tfreq2=$Larr[0]/($Larr[0]+$Larr[1]);  

  return ($tfreq1,$tfreq2);
}

#chi-square test computation for p value
sub chi_squared_test {
  my %args = @_;
  my $observed = delete $args{observed} // croak q(Argument "observed" required);
  my $expected = delete $args{expected} // croak q(Argument "expected" required);
  @$observed == @$expected or croak q(Input arrays must have same length);

  my $chi_squared = sum map {
    ($observed->[$_] - $expected->[$_])**2 / $expected->[$_];
  } 0 .. $#$observed;
  my $degrees_of_freedom = @$observed - 1;
  my $probability = chisqrprob($degrees_of_freedom, $chi_squared);
  return $probability;
}

sub trim()
{
	my $string=shift;
	$string=~s/^\s+//;
	$string=~s/\s+$//;
	return $string;
}
