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
#modified at 20180423 to add new function that calculates all the statistics within window spanning similar number of markers 
#rather than similar genomic length

use strict;
use warnings;
use Carp qw< croak >;
use List::Util qw< sum >;
use Statistics::Distributions qw< chisqrprob >;

#require lilinperl;

&main(); # the main interface that the program starts from

sub main(){
	#input parameters
	my ($vcf,$rdepth,$wsize,$outputfile)=@ARGV;
	unless($outputfile)
	{
		print "Syntax:\n";
		print "./QTG_Parser.pl <VCF file> <Read Cov> <Window Size> <output file> \n";
		print "Note:\n";
		print "to parse VCF data and fine-map QTG.\n";
		exit(0);
	}
	my @tarr=getDatafromFile($vcf);
	my $summaryfile=$outputfile."_res";
	open OUT,">$outputfile" or die "Cannot create the output file $!";
	open SUMOUT,">$summaryfile" or die "Cannot create result summary file $!";
	print SUMOUT "chr,position,#SNP,Afreq_1stPool,afreq_2ndPool,Afreq_1stPool(median),afreq_2ndPool(median),P_1stPool,P_2ndPool,P_1stPool(median),P_2ndPool(median),EuclideanDist,SNPindex,SNPindex(median),-log10(P1*P2)\n";
	
	my $i=0;
	my $spos=0;
	my $slidingnum=0;
	my $epos=$spos+$wsize*$slidingnum;
	while($i<@tarr){
		if($tarr[$i] =~/^#/){
			$i++;
		}else{
			$slidingnum=0;
			$spos=$i;
			$epos=$i+$wsize;
			my $filteredSNPfile="range_".$spos."_".$epos.".csv";
			#open SNPOUT,">$filteredSNPfile" || die "Cannot create filtered SNP file:$!";
			#print out the header
			#print SNPOUT "chr,position,A_1stPool,a_1stPool,A_2ndPool,a_2ndPool\n";
			print OUT "chr,position,A_1stPool,a_1stPool,A_2ndPool,a_2ndPool,P_1stPool,P_2ndPool,Afreq_1stPool,afreq_1stPool,Afreq_2ndPool,afreq_2ndPool\n";
			my @trow=split("\t",trim($tarr[$i]));
			if($slidingnum==0){
				my $SNPstr=$trow[0].",".$trow[1];
				my $RESstr=$trow[0].",".$trow[1];
				my $flag=0;
				my $flag2=0;
				#A and a allele in High pool and Low pool
				my $HA=-1;
				my $Ha=-1;
				my $LA=-1;
				my $La=-1;
				my $squarefreq=0;
				my @LAfreqarr=();
				my @Hafreqarr=();
				my @Lchiparr=();
				my @Hchiparr=();
				my $tnum=0;
				
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
						$SNPstr.=",".$HA.",".$Ha.",".$LA.",".$La.",".$p1.",".$p2.",".$HA/($HA+$Ha).",".$Ha/($HA+$Ha).",".$LA/($LA+$La).",".$La/($LA+$La);
						#print SNPOUT $SNPstr,"\n";
						print OUT $SNPstr,"\n";	
						push(@LAfreqarr,$HA/($HA+$Ha));
						push(@Hafreqarr,$La/($LA+$La));
						push(@Lchiparr,$p1);
						push(@Hchiparr,$p2);
						$squarefreq+=($HA/($HA+$Ha) - $LA/($LA+$La)) * ($HA/($HA+$Ha) - $LA/($LA+$La));
					}elsif( ($p1<0.05 && $p2<0.05 && (($HA/($HA+$Ha)>0.5 && $LA/($LA+$La)<0.5) || ($HA/($HA+$Ha)<0.5 && $LA/($LA+$La)>0.5)))){
						if($HA<$Ha){
							$SNPstr.=",".$Ha.",".$HA.",".$La.",".$LA.",".$p1.",".$p2.",".$Ha/($HA+$Ha).",".$HA/($HA+$Ha).",".$La/($LA+$La).",".$LA/($LA+$La);
							push(@LAfreqarr,$Ha/($HA+$Ha));
							push(@Hafreqarr,$LA/($LA+$La));
							push(@Lchiparr,$p1);
							push(@Hchiparr,$p2);
							$squarefreq+=($Ha/($HA+$Ha) - $LA/($LA+$La)) * ($Ha/($HA+$Ha) - $LA/($LA+$La));							
						}else{
							$SNPstr.=",".$HA.",".$Ha.",".$LA.",".$La.",".$p1.",".$p2.",".$HA/($HA+$Ha).",".$Ha/($HA+$Ha).",".$LA/($LA+$La).",".$La/($LA+$La);
							push(@LAfreqarr,$HA/($HA+$Ha));
							push(@Hafreqarr,$La/($LA+$La));
							push(@Lchiparr,$p1);
							push(@Hchiparr,$p2);
							$squarefreq+=($HA/($HA+$Ha) - $LA/($LA+$La)) * ($HA/($HA+$Ha) - $LA/($LA+$La));							
						}
						#print SNPOUT $SNPstr,"\n";
						print OUT $SNPstr,"\n";						
					}
					$tnum++;
					$slidingnum++;
				}
				my $j=$i+1;
				while($j<@tarr){
					@trow=split("\t",trim($tarr[$j]));
					$flag=0;
					$flag2=0;
					my $HA=-1;
					my $Ha=-1;
					my $LA=-1;
					my $La=-1;					
					if($slidingnum<$wsize){
						$SNPstr=$trow[0].",".$trow[1];
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
								$SNPstr.=",".$HA.",".$Ha.",".$LA.",".$La.",".$p1.",".$p2.",".$HA/($HA+$Ha).",".$Ha/($HA+$Ha).",".$LA/($LA+$La).",".$La/($LA+$La);
								#print SNPOUT $SNPstr,"\n";
								print OUT $SNPstr,"\n";	
								push(@LAfreqarr,$HA/($HA+$Ha));
								push(@Hafreqarr,$La/($LA+$La));
								push(@Lchiparr,$p1);
								push(@Hchiparr,$p2);
								$squarefreq+=($HA/($HA+$Ha) - $LA/($LA+$La)) * ($HA/($HA+$Ha) - $LA/($LA+$La));							
							}elsif( ($p1<0.05 && $p2<0.05 && (($HA/($HA+$Ha)>0.5 && $LA/($LA+$La)<0.5) || ($HA/($HA+$Ha)<0.5 && $LA/($LA+$La)>0.5)))){
								if($HA<$Ha){
									$SNPstr.=",".$Ha.",".$HA.",".$La.",".$LA.",".$p1.",".$p2.",".$Ha/($HA+$Ha).",".$HA/($HA+$Ha).",".$La/($LA+$La).",".$LA/($LA+$La);
									push(@LAfreqarr,$Ha/($HA+$Ha));
									push(@Hafreqarr,$LA/($LA+$La));
									push(@Lchiparr,$p1);
									push(@Hchiparr,$p2);
									$squarefreq+=($Ha/($HA+$Ha) - $LA/($LA+$La)) * ($Ha/($HA+$Ha) - $LA/($LA+$La));							

								}else{
									$SNPstr.=",".$HA.",".$Ha.",".$LA.",".$La.",".$p1.",".$p2.",".$HA/($HA+$Ha).",".$Ha/($HA+$Ha).",".$LA/($LA+$La).",".$La/($LA+$La);
									push(@LAfreqarr,$HA/($HA+$Ha));
									push(@Hafreqarr,$La/($LA+$La));
									push(@Lchiparr,$p1);
									push(@Hchiparr,$p2);
									$squarefreq+=($HA/($HA+$Ha) - $LA/($LA+$La)) * ($HA/($HA+$Ha) - $LA/($LA+$La));							
								}
								#print SNPOUT $SNPstr,"\n";
								print OUT $SNPstr,"\n";						
							}
							$tnum++;
							$slidingnum++;
						}
						$j++;
					}else{
						last;
					}
				}
				#print length(@LAfreqarr),"\t",$tnum,"\n";
				if(@LAfreqarr>=1){
				#the number of SNP in a specific block will be subject to change
					my $SNPindex=ave(@LAfreqarr)+ave(@Hafreqarr);
					my $SNPindex2=median(@LAfreqarr)+median(@Hafreqarr);
					my $log10P1P2=-log10(ave(@Lchiparr))-log10(ave(@Hchiparr));
					$RESstr.=",".$tnum.",".ave(@LAfreqarr).",".ave(@Hafreqarr).",".median(@LAfreqarr).",".median(@Hafreqarr).",".ave(@Lchiparr).",".ave(@Hchiparr).",".median(@Lchiparr).",".median(@Hchiparr).",".sqrt($squarefreq).",".$SNPindex.",".$SNPindex2.",".$log10P1P2;
					print SUMOUT $RESstr,"\n";
				}
				$i=$j;
			}else{
				$i++;
			}
			# recall the R script for the calculation of Bayesian proportion
			
			# the end of recalling of R scripts
			
			
		}
	}

	close(OUT);
	close(SUMOUT);
}
sub log10 {
	my $n = shift;
	return log($n)/log(10);
}
#----------------------------------------------------------------------
#获取一个数组的平均数
sub ave{
    return sum(@_)/@_;
}
#----------------------------------------------------------------------
#获取一个数组的中位数
#----------------------------------------------------------------------
sub median
{
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
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
