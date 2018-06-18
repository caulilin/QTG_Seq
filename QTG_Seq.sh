########################################################################################
########################QTG-Seq bioinformatic pipeline##################################
########################Author:Dr. Lin Li ##############################################
########################input files: VCF file###########################################
########################with genomic variants ofboth high and low pools#################
########################GFF file########################################################
########################Coverage Threshold##############################################
########################Window Size (Integrater)########################################
########################Date:2018-06-16#################################################

#!/bin/bash

echo 'Reading user inputs...'
if [ $# -ne 4 ] && [ $# -ne 5 ]
then
	echo "usage:./QTG_Seq.sh [VCF] [GFF] [Cov Threshold] [Win Size] <EuclideanDist>"
	echo "QTG-Seq provides several different statistics for analysis as follows:"
	echo "EuclideanDist (default) - Euclidean Distance"
	echo "SNPindex - Delta SNPindex"
	echo "Pvalue - Delta P(Chi-Sqrt)"
	echo "ED4 "
	exit 0
else
	if [ -f "$1" ]; then
		if [ -f "$2" ]; then
			if grep '^[[:digit:]]*$' <<< "$3"; then
				echo "running QTG_pretreatment.pl $1 $3 $1_$3"
				#exit 0
				if ! grep '^[[:digit:]]*$' <<< "$4"; then
				#else
					echo 'Please input a valid integrater value of Window Size'
					echo 'QTG-Seq terminated'
					exit 0				
				fi
				perl QTG_pretreatment.pl $1 $3 $1_$3
			else
				echo 'Please input a valid integrater value of Sequencing coverage'
				echo 'QTG-Seq terminated'
				exit 0
			fi
		else
			echo 'Cannot read GFF file'
			echo 'QTG-Seq terminated'
			exit 0
		fi
	else 
		echo 'Cannot read VCF file with genomic variants'
		echo 'QTG-Seq terminated'
		exit 0
	fi
fi
if [ $# -ne 5 ]
then
	$5='EuclideanDist'
fi

cut -f 1 $1_$3 | sort |uniq >chrname.txt; 

while read -r tchr
do  
	echo $tchr;
	grep -P "$tchr\t" $1_$3 >$1_$3_$tchr; 
	perl QTG_parser.pl $1_$3_$tchr $3 $4 res_$tchr.txt
done <chrname.txt

#extract the statistics and order them by coordinates
perl QTG_summarizer.pl chrname.txt $5 allBinres.csv

#summarize and extract the fine-mapped region of QTG by R
module load R
R CMD BATCH QTG_Seq.R
# Check the output
cat QTG_Seq.Rout

perl QTG_Miner.pl QTG_Seq_R_summaryfile.txt QTG_region.txt

rm chrname.txt
rm QTG_Seq_R_summaryfile.txt

exit 0
