# QTG_Seq

The open source code of pipeline for the QTG-Seq, which is dedicated to accelerating QTL fine-mapping through QTL partitioning and whole genome sequencing on bulked segregant samples

Source code of the QTG-Seq

1, Copyright (c) 2018

Huazhong Agricultural University, All Rights Reserved. Authors: Lin Li, Guoying Wang, Hongwei Zhang, Xi Wang

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

o Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

o Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

o Neither the name of the University of Minnesota nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LIN LI, GUOYING WANG, HONGWEI ZHANG, XI WANG (OR HUAZHONG AGRICULTURAL UNIVERSITY) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

DESCRIPTION: QTG-Seq enables rapid QTL fine-mapping using native 2nd Next Generation Sequence Data derived from the bulked high and low segregant pools. This script essentially uses the results from external alignment programs and performs a series of statistic analyses via a set of specified parameters.

CITATION: QTG-Seq can be cited as: Zhang HW, Wang X, Pan QC, Li P, Li J, Han LQ, Liu YJ, Wang PX, Li DD, Liu Y, Zhang YM, Wang GY, Li L: QTG-seq accelerates QTL fine-mapping through QTL partitioning and whole genome sequencing on bulked segregating samples, 2018.

2, Prerequisite of external softwares

Need to install Perl Statistics module, bwa in your local environment correctly and set the running path of external programs at the begining of the pipeline

3, Usage

usage:./QTG_Seq.sh [VCF] [GFF] [Cov Threshold] [Win Size] \<EuclideanDist\>

QTG-Seq provides several different statistics for analysis as follows:

	EuclideanDist (default) - Euclidean Distance
	
	SNPindex - Delta SNPindex
	
	Pvalue - Delta P(Chi-Sqrt)
	
	ED4 
