------------------------------------------------------------
HKAdirect version beta 20150716
Sebastian E. Ramos-Onsins, Emanuelle Raineri and Luca Ferretti
------------------------------------------------------------

This program computes the HKA from a dataset table for a sampled population 
and a single outgroup.

The program computes the expected polymorphism and divergence as well as the 
theta values per nucleotide, the Time to the ancestor, the partial HKA for each
locus (window), the Chi-square and the P-value. 

NOTE: Check the final Chi-square result by simulation methods. The final result follows approximately a Chi-square distribution.

COMPILE the code: 
Go to source folder and type

gcc *.c -lm -o HKAdirect -lm -Wall -pedantic

RUN HKAdirect:

Two options:
1. Double Click on it
2. From the command line, type
 ./HKAdirect [input_file] > [output_file]
 
INPUT FILE:

The input file must be the following format

1st line: Title
2nd line: nloci 
3er line: header: IDlocus nsam SegSites Divergence length_pol length_div factor_chrn [%missing]
4th and rest: name and values for each locus 
Note: [] is  an optional value.

EXAMPLE INPUT FILE:

Checking HKAdirect: Hudson et al. 1987
nloci 2
IDloci	nsam	S	D	Lp	Ld	fchrn
Adh5'	81	9	210	414	4052	1
Adh	81	8	18	79	324	1


EXAMPLE FILES:

In the folder example two input files with two loci (or windows) are included. 
Note that the factor_chr can be 1.0 (Autosome), 0.75 (X-linked loci), 0.25 (Y-linked loci)
REMEMBER that the Divergence is NOT only the number of fixed positions but the number of variants
between any sample and the outgroup.

OUTPUT RESULTS

The program computes the expected polymorphism and divergence as well as the 
theta values per nucleotide, the Time to the ancestor, the partial HKA for each
locus (window), the Chi-square and the P-value.

EXAMPLE OUTPUT FILE: date Tue May 15 14:54:08 2012 

HKAdirect beta(20120521).
Sebastian E. Ramos-Onsins.

OUTPUT FILE: date Mon May 21 18:08:35 2012 

Input data from the file: ./checkHKA.txt

Title: Checking HKAdirect: Hudson et al. 1987 
nloci: 2

#IDloci	nsam	obs_S	obs_div	length_pol	length_div	factor_chrn	expHKA_S	expVar_S	expHKA_div	expVar_D	expHKA_theta	partialHKA
Adh5'	81	9	210.00	414.00	4052.00	1.00	13.48	25.51	205.52	911.61	0.00656	0.809	
Adh	81	8	18.00	79.00	324.00	1.00	3.52	4.34	22.48	30.93	0.00897	5.276	

Time to the ancestor (in 2N generations): 6.734	Chi-square value: 6.085	

Significance of Chi-square: Significant, P(dgf=1) = 1.363e-02


Good luck.

Sebas

------------------------------------------------------------
Sebastian E. Ramos-Onsins, PhD
Ramon y Cajal Research Position

Consorci CSIC-IRTA-UAB (CRAG)
Centre for Research in Agricultural Genomics
Dept. Animal Genetics
Despatx 307, Edifici CRAG, Campus UAB
08193 Bellaterra, SPAIN

Phone: +34 93 563 6600 Ext 3348
Fax: +34 93 563 66 01

email: sebastian.ramos@cragenomica.es
http://bioinformatics.cragenomica.es/numgenomics/people/sebas
-----------------------------------------------------------
