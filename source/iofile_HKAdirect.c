/*
 *  iofile_HKA.c
 *  mlfreqS
 *
 *  Created by Sebastian E. Ramos-Onsins.
 *
 */

#include "hka.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUMBER_COLS 9

int input_datag(FILE *file_in,struct parametersg *data)
{
	/****************************************************************************************/
	/* INPUT FILE WITH THE STRUCTURE:                                                       */
	/* 1st line: Title                                                                      */
	/* 2nd line: nloci                                                                      */
	/* 3er line: header: IDlocus nsam SegSites Divergence length_pol length_div factor_chrn */
	/*	         [factor_an factor_bn theta_square]                                         */ 
	/* 4th and rest: name and values for each locus                                         */
	/****************************************************************************************/
	
	char *f;
    char *line;
    char *number;
    long int c,i,j,k;
	long int nline;
        
    if(!(line = (char *)calloc(10000,sizeof(char)))) {
        perror("calloc error.0");
        return 1;
    }
    if(!(number = (char *)calloc(128,sizeof(char)))) {
        perror("calloc error.0");
        return 1;
    }
    if(!(f = (char *)malloc(BUFSIZ))) {
        puts("\nError: Not enough memory to read  the input file.\n");
        return 1;
    }
	setbuf(file_in,f);
    
	nline = 0;
	do {
		j = 0;
		while(((c=fgetc(file_in)) != 10 && c!=13 && c!=-1 && c!=0) || j==0)
			line[j++] = c;
		if(c==-1 || c==0) {
			c=0;
			break;
		}
		line[j] = '\0';
		if(nline==0) {
			 strcpy(data[0].title,line);
		}
		if(nline == 1) {
			j=0;
			while(line[j] < '0' || line[j] >'9') j++;
			k = 0;
			while(line[j] != 32 && line[j] != 9 && line[j] !=0) {
				number[k] = line[j];
				k++;
				j++;
			}
			number[k] = '\0';
			data[0].nloci = atol(number);
			
			if(data[0].nloci < 2) {
				printf("\nError: nloci must be higher than 1. Exiting.\n");
				exit(1);
			}
			if(!(data[0].idloci = (char **)calloc(data[0].nloci,sizeof(char *)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			for(i=0;i<data[0].nloci;i++) {
				if(!(data[0].idloci[i] = (char *)calloc(128,sizeof(char)))) {
					puts("\nError: Not enough memory.\n");
					return 1;
				}
			}
			if(!(data[0].nsamples = (long int *)calloc(data[0].nloci,sizeof(long int)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].Ssites = (long int *)calloc(data[0].nloci,sizeof(long int)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].divergence = (float *)calloc(data[0].nloci,sizeof(float)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].length_pol = (float *)calloc(data[0].nloci,sizeof(float)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].length_div = (float *)calloc(data[0].nloci,sizeof(float)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].factor_chrn = (float *)calloc(data[0].nloci,sizeof(float)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].factor_pool_an = (float *)calloc(data[0].nloci,sizeof(float)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].factor_pool_bn = (float *)calloc(data[0].nloci,sizeof(float)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
			if(!(data[0].pmissing = (float *)calloc(data[0].nloci,sizeof(float)))) {
				puts("\nError: Not enough memory.\n");
				return 1;
			}
		} 
		if(nline > 2 && nline-3 < data[0].nloci) {
			j = 0;
			i = 0;
			while(line[i]!='\0' && j < NUMBER_COLS) {
				k = 0;
				while(line[i] != 32 && line[i] != 9 && line[i] !=0 && k<128) {
					number[k] = line[i];
					k++;
					i++;
				}
				if(j>=7 && k==0) {
					j++;
					continue;
				}
				number[k] = '\0';
				switch(j) {
					case 0: 
						memcpy(data[0].idloci[nline-3],number,k);
						break;
					case 1: 
						data[0].nsamples[nline-3] = atol(number);
						if(data[0].nsamples[nline-3] < 2) {
							printf("\nError in locus %ld: nsam must be always > 1. Sorry.\n",nline-2);
							exit(1);
						}
						break;
					case 2: 
						data[0].Ssites[nline-3] = atol(number);
						break;
					case 3: 
						data[0].divergence[nline-3] = atof(number);
						break;
					case 4: 
						data[0].length_pol[nline-3] = atof(number);
						if(data[0].length_pol[nline-3] < 2) {
							printf("\nError in locus %ld: length must be always > 0. Sorry.\n",nline-2);
							exit(1);
						}
						break;
					case 5: 
						data[0].length_div[nline-3] = atof(number);
						if(data[0].length_div[nline-3] < 2) {
							printf("\nErrorin locus %ld: length must be always > 0. Sorry.\n",nline-2);
							exit(1);
						}
						break;
					case 6: 
						data[0].factor_chrn[nline-3] = atof(number);
						break;
					case 7: 
						data[0].factor_pool_an[nline-3] = atof(number);
						data[0].pmissing[nline-3] = atof(number);
						data[0].take_anbn = 2;
						break;
					case 8: 
						data[0].factor_pool_bn[nline-3] = atof(number);
						data[0].take_anbn = 1;
						break;
					default:
						break;
				}
				j++;
				while((line[i] < '0' || line[i] >'9') && (line[i] !=0)) 
					i++;
			}
		}
		nline++;
	} while(nline-3 < data[0].nloci);
	if(c==0 && nline-3 < data[0].nloci) {
		printf("\nFormat Error: Less rows (%ld) than loci (%ld)?. Sorry.\n",nline-2,data[0].nloci);
		exit(1);
	}
	    
    return 0;
}

int input_Cij_matrix(float ****Cij_matrix,int *nx)
{
	/****************************************************************************************/
	/* sum of Cov theta values for mising values: filename: "Covij_sumcov_matrix_theta2.txt */
	/* # means comment. skip                                                                */
	/* Each matrix lower diagonal contains in rows the values of nxy from 0 to ny           */
	/* in columns represent the value of nx (fixed) and ny from 0 to nx                     */
	/* Next matrix following with nx inceasing                                              */ 
	/****************************************************************************************/
	
	char *f;
    char *line;
    char *number;
    long int c,j,i,k;
	FILE *file_input_Cij;
	int ny,yy,xy;
	

	if (!(file_input_Cij = fopen ("Covij_sumcov_matrix_theta2.txt","r"))) {
        puts("\nError opening the file 'Covij_sumcov_matrix_theta2.txt'\n");
        return(1);
    }
	
	if(!(line = (char *)calloc(10000,sizeof(char)))) {
        perror("\nError:calloc error 222.\n");
        return 1;
    }
    if(!(number = (char *)calloc(128,sizeof(char)))) {
        perror("\nError:calloc error 226.\n");
        return 1;
    }
    if(!(f = (char *)malloc(BUFSIZ))) {
        puts("\nError:calloc error 230.\n");
        return 1;
    }
	setbuf(file_input_Cij,f);
    
	*nx = 0;
	ny  = 0;
    Cij_matrix[0] = (float ***)calloc(*nx,sizeof(float **));
    Cij_matrix[0][0] = (float **)calloc(1,sizeof(float *));
    Cij_matrix[0][0][0] = (float *)calloc(1,sizeof(float));
	
	do {
		j = 0;
		while(((c=fgetc(file_input_Cij)) != 10 && c!=13 && c!=-1 && c!=0) || j==0)
			line[j++] = c;
		if(c==-1 || c==0) {
			c=0;
			break;
		}
		line[j] = '\0';
		if(line[0]=='#') {
			continue;
		}
		else {
			j = 0;
			i = 0;
			while(line[i]!='\0') {
				k = 0;
				while(line[i] != 32 && line[i] != 9 && line[i] !=0 && line[i] !=-1 && k<128) {
					number[k] = line[i];
					k++;
					i++;
				}
				number[k] = '\0';
				/**/
				Cij_matrix[0][*nx][ny][j] = atof(number);
				/**/
				j++;
				while(line[i] == 32 || line[i] == 9) 
					i++;
			}
			ny++;
			if(ny > *nx) {
				*nx += 1;
				ny = 0;
				Cij_matrix[0] = (float ***)realloc(Cij_matrix[0],(*nx+1)*sizeof(float **));
				for(yy=*nx;yy<*nx+1;yy++) {
					Cij_matrix[0][yy] = (float **)calloc((*nx+1),sizeof(float *));
					for(xy=0;xy<=yy;xy++) {
						Cij_matrix[0][yy][xy] = (float *)calloc(xy+1,sizeof(float));
					}
				}
			}
		}
	} while((c=!0) && (c!=-1));
	
	return 0;
}


