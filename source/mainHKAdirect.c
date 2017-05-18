/*
 *  Program created to calculate directly HKA (functions from MuLoNe program)
 *
 *  Created by Sebastian E. Ramos-Onsins.
 *
 */

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "hka.h"

int main (int argc,char *argv[]) {

    FILE *file_input;
    FILE *file_output;
    char *file_in;
    char *file_out;
    
	struct parametersg *data;
    struct statistics *matrix;
    struct statmulo matrixml;
    
    int input_datag(FILE *,struct parametersg *);
	int input_Cij_matrix(float ****,int *);
	int dohka_obs(struct statistics *,struct statmulo,int , long int , long int *, float ***);
    void print_prob_obshka(struct statmulo, long int, FILE *);
    
    long int i, nlocihka;
    int outgroup = 1;
	
	float ***Cij_matrix=0;
	int maxCijnsam,maxnsam;
	/*int xx,yy;*/
    
    time_t start;
    time_t end;
    struct tm *date;
    char s[80];
    int hour,min,sec;
    
    if(!(data = (struct parametersg *)calloc(1,sizeof(struct parametersg)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(file_in = (char *)calloc(1024,sizeof(char)))) {
        puts("calloc error.main.00");
        exit(1);
    }
    if(!(file_out = (char *)calloc(1024,sizeof(char)))) {
        puts("calloc error.main.00");
        exit(1);
    }

    if(argc != 1 && argc != 2) {
        puts("usage: HKAdirect input_file > output_file");
        exit(1);
    }

    if(argc == 2) {
        strcpy(file_in,argv[1]);
		/*
		if (!(file_input = fopen (file_in,"r"))) {
			puts("Error in input/output");
			exit(1);
		} 
		*/
    }    
    if(argc == 1) { 
		/*file_input = stdin;*/
		/**/
		printf(HKAdirect);
		printf("\nINTRODUCE INPUT FILE:\n");
		printf("1st line: Title\n");
		printf("2nd line: nloci \n");
		printf("3er line: header: IDlocus nsam SegSites Divergence length_pol length_div factor_chrn [pmissing]\n");
		printf("4th and rest: name and values for each locus \n");
		printf("Note: [optional parameters] \n");

        puts("\nInput file?");
        scanf("%s",file_in);
		/**/
    }
	/**/
    if (!(file_input = fopen (file_in,"r"))) {
        puts("Error in input/output");
        exit(1);
    }    
	/**/
	/******** INTRODUCE THE DATA *********/
    if((i=input_datag(file_input,data)) != 0)
        exit(1);
  
    fclose(file_input);

    if(argc == 1) {
        puts("\nOutput file?");
        scanf("%s",file_out);	
		if (!(file_output = fopen (file_out,"w"))) {
			puts("Error in input/output");
			exit(1);
		}    
    }
	else {
		file_output = stdout;
	}
    /*Starting time*/
    time(&start);
    date = localtime(&start);
    strftime(s,80,"%c",date);

    /*introduce in matrix all the values that are in data[0]. Define matrixml*/
    if(!(matrixml.Shka = (float *)calloc(1,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.hka_T = (float *)calloc(1,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.signif_hka= (float *)calloc(1,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.hka_theta = (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.varS0 = (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.varD = (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.varSmiss = (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.hka= (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.Sobshka= (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.Sexphka= (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.Dobshka= (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrixml.Dexphka= (float *)calloc(data[0].nloci,sizeof(float)))) {
        perror("calloc error.0");
        exit(1);
    }
    if(!(matrix = (struct statistics *)calloc(data[0].nloci,sizeof(struct statistics )))) {
        perror("calloc error.0");
        exit(1);
    }
    
    for(i=0;i<data[0].nloci;i++) {
        matrix[i].nsamples = data[0].nsamples[i];
		matrix[i].biallsitesn = data[0].Ssites[i];
        matrix[i].ndivergence = data[0].divergence[i];
        matrix[i].length_pol = data[0].length_pol[i];
        matrix[i].length_div = data[0].length_div[i];
        matrix[i].factor_chrn = data[0].factor_chrn[i];
		matrix[i].take_anbn = data[0].take_anbn;
		matrix[i].factor_pool_an = data[0].factor_pool_an[i];
		matrix[i].factor_pool_bn = data[0].factor_pool_bn[i];
		matrix[i].pmissing = data[0].pmissing[i];
    }
	
	/* if missing, include the matrix Cij. Look for maxnsam. If maxnsam of Cij < maxnsam data then Cij=0*/
	
	maxnsam = 0;
	for(i=0;i<data[0].nloci;i++) {
		if(matrix[i].pmissing > 0.0) {
			if(matrix[i].nsamples > maxnsam)
				maxnsam = matrix[i].nsamples;
		}
	}
	if((i=input_Cij_matrix(&Cij_matrix,&maxCijnsam)) != 0) {
		/*in case o file or error, not use*//*
		for(xx=0;xx<maxCijnsam+1;xx++) {
			for(yy=0;yy<xx+1;yy) {
				free(Cij_matrix[xx][yy]);
			}
			free(Cij_matrix[xx]);
		}
		free(Cij_matrix);*/		
		Cij_matrix = 0;
	}
	if(maxnsam > maxCijnsam){
		/*in case matrix is smaller than data, not use*//*
		for(xx=0;xx<maxCijnsam+1;xx++) {
			for(yy=0;yy<xx+1;yy) {
				free(Cij_matrix[xx][yy]);
			}
			free(Cij_matrix[xx]);
		}
		free(Cij_matrix);*/
		Cij_matrix = 0;
	}
	
    /******************* calculate HKA ***********************/
	if(dohka_obs(matrix,matrixml,outgroup,data[0].nloci, &nlocihka, Cij_matrix) == 0) {
		fprintf(file_output,HKAdirect);
		fprintf(file_output,"OUTPUT FILE: date %s \n\nInput data from the file: %s\n\n",s,file_in);
        fprintf(file_output,"ERROR: Results are not available. Time obtained is negative.\n\n");
		return 1;    
	}
    /*print HKA results and probaility using chi-square*/
    /*output  file*/
	fprintf(file_output,HKAdirect);
    fprintf(file_output,"OUTPUT FILE: date %s \n\nInput data from the file: %s\n\n",s,file_in);
	
    /*print names and all values*/
    fprintf(file_output,"Title: %s ",data[0].title);
    fprintf(file_output,"\nnloci: %ld\n",data[0].nloci);
    if(data[0].factor_pool_an[0] == 0.)
		fprintf(file_output,"\n#IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\texpHKA_S\texpVar_S\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA");
    else {
		if(data[0].take_anbn == 1) 
			fprintf(file_output,"\n#IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\tfactor_an\tfactor_bn\texpHKA_S\texpHKA_div\texpHKA_theta\tpartialHKA");
		if(data[0].take_anbn == 2) 
			fprintf(file_output,"\n#IDloci\tnsam\tobs_S\tobs_div\tlength_pol\tlength_div\tfactor_chrn\tpmissing\texpHKA_S\texpVar_S0\texpVar_Smiss\texpHKA_div\texpVar_D\texpHKA_theta\tpartialHKA");
	}
	/*header in output_datag*/
    for(i=0;i<data[0].nloci;i++) {
        fprintf(file_output,"\n%s\t",data[0].idloci[i]);
        fprintf(file_output,"%ld\t",data[0].nsamples[i]);
        fprintf(file_output,"%ld\t",data[0].Ssites[i]);
        fprintf(file_output,"%.2f\t",data[0].divergence[i]);
        fprintf(file_output,"%.2f\t",data[0].length_pol[i]);
		fprintf(file_output,"%.2f\t",data[0].length_div[i]);
		fprintf(file_output,"%.2f\t",data[0].factor_chrn[i]);
		
		if(data[0].factor_pool_an[0] != 0.) {
			if(data[0].take_anbn == 1) {
				fprintf(file_output,"%.3f\t",data[0].factor_pool_an[i]);
				fprintf(file_output,"%.3f\t",data[0].factor_pool_bn[i]);
			}
			if(data[0].take_anbn == 2) {
				fprintf(file_output,"%.3f\t",data[0].pmissing[i]);
			}
		}
        fprintf(file_output,"%.2f\t",matrixml.Sexphka[i]);/*S*/
		if(data[0].take_anbn == 0) 
			fprintf(file_output,"%.2f\t",matrixml.varS0[i]);
		if(data[0].take_anbn == 2) {
			fprintf(file_output,"%.2f\t",matrixml.varS0[i]);
			fprintf(file_output,"%.2f\t",matrixml.varSmiss[i]);
		}
        fprintf(file_output,"%.2f\t",matrixml.Dexphka[i]);/*Div*/
		if(data[0].take_anbn == 0 || data[0].take_anbn == 2) fprintf(file_output,"%.2f\t",matrixml.varD[i]);
        fprintf(file_output,"%.5f\t",matrixml.hka_theta[i]);/*theta/nt*/
        fprintf(file_output,"%.3f\t",matrixml.hka[i]);/*partial HKA*/
    }
    fputs("\n\n",file_output);
    fprintf(file_output,"Time to the ancestor (in 2N generations): %.3f\t",matrixml.hka_T[0]);
    fprintf(file_output,"Chi-square value: %.3f\t",matrixml.Shka[0]);
    fputs("\n\n",file_output);
    print_prob_obshka(matrixml,nlocihka,file_output);
    
    /*time and duration*/
    time(&end);
    date = localtime(&end);
    strftime(s,80,"%c",date);
    fprintf(file_output,"\n\nDate of completion: %s",s);
    fprintf(file_output,"\nDuration of the process: %.3f seconds.",difftime(end,start));
    hour =  difftime(end,start)/3600.;
    min  = (difftime(end,start) - (float)hour*3600.)/60.;
    sec  =  difftime(end,start) - (float)hour*3600. - (float)min*60.;
    fprintf(file_output," (%dh:%dm:%ds)\n",hour,min,sec);
	
    fclose(file_output);
/*
	if(Cij_matrix != 0) {
		for(xx=0;xx<maxCijnsam+1;xx++) {
			for(yy=0;yy<xx+1;yy) {
				free(Cij_matrix[xx][yy]);
			}
			free(Cij_matrix[xx]);
		}
		free(Cij_matrix);		
	}
*/
    return(0);
}
