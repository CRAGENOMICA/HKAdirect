/*
 *  coaltrees.c
 *  HKAdirect
 *
 *  Created by Sebastian E. Ramos Onsins on 12/03/2013.
 *  Copyright 2013 CRAG. All rights reserved.
 *
 */

#include "coaltrees.h"

float calcTi(int i) 
{
	double ran1(void);
	float sTi;
	sTi = -log(1-ran1())/((float)i*(float)(i-1));
	return(sTi);
}

void dotree(int nt,int **node)
{
	int i,j,k;
	int res[2];
	float *prob;
	int *dat;
	long int len;
	int repl,nrep;

	dat  = (int *)calloc(nt,sizeof(int *));
	prob = (float *)calloc(nt,sizeof(float *));
	for(i=0;i<nt;i++) dat[i]  = i+1;
	for(i=nt;i>1;i--) for(j=0;j<i;j++) node[i-1][j] = 0;
	
	/*choose a couple of numbers that will coalesce. Less than the size-1 each coalescent*/
	/*the value indicate the connection of samples at the next step.*/
	for(i=nt;i>1;i--) {
		for(j=0;j<i-0;j++) prob[j] = 1.0;
		sample(nrep=1,len=i-0,dat,repl=0,prob,res+0);
		for(j=0;j<i-1;j++) prob[j] = 1.0;
		sample(nrep=1,len=i-1,dat,repl=0,prob,res+1);
		node[i-1][(int)res[0]-1] = (int)res[1];
		k = 1;
		for(j=0;j<i;j++) {
			if(node[i-1][j] == 0) {
				node[i-1][j] = k;
				k++;
			}
			/*printf("%d ",node[i-1][j]);*/
		}
		/*printf("\n");*/
	}
	free(dat);
	free(prob);
	return;
}

float sumTi_nx(int nt, int nxi, int **node, float *Ti)
{
	int i,j,lensam0,lensam1;
	float sumTi;
	int *res0,*res1;
	float *prob;
	
	prob = (float *)calloc(nt,sizeof(float *));
	for(i=0;i<nt;i++) prob[i] = 1.0;
	res0  = (int *)calloc(nxi,sizeof(int *));
	res1  = (int *)calloc(nxi,sizeof(int *));
	
	/*collect a suubsample nx from nt*/
	sample(nxi,nt,node[nt-1],0,prob,res0);
	sumTi = 0.0;
	i = nt-1;
	lensam0 = nxi;
	/*sum the number of branches in each nt step but for the subsample nx*/
	/*the coalescence is given in node*/
	while(lensam0 > 1) {
		sumTi += lensam0 * Ti[i];
		lensam1 = unique(res0,lensam0,res1);
		for(j=0;j<lensam1;j++) 
			res0[j] = node[i-1][res1[j]-1];
		lensam0 = lensam1;
		i--;
	}
	
	free(prob);
	free(res0);
	free(res1);
	return(sumTi);
}

long int calculateS(double theta, float sumTi) 
{
	long int S;
	double poissondist(double);
	
	S = poissondist(theta*sumTi);
	return(S);
}

float varScoalnx(double theta, int nt, int *nx, long int L, long int niter)
{
	float meanS,varS;
	float *Ti,sumTi;
	long int i,j,*Stot;
	int **node;
	
	Ti = (float *)calloc(nt,sizeof(float));
	node = (int **)calloc(nt,sizeof(int *));
	for(i=0;i<nt;i++)
		node[i] = (int *)calloc(nt,sizeof(int));
	Stot = (long int *)calloc(niter,sizeof(long int));
	
	for(i=0;i<niter;i++) {
		for(j=nt;j>=2;j--)
			Ti[j-1] = calcTi(j);
		dotree(nt,node);
		sumTi = 0.0;
		for(j=0;j<L;j++)
			sumTi += sumTi_nx(nt,nx[j],node,Ti);
		Stot[i] = calculateS(theta,sumTi);
	}
	meanS = li_mean(Stot,niter);
	varS  = li_var(Stot,niter);
	
	free(Ti);
	for(i=0;i<nt;i++) free(node[i]);
	free(node);
	free(Stot);
	return varS;
}

