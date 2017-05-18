/*
 *  hka.c
 *  MuLoNeTests
 *
 *  Created by sonsins on Wed Mar 19 2003. 
 *  Thanks to J. Schumacher (BGC-MPI Jena, Germany) to explain me how to do the algorithm.
 *  Biometry functions included
 */

#include "hka.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "missing_freqs.h"

#define JCCorrection 0
#define NITER 100

int dohka_obs(struct statistics *matrix,struct statmulo matrixml,int outgroup, long int n_loci, long int *nlocihka, float ***Cij_matrix)
{
    long int i,j,k,n,ns,ns2;
	int *nx;
    double *S,r;
    double *D,*a,*theta,T,*b,*f,*Lp,*Ld,*theta2/*,*a2,*b2*/,an,bn,dbinom;
    int hka_par(double *,double *, double *,double *, double *, double *,long int, double **, double *);
	void init_coef(double *,int);
	long int rbinom(double,double,long int);
	float lnbinomial(int, int);
	double sum2bin;
	float varScoalnx(double, int, int *, long int , long int);
	float calculate_covarxy(long int,int *, int, float ***);
	float varS;
	float theta_rel;
	long int length_rel;
	float covarxy,sumanx;
	/*long int niter;*/
	
    if(outgroup == 1 && n_loci > 1) {
        matrixml.Shka[0] = (double)0.0;
        matrixml.hka_T[0] = (double)-10000;
    
        if((S = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((D = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((a = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((b = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((theta = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((theta2 = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((f = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((Lp = (double *) calloc(n_loci,sizeof(double))) == 0) {
            return 0;
        }
        if((Ld = (double *) calloc(n_loci,sizeof(double))) == 0) {
           return 0;
        }
		
		init_seed1((unsigned)time(NULL)/2);
		/*srand((unsigned)time(NULL)/2);*/
		ns2 = 0;
		
        for(i=0,k=0;i<n_loci;i++) {
            f[k] = matrix[i].factor_chrn;
			Lp[k] = matrix[i].length_pol;
			Ld[k] = matrix[i].length_div;
            /*with no correction... debugging.*/
			#if JCCorrection == 0
            S[k] = (double)matrix[i].biallsitesn;
            D[k] = matrix[i].ndivergence;
			#else
            /*or with correction of Jukes and Cantor.*/
            if(matrix[i].ndivergence/(double)matrix[i].length_div < (double)0.75) {
				if(matrix[i].ndivergence == (double)0) D[k] = (double)0;
				else D[k] = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * matrix[i].ndivergence/(double)matrix[i].length_div) * ((double)matrix[i].length_div + (double)matrix[i].nmhits);
                
				if(matrix[i].biallsitesn == 0) S[k] = (double)0;
				else { 
					/*S[k] = -(double)0.75 * (double)log((double)1. - (double)4./(double)3. * ((double)matrix[i].biallsitesn/coef[0])/(double)matrix[i].length_div) * (double)matrix[i].length_div * coef[0];*/
					S[k] = D[k]/matrix[i].ndivergence * (double)matrix[i].biallsitesn;
                }
            }
			#endif
            if(matrix[i].take_anbn == 0 || matrix[i].pmissing == 0.0) {
				a[k] = b[k] = (double)0.;
				for(j=1;j<matrix[i].nsamples;j++) {
					a[k] += (double)1./(double)j;
					b[k] += (double)1./((double)j*(double)j);
				}
			}
			else {
				if(matrix[i].take_anbn == 1) { 
					a[k] = matrix[i].factor_pool_an;
					b[k] = matrix[i].factor_pool_bn; /*INVALID*/
				}
				if(matrix[i].take_anbn == 2) {  /*assume the missing values are uniformly distributed across segsites*/
					/*include a calculation using a binomial distribution p missing values to have the number of samples per S*/
					/*CONDITIONAL ON BINOMIAL: the position has at least two samples !!!!*/
					/*
					a[k] = b[k] = (double)0.;
					for(j=0;j<S[k];j++) {
						do {
							r = rand()/(double)RAND_MAX;
							ns = rbinom(r,(1.0-matrix[i].pmissing),matrix[i].nsamples);
						} while(ns < 2);
						for(n=1;n<ns;n++) {
							a[k] += (double)1./(double)n;
							b[k] += (double)1./((double)n*(double)n);
						}
					}
					a[k] /= (double)S[k];
					b[k] /= (double)S[k];
					*/
					/*calculate the expected an and bn value conditioned on having at least two samples*/
					/*Assume uniform distribution of missing values, each position is a binomial*/
					/**/
					a[k] = b[k] = (double)0.;
					sum2bin = 0.;
					for(j=2;j<=matrix[i].nsamples;j++) 
						sum2bin += exp(lnbinomial((int)matrix[i].nsamples,(int)j)+(float)j*log((1.0-matrix[i].pmissing))+(float)(matrix[i].nsamples-j)*log(matrix[i].pmissing));
					for(j=2;j<=matrix[i].nsamples;j++) {
						an = bn = 0.0;					
						for(n=1;n<j;n++) {
							an += (double)1./(double)n;
							bn += (double)1./((double)n*(double)n);
						}
						dbinom = exp(lnbinomial((int)matrix[i].nsamples,(int)j)+(float)j*log((1.0-matrix[i].pmissing))+(float)(matrix[i].nsamples-j)*log(matrix[i].pmissing));
						a[k] += an * dbinom/sum2bin;
						b[k] += bn * dbinom/sum2bin; /*INVALID*/
					}
					/*
					printf("an[%02ld]=%f\n",k,a[k]);
					*/
				}
			}
			/*in case S and D data are zero or an is zero, not counting*/
            if((matrix[i].biallsitesn > 0 || matrix[i].ndivergence > (double)0) /**/&&
               (matrix[i].ndivergence/(double)matrix[i].length_div < (double)0.75)/**/ &&
			   a[k] > 0.0){
                k++;
            }
        }        
        *nlocihka = k;/*subtracting loci with no variation or an=0*/  
        
        if(*nlocihka) {
            if(!(hka_par(S,D,a,f,Lp,Ld,*nlocihka,&theta,&T))) {
				return 0;
			}
        }else {
			return 0;
        }
        for(i=0,k=0;i<n_loci;i++) {
            if((matrix[i].biallsitesn > 0 || matrix[i].ndivergence > (double)0) /**/&&
               (matrix[i].ndivergence/(double)matrix[i].length_div < (double)0.75)/**/){
				if(matrix[i].pmissing > 0.0) {
					if(Cij_matrix != 0) {
						/*calculate varS if missing*/
						/*to do "faster" we modify theta to have a fix number of positions L=100 (if theta is <= 0.1)*/
						if(Lp[k] < 100) 
							length_rel = Lp[k];
						else {
							if(theta[k]*f[k]*Lp[k] / 0.1 <= 100.0) length_rel = 100;
							else length_rel = (long int)floor(theta[k]*f[k]*Lp[k] / 0.1);
						}
						/*calculate by a binomial the number of samples at each position, conditioned on n>=2*/
						nx = (int *)calloc(length_rel,sizeof(int));
						sumanx = 0.;
						for(j=0;j<length_rel;j++) {
							do {
								r = ran1();
								ns = rbinom(r,(1.0-matrix[i].pmissing),matrix[i].nsamples);
							} while(ns < 2);
							nx[j] = (int)ns;
							an = 0;
							for(n=1;n<ns;n++) {
								an += (double)1./(double)ns;
							}
							sumanx += an;
						}
						/*todo*/
						theta_rel = theta[k]*f[k]*Lp[k]/(float)length_rel;
						covarxy = calculate_covarxy(length_rel,nx,matrix[i].nsamples,Cij_matrix);
						varS = theta_rel*f[k]*length_rel*sumanx + (S[k]*S[k] - S[k]) * covarxy / (sumanx * sumanx + covarxy);

						matrixml.varSmiss[i] = varS;
						matrixml.varS0[i] = theta[k]*Lp[k]*f[k]*a[k] + theta[k]*theta[k]*Lp[k]*f[k]*Lp[k]*f[k]*b[k];
						free(nx);
					}
					else {
						/*calculate the variance by simulation (it is easier...)*/
						/*to do "faster" we modify theta to have a fix number of positions L=100 (if theta is <= 0.1)*/
						if(Lp[k] < 100) 
							length_rel = Lp[k];
						else {
							if(theta[k]*f[k]*Lp[k] / 0.1 <= 100.0) length_rel = 100;
							else length_rel = (long int)floor(theta[k]*f[k]*Lp[k] / 0.1);
						}
						/*calculate by a binomial the number of samples at each position, conditioned on n>=2*/
						nx = (int *)calloc(length_rel,sizeof(int));
						for(j=0;j<length_rel;j++) {
							do {
								r = ran1();
								ns = rbinom(r,(1.0-matrix[i].pmissing),matrix[i].nsamples);
							} while(ns < 2);
							nx[j] = (int)ns;
						}
						theta_rel = theta[k]*f[k]*Lp[k]/(float)length_rel;
						varS = varScoalnx(theta_rel,matrix[i].nsamples,nx,length_rel,NITER);
						matrixml.varSmiss[i] = varS;
						matrixml.varS0[i] = theta[k]*Lp[k]*f[k]*a[k] + theta[k]*theta[k]*Lp[k]*f[k]*Lp[k]*f[k]*b[k];
						free(nx);
					}
				}
				else {
					varS = theta[k]*Lp[k]*f[k]*a[k] + theta[k]*theta[k]*Lp[k]*f[k]*Lp[k]*f[k]*b[k]; /*no missing*/
					matrixml.varS0[i] = varS;
				}
				matrixml.hka[i] = 
                    (S[k] - theta[k]*Lp[k]*f[k]*a[k]) * (S[k] - theta[k]*Lp[k]*f[k]*a[k]) / (varS) /*in case missing data the variance is modified*/
                    +
                    (D[k] - theta[k]*Ld[k]*(T+(double)1*f[k])) * (D[k] - theta[k]*Ld[k]*(T+(double)1*f[k])) / (theta[k]*Ld[k]*(T+(double)1*f[k]) + theta[k]*theta[k]*Ld[k]*Ld[k]*((double)1*f[k])*((double)1*f[k]));
				/*keep resuts and data*/
				matrixml.Sexphka[i] = theta[k]*Lp[k]*f[k]*a[k];
				matrixml.Dexphka[i] = theta[k]*Ld[k]*(T+(double)1*f[k]);
				matrixml.varD[i] = (theta[k]*Ld[k]*(T+(double)1*f[k]) + theta[k]*theta[k]*Ld[k]*Ld[k]*((double)1*f[k])*((double)1*f[k]));
				matrixml.Sobshka[i] = S[k];
				matrixml.Dobshka[i] = D[k];
				matrixml.Shka[0] += matrixml.hka[i];
				matrixml.hka_theta[i] = theta[k];
                k++;
            }
            /*in case S[i]+D[i] is zero, we decide to include in the analysis, and Chi-square at this locus is zero.*/
            else {
				matrixml.hka[i] = (double)0.;
				matrixml.hka_theta[i] = (double)0;
			}            
        }
        matrixml.hka_T[0] = T;
                
		free(S);
        free(D);
        free(a);
		free(b);
        free(theta);
        free(theta2);
        free(f);
		free(Lp);
		free(Ld);
        
        return 1;
    }
    else return 0;
    
    return 0;
}

int hka_par(double *S,double *D, double *a,double *f,double *Lp, double *Ld,long int nloci,double **theta,double *T)
{
    long int i;
    double x1,x2;
    double xacc;
    double functionT(double,double *,double *,double *, double *,double *, double *,long int);
    int zbrac(double *,double *,double *,double *,double *, double *,double *, double *,long int);
    double zriddr(double,double,double,double *,double *,double *, double *,double *, double *,long int);
    
    if(functionT(-0.5,S,D,a,f,Lp,Ld,nloci) < (double)0) { /*check if negative at T=-0.5*/
        /*calculate the rangs*/
        x1 = (double)0.;
        x2 = 10.;
        if(zbrac(&x1,&x2,S,D,a,f,Lp,Ld,nloci) == 1) {
            /*estimate the value of T*/
            xacc = (double)1e-6*x2; /* accuracy of five decimals*/
            *T = (double)zriddr(x1,x2,xacc,S,D,a,f,Lp,Ld,nloci);
            if(*T == -1.11e30) return 0;
            /*estimate all thetas/nt*/
            for(i=0;i<nloci;i++)
                theta[0][i] = (S[i] + D[i]) / (((*T+(double)1.0*f[i])*Ld[i]+a[i]*Lp[i]*f[i]));
        }
        else return 0;/*if not, that means that T is less than -0.5*/
    }
    else return 0;

    return 1;
}

double functionT(double T,double *S,double *D, double *a,double *f,double *Lp, double *Ld,long int nloci) 
{
	/*solution for T is obtained when the result is 0*/
	/*isolate theta_L for sum(Si) equation*/
	/*isolate theta_L for sum(Di) equation*/
	/*isolate theta_i for Si+Di equation*/
	/*equality between eq1 and eq2 and left only T by including eq3 in both*/
    long int i;
    double fdT;
    double s1,s2,s3,s4;
    
    s1 = s2 = s3 = s4 = (double)0.0;
    for(i=0;i<nloci-1;i++) {
        s1 += S[i]/(a[nloci-1]*Lp[nloci-1]*f[nloci-1]);
        s2 += (S[i] + D[i])/(a[i]*Lp[i]*f[i]+(T+(double)1*f[i])*Ld[i]) * a[i]*Lp[i]*f[i];
        s3 += D[i]/(T+(double)1*f[i]);
		s4 += (S[i] + D[i])/(a[i]*Lp[i]*f[i]+(T+(double)1*f[i])*Ld[i]) * Ld[i];
    }
    s1 += S[nloci-1]/(a[nloci-1]*Lp[nloci-1]*f[nloci-1]);
	s2  = s2/(a[nloci-1]*Lp[nloci-1]*f[nloci-1]);
	s3 += D[nloci-1]/(T+(double)1*f[nloci-1]);
	s3  = (s3 - s4)/(Ld[nloci-1]);
	
    fdT = s1 - s2 - s3;

    return fdT;
}

int zbrac(double *x1,double *x2,double *S,double *D, double *a,double *f,double *Lp, double *Ld,long int nloci)
{
	/* Based on Numerical Recipes in C. Press et al. 1992. 
    We need *x1 and *x2 be the range where a root is within them. We expand geometrically the range until finding
	(one a positive and one a negative value). If not, return 0.
	*/

    double f1,f2;
    int k=60;
    double functionT(double,double *,double *,double *,double *,double *, double *,long int);
    
    if(*x1 == *x2) return 0;

    f1 = functionT(*x1,S,D,a,f,Lp,Ld,nloci);
    f2 = functionT(*x2,S,D,a,f,Lp,Ld,nloci);
	
	if(f1*f2 < (double)0) return 1;

    while(k--) {
        if(fabs(f1) < fabs(f2)) {
            *x1 += (double)1.5 * (*x1 - *x2);
            f1 = functionT(*x1,S,D,a,f,Lp,Ld,nloci);
        }
        else {
            *x2 += (double)1.5 * (*x2 - *x1);
            f2 = functionT(*x2,S,D,a,f,Lp,Ld,nloci);
        }
        if(f1*f2 < (double)0) return 1;
    }
    return 0;
}

double zriddr(double xlow,double xhigh,double xacc,double *S,double *D, double *a,double *f,double *Lp, double *Ld,long int nloci)
{
	/* Based on Numerical Recipes in C. Press et al. 1992., p. 358 an on
	Ridders, 1979, IEEE Transactions on Circuits and systems, Vol. Cas-26, No. 11, pp. 979-980.
	*/
    int k=60;
	double flow,fhigh;
	double f1,f2,f3,f4;
	double x1,x2,x3,x4;
	double den,num,nsign;
    double functionT(double,double *,double *,double *, double *,double *, double *,long int);
	

    flow  = functionT(xlow,S,D,a,f,Lp,Ld,nloci);
    fhigh = functionT(xhigh,S,D,a,f,Lp,Ld,nloci);

	if(flow  == (double)0) return xlow;
	if(fhigh == (double)0) return xhigh;
	if(flow*fhigh > (double)0) 
		return (double)-1e32;
	
	x1 = xlow;
	x2 = xhigh;
	f1 = flow;
	f2 = fhigh;
		
	while(k--) {
		x3 = (x1+x2)/(double)2;
		f3 = functionT(x3,S,D,a,f,Lp,Ld,nloci);
		if(f1 - f2 < (double)0) nsign = (double)-1;
		else nsign = (double)1;
		num = (x3-x1) * f3 * nsign;
		den = (double)sqrt((double)f3*(double)f3 - (double)f1*(double)f2);
		if(den <= xacc && -den <= xacc) return x3;
		x4 = x3 + num/den;
		f4 = functionT(x4,S,D,a,f,Lp,Ld,nloci);
		if(f4 <= xacc && -f4 <= xacc) return x4;
		if(f3*f4<(double)0) {
			x1 = x3;
			f1 = f3;
			x2 = x4;
			f2 = f4;
		}
		else {
			if(f1*f4<(double)0) {
				x2 = x4;
				f2 = f4;
			}
			else {
				if(f2*f4<(double)0) {
					x1 = x4;
					f1 = f4;
				}
			}
		}
		if(fabs(x1-x2) <= xacc) return x1;
	}	
	return (double)-1e32;
}

void print_prob_obshka(struct statmulo matrixml, long int n_loci, FILE *file_output)
{

    double beta;
    double probQ_chisquare(int,double);
    
    beta = probQ_chisquare(n_loci-1,matrixml.Shka[0]);
    
    if(beta > (double)0.05) {
        if(!file_output) printf("Significance of Chi-square: Non-Significant, P(dgf=%ld) = %.3e\n",n_loci-1,beta);
        if(file_output) 
            fprintf(file_output,"Significance of Chi-square: Non-Significant, P(dgf=%ld) = %.3e\n",n_loci-1,beta);
    }
    else {
        if(!file_output) printf("Significance of Chi-square: Significant, P(dgf=%ld) = %.3e\n",n_loci-1,beta);
        if(file_output)
            fprintf(file_output,"Significance of Chi-square: Significant, P(dgf=%ld) = %.3e\n",n_loci-1,beta);
    }

    return;
}

float calculate_covarxy(long int L ,int *nx, int nsam,float ***Cij) {
	float cov;
	long int x,y,nxy;
	int xx,yy;
	float rhyper(int,int,int);
	
	cov = 0.0;
	for(x=0;x<L-1;x++) {
		for(y=x+1;y<L;y++) {
			nxy = rhyper(nx[x],nx[y],nsam);
			xx = max(nx[x],nx[y]);
			yy = min(nx[x],nx[y]);
			cov += Cij[xx][yy][nxy];
		}
	}
	
	return cov;
}


