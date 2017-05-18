#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "missing_freqs.h"

#define PI 3.14159265359

/* TABLE CHI-SQUARE. UP TO 30 dgf. */

double table_chisquare(int dgf,double alfa) 
{
    /*table chi-square for 0.10, 0.05, 0.025, 0.01, 0.005 and 0.001. first 30 degrees of freedom. 1 tail*/
    double chsq[31][6] = { /*row[0] indicate the probabilities.*/
        {0.100,		0.050,	0.025,	0.010,	0.005, 	0.001},
        {2.706,		3.841,	5.024,	6.635,	7.879,	10.83},
        {4.605,		5.991,	7.378,	9.210,	10.597,	13.82},
        {6.251,		7.815,	9.348,	11.345,	12.838,	16.27}, 
        {7.779,		9.488,	11.143,	13.277,	14.860,	18.47}, 
        {9.236,		11.070,	12.833,	15.086,	16.750,	20.52}, 
        {10.645,	12.592,	14.449,	16.812,	18.548,	22.46}, 
        {12.017,	14.067,	16.013,	18.475,	20.278,	24.32}, 
        {13.362,	15.507,	17.535,	20.090,	21.955,	26.13}, 
        {14.684,	16.919,	19.023,	21.666,	23.589,	27.88}, 
        {15.987,	18.307,	20.483,	23.209,	25.188,	29.59}, 
        {17.275,	19.675,	21.920,	24.725,	26.757,	31.26}, 
        {18.549,	21.026,	23.337,	26.217,	28.300,	32.91}, 
        {19.812,	22.362,	24.736,	27.688,	29.819,	34.53}, 
        {21.064,	23.685,	26.119,	29.141,	31.319,	36.12}, 
        {22.307,	24.996,	27.488,	30.578,	32.801,	37.70}, 
        {23.542,	26.296,	28.845,	32.000,	34.267,	39.25}, 
        {24.769,	27.587,	30.191,	33.409,	35.718,	40.79}, 
        {25.989,	28.869,	31.526,	34.805,	37.156,	42.31}, 
        {27.204,	30.144,	32.852,	36.191,	38.582,	43.82}, 
        {28.412,	31.410,	34.170,	37.566,	39.997,	45.32}, 
        {29.615,	32.671,	35.479,	38.932,	41.401,	46.80}, 
        {30.813,	33.924,	36.781,	40.289,	42.796,	48.27}, 
        {32.007,	35.172,	38.076,	41.638,	44.181,	49.73}, 
        {33.196,	36.415,	39.364,	42.980,	45.559,	51.18}, 
        {34.382,	37.652,	40.646,	44.314,	46.928,	52.62}, 
        {35.563,	38.885,	41.923,	45.642,	48.290,	54.05}, 
        {36.741,	40.113,	43.195,	46.963,	49.645,	55.48}, 
        {37.916,	41.337,	44.461,	48.278,	50.993,	56.89}, 
        {39.087,	42.557,	45.722,	49.588,	52.336,	58.30}, 
        {40.256,	43.773,	46.979,	50.892,	53.672,	59.70},
    };

    if(dgf > 0 && dgf <= 30) {
        if(alfa == (double)chsq[0][0]) return (double)chsq[dgf][0];
        if(alfa == (double)chsq[0][1]) return (double)chsq[dgf][1];
        if(alfa == (double)chsq[0][2]) return (double)chsq[dgf][2];
        if(alfa == (double)chsq[0][3]) return (double)chsq[dgf][3];
        if(alfa == (double)chsq[0][4]) return (double)chsq[dgf][4];
        if(alfa == (double)chsq[0][5]) return (double)chsq[dgf][5];
        return -1.;
    }
    return -1.;
}

/* TABLE WILCOXON'S SIGNED RANK TEST. UP TO 50 VALUES. */

double table_wilcoxonsr(int n,double alfa) 
{
    /*table Wilcoxon's signed rank test for 0.10, 0.05, 0.02 and  0.01. first 50 n. 2 tails.*/
    double wsrt[51][4] = { /*row[0] indicate the probabilities.*/
        {(double)0.1,	(double)0.05,	(double)0.02,	(double)0.01},
        {1,	-1,	-1,	-1},
        {2,	1,	-1,	-1},
        {4,	2,	0,	-1},
        {6,	4,	2,	0},
        {8,	6,	3,	2},
        {11,	8,	5,	3},
        {14,	11,	7,	5},
        {17,	14,	10,	7},
        {21,	17,	13,	10},
        {26,	21,	16,	13},
        {30,	25,	20,	16},
        {36,	30,	24,	19},
        {41,	35,	28,	23},
        {47,	40,	33,	28},
        {54,	46,	38,	32},
        {60,	52,	43,	37},
        {68,	59,	49,	43},
        {75,	66,	56,	49},
        {83,	73,	62,	55},
        {92,	81,	69,	61},
        {101,	90,	77,	68},
        {110,	98,	85,	76},
        {120,	107,	93,	84},
        {130,	117,	102,	92},
        {141,	127,	111,	100},
        {152,	137,	120,	109},
        {163,	148,	130,	118},
        {175,	159,	141,	128},
        {188,	171,	151,	138},
        {201,	183,	162,	149},
        {214,	195,	174,	160},
        {228,	208,	186,	171},
        {242,	222,	198,	183},
        {256,	235,	211,	195},
        {271,	250,	224,	208},
        {287,	264,	238,	221},
        {303,	279,	252,	234},
        {319,	295,	267,	248},
        {336,	311,	281,	262},
        {353,	327,	297,	277},
        {371,	344,	313,	292},
        {389,	361,	329,	307},
        {408,	379,	345,	323},
        {427,	397,	362,	339},
        {446,	415,	380,	356},
        {466,	434,	398,	373},
    };

    if(n > 4 && n <= 50) {
        if(alfa == wsrt[0][0]) return wsrt[n-4][0];
        if(alfa == wsrt[0][1]) return wsrt[n-4][1];
        if(alfa == wsrt[0][2]) return wsrt[n-4][2];
        if(alfa == wsrt[0][3]) return wsrt[n-4][3];
        return (double)-1.;
    }
    return (double)-1.;
}

/* TABLE STANDARD NORMAL */

double table_standardnormal(double alfa) 
{
    double normalt[2][7] = {
        {0.10,	0.05,	0.025,	0.01,	0.05,	0.025,	0.001},
        {1.63,	1.96,	2.24,	2.58,	2.81,	3.04,	3.32},
    };
    
	if(alfa == (double)normalt[0][0]) return (double)normalt[1][0];
	if(alfa == (double)normalt[0][1]) return (double)normalt[1][1];
	if(alfa == (double)normalt[0][2]) return (double)normalt[1][2];
	if(alfa == (double)normalt[0][3]) return (double)normalt[1][3];
	if(alfa == (double)normalt[0][4]) return (double)normalt[1][4];
	if(alfa == (double)normalt[0][5]) return (double)normalt[1][5];
	if(alfa == (double)normalt[0][6]) return (double)normalt[1][6];
	return (double)-1.;
}

/*FUNCTION CHI TO NORMAL */

double chitonormal(int dgf,double chi)
{
    double z;
    if(dgf > 30) {
        z = ((double)pow(chi/(double)dgf,(double)1./(double)3.) - ((double)1.0 - (double)2.0/((double)9.0*(double)dgf)))/(double)sqrt((double)2.0/((double)9.0*(double)dgf));
        return z;
    }
    
    return (double)-1.;
}

/*GAMMA and Chi-square FUNCTIONS */

double gammln(double zz)
{
	/*Based on Numerical Recipes in C, Press et al. 1992. p. 213. and on 
	Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*gamma distribution for a z integer*/
	double loggammaz;
	double z,logg,h,sumc;
	static double gamma = 5.0;
	static double c0 =  1.000000000178;
	static double c1 = 76.180091729406;
	static double c2 = 86.505320327112;
	static double c3 = 24.014098222230;
	static double c4 =  1.231739516140;
	static double c5 =  0.001208580030;
	static double c6 =  0.000005363820;
	
	if(zz <= 0.) {
		puts("Error gamma");
		return (double)-10000.;
	}
	
	z = (double)zz;
	h = (double)sqrt(2. * PI);
	sumc = c0 + c1/(z+1.) - c2/(z+2.) + c3/(z+3.)  - c4/(z+4.) + c5/(z+5.) - c6/(z+6.);
	logg = (z + 0.5)*(double)log((double)(z + gamma + 0.5)) - (z + gamma + 0.5);
	loggammaz = log((double)h);
	loggammaz += logg + log((double)sumc);
	loggammaz -= log((double)z);
	
	return (double)loggammaz;
}

#define ITMAX 1000
#define EPS 3.0e-7
void gser(double *gamser,double a, double x, double *gln)
{
    /*Returns the incomplete gamma function P(a,x) evaluated by its series
    representation as gamser. Also returns ln T(a) as gln. numerical recipes in C. 2nd ed. p. 218
    */
    double gammln(double xx);
    int n;
    double sum,del,ap;
    
    *gln = gammln(a);
    if(x <= (double)0.0) {
        if(x < (double)0.0) puts("\n\nError in gser. Gamma is not well calculated.");
        *gamser = (double)0.0;
        return;
    }
    else {
        ap = a;
        del = sum = (double)1.0/a;
        for(n=1;n<=ITMAX;n++) {
            ++ap;
            del *= x/ap;
            sum += del;
            if(fabs(del) < fabs(sum)*EPS) {
                *gamser = sum*(double)exp(-x+a*(double)log(x)-(*gln));
                return;
            }
        }
        puts("\n\nError in gser. Gamma is not well calculated. ITMAX too small.");
		*gamser = (double)-10000;
        return;
    }
}

#define FPMIN 1.0e-30

void gcf(double *gammcf, double a, double x, double *gln)
{
    /*
        Returns the incomplete gamma function Q(a,x) evaluated by its continued fraction
        representation as gammcf. also returns ln T(a) as gln. numerical recipes in C. 2nd ed. p. 218
    */
    double gammln(double xx);
    int i;
    double an,b,c,d,del,h;
    
    *gln = gammln(a);
    b = x + (double)1.0 - a;
    c = (double)1.0/(double)FPMIN;
    d = (double)1.0/b;
    h = d;
    
    for(i=1;i<=ITMAX;i++) {
        an = - i*(i-a);
        b += (double)2.0;
        d = an*d+b;
        if(fabs(d) < FPMIN) d = (double)FPMIN;
        c = b + an/c;
        if(fabs(c) < FPMIN) c = (double)FPMIN;
        d = (double)1.0/d;
        del = d*c;
        h *= del;
        if(fabs(del - (double)1.0) < EPS) break;
    }
    if(i > ITMAX) {
		puts("Error in gcf. Gamma is not well calculated. ITMAX too small.");
		*gammcf = (double)-10000;
	}
    *gammcf = (double)exp(-x+a*(double)log(x)-(*gln))*h;
}

double gammq( double a, double x)
{
    /*Returns the incomplete gamma function Q(a,x) == 1 - P(a,x).*/
    /*numerical recipes in C. 2nd ed. p. 218*/
    void gcf(double *gammcf, double a, double x, double *gln);
    void gser(double *gamser, double a, double x, double *gln);
    double gamser,gammcf,gln;
    
    if(x < (double)0.0 || a <= (double)0.0) return (double)-1.0;
    if(x < (a+1.0)) {
        gser(&gamser,a,x,&gln);
        return (double)1.0-gamser;
    }
    else {
        gcf(&gammcf,a,x,&gln);
        return gammcf;
    }
}

double probQ_chisquare(int dgf,  double chisq)
{
    double gammq( double, double);
	/* probability of a chi square Q(a,x). Counting from the value of chi to infinite.*/
    return gammq((double)dgf/(double)2.0,chisq/(double)2.0);
}

/* CUMMULATIVE BINOMIAL DISTRIBUTION using the incomplete beta function. */

#define MAXIT 100

double betacf(double a, double b, double x)
{
    /*
    used by betai: Evaluates continued fraction for incomplete beta function by modified
    Lentz's method. numerical recipes in C. 2nd ed. p. 227.
    */
    int m,m2;
    double aa,c,d,del,h,qab,qam,qap;
    
    qab = a + b;
    qap = a + (double)1.0;
    qam = a - (double)1.0;
    c = (double)1.0;
    d = (double)1.0 - qab*x/qap;
    if(fabs(d) < FPMIN) d = (double)FPMIN;
    d = (double)1.0/d;
    h = d;
    for(m = 1;m <= MAXIT; m++) {
        m2 = 2*m;
        aa = m*(b-m)*x/((qam+m2)*(a+m2));
        d = (double)1.0 + aa*d;
        if((double)fabs(d) < FPMIN) d = (double)FPMIN;
        c = (double)1.0 + aa/c;
        if((double)fabs(c) < FPMIN) c = (double)FPMIN;
        d = (double)1.0/d;
        h *= d*c;
        aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
        d = (double)1.0 + aa*d;
        if((double)fabs(d) < FPMIN) d = (double)FPMIN;
        c = (double)1.0 + aa/c;
        if((double)fabs(c) < FPMIN) c = (double)FPMIN;
        d = (double)1.0/d;
        del = d*c;
        h *= del;
        if((double)fabs(del-(double)1.0) < EPS) break;
    }
    if(m > MAXIT)  {
		puts("Error in betacf. MAXIT is too small.");
		return (double)-10000;
	}
    return h;
}

double betai(double a, double b,double x)
{
    /*
    Returns the incomplete beta function Ix(a,b). numerical recipes in C. 2nd ed. p. 227.
    */
    double betacf(double a,double b,double x);
    double gammln(double xx);
    double bt;
    
    if(x < (double)0.0 || x > (double)1.0) {
        puts("Error in betai.");
        return (double)-1.0;
    }
    if(x == (double)0.0 || x == (double)1.0) bt = (double)0.0;
    else bt = (double)exp(gammln(a+b)-gammln(a)-gammln(b)+a*(double)log(x)+b*(double)log((double)1.0-x));
    
    if(x < (a + (double)1.0)/(a + b + (double)2.0)) return bt*betacf(a,b,x)/a;
    else return (double)1.0 - bt*betacf(b,a,(double)1.0-x)/b;
}

double prob_cumbinomial(int n, int k,double p)
{
    double prob;
	double betai(double,double,double);
	
	if(n==0) return -10000;
	/*modified to give allways values between 0 and 0.5*/
	if(n-k > k) k = n-k;
	/*Acumulated probability for a binomial prob p and n ind with k cases.*/
    prob = betai((double)k,(double)(n-k+(double)1.0),p);
	if(n-k != k) prob = (double)2.*prob;
	return prob;
}

/*modified from Emanuele Raineri:*//*
float lnbinomial(int n, int k)
{
	
	static int max_n=0;
	static float **table;
	double factln(long int n);
	int x;
	float result;
	
	if(max_n == 0) {
		max_n = 100;
		table = (float **)calloc(max_n,sizeof(float *));
		for(x=0;x<max_n;x++) {
			table[x] = (float *)calloc(max_n,sizeof(float));
		}
	}
	if(n+1 > max_n) {
		table = (float **)realloc(table,(n+1)*sizeof(float *));
		for(x=max_n;x<n+1;x++) {
			table[x] = (float *)calloc(n+1,sizeof(float));
		}
		for(x=0;x<max_n;x++) {
			table[x] = (float *)realloc(table[x],(n+1)*sizeof(float));
		}
		max_n = n+1;
	}
	if (n<0 || k<0)  {return -1e+6;}
	if (k>n) {return -1e+6;}
	result = (table[n][k]>0)?table[n][k]:(table[n][k]=factln((long int)n)-factln((long int)k)-factln((long int)(n-k)));
	return (result); 
}
*/
long int rbinom(double r,double p, long int ntot) 
{
	long int i,n;	
	float lnbinomial(int, int);
	static int max_n=0;
	static float **pnk;
	float sump;
	
	if(max_n==0) {
		max_n = ntot+1;
		pnk = (float **)calloc(max_n,sizeof(float *));
		for(i=0;i<max_n;i++) {
			pnk[i] = (float *)calloc(max_n,sizeof(float));
		}
	}
	if(ntot+1 > max_n) {
		pnk = (float **)realloc(pnk,(ntot+1)*sizeof(float *));
		for(i=max_n;i<ntot+1;i++) {
			pnk[i] = (float *)calloc(ntot+1,sizeof(float));
		}
		for(i=0;i<max_n;i++) {
			pnk[i] = (float *)realloc(pnk[i],(ntot+1)*sizeof(float));
		}
		max_n = ntot+1;
	}
	
	sump = 0.0;
	for(n=0;n<ntot;n++) {
		if(pnk[n][ntot] == 0.0) 
			pnk[n][ntot] = exp(lnbinomial((int)ntot,(int)n)+(float)n*log(p)+(float)(ntot-n)*log(1.0-p));
		if(pnk[n][ntot] >= 0 && pnk[n][ntot] < 1.0) {
			if(sump+pnk[n][ntot] >= r)
				break;
			sump += pnk[n][ntot];
		}
	}
	
	return n;
}

double factln(long int x)
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	 Lanczos 1964 J. SIAM Numer. Anal. Ser. B, Vol. 1 pp 86-96.*/
	
	/*G(n+1) = n!*/
	/*do the log(n!)*/
	double gammln(double);
	static double factlog[120];
	
	/*
	 if(x < 0) { 
	 puts("Error factln");
	 return (double)-10000.;
	 }
	 */
	if(x == 0) return 0.;
	
	if(x < 120) {
		if(factlog[x] == (double)0) {
			factlog[x] = gammln((double)x+(double)1.0);
			return factlog[x];
		}
		else return factlog[x];
	}
	return (gammln((double)x+(double)1.0));
}

double poissondist(double lambda) 
{
	/*Based on Atkinson 1979 Appl. Statist. 28: No. 1, pp, 29-35.*/
	double ran1(void);	
	double r,s;
	int N;
	
	double factln(long int);
	double alfa,beta,k,X;
	double rand1,rand2;
	static double c = (double)0.6;
	
	if(lambda < (double)0) {
		puts("Error poissondist");
		return (double)-10000.;
	}
	
	if(lambda == (double)0) return (double)0;
	if(lambda <= (double)20) {
		/*included for having not biased small values with mhits=0...*/
		/*if(lambda < 0.25) {
		 r = 1. - lambda;
		 s = ran1();
		 if(s >= r) N = (int)1;
		 else N = (int)0;
		 }
		 else {
		 */	r = (double)exp(-(double)lambda);
		N = (int)0;
		s = (double)1;
		do {
			s *= ran1();
			if(s >= r) N += (int)1;
			else break;
		}while(1);
		/*}*/
	}
	else {
		beta = (double)PI * (double)1./(double)sqrt((double)3*(double)lambda);
		alfa = beta * (double)lambda;
		k = (double)log((double)c) - (double)lambda - (double)log((double)beta);
		do{
			rand1 = ran1();
			X = (alfa-(double)log(((double)1-(double)rand1)/(double)rand1))/beta;
			if(X >= (double)-0.5) {
				N = (int)(X + (double)0.5);
				rand2 = ran1();
				if(alfa - beta*X +
				   (double)log((double)(rand2/(((double)1+(double)exp((double)(alfa-beta*X)))*((double)1+(double)exp((double)(alfa-beta*X)))))) 
				   <= k + (double)N*(double)log((double)lambda) - (double)factln((long int)N)) 
					break;
			}
		}while(1);
	}
	return (double)N;
}

double largebinomialdist(double pp, double n) 
{
	/*Based on Numerical Recipes in C, Press et al. 1992 and on
	 Fishman 1979 J. American Statistical Association Vol. 74, No. 366, pp 418-423*/
	
	double ran1(void);	
	double p,np;
	int N;
	double A,B,C,D,V,s;
	double g,plog,pclog,sq,angle,y,em,tt;
	double gammln(double);	
	
	if(pp > 0.5) p = (double)1.-pp;
	else p = pp;
	
	np = n * p;
	
	if(n==0) {
		puts("Error bindist");
		return (double)-10000.;
	}
	if(p==(double)0) {
		if(pp > 0.5) return (double)n;
		return (double)0;
	}
	
	if(np < (double)10) {
		/*Rejection Method: BI Algorithm*/
		s = (double)1- p;
		A = (double)1;
		B = p/s;
		C = ((double)n+(double)1)*B;
		D = A;
		N = 0;
		V = ran1()/(double)pow(s,(double)n);
		while(V > A) {
			N++;
			D *= (C/(double)N - B);
			A += D;
			if(N > n) break;
		}
	}
	else { /*Rejection method with a Lorentzian comparison distribution*/
		g = gammln(n+1.);
		plog  = (double)log(p);
		pclog = (double)log((1.0 - p));
		sq = (double)sqrt(2.0*np*(1.0 - p));
		do {
			do {
				angle = PI*ran1();
				y = tan(angle);
				em=sq*y+np;
			} while(em < 0.0 || em >= (n + 1.0));
			em = floor(em);
			tt = 1.2*sq*(1.0+y*y)*exp(g-gammln(em+1.0) - gammln(n-em+1.0)+em*plog+(n-em)*pclog);
		} while(ran1() > tt);
		N = (int)em;
	}
	
	if(pp > 0.5) N = (int)n - N;
	return (double)N;
}

int sample(long int nrep, long int len, int *data, int replace, float *prob, int *result)
{
	double ran1(void);
	long int i,j,k,lenr;
	float sump,sampl,sumsampl;
	
	if(replace==0 && len < nrep) {
		printf("\nError sampling,\n");
		return 0;
	}
	if(prob == 0 || data == 0 || result == 0 || len < 1 || nrep < 1 || (replace != 1 && replace != 0) ) {
		printf("\nError sampling,\n");
		return 0;
	}
	/*define sump*/
	sump = 0;
	for(i=0;i<len;i++) sump += prob[i];
		
	/*sample*/
	lenr = len;
	for(j=0;j<nrep;j++) {
		sampl = (long int)floor(ran1()*sump);
		sumsampl = 0.0;
		k=0; while(prob[k]==0.0) {sumsampl += prob[k]; k++;}
		while(sumsampl < sampl) {
			k++; sumsampl += prob[k];
		}
		/*result*/
		result[j] = data[k];
		if(replace==0) {
			sump -= prob[k];
			prob[k] = 0.0;
		}
	}
	
	return 1;
}

int unique(int *data, long int lendata, int *result)
{
	long int lenres,i,j;
	
	lenres = 0;
	for(i=0;i<lendata;i++) {
		for(j=i+1;j<lendata;j++) {
			if(data[i] == data[j]) break;
		}
		if(j==lendata) {
			result[lenres] = data[i];
			lenres++;
		}
	}
	return(lenres);
}

float li_mean(long int *S, long int len)
{
	long int i,sumS;
	float mean;
	sumS = 0;
	for(i=0;i<len;i++) sumS += S[i];
	mean = (float)sumS/(float)len;
	return(mean);
}
float li_var(long int *S, long int len)
{
	long int i;
	float meanS,varS;
	float li_mean(long int *, long int);
	
	meanS = li_mean(S,len);
	varS = 0.0;
	for(i=0;i<len;i++) 
		varS += (S[i]-meanS)*(S[i]-meanS);
	varS *= 1.0/((float)len-1);
	return(varS);
}

float rhyper(int nx,int ny,int nsam) {
	float lnbinomial(int, int);
	double ran1(void);
	int nxy;
	float pk,p;
	
	pk = ran1();
	p = 0.0;
	for(nxy=0;nxy<=min(nx,ny);nxy++) {
		p += exp(lnbinomial(nx,nxy) + lnbinomial(nsam-nx,ny-nxy) - lnbinomial(nsam,ny));
		if(p>pk) 
			break;
	}
	return(nxy);
}


