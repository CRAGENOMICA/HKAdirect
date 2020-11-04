
#define HKAdirect "HKAdirect 0.70beta(20201104beta).\nSebastian E. Ramos-Onsins, Emanuelle Raineri and Luca Ferretti\n\n"

struct parametersg {
    char title[1024];
    /*observed data*/
    long int nloci;
    char **idloci;
    long int *nsamples;
    long int *Ssites;
    float *divergence;
	float *length_pol;
	float *length_div;
	float *factor_chrn;
	int take_anbn;
	float *factor_pool_an;
	float *factor_pool_bn;
	float *pmissing;
};

struct statistics /*ONE VECTOR FOR EACH LOCUS.*/
{
    long int nsamples;
    long int biallsitesn;
    float ndivergence;
    double length_pol;
	double length_div;
	float factor_chrn;
	int take_anbn;
	float factor_pool_an;
	float factor_pool_bn;
	float pmissing;

	int nmhits;
};

struct statmulo /*HKA results*/
{
    long int nloci;
   
    float *Shka;
    float *signif_hka;
    float *hka_theta;
    float *hka_T;
	float *varS0;
	float *varD;
	float *varSmiss;
    float *hka;
	
	float *Sexphka;
	float *Dexphka;
	float *Sobshka;
	float *Dobshka;
};


void init_seed1(long int seed);
double ran1(void);
