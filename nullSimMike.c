/**************************************************************************/
/* This program simulates and analyzes datasets of inbred individuals.    */
/* It was written during Spring 2008*/
/* This is a modified version of simulation.c that includes the algorithm */
/* for null allele frequency estimation.*/
/* */
/* The null allele algorithm is slow, so I have cut some stuff out.       */
/* Specifically, I cut out CIs for all algorithms */
/* */
/* Also, I have added Vogl's algorithm. We run it twice. vogl2 denotes the */
/* results from a run with prior on the inbreeding coefficients of */
/* alpha = beta = 0.001, as in vogl's paper (seems to put a lot of prior */
/* weight on values of f near 0 or 1). */
/* vogl1 or just vogl denotes vogl run with a prior of alpha = beta = 1 */
/**************************************************************************/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define numNonXLinked 18
#define simnumber 0 /* I will use this to keep track of which set of */
                    /* simulations I am running */
//#define numnmark 6 /* the number of different values of nmark we will consider*/
#define numnmark 1 /* the number of different values of nmark we will consider*/
#define numnumallele 1 /* number of different values of numalleles we will use*/
#define nloops 400 /* number of sets of individuals to be generated. Each */
                  /* set of individuals includes one each of individuals with*/
                  /* f = 0, 0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, and 0.9*/
#define numfs 10 /* as in the previous line, we have 10 values of f*/
#define numnumind 1 /* number of different values for numind we will consider*/
/* from generate.c */
#define alleletype  9 /* 1 Dirichlet, 2 equally frequent, 3 common allele,  */
                     /* 4 set microsatellite values (must have 10 alleles) */
                     /* 5 biallelic, with MAF specified below */
                     /* 6 means to use a triangle distribution */
                     /* 7 means allele frequencies will be imported */
                     /* 8 means to use special frequencies (numalleles = 10)*/
                     /*   and have allele 9 be a null allele*/
                     /* 9 special frequencies: last allele is null, remaining*/
                     /*   probabilities are triangled out to t'other alleles*/
#define nullfreq1 0.00/* null allele frequencies are only used if alleletype*/
#define nullfreq2 0.10/* is 8*/
#define nullfreq3 0.20
#define MAF 0.4 /* minor allele frequency */
#define missrate 0.05 /* missing data rate */
#define linked 0 /* 0 for unlinked loci, 1 for linked loci, 2 for linked, different input file style*/

/* values having to do with EM algorithm */
//#define NMAX 100 /* maximum number of iterations allowed through EM algorithm*/
#define tol 1E-6		//tolerance level for EM for joint algorithm
#define tolerance 1E-6		/*tolerance level for betaEM */
#define nathantol 0.0001      //tolerance level for simple EM algorithm
#define liktol .1
#define ktol 1e-6 /* tolerance for Kalinowski and Taper's algorithm */
#define numstarts 1 /*number of different starting values for basic EM (hardcoded as 5, actually*/
#define nstart 1 /* number of different starting values for joint EM (should be, oh, 20 or more)*/
#define nullstarts 1 /* number of different starting values for null allele EM*/
#define M2nullstarts 1 /*number of different starting values for the null EM w/differing betas for the different markers.*/
#define kstarts 1 /* number of different starting values for null allele EM (Kalinowski and Taper)*/
#define MAXNUMTIMES 10000 /* maximum number of times to go through EM algorithm*/
#define alleletol 0.005
/* things to do with vogl's algorithm*/
#define burnin 100
#define numgibbs 100
#define voglbeta 0.001
#define voglalpha 0.001
int main()
{
  /* the next several functions include algorithms to estimate inbreeding */
  /* coefficients and a function that evaluates the likelihood of one of */
  /* the models.*/
  void EMalg(int*** genos, double** allelefreqs, int nmark, int numind, double* finblist);/* Nathan's algorithm */
  void finbreeding(int*** genos, double** allelefreq, int* numalleles, double* indf, int numind, int nmark);/* Daisy's algorithm--jointly estimates inbreeding coef and allefreqs*/
  double loglike(int*** genos, double** allelefreq, double* indf, int numind, int nmark);/* loglikelihood for Daisy's model*/
void nullEMalg(int*** myind1, double** allelefreqs, int* numalleles, double* finb,double* beta, int numind, int nmark);
  double nulllikelihood(int*** myind1,double** allelefreqs, double* beta, double* finb, int* numalleles, int numind, int nmark);
  void betaLaina(int*** myind1, double** allelefreqs, int* numalleles, double* finb, double* beta, int nmark, int numind);
  double betalikelihood(int*** myind1,double** allelefreqs, double* beta, double* finb, int* numalleles, int nmark, int numind);

  void Kalinowski(int*** myind1, double** kallelefreqs, int* numalleles, double* kbeta, int nmark, int numind);
  double klikelihood(int*** myind1,double** allelefreqs, double* beta, int* numalleles, int nmark, int numind);
  void simplecalc(int*** genos, double** allelefreqs, int nmark, int numind, double* simpfinblist, int* numalleles);
  void RitlandMME(int*** genos, double** allelefreqs, int nmark, int* numalleles, int numind, double* ritfinblist);
  /* these are the functions we use to set the starting values for the EM */
  /* algorithm for the joint estimation.*/
  void startallele(double **randallelefreq, int *numalleles, int nmark);
  void startindf(double *randindf, int numind);

  FILE *foutout, *foutlog;
  FILE *foutfreqs, *foutmarker;
  FILE *finalleles, *finmark;
  int generate(int toprint, int nmark, int* numalleles, int*** genos, double** truefreqs, double *recomb, double* finbtrue, int numind);
  /* functions that have to do with generating random allele frequencies,*/
  /* if applicable.*/
  double rbeta(double alpha, double beta);
  double rgamma(double alpha);
  double rexp(double lambda);
  void rdirichlet(double *alpha, double *p, int mynumalleles);
  void vogl(int*** genos, double** allelefreqs, double* voglfinblist, int* numalleles, int mynumind, int mynmark, double myalpha, double mybeta);/* Vogl's (2002) algorithm */
  /* names of files */
  char *logfile, *outputfile, *allelefile, *markerfile;
  /* index stuff: */
  int i, j, k, n, index;
  int markindex, nmarkindex, numalleleindex, numindindex, loopindex;/* which marker*/

  /*********************************************************************/
  /* having to do with assessing running time                          */
  /*********************************************************************/
  time_t t1, t2, t3;

  /* having to do with describing the families, individuals, and their genos:*/
  int ***genos;

  /*********************************************************************/
  /* having to do with describing the markers                          */
  /*********************************************************************/
  int nmarkarray[numnmark] = {757};
  //int nmarkarray[numnmark] = {100, 500};
  int numallelearray[numnumallele] = {3};
  int nmark;
  int maxnumalleles;
  int *numalleles;
  int *tempnumalleles;
  double* asum;
  double **allelefreqs;
  double **M2allelefreqs;
  double **kallelefreqs;
  double **voglfreqs;/* vogl: the allele frequencies [marker][allele] */
  double **vogl2freqs;/* vogl: the allele frequencies [marker][allele] */
  double **naivefreqs;
  //double **randallelefreq;        /* the random starting points for the algorithm */
  double **truefreqs;
  double **nullfreqs;
  double **alpha;/* Dirichlet parameters used for setting true allele freqs*/
                 /* if alleletype == 1*/
  double markdist, length, chromcov;/* distance in cM between markers,*/
  /* length in cM of chromosome, prop. of chromosome covered by our markers*/
  int nchrom; /* number of chromosomes*/
  int* nmarkchrom; /* number of markers per chromosome*/
  double* recomb; /* vector of inter-marker recombination frequencies*/


  /*********************************************************************/
  /* describing inbreeding.  */
  /*********************************************************************/

  int numind;
  //double* finb;/* a vector that stores the inbreeding coefficients for the */
  /* random single individuals in our simulations. */
  double *finbtrue;

  /*********************************************************************/
  /* miscellaneous stuff                                               */
  /*********************************************************************/
  double tempdoub, tempdoub1, tempdoub2; /* used if I need a temporary double for a couple of lines*/
  int count;
  int maxobsnamesize;/* maximum observed name size in the dataset */
  int toprint = 1; /* 0 if I do not want to print out a lot of extra */
                      /* information about how the MLE maximization is */
                      /* working. */
  double sum; /* Generally used for summing things (e.g. to check deltas */
              /* sum to 1, or to fine sum-squared difference between */
             /* deltaest and deltatrue.*/
  //double temp;/* not used?*/
  int tempint, tempint1, tempint2;/* temporary integers */
  /*********************************************************************/
  /* arrays to store results of different algorithms.                  */
  /*********************************************************************/
  double *finblist;	//stores convergence value for each individual
  double *nfinblist;	//stores convergence value for each individual
  double *simpfinblist; //stores simple finb values
  double *nsimpfinblist; //stores simple finb values
  double *modsimpfinblist;
  double *nmodsimpfinblist;
  double *ritfinblist;
  double *nritfinblist;
  double *modritfinblist;
  double *nmodritfinblist;
  double *jointfinblist;
  double *nullfinblist;
  double *M2finblist;
  double *beta;
  double *kbeta;
  double *M2beta;
  double *voglfinblist;
  double *vogl2finblist;
  /*********************************************************************/
  /* Dealing with different starting values.                  */
  /*********************************************************************/
 double logl, oldlogl;           /* stores calculated loglikelihood*/

  /*********************************************************************/
  /* General information about the markers.                            */
  /*********************************************************************/
 int* countmiss;
 int* counthet;
 int* counthom;


 (void) time(&t1);
 //srand48((long)t1);
 printf("\n not a random seed (168)");
 srand48(1254150884);

 maxnumalleles = numallelearray[numnumallele-1];

 logfile = (char*)malloc(14*sizeof(char));
 if(logfile == NULL)
   {
     printf("\n(181) Error allocating memory for logfilename.Exiting.\n");exit(0);
   }
 outputfile = (char*)malloc(14*sizeof(char));
 if(outputfile == NULL)
   {
     printf("\n(181) Error allocating memory for outputfilename.Exiting.\n");exit(0);
   }
 allelefile = (char*)malloc(15*sizeof(char));
 markerfile = (char*)malloc(15*sizeof(char));
 logfile = "simresults.log";
 allelefile = "simresults.freq";
 markerfile = "simresults.mark";
 outputfile = "simresults.out";

 if(simnumber == 0)
   {
     foutlog = fopen(logfile, "w");
     foutfreqs = fopen(allelefile, "w");
     foutmarker = fopen(markerfile, "w");
     foutout = fopen(outputfile, "w");
   }
 else{
   foutlog = fopen(logfile, "a");
   foutfreqs = fopen(allelefile, "a");
   foutmarker = fopen(markerfile, "a");
   foutout = fopen(outputfile, "a");
 }

 beta = (double*)malloc(1*sizeof(double));
 if(beta == NULL)
   {
     printf("\n(200) Error allocating memory for beta.Exiting.\n");exit(0);
   }
 fprintf(foutlog, "\n simulation %d; seed = %d", simnumber, t1);
 fprintf(foutlog, "\n\t linked = %d;", linked);
 fprintf(foutlog, "\n\t alleletype = %d; missrate = %4.4f", alleletype, missrate);
 fprintf(foutlog, "\n\t tol = %4.4e; nathantol = %4.4f; liktol = %4.4e", tol, nathantol, liktol);
 fflush(foutlog);
 /* beginning simulation loop */
 for(numindindex = 0; numindindex < numnumind; numindindex++)
   {
     if(numnumind == 4)
       {
	 if(numindindex == 0){numind = 20;}
	 else if(numindindex == 1){numind = 50;}
	 else if(numindindex == 2){numind = 100;}
	 else if(numindindex == 3){numind = 200;}
       }
     else if(numnumind == 1)
       {
	 numind = 94;
       }
	 else
	   {
	     printf("\n(196) Actually, this part of the program is hardcoded to have numnumind == 4.\n");exit(0);
	   }

     finbtrue = (double*)malloc(numind*sizeof(double));

     for(i = 0; i < (int)numind/numfs; i++)
       {
	 finbtrue[i*numfs+0] = 0.0;
	 finbtrue[i*numfs+1] = 0.02;
	 finbtrue[i*numfs+2] = 0.05;
	 finbtrue[i*numfs+3] = 0.1;
	 finbtrue[i*numfs+4] = 0.2;
	 finbtrue[i*numfs+5] = 0.3;
	 finbtrue[i*numfs+6] = 0.4;
	 finbtrue[i*numfs+7] = 0.5;
	 finbtrue[i*numfs+8] = 0.7;
	 finbtrue[i*numfs+9] = 0.9;
       }
     /*****************************************/
     /* Allocating Memory for results of various algorithms.*/
     /*****************************************/
     simpfinblist=(double*)malloc(numind*sizeof(double));
     modsimpfinblist=(double*)malloc(numind*sizeof(double));
     ritfinblist=(double*)malloc(numind*sizeof(double));
     modritfinblist=(double*)malloc(numind*sizeof(double));
     finblist=(double*)malloc(numind*sizeof(double));
     jointfinblist=(double*)malloc(numind*sizeof(double));
     nullfinblist=(double*)malloc(numind*sizeof(double));
     M2finblist=(double*)malloc(numind*sizeof(double));
     voglfinblist=(double*)malloc(numind*sizeof(double));
     vogl2finblist=(double*)malloc(numind*sizeof(double));
     /* next, the same stuff, but new ones to indicate that they used naive allelefreqs*/
     nsimpfinblist=(double*)malloc(numind*sizeof(double));
     nmodsimpfinblist=(double*)malloc(numind*sizeof(double));
     nritfinblist=(double*)malloc(numind*sizeof(double));
     nmodritfinblist=(double*)malloc(numind*sizeof(double));
     nfinblist=(double*)malloc(numind*sizeof(double));

     /* here we start looping*/
     if(simnumber == 0)
       {
	 fprintf(foutout, "simset simnumber nmark numalleles numind truef MLE nMLE Simple nSimple ModSimple nModSimple Ritland nRitland ModRitland nModRitland jMLE M1MLE M2MLE Vogl1 Vogl2\n");
	 fprintf(foutmarker, "simset simnumber nmark numalleles numind missing hom het betatrue M1beta M2beta KTbeta nulltrue M1pnull M2pnull KTpnull V1pnull V2pnull\n");
	 fprintf(foutfreqs, "simset simnumber nmark numalleles numind nulltrue L1pnull L2pnull kpnull V1pnull V2pnull barrier1");
	 for(i = 1; i < maxnumalleles; i++)
	   {
	     fprintf(foutfreqs, " jointfreq%d", i);
	   }
	 fprintf(foutfreqs, " barrier2");
	 for(i = 1; i < maxnumalleles; i++)
	   {
	     fprintf(foutfreqs, " M1freq%d", i);
	   }
	 fprintf(foutfreqs, " barrier3");
	 for(i = 1; i < maxnumalleles; i++)
	   {
	     fprintf(foutfreqs, " M2freq%d", i);
	   }
	 fprintf(foutfreqs, " barrier4");
	 for(i = 1; i < maxnumalleles; i++)
	   {
	     fprintf(foutfreqs, " kfreq%d", i);
	   }
	 fprintf(foutfreqs, " barrier5");
	 for(i = 1; i < maxnumalleles; i++)
	   {
	     fprintf(foutfreqs, " V1freq%d", i);
	   }
	 fprintf(foutfreqs, " barrier6");
	 for(i = 1; i < maxnumalleles; i++)
	   {
	     fprintf(foutfreqs, " V2freq%d", i);
	   }
	 fprintf(foutfreqs, "\n");
   }

     for(nmarkindex = 0; nmarkindex < numnmark; nmarkindex++)
       {
	 nmark = nmarkarray[nmarkindex];
	 kbeta = (double*)malloc(nmark*sizeof(double));
	 if(kbeta == NULL)
	   {
	     printf("\n(201) Error allocating memory for kbeta.Exiting.\n");exit(0);
	   }
	 M2beta = (double*)malloc(nmark*sizeof(double));
	 if(M2beta == NULL)
	   {
	     printf("\n(306) Error allocating memory for M2beta.Exiting.\n");exit(0);
	   }
		   countmiss = (int*)malloc(nmark*sizeof(int));
		   counthom = (int*)malloc(nmark*sizeof(int));
		   counthet = (int*)malloc(nmark*sizeof(int));
	 for(numalleleindex = 0; numalleleindex < numnumallele; numalleleindex++)
	   {
	     printf("\n(246)");fflush(stdout);
	     numalleles = (int*)malloc(nmark*sizeof(int));
	     tempnumalleles = (int*)malloc(nmark*sizeof(int));

	     for(markindex = 0; markindex < nmark; markindex++)
	       {
		 numalleles[markindex] = numallelearray[numalleleindex];
	       }

	     /* allocating memory and initializing */
	     asum = (double*)malloc(numallelearray[numalleleindex]*sizeof(double));
	     /* The allele frequency arrays have dimensions: 3Xmaxnumalleles */
	     /* The number of rows is the number of different null allele */
	     /* frequencies we are looking at. Then, we need a column for */
	     /* each allele. We went with maxnumalleles instead of the exact*/
	     /* number of alleles for these markers so that we could print*/
	     /* exactly the same regardless of the value of numalleles (and*/
	     /* so the output file looks like a table that can be imported */
	     /* into R.*/
	     genos = (int***)malloc(numind*sizeof(int**));
	     for(j = 0; j < numind; j++)
	       {
		 genos[j] = (int**)malloc(nmark*sizeof(int*));
		 for(k = 0; k < nmark; k++)
		   {
		     genos[j][k] = (int*)malloc(2*sizeof(int));
		     genos[j][k][0] = 0;
		     genos[j][k][1] = 0;
		   }
	       }
	     allelefreqs = (double**)malloc(nmark*sizeof(double*));
	     M2allelefreqs = (double**)malloc(nmark*sizeof(double*));
	     kallelefreqs = (double**)malloc(nmark*sizeof(double*));
	     voglfreqs = (double**)malloc(nmark*sizeof(double*));
	     vogl2freqs = (double**)malloc(nmark*sizeof(double*));
	     naivefreqs = (double**)malloc(nmark*sizeof(double*));
	     //randallelefreq = (double**)malloc(nmark*sizeof(double*));
	     truefreqs = (double**)malloc(nmark*sizeof(double*));
	     nullfreqs = (double**)malloc(nmark*sizeof(double*));
	     for(i = 0; i < nmark; i++)
	       {
		 allelefreqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 M2allelefreqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 kallelefreqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 voglfreqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 vogl2freqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 naivefreqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 //randallelefreq[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 truefreqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
		 nullfreqs[i] = (double*)malloc(numalleles[i]*sizeof(double));
	       }

	     for(loopindex = 0; loopindex < nloops; loopindex++)
	       {
		 //printf("\n(287, loopindex = %d)", loopindex);fflush(stdout);
		 /* set allele frequencies */
		 if(alleletype == 1)
		   {
		     alpha = (double**)malloc(nmark*sizeof(double*));
		     for(markindex = 0; markindex < nmark; markindex++)
		       {
			 alpha[markindex] = (double*)malloc(numalleles[markindex]*sizeof(double));
		       }
		     //if(toprint == 1)printf("\n Dirichlet allele frequencies.");
		     for(markindex = 0; markindex < nmark; markindex++)
		       {
			 for(i = 0; i < numalleles[markindex]; i++)
			   {
			     alpha[markindex][i] = 1;/* 1 is the parameter used in Milligan's paper */
			   }

			 rdirichlet(alpha[markindex], truefreqs[markindex], numalleles[markindex]);
		       }
		     for(markindex = 0; markindex < nmark; markindex++){free(alpha[markindex]);}
		     free(alpha);
		   }
		 else if(alleletype == 2)
		   {
		     for(markindex = 0; markindex < nmark; markindex++)
		       {
			 for(i = 0; i < numalleles[markindex]; i++)
			   {
			     truefreqs[markindex][i] = (double) 1.0/numalleles[markindex];
			   }
		       }
		   }
		 else if(alleletype == 3)
		   {
		     if(toprint == 1)printf("\n One common allele");
		     for(markindex = 0; markindex < nmark; markindex++)
		       {
			 truefreqs[markindex][0] = 0.8;
			 for(i = 1; i < numalleles[markindex]; i++)
			   {
			     truefreqs[markindex][i] = (double) 0.2/(numalleles[markindex]-1);
			   }
		       }
		   }
		 else if(alleletype == 4)
		   {
		     for(markindex = 0; markindex < nmark; markindex++)
		       {
			 if(numalleles[0] != 10)
			   {
			     printf("\n Set allele frequencies are only available for the case in which there are ten alleles. Please reset numalleles and start again.\n");
			     exit(0);
			   }
		       }
		     for(j = 0; j < nmark; j++)
		       {
			 truefreqs[j][0] = 0.34;
			 truefreqs[j][1] = 0.24;
			 truefreqs[j][2] = 0.16;
			 truefreqs[j][3] = 0.11;
			 truefreqs[j][4] = 0.06;
			 truefreqs[j][5] = 0.04;
			 truefreqs[j][6] = 0.02;
			 truefreqs[j][7] = 0.01;
			 truefreqs[j][8] = 0.01;
			 truefreqs[j][9] = 0.01;
		       }
		   }
		 else if(alleletype == 5)
		   {
		     for(markindex = 0; markindex < nmark; markindex++)
		       {
			 if(numalleles[markindex] != 2)
			   {
			     printf("\n alleletype 5 requires that the markers be diallelic. Please reset either your alleletype or your number of alleles.\n");
			     exit(0);
			   }
		       }
		     for(j = 0; j < nmark; j++)
		       {
			 truefreqs[j][0] = MAF;
			 truefreqs[j][1] = 1.0 - MAF;
		       }
		   }
		 else if(alleletype == 6)/* triangle distribution */
		   {
		     for(j = 0;j < nmark; j++)
		       {
			 sum = 0.0;
			 for(i = 1; i <= numalleles[j]; i++)
			   {
			     sum += i;
			   }

			 for(i = 0; i < numalleles[j]; i++)
			   {
			     truefreqs[j][i] = (numalleles[j] - i)/sum;
			   }
		       }
		   }
		 else if(alleletype == 7)
		   {  /* next,  read in allele frequencies from an input file */
		     finalleles = fopen("allelefile.dat", "r");
		     if(finalleles == NULL)
		       {
			 printf("\nUnable to locate input file (allelefile.dat)\n");

			 exit(1);
		       }

		     printf("\n(316) I'm cheating here--the number of alleles per marker are not being read in from the input file.\n");
		     for(markindex = 0; markindex < nmark; markindex++)
		       {
			 fscanf(finalleles, " %*d");
			 for(i = 0; i < numalleles[markindex]; i++)
			   {
			     fscanf(finalleles, " %lf", &truefreqs[markindex][i]);
			   }
		       }
		     fclose(finalleles);
		   }
		 else if(alleletype == 8)
		   {
		     if(numalleles[0] <= 2)
		       {
			 printf("\n (line 167) You might want to rethink running the program with numalleles = %d. With one allele null, that doesn't give you observable variation at this marker.\n");
			 exit(0);
		       }

		     i = 0;
		     for(j = 0; j < (int)nmark/3+1; j++)
		       {
			 truefreqs[i][numalleles[i]-1] = nullfreq1;
			 for(k = 0; k < numalleles[i]-1; k++)
			   {
			     truefreqs[i][k] = (1.0 - nullfreq1)/(numalleles[i]-1);
			   }
			 i++;
			 if(i < nmark)
			   {
			     truefreqs[i][numalleles[i]-1] = nullfreq2;
			     for(k = 0; k < numalleles[i]-1; k++)
			       {
				 truefreqs[i][k] = (1.0 - nullfreq2)/(numalleles[i]-1);
			       }
			     i++;
			   }
			 if(i < nmark)
			   {
			     truefreqs[i][numalleles[i]-1] = nullfreq3;
			     for(k = 0; k < numalleles[i]-1; k++)
			       {
				 truefreqs[i][k] = (1.0 - nullfreq3)/(numalleles[i]-1);
			       }
			   }
			 i++;
		       }
		   }
		 else if(alleletype == 9)
		   {
		     if(numalleles[0] <= 2)
		       {
			 printf("\n (line 439) You might want to rethink running the program with numalleles = %d. With one allele null, that doesn't give you observable variation at this marker.\n");
			 exit(0);
		       }

		     i = 0;
		     for(j = 0; j < (int)nmark/3+1; j++)
		       {
			 truefreqs[i][numalleles[i]-1] = nullfreq1;
			 for(k = 0; k < numalleles[i]-1; k++)
			   {
			     truefreqs[i][k] = 2*(k+1)*(1.0 - nullfreq1)/((numalleles[i])*(numalleles[i]-1));
			   }
			 i++;
			 if(i < nmark)
			   {
			     truefreqs[i][numalleles[i]-1] = nullfreq2;
			     for(k = 0; k < numalleles[i]-1; k++)
			       {
				 truefreqs[i][k] = 2*(k+1)*(1.0 - nullfreq2)/((numalleles[i])*(numalleles[i]-1));
			       }
			     i++;
			   }
			 if(i < nmark)
			   {
			     truefreqs[i][numalleles[i]-1] = nullfreq3;
			     for(k = 0; k < numalleles[i]-1; k++)
			       {
				 truefreqs[i][k] = 2*(k+1)*(1.0 - nullfreq3)/((numalleles[i])*(numalleles[i]-1));
			       }
			   }
		     i++;
		       }
		   }
		 else
		   {
		     printf("\n Undefined allele type. Exiting. ");
		     exit(1);
		   }

		 /* done reading in allele frequencies*/

		 /*for(markindex = 0; markindex < nmark; markindex++)
		   {
		   printf("\n %d", markindex);
		   for(j = 0; j < numalleles[markindex]; j++)
		   {
		   printf(" %4.5f", truefreqs[markindex][j]);
		   }
		   }
		   printf("\n(482)");exit(0);*/
		 /* now, one way or the other, we have the allele frequencies*/

		 /* Next, our marker spacing */
		 /*if(linked == 0){printf("\n Simulating %d unlinked loci.", nmark);}
		   else{printf("\n Simulating %d linked loci.", nmark);}
		   fflush(stdout);*/

		 recomb = (double*)malloc((nmark-1)*sizeof(double));
		 if(linked == 0)
		   {/* unlinked*/
		     for(markindex = 0; markindex < nmark-1; markindex++)
		       {
			 recomb[markindex] = 0.5;
		       }
		   }
		 else if(linked == 1)
		   {
		     finmark = fopen("marker.dat", "r");
		     if(finmark == NULL)
		       {
			 printf("\nUnable to locate input file (markfile.dat)\n");
			 exit(1);
		       }
		     fscanf(finmark, " nchrom %d", &nchrom);
		     nmarkchrom = (int*)malloc(nchrom*sizeof(int));
		     for(i = 0; i < nchrom; i++)
		       {
			 fscanf(finmark, " %d", &nmarkchrom[i]);
		       }

		     markindex = 0;
		     for(i = 0; i < nchrom; i++)
		       {
			 fscanf(finmark, " %*d %*d %lf %*f %lf", &length, &chromcov);
			 /* That length is in cM. We divide by 100 to get Morgans.*/
		  /* We want our markers to cover the same distance as the */
		  /* mouse markers: namely chromcov*lengthofchrom*/
		  /* finally, we divide that length evenly among our markers */
		  markdist = length*chromcov/(100*nmarkchrom[i]);
			  /* convert to recombination frequency using Haldane's map funct.*/
		    tempdoub = 0.5*(1-exp(-2*markdist));
			  for(j = 0; j < nmarkchrom[i]-1; j++)
			    {/* use Haldane's map function to go from distance to recomb*/
			      recomb[markindex] = tempdoub;
			      markindex++;
			    }
			  if(i < nchrom-1)
			    {
			      recomb[markindex] = 0.5;/* last marker on this chromosome is*/
			      /* unlinked to first marker on next chromosome*/
			      markindex++;
			    }
			}
		      fclose(finmark);
		    }
		  else if(linked == 2)
		    {/* different style of input file than if linked ==1*/
		      finmark = fopen("marker.dat", "r");
		      fscanf(finmark, " nchrom %d", &nchrom);
		      nmarkchrom = (int*)malloc(nchrom*sizeof(int));
		      for(i = 0; i < nchrom; i++)
			{
			  fscanf(finmark, " %d", &nmarkchrom[i]);
			}
		      fscanf(finmark, " chr pos.bp marker.id pos.cM");
		      markindex = 0;
		      for(i = 0; i < nchrom; i++)
			{
			  fscanf(finmark, "%d %*d rs%*d %lf", &tempint1, &tempdoub1);
			  if(i > 0)
			    {
			      recomb[markindex] = 0.5; /* first marker on any chrom has a */
			      /* recomb freq of 0.5 with last marker on previous chrom.*/
			      markindex++;
			    }
			  if(tempint1 != (i+1))
			    {
			      printf("\n something wrong with input marker file.\n");exit(0);
			    }
			  for(j = 1; j < nmarkchrom[i]; j++)
			    {
			      fscanf(finmark, "%*d %*d rs%*d %lf", &tempdoub2);
			      markdist = (tempdoub2 - tempdoub1)/100;/* divide by 100 since these are measured in cM*/
			      tempdoub = 0.5*(1-exp(-2*markdist));/* Haldane*/
			      recomb[markindex] = tempdoub;
			      markindex++;
			      tempdoub1 = tempdoub2;
			    }
			}
		      fclose(finmark);
		    }
		  /*for(markindex = 0; markindex < nmark-1; markindex++)
		    {
		    printf("\n %3d %3d: %4.4e", markindex+1, markindex+2, recomb[markindex]);
		    }*/
		  /****************************************************************/
		  /* At this point, all preliminary stuff has been taken care of. */
		  /* We now begin running the simulations.                        */
		  /****************************************************************/
		  /* now that we have the allele frequencies, we can generate the data */
		  /* call generate here */

		  tempint = generate(toprint, nmark, numalleles, genos, truefreqs, recomb, finbtrue, numind);

		  /* truefreqs are the true allele frequencies. naivefreqs are the*/
		  /* naively estimated allele freqs.*/

		  for(markindex = 0; markindex < nmark; markindex++)
		    {

		      tempint = 0;/* in this loop, tempint keeps track of how many*/
		      /* non-missing alleles we have seen */
		      for(j = 0; j < numalleles[markindex]; j++)
			{
			  asum[j] = 0;
			}
		      for(i = 0; i < numind; i++)
			{
			  if(genos[i][markindex][0] != 99)
			    {
			      asum[genos[i][markindex][0]]++;
			      asum[genos[i][markindex][1]]++;
			      tempint = tempint+2;
			    }
			}

		      for(j = 0; j < numalleles[markindex]; j++)
			{
			  naivefreqs[markindex][j] = asum[j]/(tempint);
			}

		    }/* end for(markindex...)*/

		  /*************************************************************/
		  /* Now the data is done being generated.                   */
		  /*************************************************************/
		  /* Note, we are making the algorithms just analyze the first numfs */
		  /* individuals--one ind. for each f-value for each dataset */

		  printf("\n(725) starting new dataset. simnumber = %d, nmark = %d, numalleles = %d, numind = %d.", loopindex, nmark, numalleles[0], numind);fflush(stdout);
		  fflush(foutlog);fflush(foutfreqs);fflush(foutmarker);fflush(foutout);
		  EMalg(genos, truefreqs, nmark, numfs, finblist);
		  EMalg(genos, naivefreqs, nmark, numfs, nfinblist);

		  simplecalc(genos, truefreqs, nmark, numfs, simpfinblist, numalleles);
		  simplecalc(genos, naivefreqs, nmark, numfs, nsimpfinblist, numalleles);

		  for(i = 0; i < numfs; i++)
		    {
		      if((simpfinblist[i] > 0)&&(simpfinblist[i] < 1.0))
			{
			  modsimpfinblist[i] = simpfinblist[i];
			}
		      else if(simpfinblist[i] <= 0){modsimpfinblist[i] = 0.0;}
		      else if(simpfinblist[i] >= 1){modsimpfinblist[i] = 1.0;}

		      if((nsimpfinblist[i] > 0)&&(nsimpfinblist[i] < 1.0))
			{
			  nmodsimpfinblist[i] = nsimpfinblist[i];
			}
		      else if(nsimpfinblist[i] <= 0){nmodsimpfinblist[i] = 0.0;}
		      else if(nsimpfinblist[i] >= 1){nmodsimpfinblist[i] = 1.0;}
		    }

		  RitlandMME(genos, truefreqs, nmark, numalleles, numfs, ritfinblist);
		  RitlandMME(genos, naivefreqs, nmark, numalleles, numfs, nritfinblist);

		  for(i = 0; i < numfs; i++)
		    {
		      if((ritfinblist[i] > 0)&&(ritfinblist[i] < 1.0))
			{
			  modritfinblist[i] = ritfinblist[i];
			}
		      else if(ritfinblist[i] <= 0){modritfinblist[i] = 0.0;}
		      else if(ritfinblist[i] >= 1){modritfinblist[i] = 1.0;}
		    }

		  for(i = 0; i < numfs; i++)
		    {
		      if((nritfinblist[i] > 0)&&(nritfinblist[i] < 1.0))
			{
			  nmodritfinblist[i] = nritfinblist[i];
			}
		      else if(nritfinblist[i] <= 0){nmodritfinblist[i] = 0.0;}
		      else if(nritfinblist[i] >= 1){nmodritfinblist[i] = 1.0;}
		    }
		  /* Starting Daisy's algorithm for joint estimation of inbreeding*/
		  /* and allele frequencies. (579)*/

		  for(n=0;n<nstart;n++){

		    for(j = 0; j < nmark; j++)
		      {
			tempnumalleles[j] = numalleles[j] - 1;
		      }
		    startallele(allelefreqs, tempnumalleles, nmark);
		    startindf(jointfinblist, numind);

		    finbreeding(genos, allelefreqs, numalleles, jointfinblist, numind, nmark);

		    if(n==0){logl= loglike(genos, allelefreqs, jointfinblist, numind, nmark);}// if n == 0
		    /**************************************************************************/
		    /*This entire long section basically tests to see if we've ended up in a new place, and if we have then we print*/
		    /**************************************************************************/

		    if(n>0){
		      oldlogl = logl;
		      logl = loglike(genos, allelefreqs, jointfinblist, numind, nmark);
		      if(fabs(oldlogl-logl)>1000*tol){
			printf("\n(671) Joint algorithm: The log likelihoods differ\n  Old: %f new: %f\n", oldlogl, logl);
			  if(nmark < 20)
			    {
			      printf("\n\t");
			      for(i = 0; i < nmark; i++)
				{
				  printf(" %d%d", genos[j][i][0], genos[j][i][1]);
				}
			      printf("\n");
			    }
			  exit(1);

		      }//if fabs
		    }//if n>0
		  }
		  /* end of Daisy's algorithm*/
		  /**************************************/
		  /* Laina's algorithm for null alleles */
		  /**************************************/

		  logl = 0.0; oldlogl = 0.0;
		  for(n = 0; n < nullstarts; n++)
		    {
		      //printf("\n (305) Laina's MLE startconfig %d nmark = %d, numallele = %d, numind = %d, dataset = %d", n, nmark, numalleles[0], numind, loopindex);fflush(stdout);
		      /****** Set starting values for EM algorithm*************/
		      startallele(nullfreqs, numalleles, nmark);
		      startindf(nullfinblist, numind);

		      beta[0] = 0.2*drand48();/* uniform(0, 0.2)*/

		      nullEMalg(genos, nullfreqs, numalleles, nullfinblist, beta, numind, nmark);
		      if(n==0){logl= nulllikelihood(genos, nullfreqs, beta, nullfinblist, numalleles, numind, nmark);}// if n == 0
		      /**************************************************************************/
		      /*This section basically tests to see if we've ended up in a new place*/
  /**************************************************************************/

		      oldlogl = logl;
		      logl = nulllikelihood(genos, nullfreqs, beta, nullfinblist, numalleles, numind, nmark);


		      if(fabs(oldlogl - logl) > 100*tol){
			printf("\n (341) n = %d. Two starting configs gave different values. seed = %d\n\t oldlogl = %4.6f\n\tnewlogl = %4.6f\nExiting.\n", n, t1, oldlogl, logl);
			exit(1);
		      }/*}if fabs*///}
		    }/* end loop through starting values*/
		  /* end of Laina's algorithm*/
      /*********************************************************************/
  /* beta algorithm: EM for inbreeding coefficients allowing for null  */
  /* alleles, assuming varying missing data rates among the markers    */
  /*********************************************************************/

  for(n = 0; n < M2nullstarts; n++)
    {
      //printf("\n (305) startconfig %d", n);fflush(stdout);
      /****** Set starting values for EM algorithm*************/
      startallele(M2allelefreqs, numalleles, nmark);
      startindf(M2finblist, numind);
      for(j = 0; j<nmark; j++){M2beta[j] = .2*drand48();}

      /*** run EM algorithm to find potential MLE ****/
      betaLaina(genos, M2allelefreqs, numalleles, M2finblist, M2beta, nmark, numind);

      /**/
      if(n==0){logl= betalikelihood(genos, M2allelefreqs, M2beta, M2finblist, numalleles, nmark, numind);}// if n == 0
	//if(n>0){
      oldlogl = logl;
      logl = betalikelihood(genos, M2allelefreqs, M2beta, M2finblist, numalleles, nmark, numind);


      if(fabs(oldlogl - logl) > 0.2){
	printf("\n (377) Two starting configs gave different values. seed = %d\n\t oldlogl = %4.6f\n\tnewlogl = %4.6f\nExiting.\n", t1, oldlogl, logl);
	exit(1);
      }/*}if fabs*///}//if n>0*/
      /*** decide if this new starting configuration gives the best results**/
    }

		  /*****************************************************/
		  /* Kalinowski and Taper's algorithm for null alleles */
		  /*****************************************************/
	           for(n = 0; n < kstarts; n++)
		     {
		       //printf("\n (810) KT startconfig %d", n);fflush(stdout);
		       /****** Set starting values for EM algorithm*************/
		       startallele(kallelefreqs, numalleles, nmark);
		       for(j = 0; j < nmark; j++)
			 {
			   kbeta[j] = .2*drand48();
			 }
		       /*** run EM algorithm to find potential MLE ****/
		       Kalinowski(genos, kallelefreqs, numalleles, kbeta, nmark, numind);

		       /**/
		       if(n==0){logl= klikelihood(genos,kallelefreqs, kbeta, numalleles, nmark, numind);}// if n == 0
		       /**************************************************************************/
		       /*This entire long section basically tests to see if we've ended up in a new place.*/
		       /**************************************************************************/

		       oldlogl = logl;
		       logl = klikelihood(genos, kallelefreqs, kbeta, numalleles, nmark, numind);

		       if(fabs(oldlogl - logl) > 0.2){
			 printf("\n (341) Two starting configs gave different values. seed = %d\n\t oldlogl = %4.6f\n\tnewlogl = %4.6f\nExiting.\n", t1, oldlogl, logl);
			 exit(1);
		       }/*if fabs*/
		       /*** decide if this new starting configuration gives the best results**/
		     }/* end loop through different starting points for*/
		   /* Kalinowski and Taper*/

		   /*********************************************************************/
		   /* Begin Vogl's algorithm */
		   /*********************************************************************/

		   //printf("\n(493) starting Vogl1 (prior alpha = %4.3f, beta = %4.3f)", 1.0, 1.0);
		   vogl(genos, voglfreqs, voglfinblist, numalleles, numind, nmark, 1.0, 1.0);
		   //printf("\n(501) Starting Vogl2 (prior alpha = %4.3f, beta = %4.3f)", voglalpha, voglbeta);
		   vogl(genos, vogl2freqs, vogl2finblist, numalleles, numind, nmark, voglalpha, voglbeta);
		   //printf("\n(502) Finished Vogl2");

		   /*********************************************************************/
		   /* End Vogl's algorithm */

		   /*****************************************************/
		   /* Now we find a few marker facts for printout. */
		   /*****************************************************/

		   for(j = 0; j < nmark; j++)
		     {
		       countmiss[j] = 0;
		       counthom[j] = 0;
		       counthet[j] = 0;
		       for(i = 0; i < numind; i++)
			 {
			   if(genos[i][j][0] == 99){countmiss[j]++;}
			   else if(genos[i][j][0] == genos[i][j][1]){counthom[j]++;}
			   else{counthet[j]++;}
			 }
		     }
		   /*****************************************/
		   /*	printing to file               	 */
		   /*****************************************/

		   for(i = 0; i < numfs; i++)
		    {
		      fprintf(foutout, "%d %d %d %d %d % 4.4f  % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f % 4.4f\n", simnumber, loopindex, nmark, numalleles[0], numind, finbtrue[i], finblist[i], nfinblist[i], simpfinblist[i], nsimpfinblist[i], modsimpfinblist[i], nmodsimpfinblist[i], ritfinblist[i], nritfinblist[i], modritfinblist[i], nmodritfinblist[i], jointfinblist[i], nullfinblist[i], M2finblist[i], voglfinblist[i], vogl2finblist[i]);
		    }
		  for(i = 0; i < 3; i++)
		    {
			fprintf(foutfreqs, "%3d %4d %3d %2d %3d ", simnumber, loopindex, nmark, numalleles[0]-1, numind);

		      fprintf(foutfreqs, " %4.4f %4.4f %4.4f %4.4f %4.4f %4.4f", truefreqs[i][numalleles[i]-1], nullfreqs[i][numalleles[i]-1], M2allelefreqs[i][numalleles[i]-1], kallelefreqs[i][numalleles[i]-1], voglfreqs[i][numalleles[i]-1], vogl2freqs[i][numalleles[i]-1]);

		       fprintf(foutmarker, "%4d %3d %3d %2d %3d %4d %4d %4d %4.3f %4.3f %4.3f %4.3f %4.3f", simnumber, loopindex, nmark, numalleles[i], numind, countmiss[i], counthom[i], counthet[i], missrate, beta[0], M2beta[i], kbeta[i], truefreqs[i][numalleles[i]-1]);
		      fprintf(foutmarker, " %4.4f %4.4f %4.4f %4.3f %4.3f", nullfreqs[i][numalleles[i]-1], M2allelefreqs[i][numalleles[i]-1], kallelefreqs[i][numalleles[i]-1], voglfreqs[i][numalleles[i]-1], vogl2freqs[i][numalleles[i]-1]);
		      fprintf(foutfreqs, " 967");
		      for(k = 0; k < numalleles[i]-1; k++)
			{
			  fprintf(foutfreqs, " %4.4f", allelefreqs[i][k]);
			}
		      for(k = numalleles[i]; k < maxnumalleles; k++)
			{ fprintf(foutfreqs, " 0.0000");}
		      /* null MLE estimated allelefreqs*/
		      fprintf(foutfreqs, " 999");
		      for(k = 0; k < numalleles[i]-1; k++)
			{
			  fprintf(foutfreqs, " %4.4f", nullfreqs[i][k]);
			}
		      for(k = numalleles[i]; k < maxnumalleles; k++)
			{ fprintf(foutfreqs, " 0.0000");}
		      /* model 2 null allele/ differing beta's output*/
		      fprintf(foutfreqs, " 999");
		      for(k = 0; k < numalleles[i]-1; k++)
			{
			  fprintf(foutfreqs, " %4.4f", M2allelefreqs[i][k]);
			}
		      for(k = numalleles[i]; k < maxnumalleles; k++)
			{ fprintf(foutfreqs, " 0.0000");}
		      /* kalinowski and taper estimated allelefreqs*/
		      fprintf(foutfreqs, " 999");
		      for(k = 0; k < numalleles[i]-1; k++)
			{
			  fprintf(foutfreqs, " %4.4f", kallelefreqs[i][k]);
			}
		      for(k = numalleles[i]; k < maxnumalleles; k++)
			{ fprintf(foutfreqs, " 0.0000");}
		      /* Vogl1 estimated allelefreqs*/
		      fprintf(foutfreqs, " 999");
		      for(k = 0; k < numalleles[i]-1; k++)
			{
			  fprintf(foutfreqs, " %4.4f", voglfreqs[i][k]);
			}
		      for(k = numalleles[i]; k < maxnumalleles; k++)
			{ fprintf(foutfreqs, " 0.0000");}
		      /* Vogl2 estimated allelefreqs*/
		      fprintf(foutfreqs, " 999");
		      for(k = 0; k < numalleles[i]-1; k++)
			{
			  fprintf(foutfreqs, " %4.4f", vogl2freqs[i][k]);
			}
		      for(k = numalleles[i]; k < maxnumalleles; k++)
			{ fprintf(foutfreqs, " 0.0000");}
		      /* */
		      fprintf(foutfreqs, "\n");
		      fprintf(foutmarker, "\n");
		    }
		  fflush(foutfreqs);
		  fflush(foutmarker);
		  free(recomb);
		  if(linked > 0){free(nmarkchrom);}

		}/* end loopindex */
	  /*****************************************************/
	  /* Freeing things that were "malloc'd" in this loop. */
	  /*****************************************************/
	  free(numalleles);
	  free(tempnumalleles);
	  for(j = 0; j < numind; j++)
	    {
	      for(k = 0; k < nmark; k++)
		{
		  free(genos[j][k]);
		}
	      free(genos[j]);
	    }
	  free(genos);

	  printf("(886)");fflush(stdout);

	  for(k = 0; k < nmark; k++)
	    {
	      free(naivefreqs[k]);
	      free(allelefreqs[k]);
	      free(kallelefreqs[k]);
	      free(M2allelefreqs[k]);
	      free(truefreqs[k]);
	      free(nullfreqs[k]);
	    }
	  free(allelefreqs);
	  free(kallelefreqs);
	  free(naivefreqs);
	  free(M2allelefreqs);
	  free(truefreqs);
	  free(nullfreqs);
	  free(asum);
	  if(numalleleindex == 0){t3 = t1;}
	  else
	    {
	      t3 = t2;
	    }
	  time(&t2);
	  //printf("\n(736) numalleleindex = %d, t2 = %d", numalleleindex, t2);
	  //printf("\n t2 = %d; t3 = %d; t1 = %d",t2, t3, t1);fflush(stdout);
	  printf("\n Finished with numallele loop. numind = %4d, nmark = %4d, numallele = %3d, Time for loop: %4.3f min", numind, nmark, numallelearray[numalleleindex], ((double)t2-(double)t3)/60.0);fflush(stdout);
	}/* end for(numalleleindex = 0; ...*/
	 printf("(973)");fflush(stdout);
	 free(kbeta);
	 free(M2beta);
	 free(counthom);
	 free(counthet);
	 free(countmiss);
    }/* end for(nmarkindex = 0; nmarkindex < numnmark; nmarkindex++)...*/

      free(simpfinblist);
      free(nsimpfinblist);
      free(modsimpfinblist);
      free(nmodsimpfinblist);
      free(ritfinblist);
      free(nritfinblist);
      free(modritfinblist);
      free(nmodritfinblist);
      free(finblist);
      free(nfinblist);
      free(nullfinblist);
      free(M2finblist);
      free(finbtrue);
      free(jointfinblist);
    }/* end loop through different numind values */
  /***** Supplementary output *****************************/

  (void) time(&t2);
  fprintf(foutlog, "\n Running time %4.3f minutes\n", (double) (t2-t1)/60.0);

  fprintf(foutout, "\n");

      printf("\n\n Running time %4.3f minutes\n", (double) (t2-t1)/60.0);

  fflush(foutout);
  fflush(foutlog);
  fclose(foutout);
  fclose(foutlog);
  printf("(620)");fflush(stdout);

  free(beta);
} /* end main */


/***********************************************************************/

int generate(int toprint, int nmark, int* numalleles, int*** genos, double** truefreqs, double *recomb, double* finbtrue, int numind)
{
  int i,j,k, markindex;
  double sum;

  int happoint; /* points to which haplotype (maternal or paternal) we are */
                /* copying.*/
  int count;
  double x;

  /* initializing genotypes */
  for(i = 0; i < numind; i++)
    {
      for(markindex = 0; markindex < nmark; markindex++)
	{
	  genos[i][markindex][0] = 0;
	  genos[i][markindex][1] = 0;
	}
    }
  /******************************************************************/
  /* And then we set their genotypes.                               */
  /*******************************************************************/

 for(i = 0; i < numind; i++)
    {
      for(markindex = 0; markindex < nmark; markindex++)
	{
	  /* set first allele*/
	  sum = truefreqs[markindex][0];
	  count = 0;
	  x = drand48();
	  while((sum <= 1.001) && (count <= numalleles[markindex]))
	    {
	      if(x < sum)
		{
		  genos[i][markindex][0] = count;
		  sum = 1.2; /* signals loop to end */
		}
	      else
		{
		  count++;
		  sum += truefreqs[markindex][count];
		}
	    }/* end of while loop */

	  /* set second allele*/
	  x = drand48();
	  if(x < finbtrue[i])/* if the two alleles are IBD*/
	    {
	      genos[i][markindex][1] = genos[i][markindex][0];
	    }
	  else /* if not IBD*/
	    {
	  sum = truefreqs[markindex][0];
	  count = 0;
	  x = drand48();
	  while((sum <= 1.001) && (count <= numalleles[markindex]))
	    {
	      if(x < sum)
		{
		  genos[i][markindex][1] = count;
		  sum = 1.2; /* signals loop to end */
		}
	      else
		{
		  count++;
		  sum += truefreqs[markindex][count];
		}
	    }/* end of while loop */

	    }
	}/* end for each marker */
    }/* end loop through all individuals*/

  /******************************************************************/
  /* Now we delete all missing genotypes and deal with null alleles*/
  /******************************************************************/
  for(i = 0; i < numind; i++)
    {
      for(markindex = 0; markindex < nmark; markindex++)
	{
	  x = drand48();
	  if(x < missrate)
	    {
	      genos[i][markindex][0] = 99;
	      genos[i][markindex][1] = 99;
	    }
	  if((alleletype == 8)||(alleletype == 9))
	    {
	      if(genos[i][markindex][0] == numalleles[markindex]-1)
		{
		  if(genos[i][markindex][1] == numalleles[markindex]-1)
		    {/* if both alleles are null, then it shows up as missing*/
		      genos[i][markindex][0] = 99;
		      genos[i][markindex][1] = 99;
		    }
		  else /* if only the 0 allele is null */
		    {/* we read as if it were homozygous for the non-null allele*/
		      genos[i][markindex][0] = genos[i][markindex][1];
		    }
		}
	      else if(genos[i][markindex][1] == numalleles[markindex]-1)
		{/* if only the second allele is null, we read both alleles as being the first allele */
		  genos[i][markindex][1] = genos[i][markindex][0];
		}
	    }
	}
    }
  /* Now we have the genotypes. Let's look at them */
  /*printf("\n (1823)At the end of generate:");
    for(i = 0; i < numind; i++)
    {
    printf("\n %2d:", i+1);
    for(markindex = 0; markindex < nmark; markindex++)
    {
    printf(" %d%d", genos[i][markindex][0], genos[i][markindex][1]);
    }
    }*/

  return 0;

}
/***********************************************************************/
double rexp(double lambda)
     /*****************************************/
     /* Generate exponential random variable */
     /***************************************/
     /* double lambda;    Rate */
{
  return -log(drand48())/lambda;
  }

/*---------------------------------------*/
double rgamma(double alpha)
     /***********************************************************************/
     /* Generate gamma random variable--the gamma with shape parameter alpha*/
     /* and scale parameter 1.0: with pdf:                                  */
     /* f(alpha) = x^(alpha - 1)e^(-x)/Gamma(alpha).                        */
     /***********************************************************************/
     /* double alpha;   Shape parameter */
{
  double r1,r2,aa,x,w,c1,c2,c3,c4,c5;
  /* new line:*/
  double EE = 2.718281828459045;
  if (alpha<=0.) return 0.;
  if (alpha == 1.) return rexp(1.);
  if (alpha<1)
    {
      aa=(alpha+EE)/EE;
      //printf("\n Ouch! there's something wrong with my gamma generator. I can't do gamma's with parameter values between zero and one...long story. It has to do with not being able to find a library that my code wants to call and then (related) not knowing what a certain symbol (EE) means. I think it's probably the exponential E, but I'm not sure. Exiting.");

      do
	{
	  r1=drand48();
	  r2=drand48();
	  if (r1>1./aa)
	    {
	      x = -log(aa*(1.-r1)/alpha);
	      if (r2<pow(x,(alpha-1.))) return x;
	    }
	  else
	    {
	      x = pow((aa*r1),(1./alpha));
	      if (r2<exp(-x)) return x;
	    }
	}
      while (r2<2);
    }
  else
    {
      c1 = alpha-1;
      c2 = (alpha-1./(6.*alpha))/c1;
      c3 = 2./c1;
      c4 = c3+2.;
      c5 = 1./sqrt(alpha);
      do
	{
         do
	   {
	     r1=drand48();
	      r2=drand48();
	     if (alpha>2.5) r1=r2+c5*(1.-1.86*r1);
	   }
	 while (r1<=0 || r1 >= 1);
         w = c2*r2/r1;
         if (c3*r1+w+1/w <= c4) return c1*w;
         if (c3*log(r1)-log(w)+w<1) return c1*w;
	}
      while (r2<2);
    }

/* This is to remove the warning statement */
  return 0.;
}

double rbeta(double alpha, double beta)
    /**********************************/
     /* Generate beta random variable */
     /********************************/
     /* double alpha, beta;   Shape parameters */
{
  double r1;
  if (alpha <=0. || beta <= 0.) return 0.;
  r1 = rgamma(alpha);
  return r1/(r1+rgamma(beta));
  }
/*---------------------------------------*/

/*---------------------------------------*/

void rdirichlet(double *alpha,  double *p, int mynumalleles)
     /***************************************/
     /* Generate dirichlet random variable */
     /*************************************/
     /* double *alpha;   Shape parameters */
     /* int k;           Dimension */
{
  int i, k;
  double r[mynumalleles];
  double total, rgamma(double);

  k = mynumalleles;

  total = 0.0;
  for (i=0; i<k; i++)
    {
      r[i] = rgamma(alpha[i]);
      total += r[i];
    }
  for (i=0; i<k; i++)
    {
      p[i] = r[i]/total;
    }
  //printf("\n In the function, allele freqs are %f %f %f %f\n ", p[0], p[1], p[2], p[3]);

  return;
}
/************************************************************************/
/************* Nathan's algorithm: assumes allele freqs known and *******/
/******** calculates the inbreeding coefficient for each person individually*/
//Apply EM for all ind
void EMalg(int*** genos, double** allelefreqs, int nmark, int numind, double* finblist)
{
  double valloglik(int** myind, double** allelefreqs, int nmark, double finb);

  int i,j,k,m, n;
  double liktemp;
  double bestfinb=0.0;
  double maxliktemp=0.0;
  double finb=.5;
  double finbold=0;
  double finbsum;
  double startvals[numstarts];
  double *probIBD; //prob IBD at marker M (E[z|data])

  for(j = 0; j < numstarts; j++)
    {
      startvals[j] = drand48();
    }

	probIBD=(double*)malloc(nmark*sizeof(double));

	for (j=0; j<numind; j++)
	  {
	    for (m=0; m<numstarts; m++)
	      {
		//finb=.5;
		finb=startvals[m];
		finbold=0;
		n = 0;
		while((fabs(finbold-finb)>nathantol) && (n < MAXNUMTIMES))
		  {
		    n++;
		    finbsum=0.0;
		    finbold=finb;
		    k=0;
		    for (i=0; i<nmark; i++){
		      if ((genos[j][i][0]==genos[j][i][1]) && (genos[j][i][0]!=99))
			probIBD[i]=((finb)*(allelefreqs[i][genos[j][i][0]]))/((finb)*(allelefreqs[i][genos[j][i][0]])+(1-finb)*pow(allelefreqs[i][genos[j][i][0]],2));
		      else if ((genos[j][i][0]==99) || (genos[j][i][1]==99))
			{
			  probIBD[i]=0;
			  k++;
			}
		      else
			{
			  probIBD[i]=0;
			}
		    }

		    for (i=0; i<nmark; i++)
		      {
			finbsum+=probIBD[i];
		      }

		    finb=finbsum/(nmark-k);			//calculate new value of finb
		  }	//end of while loop

		if(n > 2500){printf("\n\t (1065)ind = %d, number EM iterations = %d", j+1, n);}
		if(n >= MAXNUMTIMES)
		  {
		    printf("\n(937) Failure for EM algorithm to converge. ind = %d, n = %d. exiting.\n", j+1, n);
		    exit(0);
		  }
		liktemp=valloglik(genos[j], allelefreqs, nmark, finb);

		//compare likelihood values and keep best one, recording associated finb
		if (m==0)
		  {
		    maxliktemp=liktemp;
		    bestfinb=finb;
		  }
		else
		  {
		    if ((fabs(liktemp-maxliktemp)>liktol) && (fabs(bestfinb - finb) > nathantol))
		      {
			printf("\n(954) Nathan's program. Different likelihoods for different finb start values. (Individual %d, m = %d)\nQuitting. . .\n", (j+1), m);
			exit(0);
		      }
		    if(liktemp>maxliktemp)
		      {
			maxliktemp=liktemp;
			bestfinb=finb;
		      }
		  }

		//printf("numind: %d startval: %lf finb: %lf likelihood: %lf \n\t maxlik: %lf bestfinb: %lf\n", j, startvals[m], finb, liktemp, maxliktemp, bestfinb); //test output
	      }
	    finblist[j]=bestfinb;
	  } //numind for loop

	free(probIBD);

	return;
}

/*****************************************/
/*	  Nathan's second program, valloglik			 */
/*	  calculates loglikelihood values    */
/*	  for different convergence values for the basic algorithm.  */
/*****************************************/

double valloglik(int** myind, double** allelefreqs, int nmark, double finb)
{
	int i,j,k;
	double loglik=0;

	for (i=0; i<nmark; i++)
	{
		if (myind[i][0]!=99)
		{
			if (myind[i][0]==myind[i][1])
			{
				loglik+=log(finb*allelefreqs[i][myind[i][0]]+(1-finb)*pow(allelefreqs[i][myind[i][0]],2));
			}
			else
			{
				loglik+=log(2*(1-finb)*allelefreqs[i][myind[i][0]]*allelefreqs[i][myind[i][1]]);
			}
		}
	}
	return loglik;
}

/*****************************************/
/*	  Third program, simplecalc			 */
/*	  calculates "simple" value		     */
/*	  for finb given hetero vs homo      */
/*****************************************/

void simplecalc(int*** genos, double** allelefreqs, int nmark, int numind, double* simpfinblist, int* numalleles)
{
	int i,j,k,m;
	double numhetero;	//total number of hetero for each ind
	double exphetero;	//expected hetero at marker i
	double sumexphetero;

	for (j=0; j<numind; j++)
	{
	  k=0;
	  numhetero=0;
	  sumexphetero=0;
	  for (i=0; i<nmark; i++)
	    {
	      exphetero=0;
	      if ((genos[j][i][0]!=genos[j][i][1]))
		{
		  numhetero++;
		}
	      else if ((genos[j][i][0]==99) || (genos[j][i][1]==99))
		{
		  k++;
		}
	      if ((genos[j][i][0]!=99))
		{
		  for (m=0; m<numalleles[i]; m++)
		    {
		      exphetero+=pow(allelefreqs[i][m],2);
		    }
		  sumexphetero+=(1-exphetero);
		}
	    }
	  simpfinblist[j]=(1-((numhetero/(nmark-k))/(sumexphetero/(nmark-k))));
	}
	return;
}


/*****************************************/
/*	  Fifth program, RitlandMME			 */
/*	  calculates value for finb		     */
/*	  based on Ritland article.     	 */
/*****************************************/
void RitlandMME(int*** genos, double** allelefreqs, int nmark, int* numalleles, int numind, double* ritfinblist)
{
	int i,j,k,l;
	double S;
	double numeratorsum;
	double denominatorsum;

	for (k=0; k<numind; k++)
	{
	  numeratorsum=0.0;
	  denominatorsum=0.0;
	  for (l=0; l<nmark; l++)
	    {
	      if(genos[k][l][0]!=99)
		{
		  for (i=0; i<numalleles[l]; i++)
		    {
		      if(allelefreqs[l][i]>0)
			{
			  S=0.0;
			  if ((genos[k][l][0]==genos[k][l][1])&&(genos[k][l][0]==i))
			    {
			      S=1.0;
			    }
			  else
			    {
			      S=0.0;
			    }
			  numeratorsum+=(S-pow(allelefreqs[l][i],2))/allelefreqs[l][i];
			}
		    }
		  denominatorsum+=(numalleles[l]-1.0);
		}
	    }
	  if(denominatorsum < 0.0000001)
	    {
	      printf("\n(1235) Ritland denominator 0. k = %d, numeratorsum = %4.4f. exiting\n", k, numeratorsum);
	      exit(0);
	    }
	  ritfinblist[k]=(numeratorsum/denominatorsum);
	  //printf("\n k = %d, ritfinblist[k] = %4.4f numer = %4.4f, denom = %4.4f", k, ritfinblist[k], numeratorsum, denominatorsum);
	}
	return;
}
/***************************************************************************/
/*Now we call a different program to calculate the log likelihood of our   */
/*data set given our indf and allele frequencies. This is the loglikelihood*/
/* for Daisy's program.                                                    */
/***************************************************************************/
  double loglike(int*** genos, double** allelefreq, double* indf, int numind, int nmark){
	  double p,logl;
	  int i,j,m,n;
	  logl = 0.0;
	  for(i=0;i<numind;i++){
		  for(j=0;j<nmark;j++){
			  if(genos[i][j][0] != 99){
			 m=genos[i][j][0];
			 n=genos[i][j][1];
			 p=(1-indf[i])*2*allelefreq[j][m]*allelefreq[j][n];
			 if (m==n){
				p=allelefreq[j][m]*indf[i]+(1-indf[i])*pow(allelefreq[j][m],2);
			  }//if
			  logl=logl+log(p);
			  }//if
		  }//for j
		 }//for i
//printf(" loglikelihood = %4.8f", logl);
return logl;
  }// loglike function


/************************************************************************/
/**** Daisy's finbreeding algorithm that runs an EM algorithm to jointly*/
/**** estimate inbreeding coefficients and allele frequencies. **********/
/************************************************************************/
void finbreeding(int*** mydata, double** allelefreq, int* numalleles, double* indf, int numind, int nmark){
  int tempnmark = nmark;
  double fold = 0;
  int i,j,k,l;
  double** finb;
  double numerator,denominator,p,q,m;
  double logl;
  double oldlogl;
  int ntimes;

  finb=(double**)malloc(numind*sizeof(double*));
  for(i=0;i<numind;i++){finb[i]=(double*)malloc(nmark*sizeof(double));}

  /*now start the real work...find array of X's, (p(IBD at each marker))*/
  /*we set it so that X is 0.0 at every marker*/
  for(i=0;i<numind;i++){
    for(j=0;j<nmark;j++){finb[i][j]=0.0;}//for j
  }//for i

  //begin while loop, which keeps going until logl gets below a tolerance, or until we've gone through a ton of iterations
  logl=2.0;
  oldlogl=1.0;
  ntimes=0;
  while((fabs(logl-oldlogl)>tol)&&(ntimes<MAXNUMTIMES)){
    //this part of the loops uses allelefreqs, mydata, and individual inbreeding coefficients to calculate X matrix
    //printf("\n Daisy iteration %d", ntimes);fflush(stdout);
    for(i=0;i<numind;i++){
      for(j=0;j<nmark;j++){
	if(mydata[i][j][0] != 99){
	  k=mydata[i][j][0];
	  l=mydata[i][j][1];
	  if(k == l){
	    m=indf[i];
	    finb[i][j]=((m*allelefreq[j][k])/(m*allelefreq[j][k]+(1-m)*pow(allelefreq[j][k],2)));

	    if((finb[i][j] >= 1.0000001) || (finb[i][j] < 0))
	      {
		printf("\n(292) Individual %d marker %d geno %d %d; finb = %4.6f\n", i, j, k, l, finb[i][j]);
		exit(0);
	      }
	  }//if
	}//if
      }//for j
    }//for i

    //this part uses mydata and X matrix to calculate allelefreqs
    for(j=0;j<nmark;j++){
      denominator=0;
      for(i=0;i<numind;i++){
	if(mydata[i][j][0] != 99){
	  q=2;
	  if(mydata[i][j][0]==mydata[i][j][1]){
	    q=2-finb[i][j];
	  }//if2
	  denominator+=q;
	}//if1
      }//for i
      for(k=0;k<numalleles[j];k++){
	numerator=0;
	for(i=0;i<numind;i++){
	  p=0;
	  if(mydata[i][j][0] != 99){
	    if((mydata[i][j][0]==(k))||(mydata[i][j][1]==(k))){
	      p=1;
	      if(mydata[i][j][0]==mydata[i][j][1]){
		p=2-finb[i][j];
	      }//if1
	    }//if2
	  }//if3
	  numerator=numerator+p;
	}// for i
	allelefreq[j][k]=numerator/denominator;
	if((allelefreq[j][k] > 1)|| (allelefreq[j][k] < 0))
	  {
	    printf("\n(329) allelefreq[%d][%d] = %4.4e\n exiting.", j, k, allelefreq[j][k]);exit(0);
	  }
      }//for k
    }// for j

    //this part calculates indf - the individual inbreeding coefficients - by summing up entries of X for each individual.

    for(i=0;i<numind;i++){
      m=0;
      k=0;
      for(j=0;j<nmark;j++){
	if(mydata[i][j][0] != 99){
	  k=k+1;
	  m=m+finb[i][j];
	}// if
      }// for j
      indf[i]=m/(double)k;

    }//for i
    ntimes++;
    oldlogl=logl;
    logl= loglike(mydata, allelefreq, indf, numind, nmark);
    m=fabs(oldlogl-logl);

    if(ntimes>=MAXNUMTIMES){
      printf("(381) failure to converge.  Difference between subsequent log likelihoods is %f\n", m); exit(1);
    }
  }//for while loop

  for(i=0;i<numind;i++){free(finb[i]);}//for i
  free(finb);
  return;

}

/****************************************************************************/
/**************************************************************************/
/* the next two functions set the starting values for the EM algorithm for */
/* joint estimation of inbreeding coefficients and allele frequencies.*/
/**************************************************************************/
  void startindf(double* randindf, int numind)
  {
	  int i;
	  double r1;

for(i=0;i<numind;i++)
{
	r1 = (double)((rand()/(RAND_MAX + 1.0)));
		randindf[i] = (double)((rand()/(RAND_MAX + 1.0)));
		//printf("\n%.4f", randindf[i]);
	}//for i

return;
  }
/**************************************************************************************************************/
void  startallele(double** randallelefreq, int* numalleles, int nmark)
{
  int j, i, k;
  double* alpha;
  int good = 1;

  alpha = (double*)malloc(numalleles[0]*sizeof(double));
  for(j=0;j<nmark;j++)
    {
      for(i=0;i<numalleles[j];i++)
	{
	  randallelefreq[j][i] = 0.0;
	}
    }

  for(j=0;j<nmark;j++)
    {
      for(i=0;i<numalleles[j];i++)
	{
	  alpha[i] = 1.0;
	}
      k = numalleles[j];

      good = 1;
      while(good > 0)
	{
	  rdirichlet(alpha, randallelefreq[j], k);
	  good = 0;
	  for(i=0;i<numalleles[j];i++){
	    if(randallelefreq[j][i] < alleletol)
	    {
	      good = 1;
	    }
	}
      }
  }//for j
  free(alpha);
  return;
}
/***************************************************************************/
/********* likelihood function for kalinowski and taper**************/
/***************************************************************************/
double klikelihood(int*** myind1, double** allelefreqs, double* beta, int* numalleles, int nmark, int numind)
{
  int i, j, l, k;
  double pnull, pmissing, LogL;

  LogL=0.0;
	for(i=0;i<numind;i++)
	{
		for(j=0;j<nmark;j++)
		{
		  k=myind1[i][j][0];
		  l=myind1[i][j][1];
		  pnull=allelefreqs[j][numalleles[j]-1];
		  if(k==99)/* if the genotype is missing */
		    {
		      LogL = LogL + log(beta[j]+(1-beta[j])*pow(pnull,2));
		    }
		  else if(k==l)/* if it is homozygous (nonmissing)*/
		    {
		      LogL = LogL + log((1-beta[j])*(pow(allelefreqs[j][k],2)+2*allelefreqs[j][k]*pnull));
		    }
		  else/* if it is heterozygous (nonmissing)*/
		    {
		      LogL = LogL + log((1-beta[j])*2*allelefreqs[j][k]*allelefreqs[j][l]);
		    }
		}
	}
return LogL;

}
/**********************************************************************/ /********* likelihood function--gives the loglikelihood value for ***********/
/* the nullalleleEMalgorithm*****************************/
double nulllikelihood(int*** myind1, double** allelefreqs, double* beta, double* finb, int* numalleles, int numind, int nmark)
{
  int i, j, l, a, k, m;
  int markindex, loopindex, index, index2, index3, index4;
  double sum;
  double pnull, pmissing, LogL;

  LogL=0.0;
  for(i=0;i<numind;i++)
    {
      for(j=0;j<nmark;j++)
	{
	  k=myind1[i][j][0];
	  l=myind1[i][j][1];

	  pnull=allelefreqs[j][numalleles[j]-1];
	  if(k==99)
	    {
	      LogL = LogL + log(beta[0]+(1-beta[0])*(finb[i]*pnull+(1-finb[i])*pow(pnull,2)));
	    }
	  else if(k==l)/**** I did this so it would compile, I couldn't get the elseif to work****/
	    {
	      LogL = LogL + log((1-beta[0])*(finb[i]*allelefreqs[j][k]+(1-finb[i])*pow(allelefreqs[j][k],2)+(1-finb[i])*2*allelefreqs[j][k]*pnull));
	    }
	  else
	    {
	      LogL = LogL + log((1-beta[0])*(1.0 - finb[i])*2*allelefreqs[j][k]*allelefreqs[j][l]);
	    }
	}
    }
return LogL;

}
/************************************************************************/
/********* likelihood function for nullEM with a separate beta for each */
/* markers. (That is, the rate at which genotypes are missing at random */
/* varies among the different markers.*/
/***************************************************************************/
double betalikelihood(int*** myind1, double** allelefreqs, double* beta, double* finb, int* numalleles, int nmark, int numind)
{
  int i, j, l, a, k;
  double sum;
  double pnull, LogL;

  LogL=0.0;
  for(i=0;i<numind;i++)
    {
      for(j=0;j<nmark;j++)
	{
	  k=myind1[i][j][0];
	  l=myind1[i][j][1];
	  pnull=allelefreqs[j][numalleles[j]-1];
    if(j >= numNonXLinked && sex[i]==1)
      {
      if(k==99)
        {
        LogL = LogL + log(beta[j]+(1-beta[j])*pnull);
        }
      else
        {
        LogL = LogL + log((1-beta[j])*(allelefreqs[j][k]));
        }
      }
	  else if(k==99)
	    {
	      LogL = LogL + log(beta[j]+(1-beta[j])*(finb[i]*pnull+(1-finb[i])*pow(pnull,2)));
	    }
	  else if(k==l)/**** I did this so it would compile, I couldn't get the elseif to work****/
		    {
		      LogL = LogL + log((1-beta[j])*(finb[i]*allelefreqs[j][k]+(1-finb[i])*pow(allelefreqs[j][k],2)+(1-finb[i])*2*allelefreqs[j][k]*pnull));
		    }
		  else
		    {
		      LogL = LogL + log((1-beta[j])*(1.0 - finb[i])*2*allelefreqs[j][k]*allelefreqs[j][l]);
		    }
		}
	}
return LogL;

}

/****************************************************************/

void nullEMalg(int*** myind1, double** allelefreqs, int* numalleles, double* finb, double* beta, int numind, int nmark)
{
  int i, j, l, a, k, m, markindex, loopindex, index, index2, index3, index4;

  double sum, sum0, sum1, sumtop, sumbottom, likval, oldlikval;
  double** ExpX; /*Exp[numind][nmark]*/
  /* probgenos is the prob that an individual has each genotype at each */
  /* marker. We need this to be conditional on whether the two alleles are*/
  /* IBD or not, hence the 0 (not IBD) and the 1 (IBD). So */
  /* Probgenos0 = Pr(genotype, xij|obs data, current param values)*/
  double**** Probgenos0; /*Probgenos[numind][nmark][numalleles][numalleles]*/
  double**** Probgenos1; /*Probgenos[numind][nmark][numalleles][numalleles]*/
  double** ExpB;/*ExpB[numind][nmark]*/

  /* pmissing = Pr(observed = missing|data). pmissing0 = Pr(obs = missing and notIBD|data), pmissing1 = Pr(obs = missing and IBD| data) */
  double pnull, pmissing;
  //double pmissing0, pmissing1;
  double pmissing1notf, pmissing0notf; /* pmissing1 is a product of f times something. pmissing1notf is that something.*/

  /* allocating memory *****/
  ExpX = (double**)malloc(numind*sizeof(double*));
  ExpB = (double**)malloc(numind*sizeof(double*));
  for(i = 0; i < numind; i++)
    {
      ExpX[i] = (double*)malloc(nmark*sizeof(double));
      ExpB[i] = (double*)malloc(nmark*sizeof(double));
    }
  Probgenos0 = (double****)malloc(numind*sizeof(double***));
  Probgenos1 = (double****)malloc(numind*sizeof(double***));
  if(Probgenos0 == NULL)
    {
      printf("\n(488) Problem allocating memory for Probgenos0. Exiting.\n");exit(0);
    };
  if(Probgenos1 == NULL)
    {
      printf("\n(497) Problem allocating memory for Probgenos1. Exiting.\n");exit(0);
    }

  for(i = 0; i < numind; i++)
    {
      Probgenos0[i] = (double***)malloc(nmark*sizeof(double**));
      Probgenos1[i] = (double***)malloc(nmark*sizeof(double**));
      if(Probgenos0[i] == NULL)
	{
	  printf("\n(506) Problem allocating memory for Probgenos0[i = %d]. Exiting.\n", i);exit(0);
	}
      if(Probgenos0[i] == NULL)
	{
	  printf("\n(510) Problem allocating memory for Probgenos1[i = %d]. Exiting.\n", i);exit(0);
	}
      for(markindex = 0; markindex < nmark; markindex++)
        {
          Probgenos0[i][markindex] = (double**)malloc(numalleles[markindex]*sizeof(double*));
          Probgenos1[i][markindex] = (double**)malloc(numalleles[markindex]*sizeof(double*));
	  if(Probgenos0[i][markindex] == NULL)
	    {
	      printf("\n(1801) Problem allocating memory for Probgenos0[i = %d][m = %d]. Exiting.\n", i, markindex);exit(0);
	    }
	  if(Probgenos1[i][markindex] == NULL)
	    {
	      printf("\n(1805) Problem allocating memory for Probgenos1[i = %d][m = %d]. Exiting.\n", i, markindex);exit(0);
	    }
          for(j = 0; j < numalleles[markindex]; j++)
            {
              Probgenos0[i][markindex][j] = (double*)malloc(numalleles[markindex]*sizeof(double));
              Probgenos1[i][markindex][j] = (double*)malloc(numalleles[markindex]*sizeof(double));
	      if(Probgenos0[i][markindex][j] == NULL)
		{
		  printf("\n(530) Problem allocating memory for Probgenos0[i = %d][m = %d][j = %d]. Exiting.\n", i, markindex, j);exit(0);
		}
	      if(Probgenos1[i][markindex][j] == NULL)
		{
		  printf("\n(534) Problem allocating memory for Probgenos1[i = %d][m = %d][j = %d]. Exiting.\n", i, markindex, j);exit(0);
		}
              for(k = 0; k < numalleles[markindex]; k++)
                {
                  Probgenos0[i][markindex][j][k] = 0.0;
                  Probgenos1[i][markindex][j][k] = 0.0;
                }
            }
        }
    }

  /*printf("\n(498) WARNING: I'm resetting the starting values as a debugging\nexercise\n");fflush(stdout);
    for(j = 0; j < nmark; j++)
    {
    allelefreqs[j][0] = 0.5;
    allelefreqs[j][1] = 0.3;
    allelefreqs[j][2] = 0.2;
    }
    for(i = 0; i < numind; i++)
    {
    finb[i] = 0.2;
    }
    beta[0] = 0.1;*/

  likval=1.0;
  oldlikval=2.0;

  /* Next there will be a giant while loop that runs the EM algorithm. */
  /* When you are first writing this, make it a for loop*/
  //for(loopindex = 0; loopindex < 10; loopindex++)
  loopindex = 0;
  while(fabs(likval-oldlikval)>tol)
  {
    loopindex = loopindex + 1;
    oldlikval = likval;
    //printf("\n \t(359) Start of %dth iteration through EM algorithm loglik = %4.8f", loopindex, oldlikval);fflush(stdout);

    for(i=0;i<numind;i++)
      {
	for(j=0;j<nmark;j++)
	  {
	    k=myind1[i][j][0];
	    l=myind1[i][j][1];
	    pnull=allelefreqs[j][numalleles[j]-1];
	    pmissing=beta[0] + (1-beta[0])*(finb[i]*pnull+(1-finb[i])*pow(pnull,2));
	    /* pmissing0 and pmissing1 are Pr(geno = M and notIBD) and Pr(geno = M and IBD), respectively*/
	    //pmissing0 = beta[0]*(1-finb[i]) + (1-beta[0])*(1-finb[i])*pow(pnull,2);
	    //pmissing1 = beta[0]*(finb[i]) + (1-beta[0])*finb[i]*pnull;
	    /* the probability that an individual's genotype is missing, */
	    /* given that their alleles are IBD (ie pmissing1) has a factor*/
	    /* of finb[i] in it. When we use pmissing1 to find values in */
	    /*Probgenos1, the probgenos value has the form:                */
	    /* finb[i]*something/pmissing1. We are having problems with this*/
	    /* when finb[i] is converging to zero. The solution is to cancel*/
	    /* the finb[i] term out of top and bottom of the ratio. Hence,*/
	    /* instead of pmissing1, we want pmissing1 without the finb[i].*/
	    /* that is, we want pmissing1notf.*/
	    pmissing1notf = beta[0] + (1-beta[0])*pnull;
	    pmissing0notf =  beta[0] + (1-beta[0])*pow(pnull,2);
	    /*if((i == 1) && (j == 2) && (nmark == 5) && (numalleles[j] == 6) && (numind == 20))
	      {
	      printf("\n(1666) beta = %4.4e, finb = %4.4e, pnull = %4.ef, pmissing1notf = %4.3e\n\t", beta[0], finb[i], pnull, pmissing1notf);
	      for(index = 0; index < numalleles[j]; index++)
	      {
	      printf(" %4.3f", allelefreqs[j][index]);
	      }
	      }*/
	    if(k==l)
	      {
		if(k==99)
		  {
		    for(index=0;index<numalleles[j]-1;index++)
		      {
			Probgenos1[i][j][index][index] = (beta[0]*allelefreqs[j][index])/pmissing1notf;
			Probgenos0[i][j][index][index] = (beta[0]*pow(allelefreqs[j][index],2))/pmissing0notf;
			for(index2=index+1;index2<numalleles[j];index2++)
			  {
			    Probgenos0[i][j][index][index2]=(beta[0]*2*allelefreqs[j][index]*allelefreqs[j][index2])/pmissing0notf;
			    Probgenos1[i][j][index][index2]=0.0;
			  }//index2 loop
		      }//index loop
		    Probgenos0[i][j][numalleles[j]-1][numalleles[j]-1]=pow(pnull,2)/pmissing0notf; //prob of null given missing obs.
		    Probgenos1[i][j][numalleles[j]-1][numalleles[j]-1]=pnull/pmissing1notf; //prob of null given missing obs.
		    /*if((i == 17) && (j == 3) && (nmark == 5) && (numalleles[j] == 6) && (numind == 20))
		      {
		      printf("\n(1685) Probgenos1[i][j][null][null] = %4.3eX%4.3e/%4.3e", finb[i],pnull,pmissing1);
		      }*/
		    //Probgenos[i][j][numalleles[j]-1][numalleles[j]-1]=(finb[i]*pnull+(1-finb[i])*pow(pnull,2))/pmissing; //prob of null given null obs.

		    ExpX[i][j] = (finb[i]*(beta[0]+(1-beta[0])*pnull))/(beta[0]+(1-beta[0])*finb[i]*pnull+(1-beta[0])*(1-finb[i])*pow(pnull,2));

		    ExpB[i][j] = beta[0]/pmissing;

		  }//data is missing
		else
		  {
		    Probgenos0[i][j][k][k] = allelefreqs[j][k]/(allelefreqs[j][k]+2*pnull);
		    Probgenos1[i][j][k][k] = 1.0;
		    //Probgenos[i][j][k][k] = (finb[i]+(1-finb[i])*allelefreqs[j][k])/(finb[i]+(1-finb[i])*allelefreqs[j][k]+(1-finb[i])*2*pnull);
		    Probgenos0[i][j][k][numalleles[j]-1] = 2*pnull/(allelefreqs[j][k]+2*pnull);
		    Probgenos1[i][j][k][numalleles[j]-1] = 0.0; /* you can't be A_kA_n if your alleles are IBD */

		    ExpX[i][j] = (finb[i]*allelefreqs[j][k])/(finb[i]*allelefreqs[j][k]+(1-finb[i])*(pow(allelefreqs[j][k],2)+2*pnull*allelefreqs[j][k]));//allelefreqs[j][k?]

		    ExpB[i][j] = 0.0;
		    //if((i == 0) && (j < 10)) printf("\n (406) ind %d, marker %d geno %d%d, ExpX = %4.3f\n\t ExpX = %4.3f*%4.3f*%4.3f/%4.3f", i, j, myind1[i][j][0], myind1[i][j][1], ExpX[i][j],     (1-beta[0]),finb[i],allelefreqs[j][k],    ((1-beta[0])*finb[i]*allelefreqs[i][k]+(1-beta[0])*(1-finb[i])*(pow(allelefreqs[i][k],2)+2*pnull*allelefreqs[i][k])));
		  }//homozygous not missing
	      }
	    else
	      {
		Probgenos0[i][j][k][l]=1.0;
		Probgenos1[i][j][k][l]=0.0;
		Probgenos0[i][j][l][k]=1.0;
		Probgenos1[i][j][l][k]=0.0;

		ExpB[i][j] = 0.0;
		ExpX[i][j] = 0.0;

	      }//heterozygous
	    for(index=0;index<numalleles[j];index++)
	      {
		for(index2=index+1;index2<numalleles[j];index2++)
		  {
		    Probgenos0[i][j][index2][index]=Probgenos0[i][j][index][index2];
		    Probgenos1[i][j][index2][index]=Probgenos1[i][j][index][index2];
		  }
	      }//reflecting the matrix

	    /****debugging lines *****/
	    /*printf("\n(439) ind = %d, marker = %d; expX = %4.3f, expB = %4.3f", i, j, ExpX[i][j], ExpB[i][j]);
	      printf("\n Printing Probgenos:");

	      printf("\n\n (448) ind %d, marker %d, geno %d%d", i, j, k, l);
	      for(index = 0; index < numalleles[j]; index++)
	      {
	      printf("\n");
	      for(index2 = 0; index2 < numalleles[j]; index2++)
	      {
	      //sum += Probgenos0[i][j][index][index2];
	      printf(" %4.3f", Probgenos0[i][j][index][index2]);
	      }
	      }*/
	    /****end of debugging lines *****/
	    if((ExpX[i][j] < -0.0000000)|| (ExpX[i][j] > 1.000000001))
	      {
		printf("\n (662) i = %d, j = %d, ExpX = %4.6f (should be between 0 and 1). Exiting.\n", i, j, ExpX[i][j]);exit(0);
	      }
	    /* the next few lines check that Probgenos sums to 1.0*/
	     sum0= 0.0;
	     sum1= 0.0;
	     for(index = 0; index < numalleles[j]; index++)
	       {
		 for(index2 = index; index2 < numalleles[j]; index2++)
		   {
		     sum0 += Probgenos0[i][j][index][index2];
		     sum1 += Probgenos1[i][j][index][index2];
		   }
	       }//sum loop
	     if((sum0 < 0.99)|| (sum0 > 1.01))
	       {
		 printf("\n i = %d, j = %d, geno = %d%d, Probgenos0 sum = %4.6f (should be 1 but isn't)",i, j, myind1[i][j][0], myind1[i][j][1],sum0);
		 printf("\n (1729)Printing Probgenos0: geno is %d %d, numalleles = %d", myind1[i][j][0], myind1[i][j][1], numalleles[j]);
		 for(index = 0; index < numalleles[j]; index++)
		   {
		     printf("\n(466)");
		     for(index2 = 0; index2 < numalleles[j]; index2++)
		       {
			 sum += Probgenos0[i][j][index][index2];
			 printf(" %4.3f", Probgenos0[i][j][index][index2]);
		       }
		   }
		 printf("\n(595) Exiting\n");exit(0);
	       }//end if sum isn't 1
	     if((sum1 < 0.99)|| (sum1 > 1.01))
	       {
		 if(k == l)/* if heterozygous, p(x_ij = 1) = 0, so Probgenos1 is irrelevant. Hence, only look at homozygous case */
		   {
		     printf("\n (1745) i = %d, j = %d, geno = %d%d, Probgenos1 sum = %4.6f (should be 1 but isn't)\n\t numalleles = %d",i, j, myind1[i][j][0], myind1[i][j][1],sum1, numalleles[j]);
		     printf("\n Printing Probgenos1:");
		     for(index = 0; index < numalleles[j]; index++)
		       {
			 printf("\n(466)");
			 for(index2 = 0; index2 < numalleles[j]; index2++)
			   {
			     sum += Probgenos1[i][j][index][index2];
			     printf(" %4.3f", Probgenos1[i][j][index][index2]);
			   }
		       }
		 printf("\n(1756) Exiting \n");exit(0);
		   }
	       }//end if sum isn't 1
	     /* end of lines that check that Probgenos sums to 1.0*/

	  }/*for j*/
      }/*for i*/

    //printf("\n\t (429) done updating expected values");
    for(i=0;i<numind;i++)
      {
	sum = 0.0;
	for(j=0;j<nmark;j++)
	  {
	    sum = sum + ExpX[i][j];
	  }//for j

	finb[i] = sum/nmark;
      }//updating finb
    /*printf("\n(667) finb");
      for(i = 0; i < numind; i++)
      {
      printf(" %4.3f", finb[i]);
      }*/

	/***update beta****/

    sum = 0.0;
    for(i=0;i<numind;i++)
      {
	for(j=0;j<nmark;j++)
	  {
	    sum = sum + ExpB[i][j];
	  }
      }
    beta[0]= sum/(nmark*numind);
    //printf("\n (684) beta = %4.6f", beta[0]);

    /****updating allelefreqs****/

    for(j=0;j<nmark;j++)
      {
	sumbottom=0.0;
	for(i = 0; i < numind; i++)
	  {
	    for(k = 0; k < numalleles[j]; k++)
	      {
		sumbottom= sumbottom +Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j];
		for(index2 = k+1; index2 < numalleles[j]; index2++)
		  {
		    sumbottom = sumbottom + 2*Probgenos0[i][j][k][index2]*(1-ExpX[i][j]);
		  }
	      }
	  }

	for(k=0;k<numalleles[j];k++)
	  {
	    sumtop=0.0;
	    for(i=0;i<numind;i++)
	      {
		sumtop = sumtop + Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j];
		//if((i < 11) &&(Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j]>0.001) )printf("\n(794) i = %d, k = %d, index2 = %d, adding %4.3f", i, k, k, Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j]);
		for(index=0;index<numalleles[j];index++)
		  {
		    if(index!=k)
		      {
			sumtop = sumtop + Probgenos0[i][j][k][index]*(1-ExpX[i][j]);
			//if((i < 11)&& (Probgenos0[i][j][k][index]*(1-ExpX[i][j])>0.001))printf("\n(800) i = %d, k = %d, index2 = %d, adding %4.3f", i, k, index,  Probgenos0[i][j][k][index]*(1-ExpX[i][j]));
		      }
		  }//for index

		//printf("\n (534) ind %d. marker %d allele k = %d. (kk) %4.3f; (%4.3f, %4.3f, %4.3f)", i, j, k, Probgenos[i][j][k][k]*(2*(1-ExpX[i][j])+ExpX[i][j]),Probgenos[i][j][k][0], Probgenos[i][j][k][1],Probgenos[i][j][k][2]);
		//printf("\n\t (468) ExpX[i][j] = %4.3f", ExpX[i][j]);

	      }//for i
	    allelefreqs[j][k] = sumtop/sumbottom;
	  }//for k
	/*printf("\n (544)marker %d", j);
	  for(k = 0; k < numalleles[j]; k++)
	  {
	  printf(" %4.3f", allelefreqs[j][k]);
	  }*/
	sum = 0.0;
	for(k=0;k<numalleles[j];k++)
	  {
	    sum= sum + allelefreqs[j][k];
	  }
	if((sum < 0.99) || (sum > 1.01))
	  {
	    printf("\n(483) j = %d sum=%4.3f\n\tsumtop = %4.3f, sumbottom = %4.3f", j, sum, sumtop, sumbottom);
	    printf("\n(819)\n");exit(0);
	  }
      }//for j (loop through markers in which we update the allele freqs.*/
    /*for(i = 0; i < nmark; i++)
      {
      printf("\n(1884) marker %d: ", i);
      for(j =0; j < numalleles[i]; j++)
      {
      printf(" %4.4f", allelefreqs[i][j]);
      }
      }*/

    likval = nulllikelihood(myind1, allelefreqs, beta, finb, numalleles, numind, nmark);
  } /* end while loop (for for now) */


  /************* end of new stuff ***************/

  /* freeing memory *****/
  for(i = 0; i < numind; i++)
    {
      free(ExpX[i]);
      free(ExpB[i]);
    }
  free(ExpX);
  free(ExpB);

  for(i = 0; i < numind; i++)
    {
      for(markindex = 0; markindex < nmark; markindex++)
        {
          for(j = 0; j < numalleles[markindex]; j++)
            {
              free(Probgenos0[i][markindex][j]);
              free(Probgenos1[i][markindex][j]);
            }
          free(Probgenos0[i][markindex]);
          free(Probgenos1[i][markindex]);
        }
      free(Probgenos0[i]);
      free(Probgenos1[i]);
    }
    free(Probgenos0);
    free(Probgenos1);
      return;
//return myf;*/
}
/****************************************************************************/
/* Kalinowski and Taper (2006) Conservation Genetics 7:991--995.*/
/****************************************************************************/
  void Kalinowski(int*** myind1, double** kallelefreqs, int* numalleles, double* kbeta, int nmark, int numind)
{
  int i, j, k, l;
  int loopindex;
  double** oldfreqs;
  double* betaold;
  int*** nmatrix;/* holds counts of individuals with each genotype */
  double sumjk;/* the number of individuals heterozygous for allele k*/
  double specialsum;
  double pnull, nmiss;/* shorthand for null allele freq and number missing people*/
  double likval, oldlikval;

  /* allocating memory*/
  betaold = (double*)malloc(nmark*sizeof(double));
  oldfreqs = (double**)malloc(nmark*sizeof(double*));
  for(j = 0; j < nmark; j++)
    {
      oldfreqs[j] = (double*)malloc(numalleles[j]*sizeof(double));
    }
  nmatrix = (int***)malloc(nmark*sizeof(int**));
  for(j = 0; j < nmark; j++)
    {
      nmatrix[j] = (int**)malloc((numalleles[j]+1)*sizeof(int*));
      for(k = 0; k < numalleles[j]+1; k++)
	{
	  nmatrix[j][k] = (int*)malloc((numalleles[j]+1)*sizeof(int));
	  for(l = 0; l < numalleles[j]+1; l++)
	    {
	      nmatrix[j][k][l] = 0;
	    }
	}
    }

  /* calculating nmatrix*/
  for(j = 0; j < nmark; j++)
    {
      for(i = 0; i < numind; i++)
	{
	  if(myind1[i][j][0] == 99)
	    {
	      nmatrix[j][numalleles[j]][numalleles[j]]++;
	    }
	  else
	    {
	      if(myind1[i][j][0] == myind1[i][j][1])
		{
		  nmatrix[j][myind1[i][j][0]][myind1[i][j][0]]++;
		}
	      else
		{
		 nmatrix[j][myind1[i][j][0]][myind1[i][j][1]]++;
		 nmatrix[j][myind1[i][j][1]][myind1[i][j][0]]++;
		}
	    }
	}
    }
  /* initial parameter estimates */
  /*for(j = 0; j < nmark; j++)
    {
    kbeta[j] = 0.05;
    }
    for(j = 0; j < nmark; j++)
    {
    for(k = 0; k < numalleles[j]; k++)
    {
    kallelefreqs[j][k] = 1/numalleles[j];
    }
    }*/
  /*printf("\n\n(652) starting values");
    for(j = 0; j < nmark; j++)
    {
    printf("\n marker %d; beta = %4.4f; allelefreqs =", j, kbeta[j]);
    for(k = 0; k < numalleles[j]; k++)
    {
    printf(" %4.4f", kallelefreqs[j][k]);
    }
    }*/


  /* main part of function */
  likval=1.0;
  oldlikval=2.0;
loopindex = 0;
  while((fabs(likval-oldlikval)> ktol) && (loopindex < MAXNUMTIMES))
  {
    loopindex = loopindex + 1;
    oldlikval = likval;
    // printf("\n \t(359) %dth iteration through KT EM algorithm. like = %4.4f", loopindex, likval);fflush(stdout);

      /* set old values */
      for(j = 0; j < nmark; j++)
	{
	  betaold[j] = kbeta[j];
	}

      for(j = 0; j < nmark; j++)
	{
	  for(k = 0; k < numalleles[j]; k++)
	    {
	      oldfreqs[j][k] = kallelefreqs[j][k];
	    }
	}
      /* old values set.(finished)*/

      /* update */

      for(j = 0; j < nmark; j++)
	{
	  /* update allele freqs */
	  pnull = oldfreqs[j][numalleles[j]-1];
	  nmiss = nmatrix[j][numalleles[j]][numalleles[j]];
	  for(k = 0; k < numalleles[j]-1; k++)
	    {
	      /* we need the number of individuals that hare heterozygous*/
	      /* for allele k*/
	      sumjk = 0.0;
	      for(l = 0; l < numalleles[j]-1; l++)
		{
		  if(l != k){sumjk += nmatrix[j][k][l];}
		}
	      kallelefreqs[j][k] = (1/(2*(double)numind))*(2*nmatrix[j][k][k]*(oldfreqs[j][k] + pnull)/(oldfreqs[j][k] + 2*pnull) + sumjk + 2*nmiss*betaold[j]*oldfreqs[j][k]/(betaold[j] + (1-betaold[j])*pow(pnull, 2)));
	    }
	  specialsum = 0.0;
	  for(k = 0; k < numalleles[j]-1; k++)
	    {
	      specialsum += nmatrix[j][k][k]*2*pnull/(oldfreqs[j][k] + 2*pnull);
	    }

	  kallelefreqs[j][numalleles[j]-1] = (1/(2*(double)numind))*(specialsum+ 2*nmiss*(pow(pnull, 2)*(1-betaold[j]) + betaold[j]*pnull)/(betaold[j] + pow(pnull, 2)*(1-betaold[j])));
	  /* update beta */
	  kbeta[j] = (1/(double)numind)*nmiss*betaold[j]/(betaold[j] + pow(pnull, 2)*(1-betaold[j]));
	  //printf("\n(706) oldpnull = %4.4f, newpnull = (%4.4e) = %4.4e", betaold[j], (1/(2*(double)numind))*( specialsum+ 2*nmiss*(pow(pnull, 2)*(1-betaold[j]) + betaold[j]*pnull)/(betaold[j] + pow(pnull, 2)*(1-betaold[j]))), kallelefreqs[j][numalleles[j]-1]);
	}/* end loop through markers */

    likval = klikelihood(myind1, kallelefreqs, kbeta, numalleles, nmark, numind);
    }/* end main EM loop */


  /*for(j = 0; j < nmark; j++)
    {
    printf("\n marker %d; beta = %4.4f; allelefreqs =", j, kbeta[j]);
    for(k = 0; k < numalleles[j]; k++)
    {
    printf(" %4.4f", kallelefreqs[j][k]);
    }
    }
    printf("\n(710)\n");exit(0);*/

  /* freeing */
  for(j = 0; j < nmark; j++)
    {
      free(oldfreqs[j]);
    }
  free(oldfreqs);
  free(betaold);
  return;
}
/****************************************************************/

/*************************************************************************/
/* An EM algorithm that jointly estimates allele frequencies (including) */
/* null allele frequencies) and individual inbreeding coefficients.      */
/* in addition to missingness due to homozygosity for null alleles, there*/
/* is some probability that any genotype can be simply randomly missing. */
/* In this model, each marker has its own probability that its genotypes */
/* will be randomly missing.                                             */
/*************************************************************************/
void betaLaina(int*** myind1, double** allelefreqs, int* numalleles, double* finb, double* beta, int nmark, int numind, int* sex)
{
  int i, j, l, a, k, m, markindex, loopindex, index, index2, index3, index4;

  double sum, sum0, sum1, sumtop, sumbottom, likval, oldlikval;
  double** ExpX; /*Exp[numind][nmark]*/
  /* probgenos is the prob that an individual has each genotype at each */
  /* marker. We need this to be conditional on whether the two alleles are*/
  /* IBD or not, hence the 0 (not IBD) and the 1 (IBD). So */
  /* Probgenos0 = Pr(genotype, xij|obs data, current param values)*/
  double**** Probgenos0; /*Probgenos[numind][nmark][numalleles][numalleles]*/
  double**** Probgenos1; /*Probgenos[numind][nmark][numalleles][numalleles]*/
  double** ExpB;/*ExpB[numind][nmark]*/

  /* pmissing = Pr(observed = missing|data). pmissing0 = Pr(obs = missing and notIBD|data), pmissing1 = Pr(obs = missing and IBD| data) */
  double pnull, pmissing;
  //double pmissing0, pmissing1;
  double pmissing1notf, pmissing0notf; /* pmissing1 is a product of f times something. pmissing1notf is that something.*/

  /* allocating memory *****/

  ExpX = (double**)malloc(numind*sizeof(double*));
  ExpB = (double**)malloc(numind*sizeof(double*));
  for(i = 0; i < numind; i++)
    {
      ExpX[i] = (double*)malloc(nmark*sizeof(double));
      ExpB[i] = (double*)malloc(nmark*sizeof(double));
    }
  Probgenos0 = (double****)malloc(numind*sizeof(double***));
  Probgenos1 = (double****)malloc(numind*sizeof(double***));
  if(Probgenos0 == NULL)
    {
      printf("\n(488) Problem allocating memory for Probgenos0. Exiting.\n");exit(0);
    };
  if(Probgenos1 == NULL)
    {
      printf("\n(497) Problem allocating memory for Probgenos1. Exiting.\n");exit(0);
    }

  for(i = 0; i < numind; i++)
    {
      Probgenos0[i] = (double***)malloc(nmark*sizeof(double**));
      Probgenos1[i] = (double***)malloc(nmark*sizeof(double**));
      if(Probgenos0[i] == NULL)
	{
	  printf("\n(506) Problem allocating memory for Probgenos0[i = %d]. Exiting.\n", i);exit(0);
	}
      if(Probgenos0[i] == NULL)
	{
	  printf("\n(510) Problem allocating memory for Probgenos1[i = %d]. Exiting.\n", i);exit(0);
	}
      for(markindex = 0; markindex < nmark; markindex++)
        {
          Probgenos0[i][markindex] = (double**)malloc(numalleles[markindex]*sizeof(double*));
          Probgenos1[i][markindex] = (double**)malloc(numalleles[markindex]*sizeof(double*));
	  if(Probgenos0[i][markindex] == NULL)
	    {
	      printf("\n(518) Problem allocating memory for Probgenos0[i = %d][m = %d]. Exiting.\n", i, markindex);exit(0);
	    }
	  if(Probgenos1[i][markindex] == NULL)
	    {
	      printf("\n(518) Problem allocating memory for Probgenos1[i = %d][m = %d]. Exiting.\n", i, markindex);exit(0);
	    }
          for(j = 0; j < numalleles[markindex]; j++)
            {
              Probgenos0[i][markindex][j] = (double*)malloc(numalleles[markindex]*sizeof(double));
              Probgenos1[i][markindex][j] = (double*)malloc(numalleles[markindex]*sizeof(double));
	      if(Probgenos0[i][markindex][j] == NULL)
		{
		  printf("\n(530) Problem allocating memory for Probgenos0[i = %d][m = %d][j = %d]. Exiting.\n", i, markindex, j);exit(0);
		}
	      if(Probgenos1[i][markindex][j] == NULL)
		{
		  printf("\n(534) Problem allocating memory for Probgenos1[i = %d][m = %d][j = %d]. Exiting.\n", i, markindex, j);exit(0);
		}
              for(k = 0; k < numalleles[markindex]; k++)
                {
                  Probgenos0[i][markindex][j][k] = 0.0;
                  Probgenos1[i][markindex][j][k] = 0.0;
                }
            }
        }
    }

  /*printf("\n(498) WARNING: I'm resetting the starting values as a debugging\nexercise\n");fflush(stdout);
    for(j = 0; j < nmark; j++)
    {
    allelefreqs[j][0] = 0.5;
    allelefreqs[j][1] = 0.3;
    allelefreqs[j][2] = 0.2;
    }
    for(i = 0; i < numind; i++)
    {
    finb[i] = 0.2;
    }
    beta[0] = 0.1;*/

  likval=1.0;
  oldlikval=2.0;

  /* Next there will be a giant while loop that runs the EM algorithm. */
  /* When you are first writing this, make it a for loop*/
  //for(loopindex = 0; loopindex < 10; loopindex++)
  loopindex = 0;
  while((fabs(likval-oldlikval)>tolerance) && (loopindex < MAXNUMTIMES))
  {
    loopindex = loopindex + 1;
    oldlikval = likval;

    for(i=0;i<numind;i++)
      {
	for(j=0;j<nmark;j++)
	  {
	    k=myind1[i][j][0];
	    l=myind1[i][j][1];
	    pnull=allelefreqs[j][numalleles[j]-1];
	    pmissing=beta[j] + (1-beta[j])*(finb[i]*pnull+(1-finb[i])*pow(pnull,2));
	    /* pmissing0 and pmissing1 are Pr(geno = M and notIBD) and Pr(geno = M and IBD), respectively*/
	    //pmissing0 = beta[j]*(1-finb[i]) + (1-beta[j])*(1-finb[i])*pow(pnull,2);
	    //pmissing1 = beta[j]*(finb[i]) + (1-beta[j])*finb[i]*pnull;
	    /* the probability that an individual's genotype is missing, */
	    /* given that their alleles are IBD (ie pmissing1) has a factor*/
	    /* of finb[i] in it. When we use pmissing1 to find values in */
	    /*Probgenos1, the probgenos value has the form:                */
	    /* finb[i]*something/pmissing1. We are having problems with this*/
	    /* when finb[i] is converging to zero. The solution is to cancel*/
	    /* the finb[i] term out of top and bottom of the ratio. Hence,*/
	    /* instead of pmissing1, we want pmissing1 without the finb[i].*/
	    /* that is, we want pmissing1notf.*/
	    pmissing1notf = beta[j] + (1-beta[j])*pnull;
	    pmissing0notf =  beta[j] + (1-beta[j])*pow(pnull,2);
	    /*if((i == 1) && (j == 2) && (nmark == 5) && (numalleles[j] == 6) && (numind == 20))
	      {
	      printf("\n(1666) beta = %4.4e, finb = %4.4e, pnull = %4.ef, pmissing1notf = %4.3e\n\t", beta[j], finb[i], pnull, pmissing1notf);
	      for(index = 0; index < numalleles[j]; index++)
	      {
	      printf(" %4.3f", allelefreqs[j][index]);
	      }
	      }*/
        if(j >= numNonXLinked && sex[i]==1)
            {
              if(k==99)
                {
                  for(index=0;index<numalleles[j]-1;index++)
          		      {
          			        Probgenos1[i][j][index][index] = (beta[j]*allelefreqs[j][index])/pmissing1notf;
                        Probgenos0[i][j][index][index] = (beta[j]*allelefreqs[j][index])/pmissing1notf;
          		      }//index loop
                  Probgenos1[i][j][numalleles[j]-1][numalleles[j]-1]=pnull/pmissing1notf;
                  Probgenos0[i][j][numalleles[j]-1][numalleles[j]-1]=pnull/pmissing1notf;
                  ExpX[i][j] = 1.0;
                  ExpB[i][j] = beta[j]/pmissing1notf;
                }
              else
                {
                  Probgenos1[i][j][k][k] = 1.0;
                  Probgenos0[i][j][k][k] = 1.0;
                  ExpX[i][j] = 1.0;
                  ExpB[i][j] = 0;
                }
            }
        else
            {
	    if(k==l)
	      {
		if(k==99)
		  {
		    for(index=0;index<numalleles[j]-1;index++)
		      {
			Probgenos1[i][j][index][index] = (beta[j]*allelefreqs[j][index])/pmissing1notf;
			Probgenos0[i][j][index][index] = (beta[j]*pow(allelefreqs[j][index],2))/pmissing0notf;
			for(index2=index+1;index2<numalleles[j];index2++)
			  {
			    Probgenos0[i][j][index][index2]=(beta[j]*2*allelefreqs[j][index]*allelefreqs[j][index2])/pmissing0notf;
			    Probgenos1[i][j][index][index2]=0.0;
			  }//index2 loop
		      }//index loop
		    Probgenos0[i][j][numalleles[j]-1][numalleles[j]-1]=pow(pnull,2)/pmissing0notf; //prob of null given missing obs.
		    Probgenos1[i][j][numalleles[j]-1][numalleles[j]-1]=pnull/pmissing1notf; //prob of null given missing obs.
		    /*if((i == 17) && (j == 3) && (nmark == 5) && (numalleles[j] == 6) && (numind == 20))
		      {
		      printf("\n(1685) Probgenos1[i][j][null][null] = %4.3eX%4.3e/%4.3e", finb[i],pnull,pmissing1);
		      }*/
		    //Probgenos[i][j][numalleles[j]-1][numalleles[j]-1]=(finb[i]*pnull+(1-finb[i])*pow(pnull,2))/pmissing; //prob of null given null obs.

		    ExpX[i][j] = (finb[i]*(beta[j]+(1-beta[j])*pnull))/(beta[j]+(1-beta[j])*finb[i]*pnull+(1-beta[j])*(1-finb[i])*pow(pnull,2));

		    ExpB[i][j] = beta[j]/pmissing;

		  }//data is missing
		else
		  {
		    Probgenos0[i][j][k][k] = allelefreqs[j][k]/(allelefreqs[j][k]+2*pnull);
		    Probgenos1[i][j][k][k] = 1.0;
		    //Probgenos[i][j][k][k] = (finb[i]+(1-finb[i])*allelefreqs[j][k])/(finb[i]+(1-finb[i])*allelefreqs[j][k]+(1-finb[i])*2*pnull);
		    Probgenos0[i][j][k][numalleles[j]-1] = 2*pnull/(allelefreqs[j][k]+2*pnull);
		    Probgenos1[i][j][k][numalleles[j]-1] = 0.0; /* you can't be A_kA_n if your alleles are IBD */

		    ExpX[i][j] = (finb[i]*allelefreqs[j][k])/(finb[i]*allelefreqs[j][k]+(1-finb[i])*(pow(allelefreqs[j][k],2)+2*pnull*allelefreqs[j][k]));//allelefreqs[j][k?]

		    ExpB[i][j] = 0.0;
		    //if((i == 0) && (j < 10)) printf("\n (406) ind %d, marker %d geno %d%d, ExpX = %4.3f\n\t ExpX = %4.3f*%4.3f*%4.3f/%4.3f", i, j, myind1[i][j][0], myind1[i][j][1], ExpX[i][j],     (1-beta[0]),finb[i],allelefreqs[j][k],    ((1-beta[0])*finb[i]*allelefreqs[i][k]+(1-beta[0])*(1-finb[i])*(pow(allelefreqs[i][k],2)+2*pnull*allelefreqs[i][k])));
		  }//homozygous not missing
	      }
	    else
	      {
		Probgenos0[i][j][k][l]=1.0;
		Probgenos1[i][j][k][l]=0.0;
		Probgenos0[i][j][l][k]=1.0;
		Probgenos1[i][j][l][k]=0.0;

		ExpB[i][j] = 0.0;
		ExpX[i][j] = 0.0;

	      }//heterozygous
      }
	    for(index=0;index<numalleles[j];index++)
	      {
		for(index2=index+1;index2<numalleles[j];index2++)
		  {
		    Probgenos0[i][j][index2][index]=Probgenos0[i][j][index][index2];
		    Probgenos1[i][j][index2][index]=Probgenos1[i][j][index][index2];
		  }
	      }//reflecting the matrix


	    /****debugging lines *****/
	    /*printf("\n(439) ind = %d, marker = %d; expX = %4.3f, expB = %4.3f", i, j, ExpX[i][j], ExpB[i][j]);
	      printf("\n Printing Probgenos:");

	      printf("\n\n (448) ind %d, marker %d, geno %d%d", i, j, k, l);
	      for(index = 0; index < numalleles[j]; index++)
	      {
	      printf("\n");
	      for(index2 = 0; index2 < numalleles[j]; index2++)
	      {
	      //sum += Probgenos0[i][j][index][index2];
	      printf(" %4.3f", Probgenos0[i][j][index][index2]);
	      }
	      }*/
	    /****end of debugging lines *****/
	    if((ExpX[i][j] < -0.0000000)|| (ExpX[i][j] > 1.000000001))
	      {
		printf("\n (662) i = %d, j = %d, ExpX = %4.6f (should be between 0 and 1). Exiting.\n", i, j, ExpX[i][j]);exit(0);
	      }
	    /* the next few lines check that Probgenos sums to 1.0*/
	     sum0= 0.0;
	     sum1= 0.0;
	     for(index = 0; index < numalleles[j]; index++)
	       {
		 for(index2 = index; index2 < numalleles[j]; index2++)
		   {
		     sum0 += Probgenos0[i][j][index][index2];
		     sum1 += Probgenos1[i][j][index][index2];
		   }
	       }//sum loop
	     if((sum0 < 0.99)|| (sum0 > 1.01))
	       {
		 printf("\n i = %d, j = %d, geno = %d%d, Probgenos0 sum = %4.6f (should be 1 but isn't)",i, j, myind1[i][j][0], myind1[i][j][1],sum0);
		 printf("\n (1729)Printing Probgenos0: geno is %d %d, numalleles = %d", myind1[i][j][0], myind1[i][j][1], numalleles[j]);
		 for(index = 0; index < numalleles[j]; index++)
		   {
		     printf("\n(466)");
		     for(index2 = 0; index2 < numalleles[j]; index2++)
		       {
			 sum += Probgenos0[i][j][index][index2];
			 printf(" %4.3f", Probgenos0[i][j][index][index2]);
		       }
		   }
		 printf("\n(595) Exiting\n");exit(0);
	       }//end if sum isn't 1
	     if((sum1 < 0.99)|| (sum1 > 1.01))
	       {
		 if(k == l)/* if heterozygous, p(x_ij = 1) = 0, so Probgenos1 is irrelevant. Hence, only look at homozygous case */
		   {
		     printf("\n (1745) i = %d, j = %d, geno = %d%d, Probgenos1 sum = %4.6f (should be 1 but isn't)\n\t numalleles = %d",i, j, myind1[i][j][0], myind1[i][j][1],sum1, numalleles[j]);
		     printf("\n Printing Probgenos1:");
		     for(index = 0; index < numalleles[j]; index++)
		       {
			 printf("\n(466)");
			 for(index2 = 0; index2 < numalleles[j]; index2++)
			   {
			     sum += Probgenos1[i][j][index][index2];
			     printf(" %4.3f", Probgenos1[i][j][index][index2]);
			   }
		       }
		 printf("\n(1756) Exiting \n");exit(0);
		   }
	       }//end if sum isn't 1
	     /* end of lines that check that Probgenos sums to 1.0*/

	  }/*for j*/
      }/*for i*/

    //printf("\n\t (429) done updating expected values");
    for(i=0;i<numind;i++)
      {
	sum = 0.0;
	for(j=0;j<nmark;j++)
	  {
      if(j >= numNonXLinked && sex[i]==1)
          {}
      else
        {
	         sum = sum + ExpX[i][j];
        }
	  }//for j
    if(sex[i]==1){
      finb[i] = sum/numNonXLinked;
    }
    else{
	     finb[i] = sum/nmark;
     }
      }//updating finb
    /*printf("\n(667) finb");
      for(i = 0; i < numind; i++)
      {
      printf(" %4.3f", finb[i]);
      }*/

	/***update beta****/

    for(j=0;j<nmark;j++)
      {
	sum = 0.0;
	for(i=0;i<numind;i++)
	  {
	    sum = sum + ExpB[i][j];
	  }
	beta[j]= sum/(numind);
      }
    //printf("\n (684) beta = %4.6f", beta[0]);

    /****updating allelefreqs****/

    for(j=0;j<nmark;j++)
      {
	sumbottom=0.0;
	for(i = 0; i < numind; i++)
	  {
	    for(k = 0; k < numalleles[j]; k++)
	      {
		sumbottom= sumbottom +Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j];
		for(index2 = k+1; index2 < numalleles[j]; index2++)
		  {
		    sumbottom = sumbottom + 2*Probgenos0[i][j][k][index2]*(1-ExpX[i][j]);
		  }
	      }
	  }

	for(k=0;k<numalleles[j];k++)
	  {
	    sumtop=0.0;
	    for(i=0;i<numind;i++)
	      {
		sumtop = sumtop + Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j];
		//if((i < 11) &&(Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j]>0.001) )printf("\n(794) i = %d, k = %d, index2 = %d, adding %4.3f", i, k, k, Probgenos0[i][j][k][k]*2*(1-ExpX[i][j]) + Probgenos1[i][j][k][k]*ExpX[i][j]);
		for(index=0;index<numalleles[j];index++)
		  {
		    if(index!=k)
		      {
			sumtop = sumtop + Probgenos0[i][j][k][index]*(1-ExpX[i][j]);
			//if((i < 11)&& (Probgenos0[i][j][k][index]*(1-ExpX[i][j])>0.001))printf("\n(800) i = %d, k = %d, index2 = %d, adding %4.3f", i, k, index,  Probgenos0[i][j][k][index]*(1-ExpX[i][j]));
		      }
		  }//for index

		//printf("\n (534) ind %d. marker %d allele k = %d. (kk) %4.3f; (%4.3f, %4.3f, %4.3f)", i, j, k, Probgenos[i][j][k][k]*(2*(1-ExpX[i][j])+ExpX[i][j]),Probgenos[i][j][k][0], Probgenos[i][j][k][1],Probgenos[i][j][k][2]);
		//printf("\n\t (468) ExpX[i][j] = %4.3f", ExpX[i][j]);

	      }//for i
	    allelefreqs[j][k] = sumtop/sumbottom;
	  }//for k
	/*printf("\n (544)marker %d", j);
	  for(k = 0; k < numalleles[j]; k++)
	  {
	  printf(" %4.3f", allelefreqs[j][k]);
	  }*/
	sum = 0.0;
	for(k=0;k<numalleles[j];k++)
	  {
	    sum= sum + allelefreqs[j][k];
	  }
	if((sum < 0.99) || (sum > 1.01))
	  {
	    printf("\n(483) j = %d sum=%4.3f\n\tsumtop = %4.3f, sumbottom = %4.3f", j, sum, sumtop, sumbottom);
	    printf("\n(819)\n");exit(0);
	  }
      }//for j (loop through markers in which we update the allele freqs.*/
    /*for(i = 0; i < nmark; i++)
      {
      printf("\n(1884) marker %d: ", i);
      for(j =0; j < numalleles[i]; j++)
      {
      printf(" %4.4f", allelefreqs[i][j]);
      }
      }*/

    likval = betalikelihood(myind1, allelefreqs, beta, finb, numalleles, nmark, numind);
  } /* end while loop (for for now) */

  /************* end of new stuff ***************/

  /* freeing memory *****/
  for(i = 0; i < numind; i++)
    {
      free(ExpX[i]);
      free(ExpB[i]);
    }
  free(ExpX);
  free(ExpB);

  for(i = 0; i < numind; i++)
    {
      for(markindex = 0; markindex < nmark; markindex++)
        {
          for(j = 0; j < numalleles[markindex]; j++)
            {
              free(Probgenos0[i][markindex][j]);
              free(Probgenos1[i][markindex][j]);
            }
          free(Probgenos0[i][markindex]);
          free(Probgenos1[i][markindex]);
        }
      free(Probgenos0[i]);
      free(Probgenos1[i]);
    }
    free(Probgenos0);
    free(Probgenos1);
      return;
//return myf;*/
}
/*************************************************************************/
/****************************************************/
/* The following is vogl's algorithm for jointly    */
/* estimating inbreeding coefficients and allele */
/* frequencies (including null alleles).*/
/* vogl et al. (2002) J Evol Bio 15 433--439.*/
/*****************************************************/
void vogl(int*** genos, double** allelefreqs, double* voglfinblist, int* numalleles, int numind, int nmark, double myalpha, double mybeta)/* Vogl's (2002) algorithm */
{
  void voglupdatefinb(int*** genos, int*** tempgenos, double*** tempfreqs, double** tempfinb, int gibbsindex, int** xij, double myalpha, double mybeta, int* numalleles, int numind, int nmark);
  int i, j, k, gibbsindex, markindex;
  int** xij;
  double x;
  double **tempfinb;
  double ***tempfreqs;
  double *alpha;
  double sumx, sumy;
  double *sumZ, *sumfinb, *allelesum;
  int *nummiss;
  int ***tempgenos;
  double sum;
  double tempZ;

  //tempout = fopen("tempout.txt", "w");
  /* memory allocation.*/
  tempgenos = (int***)malloc(numind*sizeof(int**));
  for(i = 0; i < numind; i++)
    {
      tempgenos[i] = (int**)malloc(nmark*sizeof(int*));
      for(j = 0; j < nmark; j++)
        {
          tempgenos[i][j] = (int*)malloc(2*sizeof(int));
          tempgenos[i][j][0] = genos[i][j][0];
          tempgenos[i][j][1] = genos[i][j][1];
        }
    }
  /* nummiss records how many markers each individual has a missing genotype for*/
  nummiss = (int*)malloc(numind*sizeof(int));
  for(i = 0; i < numind; i++)
    {
      nummiss[i] = 0;
      for(j = 0; j < nmark; j++)
	{
	  if(genos[i][j][0] == 99){nummiss[i]++;}
	}
    }

  sumfinb = (double*)malloc(numind*sizeof(double));
  /*xij are indicators: 1 if individual i's alleles are IBD at marker j, 0 else*/
  xij = (int**)malloc(numind*sizeof(int*));
  for(i = 0; i < numind; i++)
    {
      xij[i] = (int*)malloc(nmark*sizeof(int));
      for(j = 0; j < nmark; j++)
	{
	  xij[i][j] = 0;
	}
    }
  /* tempfinb keeps track of all the finb values from our Gibbs sampler*/
  tempfinb = (double**)malloc(numgibbs*sizeof(double*));
  for(gibbsindex = 0; gibbsindex <= numgibbs; gibbsindex++)
    {
      tempfinb[gibbsindex] = (double*)malloc(numind*sizeof(double));
    }
  /* initializing tempfinb (0th iteration)*/
  for(i = 0; i < numind; i++)
    {
      x = drand48();
      tempfinb[0][i] = x/10.0;
    }

  /* tempfreqs keeps track of all our allele frequency estimates */
   tempfreqs = (double***)malloc(numgibbs*sizeof(double**));
  for(gibbsindex = 0; gibbsindex <= numgibbs; gibbsindex++)
    {
      tempfreqs[gibbsindex] = (double**)malloc(nmark*sizeof(double*));
      for(j = 0; j < nmark; j++)
	{
	  tempfreqs[gibbsindex][j] = (double*)malloc(numalleles[j]*sizeof(double));
	  for(k = 0; k < numalleles[j]; k++)
	    {
	      tempfreqs[gibbsindex][j][k] = 0.0;
	    }
	}
    }

  /* initializing tempfreqs (just values for 0th iteration)*/
  for(j = 0; j < nmark; j++)
    {
      alpha =(double*)malloc(numalleles[j]*sizeof(double));
      for(k = 0; k < numalleles[j]; k++)
	{
	  alpha[k] = 1.0;
	}
      rdirichlet(alpha, tempfreqs[0][j], numalleles[j]);
      free(alpha);
    }

  /*for(j = 0; j < nmark; j++)
    {
    printf("\n(556) marker %d: ", j);
    for(k = 0; k < numalleles[j]; k++)
    {
    printf(" %4.4f", tempfreqs[0][j][k]);
    }
    }*/
  /* done with mallocing and initialization */

  if(burnin > gibbsindex)
    {
      printf("\n the program requires that the number of burnin runs be less than the number of cycles through the gibbs sampler. Exiting.\n");exit(0);
    }

  /* burnin loop:*/
  for(gibbsindex = 0; gibbsindex < burnin; gibbsindex++)
    {
      /* updating xij*/
      voglupdatefinb(genos, tempgenos, tempfreqs, tempfinb, gibbsindex, xij, myalpha, mybeta, numalleles, numind, nmark);/* also sets the xij values and tempgenos*/

      /* updating allele frequencies */
      for(j = 0; j < nmark; j++)
	{
	  sumZ = (double*)malloc(numalleles[j]*sizeof(double));
	  for(k = 0; k < numalleles[j]; k++)
	    {
	      sumZ[k] = 0.0;
	    }
	  /* sumZ records how many of each allele we see in the data */
	  for(i = 0; i < numind; i++)
	    {
	      if(xij[i][j] == 1)
		{
		  sumZ[tempgenos[i][j][0]]++; /* two IBD alleles, but count as 1*/
		}
	      else
		{
		  sumZ[tempgenos[i][j][0]]++;
		  sumZ[tempgenos[i][j][1]]++;
		}
	    }
	  /* the parameters of the dirichlet distributions are 1 more than the counts */
	  for(i = 0; i < numalleles[j]; i++)
	    {
	      sumZ[i] += 1;
	    }
	  /* the new allele frequencies will be a random draw from */
	  /* a dirichlet distribution with a mean of the old frequencies.*/
	  rdirichlet(sumZ, tempfreqs[gibbsindex+1][j], numalleles[j]);
	  /*printf("\n(628) marker %d: ", j);
	    for(k = 0; k < numalleles[j]; k++)
	    {
	    printf(" %4.3f",sumZ[k]);
	    }
	    printf("\n(633) marker %d: ", j);
	    for(k = 0; k < numalleles[j]; k++)
	    {
	    printf(" %4.3f",tempfreqs[gibbsindex+1][j][k]);
	    }*/
	  free(sumZ);
	}/* end loop through markers*/
      //printf("\n(1884)\n");exit(0);
    }/* done with burnin samples */

  /* the last of the burnin samples will be the start of the real samples*/
  for(j = 0; j < nmark; j++)
    {
      for(k = 0; k < numalleles[j]; k++)
	{
	  tempfreqs[0][j][k] = tempfreqs[burnin][j][k];
	}
    }
  for(i = 0; i < numind; i++)
    {
      tempfinb[0][i] = tempfinb[burnin-1][i];
    }

  /* main loop:*/
  for(gibbsindex = 0; gibbsindex < numgibbs; gibbsindex++)
    {
      /* updating xij*/
      voglupdatefinb(genos, tempgenos, tempfreqs, tempfinb, gibbsindex, xij, myalpha, mybeta, numalleles, numind, nmark);/* also sets the xij values and tempgenos*/

   /* updating allele frequencies */
      for(j = 0; j < nmark; j++)
	{
	  sumZ = (double*)malloc(numalleles[j]*sizeof(double));
	  for(k = 0; k < numalleles[j]; k++)
	    {
	      sumZ[k] = 0.0;
	    }
	  /* sumZ records how many of each allele we see in the data */
	  for(i = 0; i < numind; i++)
	    {
	      if(xij[i][j] == 1)
		{
		  sumZ[tempgenos[i][j][0]]= sumZ[tempgenos[i][j][0]]+1.0; /* two IBD alleles, but count as 1*/
		}
	      else
		{
		  sumZ[tempgenos[i][j][0]]++;
		  sumZ[tempgenos[i][j][1]]++;
		}
	      /*if(j == 99)
		{
		printf("\n(2130) xij[%d][99] = %d. tempgenos = %d%d. sumZ[9] = %4.0f", i, xij[i][j], tempgenos[i][j][0], tempgenos[i][j][1], sumZ[9]);
		}*/
	    }
	    tempZ = sumZ[numalleles[j]-1];
	  for(i = 0; i < numalleles[j]; i++)
	    {
	      sumZ[i] += 1;
	    }
	  rdirichlet(sumZ, tempfreqs[gibbsindex+1][j], numalleles[j]);
	  free(sumZ);
	}/* end loop through markers*/
    }/* done with main loop */

  /* Now we have our gibbs sample from the unconditional distribution. We */
  /* need to summarize it. First, each person's estimated inbreeding */
  /* coefficient is the mean of the estimates from each sample.*/
  for(i = 0; i < numind; i++){sumfinb[i] = 0.0;}
    for(gibbsindex = 0; gibbsindex < numgibbs; gibbsindex++)
    {
      for(i = 0; i < numind; i++)
	{
	  sumfinb[i] += tempfinb[gibbsindex][i];
	}
    }
    for(i = 0; i < numind; i++)
      {
	voglfinblist[i] = sumfinb[i]/numgibbs;
      }
    /* next, each allele has a frequency that is the average of its frequencies over all samples */
    for(j = 0; j < nmark; j++)
      {
	allelesum = (double*)malloc(numalleles[j]*sizeof(double));
	for(k = 0; k < numalleles[j]; k++)
	  {
	    allelesum[k] = 0;
	    for(gibbsindex = 0; gibbsindex < numgibbs; gibbsindex++)
	      {
		allelesum[k]+= tempfreqs[gibbsindex][j][k];
	      }
	    allelefreqs[j][k] = allelesum[k]/numgibbs;
	  }
	if(allelefreqs[j][numalleles[j]-1] < 0)
	  {
	    printf("\n(3097) marker %d\n", j);exit(0);
	  }
	    free(allelesum);
      }
    //printf("\n(2036) vogl %4.4f, myalpha = %4.4f", voglfinblist[0], myalpha);fflush(stdout);
  /* freeing stuff*/
    free(sumfinb);
  for(gibbsindex = 0; gibbsindex < numgibbs; gibbsindex++)
    {
      for(j = 0; j < nmark; j++)
	{
	  free(tempfreqs[gibbsindex][j]);
	}
      free(tempfreqs[gibbsindex]);
    }
  free(tempfreqs);

  for(gibbsindex = 0; gibbsindex < numgibbs; gibbsindex++)
    {
      free(tempfinb[gibbsindex]);
    }
  free(tempfinb);

  for(i = 0; i < numind; i++)
    {
      free(xij[i]);
    }
  free(xij);

  return;
}
/*******************************************/
/* Updates xij for the Vogl algorithm.  ****/
/*******************************************/
void voglupdatefinb(int*** genos, int*** tempgenos, double*** tempfreqs, double** tempfinb, int gibbsindex, int** xij, double myalpha, double mybeta, int* numalleles, int numind, int nmark)
{
  int i, j, k;
  double F;
  double sumx, sumy;
  double x;

  for(i = 0; i < numind; i++)
    {
      sumx = 0;/* counts number of IBD genotypes*/
      sumy = 0;/* counts number of nonIBD genotypes*/
      for(j = 0; j < nmark; j++)
	{
	  //printf("\n(580) i = %d, j = %d", i, j);fflush(stdout);
	  if(genos[i][j][0] == 99)/* missing genotypes are assumed to be homozygous for the null allele*/
	    {
	      /* automatically has genotype Null,Null*/
	      x = drand48();/* decide if they are IBD alleles*/
	      if(x < (1 - tempfinb[gibbsindex][i])*tempfreqs[gibbsindex][j][numalleles[j]-1]/((1 - tempfinb[gibbsindex][i])*tempfreqs[gibbsindex][j][numalleles[j]-1] + tempfinb[gibbsindex][i]))/* formula (1-f)p_n/((1-f)p_n + f) verified*/
		{/* not IBD*/
		  xij[i][j] = 0;
		  sumy++;
		}
	      else
		{/* IBD*/
		  xij[i][j] = 1;
		  sumx++;
		}
	      /* IBD or not, genotype for missing is AnAn*/
	      tempgenos[i][j][0] = numalleles[j]-1;
	      tempgenos[i][j][1] = numalleles[j]-1;
	    }
	  else if(genos[i][j][0] == genos[i][j][1])/* observed homozygous, not missing*/
	    {
	      F =tempfinb[gibbsindex][i];
	      x = drand48();
	      if(x < (F + (1.0 - F)*tempfreqs[gibbsindex][j][genos[i][j][0]])/(F + (1-F)*(tempfreqs[gibbsindex][j][genos[i][j][0]] + 2*tempfreqs[gibbsindex][j][numalleles[j]-1])))/*With this prob we are really homozygous for these observed alleles*/
		{
		  x = drand48();
		  if(x < (1 - tempfinb[gibbsindex][i])*tempfreqs[gibbsindex][j][genos[i][j][0]]/((1 - tempfinb[gibbsindex][i])*tempfreqs[gibbsindex][j][genos[i][j][0]] + tempfinb[gibbsindex][i]))/* given that we are homozygous for these alleles, this makes the alleles not IBD with the correct probability*/
		    {
		      xij[i][j] = 0;
		      sumy++;
		    }
		  else/* not ibd*/
		    {
		      xij[i][j] = 1;
		      sumx++;
		    }

		  /* setting tempgenos*/
		    /* homozygous for the observed alleles*/
		  tempgenos[i][j][0] = genos[i][j][0];
		  tempgenos[i][j][1] = genos[i][j][1];

		}
	      else /* we really are heterozygous with the null allele*/
		{
		  xij[i][j] = 0;/* can't be IBD if we are heterozygous*/
		  sumy++;
		  /* setting tempgenos*/
		  /* the genotype has been observed to be homozygous, so we want to make it heterozygous with the null allele.*/

		      tempgenos[i][j][0] = genos[i][j][0];
		      tempgenos[i][j][1] = numalleles[j]-1;
		}
	    }
	  else
	    {/* we are visibly heterozygous, and so can't be IBD*/
	      xij[i][j] = 0;
	      sumy++;
	      tempgenos[i][j][0] = genos[i][j][0];
	      tempgenos[i][j][1] = genos[i][j][1];
	    }
	}
      /* inherent assumption here that our prior for finb was beta with alpha = myalpha, beta = mybeta*/
      tempfinb[gibbsindex+1][i] = rbeta(sumx+myalpha, sumy+mybeta);
    }
  return;
}
