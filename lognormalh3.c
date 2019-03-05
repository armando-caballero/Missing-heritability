
/* lognormalh3.c (data from lognormal to check heritability)*/

/* *********************************************** */

#include "libhdr"
#define NN 661  /* max number of loci */
#define MM 1000001  /* max number of NIND */
#define MMM 23  /* max NCRO */

int NIND, NCRO, NLOCI, REP;
int RM[NN], i, j, k, l, r;
int gm[MM][MMM][2], sm[MM][MMM][2];
int gmDZ[MM][MMM][2], EPIS;
int num_homo_whole, homo[MM][31];

double genvalm_a[MM], pm_a[MM];
double aLN[1000001], qLN[1000001], a[MMM][31], ha[MMM][31], q[MMM][31];
double d_s, d_a, alfa_a, alfa_s, VA, VD;
double h_a, VE;
double gmean_a, gvar_a;
double genvalm_a_DZ[MM], pm_a_DZ[MM], genvalm_a_MZ[MM], pm_a_MZ[MM];
double tMZ, tDZ, slogN[1000001], w, h2op;
double VAgwas, VPgwas, alfa_gwas, var_pm;

struct acc AVE_VG, AVE_VA, AVE_VD, AVE_H2, AVE_h2, AVE_d2, AVE_tMZ, AVE_tDZ, AVE_h2op, AVE_VAgwas, AVE_VPgwas;


FILE *fptr, *flogN;


/* ***************************************************************** */
     

main()
{
	fptr = fopen ("dfilename","w");

	flogN = fopen ("sq.txt","r");
	for (i=1; i<=100000; i++)
	{
		fscanf(flogN,"%lf", &w);
		aLN[i] = 2.0*w;
		fscanf(flogN,"%lf", &w);
		qLN[i] = w;
		if (qLN[i] > 1.0)	qLN[i]=1.0;

		if (tracelevel!=0)	fprintf(fptr,"aLN=%f  qLN=%f\n", aLN[i], qLN[i]);
	}
	fclose(flogN);

	getinputs();
	recombination_masks();

	for (r=1; r<=REP; r++)
	{
		initialize();
		selected_genes();
		genotypic_values();
		if (EPIS==1) genotypic_values_epistasis_whole();
		phenotypeB();
		mating_MZDZ();
		ANOVA_DZ();
		ANOVA_MZ();
		OP_regression();
		GWAS();
		averages();
	}

	printout();
	writeseed();
}


/* ***************************************************************** */


getinputs()
{
	tracestart();
	getseed();
	getintandskip("Number of individuals (max 2000000) :",&NIND,2,1000000);
	getintandskip("Number of pieces of chromosome (min 1, max 50) :",&NCRO,1,50);
	getintandskip("Number of loci per piece (first is neutral) (min 2, max 30) :",&NLOCI,2,30);
	getintandskip("EPIS (1: yes; 0:no) :",&EPIS,0,1);
	getrealandskip("h_a :", &h_a,0.0,(double)infinity);
	getrealandskip("VE :", &VE,0.0,(double)infinity);
	getintandskip("REP :",&REP,1,10000);
}


/* **************************************************************** */


void recombination_masks()
{
	for (l=0; l<NLOCI; l++)   RM[l]=pow(2.0,(double)l);
}


/* ***************************************************************** */


void initialize ()
{
    int rnd;

    for (k=0; k<NCRO; k++)
    for (l=0; l<NLOCI; l++)
    {
		rnd = (int)(uniform()*100000)+1;
		a[k][l] = aLN[rnd];
		q[k][l] = qLN[rnd];
		ha[k][l] = h_a;
		if (tracelevel!=0)	fprintf(fptr,"a=%f  h=%f\n", a[k][l], ha[k][l]);
    }

    for (k=0; k<NCRO; k++)
    for (l=0; l<NLOCI; l++)   
    for (i=0; i<NIND; i++)
    {
	    	gm[i][k][0]=0;
	    	gm[i][k][1]=0;
    }

    for (k=0; k<NCRO; k++)
    for (l=0; l<NLOCI; l++)   
    {
	    for (i=0; i<NIND; i++)
	    {
	    	if (uniform() < q[k][l])	gm[i][k][0]=(gm[i][k][0] | RM[l]);
	    	if (uniform() < q[k][l])	gm[i][k][1]=(gm[i][k][1] | RM[l]);
	    }
    }

    if (tracelevel!=0)	for (i=0; i<NIND; i++)   fprintf(fptr,"%d  %d\n",gm[i][0][0],gm[i][0][1]);
}


/* **************************************************************** */


void selected_genes ()
{
    double AA, Aa, aa;

    VA=0.0; VD=0.0;

    if (tracelevel!=0)	fprintf(fptr,"\n selected genes\npiece\tlocus\tAA\tAa\taa\t q\t d_a\t alfa_a\t VA\n");

    for (k=0; k<NCRO; k++)
    {
	for (l=0; l<NLOCI; l++)
	{
	    AA=0.0; Aa=0.0; aa=0.0;

	    for (i=0; i<NIND; i++)
	    {
	         if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))           aa+=1.0;
	         else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))   AA+=1.0;
	         else							       Aa+=1.0;
	    }

	    q[k][l] = (aa/(double)NIND)+(Aa/(2.0*(double)NIND));

		d_a = (a[k][l]/2.0) * (2.0*ha[k][l] - 1.0);
		if (a[k][l] >= 0.0)		alfa_a = (a[k][l]/2.0) + ( d_a * (1.0 - 2.0*q[k][l]) );
		else			alfa_a = (-a[k][l]/2.0) + ( d_a * (2.0*q[k][l] - 1.0) );
		
		VA += 2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]);
		VD +=  pow(2.0 * d_a * q[k][l] * (1.0 - q[k][l]), 2.0);

	    if (tracelevel!=0)    fprintf(fptr,"%d\t%d\t%1.0f\t%1.0f\t%1.0f\t%f\t%f\t%f\t%f\n",k,l,AA,Aa,aa,q[k][l],d_a,alfa_a,2.0 * alfa_a * alfa_a * q[k][l] * (1.0 - q[k][l]));
	}
    }
}


/* ***************************************************************** */


void genotypic_values()
{
	for (i=0; i<NIND; i++)
	{
		genvalm_a[i]=0.0;

		for (k=0; k<NCRO; k++)
		{
		   for (l=0; l<NLOCI; l++)
		   {
	    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				genvalm_a[i] += a[k][l];
			}
			else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l])) 								/* AA */;
			else
			{
				genvalm_a[i] += a[k][l]*ha[k][l];
			}
		   }
		}
		if (tracelevel!=0)	fprintf(fptr, " %d    genvalm_a = %f\n", i, genvalm_a[i]);
	}
}


/* **************************************************************** */


genotypic_values_epistasis_whole()
{
	for (i=0; i<NIND; i++)
	{
		for (k=0; k<NCRO; k++)
		for (l=0; l<NLOCI; l++)	homo[k][l] = 0;

		num_homo_whole = 0;

		for (k=0; k<NCRO; k++)
		{
			for (l=0; l<NLOCI; l++)
			{
		    		if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l])) 
				{
					homo[k][l] = 1;
					num_homo_whole ++;
				}
			}
		}

		if (num_homo_whole >= 2)
		{
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)
			{
				if (homo[k][l] == 1)
				{
					genvalm_a[i] += a[k][l];
				}
			}
		}
	}
}


/* ***************************************************** */


void phenotypeB ()
{
	double gsum_a=0.0, gsum2_a=0.0;

	for (i=0; i<NIND; i++)
	{	
	    pm_a[i] = genvalm_a[i] + normal(0.0, sqrt(VE));
	    gsum_a += genvalm_a[i];
	    gsum2_a += (genvalm_a[i]*genvalm_a[i]);
	}
	gmean_a = gsum_a/(double)NIND;
	gvar_a = (gsum2_a - (gsum_a*gsum_a / (double)NIND)) / (double)NIND;

	if (tracelevel!=0)	fprintf(fptr,"\ngmean_a = %f  gvar_a = %f\n", gsum_a/(double)NIND, (gsum2_a - (gsum_a*gsum_a / (double)NIND)) / ((double)NIND-1.0));
}


/* **************************************************************** */


void mating_MZDZ()
{
	int EE[MMM], FF[MMM], p1, p2;
		
	for (i=0; i<NIND; i++)
	{
		for (k=0; k<NCRO; k++)
		{
			sm[i][k][0]=gm[i][k][0];
			sm[i][k][1]=gm[i][k][1];
		}
	}

	if (tracelevel!=0)	fprintf(fptr, "\n Parents MZ and DZ\n");

	for (i=0; i<NIND; i++)
	{
		if (i%2 == 0)
		{
			p1 = i;
			p2 = p1 + 1;
		}
		else
		{
			p1 = i - 1;
			p2 = p1 + 1;
		}

		if (tracelevel!=0)	fprintf(fptr, "i=%d p1=%d\tp2=%d\n", i, p1, p2);

		for (k=0; k<NCRO; k++)
		{
	   		EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
	   		FF[k] = ~EE[k];
	   		gmDZ[i][k][0]=((EE[k]&sm[p1][k][0])|(FF[k]&sm[p1][k][1]));
		}
		for (k=0; k<NCRO; k++)
		{
	   		EE[k] = (int)(uniform()*(pow(2.0,(double)NLOCI)));
	   		FF[k] = ~EE[k];
	   		gmDZ[i][k][1]=((EE[k]&sm[p2][k][0])|(FF[k]&sm[p2][k][1]));
		}


		genvalm_a_DZ[i]=0.0;

		for (k=0; k<NCRO; k++)
		{
		   for (l=0; l<NLOCI; l++)
		   {
	    		if (((gmDZ[i][k][0] & RM[l])==RM[l])&&((gmDZ[i][k][1] & RM[l])==RM[l]))  			
			{
				genvalm_a_DZ[i] += a[k][l];
			}
			else if (((gmDZ[i][k][0] & RM[l])!=RM[l])&&((gmDZ[i][k][1] & RM[l])!=RM[l])) /* AA */;
			else
			{
				genvalm_a_DZ[i] += a[k][l]*ha[k][l];
			}
		   }
		}
		if (tracelevel!=0)	fprintf(fptr, " %d    genvalm_a_DZ = %f\n", i, genvalm_a_DZ[i]);

		pm_a_DZ[i] = genvalm_a_DZ[i] + normal(0.0, sqrt(VE));

		if (EPIS==1) 
		{
			for (k=0; k<NCRO; k++)
			for (l=0; l<NLOCI; l++)	homo[k][l] = 0;

			num_homo_whole = 0;

			for (k=0; k<NCRO; k++)
			{
				for (l=0; l<NLOCI; l++)
				{
				    	if (((gmDZ[i][k][0] & RM[l])==RM[l])&&((gmDZ[i][k][1] & RM[l])==RM[l])) 
					{
						homo[k][l] = 1;
						num_homo_whole ++;
					}
				}
			}

			if (num_homo_whole >= 2)
			{
				for (k=0; k<NCRO; k++)
				for (l=0; l<NLOCI; l++)
				{
					if (homo[k][l] == 1)
					{
						genvalm_a_DZ[i] += a[k][l];
					}
				}
			}
			pm_a_DZ[i] = genvalm_a_DZ[i] + normal(0.0, sqrt(VE));

			if (tracelevel!=0)	fprintf(fptr, "\nEpistatic\n");
			if (tracelevel!=0)	fprintf(fptr, "i=%d genval_a_DZ=%f pm_a_DZ=%f\n", i, genvalm_a_DZ[i], pm_a_DZ[i]);
		}
	}
}


/* **************************************************************** */


ANOVA_DZ()
{
	int n = NIND/2;

	double Xip2=0.0, Xpp=0.0, Xpp2=0.0, Xij2=0.0;
	double SCG=0.0, SCE=0.0, SCT=0.0, MCG=0.0, MCE=0.0, Vb=0.0, Vw=0.0;

	tDZ=0.0;

	if (tracelevel!=0)	fprintf(fptr, "\nAnova DZ\n");

	for (i=0; i<(2*n); i++)
	{
		if (tracelevel!=0)	fprintf(fptr, "i=%d  genvalm_a_DZ=%f  pm_a_DZ=%f\n", i, genvalm_a_DZ[i], pm_a_DZ[i]);

		if (i%2 == 0)   Xip2 += pow((pm_a_DZ[i] + pm_a_DZ[i+1]), 2.0);
		Xpp += pm_a_DZ[i];
		Xij2 += (pm_a_DZ[i]*pm_a_DZ[i]);
	}
	Xpp2 = Xpp * Xpp;

	if (tracelevel!=0)	fprintf(fptr, "\n Xij2=%f  Xpp2=%f  Xip2=%f\n", Xij2, Xpp2, Xip2);

	SCG = (Xip2/2.0) - (Xpp2/(2.0*n));
	SCE = Xij2 - (Xip2/2.0);
	SCT = Xij2 - (Xpp2/(2.0*n));

	MCG = SCG / (n - 1.0);
	MCE = SCE / n;

	if (tracelevel!=0)	fprintf(fptr, "\n SCG=%f  SCE=%f  SCT=%f  MCG=%f  MCE=%f\n", SCG, SCE, SCT, MCG, MCE);

	Vw = MCE;
	Vb = (MCG - MCE) / 2.0;

	tDZ = Vb / (Vb + Vw);

	if (tracelevel!=0)	fprintf(fptr, "\n Vb=%f  Vw=%f  tDZ=%f  h2+d2/2=%f\n", Vb, Vw, tDZ, 2.0*tDZ);
}


/* ***************************************************************** */


ANOVA_MZ()
{
	int n = NIND/2;

	double Xip2=0.0, Xpp=0.0, Xpp2=0.0, Xij2=0.0;
	double SCG=0.0, SCE=0.0, SCT=0.0, MCG=0.0, MCE=0.0, Vb=0.0, Vw=0.0;

	tMZ=0.0;

	if (tracelevel!=0)	fprintf(fptr, "\nAnova MZ\n");

	for (i=0; i<NIND; i++)
	{
	    genvalm_a_MZ[i] = genvalm_a_DZ[i];
	    pm_a_MZ[i] = pm_a_DZ[i];
	}

	for (i=1; i<NIND; i+=2)
	{
	    genvalm_a_MZ[i] = genvalm_a_MZ[i-1];
	    pm_a_MZ[i] = genvalm_a_MZ[i] + normal(0.0, sqrt(VE));
	}

	for (i=0; i<(2*n); i++)
	{
		if (tracelevel!=0)	fprintf(fptr, "i=%d  genvalm_a_MZ=%f  pm_a_MZ=%f\n", i, genvalm_a_MZ[i], pm_a_MZ[i]);

		if (i%2 == 0)   Xip2 += pow((pm_a_MZ[i] + pm_a_MZ[i+1]), 2.0);
		Xpp += pm_a_MZ[i];
		Xij2 += (pm_a_MZ[i]*pm_a_MZ[i]);
	}
	Xpp2 = Xpp * Xpp;

	if (tracelevel!=0)	fprintf(fptr, "\n Xij2=%f  Xpp2=%f  Xip2=%f\n", Xij2, Xpp2, Xip2);

	SCG = (Xip2/2.0) - (Xpp2/(2.0*n));
	SCE = Xij2 - (Xip2/2.0);
	SCT = Xij2 - (Xpp2/(2.0*n));

	MCG = SCG / (n - 1.0);
	MCE = SCE / n;

	if (tracelevel!=0)	fprintf(fptr, "\n SCG=%f  SCE=%f  SCT=%f  MCG=%f  MCE=%f\n", SCG, SCE, SCT, MCG, MCE);

	Vw = MCE;
	Vb = (MCG - MCE) / 2.0;

	tMZ = Vb / (Vb + Vw);

	if (tracelevel!=0)	fprintf(fptr, "\n Vb=%f  Vw=%f  tMZ=%f  h2+d2=%f\n", Vb, Vw, tMZ, tMZ);

	if (tracelevel!=0)	fprintf(fptr, "\n 2(tMZ-tDZ)=h2+(3/2)d2=%f\n", 2.0*(tMZ-tDZ));
}


/* ***************************************************************** */


OP_regression()
{
	int n = NIND/2;
	
	double sumP=0.0, sumO=0.0, sumP2=0.0, sumOP=0.0;
	double varP, covOP;

	if (tracelevel!=0)	fprintf(fptr, "\nOP regression\n");

	for (i=0; i<(2*n); i++)	
	{
		if (tracelevel!=0)	fprintf(fptr, "i=%d %f %f\n", i, pm_a[i], pm_a_DZ[i]);
	}

	for (i=0; i<(2*n); i++)
	if (i%2 == 0)
	{
		sumP += (pm_a[i] + pm_a[i+1]) / 2.0;
		sumO += (pm_a_DZ[i] + pm_a_DZ[i+1]) / 2.0;
		sumP2 += ((pm_a[i] + pm_a[i+1]) / 2.0) * ((pm_a[i] + pm_a[i+1]) / 2.0);
		sumOP += ((pm_a[i] + pm_a[i+1]) / 2.0) * ((pm_a_DZ[i] + pm_a_DZ[i+1]) / 2.0);
	}
	varP = ( sumP2 - (sumP*sumP/n) ) / (n-1.0);
	covOP = ( sumOP - (sumP*sumO/n) ) / (n-1.0);
	h2op = covOP / varP;

	if (tracelevel!=0)	fprintf(fptr, "\n sumP=%f  sumO=%f  sumOP=%f\n", sumP, sumO, sumOP);
	if (tracelevel!=0)	fprintf(fptr, "\n varP=%f  covOP=%f  h2op=%f\n", varP, covOP, h2op);

}


/* ***************************************************************** */


GWAS()
{
	double sum_cop, sum_pm, sum_gm, sum_pm2, sum_cop2, sum_pm_cop, sum_gm_cop, copies[MM];
	double var_cop, cov_pm_cop, cov_gm_cop, alfa_gwas;

	VAgwas = 0.0;

	for (k=0; k<NCRO; k++)
	for (l=0; l<NLOCI; l++)
	{
		sum_cop = 0.0;
		sum_pm = 0.0;
		sum_gm = 0.0;
		sum_pm2 = 0.0;
		sum_cop2 = 0.0;
		sum_pm_cop = 0.0;
		sum_gm_cop = 0.0;

		for (i=0; i<NIND; i++)
		{
			if (((gm[i][k][0] & RM[l])==RM[l])&&((gm[i][k][1] & RM[l])==RM[l]))
			{
				copies[i]=2.0;
			}
			else    if (((gm[i][k][0] & RM[l])!=RM[l])&&((gm[i][k][1] & RM[l])!=RM[l]))
			{
				copies[i]=0.0;
			}
			else
			{
				copies[i]=1.0;
			}

			sum_cop += copies[i];
			sum_pm += pm_a[i];
			sum_pm2 += pm_a[i]*pm_a[i];
//			sum_gm += genvalm_a[i];
			sum_cop2 += copies[i] * copies[i];
			sum_pm_cop += copies[i] * pm_a[i];
//			sum_gm_cop += copies[i] * genvalm_a[i];

			if ((tracelevel!=0) && (k==0))	fprintf(fptr,"i=%d  k=%d  l=%d  genvalm_a=%f  pm_a=%f  copies=%f\n", i, k, l, genvalm_a[i], pm_a[i], copies[i]);
		}

		var_pm = ( sum_pm2 - (sum_pm*sum_pm/(double)NIND) ) / ((double)NIND-1.0);

		//	CALCULATION of VAgwas

		var_cop = ( sum_cop2 - (sum_cop*sum_cop/(double)NIND) ) / ((double)NIND-1.0);

		cov_pm_cop = ( sum_pm_cop - (sum_pm*sum_cop/(double)NIND) ) / ((double)NIND-1.0);
		if (var_cop != 0.0)	alfa_gwas = cov_pm_cop / var_cop;

		if (tracelevel!=0)	fprintf(fptr,"k=%d  l=%d  cov_pm_cop=%f var_cop=%f\n", k, l, cov_pm_cop, var_cop);
		if (tracelevel!=0)	fprintf(fptr,"k=%d  l=%d  q=%f  a=%f  alfa_gwas=%f\n", k, l, q[k][l], a[k][l], alfa_gwas);

		VAgwas += (2.0 * alfa_gwas * alfa_gwas * q[k][l] * (1.0 - q[k][l]));

//		if ((k==0)&&(l==0))	printf("a=%f  alfa_gwas=%f\n", q[k][l], alfa_gwas);
	}
}


/* ***************************************************************** */


averages()
{
	accum (&AVE_VG, gvar_a);
	accum (&AVE_VA, VA);
	accum (&AVE_VD, VD);
	accum (&AVE_H2, gvar_a / (gvar_a + VE));
	accum (&AVE_h2, VA / (gvar_a+VE));
	accum (&AVE_d2, VD / (gvar_a+VE));
	accum (&AVE_tMZ, tMZ);
	accum (&AVE_tDZ, tDZ);
	accum (&AVE_h2op, h2op);
	accum (&AVE_VAgwas, VAgwas);
	accum (&AVE_VPgwas, var_pm);
}


/* ***************************************************************** */


printout()
{
	printf("\n\nLOCI=%d    h=%f    EPIS=%d\n", NCRO*NLOCI, h_a, EPIS); 

	printf("VG=%6.4f+-%6.4f    VA=%6.4f+-%6.4f    VD=%6.4f+-%6.4f\n",
		accmean(&AVE_VG), se(&AVE_VG), accmean(&AVE_VA), se(&AVE_VA), accmean(&AVE_VD), se(&AVE_VD));
	printf("t_MZ=%6.4f+-%6.4f    t_DZ=%6.4f+-%6.4f     tMZ-2tDZ=%6.4f\n",
		accmean(&AVE_tMZ), se(&AVE_tMZ), accmean(&AVE_tDZ), se(&AVE_tDZ), accmean(&AVE_tMZ)-2.0*accmean(&AVE_tDZ)); 
	printf("VAgwas=%6.4f+-%6.4f    VPgwas=%6.4f+-%6.4f\n",
		accmean(&AVE_VAgwas), se(&AVE_VAgwas), accmean(&AVE_VPgwas), se(&AVE_VPgwas)); 
	printf("H2=%6.4f+-%6.4f    h2=%6.4f+-%6.4f    d2=%6.4f+-%6.4f\n",
		accmean(&AVE_H2), se(&AVE_H2), accmean(&AVE_h2), se(&AVE_h2), accmean(&AVE_d2), se(&AVE_d2));
	printf("h2tMZtDZ=%6.4f   h2OP=%6.4f+-%6.4f    h2gwas=%6.4f\n\n",
		2.0*(accmean(&AVE_tMZ)-accmean(&AVE_tDZ)), accmean(&AVE_h2op), se(&AVE_h2op), accmean(&AVE_VAgwas)/(gvar_a+VE)); 
}


/* ************************************************************ */
