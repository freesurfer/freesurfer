
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "cma.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "gca.h"
#include "flash.h"

static char vcid[] = "$Id: mri_gca_ambiguous.c,v 1.1 2003/10/08 21:23:37 fischl Exp $";

static double scale = 5*5e3 ;

int main(int argc, char *argv[]) ;

int GCAscale(GCA *gca_flash, double min_val, double max_val) ;
static  MRI *GCAcomputeAmbiguity(GCA *gca, MRI *mri, double *pamb) ;
static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int optimize = 0 ;

static double TR = 20;
static double TE  = 3  ;
static double MIN_FA1 = 1  ;
static  double MAX_FA1 = 40 ;
static double MIN_FA2 = 1  ;
static  double MAX_FA2 = 40 ;
static  double FA_STEP = 1 ;
static int append = 0 ;
static char *fname = "amb.log" ; 

int
main(int argc, char *argv[])
{
  char               **av, *out_name, *gca_name ;
  int                ac, nargs ;
	GCA                *gca  ;
	MRI                *mri = NULL ;
	double             amb ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_gca_ambiguous.c,v 1.1 2003/10/08 21:23:37 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    print_help() ;

  gca_name = argv[1] ;
  out_name = argv[2] ;

	printf("reading gca  from %s....\n", gca_name) ;
  gca = GCAread(gca_name) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read gca file %s", Progname, gca_name) ;


	if (optimize)
	{
		double fa1,  fa2, min_fa1,  min_fa2  ;
		GCA    *gca_flash ;
		double TEs[2],  TRs[2],  FAs[2], min_amb, scale ;
		FILE   *fp ;

		fp = fopen(fname, append ? "a" : "w") ;

		TRs[0]  = TR ; TRs[1] = TR ;
		TEs[0]  = TE ; TEs[1] = TE ;
		min_amb = 10000000  ; min_fa1 = min_fa2 = -1 ; 
		scale = 1.0/(gca->prior_width * gca->prior_height * gca->prior_depth);
		for (fa1 = MIN_FA1  ; fa1 <= MAX_FA1 ; fa1 += FA_STEP)
		{
			FAs[0]  = RADIANS(fa1)  ;
			for (fa2 = MIN_FA2 ; fa2 <= MAX_FA2 ; fa2 += FA_STEP)
			{
				printf("testing flip angles %2.3f, %2.3f\n", fa1, fa2) ;
				FAs[1] = RADIANS(fa2) ;
				gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 2);
				/*				GCAscale(gca_flash, 0, 75) ;*/
				mri =  GCAcomputeAmbiguity(gca_flash, mri, &amb) ;
				amb *= scale ;
				printf("\tflip angles %2.3f, %f: ambiguity == %f\n", fa1, fa2, amb) ;
				if (amb  < min_amb || min_fa1  < 0)
				{
					min_amb = amb ;
					min_fa1 = fa1 ;
					min_fa2 = fa2 ;
					printf("minimimum ambiguity at fas %2.3f, %2.3f (%f)\n", min_fa1,  min_fa2, min_amb) ;
				}
				fprintf(fp, "%f %f %f\n", fa1,fa2, amb) ;
				fflush(fp) ;
				GCAfree(&gca_flash)  ;
			}
		}

		MRIfree(&mri) ;
		fclose(fp) ;
		printf("minimimum ambiguity at  fas  %2.3f,  %2.3f (%2.5f)\n", min_fa1,  min_fa2, min_amb) ;
		FAs[0] = min_fa1 ; FAs[1] = min_fa2 ;
		gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 2);
		mri =  GCAcomputeAmbiguity(gca_flash, NULL, &amb) ;
		MRIwrite(mri, out_name) ;
		MRIfree(&mri) ;
		mri = MRIallocSequence(gca_flash->prior_width, gca_flash->prior_height,
									 gca_flash->prior_depth, MRI_FLOAT, gca_flash->ninputs) ;
		GCAbuildMostLikelyVolume(gca_flash,mri) ;
		MRIwrite(mri, "mri_gca.mgh")  ;
		GCAfree(&gca) ;
		GCAfree(&gca_flash) ;
		MRIfree(&mri) ;
	}
	else
	{
		mri =  GCAcomputeAmbiguity(gca, NULL, &amb) ;
		GCAbuildMostLikelyVolume(gca,mri) ;
		MRIwrite(mri, "mri_gca.mgh")  ;
	}


  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "debug_voxel"))
	{
		Gx = atoi(argv[2]) ;
		Gy = atoi(argv[3]) ;
		Gz = atoi(argv[4]) ;
		printf("debugging  voxel (%d, %d,  %d)\n", Gx, Gy,  Gz) ;
		nargs = 3 ;
	}
  else if (!stricmp(option, "FA1"))
	{
		MIN_FA1 = atof(argv[2]) ;
		MAX_FA1 = atof(argv[3]) ;
		printf("covering range for 1st flip angle %2.3f-->%2.3f\n", MIN_FA1, MAX_FA1) ;
		nargs = 2 ;
	}
  else if (!stricmp(option, "FA2"))
	{
		MIN_FA2 = atof(argv[2]) ;
		MAX_FA2 = atof(argv[3]) ;
		printf("covering range for 2nd flip angle %2.3f-->%2.3f\n", MIN_FA2, MAX_FA2) ;
		nargs = 2 ;
	}
  else switch (toupper(*option))
  {
	case 'A':
		append = 1 ;
		nargs = 1 ;
		fname = argv[2] ;
		printf("appending output to %s...\n", fname) ;
		break ;
	case 'S':
		FA_STEP = atof(argv[2]) ;
		printf("using step size %2.3f\n", FA_STEP) ;
		nargs = 1 ;
		break ;
	case 'O':
		optimize = atoi(argv[2])  ;
		nargs = 1  ;
		break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  default:
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    exit(1) ;
    break ;
  }

  return(nargs) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <gca file> <output volume>\n",
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
					"\nThis program computes an ambiguity measure across a  GCA and output an MR image of it\n") ;
  exit(1) ;
}


static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;

}
static  MRI *
GCAcomputeAmbiguity(GCA *gca, MRI *mri, double *pamb)
{
	int        xp, yp, zp, l1, l2, i1,  i2,  i, label_count[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1], label1_count[MAX_CMA_LABEL+1] ;
	GCA_PRIOR  *gcap  ;
	double     ambiguity, atotal, p1,  p2, amax, total_ambiguity, min_I, max_I, std,  Istep,  I ;
	float      vals[MAX_GCA_INPUTS];
	double     label_ambiguity[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1];
	GC1D       *gc1, *gc2 ;

	memset(label_ambiguity, 0, sizeof(label_ambiguity)) ;
	memset(label_count, 0, sizeof(label_count)) ;
	memset(label1_count, 0, sizeof(label1_count)) ;
	if (!mri)
		mri =  MRIalloc(gca->prior_width,gca->prior_height, gca->prior_depth, MRI_FLOAT)  ;
	mri->xsize = mri->ysize = mri->zsize = gca->prior_spacing ;

	for (total_ambiguity = xp = 0  ; xp  < gca->prior_width ; xp++)
	{
		for  (yp  = 0 ;  yp < gca->prior_height ; yp++)
		{
			for (zp = 0 ;  zp  < gca->prior_depth  ; zp++)
			{
				gcap =  &gca->priors[xp][yp][zp] ;
				if (xp == Gx && yp == Gy && zp == Gz)
					DiagBreak() ;
				if (gcap->nlabels <= 1)
					continue ;

				min_I  = 10000000 ; max_I = 0 ;
				for (i1 = 0 ; i1  <  gcap->nlabels ; i1++)
				{
					l1 = gcap->labels[i1] ;
					gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
					if (!gc1)
						continue ;
					std = sqrt(covariance_determinant(gc1, gca->ninputs)) ;
					for (i = 0  ; i < gca->ninputs ; i++)
					{
#define NSTDS 1
						if (gc1->means[i]+NSTDS*std > max_I)
							max_I = gc1->means[i]+NSTDS*std ;
						if (gc1->means[i]-NSTDS*std < min_I)
							min_I = gc1->means[i]-NSTDS*std ;
					}
				}

				for (amax =  atotal = 0.0, i1 = 0 ;  i1 < gcap->nlabels ; i1++)
				{
					for (i1 = 0 ;  i1 < gcap->nlabels ; i1++)
					{
						l1 = gcap->labels[i1] ;
						gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
						if (!gc1)
							continue ;
						label1_count[l1]++ ;
						for (i2 = 0 ;  i2 < i1 ; i2++)
						{
							l2 = gcap->labels[i2] ;
							gc2 = GCAfindPriorGC(gca, xp, yp, zp,  l2) ;
							if (!gc2)
								continue ;

							if (l1 == 0 || l2 == 0)
								continue ;
							if  ((l1 == 17 &&  l2  == 7) ||(l2 == 17 &&  l1  == 7))
 								DiagBreak() ;

#define ISTEPS 100

							/* marginalize over intensity */
							Istep = (max_I - min_I) / ISTEPS ;
							ambiguity = 0.0 ;
							if (gca->ninputs == 1)
							{
								for (ambiguity = 0,  I = min_I;  I <= max_I  ; I += Istep)
								{
									vals[0] =  I ;
									p1 = scale*GCAcomputeConditionalDensity(gc1, vals, gca->ninputs, l1) ;
									p2 = scale*GCAcomputeConditionalDensity(gc2, vals, gca->ninputs, l2) ;
									ambiguity += (p1*p2) ;
								}
							}
							else if (gca->ninputs == 2)
							{
								double I1, I2, Istep1, Istep2, std1, std2, scale1, scale2, dsq ;
								MATRIX *m1_cov_inv, *m2_cov_inv ;
								VECTOR *v1, *v2, *vtmp ;

								m1_cov_inv = load_inverse_covariance_matrix(gc1, NULL, gca->ninputs) ;
								m2_cov_inv = load_inverse_covariance_matrix(gc2, NULL, gca->ninputs) ;
								v1 = VectorAlloc(2, MATRIX_REAL) ;
								v2 = VectorAlloc(2, MATRIX_REAL) ;
								vtmp = VectorAlloc(2, MATRIX_REAL) ;
								scale1 = scale*(1.0 / (pow(2*M_PI,gca->ninputs/2.0)*sqrt(covariance_determinant(gc1, gca->ninputs))));
								scale2 = scale*(1.0 / (pow(2*M_PI,gca->ninputs/2.0)*sqrt(covariance_determinant(gc2, gca->ninputs))));

								if (l1 == Left_Hippocampus && l2 == Left_Amygdala)
									DiagBreak() ;
								std1 = sqrt(gc1->covars[0]) ; std2 = sqrt(gc1->covars[2]) ; 
#undef ISTEPS
#define ISTEPS 20
								Istep1 = std1/(ISTEPS/4) ; Istep2 = std2/(ISTEPS/4) ;
								for (ambiguity = 0,  I1 = gc1->means[0]-(ISTEPS/2)*Istep1;  I1 <= gc1->means[0]+(ISTEPS/2)*Istep1  ; I1 += Istep1)
								{
									VECTOR_ELT(v1, 1) = I1-gc1->means[0] ; VECTOR_ELT(v2, 1) = I1-gc2->means[0] ;
									for (I2 = gc1->means[1]-(ISTEPS/2)*Istep2;  I2 <= gc1->means[1]+(ISTEPS/2)*Istep2  ; I2 += Istep2)
									{
										VECTOR_ELT(v1, 2) = I2-gc1->means[1] ;
										MatrixMultiply(m1_cov_inv, v1, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
										dsq = VectorDot(v1, vtmp) ;
										p1 = scale1*exp(-0.5*dsq) ;

										VECTOR_ELT(v2, 2) = I2-gc2->means[1] ;
										MatrixMultiply(m2_cov_inv, v2, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
										dsq = VectorDot(v2, vtmp) ;
										p2 = scale2*exp(-0.5*dsq) ;
										ambiguity += (p1*p2) ;
									}
								}
								std1 = sqrt(gc2->covars[0]) ; std2 = sqrt(gc2->covars[2]) ; 
								Istep1 = std1/(ISTEPS/4) ; Istep2 = std2/(ISTEPS/4) ;
								for (I1 = gc2->means[0]-(ISTEPS/2)*Istep1;  I1 <= gc2->means[0]+(ISTEPS/2)*Istep1  ; I1 += Istep1)
								{
									VECTOR_ELT(v1, 1) = I1-gc1->means[0] ; VECTOR_ELT(v2, 1) = I1-gc2->means[0] ;
									for (I2 = gc2->means[1]-(ISTEPS/2)*Istep2;  I2 <= gc2->means[1]+(ISTEPS/2)*Istep2  ; I2 += Istep2)
									{
										VECTOR_ELT(v1, 2) = I2-gc1->means[1] ;
										MatrixMultiply(m1_cov_inv, v1, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
										dsq = VectorDot(v1, vtmp) ;
										p1 = scale1*exp(-0.5*dsq) ;

										VECTOR_ELT(v2, 2) = I2-gc2->means[1] ;
										MatrixMultiply(m2_cov_inv, v2, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
										dsq = VectorDot(v2, vtmp) ;
										p2 = scale2*exp(-0.5*dsq) ;
										ambiguity += (p1*p2) ;
									}
								}
								MatrixFree(&m1_cov_inv) ; MatrixFree(&m2_cov_inv) ;
								VectorFree(&v1) ; VectorFree(&v2) ; VectorFree(&vtmp) ;
							}
							else
								ErrorExit(ERROR_UNSUPPORTED, "%s: ambiguity  only supported  for  ninputs=1 or 2", Progname)  ;

							label_ambiguity[l1][l2] += ambiguity ;
							label_ambiguity[l2][l1] += ambiguity ;
							label_count[l1][l2]++;
							label_count[l2][l1]++;
							ambiguity *= ambiguity * (gcap->priors[i1]+gcap->priors[i2]) ;
							atotal += ambiguity ;
							if (ambiguity > amax)
								amax  =  ambiguity ;
						}
					}
				}
				total_ambiguity += atotal ;
 				if (amax > 0)
					MRIFvox(mri, xp, yp,  zp) = amax  ;
			}
		}
	}

	{
		FILE *fp ;
		MRI  *mri2;
		float norm;

		mri2 = MRIalloc(MAX_CMA_LABEL+1, MAX_CMA_LABEL+1, 2, MRI_FLOAT) ;
		fp = fopen("label_amb.dat", "w") ;
		for  (l1 = 0 ; l1 <= MAX_CMA_LABEL ; l1++)
		{
			for (l2 =  0  ; l2  <=  MAX_CMA_LABEL  ; l2++)
			{
#if  1
				norm = label_count[l1][l2];
#else
				norm = (label1_count[l1]+label1_count[l2])/2 ;
#endif
				if  (norm> 0)
					label_ambiguity[l1][l2] /= (float)norm ;
				fprintf(fp, "%f   ", label_ambiguity[l1][l2])  ;
				MRIFvox(mri2, l1, l2, 0) = label_ambiguity[l1][l2] ;
			}
			fprintf(fp, "\n") ;
		}
		fclose(fp)  ;
		MRIwrite(mri2, "label_amb.mgh") ;
		MRIfree(&mri2);
	}

	*pamb = total_ambiguity;
	return(mri) ;
}


int
GCAscale(GCA *gca, double min_val, double max_val)
{
	int        xp, yp, zp,i, n, r, c, k ;
	double     scale, umin[MAX_GCA_INPUTS], umax[MAX_GCA_INPUTS], range ;
	GCA_PRIOR  *gcap  ;
	GC1D       *gc ;

	for (xp = 0  ; xp  < gca->prior_width ; xp++)
	{
		for  (yp  = 0 ;  yp < gca->prior_height ; yp++)
		{
			for (zp = 0 ;  zp  < gca->prior_depth  ; zp++)
			{
				gcap =  &gca->priors[xp][yp][zp] ;
				if (xp == Gx && yp == Gy && zp == Gz)
					DiagBreak() ;
				if (gcap->nlabels <= 1)
					continue ;

				for (i = 0 ; i  <  gcap->nlabels ; i++)
				{
					gc = GCAfindPriorGC(gca, xp, yp, zp,  gcap->labels[i]) ;
					for (n = 0 ; n < gca->ninputs ; n++)
					{
						if (i == 0 && xp == 0 && yp == 0 && zp == 0)
						{
							umin[n] = gc->means[n] ; umax[n] = gc->means[n] ;
						}
						else 
						{
							if (umin[n] > gc->means[n])
								umin[n] = gc->means[n] ;
							if (umax[n] < gc->means[n])
								umax[n] = gc->means[n] ;
						}
					}
				}
			}
		}
	}

	printf("ranges:\n") ; range = 0 ;
	for (n = 0 ; n < gca->ninputs ; n++)
	{
		if ((umax[n] - umin[n]) > range)
			range = umax[n] - umin[n] ;
		printf("%d: %f --> %f\n", n, umin[n], umax[n]) ;
	}
	scale = (max_val - min_val) / range ;
	printf("scaling by %2.3f\n", scale) ;

	for (xp = 0  ; xp  < gca->prior_width ; xp++)
	{
		for  (yp  = 0 ;  yp < gca->prior_height ; yp++)
		{
			for (zp = 0 ;  zp  < gca->prior_depth  ; zp++)
			{
				gcap =  &gca->priors[xp][yp][zp] ;
				if (xp == Gx && yp == Gy && zp == Gz)
					DiagBreak() ;
				if (gcap->nlabels <= 1)
					continue ;

				for (i = 0 ; i  <  gcap->nlabels ; i++)
				{
					gc = GCAfindPriorGC(gca, xp, yp, zp,  gcap->labels[i]) ;
					for (n = 0 ; n < gca->ninputs ; n++)
					{
						gc->means[n] = (gc->means[n] - umin[n])*scale ;
						for (r = k = 0 ; r < gca->ninputs ; r++)
							for (c = r ;  c < gca->ninputs ; c++, k++)
								gc->covars[k] *= (scale*scale) ;
					}
				}
			}
		}
	}


	return(NO_ERROR) ;
}

