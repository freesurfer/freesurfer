/**
 * @brief computes the mi between two input volumes
 *
 * computes the mutual information (mi) between two input volumes
 * which must all have the same geometry/ras coords.
 */

/*
 * Original Author: Lilla Zollei
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer Software License Agreement' contained
 * in the file 'LICENSE' found in the FreeSurfer distribution, and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferSoftwareLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

#include <stdio.h>
#include <sys/stat.h>
#include <unistd.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include "const.h"
#include "utils.h"
#include "mri.h"
#include "mri2.h"
#include "error.h"
#include "diag.h"
#include "version.h"
#include "macros.h"
#include "mri_mi.h"

static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;

static int get_option(int argc, char *argv[]) ;

const char *Progname = NULL;
int silent_mode = 0;

static int bins_1 = 64;
static int bins_2 = 64;

/***-------------------------------------------------------****/
int main(int argc, char *argv[])
{
  int  nargs, index, ac, nvolumes;
  char **av ;
  MRI  *mri_1, *mri_2 ;

  // printf("mri_mi 1\n") ;

  nargs = handleVersionOption(argc, argv, "mri_mi");
  if (nargs && argc - nargs == 1)
    exit (0);
  Progname = argv[0] ;
  argc -= nargs;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  //printf("mri_mi 2\n") ;

  nvolumes = argc-2 ;
  if (nvolumes <= 0)
    usage_exit() ;
  if(!silent_mode)
    printf("processing %d input files\n", nvolumes) ;
  
  //printf("mri_mi 3\n") ;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  
  //printf("mri_mi 4\n") ;

  // Input volume #1
  index = 0;
  char *fname = argv[index+1] ;
  if(!silent_mode)
    printf("processing input volume %d of %d: %s\n",
	   index+1, nvolumes, fname) ;
  mri_1 = MRIread(fname) ;
  if(mri_1 == NULL) exit(1);
  //MRIwrite(mri_1,"/tmp/mri_1.mgz");

  //printf("mri_mi 5\n") ;

  // Input volume #2
  index++;
  fname = argv[index+1] ;
  if(!silent_mode)
    printf("processing input volume %d of %d: %s\n",
	   index+1, nvolumes, fname) ;
  mri_2 = MRIread(fname) ;
  if(mri_2 == NULL) exit(1);
  //MRIwrite(mri_2,"/tmp/mri_2.mgz");

  //printf("mri_mi 6\n") ;

  float min1 = 0, max1 = 0, min2 = 0, max2 = 0;
  mri_minmax(mri_1, &min1, &max1);
  mri_minmax(mri_2, &min2, &max2);
  HISTOGRAM* histo1 = HISTOalloc(bins_1);
  HISTOGRAM* histo2 = HISTOalloc(bins_2);
  HISTOinit(histo1, bins_1, min1, max1);  // grab min-max from image!
  HISTOinit(histo2, bins_2, min2, max2);  // grab min-max from image!

  //printf("mri_mi 7\n") ;

  int width1 = mri_1->width;
  int height1 = mri_1->height ;
  int depth1 = mri_1->depth;
  // int N = width1*height1*depth1;
  int N = depth1;
  //printf("mri_mi 1: N  = %d \n", N) ;
  double samples1[N];
  int i, j, k;
  for(i = 0; i < width1; i++)
    for(j = 0; j < height1; j++)
      {
	for(k = 0; k < depth1; k++)
	  {
	    //printf("mri_mi 1: (i, j, k)  = (%d, %d, %d) \n", i, j, k) ;
	    //MRIsampleVolume(mri_1, i, j, k, &samples1[i*(height1*depth1)+j*depth1+(k)]); 
	    MRIsampleVolume(mri_1, i, j, k, &samples1[k]); 
	  }
	HISTOcount(histo1, samples1, N);
      }
  
  //printf("mri_mi 8\n") ;

  int width2 = mri_2->width;
  int height2 = mri_2->height ;
  int depth2 = mri_2->depth;
  // N = width2*height2*depth2;
  N = depth2;
  // printf("mri_mi 2: N  = %d \n", N) ;
  double samples2[N];
  for(i = 0; i < width2; i++)
    for(j = 0; j < height2; j++)
      {
	for(k = 0; k < depth2; k++)
	  {
	    //printf("mri_mi 2: (i, j, k)  = (%d, %d, %d) \n", i, j, k) ;
	    //MRIsampleVolume(mri_2, i, j, k, &samples2[i*(height2*depth2)+j*depth2+(k)]);
	    MRIsampleVolume(mri_2, i, j, k, &samples2[k]);
	  }
	HISTOcount(histo2, samples2, N);
      }

  //printf("mri_mi 9\n") ;

  JOINT_HISTOGRAM* jhisto = JHISTOalloc(bins_1, bins_2);
  JHISTOfill(mri_1, mri_2, jhisto);

  //printf("mri_mi 10\n") ;

  double marginalentropy1 =   HISTOgetEntropy(histo1);
  double marginalentropy2 =   HISTOgetEntropy(histo2);
  double jointentropy     =   JHISTOgetEntropy(jhisto);
  double mi_score = marginalentropy1 + marginalentropy2 - jointentropy;
  //printf("MI SCORE: %f\n", mi_score) ;

  MRIfree(&mri_1) ;
  MRIfree(&mri_2) ;
  HISTOfree(&histo1);
  HISTOfree(&histo2);
  JHISTOfree(&jhisto);

  //printf("mri_mi 11\n") ;

  if(silent_mode)
    printf("%f\n", mi_score) ;
  else
    printf("The mutual information between the input volumes: %f\n", mi_score) ;
  exit(0);

} /* end main() */


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
  {
    print_help() ;
  }
  else if (!stricmp(option, "-bins"))
  { 
    bins_1 = atoi(argv[2]) ;
    bins_2 = atoi(argv[3]) ;
    nargs = 2 ;
    // printf("Histogram bins = (%d, %d)\n", bins_1, bins_2);
  }
  else if (!stricmp(option, "-silent"))
    { 
      silent_mode = 1;
    }
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      nargs = 0 ;
      print_usage() ;
      exit(1) ;
      break ;
    case 'V':
      print_version() ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s  <options> fname1 fname2 \n",Progname) ;
  printf("\n");
  printf("Options:\n\n") ;
  printf("  --bins bin1 bin2 (default = 64x64)\n");
  printf("  --silent: only write out final MI result (default = 0) \n");
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
    "\n"
    "Computes mutual information (mi) between two input volumes.\n"
  );
  exit(1) ;
}

/* --------------------------------------------- */
static void print_version(void)
{
  std::cout << getVersion() << std::endl;
  exit(1) ;
}

/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
