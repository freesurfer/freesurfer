/**
 * @brief performs a logical and at each voxel on a series of volumes
 *
 * performs a logical and at each voxel on a series of volumes
 * which must all have the same geometry/ras coords.
 */
/*
 * Original Author: Bruce Fischl
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
#include "error.h"
#include "diag.h"
#include "version.h"
#include "macros.h"

static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;

static int get_option(int argc, char *argv[]) ;

const char *Progname ;

static float thresh = -1 ;

/***-------------------------------------------------------****/
int main(int argc, char *argv[])
{
  int  nargs, index, ac, nvolumes;
  char **av ;
  MRI  *mri_and = NULL, *mri ;

  nargs = handleVersionOption(argc, argv, "mri_and");
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

  nvolumes = argc-2 ;
  if (nvolumes <= 0)
    usage_exit() ;
  printf("processing %d input files\n", nvolumes) ;

  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  for (index = 0 ; index < nvolumes ; index++)
  {
    char *fname = argv[index+1] ;
    printf("processing input volume %d of %d: %s\n",
           index+1, nvolumes, fname) ;
    mri = MRIread(fname) ;
    if (thresh >= 0)
    { 
      MRI *mri_tmp ;
      mri_tmp = MRIthreshold(mri, NULL, thresh) ;
      MRIfree(&mri) ;
      mri = mri_tmp ;
    }
    if (index == 0)
      mri_and = MRIcopy(mri, NULL) ;
    else
      MRIand(mri, mri_and, mri_and, 0) ;

    MRIfree(&mri) ;
  }

  printf("writing output to %s\n", argv[argc-1]) ;
  MRIwrite(mri_and, argv[argc-1]) ;
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
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      nargs = 0 ;
      print_usage() ;
      exit(1) ;
      break ;
    case 'T':
      thresh = atof(argv[2]) ;
      printf("applying threshold %2.3f to input data prior to anding\n", thresh);
      nargs = 1 ;
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
  printf("USAGE: %s  <options> fname1 fname2 .. \n",Progname) ;
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
    "\n"
    "Performs a logical voxel-wise AND on a series of volumes\n"
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


