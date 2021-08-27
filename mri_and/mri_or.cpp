/**
 * @brief performs a logical or at each voxel on a series of volumes
 *
 * performs a logical or at each voxel on a series of volumes
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
int use_orig_value = 0;

/***-------------------------------------------------------****/
int main(int argc, char *argv[])
{
  int  nargs, index, ac, nvolumes;
  char **av ;
  MRI  *mri_or = NULL, *mri ;

  nargs = handleVersionOption(argc, argv, "mri_or");
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
    if (index == 0){
      mri_or = MRIcopy(mri, NULL) ;
    // if nvolumes == 1 binarize the volume! LZ: MRIbinarize(MRI *mri_src, MRI *mri_dst, float threshold, float low_val,float hi_val)
      if (nvolumes == 1) {
	if(use_orig_value)
	  MRIorVal(mri, mri_or, mri_or, 0) ;
	else
	  MRIor(mri, mri_or, mri_or, 0) ;
      }
    } 
    else {
      if(use_orig_value)
	MRIorVal(mri, mri_or, mri_or, 0) ;
      else
	MRIor(mri, mri_or, mri_or, 0) ;
    }


    MRIfree(&mri) ;
  }

  printf("writing output to %s\n", argv[argc-1]) ;
  MRIwrite(mri_or, argv[argc-1]) ;
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
    case 'V':
      print_version() ;
      break ;
    case 'O':
      use_orig_value = 1;
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
	 "Performs a logical voxel-wise OR on a series of volumes\n"
	 "\n"
	 );
  printf("Options:\n\n") ;
  printf("  -o : The logical OR operation will keep the original label values in the input files when creating the output.\n");
  printf("  -help : print usage\n");
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


