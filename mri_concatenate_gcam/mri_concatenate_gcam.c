/**
 * @file  mri_concatenate_gcam.c
 * @brief concatenates 2 gcam (nonlinear) morphs
 *
 */
/*
 * Original Author: Lilla Zollei
 * CVS Revision Info:
 *    $Author: lzollei $
 *    $Date: 2011/11/01 14:24:39 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "gcamorph.h"
#include "mri.h"
#include "matrix.h"
#include "proto.h"
#include "macros.h"
#include "error.h"
#include "timer.h"
#include "diag.h"
#include "mrimorph.h"
#include "utils.h"
#include "gca.h"
#include "cma.h"
#include "version.h"
#include "transform.h"
#include "fastmarching.h"
#include "voxlist.h"
#include "mrisurf.h"


static void  usage_exit(int ecode) ;
static int get_option(int argc, char *argv[]) ;

char *Progname ;
int DoWriteInverse = 0;
int InvertGcam1 = 0;
int InvertGcam2 = 0;

int
main(int argc, char *argv[])
{
  char         *gcam1_fname, *gcam2_fname, *gcam_composed_fname, **av;
  GCA_MORPH    *gcam1, *gcam2, *gcam_composed = NULL;
  int          ac, nargs;

  Progname = argv[0] ;
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
    {
      nargs = get_option(argc, argv) ;
      argc -= nargs ;
      argv += nargs ;
    }
  
  if (argc < 4)
    usage_exit(1) ;
  
  gcam1_fname         = argv[1] ;
  gcam2_fname         = argv[2] ;
  gcam_composed_fname = argv[3] ;

  if(InvertGcam1)
    {
      gcam1 = GCAMreadAndInvert(gcam1_fname); // TODO: DON'T use that. Not general enough! -- assumes name talairach.m3z
    }
  else
    gcam1 = GCAMread(gcam1_fname); 

  if (!gcam1)
    ErrorExit(ERROR_NOFILE, "%s: could not read 1st morph %s",
	      Progname, gcam1_fname) ;

  if(InvertGcam2)
    {
      gcam2 = GCAMreadAndInvert(gcam2_fname);
    }
  else
    gcam2 = GCAMread(gcam2_fname); 
  if (!gcam2)
    ErrorExit(ERROR_NOFILE, "%s: could not read 2nd morph %s",
	      Progname, gcam2_fname) ;

  gcam_composed = GCAMalloc(gcam1->width, gcam1->height, gcam1->depth) ;
  gcam_composed->image = gcam1->image; 
  gcam_composed->atlas = gcam2->atlas; 
  GCAMconcatenate(gcam1, gcam2, gcam_composed);

  if(DoWriteInverse)
    {
      printf("Writing inverse: not yet implemented\n") ;
    }
  else
    {
      printf("Writing the concatenated warp to %s\n", gcam_composed_fname) ;
      GCAMwrite(gcam_composed, gcam_composed_fname) ;
      //printf("Done writing\n") ;
    }

  printf("Done in %s\n", Progname) ;
  exit(0) ;
  return(0) ;
}


static void 
usage_exit(int ecode)
{
  printf("usage: %s [options] <gcam1> <gcam2> <gcam_combined>\n",
	 Progname) ;
  printf("\n") ;
  printf("Options  -- PLACEHOLDER -- NOT YET IMPLEMENTED:\n\n") ;
  printf("  -i1 (default = 0 / FALSE)\n");
  printf("  -i2 (default = 0 / FALSE)\n");
  printf("  -write_invert (default = 0 / FALSE)\n");

  printf("  -?			: print usage\n");
  printf("  -help         	: print usage\n");
  printf("  -U			: print usage \n");

  exit(ecode) ;
}

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  StrLower(option) ;
  if (!stricmp(option, "write_invert"))
    {
      DoWriteInverse = 1 ;
      if(DoWriteInverse)
	printf("Will write the inverse of the concatenated morphs -- PLACEHOLDER -- NOT YET IMPLEMENTED.\n") ;
    }
  else if (!stricmp(option, "i1"))
    {
      InvertGcam1 = 1 ;
      if(InvertGcam1)
	printf("Will invert 1st morph before concatenating with 2nd morph -- PLACEHOLDER -- NOT YET IMPLEMENTED.\n") ;
    }
  else if (!stricmp(option, "i2"))
    {
      InvertGcam2 = 1 ;
      if(InvertGcam2)
	printf("Will invert 2nd morph before concatenating with 1st morph -- PLACEHOLDER -- NOT YET IMPLEMENTED.\n") ;
    }
  else if (!stricmp(option, "help"))
    {
      usage_exit(1);
    }
  else switch (*option)
    {
    case '?':
    case 'u':
      usage_exit(1);
      break ;
    default:
      printf("Unknown option %s\n", argv[1]) ;
      usage_exit(1) ;
      break ;
    }
  return(nargs) ;
}
