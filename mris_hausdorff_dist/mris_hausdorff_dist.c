/**
 * @file  mris_distance_transform.c
 * @brief program for computing a distance transform on the surface
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:32 $
 *    $Revision: 1.3 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "label.h"
#include "MARS_DT_Boundary.h"

static char vcid[] = 
"$Id: mris_hausdorff_dist.c,v 1.3 2011/03/02 00:04:32 nicks Exp $";

int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;
static char *annot_name = NULL ;
static float binarize = 0 ;

int
main(int argc, char *argv[])
{
  char          **av ;
  int           ac, nargs, msec, mode = HDIST_MODE_SYMMETRIC_MEAN, a1, a2 ;
  LABEL         *area1, *area2 ;
  MRI_SURFACE   *mris ;
  struct timeb  then ;
  double        hdist ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mris_hausdorff_dist.c,v 1.3 2011/03/02 00:04:32 nicks Exp $", 
     "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Gdiag |= DIAG_SHOW ;
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
    usage_exit() ;

  TimerStart(&then) ;
  mris = MRISread(argv[1]) ;
  if (mris == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s",
              Progname, argv[1]) ;
  if (annot_name)
  {
    if (MRISreadAnnotation(mris, annot_name) != NO_ERROR)
      ErrorExit(ERROR_BADFILE, "%s:could not open annotation %s", Progname, annot_name);
    if (mris->ct == NULL)
      ErrorExit(ERROR_BADFILE, "Annotation %s does not contain a color table", annot_name);
    for (a1 = 0 ; a1 < mris->ct->nentries ; a1++)
    {
      area1 = MRISannotation_to_label(mris, a1) ;
      if (area1 == NULL)
        continue ;
      for (a2 = 0 ; a2 < mris->ct->nentries ; a2++)
      {
        if (a1 == a2)
          continue ;
        area2 = MRISannotation_to_label(mris, a2) ;
        if (area2 == NULL)
          continue ;
        MRISsetValues(mris, 0) ;
        MRIScopyValToVal2(mris) ;
        MRIScomputeMetricProperties(mris) ;
        MRISdistanceTransform(mris, area1, DTRANS_MODE_SIGNED) ;
        MRIScopyValToVal2(mris) ;
    
        MRIScomputeMetricProperties(mris) ;
        MRISdistanceTransform(mris, area2, DTRANS_MODE_SIGNED) ;
        hdist = MRIScomputeHausdorffDistance(mris, mode) ;
        printf("%s %s %d %d %2.3f\n", 
               mris->ct->entries[a1]->name,
               mris->ct->entries[a2]->name,
               a1,
               a2,
               hdist) ;
        fflush(stdout) ;
        LabelFree(&area2) ;
      }
      LabelFree(&area1) ; 
    }
  }
  else
  {
    if (argc < 4)
      usage_exit() ;
    area1 = LabelRead(NULL, argv[2]) ;
    if (area1 == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s",
                Progname, argv[2]) ;
    
    area2 = LabelRead(NULL, argv[3]) ;
    if (area2 == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s",
                Progname, argv[3]) ;
    

    if (binarize > 0)
    {
      LabelThreshold(area1, binarize) ;
      LabelThreshold(area2, binarize) ;
    }
    MRIScomputeMetricProperties(mris) ;
    MRISdistanceTransform(mris, area1, DTRANS_MODE_SIGNED) ;
    MRIScopyValToVal2(mris) ;
    
    MRIScomputeMetricProperties(mris) ;
    MRISdistanceTransform(mris, area2, DTRANS_MODE_SIGNED) ;
    hdist = MRIScomputeHausdorffDistance(mris, mode) ;
    printf("%2.3f\n", hdist) ;
  }

  msec = TimerStop(&then) ;
  fprintf(stderr,"hausdorff transform took %2.1f minutes\n", (float)msec/(60*1000.0f));

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
  else switch (toupper(*option))
    {
    case 'B':
      binarize = atof(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'A':
      annot_name = argv[2] ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      print_help() ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void)
{
  printf("usage: %s [options] <surface> <label1> <label2>\n",
         Progname) ;
  printf("\t-a <annot name>    compute pairwise HD between all annotations\n");
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr,
          "This program computes the Hausdorff distance between two labels on "
          "the surface\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "(see the source code!)\n") ;
  exit(1) ;
}


static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

