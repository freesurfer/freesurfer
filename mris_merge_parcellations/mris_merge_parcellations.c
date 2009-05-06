/**
 * @file  mris_merge_parcellations
 * @brief program for merging two parcellations into one, taking some divisions from one and some from the other
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2009/05/06 01:35:12 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2004-2007,
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
"$Id: mris_merge_parcellations.c,v 1.2 2009/05/06 01:35:12 fischl Exp $";

int main(int argc, char *argv[]) ;
static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int merge_annotations(COLOR_TABLE *ct, MRI_SURFACE *mris1, MRI_SURFACE *mris2, MRI_SURFACE *mris) ;

const char *Progname ;
static char          fsdir[STRLEN] = "" ;

int
main(int argc, char *argv[])
{
  char          **av, *parc1, *parc2, *oname, surf_name[STRLEN], path[STRLEN],hemi[STRLEN],*cp,fname[STRLEN];
  int           ac, nargs, msec ;
  MRI_SURFACE   *mris1, *mris2 ;
  struct timeb  then ;
  COLOR_TABLE   *ct ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option 
    (argc, argv, 
     "$Id: mris_merge_parcellations.c,v 1.2 2009/05/06 01:35:12 fischl Exp $", 
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

  if (argc < 4)
    usage_exit() ;
  TimerStart(&then) ;

  if (strlen(fsdir) == 0)
  {
    cp = getenv("FREESURFER_HOME") ;
    if (cp == NULL)
      ErrorExit(ERROR_BADPARM, "FRESURFER_HOME must be defined in the environment") ;
    strcpy(fsdir, cp) ;
  }
  sprintf(fname, "%s/FreeSurferColorLUT.txt", fsdir) ;
  ct = CTABreadASCII(fname) ;
  if (ct == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read color table from %s", Progname, fname) ;

  parc1 = argv[1] ;
  parc2 = argv[2] ;
  oname = argv[3] ;
  FileNamePath(parc1, path);
  FileNameOnly(parc1, fname) ;
  cp = strstr(fname, "h.") ;
  if (cp == NULL)
    ErrorExit(ERROR_UNSUPPORTED, "%s: could not scan hemisphere from fname %s", Progname, fname) ;
  strncpy(hemi, cp-1, 2) ; hemi[2] = 0 ;
  sprintf(surf_name, "%s/../surf/%s.orig", path, hemi) ;
  mris1 = MRISread(surf_name) ;
  if (mris1 == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, surf_name) ;
  if (MRISreadAnnotation(mris1, parc1) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s:could not open annotation %s", Progname, parc1);
  if (mris1->ct == NULL)
    ErrorExit(ERROR_BADFILE, "Annotation %s does not contain a color table", parc1);

  mris2 = MRISread(surf_name) ;
  if (mris2 == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", Progname, surf_name) ;
  if (MRISreadAnnotation(mris2, parc2) != NO_ERROR)
    ErrorExit(ERROR_BADFILE, "%s:could not open annotation %s", Progname, parc2);
  if (mris2->ct == NULL)
    ErrorExit(ERROR_BADFILE, "Annotation %s does not contain a color table", parc2);

  merge_annotations(ct, mris1, mris2, mris2) ;

  printf("writing output parcellation to %s\n", oname) ;
  MRISwriteAnnotation(mris2, oname) ;
  msec = TimerStop(&then) ;
  printf("parcellation merging took %2.1f minutes\n", (float)msec/(60*1000.0f));

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
  else if (!stricmp(option, "fsdir"))
  {
    strcpy(fsdir, argv[2]) ;
    printf("using %s as FREESURFER_HOME", fsdir) ;
  }
  else switch (toupper(*option))
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
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

/*
  these are actually the same surface, the 1st one has Rahul's parcellation, the
  2nd one Christophe's. Grab the cingulate subdivisions from Rahul's and put them in Christophe's.
*/
static int
merge_annotations(COLOR_TABLE *ct, MRI_SURFACE *mris1, MRI_SURFACE *mris2, MRI_SURFACE *mris)
{
  int       vno, annot, s_cingulate, g_cingulate, s_pericallosal ;
  int       caudal_acc, posterior_cingulate, rostral_acc ;
  int       s_caudal_acc, s_posterior_cingulate, s_rostral_acc ;
  int       g_caudal_acc, g_posterior_cingulate, g_rostral_acc ;
  VERTEX    *v ;

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    CTABprintASCII(ct, stdout) ;
  if (mris2->hemisphere == RIGHT_HEMISPHERE)
  {
    caudal_acc = CTABentryNameToAnnotation("ctx-rh-caudalanteriorcingulate", ct) ;
    posterior_cingulate = CTABentryNameToAnnotation("ctx-rh-posteriorcingulate", ct) ;
    rostral_acc = CTABentryNameToAnnotation("ctx-rh-rostralanteriorcingulate", ct) ;

    s_caudal_acc = CTABentryNameToAnnotation("ctx-rh-S_cingulate-caudal_ACC", ct) ;
    s_rostral_acc = CTABentryNameToAnnotation("ctx-rh-S_cingulate-rostral_ACC", ct);
    s_posterior_cingulate = CTABentryNameToAnnotation("ctx-rh-S_cingulate-posterior", ct) ;
    g_caudal_acc = CTABentryNameToAnnotation("ctx-rh-G_cingulate-caudal_ACC", ct) ;
    g_rostral_acc = CTABentryNameToAnnotation("ctx-rh-G_cingulate-rostral_ACC", ct);
    g_posterior_cingulate = CTABentryNameToAnnotation("ctx-rh-G_cingulate-posterior", ct) ;


    s_cingulate = CTABentryNameToAnnotation("ctx-lh-S_cingulate-Main_part_and_Intracingulate", ct) ;
    g_cingulate = CTABentryNameToAnnotation("ctx-rh-G_cingulate-Main_part", ct) ;
    s_pericallosal = CTABentryNameToAnnotation("ctx-rh-S_pericallosal", ct) ;
  }
  else // left hemi
  {
    caudal_acc = CTABentryNameToAnnotation("ctx-lh-caudalanteriorcingulate", ct) ;
    posterior_cingulate = CTABentryNameToAnnotation("ctx-lh-posteriorcingulate", ct) ;
    rostral_acc = CTABentryNameToAnnotation("ctx-lh-rostralanteriorcingulate", ct) ;

    s_cingulate = CTABentryNameToAnnotation("ctx-lh-S_cingulate-Main_part_and_Intracingulate", ct) ;
    g_cingulate = CTABentryNameToAnnotation("ctx-lh-G_cingulate-Main_part", ct) ;
    s_pericallosal = CTABentryNameToAnnotation("ctx-lh-S_pericallosal", ct) ;

    s_caudal_acc = CTABentryNameToAnnotation("ctx-lh-S_cingulate-caudal_ACC", ct) ;
    s_rostral_acc = CTABentryNameToAnnotation("ctx-lh-S_cingulate-rostral_ACC", ct);
    s_posterior_cingulate = CTABentryNameToAnnotation("ctx-lh-S_cingulate-posterior", ct) ;
    g_caudal_acc = CTABentryNameToAnnotation("ctx-lh-G_cingulate-caudal_ACC", ct) ;
    g_posterior_cingulate = CTABentryNameToAnnotation("ctx-lh-G_cingulate-posterior", ct);
    g_rostral_acc = CTABentryNameToAnnotation("ctx-lh-G_cingulate-rostral_ACC", ct) ;
  }

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    annot = mris2->vertices[vno].annotation ;
    if  (mris1->vertices[vno].annotation == caudal_acc)
    {
      if (annot == s_cingulate /* || annot == s_pericallosal*/)
        annot = s_caudal_acc ; 
      else if (annot == g_cingulate)
        annot = g_caudal_acc ; 
    }
    else if  (mris1->vertices[vno].annotation == rostral_acc)
    {
      if (annot == s_cingulate /* || annot == s_pericallosal*/)
        annot = s_rostral_acc ;
      else if (annot == g_cingulate)
        annot = g_rostral_acc ; 
    }
    else if (mris1->vertices[vno].annotation == posterior_cingulate)
    {
      if (annot == s_cingulate /* || annot == s_pericallosal*/)
        annot = s_posterior_cingulate ; 
      else if (annot == g_cingulate)
        annot = g_posterior_cingulate ; 
    }
    if (annot < 0)
      DiagBreak() ;
    v->annotation = annot ;
  }

  mris->ct = ct ;
  return(NO_ERROR) ;
}

