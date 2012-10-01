/**
 * @file  mris_sample_label.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2012/10/01 18:59:32 $
 *    $Revision: 1.5 $
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
#include "mrisurf.h"
#include "macros.h"
#include "fio.h"
#include "label.h"

static char vcid[] = "$Id: mris_sample_label.c,v 1.5 2012/10/01 18:59:32 fischl Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

/*-------------------------------- DATA ----------------------------*/

char *Progname ;
static int voxval = 1 ;

/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  MRI_SURFACE  *mris ;
  char         **av, *in_label_fname, *out_label_fname, *surf_fname, ext[STRLEN] ; ;
  int          ac, nargs ;
  LABEL        *label, *label_out ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 4)
    usage_exit() ;

  in_label_fname = argv[1] ;
  surf_fname = argv[2] ;
  out_label_fname = argv[3] ;

  printf("reading label from %s...\n", in_label_fname) ;
  if (!strcmp(FileNameExtension(in_label_fname, ext), "mgz"))
  {
    MRI *mri = MRIread(in_label_fname) ;
    printf("creating label from volumetric inputs with voxval = %d\n", voxval) ;
    if (mri == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read input volume from %s", Progname, in_label_fname);
    label = LabelfromASeg(mri, voxval) ;
    MRIfree(&mri) ;
  }
  else
  {
    label = LabelRead(NULL, in_label_fname) ;
    if (!label)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s", Progname, in_label_fname) ;
  }
  printf("reading surface from %s...\n", surf_fname) ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",Progname, surf_fname) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;

#if 0
  LabelFillUnassignedVertices(mris, label) ;
#else
  label_out = LabelFillHoles(label, mris, ORIGINAL_VERTICES) ;
#endif
  printf("writing sampled label to %s...\n", out_label_fname) ;
  LabelWrite(label_out, out_label_fname) ;
  MRISfree(&mris) ;
  LabelFree(&label) ;

  exit(0) ;
  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option)) {
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <input label file> <input surface file> <output label file>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will sample a label onto a surface model.\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
