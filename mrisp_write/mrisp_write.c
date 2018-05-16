/**
 * @file  mrisp_write.c
 * @brief extracts an array ("a variable") from surface-registration template
 *
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2016/03/22 14:47:57 $
 *    $Revision: 1.12 $
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
#include "mri.h"
#include "macros.h"
#include "version.h"

static char vcid[] = "$Id: mrisp_write.c,v 1.12 2016/03/22 14:47:57 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;

static int normalize = 0 ;
static int navgs = 0 ;


static char subjects_dir[STRLEN] ;

static float scale = 1 ;
static int nframes = 1 ;

int
main(int argc, char *argv[])
{
  char         **av, *out_fname;
  int          ac, nargs ;
  char         *in_surf, *in_overlay;
  MRI_SURFACE  *mris;
  MRI_SP       *mrisp ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mrisp_write.c,v 1.12 2016/03/22 14:47:57 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
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

  if (argc < 4)
  {
    usage_exit() ;
  }

  in_surf = argv[1] ;
  in_overlay = argv[2] ;
  out_fname = argv[3] ;

  mrisp = MRISPalloc(scale, 1) ;
  fprintf(stderr, "reading surface from %s...\n", in_surf) ;
  mris = MRISread(in_surf) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s", Progname, in_surf) ;

  printf("reading overlay from %s\n", in_overlay) ;
  if (MRISreadCurvatureFile(mris, in_overlay) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read input overlay %s", Progname, in_overlay) ;

  if (normalize)
    MRISnormalizeCurvature(mris, NORM_MEAN);

  MRISaverageCurvatures(mris, navgs) ;

  MRIStoParameterization(mris, mrisp, scale, 0) ;

  printf("writing output file to %s\n", out_fname) ;

  MRISPwrite(mrisp, out_fname) ;

  MRISPfree(&mrisp) ;
  MRISfree(&mris) ;
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
  int    nargs = 0 ;
  char   *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    print_help() ;
  }
  else if (!stricmp(option, "-version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "SDIR"))
  {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  }
  else if (!stricmp(option, "NFRAMES")) // not implemented yet
  {
    nframes = atoi(argv[2]) ;
    nargs = 1 ;
    printf("writing out %d frames - NOT IMPLEMENTED YET\n", nframes) ;
    exit(1) ;
  }
  else if (!stricmp(option, "scale"))
  {
    scale = atof(argv[2]);
    printf("scaling width/height MRISP by %2.1f\n", scale) ;
    nargs=1;
  }
  else switch (toupper(*option))
    {
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "averaging overlay %d times...\n", navgs) ;
      break ;
    case 'N':
      normalize = 1 ;
      fprintf(stderr, "normalizing curvature by variance.\n") ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
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
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

#include "mrisp_write.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mrisp_write_help_xml,
                mrisp_write_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

