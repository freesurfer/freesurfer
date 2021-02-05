/**
 * @brief extracts an array ("a variable") from surface-registration template
 *
 */
/*
 * Original Author: Bruce Fischl
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


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static int coords = -1 ;
static int normalize = 0 ;
static int variance =  0 ;
static int navgs = 0 ;
static int sqrt_flag = 0 ;

const static char *surface_names[] =
{
  "inflated",
  "smoothwm",
  "smoothwm"
} ;

static int field_no = -1;
static char subjects_dir[STRLEN] ;
static char *hemi=NULL;
static char *subject_name;

static int frame_number = 0  ;
static int nframes = 1 ;

int
main(int argc, char *argv[])
{
  char         **av, *surf_fname, *template_fname, *out_fname, *cp;
  int          n,ac, nargs ;
  float        sse,var;
  char fname[STRLEN];
  MRI_SURFACE  *mris ,*mris_var;
  MRI_SP       *mrisp ;
  VERTEX *v;

  nargs = handleVersionOption(argc, argv, "mrisp_paint");
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

  template_fname = argv[1] ;
  surf_fname = argv[2] ;
  out_fname = argv[3] ;

  fprintf(stderr, "reading surface from %s...\n", surf_fname) ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, surf_fname) ;

  if (normalize)
  {
    cp = strchr(template_fname, '#') ;
    if (cp)   /* # explicitly given */
    {
      frame_number = atoi(cp+1) ;
      *cp = 0 ;
    }
  }
  else
  {
    cp = strchr(template_fname, '#') ;
    if (cp)   /* # explicitly given */
    {
      frame_number = atoi(cp+1) ;
      *cp = 0 ;
    }
  }
  fprintf(stderr, "reading template parameterization from %s...\n",
          template_fname) ;
  mrisp = MRISPread(template_fname) ;
  if (!mrisp)
    ErrorExit(ERROR_NOFILE, "%s: could not open template file %s",
              Progname, template_fname) ;

  if (coords >= 0)
  {
#if 1
    MRIScoordsFromParameterizationBarycentric(mris, mrisp, coords) ;
#else
    int vno ;
    MRISfromParameterization(mrisp, mris, 0) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
      mris->vertices[vno].whitex = mris->vertices[vno].curv ;
    MRISfromParameterization(mrisp, mris, 1) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
      mris->vertices[vno].whitey = mris->vertices[vno].curv ;
    MRISfromParameterization(mrisp, mris, 2) ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
      mris->vertices[vno].whitez = mris->vertices[vno].curv ;
#endif
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
    printf("writing surface to %s\n", out_fname);
    MRISwrite(mris,out_fname) ;
    exit(0) ;
  }
  if (normalize)
  {
    MRISnormalizeFromParameterization(mrisp, mris, frame_number) ;
  }
  else
  {
    MRISfromParameterization(mrisp, mris, frame_number) ;
  }

  MRISaverageCurvatures(mris, navgs) ;

  if (variance)
  {
    /* check if SUBJECTS_DIR is set up */
    if (!strlen(subjects_dir)) /* hasn't been set on command line */
    {
      cp = getenv("SUBJECTS_DIR") ;
      if (!cp)
        ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                  Progname);
      strcpy(subjects_dir, cp) ;
    }
    mris_var=MRISclone(mris);
    /* reading or generating the field */
    if (ReturnFieldName(field_no))    /* read in precomputed curvature file */
    {
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", subjects_dir,subject_name,hemi,ReturnFieldName(field_no)) ;
      if (req >= STRLEN) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (MRISreadCurvatureFile(mris_var, fname) != NO_ERROR)
      {
        MRISfree(&mris_var);
        ErrorExit(ERROR_BADPARM,"%s: could not load file %s\n",fname);
      }
    }
    else                           /* compute curvature of surface */
    {
      int req = snprintf(fname, STRLEN, "%s/%s/surf/%s.%s", subjects_dir,subject_name,hemi,surface_names[field_no]) ;
      if (req >= STRLEN) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      /*    if(parms->fields[n].field==0) */
      /*     sprintf(fname, "inflated") ; */
      /*    else */
      /*     sprintf(fname, "smoothwm") ; */
      MRISsaveVertexPositions(mris_var, CANONICAL_VERTICES) ;
      if (MRISreadVertexPositions(mris_var, fname) != NO_ERROR)
      {
        MRISfree(&mris_var);
        ErrorExit(ERROR_BADPARM,"%s: could not load file %s\n",fname);
      }
      MRISsetNeighborhoodSizeAndDist(mris_var, -1) ;  /* back to max */
      MRIScomputeMetricProperties(mris_var) ;
      MRIScomputeSecondFundamentalForm(mris_var) ;
      MRISuseMeanCurvature(mris_var) ;
      MRISaverageCurvatures(mris_var, navgs) ;
      MRISresetNeighborhoodSize(mris_var,1);/*only use nearest neighbor distances*/
      MRISrestoreVertexPositions(mris_var, CANONICAL_VERTICES) ;
    }
    MRISnormalizeField(mris_var,IsDistanceField(field_no), NORM_MEAN) ;
    MRISnormalizeField(mris,IsDistanceField(field_no), NORM_MEAN);
    /* save curv into curvbak*/
    for (n=0 ; n < mris->nvertices; n++)
    {
      v=&mris_var->vertices[n];
      v->curvbak=v->curv;
    }
    /* computing variance */
    MRISfromParameterization(mrisp, mris_var, frame_number+1) ;
    for (sse=0.0f,n=0 ; n < mris->nvertices; n++)
    {
      v=&mris_var->vertices[n];
      var=MAX(0.01,v->curv);
      v->curv = SQR(v->curvbak-mris->vertices[n].curv)/var;
      sse += v->curv;
    }
    fprintf(stderr,"XXXXXXXXXXXXXXXXXXX\n\n");
    fprintf(stderr,"The SSE for this field is %f\n",sqrt(sse/(float)mris->nvertices));
    fprintf(stderr,"\nXXXXXXXXXXXXXXXXXXX\n\n");

    MRISfree(&mris_var);
  }

  fprintf(stderr, "writing curvature file to %s...\n", out_fname) ;
  if (sqrt_flag)
  {
    MRIScopyCurvatureToValues(mris) ;
    MRISsqrtVal(mris) ;
    MRIScopyCurvatureFromValues(mris) ;
  }

  MRISwriteCurvature(mris, out_fname) ;

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
  else if (!stricmp(option, "COORDS"))
  {
    if (!stricmp(argv[2], "white"))
      coords = WHITE_VERTICES ;
    else if (!stricmp(argv[2], "pial"))
      coords = PIAL_VERTICES ;
    else
      ErrorExit(ERROR_UNSUPPORTED, "Unknown coords value %s", argv[2]) ;
    nargs = 1 ;
    printf("writing coords %s (%d) into parameterizaton\n", argv[2], coords) ;
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
  else if (!stricmp(option, "variance"))
  {
    variance=1;
    hemi = argv[3];
    subject_name = argv[2];
    field_no = atoi(argv[4]);
    if (field_no<0)
    {
      fprintf(stderr,"Incorrect Field Number\n");
      exit(-1);
    }
    fprintf(stderr,"generating variance map for with field %d (%s/surf/%s.%s) \n",field_no,subject_name,hemi,ReturnFieldName(field_no));
    nargs=3;
  }
  else switch (toupper(*option))
    {
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "averaging curvature patterns %d times...\n", navgs) ;
      break ;
    case 'F':
      frame_number = atoi(argv[2]) ;
      nargs = 1 ;
      printf("writing out frame %d\n", frame_number) ;
      break ;
    case 'N':
      normalize = 1 ;
      fprintf(stderr, "normalizing curvature by variance.\n") ;
      break ;
    case 'W':
      Gdiag |= DIAG_WRITE ;
      break ;
    case 'S':
      printf("taking sqrt before writing...\n") ;
      sqrt_flag = 1 ;
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

#include "mrisp_paint.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mrisp_paint_help_xml,
                mrisp_paint_help_xml_len);
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
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

