
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
#include "annotation.h"

static char vcid[] = "$Id: mris_sample_parc.c,v 1.2 2002/01/24 17:48:12 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static int  translate_indices_to_annotations(MRI_SURFACE *mris, char *translation_fname) ;


char *Progname ;
static int mode_filter = 0 ;
static char *surf_name = WHITE_MATTER_NAME ;
static char *thickness_name = "thickness" ;
static char sdir[STRLEN] ;
static char *translation_fname = "cma_parcellation_colors.txt" ;
static int wsize = 7 ;

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *subject_name, *cp, fname[STRLEN], *parc_name, *annot_name ;
  int           ac, nargs, vno ;
  MRI_SURFACE   *mris ;
  MRI           *mri_parc ;
  VERTEX        *v ;
  double        d ;
  Real          x, y, z, xw, yw, zw ;

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

  subject_name = argv[1] ;
  hemi = argv[2] ; 
  parc_name = argv[3] ;
  annot_name = argv[4] ;

  if (strlen(sdir) == 0)  /* if not specified explicitly as option */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, 
                "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
    strcpy(sdir, cp) ;
  }

  sprintf(fname, "%s/%s/mri/%s", sdir, subject_name, parc_name) ;
  printf("reading parcellation volume from %s...\n", fname) ;
  mri_parc = MRIread(fname) ;
  if (!mri_parc)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s",
              Progname, fname) ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, subject_name, hemi, surf_name) ;
  printf("reading input surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  if (MRISreadCurvatureFile(mris, thickness_name) != NO_ERROR)
    ErrorExit(ERROR_NOFILE, "%s: could not read thickness file %s",
              Progname, thickness_name) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    d = v->curv*.5 ;  /* halfway out */
    x = v->x+d*v->nx ; y = v->y+d*v->ny ; z = v->z+d*v->nz ;
    MRIworldToVoxel(mri_parc, x, y, z, &xw, &yw, &zw) ;
    v->annotation = v->val = 
      MRIfindNearestNonzero(mri_parc, wsize, xw, yw, zw) ;
#if 0
    v->val = v->annotation = MRIvox(mri_parc, nint(xw), nint(yw), nint(zw)) ;
#endif
  }
  if (mode_filter)
  {
    printf("mode filtering sample labels...\n") ;
#if 0
    MRISmodeFilterZeroVals(mris) ;
#else
    MRISmodeFilterVals(mris, mode_filter) ;
#endif
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      v->annotation = v->val ;
    }
  }

  /* this will fill in the v->annotation field from the v->val ones */
  translate_indices_to_annotations(mris, translation_fname) ;

  sprintf(fname, "%s-%s.annot", hemi, annot_name) ;
  MRISwriteAnnotation(mris, fname) ; 
  MRISreadAnnotation(mris, fname) ;
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
  else if (!stricmp(option, "sdir"))
  {
    strcpy(sdir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
  }
  else switch (toupper(*option))
  {
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'F':
    mode_filter = atoi(argv[2]) ;
    nargs = 1 ;
    printf("applying mode filter %d times to parcellation\n", mode_filter) ; 
    break ;
  case 'W':
    wsize = atoi(argv[2]) ;
    nargs = 1 ;
    printf("using window size=%d for sampling\n", wsize) ;
    break ;
  case 'T':
    thickness_name = argv[2] ;
    nargs = 1 ;
    printf("using thickness file %s\n", thickness_name) ;
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
usage_exit(void)
{
  print_help() ;
  exit(1) ;
}

static void
print_usage(void)
{
  fprintf(stderr, 
          "usage: %s [options] <subject name> <hemi> <parc name> <output annot>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
          "\nThis program samples a volumetric parcellation onto a surface\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  fprintf(stderr,
          "-t <thickness file>          - use specified file for computing "
          "thickness statistics\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

#if 0
static int
translate_indices_to_annotations(MRI_SURFACE *mris, char *translation_fname)
{
  char   fname[STRLEN], *cp, line[STRLEN], name[STRLEN], **names ;
  FILE   *fp ;
  int    nlines, i, vno, *r, *g, *b, *indices ;
  VERTEX *v ;

  cp = getenv("CSURF_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: CSURF_DIR not defined in environment", Progname) ;
  sprintf(fname, "%s/%s", cp, translation_fname) ;
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not read translation file %s", Progname, fname) ;

  nlines = 0 ;
  cp = fgetl(line, STRLEN-1, fp) ;
  while (cp != NULL)
  {
    nlines++ ;
    cp = fgetl(line, STRLEN-1, fp) ;
  }
  printf("%d lines found in file\n", nlines) ;
  rewind(fp) ;

  r = (int *)calloc(nlines, sizeof(int)) ;
  g = (int *)calloc(nlines, sizeof(int)) ;
  b = (int *)calloc(nlines, sizeof(int)) ;
  indices = (int *)calloc(nlines, sizeof(int)) ;
  names = (char **)calloc(nlines, sizeof(char *)) ;
  if (!r || !g || !b || !indices || !names)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate %d-len internal buffers",
              Progname, nlines) ;

  for (i = 0 ; i < nlines ; i++)
  {
    cp = fgetl(line, STRLEN-1, fp) ;
    sscanf(cp, "%d %s %d %d %d %*d",
           &indices[i], name, &r[i], &g[i], &b[i]) ;
    names[i] = (char *)calloc(strlen(name)+1, sizeof(char)) ;
    strcpy(names[i], name) ;
    printf("parsing parcellation unit %s...\n", name) ;
  }
  fclose(fp) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
  }
  return(NO_ERROR) ;
}
#else
static int
translate_indices_to_annotations(MRI_SURFACE *mris, char *translation_fname)
{
  int    vno ;
  VERTEX *v ;

  read_named_annotation_table(translation_fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->annotation = index_to_annotation((int)v->val) ;
  }
  return(NO_ERROR) ;
}
#endif
