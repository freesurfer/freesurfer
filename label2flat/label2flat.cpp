/*
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
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "minc.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "mrisurf.h"
#include "label.h"
#include "version.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

static Transform *load_transform(char *subject_name, General_transform *xform);

static void print_usage(void) ;
void print_help(void) ;

const char *Progname ;

static int verbose = 0 ;

static int no_talairach = 0 ;  /* use talairach coords by default */
static int ndilate = 0 ;
static int nerode = 0 ;
static int nclose = 0 ;

#define MAX_AREAS     100
#define MAX_SUBJECTS  1000
#define NAME_LEN      100

static char subjects_dir[NAME_LEN] = "" ;

static General_transform    transform ;
static Transform            *linear_transform ;

static char *output_subject = NULL ;

static char *canon_name = NULL ;

#if 0
typedef struct {
  int    n_points ;           /* # of points in area */
  char   name[100] ;          /* name of label file */
  char   subject_name[100] ;  /* name of subject */
  int    *vno ;       /* indices into surface tessellation */
  float  *x ;
  float  *y ;
  float  *z ;
}
LABEL ;


int LabelFree(LABEL **parea) ;
int LabelDump(FILE *fp, LABEL *area) ;
LABEL *LabelRead(char *subject_name, char *label_name) ;
int LabelWrite(LABEL *area, char *fname) ;
int LabelToCanonical(LABEL *area, MRI_SURFACE *mris) ;
int LabelFromCanonical(LABEL *area, MRI_SURFACE *mris) ;
int LabelToFlat(LABEL *area, MRI_SURFACE *mris) ;
#endif

int
main(int argc, char *argv[]) {
  char         **av ;
  int          ac, nargs ;
  char         *cp, label_fname[STRLEN], *subject_name, *label_name,
  *out_fname, *patch_name, surf_fname[STRLEN], hemi[10] ;
  LABEL        *area ;
  MRI_SURFACE  *mris ;

  nargs = handleVersionOption(argc, argv, "label2flat");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  /* read in command-line options */
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 5)
    print_usage() ;

  subject_name = argv[1] ;
  label_name = argv[2] ;
  patch_name = argv[3] ;
  out_fname = argv[4] ;
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "no subjects directory in environment.\n") ;
  strcpy(subjects_dir, cp) ;

  int req = snprintf(label_fname, STRLEN,
		     "%s/%s/label/%s.label",
		     subjects_dir, subject_name, label_name) ;
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }

  linear_transform = load_transform(subject_name, &transform) ;

  cp = strrchr(patch_name, '.') ;
  if (!cp) {
    strcpy(hemi, "lh") ;
  } else {
    strncpy(hemi, cp-2, 2) ;
  }
  hemi[2] = 0 ;
  req = snprintf(surf_fname, STRLEN,
		 "%s/%s/surf/%s.orig", subjects_dir, subject_name, hemi);
  if( req >= STRLEN ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  fprintf(stderr, "reading surface %s...\n", surf_fname) ;
  mris = MRISread(surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_BADFILE, "%s: could not read surface file %s.",
              Progname, surf_fname) ;
  if (MRISreadPatch(mris, patch_name) != NO_ERROR)
    exit(Gerror) ;

  area = LabelRead(subject_name, label_name) ;
  if (canon_name)   /* put it onto a canonical surface */
  {
    cp = strrchr(canon_name, '.') ;
    if (cp)  { /* hemisphere specified explicitly */
      int req = snprintf(surf_fname, STRLEN,
			 "%s/%s/surf/%s", 
			 subjects_dir, subject_name,
			 canon_name) ;
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    } else {
      int req = snprintf(surf_fname, STRLEN,
			 "%s/%s/surf/%s.%s", 
			 subjects_dir, subject_name,
			 hemi, canon_name) ;
      if( req >= STRLEN ) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
    }
    MRISreadCanonicalCoordinates(mris, surf_fname) ;
    LabelToCanonical(area, mris) ;
  }
  if (ndilate > 0)
    LabelDilate(area, mris, ndilate, CURRENT_VERTICES) ;
  if (nerode > 0)
    LabelErode(area, mris, nerode) ;
  if (output_subject)   /* write onto a different subject's flat map */
  {
    MRISfree(&mris) ;
    int req = snprintf(surf_fname, STRLEN,
		       "%s/%s/surf/%s.orig", subjects_dir, subject_name,hemi);
    if( req >= STRLEN ) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    fprintf(stderr, "reading surface %s...\n", surf_fname) ;
    mris = MRISread(surf_fname) ;
    if (!mris)
      ErrorExit(ERROR_BADFILE, "%s: could not read surface file %s.",
                Progname, surf_fname) ;
    if (MRISreadPatch(mris, patch_name) != NO_ERROR)
      exit(Gerror) ;
    if (canon_name)   /* put it onto a canonical surface */
    {
      cp = strrchr(canon_name, '.') ;
      if (cp)  { /* hemisphere specified explicitly */
        int req = snprintf(surf_fname, STRLEN,
			   "%s/%s/surf/%s",
			   subjects_dir, subject_name,
			   canon_name) ;
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      } else {
        int req = snprintf(surf_fname, STRLEN,
			   "%s/%s/surf/%s.%s",
			   subjects_dir, subject_name,
			   hemi, canon_name) ;
	if( req >= STRLEN ) {
	  std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	}
      }
      MRISreadCanonicalCoordinates(mris, surf_fname) ;
    } else {
      MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    }

    LabelFromCanonical(area, mris) ;
  }
  LabelToFlat(area, mris) ;
  /*  LabelDump(stdout, area) ;*/
  LabelWrite(area, out_fname) ;
  LabelFree(&area) ;
  MRISfree(&mris) ;
  if (verbose)
    fprintf(stderr, "done.\n") ;
  exit(0) ;
  return(0) ;
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
  if (!stricmp(option, "dilate"))
  {
    ndilate = atoi(argv[2]) ;
    printf("dilating label %d times before generating patch\n", ndilate) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "erode"))
  {
    nerode = atoi(argv[2]) ;
    printf("eroding label %d times before generating patch\n", nerode) ;
    nargs = 1 ;
  }
  else if (!stricmp(option, "close"))
  {
    nclose = atoi(argv[2]) ;
    printf("closing label %d times before generating patch\n", nclose) ;
    ndilate = nclose ;
    nerode = nclose ;
    nargs = 1 ;
  }
  else switch (toupper(*option)) {
  case 'C':
    canon_name = argv[2] ;
    fprintf(stderr, "using surface %s as canonical coordinate system.\n",
            canon_name) ;
    nargs = 1 ;
    no_talairach = 1 ;
    break ;
  case 'V':
    verbose = !verbose ;
    break ;
  case 'N':
    no_talairach = 1 ;
    break ;
  case 'O':
    output_subject = argv[2] ;
    fprintf(stderr, "generating map on subject %s\n", output_subject) ;
    nargs = 1 ;
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

static Transform *
load_transform(char *subject_name, General_transform *transform) {
  char xform_fname[100] ;

  int req = snprintf(xform_fname, 100,
		     "%s/%s/mri/transforms/talairach.xfm",
		     subjects_dir, subject_name) ;
  if( req >= 100 ) {
    std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
  }
  if (input_transform_file(xform_fname, transform) != OK)
    ErrorExit(ERROR_NOFILE, "%s: could not load transform file '%s'",
              Progname, xform_fname) ;

  if (verbose == 2)
    fprintf(stderr, "transform read successfully from %s\n", xform_fname) ;
  return(get_linear_transform_ptr(transform)) ;
}

static void
print_usage(void) {
  printf("usage: %s [options] <subject name> <label file> <patch file> "
         "<output file>\n",
         Progname) ;
  exit(1) ;
}


void
print_help(void) {
  printf("usage: %s [options] <subject name> <label file> <patch file>"
         " <output file>\n", Progname) ;
  exit(1) ;
}

#if 0
LABEL *
LabelRead(char *subject_name, char *label_name) {
  LABEL  *area ;
  char   fname[100], *cp, line[200] ;
  FILE   *fp ;
  int    vno, nlines ;
  float  x, y, z ;

  area = (LABEL *)calloc(1, sizeof(LABEL)) ;
  if (!area)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate LABEL struct.",Progname);
  sprintf(fname, "%s/%s/label/%s.label", subjects_dir,subject_name,label_name);

  strcpy(area->name, label_name) ;
  strcpy(area->subject_name, subject_name) ;

  /* read in the file */
  fp = fopen(fname, "r") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open label file %s",
              Progname, fname) ;

  cp = fgetl(line, 199, fp) ;
  if (!cp)
    ErrorExit(ERROR_BADFILE, "%s: empty label file %s", Progname, fname) ;
  if (!sscanf(cp, "%d", &area->n_points))
    ErrorExit(ERROR_BADFILE, "%s: could not scan # of lines from %s",
              Progname, fname) ;
  area->x = (float *)calloc(area->n_points, sizeof(float)) ;
  if (!area->x)
    ErrorExit(ERROR_NOMEMORY,
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;
  area->y = (float *)calloc(area->n_points, sizeof(float)) ;
  if (!area->y)
    ErrorExit(ERROR_NOMEMORY,
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;
  area->z = (float *)calloc(area->n_points, sizeof(float)) ;
  if (!area->z)
    ErrorExit(ERROR_NOMEMORY,
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;
  area->vno = (int *)calloc(area->n_points, sizeof(int)) ;
  if (!area->vno)
    ErrorExit(ERROR_NOMEMORY,
              "%s: LabelRead(%s) could not allocate %d-sized vector",
              Progname, label_name, area->n_points) ;

  nlines = 0 ;
  while ((cp = fgetl(line, 199, fp)) != NULL) {
    if (sscanf(cp, "%d %f %f %f", &vno, &x, &y, &z) != 4)
      ErrorExit(ERROR_BADFILE, "%s: could not parse %dth line in %s",
                Progname, area->n_points, fname) ;
    area->x[nlines] = x ;
    area->y[nlines] = y ;
    area->z[nlines] = z ;
    area->vno[nlines] = vno ;
    nlines++ ;
  }

  fclose(fp) ;
  if (!nlines)
    ErrorExit(ERROR_BADFILE, "%s: no data in label file %s", Progname, fname);
  return(area) ;
}
int
LabelDump(FILE *fp, LABEL *area) {
  int  n ;

  fprintf(fp, "label %s, from subject %s\n", area->name, area->subject_name) ;
  for (n = 0 ; n < area->n_points ; n++)
    fprintf(fp, "%d  %2.3f  %2.3f  %2.3f\n", area->vno[n], area->x[n],
            area->y[n], area->z[n]) ;
  return(NO_ERROR) ;
}
int
LabelFree(LABEL **parea) {
  LABEL *area ;

  area = *parea ;
  *parea = NULL ;
  free(area->x) ;
  free(area->y) ;
  free(area->z) ;
  free(area->vno) ;
  free(area) ;
  return(NO_ERROR) ;
}

int
LabelToCanonical(LABEL *area, MRI_SURFACE *mris) {
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++) {
    vno = area->vno[n] ;
    v = &mris->vertices[vno] ;
    area->vno[n] = -1 ;      /* not associated with a vertex anymore */
    area->x[n] = v->cx ;
    area->y[n] = v->cy ;
    area->z[n] = v->cz ;
  }
  return(NO_ERROR) ;
}
int
LabelFromCanonical(LABEL *area, MRI_SURFACE *mris) {
  int     n, vno ;
  VERTEX  *v ;

  for (n = 0 ; n < area->n_points ; n++) {
    vno =MRISfindClosestCanonicalVertex(mris,area->x[n],area->y[n],area->z[n]);
    v = &mris->vertices[vno] ;
    area->vno[n] = vno ;
    area->x[n] = v->cx ;
    area->y[n] = v->cy ;
    area->z[n] = v->cz ;
  }
  return(NO_ERROR) ;
}
int
LabelToFlat(LABEL *area, MRI_SURFACE *mris) {
  int     n, vno ;
  VERTEX  *v ;
  float dmin;

  for (n = 0 ; n < area->n_points ; n++) {
    vno = area->vno[n] ;
    if (vno >= 0)   /* already have associated vertex */
    {
      v = &mris->vertices[vno] ;
      area->x[n] = v->x ;
      area->y[n] = v->y ;
      area->z[n] = v->z ;
    } else    /* in canonical coordinate system - find closest vertex */
    {
      vno = MRISfindClosestVertex(mris, area->x[n], area->y[n], area->z[n],&dmin) ;
      v = &mris->vertices[vno] ;
      area->vno[n] = vno ;
      area->x[n] = v->x ;
      area->y[n] = v->y ;
      area->z[n] = v->z ;
    }
  }
  return(NO_ERROR) ;
}
int
LabelWrite(LABEL *area, char *fname) {
  FILE  *fp ;
  int  n ;

  fp = fopen(fname, "w") ;
  if (!fp)
    ErrorExit(ERROR_NOFILE, "%s: could not open label file %s",
              Progname, fname) ;

#if 0
  fprintf(fp, "# label %s, from subject %s\n", area->name, area->subject_name);
#endif
  for (n = 0 ; n < area->n_points ; n++)
    fprintf(fp, "%d  %2.3f  %2.3f  %2.3f\n", area->vno[n], area->x[n],
            area->y[n], area->z[n]) ;
  fclose(fp) ;
  return(NO_ERROR) ;
}

#endif
