
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
#include "icosahedron.h"

static char vcid[] = "$Id: mris_make_average_surface.c,v 1.2 2002/01/30 19:36:42 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static char *orig_name = "orig" ;

static int ico_no = 6 ;

char *Progname ;

int
main(int argc, char *argv[])
{
  char         **av, *avg_surf_name, *canon_surf_name, fname[STRLEN], *sdir, 
               *mdir, ico_fname[STRLEN], *hemi, *out_sname ;
  int          ac, nargs, i, vno, n ;
  VERTEX       *v ;
  MRI_SURFACE  *mris, *mris_ico ;
  MRI_SP       *mrisp, *mrisp_total ;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  mdir = getenv("MRI_DIR") ;
  if (!mdir)
    ErrorExit(ERROR_BADPARM, "%s: no MRI_DIR in envoronment.\n",Progname);
  sdir = getenv("SUBJECTS_DIR") ;
  if (!sdir)
    ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in envoronment.\n",Progname);
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 6)
    usage_exit() ;

  hemi = argv[1] ;
  avg_surf_name = argv[2] ;
  canon_surf_name = argv[3] ;
  out_sname = argv[4] ;

#define SCALE 1
  mrisp_total = MRISPalloc(SCALE, 3) ;
  for (n = 0, i = 5 ; i < argc ; i++)
  {
    fprintf(stderr, "processing subject %s...\n", argv[i]) ;
    sprintf(fname, "%s/%s/surf/%s.%s", sdir, argv[i], hemi, canon_surf_name) ;
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, fname) ;
    if (MRISreadOriginalProperties(mris, orig_name) != NO_ERROR)
      ErrorExit(ERROR_BADFILE,"%s: could not read orig file for %s.\n",
                Progname, argv[1]);
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIStalairachTransform(mris, mris) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    mrisp = MRIScoordsToParameterization(mris, NULL, SCALE) ;
    MRISPaccumulate(mrisp, mrisp_total, 0) ;
    MRISPaccumulate(mrisp, mrisp_total, 1) ;
    MRISPaccumulate(mrisp, mrisp_total, 2) ;
    MRISPfree(&mrisp) ;
    MRISfree(&mris) ;
    n++ ;
  }

  sprintf(ico_fname, "%s/lib/bem/ic%d.tri", mdir, ico_no) ;
  fprintf(stderr, "reading icosahedron from %s...\n", ico_fname) ;
  mris_ico = ICOread(ico_fname) ;
  if (!mris_ico)
    ErrorExit(ERROR_NOFILE, "%s: could not read icosahedron file %s\n",
              Progname,ico_fname) ;
  MRISscaleBrain(mris_ico, mris_ico, 
                 DEFAULT_RADIUS/MRISaverageRadius(mris_ico)) ;
  MRISsaveVertexPositions(mris_ico, CANONICAL_VERTICES) ;
  MRIScoordsFromParameterization(mrisp_total, mris_ico) ;

  if (Gdiag_no >= 0)
  {
    int n ;
    VERTEX *vn ;

    v = &mris_ico->vertices[Gdiag_no] ;
    fprintf(stderr, "v %d: x = (%2.2f, %2.2f, %2.2f)\n",
            Gdiag_no, v->origx, v->origy, v->origz) ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris_ico->vertices[v->v[n]] ;
      fprintf(stderr, "v %d: x = (%2.2f, %2.2f, %2.2f)\n",
              v->v[n], vn->origx, vn->origy, vn->origz) ;
    }
  }
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, out_sname, hemi, canon_surf_name) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,"writing average canonical surface to to %s\n", fname);
  MRISwrite(mris_ico, fname) ;

  MRISrestoreVertexPositions(mris_ico, ORIG_VERTICES);
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    v->x /= (float)n ;
    v->y /= (float)n ;
    v->z /= (float)n ;
  }
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, out_sname, hemi, avg_surf_name) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,"writing average orig surface to to %s\n", fname);
  MRISwrite(mris_ico,  fname) ;

  MRISfree(&mris_ico) ;
  MRISPfree(&mrisp_total) ;
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
  case 'I':
    ico_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'O':
    orig_name = argv[2] ;
    printf("reading vertex positions from %s...\n", orig_name) ;
    nargs = 1 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
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

static void
print_usage(void)
{
  printf(
         "usage: %s [options] <hemi> <surf name> <canon surface>\n\t<output subject name> <subject> ... "
          " <output curv file >\n", Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program will average a set of surface coordinates and genareate an average\nsurface (using Talairach coords and spherical transform).\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}

