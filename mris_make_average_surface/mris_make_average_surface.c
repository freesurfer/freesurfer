//
//  mris_make_average_surface.c
//
//
// Warning: Do not edit the following four lines.  CVS maintains them.
// Revision Author: $Author: fischl $
// Revision Date  : $Date: 2005/10/03 19:31:10 $
// Revision       : $Revision: 1.17 $
//
////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mri.h"
#include "mrisurf.h"
#include "macros.h"
#include "icosahedron.h"
#include "transform.h"
#include "version.h"
#include "fio.h"

static char vcid[] = "$Id: mris_make_average_surface.c,v 1.17 2005/10/03 19:31:10 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int normalize_area = 1 ;

static char *orig_name = "orig" ;
static char *xform_name = "talairach.xfm" ;

static int ico_no = 6 ;

char *Progname ;
static char sdir[STRLEN];

int
main(int argc, char *argv[])
{
  char         **av, *avg_surf_name, *canon_surf_name, fname[STRLEN], 
    *mdir, ico_fname[STRLEN], *hemi, *out_sname ;
  int          ac, nargs, i, vno, n ;
  VERTEX       *v ;
  MRI_SURFACE  *mris, *mris_ico ;
  MRI_SP       *mrisp, *mrisp_total ;
  LTA          *lta ;
  MRI          *mri ;
  VOL_GEOM      vg;
	float        average_surface_area = 0.0 ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mris_make_average_surface.c,v 1.17 2005/10/03 19:31:10 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  mdir = getenv("FREESURFER_HOME") ;
  if (!mdir)
    ErrorExit(ERROR_BADPARM, "%s: no FREESURFER_HOME in envoronment.\n",Progname);
  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }
  if (!strlen(sdir))
  {
    strcpy(sdir, getenv("SUBJECTS_DIR")) ;
    if (!sdir)
      ErrorExit(ERROR_BADPARM, "%s: no SUBJECTS_DIR in envoronment.\n",Progname);
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
    // read sphere.reg
    mris = MRISread(fname) ;
    if (!mris)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
								Progname, fname) ;
    // get "pial" surface vertex into ->origx, origy, origz
    if (MRISreadOriginalProperties(mris, orig_name) != NO_ERROR)
      ErrorExit(ERROR_BADFILE,"%s: could not read orig file for %s.\n",
                Progname, argv[1]);
    // read transform
    sprintf(fname, "%s/%s/mri/transforms/%s", sdir, argv[i], xform_name) ;
    lta = LTAreadEx(fname) ;
    if (!lta)
      ErrorExit(ERROR_BADPARM, "%s: could not read transform from %s", Progname, fname) ;

    // read T1 volume
    sprintf(fname, "%s/%s/mri/T1.mgz", sdir, argv[i]) ;
    if(fio_FileExistsReadable(fname)) mri = MRIreadHeader(fname,MRI_MGH_FILE);
    else{
      sprintf(fname, "%s/%s/mri/T1", sdir, argv[i]) ;
      mri = MRIreadHeader(fname, MRI_UCHAR); // MRI_CORONAL_SLICE_DIRECTORY) ;
    }

    if (!mri)
      ErrorExit(ERROR_BADPARM, "%s: could not read reference MRI volume from %s", Progname, fname) ;

    // save current vertex position into ->cx
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    // get the vertex position from ->origx, ... (get the "pial" vertex position)
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
		MRIScomputeMetricProperties(mris) ;
		printf("surface area: %2.1f cm^2\n", mris->total_area/100) ;
		average_surface_area += mris->total_area ;

    // this means that we transform "pial" surface
#if 1
    //MRIStalairachTransform(mris, mris, lta) ;
    {/*-----------------------------------------------------------------*/
      MATRIX *XFM, *sras, *tras;

      XFM = DevolveXFMWithSubjectsDir(argv[i], NULL, "talairach.xfm", sdir);
      if(XFM == NULL) exit(1);
      
      sras = MatrixAlloc(4,1,MATRIX_REAL);
      sras->rptr[4][1] = 1;
      tras = MatrixAlloc(4,1,MATRIX_REAL);
      tras->rptr[4][1] = 1;
      
			if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
				printf("Applying transform.\n");  
      for(vno=0; vno < mris->nvertices; vno++){
				v = &mris->vertices[vno] ;
				if (v->ripflag) continue ;
				sras->rptr[1][1] = mris->vertices[vno].x;
				sras->rptr[2][1] = mris->vertices[vno].y;
				sras->rptr[3][1] = mris->vertices[vno].z;
				tras = MatrixMultiply(XFM,sras,tras);
				mris->vertices[vno].x = tras->rptr[1][1];
				mris->vertices[vno].y = tras->rptr[2][1];
				mris->vertices[vno].z = tras->rptr[3][1];
				if (Gdiag_no == vno)
					printf(" v %d: (%2.1f, %2.1f, %2.1f) --> (%2.1f, %2.1f, %2.1f)\n",
								 vno, sras->rptr[1][1], sras->rptr[2][1], sras->rptr[3][1],
								 tras->rptr[1][1], tras->rptr[2][1], tras->rptr[3][1]) ;
      }
      //mrisComputeSurfaceDimensions(mris) ;

    }/*-----------------------------------------------------------------*/
#else
    MRIStransform(mris, mri, lta, 0) ;
    // copy volume geometry of the transformed volume
    memcpy((void *) &vg, (void *) &(mris->vg), sizeof(VOL_GEOM)); 
#endif
    // save transformed position in ->orig (store "pial" vertices position in orig)
		MRIScomputeMetricProperties(mris) ;
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    // get the vertex position from ->cx (note that this is not transformed)  sphere.reg vertices
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
    // mris contains sphere.reg in vertex and pial vertices in orig
    // map to a theta-phi space and accumulate values
    mrisp = MRIScoordsToParameterization(mris, NULL, SCALE) ;
    MRISPaccumulate(mrisp, mrisp_total, 0) ;
    MRISPaccumulate(mrisp, mrisp_total, 1) ;
    MRISPaccumulate(mrisp, mrisp_total, 2) ;
    MRISPfree(&mrisp) ; MRISfree(&mris) ; LTAfree(&lta) ; MRIfree(&mri) ;
    n++ ;
  }
	average_surface_area /= (float)n ;

  // mrisp_total lost info on the modified surface
  sprintf(ico_fname, "%s/lib/bem/ic%d.tri", mdir, ico_no) ;
  fprintf(stderr, "reading icosahedron from %s...\n", ico_fname) ;
  mris_ico = ICOread(ico_fname) ;
  if (!mris_ico)
    ErrorExit(ERROR_NOFILE, "%s: could not read icosahedron file %s\n",
              Progname,ico_fname) ;
  MRISscaleBrain(mris_ico, mris_ico, 
                 DEFAULT_RADIUS/MRISaverageRadius(mris_ico)) ;
  // save current ico position to ->cx, cy, cz
  MRISsaveVertexPositions(mris_ico, CANONICAL_VERTICES) ;
  // using mrisp_total to calculate position into ->origx, origy, origz (orig is the "pial" vertices)
  MRIScoordsFromParameterization(mrisp_total, mris_ico) ;
  // copy geometry info
  memcpy((void *) &mris_ico->vg, (void *) &vg, sizeof (VOL_GEOM));

  if (Gdiag_no >= 0 && Gdiag_no < mris_ico->nvertices)
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
  // write *h.sphere.reg
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, out_sname, hemi, canon_surf_name) ;
  if (Gdiag & DIAG_SHOW)
    fprintf(stderr,"writing average canonical surface to %s\n", fname);
  MRISwrite(mris_ico, fname) ;

  // get "pial vertices" from orig
  MRISrestoreVertexPositions(mris_ico, ORIG_VERTICES);
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    // n = number of subjects
    v->x /= (float)n ;
    v->y /= (float)n ;
    v->z /= (float)n ;
  }
	if (normalize_area)
	{
		MRIScomputeMetricProperties(mris_ico) ;
		printf("setting group surface area to be %2.1f cm^2 (scale=%2.2f)\n",
					 average_surface_area/100.0,
					 sqrt(average_surface_area/mris_ico->total_area)) ;

#if 0					 
		MRISscaleBrain(mris_ico, mris_ico, 
									 sqrt(average_surface_area/mris_ico->total_area)) ;
#else
		mris_ico->group_avg_surface_area = average_surface_area ;
#endif
		MRIScomputeMetricProperties(mris_ico) ;
	}
		
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, out_sname, hemi, avg_surf_name) ;
  printf("writing average %s surface to %s\n", avg_surf_name, fname);
  MRISwrite(mris_ico,  fname) ;
  {
    char path[STRLEN] ;
    LTA  *lta ;
    
    FileNamePath(fname, path) ;
    lta = LTAalloc(1, NULL) ;
    // write to a different location
    sprintf(fname, "%s/../mri/transforms/%s", path,xform_name) ;
    LTAwriteEx(lta, fname) ;
    LTAfree(&lta) ;
  }

  MRISfree(&mris_ico) ;
  MRISPfree(&mrisp_total) ;

  printf("mris_make_average_surface done\n");

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
  }
  else if (!stricmp(option, "nonorm"))
  {
		normalize_area = 0 ;
		printf("not normalizing surface area\n") ;
  }
  else switch (toupper(*option))
  {
  case 'I':
    ico_no = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'X':
    xform_name = argv[2] ;
    nargs = 1 ;
    printf("using xform %s...\n", xform_name) ;
    break ;
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'S':
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
         "usage: %s [options] <hemi> <output surf name> <canon surface>\n\t<output subject name> <subject> ... "
          "\n", Progname) ;
	printf("this program will generate an average of the orig surfaces of all the subjects\n"
				 "specified (unless the -s <surface name> flag is used)\n") ;
	printf("\tthe transform defaults to %s in the subject's mri/transforms directory, but can\n",
				 xform_name) ;
	printf("\tbe changed using the -x <xform name> switch\n") ;
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

