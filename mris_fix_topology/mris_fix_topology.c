
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
#include "icosohedron.h"
#include "mrishash.h"

static char vcid[] = "$Id: mris_fix_topology.c,v 1.1 1998/11/16 20:30:57 fischl Exp $";

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

int MRISmergeIcosahedrons(MRI_SURFACE *mri_src, MRI_SURFACE *mri_dst) ;
int MRISinverseSphericalMap(MRI_SURFACE *mris, MRI_SURFACE *mris_ico) ;

char *Progname ;

static int small = 0 ;
static int huge = 1 ;

static int avgs = 10 ;

#define MAX_VERTICES  0
#define MAX_FACES     0

int
main(int argc, char *argv[])
{
  char          **av, *hemi, *sname, sdir[400], *cp, fname[500] ;
  int           ac, nargs ;
  MRI_SURFACE   *mris, *mris_ico ;
  int           msec ;
  struct timeb  then ;
  int           vno ;
  float         radius ;

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
  sname = argv[1] ;
  hemi = argv[2] ;
  cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, 
              "%s: SUBJECTS_DIR not defined in environment.\n", Progname) ;
  strcpy(sdir, cp) ;

  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, SPHERE_NAME) ;
  fprintf(stderr, "reading input surface %s...\n", fname) ;
  mris = MRISread(fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read input surface %s",
              Progname, fname) ;
  MRIScomputeMetricProperties(mris) ;
  MRISreadOriginalProperties(mris, "orig") ;
/*
  MRISreadCanonicalCoordinates(mris, "inflated") ;
*/
/*
  MRISreadCanonicalCoordinates(mris, "sphere") ;
*/
  fprintf(stderr, "creating icosahedral representation...\n") ;

  if (huge)
    mris_ico = ic163842_make_surface(MAX_VERTICES, MAX_FACES) ;
  else
  {
    if (small)
      mris_ico = ic10242_make_surface(MAX_VERTICES, MAX_FACES) ;
    else
      mris_ico = ic40962_make_surface(MAX_VERTICES, MAX_FACES) ;
  }
  if (!mris_ico)
    ErrorExit(ERROR_NOFILE, "%s: could not create/read MR surface.",Progname) ;

  /* save orig. vertex coords */
  radius = MRISaverageRadius(mris_ico) ;
  MRISscaleBrain(mris_ico, mris_ico, 100.0f/radius) ;
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    mris_ico->vertices[vno].origx = mris_ico->vertices[vno].x;
    mris_ico->vertices[vno].origy = mris_ico->vertices[vno].y;
    mris_ico->vertices[vno].origz = mris_ico->vertices[vno].z;
  }

  fprintf(stderr, "mapping icosahedron to cortical surface...\n") ;

  MRISreadCanonicalCoordinates(mris, "inflated") ;
  MRISinverseSphericalMap(mris, mris_ico) ;
  MRISaverageVertexPositions(mris, avgs) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "inflated_ico") ;
  MRISwrite(mris_ico, fname) ;

  MRISreadCanonicalCoordinates(mris, "orig") ;
  MRISinverseSphericalMap(mris, mris_ico) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "orig_ico") ;
  MRISwrite(mris_ico, fname) ;

  MRISreadCanonicalCoordinates(mris, "smoothwm") ;
  MRISinverseSphericalMap(mris, mris_ico) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "smoothwm_ico") ;
  MRISwrite(mris_ico, fname) ;

  MRISreadCanonicalCoordinates(mris, "white") ;
  MRISinverseSphericalMap(mris, mris_ico) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "white_ico") ;
  MRISwrite(mris_ico, fname) ;

  MRISreadCanonicalCoordinates(mris, "sphere") ;
  MRISinverseSphericalMap(mris, mris_ico) ;
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "sphere_ico") ;
  MRISwrite(mris_ico, fname) ;

  MRISwriteCurvature(mris, "neg") ;
  sprintf(fname, "%s/%s/surf/%s-%s", sdir, sname, hemi, "ico.neg") ;
  MRISwriteCurvature(mris_ico, fname) ;

/*
  sprintf(fname, "%s/%s/surf/%s.%s", sdir, sname, hemi, "ico_geo") ;
  fprintf(stderr, "writing output surface to %s...\n", fname) ;
  MRISwrite(mris_ico, fname) ;
*/

  msec = TimerStop(&then) ;
  fprintf(stderr,"topology fixing took %2.1f minutes\n", 
          (float)msec/(60*1000.0f));
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
  case '?':
  case 'U':
    print_usage() ;
    exit(1) ;
    break ;
  case 'A':
    avgs = atoi(argv[2]) ;
    nargs = 1 ;
    break ;
  case 'S':
    small = 1 ;
    break ;
  case 'H':
    huge = 1 ;
    break ;
  case 'B':
    small = 0 ;
    break ;
  case 'V':
    Gdiag_no = atoi(argv[2]) ;
    nargs = 1 ;
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
  fprintf(stderr, "usage: %s [options] <subject name> <hemisphere>\n", 
          Progname) ;
}

static void
print_help(void)
{
  print_usage() ;
  fprintf(stderr, 
       "\nThis program positions the tessellation of the cortical surface\n"
          "at the white matter surface, then the gray matter surface\n"
          "and generate surface files for these surfaces as well as a\n"
          "'curvature' file for the cortical thickness, and a surface file\n"
          "which approximates layer IV of the cortical sheet.\n");
  fprintf(stderr, "\nvalid options are:\n\n") ;
  /*  fprintf(stderr, "-n    normalize output curvatures.\n") ;*/
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}



static int mrisCalculateFaceCentroid(MRI_SURFACE *mris, int fno, 
                                     float *px, float *py, float *pz) ;
static int mrisFindCommonFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2);
static int mrisCalculateOriginalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                     float *px, float *py, float *pz) ;
static int mrisCalculateCanonicalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                     float *px, float *py, float *pz) ;


int
MRISmergeIcosahedrons(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int    vno, n, unmarked = 0, fno, vlist[100], num ;
  VERTEX *vsrc, *vdst, *vn ;
  float  x, y, z, dx, dy, dz, val, ox, oy, oz ;

  /* first copy all vertices to the destination that exist in the source */
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vsrc = &mris_src->vertices[vno] ;
    vdst = &mris_dst->vertices[vno] ;
    vdst->marked = 1 ;
    vdst->x = vsrc->x ; vdst->y = vsrc->y ; vdst->z = vsrc->z ;
    vdst->origx = vsrc->origx ; vdst->origy = vsrc->origy ; 
    vdst->origz = vsrc->origz ;
    vdst->val = vsrc->val ;
    vdst->odx = vsrc->odx ; vdst->ody = vsrc->ody ; vdst->odz = vsrc->odz ;
    vdst->origarea = vsrc->origarea ;
  }

  /* now interpolate the newer vertices */
  for ( ; vno < mris_dst->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    num = 0 ; val = x = y = z = dx = dy = dz = 0.0f ;
    for (n = 0 ; n < vdst->vnum ; n++)
    {
      vn = &mris_dst->vertices[vdst->v[n]] ;
      if (!vn->marked)   /* a neighbor that hasn't been marked */
        continue ;

      dx += vn->odx ; dy += vn->ody ; dz += vn->odz ; 
      val += vn->val ;
      vlist[num++] = vdst->v[n] ;
      x += vn->x ; y += vn->y ; z += vn->z ; 
      ox += vn->origx ; oy += vn->origy ; oz += vn->origz ; 
    }
    /*    vdst->marked = 1 ;*/
    vdst->odx = dx / num ; vdst->ody = dy / num ; vdst->odz = dz / num ;
    vdst->val = val / num ;
    switch (num)
    {
    case 2:    /* put it at midpoint */
      vdst->x = x / num ; vdst->y = y / num ; vdst->z = z / num ; 
      break ;
    case 3:   /* put it at face center */
      fno = mrisFindCommonFace(mris_src, vlist[0], vlist[1], vlist[2]) ;
      if (fno >= 0)
      {
        mrisCalculateFaceCentroid(mris_src, fno, &x, &y, &z) ;
        mrisCalculateOriginalFaceCentroid(mris_src, fno, &ox, &oy, &oz) ;
        vdst->x = x ; vdst->y = y ; vdst->z = z ; 
        vdst->origx = ox ; vdst->origy = oy ; vdst->origz = oz ; 
      }
      else
        fprintf(stderr, 
                "vertices (%d, %d, %d) do not make up face in icomerge!\n",
                vlist[0], vlist[1], vlist[2]) ;
      break ;
    default:
      unmarked++ ;
      break ;
    }
  }

  mris_dst->orig_area = mris_src->orig_area ;
  if (unmarked)
    fprintf(stderr, "%d vertices left unmarked in merge!!!!\n", unmarked) ;

  MRIScomputeMetricProperties(mris_dst) ;
  return(NO_ERROR) ;
}
static int
mrisFindCommonFace(MRI_SURFACE *mris, int vno0, int vno1, int vno2)
{
  VERTEX  *v0, *v1, *v2 ;
  int     fno, n, fv ;
  FACE    *f ;

  v0 = &mris->vertices[vno0] ;
  v1 = &mris->vertices[vno1] ;
  v2 = &mris->vertices[vno2] ;

  for (n = 0 ; n < v0->num ; n++)
  {
    /* test to see if this face has all three vertices in it */
    fno = v0->f[n] ; f = &mris->faces[fno] ;
    for (fv = 0 ; fv < VERTICES_PER_FACE ; fv++)
      if (f->v[fv] != vno0 && f->v[fv] != vno1 && f->v[fv] != vno2)
        break ;
    if (fv >= VERTICES_PER_FACE)   /* they all matched one of the vertices */
      return(fno) ;
  }

  return(ERROR_BADPARM) ;
}



static int mrisChooseFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)  ;
static int mrisFindUnambiguousFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v) ;
int
MRISinverseSphericalMap(MRI_SURFACE *mris, MRI_SURFACE *mris_ico)
{
  int    vno, fno, n, num ;
  VERTEX *v, *vn ;
  MHT    *mht ;
  double r ;
  float  x, y, z ;

  float orient ;
  static char ***flagvol, ***numvol ; /* volume of flags for neg area test */
        static float ***xvol, ***yvol, ***zvol ;
        static int allocated = 0;
  float flagvolres = 2.0, flagvolfov = 300;
  int flagvoldim,i,j;
        int ix, iy, iz ;

  MRISclearCurvature(mris) ;   /* curvature will be used to calculate sulc */
  r = MRISaverageRadius(mris) ; MRISscaleBrain(mris, mris, 100.0f/r) ;
/*
  r = MRISaverageRadius(mris_ico) ; MRISscaleBrain(mris_ico,mris_ico,100.0f/r);
*/

  /*
    orig       positions are on orig
    cx,cy,cz   positions are on inflated surface
    current    positions are on sphere.
  */
/*
  printf("begin MHTfillTableWithFaces\n");

  mht = MHTfillTableWithFaces(mris, NULL) ;

  printf("end MHTfillTableWithFaces\n");
*/

  flagvoldim = ceil(flagvolfov/flagvolres);
  if (!allocated)
  {
    fprintf(stderr, "allocating flagvol...\n") ;
    flagvol = (char ***)calloc(flagvoldim, sizeof(char **));
    numvol = (char ***)calloc(flagvoldim, sizeof(char **));
    xvol = (float ***)calloc(flagvoldim, sizeof(float **));
    yvol = (float ***)calloc(flagvoldim, sizeof(float **));
    zvol = (float ***)calloc(flagvoldim, sizeof(float **));
    for (i=0;i<flagvoldim;i++)
    {
      flagvol[i] = (char **)calloc(flagvoldim, sizeof(char *));
      numvol[i] = (char **)calloc(flagvoldim, sizeof(char *));
      xvol[i] = (float **)calloc(flagvoldim, sizeof(float *));
      yvol[i] = (float **)calloc(flagvoldim, sizeof(float *));
      zvol[i] = (float **)calloc(flagvoldim, sizeof(float *));
      for (j=0;j<flagvoldim;j++)
      {
        flagvol[i][j] = (char *)calloc(flagvoldim, sizeof(char));
        numvol[i][j] = (char *)calloc(flagvoldim, sizeof(char));
        xvol[i][j] = (float *)calloc(flagvoldim, sizeof(float));
        yvol[i][j] = (float *)calloc(flagvoldim, sizeof(float));
        zvol[i][j] = (float *)calloc(flagvoldim, sizeof(float));
      }
    }
    allocated = 1;
  }
  
  for (ix=0;ix<flagvoldim;ix++)
    for (iy=0;iy<flagvoldim;iy++)
      for (iz=0;iz<flagvoldim;iz++)
        numvol[ix][iy][iz] = xvol[ix][iy][iz] = 
          yvol[ix][iy][iz] = zvol[ix][iy][iz] = flagvol[ix][iy][iz] = 0;
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    orient = mris->vertices[vno].x*mris->vertices[vno].nx+
      mris->vertices[vno].y*mris->vertices[vno].ny+
      mris->vertices[vno].z*mris->vertices[vno].nz ;
    if (orient < 0)
    {
      printf("vertex %d inside out (orient = %f)\n",vno,orient);
      /*
        printf("x = (%3.1f, %3.1f,%3.1f) n = (%3.1f,%3.1f,%3.1f)\n",
        mris->vertices[vno].x,mris->vertices[vno].y,mris->vertices[vno].z,
        mris->vertices[vno].nx,mris->vertices[vno].ny,mris->vertices[vno].nz);
        */
      mris->vertices[vno].curv = 1;
    }
    else
      mris->vertices[vno].curv = 0;
  }
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    ix = floor(0.5+(mris->vertices[vno].x+flagvolfov/2)/flagvolres);
    iy = floor(0.5+(mris->vertices[vno].y+flagvolfov/2)/flagvolres);
    iz = floor(0.5+(mris->vertices[vno].z+flagvolfov/2)/flagvolres);
    numvol[ix][iy][iz]++;
    xvol[ix][iy][iz] += mris->vertices[vno].cx;
    yvol[ix][iy][iz] += mris->vertices[vno].cy;
    zvol[ix][iy][iz] += mris->vertices[vno].cz;
    if (mris->vertices[vno].curv != 0) /* inverted vertex? */
      flagvol[ix][iy][iz] = mris->vertices[vno].curv;
  }
  
  for (ix=0;ix<flagvoldim;ix++)
    for (iy=0;iy<flagvoldim;iy++)
      for (iz=0;iz<flagvoldim;iz++)
        if (numvol[ix][iy][iz]>0)
        {
          xvol[ix][iy][iz] /= numvol[ix][iy][iz];
          yvol[ix][iy][iz] /= numvol[ix][iy][iz];
          zvol[ix][iy][iz] /= numvol[ix][iy][iz];
        }
  
  /* restore orig. vertex coords */
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    mris_ico->vertices[vno].x = mris_ico->vertices[vno].origx;
    mris_ico->vertices[vno].y = mris_ico->vertices[vno].origy;
    mris_ico->vertices[vno].z = mris_ico->vertices[vno].origz;
  }
  
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    ix = floor(0.5+(mris_ico->vertices[vno].x+flagvolfov/2)/flagvolres);
    iy = floor(0.5+(mris_ico->vertices[vno].y+flagvolfov/2)/flagvolres);
    iz = floor(0.5+(mris_ico->vertices[vno].z+flagvolfov/2)/flagvolres);
    mris_ico->vertices[vno].curv = flagvol[ix][iy][iz];
    if (mris_ico->vertices[vno].curv != 0)
    {
      mris_ico->vertices[vno].marked = 0;
      printf("ambiguous ico vertex %d\n",vno);
    }
    else
      mris_ico->vertices[vno].marked = 1;
    if (numvol[ix][iy][iz]>0)
    {
      mris_ico->vertices[vno].x = xvol[ix][iy][iz];
      mris_ico->vertices[vno].y = yvol[ix][iy][iz];
      mris_ico->vertices[vno].z = zvol[ix][iy][iz];
    }
    else
    {
      printf("### ico vertex %d missed volume\n",vno);
      mris_ico->vertices[vno].marked = 0;
    }
  }
  
  
  MRISsoapBubbleVertexPositions(mris_ico, 100) ;

/*
  for (vno = 0 ; vno < mris_ico->nvertices ; vno++)
  {
    v = &mris_ico->vertices[vno] ;
    if (v->marked || v->ripflag)
      continue ;
    fno = mrisChooseFace(mris, mht, v) ;
    mrisCalculateCanonicalFaceCentroid(mris, fno, &v->x, &v->y, &v->z) ;
  }
  MHTfree(&mht) ;
*/
  return(NO_ERROR) ;
}
static float mrisComputeFaceStretch(MRI_SURFACE *mris, int fno) ;
static int mrisCalculateCanonicalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                     float *px, float *py, float *pz) ;

#define CROSS(dest,v1,v2)                      \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define VERTEX_DIF(leg, v0, v1)   leg[0] = v1->x-v0->x, \
                                  leg[1] = v1->y-v0->y,\
                                  leg[2] = v1->z-v0->z ;

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])
/*
  mris is the original (hi-res) surface (from smoothwm).
  v is a vertex on the icosahedron.
  */
static int
mrisChooseFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int    i, fno, flist[10000], nfaces, min_fno = 0, l ;
  MHBT   *bucket, bucket_bak  ;
  MHB    *bin ;
  FACE   *f ;
  VERTEX *v1, *v2, *v3 ;
  float  ut, vt, stretch, n[3], leg[3], p[3], d[3], dot, ldot, leg2[3] ;

  bucket = MHTgetBucket(mht, v->x, v->y, v->z) ;
  bucket_bak = *bucket ;
  bin = bucket->bins ;
  for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)  /* check each face */
  {
    n[0] = n[1] = n[2] = 0.0f ;
    f = &mris->faces[bin->fno] ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      n[0] += v1->nx ; n[1] += v1->ny ; n[2] += v1->nz ; 
    }
    n[0] /= VERTICES_PER_FACE ; n[1] /= VERTICES_PER_FACE ; 
    n[2] /= VERTICES_PER_FACE ; 
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      switch (l)
      {
      default:
      case 0:
        v2 = &mris->vertices[f->v[1]] ;
        v3 = &mris->vertices[f->v[2]] ;
        break ;
      case 1:
        v2 = &mris->vertices[f->v[2]] ;
        v3 = &mris->vertices[f->v[0]] ;
        break ;
      case 2:
        v2 = &mris->vertices[f->v[0]] ;
        v3 = &mris->vertices[f->v[1]] ;
        break ;
      }

      VERTEX_DIF(leg, v1, v2) ;   /* leg of triangle */
      VERTEX_DIF(leg2, v1, v3) ;  /* leg of triangle */

      /* express p as point in triangle plane */
      VERTEX_DIF(p, v1, v) ;     /* vector from vertex to point in question */
      dot = DOT(p,n) ;
      p[0] -= dot*n[0] ; p[1] -= dot*n[1] ; p[2] -= dot*n[2] ; 
#if 0
      p[0] = ut*leg[0] + vt*leg2[0] ;
      p[1] = ut*leg[1] + vt*leg2[1] ;
      p[2] = ut*leg[2] + vt*leg2[2] ;
      ut = DOT(p, leg) ; vt = DOT(p, leg2) ;
#endif

      CROSS(d, leg, n) ;
      dot = DOT(d, p) ; ldot = DOT(d, leg2) ;

      /* on different side of leg from 3rd vertex */
      if (!FZERO(ldot) && !FZERO(dot) && dot*ldot < 0)
        break ;
    }
    if (l >= VERTICES_PER_FACE)
      flist[nfaces++] = bin->fno ;
  }
  if (!nfaces)  /* something went wrong, but Anders will fix it */
  {
    float dist, min_dist ;

    fprintf(stderr, "no faces found on sphere!\n") ;
    bin = bucket->bins ;
    min_dist = 1000000.0f ; min_fno = 0 ;
    for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      f = &mris->faces[bin->fno] ;
      v1 = &mris->vertices[f->v[0]] ;
      v2 = &mris->vertices[f->v[1]] ;
      v3 = &mris->vertices[f->v[2]] ;
#define VDIST(v1,v2) (sqrt(SQR(v1->x-v2->x)+SQR(v1->y-v2->y)+SQR(v1->z-v2->z)))
      dist = VDIST(v1,v) + VDIST(v2, v) + VDIST(v3,v) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = bin->fno ;
      }
    }
  }
  else   /* pix the face that is closest to the soap bubble location */
  {
    float min_dist, x, y, z, dist, curv ;
    int   vno ;

    if (nfaces > 1)
    {
      DiagBreak() ;
      curv = 1.0f ;
    }
    else
      curv = 0.0f ;
    for ( i = 0 ; i < nfaces ; i++)/* check each face */
    {
      f = &mris->faces[flist[i]] ;
      for (l = 0 ; l < VERTICES_PER_FACE ; l++)
        mris->vertices[f->v[l]].curv = curv ;
    }

    min_dist = 100000.0f ;
    for (i = 0 ; i < nfaces ; i++)
    {
      mrisCalculateCanonicalFaceCentroid(mris, flist[i], &x, &y, &z) ;
      dist = sqrt(SQR(v->x-x)+SQR(v->y-y)+SQR(v->z-z)) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = flist[i] ;
      }
    }
  }
  return(min_fno) ;
}

static int
mrisFindUnambiguousFace(MRI_SURFACE *mris, MHT *mht, VERTEX *v)
{
  int    i, fno, flist[10000], nfaces, min_fno = 0, l ;
  MHBT   *bucket ;
  MHB    *bin ;
  FACE   *f ;
  VERTEX *v1, *v2, *v3 ;
  float  ut, vt, stretch, n[3], leg[3], p[3], d[3], dot, ldot, leg2[3] ;

  bucket = MHTgetBucket(mht, v->x, v->y, v->z) ;
  bin = bucket->bins ;
  for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)  /* check each face */
  {
    n[0] = n[1] = n[2] = 0.0f ;
    f = &mris->faces[bin->fno] ;
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      n[0] += v1->nx ; n[1] += v1->ny ; n[2] += v1->nz ; 
    }
    n[0] /= VERTICES_PER_FACE ; n[1] /= VERTICES_PER_FACE ; 
    n[2] /= VERTICES_PER_FACE ; 
    for (l = 0 ; l < VERTICES_PER_FACE ; l++)
    {
      v1 = &mris->vertices[f->v[l]] ;
      switch (l)
      {
      default:
      case 0:
        v2 = &mris->vertices[f->v[1]] ;
        v3 = &mris->vertices[f->v[2]] ;
        break ;
      case 1:
        v2 = &mris->vertices[f->v[2]] ;
        v3 = &mris->vertices[f->v[0]] ;
        break ;
      case 2:
        v2 = &mris->vertices[f->v[0]] ;
        v3 = &mris->vertices[f->v[1]] ;
        break ;
      }

      VERTEX_DIF(leg, v1, v2) ;   /* leg of triangle */
      VERTEX_DIF(leg2, v1, v3) ;  /* leg of triangle */

      /* express p as point in triangle plane */
      VERTEX_DIF(p, v1, v) ;     /* vector from vertex to point in question */
      dot = DOT(p,n) ;
      p[0] -= dot*n[0] ; p[1] -= dot*n[1] ; p[2] -= dot*n[2] ; 
#if 0
      p[0] = ut*leg[0] + vt*leg2[0] ;
      p[1] = ut*leg[1] + vt*leg2[1] ;
      p[2] = ut*leg[2] + vt*leg2[2] ;
      ut = DOT(p, leg) ; vt = DOT(p, leg2) ;
#endif

      CROSS(d, leg, n) ;
      dot = DOT(d, p) ; ldot = DOT(d, leg2) ;

      /* on different side of leg from 3rd vertex */
      if (!FZERO(ldot) && !FZERO(dot) && dot*ldot < 0)
        break ;
    }
    if (l >= VERTICES_PER_FACE)
      flist[nfaces++] = bin->fno ;
  }
  if (!nfaces)
  {
    float dist, min_dist ;

    fprintf(stderr, "no faces found on sphere!\n") ;
    bin = bucket->bins ;
    min_dist = 1000000.0f ; min_fno = 0 ;
    for (nfaces = i = 0 ; i < bucket->nused ; i++, bin++)
    {
      f = &mris->faces[bin->fno] ;
      v1 = &mris->vertices[f->v[0]] ;
      v2 = &mris->vertices[f->v[1]] ;
      v3 = &mris->vertices[f->v[2]] ;
#define VDIST(v1,v2) (sqrt(SQR(v1->x-v2->x)+SQR(v1->y-v2->y)+SQR(v1->z-v2->z)))
      dist = VDIST(v1,v) + VDIST(v2, v) + VDIST(v3,v) ;
      if (dist < min_dist)
      {
        min_dist = dist ;
        min_fno = bin->fno ;
      }
    }
                printf("min_dist = %f (min_fno=%d)\n",min_dist,min_fno);
                
  }
  else
  {
    float min_stretch, curv ;
    int   vno ;

    if (nfaces > 1)
    {
      DiagBreak() ;
      curv = 1.0f ;
    }
    else
      curv = 0.0f ;
    for ( i = 0 ; i < nfaces ; i++)/* check each face */
    {
      f = &mris->faces[flist[i]] ;
      for (l = 0 ; l < VERTICES_PER_FACE ; l++)
        mris->vertices[f->v[l]].curv = curv ;
    }

    min_stretch = 100000.0f ;
    for (i = 0 ; i < nfaces ; i++)
    {
      stretch = mrisComputeFaceStretch(mris, flist[i]) ;
      if (stretch < min_stretch)
      {
        min_stretch = stretch ;
        min_fno = flist[i] ;
      }
    }
  }
  if (nfaces <= 1)
    return(min_fno) ;
  else
    return(-1) ;
}


static float
mrisComputeFaceStretch(MRI_SURFACE *mris, int fno)
{
  int    fv, n0, n1 ;
  float  stretch, max_stretch, dist, inflated_dist ;
  FACE   *f ;
  VERTEX *v0, *v1 ;

  f = &mris->faces[fno] ;
  max_stretch = -1.0f ;
  for (fv = 0 ; fv < VERTICES_PER_FACE ; fv++)
  {
    n0 = f->v[fv] ;
    n1 = fv < VERTICES_PER_FACE - 1 ? f->v[fv+1] : f->v[0] ;
    v0 = &mris->vertices[n0] ; v1 = &mris->vertices[n1] ;
    inflated_dist = 
      SQR(v0->cx-v1->cx) + SQR(v0->cy-v1->cy) + SQR(v0->cz-v1->cz);
    dist = 
      SQR(v0->origx-v1->origx) + 
      SQR(v0->origy-v1->origy) + SQR(v0->origz-v1->origz);
    if (!FZERO(dist))
    {
      stretch = inflated_dist / dist ;
      if (stretch > max_stretch)
        max_stretch = stretch ;
    }
  }
  return(max_stretch) ;
}

static int
mrisCalculateFaceCentroid(MRI_SURFACE *mris, int fno, float *px, float *py, 
                          float *pz)
{
  float  x, y, z ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *f ;

  f = &mris->faces[fno] ;
  v0 = &mris->vertices[f->v[0]] ; v1 = &mris->vertices[f->v[1]] ; 
  v2 = &mris->vertices[f->v[2]] ;

  /* first bisect v1->v2 line */

  x = (v1->x + v2->x) / 2.0f ;
  y = (v1->y + v2->y) / 2.0f ;
  z = (v1->z + v2->z) / 2.0f ;
  
  /* now bisect v0->bisector line */
  *px = (v0->x + x) / 2.0f ;
  *py = (v0->y + y) / 2.0f ;
  *pz = (v0->z + z) / 2.0f ;
  return(NO_ERROR) ;
}
static int
mrisCalculateOriginalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                  float *px, float *py, float *pz)
{
  float  x, y, z ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *f ;

  f = &mris->faces[fno] ;
  v0 = &mris->vertices[f->v[0]] ; v1 = &mris->vertices[f->v[1]] ; 
  v2 = &mris->vertices[f->v[2]] ;

  /* first bisect v1->v2 line */

  x = (v1->origx + v2->origx) / 2.0f ;
  y = (v1->origy + v2->origy) / 2.0f ;
  z = (v1->origz + v2->origz) / 2.0f ;
  
  /* now bisect v0->bisector line */
  *px = (v0->origx + x) / 2.0f ;
  *py = (v0->origy + y) / 2.0f ;
  *pz = (v0->origz + z) / 2.0f ;
  return(NO_ERROR) ;
}
static int
mrisCalculateCanonicalFaceCentroid(MRI_SURFACE *mris, int fno, 
                                  float *px, float *py, float *pz)
{
  float  x, y, z ;
  VERTEX *v0, *v1, *v2 ;
  FACE   *f ;

  f = &mris->faces[fno] ;
  v0 = &mris->vertices[f->v[0]] ; v1 = &mris->vertices[f->v[1]] ; 
  v2 = &mris->vertices[f->v[2]] ;

  /* first bisect v1->v2 line */

  x = (v1->cx + v2->cx) / 2.0f ;
  y = (v1->cy + v2->cy) / 2.0f ;
  z = (v1->cz + v2->cz) / 2.0f ;
  
  /* now bisect v0->bisector line */
  *px = (v0->cx + x) / 2.0f ;
  *py = (v0->cy + y) / 2.0f ;
  *pz = (v0->cz + z) / 2.0f ;
  return(NO_ERROR) ;
}
