/**
 * @brief register surface2 to surface1 
 *
 * This function will register surface2 to surface1 using ICP or input linear
 * transform (can be obtained through volume registration). It then maps
 * the coordinate function (spherical morphing for surface1) to surface2
 */
/*
 * Original Author: X. Han
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

#include <iostream>
#include <fstream>
#include "ANN.h"

extern "C"
{
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include "timer.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "mrishash.h"
#include "mri_identify.h"
#include "icosahedron.h"
#include "version.h"
#include "fmarching3dnband.h"
}

#define MAX_DATA_NUMBERS 200
#define DEBUG 0

typedef struct _double_3d
{
  double x;
  double y;
  double z;
}
double3d ;

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
 a[k][l]=h+s*(g-h*tau);
#define TINY 1.0e-20;


void jacobi(float **a, int n, float *d, float** v,int * nrot);
void ludcmp(double** a,int n,int* indx,double* d);
void lubksb(double** a,int n,int* indx,double* b);

void MRISsampleTemplateMappingToSource(MRI_SURFACE *mris, MRI_SURFACE *mris_template);


int main(int argc, char *argv[]) ;

int framesave = 0;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;


char *out_name = NULL;

int debugflag = 0;
int debugvtx = 0;
int pathflag = 0;

int normflag = 0;

int register_flag = 0;

const char *Progname ;

double transformS(double3d *V1a, double3d *V2a, int N, double TR[3][3],double shift[3]);
void FindClosest(MRI_SURFACE *TrueMesh, ANNkd_tree *annkdTree, MRI_SURFACE *EstMesh, double3d *closest);
double v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug, double *sopt, double *topt);
void register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2);
static int mrisSetVertexFaceIndex(MRI_SURFACE *mris, int vno, int fno);

/* the following two are used when applying a lta transform to the surface */
MRI          *mri = 0;
MRI          *mri_dst = 0;

static int invert = 0 ;
static char *xform_fname = NULL;

static char *template_fname = NULL;
static char *template_map_name = "sphere.reg"; //Default

using namespace std;

int main(int argc, char *argv[])
{
  char **av, *in_surf_fname, *out_fname;
  int fno, vno0, vno1, vno2;
  int nargs, ac, msec;
  Timer then;

  MRI_SURFACE *mris, *mris_template, *mris_template_map;
  FACE *face;

  LTA          *lta = 0;
  int          transform_type;

  nargs = handleVersionOption(argc, argv, "mris_indirect_morph");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  then.reset() ;
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

  if (argc != 3)
  {
    printf("Incorrect number of arguments, argc = %d\n", argc);
    usage_exit();
  }

  if (!template_fname)
  {
    fprintf(stderr, "A template surface file name must be specified using the -template option. Exit\n");
    exit(1);
  }

  in_surf_fname = argv[1] ;
  out_fname = argv[2] ;

  printf("Reading input surface file\n");
  mris = MRISread(in_surf_fname);
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s:could not read surface %s", Progname, in_surf_fname);

  if (template_fname != NULL)
  {
    printf("Reading template surface file\n");
    mris_template = MRISread(template_fname);
    if (!mris_template)
      ErrorExit(ERROR_NOFILE, "%s:could not read surface %s", Progname, template_fname);
    printf("Reading template map file (default is sphere.reg)\n");
    /* Save the original true surface cooridnates */
    MRISsaveVertexPositions(mris_template, TMP_VERTICES) ;

    //   if(keep_map_vol_geometry){
    /* Need to keep the volume geometry of the original map, right ? */
    /* Or not?? Maybe not needed if all surfaces share same RAS. Seems so */
    // need full path in order for the following "read" to work
    //  mris_template_map = MRISread(template_map_name);
    // }

    if (MRISreadVertexPositions(mris_template, template_map_name) != NO_ERROR)
      return(Gerror) ;
    /* Store the coordinate-function (mapping) to ORIGINAL_VERTICES */
    MRISsaveVertexPositions(mris_template, ORIGINAL_VERTICES) ;
    /* Restore the true surface coordinates for closest point mapping*/
    MRISrestoreVertexPositions(mris_template, TMP_VERTICES);
  }

  /* added for allowing using faces; maybe unnecessary */
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    vno0 = face->v[0] ;
    vno1 = face->v[1] ;
    vno2 = face->v[2] ;
    mrisSetVertexFaceIndex(mris, vno0, fno) ;
    mrisSetVertexFaceIndex(mris, vno1, fno) ;
    mrisSetVertexFaceIndex(mris, vno2, fno) ;
  }
  for (fno = 0 ; fno < mris_template->nfaces ; fno++)
  {
    face = &mris_template->faces[fno] ;
    vno0 = face->v[0] ;
    vno1 = face->v[1] ;
    vno2 = face->v[2] ;
    mrisSetVertexFaceIndex(mris_template, vno0, fno) ;
    mrisSetVertexFaceIndex(mris_template, vno1, fno) ;
    mrisSetVertexFaceIndex(mris_template, vno2, fno) ;
  }

  if (register_flag && (xform_fname != NULL))
  {
    printf("A LTA is specified, so only apply the LTA, and no ICP will be computed\n");
    register_flag = 0;
  }

  if (register_flag)
  {
    printf("Register input surface to template surface (using ICP)\n");
    register2to1(mris_template, mris);
  }

  if (xform_fname != NULL)
  {
    printf("Apply the given LTA xfrom to align input surface to template surface\n ...");
    // read transform
    transform_type =  TransformFileNameType(xform_fname);

    if (transform_type == MNI_TRANSFORM_TYPE ||
        transform_type == TRANSFORM_ARRAY_TYPE ||
        transform_type == REGISTER_DAT ||
        transform_type == FSLREG_TYPE
       )
    {
      printf("Reading transform ...\n");
      lta = LTAreadEx(xform_fname) ;
      if (!lta)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
                  Progname, xform_fname) ;

      if (transform_type == FSLREG_TYPE)
      {
        if (mri == 0 || mri_dst == 0)
        {
          fprintf(stderr, "ERROR: fslmat does not have information on the src and dst volumes\n");
          fprintf(stderr, "ERROR: you must give options '-src' and '-dst' to specify the src and dst volume infos\n");
        }


        LTAmodifySrcDstGeom(lta, mri, mri_dst); // add src and dst information
        LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      }

      if (lta->xforms[0].src.valid == 0)
      {
        if (mri == 0)
        {
          fprintf(stderr, "The transform does not have the valid src volume info.\n");
          fprintf(stderr, "Either you give src volume info by option --src or\n");
          fprintf(stderr, "make the transform to have the valid src info.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }
        else
        {
          LTAmodifySrcDstGeom(lta, mri, NULL); // add src information
          //   getVolGeom(mri, &lt->src);
        }
      }
      if (lta->xforms[0].dst.valid == 0)
      {
        if (mri_dst == 0)
        {
          fprintf(stderr, "The transform does not have the valid dst volume info.\n");
          fprintf(stderr, "Either you give src volume info by option --dst or\n");
          fprintf(stderr, "make the transform to have the valid dst info.\n");
          fprintf(stderr, "If the dst was average_305, then you can set\n");
          fprintf(stderr, "environmental variable USE_AVERAGE305 true\n");
          fprintf(stderr, "without giving the dst volume for RAS-to-RAS transform.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }
        else
        {
          LTAmodifySrcDstGeom(lta, NULL, mri_dst); // add  dst information
        }
      }
    }
    else
    {
      ErrorExit(ERROR_BADPARM, "transform is not of MNI, nor Register.dat, nor FSLMAT type");
    }

    if (invert)
    {
      VOL_GEOM vgtmp;
      LT *lt;
      MATRIX *m_tmp = lta->xforms[0].m_L ;
      lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
      MatrixFree(&m_tmp) ;
      lt = &lta->xforms[0];
      if (lt->dst.valid == 0 || lt->src.valid == 0)
      {
        fprintf(stderr, "WARNING:***************************************************************\n");
        fprintf(stderr, "WARNING:dst volume infor is invalid.  Most likely produce wrong inverse.\n");
        fprintf(stderr, "WARNING:***************************************************************\n");
      }
      copyVolGeom(&lt->dst, &vgtmp);
      copyVolGeom(&lt->src, &lt->dst);
      copyVolGeom(&vgtmp, &lt->src);
    }

    MRIStransform(mris, mri, lta, mri_dst) ;

    if (mri)
      MRIfree(&mri);
    if (mri_dst)
      MRIfree(&mri_dst);

    if (lta)
      LTAfree(&lta);
  }   /* if (xform_fname != NULL) */

  printf("Now mapping the spherical map from the template surface to the input surface\n");

  /* The sampled coordinate function is stored at original-vertex */
  MRISsampleTemplateMappingToSource(mris, mris_template);

  /* Move the cooridnate-function to active vertex */
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISwrite(mris, out_fname) ;

  /* Free memories */
  MRISfree(&mris);
  MRISfree(&mris_template);

  msec = then.milliseconds() ;
  fprintf(stderr, "indirect spherical mapping or morphing took %2.2f hours\n",
          (float)msec/(1000.0f*60.0f*60.0f));

  return 0;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

/* --------------------------------------------- */
static void print_usage(void)
{
  fprintf(stdout, "USAGE: %s [options] input_surface surface_map_name \n",Progname) ;
  fprintf(stdout, "\n");
  fprintf(stdout, "Options:\n");
  fprintf(stdout, "   -template  %%s template surface name\n");
  fprintf(stdout, "   -template_map  %%s template map name (sphere or sphere.reg)\n");
  fprintf(stdout, "   -register  use ICP registration to align surface to template\n");
  fprintf(stdout, "   -xform %%s apply LTA transform to align input surface to template\n");
  fprintf(stdout, "   -invert  reversely apply -xform \n");
  fprintf(stdout, "   -src %%s  source volume for -xform \n");
  fprintf(stdout, "   -dst %%s  target volume for -xform \n");
  fprintf(stdout, "\n");
  std::cout << getVersion() << std::endl;
  printf("\n");

}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
    "This program aligns the input surface to the template surface through ICP\n"
    "registration or a given LTA transform, and then maps the spherical morphing\n"
    "originally computed at the template surface to the input surface. \n"
    "\n"
    "OPTIONS\n"
    "\n"
    " -template surface_file_name\n"
    "\n"
    "  template surface to which the input surface will be mapped\n"
    "\n"
    " -template_map string (default = sphere.reg)\n"
    "\n"
    " Spherical morphing of the template that is to be sampled to the input surface\n"
    "\n"
    "  -xform LTA_file_name\n"
    "\n"
    "    LTA transform to align input surface to the template surface, \n"
    "    paint or w (for FreeSurfer paint files)."
    "\n"
    " -invert\n"
    "\n"
    "   Inversely apply the transform given by -xform option \n"
    "\n"
    " -src MRI_volume_file_name\n"
    "\n"
    "   source volume for the LTA transform given by -xform\n"
    "\n"
    " -dst MRI_volume_file_name\n"
    "\n"
    "   target volume for the LTA transform given by -xform\n"
    "\n"
    "  -register\n"
    "\n"
    "  Perform ICP registration to align input surface to template instead of using given LTA transform."
    "\n");


  exit(1) ;
}


/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stdout, "%s\n", getVersion().c_str()) ;
  exit(1) ;
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
  if (!stricmp(option, "help"))
    print_help() ;
  else if (!stricmp(option, "version"))
    print_version() ;
  else if (!stricmp(option, "template"))
  {
    template_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "template surface file is %s\n", template_fname);
  }
  else if (!stricmp(option, "template_map"))
  {
    template_map_name = argv[2];
    nargs = 1;
    fprintf(stderr, "template map name is %s\n", template_map_name);
  }
  else if (!stricmp(option, "register"))
  {
    register_flag = 1;
    printf("Perform ICP rigid registration of the input surface to the template\n");
  }
  else if (!stricmp(option, "debug"))
  {
    debugflag = 1;
    debugvtx = atoi(argv[2]);
    nargs = 1;
  }
  else if (!stricmp(option, "norm"))
  {
    normflag = 1;
    printf("Normalize the sampled mapping to same length\n");
  }
  else if (!stricmp(option, "xform"))
  {
    xform_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "transform file name is %s\n", xform_fname);
  }
  else if (!stricmp(option, "invert"))
  {
    invert = 1;
    fprintf(stderr, "Inversely apply the given LTA transform\n");
  }
  else if (!stricmp(option, "src"))
  {
    fprintf(stderr, "src volume for the given transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the src volume...\n");
    mri = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else if (!stricmp(option, "dst"))
  {
    fprintf(stderr, "dst volume for the transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the dst volume...\n");
    mri_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_dst)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else
  {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    print_help() ;
    exit(1) ;
  }

  return(nargs) ;
}


double v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug, double *sopt, double *topt)
{
  /* Compute the distance from point P0 to a face of the surface mesh */
  /* (sopt, topt) determines the closest point inside the triangle from P0 */

  double a, b, c, d, e, f, det, s, t, invDet;
  double numer, denom, tmp0, tmp1;
  VERTEX *V1, *V2, *V3;
  FACE *face;

  VERTEX E0, E1, D;

  face = &mri_surf->faces[face_number];
  V1 = &mri_surf->vertices[face->v[0]];
  V2 = &mri_surf->vertices[face->v[1]];
  V3 = &mri_surf->vertices[face->v[2]];

  E0.x = V2->x - V1->x;
  E0.y = V2->y - V1->y;
  E0.z = V2->z - V1->z;
  E1.x = V3->x - V1->x;
  E1.y = V3->y - V1->y;
  E1.z = V3->z - V1->z;
  D.x = V1->x - P0->x;
  D.y = V1->y - P0->y;
  D.z = V1->z - P0->z;

  a = E0.x *E0.x + E0.y * E0.y + E0.z *E0.z;
  b = E0.x *E1.x + E0.y * E1.y + E0.z *E1.z;
  c = E1.x *E1.x + E1.y * E1.y + E1.z *E1.z;
  d = E0.x *D.x + E0.y * D.y + E0.z *D.z;
  e = E1.x *D.x + E1.y * D.y + E1.z *D.z;
  f = D.x *D.x + D.y * D.y + D.z *D.z;

  det = a*c - b*b;
  s = b*e - c*d;
  t = b*d - a*e;

  if (debug) printf("det = %g\n", det);
  if (s + t <= det)
  {
    if (s < 0)
    {
      if (t<0)
      {
        /* Region 4 */
        if ( d < 0)
        { /* Minimum on edge t = 0 */
          s = (-d >= a ? 1 : -d/a);
          t = 0;
        }
        else if (e < 0)
        { /* Minimum on edge s = 0 */
          t = (-e >= c ? 1 : -e/c);
          s = 0;
        }
        else
        {
          s = 0;
          t = 0;
        }
        if (debug) printf("region 4,  d= %g, e = %g, s =%g, t =%g\n", d, e, s, t);
      }
      else
      {
        /* Region 3 */
        s = 0;
        t = ( e >= 0 ? 0 : (-e >= c ? 1 : (-e/c)));
        if (debug) printf("region 3, s =%g, t =%g\n", s, t);
      }
    }
    else if (t < 0)
    {
      /* Region 5 */
      t = 0;
      s = (d >= 0 ? 0 :(-d >= a ? 1 : (-d/a)));
      if (debug) printf("region 5, s =%g, t =%g\n", s, t);
    }
    else
    {
      /* Region 0 */
      invDet = 1/det;
      s *= invDet;
      t *= invDet;
      if (debug) printf("region 0, s =%g, t =%g\n", s, t);
    }

  }
  else
  {
    if ( s < 0 )
    {
      /* Region 2 */
      tmp0 = b + d;
      tmp1 = c + e;
      if (tmp1 > tmp0)
      {
        numer = tmp1 - tmp0;
        denom = a - b - b +c;
        s = (numer >= denom ? 1 : numer/denom);
        t = 1-s;

      }
      else
      {
        s = 0;
        /* t = (e >= 0 ? 0 : (-e >= c ? 0 > = c + e = tmp1) */
        t = (tmp1 <= 0 ? 1 : (e >= 0 ? 0 : -e/c));
      }
      if (debug) printf("region 2, s =%g, t =%g\n", s, t);
    }
    else if ( t < 0)
    {
      /* Region 6 */
      tmp0 = b + e;
      tmp1 = a + d;
      if (tmp1 > tmp0)
      { /* Minimum at line s + t = 1 */
        numer = tmp1 - tmp0; /* Positive */
        denom = a + c - b -b;
        t = (numer >= denom ? 1 : (numer/denom));
        s = 1 - t;
      }
      else
      { /* Minimum at line t = 0 */
        s = (tmp1 <= 0 ? 1 : (d >= 0 ? 0 : -d/a));
        t = 0;
      }
      if (debug) printf("region 6, s =%g, t =%g\n", s, t);
    }
    else
    {
      /* Region 1 */
      numer = c + e - b - d;
      if (numer <= 0)
      {
        s = 0;
      }
      else
      {
        denom = a + c - b - b; /* denom is positive */
        s = (numer >= denom ? 1 : (numer/denom));
      }
      t = 1-s;
      if (debug) printf("region 1, s =%g, t =%g\n", s, t);
    }
  }

  if ( s < 0  || s > 1 || t < 0 || t > 1)
  {
    printf("Error in computing s and t \n");
  }

  *sopt = s;
  *topt = t;
  /* return (sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f)); */
  /* The square-root will be taken later to save time */
  return (a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);

}

void MRISsampleTemplateMappingToSource(MRI_SURFACE *mris, MRI_SURFACE *mris_template)
{
  /* this function assumes the mapping function is stored as the "INFLATED vertex" of mris_template */

  int index, k, facenumber, closestface;
  double sopt, topt, tmps, tmpt; /* triangle parametrization parameters */
  double value, distance;
  double Radius, length, scale;
  double sumx, sumy, sumz, sumweight, weight;
  VERTEX *vertex, *V1, *V2, *V3;
  FACE *face;
  ANNpointArray QueryPt;

  ANNpointArray pa = annAllocPts(mris_template->nvertices, 3);

  for (index = 0; index < mris_template->nvertices; index++)
  {
    pa[index][0] = mris_template->vertices[index].x;
    pa[index][1] = mris_template->vertices[index].y;
    pa[index][2] = mris_template->vertices[index].z;
  }

  // construct and initialize the tree
  ANNkd_tree *annkdTree = new ANNkd_tree(pa, mris_template->nvertices, 3);
  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];


  QueryPt = annAllocPts(1,3);

  if (normflag)
  {
    /* Compute the radius of the template mapping */
    vertex = &mris_template->vertices[0];
    Radius = sqrt(vertex->origx *vertex->origx + vertex->origy*vertex->origy
                  + vertex->origz*vertex->origz);
    printf("Radius = %g\n", Radius);
  }

  for (index = 0; index < mris->nvertices; index++)
  {
    vertex = &mris->vertices[index];
    QueryPt[0][0] = vertex->x;
    QueryPt[0][1] = vertex->y;
    QueryPt[0][2] = vertex->z;

    annkdTree->annkSearch( // search
      QueryPt[0],       // query point
      1,   // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,  // distance (returned)
      0);   // error bound

    /* annIndex gives the closest vertex on mris_template to the vertex #index of mris */
    /* Now need to find the closest face in order to perform linear interpolation */
    distance = 1000.0;
    for (k=0; k < mris_template->vertices[annIndex[0]].num; k++)
    {

      facenumber =  mris_template->vertices[annIndex[0]].f[k]; /* index of the k-th face */
      if (facenumber < 0 || facenumber >= mris_template->nfaces) continue;
      value = v_to_f_distance(vertex, mris_template, facenumber, 0, &tmps, &tmpt);

      if (distance > value)
      {
        distance = value;
        closestface = facenumber;
        sopt = tmps;
        topt = tmpt;
      }
    }

    //    distance = sqrt(distance); /* distance is no use here */
    /* Now interpolate the functions values to point (sopt,topt) */
    /* Linear interpolation is quite simple, just equal to
       sopt*vertex3 + topt*vertex2 + (1-s-t)*vertex1. check the function
       v_to_f_distance() to see how V1, V2, V3 is defined
     */
    /* Use the orig_vertex to store the interpolated cooridnate function */
    face = &mris_template->faces[closestface];
    V1 = &mris_template->vertices[face->v[0]];
    V2 = &mris_template->vertices[face->v[1]];
    V3 = &mris_template->vertices[face->v[2]];
    sumx = 0.0;
    sumy = 0.0;
    sumz = 0.0;
    sumweight = 0.0;
    weight = 1.0/(1e-20 + (V1->x - vertex->x)*(V1->x - vertex->x) + (V1->y - vertex->y)*(V1->y - vertex->y) + (V1->z - vertex->z)*(V1->z - vertex->z));
    sumx += weight*V1->origx;
    sumy += weight*V1->origy;
    sumz += weight*V1->origz;
    sumweight += weight;
    weight = 1.0/(1e-20 + (V2->x - vertex->x)*(V2->x - vertex->x) + (V2->y - vertex->y)*(V2->y - vertex->y) + (V2->z - vertex->z)*(V2->z - vertex->z));
    sumx += weight*V2->origx;
    sumy += weight*V2->origy;
    sumz += weight*V2->origz;
    sumweight += weight;
    weight = 1.0/(1e-20 + (V3->x - vertex->x)*(V3->x - vertex->x) + (V3->y - vertex->y)*(V3->y - vertex->y) + (V3->z - vertex->z)*(V3->z - vertex->z));
    sumx += weight*V3->origx;
    sumy += weight*V3->origy;
    sumz += weight*V3->origz;
    sumweight += weight;

    vertex->origx = sumx /(sumweight + 1e-30);
    vertex->origy = sumy /(sumweight + 1e-30);
    vertex->origz = sumz /(sumweight + 1e-30);

    //    vertex->origx = sopt*V3->origx + topt*V2->origx + (1.0-sopt-topt)*V1->origx;
    //vertex->origy = sopt*V3->origy + topt*V2->origy + (1.0-sopt-topt)*V1->origy;
    //vertex->origz = sopt*V3->origz + topt*V2->origz + (1.0-sopt-topt)*V1->origz;
    if (normflag)
    {
      length = sqrt(vertex->origx *vertex->origx + vertex->origy*vertex->origy
                    + vertex->origz*vertex->origz);
      scale = Radius/(length + 1e-20);
      vertex->origx *= scale;
      vertex->origy *= scale;
      vertex->origz *= scale;
    }

  }



  if (annkdTree) delete annkdTree;
  if (annIndex) delete annIndex;
  if (annDist) delete annDist;
  if (QueryPt) delete QueryPt;

  return;
}

void register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2)
{

  double error_old, error_new;
  int iter = 0;
  double TR[3][3];
  double shift[3];
  double3d *mesh2avtx, *closevtx2a;
  int index, k;

  /* Build the ANN tree repeatedly is time consuming, so
   *  move it outside of the function
   */
  ANNpointArray pa = annAllocPts(Surf1->nvertices, 3);

  for (index = 0; index < Surf1->nvertices; index++)
  {
    pa[index][0] = Surf1->vertices[index].x;
    pa[index][1] = Surf1->vertices[index].y;
    pa[index][2] = Surf1->vertices[index].z;
  }

  ANNkd_tree *annkdTree = new ANNkd_tree(pa, Surf1->nvertices, 3);


  /* This initialization is necessary */
  TR[0][0] = TR[1][1] = TR[2][2] = 1;
  TR[0][1] = TR[0][2] = TR[1][0] = TR[1][2] = TR[2][0] = TR[2][1] = 0;

  shift[0] = 0;
  shift[1] = 0;
  shift[2]= 0;

  mesh2avtx = (double3d*) malloc(Surf2->nvertices*sizeof(double3d));
  closevtx2a = (double3d*) malloc(Surf2->nvertices*sizeof(double3d));

  /* Initialize */
  for (index = 0; index < Surf2->nvertices; index++)
  {
    mesh2avtx[index].x = Surf2->vertices[index].x;
    mesh2avtx[index].y = Surf2->vertices[index].y;
    mesh2avtx[index].z = Surf2->vertices[index].z;
  }

  error_old = 1000.0;
  error_new = 900.0;
  while ((error_old - error_new) > 0.0001 && iter < 50)
  {
    error_old = error_new;
    /* For each vertex in Surf2, find its closest vertex in Surf1,
     * and record the coordinates in closevtx2a
     */
    FindClosest(Surf1, annkdTree, Surf2, closevtx2a);

    // Find the rigid transformation
    error_new = transformS(mesh2avtx, closevtx2a, Surf2->nvertices, TR, shift);

    for (k = 0; k < Surf2->nvertices; k++)
    {
      Surf2->vertices[k].x = mesh2avtx[k].x;
      Surf2->vertices[k].y  = mesh2avtx[k].y;
      Surf2->vertices[k].z  = mesh2avtx[k].z;
    }

    iter ++;
    if (DEBUG) printf(" iteration %d, error = %15.4f\n", iter, error_new);
  }

  // free memory
  delete pa;
  if (annkdTree) delete annkdTree;

  free(mesh2avtx);
  free(closevtx2a);

  return;
}

void FindClosest(MRI_SURFACE *TrueMesh, ANNkd_tree *annkdTree, MRI_SURFACE *EstMesh, double3d *closest)
{
  int index;

  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];

  ANNpoint Qpt;

  Qpt = (ANNpoint)malloc(3*sizeof(ANNcoord));

  for (index =0; index < EstMesh->nvertices; index++)
  {
    // this is a duplicate, lame, but....ugh, to get in and out of libraries...
    Qpt[0] = EstMesh->vertices[index].x;
    Qpt[1] = EstMesh->vertices[index].y;
    Qpt[2] = EstMesh->vertices[index].z;

    annkdTree->annkSearch( // search
      Qpt,  // query point
      1,    // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,   // distance (returned)
      0);      // error bound

    closest[index].x = TrueMesh->vertices[annIndex[0]].x;
    closest[index].y = TrueMesh->vertices[annIndex[0]].y;
    closest[index].z = TrueMesh->vertices[annIndex[0]].z;
  }

  if (annIndex) delete annIndex;
  if (annDist) delete annDist;
  free(Qpt);

  return;
}


double transformS(double3d *V1a, double3d *V2a, int N, double TR[3][3],double shift[3])
// transform V1 to fit V2
// V1 is equavilent to left frame in Horn paper
{
  double3d centroid1a, centroid2a;
  double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
  float **M, **v, *d;
  float *n;
  double R[3][3],x,y,z;
  double scale1,scale2;
  float dummy;
  double error = 0;

  int k,l, nrot;
  int count;

  centroid1a.x = 0;
  centroid1a.y = 0;
  centroid1a.z = 0;
  centroid2a.x = 0;
  centroid2a.y = 0;
  centroid2a.z = 0;

  count = 0;
  for (k = 0; k < N; k ++)
  {
    //   if (label[k]){
    centroid1a.x += V1a[k].x;
    centroid1a.y += V1a[k].y;
    centroid1a.z += V1a[k].z;
    centroid2a.x += V2a[k].x;
    centroid2a.y += V2a[k].y;
    centroid2a.z += V2a[k].z;
    ++count;
  }

  /* Compute the centroid of each point set */
  centroid1a.x /= (double)count;
  centroid1a.y /= (double)count;
  centroid1a.z /= (double)count;
  centroid2a.x /= (double)count;
  centroid2a.y /= (double)count;
  centroid2a.z /= (double)count;

  Sxx = 0;
  Sxy = 0;
  Sxz = 0;
  Syx = 0;
  Syy = 0;
  Syz = 0;
  Szx = 0;
  Szy = 0;
  Szz = 0;

  /* Centralize respective point data sets */
  scale1 = 0;
  scale2 = 0;
  for (k = 0; k < N; k ++)
  {
    V1a[k].x -= centroid1a.x;
    V1a[k].y -= centroid1a.y;
    V1a[k].z -= centroid1a.z;

    V2a[k].x -= centroid2a.x;
    V2a[k].y -= centroid2a.y;
    V2a[k].z -= centroid2a.z;
  }
  for (k = 0; k < N; k++)
  {
    /* if (label[k]){ */
    scale1+=(V1a[k].x * V1a[k].x +  V1a[k].y * V1a[k].y + V1a[k].z * V1a[k].z);
    Sxx += V1a[k].x * V2a[k].x;
    Sxy += V1a[k].x * V2a[k].y;
    Sxz += V1a[k].x * V2a[k].z;

    Syx += V1a[k].y * V2a[k].x;
    Syy += V1a[k].y * V2a[k].y;
    Syz += V1a[k].y * V2a[k].z;

    Szx += V1a[k].z * V2a[k].x;
    Szy += V1a[k].z * V2a[k].y;
    Szz += V1a[k].z * V2a[k].z;
    // }
  }

  M = (float**)malloc(4*sizeof(float*));
  M --;
  for (k = 1; k <= 4; k++)
  {
    n = (float*)malloc(4*sizeof(float));
    M[k] = n-1;
  }

  v = (float**)malloc(4*sizeof(float*));
  v --;
  for (k = 1; k <= 4; k++)
  {
    n = (float*)malloc(4*sizeof(float));
    v[k] = n-1;
  }

  d = (float*)malloc(4*sizeof(float));
  d --;

  M[1][1] = Sxx+Syy+Szz;
  M[1][2] = Syz-Szy;
  M[1][3] = Szx-Sxz;
  M[1][4] = Sxy-Syx;
  M[2][1] = Syz-Szy;
  M[2][2] = Sxx-Syy-Szz;
  M[2][3] = Sxy+Syx;
  M[2][4] = Sxz+Szx;
  M[3][1] = Szx-Sxz;
  M[3][2] = Sxy+Sxy;
  M[3][3] = -Sxx+Syy-Szz;
  M[3][4] = Syz+Szy;
  M[4][1] = Sxy-Syx;
  M[4][2] = Sxz+Szx;
  M[4][3] = Szy+Syz;
  M[4][4] = -Sxx-Syy+Szz;

  for (k = 1; k <= 4; k++)
    for (l = 1; l <= 4; l++)
      M[k][l] /= (double)(N);

  /* printf("\nThe Matrix = \n");
     printf("\t %15.9f %15.9f %15.9f %15.9f\n", M[1][1],M[1][2],M[1][3], M[1][4]);
     printf("\t %15.9f %15.9f %15.9f %15.9f\n", M[2][1],M[2][2],M[2][3], M[2][4]);
     printf("\t %15.9f %15.9f %15.9f %15.9f\n", M[3][1],M[3][2],M[3][3], M[3][4]);
     printf("\t %15.9f %15.9f %15.9f %15.9f\n", M[4][1],M[4][2],M[4][3], M[4][4]);*/

  jacobi(M,4,d,v,&nrot);
  dummy = d[1];
  l = 1;
  for (k = 2; k <= 4; k++)
  {
    if (dummy < d[k])
    {
      dummy = d[k];
      l = k;
    }
  }
  for (k = 1; k <= 4; k++)
    d[k] = v[l][k];

  // printf("\nThe unit quaternion = [%f %f %f %f]\n", d[1], d[2], d[3], d[4]);
  /* R is not symmetric, because it's a rotation around an arbitrary axis, not just the origin */
  R[0][0] = d[1]*d[1] + d[2]*d[2] - d[3]*d[3] - d[4]*d[4];
  R[0][1] = 2*(d[2]*d[3] - d[1]*d[4]);
  R[0][2] = 2*(d[2]*d[4] + d[1]*d[3]);
  R[1][0] = 2*(d[2]*d[3] + d[1]*d[4]);
  R[1][1] = d[1]*d[1] - d[2]*d[2] + d[3]*d[3] - d[4]*d[4];
  R[1][2] = 2*(d[3]*d[4] - d[1]*d[2]);
  R[2][0] = 2*(d[2]*d[4] - d[1]*d[3]);
  R[2][1] = 2*(d[3]*d[4] + d[1]*d[2]);
  R[2][2] = d[1]*d[1] - d[2]*d[2] - d[3]*d[3] + d[4]*d[4];

  /* printf("\nRotation matrix R = \n");
     printf("\t %15.11f %15.11f %15.11f\n", R[0][0], R[1][0], R[2][0]);
     printf("\t %15.11f %15.11f %15.11f\n", R[0][1], R[1][1], R[2][1]);
     printf("\t %15.11f %15.11f %15.11f\n", R[0][2], R[1][2], R[2][2]);*/

  for (k = 0; k < N; k ++)
  {
    x = R[0][0] * V1a[k].x + R[1][0] * V1a[k].y + R[2][0] * V1a[k].z;
    y = R[0][1] * V1a[k].x + R[1][1] * V1a[k].y + R[2][1] * V1a[k].z;
    z = R[0][2] * V1a[k].x + R[1][2] * V1a[k].y + R[2][2] * V1a[k].z;

    V1a[k].x = x;
    V1a[k].y = y;
    V1a[k].z = z;
    // if (label[k])
    scale2+=(V1a[k].x * V2a[k].x +  V1a[k].y * V2a[k].y + V1a[k].z * V2a[k].z);
  }

  scale1 = scale2/scale1;
  //  printf ("Scaling factor: %15.4f\n", scale1);

  for (k = 0; k < N; k ++)
  {
    V1a[k].x *= scale1;
    V1a[k].y *= scale1;
    V1a[k].z *= scale1;

    //if (label[k])
    error += ((V1a[k].x-V2a[k].x)*(V1a[k].x-V2a[k].x) + (V1a[k].y-V2a[k].y)*(V1a[k].y-V2a[k].y) + (V1a[k].z-V2a[k].z)*(V1a[k].z-V2a[k].z));

    V1a[k].x += centroid2a.x;
    V1a[k].y += centroid2a.y;
    V1a[k].z += centroid2a.z;
  }

  double temp[3][3];
  /* Stores the previous transformation matrix */
  temp[0][0]=TR[0][0];
  temp[0][1]=TR[0][1];
  temp[0][2]=TR[0][2];
  temp[1][0]=TR[1][0];
  temp[1][1]=TR[1][1];
  temp[1][2]=TR[1][2];
  temp[2][0]=TR[2][0];
  temp[2][1]=TR[2][1];
  temp[2][2]=TR[2][2];

  /* Update the overall scaled rotation */
  TR[0][0]=scale1*(temp[0][0]*R[0][0]+temp[0][1]*R[1][0]+temp[0][2]*R[2][0]);
  TR[0][1]=scale1*(temp[0][0]*R[0][1]+temp[0][1]*R[1][1]+temp[0][2]*R[2][1]);
  TR[0][2]=scale1*(temp[0][0]*R[0][2]+temp[0][1]*R[1][2]+temp[0][2]*R[2][2]);

  TR[1][0]=scale1*(temp[1][0]*R[0][0]+temp[1][1]*R[1][0]+temp[1][2]*R[2][0]);
  TR[1][1]=scale1*(temp[1][0]*R[0][1]+temp[1][1]*R[1][1]+temp[1][2]*R[2][1]);
  TR[1][2]=scale1*(temp[1][0]*R[0][2]+temp[1][1]*R[1][2]+temp[1][2]*R[2][2]);

  TR[2][0]=scale1*(temp[2][0]*R[0][0]+temp[2][1]*R[1][0]+temp[2][2]*R[2][0]);
  TR[2][1]=scale1*(temp[2][0]*R[0][1]+temp[2][1]*R[1][1]+temp[2][2]*R[2][1]);
  TR[2][2]=scale1*(temp[2][0]*R[0][2]+temp[2][1]*R[1][2]+temp[2][2]*R[2][2]);

  /* The following is just the current-step transformation matrix */
  /* TR[0][0] = scale1*R[0][0];
     TR[0][1] = scale1*R[0][1];
     TR[0][2] = scale1*R[0][2];
     TR[1][0] = scale1*R[1][0];
     TR[1][1] = scale1*R[1][1];
     TR[1][2] = scale1*R[1][2];
     TR[2][0] = scale1*R[2][0];
     TR[2][1] = scale1*R[2][1];
     TR[2][2] = scale1*R[2][2]; */


  /* Update the overall shift */
  temp[0][0]=shift[0];
  temp[0][1]=shift[1];
  temp[0][2]=shift[2];
  shift[0]=scale1*(R[0][0]*(temp[0][0]-centroid1a.x)+R[1][0]*(temp[0][1]-centroid1a.y)+R[2][0]*(temp[0][2]-centroid1a.z))+centroid2a.x;
  shift[1]=scale1*(R[0][1]*(temp[0][0]-centroid1a.x)+R[1][1]*(temp[0][1]-centroid1a.y)+R[2][1]*(temp[0][2]-centroid1a.z))+centroid2a.y;
  shift[2]=scale1*(R[0][2]*(temp[0][0]-centroid1a.x)+R[1][2]*(temp[0][1]-centroid1a.y)+R[2][2]*(temp[0][2]-centroid1a.z))+centroid2a.z;

  /* The following is just the shift at the current step.
   * Note the first point data set is constantly updated/transformed every iteration
   */
  /* shift[0][0]=scale1*(R[0][0]*(-centroid1a.x)+R[1][0]*(-centroid1a.y)+R[2][0]*(-centroid1a.z))+centroid2a.x;
     shift[0][1]=scale1*(R[0][1]*(-centroid1a.x)+R[1][1]*(-centroid1a.y)+R[2][1]*(-centroid1a.z))+centroid2a.y;
     shift[0][2]=scale1*(R[0][2]*(-centroid1a.x)+R[1][2]*(-centroid1a.y)+R[2][2]*(-centroid1a.z))+centroid2a.z; */


  return(sqrt(error)/(N+1e-10));
}

void ludcmp(double** a,int n,int *indx,double* d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;

  //vv=vector(1,n);
  vv = (double*)malloc(2*n*sizeof(double));

  *d=1.0;
  for (i=1;i<=n;i++)
  {
    big=0.0;
    for (j=1;j<=n;j++)
      if ((temp=fabs(a[i][j])) > big) big=temp;
    if (big == 0.0) printf("Singular matrix in routine LUDCMP\n");
    vv[i]=1.0/big;
  }
  for (j=1;j<=n;j++)
  {
    for (i=1;i<j;i++)
    {
      sum=a[i][j];
      for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
    }
    big=0.0;
    for (i=j;i<=n;i++)
    {
      sum=a[i][j];
      for (k=1;k<j;k++)
        sum -= a[i][k]*a[k][j];
      a[i][j]=sum;
      if ( (dum=vv[i]*fabs(sum)) >= big)
      {
        big=dum;
        imax=i;
      }
    }
    if (j != imax)
    {
      for (k=1;k<=n;k++)
      {
        dum=a[imax][k];
        a[imax][k]=a[j][k];
        a[j][k]=dum;
      }
      *d = -(*d);
      vv[imax]=vv[j];
    }
    indx[j]=imax;
    if (a[j][j] == 0.0) a[j][j]=TINY;
    if (j != n)
    {
      dum=1.0/(a[j][j]);
      for (i=j+1;i<=n;i++) a[i][j] *= dum;
    }
  }
  free(vv);
}

void lubksb(double** a,int n,int* indx,double* b)
{
  int i,ii=0,ip,j;
  double sum;

  for (i=1;i<=n;i++)
  {
    ip=indx[i];
    sum=b[ip];
    b[ip]=b[i];
    if (ii)
      for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
    else if (sum) ii=i;
    b[i]=sum;
  }
  for (i=n;i>=1;i--)
  {
    sum=b[i];
    for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
    b[i]=sum/a[i][i];
  }
}


void jacobi(float **a, int n,float *d, float **v, int *nrot)

{
  int j,iq,ip,i;
  float tresh,theta,tau,t,sm,s,h,g,c,*b,*z;

  //b=vector(1,n);
  b = (float*)malloc((n+1)*sizeof(float));

  //z=vector(1,n);
  z = (float*)malloc((n+1)*sizeof(float));

  for (ip=1;ip<=n;ip++)
  {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
  for (ip=1;ip<=n;ip++)
  {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=1;i<=50;i++)
  {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++)
    {
      for (iq=ip+1;iq<=n;iq++)
        sm += fabs(a[ip][iq]);
    }
    if (sm == 0.0)
    {
      // free(z+1);
      // free(b+1);
      return;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=1;ip<=n-1;ip++)
    {
      for (iq=ip+1;iq<=n;iq++)
      {
        g=100.0*fabs(a[ip][iq]);
        if (i > 4 && fabs(d[ip])+g == fabs(d[ip])
            && fabs(d[iq])+g == fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh)
        {
          h=d[iq]-d[ip];
          if (fabs(h)+g == fabs(h))
            t=(a[ip][iq])/h;
          else
          {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=1;j<=ip-1;j++)
          {
            ROTATE(a,j,ip,j,iq)
          }
          for (j=ip+1;j<=iq-1;j++)
          {
            ROTATE(a,ip,j,j,iq)
          }
          for (j=iq+1;j<=n;j++)
          {
            ROTATE(a,ip,j,iq,j)
          }
          for (j=1;j<=n;j++)
          {
            ROTATE(v,j,ip,j,iq)
          }
          ++(*nrot);
        }
      }
    }
    for (ip=1;ip<=n;ip++)
    {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  printf("Too many iterations in routine JACOBI\n");
}
/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description
          Search the face for vno and set the v->n[] field
          appropriately.
------------------------------------------------------*/
static int
mrisSetVertexFaceIndex(MRI_SURFACE *mris, int vno, int fno)
{
  VERTEX  *v ;
  FACE    *f ;
  int     n, i ;

  v = &mris->vertices[vno] ;
  f = &mris->faces[fno] ;

  for (n = 0 ; n < VERTICES_PER_FACE ; n++)
  {
    if (f->v[n] == vno)
      break ;
  }
  if (n >= VERTICES_PER_FACE)
    return(ERROR_BADPARM) ;

  for (i = 0 ; i < v->num ; i++)
    if (v->f[i] == fno)
      v->n[i] = n ;

  return(n) ;
}


#undef ROTATE
