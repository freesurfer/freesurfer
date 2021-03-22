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


/* Given a white and pial surface, compute the thickness are every point
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


static float max_thickness = 10.0 ;


int main(int argc, char *argv[]) ;

int framesave = 0;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static int mrisSetVertexFaceIndex(MRI_SURFACE *mris, int vno, int fno);

double v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug);

const char *trgtypestring = "paint";
int trgtype = MRI_VOLUME_TYPE_UNKNOWN;

int debugflag = 0;
int debugvtx = 0;

char *middle_name = NULL; /* compute a middle surface */

int tflag = 0; /* if 1, using vertex-to-face distance to measure thickness */

const char *Progname ;

MRI *ComputeThickness(MRI_SURFACE *Mesh1, MRI_SURFACE *Mesh2, MRI *mri_res);

using namespace std;

int main(int argc, char *argv[])
{
  char **av, *surf1_name, *surf2_name;
  char *out_name;

  int nargs, ac;
  int total, index, fno, vno0, vno1, vno2;

  double scalar, std, maxV, minV, meanV, absMean;

  MRI *resVal1, *resVal;

  MRI_SURFACE *Surf1, *Surf2, *mris;

  VERTEX *vertex;
  FACE *face;

  nargs = handleVersionOption(argc, argv, "mris_compute_thickness");
  if (nargs && argc - nargs == 1)
    exit (0);
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

  /* command line: <surf1>  <surf2> <output> */

  if (argc != 4)
  {
    printf("Incorrect number of arguments, argc = %d\n", argc);
    usage_exit();
  }

  surf1_name = argv[1];
  surf2_name = argv[2];
  out_name = argv[3];

  if ( trgtypestring == NULL)
  {
    printf("Please specify output data type!\n");
    usage_exit();
  }

  printf("Reading first surface file\n");
  Surf1 = MRISread(surf1_name);
  if (!Surf1)
    ErrorExit(ERROR_NOFILE, "%s:could not read surface %s", Progname, surf1_name);

  printf("Surface1 has %d vertices\n", Surf1->nvertices);

  printf("Reading second surface file\n");
  Surf2 = MRISread(surf2_name);
  if (!Surf2)
    ErrorExit(ERROR_NOFILE, "%s:could not read surface %s", Progname, surf2_name);

  printf("Surface2 has %d vertices\n", Surf2->nvertices);

  if ((int)Surf1->nvertices != (int) Surf2->nvertices)
    ErrorExit(ERROR_NOFILE, "%s:the two surfaces have different numbers of vertices", Progname);

  /* added for allowing using faces; maybe unnecessary */
  mris = Surf1;
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
  mris = Surf2;
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


  /* one way */
  resVal1 = ComputeThickness(Surf1, Surf2, NULL);
  /* The other way */
  resVal = ComputeThickness(Surf2, Surf1, NULL);

  printf("Compute statistics \n");
  maxV = -1000.0;
  minV = 1000.0;
  meanV=0.0;
  absMean = 0.0;
  total = 0;
  for (index = 0; index < Surf1->nvertices; index++)
  {
    scalar  = MRIgetVoxVal(resVal1,index,0,0,0);
    scalar = 0.5*(MRIgetVoxVal(resVal,index,0,0,0) + scalar);

    MRIsetVoxVal(resVal,index,0,0,0, scalar);

    vertex = &Surf1->vertices[index];
    if (vertex->border == 1) continue;

    total++;


    if (maxV < scalar) maxV = scalar;
    if (minV > scalar) minV = scalar;
    meanV += scalar;

    absMean += (scalar > 0 ? scalar : -scalar);
  }

  printf("total %d vertices involved in thickness computation \n", total);
  meanV /= (total + 1e-10);
  absMean /= (total + 1e-10);

  std = 0.0;
  for (index = 0; index < Surf1->nvertices; index++)
  {
    if (Surf1->vertices[index].border == 1) continue;
    scalar = MRIgetVoxVal(resVal,index,0,0,0) - meanV;

    std += scalar*scalar;
  }

  std /= (total + 1e-10);

  std = sqrt(std);

  printf("thickness stats: max = %g, min = %g, mean = %g, abs_mean = %g, std = %g\n", maxV, minV, meanV, absMean, std);


  if (out_name)
  {
    if (!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w"))
    {

      /* This function will remove a zero-valued vertices */
      /* Make sense, since default value is considered as zero */
      /* But it will confuse the processing with matlab! */
      /* So I copy the data to the curv field to force every value is
       *  written out
       */
      /* MRIScopyMRI(BaseSurf, AvgVals, framesave, "val");*/
      /* MRISwriteValues(BaseSurf,fname); */
      MRIScopyMRI(Surf1, resVal, framesave, (char*)"curv");
      MRISwriteCurvatureToWFile(Surf1,out_name);

    }
    else if (!strcmp(trgtypestring,"curv"))
    {
      MRIScopyMRI(Surf1, resVal, framesave, (char*)"curv");
      MRISwriteCurvature(Surf1,out_name);
    }
    else
    {
      fprintf(stderr, "ERROR unknown output file format.\n");
    }
  }

  if (middle_name)
  {
    for (index = 0; index < Surf1->nvertices; index++)
    {
      Surf1->vertices[index].x = (Surf1->vertices[index].x + Surf2->vertices[index].x) *0.5;
      Surf1->vertices[index].y = (Surf1->vertices[index].y + Surf2->vertices[index].y) *0.5;
      Surf1->vertices[index].z = (Surf1->vertices[index].z + Surf2->vertices[index].z) *0.5;

    }

    MRISwrite(Surf1, middle_name);
  }

  /* Free memories */
  MRISfree(&Surf1);
  MRISfree(&Surf2);
  MRIfree(&resVal);
  MRIfree(&resVal1);

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
  fprintf(stdout, "USAGE: %s [options] surface1 surface2 output_thickness \n",Progname) ;
  fprintf(stdout, "\n");
  fprintf(stdout, "Options:\n");
  fprintf(stdout, "   -trg_type  %%s output format\n");
  fprintf(stdout, "\n");
  std::cout << getVersion() << std::endl;
  printf("\n");

}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf(
    "This program computes the thickness using white and pial surfaces \n"
    "\n"
    "OPTIONS\n"
    "\n"
    "  -trg_type typestring\n"
    "\n"
    "    Format type string. Can be paint or w (for FreeSurfer paint files)."
    "  -out filename\n"
    "  -tflag\n"
    "\n"
    "    Using vertex-to-face distance to compute thickness."
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "trg_type"))
  {
    trgtypestring = argv[2];
    trgtype = string_to_type(trgtypestring);
    nargs = 1 ;
  }
  else if (!stricmp(option, "debug"))
  {
    debugflag = 1;
    debugvtx = atoi(argv[2]);
    nargs = 1;
  }
  else if (!stricmp(option, "face"))
  {
    tflag = 1;
    printf("Using vertex-to-face distance to compute thickness\n");
  }
  else if (!stricmp(option, "middle"))
  {
    middle_name = argv[2];
    nargs  = 1;
    printf("Output a middle surface to file %s\n", middle_name);
  }
  else
  {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    print_help() ;
    exit(1) ;
  }

  return(nargs) ;
}

MRI *ComputeThickness(MRI_SURFACE *Mesh1, MRI_SURFACE *Mesh2, MRI *mri_res)
{
  int index, k, facenumber;
  double value, distance;
  VERTEX vertex;

  ANNpointArray pa = annAllocPts(Mesh2->nvertices, 3);

  for (index = 0; index < Mesh2->nvertices; index++)
  {
    pa[index][0] = Mesh2->vertices[index].x;
    pa[index][1] = Mesh2->vertices[index].y;
    pa[index][2] = Mesh2->vertices[index].z;
  }

  // construct and initialize the tree
  ANNkd_tree *annkdTree = new ANNkd_tree(pa, Mesh2->nvertices, 3);
  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];
  //  ANNpoint query_pt = annAllocPt(3);
  ANNpointArray QueryPt;

  if (mri_res == NULL)
  {
    mri_res =  MRIallocSequence(Mesh1->nvertices, 1, 1,MRI_FLOAT,1);
    if (mri_res ==NULL)
    {
      printf("ERROR: compute thickness: could not alloc\n");
      exit(0);
    }
  }

  QueryPt = annAllocPts(1,3);

  for (index = 0; index < Mesh1->nvertices; index++)
  {
    if (Mesh1->vertices[index].border == 1) continue;

    QueryPt[0][0] = Mesh1->vertices[index].x;
    QueryPt[0][1] = Mesh1->vertices[index].y;
    QueryPt[0][2] = Mesh1->vertices[index].z;

    vertex.x = QueryPt[0][0];
    vertex.y = QueryPt[0][1];
    vertex.z = QueryPt[0][2];

    annkdTree->annkSearch( // search
      QueryPt[0],       // query point
      1,   // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,  // distance (returned)
      0);   // error bound

    if (tflag == 0)
    { /* compute vertex-to-vertex distance */
      value = QueryPt[0][0] - Mesh2->vertices[annIndex[0]].x;
      distance = value*value;
      value = QueryPt[0][1] - Mesh2->vertices[annIndex[0]].y;
      distance += value*value;
      value = QueryPt[0][2] - Mesh2->vertices[annIndex[0]].z;
      distance += value*value;
      distance = sqrt(distance);
    }
    else
    { /* compute vertex-to-face distance */
      distance = 1000.0;
      for (k=0; k < Mesh2->vertices[annIndex[0]].num; k++)
      {
        facenumber =  Mesh2->vertices[annIndex[0]].f[k]; /* index of the k-th face */
        if (facenumber < 0 || facenumber >= Mesh2->nfaces) continue;
        value = v_to_f_distance(&vertex, Mesh2, facenumber, 0);

        /* if(value < 0)
           printf("value at vertex %d is %g\n", index, value); */

        if (distance > value) distance = value;
      }
      if (distance <0) distance = 0;
      distance = sqrt(distance);

    }

    if (distance > max_thickness) distance = max_thickness;

    MRIsetVoxVal(mri_res,index, 0, 0, 0, (float)distance);
  }

  if (annkdTree) delete annkdTree;
  if (annIndex) delete annIndex;
  if (annDist) delete annDist;
  if (QueryPt) delete QueryPt;

  return (mri_res);
}

double v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug)
{
  /* Compute the distance from point P0 to a face of the surface mesh */

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


  /* return (sqrt(a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f)); */
  /* The square-root will be taken later to save time */
  return (a*s*s + 2*b*s*t + c*t*t + 2*d*s + 2*e*t + f);

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
