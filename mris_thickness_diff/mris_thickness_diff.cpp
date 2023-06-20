/**
 * @brief  computes the difference of two surface data sets
 *
 * Given two surfaces with two functions defined on each of them,
 * compute the difference of the function values by first finding the closest
 * vertex on surface2 to a vertex on surface1 and then take the
 * difference of the function-values at the pair of vertices.
 * This function aims to compute thickness difference without mapping to the
 * atlas. That is, to separate true thickness difference from noise caused
 * by sphere-morphing.
 */
/*
 * Original Author: Xaio Han
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
#include "annotation.h"
#include "icosahedron.h"
#include "version.h"

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

static void jacobi(float **a, int n, float *d, float** v,int * nrot);

static char *log_fname = NULL ;
static  char  *subject_name = NULL ;


int main(int argc, char *argv[]) ;

int framesave = 0;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

static char *srctypestring = NULL;
static int srctype = MRI_VOLUME_TYPE_UNKNOWN;
static char *trgtypestring = NULL;
static int trgtype = MRI_VOLUME_TYPE_UNKNOWN;

static char *out_name = NULL;
static char *out_resampled_name = NULL;
static char *maplike_fname = NULL; /*Generate maps to indicate 
                                     thickness difference */
static char *mapout_fname = NULL; /* must be used together with above */

static int debugflag = 0;
static int debugvtx = 0;
static int abs_flag = 0;
static int percentage = 0;
static int register_flag = 0;

static int annotation_flag = 0;
static char *annotation_fname = NULL;
static char *annotation_name = NULL;

static char *label_name = NULL ;

static int compute_distance = 0;

static int nSmoothSteps = 0;

const char *Progname ;

MRI *ComputeDifferenceNew(MRI_SURFACE *Mesh1,
                          MRI *mri_data1,
                          MRI_SURFACE *Mesh2,
                          MRI *mri_data2,
                          MRI *mri_res);
double transformS(double3d *V1a,
                  double3d *V2a,
                  int N,
                  double TR[3][3],
                  double shift[3]);
void FindClosest(MRI_SURFACE *TrueMesh,
                 MHT *srcHash,
                 MRI_SURFACE *EstMesh,
                 double3d *closest);
void register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2);

/* the following two are used when applying a lta transform to the surface */
static MRI *mri = 0;
static MRI *mri_dst = 0;
static int invert = 0 ;
static char *xform_fname = NULL;

double v_to_f_distance(VERTEX *P0,
                       MRI_SURFACE *mri_surf,
                       int face_number,
                       int debug,
                       double *sopt,
                       double *topt);


int main(int argc, char *argv[])
{
  char **av, *surf1_name, *surf2_name;
  char *data1_name, *data2_name;

  int nargs, ac;
  int total, index, label, x, y, z;

  double scalar, std, maxV, minV, meanV, absMean, tmpval;
  int annotation;

  MRI *SrcVal1, *SrcVal2, *resVal;

  MRI_SURFACE *Surf1, *Surf2;
  FILE *log_fp;
  FACE *face;
  VERTEX *vertex;
  double cx, cy, cz;
  double vx, vy, vz;

  MRI *mri_map = NULL;
  int fno, vno0, vno1, vno2;

  LTA          *lta = 0;
  int          transform_type;

  label = 0;
  annotation = 0;

  nargs = handleVersionOption(argc, argv, "mris_thickness_diff");
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

  /* command line: <surf1> <datafile 1> <surf2> <datafile 2> */

  if (argc != 5)
  {
    printf("Incorrect number of arguments, argc = %d\n", argc);
    usage_exit();
  }

  surf1_name = argv[1];
  data1_name = argv[2];
  surf2_name = argv[3];
  data2_name = argv[4];

  if (srctypestring == NULL || (out_name != NULL && trgtypestring == NULL))
  {
    printf("Please specify input and output data type!\n");
    usage_exit();
  }

  printf("Reading first surface file %s\n",surf1_name);
  Surf1 = MRISread(surf1_name);
  if (!Surf1)
    ErrorExit
    (ERROR_NOFILE, "%s:could not read surface %s", Progname, surf1_name);

  printf("Surface1 has %d vertices\n", Surf1->nvertices);

  printf("Reading in first data file %s\n",data1_name);
  /* Read in the first data file */
  if (!strcmp(srctypestring,"curv"))
  { /* curvature file */
    if (MRISreadCurvatureFile(Surf1, data1_name) != 0)
    {
      printf("ERROR: reading curvature file\n");
      exit(1);
    }
    SrcVal1 = MRIcopyMRIS(NULL, Surf1, 0, "curv");
  }
  else if (!strcmp(srctypestring,"paint") || !strcmp(srctypestring,"w"))
  {
    MRISreadValues(Surf1,data1_name);
    SrcVal1 = MRIcopyMRIS(NULL, Surf1, 0, "val");
  }
  else
  {
    SrcVal1=MRIread(data1_name);
    if (NULL==SrcVal1) 
    {
      printf("ERROR: reading file %s as surface file\n",data1_name);
      exit(1);
    }
  }

  if (SrcVal1 == NULL)
  {
    fprintf(stderr, "ERROR loading data values from %s\n", data1_name);
  }

  if (nSmoothSteps > 0)
  {
    printf("Smooth input data 1 by %d steps\n", nSmoothSteps);
    MRISsmoothMRI(Surf1, SrcVal1, nSmoothSteps, NULL,SrcVal1);

  }

  printf("Reading second surface file %s\n",surf2_name);
  Surf2 = MRISread(surf2_name);
  if (!Surf2)
    ErrorExit
    (ERROR_NOFILE, "%s:could not read surface %s", Progname, surf2_name);

  printf("Surface2 has %d vertices\n", Surf2->nvertices);

  printf("Reading in second data file %s\n",data2_name);
  /* Read in the second data file */
  /* only two data types are supported */
  if (!strcmp(srctypestring,"curv"))
  { /* curvature file */
    if (MRISreadCurvatureFile(Surf2, data2_name) != 0)
    {
      printf("ERROR: reading curvature file\n");
      exit(1);
    }
    SrcVal2 = MRIcopyMRIS(NULL, Surf2, 0, "curv");
  }
  else if (!strcmp(srctypestring,"paint") || !strcmp(srctypestring,"w"))
  {
    MRISreadValues(Surf2,data2_name);
    SrcVal2 = MRIcopyMRIS(NULL, Surf2, 0, "val");
  }
  else
  {
    SrcVal2=MRIread(data2_name);
    if (NULL==SrcVal2) 
    {
      printf("ERROR: reading file %s as surface file\n",data2_name);
      exit(1);
    }
  }

  if (SrcVal2 == NULL)
  {
    fprintf(stderr, "ERROR loading data values from %s\n", data2_name);
  }

  if (nSmoothSteps > 0)
  {
    printf("Smooth input data 2 by %d steps\n", nSmoothSteps);
    MRISsmoothMRI(Surf2, SrcVal2, nSmoothSteps, NULL, SrcVal2);
  }

  /* added for allowing using faces; maybe unnecessary */
  for (fno = 0 ; fno < Surf1->nfaces ; fno++)
  {
    face = &Surf1->faces[fno] ;
    vno0 = face->v[0] ;
    vno1 = face->v[1] ;
    vno2 = face->v[2] ;
    mrisSetVertexFaceIndex(Surf1, vno0, fno) ;
    mrisSetVertexFaceIndex(Surf1, vno1, fno) ;
    mrisSetVertexFaceIndex(Surf1, vno2, fno) ;
  }
  for (fno = 0 ; fno < Surf2->nfaces ; fno++)
  {
    face = &Surf2->faces[fno] ;
    vno0 = face->v[0] ;
    vno1 = face->v[1] ;
    vno2 = face->v[2] ;
    mrisSetVertexFaceIndex(Surf2, vno0, fno) ;
    mrisSetVertexFaceIndex(Surf2, vno1, fno) ;
    mrisSetVertexFaceIndex(Surf2, vno2, fno) ;
  }

  if (debugflag)
  {
    printf("Data%d at vertex %d has value %g\n",
           1, debugvtx, MRIFseq_vox(SrcVal1, debugvtx, 0, 0, 0));
  }

  if (register_flag && (xform_fname != NULL))
  {
    printf("A LTA is specified, so only apply the LTA, "
           "and no ICP will be computed\n");
    register_flag = 0;
  }

  if (register_flag)
  {
    printf("Register surface 2 to surface 1 (rigid mapping using ICP)\n");
    register2to1(Surf1, Surf2);
  }


  if (xform_fname != NULL)
  {
    printf("Apply the given LTA xform to align input "
           "surface 1 to surface 2 \n ...");
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
          fprintf
          (stderr,
           "ERROR: fslmat does not have information on the "
           "src and dst volumes\n");
          fprintf
          (stderr,
           "ERROR: you must give options '-src' and '-dst' "
           "to specify the src and dst volume infos\n");
        }

        LTAmodifySrcDstGeom(lta, mri, mri_dst); // add src and dst information
        LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      }

      if (lta->xforms[0].src.valid == 0)
      {
        if (mri == 0)
        {
          fprintf(stderr,
                  "The transform does not have the valid src volume info.\n");
          fprintf(stderr,
                  "Either you give src volume info by option --src or\n");
          fprintf(stderr,
                  "make the transform to have the valid src info.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }
        else
        {
          LTAmodifySrcDstGeom(lta, mri, NULL); // add src information
          //      getVolGeom(mri, &lt->src);
        }
      }
      if (lta->xforms[0].dst.valid == 0)
      {
        if (mri_dst == 0)
        {
          fprintf(stderr,
                  "The transform does not have the valid dst volume info.\n");
          fprintf(stderr,
                  "Either you give src volume info by option --dst or\n");
          fprintf(stderr,
                  "make the transform to have the valid dst info.\n");
          fprintf(stderr,
                  "If the dst was average_305, then you can set\n");
          fprintf(stderr,
                  "environmental variable USE_AVERAGE305 true\n");
          fprintf(stderr,
                  "without giving the dst volume for RAS-to-RAS transform.\n");
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
      ErrorExit
      (ERROR_BADPARM,
       "transform is not of MNI, nor Register.dat, nor FSLMAT type");
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
        fprintf(stderr,
                "WARNING:*************************************************\n");
        fprintf(stderr,
                "WARNING:dst volume infor is invalid. "
                " Most likely produce wrong inverse.\n");
        fprintf(stderr,
                "WARNING:*************************************************\n");
      }
      //copyVolGeom(&lt->dst, &vgtmp);
      vgtmp = lt->dst;
      //copyVolGeom(&lt->src, &lt->dst);
      lt->dst = lt->src;
      //copyVolGeom(&vgtmp, &lt->src);
      lt->src = vgtmp;
    }

    /* save the original cooridnates into orig */
    MRISsaveVertexPositions(Surf1, ORIGINAL_VERTICES);
    {
      TRANSFORM transform ;

      transform.type = lta->type ;
      transform.xform = (void *)lta ;
      MRIStransform(Surf1, mri, &transform, mri_dst) ;
    }

    if (mri)
      MRIfree(&mri);
    if (mri_dst)
      MRIfree(&mri_dst);

    if (lta)
      LTAfree(&lta);
  }   /* if (xform_fname != NULL) */


  if (annotation_flag)
  {
    printf("read annotation file for surface #1\n");
    MRISreadAnnotation(Surf1, annotation_fname) ;
    printf("mask out vertices not have annotation %d\n", annotation);
    //CTABindexToAnnotation(Surf1->ct, annotation_index, &annotation);
    CTABfindName(Surf1->ct, annotation_name, &annotation);
    for (index = 0; index < Surf1->nvertices; index++)
    {
      vertex = &Surf1->vertices[index];
      if (vertex->annotation != annotation)
        vertex->border = 1;
    }
  }
	
  if (label_name)
  {
    LABEL  *area ;
    area = LabelRead(NULL, label_name) ;
    if (!area)
      ErrorExit(ERROR_NOFILE, "%s: could not read label file %s\n", label_name) ;
    MRISsetMarks(Surf1, -1) ;
    LabelMark(area, Surf1) ;  // set the vertices in the label to marked=1
    LabelFree(&area) ;
  }
	else MRISsetMarks(Surf1, 1);

  resVal = ComputeDifferenceNew(Surf1, SrcVal1, Surf2, SrcVal2, NULL);

  printf("Compute difference\n");
  maxV = -1000.0;
  minV = 1000.0;
  meanV=0.0;
  absMean = 0.0;
  total = 0;
  for (index = 0; index < Surf1->nvertices; index++)
  {
    vertex = &Surf1->vertices[index];
    if (vertex->border == 1) continue;
    if (vertex->marked != 1) continue;

    total++;
    scalar = MRIgetVoxVal(resVal,index,0,0,0);

    if (maxV < scalar) maxV = scalar;
    if (minV > scalar) minV = scalar;
    meanV += scalar;

    if (scalar < 0) scalar = -scalar;
    absMean += scalar;

    //    if(abs_flag) MRIsetVoxVal(resVal, index, 0, 0, 0, scalar);
    //    absMean += (scalar > 0 ? scalar : -scalar);
  }

  printf("total %d vertices involved in thickness comparison \n", total);
  meanV /= (total + 1e-20);
  absMean /= (total + 1e-20);

  if (maplike_fname && mapout_fname)
  {
    mri_map = MRIread(maplike_fname);
    for (z = 0; z < mri_map->depth; z++)
      for (y=0; y < mri_map->height; y++)
        for (x = 0; x < mri_map->width; x++)
        {
          MRIsetVoxVal(mri_map, x, y, z, 0, 0);
        }
  }

  if (abs_flag)
  {
    printf("Compute stdev of absolute-valued thickness difference\n");
  }

  std = 0.0;
  for (index = 0; index < Surf1->nvertices; index++)
  {
    if (Surf1->vertices[index].border == 1) continue;
    if (Surf1->vertices[index].marked != 1) continue;
    tmpval = MRIgetVoxVal(resVal,index,0,0,0);
    if (abs_flag)
    { //compute stat of absolute-value-thickness-differences
      if (tmpval < 0) tmpval = -tmpval;
      scalar = tmpval - absMean;
    }
    else
      scalar = tmpval - meanV;

    std += scalar*scalar;

    if (maplike_fname && mapout_fname)
    {
      vertex = &Surf1->vertices[index];
      cx = vertex->origx;
      cy = vertex->origy;
      cz = vertex->origz;

      MRIsurfaceRASToVoxel(mri_map, cx, cy, cz, &vx, &vy, &vz);
      /* Nearest neighbor */
      x = (int) (vx + 0.5);
      y = (int)(vy+0.5);
      z = (int)(vz + 0.5);

      if (x < 0 || x >= mri_map->width ||
          y < 0 || y >= mri_map->height ||
          z < 0 || z >= mri_map->depth) continue;

      if (tmpval < 0) tmpval = -tmpval;

      /* index = 69894 & 69916
         if(x==162 && y==97 && z == 124){
         printf("index = %d, diff = %g\n", index, tmpval);
         } */

      if (DEBUG)
      {
        if (tmpval > 1.0) tmpval = 4;
      }
      else
        tmpval *= 10;

      MRIsetVoxVal(mri_map, x, y, z, 0, tmpval);
    }
  }

  std /= (total + 1e-20);

  std = sqrt(std);

  printf("difference stats: max = %g, min = %g, "
         "mean = %g, abs_mean = %g, std = %g\n",
         maxV, minV, meanV, absMean, std);

  if (log_fname)
  {
    char fname[STRLEN] ;

    sprintf(fname, log_fname, label) ;
    printf("logging to %s...\n", fname) ;
    log_fp = fopen(fname, "a+") ;
    if (!log_fp)
      ErrorExit(ERROR_BADFILE, "%s: could not open %s for writing",
                Progname, fname) ;
  }
  else
    log_fp = NULL ;

  if (log_fp)
  {
    fprintf(log_fp, "%s  ", subject_name) ;
    fprintf(log_fp, "%d %8.4f %8.4f %8.4f\n", total, meanV, absMean, std);
    fclose(log_fp) ;
  }

  if (maplike_fname && mapout_fname)
  {
    MRIwrite(mri_map, mapout_fname);
    MRIfree(&mri_map);
  }

  if (out_name)
  {
    MRIScopyMRI(Surf1, resVal, framesave, "curv");
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
      MRISwriteCurvatureToWFile(Surf1,out_name);
    }
    else if (!strcmp(trgtypestring,"curv"))
    {
      MRISwriteCurvature(Surf1, out_name);
    }
    else
    {
      if (MRIwrite(resVal, out_name)) 
      {
        fprintf(stderr,"ERROR: failed MRIwrite file %s\n",out_name);
        exit(1);
      }
    }
  }

  /* Free memories */
  MRISfree(&Surf1);
  MRISfree(&Surf2);
  MRIfree(&resVal);
  MRIfree(&SrcVal1);
  MRIfree(&SrcVal2);

  return 0;
}


static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

/* --------------------------------------------- */
static void print_usage(void)
{
  fprintf(stdout, "USAGE: %s [options] surface1 data1 surface2 data2 \n",
          Progname) ;
  fprintf(stdout, "\n");
  fprintf(stdout, "Options:\n");
  fprintf(stdout, "   -src_type %%s       input surface data format "
          "(curv, paint or w)\n");
  fprintf(stdout, "   -trg_type %%s       output format\n");
  fprintf(stdout, "   -out %%s            output file name\n");
  fprintf(stdout, "   -out_resampled %%s  output resampled thickness\n");
  fprintf(stdout, "   -nsmooth %%d        number of smoothing steps\n");
  fprintf(stdout, "   -register          force a rigid registration of "
          "surf2 to surf1\n");
  fprintf(stdout, "   -xform %%s          apply LTA transform to align input "
          "surf1 to surf2\n");
  fprintf(stdout, "   -invert            reversely apply -xform \n");
  fprintf(stdout, "   -src %%s            source volume for -xform \n");
  fprintf(stdout, "   -dst %%s            target volume for -xform \n");
  fprintf(stdout, "   -abs               compute the std of abs-thickness-diff \n");
  fprintf(stdout, "   -L %%s              log_file name \n");
  fprintf(stdout, "   -S %%s              subject name \n");
  fprintf(stdout, "   --help             more help \n");
  fprintf(stdout, "\n");
  std::cout << getVersion() << std::endl;
  printf("\n");

}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf
  (
    "mris_thickness_diff computes the difference of two surface data\n"
    "sets defined on the two surface mesh. Result = data2 - data1\n"
    "(in the sense of closest vertex).\n\n"

    "mris_thickness_diff uses closest Euclidean distance to \n"
    "define correspondence across the two surfaces. It does \n"
    "not assume that the correspondence is given by the vertex IDs.\n"
    "But if one has two surfaces that are the 'same' but not aligned \n"
    "in space yet, the trick is to just use one of them as both the \n"
    "src and trgt surface, and mris_thickness_diff will then associate \n"
    "data by their corresponding vertex ID. \n"
    "\n"
    "In contrast, mri_sur2surf builds the surface correspondence \n"
    "through the spherical map. This 'nonlinear' way of building\n"
    "surface correspondence is less reliable than using linear \n"
    "registration and building correspondence by surface distance \n"
    "in the longitudinal case. Overall, the key is that mri_surf2surf \n"
    "is better to align surfaces across different subjects while ICP \n"
    "or linear volume registration is better for aligning surfaces \n"
    "from the same subject across time.  This means that \n"
    "mris_thickness_diff is more useful for longitudinal studies.\n"
    "\n"
    "The -register flag to mris_thickness_diff will run ICP internally \n"
    "to align the two input surfaces if they have not been aligned yet. \n"
    "But ICP (iterative-closest-point) has its limited capture range. \n"
    "For example, if one surface is the other one rotated by more than \n"
    "90 degree, ICP may not align them.\n"

    "\n"
    "OPTIONS\n"
    "\n"
    "  -src_type typestring\n"
    "    Format type string. Can be either curv (for FreeSurfer\n"
    "    curvature file), paint or w (for FreeSurfer paint files).\n"
    "\n"
    "  -trg_type typestring\n"
    "    Format type string. Can be paint or w (for FreeSurfer\n"
    "    paint files).\n"
    "\n"
    "  -out filename\n"
    "    Output the difference to the paint file. \n"
    "    Note the # of entries is same as VN of surface1.\n"
    "\n"
    "  -nsmooth #\n"
    "    Perform # of iterations (smoothing) before computing\n"
    "    the difference.\n"
    "\n"
    "  -annot annotationfile annotation_name\n"
    "    Limit thickness comparison to region with\n"
    "    annotation = annotation_name.\n"
    "\n"
    "  -label labelfile\n"
    "    Limit thickness comparison to region inside label.\n"
    "\n"
    "  -s subject_name (to be recorded in logfile).\n"
    "\n"
    "  -register\n"
    "    Perform ICP rigid registration before computing\n"
    "    closet-vertex difference.\n"
    "\n"
    "  -xform xformfilename\n    Apply LTA transform to align "
    "surface1 to surface2\n\n"
    "  -invert\n    Reversely apply -xform \n\n"
    "  -src filename\n    Source volume for -xform \n\n"
    "  -dst filename\n    Target volume for -xform \n\n"
    "\n"
    "EXAMPLE\n\n"
    "mris_thickness_diff \\ \n"
    "  -out $SUBJECTS_DIR/subj_tp1/surf/lh.thickness_diff \\ \n"
    "  -trg_type curv \\ \n"
    "  -register \\ \n"
    "  $SUBJECTS_DIR/subj_tp1/surf/lh.white \\ \n"
    "  $SUBJECTS_DIR/subj_tp1/surf/lh.thickness \\ \n"
    "  $SUBJECTS_DIR/subj_tp2/surf/lh.white \\ \n"
    "  $SUBJECTS_DIR/subj_tp2/surf/lh.thickness\n\n"
    "tksurfer subj_tp1 lh inflated \\ \n"
    "  -overlay $SUBJECTS_DIR/subj_tp1/surf/lh.thickness_diff\n"
    "Select menu View->Configure->Overlap to adjust threshold.\n\n"
  );

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
static int get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "src_type"))
  {
    srctypestring = argv[2];
    srctype = string_to_type(srctypestring);
    nargs = 1 ;
  }
  else if (!stricmp(option, "register"))
  {
    register_flag = 1;
    printf("Do a rigid alignment of two surfaces before "
           "mapping to each other\n");
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
    fprintf(stderr, "src volume for the given transform "
            "(given by -xform) is %s\n",argv[2]);
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
    fprintf(stderr, "dst volume for the transform "
            "(given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the dst volume...\n");
    mri_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!mri_dst)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else if (!stricmp(option, "out") ||
           !stricmp(option, "out_file") ||
           !stricmp(option, "out_name"))
  {
    out_name = argv[2];
    printf("Output differences to file %s\n", out_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "out_resampled"))
  {
    out_resampled_name = argv[2];
    printf("Output TP2 thickness resampled to TP1 as file %s\n", out_resampled_name);
    nargs = 1 ;
  }
  else if (!stricmp(option, "nsmooth"))
  {
    nSmoothSteps = atoi(argv[2]);
    nargs = 1;
    printf("Perform %d steps of smoothing of input data\n", nSmoothSteps);
  }
  else if (!stricmp(option, "trg_type"))
  {
    trgtypestring = argv[2];
    trgtype = string_to_type(trgtypestring);
    nargs = 1 ;
  }
  else if (!stricmp(option, "annot") ||
           !stricmp(option, "annotation"))
  {
    annotation_fname = argv[2];
    annotation_flag = 1;
    annotation_name = argv[3];
    nargs = 2;
    printf("only compute thickness stats for region with label %s\n",
           annotation_name);
  }
  else if (!stricmp(option, "label") )
  {
      label_name = argv[2] ;
      nargs = 1 ;
      printf("limiting computations to label %s.\n", label_name) ;
  }
  else if (!stricmp(option, "distance"))
  {
    compute_distance = 1;
  }
  else if (!stricmp(option, "abs"))
  {
    abs_flag = 1;
  }
  else if (!stricmp(option, "L") ||
           !stricmp(option, "Log")
          )
  {
    log_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "logging results to %s\n", log_fname) ;
  }
  else if (!stricmp(option, "S") ||
           !stricmp(option, "subj")
          )
  {
    subject_name = argv[2];
    nargs = 1;
  }
  else if (!stricmp(option, "debug"))
  {
    debugflag = 1;
    debugvtx = atoi(argv[2]);
    nargs = 1;
  }
  else if (!stricmp(option, "percentage"))
  {
    percentage = 1;
    printf("Compute percentage thickness-difference\n");
  }
  else if (!stricmp(option, "map_like"))
  {
    maplike_fname = argv[2];
    printf("Create a volume indicating thickness diff\n");
    nargs = 1;
  }
  else if (!stricmp(option, "map_out"))
  {
    mapout_fname = argv[2];
    printf("Output thickness volume map to %s\n", mapout_fname);
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


double v_to_f_distance(VERTEX *P0,
                       MRI_SURFACE *mri_surf,
                       int face_number,
                       int debug,
                       double *sopt,
                       double *topt)
{
  /* Compute the distance from point P0 to a face of the surface mesh */
  /* (sopt, topt) determines the closest point inside the triangle from P0 */

  double a, b, c, d, e, f, det, s, t, invDet;
  double numer, denom, tmp0, tmp1;
  VERTEX *V1, *V2, *V3;
  FACE *face;

  struct { float x,y,z; } E0, E1, D;

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
        if (debug)
          printf("region 4,  d= %g, e = %g, s =%g, t =%g\n", d, e, s, t);
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


MRI *ComputeDifferenceNew(MRI_SURFACE *Mesh1,
                          MRI *mri_data1,
                          MRI_SURFACE *Mesh2,
                          MRI *mri_data2,
                          MRI *mri_res)
{
  /* This one will do a more accurate interpolation */
  int index, k, facenumber, closestface, nnindex;
  VERTEX *vertex;
  double sumcurv, sumweight, weight;
  VERTEX *V1, *V2, *V3;
  FACE *face;
  double value, distance, tmps, tmpt;
  double total_distance, tmp_distance;
  double max_distance, std_dist;
  float dmin;
  MHT *SrcHash;
	MRI *mri_resampled = MRIclone(mri_data1, NULL);

  SrcHash = MHTcreateVertexTable_Resolution(Mesh2, CURRENT_VERTICES,16);

  max_distance = 0;

  if (mri_res == NULL)
    mri_res = MRIclone(mri_data1, NULL);


  total_distance = 0.0;
  std_dist = 0;

  for (index = 0; index < Mesh1->nvertices; index++)
  {
    if (Mesh1->vertices[index].border == 1) continue;
    if (Mesh1->vertices[index].marked != 1) continue;
    vertex = &Mesh1->vertices[index];
    nnindex = MHTfindClosestVertexNo2(SrcHash,Mesh2,Mesh1,vertex,&dmin);

    /* nnIndex gives the closest vertex on
       mris_template to the vertex #index of mris 
       sanity-check it: */
    if (nnindex >= Mesh2->nvertices)
    {
      printf("\nERROR: ComputeDifferenceNew: MHTfindClosestVertexNo "
             "vidx %d returned nnindex=%d, which is >= nvertices=%d\n", 
             index,nnindex, Mesh2->nvertices);
      exit(1);
    }
    if (nnindex < 0)
    {
      printf("\nERROR: ComputeDifferenceNew: MHTfindClosestVertexNo "
             "vidx %d returned nnindex=%d, which is invalid\n", index, nnindex);
      exit(1);
    }

    /* Now need to find the closest face in order
       to perform linear interpolation,
       but first, a sanity-check: */
    if (Mesh2->vertices_topology[nnindex].num <= 0) 
    {
      printf("\nERROR: ComputeDifferenceNew: \n"
             "(Mesh2->vertices[nnindex].num=%d <= 0)\n", 
             Mesh2->vertices_topology[nnindex].num);
      exit(1);
    }
    distance = 1e30;
    closestface = 0;
    for (k=0; k < Mesh2->vertices_topology[nnindex].num; k++)
    {

      facenumber =  Mesh2->vertices_topology[nnindex].f[k]; /* index of the k-th face */
      if (facenumber < 0 || facenumber >= Mesh2->nfaces) continue;
      value = v_to_f_distance(vertex, Mesh2, facenumber, 0, &tmps, &tmpt);
      /* Here it would be better to compute the correct distance together
			   with the closest point and then use barycentric coordinates below (mr) */

      if (distance > value)
      {
        distance = value;
        closestface = facenumber;
      }
    }

    if (compute_distance)
    {
      tmp_distance =  sqrt(distance);
      if (max_distance < tmp_distance) max_distance = tmp_distance;
      total_distance += (tmp_distance);
      std_dist += distance;
    }

    face = &Mesh2->faces[closestface];
    V1 = &Mesh2->vertices[face->v[0]];
    V2 = &Mesh2->vertices[face->v[1]];
    V3 = &Mesh2->vertices[face->v[2]];
    sumcurv = 0.0;
    sumweight = 0.0;
    weight = 1.0/(1e-20 +
                  (V1->x - vertex->x)*(V1->x - vertex->x) +
                  (V1->y - vertex->y)*(V1->y - vertex->y) +
                  (V1->z - vertex->z)*(V1->z - vertex->z));
    sumcurv += weight*MRIgetVoxVal(mri_data2,face->v[0],0,0,0);
    sumweight += weight;

    weight = 1.0/(1e-20 +
                  (V2->x - vertex->x)*(V2->x - vertex->x) +
                  (V2->y - vertex->y)*(V2->y - vertex->y) +
                  (V2->z - vertex->z)*(V2->z - vertex->z));
    sumcurv += weight*MRIgetVoxVal(mri_data2,face->v[1],0,0,0);
    sumweight += weight;

    weight = 1.0/(1e-20 +
                  (V3->x - vertex->x)*(V3->x - vertex->x) +
                  (V3->y - vertex->y)*(V3->y - vertex->y) +
                  (V3->z - vertex->z)*(V3->z - vertex->z));
    sumcurv += weight*MRIgetVoxVal(mri_data2,face->v[2],0,0,0);
    sumweight += weight;

    sumcurv /= (sumweight + 1e-30);
		
    MRIsetVoxVal(mri_resampled,index, 0, 0, 0, sumcurv);
		
		
    value =  sumcurv -
             MRIgetVoxVal(mri_data1,index,0,0,0);

    if (percentage)
    {
      value /= 0.5*(sumcurv + MRIgetVoxVal(mri_data1,index,0,0,0) + 1e-15);
    }

    MRIsetVoxVal(mri_res,index, 0, 0, 0, value);
  }

  if (compute_distance)
  {
    std_dist /= (Mesh1->nvertices + 1e-30);
    total_distance /= (Mesh1->nvertices + 1e-30);
    std_dist = sqrt(std_dist - total_distance* total_distance);
    printf("Max, average, and std of vertex-to-vertex distances of "
           "the two surfaces are %g, %g and %g\n",
           max_distance, total_distance, std_dist);
  }
	
	if (out_resampled_name) // output thickness TP2 resampled to TP1
	{
	  MRI_SURFACE *Surfresampled = MRISclone(Mesh1);
    MRIScopyMRI(Surfresampled, mri_resampled, framesave, "curv");
    if (!strcmp(trgtypestring,"paint") || !strcmp(trgtypestring,"w"))
    {

      /* This function will remove a zero-valued vertices */
      /* Make sense, since default value is considered as zero */
      /* But it will confuse the processing with matlab! */
      /* So I copy the data to the curv field to force every value is
       *  written out
       */
      MRISwriteCurvatureToWFile(Surfresampled,out_resampled_name);
    }
    else if (!strcmp(trgtypestring,"curv"))
    {
      MRISwriteCurvature(Surfresampled, out_resampled_name);
    }
    else
    {
      if (MRIwrite(mri_resampled, out_resampled_name)) 
      {
        fprintf(stderr,"ERROR: failed MRIwrite file %s\n",out_resampled_name);
        exit(1);
      }
    }
	   
	
	}
	
	MRIfree(&mri_resampled);

  return (mri_res);
}


void register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2)
{
  MHT* SrcHash = MHTcreateVertexTable_Resolution(Surf1, CURRENT_VERTICES,16);

  /* This initialization is necessary */
  double TR[3][3];
  TR[0][0] = TR[1][1] = TR[2][2] = 1;
  TR[0][1] = TR[0][2] = TR[1][0] = TR[1][2] = TR[2][0] = TR[2][1] = 0;

  double shift[3];
  shift[0] = 0;
  shift[1] = 0;
  shift[2] = 0;

  double3d* mesh2avtx  = (double3d*) malloc(Surf2->nvertices*sizeof(double3d));
  double3d* closevtx2a = (double3d*) malloc(Surf2->nvertices*sizeof(double3d));

  /* Initialize */
  int index;
  for (index = 0; index < Surf2->nvertices; index++)
  {
    mesh2avtx[index].x = Surf2->vertices[index].x;
    mesh2avtx[index].y = Surf2->vertices[index].y;
    mesh2avtx[index].z = Surf2->vertices[index].z;
  }

  double error_old = 1000.0;
  double error_new = 900.0;

  int iter = 0;
  while ((error_old - error_new) > 0.0001 && iter < 50)
  {
    error_old = error_new;
    
    /* For each vertex in Surf2, find its closest vertex in Surf1,
     * and record the coordinates in closevtx2a
     */
    FindClosest(Surf1, SrcHash, Surf2, closevtx2a);

    // Find the rigid transformation
    error_new = transformS(mesh2avtx, closevtx2a, Surf2->nvertices, TR, shift);

    int k;
    for (k = 0; k < Surf2->nvertices; k++)
    {
      MRISsetXYZ(Surf2,k,
        mesh2avtx[k].x,
        mesh2avtx[k].y,
        mesh2avtx[k].z);
    }

    iter ++;
    if (debugflag)  printf(" iteration %d, error = %15.4f\n", iter, error_new);
  }

  MHTfree(&SrcHash);

  free(mesh2avtx);
  free(closevtx2a);

  return;
}


void FindClosest(MRI_SURFACE *TrueMesh,
                 MHT *SrcHash,
                 MRI_SURFACE *EstMesh,
                 double3d *closest)
{
  int index;
  float dmin;
  int annIndex;
  VERTEX *v;

  for (index =0; index < EstMesh->nvertices; index++)
  {
    // this is a duplicate, lame, but....ugh, to get in and out of libraries...
    v = &(EstMesh->vertices[index]);

    annIndex = MHTfindClosestVertexNo2(SrcHash, TrueMesh,EstMesh,v,&dmin);

    if (annIndex < 0)
    {
      printf("\nERROR:  FindClosest: MHTfindClosestVertexNo "
             "returned annIndex=%d, which is invalid!\n", annIndex);
      exit(1);
    }

    closest[index].x = TrueMesh->vertices[annIndex].x;
    closest[index].y = TrueMesh->vertices[annIndex].y;
    closest[index].z = TrueMesh->vertices[annIndex].z;
  }

  return;
}


double transformS(double3d *V1a,
                  double3d *V2a,
                  int N,
                  double TR[3][3],
                  double shift[3])
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
  double temp[3][3];
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
    scale1+=(V1a[k].x * V1a[k].x + V1a[k].y * V1a[k].y + V1a[k].z * V1a[k].z);
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
     printf("\t %15.9f %15.9f %15.9f %15.9f\n",
     M[1][1],M[1][2],M[1][3], M[1][4]);
     printf("\t %15.9f %15.9f %15.9f %15.9f\n",
     M[2][1],M[2][2],M[2][3], M[2][4]);
     printf("\t %15.9f %15.9f %15.9f %15.9f\n",
     M[3][1],M[3][2],M[3][3], M[3][4]);
     printf("\t %15.9f %15.9f %15.9f %15.9f\n",
     M[4][1],M[4][2],M[4][3], M[4][4]);*/

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
  /* R is not symmetric, because it's a rotation
     around an arbitrary axis, not just the origin */
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
    error += ((V1a[k].x-V2a[k].x)*(V1a[k].x-V2a[k].x) +
              (V1a[k].y-V2a[k].y)*(V1a[k].y-V2a[k].y) +
              (V1a[k].z-V2a[k].z)*(V1a[k].z-V2a[k].z));

    V1a[k].x += centroid2a.x;
    V1a[k].y += centroid2a.y;
    V1a[k].z += centroid2a.z;
  }

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
  shift[0]=scale1*(R[0][0]*(temp[0][0]-centroid1a.x)+
                   R[1][0]*(temp[0][1]-centroid1a.y)+
                   R[2][0]*(temp[0][2]-centroid1a.z))+centroid2a.x;
  shift[1]=scale1*(R[0][1]*(temp[0][0]-centroid1a.x)+
                   R[1][1]*(temp[0][1]-centroid1a.y)+
                   R[2][1]*(temp[0][2]-centroid1a.z))+centroid2a.y;
  shift[2]=scale1*(R[0][2]*(temp[0][0]-centroid1a.x)+
                   R[1][2]*(temp[0][1]-centroid1a.y)+
                   R[2][2]*(temp[0][2]-centroid1a.z))+centroid2a.z;

  /* The following is just the shift at the current step.
   * Note the first point data set is constantly
   updated/transformed every iteration
   */
  /* shift[0][0]=scale1*(R[0][0]*(-centroid1a.x)+
     R[1][0]*(-centroid1a.y)+R[2][0]*(-centroid1a.z))+centroid2a.x;
     shift[0][1]=scale1*(R[0][1]*(-centroid1a.x)+
     R[1][1]*(-centroid1a.y)+R[2][1]*(-centroid1a.z))+centroid2a.y;
     shift[0][2]=scale1*(R[0][2]*(-centroid1a.x)+
     R[1][2]*(-centroid1a.y)+R[2][2]*(-centroid1a.z))+centroid2a.z; */

  return(sqrt(error)/(N+1e-10));
}


static void jacobi(float **a, int n,float *d, float **v, int *nrot)
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


#undef ROTATE
