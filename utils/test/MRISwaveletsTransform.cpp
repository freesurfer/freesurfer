/**
 * @brief wavelet utils
 *
 */
/*
 * Original Author: Peng Yu
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

#include "ANN.h"

extern "C"
{
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <ctype.h>
#include "mri.h"
#include "mrisurf.h"
#include "icosahedron.h"
#include "const.h"
#include "diag.h"
#include "error.h"
#include "macros.h"
#include "proto.h"
#include "timer.h"
#include "mrinorm.h"
#include "cma.h"
#include "version.h"
#include "error.h"
#include "matrix.h"
}

#define VERTEX_EDGE(vec, v0, v1)   VECTOR_LOAD(vec,v1->x-v0->x,v1->y-v0->y, v1->z-v0->z)
typedef struct _double_3d
{
  double x;
  double y;
  double z;
}
double3d ;

#define TINY 1.0e-20;
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
 a[k][l]=h+s*(g-h*tau);

int             main(int argc, char *argv[]) ;
static int      get_option(int argc, char *argv[]) ;
const char            *Progname ;
static MRI_SURFACE *center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst);
static MRI_SURFACE *sample_origposition(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
static MRI_SURFACE *sample_origcurvature(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst) ;
static double   v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug) ;
static double   brain_volume(MRI_SURFACE *mris);
//static int      sort(double **array, int order, int number, int total_number);
/* ICP functions */
#if 0
#if 0
static void     register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2);
static void     FindClosest(MRI_SURFACE *TrueMesh, ANNkd_tree *annkdTree, MRI_SURFACE *EstMesh, double3d *closest);
static double   transformS(double3d *V1a, double3d *V2a, int N, double TR[3][3], double shift[3]);
static void     jacobi(float **a, int n, float *d, float** v,int * nrot);
static int      DEBUG=0; //warning
#else
static void     register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2, double T[5][5]);
static void     FindClosest(MRI_SURFACE *TrueMesh, ANNkd_tree *annkdTree, MRI_SURFACE *EstMesh, double3d *closest);
static double   affineS(double3d *V1, double3d *V2, int N, double A[5][5]);
static void     ludcmp(double** a,int n,int* indx,double* d);
static void     lubksb(double** a,int n,int* indx,double* b);
#endif
#endif

static MRI_SURFACE  *wavelet_analysis_curv(MRI_SURFACE *mris_out, int order) ;
static MRI_SURFACE  *wavelet_analysis_vec(MRI_SURFACE *mris_out, int order);
static MRI_SURFACE  *wavelet_synthesis_curv(MRI_SURFACE *mris_out, int order) ;
static MRI_SURFACE  *wavelet_synthesis_vec(MRI_SURFACE *mris_out, int order);

static int      ANALYSIS=0;
static int      SYNTHESIS=0;
static int      COMPARE=0;
static int      CURV=0;
static int      RADIUS=0;
static char     *outs_fname;
static char     *inc_fname;
static int      SAMPLE_OUT=0;
static float    threshold=0;
static float    shrink = 0;

int
main(int argc, char *argv[])
{
  int           nargs, msec, order, i, number;
  Timer then ;
  MRIS          *mris_in, *mris_out;
  double        volume; // T[5][5] the transformation matrx (using index 1-4)

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    ErrorExit(ERROR_BADPARM,
              "usage: %s <input sphere> <input surface> <int which> <output wavelets> ", Progname);

  then.reset() ;

  //order = atoi (argv[3]);
  order = 7;
  //fprintf(stdout, "Set %s as the finest scale level\n", argv[3]);
  if (order > 8)
    ErrorExit(ERROR_BADPARM, "the highest order is 7\n");

  if (atoi(argv[3])  == 0)
  {
    ANALYSIS=1;
  }
  else if (atoi(argv[3]) == 1)
  {
    ANALYSIS = 1;
    CURV=1;
  }
  else if (atoi(argv[3]) == 2)
  {
    SYNTHESIS = 1;
    CURV = 0;
  }
  else if (atoi(argv[3]) == 3)
  {
    SYNTHESIS = 1;
    CURV = 1;
  }
  else
  {
    fprintf(stdout, "Set mode to 1\n");
    ANALYSIS=1;
  }

  /*Spherical Wavelet Analysis*/

  if ( ANALYSIS && !CURV && !RADIUS)
  {
    //MRI   *mri_out;
    char  *cp=0, afname[STRLEN], meas[STRLEN], path[STRLEN], sname[STRLEN];

    /* using talairach transformation */
    fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);
    mris_in = MRISread(argv[1]) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    MRISreadOriginalProperties(mris_in, argv[2]) ;
    fprintf(stdout, "Reading original surface from %s\n", argv[2]);
    MRISsaveVertexPositions(mris_in, TMP_VERTICES);
    MRISrestoreVertexPositions(mris_in, ORIGINAL_VERTICES);

    strcpy(afname,argv[2]);
    cp=strchr(afname, '.');
    strcpy(meas,cp-2);
    FileNamePath(afname,path);
    FileNamePath(path,afname);
    cp=strrchr(afname, '/');
    strcpy(sname,cp);
    FileNamePath(afname,path);

    if (0) //comment out this talairach transformation
    {
      MRI          *mri = 0, *mri_dst = 0;
      LTA          *lta=0;
      int          transform_type;

      sprintf(afname, "%s/%s/mri/transforms/talairach.lta",path,sname);
      printf("Apply the given LTA xfrom to align input surface\n ...");
      printf("Reading transform from %s\n", afname);
      // read transform
      transform_type =  TransformFileNameType(afname);

      if (transform_type == MNI_TRANSFORM_TYPE ||
          transform_type == TRANSFORM_ARRAY_TYPE || \
          transform_type == REGISTER_DAT ||
          transform_type == FSLREG_TYPE
         )
      {
        lta = LTAreadEx(afname) ;
        if (!lta)
          ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
                    Progname, afname) ;

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

      MRIStransform(mris_in, mri, lta, mri_dst) ;

      if (mri)
        MRIfree(&mri);
      if (mri_dst)
        MRIfree(&mri_dst);

      if (lta)
        LTAfree(&lta);
    }   /* if (xform_fname != NULL) */

    /* save the transformed original coordinates */
    center_brain(mris_in, mris_in);
    MRISsaveVertexPositions(mris_in, ORIGINAL_VERTICES);
    MRISrestoreVertexPositions(mris_in, TMP_VERTICES);
    /* sample the surface onto ic7 */
    mris_out = ReadIcoByOrder(order, 100);
    sample_origposition(mris_in, mris_out) ;
    //MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    //MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES);
    //if (SAMPLE_OUT)
    //  {MRISwrite(mris_out,outs_fname) ;
    //  fprintf(stdout,"Write sampled surface to %s\n",outs_fname);}
    //MRISsaveVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    //MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;

#if 0
    /* method 1: ICP registration to an average surface*/
    MRISfree(&mris_in);
    sprintf(afname, "%s/average/surf/%s",path,meas);
    mris_in = MRISread(afname) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read registration surface %s",
                Progname, afname) ;
    MRISsaveVertexPositions(mris_in, ORIGINAL_VERTICES);
    fprintf(stdout, "Register %s to surface %s (rigid mapping using ICP)\n", argv[2],afname);
    sprintf(afname,"sphere.reg");
    if (MRISreadVertexPositions(mris_in, afname) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s:could not read sphere.reg from %s", Progname, afname);
    register2to1(mris_in, mris_out, T);  //rigid registration
#endif

    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    center_brain(mris_out, mris_out);
    volume=brain_volume(mris_out);
    //MRISscaleBrain(mris_out, mris_out, cbrt(300000/volume)) ;
    if (SAMPLE_OUT)
    {
      MRISwrite(mris_out,outs_fname) ;
      fprintf(stdout,"Write sampled surface to %s\n",outs_fname);
    }

    volume=brain_volume(mris_out);
    fprintf(stderr, "brain volume becomes %f\n", volume) ;
    MRISsaveVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    MRISupdateSurface(mris_out);

    /* wavelet analysis of origx, origy, origz */
    wavelet_analysis_vec(mris_out,order);

    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
#if 0
    for (i=0;i<mris_out->nvertices;i++)
      if (mris_out->vertices[i].z>6)
        fprintf(stdout, "%d %f %f %f\n", i,mris_out->vertices[i].x, mris_out->vertices[i].y, mris_out->vertices[i].z);
#endif
    fprintf(stdout, "Writing wavelets coefficient of original surface to %s\n", argv[4]);
    MRISwrite(mris_out,argv[4] ) ;
#if 0
    mri_out = MRIallocSequence(mris_out->nvertices,1,1,MRI_FLOAT,3);
    for ( i=0;i<mris_out->nvertices;i++)
    {
      MRIFseq_vox(mri_out, i, 0, 0, 0) = mris_out->vertices[i].x;
      MRIFseq_vox(mri_out, i, 0, 0, 1) = mris_out->vertices[i].y;
      MRIFseq_vox(mri_out, i, 0, 0, 2) = mris_out->vertices[i].z;
    }
    MRIwrite(mri_out,argv[5]);
    MRIfree(&mri_out);
#endif
    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
    MRISfree(&mris_in) ;
    /*End of Analysis*/
  }
  else if ( ANALYSIS && CURV)
  {
    //MRI   *mri_out;

    mris_in = MRISread(argv[1]) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);

    MRISreadCurvatureFile(mris_in, argv[2]) ;
    fprintf(stdout, "Reading input from %s\n", argv[2]);

    mris_out = ReadIcoByOrder(order, 100);
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    sample_origcurvature(mris_in, mris_out) ;
    if (SAMPLE_OUT)
      MRISwriteCurvature(mris_out,outs_fname);

    /* wavelet analysis of origx, origy, origz */
    wavelet_analysis_curv(mris_out,order);

    //fprintf(stdout, "Writing wavelets coefficient of original surface to %s\n", argv[4]);
    MRISwriteCurvature(mris_out,argv[4] ) ;
#if 0
    mri_out = MRIallocSequence(mris_out->nvertices,1,1,MRI_FLOAT,1);
    for ( i=0;i<mris_out->nvertices;i++)
    {
      MRIFseq_vox(mri_out, i, 0, 0, 0) = mris_out->vertices[i].curv;
    }
    MRIwrite(mri_out,argv[5]);
    MRIfree(&mri_out);
#endif
    MRISfree(&mris_in) ;
    /*End of Analysis*/
  }
  else if ( ANALYSIS && RADIUS )
  {
    MRI   *mri_out;
    float rad;

#if 1
    mris_in = MRISread(argv[1]) ;
    if (!mris_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    fprintf(stdout, "Reading input spherical surface from %s\n", argv[1]);

    MRISreadOriginalProperties(mris_in, argv[2]) ;
    mris_out = ReadIcoByOrder(order, 100);

    sample_origposition(mris_in, mris_out) ;

    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    center_brain(mris_out, mris_out);
    volume = brain_volume(mris_out);
    MRISscaleBrain(mris_out, mris_out, cbrt(300000/volume)) ;
#else

    MRIS *mris_high;
    char  *cp, afname[STRLEN], meas[10];

    mris_high = ReadIcoByOrder(order, 100);
    fprintf(stdout, "Reading input surface from %s\n", argv[1]);
    mris_out = MRISread(argv[1]) ;
    if (!mris_out)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, argv[1]) ;
    MRISreadOriginalProperties(mris_out, argv[2]) ;
    sample_origposition(mris_out, mris_high) ;
    MRISfree(&mris_out);
    mris_out = MRISclone(mris_high);
    for (i = 0; i<mris_out->nvertices; i++)
    {
      mris_out->vertices[i].origx=mris_high->vertices[i].origx;
      mris_out->vertices[i].origy=mris_high->vertices[i].origy;
      mris_out->vertices[i].origz=mris_high->vertices[i].origz;
    }

    cp=strchr(argv[2], '.');
    if (strstr(cp-2,"rh"))
      mris_in = MRISread("/space/birn/42/users/pengyu/AD/average7/surf/rh.sphere.reg") ;
    else
      mris_in = MRISread("/space/birn/42/users/pengyu/AD/average7/surf/lh.sphere.reg") ;
    strcpy(meas,cp-2);
    sprintf(afname, "/space/birn/42/users/pengyu/AD/average7/surf/%s",meas);
    //sprintf(afname, "/space/birn/42/users/pengyu/AD/average7/surf/lh.white");
    MRISreadOriginalProperties(mris_in, afname) ;
    sample_origposition(mris_in, mris_high) ;
    MRISfree(&mris_in);
    mris_in = MRISclone(mris_high);
    for (i = 0; i<mris_out->nvertices; i++)
    {
      mris_in->vertices[i].origx=mris_high->vertices[i].origx;
      mris_in->vertices[i].origy=mris_high->vertices[i].origy;
      mris_in->vertices[i].origz=mris_high->vertices[i].origz;
    }
    MRISfree(&mris_high);

    MRISsaveVertexPositions(mris_in, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_in, ORIGINAL_VERTICES) ;
    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;

    center_brain(mris_out, mris_out);
    center_brain(mris_in, mris_in);
    fprintf(stdout, "Register %s to surface %s (rigid mapping using ICP)\n", argv[2],afname);
    register2to1(mris_in, mris_out, T);
    center_brain(mris_out, mris_out);
    if (SAMPLE_OUT)
    {
      MRISwrite(mris_out,outs_fname) ;
      fprintf(stdout,"Write sampled surface to %s\n",outs_fname);
    }
#endif
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    MRISupdateSurface(mris_out);
    volume=brain_volume(mris_out);
    fprintf(stderr, "brain volume becomes %f\n", volume) ;


    for (i=0;i<mris_out->nvertices;i++)
    {
      rad = (mris_out->vertices[i].x*mris_out->vertices[i].x+mris_out->vertices[i].y*mris_out->vertices[i].y+mris_out->vertices[i].z*mris_out->vertices[i].z);
      mris_out->vertices[i].curv = sqrt(rad);
    }

    if (SAMPLE_OUT)
      MRISwriteCurvature(mris_out,outs_fname);

    /* wavelet analysis of radius */
    wavelet_analysis_curv(mris_out,order);

    fprintf(stdout, "Writing wavelets coefficient of original surface to %s\n", argv[4]);
    MRISwriteCurvature(mris_out,argv[4] ) ;
#if 0
    mri_out = MRIallocSequence(mris_out->nvertices,1,1,MRI_FLOAT,1);
    for ( i=0;i<mris_out->nvertices;i++)
    {
      MRIFseq_vox(mri_out, i, 0, 0, 0) = mris_out->vertices[i].curv;
    }
    MRIwrite(mri_out,argv[5]);
#endif
    MRIfree(&mri_out);
    MRISfree(&mris_in) ;
    /*End of Analysis*/
  }
  else if ( SYNTHESIS && !CURV && !RADIUS ) /*Spherical Wavelet Synthesis*/
  {
    mris_out = ReadIcoByOrder(order, 100); //higher order surface
    fprintf(stdout, "Creating a %d order spherical surface\n", order);
#if 0
    MRISreadOriginalProperties(mris_out, argv[1]) ;
#else
    mris_in = MRISread(argv[1]);
    for (i = 0; i<mris_out->nvertices; i++)
    {
      mris_out->vertices[i].origx = mris_in->vertices[i].x ;
      mris_out->vertices[i].origy = mris_in->vertices[i].y ;
      mris_out->vertices[i].origz = mris_in->vertices[i].z ;
    }
    MRISfree(&mris_in);
#endif
    fprintf(stdout, "Reading wavelet coefficients from %s\n", argv[1]);
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_out, 3) ;

    fprintf(stdout, "Recover the surface using %s order coefficients\n",argv[2]);
    if (atoi(argv[2])==-1) number = 0;
    else  number = IcoNVtxsFromOrder(atoi(argv[2]));
    for (i = number; i<mris_out->nvertices; i++)
    {
      mris_out->vertices[i].origx = 0;
      mris_out->vertices[i].origy = 0;
      mris_out->vertices[i].origz = 0;
    }

#if 0 //for SPHARM comparison
    if (COMPARE)
    {
      int count = 0, k;
      double **array;

      mris_in = MRISread(inc_fname);
      array = (double **)calloc(2,sizeof(char *));
      array[0] = (double *)calloc(IcoNVtxsFromOrder(order)*3, sizeof(double)) ;
      array[1] = (double *)calloc(IcoNVtxsFromOrder(order)*3, sizeof(double)) ;
      for (i=0; i<IcoNVtxsFromOrder(order); i++)
      {
        //if (mris_out->vertices[i].origx==0)
        area =  fabs(mris_out->vertices[i].origx-mris_in->vertices[i].x);
        //else area = fabs((mris_out->vertices[i].origx-mris_in->vertices[i].x)/mris_out->vertices[i].origx);
        array[0][i*3]=area;
        array[1][i*3]=i*3;
        //if ( area>threshold )
        //{mris_out->vertices[i].origx = mris_in->vertices[i].x ;
        //fprintf(stdout, "%d %f\n", i, area); count++;}
        //if (mris_out->vertices[i].origy==0)
        area =  fabs(mris_out->vertices[i].origy-mris_in->vertices[i].y);
        //else area = fabs((mris_out->vertices[i].origy-mris_in->vertices[i].y)/mris_out->vertices[i].origy);
        array[0][i*3+1]=area;
        array[1][i*3+1]=i*3+1;
        //if ( area>threshold )
        //{mris_out->vertices[i].origy = mris_in->vertices[i].y ;
        //fprintf(stdout, "%d %f\n", i, area);  count++;}
        //if (mris_out->vertices[i].origz==0)
        area =  fabs(mris_out->vertices[i].origz-mris_in->vertices[i].z);
        //else area = fabs((mris_out->vertices[i].origz-mris_in->vertices[i].z)/mris_out->vertices[i].origz);
        array[0][i*3+2]=area;
        array[1][i*3+2]=i*3+2;
        //if ( area>threshold )
        //{mris_out->vertices[i].origz = mris_in->vertices[i].z ;
        //fprintf(stdout, "%d %f\n", i, area);  count++;}
      }
      sort(array,0,1000,IcoNVtxsFromOrder(order)*3);
      for (count=0; count<threshold; count++)
      {
        i=(int)floor(array[1][count]/3);
        k=(int)array[1][count]-i*3;
        if (k==0)
        {
          mris_out->vertices[i].origx = mris_in->vertices[i].x ;
          fprintf(stdout, "%d %d %f\n", i, k, array[0][count]);
        }
        else if (k==1)
        {
          mris_out->vertices[i].origy = mris_in->vertices[i].y ;
          fprintf(stdout, "%d %d %f\n", i, k, array[0][count]);
        }
        else if (k==2)
        {
          mris_out->vertices[i].origz = mris_in->vertices[i].z ;
          fprintf(stdout, "%d %d %f\n", i, k, array[0][count]);
        }
        else fprintf(stdout, "%d has an error\n", count);
      }
      fprintf(stdout, "%d coefficients used in total\n",count);
      MRISfree(&mris_in);
    }
#else
    if (COMPARE)
    {
      mris_in = MRISread(inc_fname);
      for (i = 0; i<mris_out->nvertices; i++)
      {
        mris_out->vertices[i].origx += mris_in->vertices[i].x ;
        mris_out->vertices[i].origy += mris_in->vertices[i].y ;
        mris_out->vertices[i].origz += mris_in->vertices[i].z ;
      }
      MRISfree(&mris_in);
    }
#endif

    /* wavelet synthesis of origx, origy, origz */
    wavelet_synthesis_vec(mris_out,order);
#if 1
    {   MRI_SURFACE *mris_high;
      mris_high = ReadIcoByOrder(1, 100);
      for (i=0;i<mris_high->nvertices;i++)
      {
        mris_high->vertices[i].x=mris_out->vertices[i].origx;
        mris_high->vertices[i].y=mris_out->vertices[i].origy;
        mris_high->vertices[i].z=mris_out->vertices[i].origz;
      }
      MRISwrite(mris_high, "/space/birn/42/users/pengyu/simulation/average/surf/lh.wavelet1") ;
      MRISfree(&mris_high);
    }
#endif
    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    fprintf(stdout, "Writing recovered surface to %s\n", argv[4]);
    MRISwrite(mris_out, argv[4]) ;

    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
    /*End of Synthesis*/
  }
  else if ( SYNTHESIS && CURV ) /*Spherical Wavelet Synthesis*/
  {

    mris_out = ReadIcoByOrder(order, 100); //higher order surface
    fprintf(stdout, "Creating a %d order spherical surface\n", order);
    MRISreadCurvatureFile(mris_out, argv[1]) ;
    fprintf(stdout, "Reading wavelet coefficients from %s\n", argv[1]);
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_out, 3) ;


    fprintf(stdout, "Recover the surface using %s order coefficients\n",argv[2]);
    number = IcoNVtxsFromOrder(atoi(argv[2]));
    for (i = number; i<mris_out->nvertices; i++)
    {
      mris_out->vertices[i].curv = 0;
    }

    /* wavelet analysis of curvature */
    wavelet_synthesis_curv(mris_out,order);

    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    fprintf(stdout, "Writing recovered curvature file to %s\n", argv[4]);
    MRISwriteCurvature(mris_out,argv[4]) ;
#if 0
    {  MRI_SURFACE *mris_high;
      mris_high = ReadIcoByOrder(4, 100);
      for (i=0;i<mris_high->nvertices;i++)
      {
        mris_high->vertices[i].x=mris_out->vertices[i].x;
        mris_high->vertices[i].y=mris_out->vertices[i].y;
        mris_high->vertices[i].z=mris_out->vertices[i].z;
      }
      MRISwrite(mris_high, "/space/xrt/1/users/btquinn/buckner_paper/010223_61223/surf/lh.wavelet.recon") ;
      MRISfree(&mris_high);
    }
#endif
    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
    /*End of Synthesis*/
  }
  else if ( SYNTHESIS && RADIUS ) /*Spherical Wavelet Synthesis*/
  {
    char *cp, path[STRLEN], afname[STRLEN];
    float rad;

    mris_out = ReadIcoByOrder(order, 100); //higher order surface
    fprintf(stdout, "Creating a %d order spherical surface\n", order);
    MRISreadCurvatureFile(mris_out, argv[1]) ;
    fprintf(stdout, "Reading wavelet coefficients from %s\n", argv[1]);
    for (i = 0; i<mris_out->nvertices; i++)
      mris_out->vertices[i].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_out, 3) ;

    fprintf(stdout, "Recover the surface using %s order coefficients\n",argv[2]);
    number = IcoNVtxsFromOrder(atoi(argv[2]));
    for (i = number; i<mris_out->nvertices; i++)
    {
      mris_out->vertices[i].curv = 0;
    }

    if (shrink > 0)
      for (i = 0; i<mris_out->nvertices; i++)
        if (fabs(mris_out->vertices[i].curv)<shrink)
          mris_out->vertices[i].curv = 0;

    if (COMPARE)
    {
      mris_in= ReadIcoByOrder(order, 100);
      MRISreadCurvatureFile(mris_in, inc_fname ) ;

      for (i = 0; i<mris_out->nvertices; i++)
        mris_out->vertices[i].curv += mris_in->vertices[i].curv ;
    }

    /* wavelet analysis of curvature */
    wavelet_synthesis_curv(mris_out,order);

    // Read in a template
    FileNamePath(argv[1], path) ;

    cp=strchr(argv[1], '.');
    if (strstr(cp-2,"rh"))
    {
      sprintf(afname, "%s/rh.sphere.reg",path);
      mris_in = MRISread(afname) ;
    }
    else
    {
      sprintf(afname, "%s/lh.sphere.reg",path);
      mris_in = MRISread(afname) ;
    }
    if (strstr(cp-2,"rh.white"))
      sprintf(afname, "%s/rh.white",path);
    else if (strstr(cp-2,"lh.white"))
      sprintf(afname, "%s/lh.white",path);
    else if (strstr(cp-2,"lh.pial"))
      sprintf(afname, "%s/lh.pial",path);
    else if (strstr(cp-2,"rh.pial"))
      sprintf(afname, "%s/rh.pial",path);
    else if (strstr(cp-2,"rh.smoothwm"))
      sprintf(afname, "%s/rh.smoothwm",path);
    else if (strstr(cp-2,"lh.smoothwm"))
      sprintf(afname, "%s/lh.smoothwm",path);
    else if (strstr(cp-2,"rh.inflated"))
      sprintf(afname, "%s/rh.inflated",path);
    else if (strstr(cp-2,"lh.inflated"))
      sprintf(afname, "%s/lh.inflated",path);
    else
      fprintf(stdout, "unrecognized file name %s\n", argv[1]);

    fprintf(stdout, "reading template from %s\n", afname);
    MRISreadOriginalProperties(mris_in, afname) ;
    sample_origposition(mris_in, mris_out) ;

    MRISsaveVertexPositions(mris_out, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris_out, ORIGINAL_VERTICES) ;
    center_brain(mris_out, mris_out);
    volume = brain_volume(mris_out);
    MRISscaleBrain(mris_out, mris_out, cbrt(300000/volume)) ;

    for (i=0;i<mris_out->nvertices;i++)
    {
      rad = (mris_out->vertices[i].x*mris_out->vertices[i].x+mris_out->vertices[i].y*mris_out->vertices[i].y+mris_out->vertices[i].z*mris_out->vertices[i].z);
      rad=sqrt(rad);
      fprintf(stdout, "%f %f\n", rad, mris_out->vertices[i].curv);
      mris_out->vertices[i].x = mris_out->vertices[i].x*mris_out->vertices[i].curv/rad;
      mris_out->vertices[i].y = mris_out->vertices[i].y*mris_out->vertices[i].curv/rad;
      mris_out->vertices[i].z = mris_out->vertices[i].z*mris_out->vertices[i].curv/rad;
    }


    /*write out*/
    volume = brain_volume(mris_out);
    fprintf(stdout, "Writing recovered surface to %s volume %f\n", argv[4],volume);
    MRISwrite(mris_out,argv[4]) ;
    MRISsetNeighborhoodSizeAndDist(mris_out, 2) ;
    MRIScomputeSecondFundamentalForm(mris_out) ;
    MRISuseMeanCurvature(mris_out) ;
    sprintf(afname,"%s.curv", argv[4]);
    MRISwriteCurvature(mris_out, afname) ;

    MRISrestoreVertexPositions(mris_out, TMP_VERTICES) ;
    /*End of Synthesis*/
  }

  MRISfree(&mris_out);
  msec = then.milliseconds() ;
  fprintf(stdout, "spherical wavelet took %2.1f minutes\n", (float)msec/(1000.0f*60.0f));
  exit(0) ;
  return(0) ;
}

static MRI_SURFACE *
wavelet_analysis_curv(MRI_SURFACE *mris_out, int order)
{
  int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0;
  MRIS          *mris_high;
  VERTEX        *vm_out, *vm_high, *v;
  double         s_jkm;

  /* Initialize Ij,k*/
  for (vno = 0 ; vno<mris_out->nvertices; vno++)
  {
    vm_out = &mris_out->vertices[vno];
    vm_out->val = 1;
  }

  /*Iteratively compute Ij,k*/
  for (i=order;i>0;i--)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
    number = IcoNVtxsFromOrder(i-1); //the start of m vertices
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          v->val += 0.5*vm_out->val ;
        }
      for (; nnum<vm_high->v2num; nnum++)
        if ( vm_high->v[nnum]<number ) //B(j,m)
        {
          k = vm_high->v[nnum];
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          v->val += 0.125*vm_out->val ;
        }
      for (; nnum<vm_high->v3num; nnum++)
        if ( vm_high->v[nnum]<number ) //C(j,m)
        {
          k = vm_high->v[nnum];
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            v->val -= 0.0625*vm_out->val ;
          }
        }
    }
  }

  /*Analysis Stage I:*/
  for (i=order;i>0;i--)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;

    number = IcoNVtxsFromOrder(i-1); //the start of m vertices
    /* compute Yj,m for each m vertices */
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
        {
          k = vm_high->v[nnum] ;
          v = &mris_out->vertices[k];
          vm_out->curv -= 0.5*v->curv;
          //if(m==67770) fprintf(stdout, "%f, %d, %f\n", v->curv, k, vm_out->curv);
        }
      for (; nnum<vm_high->v2num; nnum++) //second order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
        {
          k = vm_high->v[nnum] ;
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          vm_out->curv -= 0.125*v->curv;
          //if(m==67770) fprintf(stdout, "%f, %d, %f\n", v->curv, k, vm_out->curv);
        }
      for (; nnum<vm_high->v3num; nnum++)
        if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
        {
          k = vm_high->v[nnum] ;
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            vm_out->curv += 0.0625*v->curv;
            //if(m==67770) fprintf(stdout, "%f, %d, %f\n", v->curv, k, vm_out->curv);
          }
        }
    }


    /*Analysis Stage II: */
    /*Compute Lamda(j,k) using the Yita(j,m)*/
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          s_jkm = vm_out->val/2/v->val;
          //if(k==6642) fprintf(stdout, "%f, %f, %f, %f, %d\n", vm_out->val, v->val, s_jkm, vm_out->curv, m);
          v->curv += s_jkm*vm_out->curv;
        }

    }
  }
  MRISfree(&mris_high) ;
  return(mris_out);
}

static MRI_SURFACE *
wavelet_analysis_vec(MRI_SURFACE *mris_out, int order)
{
  int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0;
  MRIS          *mris_high;
  VERTEX        *vm_out, *vm_high, *v;
  double         s_jkm;

  /* Initialize Ij,k*/
  for (vno = 0 ; vno<mris_out->nvertices; vno++)
  {
    vm_out = &mris_out->vertices[vno];
    vm_out->val = 1;
  }

  /*Iteratively compute Ij,k*/
  for (i=order;i>0;i--)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
    number = IcoNVtxsFromOrder(i-1); //the start of m vertices
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          v->val += 0.5*vm_out->val ;
        }
      for (; nnum<vm_high->v2num; nnum++)
        if ( vm_high->v[nnum]<number ) //B(j,m)
        {
          k = vm_high->v[nnum];
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          v->val += 0.125*vm_out->val ;
        }
      for (; nnum<vm_high->v3num; nnum++)
        if ( vm_high->v[nnum]<number ) //C(j,m)
        {
          k = vm_high->v[nnum];
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            v->val -= 0.0625*vm_out->val ;
          }
        }
    }
  }


  /*Analysis Stage I:*/
  for (i=order;i>0;i--)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;

    number = IcoNVtxsFromOrder(i-1); //the start of m vertices
    /* compute Yj,m for each m vertices */
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
        {
          k = vm_high->v[nnum] ;
          v = &mris_out->vertices[k];
          vm_out->origx -= 0.5*v->origx;
          vm_out->origy -= 0.5*v->origy;
          vm_out->origz -= 0.5*v->origz;
        }
      for (; nnum<vm_high->v2num; nnum++) //second order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
        {
          k = vm_high->v[nnum] ;
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          vm_out->origx -= 0.125*v->origx;
          vm_out->origy -= 0.125*v->origy;
          vm_out->origz -= 0.125*v->origz;
        }
      for (; nnum<vm_high->v3num; nnum++)
        if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
        {
          k = vm_high->v[nnum] ;
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            vm_out->origx += 0.0625*v->origx;
            vm_out->origy += 0.0625*v->origy;
            vm_out->origz += 0.0625*v->origz;
          }
        }
    }


    /*Analysis Stage II: */
    /*Compute Lamda(j,k) using the Yita(j,m)*/
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          s_jkm = vm_out->val/2/v->val;
          v->origx += s_jkm*vm_out->origx;
          v->origy += s_jkm*vm_out->origy;
          v->origz += s_jkm*vm_out->origz;
        }

    }

  }
  MRISfree(&mris_high) ;
  return(mris_out);
}

static MRI_SURFACE *
wavelet_synthesis_curv(MRI_SURFACE *mris_out, int order)
{
  int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0;
  MRIS          *mris_high;
  VERTEX        *vm_out, *vm_high, *v;
  double         s_jkm;

  /*Initialize Ij,k*/
  for (vno = 0; vno<mris_out->nvertices; vno++)
  {
    vm_out = &mris_out->vertices[vno];
    vm_out->val = 1;
  }

  /*Iteratively compute Ij,k*/
  for (i=order;i>0;i--)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
    number = IcoNVtxsFromOrder(i-1); //the start of m vertices
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          v->val += 0.5*vm_out->val ;
        }
      for (; nnum<vm_high->v2num; nnum++)
        if ( vm_high->v[nnum]<number ) //B(j,m)
        {
          k = vm_high->v[nnum];
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          v->val += 0.125*vm_out->val ;
        }
      for (; nnum<vm_high->v3num; nnum++)
        if ( vm_high->v[nnum]<number ) //C(j,m)
        {
          k = vm_high->v[nnum];
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            v->val -= 0.0625*vm_out->val ;
          }
        }
    }
  }


  for (i=1;i<=order;i++)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
    number = IcoNVtxsFromOrder(i-1); //the start of m vertices

    /* Synthesis Stage I */
    /* Compute Lamda(j+1,k) using the Yita(j,m) */
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          s_jkm = vm_out->val/2/v->val;
          v->curv -= s_jkm*vm_out->curv;
        }
    }

    /* compute Lamda(j+1,m) for each m vertices */
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
        {
          k = vm_high->v[nnum] ;
          v = &mris_out->vertices[k];
          vm_out->curv += 0.5*v->curv;
        }
      for (; nnum<vm_high->v2num; nnum++) //second order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
        {
          k = vm_high->v[nnum] ;
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          vm_out->curv += 0.125*v->curv;
        }
      for (; nnum<vm_high->v3num; nnum++) //third order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
        {
          k = vm_high->v[nnum] ;
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            vm_out->curv -= 0.0625*v->curv;
          }
        }
    }
  }
  MRISfree(&mris_high) ;
  return(mris_out);
}

static MRI_SURFACE *
wavelet_synthesis_vec(MRI_SURFACE *mris_out, int order)
{
  int           i, number, vno, nnum, m, k, b1=0, b2=0, cno, flag=0;
  MRIS          *mris_high;
  VERTEX        *vm_out, *vm_high, *v;
  double         s_jkm;


  /*Initialize Ij,k*/
  for (vno = 0; vno<mris_out->nvertices; vno++)
  {
    vm_out = &mris_out->vertices[vno];
    vm_out->val = 1;
  }

  /*Iteratively compute Ij,k*/
  for (i=order;i>0;i--)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
    number = IcoNVtxsFromOrder(i-1); //the start of m vertices
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          v->val += 0.5*vm_out->val ;
        }
      for (; nnum<vm_high->v2num; nnum++)
        if ( vm_high->v[nnum]<number ) //B(j,m)
        {
          k = vm_high->v[nnum];
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          v->val += 0.125*vm_out->val ;
        }
      for (; nnum<vm_high->v3num; nnum++)
        if ( vm_high->v[nnum]<number ) //C(j,m)
        {
          k = vm_high->v[nnum];
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            v->val -= 0.0625*vm_out->val ;
          }
        }
    }
  }


  for (i=1;i<=order;i++)
  {
    mris_high = ReadIcoByOrder(i, 100); //higher order surface
    for (m = 0; m<mris_high->nvertices; m++)
      mris_high->vertices[m].nsize=1;
    MRISsetNeighborhoodSizeAndDist(mris_high, 3) ;
    number = IcoNVtxsFromOrder(i-1); //the start of m vertices

    /* Synthesis Stage I */
    /* Compute Lamda(j+1,k) using the Yita(j,m) */
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      for (nnum=0; nnum<vm_high->vnum; nnum++)
        if ( vm_high->v[nnum]<number ) //A(j,m)
        {
          k = vm_high->v[nnum];
          v = &mris_out->vertices[k];
          s_jkm = vm_out->val/2/v->val;
          v->origx -= s_jkm*vm_out->origx;
          v->origy -= s_jkm*vm_out->origy;
          v->origz -= s_jkm*vm_out->origz;
        }
    }

    /* compute Lamda(j+1,m) for each m vertices */
    for (m = number; m<mris_high->nvertices; m++)
    {
      vm_out = &mris_out->vertices[m];
      vm_high = &mris_high->vertices[m];
      flag=0;
      for (nnum=0; nnum<vm_high->vnum; nnum++)  //first order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor A(j,m)
        {
          k = vm_high->v[nnum] ;
          v = &mris_out->vertices[k];
          vm_out->origx += 0.5*v->origx;
          vm_out->origy += 0.5*v->origy;
          vm_out->origz += 0.5*v->origz;
        }
      for (; nnum<vm_high->v2num; nnum++) //second order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor B(j,m)
        {
          k = vm_high->v[nnum] ;
          if (flag==0) b1=k;
          else b2=k;
          flag++;
          v = &mris_out->vertices[k];
          vm_out->origx += 0.125*v->origx;
          vm_out->origy += 0.125*v->origy;
          vm_out->origz += 0.125*v->origz;
        }
      for (; nnum<vm_high->v3num; nnum++) //third order neighborhood
        if ( vm_high->v[nnum]<number ) //neighbor C(j,m)
        {
          k = vm_high->v[nnum] ;
          flag=0; //C has to be a second-order neighbor of B
          for (cno=mris_high->vertices[b1].vnum; cno<mris_high->vertices[b1].v2num;cno++)
            if (mris_high->vertices[b1].v[cno]==k) flag=1;
          for (cno=mris_high->vertices[b2].vnum; cno<mris_high->vertices[b2].v2num;cno++)
            if (mris_high->vertices[b2].v[cno]==k) flag=1;
          if (flag)
          {
            v = &mris_out->vertices[k];
            vm_out->origx -= 0.0625*v->origx;
            vm_out->origy -= 0.0625*v->origy;
            vm_out->origz -= 0.0625*v->origz;
          }
        }
    }
  }
  MRISfree(&mris_high) ;
  return(mris_out);
}

double
brain_volume(MRI_SURFACE *mris)
{
  int fno;
  FACE *face;
  double total_volume, face_area;
  VECTOR *v_a, *v_b, *v_n, *v_cen;
  VERTEX  *v0, *v1, *v2;

  v_a = VectorAlloc(3, MATRIX_REAL) ;
  v_b = VectorAlloc(3, MATRIX_REAL) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;       /* normal vector */
  v_cen = VectorAlloc(3, MATRIX_REAL) ;     /* centroid vector */

  total_volume = 0;
  for (fno = 0 ; fno < mris->nfaces ; fno++)
  {
    face = &mris->faces[fno] ;
    if (face->ripflag)
      continue ;

    v0 = &mris->vertices[face->v[0]] ;
    v1 = &mris->vertices[face->v[1]] ;
    v2 = &mris->vertices[face->v[2]] ;

    VERTEX_EDGE(v_a, v0, v1) ;
    VERTEX_EDGE(v_b, v0, v2) ;

    /* face normal vector */
    V3_CROSS_PRODUCT(v_a, v_b, v_n) ;
    face_area = V3_LEN(v_n) * 0.5f ;

    V3_NORMALIZE(v_n, v_n) ;             /* make it a unit vector */

    /* compute face centroid */
    V3_X(v_cen) = (v0->x + v1->x + v2->x)/3.0;
    V3_Y(v_cen) = (v0->y + v1->y + v2->y)/3.0;
    V3_Z(v_cen) = (v0->z + v1->z + v2->z)/3.0;

    total_volume += V3_DOT(v_cen, v_n)*face_area;
  }

  total_volume /= 3.0;

  return(total_volume);
}


static MRI_SURFACE *
center_brain(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int         fno, vno ;
  FACE        *face ;
  VERTEX      *vdst;
  double       x, y, z, x0, y0, z0 ;

  if (!mris_dst)
    mris_dst = MRISclone(mris_src) ;

  x0 = y0 = z0 = 0 ;   /* silly compiler warning */
  //MRISupdateSurface(mris_dst);
#if 1 //area weighting
  for (fno = 0 ; fno < mris_src->nfaces ; fno++)
  {
    face = &mris_dst->faces[fno] ;
    if (face->ripflag)
      continue ;
    x = mris_dst->vertices[face->v[0]].x;
    y = mris_dst->vertices[face->v[0]].y;
    z = mris_dst->vertices[face->v[0]].z;
    x += mris_dst->vertices[face->v[1]].x;
    y += mris_dst->vertices[face->v[1]].y;
    z += mris_dst->vertices[face->v[1]].z;
    x += mris_dst->vertices[face->v[2]].x;
    y += mris_dst->vertices[face->v[2]].y;
    z += mris_dst->vertices[face->v[2]].z;
    x = x*face->area/3;
    y = y*face->area/3;
    z = z*face->area/3;
    x0 += x/mris_dst->total_area;
    y0 += y/mris_dst->total_area;
    z0 += z/mris_dst->total_area;
    if (face->area == 0) fprintf(stdout, "%d %f\n", fno, face->area);
  }
#else
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno];
    x = vdst->x;
    y = vdst->y;
    z = vdst->z;
    x0 += x/mris_dst->nvertices;
    y0 += y/mris_dst->nvertices;
    z0 += z/mris_dst->nvertices;
  }
#endif
  fprintf(stdout, "brain center is found at %f %f %f\n", x0, y0, z0);
  for (vno = 0 ; vno < mris_src->nvertices ; vno++)
  {
    vdst = &mris_dst->vertices[vno] ;
    if (vdst->ripflag)
      continue ;
    vdst->x -= x0 ;
    vdst->y -= y0 ;
    vdst->z -= z0 ;
  }

  mris_dst->xctr = mris_dst->yctr = mris_dst->zctr = 0 ;
  return(mris_dst) ;
}


static MRI_SURFACE *
sample_origposition(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int index, fno, fnum=0, i;
  VERTEX *vertex;
  double nearest, dist, r, s, t;
  double a, b, c, p;
  ANNpointArray pa = annAllocPts(mris_src->nvertices, 3);

  for (index = 0; index < mris_src->nvertices; index++)
  {
    pa[index][0] = mris_src->vertices[index].x;
    pa[index][1] = mris_src->vertices[index].y;
    pa[index][2] = mris_src->vertices[index].z;
  }

  ANNkd_tree *annkdTree = new ANNkd_tree(pa, mris_src->nvertices, 3);
  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];
  //  ANNpoint query_pt = annAllocPt(3);
  ANNpointArray QueryPt;

  //if(mris_dst == NULL)
  // mris_dst = MRIScopy(mris_src, NULL);

  QueryPt = annAllocPts(1,3);

  for (index = 0; index < mris_dst->nvertices; index++)
  {
    if (mris_dst->vertices[index].border == 1) continue;

    QueryPt[0][0] = mris_dst->vertices[index].x;
    QueryPt[0][1] = mris_dst->vertices[index].y;
    QueryPt[0][2] = mris_dst->vertices[index].z;

    annkdTree->annkSearch( // search
      QueryPt[0],       // query point
      1,   // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,  // distance (returned)
      0);   // error bound

#if 1
    vertex = &mris_src->vertices[annIndex[0]];
    nearest = 100000;
    for (i=0; i<vertex->num; i++)
    {
      fno = vertex->f[i];
      dist = v_to_f_distance(&mris_dst->vertices[index], mris_src, fno, 0);
      if (dist<nearest)
      {
        nearest = dist;
        fnum = fno;
      }
    }

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    r = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[0]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[0]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[0]].z),2));
    p = (a+b+c)/2;
    s = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    t = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));
    p= (r+s+t);
    r = r/p;
    s = s/p;
    t = t/p;

    mris_dst->vertices[index].origx = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origx + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origx + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origx;
    mris_dst->vertices[index].origy = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origy + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origy + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origy;
    mris_dst->vertices[index].origz = r*mris_src->vertices[mris_src->faces[fnum].v[0]].origz + s*mris_src->vertices[mris_src->faces[fnum].v[1]].origz + t*mris_src->vertices[mris_src->faces[fnum].v[2]].origz;

#else
    mris_dst->vertices[index].origx = mris_src->vertices[annIndex[0]].origx;
    mris_dst->vertices[index].origy = mris_src->vertices[annIndex[0]].origy;
    mris_dst->vertices[index].origz = mris_src->vertices[annIndex[0]].origz;
#endif
    if (annIndex[0] == 81426)
      printf("src index %d dst index %d face %d r %f s %f t %f\n", index, annIndex[0], fnum, r, s, t);
  }
  return(mris_dst) ;
}

static MRI_SURFACE *
sample_origcurvature(MRI_SURFACE *mris_src, MRI_SURFACE *mris_dst)
{
  int index, fno, fnum=0, i;
  VERTEX *vertex;
  double nearest, dist, r, s, t;
  double a, b, c, p;
  ANNpointArray pa = annAllocPts(mris_src->nvertices, 3);

  for (index = 0; index < mris_src->nvertices; index++)
  {
    pa[index][0] = mris_src->vertices[index].x;
    pa[index][1] = mris_src->vertices[index].y;
    pa[index][2] = mris_src->vertices[index].z;
  }

  ANNkd_tree *annkdTree = new ANNkd_tree(pa, mris_src->nvertices, 3);
  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];
  //  ANNpoint query_pt = annAllocPt(3);
  ANNpointArray QueryPt;

  //if(mris_dst == NULL)
  // mris_dst = MRIScopy(mris_src, NULL);

  QueryPt = annAllocPts(1,3);

  for (index = 0; index < mris_dst->nvertices; index++)
  {
    if (mris_dst->vertices[index].border == 1) continue;

    QueryPt[0][0] = mris_dst->vertices[index].x;
    QueryPt[0][1] = mris_dst->vertices[index].y;
    QueryPt[0][2] = mris_dst->vertices[index].z;

    annkdTree->annkSearch( // search
      QueryPt[0],       // query point
      1,   // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,  // distance (returned)
      0);   // error bound

#if 1
    vertex = &mris_src->vertices[annIndex[0]];
    nearest = 100000;
    for (i=0; i<vertex->num; i++)
    {
      fno = vertex->f[i];
      dist = v_to_f_distance(&mris_dst->vertices[index], mris_src, fno, 0);
      if (dist<nearest)
      {
        nearest = dist;
        fnum = fno;
      }
    }

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    r = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[2]].x - mris_src->vertices[mris_src->faces[fnum].v[0]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].y - mris_src->vertices[mris_src->faces[fnum].v[0]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[2]].z - mris_src->vertices[mris_src->faces[fnum].v[0]].z),2));
    p = (a+b+c)/2;
    s = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    a = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[1]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[1]].z - mris_dst->vertices[index].z),2));
    b = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_dst->vertices[index].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_dst->vertices[index].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_dst->vertices[index].z),2));
    c = sqrt( pow((mris_src->vertices[mris_src->faces[fnum].v[0]].x - mris_src->vertices[mris_src->faces[fnum].v[1]].x),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].y - mris_src->vertices[mris_src->faces[fnum].v[1]].y),2)+
              pow((mris_src->vertices[mris_src->faces[fnum].v[0]].z - mris_src->vertices[mris_src->faces[fnum].v[1]].z),2));
    p = (a+b+c)/2;
    t = sqrt(fabs(p*(p-a)*(p-b)*(p-c)));

    p= (r+s+t);
    r = r/p;
    s = s/p;
    t = t/p;

    mris_dst->vertices[index].curv = r*mris_src->vertices[mris_src->faces[fnum].v[0]].curv + s*mris_src->vertices[mris_src->faces[fnum].v[1]].curv + t*mris_src->vertices[mris_src->faces[fnum].v[2]].curv;
#else

    mris_dst->vertices[index].curv = mris_src->vertices[annIndex[0]].curv;
    if (index == 67770)
      printf("src index %d dst index %d face %d r %f s %f t %f\n", index, annIndex[0], fnum, r, s, t);
#endif

  }
  return(mris_dst) ;
}
#if 0
#if 0
static void
register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2)
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

static void
FindClosest(MRI_SURFACE *TrueMesh, ANNkd_tree *annkdTree, MRI_SURFACE *EstMesh, double3d *closest)
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

static double
transformS(double3d *V1a, double3d *V2a, int N, double TR[3][3],double shift[3])
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
static void
jacobi(float **a, int n,float *d, float **v, int *nrot)
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

#else
static void
register2to1(MRI_SURFACE *Surf1, MRI_SURFACE *Surf2, double T[5][5])
{
  //transform Surf2 to match Surf1
  //This version uses (x,y,z) to build correspondence, and (origx, origy, origz)to compute coordinate transformation
  double3d *mesh2avtx, *closevtx2a;
  int index, k;

  /* Build the ANN tree repeatedly is time consuming, so
  *  move it outside of the function
  */
  ANNpointArray pa = annAllocPts(Surf1->nvertices, 3);

  //note that using the sphere.reg coordinates to build correspondence
  // but use original coordinates to compute coordiante mapping
  for (index = 0; index < Surf1->nvertices; index++)
  {
    pa[index][0] = Surf1->vertices[index].x;
    pa[index][1] = Surf1->vertices[index].y;
    pa[index][2] = Surf1->vertices[index].z;
  }

  ANNkd_tree *annkdTree = new ANNkd_tree(pa, Surf1->nvertices, 3);

  mesh2avtx = (double3d*) malloc(Surf2->nvertices*sizeof(double3d));
  closevtx2a = (double3d*) malloc(Surf2->nvertices*sizeof(double3d));

  /* Initialize */
  for (index = 0; index < Surf2->nvertices; index++)
  {
    mesh2avtx[index].x = Surf2->vertices[index].origx;
    mesh2avtx[index].y = Surf2->vertices[index].origy;
    mesh2avtx[index].z = Surf2->vertices[index].origz;
  }

  {
    /* For each vertex in Surf2, find its closest vertex in Surf1
    * (based on their spherical cooridnates),
    * and record the original coordinates in closevtx2a
    */
    FindClosest(Surf1, annkdTree, Surf2, closevtx2a);


    // Find the affine transformation
    // transformS(mesh2avtx, closevtx2a, Surf2->nvertices, TR, shift);
    affineS(mesh2avtx, closevtx2a,Surf2->nvertices, T);

    //Surf2 transformed
    for (k = 0; k < Surf2->nvertices; k++)
    {
      Surf2->vertices[k].origx = mesh2avtx[k].x;
      Surf2->vertices[k].origy  = mesh2avtx[k].y;
      Surf2->vertices[k].origz  = mesh2avtx[k].z;
    }
  }

  // free memory
  delete pa;
  if (annkdTree) delete annkdTree;

  free(mesh2avtx);
  free(closevtx2a);

  return;
}

static void
FindClosest(MRI_SURFACE *TrueMesh, ANNkd_tree *annkdTree, MRI_SURFACE *EstMesh, double3d *closest)
{
  int index;

  ANNidxArray annIndex = new ANNidx[1];
  ANNdistArray annDist = new ANNdist[1];

  ANNpoint Qpt;

  Qpt = (ANNpoint)malloc(3*sizeof(ANNcoord));

  for (index =0; index < EstMesh->nvertices; index++)
  {
    Qpt[0] = EstMesh->vertices[index].x;
    Qpt[1] = EstMesh->vertices[index].y;
    Qpt[2] = EstMesh->vertices[index].z;

    annkdTree->annkSearch( // search
      Qpt,  // query point
      1,    // number of near neighbors
      annIndex,  // nearest neighbors (returned)
      annDist,   // distance (returned)
      0);      // error bound

    //NOTE: ORIG COORD are used to compute transformation
    closest[index].x = TrueMesh->vertices[annIndex[0]].origx;
    closest[index].y = TrueMesh->vertices[annIndex[0]].origy;
    closest[index].z = TrueMesh->vertices[annIndex[0]].origz;
  }

  if (annIndex) delete annIndex;
  if (annDist) delete annDist;
  free(Qpt);

  return;
}

static double
affineS(double3d *V1, double3d *V2, int N, double T[5][5])
//transform V1 to fit V2. Here V1 and V2 are 4XN matrix, the 4th row of them are all 1's.
// T = V2*Trans(V1)*Inv(V1*Trans(V1)), B = V2*Trans(V1), A = V1*Trans(V1)
// T is a 4X4 matrix, which takes care of everything: translation, rotation, and scaling
{
  double *A[5], *B[5], C[5][5]; // C = inv(A); A is symmetric.
  int k,l,m,n;
  double col[5], d;
  int *indx;
  double3d v;
  //  FILE *fp;
  double error;
  // int count;

  indx  = (int*)malloc(10*sizeof(int));
  for (k = 1; k <= 4; k++)
  {
    A[k] = (double*)malloc(5*sizeof(double));
    B[k] = (double*)malloc(5*sizeof(double));
    for (l = 1; l <= 4; l ++)
    {
      A[k][l] = 0;
      B[k][l] = 0;
    }
  }

  for (n = 0; n < N; n++)
    if (1)
    {
      A[1][1] += V1[n].x * V1[n].x;
      A[1][2] += V1[n].x * V1[n].y;
      A[1][3] += V1[n].x * V1[n].z;
      A[1][4] += V1[n].x;

      A[2][1] += V1[n].y * V1[n].x;
      A[2][2] += V1[n].y * V1[n].y;
      A[2][3] += V1[n].y * V1[n].z;
      A[2][4] += V1[n].y;

      A[3][1] += V1[n].z * V1[n].x;
      A[3][2] += V1[n].z * V1[n].y;
      A[3][3] += V1[n].z * V1[n].z;
      A[3][4] += V1[n].z;

      A[4][1] += V1[n].x;
      A[4][2] += V1[n].y;
      A[4][3] += V1[n].z;
      A[4][4] += 1;

      B[1][1] += V2[n].x * V1[n].x;
      B[1][2] += V2[n].x * V1[n].y;
      B[1][3] += V2[n].x * V1[n].z;
      B[1][4] += V2[n].x;

      B[2][1] += V2[n].y * V1[n].x;
      B[2][2] += V2[n].y * V1[n].y;
      B[2][3] += V2[n].y * V1[n].z;
      B[2][4] += V2[n].y;

      B[3][1] += V2[n].z * V1[n].x;
      B[3][2] += V2[n].z * V1[n].y;
      B[3][3] += V2[n].z * V1[n].z;
      B[3][4] += V2[n].z;

      B[4][1] += V1[n].x;
      B[4][2] += V1[n].y;
      B[4][3] += V1[n].z;
      B[4][4] += 1;
    }

#if 0
  fp = fopen("m.dat","a");
  for (k = 1; k <= 4; k++)
  {
    for (l = 1; l <= 4; l ++)
      fprintf(fp, "%g ", A[k][l]);
    fprintf(fp,"\n");
  }
  for (k = 1; k <= 4; k++)
  {
    for (l = 1; l <= 4; l ++)
      fprintf(fp, "%g ", B[k][l]);
    fprintf(fp,"\n");
  }
  fclose(fp);
#endif

  ludcmp(A,4,indx,&d);
  for (k = 1; k <= 4; k ++)
  {
    for (l = 1; l <= 4; l++) col[l] = 0;
    col[k] = 1.0;
    lubksb(A,4,indx,col);
    for (l = 1; l <= 4; l++)
      C[l][k] = col[l];
  }

#if 0
  fp = fopen("m.dat","a");
  for (k = 1; k <= 4; k++)
  {
    for (l = 1; l <= 4; l ++)
      fprintf(fp, "%g ", C[k][l]);
    fprintf(fp,"\n");
  }
  fclose(fp);

#endif

  for (k = 1; k <= 4; k++)
    for (l = 1; l <= 4; l ++)
    {
      A[k][l] = 0;
      for (m = 1; m <= 4; m++)
        A[k][l] += B[k][m]*C[m][l];
    }

#if 0
  fp = fopen("m.dat","a");
  for (k = 1; k <= 4; k++)
  {
    for (l = 1; l <= 4; l ++)
      fprintf(fp, "%g ", A[k][l]);
    fprintf(fp,"\n");
  }
  fclose(fp);
#endif

  printf("The transformation matrix is:\n");
  for (k = 1; k <= 4; k++)
  {
    printf("\t");
    for (l = 1; l <= 4; l ++)
    {
      printf("%13.9f ",A[k][l]);
      T[k][l] = A[k][l];
    }
    printf("\n");
  }

  error = 0;
  for (k = 0; k < N; k++)
  {
    v.x = A[1][1] * V1[k].x + A[1][2] * V1[k].y + A[1][3] * V1[k].z + A[1][4];
    v.y = A[2][1] * V1[k].x + A[2][2] * V1[k].y + A[2][3] * V1[k].z + A[2][4];
    v.z = A[3][1] * V1[k].x + A[3][2] * V1[k].y + A[3][3] * V1[k].z + A[3][4];

    error += (v.x-V2[k].x)*(v.x-V2[k].x) + (v.y-V2[k].y)*(v.y-V2[k].y) + (v.z-V2[k].z)*(v.z-V2[k].z);
    V1[k].x = v.x;
    V1[k].y = v.y;
    V1[k].z = v.z;
  }

  for (k = 1; k <= 4; k++)
  {
    free(A[k]);
    free(B[k]);
  }

  return(sqrt(error));
}
static void
ludcmp(double** a,int n,int *indx,double* d)
{
  int i,imax,j,k;
  double big,dum,sum,temp;
  double *vv;

  imax = 0;

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

static void
lubksb(double** a,int n,int* indx,double* b)
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

#endif
#endif
static double
v_to_f_distance(VERTEX *P0, MRI_SURFACE *mri_surf, int face_number, int debug)
{
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
        if (debug) printf("region 4, s =%g, t =%g\n", s, t);
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


#if 0
static int
sort(double **array, int order, int number, int total_number)
{
  int i, j;
  double temp, index;

  if (order)
  {
    for (j=0; j<number; j++)
      for (i=total_number-1; i>j ; i--)
        if (array[0][i]<array[0][i-1])
        {
          temp=array[0][i];
          array[0][i]=array[0][i-1];
          array[0][i-1]=temp;
          index=array[1][i];
          array[1][i]=array[1][i-1];
          array[1][i-1]=index;
        }
  }
  else
  {
    for (j=0; j<number; j++)
      for (i=total_number-1; i>j ; i--)
        if (array[0][i]>array[0][i-1])
        {
          temp=array[0][i];
          array[0][i]=array[0][i-1];
          array[0][i-1]=temp;
          index=array[1][i];
          array[1][i]=array[1][i-1];
          array[1][i-1]=index;
        }
  }
  return(NO_ERROR);
}
#endif
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

  if (!stricmp(option, "A"))
  {
    ANALYSIS = 1 ;
    fprintf(stdout,"spherical wavelet analysis\n");
  }
  else if (!stricmp(option, "S"))
  {
    SYNTHESIS = 1 ;
    fprintf(stdout,"Spherical wavelet synthesis\n");
  }
  else if (!stricmp(option, "C"))
  {
    COMPARE = 1;
    inc_fname = argv[2] ;
    threshold = atof(argv[3]);
    fprintf(stdout,"Reconstruct the surface using only coefs changes %f\n", threshold);
    nargs = 2;
  }
  else if (!stricmp(option, "sampled"))
  {
    SAMPLE_OUT = 1;
    outs_fname = argv[2] ;
    fprintf(stdout,"Write sampled surface to %s\n",outs_fname);
    nargs = 1;
  }
  else if (!stricmp(option, "CURV"))
  {
    CURV = 1;
    fprintf(stdout,"Read in the curvature file and decompose it\n");
  }
  else if (!stricmp(option, "RADIUS"))
  {
    RADIUS = 1;
    fprintf(stdout,"Decompose the radius\n");
  }
  else if (!stricmp(option, "WS"))
  {
    shrink = atof(argv[2]);
    fprintf(stdout, "Wavelet shrinkage using threshold %f\n", shrink);
    nargs = 1;
  }
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      fprintf(stdout,
              "usage: %s <input surface> <orig surface> <finest order> <output surface> <output volume for glm>\n",
              Progname) ;
      fprintf(stdout, "-Option you must have: \n");
      fprintf(stdout, "        -S (synthesis) or -A (Analysis) \n");
      fprintf(stdout, "-Other options: \n");
      fprintf(stdout, "        -CURV   Handle the orig surface file as a curv file) \n");
      fprintf(stdout, "        -C surface_name threshold  Compare orig surface with another surface \
              and synthesize a new surface using part of the new coefficients\n");
      fprintf(stdout, "        -sampled   Write out a sampled version of orig surface \n");
      exit(1) ;
      break ;
    default:
      fprintf(stdout, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}






