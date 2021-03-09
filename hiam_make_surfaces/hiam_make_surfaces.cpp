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



///////////////////////////////////////////
// hiam_make_surfaces.c
//
// written by Peng Yu
// date: 01/27/04
//
// Warning: Do not edit the following four lines.  CVS maintains them.
////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"

#include "mri.h"
#include "mrimorph.h"
#include "mrinorm.h"

#include "mrisurf.h"
#include "mrisurf_project.h"
#include "mrishash_internals.h"

#include "error.h"
#include "diag.h"
#include "proto.h"
#include "timer.h"
#include "matrix.h"

#define MAX_4_NEIGHBORS     100
#define MAX_3_NEIGHBORS     70
#define MAX_2_NEIGHBORS     20
#define MAX_1_NEIGHBORS     8
#define REPULSE_K           1.0

#undef  REPULSE_E           // differs from mrisurf_base.h 0.25
#define REPULSE_E           0.5

#define MAX_MOMENTUM_MM     1

static int    get_option(int argc, char *argv[]) ;
static void   usage_exit(void) ;
static void   print_usage(void) ;
static int    extractlabelvolume(MRI *mri_label[5], MRI *mri_orig);
static int    mrisFindneighborlabel(MRI_SURFACE *mris, char surftype[10], MRI *mri_label[5], MRI *mri_orig);
static int    mrisExaminemovelength(MRI_SURFACE *mris);
static int    mrisClearGradient(MRI_SURFACE *mris);
static int    mrisClearMomentum(MRI_SURFACE *mris);

static int    mrisComputeLabelTerm1(MRI_SURFACE *mris, double weight_label,  MRI *mri_smooth[5],
                                    MRI *mri_label[5], MRI *mri_orig);
static double mrisComputeLabelEnergy(MRI_SURFACE *mris, MRI *mri_smooth[5],
                                     MRI *mri_label[5], MRI *mri_orig);
static int    mrisComputeRepulsiveTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht);
static double mrisComputeRepulsiveEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht);
static int    mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring);
static double mrisComputeNormalSpringEnergy(MRI_SURFACE *mris);
static int    mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring);
static double mrisComputeTangentialSpringEnergy(MRI_SURFACE *mris);
static int    mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv);
static double mrisComputeQuadraticCurvatureEnergy(MRI_SURFACE *mris);
static int    mriSspringTermWithGaussianCurvature(MRI_SURFACE *mris, double gaussian_norm, double l_spring);
static double mrisComputeGaussianCurvatureSpringEnergy(MRI_SURFACE *mris, double gaussian_norm);

static double mrismomentumTimeStep(MRI_SURFACE *mris, float momentum, float dt, float tol,
                                   float n_averages);
static int    my_mrisProjectSurface(MRI_SURFACE *mris);
static int    my_mrisComputeTangentPlanes(MRI_SURFACE *mris);
static int    FindSpikes(MRI_SURFACE *mris, int iter);
static int    SmoothSpikes(MRI_SURFACE *mris, int niter);

static int      all_flag = 0 ;
static const char     *suffix = "hippocampus" ;
static const char     *labelvolume = "mri/aseg" ;
//static char     *labelvolume = "mri/aseg_new.mgh" ;
static const char     *orig_name =  "hippocampus.orig" ;
static float    weight_Gspring = 0.0 ;
static double   gaussian_norm = 2.0;
static int      write_iterations = 0 ;
static int      niteration = 0 ;
static int      smooth_spikes = 500;
static int      nbrs = 3 ;
//static float    weight_quadcur = 0.0, weight_label = 1.0, weight_repulse = 0.0, weight_Nspring = 0.1, weight_Tspring = 0.1;
static float    weight_quadcur = 1.2, weight_label = 1.2, weight_repulse = 3.0, weight_Nspring = 0.5, weight_Tspring = 0.5;
const char            *Progname ;
int             t=0;
int             table[2][20000];
char            surf[5][10]= {"lh","rh","LA","RA","OTHER"
                             };
float           sigma = 2.0f ;
float           MAX_mag = 2.0, threshold = 0.5 ;
float           rmax = 10, rmin = 3.3 ;
float           step_size = 1;


int main(int argc, char *argv[]) ;

int
main(int argc, char *argv[]) {
  char          data_dir[400], *cp, ifname[200], ofname[200], labelfilename[200], surftype[10];
  int           nargs, s, counter=0, i, spikes=1;
  float         ratio = 1, energy_new = 0, energy_old = 0 ;
  // float         weight_quadcur = 0.2, weight_label = 0.5, weight_repulse = 0.0,
  //            weight_Nspring = 0.3, weight_Tspring = 0.3;
  MRI_SURFACE   *mris;
  MRI           *mri_orig, *mri_label[5];
  MHT           *mht_v_current = NULL ;
  double        energy_quadcur=0, energy_label=0, energy_repulse = 0,
                               energy_Nspring = 0 , energy_Tspring = 0, energy_Gspring = 0;
  MRI *mri_e;
  int xi, yi, zi;

  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;
  Progname = argv[0] ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    usage_exit() ;

  cp = getenv("SUBJECTS_DIR");
  if (cp==NULL) {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(1);
  }

  strcpy(data_dir, cp) ;
  strcpy(surftype, argv[2]);

  //***** Extract volume for each label  *****//
  sprintf(labelfilename,"%s/%s/%s",data_dir,argv[1],labelvolume);
  fprintf(stderr, "reading segmentation volume from %s...\n", labelfilename) ;
  mri_orig = MRIread(labelfilename) ;
  extractlabelvolume(mri_label, mri_orig);
  /* added temporarily */
  mri_e = MRIcopy(mri_orig, NULL);
  for (xi=0; xi<mri_e->depth; xi++)
    for (yi=0; yi<mri_e->height; yi++)
      for (zi=0; zi<mri_e->width; zi++) {
        if ( ((xi-127)*(xi-127)/10/10)+((yi-127)*(yi-127)/3/3)+((zi-127)*(zi-127)/3/3) <= 1 )
          mri_e->slices[xi][yi][zi] = 17;
        else
          mri_e->slices[xi][yi][zi] = 0;
      }
  MRIwrite(mri_e, "/autofs/space/dijon_004/ksong/BIRN_processed_data/prospective/buckner/CORTICAL_THINNING/RECONS/001009_vc5398/surf/ellipsoid.mgh");
  MRIfree(&mri_e);
  /******************************/
  //******** read original tesselation ********//
  sprintf(ifname,"%s/%s/surf/%s.%s",data_dir,argv[1],argv[2],orig_name);
  fprintf(stderr, "reading original surface position from %s...\n", ifname) ;
  mris = MRISread(ifname) ;

  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, ofname) ;

  //*** Set the neighborhoodsize to be 3 before calculate curvature term the first time ****//
  //  MRISaverageVertexPositions(mris, 15);
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
#if 0
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  sprintf(ifname,"%s/%s/surf/%s.hippocampus.orig",data_dir,argv[1],argv[2]);
  if (MRISreadVertexPositions(mris, ifname) != NO_ERROR)
    ErrorReturn(ERROR_BADFILE,
                (ERROR_BADFILE,
                 "MRISreadOriginalProperties: could not read surface file %s",
                 ifname)) ;
  MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
#endif
  MRISsetNeighborhoodSizeAndDist(mris, nbrs) ;
  MRISupdateSurface(mris);
  
  MHTfree(&mht_v_current);
  mht_v_current = MHTcreateVertexTable_Resolution(mris, CURRENT_VERTICES, 1.0f);

  //*** Find outside label and inside label for each surface vertex  ****//
  mrisFindneighborlabel(mris,surftype, mri_label, mri_orig);

  MRISuseMeanCurvature(mris);

  //*** Iterate four times using sigma 2, 1.5, 1.0, 0.5 ************//
  for ( s=0; s<4; s++) {
    MRI *mri_kernel;
    MRI *mri_smooth[5];

    if (s==0) weight_repulse = 0.0;

    mri_smooth[0] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
    mri_smooth[1] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
    mri_smooth[2] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
    mri_smooth[3] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
    mri_smooth[4] = MRIalloc(256, 256, 256, MRI_FLOAT) ;

    t = 0;
    ratio = 1;
    energy_old = energy_new = 0;
    mri_kernel = MRIgaussian1d(sigma, 100) ;
    fprintf(stderr, "smoothing label volume with sigma = %2.3f\n", sigma) ;

    mri_smooth[0] = MRIconvolveGaussian(mri_label[0], NULL, mri_kernel) ;
    mri_smooth[1] = MRIconvolveGaussian(mri_label[1], NULL, mri_kernel) ;
    mri_smooth[2] = MRIconvolveGaussian(mri_label[2], NULL, mri_kernel) ;
    mri_smooth[3] = MRIconvolveGaussian(mri_label[3], NULL, mri_kernel) ;
    mri_smooth[4] = MRIconvolveGaussian(mri_label[4], NULL, mri_kernel) ;
    MRIfree(&mri_kernel) ;


    energy_label = mrisComputeLabelEnergy(mris, mri_smooth, mri_label, mri_orig);
    energy_quadcur = mrisComputeQuadraticCurvatureEnergy(mris);
    energy_repulse = mrisComputeRepulsiveEnergy(mris, weight_repulse, mht_v_current);
    energy_Tspring = mrisComputeTangentialSpringEnergy(mris) ;
    energy_Nspring = mrisComputeNormalSpringEnergy(mris);
    energy_Gspring = mrisComputeGaussianCurvatureSpringEnergy(mris, gaussian_norm);

    //****** Total Energy  ****//

    energy_old = weight_quadcur*energy_quadcur + weight_label*energy_label + \
                 + weight_repulse*energy_repulse + weight_Nspring*energy_Nspring \
                 + weight_Tspring*energy_Tspring +  weight_Gspring*energy_Gspring;

    mrisClearMomentum(mris);
    while ( !(ratio <= 0.0001 || t > 5000)) {
      /***reset the parameters***/
      energy_quadcur = 0;
      energy_label = 0;
      energy_Tspring = 0;
      energy_Nspring = 0;
      energy_repulse = 0;
      energy_Gspring = 0;
      mrisClearGradient(mris);
#if 1
      if (t <= 100)  step_size = 0.1;
      else if (t <= 200) step_size = 0.2;
      else if (t <= 1000) step_size = t/1000.0;
      else step_size = 1;
#else
      if (counter <= 500)  step_size = 0.1;
      else if (counter <= 1000) step_size = 0.2;
      else if (counter <= 5000) step_size = counter/5000.0;
      else step_size = 1;
#endif
      MHTfree(&mht_v_current);
      mht_v_current = MHTcreateVertexTable_Resolution(mris,CURRENT_VERTICES, 1.0f);
      mrisComputeQuadraticCurvatureTerm(mris, weight_quadcur);
      mrisComputeLabelTerm1(mris,weight_label,mri_smooth,mri_label,mri_orig);
      mrisComputeRepulsiveTerm(mris,weight_repulse,mht_v_current) ;
      mrisComputeNormalSpringTerm(mris,weight_Nspring);
      mrisComputeTangentialSpringTerm(mris, weight_Tspring);
      mriSspringTermWithGaussianCurvature(mris, gaussian_norm, weight_Gspring);

      //*****  Examine if movement is too large ****//
      mrisExaminemovelength(mris);

      MRISupdateSurface(mris);

      //******** Compute Energy  *****//
      MRISuseMeanCurvature(mris);
      energy_label = mrisComputeLabelEnergy(mris,mri_smooth,mri_label,mri_orig);
      energy_quadcur = mrisComputeQuadraticCurvatureEnergy(mris);
      energy_repulse = mrisComputeRepulsiveEnergy(mris, weight_repulse, mht_v_current);
      energy_Tspring = mrisComputeTangentialSpringEnergy(mris) ;
      energy_Nspring = mrisComputeNormalSpringEnergy(mris);
      energy_Gspring = mrisComputeGaussianCurvatureSpringEnergy(mris, gaussian_norm);

      //****** Total Energy  ****//
      //energy_new = weight_quadcur*energy_quadcur + weight_label*energy_label + energy_repulse
      //                + weight_Nspring*energy_Nspring + weight_Tspring*energy_Tspring;

      energy_new = weight_quadcur*energy_quadcur + weight_label*energy_label \
                   + weight_repulse*energy_repulse +  weight_Nspring*energy_Nspring \
                   + weight_Tspring*energy_Tspring +  weight_Gspring*energy_Gspring;

      if (energy_new == 0) ratio = 0;
      else ratio = fabs(energy_new - energy_old) / energy_old;

      if (write_iterations > 0) {
        if (((++counter) % write_iterations) == 0) {
          sprintf(ofname, "%s/%s/surf/movie/%s.firstrefined%3.3d", data_dir, argv[1], argv[2],counter/write_iterations);
          MRISwrite(mris, ofname) ;
        }
      }

      /*fprintf(stderr, "%dth iteration: ratio = %2.5f, energy_quadcur = %5.2f, energy_label = %5.2f,"\
          " energy_repulse = %5.2f, energy_Nspring = %5.2f, energy_Tspring = %5.2f\n",\
          t,ratio,energy_quadcur,energy_label,energy_repulse,energy_Nspring,energy_Tspring) ; */
      t++;
      energy_old = energy_new ;

    }

    MRIfree(&mri_smooth[0]);
    MRIfree(&mri_smooth[1]);
    MRIfree(&mri_smooth[2]);
    MRIfree(&mri_smooth[3]);
    MRIfree(&mri_smooth[4]);

    sigma -= 0.5;

  }


  if (smooth_spikes>0) {
    fprintf(stderr, "Smooth the spikes for %d times..............\n", smooth_spikes);
    i=0;
    while ( i<smooth_spikes ) {
      spikes = FindSpikes(mris,i);
      if (write_iterations > 0) {
        if (((++counter) % write_iterations) == 0) {
          sprintf(ofname, "%s/%s/surf/movie/%s.refined%3.3d", data_dir, argv[1], argv[2],counter/write_iterations);
          fprintf(stderr, "writing out reconstructed surface after %d iteration to %3.3d \n", i, counter/write_iterations );
          MRISwrite(mris, ofname) ;
          //measure the volume this new surface encloses and compare with the original//
        }
      }
      if (i<0) SmoothSpikes(mris, 3);
      else SmoothSpikes(mris, 2);
      fprintf(stderr, "Find %d spikes in %d iteration \n", spikes, i );
      i++;
    }
  }

  fprintf(stderr, "Gaussian smoothing for %d times...............\n", niteration);
  if (niteration>0) {
    MRISresetNeighborhoodSize(mris, 2) ;
    mrisClearMomentum(mris);
    mrisClearGradient(mris);
    MRIScomputeMetricProperties(mris) ;
    MRISstoreMetricProperties(mris) ;
    for ( i=0; i<niteration; i++ ) {
      MRIScomputeSecondFundamentalForm(mris) ;
      mriSspringTermWithGaussianCurvature(mris, gaussian_norm, 1) ;
      mrismomentumTimeStep(mris, 0.5, 1, 1, 0) ;
      mrisClearGradient(mris);
      if (write_iterations > 0) {
        if (((++counter) % write_iterations) == 0) {
          sprintf(ofname, "%s/%s/surf/movie/%s.refined%3.3d", data_dir, argv[1], argv[2],counter/write_iterations);
          MRISwrite(mris, ofname) ;
        }
      }
    }
  }

  fprintf(stderr, "Average Vertex Position for 1 times...........\n");
  MRISaverageVertexPositions(mris, 1);
  //write out final smoothing  result//
  if (write_iterations > 0) {
    counter = floor(counter/write_iterations)+1 ;
    sprintf(ofname, "%s/%s/surf/movie/%s.refined%3.3d", data_dir, argv[1], argv[2],counter);
    MRISwrite(mris, ofname) ;
  }

  sprintf(ofname, "%s/%s/surf/%s.%s", data_dir, argv[1], argv[2],suffix) ;
  fprintf (stderr,"writing refined surface to %s\n", ofname);
  MRISwrite(mris, ofname) ;
  MRISuseMeanCurvature(mris);
  sprintf(ofname, "%s/%s/surf/%s.%s.curv", data_dir, argv[1], argv[2],suffix) ;
  MRISwriteCurvature(mris, ofname) ;
  MRISfree(&mris);
  MRIfree(&mri_label[0]);
  MRIfree(&mri_label[1]);
  MRIfree(&mri_orig);
  MRIfree(&mri_label[2]);
  MRIfree(&mri_label[3]);
  MRIfree(&mri_label[4]);
  exit(0) ;
  return(0) ;
}

#if 0
static int
extractlabelvolume(MRI *mri_label[5], MRI *mri_orig) {
  int  xi, yi, zi;

  mri_label[0] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[1] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[2] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[3] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[4] = MRIalloc(256, 256, 256, MRI_FLOAT) ;

  for (xi=0; xi<mri_orig->depth; xi++)
    for (yi=0; yi<mri_orig->height; yi++)
      for (zi=0; zi<mri_orig->width; zi++) {
        if ( MRIvox(mri_orig,zi,yi,xi) >=1 && MRIvox(mri_orig,zi,yi,xi) <=4 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else if ( MRIvox(mri_orig,zi,yi,xi) >=5 && MRIvox(mri_orig,zi,yi,xi) <=8 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else if ( MRIvox(mri_orig,zi,yi,xi) == 13 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else if ( MRIvox(mri_orig,zi,yi,xi) == 14 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 1.0;
        }
      }
  return(NO_ERROR);

}

static int
mrisFindneighborlabel(MRI_SURFACE *mris, char surftype[10], MRI *mri_label[5], MRI *mri_orig) {
  int           vno, type=0, tt;
  VERTEX        *v;
  float         x, y, z, step;
  double        xw, yw, zw, val=0;

  for (tt=0; tt<4; tt++) {
    if ( !strcmp(surftype,surf[tt]) ) type = tt;
  }

  for ( vno=0; vno < mris->nvertices; vno++ ) {
    v = &mris->vertices[vno] ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    step = 0.5;
    table[0][vno] = type;

    while ( table[0][vno] == type && step <= 2 ) {
      MRIsurfaceRASToVoxel(mri_orig, v->x+step*v->nx, v->y+step*v->ny, v->z+step*v->nz, &xw, &yw, &zw) ;
      MRIsampleVolumeType(mri_orig, xw, yw, zw, &val, SAMPLE_NEAREST);
      if ( val >=1 && val <=4 )
        table[0][vno] = 0;
      else if ( val >=5 && val <=8 )
        table[0][vno] = 1;
      else if ( val == 13 )
        table[0][vno] = 2;
      else if ( val == 14 )
        table[0][vno] = 3;
      else
        table[0][vno] = 4;
      step +=0.25;
    }

#if 0
    MRIsurfaceRASToVoxel(mri_orig, v->x-0.5*v->nx, v->y-0.5*v->ny, v->z-0.5*v->nz, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_orig, xw, yw, zw, &val, SAMPLE_NEAREST);
    if ( val >=1 && val <=4 )
      table[1][vno] = 0;
    else if ( val >=5 && val <=8 )
      table[1][vno] = 1;
    else if ( val == 13 )
      table[1][vno] = 2;
    else if ( val == 14 )
      table[1][vno] = 3;
    else
      table[1][vno] = 4;
#else
    table[1][vno] = type;
#endif
  }
  return(NO_ERROR);
}


#else

static int
extractlabelvolume(MRI *mri_label[5], MRI *mri_orig) {
  int  xi, yi, zi;

  mri_label[0] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[1] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[2] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[3] = MRIalloc(256, 256, 256, MRI_FLOAT) ;
  mri_label[4] = MRIalloc(256, 256, 256, MRI_FLOAT) ;

  for (xi=0; xi<mri_orig->depth; xi++)
    for (yi=0; yi<mri_orig->height; yi++)
      for (zi=0; zi<mri_orig->width; zi++) {
        if ( MRIvox(mri_orig,zi,yi,xi) == 17 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else if ( MRIvox(mri_orig,zi,yi,xi) == 53 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else if ( MRIvox(mri_orig,zi,yi,xi) == 18 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else if ( MRIvox(mri_orig,zi,yi,xi) == 54 ) {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 1.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 0.0;
        } else {
          MRIFvox(mri_label[0],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[1],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[2],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[3],zi,yi,xi) = 0.0;
          MRIFvox(mri_label[4],zi,yi,xi) = 1.0;
        }
      }
  return(NO_ERROR);

}


static int
mrisFindneighborlabel(MRI_SURFACE *mris, char surftype[10], MRI *mri_label[5], MRI *mri_orig) {
  int           vno, type=0, tt;
  VERTEX        *v;
  float         x, y, z, step;
  double        xw, yw, zw, val=0;

  for (tt=0; tt<4; tt++) {
    if ( !strcmp(surftype,surf[tt]) ) type = tt;
  }

  for ( vno=0; vno < mris->nvertices; vno++ ) {
    v = &mris->vertices[vno] ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    step = 0.5;
    table[0][vno] = type;

    while ( table[0][vno] == type && step <= 2 ) {
      MRIsurfaceRASToVoxel(mri_orig, v->x+step*v->nx, v->y+step*v->ny, v->z+step*v->nz, &xw, &yw, &zw) ;
      MRIsampleVolumeType(mri_orig, xw, yw, zw, &val, SAMPLE_NEAREST);
      if ( val == 17 )
        table[0][vno] = 0;
      else if ( val == 53 )
        table[0][vno] = 1;
      else if ( val == 18 )
        table[0][vno] = 2;
      else if ( val == 54 )
        table[0][vno] = 3;
      else
        table[0][vno] = 4;
      step +=0.25;
    }

#if 1
    MRIsurfaceRASToVoxel(mri_orig, v->x-0.5*v->nx, v->y-0.5*v->ny, v->z-0.5*v->nz, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_orig, xw, yw, zw, &val, SAMPLE_NEAREST);
    if ( val ==17 )
      table[1][vno] = 0;
    else if ( val == 53 )
      table[1][vno] = 1;
    else if ( val == 18 )
      table[1][vno] = 2;
    else if ( val == 54 )
      table[1][vno] = 3;
    else
      table[1][vno] = 4;
#else
table[1][vno] = type;
#endif
  }
  return(NO_ERROR);
}

#endif



//////************ Used to calcuate the label term*************/////

static int
mrisComputeLabelTerm1(MRI_SURFACE *mris, double weight_label, MRI *mri_smooth[5],MRI *mri_label[5], MRI *mri_orig) {
  int     vno;
  VERTEX  *v ;
  float   x, y, z, dx=0, dy=0, dz=0, nx=0, ny=0, nz=0;
  double  xw, yw, zw, dn, xw1, yw1, zw1, outlabel=0, inlabel=0;

  if (FZERO(weight_label))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
#if 0
    x = v->x+1.5*v->nx ;
    y = v->y+1.5*v->ny ;
    z = v->z+1.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw1, &yw1, &zw1) ;
    x = v->x+0.5*v->nx ;
    y = v->y+0.5*v->ny ;
    z = v->z+0.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw, &yw, &zw) ;
    nx = xw1-xw ;
    ny = yw1-yw ;
    nz = zw1-zw ;
    MRIsampleVolumeDerivative(mri_smooth[table[0][vno]], xw, yw, zw, nx, ny, nz, &dn) ;

    MRIsurfaceRASToVoxel(mri_smooth[table[0][vno]], v->x+0.5*v->nx, v->y+0.5*v->ny, v->z+0.5*v->nz, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_smooth[table[0][vno]], xw, yw, zw, &outlabel, SAMPLE_TRILINEAR);

    dx = (1-outlabel)* v->nx * dn;
    dy = (1-outlabel)* v->ny * dn;
    dz = (1-outlabel)* v->nz * dn;

    x = v->x+0.5*v->nx ;
    y = v->y+0.5*v->ny ;
    z = v->z+0.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw1, &yw1, &zw1) ;
    x = v->x-0.5*v->nx ;
    y = v->y-0.5*v->ny ;
    z = v->z-0.5*v->nz ;
    MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw, &yw, &zw) ;
    nx = xw1-xw ;
    ny = yw1-yw ;
    nz = zw1-zw ;
    MRIsampleVolumeDerivative(mri_smooth[table[1][vno]], xw, yw, zw, nx, ny, nz, &dn) ;

    MRIsurfaceRASToVoxel(mri_smooth[table[1][vno]], v->x-0.5*v->nx, v->y-0.5*v->ny, v->z-0.5*v->nz, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_smooth[table[1][vno]], xw, yw, zw, &inlabel, SAMPLE_TRILINEAR);

    dx += (1-inlabel)* v->nx * dn;
    dy += (1-inlabel)* v->ny * dn;
    dz += (1-inlabel)* v->nz * dn;
#else
    x = v->x+v->nx ;
    y = v->y+v->ny ;
    z = v->z+v->nz ;
    MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw1, &yw1, &zw1) ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw, &yw, &zw) ;
    nx = xw1-xw ;
    ny = yw1-yw ;
    nz = zw1-zw ;

    MRIsampleVolumeDerivative(mri_smooth[table[0][vno]], xw, yw, zw, nx, ny, nz, &dn) ;
    MRIsurfaceRASToVoxel(mri_orig, x+0.5*v->nx, y+0.5*v->ny, z+0.5*v->nz, &xw, &yw, &zw);
    MRIsampleVolumeType(mri_label[table[0][vno]], xw, yw, zw, &outlabel, SAMPLE_NEAREST);

    dx = (1 - outlabel)* v->nx * dn;
    dy = (1 - outlabel)* v->ny * dn;
    dz = (1 - outlabel)* v->nz * dn;

    MRIsampleVolumeDerivative(mri_smooth[table[1][vno]], xw, yw, zw, nx, ny, nz, &dn) ;
    MRIsurfaceRASToVoxel(mri_orig, x-0.5*v->nx, y-0.5*v->ny, z-0.5*v->nz, &xw, &yw, &zw);
    MRIsampleVolumeType(mri_label[table[1][vno]], xw, yw, zw, &inlabel, SAMPLE_NEAREST);

    dx += (1 - inlabel)* v->nx * dn;
    dy += (1 - inlabel)* v->ny * dn;
    dz += (1 - inlabel)* v->nz * dn;

#endif

    v->dx += weight_label *dx;
    v->dy += weight_label *dy;
    v->dz += weight_label *dz;
  }
  return(NO_ERROR) ;
}

static int
mrisComputeNormalSpringTerm(MRI_SURFACE *mris, double l_spring) {
  int     vno, n, m ;
  float   sx, sy, sz, nx, ny, nz, nc, x, y, z ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    nx = vertex->nx ;
    ny = vertex->ny ;
    nz = vertex->nz ;
    x = vertex->x ;
    y = vertex->y ;
    z = vertex->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < vertext->vnum ; m++) {
      VERTEX const * const vn = &mris->vertices[vertext->v[m]] ;
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n>0) {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
    nc = sx*nx+sy*ny+sz*nz;   /* projection onto normal */
    sx = l_spring*nc*nx ;              /* move in normal direction */
    sy = l_spring*nc*ny ;
    sz = l_spring*nc*nz ;

    vertex->dx += sx ;
    vertex->dy += sy ;
    vertex->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring normal term:  (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }


  return(NO_ERROR) ;
}

static double
mrisComputeNormalSpringEnergy(MRI_SURFACE *mris) {
  int     vno, n ;
  double  area_scale, sse_spring, v_sse ;
  float   dx, dy, dz, x, y, z, nc, dist_sq ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  for (sse_spring = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;

    x = v->x ;
    y = v->y ;
    z = v->z ;

    for (v_sse = 0.0, n = 0 ; n < vt->vnum ; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      dx = vn->x - x ;
      dy = vn->y - y ;
      dz = vn->z - z ;
      nc = dx * v->nx + dy*v->ny + dz*v->nz ;
      dist_sq = nc*nc ;
      v_sse += dist_sq ;
    }
    sse_spring += area_scale * v_sse ;
  }
  return(sse_spring) ;
}

static int
mrisComputeTangentialSpringTerm(MRI_SURFACE *mris, double l_spring) {
  int     vno, n, m ;
  float   sx, sy, sz, x, y, z, nc ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    if (v->border && !v->neg)
      continue ;

    x = v->x ;
    y = v->y ;
    z = v->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < vt->vnum ; m++) {
      VERTEX const * const vn = &mris->vertices[vt->v[m]] ;
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
#if 0
    n = 4 ;  /* avg # of nearest neighbors */
#endif
    if (n>0) {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }

    nc = sx*v->nx+sy*v->ny+sz*v->nz;   /* projection onto normal */
    sx -= l_spring*nc*v->nx ;                   /* remove  normal component */
    sy -= l_spring*nc*v->ny ;
    sz -= l_spring*nc*v->nz;

    v->dx += sx ;
    v->dy += sy ;
    v->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring tangent term: (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }


  return(NO_ERROR) ;
}


static double
mrisComputeTangentialSpringEnergy(MRI_SURFACE *mris) {
  int     vno, n ;
  double  area_scale, sse_spring, v_sse ;
  float   dx, dy, dz, x, y, z, nc, dist_sq ;

#if METRIC_SCALE
  if (mris->patch)
    area_scale = 1.0 ;
  else
    area_scale = mris->orig_area / mris->total_area ;
#else
  area_scale = 1.0 ;
#endif

  for (sse_spring = 0.0, vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;

    x = v->x ;
    y = v->y ;
    z = v->z ;

    for (v_sse = 0.0, n = 0 ; n < vt->vnum ; n++) {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      dx = vn->x - x ;
      dy = vn->y - y ;
      dz = vn->z - z ;
      nc = dx * v->nx + dy*v->ny + dz*v->nz ;
      dx -= nc*v->nx ;
      dy -= nc*v->ny ;
      dz -= nc*v->nz ;
      dist_sq = dx*dx+dy*dy+dz*dz ;
      v_sse += dist_sq ;
    }
    sse_spring += area_scale * v_sse ;
  }
  return(sse_spring) ;
}


//************* For Calculating the Label Term Energy and curvature only *****//

static double
mrisComputeLabelEnergy(MRI_SURFACE *mris, MRI *mri_smooth[5],MRI *mri_label[5], MRI *mri_orig) {
  int           vno ;
  double        xw, yw, zw, inval=0, outval=0;
  float         x, y, z;
  //float         target_I = 0.5;
  VERTEX        *v;
  double        energy = 0;

  for (vno=0;vno<mris->nvertices;vno++) {
    v = &mris->vertices[vno];
    x = v->x;
    y = v->y;
    z = v->z;
#if 1
    MRIsurfaceRASToVoxel(mri_orig, x+0.5*v->nx, y+0.5*v->ny, z+0.5*v->nz, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_label[table[0][vno]], xw, yw, zw, &outval, SAMPLE_NEAREST);
    MRIsurfaceRASToVoxel(mri_orig, x-0.5*v->nx, y-0.5*v->ny, z-0.5*v->nz, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_label[table[1][vno]], xw, yw, zw, &inval, SAMPLE_NEAREST);
    //if ( outval!=1 || inval!=1)
    //fprintf(stderr, "vertex %d is not correct: outval=%f inval=%f \n", vno, outval, inval) ;
    energy += (1 - outval) * (1 - outval) + (1 - inval) * (1 - inval) ;
#else
    MRIsurfaceRASToVoxel(mri_orig, x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_label[table[0][vno]], xw, yw, zw, &outval, SAMPLE_NEAREST);
    //MRIsurfaceRASToVoxel(mri_smooth[table[1][vno]], x, y, z, &xw, &yw, &zw) ;
    MRIsampleVolumeType(mri_label[table[1][vno]], xw, yw, zw, &inval, SAMPLE_NEAREST);
    if ( outval!=target_O || inval!=target_I)
      fprintf(stderr, "vertex %d is not correct: outval=%f inval=%f \n", vno, outval, inval) ;
    energy += (target_O - outval) * (target_O - outval) + (target_I - inval) * (target_I - inval) ;
#endif
  }
  return(energy) ;
}

#if 0
//////*********** Revised a little bit to fit this problem  ****************///
static int
mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv) {
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n;
  VERTEX   *v, *vn ;
  float    ui, vi, rsqx, rsqy, a, b, c;

  if (FZERO(l_curv))
    return(NO_ERROR) ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 3, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;

    /*For each vertex, compute the movement it should make*/

    for (n = 0 ; n < v->vtotal ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ;
      vi = V3_DOT(v_e2, v_nbr) ;
      rsqx = ui*ui;
      rsqy = vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsqx ;
      *MATRIX_RELT(m_R, n+1, 2) = rsqy ;
      *MATRIX_RELT(m_R, n+1, 3) = 1 ;
    }
    //m_R_inv = MatrixSVDInverse(m_R, NULL) ;
    m_R_inv = MatrixPseudoInverse(m_R, NULL);
    if (!m_R_inv) {
      MatrixFree(&m_R) ;
      VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
    c = VECTOR_ELT(v_A, 3) ;

    c *= l_curv ;
    v->dx = +c * v->nx ;
    v->dy = +c * v->ny ;
    v->dz = +c * v->nz ;

    if (vno == Gdiag_no)
      fprintf(stdout, "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "a=%2.2f, b=%2.2f, c=%2.2f\n",
              vno, c*v->nx, c*v->ny, c*v->nz, a, b, c) ;

    MatrixFree(&m_R) ;
    VectorFree(&v_Y) ;
    MatrixFree(&m_R_inv) ;

  }

  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_nbr) ;
  VectorFree(&v_A) ;
  return(NO_ERROR) ;
}




//////*********** Revised a little bit to fit this problem  ****************///
static double
mrisComputeQuadraticCurvatureEnergy(MRI_SURFACE *mris) {
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n;
  VERTEX   *v, *vn ;
  float    ui, vi, rsqx, rsqy, a, b, c;
  double   energy = 0;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(v->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(v->vtotal, 3, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;

    /*For each vertex, compute the movement it should make*/

    for (n = 0 ; n < v->vtotal ; n++)  /* build data matrices */
    {
      vn = &mris->vertices[v->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ;
      vi = V3_DOT(v_e2, v_nbr) ;
      rsqx = ui*ui;
      rsqy = vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsqx ;
      *MATRIX_RELT(m_R, n+1, 2) = rsqy ;
      *MATRIX_RELT(m_R, n+1, 3) = 1 ;
    }
    //m_R_inv = MatrixSVDInverse(m_R, NULL) ;
    m_R_inv = MatrixPseudoInverse(m_R, NULL);
    if (!m_R_inv) {
      MatrixFree(&m_R) ;
      VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
    c = VECTOR_ELT(v_A, 3) ;
    energy += c*c;

    if (vno == Gdiag_no)
      fprintf(stdout, "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "a=%2.2f, b=%2.2f, c=%2.2f\n",
              vno, c*v->nx, c*v->ny, c*v->nz, a, b, c) ;

    MatrixFree(&m_R) ;
    VectorFree(&v_Y) ;
    MatrixFree(&m_R_inv) ;

  }

  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_nbr) ;
  VectorFree(&v_A) ;
  return(energy) ;
}

#else

/*-----------------------------------------------------
Parameters:

Returns value:

Description
Fit a 1-d quadratic to the surface locally and move the
vertex in the normal direction to improve the fit.
------------------------------------------------------*/
static int
mrisComputeQuadraticCurvatureTerm(MRI_SURFACE *mris, double l_curv) {
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n ;
  float    ui, vi, rsq, a, b ;

  if (FZERO(l_curv))
    return(NO_ERROR) ;

  my_mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(2, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(vt->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(vt->vtotal, 2, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < vt->vtotal ; n++)  /* build data matrices */
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ;
      vi = V3_DOT(v_e2, v_nbr) ;
      rsq = ui*ui + vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsq ;
      *MATRIX_RELT(m_R, n+1, 2) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    if (!m_R_inv) {
      MatrixFree(&m_R) ;
      VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
    b *= l_curv ;
    v->dx += b * v->nx ;
    v->dy += b * v->ny ;
    v->dz += b * v->nz ;

    if (vno == Gdiag_no)
      fprintf(stdout, "v %d curvature term:      (%2.3f, %2.3f, %2.3f), "
              "a=%2.2f, b=%2.1f\n",
              vno, b*v->nx, b*v->ny, b*v->nz, a, b) ;
    MatrixFree(&m_R) ;
    VectorFree(&v_Y) ;
    MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_nbr) ;
  VectorFree(&v_A) ;
  return(NO_ERROR) ;
}
/*-----------------------------------------------------
Parameters:

Returns value:

Description
Fit a 1-d quadratic to the surface locally and move the
vertex in the normal direction to improve the fit.
------------------------------------------------------*/
static double
mrisComputeQuadraticCurvatureEnergy(MRI_SURFACE *mris) {
  MATRIX   *m_R, *m_R_inv ;
  VECTOR   *v_Y, *v_A, *v_n, *v_e1, *v_e2, *v_nbr ;
  int      vno, n ;
  float    ui, vi, rsq, a, b ;
  double   sse = 0.0 ;


  my_mrisComputeTangentPlanes(mris) ;
  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_A = VectorAlloc(2, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v_nbr = VectorAlloc(3, MATRIX_REAL) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    v_Y = VectorAlloc(vt->vtotal, MATRIX_REAL) ;    /* heights above TpS */
    m_R = MatrixAlloc(vt->vtotal, 2, MATRIX_REAL) ; /* radial distances */
    VECTOR_LOAD(v_n, v->nx, v->ny, v->nz) ;
    VECTOR_LOAD(v_e1, v->e1x, v->e1y, v->e1z) ;
    VECTOR_LOAD(v_e2, v->e2x, v->e2y, v->e2z) ;
    for (n = 0 ; n < vt->vtotal ; n++)  /* build data matrices */
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      VERTEX_EDGE(v_nbr, v, vn) ;
      VECTOR_ELT(v_Y, n+1) = V3_DOT(v_nbr, v_n) ;
      ui = V3_DOT(v_e1, v_nbr) ;
      vi = V3_DOT(v_e2, v_nbr) ;
      rsq = ui*ui + vi*vi ;
      *MATRIX_RELT(m_R, n+1, 1) = rsq ;
      *MATRIX_RELT(m_R, n+1, 2) = 1 ;
    }
    m_R_inv = MatrixPseudoInverse(m_R, NULL) ;
    if (!m_R_inv) {
      MatrixFree(&m_R) ;
      VectorFree(&v_Y) ;
      continue ;
    }
    v_A = MatrixMultiply(m_R_inv, v_Y, v_A) ;
    a = VECTOR_ELT(v_A, 1) ;
    b = VECTOR_ELT(v_A, 2) ;
    sse += b*b ;
    if (vno == Gdiag_no)
      printf("v %d: curvature sse %2.2f\n", vno, b*b) ;
    MatrixFree(&m_R) ;
    VectorFree(&v_Y) ;
    MatrixFree(&m_R_inv) ;
  }

  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  VectorFree(&v_nbr) ;
  VectorFree(&v_A) ;
  return(sse) ;
}

#endif

static int
mrisExaminemovelength(MRI_SURFACE *mris) {
  int           k, overnum=0, over_num=0;
  //float         dx, dy, dz, x, y, z, A, B, C, f;
  float         mag, th ;

  for (k=0;k<mris->nvertices;k++) {
    VERTEX* v = &mris->vertices[k];

    /********** First, restrict vertex movement in each step by threshold *******/
    mag = sqrt ((v->dx)*(v->dx) + (v->dy)*(v->dy)+ (v->dz)*(v->dz));
    th = threshold /step_size ;
    if ( mag > th ) {
      //fprintf(stdout, "%dth :  (%2.3f, %2.3f, %2.3f), movement:"
      //   "dx=%2.2f, dy=%2.2f, dz=%2.2f\n",
      //   k, v->x, v->y, v->z, v->dx, v->dy, v->dz) ;
      over_num++;
      v->dx = v->dx * th / mag;
      v->dy = v->dy * th / mag;
      v->dz = v->dz * th / mag;
    }
#if 0
    /********** Then, Constrain total movement to be within 1 mm ******/
    x = v->x - v->origx;
    y = v->y - v->origy;
    z = v->z - v->origz;   //the movement it has made in past
    dx = step_size*v->dx;
    dy = step_size*v->dy;
    dz = step_size*v->dz; // the movement it's going to make
    mag = MAX_mag*MAX_mag;
    if ( ((x+dx)*(x+dx) + (y+dy)*(y+dy) + (z+dz)*(z+dz)) > mag ) // if total movement is larger than 1mm
    {
      A = dx*dx + dy*dy + dz*dz;
      B = 2*x*dx + 2*y*dy + 2*z*dz;
      C =  x*x + y*y + z*z - mag;
      if ( C > 0) C = 0;
      f = ( sqrt( B*B-4*A*C ) - B) / 2 / A;
      dx *= f;
      dy *= f;
      dz *= f;
      overnum++;
    }
    v->x += dx ;
    v->y += dy ;
    v->z += dz ;                  //then make it to 1mm
#else
    //v->dx = step_size*v->dx; v->dy = step_size*v->dy; v->dz = step_size*v->dz;
    v->odx = step_size*(v->dx + 0.2*v->odx) ;
    v->ody = step_size*(v->dy + 0.2*v->ody) ;
    v->odz = step_size*(v->dz + 0.2*v->odz) ;
    mag =
      sqrt(v->odx*v->odx +
           v->ody*v->ody +
           v->odz*v->odz) ;
    if (mag > MAX_MOMENTUM_MM) /* don't let step get too big */
    {
      mag =  MAX_MOMENTUM_MM/ mag ;
      v->odx *= mag ;
      v->ody *= mag ;
      v->odz *= mag ;
      overnum++;
      //fprintf(stdout, "%dth vertex :  (%2.3f, %2.3f, %2.3f), movement:"
      //   "dx=%2.2f, dy=%2.2f, dz=%2.2f\n",
      //   k,v->dx, v->dy, v->dz, v->odx, v->ody, v->odz ) ;
    }
    MRISsetXYZ(mris,k,
      v->x + v->odx,
      v->y + v->ody,
      v->z + v->odz);
#endif
  }

  //fprintf(stderr,"There are %d vertex make large movement in a single step\n", over_num);
  //fprintf(stderr,"%d vertex move out of region in %dth iteration\n", overnum, t+1);
  return (NO_ERROR);
}

static int
mrisComputeRepulsiveTerm(MRI_SURFACE *mris, double l_repulse, MHT *mht) {
  int     vno, num, min_vno, i, n ;
  float   dist, dx, dy, dz, x, y, z, sx, sy, sz, min_d, min_scale, norm ;
  double  scale=0 ;

  if (FZERO(l_repulse))
    return(NO_ERROR) ;

  min_d = 100000.0 ;
  min_scale = 1.0 ;
  min_vno = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    MHBT *bucket = MHTacqBucket(mht, x, y, z) ;
    if (!bucket)
      continue ;
    sx = sy = sz = 0.0 ;
    MHB     *bin ;
    for (bin = bucket->bins, num = i = 0 ; i < bucket->nused ; i++, bin++) {
      if (bin->fno == vno)
        continue ;  /* don't be repelled by myself */
      for (n = 0 ; n < vt->vtotal ; n++)
        if (vt->v[n] == bin->fno)
          break ;
      if (n < vt->vtotal)   /* don't be repelled by a neighbor */
        continue ;
      VERTEX const * const vn = &mris->vertices[bin->fno] ;
      if (!vn->ripflag) {
        dx = vn->x - x ;
        dy = vn->y - y ;
        dz = vn->z - z ;
        dist = sqrt(dx*dx+dy*dy+dz*dz) + REPULSE_E ;
        scale = -4*REPULSE_K / (dist*dist*dist*dist*dist*dist*dist) ;  /* ^-7 */
        if (vno == Gdiag_no) {
          if (dist-REPULSE_E < min_d) {
            min_vno = bin->fno ;
            min_d = dist-REPULSE_E ;
            min_scale = scale ;
          }
        }
        norm = sqrt(dx*dx+dy*dy+dz*dz) ;
        dx /= norm ;
        dy /= norm ;
        dz /= norm ;
        sx += scale * dx ;
        sy += scale * dy ;
        sz += scale*dz ;
        num++ ;
      }
    }
    if (num) {
      scale = l_repulse / (double)num ;
      sx *= scale ;
      sy *= scale ;
      sz *= scale ;
    }
    v->dx += sx ;
    v->dy += sy ;
    v->dz += sz ;
    if ((vno == Gdiag_no) && min_d < 1000) {
      VERTEX const * const vn = &mris->vertices[min_vno] ;
      dx = x - vn->x ;
      dy = y - vn->y ;
      dz = z - vn->z ;

      fprintf(stdout, "v %d self repulse term:   (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
      fprintf(stdout, "min_dist @ %d = %2.2f, scale = %2.1f\n",
              min_vno, min_d, min_scale) ;
    }
    MHTrelBucket(&bucket);
  }
  return(NO_ERROR) ;
}


static double
mrisComputeRepulsiveEnergy(MRI_SURFACE *mris, double l_repulse, MHT *mht) {
  int     vno, num, min_vno, i, n ;
  float   dist, dx, dy, dz, x, y, z, min_d ;
  double  sse_repulse, v_sse ;

  if (FZERO(l_repulse))
    return(NO_ERROR) ;

  min_d = 1000.0 ;
  min_vno = 0 ;
  for (sse_repulse = vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX          const * const v  = &mris->vertices         [vno];
    if (v->ripflag)
      continue ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    MHBT* bucket = MHTacqBucket(mht, x, y, z) ;
    if (!bucket)
      continue ;
    MHB *bin ;
    for (v_sse = 0.0, bin = bucket->bins, num = i = 0 ; i < bucket->nused ; i++, bin++) {
      if (bin->fno == vno)
        continue ;  /* don't be repelled by myself */
      for (n = 0 ; n < vt->vtotal ; n++)
        if (vt->v[n] == bin->fno)
          break ;
      if (n < vt->vtotal)   /* don't be repelled by a neighbor */
        continue ;
      VERTEX const * const vn = &mris->vertices[bin->fno] ;
      if (!vn->ripflag) {
        dx = vn->x - x ;
        dy = vn->y - y ;
        dz = vn->z - z ;
        dist = sqrt(dx*dx+dy*dy+dz*dz) + REPULSE_E ;
        if (vno == Gdiag_no) {
          if (dist-REPULSE_E < min_d) {
            min_vno = bin->fno ;
            min_d = dist-REPULSE_E ;
          }
        }
        dist = dist*dist*dist ;
        dist *= dist ; /* dist^6 */
        v_sse += REPULSE_K / dist ;
      }
    }
    sse_repulse += v_sse ;

    if (vno == Gdiag_no && !FZERO(v_sse)) {
      printf("v %d: repulse sse:    min_dist=%2.4f, v_sse %2.4f\n", vno,
             min_d, v_sse) ;
    }
    MHTrelBucket(&bucket);
  }
  return(sse_repulse) ;
}

static int
mriSspringTermWithGaussianCurvature(MRI_SURFACE *mris, double gaussian_norm, double l_spring) {
  int     vno, n, m ;
  float   sx, sy, sz, x, y, z, scale ;

  if (FZERO(l_spring))
    return(NO_ERROR) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX                * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = vertex->x ;
    y = vertex->y ;
    z = vertex->z ;

    sx = sy = sz = 0.0 ;
    n=0;
    for (m = 0 ; m < vertext->vnum ; m++) {
      VERTEX const * const vn = &mris->vertices[vertext->v[m]] ;
      if (!vn->ripflag) {
        sx += vn->x - x;
        sy += vn->y - y;
        sz += vn->z - z;
        n++;
      }
    }
    if (n>0) {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
    }
    scale = pow(vertex->K, gaussian_norm) ;
    if (scale > 1)
      scale = 1 ;
    scale *= l_spring ;
    sx *= scale ;              /* move in normal direction */
    sy *= scale ;
    sz *= scale ;

    vertex->dx += sx ;
    vertex->dy += sy ;
    vertex->dz += sz ;
    if (vno == Gdiag_no)
      fprintf(stdout, "v %d spring normal term:  (%2.3f, %2.3f, %2.3f)\n",
              vno, sx, sy, sz) ;
  }

  return(NO_ERROR) ;
}

static double
mrisComputeGaussianCurvatureSpringEnergy(MRI_SURFACE *mris, double gaussian_norm ) {
  int     vno, m ;
  float   sx, sy, sz, x, y, z, scale ;
  double  sse_spring, v_sse;


  for (sse_spring = 0, vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
    VERTEX          const * const vertex  = &mris->vertices         [vno];
    if (vertex->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    x = vertex->x ;
    y = vertex->y ;
    z = vertex->z ;

    sx = sy = sz = 0.0 ;

    for (v_sse = 0,m = 0 ; m < vertext->vnum ; m++) {
      VERTEX const * const vn = &mris->vertices[vertext->v[m]] ;
      if (!vn->ripflag) {
        sx = vn->x - x;
        sy = vn->y - y;
        sz = vn->z - z;
        v_sse += sx*sx + sy*sy + sz*sz ;
      }
    }
    scale = pow(vertex->K, gaussian_norm) ;
    if (scale > 1)
      scale = 1 ;

    v_sse *= scale ;
    sse_spring += v_sse;

  }

  return(sse_spring) ;
}

static int
mrisClearMomentum(MRI_SURFACE *mris) {
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->odx = 0 ;
    v->ody = 0 ;
    v->odz = 0 ;
  }
  return(NO_ERROR) ;
}

static int
mrisClearGradient(MRI_SURFACE *mris) {
  int     vno, nvertices ;
  VERTEX  *v ;

  nvertices = mris->nvertices ;
  for (vno = 0 ; vno < nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->dx = 0 ;
    v->dy = 0 ;
    v->dz = 0 ;
  }
  return(NO_ERROR) ;
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
  if (!stricmp(option, "an option")) {}
  else switch (toupper(*option)) {
    case 'A':
      all_flag = 1 ;
      printf("tessellating the surface of all voxels with different labels\n") ;
      break ;
    case 'G':
      gaussian_norm = atof(argv[2]) ;
      printf("using Gaussian curvature smoothing with norm %2.2f\n", gaussian_norm) ;
      nargs = 1 ;
      break ;
    case 'W':
      write_iterations = atoi(argv[2]) ;
      printf("writing out snapshots every %d iterations\n", write_iterations) ;
      nargs = 1 ;
      break ;
    case 'N':
      niteration = atoi(argv[2]) ;
      printf("Gaussian Curvature Smoothing %d iterations\n", niteration) ;
      nargs = 1 ;
      break ;
    case 'C':
      weight_Gspring = atof (argv[2]) ;
      printf("changing weighting of Gaussian Curvature Spring term to %2.2f\n", weight_Gspring) ;
      nargs = 1 ;
      break ;
    case 'S':
      suffix = argv[2] ;
      printf("changing output surface file name to %s\n", suffix) ;
      nargs = 1 ;
      break ;
    case 'L':
      labelvolume = argv[2] ;
      printf("changing input segmentation volume to %s\n", labelvolume) ;
      nargs = 1 ;
      break ;
    case 'O':
      orig_name = argv[2] ;
      printf("changing original input surface to %s\n", orig_name) ;
      nargs = 1 ;
      break ;
    case 'P':
      weight_quadcur = atof(argv[2]);
      weight_label = atof(argv[3]);
      weight_repulse = atof(argv[4]);
      weight_Nspring = atof(argv[5]);
      weight_Tspring = atof(argv[6]);
      printf("changing term weighting to %2.2f %2.2f %2.2f %2.2f %2.2f\n", weight_quadcur, weight_label, weight_repulse,weight_Nspring, weight_Tspring) ;
      nargs = 5 ;
      break ;
    case 'M':
      smooth_spikes = atoi(argv[2]) ;
      printf("Spike Smoothing %d iterations\n", smooth_spikes) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      printf("Usage: %s <input Subject> <label> \n",Progname);
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr, "usage: %s [options] <subject name> <structure: RA LA RH LH>\n",
          Progname) ;
}


static double
mrismomentumTimeStep(MRI_SURFACE *mris, float momentum, float dt, float tol,
                     float n_averages) {
  double  delta_t, mag ;
  int     vno ;
  VERTEX  *v ;
#if 0
  double  max_delta ;
  float   dx, dy, dz ;
#endif

  delta_t = dt * sqrt((double)n_averages+1.0) ;
  ;


#if 0
  /* find the largest delta, and scale the gradient by it */
  max_delta = 0.0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dx = v->dx ;
    dy = v->dy ;
    dz = v->dz ;
    mag = sqrt(dx*dx+dy*dy+dz*dz) ;
    if (mag > max_delta)
      max_delta = mag ;
  }
  if (FZERO(max_delta))
    max_delta = tol ;

  if (delta_t > MAX_MOMENTUM_MM / max_delta)   /* no bigger than 1mm */
    delta_t = MAX_MOMENTUM_MM / max_delta ;
#endif

  /* take a step in the gradient direction modulated by momentum */
  if (mris->status == MRIS_RIGID_BODY) {
  
    mris->da = delta_t * mris->alpha + momentum * mris->da ;
    mris->db = delta_t * mris->beta  + momentum * mris->db ;
    mris->dg = delta_t * mris->gamma + momentum * mris->dg ;
    MRISrotate(mris, mris, mris->da, mris->db, mris->dg) ;
    
  } else {
    for (vno = 0 ; vno < mris->nvertices ; vno++) {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (vno == Gdiag_no)
        DiagBreak() ;
      v->odx = delta_t * v->dx + momentum*v->odx ;
      v->ody = delta_t * v->dy + momentum*v->ody ;
      v->odz = delta_t * v->dz + momentum*v->odz ;
      mag =
        sqrt(v->odx*v->odx +
             v->ody*v->ody +
             v->odz*v->odz) ;
      if (mag > MAX_MOMENTUM_MM) /* don't let step get too big */
      {
        mag = MAX_MOMENTUM_MM / mag ;
        v->odx *= mag ;
        v->ody *= mag ;
        v->odz *= mag ;
      }
      if (vno == Gdiag_no) {
        float dist, dot, dx, dy, dz ;

        dx = v->x - v->origx ;
        dy = v->y - v->origy ;
        dz = v->z - v->origz ;
        dist = sqrt(dx*dx+dy*dy+dz*dz) ;
        dot = dx*v->nx + dy*v->ny + dz*v->nz ;
        fprintf(stdout, "moving v %d by (%2.2f, %2.2f, %2.2f) dot=%2.2f-->"
                "(%2.1f, %2.1f, %2.1f)\n", vno, v->odx, v->ody, v->odz,
                v->odx*v->nx+v->ody*v->ny+v->odz*v->nz,
                v->x, v->y, v->z) ;
        fprintf(stdout, "n = (%2.1f,%2.1f,%2.1f), total dist=%2.3f, total dot = %2.3f\n",
                v->nx, v->ny, v->nz, dist, dot) ;
      }
      MRISsetXYZ(mris,vno,
        v->x + v->odx,
        v->y + v->ody,
        v->z + v->odz);
    }
  }
  
  my_mrisProjectSurface(mris) ;
  return(delta_t) ;
}

static int
my_mrisProjectSurface(MRI_SURFACE *mris) {

  /*  MRISupdateSurface(mris) ;*/
  switch (mris->status) {
  case MRIS_PARAMETERIZED_SPHERE:
    MRISprojectOntoSphere(mris, mris, mris->radius) ;
    break ;
  case MRIS_SPHERE:
    MRISprojectOntoSphere(mris, mris, mris->radius) ;
    break ;
  case MRIS_ELLIPSOID:
    MRISprojectOntoEllipsoid(mris, mris, 0.0f, 0.0f, 0.0f);
    break ;
  case PROJECT_PLANE:
    /*    mrisOrientPlane(mris) ;*/
    break ;
  case MRIS_RIGID_BODY:
    /*    MRISprojectOntoSphere(mris, mris, mris->radius) ;*/
    mris->status = MRIS_RIGID_BODY ;
    break ;
  default:
    break ;
  }
  return(NO_ERROR) ;
}


static int
my_mrisComputeTangentPlanes(MRI_SURFACE *mris) {
  VECTOR  *v_n, *v_e1, *v_e2, *v ;
  int     vno ;
  VERTEX  *vertex ;

  v_n = VectorAlloc(3, MATRIX_REAL) ;
  v_e1 = VectorAlloc(3, MATRIX_REAL) ;
  v_e2 = VectorAlloc(3, MATRIX_REAL) ;
  v = VectorAlloc(3, MATRIX_REAL) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    vertex = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    VECTOR_LOAD(v_n, vertex->nx, vertex->ny, vertex->nz) ;
    /* now find some other non-parallel vector */
#if 0
    if (!FZERO(vertex->nx) || !FZERO(vertex->ny)) {
      VECTOR_LOAD(v, 0.0, 0.0, 1.0) ;
    } else {
      VECTOR_LOAD(v, 0.0, 1.0, 0.0) ;
    }
#else
    VECTOR_LOAD(v, vertex->ny, vertex->nz, vertex->nx) ;
#endif
    V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    if (FZERO(V3_LEN(v_e1)))  /* happened to pick a parallel vector */
    {
      VECTOR_LOAD(v, vertex->ny, -vertex->nz, vertex->nx) ;
      V3_CROSS_PRODUCT(v_n, v, v_e1) ;
    }

    if (FZERO(V3_LEN(v_e1)) && DIAG_VERBOSE_ON)  /* happened to pick a parallel vector */
      fprintf(stderr, "vertex %d: degenerate tangent plane\n", vno) ;
    V3_CROSS_PRODUCT(v_n, v_e1, v_e2) ;
    V3_NORMALIZE(v_e1, v_e1) ;
    V3_NORMALIZE(v_e2, v_e2) ;
    vertex->e1x = V3_X(v_e1) ;
    vertex->e2x = V3_X(v_e2) ;
    vertex->e1y = V3_Y(v_e1) ;
    vertex->e2y = V3_Y(v_e2) ;
    vertex->e1z = V3_Z(v_e1) ;
    vertex->e2z = V3_Z(v_e2) ;
  }

  VectorFree(&v) ;
  VectorFree(&v_n) ;
  VectorFree(&v_e1) ;
  VectorFree(&v_e2) ;
  return(NO_ERROR) ;
}


static int
FindSpikes(MRI_SURFACE *mris, int iter) {
  int      alarm=0, step=4;

  MRISupdateSurface(mris);
  MRISuseMeanCurvature(mris);

  if ( iter<= 50 && iter%5==0) MRISresetNeighborhoodSize(mris, 2) ;
  else MRISresetNeighborhoodSize(mris, 1) ;

  if (iter<=100) step = 4;
  else step =  9;

  MRISclearMarks(mris);

  int counter=0 ;

  int vno;
  for (vno = 0 ; vno < mris->nvertices ; vno++) {
    VERTEX* vertex = &mris->vertices[vno] ;

    if ( (fabs(vertex->curv)>=3) || (vertex->K*vertex->K>=4) || fabs(vertex->k1)>=3 || fabs(vertex->K)<0.01 ) {
      if (vertex->marked == 0) {
        counter++;
        vertex->marked = 1;
        if (iter==0) {
          MRISsetXYZ(mris, vno, vertex->origx, vertex->origy, vertex->origz);
        }
      }
    }

    if ( iter<= smooth_spikes/1.25 && iter%5==0) {
      if ( (fabs(vertex->curv)>=4) || (vertex->K*vertex->K>=step) ||fabs(vertex->k1)>=4 || fabs(vertex->K)<0.01 ) {
        float weight = vertex->k1/100;
        if (vertex->k1 > 0 && vertex->K < 0) weight = -weight;
        if ( fabs(weight)>0.4) weight = weight/fabs(weight)/2.5;
        
        MRISsetXYZ(mris, vno,
          vertex->x - weight*vertex->nx,
          vertex->y - weight*vertex->ny,
          vertex->z - weight*vertex->nz);
          
        alarm++;
      }
    }
  }
  
  MRISupdateSurface(mris);
  MRISuseMeanCurvature(mris);
  
  fprintf(stderr, "Marked %d vertices after %d iterations \n", counter, iter);
  return(alarm);
}

static int
SmoothSpikes(MRI_SURFACE *mris, int niter) {

  if (FZERO(niter))
    return(NO_ERROR) ;

  int const nvertices = mris->nvertices;
  
  float *tx, *ty, *tz;
  MRISmemalignNFloats(nvertices, &tx, &ty, &tz);

  float* px, *py, *pz;
  MRISexportXYZ(mris, &px,&py,&pz);

  int i;
  for (i=0; i<niter; i++) {

    int  vno;
    for (vno = 0 ; vno < nvertices ; vno++) {
      VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
      VERTEX                * const vertex  = &mris->vertices         [vno];

      if (vertex->marked != 1) continue;
       
      float x = px[vno];
      float y = py[vno];
      float z = pz[vno];
      
      int num = 1;  /* account for central vertex */
      int n;
      for (n = 0 ; n < vertext->vtotal ; n++) {
        int const vno2 = vertext->v[n];
        VERTEX const * const vn = &mris->vertices[vno2] ;
        if (vn->ripflag)
          continue ;
        num++ ;
        x += px[vno2] ;
        y += py[vno2] ;
        z += pz[vno2] ;
      }
      
      tx[vno] = x / (float)(num) ;
      ty[vno] = y / (float)(num) ;
      tz[vno] = z / (float)(num) ;
    }

    for (vno = 0 ; vno < nvertices ; vno++) {
      VERTEX * const vertex = &mris->vertices[vno] ;
      if (vertex->marked == 1) {
        px[vno] = tx[vno] ;
        py[vno] = ty[vno] ;
        pz[vno] = tz[vno] ;
      }
    }

    int j;
    for (j=0; j<2; j++) {
      for (vno = 0 ; vno < nvertices ; vno++) {
        VERTEX_TOPOLOGY const * const vertext = &mris->vertices_topology[vno];
        VERTEX                * const vertex  = &mris->vertices         [vno];
        
        if ( (fabs(vertex->curv)>=4) || (vertex->K*vertex->K>=4) ||fabs(vertex->k1)>=4 || fabs(vertex->K)<0.01 ) {
          int num = 1;  /* account for central vertex */
          float x = px[vno] ;
          float y = py[vno] ;
          float z = pz[vno] ;
          
          int n;
          for (n = 0 ; n < vertext->vtotal ; n++) {
            int const vno2 = vertext->v[n];
            VERTEX const * const vn = &mris->vertices[vno2] ;
            if (vn->ripflag)
              continue ;
            num++ ;
            x += px[vno2] ;
            y += py[vno2] ;
            z += pz[vno2] ;
          }
          
          tx[vno] = x / (float)(num) ;
          ty[vno] = y / (float)(num) ;
          tz[vno] = z / (float)(num) ;
        }
      }
      
      for (vno = 0 ; vno < nvertices ; vno++) {
        VERTEX * const vertex = &mris->vertices[vno] ;
        if ( (fabs(vertex->curv)>=4) || (vertex->K*vertex->K>=4) ||fabs(vertex->k1)>=4 || fabs(vertex->K)<0.01 ) {
          px[vno] = tx[vno] ;
          py[vno] = ty[vno] ;
          pz[vno] = tz[vno] ;
        }
      }
    }
  }

  MRISimportXYZ(mris,  px, py, pz);

  freeAndNULL(px);
  freeAndNULL(py);
  freeAndNULL(pz);

  freeAndNULL(tx);
  freeAndNULL(ty);
  freeAndNULL(tz);

  return(NO_ERROR);
}
