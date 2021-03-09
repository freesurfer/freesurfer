/**
 * @brief computes distance maps for subcortical structures
 *
 * compute distance maps for amygdala, hippocampus, pallidum, putamen,
 * caudate, lateral ventricle, and layer IV gray
 */
/*
 * Original Author: Bruce Fischl
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
#include "version.h"
#include "timer.h"
#include "fastmarching.h"


static char *aseg_fname=NULL;

static const char *FRAME_FIELD_NAMES[]=  /* order correspond to 
                                            macros defined in mrisurf.h */
{
  NULL,
  "sulc",
  NULL, /* curvature directly computed */
  GRAYMID_NAME,
  T1MID_NAME,
  T2MID_NAME,
  PDMID_NAME,
  AMYGDALA_DIST_NAME,
  HIPPOCAMPUS_DIST_NAME,
  PALLIDUM_DIST_NAME,
  PUTAMEN_DIST_NAME,
  CAUDATE_DIST_NAME,
  LAT_VENTRICLE_DIST_NAME,
  INF_LAT_VENTRICLE_DIST_NAME,
};

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;

static int labels[50];
static int nlabels=0;
static int navgs=0;
static float fdistance=20.0f;
static int mode = 1 ;

static char subjects_dir[STRLEN] ;

static void mrisExtractMRIvalues(MRIS * mris,
                                 MRI * mri,
                                 MRI *mri_distance,
                                 float distance,
                                 int mode) {
  int n;
  double xw,yw,zw,xv,yv,zv,val;
  VERTEX *v;

  MRISclearCurvature(mris);
  for (n=0;n<mris->nvertices;n++) {
    v = &mris->vertices[n] ;
    xw = v->x ;
    yw = v->y ;
    zw = v->z ;
    MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    MRIsampleVolumeType(mri_distance, xv, yv, zv, &val, SAMPLE_NEAREST) ;
    val=MIN(distance,MAX(-distance,val));
    if (mode==1)
      val=MAX(0,val);
    if (mode==2)
      val=MIN(val,0);
    v->curv=val;
  }
}

#define MAXIMUM_DISTANCE 10.0

static void mrisProcessDistanceValues(MRIS *mris) {
  int n;
  VERTEX *v;

  for (n=0;n<mris->nvertices;n++) {
    v = &mris->vertices[n] ;

    v->curv=MIN(MAXIMUM_DISTANCE,
                MAX(0,(MAXIMUM_DISTANCE-v->curv)))/MAXIMUM_DISTANCE;
  }
}


static void mrisExtractMidGrayValues(MRIS *mris, MRI *mri) {

  int n;
  float th;
  double xw,yw,zw,xv,yv,zv,val;
  VERTEX *v;

  for (n=0;n<mris->nvertices;n++) {
    v = &mris->vertices[n] ;
    v->val=0;
  }

  MRIScomputeMetricProperties(mris);
  for (n=0;n<mris->nvertices;n++) {
    v = &mris->vertices[n] ;
    th=v->curv/2.0f;
    xw = v->x + th*v->nx;
    yw = v->y + th*v->ny;
    zw = v->z + th*v->nz;
    MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv) ;
    MRIsampleVolume(mri, xv, yv, zv, &val) ;
    v->val=val;
  }

  for (n=0;n<mris->nvertices;n++) {
    v = &mris->vertices[n] ;
    v->curv = v->val;
  }
}


/* see definition in mrisurf.h */
static int findSurfaceReference(int label) {

  switch (label) {
  case 18:
  case 54:
    return 7;
  case 17:
  case 53:
    return 8;
  case 13:
  case 52:
    return 9;
  case 12:
  case 51:
    return 10;
  case 11:
  case 50:
    return 11;
  case 4:
  case 43:
    return 12;
  case 5:
  case 44:
    return 13;
  default:
    return 0;
  }
}

int main(int argc, char *argv[]) {
  char *subject_fname,*subjects_fname[STRLEN],fname[STRLEN],*cp,*hemi;
  int  nargs,n , m,surface_reference,nsubjects;
  MRI_SURFACE  *mris;
  MRI *mri,*mri_distance, *mri_orig;

  int msec, minutes, seconds ;
  Timer start;

  nargs = handleVersionOption(argc, argv, "mris_distance_to_label");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (!strlen(subjects_dir)) /* hasn't been set on command line */
  {
    cp = getenv("SUBJECTS_DIR") ;
    if (!cp)
      ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
                Progname);
    strcpy(subjects_dir, cp) ;
  }

  if (argc < 3)
    usage_exit() ;

  /* hemisphere information */
  hemi = argv[1];
  for (nsubjects=0 , n = 2 ; n < argc ; n++)
    subjects_fname[nsubjects++]=argv[n];

  if (nlabels==0) {
    fprintf(stderr,"using default option\n");
    fprintf(stderr,"computing distance maps for :\n");
    fprintf(stderr,"      amygdala\n");
    fprintf(stderr,"      hippocampus\n");
    fprintf(stderr,"      pallidum\n");
    fprintf(stderr,"      putamen\n");
    fprintf(stderr,"      caudate\n");
    fprintf(stderr,"      lateral ventricle\n");
    //  fprintf(stderr,"      inferior lateral ventricle\n");
    fprintf(stderr,"      layer IV gray\n");
    nlabels=8;
    if (!stricmp(hemi,(char*)"rh")) { /* right hemisphere */
      labels[0]=54;
      labels[1]=53;
      labels[2]=52;
      labels[3]=51;
      labels[4]=50;
      labels[5]=43;
      labels[6]=44;
      labels[7]=-1;
    } else {
      labels[0]=18;
      labels[1]=17;
      labels[2]=13;
      labels[3]=12;
      labels[4]=11;
      labels[5]=4;
      labels[6]=5;
      labels[7]=-1;
    }
  }

  for ( m = 0 ; m < nsubjects ; m++) {
    subject_fname=subjects_fname[m];

    fprintf(stderr,"\n\nPROCESSING SUBJECT '%s' \n",subject_fname);

    sprintf(fname,"%s/%s/surf/%s.white", subjects_dir,subject_fname,hemi);
    fprintf(stderr, "reading surface from %s...\n", fname) ;
    mris=MRISread(fname);

    if (aseg_fname)
      sprintf(fname,"%s/%s/mri/%s", subjects_dir,subject_fname,aseg_fname);
    else
      sprintf(fname,"%s/%s/mri/aseg.mgz", subjects_dir,subject_fname);

    fprintf(stderr, "reading mri segmentation from %s...\n", fname) ;
    mri=MRIread(fname);

    fprintf(stderr, "allocating distance map\n") ;
    mri_distance=MRIalloc(mri->width,mri->height,mri->depth,MRI_FLOAT);

    for (n=0 ; n < nlabels ; n++) {

      if (labels[n]>=0) {
        fprintf(stderr, "generating distance map for label %d\n", labels[n]) ;
        MRIextractDistanceMap(mri,mri_distance,labels[n],fdistance,mode,NULL);

        fprintf(stderr,
                "extracting distance values for label %d\n", labels[n]) ;
        mrisExtractMRIvalues(mris,mri,mri_distance,fdistance,mode);

        mrisProcessDistanceValues(mris);

        surface_reference=findSurfaceReference(labels[n]);
        if (surface_reference>=3 and surface_reference<=14)
          sprintf(fname,"%s/%s/surf/%s.%s",
                  subjects_dir,subject_fname,hemi,
                  FRAME_FIELD_NAMES[surface_reference]);
        else
          sprintf(fname,"%s/%s/surf/%s.dist_%d",
                  subjects_dir,subject_fname,hemi,labels[n]);

        fprintf(stderr,
                "writing out surface distance file for label %d in %s...\n",
                labels[n],fname) ;
        MRISaverageCurvatures(mris,navgs);
        MRISwriteCurvature(mris,fname);
      } else { /* extract layer IV */
        sprintf(fname,"%s/%s/surf/%s.thickness",
                subjects_dir,subject_fname,hemi);
        fprintf(stderr, "reading curvature from %s...\n", fname) ;
        MRISreadCurvature(mris,fname);

        sprintf(fname,"%s/%s/mri/T1.mgz", subjects_dir,subject_fname);
        fprintf(stderr, "reading orig mri segmentation from %s...\n", fname) ;
        mri_orig=MRIread(fname);
        mrisExtractMidGrayValues(mris,mri_orig);
        MRIfree(&mri_orig);

        surface_reference=3;
        sprintf(fname,"%s/%s/surf/%s.%s",
                subjects_dir,subject_fname,hemi,
                FRAME_FIELD_NAMES[surface_reference]);
        fprintf(stderr,
                "writing out surface distance file for label %d in %s...\n",
                labels[n],fname) ;
        MRISaverageCurvatures(mris,navgs);
        MRISwriteCurvature(mris,fname);
      }
    }

    MRIfree(&mri_distance);
    MRIfree(&mri);
    MRISfree(&mris);
  }

  msec = start.milliseconds() ;
  seconds = (int)((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("mris_distance_to_label took %d minutes and %d seconds.\n",
         minutes, seconds) ;
  exit(0) ;
  return(0) ;  /* for ansi */
}


/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[]) {
  int    nargs=0;
  char   *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, (char*)"SDIR")) {
    strcpy(subjects_dir, argv[2]) ;
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else if (!stricmp(option, (char*)"aseg")) {
    aseg_fname=argv[2];
    fprintf(stderr,"using %s for the aseg volume\n",aseg_fname);
    nargs = 1 ;
  } else if (!stricmp(option,(char*) "mode")) {
    mode=atoi(argv[2]);
    fprintf(stderr,"mode %d : (1 == outside ; 2 == inside ; 3 == both) \n",
            mode);
    nargs = 1 ;
  } else if (!stricmp(option, (char*)"distance")) {
    fdistance=atof(argv[2]);
    fprintf(stderr,"computing distance map for distances smaller than %f\n",
            fdistance);
    nargs = 1 ;
    printf("using %s as subjects directory\n", subjects_dir) ;
  } else if (!stricmp(option,(char*) "-help"))
    print_help() ;
  else if (!stricmp(option, (char*)"-version"))
    print_version() ;
  else if (!stricmp(option,(char*) "navgs")) {
    navgs=atoi(argv[2]);
    fprintf(stderr,"smoothing curv for %d iterations\n",navgs);
    nargs=1;
  } else switch (toupper(*option)) {
  case 'L':
    labels[nlabels++]=atoi(argv[2]);
    fprintf(stderr,"computing distance map for label %d\n",labels[nlabels-1]);
    nargs=1;
    break;
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
usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <hemisphere> <subjects_1>\n",Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

