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
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "label.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static char *log_fname = NULL ;
static void usage_exit(int code) ;

static int  segmentation_flag = -1 ;
static char *annot_prefix = NULL ;
static int quiet = 0 ;
static int scaleup_flag = 0 ;
static int cras =0; // 0 is false.  1 is true
static int cras_not_set = 1 ;
static char *surface_dir = NULL ;
static char *hemi ;
static int erode = 0 ;
static int coords = 0 ;

int
main(int argc, char *argv[]) {
  char   **av, *label_name, *vol_name, *out_name ;
  int    ac, nargs ;
  int    msec, minutes, seconds,  i  ;
  LABEL  *area ;
  Timer start ;
  MRI    *mri,  *mri_seg ;
  double xw, yw, zw, xv, yv, zv, val;

  nargs = handleVersionOption(argc, argv, "mri_label_vals");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  start.reset() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
    usage_exit(1) ;

  vol_name = argv[1] ;
  label_name = argv[2]  ;
  out_name = argv[3] ;

  mri = MRIread(vol_name) ;
  if (!mri)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname, vol_name) ;
  if (scaleup_flag) {
    float scale, fov_x, fov_y, fov_z  ;

    scale = 1.0/MIN(MIN(mri->xsize, mri->ysize),mri->zsize) ;
    fprintf(stderr, "scaling voxel sizes up by %2.2f\n", scale) ;
    mri->xsize *= scale ;
    mri->ysize *= scale ;
    mri->zsize *= scale ;
    fov_x = mri->xsize * mri->width;
    fov_y = mri->ysize * mri->height;
    fov_z = mri->zsize * mri->depth;
    mri->xend = fov_x / 2.0;
    mri->xstart = -mri->xend;
    mri->yend = fov_y / 2.0;
    mri->ystart = -mri->yend;
    mri->zend = fov_z / 2.0;
    mri->zstart = -mri->zend;

    mri->fov = (fov_x > fov_y ? (fov_x > fov_z ? fov_x : fov_z) : (fov_y > fov_z ? fov_y : fov_z) );
  }

  if (segmentation_flag >= 0) {
    int x, y, z  ;
    VECTOR *v_seg, *v_mri ;
    MATRIX *m_seg_to_mri ;

    v_seg = VectorAlloc(4, MATRIX_REAL) ;
    v_mri = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v_seg, 4) = 1.0 ;
    VECTOR_ELT(v_mri, 4) = 1.0 ;

    mri_seg = MRIread(argv[2]) ;
    if (!mri_seg)
      ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s",Progname, argv[2]) ;
    if (erode) {
      MRI *mri_tmp ;

      mri_tmp = MRIclone(mri_seg, NULL) ;
      MRIcopyLabel(mri_seg, mri_tmp, segmentation_flag) ;
      while (erode-- > 0)
        MRIerode(mri_tmp, mri_tmp) ;
      MRIcopy(mri_tmp, mri_seg) ;
      MRIfree(&mri_tmp) ;
    }

    m_seg_to_mri = MRIgetVoxelToVoxelXform(mri_seg, mri) ;
    for (x = 0  ; x  < mri_seg->width ; x++) {
      V3_X(v_seg) = x ;
      for (y = 0  ; y  < mri_seg->height ; y++) {
        V3_Y(v_seg) = y ;
        for (z = 0  ; z  < mri_seg->depth ; z++) {
          V3_Z(v_seg) = z ;
          if (MRIvox(mri_seg, x, y,  z) == segmentation_flag) {
            MatrixMultiply(m_seg_to_mri, v_seg, v_mri) ;
            xv = V3_X(v_mri) ;
            yv = V3_Y(v_mri) ;
            zv = V3_Z(v_mri) ;
            MRIsampleVolumeType(mri, xv,  yv, zv, &val, SAMPLE_NEAREST);
#if  0
            if (val < .000001) {
              val *= 1000000;
              printf("%f*0.000001\n", val);
            } else
#endif
              if (coords)
                printf("%2.1f %2.1f %2.1f %f\n", xv, yv, zv, val);
              else
                printf("%f\n", val);
          }
        }
      }
    }
    MatrixFree(&m_seg_to_mri) ;
    VectorFree(&v_seg) ;
    VectorFree(&v_mri) ;
  } else {
    if (cras == 1)
      fprintf(stderr,"using the label coordinates to be c_(r,a,s) != 0.\n");

    if (surface_dir) {
      MRI_SURFACE *mris ;
      char fname[STRLEN] ;
      sprintf(fname, "%s/%s.white", surface_dir, hemi) ;
      mris = MRISread(fname) ;
      if (mris == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s...\n", Progname,fname) ;
      sprintf(fname, "%s/%s.thickness", surface_dir, hemi) ;
      if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
        ErrorExit(ERROR_BADPARM, "%s: could not read thickness file %s...\n", Progname,fname) ;

      if (annot_prefix)   /* read an annotation in and print vals in it */
      {
#define MAX_ANNOT 10000
        int   vno, annot_counts[MAX_ANNOT], index ;
        VERTEX *v ;
        double xw, yw, zw, xv, yv, zv, val ;
        float annot_means[MAX_ANNOT] ;
        FILE  *fp ;

        memset(annot_means, 0, sizeof(annot_means)) ;
        memset(annot_counts, 0, sizeof(annot_counts)) ;
        if (MRISreadAnnotation(mris, label_name) != NO_ERROR)
          ErrorExit(ERROR_BADPARM, "%s: could not read annotation file %s...\n", Progname,fname) ;
        if (mris->ct == NULL)
          ErrorExit(ERROR_BADPARM, "%s: annot file does not contain a color table, specifiy one with -t ", Progname);
        for (vno = 0 ; vno < mris->nvertices ; vno++) {
          v = &mris->vertices[vno] ;
          if (v->ripflag)
            continue ;
          CTABfindAnnotation(mris->ct, v->annotation, &index) ;
          if (index >= 0 && index < mris->ct->nentries) {
            annot_counts[index]++ ;
            xw = v->x + v->curv*.5*v->nx ;
            yw = v->y + v->curv*.5*v->ny ;
            zw = v->z + v->curv*.5*v->nz ;
            if (cras == 1)
              MRIworldToVoxel(mri, xw, yw,  zw, &xv, &yv, &zv) ;
            else
              MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv);
            MRIsampleVolume(mri, xv, yv, zv, &val) ;
            annot_means[index] += val ;
            int req = snprintf(fname, STRLEN, "%s-%s-%s.dat",
			       annot_prefix, hemi, mris->ct->entries[index]->name) ;
	    if( req >= STRLEN ) {
	      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
	    }
            fp = fopen(fname, "a") ;
            fprintf(fp, "%f\n", val) ;
            fclose(fp) ;
          }
        }
      }
      else  /* read label in and print vals in it */
      {
        area = LabelRead(NULL, label_name) ;
        if (!area)
          ErrorExit(ERROR_NOFILE, "%s: could not read label from %s",Progname, label_name) ;

      }
    } else {
      area = LabelRead(NULL, label_name) ;
      if (!area)
        ErrorExit(ERROR_NOFILE, "%s: could not read label from %s",Progname, label_name) ;

      for (i = 0 ;  i  < area->n_points ; i++) {

        xw =  area->lv[i].x ;
        yw =  area->lv[i].y ;
        zw =  area->lv[i].z ;
        if (cras == 1 || (cras_not_set && area->coords == LABEL_COORDS_SCANNER_RAS))
          MRIworldToVoxel(mri, xw, yw,  zw, &xv, &yv, &zv) ;
        else
          MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv);
        MRIsampleVolumeType(mri, xv,  yv, zv, &val, SAMPLE_NEAREST);
#if 0
        if (val < .000001) {
          val *= 1000000;
          printf("%f*0.000001\n", val);
        } else
#endif
          printf("%f\n", val);
      }
    }
  }
  msec = start.milliseconds() ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;

  if (DIAG_VERBOSE_ON)
    fprintf(stderr, "label value extractiong took %d minutes and %d seconds.\n",
            minutes, seconds) ;

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
  if (strcmp("cras", option) == 0)
  {
    cras = 1;
    cras_not_set = 0 ;
  }
  else if (strcmp("scaleup", option) == 0)
    scaleup_flag = 1 ;
  else if (strcmp("segmentation", option) == 0) {
    segmentation_flag = atoi(argv[2]) ;
    nargs =  1  ;
    fprintf(stderr,"using segmentation %d (%s) as label...\n", segmentation_flag,
            cma_label_to_name(segmentation_flag))  ;
  } else if (strcmp("coords", option) == 0) {
    coords = 1 ;
    fprintf(stderr,"printing out coordinates with values\n") ;
  } else if (strcmp("erode", option) == 0) {
    erode = atoi(argv[2]) ;
    nargs = 1;
    fprintf(stderr,"eroding segmentation label %d times (must specify label with -segmentaiton)\n", erode);
  } else {
    switch (toupper(*option)) {
    case 'S':
        surface_dir = argv[2] ;
      hemi = argv[3] ;
      printf("sampling from midpoint of cortical ribbon in %s (hemi=%s)\n", surface_dir, hemi) ;
      nargs = 2 ;
      break ;
    case 'Q':
      quiet = 1 ;
      break ;
    case 'A':
      annot_prefix = argv[2] ;
      nargs =  1  ;
      fprintf(stderr,"reading annotation file, and outputting with prefix %s ...\n", annot_prefix)  ;
      break ;
    case 'L':
      log_fname = argv[2] ;
      nargs = 1 ;
      fprintf(stderr, "logging results to %s\n", log_fname) ;
      break ;
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }
  }
  return(nargs) ;
}
/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static void
usage_exit(int code) {
  printf("usage: %s [options] <volume> <label file>\n",  Progname) ;
  printf("where options are\n");
  printf("   -cras   label created in the coordinates where c_(r,a,s) != 0\n");
  printf("             if it did not work, try using this option.\n");
  // printf("   -q      quiet\n");
  // printf("   -a      all\n");
  // printf("   -l file log results to a file.\n");
  printf("   -u      print this help.\n");
  exit(code) ;
}
