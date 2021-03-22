/*
 * Original Author: Florent Segonne 
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
#include "cma.h"


#include "fastmarching.h"

const char *Progname ;

static float remove_csf_from_paths(MRI *mri_distance, MRI *mri_area, 
                                   MRI *mri_csf) ;
static int get_option(int argc, char *argv[]) ;
static int MRIexpandDistancesIntoGrayMatter(MRI *mri_distance, MRI *mri_aseg, MRI *mri_white, int target_label) ;

static MRI *mri_area = NULL ;
static MRI *mri_white = NULL, *mri_aseg = NULL ;
static char *surf_name = NULL ;
static char *csf_name = NULL ;
static float normalize = -1.0 ;
static float wthresh = -1 ;
static MRI *mri_thresh = NULL ;
static float anterior_dist = -1 ;
static float posterior_dist = -1 ;

static float binarize = 0.0 ;
static int percent = 0;

static int ndilations = 0 ;
MRI *MRIthresholdPosterior(MRI *mri_src, MRI *mri_dst, float posterior_dist) ;
MRI *MRIthresholdAnterior(MRI *mri_src, MRI *mri_dst, float anterior_dist) ;
MRI *MRIscaleDistanceTransformToPercentMax(MRI *mri_in, MRI *mri_out);

int main(int argc, char *argv[]) {
  MRI         *mri,*mri_distance;
  int         label,mode;
  float       max_distance;
  int         nargs ;
  MRI_SURFACE *mris = NULL ;
  MRI         *mri_csf = NULL ;

  max_distance=10;
  mode=1;

  Progname=argv[0];

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  fprintf(stderr,"mri_distance_transform <input volume> <label> <max_distance> <mode[=1]> <output volume>\n");
  fprintf(stderr,"mode : 1 = outside , mode : 2 = inside , mode : 3 = both, mode : 4 = both unsigned \n");

  if (argc < 5)
    exit(0) ;

  mri=MRIread(argv[1]);
  if (mri == NULL)
    ErrorExit(ERROR_NOFILE, "%s: could not read volume from %s", Progname, argv[1]) ;
  label=atoi(argv[2]);
  max_distance=atof(argv[3]);
  mode=atoi(argv[4]);

  if (binarize > 0)
    MRIbinarize(mri, mri, binarize, 0, label) ;
  if (ndilations > 0)
      MRIdilateLabel(mri, mri, label, ndilations) ;
  if (posterior_dist > 0)
    MRIthresholdPosterior(mri, mri, posterior_dist) ;
  if (anterior_dist > 0)
    MRIthresholdAnterior(mri, mri, anterior_dist) ;
  if (mri->type != MRI_FLOAT)
    {
      MRI *mri_tmp;
      mri_tmp = MRIchangeType(mri, MRI_FLOAT, 0, 1, 1) ;
      MRIfree(&mri) ;
      mri = mri_tmp ;
    }
  fprintf(stderr,"label=%d distance=%.0f mode=%d\n",label,max_distance,mode);

  mri_distance=MRIalloc(mri->width,mri->height,mri->depth,MRI_FLOAT);
  MRIcopyHeader(mri, mri_distance) ;

  if (surf_name)
    {
      int         x, y, z ;

      mris = MRISread(surf_name) ;
      if (mris == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read surface from %s", Progname, surf_name) ;
      mri_white = MRIclone(mri, NULL) ;
      MRISfillInterior(mris, mri_white->xsize, mri_white) ;

      if (mri_aseg)  // turn off non-wm/cortex labels
      {
        int label ;
        for (x = 0 ; x < mri_aseg->width ; x++)
          for (y = 0 ; y < mri_aseg->height ; y++)
            for (z = 0 ; z < mri_aseg->depth ; z++)
            {
              if (nint(MRIgetVoxVal(mri_white, x, y, z, 0)) == 0)
                continue ;
              label = nint(MRIgetVoxVal(mri_aseg, x, y, z, 0)) ;
              if (!IS_CORTEX(label) && !IS_WHITE_CLASS(label))
                MRIsetVoxVal(mri_white, x, y, z, 0, 0) ; // turn it off
            }
      }

      if (normalize > 0)
      {
#if 0
        normalize = MRItotalVoxelsOn(mri_white, 1) / (mri_white->xsize*mri_white->ysize*mri_white->zsize); ;
        printf("normalizing distances by pow(%2.0f, 1/3) = %2.1f\n", normalize, pow(normalize, 1/3.0)) ;
        normalize = pow(normalize, 1.0/3.0) ;
#endif

        normalize = sqrt(mris->total_area) ;
        printf("normalizing surface distances by sqrt(%2.1f) = %2.1f\n", 
               mris->total_area,normalize) ;
        
      }
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_white, "w.mgz") ;
      if (mri_thresh != NULL)  // remove voxels that are in the ventricles
      {
        int  csf_vox = 0 ;

        mri_csf = MRIclone(mri_thresh, NULL) ;
        for (x = 0 ; x < mri_thresh->width ; x++)
          for (y = 0 ; y < mri_thresh->height ; y++)
            for (z = 0 ; z < mri_thresh->depth ; z++)
              {
                if (MRIgetVoxVal(mri_thresh, x, y, z, 0) < wthresh)
                  {
                    if ((int)MRIgetVoxVal(mri_white, x, y, z, 0) > 0)
                    {
                      MRIsetVoxVal(mri_csf, x, y, z, 0, 1) ;
                      MRIsetVoxVal(mri_white, x, y, z, 0, 0) ;
                      csf_vox++ ;
                    }
                  }
              }
        if (csf_name)
          {
            FILE *fp ;

            printf("writing csf vox (%d) to %s\n", csf_vox, csf_name) ;
            fp = fopen(csf_name, "w") ;
            fprintf(fp, "%2.1f\n", (float)csf_vox) ;
            fclose(fp) ;
          }
      }

      mri_aseg = MRIdilate(mri_white, NULL) ;
      for (x = 0 ; x < mri_aseg->width ; x++)
        for (y = 0 ; y < mri_aseg->height ; y++)
          for (z = 0 ; z < mri_aseg->depth ; z++)
          {

            if (MRIgetVoxVal(mri_white, x, y, z, 0) != 0)
              MRIsetVoxVal(mri_aseg, x, y, z, 0, Left_Cerebral_White_Matter) ;
            else if ((MRIgetVoxVal(mri_aseg, x, y, z, 0) != 0)) 
              MRIsetVoxVal(mri_aseg, x, y, z, 0, Left_Cerebral_Cortex) ;
          }
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_aseg, "a.mgz") ;
    }

  mri_distance=MRIextractDistanceMap(mri,mri_distance,label, max_distance, mode, mri_white);

  if (mri_aseg)
    {
      MRIexpandDistancesIntoGrayMatter(mri_distance, mri_aseg, mri_white, label) ;
    }

  if (normalize > 0)
    MRIscalarMul(mri_distance, mri_distance, 1.0/normalize) ;
  if (percent)
    MRIscaleDistanceTransformToPercentMax(mri_distance,mri_distance) ;
  if (mri_area)
    {
      float csf_vox ;

      if (!mri_thresh || !mris)
        ErrorExit(ERROR_BADPARM, 
                  "%s: must specify -wsurf, -wthresh and -csf with -label",
                  Progname) ;
      MRIbinarize(mri_area, mri_area, 1, 0, 1) ;
      MRIwrite(mri_csf, "csf.mgz") ;
      csf_vox = remove_csf_from_paths(mri_distance, mri_area, mri_csf) ;
      if (csf_name)
        {
#if 0
          int      x, y, z ;
          double   label_dist, dist ;
#endif
          FILE *fp ;

#if 0
          printf("checking for csf in regions closer than %2.1fmm\n", label_dist);
          
          /*
            count the # of voxels that are csf in the interior of the white
            matter, that are also closer to the source of the distance
            transform than a specified label (to ignore distant csf)
          */
          label_dist = MRImeanInLabel(mri_distance, mri_area, 1) ;
          csf_vox = 0 ;
          for (x = 0 ; x < mri_thresh->width ; x++)
            for (y = 0 ; y < mri_thresh->height ; y++)
              for (z = 0 ; z < mri_thresh->depth ; z++)
                {
                  if ((int)MRIgetVoxVal(mri_orig_white, x, y, z, 0) == 0)
                    continue ;  // not in interior of white surface
                  if (MRIgetVoxVal(mri_thresh, x, y, z, 0) < wthresh)
                    {
                      dist = MRIvoxelMin(mri_distance, x, y, z, 3) ;
                      if (dist < label_dist)
                        csf_vox++ ;
                    }
                }
#endif
          printf("writing csf vox (%2.0f) to %s\n", csf_vox, csf_name) ;
          fp = fopen(csf_name, "w") ;
          fprintf(fp, "%2.1f\n", csf_vox) ;
          fclose(fp) ;
        }
    }
  
  MRIwrite(mri_distance,argv[5]);
  
  if (mri_thresh)
    MRIfree(&mri_thresh) ;
  if (mris)
    MRISfree(&mris) ;
  MRIfree(&mri);
  MRIfree(&mri_distance);

  return 0;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static int wm_labels[] =
  {
    Left_Cerebral_White_Matter,
    Right_Cerebral_White_Matter,
    Left_Cerebellum_White_Matter,
    Right_Cerebellum_White_Matter,
    Left_Thalamus_Proper,
    Right_Thalamus_Proper,
    Brain_Stem,
    Left_VentralDC,
    Right_VentralDC
  } ;
#define NWM_LABELS (int)((sizeof(wm_labels) / sizeof(wm_labels[0])))

static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  /*  StrUpper(option) ;*/
  if (!stricmp(option, "wm"))
    {
      int i ;
      nargs = 1 ;
      mri_aseg = MRIread(argv[2]) ;
      if (mri_aseg == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s", Progname, argv[2]) ;
      mri_white = MRIclone(mri_aseg, NULL);
      for (i = 0 ; i < NWM_LABELS ; i++)
        MRIcopyLabel(mri_aseg, mri_white, wm_labels[i]) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        MRIwrite(mri_white, "white.mgz");
    }
  else if (!stricmp(option, "anterior"))
  {
    anterior_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("using anterior-most %2.3f mm of label only\n", anterior_dist) ;
  }
  else if (!stricmp(option, "aseg"))
  {
    mri_aseg = MRIread(argv[2]) ;
    if (mri_aseg == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read aseg volume from %s", Progname, argv[2]) ;
    nargs = 1 ;
    printf("using aseg volume %s\n", argv[2]) ;
  }
  else if (!stricmp(option, "label"))
  {
    mri_area = MRIread(argv[2]) ;
    if (mri_area == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label volume from %s", 
                Progname,
                argv[2]) ;
    nargs = 1 ;
    printf("computing CSF volume bordering regions closer than label %s\n",
           argv[2]);
  }
  else if (!stricmp(option, "posterior"))
  {
    posterior_dist = atof(argv[2]) ;
    nargs = 1 ;
    printf("using posterior-most %2.3f mm of label only\n", posterior_dist) ;
  }
  else if (!stricmp(option, "wsurf"))
    {
      nargs = 1 ;
      surf_name = argv[2] ;
    }
  else if (!stricmp(option, "csf"))
    {
      nargs = 1 ;
      csf_name = argv[2] ;
    }
  else if (!stricmp(option, "normalize"))
    {
      normalize = 1 ;
      printf("normalizing distances by sqrt(surface area)\n") ;
    }
  else if (!stricmp(option, "wthresh"))
    {
      nargs = 2 ;
      mri_thresh = MRIread(argv[2]) ;
      if (mri_thresh == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read threshold volume %s",Progname, argv[2]) ;
      wthresh = atof(argv[3]) ;
      printf("threshing interior volume %s at %2.0f to remove ventricles from mask\n", argv[2], wthresh) ;
    }
  else if (!stricmp(option, "dilate"))
    {
      ndilations = atoi(argv[2]) ;
      nargs = 1;
      printf("performing %d dilations on labeled volume before computing distance transform\n", 
             ndilations) ;
    }
  else if (!stricmp(option, "b"))
    {
      binarize = atof(argv[2]) ;
      nargs = 1;
      printf("binarizing input data with thresh = %2.1f\n", binarize) ;
    }
  else if (!stricmp(option, "p"))
    {
      percent=1;
      printf("scaling distances to be percent of max\n");
    }
  return(nargs) ;
}

static int
MRIexpandDistancesIntoGrayMatter(MRI *mri_distance, MRI *mri_aseg, MRI *mri_white, int target_label)
{
  int   x, y, z, xi, yi, zi, xk, yk, zk, label ;
  float dist, min_dist, vdist ;

  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (label == target_label)
          MRIsetVoxVal(mri_distance, x, y, z, 0, 0) ;
      }

  // sample from wm into adjacent gm to prevent jumping the banks of a sulcus
  for (x = 0 ; x < mri_aseg->width ; x++)
    for (y = 0 ; y < mri_aseg->height ; y++)
      for (z = 0 ; z < mri_aseg->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0) ;
        if (IS_GM(label) == 0)
          continue ;
        min_dist = MRIgetVoxVal(mri_distance, x, y, z, 0) ;
        for (xk = -1 ; xk <= 1 ; xk++)
        {
          xi = mri_aseg->xi[x+xk] ;
          for (yk = -1 ; yk <= 1 ; yk++)
          {
            yi = mri_aseg->yi[y+yk] ;
            for (zk = -1 ; zk <= 1 ; zk++)
            {
              zi = mri_aseg->zi[z+zk] ;
              /* only allow distances to propagate through the wm.
                 This will prevent the distances from jumping
                 across the banks of a sulcus (which is the whole point)
              */
              if (MRIgetVoxVal(mri_white, xi, yi, zi, 0) == 0)
                continue ;  
              vdist = sqrt(xk*xk + yk*yk + zk*zk);
              zi = mri_aseg->zi[z+zk] ;
              dist = MRIgetVoxVal(mri_distance, xi, yi, zi, 0) + vdist ;
              if (dist < min_dist)
                min_dist = dist ;
            }
          }
        }
        MRIsetVoxVal(mri_distance, x, y, z, 0, min_dist);
      }

  return(NO_ERROR) ;
}

MRI *
MRIthresholdAnterior(MRI *mri_src, MRI *mri_dst, float anterior_dist)
{
  MATRIX  *m_vox2ras ;
  VECTOR  *v_vox, *v_ras ;
  int     x, y, z, f ;
  float  amax ;

  m_vox2ras = MRIgetVoxelToRasXform(mri_src) ;
  v_vox = VectorAlloc(4, MATRIX_REAL) ;
  v_ras = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_vox,4) = VECTOR_ELT(v_ras,4) = 1.0 ;

  if (mri_dst == NULL)
    mri_dst = MRIcopy(mri_src, NULL) ;

  amax = -100000;
  for (f = 0 ; f < mri_dst->nframes ; f++)
    {
      for (x = 0 ; x < mri_dst->width ; x++)
        {
          for (y = 0 ; y < mri_dst->height ; y++)
            {
              for (z = 0 ; z < mri_dst->depth ; z++)
                {
                  if (x == Gx && y == Gy && z == Gz)
                    DiagBreak() ;
                  if ((int)MRIgetVoxVal(mri_dst, x, y, z, f) == 0)
                    continue ;
                  VECTOR_ELT(v_vox, 1) = x ;
                  VECTOR_ELT(v_vox, 2) = y ;
                  VECTOR_ELT(v_vox, 3) = z ;
                  MatrixMultiply(m_vox2ras, v_vox, v_ras) ;
                  if (VECTOR_ELT(v_ras, 2) > amax)
                    amax = VECTOR_ELT(v_ras, 2) ;
                }
            }
        }
    }
  amax -= anterior_dist ;
  for (f = 0 ; f < mri_dst->nframes ; f++)
    {
      for (x = 0 ; x < mri_dst->width ; x++)
        {
          for (y = 0 ; y < mri_dst->height ; y++)
            {
              for (z = 0 ; z < mri_dst->depth ; z++)
                {
                  if (x == Gx && y == Gy && z == Gz)
                    DiagBreak() ;
                  if ((int)MRIgetVoxVal(mri_dst, x, y, z, f) == 0)
                    continue ;
                  VECTOR_ELT(v_vox, 1) = x ;
                  VECTOR_ELT(v_vox, 2) = y ;
                  VECTOR_ELT(v_vox, 3) = z ;
                  MatrixMultiply(m_vox2ras, v_vox, v_ras) ;
                  if (VECTOR_ELT(v_ras, 2) < amax)
                    MRIsetVoxVal(mri_dst, x, y, z, f, 0) ;
                }
            }
        }
    }

  MatrixFree(&m_vox2ras) ; VectorFree(&v_vox) ; VectorFree(&v_ras) ;
  return(mri_dst) ;
}

MRI *
MRIthresholdPosterior(MRI *mri_src, MRI *mri_dst, float posterior_dist)
{
  MATRIX  *m_vox2ras ;
  VECTOR  *v_vox, *v_ras ;
  int     x, y, z, f ;
  float  amin ;

  m_vox2ras = MRIgetVoxelToRasXform(mri_src) ;
  v_vox = VectorAlloc(4, MATRIX_REAL) ;
  v_ras = VectorAlloc(4, MATRIX_REAL) ;
  VECTOR_ELT(v_vox,4) = VECTOR_ELT(v_ras,4) = 1.0 ;

  if (mri_dst == NULL)
    mri_dst = MRIcopy(mri_src, NULL) ;

  amin = 100000;
  for (f = 0 ; f < mri_dst->nframes ; f++)
    {
      for (x = 0 ; x < mri_dst->width ; x++)
        {
          for (y = 0 ; y < mri_dst->height ; y++)
            {
              for (z = 0 ; z < mri_dst->depth ; z++)
                {
                  if (x == Gx && y == Gy && z == Gz)
                    DiagBreak() ;
                  if ((int)MRIgetVoxVal(mri_dst, x, y, z, f) == 0)
                    continue ;
                  VECTOR_ELT(v_vox, 1) = x ;
                  VECTOR_ELT(v_vox, 2) = y ;
                  VECTOR_ELT(v_vox, 3) = z ;
                  MatrixMultiply(m_vox2ras, v_vox, v_ras) ;
                  if (VECTOR_ELT(v_ras, 2) < amin)
                    amin = VECTOR_ELT(v_ras, 2) ;
                }
            }
        }
    }
  amin += posterior_dist ;
  for (f = 0 ; f < mri_dst->nframes ; f++)
    {
      for (x = 0 ; x < mri_dst->width ; x++)
        {
          for (y = 0 ; y < mri_dst->height ; y++)
            {
              for (z = 0 ; z < mri_dst->depth ; z++)
                {
                  if (x == Gx && y == Gy && z == Gz)
                    DiagBreak() ;
                  if ((int)MRIgetVoxVal(mri_dst, x, y, z, f) == 0)
                    continue ;
                  VECTOR_ELT(v_vox, 1) = x ;
                  VECTOR_ELT(v_vox, 2) = y ;
                  VECTOR_ELT(v_vox, 3) = z ;
                  MatrixMultiply(m_vox2ras, v_vox, v_ras) ;
                  if (VECTOR_ELT(v_ras, 2) > amin)
                    MRIsetVoxVal(mri_dst, x, y, z, f, 0) ;
                }
            }
        }
    }

  MatrixFree(&m_vox2ras) ; VectorFree(&v_vox) ; VectorFree(&v_ras) ;
  return(mri_dst) ;
}

static float
remove_csf_from_paths(MRI *mri_distance, MRI *mri_area, MRI *mri_csf)
{
  int    x, y, z, steps, xk, yk, zk, xi, yi, zi, csf_nbr,total_label_vox;
  float  xf, yf, zf, mean_csf_len, csf_len, total_csf_len ;
  double dx, dy, dz, distance ;
  MRI   *mri_orig_csf ;

  mri_orig_csf = MRIcopy(mri_csf, NULL) ;

  total_label_vox = 0;
  total_csf_len = 0.0 ;
  for (x = 0 ; x < mri_area->width ; x++)
    for (y = 0 ; y < mri_area->height ; y++)
      for (z = 0 ; z < mri_area->depth ; z++)
      {
        if ((int)MRIgetVoxVal(mri_area, x, y, z, 0) == 0)
          continue ;
        total_label_vox++ ;
        /*
          now trace path down distance transform gradient until we get to
          the interior (dtrans <= 0) and count the # of adjacent csf
          voxels.
        */
#define STEP_SIZE 0.1
#define MAX_STEPS (1000/STEP_SIZE)
        xf = x ; yf = y ; zf = z ;
        steps = 0 ;
        csf_len = 0.0 ;
        do
        {
          csf_nbr = 0 ;
          for (xk = -1 ; xk <= 1 && !csf_nbr ; xk++)
          {
            xi = mri_csf->xi[nint(xf+xk)] ;
            for (yk = -1 ; yk <= 1 && !csf_nbr ; yk++)
            {
              yi = mri_csf->yi[nint(yf+yk)] ;
              for (zk = -1 ; zk <= 1 && !csf_nbr ; zk++)
              {
                zi = mri_csf->zi[nint(zf+zk)] ;
                if (fabs(xk)+ fabs(yk) + fabs(zk) > 1)
                  continue ; // only 6-connected nbrs
                if ((int)MRIgetVoxVal(mri_csf, xi, yi, zi, 0) != 0)
                {
                  // don't count it twice
                  //                  MRIsetVoxVal(mri_csf, xi, yi, zi, 0, 0) ; 
                  csf_nbr = 1 ;
                }
              }
            }
          }
          if (csf_nbr)
            csf_len += STEP_SIZE  ;

          // take next step
          MRIsampleVolumeGradient(mri_distance, xf, yf, zf,
                                  &dx, &dy, &dz) ;
          xf -= STEP_SIZE *dx ; yf -= STEP_SIZE *dy ; zf -= STEP_SIZE *dz ;
          MRIsampleVolume(mri_distance, xf, yf, zf, &distance) ;
        } while (distance > 0 && steps++ < MAX_STEPS) ;
        if (steps >= MAX_STEPS) // couldn't get to label??
        {
          total_label_vox-- ;
          continue ;
        }
        distance = MRIgetVoxVal(mri_distance, x, y, z, 0) ;
        if (steps*STEP_SIZE > distance*1.1)
        {
          total_label_vox-- ;
          DiagBreak();
          continue ; // too suboptimal a path
        }
        if (csf_len > 100)
          DiagBreak() ;
        // correct for suboptimal paths
        total_csf_len += (csf_len * distance / (steps*STEP_SIZE));
        //        MRIsetVoxVal(mri_distance, x, y, z, 0, distance-csf_len) ;
        /*        MRIcopy(mri_orig_csf, mri_csf) ; */ // restore csf labels
      }
  MRIfree(&mri_orig_csf) ;
  if (total_label_vox > 0)
    mean_csf_len = total_csf_len / total_label_vox ;
  else
    mean_csf_len = 0.0 ;
  return(mean_csf_len) ;
}

MRI *
MRIscaleDistanceTransformToPercentMax(MRI *mri_in, MRI *mri_out)
{
  float min_val, max_val, val ;
  int   x, y, z ;

  MRIvalRange(mri_in, &min_val, &max_val) ;
  printf("scaling distance transforms by [%2.1f, %2.1f]\n", min_val, max_val);
  mri_out = MRIcopy(mri_in, mri_out) ;

  for (x = 0 ; x < mri_out->width ; x++)
    for (y = 0 ; y < mri_out->height ; y++)
      for (z = 0 ; z < mri_out->depth ; z++)
      {
	val = MRIgetVoxVal(mri_out, x, y, z, 0) ;
	if (val > 0)
	  val /= max_val ;
	else if (min_val < 0)
	  val /= -min_val ;
	MRIsetVoxVal(mri_out, x, y,z, 0, val) ;
      }

  return(mri_out) ;
}

