/**
 * @brief use tp1's ctrl volume to normalize tp2
 *
 * Used in longitudinal processing.  Rather than calculating the control
 * points to use to normalize the volume, just use the control points from the
 * prior scan (timepoint 1, aka tp1).
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
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "mrinorm.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "fio.h"
#include "version.h"
#include "transform.h"

static int debug_flag = 0;

const char *Progname;

static char *tp1_T1_fname = NULL;
static char *tp1_ctrl_fname  = NULL;

/* The following specifies the src and dst volumes
   of the input FSL/LTA transform */
static MRI *lta_src = 0;
static MRI *lta_dst = 0;
static int invert = 0 ;
static char *xform_fname = NULL;
static float bias_sigma = 8.0 ;
static char *mask1_fname = NULL;
static char *mask2_fname = NULL;

static float noise_threshold = 1.0;

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;
static void usage(int exit_val);

int main(int argc, char *argv[]) {
  char **av;
  MRI *mri_T1, *mri_tmp, *mri_ctrl, *mri_in, *mri_out;
  MRI *mri_snr, *mri_bias;
  MRI *mri_mask1 = NULL;
  MRI *mri_mask2 = NULL;

  int ac, nargs;
  int  width, height, depth, x, y, z;
  int mask1_set = 0;
  int mask2_set = 0;
  int i, j, k, cx, cy, cz, count;
  LTA          *lta = 0;
  int          transform_type;
  double mean, std, value, src, bias, norm;
//  HISTOGRAM *h;
//  float bin_size;
//  int nbins, bin_no;
  double mean1, std1, mean2, std2, count1, count2, slope, offset;
  VOL_GEOM vgtmp;
  LT *lt = NULL;
  MATRIX *m_tmp = NULL;

  Progname = argv[0];

  nargs = handleVersionOption(argc, argv, "mri_normalize_tp2");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc  !=  3)
    usage(1);

  if (tp1_ctrl_fname == NULL  || tp1_T1_fname == NULL) {
    printf("Use options to specify ctrl volume and T1 volume for tp1\n");
    usage(1);
  }

  mri_in = MRIread(argv[1]) ;
  if (!mri_in)
    ErrorExit(ERROR_BADPARM, "%s: could not read input volume %s",
              Progname, argv[1]) ;

  mri_T1 = MRIread(tp1_T1_fname) ;
  if (!mri_T1)
    ErrorExit(ERROR_BADPARM, "%s: could not read T1 volume for tp1 %s",
              Progname, tp1_T1_fname) ;

  mri_ctrl = MRIread(tp1_ctrl_fname) ;
  if (!mri_ctrl)
    ErrorExit(ERROR_BADPARM,
              "%s: could not read control points volume for tp1 %s",
              Progname, tp1_ctrl_fname) ;

  if ((mri_in->width != mri_T1->width) ||
      (mri_in->height != mri_T1->height) ||
      (mri_in->depth != mri_T1->depth) ||
      (mri_in->width != mri_ctrl->width) ||
      (mri_in->height != mri_ctrl->height) ||
      (mri_in->depth != mri_ctrl->depth)
     ) ErrorExit
         (ERROR_BADPARM,
          "%s: three input volumes have different sizes \n", Progname);

  if (mask1_fname) {
    mri_mask1 = MRIread(mask1_fname) ;
    if (!mri_mask1)
      ErrorExit(ERROR_BADPARM,
                "%s, could not read mask volume for tp1 %s",
                Progname, mask1_fname);
    mask1_set = 1;
    if ((mri_mask1->width != mri_in->width) ||
        (mri_mask1->height != mri_in->height) ||
        (mri_mask1->depth != mri_in->depth))
      ErrorExit
      (ERROR_BADPARM,
       "%s:  mask volumes have different sizes than other volumes \n",
       Progname);
  }
  if (mask2_fname) {
    mri_mask2 = MRIread(mask2_fname) ;
    if (!mri_mask2)
      ErrorExit
      (ERROR_BADPARM,
       "%s, could not read mask volume for tp2 %s",
       Progname, mask2_fname);
    mask2_set = 1;
    if ((mri_mask2->width != mri_T1->width) ||
        (mri_mask2->height != mri_T1->height) ||
        (mri_mask2->depth != mri_T1->depth)
       )
      ErrorExit
      (ERROR_BADPARM,
       "%s:  mask volumes have different sizes than other volumes \n",
       Progname);
  }

  width = mri_in->width ;
  height = mri_in->height ;
  depth = mri_in->depth ;

  //nbins = 200;
  //h = HISTOalloc(nbins);

  mri_out = MRIclone(mri_in, NULL) ;

  /* Read LTA transform and apply it to mri_ctrl */
  if (xform_fname != NULL) {
    // read transform
    transform_type =  TransformFileNameType(xform_fname);

    if (transform_type == MNI_TRANSFORM_TYPE ||
        transform_type == TRANSFORM_ARRAY_TYPE ||
        transform_type == REGISTER_DAT ||
        transform_type == FSLREG_TYPE
       ) {
      printf("Reading transform ...\n");
      lta = LTAreadEx(xform_fname) ;
      if (!lta)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
                  Progname, xform_fname) ;

      if (transform_type == FSLREG_TYPE) {
        if (lta_src == 0 || lta_dst == 0) {
          fprintf
          (stderr,
           "ERROR: fslmat does not have information "
           "on the src and dst volumes\n");
          fprintf
          (stderr,
           "ERROR: you must give options '-lta_src' and "
           "'-lta_dst' to specify the src and dst volume infos\n");
        }

        LTAmodifySrcDstGeom
        (lta, lta_src, lta_dst); // add src and dst information
        LTAchangeType(lta, LINEAR_VOX_TO_VOX); //this is necessary
      }

      if (lta->xforms[0].src.valid == 0) {
        if (lta_src == 0) {
          fprintf
          (stderr,
           "The transform does not have the valid src volume info.\n");
          fprintf
          (stderr,
           "Either you give src volume info by option -lta_src or\n");
          fprintf(stderr, "make the transform to have the valid src info.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        } else {
          LTAmodifySrcDstGeom(lta, lta_src, NULL); // add src information
          //   getVolGeom(lta_src, &lt->src);
        }
      }
      if (lta->xforms[0].dst.valid == 0) {
        if (lta_dst == 0) {
          fprintf
          (stderr,
           "The transform does not have the valid dst volume info.\n");
          fprintf
          (stderr,
           "Either you give src volume info by option -lta_dst or\n");
          fprintf
          (stderr,
           "make the transform to have the valid dst info.\n");
          fprintf
          (stderr,
           "If the dst was average_305, then you can set\n");
          fprintf
          (stderr,
           "environmental variable USE_AVERAGE305 true\n");
          fprintf
          (stderr,
           "without giving the dst volume for RAS-to-RAS transform.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        } else {
          LTAmodifySrcDstGeom(lta, NULL, lta_dst); // add  dst information
        }
      }
    } else {
      ErrorExit
      (ERROR_BADPARM,
       "transform is not of MNI, nor Register.dat, nor FSLMAT type");
    }

    if (invert) {
      m_tmp = lta->xforms[0].m_L ;
      lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
      MatrixFree(&m_tmp) ;
      lt = &lta->xforms[0];
      if (lt->dst.valid == 0 || lt->src.valid == 0) {
        fprintf
        (stderr,
         "WARNING:***********************************************\n");
        fprintf
        (stderr,
         "WARNING: dst volume infor is invalid.  "
         "Most likely produce wrong inverse.\n");
        fprintf
        (stderr,
         "WARNING:***********************************************\n");
      }
      //copyVolGeom(&lt->dst, &vgtmp);
      vgtmp = lt->dst;
      //copyVolGeom(&lt->src, &lt->dst);
      lt->dst = lt->src;
      //copyVolGeom(&vgtmp, &lt->src);
      lt->src = vgtmp;
    }

    //    LTAchangeType(lta, LINEAR_VOX_TO_VOX);

    /* apply lta to the ctrl volume */
    mri_tmp = MRIalloc(mri_ctrl->width,
                       mri_ctrl->height,
                       mri_ctrl->depth,
                       mri_ctrl->type) ;
    MRIcopyHeader(mri_in, mri_tmp) ;
    // this function doesn't do NEAREST at all!!
    // I found the bug, in LTAtransformInterp()
    mri_tmp = LTAtransformInterp(mri_ctrl, mri_tmp, lta, SAMPLE_NEAREST);

    MRIfree(&mri_ctrl);
    mri_ctrl = mri_tmp;

    if (mask1_fname != NULL && mask2_fname == NULL) {
      printf("map mask for tp1 to get mask for tp2 ...\n");
      mri_mask2 = MRIalloc(mri_in->width,
                           mri_in->height,
                           mri_in->depth,
                           mri_mask1->type) ;
      MRIcopyHeader(mri_in, mri_mask2) ;

      mri_mask2 = LTAtransformInterp(mri_mask1,
                                     mri_mask2,
                                     lta,
                                     SAMPLE_NEAREST);
      mask2_set = 1;
      if (debug_flag)
        MRIwrite(mri_mask2, "mri_mask2.mgz");
    } else if (mask2_fname != NULL && mask1_fname == NULL) {
      printf("map mask for tp2 to get mask for tp1 ...\n");
      //need to invert lta first
      m_tmp = lta->xforms[0].m_L ;
      lta->xforms[0].m_L = MatrixInverse(lta->xforms[0].m_L, NULL) ;
      MatrixFree(&m_tmp) ;
      lt = &lta->xforms[0];

      //copyVolGeom(&lt->dst, &vgtmp);
      vgtmp = lt->dst;
      //copyVolGeom(&lt->src, &lt->dst);
      lt->dst = lt->src;
      //copyVolGeom(&vgtmp, &lt->src);
      lt->src = vgtmp;

      mri_mask1 = MRIalloc(mri_T1->width,
                           mri_T1->height,
                           mri_T1->depth,
                           mri_mask2->type) ;
      MRIcopyHeader(mri_T1, mri_mask1) ;

      mri_mask1 = LTAtransformInterp(mri_mask2,
                                     mri_mask1,
                                     lta,
                                     SAMPLE_NEAREST);
      mask1_set = 1;
      if (debug_flag)
        MRIwrite(mri_mask1, "mri_mask1.mgz");
    }

    if (lta_src)
      MRIfree(&lta_src);
    if (lta_dst)
      MRIfree(&lta_dst);

    if (lta)
      LTAfree(&lta);
  }   /* if (xform_fname != NULL) */

  if (debug_flag) {
    //    MRIwrite(mri_snr, "snr.mgz");
    MRIwrite(mri_ctrl, "ctrl.mgz");
  }

  if (mask1_set == 0) {
    //create mask1
    mri_mask1 = MRIalloc(mri_T1->width,
                         mri_T1->height,
                         mri_T1->depth,
                         MRI_UCHAR) ;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++) {
          if (MRIgetVoxVal(mri_T1, x, y, z, 0) < noise_threshold) {
            MRIvox(mri_mask1,x,y,z) = 0;
          } else
            MRIvox(mri_mask1,x,y,z) = 1;
        }
  }

  if (mask2_set == 0) {
    //create mask2
    mri_mask2 = MRIalloc(mri_in->width,
                         mri_in->height,
                         mri_in->depth,
                         MRI_UCHAR) ;
    for (z=0; z < depth; z++)
      for (y=0; y< height; y++)
        for (x=0; x < width; x++) {
          if (MRIgetVoxVal(mri_in, x, y, z, 0) < noise_threshold) {
            MRIvox(mri_mask2,x,y,z) = 0;
          } else
            MRIvox(mri_mask2,x,y,z) = 1;
        }
  }


#if 0
  /* compute the mean and std of T1 volume */
  /* Using only high SNR points */
  mri_snr = MRIalloc(mri_T1->width, mri_T1->height, mri_T1->depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_T1, mri_snr) ;

  h->bin_size = bin_size = 0.5;
  for (bin_no = 0; bin_no < nbins; bin_no++)
    h->bins[bin_no] = (bin_no)*bin_size;

  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        if (MRIgetVoxVal(mri_T1, x, y, z, 0) < noise_threshold) {
          MRIFvox(mri_snr,x,y,z) = 0;
          continue;
        }
        mean = 0;
        std = 0;
        count = 0;
        for (i=-1; i<=1; i++)
          for (j=-1; j<=1; j++)
            for (k=-1;k<=1;k++) {
              cx = x+i;
              cy = y+j, cz = z+k;

              if (cx < 0 ||
                  cx >= width ||
                  cy < 0 ||
                  cy >= height ||
                  cz < 0 ||
                  cz >= depth) continue;
              count++;
              value = MRIgetVoxVal(mri_T1, cx, cy, cz, 0);
              mean += value;
              std += value*value;
            }

        mean /= (count + 1e-30);
        std /=  (count + 1e-30);

        std = std - mean *mean;

        if (std <= 0)
          std = 0;
        value = mean/sqrt(std);

        MRIFvox(mri_snr,x,y,z) = value;
        bin_no = nint((float)value/(float)bin_size);
        if (bin_no >= nbins) bin_no = nbins - 1;
        h->counts[bin_no]++;
      }

  for (num = 0.0f, b = h->nbins - 1; b >= 1; b --) {
    num += h->counts[b];
    if (num > 20000) /* this may make me only use WM points,
                        is it good to use only WM to compute
                        scale of intensity?? */
      break;
  }

  printf("using SNR threshold %2.3f at bin %d\n", h->bins[b], b);

  mean1 = 0;
  std1 = 0;
  count1 = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        if (MRIgetVoxVal(mri_T1, x, y, z, 0) < noise_threshold) {
          continue;
        }

        value = MRIFvox(mri_snr,x,y,z);

        if (value < h->bins[b]) continue;

        value = MRIgetVoxVal(mri_T1, x, y, z, 0);
        count1++;
        mean1 += value;
        std1 += value*value;

      }

  MRIfree(&mri_snr);
#else
  printf("compute mean and std of tp1 volume within masked area...\n");
  mean1 = 0;
  std1 = 0;
  count1 = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        if (MRIgetVoxVal(mri_mask1, x, y, z, 0) <= 1e-30) {
          continue;
        }

        value = MRIgetVoxVal(mri_T1, x, y, z, 0);


        count1++;
        mean1 += value;
        std1 += value*value;

      }

#endif

  mean1 /= (count1 + 1e-30);
  std1 /= (count1 + 1e-30);
  std1 = std1 - mean1*mean1;
  if (std1 <= 0)
    printf("warning: negative std for T1 volume. \n");
  else
    printf("mean and variance for tp1 volume are %g and %g\n", mean1, std1);

  printf("now compute SNR and stats for input volume ... \n");
  mri_snr = MRIalloc(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT) ;
  MRIcopyHeader(mri_in, mri_snr) ;
  //HISTOclear(h,h);
  //h->bin_size = bin_size = 0.5;
  //for (bin_no = 0; bin_no < nbins; bin_no++)
  //  h->bins[bin_no] = (bin_no)*bin_size;

  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        if (MRIgetVoxVal(mri_in, x, y, z, 0) < noise_threshold) {
          MRIFvox(mri_snr,x,y,z) = 0;
          continue;
        }
        mean = 0;
        std = 0;
        count = 0;
        for (i=-1; i<=1; i++)
          for (j=-1; j<=1; j++)
            for (k=-1;k<=1;k++) {
              cx = x+i;
              cy = y+j, cz = z+k;

              if (cx < 0 ||
                  cx >= width ||
                  cy < 0 ||
                  cy >= height ||
                  cz < 0 ||
                  cz >= depth) continue;
              count++;
              value = MRIgetVoxVal(mri_in, cx, cy, cz, 0);
              mean += value;
              std += value*value;
            }

        mean /= (count + 1e-30);
        std /=  (count + 1e-30);

        std = std - mean *mean;

        if (std <= 0)
          std = 0;
        value = mean/sqrt(std);

        MRIFvox(mri_snr,x,y,z) = value;
        //bin_no = nint((float)value/(float)bin_size);
        //if (bin_no >= nbins) bin_no = nbins - 1;
        //h->counts[bin_no]++;
      }

#if 0
  for (num = 0.0f, b = h->nbins - 1; b >= 1; b --) {
    num += h->counts[b];
    if (num > 20000) /* this may make me only use WM points, is it good to
                        use only WM to compute scale of intensity?? */
      break;
  }

  printf("using SNR threshold %2.3f at bin %d\n", h->bins[b], b);

  mean2 = 0;
  std2 = 0;
  count2 = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        if (MRIgetVoxVal(mri_in, x, y, z, 0) < noise_threshold) {
          continue;
        }

        value = MRIFvox(mri_snr,x,y,z);

        if (value >= h->bins[b]) {
          count2++;
          mean2 += value;
          std2 += value*value;
        }
      }
#else
  printf("compute mean and std of tp2 volume within masked area\n");
  /* somehow mri_watershed seems to leave some unzero voxels around
  image border, so I will skip image boundaries
  no, that's not a problem of most recent mri_watershed;
  something wrong previously */
  mean2 = 0;
  std2 = 0;
  count2 = 0;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        if (MRIgetVoxVal(mri_mask2, x, y, z, 0) <= 1e-30) {
          continue;
        }

        value = MRIgetVoxVal(mri_in, x, y, z, 0);

        count2++;
        mean2 += value;
        std2 += value*value;
      }
#endif

  mean2 /= (count2 + 1e-30);
  std2 /= (count2 + 1e-30);
  std2 = std2 - mean2*mean2;
  if (std2 <= 0)
    printf("warning: negative std for input volume. \n");
  else
    printf("mean and variance for input tp2 volume are %g and %g\n",
           mean2, std2);

  //compute intensity scale
  slope = sqrt(std1/std2);
  offset = mean1 - slope*mean2;

  printf("scale input volume by %g x + %g\n", slope, offset);
  // first change mri_in to FLOAT type
  mri_tmp = MRIchangeType(mri_in, MRI_FLOAT, 0, 1.0, 1);
  MRIfree(&mri_in);
  mri_in = mri_tmp;
  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        value = MRIFvox(mri_in, x, y, z);
        MRIFvox(mri_in, x, y, z) = value*slope + offset;
      }


  //  printf("compute SNR map of tp2 volume\n"); //already done above
  //  mri_snr = MRIalloc(mri_ctrl->width,
  //     mri_ctrl->height, mri_ctrl->depth, MRI_FLOAT) ;

  for (z=0; z < depth; z++)
    for (y=0; y< height; y++)
      for (x=0; x < width; x++) {
        if (MRIgetVoxVal(mri_in, x, y, z, 0) < noise_threshold) {
          //    MRIFvox(mri_snr,x,y,z) = 0;
          continue;
        }

        value = MRIFvox(mri_snr,x,y,z);

        if (value < 20) MRIvox(mri_ctrl, x, y, z) = 0;
        else if (MRIvox(mri_ctrl, x, y, z) > 0) {
          MRIvox(mri_ctrl, x, y, z) = 1;
        }
      }
  if (debug_flag) {
    MRIwrite(mri_snr, "snr.mgz");
    //    MRIwrite(mri_ctrl, "ctrl.mgz");
  }

  // SNR >= 20 seems a good threshold

  // Now use ctrl points to normalize tp2
  printf("normalize tp2...\n");
  mri_bias = MRIbuildBiasImage(mri_in, mri_ctrl, NULL, bias_sigma) ;

  for (z = 0 ; z < depth ; z++) {
    for (y = 0 ; y < height ; y++) {
      for (x = 0 ; x < width ; x++) {
        src = MRIgetVoxVal(mri_in, x, y, z, 0) ;
        bias = MRIgetVoxVal(mri_bias, x, y, z, 0) ;
        if (!bias)   /* should never happen */
          norm = (float)src ;
        else
          norm = (float)src * 110.0 / (float)bias ;
        if (norm > 255.0f && mri_out->type == MRI_UCHAR)
          norm = 255.0f ;
        else if (norm < 0.0f && mri_out->type == MRI_UCHAR)
          norm = 0.0f ;
        MRIsetVoxVal(mri_out, x, y, z, 0, norm) ;
      }
    }
  }

  printf("writing normalized volume to %s...\n", argv[2]) ;
  MRIwrite(mri_out, argv[2]);

  MRIfree(&mri_in);
  MRIfree(&mri_bias);
  MRIfree(&mri_out);
  MRIfree(&mri_T1);
  MRIfree(&mri_ctrl);
  MRIfree(&mri_snr);
  //HISTOfree(&h);

  exit(0);

}  /*  end main()  */


static void usage(int exit_val) {
  FILE *fout;

  fout = (exit_val ? stderr : stdout);

  fprintf(fout, "usage: %s <input vol> <normalized vol>\n", Progname);
  fprintf(fout, "This program uses the controls points of "
          "tp1 to help normalize tp2 \n") ;
  fprintf(fout, "Options are: \n") ;
  fprintf(fout, "\t -T1 %%s: T1 volume for tp1 (actually the "
          "volume to perform mri_normalize on for tp1) \n");
  fprintf(fout, "\t -mask1 %%s: brain mask for tp1 (will be "
          "mapped to tp2 by applying the xform if -mask2 not specified) \n");
  fprintf(fout, "\t -mask2 %%s: brain mask for tp2 (will be mapped to tp1 by "
          "applying inverse-xform if -mask1 not specified) \n");
  fprintf(fout, "\t -threshold %%d: threshold for background "
          "(default = 1.0) \n");
  fprintf(fout, "\t -ctrl %%s: ctrl point volume for tp1 \n") ;
  fprintf(fout, "\t -xform %%s: LTA transform that aligns tp1 to tp2 \n");
  fprintf(fout, "\t -invert  reversely apply -xform \n");
  fprintf(fout, "\t -lta_src %%s:  source volume for -xform "
          "(if not available from the xform file) \n");
  fprintf(fout, "\t -lta_dst %%s:  target volume for -xform "
          "(if not available from the xform file) \n");

  exit(exit_val);
}  /*  end usage()  */


static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    debug_flag = 1;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)...\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "T1")) {
    tp1_T1_fname  = argv[2];
    printf("Use file %s for tp1's T1 volume \n", tp1_T1_fname) ;
    nargs = 1;
  } else if (!stricmp(option, "threshold")) {
    noise_threshold  = atof(argv[2]);
    printf("Take values < %g as background for both volumes \n",
           noise_threshold) ;
    nargs = 1;
  } else if (!stricmp(option, "mask1")) {
    mask1_fname  = (argv[2]);
    printf("Use volume %s as brain mask for tp1, and will be mapped to\n"
           "tp2 to get rough mask for tp2 in computing "
           "intensity normalization \n", mask1_fname) ;
    nargs = 1;
  } else if (!stricmp(option, "mask2")) {
    mask2_fname  = (argv[2]);
    printf("Use volume %s as brain mask for tp2, and will be mapped to \n"
           "tp1 to get rough mask for tp1 in computing "
           "intensity normalization \n", mask2_fname) ;
    nargs = 1;
  } else if (!stricmp(option, "ctrl")) {
    tp1_ctrl_fname  = argv[2];
    printf("Use file %s for control points volume for tp1\n", tp1_ctrl_fname) ;
    nargs = 1;
  } else if (!stricmp(option, "xform")) {
    xform_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "transform file name is %s\n", xform_fname);
  } else if (!stricmp(option, "invert")) {
    invert = 1;
    fprintf(stderr, "Inversely apply the given LTA transform\n");
  } else if (!stricmp(option, "lta_src") ||
             !stricmp(option, "src")
            ) {
    fprintf
    (stderr,
     "src volume for the given transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the src volume...\n");
    lta_src = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_src) {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else if (!stricmp(option, "lta_dst") ||
             !stricmp(option, "dst")
            ) {
    fprintf
    (stderr,
     "dst volume for the transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the dst volume...\n");
    lta_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_dst) {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  } else switch (toupper(*option)) {
    case '?':
    case 'U':
      usage(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
