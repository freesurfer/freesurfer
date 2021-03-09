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


/* mri_transform_to_COR.c */
/* Convert an input volume to COR and at the same time applying
 * a registration matrix. In other words, perform the registration
 * and resampling the registered volume in COR format.
 * The purpose of this program is to perform registration and resampling
 * together to reduce number of interpolations
 */
/* Note that the volume data (raw data) itself is meaningless without
 * specifying its i_to_r; especially when need to compare multiple volumes
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "transform.h"
#include "version.h"
#define MYFZERO(f)     (fabs(f) < 0.0001F)
#define SAMPLE_BSPLINE 5
//#define DBL_EPSILON 1e-10


LTA *ltaReadFileEx(const char *fname);
int MYvg_isEqual(const VOL_GEOM *vg1, const VOL_GEOM *vg2);

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(int) ;
static void print_version(void) ;
static int debug_flag = 0;
static int sinchalfwindow = 6;
static int SplineDegree = 3;
const char *Progname ;
static char *out_like_fname = NULL ;
static int invert_flag = 0 ;
//static int InterpMethod = SAMPLE_BSPLINE;
static int InterpMethod = SAMPLE_TRILINEAR;
static char *transform_fname = NULL;
static int transform_flag = 0;
static float scale = 1.0;
//static int out_type = 3; /* MRI_FLOAT */
static int out_type = 0; /* MRI_COR */
static int noscale = 0;
static int autoscale = 0; /* automatically scale histogram peak to 110 */

float thred_low = 0.0;
float thred_high = 1.0;

MRI *lta_src = NULL;
MRI *lta_dst = NULL;

double InterpolatedValue
(
  MRI *Bcoeff,     /* input B-spline array of coefficients */
  double x,      /* x coordinate where to interpolate */
  double y,      /* y coordinate where to interpolate */
  double  z,           /* z coordinate where to interpolate */
  long SplineDegree /* degree of the spline model */
);

int  SamplesToCoefficients
(
  MRI *mri_vol,  /* in-place processing */
  long SplineDegree/* degree of the spline model */
);

MRI *MRIlinearTransformInterpBSpline(MRI *mri_src, MRI *mri_dst, MATRIX *mA,
                                     int splinedegree);

int
main(int argc, char *argv[])
{
  char        **av, *in_vol, *out_vol;
  int         ac, nargs;

  MRI         *mri_in, *mri_out, *mri_tmp ;
  LTA         *lta = 0;
  MATRIX *i_to_r_src = 0; /* src geometry of the input LTA */
  MATRIX *V_to_V = 0; /* Final voxel-to-voxel transform */
  MATRIX *r_to_i_dst = 0; /* dst geometry of the input LTA */
  MATRIX *m_tmp = 0;
  MATRIX *i_to_r_reg = 0; /* i_to_r of the volume after registration */
  MATRIX *r_to_i_out = 0; /* r_to_i of the final output volume */
  VOL_GEOM vgm_in;
  int x, y, z;
  double maxV, minV, value;
  //  MATRIX *i_to_r, *r_to_i;

  nargs = handleVersionOption(argc, argv, "mri_transform_to_COR");
  if (nargs && argc - nargs == 1)
    usage_exit (0);
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

  if (argc < 3)
    usage_exit(0) ;

  in_vol = argv[1] ;
  out_vol = argv[2] ;

  printf("reading volume from %s...\n", in_vol) ;
  mri_in = MRIread(in_vol) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname,
              in_vol) ;

  /* Convert mri_in to float type */
  /* double would be more accurate */
  if (mri_in->type != MRI_FLOAT)
  {
    printf("Input volume type is %d\n", mri_in->type);
    printf("Change input volume to float type for convenience and accuracy");
    mri_tmp = MRIchangeType(mri_in, MRI_FLOAT, 0, 1.0, 1);
    MRIfree(&mri_in);
    mri_in = mri_tmp; //swap
  }

  /* Get input volume geometry, which is needed to compute i_to_r
   * and r_to_i of input volume. Note that i_to_r and r_to_i assumed
   * a certain prespecified c_r, c_a, c_s
   */
  getVolGeom(mri_in, &vgm_in);

  maxV = -10000.0;
  minV = 10000.0;
  for (z=0; z < mri_in->depth; z++)
    for (y=0; y< mri_in->height; y++)
      for (x=0; x < mri_in->width; x++)
      {
        if (MRIFvox(mri_in, x, y, z) > maxV )
          maxV = MRIFvox(mri_in, x, y,z) ;
        if (MRIFvox(mri_in, x, y, z) < minV )
          minV = MRIFvox(mri_in, x, y,z) ;
      }

  printf("Input volume has max = %g, min =%g\n", maxV, minV);

  printf("Scale input volume by %g \n", scale);

  maxV = -10000.0;
  minV = 10000.0;
  for (z=0; z < mri_in->depth; z++)
    for (y=0; y< mri_in->height; y++)
      for (x=0; x < mri_in->width; x++)
      {
        MRIFvox(mri_in, x, y, z) *= scale;
        if (MRIFvox(mri_in, x, y, z) > maxV )
          maxV = MRIFvox(mri_in, x, y,z) ;
        if (MRIFvox(mri_in, x, y, z) < minV )
          minV = MRIFvox(mri_in, x, y,z) ;
      }

  printf("Input volume after scaling has max = %g, min =%g\n", maxV, minV);

  /* Try to compute the Voxel_to_Voxel transform from the input volume
   * and the registration target/reference volume!
   * If no registration is involved, vox_to_vox is simply identity
   */
  /* Things become more complicated when allowing inverse transform */
  if (transform_flag)
  {
    int transform_type;

    printf("INFO: Applying transformation from file %s...\n",
           transform_fname);
    transform_type =  TransformFileNameType(transform_fname);

    /* Read in LTA transform file name */
    if (transform_type == MNI_TRANSFORM_TYPE ||
        transform_type == TRANSFORM_ARRAY_TYPE ||
        transform_type == REGISTER_DAT ||
        transform_type == FSLREG_TYPE
       )
    {

      printf("Reading transform ...\n");
      lta = LTAreadEx(transform_fname) ;

      if (!lta)
        ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
                  Progname, transform_fname) ;

      if (transform_type == FSLREG_TYPE)
      {
        if (lta_src == 0 || lta_dst == 0)
        {
          fprintf(stderr, "ERROR: fslmat does not have information on the src and dst volumes\n");
          fprintf(stderr, "ERROR: you must give options '-src' and '-dst' to specify the src and dst volume infos for the registration\n");
        }


        LTAmodifySrcDstGeom(lta, lta_src, lta_dst); // add src and dst information
        //The following is necessary to interpret FSLMAT correctly!!!
        LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      }
      if (lta->xforms[0].src.valid == 0)
      {
        if (lta_src == 0)
        {
          fprintf(stderr, "The transform does not have the valid src volume info.\n");
          fprintf(stderr, "Either you give src volume info by option -src or\n");
          fprintf(stderr, "make the transform to have the valid src info.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }
        else
        {
          LTAmodifySrcDstGeom(lta, lta_src, NULL); // add src information
        }
      }
      if (lta->xforms[0].dst.valid == 0)
      {
        if (lta_dst == 0)
        {
          fprintf(stderr, "The transform does not have the valid dst volume info.\n");
          fprintf(stderr, "Either you give src volume info by option -dst or\n");
          fprintf(stderr, "make the transform to have the valid dst info.\n");
          fprintf(stderr, "If the dst was average_305, then you can set\n");
          fprintf(stderr, "environmental variable USE_AVERAGE305 true\n");
          fprintf(stderr, "instead.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }
        else
        {
          LTAmodifySrcDstGeom(lta, NULL, lta_dst); // add  dst information
        }
      }


      // The following procedure aims to apply an LTA computed from COR format to a volume in non-COR format, or vice versa, as long as they share the same RAS
      // first change to LINEAR RAS_TO_RAS using old info
      if (lta->type != LINEAR_RAS_TO_RAS)
      {
        LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      }

      // now possiblly reset the src and dst
      if (lta_src != NULL)
      {
        //always trust the user
        LTAmodifySrcDstGeom(lta, lta_src, NULL);
      }
      if (lta_dst != NULL)
      {
        //always trust the user
        LTAmodifySrcDstGeom(lta, NULL, lta_dst);
      }

      if (lta->type == LINEAR_RAS_TO_RAS)
      {
        /* Convert it to VOX_TO_VOX */
        /* VOXELsrc_to_VOXELdst = R2Vdst*R2Rlta*V2Rsrc */
        /* Note whether the input should be identical to src or dst here depends
         * on whether the LTA here is the direct or inverse transform
         */
        i_to_r_src = vg_i_to_r(&lta->xforms[0].src);
        r_to_i_dst = vg_r_to_i(&lta->xforms[0].dst);

        if (!r_to_i_dst || !i_to_r_src)
          ErrorExit(ERROR_BADFILE, "%s: failed to extract volume geometries from input LTA file",Progname);
        m_tmp = MatrixMultiply(lta->xforms[0].m_L, i_to_r_src, NULL);
        V_to_V = MatrixMultiply(r_to_i_dst, m_tmp, NULL);
        MatrixFree(&m_tmp);

        MatrixFree(&i_to_r_src);
        MatrixFree(&r_to_i_dst);
      }
    }
    else
    {
      fprintf(stderr, "unknown transform type in file %s\n",
              transform_fname);
      exit(1);
    }

    if (invert_flag)
    {
      /* Geometry of input volume should match that of the dst of the LTA */
      if (MYvg_isEqual(&lta->xforms[0].dst, &vgm_in) == 0)
      {
        ErrorExit(ERROR_BADFILE, "%s: dst volume of lta doesn't match that of input volume",Progname);
      }

      i_to_r_reg = vg_i_to_r(&lta->xforms[0].src);

      if (!i_to_r_reg)
        ErrorExit(ERROR_BADFILE, "%s: failed to extract i_to_r of registered volume from LTA",Progname);

      m_tmp =  MatrixInverse(V_to_V, NULL);
      if (!m_tmp)
        ErrorExit(ERROR_BADPARM, "%s: transform is singular!", Progname);

      MatrixFree(&V_to_V);
      V_to_V = m_tmp;
    }
    else
    {
      /* Geometry of input volume should match that of the src of the LTA */
      if (MYvg_isEqual(&lta->xforms[0].src, &vgm_in) == 0)
      {
        ErrorExit(ERROR_BADFILE, "%s: src volume of lta doesn't match that of input volume",Progname);
      }

      i_to_r_reg = vg_i_to_r(&lta->xforms[0].dst);

      if (!i_to_r_reg)
        ErrorExit(ERROR_BADFILE, "%s: failed to extract i_to_r of registered volume from LTA",Progname);
    }

  }
  else
  {
    /* No registration transform need be applied */
    V_to_V = MatrixIdentity(4, NULL);
    i_to_r_reg = extract_i_to_r(mri_in);
    if (!i_to_r_reg)
      ErrorExit(ERROR_BADFILE, "%s: failed to extract i_to_r from input volume",Progname);
  }

  /* Now need to find the vox-to-vox transformation between registered volume
   * (or input volume itself if no registration involved) and the output
   * volume, either in COR format or as the out-like volume
   */
  /* Given a volume with a certain i_to_r, we need to compute the necessary
   * vox-to-voxel transform to change its i_to_r to like another volume.
   * The vox-to-vox is equal to R2V(r_to_i)_likevol*i_to_r_current_vol.
   */
  if (out_like_fname)
  {
    mri_tmp = MRIread(out_like_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read template volume from %s",out_like_fname) ;

    /* out_type = mri_tmp->type; */

    /* specify the out-type to float initially so as not to lose accuracy
     * during reslicing, will change type to correct type later.
     */
    mri_out = MRIalloc(mri_tmp->width, mri_tmp->height, mri_tmp->depth, MRI_FLOAT) ;

    MRIcopyHeader(mri_tmp, mri_out) ;
    MRIfree(&mri_tmp);
  }
  else  /* assume output is in COR format */
  {
    mri_out = MRIalloc(256, 256, 256, MRI_FLOAT) ;
    /* out_type = MRI_UCHAR; */

    /* Who says MRIlinearTransformInterp will change the header??
     * I don't think so!
     */
    //E/ set xyzc_ras to coronal ones.. - these'll get zorched
    //by MRIlinearTransformInterp() - copy again later - is there
    //any use in having them here now?  yes, so we can pass mri_out
    //to the ras2vox fns.


    mri_out->imnr0 = 1; /* what's this? */
    mri_out->imnr1 = 256; /* what's this? */
    mri_out->thick = 1.0;
    mri_out->ps = 1.0; /* what's this? */
    mri_out->xsize = mri_out->ysize = mri_out->zsize = 1.0;
    mri_out->xstart = mri_out->ystart = mri_out->zstart = -128.0;
    mri_out->xend = mri_out->yend = mri_out->zend = 128.0;
    mri_out->x_r =-1;
    mri_out->y_r = 0;
    mri_out->z_r = 0;
    mri_out->x_a = 0;
    mri_out->y_a = 0;
    mri_out->z_a = 1;
    mri_out->x_s = 0;
    mri_out->y_s =-1;
    mri_out->z_s = 0;

    /* In this case, the RAS itself is not fully determined, i.e., c_ras.
     * It's quite arbitrary, different values just change the final
     * sitting of the volume inside the RAS system.
     */
    /* NO! The C_RAS has to be set correctly, depending which target
     * volume the previous Vox_to_Vox transformation assumes!
     * When a registration is involved, the target volume is either
     * the src of LTA (direct) or the dst (inverse transform). When
     * just change format, the target volume is the input itself!!
     */
    if (transform_flag)
    {
      if (invert_flag)
      {
        mri_out->c_r = lta->xforms[0].src.c_r;
        mri_out->c_a = lta->xforms[0].src.c_a;
        mri_out->c_s = lta->xforms[0].src.c_s;

      }
      else
      {
        mri_out->c_r = lta->xforms[0].dst.c_r;
        mri_out->c_a = lta->xforms[0].dst.c_a;
        mri_out->c_s = lta->xforms[0].dst.c_s;
      }
    }
    else
    {
      mri_out->c_r = mri_in->c_r;
      mri_out->c_a = mri_in->c_a;
      mri_out->c_s = mri_in->c_s;
    }

    mri_out->ras_good_flag=1; /* What does this flag mean ? */

    /* since output is just transformed input */
    MRIcopyPulseParameters(mri_in, mri_out) ;

  }

  /* Compute the final input-to-output VOX_to_VOX transformation matrix */
  r_to_i_out = extract_r_to_i(mri_out);

  m_tmp = MatrixMultiply(r_to_i_out, i_to_r_reg, NULL);
  V_to_V = MatrixMultiply(m_tmp, V_to_V, V_to_V);
  MatrixFree(&m_tmp);

  printf("InterpMethod = %d\n", InterpMethod);

  /* Modify the MyMRIlinearTr... if I want to implement my cubic-B-spline
   * interpolation method. Otherwise, unnecessary
   */
  /* mri_out = MyMRIlinearTransformInterp(mri_in, mri_out, V_to_V, InterpMethod); */
  if (InterpMethod == SAMPLE_BSPLINE)
    mri_out = MRIlinearTransformInterpBSpline(mri_in, mri_out, V_to_V,
              SplineDegree);
  else
    mri_out = MRIlinearTransformInterp(mri_in, mri_out, V_to_V, InterpMethod);

  maxV = -10000.0;
  minV = 10000.0;
  for (z=0; z < mri_out->depth; z++)
    for (y=0; y< mri_out->height; y++)
      for (x=0; x < mri_out->width; x++)
      {
        if (MRIFvox(mri_out, x, y, z) > maxV )
          maxV = MRIFvox(mri_out, x, y,z) ;
        if (MRIFvox(mri_out, x, y, z) < minV )
          minV = MRIFvox(mri_out, x, y,z) ;
      }

  if (autoscale)
  {
    noscale = 1;

    /* compute histogram of output volume */
    HISTOGRAM *h, *hsmooth ;
    float fmin, fmax, val, peak, smooth_peak;
    int i, nbins, bin;

    fmin = minV;
    fmax = maxV;
    if (fmin < 0) fmin = 0;
    nbins = 256 ;
    h = HISTOalloc(nbins) ;
    hsmooth = HISTOcopy(h, NULL) ;
    HISTOclear(h, h) ;
    h->bin_size = (fmax-fmin)/255.0 ;

    for (i = 0 ; i < nbins ; i++)
      h->bins[i] = (i+1)*h->bin_size ;

    for (z=0; z < mri_out->depth; z++)
      for (y=0; y< mri_out->height; y++)
        for (x=0; x < mri_out->width; x++)
        {
          val = MRIFvox(mri_out, x, y, z);
          if (val <= 0) continue;

          bin = nint((val - fmin)/h->bin_size);
          if (bin >= h->nbins)
            bin = h->nbins-1;
          else if (bin < 0)
            bin = 0;

          h->counts[bin] += 1.0;
        }
    HISTOfillHoles(h) ;
    HISTOsmooth(h, hsmooth, 5)  ;
    peak =
      hsmooth->bins[HISTOfindHighestPeakInRegion(h, 1, h->nbins)] ;
    //   smooth_peak =
    //  hsmooth->bins[HISTOfindHighestPeakInRegion(hsmooth, 1, hsmooth->nbins)] ;

    smooth_peak =
      hsmooth->bins[HISTOfindLastPeak(hsmooth, 5, 0.8)] ;

    /*
      bin = nint((smooth_peak - fmin)/hsmooth->bin_size) ;

      printf("Highest peak has count = %d\n", (int)hsmooth->counts[bin]);

      bin = nint((420 - fmin)/hsmooth->bin_size) ;

      printf("bin at 420 has count = %d\n", (int)hsmooth->counts[bin]);
    */

    scale =  110.0/smooth_peak;
    printf("peak of output volume is %g, smooth-peak is %g, multiply by %g to scale it to 110\n", peak, smooth_peak, scale);
    for (z=0; z < mri_out->depth; z++)
      for (y=0; y< mri_out->height; y++)
        for (x=0; x < mri_out->width; x++)
        {
          val = MRIFvox(mri_out, x, y, z);
          MRIFvox(mri_out, x, y, z) = val*scale;
        }

  }


  printf("Output volume (before type-conversion) has max = %g, min =%g\n", maxV, minV);

  /* Finally change type to desired */
  if (mri_out->type != out_type)
  {
    printf("Change output volume to type %d\n", out_type);
    /* I need to modify the MIRchangeType function to make sure
     * it does roundoff instead of simple truncation!
     */
    /* Note if the last flag is set to 1, then it won't do scaling
       and small float numbers will become zero after convert to
       BYTE
    */
    if (out_type == 0 && noscale == 1)
    {
      //convert data to UCHAR
      mri_tmp = MRIalloc(mri_out->width, mri_out->height, mri_out->depth, out_type) ;
      MRIcopyHeader(mri_out, mri_tmp);
      for (z=0; z < mri_out->depth; z++)
        for (y=0; y< mri_out->height; y++)
          for (x=0; x < mri_out->width; x++)
          {
            value = floor(MRIgetVoxVal(mri_out, x, y, z, 0) + 0.5);
            if (value < 0 ) value = 0;
            if (value > 255) value = 255;
            MRIvox(mri_tmp,x,y,z) = (unsigned char)value;
          }
    }
    else
      mri_tmp = MRIchangeType(mri_out, out_type, thred_low, thred_high, noscale);

    MRIfree(&mri_out);
    mri_out = mri_tmp; //swap
  }

  MRIwrite(mri_out, out_vol) ;

  MRIfree(&mri_in);
  MRIfree(&mri_out);

  if (lta_src)
    MRIfree(&lta_src);
  if (lta_dst)
    MRIfree(&lta_dst);

  MatrixFree(&V_to_V);

  if (!r_to_i_out)
    MatrixFree(&r_to_i_out);

  if (!i_to_r_reg)
    MatrixFree(&i_to_r_reg);

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
  if (!stricmp(option, "debug_voxel"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    debug_flag = 1;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)...\n", Gx, Gy, Gz) ;
  }
  else if (!stricmp(option, "scaling"))
  {
    scale = atof(argv[2]);
    nargs = 1;
    printf("scale input by %g\n", scale);
  }
  else if (!stricmp(option, "low"))
  {
    thred_low = atof(argv[2]);
    nargs = 1;
    printf("Lower threshold for histogram scaling is %g\n", thred_low);
  }
  else if (!stricmp(option, "high"))
  {
    thred_high = atof(argv[2]);
    nargs = 1;
    printf("Higher threshold for histogram scaling is %g\n", thred_high);
  }
  else if (!stricmp(option, "noscale"))
  {
    noscale = 1;
    printf("do not do histogram scaling at output stage\n");
  }
  else if (!stricmp(option, "autoscale"))
  {
    autoscale = 1;
    printf("automatically scale output volume histo peak to 110 \n");
  }
  else if (!stricmp(option, "out_type"))
  {
    out_type = atoi(argv[2]);
    nargs = 1;
    printf("Output type is %d\n", out_type);
  }
  else if (!stricmp(option, "out_like") ||
           !stricmp(option, "like"))
  {
    out_like_fname = argv[2];
    nargs = 1;
    printf("Shape the output like the volume in file %s\n", out_like_fname);
  }
  else if (!stricmp(option, "transform") ||
           !stricmp(option, "at"))
  {
    transform_fname = argv[2];
    transform_flag = 1;
    nargs = 1;
    printf("Apply transformation specified by file %s\n", transform_fname);
  }
  else if (!stricmp(option, "lta_src") ||
           !stricmp(option, "src")
          )
  {
    fprintf(stderr, "src volume for the given transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the src volume...\n");
    lta_src = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_src)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else if (!stricmp(option, "lta_dst") ||
           !stricmp(option, "dst")
          )
  {
    fprintf(stderr, "dst volume for the transform (given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the dst volume...\n");

    lta_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);

    if (!lta_dst)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else if (!stricmp(option, "invert_transform") ||
           !stricmp(option, "ait"))
  {
    transform_fname = argv[2];
    transform_flag = 1;
    invert_flag = 1;
    nargs = 1;
    printf("Apply inversely the transformation specified by file %s\n", transform_fname);
  }
  /*E* Interpolation method.  Default is trilinear, other options are
    nearest, cubic, sinc.  You can say -foo or -interp foo.  For sinc,
    you can say -interp sinc 3 or -interp sinc -hw 3 or -sinc 3 or
    -sinc -hw 3.  Maybe -hw 3 should imply sinc, but right now it
    doesn't.  */

  else if (!stricmp(option, "st") ||
           !stricmp(option, "sample") ||
           !stricmp(option, "sample_type") ||
           !stricmp(option, "interp"))
  {
    nargs = 1;

    if (!stricmp(argv[2], "bspline"))
      InterpMethod = SAMPLE_BSPLINE;
    else
      InterpMethod = MRIinterpCode(argv[2]) ;

    if (InterpMethod==SAMPLE_SINC)
    {
      if ((argc<4) || !strncmp(argv[3],"-",1)) /*E* i.e. no sinchalfwindow value supplied */
      {
        printf("using sinc interpolation (default windowwidth is 6)\n");
      }
      else
      {
        sinchalfwindow = atoi(argv[3]);
        nargs = 2;
        printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
      }
    }
    else if (InterpMethod == SAMPLE_BSPLINE)
    {
      if ((argc<4) || !strncmp(argv[3],"-",1)) /* i.e. no spline-degree value supplied */
      {
        printf("using BSPline interpolation (default Bspline degree is 3)\n");
      }
      else
      {
        SplineDegree = atoi(argv[3]);
        nargs = 2;
        printf("using BSpline interpolation with degree of %d\n", SplineDegree);
      }

    }
  }
  else if (!stricmp(option, "sinc"))
  {
    InterpMethod = SAMPLE_SINC;
    if ((argc<3) || !strncmp(argv[2],"-",1)) /*E* i.e. no sinchalfwindow value supplied */
    {
      printf("using sinc interpolation (default windowwidth is 6)\n");
    }
    else
    {
      sinchalfwindow = atoi(argv[2]);
      nargs = 1;
      printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
    }
  }
  else if (!stricmp(option, "bspline"))
  {
    InterpMethod = SAMPLE_BSPLINE;
    if ((argc<3) || !strncmp(argv[2],"-",1)) /*E* i.e. no spline degree supplied */
    {
      printf("using cubic-bspline interpolation \n");
    }
    else
    {
      SplineDegree = atoi(argv[2]);
      nargs = 1;
      printf("using B-spline interpolation with degree of %d\n", SplineDegree);
    }
  }
  else if (!stricmp(option, "sinchalfwindow") ||
           !stricmp(option, "hw"))
  {
    InterpMethod = SAMPLE_SINC;
    sinchalfwindow = atoi(argv[2]);
    nargs = 1;
    printf("using sinc interpolation with windowwidth of %d\n", 2*sinchalfwindow);
  }
  else if (!stricmp(option, "trilinear"))
  {
    InterpMethod = SAMPLE_TRILINEAR;
    printf("using trilinear interpolation\n");
  }
  else if (!stricmp(option, "cubic"))
  {
    InterpMethod = SAMPLE_CUBIC;
    printf("using cubic interpolation\n");
  }
  else if (!stricmp(option, "nearest"))
  {
    InterpMethod = SAMPLE_NEAREST;
    printf("using nearest-neighbor interpolation\n");
  }
  else switch (toupper(*option))
    {
    case '?':
    case 'U':
      usage_exit(0) ;
      break ;
    default:
      printf("unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}


static void
usage_exit(int exit_val)
{
  printf("usage: %s  <input> <output> \n", Progname);
  printf("this program convert input to output while at the same time\n") ;

  printf("applying an LTA transform if available. \n") ;

  printf("Options includes:\n");
  printf("\t -at fname to apply an LTA transform\n");
  printf("\t -ait fname to inversely apply an LTA transform \n");
  printf("\t -like fname to force the output be shaped like this volume \n");
  printf("\t -interp: resample type:<trilinear,nearest,sinc,cubic,bspline> \n");
  printf("\t -scaling #: scale the input values by the number \n");
  printf("\t -autoscale #: automatically scale the output volume peak to 110 \n");
  printf("\t -noscale: donot scale output during type-conversion\n");
  printf("\t -out_type #: specify output volume type to be the number\n");
  printf("\t -high #: value near 1 to specify higher-end percentage for histogram guided float-to-byte conversion \n");
  printf("\t -low #: value near 0 for lower-end of histogram scaling \n");

  print_version();

  exit(exit_val);
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

/* Actually no need to modify this function, since I will only use the float
 * type here, so the roundoff I added will never take effect.
 * What I need to modify is the MRIchangeType function!
 */
MRI *MRIlinearTransformInterpBSpline(MRI *mri_src, MRI *mri_dst, MATRIX *mA,
                                     int splinedegree)
{
  int    y1, y2, y3, width, height, depth ;
  VECTOR *v_X, *v_Y ;   /* original and transformed coordinate systems */
  MATRIX *mAinv ;     /* inverse of mA */
  double   val, x1, x2, x3 ;
  MRI *mri_Bcoeff;

  mAinv = MatrixInverse(mA, NULL) ;      /* will sample from dst back to src */
  if (!mAinv)
    ErrorReturn(NULL, (ERROR_BADPARM,
                       "MRIlinearTransformBSpline: xform is singular")) ;

  if (!mri_dst) mri_dst = MRIclone(mri_src, NULL) ;
  else          MRIclear(mri_dst) ; /* set all values to zero */

  if (mri_src->type != MRI_FLOAT)
    mri_Bcoeff = MRIchangeType(mri_src, MRI_FLOAT, 0, 1.0, 1);
  else
    mri_Bcoeff = MRIcopy(mri_src, NULL);

  /* convert between a representation based on image samples */
  /* and a representation based on image B-spline coefficients */
  if (SamplesToCoefficients(mri_Bcoeff, splinedegree))
  {
    ErrorReturn(NULL, (ERROR_BADPARM, "Change of basis failed\n"));
  }

  printf("Direct B-spline Transform Finished. \n");

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */

  v_Y->rptr[4][1] = 1.0f ;
  for (y3 = 0 ; y3 < mri_dst->depth ; y3++)
  {
    V3_Z(v_Y) = y3 ;
    for (y2 = 0 ; y2 < mri_dst->height ; y2++)
    {
      V3_Y(v_Y) = y2 ;
      for (y1 = 0 ; y1 < mri_dst->width ; y1++)
      {
        V3_X(v_Y) = y1 ;
        MatrixMultiply(mAinv, v_Y, v_X) ;

        x1 = V3_X(v_X) ;
        x2 = V3_Y(v_X) ;
        x3 = V3_Z(v_X) ;

        if (x1 <= -0.5 || x1 >= (width - 0.5) || x2 <= -0.5 ||
            x2 >= (height - 0.5) || x3 <= -0.5 || x3 >= (depth-0.5))
          val = 0.0;
        else
          val = InterpolatedValue(mri_Bcoeff, x1, x2, x3, splinedegree);

        switch (mri_dst->type)
        {
        case MRI_UCHAR:
          if (val <-0.5) val = -0.5;
          if (val > 254.5) val = 254.5;
          MRIvox(mri_dst,y1,y2,y3) = (BUFTYPE)nint(val+0.5) ;
          break ;
        case MRI_SHORT:
          MRISvox(mri_dst,y1,y2,y3) = (short)nint(val+0.5) ;
          break ;
        case MRI_FLOAT:
          MRIFvox(mri_dst,y1,y2,y3) = (float)(val) ;
          break ;
        case MRI_INT:
          MRIIvox(mri_dst,y1,y2,y3) = nint(val+0.5) ;
          break ;
        default:
          ErrorReturn(NULL,
                      (ERROR_UNSUPPORTED,
                       "MRIlinearTransformBSpline: unsupported dst type %d",
                       mri_dst->type)) ;
          break ;
        }
      }
    }
  }

  MatrixFree(&v_X) ;
  MatrixFree(&mAinv) ;
  MatrixFree(&v_Y) ;

  MRIfree(&mri_Bcoeff);
  return(mri_dst) ;
}

LTA *ltaReadFileEx(const char *fname)
{
  FILE             *fp;
  LINEAR_TRANSFORM *lt ;
  int              i, nxforms, type ;
  char             line[STRLEN], *cp ;
  LTA              *lta ;

  fp = fopen(fname,"r");
  if (fp==NULL)
    ErrorReturn(NULL,
                (ERROR_BADFILE, "ltaReadFile(%s): can't open file",fname));
  cp = fgetl(line, 199, fp) ;
  if (cp == NULL)
  {
    fclose(fp) ;
    ErrorReturn(NULL, (ERROR_BADFILE, "ltaReadFile(%s): can't read data",fname));
  }
  sscanf(cp, "type      = %d\n", &type) ;
  cp = fgetl(line, 199, fp) ;
  sscanf(cp, "nxforms   = %d\n", &nxforms) ;
  lta = LTAalloc(nxforms, NULL) ;
  lta->type = type ;
  for (i = 0 ; i < lta->num_xforms ; i++)
  {
    lt = &lta->xforms[i] ;
    fscanf(fp, "mean      = %f %f %f\n", &lt->x0, &lt->y0, &lt->z0) ;
    fscanf(fp, "sigma     = %f\n", &lt->sigma) ;
    MatrixAsciiReadFrom(fp, lt->m_L) ;
  }
  // oh, well this is the added part
  for (i=0; i < lta->num_xforms; i++)
  {
    if (fgets(line, 199, fp))
    {
      if (strncmp(line, "src volume info", 15)==0)
      {
        char *p;
        readVolGeom(fp, &lta->xforms[i].src);
        p = fgets(line, 199, fp);
        if (strncmp(line, "dst volume info", 15)==0)
          readVolGeom(fp, &lta->xforms[i].dst);
      }
    }
  }
  fclose(fp) ;
  return(lta) ;
}

int MYvg_isEqual(const VOL_GEOM *vg1, const VOL_GEOM *vg2)
{
  if (vg1->valid == vg2->valid)
    if (vg1->width == vg2->width)
      if (vg1->height == vg2->height)
        if (vg1->depth == vg2->depth)
          if (MYFZERO(vg1->xsize - vg2->xsize))
            if (MYFZERO(vg1->ysize - vg2->ysize))
              if (MYFZERO(vg1->zsize - vg2->zsize))
                if (MYFZERO(vg1->x_r - vg2->x_r))
                  if (MYFZERO(vg1->x_a - vg2->x_a))
                    if (MYFZERO(vg1->x_s - vg2->x_s))
                      if (MYFZERO(vg1->y_r - vg2->y_r))
                        if (MYFZERO(vg1->y_a - vg2->y_a))
                          if (MYFZERO(vg1->y_s - vg2->y_s))
                            if (MYFZERO(vg1->z_r - vg2->z_r))
                              if (MYFZERO(vg1->z_a - vg2->z_a))
                                if (MYFZERO(vg1->z_s - vg2->z_s))
                                  if (MYFZERO(vg1->c_r - vg2->c_r))
                                    if (MYFZERO(vg1->c_a - vg2->c_a))
                                      if (MYFZERO(vg1->c_s - vg2->c_s))
                                        return 1;
  return 0;
}


/*************************************************************************
       The following functions are used for B-Spline interpolation
************************************************************************/
/*--------------------------------------------------------------------------*/
double InitialAntiCausalCoefficient
(
  double c[],  /* coefficients */
  long DataLength, /* number of samples or coefficients */
  double z   /* actual pole */
)

{ /* begin InitialAntiCausalCoefficient */

  /* this initialization corresponds to mirror boundaries */
  return((z / (z * z - 1.0)) * (z * c[DataLength - 2L] + c[DataLength - 1L]));
} /* end InitialAntiCausalCoefficient */

/*--------------------------------------------------------------------------*/
double InitialCausalCoefficient
(
  double c[],  /* coefficients */
  long DataLength, /* number of coefficients */
  double z,   /* actual pole */
  double Tolerance /* admissible relative error */
)

{ /* begin InitialCausalCoefficient */

  double Sum, zn, z2n, iz;
  long n, Horizon;

  /* this initialization corresponds to mirror boundaries */
  Horizon = DataLength;
  if (Tolerance > 0.0)
  {
    Horizon = (long)ceil(log(Tolerance) / log(fabs(z)));
  }
  if (Horizon < DataLength)
  {
    /* accelerated loop */
    zn = z;
    Sum = c[0];
    for (n = 1L; n < Horizon; n++)
    {
      Sum += zn * c[n];
      zn *= z;
    }
    return(Sum);
  }
  else
  {
    /* full loop */
    zn = z;
    iz = 1.0 / z;
    z2n = pow(z, (double)(DataLength - 1L));
    Sum = c[0] + z2n * c[DataLength - 1L];
    z2n *= z2n * iz;
    for (n = 1L; n <= DataLength - 2L; n++)
    {
      Sum += (zn + z2n) * c[n];
      zn *= z;
      z2n *= iz;
    }
    return(Sum / (1.0 - zn * zn));
  }
} /* end InitialCausalCoefficient */

/*--------------------------------------------------------------------------*/
static void  ConvertToInterpolationCoefficients
(
  double c[],  /* input samples --> output coefficients */
  long DataLength, /* number of samples or coefficients */
  double z[],  /* poles */
  long NbPoles, /* number of poles */
  double Tolerance /* admissible relative error */
)

{ /* begin ConvertToInterpolationCoefficients */

  double Lambda = 1.0;
  long n, k;

  /* special case required by mirror boundaries */
  if (DataLength == 1L)
  {
    return;
  }
  /* compute the overall gain */
  for (k = 0L; k < NbPoles; k++)
  {
    Lambda = Lambda * (1.0 - z[k]) * (1.0 - 1.0 / z[k]);
  }
  /* apply the gain */
  for (n = 0L; n < DataLength; n++)
  {
    c[n] *= Lambda;
  }
  /* loop over all poles */
  for (k = 0L; k < NbPoles; k++)
  {
    /* causal initialization */
    c[0] = InitialCausalCoefficient(c, DataLength, z[k], Tolerance);
    /* causal recursion */
    for (n = 1L; n < DataLength; n++)
    {
      c[n] += z[k] * c[n - 1L];
    }
    /* anticausal initialization */
    c[DataLength - 1L] = InitialAntiCausalCoefficient(c, DataLength, z[k]);
    /* anticausal recursion */
    for (n = DataLength - 2L; 0 <= n; n--)
    {
      c[n] = z[k] * (c[n + 1L] - c[n]);
    }
  }
} /* end ConvertToInterpolationCoefficients */


/*--------------------------------------------------------------------------*/
static void  GetColumn
(
  MRI   *mri_vol,  /* input image volume */
  long x,   /* x coordinate of the selected line */
  long    z,          /* Slice number of the line */
  double Line[]  /* output linear array */
)

{ /* begin GetColumn */

  long y;

  for (y = 0L; y < mri_vol->height; y++)
  {
    Line[y] = (double)MRIgetVoxVal(mri_vol, x, y, z, 0);
  }
} /* end GetColumn */

/*--------------------------------------------------------------------------*/
static void  GetRow
(
  MRI *mri_vol,  /* input image volume */
  long y,   /* y coordinate of the selected line */
  long    z,          /* Slice number of the line */
  double Line[]  /* output linear array */
)

{ /* begin GetRow */

  long x;

  for (x = 0L; x < mri_vol->width; x++)
  {
    Line[x] = (double)MRIgetVoxVal(mri_vol, x, y, z, 0);
  }
} /* end GetRow */

/*--------------------------------------------------------------------------*/
static void GetVertical
(
  MRI *mri_vol,      /* input image volume */
  long x,            /* x coordinate of the selected line */
  long y,            /* y coordinate of the selected line */
  double Line[]     /* output linear array */
)
{ /* begin GetVertical */
  long z;

  for (z = 0L; z < mri_vol->depth; z++)
  {
    Line[z] = (double) MRIgetVoxVal(mri_vol, x, y, z, 0);
  }

} /* end GetVertical */


/*--------------------------------------------------------------------------*/
static void  PutColumn
(
  MRI *mri_vol,  /* output image volume */
  long x,   /* x coordinate of the selected line */
  long    z,          /* Slice number of the line */
  double Line[]  /* input linear array */
)

{ /* begin PutColumn */

  long y;

  for (y = 0L; y < mri_vol->height; y++)
  {
    MRIsetVoxVal(mri_vol, x, y, z, 0, (float)Line[y]);
  }
} /* end PutColumn */

/*--------------------------------------------------------------------------*/
static void  PutRow
(
  MRI *mri_vol,  /* output image volume */
  long y,   /* y coordinate of the selected line */
  long    z,          /* Slice number of the line */
  double *Line  /* input linear array */
)

{ /* begin PutRow */

  long x;

  for (x = 0L; x < mri_vol->width; x++)
  {
    MRIsetVoxVal(mri_vol, x, y, z, 0, (float)Line[x]);
  }
} /* end PutRow */

/*--------------------------------------------------------------------------*/
static void PutVertical
(
  MRI *mri_vol,      /* output image volume */
  long x,            /* x coordinate of the selected line */
  long y,            /* y coordinate of the selected line */
  double Line[]     /* output linear array */
)
{ /* begin PutVertical */
  long z;
  for (z = 0L; z < mri_vol->depth; z++)
  {
    MRIsetVoxVal(mri_vol, x, y, z, 0, (float)Line[z]);
  }
  return;
} /* end PutVertical */

/*****************************************************************************
 * Definition of extern procedures
 ****************************************************************************/
/*--------------------------------------------------------------------------*/
int  SamplesToCoefficients
(
  MRI *mri_vol,  /* in-place processing */
  long SplineDegree /* degree of the spline model */
)

{ /* begin SamplesToCoefficients */

  double  *Line;
  double  Pole[2];
  long NbPoles;
  long x, y, z;
  long Width, Height, Depth;

  /* recover the poles from a lookup table */
  switch (SplineDegree)
  {
  case 2L:
    NbPoles = 1L;
    Pole[0] = sqrt(8.0) - 3.0;
    break;
  case 3L:
    NbPoles = 1L;
    Pole[0] = sqrt(3.0) - 2.0;
    break;
  case 4L:
    NbPoles = 2L;
    Pole[0] = sqrt(664.0 - sqrt(438976.0)) + sqrt(304.0) - 19.0;
    Pole[1] = sqrt(664.0 + sqrt(438976.0)) - sqrt(304.0) - 19.0;
    break;
  case 5L:
    NbPoles = 2L;
    Pole[0] = sqrt(135.0 / 2.0 - sqrt(17745.0 / 4.0)) + sqrt(105.0 / 4.0)
              - 13.0 / 2.0;
    Pole[1] = sqrt(135.0 / 2.0 + sqrt(17745.0 / 4.0)) - sqrt(105.0 / 4.0)
              - 13.0 / 2.0;
    break;
  default:
    printf("Invalid spline degree\n");
    return(1);
  }

  Width = mri_vol->width;
  Height = mri_vol->height;
  Depth = mri_vol->depth;

  /* convert the image samples into interpolation coefficients */
  /* in-place separable process, along x */
  Line = (double *)malloc((size_t)(Width * (long)sizeof(double)));
  if (Line == (double *)NULL)
  {
    printf("Row allocation failed\n");
    return(1);
  }

  for (z = 0L; z < Depth; z++)
  {
    for (y = 0L; y < Height; y++)
    {
      GetRow(mri_vol, y, z, Line);
      ConvertToInterpolationCoefficients(Line, Width, Pole, NbPoles, DBL_EPSILON);
      PutRow(mri_vol, y, z, Line);

    }
  }
  free(Line);

  /* in-place separable process, along y */
  Line = (double *)malloc((size_t)(Height * (long)sizeof(double)));
  if (Line == (double *)NULL)
  {
    printf("Column allocation failed\n");
    return(1);
  }
  for (z = 0L; z < Depth; z++)
  {
    for (x = 0L; x < Width; x++)
    {
      GetColumn(mri_vol,  x, z, Line);
      ConvertToInterpolationCoefficients(Line, Height, Pole, NbPoles, DBL_EPSILON);
      PutColumn(mri_vol, x, z, Line);
    }
  }
  free(Line);

  /* in-place separable process, along z */
  Line = (double *)malloc((size_t)(Depth * (long)sizeof(double)));
  if (Line == (double *)NULL)
  {
    printf("Column allocation failed\n");
    return(1);
  }
  for (y = 0L; y < Height; y++)
  {
    for (x = 0L; x < Width; x++)
    {
      GetVertical(mri_vol, x, y, Line);
      ConvertToInterpolationCoefficients(Line, Depth, Pole, NbPoles, DBL_EPSILON);
      PutVertical(mri_vol, x, y, Line);
    }
  }
  free(Line);

  return(0);
} /* end SamplesToCoefficients */


double InterpolatedValue
(
  MRI *Bcoeff,     /* input B-spline array of coefficients */
  double x,      /* x coordinate where to interpolate */
  double y,      /* y coordinate where to interpolate */
  double  z,           /* z coordinate where to interpolate */
  long SplineDegree /* degree of the spline model */
)

{ /* begin InterpolatedValue */

  double xWeight[6], yWeight[6], zWeight[6];
  double interpolated;
  double w, w2, w4, t, t0, t1;
  long xIndex[6], yIndex[6], zIndex[6];
  long Height, Width, Depth, Height2, Width2, Depth2;
  long i, j, d, k; /* i,j, d are indices for x,y, and z respectively; k for spline index */
  long cz, cy, cx;

  Width = Bcoeff->width;
  Height = Bcoeff->height;
  Depth = Bcoeff->depth;

  Width2 = 2L * Width - 2L;
  Height2 = 2L * Height - 2L;
  Depth2 = 2L * Depth - 2L;

  /* compute the interpolation indexes */
  if (SplineDegree & 1L)
  {
    i = (long)floor(x) - SplineDegree / 2L;
    j = (long)floor(y) - SplineDegree / 2L;
    d = (long)floor(z) - SplineDegree / 2L;
    for (k = 0L; k <= SplineDegree; k++)
    {
      xIndex[k] = i++;
      yIndex[k] = j++;
      zIndex[k] = d++;
    }
  }
  else
  {
    i = (long)floor(x + 0.5) - SplineDegree / 2L;
    j = (long)floor(y + 0.5) - SplineDegree / 2L;
    d = (long)floor(z + 0.5) - SplineDegree / 2L;
    for (k = 0L; k <= SplineDegree; k++)
    {
      xIndex[k] = i++;
      yIndex[k] = j++;
      zIndex[k] = d++;
    }
  }

  /* compute the interpolation weights */
  switch (SplineDegree)
  {
  case 2L:
    /* x */
    w = x - (double)xIndex[1];
    xWeight[1] = 3.0 / 4.0 - w * w;
    xWeight[2] = (1.0 / 2.0) * (w - xWeight[1] + 1.0);
    xWeight[0] = 1.0 - xWeight[1] - xWeight[2];
    /* y */
    w = y - (double)yIndex[1];
    yWeight[1] = 3.0 / 4.0 - w * w;
    yWeight[2] = (1.0 / 2.0) * (w - yWeight[1] + 1.0);
    yWeight[0] = 1.0 - yWeight[1] - yWeight[2];
    /* z */
    w = z - (double)zIndex[1];
    zWeight[1] = 3.0 / 4.0 - w * w;
    zWeight[2] = (1.0 / 2.0) * (w - zWeight[1] + 1.0);
    zWeight[0] = 1.0 - yWeight[1] - zWeight[2];
    break;
  case 3L:
    /* x */
    w = x - (double)xIndex[1];
    xWeight[3] = (1.0 / 6.0) * w * w * w;
    xWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - xWeight[3];
    xWeight[2] = w + xWeight[0] - 2.0 * xWeight[3];
    xWeight[1] = 1.0 - xWeight[0] - xWeight[2] - xWeight[3];
    /* y */
    w = y - (double)yIndex[1];
    yWeight[3] = (1.0 / 6.0) * w * w * w;
    yWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - yWeight[3];
    yWeight[2] = w + yWeight[0] - 2.0 * yWeight[3];
    yWeight[1] = 1.0 - yWeight[0] - yWeight[2] - yWeight[3];
    /* z */
    w = z - (double)zIndex[1];
    zWeight[3] = (1.0 / 6.0) * w * w * w;
    zWeight[0] = (1.0 / 6.0) + (1.0 / 2.0) * w * (w - 1.0) - zWeight[3];
    zWeight[2] = w + zWeight[0] - 2.0 * zWeight[3];
    zWeight[1] = 1.0 - zWeight[0] - zWeight[2] - zWeight[3];
    break;
  case 4L:
    /* x */
    w = x - (double)xIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    xWeight[0] = 1.0 / 2.0 - w;
    xWeight[0] *= xWeight[0];
    xWeight[0] *= (1.0 / 24.0) * xWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    xWeight[1] = t1 + t0;
    xWeight[3] = t1 - t0;
    xWeight[4] = xWeight[0] + t0 + (1.0 / 2.0) * w;
    xWeight[2] = 1.0 - xWeight[0] - xWeight[1] - xWeight[3] - xWeight[4];
    /* y */
    w = y - (double)yIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    yWeight[0] = 1.0 / 2.0 - w;
    yWeight[0] *= yWeight[0];
    yWeight[0] *= (1.0 / 24.0) * yWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    yWeight[1] = t1 + t0;
    yWeight[3] = t1 - t0;
    yWeight[4] = yWeight[0] + t0 + (1.0 / 2.0) * w;
    yWeight[2] = 1.0 - yWeight[0] - yWeight[1] - yWeight[3] - yWeight[4];
    /* z */
    w = z - (double)zIndex[2];
    w2 = w * w;
    t = (1.0 / 6.0) * w2;
    zWeight[0] = 1.0 / 2.0 - w;
    zWeight[0] *= zWeight[0];
    zWeight[0] *= (1.0 / 24.0) * zWeight[0];
    t0 = w * (t - 11.0 / 24.0);
    t1 = 19.0 / 96.0 + w2 * (1.0 / 4.0 - t);
    zWeight[1] = t1 + t0;
    zWeight[3] = t1 - t0;
    zWeight[4] = zWeight[0] + t0 + (1.0 / 2.0) * w;
    zWeight[2] = 1.0 - zWeight[0] - zWeight[1] - zWeight[3] - zWeight[4];
    break;
  case 5L:
    /* x */
    w = x - (double)xIndex[2];
    w2 = w * w;
    xWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    xWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - xWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    xWeight[2] = t0 + t1;
    xWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    xWeight[1] = t0 + t1;
    xWeight[4] = t0 - t1;
    /* y */
    w = y - (double)yIndex[2];
    w2 = w * w;
    yWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    yWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - yWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    yWeight[2] = t0 + t1;
    yWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    yWeight[1] = t0 + t1;
    yWeight[4] = t0 - t1;
    /* z */
    w = z - (double)zIndex[2];
    w2 = w * w;
    zWeight[5] = (1.0 / 120.0) * w * w2 * w2;
    w2 -= w;
    w4 = w2 * w2;
    w -= 1.0 / 2.0;
    t = w2 * (w2 - 3.0);
    zWeight[0] = (1.0 / 24.0) * (1.0 / 5.0 + w2 + w4) - zWeight[5];
    t0 = (1.0 / 24.0) * (w2 * (w2 - 5.0) + 46.0 / 5.0);
    t1 = (-1.0 / 12.0) * w * (t + 4.0);
    zWeight[2] = t0 + t1;
    zWeight[3] = t0 - t1;
    t0 = (1.0 / 16.0) * (9.0 / 5.0 - t);
    t1 = (1.0 / 24.0) * w * (w4 - w2 - 5.0);
    zWeight[1] = t0 + t1;
    zWeight[4] = t0 - t1;
    break;
  default:
    printf("Invalid spline degree\n");
    return(0.0);
  }

  /* apply the mirror boundary conditions */
  if (Width == 1L)
  {
    for (k = 0L; k <= SplineDegree; k++)
    {
      xIndex[k] = 0L;
    }
  }
  if (Height == 1L)
  {
    for (k = 0L; k <= SplineDegree; k++)
    {
      yIndex[k] = 0L;
    }
  }
  if (Depth == 1L)
  {
    for (k = 0L; k <= SplineDegree; k++)
    {
      zIndex[k] = 0L;
    }
  }

  for (k = 0L; k <= SplineDegree; k++)
  {
    xIndex[k] = (xIndex[k] < 0L) ? (-xIndex[k] - Width2 * ((-xIndex[k]) / Width2))
                : (xIndex[k] - Width2 * (xIndex[k] / Width2));
    if (Width <= xIndex[k])
    {
      xIndex[k] = Width2 - xIndex[k];
    }
    yIndex[k] = (yIndex[k] < 0L) ? (-yIndex[k] - Height2 * ((-yIndex[k]) / Height2))
                : (yIndex[k] - Height2 * (yIndex[k] / Height2));
    if (Height <= yIndex[k])
    {
      yIndex[k] = Height2 - yIndex[k];
    }
    zIndex[k] = (zIndex[k] < 0L) ? (-zIndex[k] - Depth2 * ((-zIndex[k]) / Depth2))
                : (zIndex[k] - Depth2 * (zIndex[k] / Depth2));
    if (Depth <= zIndex[k])
    {
      zIndex[k] = Depth2 - zIndex[k];
    }
  }

  /* perform interpolation */
  interpolated = 0.0;
  for (d = 0L; d <=SplineDegree; d++)
  {
    t = 0;
    cz = zIndex[d];
    for (j = 0L; j <= SplineDegree; j++)
    {
      w = 0.0;
      cy = yIndex[j];
      for (i = 0L; i <= SplineDegree; i++)
      {
        cx = xIndex[i];
        w += xWeight[i] * MRIgetVoxVal(Bcoeff, cx, cy, cz, 0);
      }
      t += yWeight[j] * w;
    }
    interpolated += zWeight[d] * t;
  }

  return(interpolated);
} /* end InterpolatedValue */
