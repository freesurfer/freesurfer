/**
 * @brief applies a mask volume
 *
 * Usage: mri_mask [options] <in vol> <mask vol> <out vol>
 * This program applies a mask volume (typically skull stripped).
 * Options:
 *    -xform %s    apply LTA transform to align mask to input volume
 *    -invert      reversely apply -xform
 *    -lta_src %s  source volume for -xform
 *    -lta_dst %s  target volume for -xform
 *    -T #         threshold mask volume at #
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
#include <math.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "fio.h"
#include "version.h"
#include "transform.h"
#include "region.h"
#include "cma.h"


void usage(int exit_val);

const char *Progname;

static int  get_option(int argc, char *argv[]) ;

static int   InterpMethod = SAMPLE_NEAREST;
static int dilate = 0 ;

/* The following specifies the src and dst volumes
   of the input FSL/LTA transform */
MRI          *lta_src = 0;
MRI          *lta_dst = 0;

static int no_cerebellum = 0 ;
static int DoLH = 0;
static int DoRH = 0;

static int samseg = 0 ;
static  float out_val=0;
static int invert = 0 ;
static char *xform_fname = NULL;
static float threshold = -1e10;
int ThresholdSet = 0;
static int do_transfer=0;
static float transfer_val;
static int keep_mask_deletion_edits = 0; // if 1, keep mask voxels with value=1
int DoAbs = 0;
int DoBB = 0, nPadBB[6], DoBBEq = 0;
double bgnoisescale = 0;
int bgnoisesign=0;
int InvertMask = 0;
COLOR_TABLE *ctab=NULL;
int CropToFoV[3];
int DoCropToFoV =0;

int main(int argc, char *argv[])
{
  char **av;
  MRI *mri_in, *mri_mask, *mri_out, *mri_mask_orig, *mri_tmp;
  int nargs, ac, nmask;
  int x, y, z;
  float value;
  MRI_REGION *region;

  nargs = handleVersionOption(argc, argv, "mri_mask");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs ;

  Progname = argv[0];
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

  if (argc != 4)
  {
    printf("Incorrect number of arguments, argc = %d\n", argc);
    usage(1);
  }

  mri_in = MRIread(argv[1]) ;
  if (!mri_in)
    ErrorExit(ERROR_BADPARM, "%s: could not read source volume %s",
              Progname, argv[1]) ;
  if(mri_in->ct) ctab = CTABdeepCopy(mri_in->ct);
  mri_mask = MRIread(argv[2]) ;
  if (!mri_mask)
    ErrorExit(ERROR_BADPARM, "%s: could not read mask volume %s",
              Progname, argv[2]) ;

  if(InvertMask){
    double thresh = threshold;
    if(! ThresholdSet) thresh = 0.5;
    printf("Inverting and binarizing mask thresh = %g\n",thresh);
    int c,count=0;
    for(c=0; c < mri_mask->width; c++){
      int r,s;
      double v;
      for(r=0; r < mri_mask->height; r++){
	for(s=0; s < mri_mask->depth; s++){
	  v = MRIgetVoxVal(mri_mask,c,r,s,0);
	  if(DoAbs) v = fabs(v);
	  if(v < thresh) {
	    MRIsetVoxVal(mri_mask,c,r,s,0,1);
	    count++;
	  }
	  else MRIsetVoxVal(mri_mask,c,r,s,0,0);
	}
      }
    }
    printf("count %d\n",count);
    printf("Resetting threshold to 0.5\n");
    threshold = 0.5;
    ThresholdSet = 1;
  }

  printf("DoAbs = %d\n",DoAbs);

  /* Read transform and apply it to mri_mask */
  if (xform_fname != NULL && strcmp(xform_fname, "identity.nofile"))
  {
    printf("Applying transform to the mask volume\n");
    TRANSFORM *transform = TransformRead(xform_fname);
    
    if (!transform)
      ErrorExit(ERROR_NOFILE, "%s: could not read transform file %s",
                Progname, xform_fname) ;
  
    if (transform->type != MORPH_3D_TYPE)
    {
      LTA *tmp = (LTA *)transform->xform;
      LTA *lta = LTAreduce(tmp); // Apply full array, allocation.
      LTAfree(&tmp);
      if (lta->xforms[0].src.valid == 0) // For LTAs such as FSLREG_TYPE.
      {
        if (lta_src == 0)
        {
          fprintf(stderr,
                  "The transform does not have the valid src volume info.\n");
          fprintf(stderr,
                  "Either you give src volume info by option -lta_src or\n");
          fprintf(stderr,
                  "make the transform to have the valid src info.\n");
          ErrorExit(ERROR_BAD_PARM, "Bailing out...\n");
        }
        LTAmodifySrcDstGeom(lta, lta_src, NULL);
      }
      
      if (lta->xforms[0].dst.valid == 0)
      {
        if (lta_dst == 0)
        {
          fprintf(stderr,
                  "The transform does not have the valid dst volume info.\n");
          fprintf(stderr,
                  "Either you give src volume info by option -lta_dst or\n");
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
        LTAmodifySrcDstGeom(lta, NULL, lta_dst);
      }
      
      LTAchangeType(lta, LINEAR_VOX_TO_VOX); // Support more types.
      transform->type = LINEAR_VOX_TO_VOX;
      transform->xform = (void *)lta;
    } /* if (transform->type != MORPH_3D_TYPE) */
    
    vg_isEqual_Threshold = 10e-4; // Override, include/transform.h.
    VOL_GEOM vg_mask, vg_in, vg_transform_src, vg_transform_dst;
    TransformGetSrcVolGeom(transform, &vg_transform_src);
    TransformGetDstVolGeom(transform, &vg_transform_dst);
    getVolGeom(mri_mask, &vg_mask);
    getVolGeom(mri_in, &vg_in);
    if (!vg_isEqual(&vg_mask, &vg_transform_src) ||
        !vg_isEqual(&vg_in, &vg_transform_dst))
    {
      if (!vg_isEqual(&vg_mask, &vg_transform_dst) ||
          !vg_isEqual(&vg_in, &vg_transform_src))
      {
        ErrorExit(ERROR_BADPARM, "ERROR: transform geometry does not match");
      }
      printf("Inverting transform to match MRI geometries\n");
      // Pass MRI in case path changed (see GCAMfillInverse):
      TransformInvertReplace(transform, mri_in);
    }
    
    mri_tmp = TransformApplyType(transform, mri_mask, NULL, InterpMethod);
    MRIfree(&mri_mask);
    mri_mask = mri_tmp;
    TransformFree(&transform);
  } /* if (xform_fname != NULL) */
  
  if (lta_src) MRIfree(&lta_src);
  if (lta_dst) MRIfree(&lta_dst);
  mri_mask_orig = MRIcopy(mri_mask, NULL);

  // Threshold and binarize mask. Without binarization, this fails 
  // when values are less than 1
  if(ThresholdSet){
    nmask = 0;
    for (z = 0 ; z <mri_mask->depth ; z++) {
      for (y = 0 ; y < mri_mask->height ; y++) {
	for (x = 0 ; x < mri_mask->width ; x++) {
	  value = MRIgetVoxVal(mri_mask, x, y, z, 0);
	  if(DoAbs) value = fabs(value);
	  if(value <= threshold) MRIsetVoxVal(mri_mask,x,y,z,0,0);
	  else {
	    MRIsetVoxVal(mri_mask,x,y,z,0,1);
	    nmask ++;
	  }
	}
      }
    }
    printf("Found %d voxels in mask (pct=%6.2f)\n",nmask,
	   100.0*nmask/(mri_mask->width*mri_mask->height*mri_mask->depth));
  }
  if(samseg){
    nmask = 0;
    for (z = 0 ; z <mri_mask->depth ; z++) {
      for (y = 0 ; y < mri_mask->height ; y++) {
	for (x = 0 ; x < mri_mask->width ; x++) {
	  value = MRIgetVoxVal(mri_mask, x, y, z, 0);
	  switch ((int)value)
	  {
	  case Unknown:
	  case CSF:
	  case Skull:
	  case Air:
	  case Head_ExtraCerebral:
	  case CSF_ExtraCerebral:
	    MRIsetVoxVal(mri_mask,x,y,z,0,0);
	    break;
	  default:
	    MRIsetVoxVal(mri_mask,x,y,z,0,1);
	    nmask ++;
	  }
	}
      }
    }
    printf("Found %d voxels in mask (pct=%6.2f)\n",nmask,
	   100.0*nmask/(mri_mask->width*mri_mask->height*mri_mask->depth));
  }

  if (DoLH || DoRH)  // mri_mask is an aseg - convert it to a mask
  {
    MRI *mri_aseg = mri_mask ;
    int label ;

    mri_mask = MRIclone(mri_aseg, NULL) ;

    for (z = 0 ; z <mri_mask->depth ; z++) {
      for (y = 0 ; y < mri_mask->height ; y++) {
	for (x = 0 ; x < mri_mask->width ; x++) {
	  if (x == Gx && y == Gy && z == Gz)
	    DiagBreak() ;
 	  label = (int)MRIgetVoxVal(mri_aseg, x, y, z, 0);
	  if (IS_CEREBELLUM(label) && no_cerebellum)
	    label = 0 ;
	  if (DoLH && IS_LH_CLASS(label))
	    MRIsetVoxVal(mri_mask, x, y, z, 0, 1) ;
	  else if (DoRH && IS_RH_CLASS(label))
	    MRIsetVoxVal(mri_mask, x, y, z, 0, 1) ;
	  else
	    MRIsetVoxVal(mri_mask, x, y, z, 0, 0) ;
	  
	}
      }
    }
  }

  if(DoBB || DoBBEq || DoCropToFoV){
    if(DoBB){
      printf("Computing bounding box, npad = %d, %d, %d, %d, %d, %d\n",
	     nPadBB[0],nPadBB[1],nPadBB[2],nPadBB[3],nPadBB[4],nPadBB[5]);
      region = REGIONgetBoundingBoxM(mri_mask,nPadBB);
    }
    if(DoBBEq) {
      printf("Computing bounding box with ncols=nrows, npad = %d\n",nPadBB[0]);
      region = REGIONgetBoundingBoxEqOdd(mri_mask,nPadBB[0]);
    }
    if(DoCropToFoV){
      region = MRIcropToFoV(mri_mask, CropToFoV[0], CropToFoV[1], CropToFoV[2], threshold);
    }
    REGIONprint(stdout, region);
    mri_tmp = MRIextractRegion(mri_mask, NULL, region);
    if(mri_tmp == NULL) exit(1);
    MRIfree(&mri_mask);
    mri_mask = mri_tmp;
    mri_tmp = MRIextractRegion(mri_in, NULL, region);
    if(mri_tmp == NULL) exit(1);
    MRIfree(&mri_in);
    mri_in = mri_tmp;
  }

  int mask=0;
  if (do_transfer)  {
    mask = (int)transfer_val;
    out_val = transfer_val;
  }
  if (dilate > 0)  {
    int i ;
    for (i = 0 ; i < dilate ; i++)
      MRIdilate(mri_mask, mri_mask);
  }

  printf("maskval=%d, outval=%g\n",mask,out_val);
  mri_out = MRImask(mri_in, mri_mask, NULL, mask, out_val) ;
  if (!mri_out)  {
    printf("%s: stripping failed\n", Progname) ;
    exit(1);
  }

  if (keep_mask_deletion_edits)  {
    MRI *mri_tmp ;
    mri_tmp = MRImask(mri_out, mri_mask_orig, NULL, 1, 1) ; // keep voxels = 1
    if (!mri_tmp) {
      printf("%s: stripping failed on keep_mask_deletion_edits\n", Progname) ;
      exit(1);
    }
    MRIfree(&mri_out) ; mri_out = mri_tmp ;
  }

  if(bgnoisescale > 0){
    printf("Adding background noise %g %d\n",bgnoisescale,bgnoisesign);
    for (z = 0 ; z <mri_mask->depth ; z++) {
      for (y = 0 ; y < mri_mask->height ; y++) {
	for (x = 0 ; x < mri_mask->width ; x++) {
 	  int m  = (int)MRIgetVoxVal(mri_mask, x, y, z, 0);
	  if(m) continue;
	  double v = (drand48() - bgnoisesign/2.0)*bgnoisescale;
	  MRIsetVoxVal(mri_out, x, y, z, 0, v) ;
	  
	}
      }
    }


  }

  if(ctab) mri_out->ct = CTABdeepCopy(ctab);				     

  printf("Writing masked volume to %s...", argv[3]) ;
  int err = MRIwrite(mri_out, argv[3]);
  if(err) exit(1);
  printf("done.\n") ;

  MRIfree(&mri_in);
  MRIfree(&mri_mask);
  MRIfree(&mri_out);
  MRIfree(&mri_mask_orig);

  exit(0);

}  /*  end main()  */

#include "mri_mask.help.xml.h"
void usage(int exit_val)
{
  outputHelpXml(mri_mask_help_xml,
                mri_mask_help_xml_len);
  exit(exit_val);

}  /*  end usage()  */

/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stdout, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage"))
  {
    usage(1) ;
  }
  else if (!stricmp(option, "version"))
  {
    print_version() ;
  }
  else if (!stricmp(option, "samseg"))
  {
    samseg = 1 ;
    printf("masking nonbrain using samseg segmentation\n") ;
  }
  else if (!stricmp(option, "abs"))
  {
    DoAbs = 1;
  }
  else if (!stricmp(option, "invert"))
  {
    InvertMask = 1;
  }
  else if (!stricmp(option, "no-invert"))
  {
    InvertMask = 0;
  }
  else if (!stricmp(option, "lh"))
  {
    DoLH = 1;
  }
  else if (!stricmp(option, "no_cerebellum"))
  {
    no_cerebellum = 1;
  }
  else if (!stricmp(option, "rh"))
  {
    DoRH = 1;
  }
  else if (!stricmp(option, "dilate"))
  {
    dilate = atoi(argv[2]) ;
    nargs = 1 ;
    printf("dilating mask %d times before applying\n", dilate);
  }
  else if (!stricmp(option, "xform"))
  {
    xform_fname = argv[2];
    nargs = 1;
    fprintf(stderr, "transform file name is %s\n", xform_fname);
  }
  else if (!stricmp(option, "T")
           || !stricmp(option, "threshold")
          )
  {
    threshold = (float)atof(argv[2]);
    nargs = 1;
    fprintf(stderr, "threshold mask volume at %g\n", threshold);
    ThresholdSet = 1;
  }
  else if (!stricmp(option, "crop-to-fov"))
  {
    CropToFoV[0] = (int)atoi(argv[2]);
    CropToFoV[1] = (int)atoi(argv[3]);
    CropToFoV[2] = (int)atoi(argv[4]);
    DoCropToFoV = 1;
    nargs = 3;
  }
  else if (!stricmp(option, "BB") || !stricmp(option, "crop") || !stricmp(option, "boundingbox"))
  {
    nPadBB[0] = (int)atoi(argv[2]);
    nPadBB[1] = (int)atoi(argv[2]);
    nPadBB[2] = (int)atoi(argv[2]);
    nPadBB[3] = (int)atoi(argv[2]);
    nPadBB[4] = (int)atoi(argv[2]);
    nPadBB[5] = (int)atoi(argv[2]);
    DoBB = 1;
    nargs = 1;
    printf("cropping npad = %d\n",nPadBB[0]);
  }
  else if (!stricmp(option, "BBM")|| !stricmp(option, "boundingboxm"))
  {
    nPadBB[0] = (int)atoi(argv[2]);
    nPadBB[1] = (int)atoi(argv[2]);
    nPadBB[2] = (int)atoi(argv[3]);
    nPadBB[3] = (int)atoi(argv[3]);
    nPadBB[4] = (int)atoi(argv[4]);
    nPadBB[5] = (int)atoi(argv[4]);
    DoBB = 1;
    nargs = 3;
    printf("bounding box M npad = %d, %d, %d\n",nPadBB[0],nPadBB[2],nPadBB[4]);
  }
  else if (!stricmp(option, "BBMM")|| !stricmp(option, "boundingboxmm"))
  {
    nPadBB[0] = (int)atoi(argv[2]);
    nPadBB[1] = (int)atoi(argv[3]);
    nPadBB[2] = (int)atoi(argv[4]);
    nPadBB[3] = (int)atoi(argv[5]);
    nPadBB[4] = (int)atoi(argv[6]);
    nPadBB[5] = (int)atoi(argv[7]);
    DoBB = 1;
    nargs = 6;
    printf("bounding box MM npad = %d, %d, %d, %d, %d, %d\n",
	   nPadBB[0],nPadBB[1],nPadBB[2],nPadBB[3],nPadBB[4],nPadBB[5]);
  }
  else if (!stricmp(option, "BBEQ"))
  {
    nPadBB[0] = (int)atoi(argv[2]);
    DoBBEq = 1;
    nargs = 1;
    printf("bounding box eq npad = %d\n",nPadBB[0]);
  }
  else if (!stricmp(option, "no-allow-diff-geom"))
  {
    setenv("FS_MRIMASK_ALLOW_DIFF_GEOM","0",1);
    nargs = 0;
  }
  else if (!stricmp(option, "transfer"))
  {
    do_transfer = 1;
    transfer_val = (float)atof(argv[2]);
    nargs = 1;
    fprintf(stderr, "transfer mask voxels=%g to dst vol\n", transfer_val);
  }
  else if (!stricmp(option, "invert"))
  {
    invert = 1;
    fprintf(stderr, "Inversely apply the given transform\n");
  }
  else if (!stricmp(option, "oval"))
  {
    out_val = atof(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "setting masked output voxels to %2.1f instead of 0\n", out_val) ;
  }
  else if (!stricmp(option, "lta_src") ||
           !stricmp(option, "src")
          )
  {
    fprintf(stderr, "src volume for the given transform "
            "(given by -xform) is %s\n",argv[2]);
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
    fprintf(stderr, "dst volume for the transform "
            "(given by -xform) is %s\n",argv[2]);
    fprintf(stderr, "Reading the dst volume...\n");
    lta_dst = MRIreadHeader(argv[2], MRI_VOLUME_TYPE_UNKNOWN);
    if (!lta_dst)
    {
      ErrorExit(ERROR_BADPARM, "Could not read file %s\n", argv[2]);
    }
    nargs = 1;
  }
  else if (!stricmp(option, "keep_mask_deletion_edits"))
  {
    keep_mask_deletion_edits = 1;
    fprintf(stderr, "Transferring mask edits ('1' voxels) to dst vol\n");
  }
  else if (!stricmp(option, "DEBUG_VOXEL"))
  {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging node (%d, %d, %d)\n", Gx,Gy,Gz) ;
  }
  else if (!stricmp(option, "bgnoise")){
    sscanf(argv[2],"%lf",&bgnoisescale);
    sscanf(argv[3],"%d",&bgnoisesign);
    nargs = 2;
  }
  else
  {
    fprintf(stderr, "unknown option %s\n", argv[1]) ;
    usage(1) ;
    exit(1) ;
  }

  return(nargs) ;
}

