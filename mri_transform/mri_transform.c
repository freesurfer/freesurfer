/**
 * @file  mri_transform.c
 * @brief program to transform an mri or cmat structure using a linear transform
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2013/04/21 17:10:41 $
 *    $Revision: 1.18 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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
#include "mri.h"
#include "proto.h"
#include "transform.h"
#include "version.h"
#include "gcamorph.h"
#include "cmat.h"
#include "transform.h"
#include "cma.h"

double MRIcomputeLinearTransformLabelDist(MRI *mri_src, MATRIX *mA, int label) ;
static char vcid[] = "$Id: mri_transform.c,v 1.18 2013/04/21 17:10:41 fischl Exp $";

//E/ For transformations: for case LINEAR_RAS_TO_RAS, we convert to
//vox2vox with MRIrasXformToVoxelXform() in mri.c; for case
//LINEAR_CORONAL_RAS_TO_CORONAL_RAS, this converts to vox2vox with
//MT_CoronalRasXformToVoxelXform(), below.  For LINEAR_VOX_TO_VOX, we
//don't know the output vox2ras tx, and just do the best we can
//(i.e. guess isotropic cor).

//E/ Eventually, maybe, this should be in more public use and live in mri.c:
static MATRIX *
MT_CoronalRasXformToVoxelXform(MRI *mri_in, MRI *mri_out,
                               MATRIX *m_ras2ras, MATRIX *m_vox2vox);

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

char *Progname ;
static int quiet_mode = 0 ;
static char *subject_name ;
static char *out_like_fname = NULL ;
static char *in_like_fname = NULL ;  // for cmat stuff
static int invert_flag = 0 ;
static int resample_type = SAMPLE_TRILINEAR ;
static int nlabels = 0 ;
static int labels[1000] ;

static int cmat_output_coords = LABEL_COORDS_VOXEL ;

int
main(int argc, char *argv[]) {
  char        **av, *in_vol, *out_vol, *xform_fname ;
  int         ac, nargs, i, ras_flag = 0, nxform_args ;
  MRI         *mri_in, *mri_out, *mri_tmp ;
  LTA         *lta ;
  MATRIX      *m, *m_total ;
  TRANSFORM   *transform = NULL ;

  VECTOR *mine;
  VECTOR *c;

#if 0
  fprintf(stderr, "Please use 'mri_convert -at xfm.'  mri_transform contains too many bugs.\n");
  exit(1);
#endif

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_transform.c,v 1.18 2013/04/21 17:10:41 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (((argc < 4) && (nlabels == 0)) || ((argc < 3) && (nlabels > 0)))
    usage_exit() ;

  in_vol = argv[1] ;
  out_vol = argv[argc-1] ;

  xform_fname = argv[argc-2] ;
  if ((strcmp(in_vol+strlen(in_vol)-5, ".cmat") == 0) ||
      (strcmp(in_vol+strlen(in_vol)-5, ".CMAT") == 0))
  {
    CMAT *cmat_in, *cmat_out ;
    MRI   *mri_in, *mri_out ;
    LTA   *lta ;

    if (in_like_fname == NULL)
      ErrorExit(ERROR_NOFILE, "%s: must specifiy --in_like MRI volume for cmat transform (use original conformed one)\n",
	Progname) ;

    mri_in = MRIread(in_like_fname) ;
    if (!mri_in)
      ErrorExit(ERROR_NOFILE, "%s: could not read template volume from %s",Progname, in_like_fname) ;

    if (out_like_fname == NULL)
      ErrorExit(ERROR_NOFILE, "%s: must specifiy --out_like MRI volume for cmat transform (use target volume)\n",
	Progname) ;
    mri_out = MRIread(out_like_fname) ;
    if (!mri_out)
      ErrorExit(ERROR_NOFILE, "%s: could not read template volume from %s",Progname, out_like_fname) ;

    printf("reading input CMAT files...\n") ;
    cmat_in = CMATread(in_vol) ;


#if 0
    if (stricmp(xform_fname, "voxel") == 0)
    {
      CMATtoVoxel(cmat_in, mri_out) ;
      printf("writing transformed cmat to %s\n", out_vol) ;
      CMATwrite(cmat_in, out_vol) ;
      exit(0) ;
    }
    else if (stricmp(xform_fname, "tkreg") == 0)
    {
      CMATtoTKreg(cmat_in, mri_out) ;
      printf("writing transformed cmat to %s\n", out_vol) ;
      CMATwrite(cmat_in, out_vol) ;
      exit(0) ;
    }
#endif
    transform = TransformRead(xform_fname) ;
    if (!transform)
      ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s",
                Progname, xform_fname) ;
    lta = (LTA *)(transform->xform) ;
    if (lta->type == LINEAR_COR_TO_COR)
    {
      MATRIX *M ;

      LTAsetVolGeom(lta, mri_in, mri_out) ;
      printf("src - %s\n", mri_in->fname) ;
      printf("dst - %s\n", mri_out->fname) ;
      printf("TKreg matrix\n") ;
      MatrixPrint(stdout, lta->xforms[0].m_L) ;
#if 0
      LTAchangeType(lta, LINEAR_RAS_TO_RAS) ;
#else
      M = MRItkReg2Native(mri_in, mri_out, lta->xforms[0].m_L) ;
      MatrixFree(&lta->xforms[0].m_L) ;
#if 0
      lta->xforms[0].m_L  = MatrixInverse(M, NULL) ;
      MatrixFree(&M) ;
#else
      lta->xforms[0].m_L  = M ;
#endif      
      lta->type = LINEAR_RAS_TO_RAS ;
#endif
      printf("scanner RAS matrix\n") ;
      MatrixPrint(stdout, lta->xforms[0].m_L) ;
      transform->type = lta->type ;
    }
    cmat_out = CMATtransform(cmat_in, transform, mri_in, mri_out, NULL) ;

    if (DIAG_VERBOSE_ON)
    {
      printf("writing before and after labels from %s (%d) to %s (%d)\n",
	     cma_label_to_name(cmat_in->labels[0]), cmat_in->labels[0],
	     cma_label_to_name(cmat_in->labels[3]), cmat_in->labels[3]) ;
      LabelWrite(cmat_in->splines[0][3], "before.label") ;
      LabelWrite(cmat_out->splines[0][3], "after.label") ;
    }
    switch (cmat_output_coords)
    {
    case LABEL_COORDS_VOXEL:
      if (cmat_out->coords != LABEL_COORDS_VOXEL)
	CMATtoVoxel(cmat_out, mri_out) ;
      break ;
    case LABEL_COORDS_SCANNER_RAS:
      if (cmat_out->coords != LABEL_COORDS_SCANNER_RAS)
	CMATtoScannerRAS(cmat_out, mri_out) ;
      break ;
    default:
    case LABEL_COORDS_TKREG_RAS:
      if (cmat_out->coords != LABEL_COORDS_TKREG_RAS)
	CMATtoScannerRAS(cmat_out, mri_out) ;
      break ;
    }

    printf("writing transformed cmat to %s\n", out_vol) ;
    CMATwrite(cmat_out, out_vol) ;
    MRIfree(&mri_in) ;
    MRIfree(&mri_out) ;
    exit(0) ;
  }

  fprintf(stderr, "reading volume from %s...\n", in_vol) ;
  mri_in = MRIread(in_vol) ;
  if (!mri_in)
    ErrorExit(ERROR_NOFILE, "%s: could not read MRI volume %s", Progname,
              in_vol) ;
  if (out_like_fname) //E/ maybe need an out_kinda_like_fname
  {
    mri_tmp = MRIread(out_like_fname) ;
    if (!mri_tmp)
      ErrorExit(ERROR_NOFILE, "%s: could not read template volume from %s",Progname,out_like_fname) ;
    mri_out = MRIalloc(mri_tmp->width, mri_tmp->height, mri_tmp->depth, mri_tmp->type) ;
    //E/ maybe better mri_in->type?
    MRIcopyHeader(mri_tmp, mri_out) ; //E/ reinstate this
    //E/ MRIfree(&mri_tmp) ; // keep this around for recopy later.

    //E/ Hey, MRIlinearTransformInterp() just sets
    //dst->ras_good_flag to zero!  and the x/y/zsize and stuff seems
    //to go away during e.g. mghWrite.  recopy later?
  } else  /* assume output should be like input */
  {
    mri_out = MRIclone(mri_in, NULL) ;
#if 0
    mri_out = MRIalloc(256, 256, 256, mri_in->type) ;
    //E/ set xyzc_ras to coronal ones.. - these'll get zorched
    //by MRIlinearTransformInterp() - copy again later - is there
    //any use in having them here now?  yes, so we can pass mri_out
    //to the ras2vox fns.

    //E/ is c_ras = 0,0,0 correct?
    mri_out->x_r =-1;
    mri_out->y_r = 0;
    mri_out->z_r = 0;
    mri_out->c_r =0;
    mri_out->x_a = 0;
    mri_out->y_a = 0;
    mri_out->z_a = 1;
    mri_out->c_a =0;
    mri_out->x_s = 0;
    mri_out->y_s =-1;
    mri_out->z_s = 0;
    mri_out->c_s =0;
    mri_out->ras_good_flag=1;
#endif
  }
  MRIcopyPulseParameters(mri_in, mri_out) ;
  m_total = MatrixIdentity(4, NULL) ;
  nxform_args = nlabels > 0 ? argc : argc-1 ;
  for (i = 2 ; i < nxform_args ; i++) {
    xform_fname = argv[i] ;
    if (strcmp(xform_fname, "-I") == 0) {
      invert_flag = 1 ;
      continue ;
    }
    fprintf(stderr, "reading transform %s...\n", xform_fname) ;
    transform = TransformRead(xform_fname) ;
    if (!transform)
      ErrorExit(ERROR_NOFILE, "%s: could not read transform from %s",
                Progname, xform_fname) ;
    if (out_like_fname == NULL)
      TransformSetMRIVolGeomToDst(transform, mri_out) ;

    if (transform->type != MORPH_3D_TYPE) {
      lta = (LTA *)(transform->xform) ;
      m = MatrixCopy(lta->xforms[0].m_L, NULL) ;
      //E/ mri_rigid_register writes out m as src2trg

      if (lta->type == LINEAR_RAS_TO_RAS || lta->type == LINEAR_CORONAL_RAS_TO_CORONAL_RAS)
        /* convert it to a voxel transform */
      {
        ras_flag = 1 ;
      } else if (ras_flag)
        ErrorExit(ERROR_UNSUPPORTED, "%s: transforms must be all RAS or all voxel",Progname) ;

      if (invert_flag) {
        MATRIX *m_tmp ;
        fprintf(stderr, "inverting transform...\n") ;
        m_tmp = MatrixInverse(m, NULL) ;
        if (!m_tmp)
          ErrorExit(ERROR_BADPARM, "%s: transform is singular!") ;
        MatrixFree(&m) ;
        m = m_tmp ;
        invert_flag = 0 ;
      }
      MatrixMultiply(m, m_total, m_total) ;
      if (ras_flag)  /* convert it to a voxel transform */
      {
        MATRIX *m_tmp ;
        fprintf(stderr, "converting RAS xform to voxel xform...\n") ;
        if (lta->type == LINEAR_RAS_TO_RAS)
          m_tmp = MRIrasXformToVoxelXform(mri_in, mri_out, m_total, NULL) ;
        else if (lta->type == LINEAR_CORONAL_RAS_TO_CORONAL_RAS)
          m_tmp =
            MT_CoronalRasXformToVoxelXform(mri_in, mri_out, m_total, NULL) ;
        else
          //E/ how else could ras_flag be set? prev tx a R2R/CR2CR tx?
          exit(1);

        //////////////////////////////////////////////////////////////
        MatrixPrint(stdout, m_tmp);
        c = VectorAlloc(4, MATRIX_REAL);
        c->rptr[1][1] = (mri_in->width)/2.;
        c->rptr[2][1] = (mri_in->height)/2.;
        c->rptr[3][1] = (mri_in->depth)/2.;
        c->rptr[4][1] = 1.;
        mine = MatrixMultiply(m_tmp, c, NULL);
        MatrixPrint(stdout, mine);
        fprintf(stderr, "voxel pos = %.2f, %.2f, %.2f for %.2f, %.2f, %.2f\n",
                mine->rptr[1][1], mine->rptr[2][1], mine->rptr[3][1],
                c->rptr[1][1], c->rptr[2][1], c->rptr[3][1]);
        VectorFree(&mine);
        /////////////////////////////////////////////////////////////

        MatrixFree(&m_total) ;
        m_total = m_tmp ;
      }
      LTAfree(&lta) ;
      invert_flag = 0 ;
    }
  }
  if (transform->type != MORPH_3D_TYPE) {
    if (nlabels > 0) {
      double dist ;
      int    n ;

      if (subject_name)
        printf("%s  ", subject_name) ;
      for (n = 0 ; n < nlabels ; n++) {
        dist = MRIcomputeLinearTransformLabelDist(mri_in, m_total, labels[n]) ;
        if (quiet_mode)
          printf("%f ", dist) ;
        //     printf("%d %f ", labels[n], dist) ;
        else
          printf("label %d moved on average %2.3f voxels\n", labels[n], dist) ;
      }
      printf("\n") ;
      exit(0) ;
    } else {
      printf("applying voxel transform:\n") ;
      MatrixPrint(stdout, m_total) ;
      MRIlinearTransformInterp(mri_in, mri_out, m_total, resample_type) ;
    }
  } else {
    GCAM *gcam ;
    gcam = (GCAM *)(transform->xform) ;
    if (gcam->type == GCAM_RAS)
    {
#if 0
      printf("!!! warning - no output geometry specified (should use -out_like <fname>) !!!!\n");
#endif
      GCAMrasToVox(gcam, mri_out) ;
    }
    if (invert_flag)
      mri_out = TransformApplyInverseType(transform,mri_in,NULL,resample_type);
    else
      mri_out = TransformApplyType(transform, mri_in, NULL, resample_type) ;
  }

  fprintf(stderr, "writing output to %s.\n", out_vol) ;

  //E/ reinstate what MRIlinearTransform zorched
  if (out_like_fname) {
    //E/ MRIcopyHeader doesn't stick because later
    //MRIlinearTransformInterp sets dst->ras_good_flag to zero!  and
    //the x/y/zsize and stuff seems to go away during e.g. mghWrite.
    //So we recopy it now.  'cause ras xform IS good.  Well, if
    //mri_tmp's is.

    MRIcopyHeader(mri_tmp, mri_out) ;
    MRIfree(&mri_tmp) ;
  } else {
    /* assume out_like the input volume */
#if 0
    //E/ set xyzc_ras to coronal 256/1mm^3 isotropic ones.. - they
    //got zorched by MRIlinearTransformInterp() - so copy again here

    mri_out->x_r =-1;
    mri_out->y_r = 0;
    mri_out->z_r = 0;
    mri_out->c_r =0;
    mri_out->x_a = 0;
    mri_out->y_a = 0;
    mri_out->z_a = 1;
    mri_out->c_a =0;
    mri_out->x_s = 0;
    mri_out->y_s =-1;
    mri_out->z_s = 0;
    mri_out->c_s =0;
    mri_out->ras_good_flag=1;
#endif
  }
  //E/

  MRIwrite(mri_out, out_vol) ;

  exit(0) ;
  return(0) ;  /* for ansi */
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
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "in_like") || !stricmp(option, "in-like") || !stricmp(option, "il")) {
    in_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
  } else if (!stricmp(option, "voxel")){
    cmat_output_coords  = LABEL_COORDS_VOXEL ;
    printf("transforming cmat labels to voxel coords before writing\n") ;
  } else if (!stricmp(option, "scanner")){
    cmat_output_coords  = LABEL_COORDS_SCANNER_RAS ;
    printf("transforming cmat labels to scanner ras coords before writing\n") ;
  } else if (!stricmp(option, "tkreg")){
    cmat_output_coords  = LABEL_COORDS_TKREG_RAS ;
    printf("transforming cmat labels to tkreg ras coords before writing\n") ;
  } else if (!stricmp(option, "out_like") || !stricmp(option, "out-like") || !stricmp(option, "ol")) {
    out_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
  } else switch (toupper(*option)) {
    case 'Q':
      quiet_mode = 1 ;
      break ;
    case 'S':
      subject_name = argv[2] ;
      nargs = 1  ;
      break ;
    case 'D':
      if (nlabels >= 1000)
        ErrorExit(ERROR_NOMEMORY, "%s: too many labels (%d)", Progname, nlabels) ;
      labels[nlabels] = atoi(argv[2]) ;
      nargs = 1 ;
      fprintf(stderr, "computing average distance traversed by label %d\n", labels[nlabels]) ;
      nlabels++ ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'I':
      invert_flag = 1 ;
      break ;
    case 'R':
      if (strcmp(StrLower(argv[2]), "interpolate") == 0)
        resample_type = SAMPLE_TRILINEAR;
      else if (strcmp(StrLower(argv[2]), "nearest") == 0)
        resample_type = SAMPLE_NEAREST;
      else if (strcmp(StrLower(argv[2]), "weighted") == 0)
        resample_type = SAMPLE_WEIGHTED;
      else if (strcmp(StrLower(argv[2]), "sinc") == 0)
        resample_type = SAMPLE_SINC;
      else if (strcmp(StrLower(argv[2]), "cubic") == 0)
        resample_type = SAMPLE_CUBIC;
      else {
        fprintf(stderr, "\n%s: unknown resample type \"%s\"\n", Progname, argv[2]);
        usage_exit();
        exit(1);
      }
      printf("setting resample type to %d\n", resample_type) ;
      nargs = 1 ;
      break ;
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
  //  print_usage() ; // print_help _calls print_usage
  print_help() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          //          "usage: %s [options] <input volume> <input surface> <registration file> <output .float file>\n",
          //E/ where did that come from??

          "usage: %s [options] <input volume> <lta file> <output file>\n",

          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
#if 0
  fprintf(stderr,
          "\nThis program will paint a average Talairach stats onto a surface\n");
  fprintf(stderr, "-imageoffset <image offset> - set offset to use\n") ;
  fprintf(stderr, "-S                          - paint using surface "
          "coordinates\n") ;
#else
  fprintf(stderr,
          "\nThis program will apply a linear transform to mri volume and write out the result.  The output volume is by default 256^3 1mm^3 isotropic, or you can specify an -out_like volume.  I think there's a bug in -i behavior if you're specifying multiple transforms.\n");
  fprintf(stderr, "-out_like <reference volume> - set out_volume parameters\n") ;
  fprintf(stderr, "-I                           - invert transform "
          "coordinates\n") ;


#endif



  exit(1) ;
}


static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


MATRIX *
MT_CoronalRasXformToVoxelXform
(MRI *mri_in, MRI *mri_out, MATRIX *m_corras_s2corras_t, MATRIX *m_vox_s2vox_t) {
  MATRIX *D_src, *D_trg, *D_trg_inv, *m_tmp;

  //E/ D's are the vox2coronalras matrices generated by
  //MRIxfmCRS2XYZtkreg().  Thanks, Doug.

  D_src = MRIxfmCRS2XYZtkreg(mri_in);
  D_trg = MRIxfmCRS2XYZtkreg(mri_out);
  D_trg_inv = MatrixInverse(D_trg, NULL);

  //E/ m_vox_s2vox_t = D_trg_inv * m_corras_s2corras_t * D_src

  m_tmp = MatrixMultiply(m_corras_s2corras_t, D_src, NULL);
  if (m_vox_s2vox_t == NULL)
    m_vox_s2vox_t = MatrixMultiply(D_trg_inv, m_tmp, NULL);
  else
    MatrixMultiply(D_trg_inv, m_tmp, m_vox_s2vox_t);

//#define _MTCRX2RX_DEBUG
#ifdef _MTCRX2RX_DEBUG
  fprintf(stderr, "m_corras_s2corras_t = \n");
  MatrixPrint(stderr, m_corras_s2corras_t);

  fprintf(stderr, "D_trg_inv = \n");
  MatrixPrint(stderr, D_trg_inv);
  fprintf(stderr, "D_src = \n");
  MatrixPrint(stderr, D_src);

  fprintf(stderr, "m_vox_s2vox_t = D_trg_inv * m_corras_s2corras_t * D_src = \n");
  MatrixPrint(stderr, m_vox_s2vox_t);
#endif

  MatrixFree(&m_tmp);
  MatrixFree(&D_src);
  MatrixFree(&D_trg);
  MatrixFree(&D_trg_inv);

  return(m_vox_s2vox_t);
}
double
MRIcomputeLinearTransformLabelDist(MRI *mri_src, MATRIX *mA, int label) {
  int    x1, x2, x3, width, height, depth, nvox ;
  VECTOR *v_X, *v_Y ;   /* original and transformed coordinate systems */
  Real   y1, y2, y3, total_dist, dist ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  v_X = VectorAlloc(4, MATRIX_REAL) ;  /* input (src) coordinates */
  v_Y = VectorAlloc(4, MATRIX_REAL) ;  /* transformed (dst) coordinates */

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("MRIlinearTransformInterp: Applying transform\n");

  total_dist = 0.0 ;
  v_Y->rptr[4][1] = 1.0f ;
  v_X->rptr[4][1] = 1.0f ;
  nvox = 0 ;
  for (x3 = 0 ; x3 < depth ; x3++) {
    V3_Z(v_X) = x3 ;
    for (x2 = 0 ; x2 < height ; x2++) {
      V3_Y(v_X) = x2 ;
      for (x1 = 0 ; x1 < width ; x1++) {
        if (MRIvox(mri_src, x1, x2, x3) != label)
          continue ;
        V3_X(v_X) = x1 ;
        MatrixMultiply(mA, v_X, v_Y) ;

        y1 = V3_X(v_Y) ;
        y2 = V3_Y(v_Y) ;
        y3 = V3_Z(v_Y) ;

        dist = sqrt(SQR(x1-y1)+SQR(x2-y2)+SQR(x3-y3)) ;
        total_dist += dist ;
        nvox++ ;
      }
    }
  }

  MatrixFree(&v_X) ;
  MatrixFree(&v_Y) ;

  if (nvox == 0)
    nvox = 1 ;

  return(total_dist / (float)nvox) ;
}
