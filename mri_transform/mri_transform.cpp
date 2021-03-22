/**
 * @brief program to transform an mri or cmat structure using a linear transform
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
#ifdef HAVE_OPENMP // mrisurf.c has numerous parallelized functions
#include "romp_support.h"
#endif

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
#include "mrinorm.h"
#include "fsinit.h"

double MRIcomputeLinearTransformLabelDist(MRI *mri_src, MATRIX *mA, int label) ;

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

const char *Progname ;
static int quiet_mode = 0 ;
static char *subject_name ;
static char *out_like_fname = NULL ;
static char *in_like_fname = NULL ;  // for cmat stuff
static int invert_flag = 0 ;
static int resample_type = SAMPLE_TRILINEAR ;
static int nlabels = 0 ;
static int labels[1000] ;

static int cmat_output_coords = LABEL_COORDS_VOXEL ;

static char *in_surf_name = NULL ;
static char *out_surf_name = NULL ;
MRI *MRIScreateVolumeWarpFromSurface(MRI *mri_in, MRI *mri_out, MRI_SURFACE *mris, int which_in, int which_out, double res_scale) ;
static MRI *MRIScreateVolumeWarpFromSphere(MRI *mri_in, MRI *mri_out, MRI_SURFACE *mris, int which_in, int which_out, double res_scale) ;

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

  nargs = handleVersionOption(argc, argv, "mri_transform");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  FSinit() ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (in_surf_name)
  {
    if ((argc < 3) && (in_surf_name != NULL))
      usage_exit() ;
  }
  else
  {
    if (((argc < 4) && (nlabels == 0)) || ((argc < 3) && (nlabels > 0)))
      usage_exit() ;
  }

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
  if (in_surf_name)  // create transform based on surface correspondences
  {
    MRI_SURFACE *mris ;
    MRI         *mri_warp ;
    char        fname[STRLEN] ;

    mris = MRISread(in_surf_name) ;
    if (mris == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read input surface from %s", Progname, in_surf_name) ;
    if (MRISreadCanonicalCoordinates(mris, out_surf_name) != NO_ERROR)
      ErrorExit(ERROR_NOFILE, "%s: could not read output surface from %s", Progname, out_surf_name) ;

    MRISrecenter(mris, CURRENT_VERTICES, CANONICAL_VERTICES) ;
    sprintf(fname, "./%s.sphere.centered", mris->hemisphere == LEFT_HEMISPHERE ? "lh" : "rh") ;
    MRISwriteVertexLocations(mris, fname, CANONICAL_VERTICES) ;
    
    mri_warp = MRIScreateVolumeWarpFromSphere(mri_in, mri_out, mris, CURRENT_VERTICES, CANONICAL_VERTICES,1) ;
    MRIapplyMorph(mri_in, mri_warp, mri_out, resample_type) ;
    printf("writing output volume to %s\n", out_vol) ;
    MRIwrite(mri_out, out_vol) ;
    exit(0) ;
  }

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
  else if (!stricmp(option, "debug_voxel") || !stricmp(option, "debug-voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    nargs = 3 ;
    printf("debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
  } else if (!stricmp(option, "in_like") || !stricmp(option, "in-like") || !stricmp(option, "il")) {
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
  } else if (!stricmp(option, "surf") || !stricmp(option, "surface")) {
    in_surf_name = argv[2] ;
    out_surf_name = argv[3] ;
    nargs = 2 ;
    printf("creating volume mapping from %s-->%s surface map\n", in_surf_name, out_surf_name) ;
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
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
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
  double y1, y2, y3, total_dist, dist ;

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
static int niter = 0000;
static double min_change= 0.001 ;

MRI *
MRIScreateVolumeWarpFromSurface(MRI *mri_in, MRI *mri_out, MRI_SURFACE *mris, int which_in, int which_out, double res_scale) 
{
  MRI    *mri_warp, *mri_ctrl, *mri_warp_x, *mri_warp_y, *mri_warp_z, *mri_warp_x_smooth, *mri_warp_y_smooth, *mri_warp_z_smooth ;
  int    vno ;
  VERTEX *v ;
  double  xi, yi, zi, xo, yo, zo, xiv, yiv, ziv, xov, yov, zov ;
  int    xv, yv, zv ;


  mri_ctrl = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_UCHAR, 1) ;
  mri_warp_x = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 1) ;
  mri_warp_y = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 1) ;
  mri_warp_z = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 1) ;
  mri_warp = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 3) ;
  MRIcopyHeader(mri_in, mri_warp) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    MRISvertexCoord2XYZ_double(v, which_in, &xi, &yi, &zi) ;
    MRISvertexCoord2XYZ_double(v, which_out, &xo, &yo, &zo) ;
    MRISsurfaceRASToVoxelCached(mris, mri_in, xi, yi, zi, &xiv, &yiv, &ziv) ;
    MRISsurfaceRASToVoxelCached(mris, mri_out, xo, yo, zo, &xov, &yov, &zov) ;
    xv = nint(xov) ; yv = nint(yov) ; zv = nint(zov) ;
    if (xv < 0 || yv < 0 || zv < 0 || xv >= mri_warp->width || yv >= mri_warp->height || zv >= mri_warp->depth)
      continue ;
    if (vno == Gdiag_no)
    {
      printf("v %d: (%d, %d, %d) <- (%2.1f, %2.1f, %2.1f) (Dx=%2.1f, %2.1f, %2.1f)\n", 
	     vno, xv, yv, zv, xiv, yiv, ziv, xov-xiv, yov-yiv, zov-ziv);
      if (Gx < 0)
      {
	Gx = xv ; Gy = yv ; Gz = zv ;
      }
    }
    MRIsetVoxVal(mri_warp_x, xv, yv, zv, 0, xiv-xov);
    MRIsetVoxVal(mri_warp_y, xv, yv, zv, 0, yiv-yov);
    MRIsetVoxVal(mri_warp_z, xv, yv, zv, 0, ziv-zov);
    MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, CONTROL_MARKED) ;
  }

  MRIbuildVoronoiDiagram(mri_warp_x, mri_ctrl, mri_warp_x);
  MRIbuildVoronoiDiagram(mri_warp_y, mri_ctrl, mri_warp_y);
  MRIbuildVoronoiDiagram(mri_warp_z, mri_ctrl, mri_warp_z);
  mri_warp_x_smooth = MRIsoapBubble(mri_warp_x, mri_ctrl, NULL, niter, min_change);
  mri_warp_y_smooth = MRIsoapBubble(mri_warp_y, mri_ctrl, NULL, niter, min_change);
  mri_warp_z_smooth = MRIsoapBubble(mri_warp_z, mri_ctrl, NULL, niter, min_change);
  MRIcopyFrame(mri_warp_x_smooth, mri_warp, 0, 0) ;
  MRIcopyFrame(mri_warp_y_smooth, mri_warp, 0, 1) ;
  MRIcopyFrame(mri_warp_z_smooth, mri_warp, 0, 2) ;

  MRIfree(&mri_ctrl) ; 
  MRIfree(&mri_warp_x) ; MRIfree(&mri_warp_y) ; MRIfree(&mri_warp_z) ;
  MRIfree(&mri_warp_x_smooth) ; MRIfree(&mri_warp_y_smooth) ; MRIfree(&mri_warp_z_smooth) ;

  return(mri_warp) ;
}

#define MAX_SAMPLE_DIST 13.0

static MRI *
MRIScreateVolumeWarpFromSphere(MRI *mri_in, MRI *mri_out, MRI_SURFACE *mris, int which_in, int which_out, double res_scale) 
{
  MRI    *mri_warp, *mri_ctrl, *mri_warp_x, *mri_warp_y, *mri_warp_z, *mri_warp_x_smooth, *mri_warp_y_smooth, *mri_warp_z_smooth ;
  MRI    *mri_interior, *mri_x = NULL, *mri_y = NULL, *mri_z = NULL, *mri_c = NULL, *mri_dilated = NULL ;
  int    vno ;
  VERTEX *v;
  double  xi, yi, zi, xo, yo, zo, xiv, yiv, ziv, xov, yov, zov, r, cnx, cny, cnz, sample_dist = 0.25, radius ;
  double  cx0, cy0, cz0;
  int     xv, yv, zv, inside, xvi, yvi, zvi, was_inside, was_outside, nchanged, iter, nbrs ;

  mri_interior = MRIclone(mri_in, NULL) ;
  MRISfillInterior(mris, mri_in->xsize, mri_interior) ;
  mri_ctrl = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_UCHAR, 1) ;
  mri_warp_x = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 1) ;
  mri_warp_y = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 1) ;
  mri_warp_z = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 1) ;
  mri_warp = MRIallocSequence(mri_in->width, mri_in->height, mri_in->depth, MRI_FLOAT, 3) ;
  MRIcopyHeader(mri_in, mri_warp) ;

  MRIsetValues(mri_warp_x, -1000) ; MRIsetValues(mri_warp_y, -1000) ; MRIsetValues(mri_warp_z, -1000) ;

  v = &mris->vertices[0] ;
  radius = sqrt(SQR(v->cx-mris->x0) + SQR(v->cy-mris->y0) + SQR(v->cz-mris->z0)) ;

  // traverse each neighbor vector and search outwards along the interpolated normal 
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    int     n ;
    double  dist, clen, len, cnx0, cny0, cnz0, cnx1, cny1, cnz1, dx, dy, dz, dcx, dcy, dcz, xi0, yi0, zi0 ;
    double  nx, ny, nz ;

    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];

    if (v->ripflag)
      continue ;
    if (vno == Gdiag_no)
      DiagBreak() ;

    cnx0 = v->cx-mris->x0 ; cny0 = v->cy-mris->y0 ; cnz0 = v->cz-mris->z0 ;
    r = sqrt(cnx0*cnx0 + cny0*cny0 + cnz0*cnz0) ;
    cnx0 /= r ; cny0 /= r ; cnz0 /= r ;

    for (n = 0 ; n < vt->vnum ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->ripflag)
	continue ;
      if (vt->v[n] == Gdiag_no)
	DiagBreak() ;

      dx =   vn->x - v->x ;   dy = vn->y - v->y ;    dz = vn->z - v->z ;
      dcx =  vn->cx - v->cx ; dcy = vn->cy - v->cy ; dcz = vn->cz - v->cz ;
      len = sqrt(dx*dx + dy*dy + dz*dz) ;
      clen = sqrt(dcx*dcx + dcy*dcy + dcz*dcz) ;
      cnx1 = vn->cx-mris->x0 ; cny1 = vn->cy-mris->y0 ; cnz1 = vn->cz-mris->z0 ;
      r = sqrt(cnx1*cnx1 + cny1*cny1 + cnz1*cnz1) ;
      cnx1 /= r ; cny1 /= r ; cnz1 /= r ;

      for (dist = 0 ; dist <= 1.0 ; dist += sample_dist)
      {
	// linearly interpolate normals
	cnx = dist * cnx1 + (1.0-dist)*cnx0 ;
	cny = dist * cny1 + (1.0-dist)*cny0 ;
	cnz = dist * cnz1 + (1.0-dist)*cnz0 ;
	nx = dist * vn->nx + (1.0-dist)*v->nx ;
	ny = dist * vn->ny + (1.0-dist)*v->ny ;
	nz = dist * vn->nz + (1.0-dist)*v->nz ;
	
	// compute point along edge on white and sphere
	cx0 = v->cx + dcx*dist ; 
	cy0 = v->cy + dcy*dist ; 
	cz0 = v->cz + dcz*dist ;
	xi0 = v->x+dist*dx ; 
	yi0 = v->y+dist*dz ; 
	zi0 = v->z+dist*dz ;

	// add points moving inwards
	MRISsurfaceRASToVoxelCached(mris, mri_in, xi0, yi0, zi0, &xiv, &yiv, &ziv) ;
	xvi = nint(xiv) ; yvi = nint(yiv) ; zvi = nint(ziv) ;
	inside = MRIgetVoxVal(mri_interior, xvi, yvi, zvi, 0) ;
	was_inside = inside ;
	for (r = -sample_dist ; r >= -MAX_SAMPLE_DIST ; r -= sample_dist)
	{
	  xi = xi0+r*nx ; yi = yi0+r*ny ; zi = zi0+r*nz ;
	  xo = cx0+r*cnx ;  yo = cy0+r*cny ;  zo = cz0+r*cnz ;
	  MRISsurfaceRASToVoxelCached(mris, mri_in, xi, yi, zi, &xiv, &yiv, &ziv) ;
	  MRISsurfaceRASToVoxelCached(mris, mri_out, xo, yo, zo, &xov, &yov, &zov) ;
	  xv = nint(xov) ; yv = nint(yov) ; zv = nint(zov) ;
	  if (xv == Gx && yv == Gy && zv == Gz)
	    DiagBreak() ;
	  if (xv < 0 || yv < 0 || zv < 0 || xv >= mri_warp->width || yv >= mri_warp->height || zv >= mri_warp->depth)
	    continue ;
	  xvi = nint(xiv) ; yvi = nint(yiv) ; zvi = nint(ziv) ;
	  inside = MRIgetVoxVal(mri_interior, xvi, yvi, zvi, 0) ;
	  
	  if (was_inside && inside == 0)  // was in the wm but have now left
	    break ;
	  if (was_inside == 0 && inside)
	    was_inside = 1 ; // entered interior of wm surface
	  
	  if (nint(MRIgetVoxVal(mri_ctrl, xv, yv, zv, 0)) == CONTROL_MARKED)
	    continue ;
	  if (vno == Gdiag_no)
	  {
	    printf("v %d: n %2.1f, (%d, %d, %d) <- (%2.1f, %2.1f, %2.1f) (Dx=%2.1f, %2.1f, %2.1f)\n", 
		   vno, r, xv, yv, zv, xiv, yiv, ziv, xov-xiv, yov-yiv, zov-ziv);
	    DiagBreak() ;
	  }
	  MRIsetVoxVal(mri_warp_x, xv, yv, zv, 0, xiv-xov);
	  MRIsetVoxVal(mri_warp_y, xv, yv, zv, 0, yiv-yov);
	  MRIsetVoxVal(mri_warp_z, xv, yv, zv, 0, ziv-zov);
	  MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, CONTROL_MARKED) ;
	}

	// add points moving outwards
	MRISsurfaceRASToVoxelCached(mris, mri_in, xi0, yi0, zi0, &xiv, &yiv, &ziv) ;
	xvi = nint(xiv) ; yvi = nint(yiv) ; zvi = nint(ziv) ;
	inside = MRIgetVoxVal(mri_interior, xvi, yvi, zvi, 0) ;
	was_inside = inside ;
	was_outside = !inside ;
	for (r = sample_dist ; r <= MAX_SAMPLE_DIST ; r += sample_dist)
	{
	  xi = xi0+r*nx ; yi = yi0+r*ny ; zi = zi0+r*nz ;
	  xo = cx0+r*cnx ;  yo = cy0+r*cny ;  zo = cz0+r*cnz ;
	  MRISsurfaceRASToVoxelCached(mris, mri_in, xi, yi, zi, &xiv, &yiv, &ziv) ;
	  MRISsurfaceRASToVoxelCached(mris, mri_out, xo, yo, zo, &xov, &yov, &zov) ;
	  xv = nint(xov) ; yv = nint(yov) ; zv = nint(zov) ;
	  if (xv == Gx && yv == Gy && zv == Gz)
	    DiagBreak() ;
	  if (xv < 0 || yv < 0 || zv < 0 || xv >= mri_warp->width || yv >= mri_warp->height || zv >= mri_warp->depth)
	    continue ;
	  xvi = nint(xiv) ; yvi = nint(yiv) ; zvi = nint(ziv) ;
	  inside = MRIgetVoxVal(mri_interior, xvi, yvi, zvi, 0) ;
	  if (was_outside && inside > 0) // found exterior of surface and went back in
	    break ;
	  if (was_inside && inside == 0)
	    was_outside = 1 ;            // found exterior of surface
	  
	  if (nint(MRIgetVoxVal(mri_ctrl, xv, yv, zv, 0)) == CONTROL_MARKED)
	    continue ;
	  if (vno == Gdiag_no)
	  {
	    printf("v %d: n %2.1f, (%d, %d, %d) <- (%2.1f, %2.1f, %2.1f) (Dx=%2.1f, %2.1f, %2.1f)\n", 
		   vno, r, xv, yv, zv, xiv, yiv, ziv, xov-xiv, yov-yiv, zov-ziv);
	    DiagBreak() ;
	  }
	  MRIsetVoxVal(mri_warp_x, xv, yv, zv, 0, xiv-xov);
	  MRIsetVoxVal(mri_warp_y, xv, yv, zv, 0, yiv-yov);
	  MRIsetVoxVal(mri_warp_z, xv, yv, zv, 0, ziv-zov);
	  MRIsetVoxVal(mri_ctrl, xv, yv, zv, 0, CONTROL_MARKED) ;
	}
      }
    }
  }

#if 0
  mri_dilated = MRIdilate(mri_ctrl, NULL) ;
  iter = 0 ;
  do
  {
    nchanged = 0 ;
    mri_x = MRIcopy(mri_warp_x, mri_x) ;
    mri_y = MRIcopy(mri_warp_y, mri_y) ;
    mri_z = MRIcopy(mri_warp_z, mri_z) ;
    mri_c = MRIcopy(mri_ctrl, mri_c) ;

    for (xv = 0 ; xv < mri_ctrl->width ; xv++)
      for (yv = 0 ; yv < mri_ctrl->height ; yv++)
	for (zv = 0 ; zv < mri_ctrl->depth ; zv++)
	{
	  double mean ;
	  if (xv == Gx && yv == Gy && zv == Gz)
	    DiagBreak() ;
	  if (MRIgetVoxVal(mri_ctrl, xv, yv, zv, 0) == CONTROL_MARKED)
	    continue ;
	  if (MRIgetVoxVal(mri_dilated, xv, yv, zv, 0) == 0) // doesn't border a control point
	    continue ;
	  if (MRIneighbors3x3(mri_ctrl, xv, yv, zv, CONTROL_MARKED) < floor(3*3*3.0*(2.0/4.0)))
	    continue ;

	  nchanged++ ;
	  mean = MRImeanInLabelInRegion(mri_warp_x, mri_ctrl, CONTROL_MARKED, xv, yv, zv, 1, NULL);
	  MRIsetVoxVal(mri_x, xv, yv, zv, 0, mean) ;
	  mean = MRImeanInLabelInRegion(mri_warp_y, mri_ctrl, CONTROL_MARKED, xv, yv, zv, 1, NULL);
	  MRIsetVoxVal(mri_y, xv, yv, zv, 0, mean) ;
	  mean = MRImeanInLabelInRegion(mri_warp_z, mri_ctrl, CONTROL_MARKED, xv, yv, zv, 1, NULL);
	  MRIsetVoxVal(mri_z, xv, yv, zv, 0, mean) ;
	  MRIsetVoxVal(mri_c, xv, yv, zv, 0, CONTROL_MARKED) ;
	}
    MRIcopy(mri_x, mri_warp_x) ;
    MRIcopy(mri_y, mri_warp_y) ;
    MRIcopy(mri_z, mri_warp_z) ;
    MRIcopy(mri_c, mri_ctrl) ;
    printf("iter %3.3d: nchanged = %d\n", iter, nchanged) ;
    if (nchanged == 0 || iter++ > 4)
    {
      MRIfree(&mri_x) ; MRIfree(&mri_y) ; MRIfree(&mri_z) ; MRIfree(&mri_c) ; MRIfree(&mri_dilated) ;
      nchanged = 0 ;  // force loop to terminate
    }
  } while (nchanged > 0) ;
#endif


  MRISsurfaceRASToVoxelCached(mris, mri_in, mris->x0, mris->y0, mris->z0, &cx0, &cy0, &cz0) ;
  for (nbrs = 3*3*3-1 ; nbrs > 0 ; nbrs--)
  {
    mri_dilated = MRIdilate(mri_ctrl, mri_dilated) ;
    iter = 0 ;
    do
    {
      nchanged = 0 ;
      mri_x = MRIcopy(mri_warp_x, mri_x) ;
      mri_y = MRIcopy(mri_warp_y, mri_y) ;
      mri_z = MRIcopy(mri_warp_z, mri_z) ;
      mri_c = MRIcopy(mri_ctrl, mri_c) ;
      
      for (xv = 0 ; xv < mri_ctrl->width ; xv++)
	for (yv = 0 ; yv < mri_ctrl->height ; yv++)
	  for (zv = 0 ; zv < mri_ctrl->depth ; zv++)
	  {
	    double mean, r ;

	    r = sqrt(SQR(xv-cx0) + SQR(yv-cy0) + SQR(zv-cz0)) ;

	    if (fabs(r-radius)>MAX_SAMPLE_DIST)
	      continue ;
	    if (xv == Gx && yv == Gy && zv == Gz)
	      DiagBreak() ;
	    if (MRIgetVoxVal(mri_ctrl, xv, yv, zv, 0) == CONTROL_MARKED)
	      continue ;
	    if (MRIgetVoxVal(mri_dilated, xv, yv, zv, 0) == 0) // doesn't border a control point
	      continue ;
	    if (MRIneighbors3x3(mri_ctrl, xv, yv, zv, CONTROL_MARKED) < nbrs)
	      continue ;
	    
	    nchanged++ ;
	    mean = MRImeanInLabelInRegion(mri_warp_x, mri_ctrl, CONTROL_MARKED, xv, yv, zv, 1, NULL);
	    MRIsetVoxVal(mri_x, xv, yv, zv, 0, mean) ;
	    mean = MRImeanInLabelInRegion(mri_warp_y, mri_ctrl, CONTROL_MARKED, xv, yv, zv, 1, NULL);
	    MRIsetVoxVal(mri_y, xv, yv, zv, 0, mean) ;
	    mean = MRImeanInLabelInRegion(mri_warp_z, mri_ctrl, CONTROL_MARKED, xv, yv, zv, 1, NULL);
	    MRIsetVoxVal(mri_z, xv, yv, zv, 0, mean) ;
	    MRIsetVoxVal(mri_c, xv, yv, zv, 0, CONTROL_MARKED) ;
	  }
      MRIcopy(mri_x, mri_warp_x) ;
      MRIcopy(mri_y, mri_warp_y) ;
      MRIcopy(mri_z, mri_warp_z) ;
      MRIcopy(mri_c, mri_ctrl) ;
      printf("nbrs %2.2d, iter %3.3d: nchanged = %d\n", nbrs, iter++, nchanged) ;
      if (nchanged == 0)
      {
	MRIfree(&mri_x) ; MRIfree(&mri_y) ; MRIfree(&mri_z) ; MRIfree(&mri_c) ; MRIfree(&mri_dilated) ;
	nchanged = 0 ;  // force loop to terminate
      }
    } while (nchanged > 0) ;
  }
#if 0
  MRIbuildVoronoiDiagram(mri_warp_x, mri_ctrl, mri_warp_x);
  MRIbuildVoronoiDiagram(mri_warp_y, mri_ctrl, mri_warp_y);
  MRIbuildVoronoiDiagram(mri_warp_z, mri_ctrl, mri_warp_z);
#endif
  if (niter > 0)
  {
    mri_warp_x_smooth = MRIsoapBubble(mri_warp_x, mri_ctrl, NULL, niter, min_change);
    mri_warp_y_smooth = MRIsoapBubble(mri_warp_y, mri_ctrl, NULL, niter, min_change);
    mri_warp_z_smooth = MRIsoapBubble(mri_warp_z, mri_ctrl, NULL, niter, min_change);
  }
  else
  {
    mri_warp_x_smooth = MRIcopy(mri_warp_x, NULL) ;
    mri_warp_y_smooth = MRIcopy(mri_warp_y, NULL) ;
    mri_warp_z_smooth = MRIcopy(mri_warp_z, NULL) ;
  }

  MRIcopyFrame(mri_warp_x_smooth, mri_warp, 0, 0) ;
  MRIcopyFrame(mri_warp_y_smooth, mri_warp, 0, 1) ;
  MRIcopyFrame(mri_warp_z_smooth, mri_warp, 0, 2) ;
  MRIfree(&mri_ctrl) ; MRIfree(&mri_interior) ;
  MRIfree(&mri_warp_x) ; MRIfree(&mri_warp_y) ; MRIfree(&mri_warp_z) ;
  MRIfree(&mri_warp_x_smooth) ; MRIfree(&mri_warp_y_smooth) ; MRIfree(&mri_warp_z_smooth) ;

  return(mri_warp) ;
}
