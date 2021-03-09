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
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "numerics.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "proto.h"
#include "transform.h"
#include "version.h"
#include "matrix.h"
#include "density.h"
#include "mrisegment.h"
#include "mri_circulars.h"

#define RGB_SIZE 500


static int    powell_minimize(MRI *mri_block,
                              MRI *mri_histo,
                              MRI *mri_seg,
                              DENSITY *pdf,
                              MATRIX *mat,
                              int cost_type) ;
static double compute_overlap(MRI *mri_src,
                              MRI *mri_dst,
                              MATRIX *m_total) ;

static double min_overlap = 0.8 ;

static int align = 1 ;
static float min_angle = -RADIANS(30) ;
static float max_angle = RADIANS(30) ;
static int nangles = 7 ;

static int snapshots = 0 ;
static char base[STRLEN] = "nissl" ;
static float min_scale = .5 ;
static float max_scale = 1.5;
static int nscales = 7 ;

static int probe = 0 ;
static int max_trans = 150 ;
static int ntrans = 7 ;
static int skip = 0 ;

static int nlevels = 1 ;

static double rotate_radians = 0 ;

static int probe_cost_function(MRI *mri_histo, MRI *mri_block, MRI *mri_seg, DENSITY *pdf, MATRIX *mat, int cost_type, char *out_name) ;
int MRIeraseImageBorder(MRI *mri, int width) ;
static MRI *rotate_image(MRI *mri_src, MRI *mri_dst, double rotate_radians);
static int compute_optimal_translation(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, MATRIX *m_best,
                                       DENSITY *pdf, int level,
                                       int max_trans, int nsteps, int skip, MRI *mri_fullres_src,
                                       int cost_type) ;
static int align_pyramid_level(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, MATRIX *xform, DENSITY *pdf,
                               int level, MRI *mri_fullres_src, int cost_type, int skip) ;
static MRI *mri_apply_slice_xform(MRI *mri_src, MRI *mri_dst, MATRIX *m, int slice) ;
static int write_snapshot(MRI *mri, MRI *mri_src, MRI *mri_dst, MATRIX *m, char *base, int n, int level, DENSITY *density, MRI *mri_seg) ;
int main(int argc, char *argv[]) ;
static int mriWriteImageView(MRI *mri, char *base_name, int target_size, int view,
                             int slice, MRI *mri_template) ;
static double compute_ml_alignment_error(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY *pdf,
    MATRIX *m_total, int skip, int level) ;
static double compute_cr_alignment_error(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY *pdf,
    MATRIX *m_total, int skip, int level) ;
static double compute_alignment_error(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY *pdf, MATRIX *m_total,
                                      int skip, int level, int cost_type) ;
static MRI *align_block_to_histo(MRI *mri_histo, MRI *mri_block, MRI *mri_seg, DENSITY **pdfs, int nlevels, MATRIX *mat,
                                 int cost_type) ;
static int  compute_optimal_xform(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY **pdfs, int nlevels, MATRIX *xform,
                                  int cost_type, int skip) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;
static float slice_thickness = 1.0 ;
static float in_plane = 1.0 ;

#define MAX_SLICES 10000

#define COST_MAXIMUM_LIKELIHOOD  0
#define COST_KULLABACK_LEIBLER   1
#define COST_CORRELATION_RATIO   2

static int cost_type = COST_MAXIMUM_LIKELIHOOD ;

#define MAX_LEVELS 30

int
main(int argc, char *argv[]) {
  char        **av, *out_name ;
  int         ac, nargs, max_width, max_height, n, i, min_val, max_val, max_vox, level ;
  MRI         *mri_block, *mri_histo, *mri_tmp, *mri_template, *mri_aligned, *mri_seg ;
  DENSITY     *pdfs[MAX_LEVELS], *pdf ;
  MATRIX      *mat ;
  MRI_SEGMENTATION *mriseg ;


  nargs = handleVersionOption(argc, argv, "histo_register_block");
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

  if (argc < 4)
    usage_exit() ;

  mri_block = MRIread(argv[1]) ;
  if (mri_block == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not open block face image %s...\n", argv[1]) ;

  mri_histo = MRIread(argv[2]) ;
  if (mri_histo == NULL)
    ErrorExit(ERROR_BADPARM, "%s: could not open histology image %s...\n", argv[2]) ;

  if (rotate_radians > 0) {
    MRI *mri_tmp ;

    printf("rotating block face image by %2.3f deg....\n", DEGREES(rotate_radians)) ;
    mri_tmp = rotate_image(mri_block, NULL, rotate_radians) ;
    MRIfree(&mri_block) ;
    mri_block = mri_tmp ;
  }

  for (level = 0 ; level < nlevels ; level++) {
    char fname[STRLEN] ;

    sprintf(fname, "%s_level%d.den", argv[3], level) ;
    pdfs[level] = DensityRead(fname) ;
    if (pdfs[level] == NULL)
      ErrorExit(ERROR_BADPARM, "%s: could not open joint density file %s...\n", fname) ;
  }

  pdf = pdfs[0] ;
  out_name = argv[argc-1] ;
  FileNameRemoveExtension(out_name, base) ;
  max_width = MAX(mri_block->width, mri_histo->width) ;
  max_height = MAX(mri_block->height, mri_histo->height) ;
  max_vox = max_width*max_height/4 ;

  /* max_trans = MAX(2*max_width/3, 2*max_height/3) ;*/

  /* find range to do segmentation over */
  n = ceil(pdf->max_val1 - pdf->min_val1 + 1) ;
  min_val = n ;
  max_val = 0 ;
  for (i = 1 ; i < n ; i++)
    if (pdf->valid1[i]) {
      if (i < min_val)
        min_val = i ;
      if (i > max_val)
        max_val = i ;
    }
  if (max_val == 255)
    max_val = 199 ;

  mriseg = MRImaxsegment(mri_histo, min_val, max_val) ;
  mri_seg = MRIsegmentToImage(mri_histo, NULL, mriseg, 0) ;
  MRIbinarize(mri_seg, mri_seg, 0, 0, 1) ;
#if 0
  MRIdilate(mri_seg, mri_seg) ;
  MRIdilate(mri_seg, mri_seg) ;
  MRIeraseImageBorder(mri_seg, 10) ;
  while ((i = MRIvoxelsInLabel(mri_seg, 1)) > max_vox) {
    printf("%d voxels in segmented image (%d max)\n", i, max_vox) ;
    MRIerode(mri_seg, mri_seg) ;
  }
#endif

  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    MRIwrite(mri_seg, "s.mgz") ;

#if 1
  MRIreplaceValues(mri_histo, mri_histo, 0, 1) ;  /* reserve 0 for background */
  MRIreplaceValues(mri_block, mri_histo, 0, 1) ;  /* reserve 0 for background */
  printf("embedding images in (%d x %d) matrix\n", max_width, max_height) ;
  mri_template = MRIalloc(max_width, max_height, 1, MRI_FLOAT) ;
  MRIcopyHeader(mri_histo, mri_template) ;
  mri_tmp = MRIresampleFill(mri_histo, mri_template, SAMPLE_NEAREST, 255) ;
  MRIfree(&mri_histo) ;
  mri_histo = mri_tmp ;
  mri_tmp = MRIresampleFill(mri_block, mri_template, SAMPLE_NEAREST, 255) ;
  MRIfree(&mri_block) ;
  mri_block = mri_tmp ;
  MRIfree(&mri_template) ;
#endif
  mat = MatrixIdentity(3, NULL) ;
  if (probe) {
    probe_cost_function(mri_histo, mri_block, mri_seg, pdf, mat, cost_type, out_name) ;
    exit(0) ;
  } else {
    mri_aligned = align_block_to_histo(mri_histo, mri_block, mri_seg, pdfs, nlevels, mat, cost_type) ;
  }


  printf("writing output to %s\n", out_name) ;
  MatrixWriteTxt(out_name, mat) ;
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.mgz", out_name) ;
    printf("writing aligned volume to %s\n", fname) ;
    MRIwrite(mri_aligned, fname) ;
  }
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
  else if (!stricmp(option, "out_like") || !stricmp(option, "ol")) {
#if 0
    out_like_fname = argv[2] ;
    nargs = 1 ;
    printf("shaping output to be like %s...\n", out_like_fname) ;
#endif
  } else if (!stricmp(option, "inplane")) {
    in_plane = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting in-plane resolution to be %2.3f\n", in_plane) ;
  } else if (!stricmp(option, "nlevels")) {
    nlevels = atoi(argv[2]) ;
    nargs = 1 ;
    printf("setting nlevels=%d...\n", nlevels) ;
  } else if (!stricmp(option, "overlap")) {
    min_overlap = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting minimum fov overlap to %2.3f\n", min_overlap) ;
  } else if (!stricmp(option, "noalign")) {
    align = 0 ;
    printf("turning off alignment...\n") ;
  } else if (!stricmp(option, "max_angle")) {
    max_angle = RADIANS(atof(argv[2])) ;
    min_angle = -max_angle ;
    nargs = 1 ;
    printf("searching angles from %2.1f to %2.1f\n", DEGREES(min_angle), DEGREES(max_angle)) ;
  } else if (!stricmp(option, "max_scale")) {
    max_scale = atof(argv[2]) ;
    min_scale = 1/max_scale ;
    nargs = 1 ;
    printf("searching scales from %2.3f to %2.3f\n", min_scale, max_scale) ;
  } else if (!stricmp(option, "max_trans")) {
    max_trans = atoi(argv[2]) ;
    nargs = 1 ;
    printf("searching translations from %d to %d\n", -max_trans, max_trans) ;
  } else if (!stricmp(option, "skip")) {
    skip = atoi(argv[2]) ;
    nargs = 1 ;
    printf("sampling every %dth voxel...\n", skip) ;
    if (skip < 0)
      ErrorExit(ERROR_BADPARM, "%s: skip %d must be >= 0\n", skip) ;
  } else if (!stricmp(option, "ntrans")) {
    ntrans = atoi(argv[2]) ;
    nargs = 1 ;
    printf("searching %d translations at each step...\n", ntrans) ;
  } else if (!stricmp(option, "nscales")) {
    nscales = atoi(argv[2]) ;
    nargs = 1 ;
    printf("searching %d scales at each step...\n", nscales) ;
  } else if (!stricmp(option, "nangles")) {
    nangles = atoi(argv[2]) ;
    nargs = 1 ;
    printf("searching %d angles at each step...\n", nangles) ;
  } else if (!stricmp(option, "slice")) {
    slice_thickness = atof(argv[2]) ;
    nargs = 1 ;
    printf("setting slice thickness to be %2.3f\n", slice_thickness) ;
  } else switch (toupper(*option)) {
    case 'R':
      rotate_radians = RADIANS(atof(argv[2])) ;
      printf("rotating image by %2.2f degrees...\n", DEGREES(rotate_radians)) ;
      nargs = 1 ;
      break ;
    case 'P':
      probe = 1 ;
      printf("probing cost functional...\n") ;
      break ;
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case 'C':
      if (!stricmp(argv[2], "cr")) {
        cost_type = COST_CORRELATION_RATIO ;
        printf("using correlation ratio as cost function\n") ;
      } else if (!stricmp(argv[2], "kl")) {
        cost_type = COST_KULLABACK_LEIBLER ;
        printf("using kullback leibler divergence as cost function\n") ;
      } else if (!stricmp(argv[2], "cr")) {
        cost_type = COST_MAXIMUM_LIKELIHOOD ;
        printf("using maximum likelihood as cost function\n") ;
      }
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
          "usage: %s [options] <seg time1> <seg time 2> <transform 1> <transform 2> <output file>\n",Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program align a histological slice with a block face image\n") ;
  fprintf(stderr, "-out_like <reference volume> - set out_volume parameters\n") ;
  fprintf(stderr, "-I                           - invert transform "
          "coordinates\n") ;
  exit(1) ;
}


static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}


#if 0
static MRI *
write_slices_into_volume(MRI **mri_slices, MRI *mri_aligned, MATRIX **slice_xforms, int nslices) {
  int  slice ;
  MRI  *mri_tmp ;

  if (!mri_aligned) {
    mri_tmp = mri_slices[0] ;
    MRIcopyHeader(mri_tmp, mri_aligned) ;
  }
  for (slice = 0 ; slice < nslices ; slice++) {
    mri_tmp = mri_apply_slice_xform(mri_slices[slice], NULL, slice_xforms[slice], 0) ;
    MRIextractInto(mri_tmp, mri_aligned, 0, 0, 0, mri_aligned->width, mri_aligned->height, 1, 0, 0, slice) ;
    MRIfree(&mri_tmp) ;
  }
  return(mri_aligned) ;
}
#endif

static MRI *
align_block_to_histo(MRI *mri_histo, MRI *mri_block, MRI *mri_seg, DENSITY **pdfs, int nlevels, MATRIX *mat, int cost_type) {
  MRI  *mri_aligned ;

  mri_aligned = MRIalloc(mri_histo->width, mri_histo->height, 1, mri_histo->type) ;
  MRIcopyHeader(mri_histo, mri_aligned) ;
  MRIextractInto(mri_block, mri_aligned, 0, 0, 0, mri_aligned->width, mri_aligned->height, 1, 0, 0, 0) ;
  compute_optimal_xform(mri_block, mri_histo, mri_seg, pdfs, nlevels, mat, cost_type, skip) ;

  powell_minimize(mri_block, mri_histo, mri_seg, pdfs[0], mat, cost_type) ;
  mri_apply_slice_xform(mri_block, mri_aligned, mat, 0) ;

  return(mri_aligned) ;
}

static int
compute_optimal_xform(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY **pdfs, int nlevels, MATRIX *xform, int cost_type, int skip) {
  int  level ;
  MRI  *mri_src_pyramid[MAX_LEVELS], *mri_dst_pyramid[MAX_LEVELS], *mri_seg_pyramid[MAX_LEVELS] ;
  char fname[STRLEN] ;
  DENSITY *pdf ;

  pdf = pdfs[0] ;

  mri_src_pyramid[0] = mri_src ;
  mri_dst_pyramid[0] = mri_dst ;
  mri_seg_pyramid[0] = mri_seg ;
  for (level = 1 ; level < nlevels ; level++) {
    mri_src_pyramid[level] = MRIreduce2D(mri_src_pyramid[level-1], NULL) ;
    mri_seg_pyramid[level] = MRIreduce2D(mri_seg_pyramid[level-1], NULL) ;
    mri_dst_pyramid[level] = MRIreduce2D(mri_dst_pyramid[level-1], NULL) ;
  }

  sprintf(fname, "%s_target", base) ;
  mriWriteImageView(mri_dst, fname, RGB_SIZE, MRI_CORONAL, 0, mri_dst) ;
  for (level = nlevels-1 ; level >= 0 ; level--) {
    align_pyramid_level(mri_src_pyramid[level], mri_dst_pyramid[level], mri_seg_pyramid[level], xform, pdf, level, mri_src,
                        cost_type, skip);
  }

  for (level = 1 ; level < nlevels ; level++) {
    MRIfree(&mri_seg_pyramid[level]) ;
    MRIfree(&mri_src_pyramid[level]) ;
    MRIfree(&mri_dst_pyramid[level]) ;
  }
  return(NO_ERROR) ;
}

static int
align_pyramid_level(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, MATRIX *xform, DENSITY *pdf, int level, MRI *mri_fullres_src, int cost_type, int skip) {
  double   s, scale, a, dx, dy, mina, maxa, maxt, astep, tstep, min_error, error, mins, maxs, sstep, thick, mint ;
  MATRIX   *m_trans, *m_rot, *m_total = NULL, *m = NULL, *m_best, *m_scale, *m_origin, *m_tmp = NULL, *m_origin_inv, *m_delta, *m_start ;
  int    found ;

  thick = pow(2.0,level) ;
  skip /= thick ;

  printf("***********  aligning pyramid level %d... **************\n", level) ;
  m_start = m_delta = NULL ;

  m_rot = MatrixAlloc(3, 3, MATRIX_REAL) ;
  m_trans = MatrixIdentity(3, NULL) ;
  m_scale = MatrixAlloc(3, 3, MATRIX_REAL) ;
  m_origin = MatrixIdentity(3, NULL) ;
  *MATRIX_RELT(m_origin, 1, 3) = mri_fullres_src->width/2 ;
  *MATRIX_RELT(m_origin, 2, 3) = mri_fullres_src->height/2 ;
  *MATRIX_RELT(m_rot, 3, 3) = 1.0 ;
  *MATRIX_RELT(m_scale, 3, 3) = 1.0 ;
  m_origin_inv = MatrixInverse(m_origin, NULL) ;
  m_best = MatrixCopy(xform, NULL) ;
  min_error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_best, skip, level, cost_type) ;
  printf("starting error: %2.4f\n", min_error) ;
  write_snapshot(mri_fullres_src, mri_src, mri_dst, m_best, base, snapshots++, level, pdf, mri_seg) ;
  for (scale = 1 ; scale >= 0.005 ; scale /= 2) {
    m_start = MatrixCopy(m_best, m_start) ;
    mins = 1-scale*(1-min_scale);
    maxs = 1+scale*(max_scale-1);
    mina = min_angle*scale ;
    maxa = max_angle*scale ;
    maxt = max_trans*scale/thick ;
    astep = (maxa - mina) / (nangles-1) ;
    if (astep < 0.000000001)
      astep = 1 ;
    tstep = (2*maxt) / (ntrans-1) ;
    if (tstep < 0.00000001)
      tstep = 1;
#if 0
    if (maxt < 1)
      break ;
#endif
    if (nscales == 1)
      sstep = 1.0 ;
    else
      sstep = (maxs-mins) / (nscales-1) ;
    if (sstep < 0.00000001)
      sstep = 1.0 ;
    printf("scale %2.4f, searching angles [%2.3f --> %2.3f] by %2.3f, scales [%2.2f->%2.2f] by %2.4f, and translations to %2.3f voxels by %2.3f\n",
           scale, DEGREES(mina), DEGREES(maxa), DEGREES(astep), mins, maxs, sstep, maxt, tstep) ;
    found = compute_optimal_translation(mri_src, mri_dst, mri_seg, m_best, pdf, maxt, 50, level, skip, mri_fullres_src, cost_type) ;
    found = compute_optimal_translation(mri_src, mri_dst, mri_seg, m_best, pdf, maxt/50, 10, level, skip, mri_fullres_src, cost_type) || found ;
    if (found)
      write_snapshot(mri_fullres_src, mri_src, mri_dst, m_best, base, snapshots++, level, pdf, mri_seg) ;
    min_error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_best, skip, level, cost_type) ;
    mint = -maxt ;
    maxs += 0.5*sstep ;
    maxa += 0.5*astep ;
    maxt += 0.5*tstep;
    if (sstep < 0.00000001)
      sstep = 1 ;
    do {
      found = 0 ;
      for (s = mins ; s <= maxs ; s += sstep) {
        printf("searching scale %2.4f...\n", s) ;
        *MATRIX_RELT(m_scale, 1, 1) = *MATRIX_RELT(m_scale, 2,2) = s ;
        for (a = mina ; a <= maxa ; a += astep) {
          *MATRIX_RELT(m_rot,1,1) = cos(a) ;
          *MATRIX_RELT(m_rot,2,1) = sin(a) ;
          *MATRIX_RELT(m_rot,1,2) = -sin(a) ;
          *MATRIX_RELT(m_rot,2,2) = cos(a) ;
          m = MatrixMultiply(m_scale, m_origin_inv, m) ;
          m_tmp = MatrixMultiply(m_rot, m, m_tmp) ;
          MatrixMultiply(m_origin, m_tmp, m) ;
          for (dx = mint ; dx <= maxt ; dx += tstep) {
            *MATRIX_RELT(m_trans, 1,3) = dx ;
            if (fabs(dx) < tstep)
              DiagBreak() ;
            for (dy = mint ; dy <= maxt ; dy += tstep) {
              if (FZERO(dy) && FZERO(dx) && FEQUAL(s,1.0) && FZERO(a))
                continue ;
              if (fabs(dy) < tstep && fabs(dx)<tstep)
                DiagBreak() ;
              *MATRIX_RELT(m_trans, 2,3) = dy ;
              m_delta = MatrixMultiply(m_trans, m, m_delta) ;
              m_total = MatrixMultiply(m_delta, m_start, m_total) ;
              error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, skip, level, cost_type) ;
              if (error < min_error) {
                double e ;
                MATRIX *m_id = MatrixIdentity(3, NULL) ;
                min_error = error ;
                MatrixCopy(m_total, m_best) ;
                printf("%3.3d: best alignment at (%2.3f deg, %2.3f vox, %2.3f vox, %2.5f): %2.3f\n", snapshots,DEGREES(a), dx,dy, s, min_error) ;
                fflush(stdout) ;
                write_snapshot(mri_fullres_src, mri_src, mri_dst, m_total, base, snapshots++, level, pdf, mri_seg) ;
                found = 1 ;
                e = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_id, skip, level, cost_type) ;
                MatrixFree(&m_id) ;
              }
            }
          }
        }
      }

      if (found)
        printf("new minimum found - repeating search at scale %2.3f\n", scale) ;
    } while (found > 0);
  }

  MatrixFree(&m_rot) ;
  MatrixFree(&m_trans) ;
  MatrixFree(&m) ;
  MatrixFree(&m_scale) ;
  MatrixCopy(m_best, xform) ;
  MatrixFree(&m_best) ;
  MatrixFree(&m_origin) ;
  MatrixFree(&m_tmp);
  MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_delta) ;
  MatrixFree(&m_start) ;
  return(NO_ERROR) ;
}

static double
compute_alignment_error(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY *pdf, MATRIX *m_total, int skip, int level, int cost_type) {
  if (min_overlap > 0 && (compute_overlap(mri_src, mri_dst, m_total) < min_overlap) && !probe)
    return(100000000) ;

  switch (cost_type) {
  case COST_MAXIMUM_LIKELIHOOD:
    return(compute_ml_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, skip, level)) ;
    break ;
  case COST_CORRELATION_RATIO:
    return(compute_cr_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, skip, level)) ;
    break ;
  }
  ErrorExit(ERROR_UNSUPPORTED, "cost type %d not supported", cost_type) ;
  return(-1) ;
}
static double
compute_ml_alignment_error(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY *pdf, MATRIX *m_total, int skip, int level) {
  double sse, num, error ;
  double val_src, val_dst, xs, ys, thick ;
  int    x, y, oob, nseg ;
  VECTOR *v_src, *v_dst ;
  MATRIX *m_inv ;


  skip++ ;
  m_inv = MatrixInverse(m_total, NULL) ;

  thick = pow(2.0, level) ;
  v_src = VectorAlloc(3, MATRIX_REAL) ;
  v_dst = VectorAlloc(3, MATRIX_REAL) ;

  VECTOR_ELT(v_src, 3) = 1.0  ;
  VECTOR_ELT(v_dst, 3) = 1.0/thick  ;

  for (nseg = oob = 0, num = sse = 0.0, x = 0 ; x < mri_dst->width ; x += skip) {
    VECTOR_ELT(v_dst,1) = (double)x ;
    for (y = 0 ; y < mri_dst->height ; y += skip) {
      if (x == Gx && y == Gy)
        DiagBreak() ;
      if (MRIvox(mri_seg, x, y, 0) == 0)
        continue ;
      nseg++ ;

      VECTOR_ELT(v_dst,2) = (double)y ;
      MatrixMultiply(m_total, v_dst, v_src);
      xs = VECTOR_ELT(v_src, 1) ;
      ys = VECTOR_ELT(v_src, 2) ;

      if (MRIindexNotInVolume(mri_src, xs, ys, 0) == 0) {
        val_dst = MRIvox(mri_dst, x, y, 0) ;
        MRIsampleVolumeSlice(mri_src, xs, ys, 0, &val_src, MRI_CORONAL) ;
      } else {
        val_src = val_dst = -1000000 ;
        oob++ ;
      }
#if 0
      if ((FZERO(val_src) || val_src >= 255) && (FZERO(val_dst) || val_dst >= 255))
        error = 0 ;
      else
#endif
        error = -DensityLogLikelihood(pdf, val_dst, val_src) ; /* nissl,block */
      if (error > 16)
        DiagBreak() ;
      sse += error ;
      num++ ;
    }
  }

#if 0
  if ((oob > nseg) && !probe)
    return(100000) ;
#endif

#if 1
  for (x = 0 ; x < mri_src->width ; x += skip) {
    VECTOR_ELT(v_src,1) = (double)x ;
    for (y = 0 ; y < mri_src->height ; y += skip) {
      if (x == Gx && y == Gy)
        DiagBreak() ;

      val_src = MRIvox(mri_src, x, y, 0) ;
#if 1
      if (val_src == 0)
        continue ;   /* part of background the image is embedded in - don't use it */
#endif
      VECTOR_ELT(v_src,2) = (double)y ;
      MatrixMultiply(m_inv, v_src, v_dst);
      xs = VECTOR_ELT(v_dst, 1) ;
      ys = VECTOR_ELT(v_dst, 2) ;
      if (MRIindexNotInVolume(mri_dst, xs, ys, 0) == 0) {
        MRIsampleVolumeSlice(mri_dst, xs, ys, 0, &val_dst, MRI_CORONAL) ;
      } else {
        val_src = val_dst = -1000000 ;
      }
#if 0
      if ((FZERO(val_src) || val_src >= 255) && (FZERO(val_dst) || val_dst >= 255))
        error = 0 ;
      else
#endif
        error = -DensityLogLikelihood(pdf, val_dst, val_src) ;
      sse += error ;
      if (error > 16)
        DiagBreak() ;
      num++ ;
    }
  }
#endif

  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;
  if (FZERO(num))
    num = 1 ;
  MatrixFree(&m_inv) ;
  return(sse / num) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description

------------------------------------------------------*/
static int
mriWriteImageView(MRI *mri, char *base_name, int target_size, int view,
                  int slice, MRI *mri_template) {
  char  fname[STRLEN];
  const char *prefix ;
  IMAGE *I ;
  float scale ;

  mri = MRIresampleFill(mri, mri_template, SAMPLE_NEAREST, 255) ;

  switch (view) {
  default:
  case MRI_CORONAL:
    prefix = "cor" ;
    break ;
  case MRI_SAGITTAL:
    prefix = "sag" ;
    break ;
  case MRI_HORIZONTAL:
    prefix = "hor" ;
    break ;
  }

  if (slice < 0)    /* pick middle slice */
  {
    switch (view) {
    default:
    case MRI_CORONAL:
      slice = mri->depth/2;
      break ;
    case MRI_SAGITTAL:
      slice = 6*mri->width/10 ;
      break ;
    case MRI_HORIZONTAL:
      slice = mri->height/2 ;
      break ;
    }
  }
  I = MRItoImageView(mri, NULL, slice, view, 0) ;
  if (!I)
    ErrorReturn(Gerror, (Gerror, "MRItoImageView failed")) ;

  scale = (float)target_size / (float)I->rows ;
  if (!FEQUAL(scale, 1.0f)) {
    IMAGE *Itmp ;

    Itmp = ImageRescale(I, NULL, scale) ;
    ImageFree(&I) ;
    I = Itmp ;
  }
  sprintf(fname, "%s_%s.rgb", prefix, base_name) ;
  ImageWrite(I, fname) ;
  ImageFree(&I) ;
  MRIfree(&mri) ;
  return(NO_ERROR) ;
}
static int
write_snapshot(MRI *mri, MRI *mri1, MRI *mri2, MATRIX *m, char *base, int n, int level, DENSITY *pdf, MRI *mri_seg) {
  MRI  *mri_xformed, *mri_ll ;
  char fname[STRLEN] ;


  mri = MRIresampleFill(mri, mri2, SAMPLE_NEAREST, 255) ;
  mri_xformed = mri_apply_slice_xform(mri, NULL, m, 0) ;

  sprintf(fname, "%s%3.3d", base, n) ;
  printf("writing snapshot to %s...\n", fname) ;
  mriWriteImageView(mri_xformed, fname, RGB_SIZE, MRI_CORONAL, -1, mri2) ;
  sprintf(fname, "%s%3.3d.mgz", base, n) ;
  MRIwrite(mri_xformed, fname) ;
  MRIfree(&mri_xformed) ;

  if (pdf) {
    MATRIX *m_inv ;

    mri_ll = DensityLikelihoodImage(mri1, mri2, NULL, m, pdf, mri_seg, 0) ;
    sprintf(fname, "%s_ll_%3.3d.mgz", base, n) ;
    printf("writing density map to %s...\n", fname) ;
    MRIwrite(mri_ll, fname) ;
    MRIfree(&mri_ll) ;

    m_inv = MatrixInverse(m, NULL) ;
    mri_ll = DensityLikelihoodImage(mri2, mri1, NULL, m_inv, pdf, mri_seg, 1) ;
    sprintf(fname, "%s_ll_src_%3.3d.mgz", base, n) ;
    printf("writing density map to %s...\n", fname) ;
    MRIwrite(mri_ll, fname) ;
    MRIfree(&mri_ll) ;
    MatrixFree(&m_inv) ;
  }
  return(NO_ERROR) ;
}

static MRI *
mri_apply_slice_xform(MRI *mri_src, MRI *mri_dst, MATRIX *m, int slice) {
  double val, xs, ys ;
  int    x, y ;
  VECTOR *v_src, *v_dst ;

  if (!mri_dst)
    mri_dst = MRIcopy(mri_src, NULL) ;

  v_src = VectorAlloc(3, MATRIX_REAL) ;
  v_dst = VectorAlloc(3, MATRIX_REAL) ;

  VECTOR_ELT(v_src, 3) = 1.0  ;
  VECTOR_ELT(v_dst, 3) = 1.0  ;
  mri_dst = MRIclone(mri_src, NULL) ;
  for (x = 0 ; x < mri_dst->width ; x++) {
    VECTOR_ELT(v_dst,1) = (double)x ;
    for (y = 0 ; y < mri_dst->height ; y++) {
      VECTOR_ELT(v_dst,2) = (double)y ;
      MatrixMultiply(m, v_dst, v_src);
      xs = VECTOR_ELT(v_src, 1) ;
      ys = VECTOR_ELT(v_src, 2) ;
      MRIsampleVolumeSlice(mri_src, xs, ys, slice, &val, MRI_CORONAL) ;
      MRIsetVoxVal(mri_dst, x, y, slice, 0, val) ;
    }
  }

  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;

  return(mri_dst) ;
}
static int
compute_optimal_translation(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, MATRIX *m_best, DENSITY *pdf,
                            int maxt, int ntrans, int level, int skip, MRI *mri_fullres_src, int cost_type) {
  MATRIX  *m_trans, *m_start, *m_total ;
  double  thick, tstep, dx, dy, error, min_error, mint ;
  int     found = 0 ;

  m_trans = MatrixIdentity(3, NULL) ;


  min_error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_best, skip, level, cost_type) ;

  tstep = (2*maxt) / (ntrans-1) ;
  if (tstep < 0.1)
    return(NO_ERROR) ;
  thick = pow(2.0,level) ;
  if (maxt < 1)
    return(NO_ERROR) ;

  printf("computing optimal translations from %d --> %d in %2.1f voxel steps\n", -maxt, maxt, tstep) ;
  m_start = MatrixCopy(m_best, NULL) ;
  m_total = NULL ;
  min_error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_best, skip, level, cost_type) ;
  mint = -maxt ;
  maxt += 0.5*tstep ;
  for (dx = mint ; dx <= maxt ; dx += tstep) {
    *MATRIX_RELT(m_trans, 1,3) = dx ;
    for (dy = mint ; dy <= maxt ; dy += tstep) {
      *MATRIX_RELT(m_trans, 2,3) = dy ;
      m_total = MatrixMultiply(m_trans, m_start, m_total) ;
      error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, skip, level, cost_type) ;
      if (error < min_error) {
        min_error = error ;
        MatrixCopy(m_total, m_best) ;
        printf("     best translation found at (%2.3f vox, %2.3f vox): %2.3f\n", dx,dy, min_error) ;
        fflush(stdout) ;
        found = 1 ;
        /*    write_snapshot(mri_fullres_src, m_total, base, snapshots++, level) ;*/
      }
    }
  }

  MatrixFree(&m_trans) ;
  MatrixFree(&m_total) ;
  MatrixFree(&m_start) ;
  return(found) ;
}

static double
compute_cr_alignment_error(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY *pdf,
                           MATRIX *m_total, int skip, int level) {
  double val_src, val_dst, xs, ys, thick ;
  int    x, y, oob, i ;
  VECTOR *v_src, *v_dst ;
  MATRIX *m_inv ;
  int    nvox[256], nseg ;
  double  means[256], vars[256], cr, mean, var, num /*, cr2*/ ;

  skip++ ;
  m_inv = MatrixInverse(m_total, NULL) ;

  thick = pow(2.0, level) ;
  v_src = VectorAlloc(3, MATRIX_REAL) ;
  v_dst = VectorAlloc(3, MATRIX_REAL) ;

  VECTOR_ELT(v_src, 3) = 1.0  ;
  VECTOR_ELT(v_dst, 3) = 1.0/thick  ;

  memset(means, 0, sizeof(means)) ;
  memset(vars, 0, sizeof(vars)) ;
  memset(nvox, 0, sizeof(nvox)) ;
  mean = var = 0.0 ;
  for (nseg = oob = 0,  num = x = 0 ; x < mri_dst->width ; x += skip) {
    VECTOR_ELT(v_dst,1) = (double)x ;
    for (y = 0 ; y < mri_dst->height ; y += skip) {
      if (x == Gx && y == Gy)
        DiagBreak() ;
      if (MRIvox(mri_seg, x, y, 0) == 0)
        continue ;
      nseg++ ;

      VECTOR_ELT(v_dst,2) = (double)y ;
      MatrixMultiply(m_total, v_dst, v_src);
      xs = VECTOR_ELT(v_src, 1) ;
      ys = VECTOR_ELT(v_src, 2) ;
      if (MRIindexNotInVolume(mri_src, xs, ys, 0) == 0)  /* a valid point */
      {
        MRIsampleVolumeSlice(mri_src, xs, ys, 0, &val_src, MRI_CORONAL) ;
        val_dst = MRIvox(mri_dst, x, y, 0) ;

        i = nint(val_dst) ;
        if (i == Gdiag_no)
          DiagBreak() ;
        means[i] += val_src ;
        vars[i] += val_src*val_src ;
        nvox[i]++ ;
        mean += val_src ;
        var += val_src*val_src ;
        num++ ;
      } else
        oob++ ;
    }
  }

  if (num == 0)
    return(1.0) ;

  if ((num < nseg/2) && !probe)
    return(1.0) ;


  mean /= num ;
  var = var / num - mean*mean ;
  for (i = 0 ; i < 256 ; i++)
    if (nvox[i]) {
      means[i] /= nvox[i] ;
      vars[i] = vars[i] / nvox[i] - means[i]*means[i] ;
    }

  for (cr = 0.0, i = 0 ; i < 256 ; i++) {
    cr += nvox[i]*vars[i] ;
  }
  cr /= (num*var) ;
  /* cr += nseg/num ;*/

#if 0
  mean = var = 0.0 ;
  memset(means, 0, sizeof(means)) ;
  memset(vars, 0, sizeof(vars)) ;
  memset(nvox, 0, sizeof(nvox)) ;
  for (oob = 0,  num = x = 0 ; x < mri_src->width ; x += skip) {
    VECTOR_ELT(v_src,1) = (double)x ;
    for (y = 0 ; y < mri_src->height ; y += skip) {
      if (x == Gx && y == Gy)
        DiagBreak() ;
      VECTOR_ELT(v_src,2) = (double)y ;
      MatrixMultiply(m_inv, v_src, v_dst);
      xs = VECTOR_ELT(v_dst, 1) ;
      ys = VECTOR_ELT(v_dst, 2) ;
      if (MRIindexNotInVolume(mri_dst, xs, ys, 0) == 0) {
        MRIsampleVolumeSlice(mri_src, x, y, 0, &val_src, MRI_CORONAL) ;
        MRIsampleVolumeSlice(mri_dst, xs, ys, 0, &val_dst, MRI_CORONAL) ;
        i = nint(val_dst) ;
        if (i == Gdiag_no)
          DiagBreak() ;
        means[i] += val_src ;
        vars[i] += val_src*val_src ;
        nvox[i]++ ;
        mean += val_src ;
        var += val_src*val_src ;
        num++ ;
      }
    }
  }

  mean /= num ;
  var = var / num - mean*mean ;
  for (i = 0 ; i < 256 ; i++)
    if (nvox[i]) {
      means[i] /= nvox[i] ;
      vars[i] = vars[i] / nvox[i] - means[i]*means[i] ;
    }

  for (cr2 = 0.0, i = 0 ; i < 256 ; i++) {
    cr2 += nvox[i]*vars[i] ;
  }
  cr2 /= (num*var) ;
  cr = (cr + cr2) / 2 ;
#endif

  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;
  MatrixFree(&m_inv) ;
  return(cr) ;  /* actually 1-cr */
}

static MRI *
rotate_image(MRI *mri_src, MRI *mri_dst, double rotate_radians) {
  MATRIX   *m_rot, *m_total, *m_origin, *m_origin_inv ;
  VECTOR   *v_src, *v_dst ;
  float    xs, ys ;
  double   val_src ;
  int      x, y ;

  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  v_src = VectorAlloc(3, MATRIX_REAL) ;
  v_dst = VectorAlloc(3, MATRIX_REAL) ;

  VECTOR_ELT(v_src, 3) = 1.0  ;
  VECTOR_ELT(v_dst, 3) = 1.0 ;
  m_origin = MatrixIdentity(3, NULL) ;
  *MATRIX_RELT(m_origin, 1, 3) = mri_src->width/2 ;
  *MATRIX_RELT(m_origin, 2, 3) = mri_src->height/2 ;
  m_rot = MatrixAlloc(3, 3, MATRIX_REAL) ;
  *MATRIX_RELT(m_rot,1,1) = cos(rotate_radians) ;
  *MATRIX_RELT(m_rot,2,1) = sin(rotate_radians) ;
  *MATRIX_RELT(m_rot,1,2) = -sin(rotate_radians) ;
  *MATRIX_RELT(m_rot,2,2) = cos(rotate_radians) ;
  *MATRIX_RELT(m_rot, 3, 3) = 1.0 ;
  m_origin_inv = MatrixInverse(m_origin, NULL) ;
  m_total = MatrixMultiply(m_rot, m_origin_inv, NULL) ;
  MatrixMultiply(m_origin, m_total, m_origin_inv) ;
  MatrixCopy(m_origin_inv, m_total) ;

  for (x = 0 ; x < mri_dst->width ; x++) {
    VECTOR_ELT(v_dst,1) = (double)x ;
    for (y = 0 ; y < mri_dst->height ; y++) {
      if (x == Gx && y == Gy)
        DiagBreak() ;
      VECTOR_ELT(v_dst,2) = (double)y ;
      MatrixMultiply(m_total, v_dst, v_src);
      xs = VECTOR_ELT(v_src, 1) ;
      ys = VECTOR_ELT(v_src, 2) ;
      MRIsampleVolumeSlice(mri_src, xs, ys, 0, &val_src, MRI_CORONAL) ;
      MRIsetVoxVal(mri_dst, x, y, 0, 0, val_src) ;
    }
  }


  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;
  MatrixFree(&m_rot) ;
  MatrixFree(&m_origin) ;
  MatrixFree(&m_origin_inv) ;
  MatrixFree(&m_total) ;
  return(mri_dst) ;
}

int
MRIeraseImageBorder(MRI *mri, int width) {
  int x, y, i ;

  for (x = 0 ; x < mri->width ; x++) {
    for (i = 0 ; i < width ; i++) {
      MRIsetVoxVal(mri, x, i, 0, 0, 0) ;
      MRIsetVoxVal(mri, x, mri->height-1-i, 0, 0, 0) ;
    }
  }
  for (y = 0 ; y < mri->height ; y++) {
    for (i = 0 ; i < width ; i++) {
      MRIsetVoxVal(mri, i, y, 0, 0, 0) ;
      MRIsetVoxVal(mri, mri->width-1-i, y, 0, 0, 0) ;
    }
  }
  return(NO_ERROR) ;
}
static double
compute_overlap(MRI *mri_src, MRI *mri_dst, MATRIX *m_total) {
  double xs, ys, thick ;
  int    x, y, nvox, ntotal ;
  VECTOR *v_src, *v_dst ;
  MATRIX *m_inv ;
  BUFTYPE val ;


  thick = mri_src->xsize ;
  v_src = VectorAlloc(3, MATRIX_REAL) ;
  v_dst = VectorAlloc(3, MATRIX_REAL) ;
  m_inv = MatrixInverse(m_total, NULL) ;

  VECTOR_ELT(v_src, 3) = 1.0  ;
  VECTOR_ELT(v_dst, 3) = 1.0/thick  ;

  for (ntotal = nvox = x = 0 ; x < mri_dst->width ; x++) {
    VECTOR_ELT(v_dst,1) = (double)x ;
    for (y = 0 ; y < mri_dst->height ; y++) {
      if (x == Gx && y == Gy)
        DiagBreak() ;
      val = MRIvox(mri_dst, x, y, 0) ;
      if (val == 0 || val == 255)  /* background - don't penalize for this */
        continue ;

      VECTOR_ELT(v_dst,2) = (double)y ;
      MatrixMultiply(m_total, v_dst, v_src);
      xs = VECTOR_ELT(v_src, 1) ;
      ys = VECTOR_ELT(v_src, 2) ;

      ntotal++ ;
      if (MRIindexNotInVolume(mri_src, xs, ys, 0) == 0)
        nvox++ ;   /* count # of voxels in dst that map within fov of src */
    }
  }

  for (x = 0 ; x < mri_src->width ; x++) {
    VECTOR_ELT(v_src,1) = (double)x ;
    for (y = 0 ; y < mri_src->height ; y++) {
      if (x == Gx && y == Gy)
        DiagBreak() ;

      val = MRIvox(mri_src, x, y, 0) ;
      if (val == 0 || val == 255)  /* background - don't penalize for this */
        continue ;
      VECTOR_ELT(v_src,2) = (double)y ;
      MatrixMultiply(m_inv, v_src, v_dst);
      xs = VECTOR_ELT(v_dst, 1) ;
      ys = VECTOR_ELT(v_dst, 2) ;

      ntotal++ ;
      if (MRIindexNotInVolume(mri_dst, xs, ys, 0) == 0)
        nvox++ ;   /* count # of voxels in dst that map within fov of src */
    }
  }

  MatrixFree(&m_inv) ;
  VectorFree(&v_src) ;
  VectorFree(&v_dst) ;
  return((double)nvox / (double)ntotal) ;
}
#define PROBE_SAMPLES 5000

static int
probe_cost_function(MRI *mri_src, MRI *mri_dst, MRI *mri_seg, DENSITY *pdf, MATRIX *mat, int cost_type, char *out_name) {
  double   s, a, dx, dy, maxt, astep, tstep, error, sstep ;
  MATRIX   *m_rot, *m_total = NULL, *m = NULL, *m_scale, *m_origin, *m_origin_inv ;
  FILE   *fp ;
  char   fname[STRLEN] ;

  m_rot = MatrixAlloc(3, 3, MATRIX_REAL) ;
  m_scale = MatrixAlloc(3, 3, MATRIX_REAL) ;
  m_origin = MatrixIdentity(3, NULL) ;
  *MATRIX_RELT(m_origin, 1, 3) = mri_src->width/2 ;
  *MATRIX_RELT(m_origin, 2, 3) = mri_src->height/2 ;
  *MATRIX_RELT(m_rot, 3, 3) = 1.0 ;
  *MATRIX_RELT(m_scale, 3, 3) = 1.0 ;
  m_origin_inv = MatrixInverse(m_origin, NULL) ;

  printf("probing scales...\n") ;
  sstep = (max_scale-min_scale) / (PROBE_SAMPLES-1) ;
  sprintf(fname, "%s_scale.dat", out_name) ;
  fp = fopen(fname, "w") ;
  for (s = min_scale ; s <= max_scale+.5*sstep ; s += sstep) {
    *MATRIX_RELT(m_scale, 1, 1) = *MATRIX_RELT(m_scale, 2,2) = s ;
    m = MatrixMultiply(m_scale, m_origin_inv, m) ;
    m_total = MatrixMultiply(m_origin, m, m_total) ;
    error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, 0, 0, cost_type) ;
    fprintf(fp, "%f %f\n", s, error) ;
    fflush(fp) ;
  }
  fclose(fp) ;

  printf("probing rotations...\n") ;
  sprintf(fname, "%s_rot.dat", out_name) ;
  fp = fopen(fname, "w") ;
  astep = (2*max_angle) / (PROBE_SAMPLES-1);
  for (a = min_angle ; a <= max_angle+.5*astep ; a += astep) {
    *MATRIX_RELT(m_rot,1,1) = cos(a) ;
    *MATRIX_RELT(m_rot,2,1) = sin(a) ;
    *MATRIX_RELT(m_rot,1,2) = -sin(a) ;
    *MATRIX_RELT(m_rot,2,2) = cos(a) ;
    m = MatrixMultiply(m_rot, m_origin_inv, m) ;
    MatrixMultiply(m_origin, m, m_total) ;
    error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, 0, 0, cost_type) ;
    fprintf(fp, "%f %f\n", DEGREES(a), error) ;
    fflush(fp) ;
  }
  fclose(fp) ;

  printf("probing x translations...\n") ;
  sprintf(fname, "%s_dx.dat", out_name) ;
  fp = fopen(fname, "w") ;
  MatrixIdentity(3, m_total) ;
  *MATRIX_RELT(m_total, 2,3) = 0 ;
  *MATRIX_RELT(m_total, 3,3) = 1.0 ;
  maxt = mri_src->width ;
  tstep = 1 ;
  for (dx = -max_trans ; dx <= max_trans+.5*tstep ; dx += tstep) {
    *MATRIX_RELT(m_total, 1,3) = dx ;
    error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, 0, 0, cost_type) ;
    fprintf(fp, "%f %f\n", dx, error) ;
    fflush(fp) ;
  }
  fclose(fp) ;

  printf("probing y translations...\n") ;
  sprintf(fname, "%s_dy.dat", out_name) ;
  fp = fopen(fname, "w") ;
  *MATRIX_RELT(m_total, 1,3) = 0 ;
  *MATRIX_RELT(m_total, 3,3) = 1.0 ;
  maxt = mri_src->height ;
  tstep = 1 ;
  for (dy = -max_trans ; dy <= max_trans+.5*tstep ; dy += tstep) {
    *MATRIX_RELT(m_total, 2,3) = dy ;
    error = compute_alignment_error(mri_src, mri_dst, mri_seg, pdf, m_total, 0, 0, cost_type) ;
    fprintf(fp, "%f %f\n", dy, error) ;
    fflush(fp) ;
  }
  fclose(fp) ;

  MatrixFree(&m_rot) ;
  MatrixFree(&m) ;
  MatrixFree(&m_scale) ;
  MatrixFree(&m_origin) ;
  MatrixFree(&m_origin_inv) ;
  return(NO_ERROR) ;
}
#define NPARMS (3*3)
#ifdef TOL
#undef TOL
#endif
#define TOL 1e-12

static MRI *Gmri_histo, *Gmri_block, *Gmri_seg ;
static DENSITY *Gpdf ;
static int Gcost_type ;

static float compute_powell_sse(float *p) ;
static int
powell_minimize(MRI *mri_block, MRI *mri_histo, MRI *mri_seg, DENSITY *pdf, MATRIX *mat, int cost_type) {
  float *p, **xi, fret, fstart;
  int   i, r, c, iter ;

  p = vector(1, NPARMS) ;
  xi = matrix(1, NPARMS, 1, NPARMS) ;
  for (i = r = 1 ; r <= 3 ; r++) {
    for (c = 1 ; c <= 3 ; c++) {
      p[i++] = *MATRIX_RELT(mat, r, c) ;
    }
  }

  Gmri_block = mri_block ;
  Gmri_histo = mri_histo ;
  Gmri_seg = mri_seg ;
  Gpdf = pdf ;
  Gcost_type = cost_type ;
  for (r = 1 ; r <= NPARMS ; r++) {
    for (c = 1 ; c <= NPARMS ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
  do {
    for (r = 1 ; r <= NPARMS ; r++) {
      for (c = 1 ; c <= NPARMS ; c++) {
        xi[r][c] = r == c ? 1 : 0 ;
      }
    }

    fstart = fret ;
    OpenPowell(p, xi, NPARMS, TOL, &iter, &fret, compute_powell_sse);
    for (i = r = 1 ; r <= 3 ; r++) {
      for (c = 1 ; c <= 3 ; c++) {
        *MATRIX_RELT(mat, r, c) = p[i++] ;
      }
    }
    printf("%3.3d: best alignment at after powell: %2.3f (%d steps)\n", snapshots,fret, iter) ;
    write_snapshot(mri_block, mri_block, mri_histo, mat, base, snapshots++, 0, pdf, mri_seg) ;
  } while (fret < fstart) ;

  free_matrix(xi, 1, NPARMS, 1, NPARMS) ;
  free_vector(p, 1, NPARMS) ;
  return(NO_ERROR) ;
}

static float
compute_powell_sse(float *p) {
  static MATRIX *mat = NULL ;
  float  error ;
  int    i, r, c ;

  if (mat == NULL)
    mat = MatrixAlloc(3, 3, MATRIX_REAL) ;
  for (i = r = 1 ; r <= 3 ; r++) {
    for (c = 1 ; c <= 3 ; c++) {
      *MATRIX_RELT(mat, r, c) = p[i++] ;
    }
  }
  error = compute_alignment_error(Gmri_block, Gmri_histo, Gmri_seg, Gpdf, mat, skip, 0, Gcost_type) ;
  return(error) ;
}

