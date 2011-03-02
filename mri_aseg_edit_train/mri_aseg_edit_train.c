/**
 * @file  mri_aseg_edit_train
 * @brief train a classifier based on manual corrections of asegs.
 *
 * use SVMs to learn to correct an aseg.
 */
/*
 * Original Author: Bruce Fischl
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:13 $
 *    $Revision: 1.3 $
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


#include "timer.h"
#include "mrinorm.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "const.h"
#include "proto.h"
#include "mrisurf.h"
#include "macros.h"
#include "fio.h"
#include "mrishash.h"
#include "cvector.h"
#include "svm.h"
#include "version.h"
#include "voxlist.h"
#include "cma.h"
#include "class_array.h"

static char vcid[] = "$Id: mri_aseg_edit_train.c,v 1.3 2011/03/02 00:04:13 nicks Exp $";


/*-------------------------------- CONSTANTS -----------------------------*/

/*-------------------------------- PROTOTYPES ----------------------------*/

int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static MATRIX *compute_ras_basis_vectors(MRI *mri_aseg_orig, MRI *mri_aseg_edit, int target_label, 
                                         MATRIX *mEvectors, double *means, int width, int height, int depth, int pad) ;


/*-------------------------------- DATA ----------------------------*/

char *Progname ;

static char *wfile_name ;
static char *c1_name = "Left_Hippocampus" ;
static char *c2_name = "Left_fimbria" ;
static char *sdir = NULL ;

static char *test_subject = NULL ;
static char *prefix = "" ;
static double svm_tol = DEFAULT_SVM_TOL ;
static double svm_C = DEFAULT_SVM_C ;
static double momentum = 0.0 ;
static double rbf_sigma = 0.0 ;
static double poly_d = 0.0 ;
static int svm_max_iter = 1000000 ;

static int ca_type = CA_GAUSSIAN ; //CA_SVM ;
static int ca_width = 8 ;

static int navgs = 0 ;

#define TEST_SVM 0
#if TEST_SVM
float x_[20][2] = {
                    {
                      -1.0000,    1.0000
                    },
                    {-0.5000,    1.5000},
                    {0,    2.0000},
                    {0.8633,    3.2528},
                    {-2.1805,    2.2737},
                    {-0.6757,    4.0846},
                    {-1.2589,    3.8058},
                    {-0.8013,    3.6031},
                    {-1.3194,    4.3115},
                    {-1.5635,    2.0287},
                    {1.0000,   -1.0000},
#if 0
                    {1.0000,   -2.0000},
#else
                    {1.5000,   -0.5000},
#endif
                    {2.0000,         0},
                    {2.2501,   -1.1316},
                    {3.0240,   -1.0661},
                    {1.0503,   -2.0306},
                    {2.7714,   -1.1275},
                    {2.2957,   -1.9387},
                    {3.2611,   -2.2788},
                    {3.1396,    0.1683}
                  } ;
float y[20] = {
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1,
                -1
              } ;
#endif

static char *aseg_edit_name = "aseg.mgz" ;
static char *aseg_orig_name = "aseg.auto.mgz" ;
static char norm_name[STRLEN] = "norm.mgz" ;

static float sigmas[] = { 0, .5, 1.0, 2.0 } ;
#define NSCALES (sizeof(sigmas) / sizeof(sigmas[0]))

static int target_label = Left_fimbria ;
static int source_label = Left_Hippocampus ;
static int which_inputs = CA_INPUT_D2I_S | CA_INPUT_INTENSITY ;

static VOXEL_LIST *vlst_add_border_voxels(MRI *mri1, MRI *mri2, VOXEL_LIST *vl_dif, VOXEL_LIST *vl_total, int target) ;

static int wsize = 1 ;
static int nscales = 1 ;


/*-------------------------------- FUNCTIONS ----------------------------*/

int
main(int argc, char *argv[]) {
  SVM          *svm ;
  MRI          *mri_aseg_orig, *mri_aseg_edit, *mri_norm ;
  char         **av, fname[STRLEN], *cp, *subject_name, subjects_dir[STRLEN],
               *out_fname ;
  int          ac, nargs, n, nsubjects ;
  struct timeb start ;
  int          msec, minutes, seconds ;
  VOXEL_LIST   *vl_dif, *vl_total ;
#if 0
  float        *classes, **inputs ;
  int          i ;
#endif
  double       means[3] ;
  float        evalues[3] ;
  MATRIX       *m_evectors, *m_xform ;
  CA           *ca = NULL;
  
 m_evectors = MatrixAlloc(3, 3, MATRIX_REAL) ;   /* eigenvectors of label */
  

#if TEST_SVM
  {
    float **x ;
    int   ntraining = 20, i, j, k ;
    double sum, w[2], out ;


    x = (float **)calloc(ntraining, sizeof(float)) ;
    for (i = 0 ; i < ntraining ; i++) {
      x[i] = (float *)calloc(2, sizeof(float)) ;/autofs/space/birn_036/users/czanner/Phantom_tests/test_7T_phantom

      x[i][0] = x_[i][0] ;
      x[i][1] = x_[i][1] ;
    }

    svm = SVMalloc(2, "c1", "c2") ;
    for (k = 0; k < 1 ; k++) {
      srandom(k);
      SVMtrain(svm, x, y, 20, svm_C, svm_tol, 100000) ;

      sum = SVMconstraint(svm, y, ntraining) ;
      printf("constraint = %f\n", sum) ;
      for (j = 0 ; j < 2 ; j++) {
        w[j] = 0 ;
        for (i = 0 ; i < svm->nsupport ; i++)
          w[j] += svm->asupport[i] * svm->ysupport[i] * svm->xsupport[i][j] ;
      }
#if 1
      for (i = 0 ; i < ntraining ; i++) {
        out = SVMclassify(svm, x[i]) ;
        printf("%d (%2.1f,%2.1f) (%2.4f): %f (%f) --> %s\n",
               i, x[i][0], x[i][1], svm->alpha[i],
               y[i], out, y[i]*out > 0 ? "CORRECT" : "INCORRECT") ;
      }
#endif
      printf("%d support vectors found, weights = %2.3f, %2.3f\n",
             svm->nsupport, w[0], w[1]) ;
    }
    exit(0) ;
  }
#endif

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_aseg_edit_train.c,v 1.3 2011/03/02 00:04:13 nicks Exp $", "$Name:  $");
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

  TimerStart(&start) ;

  /* subject_name $class1 : $class2 output*/
  if (argc < 3)
    usage_exit() ;

  if (sdir)
    cp = sdir ;
  else
    cp = getenv("SUBJECTS_DIR") ;
  if (!cp)
    ErrorExit(ERROR_BADPARM, "%s: SUBJECTS_DIR not defined in environment",
              Progname) ;

  strcpy(subjects_dir, cp) ;

  out_fname = argv[argc-1] ;
  printf("writing output to %s...\n", out_fname) ;

#define ARGV_OFFSET 1


  // every scale at each window location for each of the dx,dy,dz grads and the original image
  svm = SVMalloc(NINPUTS(which_inputs, wsize, nscales), c1_name, c2_name) ;

  /* first determine the number of subjects in each class */
  n = ARGV_OFFSET ;
  nsubjects = argc-ARGV_OFFSET-1 ; ;
  printf("training on %d subjects\n", nsubjects) ;
  subject_name = argv[ARGV_OFFSET] ;

  /* real all the asegas in for group1 */
  for (n = 0 ; n < nsubjects ; n++) {
    /* transform each subject's curvature into the output subject's space */
    subject_name = argv[n+ARGV_OFFSET] ;
    fprintf(stderr, "reading subject %d of %d: %s\n",
            n+1, nsubjects, subject_name) ;
    sprintf(fname, "%s/%s/mri/%s",
            subjects_dir,subject_name,aseg_orig_name);
    mri_aseg_orig = MRIread(fname) ;
    if (!mri_aseg_orig)
      ErrorExit(ERROR_NOFILE, "%s: could not read aseg orig file %s",
                Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s",
            subjects_dir,subject_name,aseg_edit_name);
    mri_aseg_edit = MRIread(fname) ;
    if (!mri_aseg_edit)
      ErrorExit(ERROR_NOFILE, "%s: could not read aseg edit file %s",
                Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s",
            subjects_dir,subject_name,norm_name);
    mri_norm = MRIread(fname) ;
    if (!mri_norm)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    MRIprincipleComponentsRange(mri_aseg_edit, m_evectors, evalues,
                                means, source_label, source_label) ;
    printf("%s mean (%2.1lf, %2.1lf, %2.1lf), evectors:\n", cma_label_to_name(target_label),means[0], means[1], means[2]) ;
    MatrixPrint(Gstdout, m_evectors) ;
    m_xform = compute_ras_basis_vectors(mri_aseg_orig, mri_aseg_edit, 
                                        source_label, m_evectors, means, 
                                        ca_width, ca_width, ca_width, 4) ;

    if (ca == NULL) // 1st time
      ca = CAalloc(ca_width, ca_width, ca_width, m_xform, ca_type,which_inputs,
                   wsize, nscales, c1_name, c2_name, sigmas) ;
    else
      MatrixCopy(m_xform, ca->m_vox2index) ;

    strcpy(ca->c1_name, c1_name) ; strcpy(ca->c2_name, c2_name) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRI *mri_aligned = MRIlinearTransformInterp(mri_norm, NULL, m_xform, SAMPLE_TRILINEAR) ;
      MRIwrite(mri_aligned, "a.mgz") ;
      MRIfree(&mri_aligned) ;
    }

    MatrixFree(&m_xform) ; MatrixFree(&m_evectors) ;

#if 0
    vl_dif = VLSTcreateFromDifference(mri_aseg_orig, mri_aseg_edit, NULL, target_label) ;
    vl_total = vlst_add_border_voxels(mri_aseg_orig, mri_aseg_edit, vl_dif, NULL, target_label) ;
#else
    vl_dif = VLSTcreate(mri_aseg_edit, target_label, target_label, NULL, 0, 0);
    vl_total = vlst_add_border_voxels(mri_aseg_orig, mri_aseg_edit, vl_dif, NULL, target_label) ;
#endif
    if (ca_type == CA_SVM)
      CAsetSVMparms(ca, svm_C, svm_tol, 1000) ;
    CAtrain(ca, vl_total, mri_norm, mri_aseg_orig, source_label, target_label);

    MRIfree(&mri_aseg_orig) ; MRIfree(&mri_aseg_edit) ; MRIfree(&mri_norm) ;
    VLSTfree(&vl_dif) ; VLSTfree(&vl_total) ;
  }

  CAcompleteTraining(ca) ;

  {
    MRI *mri_output ;

    sprintf(fname, "%s/%s/mri/%s",
            subjects_dir,subject_name,aseg_orig_name);
    mri_aseg_orig = MRIread(fname) ;
    if (!mri_aseg_orig)
      ErrorExit(ERROR_NOFILE, "%s: could not read aseg orig file %s",
                Progname, fname) ;

    sprintf(fname, "%s/%s/mri/%s",
            subjects_dir,subject_name,norm_name);
    mri_norm = MRIread(fname) ;
    if (!mri_norm)
      ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
                Progname, fname) ;
    mri_output = CAclassifyBorder(ca, mri_norm, mri_aseg_orig, NULL,
                                       4, source_label) ;
    MRIwrite(mri_output, "out.mgz") ;
  }
  printf("writing trained classifier to %s...\n", out_fname) ;
  CAwrite(ca, out_fname) ;
  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  printf("CA training took %d minutes and %d seconds.\n",
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
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help"))
    print_help() ;
  else if (!stricmp(option, "-version"))
    print_version() ;
  else if (!stricmp(option, "test")) {
    test_subject = argv[2] ;
    fprintf(stderr, "writing test.dat for subject %s\n", test_subject) ;
    nargs = 1 ;
  } else if (!stricmp(option, "width")) {
    ca_width = atoi(argv[2]) ;
    fprintf(stderr, "using (%d x %d x %d) classifier array\n",
            ca_width, ca_width, ca_width) ;
    nargs = 1 ;
  } else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    fprintf(stderr, "debugging voxel (%d, %d, %d)\n", Gx, Gy, Gz) ;
    nargs = 3 ;
  } else if (!stricmp(option, "aseg_edit")) {
    aseg_edit_name = argv[2] ;
    fprintf(stderr, "using %s as edited aseg volume\n", aseg_edit_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "aseg_orig")) {
    aseg_orig_name = argv[2] ;
    fprintf(stderr, "using %s as orig aseg volume\n", aseg_orig_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "sdir")) {
    sdir = argv[2] ;
    printf("using %s as SUBJECTS_DIR\n", sdir) ;
    nargs = 1 ;
  } else if (!stricmp(option, "c1")) {
    c1_name = argv[2] ;
    printf("using %s as name for class 1\n", c1_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "c2")) {
    c2_name = argv[2] ;
    printf("using %s as name for class 2\n", c2_name) ;
    nargs = 1 ;
  } else if (!stricmp(option, "max")) {
    svm_max_iter = atoi(argv[2]) ;
    printf("using max iter = %d...\n", svm_max_iter) ;
    nargs = 1 ;
  } else if (!stricmp(option, "rbf")) {
    rbf_sigma = atof(argv[2]) ;
    printf("using RBF SVM with sigma = %f\n", rbf_sigma) ;
    nargs = 1 ;
  } else if (!stricmp(option, "poly")) {
    poly_d = atof(argv[2]) ;
    printf("using polynomial SVM with dimension = %f\n", poly_d) ;
    nargs = 1 ;
  } else if (!stricmp(option, "tol")) {
    svm_tol = atof(argv[2]) ;
    printf("using integration tolerance = %e\n", svm_tol) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'L':
      target_label = atoi(argv[2]) ;
      fprintf(stderr, "resegmenting target label %s (%d)\n", cma_label_to_name(target_label), target_label) ;
      nargs = 1 ;
      break ;
    case 'C':
      svm_C = atof(argv[2]) ;
      printf("using C=%f for svm slack variable weight\n", svm_C) ;
      nargs = 1 ;
      break ;
    case 'A':
      navgs = atoi(argv[2]) ;
      nargs = 1 ;
      printf("averaging values %d times\n", navgs) ;
      break ;
    case 'P':
      prefix = argv[2] ;
      fprintf(stderr, "using label prefix %s\n", prefix) ;
      nargs = 1 ;
      break ;
    case 'W':
      wfile_name = argv[2] ;
      printf("writing svm to w file %s...\n", wfile_name) ;
      nargs = 1 ;
      break ;
    case 'M':
      momentum = atof(argv[2]) ;
      printf("using momentum = %f for SVM training.\n", momentum) ;
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
  print_usage() ;
  exit(1) ;
}

static void
print_usage(void) {
  fprintf(stderr,
          "usage: %s -o <output subject> [options] \n"
          "\t\n\t<subject1> <subject2>... <output file> \n", Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program will compute train an SVM classifier to correct an aseg\n") ;
  fprintf(stderr, "\nvalid options are:\n\n") ;
  exit(1) ;
}

static void
print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}


static VOXEL_LIST *
vlst_add_border_voxels(MRI *mri1, MRI *mri2, VOXEL_LIST *vl_dif, VOXEL_LIST *vl_total, int target)
{
  int  x, y, z, non, noff, i, xk, yk, zk, xi, yi, zi, total, in_set_already, 
       nbr, border, found, nvox, l1, l2, total_on ;
  MRI  *mri, *mri_dil, *mri_border ;

  non = noff = 0 ;
  for (i = 0 ; i < vl_dif->nvox ; i++)
  {
    x = vl_dif->xi[i] ;
    y = vl_dif->yi[i] ;
    z = vl_dif->zi[i] ;
    if (vl_dif->vsrc[i] == target)  // originally target -> turned off
      noff++ ;
    else if (vl_dif->vdst[i] == target) // otherwise was edited for a different reason
      non++ ;
  }

  mri = MRIclone(mri1, NULL) ;
  if (non > noff)  // add some off ones
  {
    vl_total = VLSTalloc(2*non) ;
    VLSTcopy(vl_dif, vl_total, 0, non+noff) ;
    for (i = 0 ; i < vl_dif->nvox ; i++)
    {
      x = vl_dif->xi[i] ;
      y = vl_dif->yi[i] ;
      z = vl_dif->zi[i] ;
      if (vl_dif->vsrc[i] != target)
        MRIsetVoxVal(mri, x, y, z, 0, 1) ;
    }
  }
  else if (noff > non)
  {
    vl_total = VLSTalloc(2*noff) ;
    VLSTcopy(vl_dif, vl_total, 0, non+noff) ;
    for (i = 0 ; i < vl_dif->nvox ; i++)
    {
      x = vl_dif->xi[i] ;
      y = vl_dif->yi[i] ;
      z = vl_dif->zi[i] ;
      if (vl_dif->vsrc[i] == target)
        MRIsetVoxVal(mri, x, y, z, 0, 1) ;
    }
    mri_dil = MRIdilate(mri, NULL) ;
    mri_border = MRImarkLabelBorderVoxels(mri2, NULL, target_label, 2, 1) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    {
      MRIwrite(mri_dil, "d.mgz") ;
      MRIwrite(mri_border, "b.mgz") ;
      MRIwrite(mri, "c.mgz") ;
    }

    // to even out the training set, look for voxels that are on the border, neighboring
    // one of the changed voxels that was itself not changed
    nvox = vl_dif->nvox ;
    do
    {
      total_on = 0 ;
      for (i = 0 ; i < vl_dif->nvox ; i++)
      {
        x = vl_dif->xi[i] ;
        y = vl_dif->yi[i] ;
        z = vl_dif->zi[i] ;
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        found = 0 ;
        for (xk = -1 ; !found && xk <= 1 ; xk++)
          for (yk = -1 ; !found && yk <= 1 ; yk++)
            for (zk = -1 ; !found && zk <= 1 ; zk++)
            {
              total = fabs(xk)+fabs(yk)+fabs(xk) ;
              if (total != 1)
                continue ; // only 6-connected nbrs
              xi = mri1->xi[xk+x] ;  yi = mri1->yi[yk+y] ; zi = mri1->zi[zk+z] ; 
              nbr = (int)MRIgetVoxVal(mri_dil, xi, yi, zi, 0) ;
              if (nbr == 0)
                continue ;
              in_set_already = (int)MRIgetVoxVal(mri, xi, yi, zi, 0) ;
              if (in_set_already)
                continue ;
              border = MRIgetVoxVal(mri_border, xi, yi, zi, 0) ;
              if (border == 0)
                continue ;
              l1 = (int)MRIgetVoxVal(mri1, xi, yi, zi, 0) ;
              if (l1 != source_label)
                continue ;
              found = 1 ;
              MRIsetVoxVal(mri, xi, yi, zi, 0, 2) ;
              if (xi == Gx && yi == Gy && zi == Gz) 
                DiagBreak() ;
              if (nvox < vl_total->nvox)
              {
                l2 = (int)MRIgetVoxVal(mri2, xi, yi, zi, 0) ;
                vl_total->xi[nvox] = xi ;
                vl_total->yi[nvox] = yi ;
                vl_total->zi[nvox] = zi ;
                vl_total->vsrc[nvox] = l1 ;
                vl_total->vdst[nvox] = l2 ;
                nvox++ ;
                non++ ;
                total_on++;
              }
            }
        if (found == 0)
          DiagBreak() ;
        if (non >= noff)
          break ;
      }
      if (total_on == 0)
        break;
    } while (non < noff) ;
    if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      MRIwrite(mri, "c2.mgz") ;
  }
  if (mri == NULL)
    ErrorExit(ERROR_BADPARM, "vlst_add_border_voxels: couldn't find any voxels");
  MRIfree(&mri) ; MRIfree(&mri_dil) ; MRIfree(&mri_border) ;

  for (i = non = noff = 0 ; i < vl_total->nvox ; i++)
  {
    if (vl_total->vdst[i] == target_label)
      non++ ;
    else
    {
      if (vl_total->vsrc[i] != target_label)
        DiagBreak() ;
      noff++ ;
    }
  }

  return(vl_total) ;
}
/*
  compute an RAS coordinate system that is most closely aligned with the eigenvectors in
  m_evectors.
*/
static MATRIX *
compute_ras_basis_vectors(MRI *mri_aseg_orig, MRI *mri_aseg_edit, int target_label, 
                          MATRIX *m_evectors, double *means, int width, int height, int depth, 
                          int pad)
{
  MRI_REGION box1, box2, box ;
  MRI        *mri_aligned, *mri_tmp ;
  MATRIX     *m_xform, *m_trans, *m_tmp, *m_id, *m_inv, *m_targ, *m_src, *m_tmp2 ;
  int        i, j, max_col ;
  double     dot, max_dot ;

  m_xform = MatrixIdentity(4, NULL) ;
  m_trans = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(m_trans, 1, 4) = -means[0] ; *MATRIX_RELT(m_trans, 2, 4) = -means[1] ;
  *MATRIX_RELT(m_trans, 3, 4) = -means[2] ;
  m_inv = MatrixInverse(m_evectors, NULL) ;

  // rearrange vectors to be as close to RAS as possible
  m_id = MatrixIdentity(3, NULL) ;
  m_tmp = MatrixMultiply(m_id, m_inv, NULL) ;  // matrix of dot products
  for (i = 1 ; i <= 3 ; i++)
  {
    // find col in the ithrow that is max
    max_col = 0 ; max_dot = 0 ;
    for (j = 1 ; j <= 3 ; j++)
    {
      dot = *MATRIX_RELT(m_tmp, i, j) ;
      if (fabs(dot) > fabs(max_dot))
      {
        max_dot = dot ;
        max_col = j ;
      }
    }
    MatrixCopyRegion(m_inv, m_xform, 1, max_col, 3, 1, 1, i) ;
    if (max_dot < 0) // reverse eigenvector
      for (j = 1 ; j <= 3 ; j++)
        *MATRIX_RELT(m_xform, j, i) *= -1 ;
  }
  MatrixFree(&m_tmp) ; MatrixFree(&m_inv) ;
  m_tmp = MatrixMultiply(m_xform, m_trans, NULL) ;
  *MATRIX_RELT(m_trans, 1, 4) = means[0] ; *MATRIX_RELT(m_trans, 2, 4) = means[1] ;
  *MATRIX_RELT(m_trans, 3, 4) = means[2] ;
  MatrixMultiply(m_trans, m_tmp, m_xform) ;

  // now compute a transform that takes the bounding box to the desired width/height/depth
  mri_tmp = MRIlinearTransformInterp(mri_aseg_edit, NULL, m_xform, SAMPLE_NEAREST) ;
  mri_aligned = MRIdilateLabel(mri_tmp, NULL, target_label, pad) ;
  MRIlabelBoundingBox(mri_aligned, target_label, &box1) ;
  MRIfree(&mri_aligned) ; MRIfree(&mri_tmp) ;
  mri_tmp = MRIlinearTransformInterp(mri_aseg_orig, NULL, m_xform, SAMPLE_NEAREST) ;
  mri_aligned = MRIdilateLabel(mri_tmp, NULL, target_label, pad) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
  {
    MRIwrite(mri_aligned, "a.mgz") ;
  }
  MRIlabelBoundingBox(mri_aligned, target_label, &box2) ;
  MRIfree(&mri_aligned) ; MRIfree(&mri_tmp) ;

  box.x = MIN(box1.x, box2.x) ; box.y = MIN(box1.y, box2.y) ; box.z = MIN(box1.z, box2.z) ;
  box.dx = MAX(box1.dx+box1.x, box2.dx+box2.x) - box.x ;
  box.dy = MAX(box1.dy+box1.y, box2.dy+box2.y) - box.y ;
  box.dz = MAX(box1.dz+box1.z, box2.dz+box2.z) - box.z ;

  // now compute transform that takes corners of bounding box to corners of atlas
  m_targ = MatrixAlloc(4, 5, MATRIX_REAL) ;
  m_src = MatrixAlloc(4, 5, MATRIX_REAL) ;
  *MATRIX_RELT(m_targ, 1, 1) = 0 ; *MATRIX_RELT(m_targ, 2, 1) = 0 ; *MATRIX_RELT(m_targ, 3, 1) = 0 ; 
  *MATRIX_RELT(m_src, 1, 1) = box.x ; *MATRIX_RELT(m_src, 2, 1) = box.y ; *MATRIX_RELT(m_src, 3, 1) = box.z ; 

  *MATRIX_RELT(m_targ, 1, 2) = 0 ; *MATRIX_RELT(m_targ, 2, 2) = 0 ; *MATRIX_RELT(m_targ, 3, 2) = depth ; 
  *MATRIX_RELT(m_src, 1, 2) = box.x ; *MATRIX_RELT(m_src, 2, 2) = box.y ; 
  *MATRIX_RELT(m_src, 3, 2) = box.z+box.dz ; 

  *MATRIX_RELT(m_targ, 1, 3) = 0 ; *MATRIX_RELT(m_targ, 2, 3) = height ; *MATRIX_RELT(m_targ, 3, 3) = 0 ; 
  *MATRIX_RELT(m_src, 1, 3) = box.x ; *MATRIX_RELT(m_src, 2, 3) = box.y+box.dy ; 
  *MATRIX_RELT(m_src, 3, 3) = box.z ;

  *MATRIX_RELT(m_targ, 1, 4) = width ; *MATRIX_RELT(m_targ, 2, 4) = 0 ; *MATRIX_RELT(m_targ, 3, 4) = 0 ; 
  *MATRIX_RELT(m_src, 1, 4) = box.x+box.dx ; *MATRIX_RELT(m_src, 2, 4) = box.y ; 
  *MATRIX_RELT(m_src, 3, 4) = box.z ;

  *MATRIX_RELT(m_targ, 1, 5) = width ; *MATRIX_RELT(m_targ, 2, 5) = height ; *MATRIX_RELT(m_targ, 3, 5) = depth ; 
  *MATRIX_RELT(m_src, 1, 5) = box.x+box.dx ; *MATRIX_RELT(m_src, 2, 5) = box.y+box.dy ; 
  *MATRIX_RELT(m_src, 3, 5) = box.z+box.dz ;

  for (i = 1 ; i <= m_src->cols ; i++)
  {
    *MATRIX_RELT(m_src, 4, i) = *MATRIX_RELT(m_targ, 4,i) = 1.0 ;
  }

  m_inv = MatrixSVDPseudoInverse(m_src, NULL) ;
  m_tmp = MatrixMultiply(m_targ, m_inv, NULL) ;
  m_tmp2 = MatrixMultiply(m_tmp, m_xform, NULL) ;
  MatrixCopy(m_tmp2, m_xform) ;

  MatrixFree(&m_targ) ; MatrixFree(&m_src) ; MatrixFree(&m_inv) ; MatrixFree(&m_tmp) ; MatrixFree(&m_tmp2) ;
  MatrixFree(&m_trans) ; 
  return(m_xform) ;
}
