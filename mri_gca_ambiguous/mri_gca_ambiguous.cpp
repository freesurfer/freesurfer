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

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "cma.h"
#include "proto.h"
#include "mri.h"
#include "macros.h"
#include "version.h"
#include "gca.h"
#include "flash.h"


static double scale = 1e12 ;

int main(int argc, char *argv[]) ;

int GCAscale(GCA *gca_flash, double min_val, double max_val) ;
static  MRI *GCAcomputeAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) ;
static int  get_option(int argc, char *argv[]) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;
static  MRI *GCAcompute1DAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) ;
static  MRI *GCAcompute2DAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) ;
static  MRI *GCAcompute3DAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) ;

const char *Progname ;

static int optimize = 0 ;
static double lambda = 100.0 ;

static double TR = 20;
static double TE  = 3  ;
static double MIN_FA1 = 1  ;
static  double MAX_FA1 = 40 ;
static double MIN_FA2 = 1  ;
static  double MAX_FA2 = 40 ;
static double MIN_FA3 = 1  ;
static  double MAX_FA3 = 40 ;
static  double FA_STEP = 1 ;
static int append = 0 ;
static const char *fname = "amb.log" ;
static int left = 0 ;
static int classnum = -1 ;

#define MAX_FAS 10
int
main(int argc, char *argv[]) {
  char         **av, *out_name, *gca_name ;
  int          ac, nargs ;
  GCA          *gca  ;
  MRI          *mri = NULL ;
  double       TEs[MAX_FAS],  TRs[MAX_FAS],  FAs[MAX_FAS], amb ;
  GCA          *gca_flash ;

  nargs = handleVersionOption(argc, argv, "mri_gca_ambiguous");
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

  if (argc < 2)
    print_help() ;

  gca_name = argv[1] ;
  out_name = argv[2] ;

  printf("reading gca  from %s....\n", gca_name) ;
  gca = GCAread(gca_name) ;
  if (!gca)
    ErrorExit(ERROR_NOFILE, "%s: could not read gca file %s", Progname, gca_name) ;


  if (left)
    GCAreplaceRightWithLeft(gca) ;

  switch (optimize) {
  case 1: {
    double fa1, min_fa1  ;
    double min_amb, scale ;
    FILE   *fp ;

    fp = fopen(fname, append ? "a" : "w") ;

    TRs[0]  = TR ;
    TEs[0]  = TE ;
    min_amb = 10000000  ;
    min_fa1 = -1 ;
    scale = 1.0/(gca->prior_width * gca->prior_height * gca->prior_depth);
    for (fa1 = MIN_FA1  ; fa1 <= MAX_FA1 ; fa1 += FA_STEP) {
      FAs[0]  = RADIANS(fa1)  ;
      printf("testing flip angle %2.1f\n", fa1) ;
      gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 1, lambda);
      /*    GCAscale(gca_flash, 0, 75) ;*/
      mri =  GCAcomputeAmbiguity(gca_flash, mri, &amb, classnum) ;
      amb *= scale ;
      printf("\tflip angle %2.3f: ambiguity == %f\n", fa1, amb) ;
      if (amb  < min_amb || min_fa1  < 0) {
        min_amb = amb ;
        min_fa1 = fa1 ;
        printf("minimimum ambiguity at fa %2.3f (%f)\n", min_fa1, min_amb) ;
      }
      fprintf(fp, "%f %f\n", fa1, amb) ;
      fflush(fp) ;
      GCAfree(&gca_flash)  ;
    }

    MRIfree(&mri) ;
    fclose(fp) ;
    printf("minimimum ambiguity at  fas  %2.3f (%2.5f)\n", min_fa1, min_amb) ;
    FAs[0] = min_fa1 ;
    gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 1, lambda);
    mri =  GCAcomputeAmbiguity(gca_flash, NULL, &amb, classnum) ;
    MRIwrite(mri, out_name) ;
    MRIfree(&mri) ;
#if 0
    mri = MRIallocSequence(gca_flash->prior_width, gca_flash->prior_height,
                           gca_flash->prior_depth, MRI_FLOAT, gca_flash->ninputs) ;
    MRIsetResolution(mri, gca_flash->prior_spacing, gca_flash->prior_spacing,gca_flash->prior_spacing);
#endif
    mri = GCAbuildMostLikelyVolume(gca_flash,NULL) ;
    {
      int i ;
      char fname[STRLEN] ;
      for (i = 0 ; i < mri->nframes ; i++) {
        sprintf(fname, "mri_gca%d.mgh", i) ;
        MRIwriteFrame(mri, fname, i) ;
      }
    }
    GCAfree(&gca) ;
    GCAfree(&gca_flash) ;
    MRIfree(&mri) ;
  }
  break ;
  case 2: {
    double fa1,  fa2, min_fa1,  min_fa2  ;
    double min_amb, scale ;
    FILE   *fp ;

    fp = fopen(fname, append ? "a" : "w") ;

    TRs[0]  = TR ;
    TRs[1] = TR ;
    TEs[0]  = TE ;
    TEs[1] = TE ;
    min_amb = 10000000  ;
    min_fa1 = min_fa2 = -1 ;
    scale = 1.0/(gca->prior_width * gca->prior_height * gca->prior_depth);
    for (fa1 = MIN_FA1  ; fa1 <= MAX_FA1 ; fa1 += FA_STEP) {
      FAs[0]  = RADIANS(fa1)  ;
      for (fa2 = MIN_FA2 ; fa2 <= MAX_FA2 ; fa2 += FA_STEP) {
        printf("testing flip angles %2.1f, %2.1f\n", fa1, fa2) ;
        FAs[1] = RADIANS(fa2) ;
        gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 2, lambda);
        /*    GCAscale(gca_flash, 0, 75) ;*/
        mri =  GCAcomputeAmbiguity(gca_flash, mri, &amb, classnum) ;
        amb *= scale ;
        printf("\tflip angles %2.3f, %2.3f: ambiguity == %f\n", fa1, fa2, amb) ;
        if (amb  < min_amb || min_fa1  < 0) {
          min_amb = amb ;
          min_fa1 = fa1 ;
          min_fa2 = fa2 ;
          printf("minimimum ambiguity at fas %2.3f, %2.3f (%f)\n", min_fa1,  min_fa2, min_amb) ;
        }
        fprintf(fp, "%f %f %f\n", fa1,fa2, amb) ;
        fflush(fp) ;
        GCAfree(&gca_flash)  ;
      }
    }

    MRIfree(&mri) ;
    fclose(fp) ;
    printf("minimimum ambiguity at  fas  %2.3f,  %2.3f (%2.5f)\n", min_fa1,  min_fa2, min_amb) ;
    FAs[0] = min_fa1 ;
    FAs[1] = min_fa2 ;
    gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 2, lambda);
    mri =  GCAcomputeAmbiguity(gca_flash, NULL, &amb, classnum) ;
    MRIwrite(mri, out_name) ;
    MRIfree(&mri) ;
#if 0
    mri = MRIallocSequence(gca_flash->prior_width, gca_flash->prior_height,
                           gca_flash->prior_depth, MRI_FLOAT, gca_flash->ninputs) ;
    MRIsetResolution(mri, gca_flash->prior_spacing, gca_flash->prior_spacing,gca_flash->prior_spacing);
#endif
    mri = GCAbuildMostLikelyVolume(gca_flash,NULL) ;
    {
      int i ;
      char fname[STRLEN] ;
      for (i = 0 ; i < mri->nframes ; i++) {
        sprintf(fname, "mri_gca%d.mgh", i) ;
        MRIwriteFrame(mri, fname, i) ;
      }
    }
    GCAfree(&gca) ;
    GCAfree(&gca_flash) ;
    MRIfree(&mri) ;
  }
  break ;
  case 3: {
    double fa1,  fa2, min_fa1,  min_fa2, fa3, min_fa3  ;
    double min_amb, scale ;
    FILE   *fp ;

    fp = fopen(fname, append ? "a" : "w") ;

    TRs[0]  = TR ;
    TRs[1] = TR ;
    TRs[2] = TR ;
    TEs[0]  = TE ;
    TEs[1] = TE ;
    TEs[2] = TE ;
    min_amb = 10000000  ;
    min_fa3 = min_fa1 = min_fa2 = -1 ;
    scale = 1.0/(gca->prior_width * gca->prior_height * gca->prior_depth);
    for (fa1 = MIN_FA1  ; fa1 <= MAX_FA1 ; fa1 += FA_STEP) {
      FAs[0]  = RADIANS(fa1)  ;
      for (fa2 = MIN_FA2 ; fa2 <= MAX_FA2 ; fa2 += FA_STEP) {
        FAs[1] = RADIANS(fa2) ;
        for (fa3 = MIN_FA3 ; fa3 <= MAX_FA3 ; fa3 += FA_STEP) {
          FAs[2] = RADIANS(fa3) ;
          printf("testing flip angles %2.1f, %2.1f %2.1f\n", fa1, fa2, fa3) ;
          gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 3, lambda);
          /*    GCAscale(gca_flash, 0, 75) ;*/
          mri =  GCAcomputeAmbiguity(gca_flash, mri, &amb, classnum) ;
          amb *= scale ;
          printf("\tflip angles %2.3f, %2.3f, %2.3f: ambiguity == %f\n",
                 fa1, fa2, fa3, amb) ;
          if (amb  < min_amb || min_fa1  < 0) {
            min_amb = amb ;
            min_fa1 = fa1 ;
            min_fa2 = fa2 ;
            min_fa3 = fa3 ;
            printf("minimimum ambiguity at fas %2.3f, %2.3f, %2.3f (%f)\n", min_fa1,  min_fa2, min_fa3, min_amb) ;
          }
          fprintf(fp, "%f %f %f %f\n", fa1,fa2, fa3,amb) ;
          fflush(fp) ;
          GCAfree(&gca_flash)  ;
        }
      }
    }
    MRIfree(&mri) ;
    fclose(fp) ;
    printf("minimimum ambiguity at  fas  %2.3f,  %2.3f,  %2.3f(%2.5f)\n", min_fa1,  min_fa2, min_fa3, min_amb) ;
    FAs[0] = min_fa1 ;
    FAs[1] = min_fa2 ;
    FAs[2] = min_fa3 ;
#if 0
    gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 3, lambda);
    mri =  GCAcomputeAmbiguity(gca_flash, NULL, &amb, classnum) ;
    MRIwrite(mri, out_name) ;
    MRIfree(&mri) ;
    mri = MRIallocSequence(gca_flash->prior_width, gca_flash->prior_height,
                           gca_flash->prior_depth, MRI_FLOAT, gca_flash->ninputs) ;
    MRIsetResolution(mri, gca_flash->prior_spacing, gca_flash->prior_spacing,gca_flash->prior_spacing);
    mri = GCAbuildMostLikelyVolume(gca_flash,NULL) ;
    {
      int i ;
      char fname[STRLEN] ;
      for (i = 0 ; i < mri->nframes ; i++) {
        sprintf(fname, "mri_gca%d.mgh", i) ;
        MRIwriteFrame(mri, fname, i) ;
      }
    }
    GCAfree(&gca) ;
    GCAfree(&gca_flash) ;
    MRIfree(&mri) ;
#endif
  }
  break ;
  case 0: {
    int i ;
    char fname[STRLEN], tmp[STRLEN] ;

    TRs[0]  = TR ;
    TRs[1] = TR ;
    TEs[0]  = TE ;
    TEs[1] = TE ;
    FAs[0]  = RADIANS(MIN_FA1)  ;
    FAs[1]  = RADIANS(MIN_FA2)  ;
    gca_flash = GCAcreateFlashGCAfromParameterGCA(gca, TRs, FAs, TEs, 2, lambda);
    mri =  GCAcomputeAmbiguity(gca_flash, NULL, &amb, classnum) ;
    printf("writing max ambiguity map to %s...\n", out_name) ;
    MRIwrite(mri, out_name) ;
    MRIfree(&mri) ;

    mri = GCAbuildMostLikelyVolume(gca,NULL) ;
    for (i = 0 ; i < mri->nframes ; i++) {
      sprintf(fname, "%s%d.mgh", FileNameRemoveExtension(out_name, tmp), i) ;
      printf("writing max likeliood interpretation of gca to %s...\n", fname) ;
      MRIwriteFrame(mri, fname, i) ;
    }
  }
  break ;
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
  else if (!stricmp(option, "left"))
    left = 1 ;
  else if (!stricmp(option, "debug_voxel")) {
    Gx = atoi(argv[2]) ;
    Gy = atoi(argv[3]) ;
    Gz = atoi(argv[4]) ;
    printf("debugging  voxel (%d, %d,  %d)\n", Gx, Gy,  Gz) ;
    nargs = 3 ;
  } else if (!stricmp(option, "FA1")) {
    MIN_FA1 = atof(argv[2]) ;
    MAX_FA1 = atof(argv[3]) ;
    printf("covering range for 1st flip angle %2.3f-->%2.3f\n", MIN_FA1, MAX_FA1) ;
    nargs = 2 ;
  } else if (!stricmp(option, "SCALE")) {
    scale = atof(argv[2]) ;
    printf("using scaling factor of %f\n", scale) ;
    nargs = 1 ;
  } else if (!stricmp(option, "FA2")) {
    MIN_FA2 = atof(argv[2]) ;
    MAX_FA2 = atof(argv[3]) ;
    printf("covering range for 2nd flip angle %2.3f-->%2.3f\n", MIN_FA2, MAX_FA2) ;
    nargs = 2 ;
  } else if (!stricmp(option, "FA3")) {
    MIN_FA3 = atof(argv[2]) ;
    MAX_FA3 = atof(argv[3]) ;
    printf("covering range for 3rd flip angle %2.3f-->%2.3f\n", MIN_FA3, MAX_FA3) ;
    nargs = 2 ;
  } else if (!stricmp(option, "LAMBDA")) {
    lambda = atof(argv[2]) ;
    printf("using 1/SNR lambda = %2.3f\n", lambda) ;
    nargs = 1 ;
  } else switch (toupper(*option)) {
    case 'W':
      append = 0 ;
      nargs = 1 ;
      fname = argv[2] ;
      printf("writing output to %s...\n", fname) ;
      break ;
    case 'A':
      append = 1 ;
      nargs = 1 ;
      fname = argv[2] ;
      printf("appending output to %s...\n", fname) ;
      break ;
    case 'C':
      classnum = atoi(argv[2]) ;
      printf("optimizing for classnum %s (%d)\n", cma_label_to_name(classnum), classnum) ;
      nargs = 1 ;
      break ;
    case 'S':
      FA_STEP = atof(argv[2]) ;
      printf("using step size %2.3f\n", FA_STEP) ;
      nargs = 1 ;
      break ;
    case 'O':
      optimize = atoi(argv[2])  ;
      nargs = 1  ;
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
print_usage(void) {
  fprintf(stderr,
          "usage: %s [options] <gca file> <output volume>\n",
          Progname) ;
}

static void
print_help(void) {
  print_usage() ;
  fprintf(stderr,
          "\nThis program computes an ambiguity measure across a  GCA and output an MR image of it\n") ;
  exit(1) ;
}


static void
print_version(void) {
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;

}


static  MRI *
GCAcomputeAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) {
  int        xp, yp, zp, l1, l2, i1,  i2,  i, label_count[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1], label1_count[MAX_CMA_LABEL+1] ;
  GCA_PRIOR  *gcap  ;
  double     ambiguity, atotal, p1,  p2, amax, total_ambiguity, min_I, max_I, std,  Istep,  I ;
  float      vals[MAX_GCA_INPUTS] ;
  double     label_ambiguity[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1];
  GC1D       *gc1, *gc2 ;

  memset(label_ambiguity, 0, sizeof(label_ambiguity)) ;
  memset(label_count, 0, sizeof(label_count)) ;
  memset(label1_count, 0, sizeof(label1_count)) ;
  if (!mri)
    mri =  MRIalloc(gca->prior_width,gca->prior_height, gca->prior_depth, MRI_FLOAT)  ;
  mri->xsize = mri->ysize = mri->zsize = gca->prior_spacing ;
  mri->x_r = -1;
  mri->y_r =  0;
  mri->z_r =  0;
  mri->c_r = 0.0;
  mri->x_a =  0;
  mri->y_a =  0;
  mri->z_a =  1;
  mri->c_a = 0.0;
  mri->x_s =  0;
  mri->y_s = -1;
  mri->z_s =  0;
  mri->c_s = 0.0;

  mri->ras_good_flag = 1 ;
  switch (gca->ninputs) {
  case 1:
    return(GCAcompute1DAmbiguity(gca, mri, pamb, classnum)) ;
    break ;
  case 2:
    return(GCAcompute2DAmbiguity(gca, mri, pamb, classnum)) ;
    break ;
  case 3:
    return(GCAcompute3DAmbiguity(gca, mri, pamb, classnum)) ;
    break ;
  default:
    break ;
  }

  for (total_ambiguity = xp = 0  ; xp  < gca->prior_width ; xp++) {
    for  (yp  = 0 ;  yp < gca->prior_height ; yp++) {
      for (zp = 0 ;  zp  < gca->prior_depth  ; zp++) {
        gcap =  &gca->priors[xp][yp][zp] ;
        if (xp == Gx && yp == Gy && zp == Gz)
          DiagBreak() ;
        if (gcap->nlabels <= 1)
          continue ;

        min_I  = 10000000 ;
        max_I = 0 ;
        for (i1 = 0 ; i1  <  gcap->nlabels ; i1++) {
          l1 = gcap->labels[i1] ;
          gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
          if (!gc1)
            continue ;
          std = sqrt(covariance_determinant(gc1, gca->ninputs)) ;
          for (i = 0  ; i < gca->ninputs ; i++) {
#define NSTDS 1
            if (gc1->means[i]+NSTDS*std > max_I)
              max_I = gc1->means[i]+NSTDS*std ;
            if (gc1->means[i]-NSTDS*std < min_I)
              min_I = gc1->means[i]-NSTDS*std ;
          }
        }

        for (amax =  atotal = 0.0, i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
          for (i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
            l1 = gcap->labels[i1] ;
            gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
            if (!gc1)
              continue ;
            label1_count[l1]++ ;
            for (i2 = 0 ;  i2 < i1 ; i2++) {
              l2 = gcap->labels[i2] ;
              gc2 = GCAfindPriorGC(gca, xp, yp, zp,  l2) ;
              if (!gc2)
                continue ;

              if (l1 == 0 || l2 == 0)
                continue ;
              if  ((l1 == 17 &&  l2  == 7) ||(l2 == 17 &&  l1  == 7))
                DiagBreak() ;

#define ISTEPS 100

              /* marginalize over intensity */
              Istep = (max_I - min_I) / ISTEPS ;
              ambiguity = 0.0 ;
              if (gca->ninputs == 1) {
                for (ambiguity = 0,  I = min_I;  I <= max_I  ; I += Istep) {
                  vals[0] =  I ;
                  p1 = scale*GCAcomputeConditionalDensity(gc1, vals, gca->ninputs, l1) ;
                  p2 = scale*GCAcomputeConditionalDensity(gc2, vals, gca->ninputs, l2) ;
                  ambiguity += (p1*p2*Istep) ;
                }
              } else if (gca->ninputs == 2) {
                double I1, I2, Istep1, Istep2, std1, std2, scale1, scale2, dsq, amb ;
                MATRIX *m1_cov_inv, *m2_cov_inv ;
                VECTOR *v1, *v2, *vtmp ;

                m1_cov_inv = load_inverse_covariance_matrix(gc1, NULL, gca->ninputs) ;
                m2_cov_inv = load_inverse_covariance_matrix(gc2, NULL, gca->ninputs) ;
                v1 = VectorAlloc(2, MATRIX_REAL) ;
                v2 = VectorAlloc(2, MATRIX_REAL) ;
                vtmp = VectorAlloc(2, MATRIX_REAL) ;
                scale1 = scale*(1.0 / (pow(2*M_PI,gca->ninputs/2.0)*sqrt(covariance_determinant(gc1, gca->ninputs))));
                scale2 = scale*(1.0 / (pow(2*M_PI,gca->ninputs/2.0)*sqrt(covariance_determinant(gc2, gca->ninputs))));

                if (l1 == Left_Hippocampus && l2 == Left_Amygdala)
                  DiagBreak() ;
                std1 = sqrt(gc1->covars[0]) ;
                std2 = sqrt(gc1->covars[2]) ;
#undef ISTEPS
#define ISTEPS 50
                Istep1 = std1/(ISTEPS/4) ;
                Istep2 = std2/(ISTEPS/4) ;
                for (ambiguity = 0,  I1 = gc1->means[0]-(ISTEPS/2)*Istep1;  I1 <= gc1->means[0]+(ISTEPS/2)*Istep1  ; I1 += Istep1) {
                  VECTOR_ELT(v1, 1) = I1-gc1->means[0] ;
                  VECTOR_ELT(v2, 1) = I1-gc2->means[0] ;
                  for (amb = 0.0, I2 = gc1->means[1]-(ISTEPS/2)*Istep2;  I2 <= gc1->means[1]+(ISTEPS/2)*Istep2  ; I2 += Istep2) {
                    VECTOR_ELT(v1, 2) = I2-gc1->means[1] ;
                    MatrixMultiply(m1_cov_inv, v1, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
                    dsq = VectorDot(v1, vtmp) ;
                    p1 = scale1*exp(-0.5*dsq) ;

                    VECTOR_ELT(v2, 2) = I2-gc2->means[1] ;
                    MatrixMultiply(m2_cov_inv, v2, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
                    dsq = VectorDot(v2, vtmp) ;
                    p2 = scale2*exp(-0.5*dsq) ;
                    amb += (p1*p2*Istep2) ;
                  }
                  ambiguity += (amb*Istep1) ;
                }
                std1 = sqrt(gc2->covars[0]) ;
                std2 = sqrt(gc2->covars[2]) ;
                Istep1 = std1/(ISTEPS/4) ;
                Istep2 = std2/(ISTEPS/4) ;
                for (I1 = gc2->means[0]-(ISTEPS/2)*Istep1;  I1 <= gc2->means[0]+(ISTEPS/2)*Istep1  ; I1 += Istep1) {
                  VECTOR_ELT(v1, 1) = I1-gc1->means[0] ;
                  VECTOR_ELT(v2, 1) = I1-gc2->means[0] ;
                  for (amb = 0.0, I2 = gc2->means[1]-(ISTEPS/2)*Istep2;  I2 <= gc2->means[1]+(ISTEPS/2)*Istep2  ; I2 += Istep2) {
                    VECTOR_ELT(v1, 2) = I2-gc1->means[1] ;
                    MatrixMultiply(m1_cov_inv, v1, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
                    dsq = VectorDot(v1, vtmp) ;
                    p1 = scale1*exp(-0.5*dsq) ;

                    VECTOR_ELT(v2, 2) = I2-gc2->means[1] ;
                    MatrixMultiply(m2_cov_inv, v2, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
                    dsq = VectorDot(v2, vtmp) ;
                    p2 = scale2*exp(-0.5*dsq) ;
                    amb += (p1*p2*Istep2) ;
                  }
                  ambiguity += (amb*Istep1) ;
                }
                MatrixFree(&m1_cov_inv) ;
                MatrixFree(&m2_cov_inv) ;
                VectorFree(&v1) ;
                VectorFree(&v2) ;
                VectorFree(&vtmp) ;
              } else
                ErrorExit(ERROR_UNSUPPORTED, "%s: ambiguity  only supported  for  ninputs=1 or 2", Progname)  ;

              label_ambiguity[l1][l2] += ambiguity ;
              label_ambiguity[l2][l1] += ambiguity ;
              label_count[l1][l2]++;
              label_count[l2][l1]++;
              ambiguity *= ambiguity * (gcap->priors[i1]+gcap->priors[i2]) ;
              atotal += ambiguity ;
              if (ambiguity > amax)
                amax  =  ambiguity ;
            }
          }
        }
        total_ambiguity += atotal ;
        if (amax > 0)
          MRIFvox(mri, xp, yp,  zp) = amax  ;
      }
    }
  }

#if 0
  {
    FILE *fp ;
    MRI  *mri2;
    float norm;

    mri2 = MRIalloc(MAX_CMA_LABEL+1, MAX_CMA_LABEL+1, 2, MRI_FLOAT) ;
    fp = fopen("label_amb.dat", "w") ;
    for  (l1 = 0 ; l1 <= MAX_CMA_LABEL ; l1++) {
      for (l2 =  0  ; l2  <=  MAX_CMA_LABEL  ; l2++) {
#if  1
        norm = label_count[l1][l2];
#else
        norm = (label1_count[l1]+label1_count[l2])/2 ;
#endif
        if  (norm> 0)
          label_ambiguity[l1][l2] /= (float)norm ;
        fprintf(fp, "%f   ", label_ambiguity[l1][l2])  ;
        MRIFvox(mri2, l1, l2, 0) = label_ambiguity[l1][l2] ;
      }
      fprintf(fp, "\n") ;
    }
    fclose(fp)  ;
    MRIwrite(mri2, "label_amb.mgh") ;
    MRIfree(&mri2);
  }
#endif

  *pamb = total_ambiguity;
  return(mri) ;
}

int
GCAscale(GCA *gca, double min_val, double max_val) {
  int        xp, yp, zp,i, n, r, c, k ;
  double     scale, umin[MAX_GCA_INPUTS], umax[MAX_GCA_INPUTS], range ;
  GCA_PRIOR  *gcap  ;
  GC1D       *gc ;

  for (xp = 0  ; xp  < gca->prior_width ; xp++) {
    for  (yp  = 0 ;  yp < gca->prior_height ; yp++) {
      for (zp = 0 ;  zp  < gca->prior_depth  ; zp++) {
        gcap =  &gca->priors[xp][yp][zp] ;
        if (xp == Gx && yp == Gy && zp == Gz)
          DiagBreak() ;
        if (gcap->nlabels <= 1)
          continue ;

        for (i = 0 ; i  <  gcap->nlabels ; i++) {
          gc = GCAfindPriorGC(gca, xp, yp, zp,  gcap->labels[i]) ;
          for (n = 0 ; n < gca->ninputs ; n++) {
            if (i == 0 && xp == 0 && yp == 0 && zp == 0) {
              umin[n] = gc->means[n] ;
              umax[n] = gc->means[n] ;
            } else {
              if (umin[n] > gc->means[n])
                umin[n] = gc->means[n] ;
              if (umax[n] < gc->means[n])
                umax[n] = gc->means[n] ;
            }
          }
        }
      }
    }
  }

  printf("ranges:\n") ;
  range = 0 ;
  for (n = 0 ; n < gca->ninputs ; n++) {
    if ((umax[n] - umin[n]) > range)
      range = umax[n] - umin[n] ;
    printf("%d: %f --> %f\n", n, umin[n], umax[n]) ;
  }
  scale = (max_val - min_val) / range ;
  printf("scaling by %2.3f\n", scale) ;

  for (xp = 0  ; xp  < gca->prior_width ; xp++) {
    for  (yp  = 0 ;  yp < gca->prior_height ; yp++) {
      for (zp = 0 ;  zp  < gca->prior_depth  ; zp++) {
        gcap =  &gca->priors[xp][yp][zp] ;
        if (xp == Gx && yp == Gy && zp == Gz)
          DiagBreak() ;
        if (gcap->nlabels <= 1)
          continue ;

        for (i = 0 ; i  <  gcap->nlabels ; i++) {
          gc = GCAfindPriorGC(gca, xp, yp, zp,  gcap->labels[i]) ;
          for (n = 0 ; n < gca->ninputs ; n++) {
            gc->means[n] = (gc->means[n] - umin[n])*scale ;
            for (r = k = 0 ; r < gca->ninputs ; r++)
              for (c = r ;  c < gca->ninputs ; c++, k++)
                gc->covars[k] *= (scale*scale) ;
          }
        }
      }
    }
  }


  return(NO_ERROR) ;
}

#ifdef NSTDS
#undef NSTDS
#endif
#define NSTDS             2.5     /* go out 3 stds on either side of the mean */
#define SAMPLING_DENSITY  0.1   /* sample every 1/10th of a std */

#define MAX_VALS         ((int)((2*NSTDS/SAMPLING_DENSITY)+0.9)+1)
#define MAX_VAL          (MAX_VALS-1)

#define VAL_TO_INDEX(v,mean,std)  ((int)(MAX_VAL*((v-mean)/(2*NSTDS*std))))
#define MAX_LABELS_PER_GCAN    10
static  MRI *
GCAcompute1DAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) {
  int        xp, yp, zp, l1, l2, i1,  i2, label_count[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1], label1_count[MAX_CMA_LABEL+1],
  index1, index1_c1, index1_c2, n,ntotal ;
  GCA_PRIOR  *gcap  ;
  double     ambiguity, atotal, p1,  p2, amax, total_ambiguity, Imin1, Imax1, pscale, dsq, std1,
  Istep1,  I1, Imins1[MAX_LABELS_PER_GCAN],
  Imaxs1[MAX_LABELS_PER_GCAN] ;
  float      pI_c[MAX_LABELS_PER_GCAN][MAX_VALS] ;
  double     label_ambiguity[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1];
  GC1D       *gc1, *gc2 ;
  MATRIX     *m_cov_inv = NULL ;
  VECTOR     *v = NULL, *vtmp = NULL ;

  memset(label_ambiguity, 0, sizeof(label_ambiguity)) ;
  memset(label_count, 0, sizeof(label_count)) ;
  memset(label1_count, 0, sizeof(label1_count)) ;
  if (!mri)
    mri =  MRIalloc(gca->prior_width,gca->prior_height, gca->prior_depth, MRI_FLOAT)  ;
  mri->xsize = mri->ysize = mri->zsize = gca->prior_spacing ;
  mri->x_r = -1;
  mri->y_r =  0;
  mri->z_r =  0;
  mri->c_r = 0.0;
  mri->x_a =  0;
  mri->y_a =  0;
  mri->z_a =  1;
  mri->c_a = 0.0;
  mri->x_s =  0;
  mri->y_s = -1;
  mri->z_s =  0;
  mri->c_s = 0.0;

  mri->ras_good_flag = 1 ;

  for (ntotal = 0, total_ambiguity = xp = 0  ; xp  < gca->prior_width ; xp++) {
    for  (yp  = 0 ;  yp < gca->prior_height ; yp++) {
      for (zp = 0 ;  zp  < gca->prior_depth  ; zp++) {
        gcap =  &gca->priors[xp][yp][zp] ;
        if (xp == Gx && yp == Gy && zp == Gz)
          DiagBreak() ;
        if (gcap->nlabels <= 1)   /* ambiguity is 0 */
          continue ;

        /* now fill in lookup tables with probabilities */
        for (i1 = 0 ; i1  <  gcap->nlabels ; i1++) {
          l1 = gcap->labels[i1] ;
          gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
          if (!gc1 || gc1->regularized)
            continue ;
          std1 = sqrt(gc1->covars[0]) ;
          Imin1 = gc1->means[0]-NSTDS*std1 ;
          Imax1 = gc1->means[0]+NSTDS*std1 ;
          Istep1 = (Imax1-Imin1)/MAX_VALS ;
          Imins1[i1] = Imin1 ;
          Imaxs1[i1] = Imax1 ;

          m_cov_inv = load_inverse_covariance_matrix(gc1, m_cov_inv, gca->ninputs) ;
          if (!v) {
            v = VectorAlloc(1, MATRIX_REAL) ;
            vtmp = VectorAlloc(1, MATRIX_REAL) ;
          }
          pscale = (1.0 / (pow(2*M_PI,gca->ninputs/2.0)*sqrt(covariance_determinant(gc1, gca->ninputs))));
          for (index1 = 0, I1 = Imin1 ; I1 <= Imax1 ; I1 += Istep1, index1++) {
            VECTOR_ELT(v, 1) = I1-gc1->means[0] ;
            MatrixMultiply(m_cov_inv, v, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
            dsq = VectorDot(v, vtmp) ;
            p1 = pscale*exp(-0.5*dsq) ;
            pI_c[i1][index1] = p1 ;
          }
          if (Gdiag & DIAG_VERBOSE) {
            FILE *fp ;
            char fname[STRLEN] ;
            sprintf(fname, "p%d.dat", i1) ;
            fp = fopen(fname, "w") ;
            for (I1 = Imin1 ; I1 <= Imax1 ; I1 += Istep1) {
              index1 = nint(MAX_VAL * (I1 - Imin1) / (Imax1-Imin1)) ;
              VECTOR_ELT(v, 1) = I1-gc1->means[0] ;
              fprintf(fp, "%d %f %f\n", index1, I1, pI_c[i1][index1]) ;
            }
            fclose(fp) ;
          }
        }

        for (n = 0, amax =  atotal = 0.0, i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
          for (i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
            l1 = gcap->labels[i1] ;
            gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
            if (!gc1 || gc1->regularized)
              continue ;
            label1_count[l1]++ ;
            for (i2 = 0 ;  i2 < i1 ; i2++) {
              l2 = gcap->labels[i2] ;
              gc2 = GCAfindPriorGC(gca, xp, yp, zp,  l2) ;
              if (!gc2 || gc2->regularized)
                continue ;
              if (l1 == l2)
                continue ;
              if (classnum > 0 && l1 != classnum && l2 != classnum)
                continue ;

              n++ ;
              if  ((l1 == 17 &&  l2  == 7) ||(l2 == 17 &&  l1  == 7))
                DiagBreak() ;

              /* marginalize over intensity */
              if (l1 == Left_Hippocampus && l2 == Left_Amygdala)
                DiagBreak() ;
              Imin1 = MAX(Imins1[i1], Imins1[i2]) ;
              Imax1 = MIN(Imaxs1[i1], Imaxs1[i2]) ;
              Istep1 = (Imax1-Imin1)/MAX_VALS ;
              for (ambiguity = 0,  I1 = Imin1 ; I1 < Imax1 ; I1 += Istep1) {
                index1_c1 = nint(MAX_VAL * (I1 - Imins1[i1]) / (Imaxs1[i1]-Imins1[i1])) ;
                if (index1_c1 < 0 || index1_c1 > MAX_VAL)
                  continue ;
                index1_c2 = nint(MAX_VAL * (I1 - Imins1[i2]) / (Imaxs1[i2]-Imins1[i2])) ;
                if (index1_c2 < 0 || index1_c2 > MAX_VAL)
                  continue ;

                p1 = pI_c[i1][index1_c1] ;
                p2 = pI_c[i2][index1_c2] ;
                ambiguity += (p1*p2*Istep1) ;
              }

              ambiguity *= scale*0.5*(gcap->priors[i1]+gcap->priors[i2]) ;
              atotal += ambiguity ;
              label_ambiguity[l1][l2] += ambiguity ;
              label_ambiguity[l2][l1] += ambiguity ;
              label_count[l1][l2]++;
              label_count[l2][l1]++;
              if (ambiguity > amax)
                amax  =  ambiguity ;
            }
          }
        }
        if ( n > 0)
          total_ambiguity += (atotal/n) ;
        if (amax > 0)
          MRIFvox(mri, xp, yp,  zp) = amax  ;
        ntotal += n ;
      }
    }
  }

  {
    FILE *fp ;
    MRI  *mri2;
    float norm;

    mri2 = MRIalloc(MAX_CMA_LABEL+1, MAX_CMA_LABEL+1, 2, MRI_FLOAT) ;
    fp = fopen("label_amb.dat", "w") ;
    for  (l1 = 0 ; l1 <= MAX_CMA_LABEL ; l1++) {
      for (l2 =  0  ; l2  <=  MAX_CMA_LABEL  ; l2++) {
#if 1
        norm = label_count[l1][l2];
#else
        norm = (label1_count[l1]+label1_count[l2])/2 ;
#endif
        if  (norm> 0)
          label_ambiguity[l1][l2] /= (float)norm ;
        fprintf(fp, "%f   ", label_ambiguity[l1][l2])  ;
        MRIFvox(mri2, l1, l2, 0) = label_ambiguity[l1][l2] ;
      }
      fprintf(fp, "\n") ;
    }
    fclose(fp)  ;
    MRIwrite(mri2, "label_amb.mgh") ;
    MRIfree(&mri2);
  }

  if (ntotal == 0) {
    total_ambiguity = scale ;
    printf("**** no valid nodes found ********* \n") ;
  } else {
    total_ambiguity /= ntotal ;
    printf("%d total nodes used in ambiguity calculation...\n", ntotal) ;
    MatrixFree(&m_cov_inv) ;
    VectorFree(&v) ;
    VectorFree(&vtmp) ;
  }
  *pamb = total_ambiguity;
  return(mri) ;
}
static MRI *
GCAcompute2DAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) {
  int        xp, yp, zp, l1, l2, i1,  i2, label_count[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1], label1_count[MAX_CMA_LABEL+1],
  index1, index2, index1_c1, index1_c2, index2_c1, index2_c2, n,ntotal ;
  GCA_PRIOR  *gcap  ;
  double     ambiguity, atotal, p1,  p2, amax, total_ambiguity, Imin1, Imin2, Imax1, Imax2, pscale, dsq, amb, std1, std2,
  Istep1, Istep2,  I1, I2, Imins1[MAX_LABELS_PER_GCAN], Imins2[MAX_LABELS_PER_GCAN],
  Imaxs1[MAX_LABELS_PER_GCAN], Imaxs2[MAX_LABELS_PER_GCAN] ;
  float      pI_c[MAX_LABELS_PER_GCAN][MAX_VALS][MAX_VALS] ;
  double     label_ambiguity[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1];
  GC1D       *gc1, *gc2 ;
  MATRIX     *m_cov_inv = NULL ;
  VECTOR     *v = NULL, *vtmp = NULL ;

  memset(label_ambiguity, 0, sizeof(label_ambiguity)) ;
  memset(label_count, 0, sizeof(label_count)) ;
  memset(label1_count, 0, sizeof(label1_count)) ;
  if (!mri)
    mri =  MRIalloc(gca->prior_width,gca->prior_height, gca->prior_depth, MRI_FLOAT)  ;
  mri->xsize = mri->ysize = mri->zsize = gca->prior_spacing ;
  mri->x_r = -1;
  mri->y_r =  0;
  mri->z_r =  0;
  mri->c_r = 0.0;
  mri->x_a =  0;
  mri->y_a =  0;
  mri->z_a =  1;
  mri->c_a = 0.0;
  mri->x_s =  0;
  mri->y_s = -1;
  mri->z_s =  0;
  mri->c_s = 0.0;

  mri->ras_good_flag = 1 ;

  for (ntotal = 0, total_ambiguity = xp = 0  ; xp  < gca->prior_width ; xp++) {
    for  (yp  = 0 ;  yp < gca->prior_height ; yp++) {
      for (zp = 0 ;  zp  < gca->prior_depth  ; zp++) {
        gcap =  &gca->priors[xp][yp][zp] ;
        if (xp == Gx && yp == Gy && zp == Gz)
          DiagBreak() ;
        if (gcap->nlabels <= 1)   /* ambiguity is 0 */
          continue ;

        /* now fill in lookup tables with probabilities */
        for (i1 = 0 ; i1  <  gcap->nlabels ; i1++) {
          l1 = gcap->labels[i1] ;
          gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
          if (!gc1 || gc1->regularized)
            continue ;
          std1 = sqrt(gc1->covars[0]) ;
          std2 = sqrt(gc1->covars[2]) ;
          Imin1 = gc1->means[0]-NSTDS*std1 ;
          Imax1 = gc1->means[0]+NSTDS*std1 ;
          Istep1 = (Imax1-Imin1)/MAX_VALS ;
          Imins1[i1] = Imin1 ;
          Imaxs1[i1] = Imax1 ;

          Imin2 = gc1->means[1]-NSTDS*std2 ;
          Imax2 = gc1->means[1]+NSTDS*std2 ;
          Istep2 = (Imax2-Imin2)/MAX_VALS ;
          Imins2[i1] = Imin2 ;
          Imaxs2[i1] = Imax2 ;

          m_cov_inv = load_inverse_covariance_matrix(gc1, m_cov_inv, gca->ninputs) ;
          if (!v) {
            v = VectorAlloc(2, MATRIX_REAL) ;
            vtmp = VectorAlloc(2, MATRIX_REAL) ;
          }
          pscale = (1.0 / (pow(2*M_PI,gca->ninputs/2.0)*sqrt(covariance_determinant(gc1, gca->ninputs))));
          for (index1 = 0, I1 = Imin1 ; I1 <= Imax1 ; I1 += Istep1, index1++) {
            VECTOR_ELT(v, 1) = I1-gc1->means[0] ;
            for (index2 = 0, I2 = Imin2 ; I2 <= Imax2 ; I2 += Istep2, index2++) {
              VECTOR_ELT(v, 2) = I2-gc1->means[1] ;
              MatrixMultiply(m_cov_inv, v, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
              dsq = VectorDot(v, vtmp) ;
              p1 = pscale*exp(-0.5*dsq) ;
              pI_c[i1][index1][index2] = p1 ;
            }
          }
          if (Gdiag & DIAG_VERBOSE) {
            FILE *fp ;
            char fname[STRLEN] ;
            sprintf(fname, "p%d.dat", i1) ;
            fp = fopen(fname, "w") ;
            for (I1 = Imin1 ; I1 <= Imax1 ; I1 += Istep1) {
              index1 = nint(MAX_VAL * (I1 - Imin1) / (Imax1-Imin1)) ;
              VECTOR_ELT(v, 1) = I1-gc1->means[0] ;
              for (I2 = Imin2 ; I2 < Imax2 ; I2 += Istep2) {
                index2 = nint(MAX_VAL*(I2 - Imin2) / (Imax2-Imin2)) ;
                fprintf(fp, "%d %d %f %f %f\n", index1, index2, I1, I2, pI_c[i1][index1][index2]) ;
              }
            }
            fclose(fp) ;
          }
        }

        for (n = 0, amax =  atotal = 0.0, i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
          for (i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
            l1 = gcap->labels[i1] ;
            gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
            if (!gc1 || gc1->regularized)
              continue ;
            label1_count[l1]++ ;
            for (i2 = 0 ;  i2 < i1 ; i2++) {
              l2 = gcap->labels[i2] ;
              gc2 = GCAfindPriorGC(gca, xp, yp, zp,  l2) ;
              if (!gc2 || gc2->regularized)
                continue ;
              if (l1 == l2)
                continue ;
              if (classnum > 0 && l1 != classnum && l2 != classnum)
                continue ;

              n++ ;
              if  ((l1 == 17 &&  l2  == 7) ||(l2 == 17 &&  l1  == 7))
                DiagBreak() ;

              /* marginalize over intensity */
              if (l1 == Left_Hippocampus && l2 == Left_Amygdala)
                DiagBreak() ;
              Imin1 = MAX(Imins1[i1], Imins1[i2]) ;
              Imax1 = MIN(Imaxs1[i1], Imaxs1[i2]) ;
              Imin2 = MAX(Imins2[i1], Imins2[i2]) ;
              Imax2 = MIN(Imaxs2[i1], Imaxs2[i2]) ;
              Istep1 = (Imax1-Imin1)/MAX_VALS ;
              Istep2 = (Imax2-Imin2)/MAX_VALS ;
              for (ambiguity = 0,  I1 = Imin1 ; I1 < Imax1 ; I1 += Istep1) {
                index1_c1 = nint(MAX_VAL * (I1 - Imins1[i1]) / (Imaxs1[i1]-Imins1[i1])) ;
                if (index1_c1 < 0 || index1_c1 > MAX_VAL)
                  continue ;
                index1_c2 = nint(MAX_VAL * (I1 - Imins1[i2]) / (Imaxs1[i2]-Imins1[i2])) ;
                if (index1_c2 < 0 || index1_c2 > MAX_VAL)
                  continue ;
                for (amb = 0.0, I2 = Imin2 ; I2 < Imax2 ; I2 += Istep2) {
                  index2_c1 = nint(MAX_VAL * (I2 - Imins2[i1]) / (Imaxs2[i1]-Imins2[i1])) ;
                  if (index2_c1 < 0 || index2_c1 > MAX_VAL)
                    continue ;
                  index2_c2 = nint(MAX_VAL * (I2 - Imins2[i2]) / (Imaxs2[i2]-Imins2[i2])) ;
                  if (index2_c2 < 0 || index2_c2 > MAX_VAL)
                    continue ;
                  p1 = pI_c[i1][index1_c1][index2_c1] ;
                  p2 = pI_c[i2][index1_c2][index2_c2] ;
                  amb += (p1*p2*Istep2) ;
                }
                ambiguity += (amb*Istep1) ;
              }

              ambiguity *= scale*0.5*(gcap->priors[i1]+gcap->priors[i2]) ;
              atotal += ambiguity ;
              label_ambiguity[l1][l2] += ambiguity ;
              label_ambiguity[l2][l1] += ambiguity ;
              label_count[l1][l2]++;
              label_count[l2][l1]++;
              if (ambiguity > amax)
                amax  =  ambiguity ;
            }
          }
        }
        if ( n > 0)
          total_ambiguity += (atotal/n) ;
        if (amax > 0)
          MRIFvox(mri, xp, yp,  zp) = amax  ;
        ntotal += n ;
      }
    }
  }

  {
    FILE *fp ;
    MRI  *mri2;
    float norm;

    mri2 = MRIalloc(MAX_CMA_LABEL+1, MAX_CMA_LABEL+1, 2, MRI_FLOAT) ;
    fp = fopen("label_amb.dat", "w") ;
    for  (l1 = 0 ; l1 <= MAX_CMA_LABEL ; l1++) {
      for (l2 =  0  ; l2  <=  MAX_CMA_LABEL  ; l2++) {
#if 1
        norm = label_count[l1][l2];
#else
        norm = (label1_count[l1]+label1_count[l2])/2 ;
#endif
        if  (norm> 0)
          label_ambiguity[l1][l2] /= (float)norm ;
        fprintf(fp, "%f   ", label_ambiguity[l1][l2])  ;
        MRIFvox(mri2, l1, l2, 0) = label_ambiguity[l1][l2] ;
      }
      fprintf(fp, "\n") ;
    }
    fclose(fp)  ;
    MRIwrite(mri2, "label_amb.mgh") ;
    MRIfree(&mri2);
  }

  if (ntotal == 0) {
    total_ambiguity = scale ;
    printf("**** no valid nodes found ********* \n") ;
  } else {
    total_ambiguity /= ntotal ;
    printf("%d total nodes used in ambiguity calculation...\n", ntotal) ;
    MatrixFree(&m_cov_inv) ;
    VectorFree(&v) ;
    VectorFree(&vtmp) ;
  }
  *pamb = total_ambiguity;
  return(mri) ;
}
static  MRI *
GCAcompute3DAmbiguity(GCA *gca, MRI *mri, double *pamb, int classnum) {
  int        xp, yp, zp, l1, l2, i1,  i2, label_count[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1], label1_count[MAX_CMA_LABEL+1],
  index1, index2, index3, index1_c1, index1_c2, index2_c1, index2_c2, index3_c1, index3_c2, n,ntotal ;
  GCA_PRIOR  *gcap  ;
  double     ambiguity, atotal, p1,  p2, amax, total_ambiguity, Imin1, Imin2, Imin3, Imax1, Imax2, Imax3, pscale, dsq, amb, std1, std2, std3,

  Istep1, Istep2,  Istep3, I1, I2, I3, Imins1[MAX_LABELS_PER_GCAN], Imins2[MAX_LABELS_PER_GCAN], Imins3[MAX_LABELS_PER_GCAN],
  Imaxs1[MAX_LABELS_PER_GCAN], Imaxs2[MAX_LABELS_PER_GCAN],
  Imaxs3[MAX_LABELS_PER_GCAN] ;
  float      pI_c[MAX_LABELS_PER_GCAN][MAX_VALS][MAX_VALS][MAX_VALS] ;
  double     label_ambiguity[MAX_CMA_LABEL+1][MAX_CMA_LABEL+1];
  GC1D       *gc1, *gc2 ;
  MATRIX     *m_cov_inv = NULL, *m_cov = NULL ;
  VECTOR     *v = NULL, *vtmp = NULL ;

  memset(label_ambiguity, 0, sizeof(label_ambiguity)) ;
  memset(label_count, 0, sizeof(label_count)) ;
  memset(label1_count, 0, sizeof(label1_count)) ;
  if (!mri)
    mri =  MRIalloc(gca->prior_width,gca->prior_height, gca->prior_depth, MRI_FLOAT)  ;
  mri->xsize = mri->ysize = mri->zsize = gca->prior_spacing ;
  mri->x_r = -1;
  mri->y_r =  0;
  mri->z_r =  0;
  mri->c_r = 0.0;
  mri->x_a =  0;
  mri->y_a =  0;
  mri->z_a =  1;
  mri->c_a = 0.0;
  mri->x_s =  0;
  mri->y_s = -1;
  mri->z_s =  0;
  mri->c_s = 0.0;

  mri->ras_good_flag = 1 ;

  for (ntotal = 0, total_ambiguity = xp = 0  ; xp  < gca->prior_width ; xp++) {
    for  (yp  = 0 ;  yp < gca->prior_height ; yp++) {
      for (zp = 0 ;  zp  < gca->prior_depth  ; zp++) {
        gcap =  &gca->priors[xp][yp][zp] ;
        if (xp == Gx && yp == Gy && zp == Gz)
          DiagBreak() ;
        if (gcap->nlabels <= 1)   /* ambiguity is 0 */
          continue ;

        /* now fill in lookup tables with probabilities */
        for (i1 = 0 ; i1  <  gcap->nlabels ; i1++) {

          l1 = gcap->labels[i1] ;
          gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
          if (!gc1 || gc1->regularized)
            continue ;
          m_cov = load_covariance_matrix(gc1, m_cov, gca->ninputs) ;

          std1 = sqrt(*MATRIX_RELT(m_cov,1,1)) ;
          std2 = sqrt(*MATRIX_RELT(m_cov,2,2)) ;
          std3 = sqrt(*MATRIX_RELT(m_cov,3,3)) ;
          Imin1 = gc1->means[0]-NSTDS*std1 ;
          Imax1 = gc1->means[0]+NSTDS*std1 ;
          Istep1 = (Imax1-Imin1)/MAX_VALS ;
          Imins1[i1] = Imin1 ;
          Imaxs1[i1] = Imax1 ;

          Imin2 = gc1->means[1]-NSTDS*std2 ;
          Imax2 = gc1->means[1]+NSTDS*std2 ;
          Istep2 = (Imax2-Imin2)/MAX_VALS ;
          Imins2[i1] = Imin2 ;
          Imaxs2[i1] = Imax2 ;

          Imin3 = gc1->means[2]-NSTDS*std3 ;
          Imax3 = gc1->means[2]+NSTDS*std3 ;
          Istep3 = (Imax3-Imin3)/MAX_VALS ;
          Imins3[i1] = Imin3 ;
          Imaxs3[i1] = Imax3 ;

          m_cov_inv = load_inverse_covariance_matrix(gc1, m_cov_inv, gca->ninputs) ;
          if (!v) {
            v = VectorAlloc(3, MATRIX_REAL) ;
            vtmp = VectorAlloc(3, MATRIX_REAL) ;
          }
          pscale = (1.0 / (pow(2*M_PI,gca->ninputs/2.0)*sqrt(covariance_determinant(gc1, gca->ninputs))));
          for (index1 = 0, I1 = Imin1 ; I1 <= Imax1 ; I1 += Istep1, index1++) {
            VECTOR_ELT(v, 1) = I1-gc1->means[0] ;
            for (index2 = 0, I2 = Imin2 ; I2 <= Imax2 ; I2 += Istep2, index2++) {
              VECTOR_ELT(v, 2) = I2-gc1->means[1] ;
              for (index3 = 0, I3 = Imin3 ; I3 <= Imax3 ; I3 += Istep3, index3++) {
                VECTOR_ELT(v, 3) = I3-gc1->means[2] ;
                MatrixMultiply(m_cov_inv, v, vtmp) ;  /* v_means is now inverse(cov) * v_vals */
                dsq = VectorDot(v, vtmp) ;
                p1 = pscale*exp(-0.5*dsq) ;
                pI_c[i1][index1][index2][index3] = p1 ;
              }
            }
          }
          if (Gdiag & DIAG_VERBOSE) {
            FILE *fp ;
            char fname[STRLEN] ;
            sprintf(fname, "p%d.dat", i1) ;
            fp = fopen(fname, "w") ;
            for (I1 = Imin1 ; I1 <= Imax1 ; I1 += Istep1) {
              index1 = nint(MAX_VAL * (I1 - Imin1) / (Imax1-Imin1)) ;
              VECTOR_ELT(v, 1) = I1-gc1->means[0] ;
              for (I2 = Imin2 ; I2 < Imax2 ; I2 += Istep2) {
                index2 = nint(MAX_VAL*(I2 - Imin2) / (Imax2-Imin2)) ;
                for (I3 = Imin3 ; I3 < Imax3 ; I3 += Istep3) {
                  index3 = nint(MAX_VAL*(I3 - Imin3) / (Imax3-Imin3)) ;
                  fprintf(fp, "%d %d %d %f %f %f %f\n", index1, index2, index3, I1, I2, I3, pI_c[i1][index1][index2][index3]) ;
                }
              }
            }
            fclose(fp) ;
          }
        }

        for (n = 0, amax =  atotal = 0.0, i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
          for (i1 = 0 ;  i1 < gcap->nlabels ; i1++) {
            l1 = gcap->labels[i1] ;
            gc1 = GCAfindPriorGC(gca, xp, yp, zp,  l1) ;
            if (!gc1 || gc1->regularized)
              continue ;
            label1_count[l1]++ ;
            for (i2 = 0 ;  i2 < i1 ; i2++) {
              l2 = gcap->labels[i2] ;
              gc2 = GCAfindPriorGC(gca, xp, yp, zp,  l2) ;
              if (!gc2 || gc2->regularized)
                continue ;
              if (l1 == l2)
                continue ;
              if (classnum > 0 && l1 != classnum && l2 != classnum)
                continue ;

              n++ ;
              if  ((l1 == 17 &&  l2  == 7) ||(l2 == 17 &&  l1  == 7))
                DiagBreak() ;

              /* marginalize over intensity */
              if (l1 == Left_Hippocampus && l2 == Left_Amygdala)
                DiagBreak() ;
              Imin1 = MAX(Imins1[i1], Imins1[i2]) ;
              Imax1 = MIN(Imaxs1[i1], Imaxs1[i2]) ;
              Imin2 = MAX(Imins2[i1], Imins2[i2]) ;
              Imax2 = MIN(Imaxs2[i1], Imaxs2[i2]) ;
              Imin3 = MAX(Imins3[i1], Imins3[i2]) ;
              Imax3 = MIN(Imaxs3[i1], Imaxs3[i2]) ;
              Istep1 = (Imax1-Imin1)/MAX_VALS ;
              Istep2 = (Imax2-Imin2)/MAX_VALS ;
              Istep3 = (Imax3-Imin3)/MAX_VALS ;
              for (ambiguity = 0,  I1 = Imin1 ; I1 < Imax1 ; I1 += Istep1) {
                index1_c1 = nint(MAX_VAL * (I1 - Imins1[i1]) / (Imaxs1[i1]-Imins1[i1])) ;
                if (index1_c1 < 0 || index1_c1 > MAX_VAL)
                  continue ;
                index1_c2 = nint(MAX_VAL * (I1 - Imins1[i2]) / (Imaxs1[i2]-Imins1[i2])) ;
                if (index1_c2 < 0 || index1_c2 > MAX_VAL)
                  continue ;
                for (I2 = Imin2 ; I2 < Imax2 ; I2 += Istep2) {
                  index2_c1 = nint(MAX_VAL * (I2 - Imins2[i1]) / (Imaxs2[i1]-Imins2[i1])) ;
                  if (index2_c1 < 0 || index2_c1 > MAX_VAL)
                    continue ;
                  index2_c2 = nint(MAX_VAL * (I2 - Imins2[i2]) / (Imaxs2[i2]-Imins2[i2])) ;
                  if (index2_c2 < 0 || index2_c2 > MAX_VAL)
                    continue ;
                  for (amb=0.0,I3 = Imin3 ; I3 < Imax3 ; I3 += Istep3) {
                    index3_c1 = nint(MAX_VAL * (I3 - Imins3[i1]) / (Imaxs3[i1]-Imins3[i1])) ;
                    if (index3_c1 < 0 || index3_c1 > MAX_VAL)
                      continue ;
                    index3_c2 = nint(MAX_VAL * (I3 - Imins3[i2]) / (Imaxs3[i2]-Imins3[i2])) ;
                    if (index3_c2 < 0 || index3_c2 > MAX_VAL)
                      continue ;
                    p1 = pI_c[i1][index1_c1][index2_c1][index3_c1] ;
                    p2 = pI_c[i2][index1_c2][index2_c2][index3_c2] ;
                    amb += (p1*p2*Istep3) ;
                  }
                  ambiguity += (amb*Istep1*Istep2) ;
                }
              }

              ambiguity *= scale*0.5*(gcap->priors[i1]+gcap->priors[i2]) ;
              atotal += ambiguity ;
              label_ambiguity[l1][l2] += ambiguity ;
              label_ambiguity[l2][l1] += ambiguity ;
              label_count[l1][l2]++;
              label_count[l2][l1]++;
              if (ambiguity > amax)
                amax  =  ambiguity ;
            }
          }
        }
        if ( n > 0)
          total_ambiguity += (atotal/n) ;
        if (amax > 0)
          MRIFvox(mri, xp, yp,  zp) = amax  ;
        ntotal += n ;
      }
    }
  }

  {
    FILE *fp ;
    MRI  *mri2;
    float norm;

    mri2 = MRIalloc(MAX_CMA_LABEL+1, MAX_CMA_LABEL+1, 2, MRI_FLOAT) ;
    fp = fopen("label_amb.dat", "w") ;
    for  (l1 = 0 ; l1 <= MAX_CMA_LABEL ; l1++) {
      for (l2 =  0  ; l2  <=  MAX_CMA_LABEL  ; l2++) {
#if 1
        norm = label_count[l1][l2];
#else
        norm = (label1_count[l1]+label1_count[l2])/2 ;
#endif
        if  (norm> 0)
          label_ambiguity[l1][l2] /= (float)norm ;
        fprintf(fp, "%f   ", label_ambiguity[l1][l2])  ;
        MRIFvox(mri2, l1, l2, 0) = label_ambiguity[l1][l2] ;
      }
      fprintf(fp, "\n") ;
    }
    fclose(fp)  ;
    MRIwrite(mri2, "label_amb.mgh") ;
    MRIfree(&mri2);
  }

  if (ntotal == 0) {
    total_ambiguity = scale ;
    printf("**** no valid nodes found ********* \n") ;
  } else {
    total_ambiguity /= ntotal ;
    printf("%d total nodes used in ambiguity calculation...\n", ntotal) ;
    MatrixFree(&m_cov_inv) ;
    VectorFree(&v) ;
    VectorFree(&vtmp) ;
    MatrixFree(&m_cov) ;
  }
  *pamb = total_ambiguity;
  return(mri) ;
}

