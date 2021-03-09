/**
 * @brief program for computing/optimizing registration of a surface to a volume
 *        
 *
 * Program to compute a rigid alignment between a surface and a volume by maximizing the gradient
 * magnitude across the gray/white boundary, divided by its variance
 * Now supports multiple similarity functions.
 */
/*
 * Original Author: Greg Grev
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


/*
BEGINUSAGE --------------------------------------------------------------

 mri_register_to_volume
  --reg regfile
  --mov fvol
  --surf surface   : surface to read in
  --pial pial surface name   : pial surface to read in
  --pial_only pial surface name   : pial surface to read in (don't use white in similarity)

  --median        : apply median filter
  --patch patch   :  patch  to read in
  --tx-mmd txmin txmax txdelta : translation (mm) in x
  --ty-mmd tymin tymax tydelta : translation (mm) in y
  --tz-mmd tzmin tzmax tzdelta : translation (mm) in z
  --ax-mmd axmin axmax axdelta : rotation (deg) about x
  --ay-mmd aymin aymax aydelta : rotation (deg) about y
  --az-mmd azmin azmax azdelta : rotation (deg) about z

  --cost costfile

  --interp interptype : interpolation trilinear or nearest (def is trilin)
  --no-crop: do not crop anat (crops by default)
  --profile : print out info about exec time

  --noise stddev : add noise with stddev to input for testing sensitivity
  --seed randseed : for use with --noise
  --skip min max  : # of vertices to skip in similarity function (for speed)
  --sigma min max  : size of blurrin kernels to use
  --CNR           : use CNR-based similarity function
  --border border : size of border region to ignore

  --out-reg outreg : reg at lowest cost (updated continuously)

ENDUSAGE ---------------------------------------------------------------
*/

/*
BEGINHELP --------------------------------------------------------------

FORMATS

Data file format can be specified implicitly (through the path name)
or explicitly. All formats accepted by mri_convert can be used.

BUGS

sinc interpolation is broken except for maybe COR to COR.


BUG REPORTING

Report bugs to analysis-bugs@nmr.mgh.harvard.edu. Include the following
formatted as a list as follows: (1) command-line, (2) directory where
the program was run (for those in the MGH-NMR Center), (3) version,
(4) text output, (5) description of the problem.

SEE ALSO

mri_vol2vol mri_convert, tkregister2 mri_segreg


ENDHELP --------------------------------------------------------------

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"

#include "matrix.h"
#include "mri.h"
#include "version.h"
#include "mri2.h"
#include "mri_identify.h"
#include "MRIio_old.h"
#include "registerio.h"
#include "resample.h"
#include "gca.h"
#include "gcamorph.h"
#include "fio.h"
#include "cmdargs.h"
#include "pdf.h"
#include "timer.h"
#include "numerics.h"
#include "mri_circulars.h"

#ifdef X
#undef X
#endif

#define MIN_RES .1 // reduce until mri->xsize is less than this in mm

static int write_lta(MATRIX *m, char *fname, MRI_SURFACE *mris, MRI *mri_reg) ;

  //static int write_register_dat(MATRIX *m, char *fname, MRI_SURFACE *mris, MRI *mri, char *subject) ;
double *GetCosts(MRI *mri_reg, MRI *seg, MATRIX *R0, MATRIX *R, double *p, double *costs);
int Min1D(MRI *mri_reg, MRI_SURFACE *mris, MATRIX *R, double *p, char *costfile, double *costs);

// For some reason, this does not seemed to be defined in math.h
double round(double x);
static int write_snapshot(MRI_SURFACE *mris, MATRIX *R0, char *fname, int n) ;

static float compute_powell_rigid_sse(float *p) ;
static int powell_minimize_rigid(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *mat, int skip,
                                 double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                           MRI *mri_mask, MATRIX *m, int skip, double scale, int diag));

static MRI *mri_surface_edges(MRI *mri, MRI *mri_grad, float ndist, float tdist, int dir, MRI *mri_edge) ;;
static int  parse_commandline(int argc, char **argv);
static double mrisRegistrationCNRSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag);
static double mrisRegistrationOverlapSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag);
static double mrisRegistrationDistanceSimilarity(MRI_SURFACE *mris, MRI *mri_dist, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag);
static double mrisRegistrationGradientSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag);
static double mrisRegistrationGradientSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag);
static double mrisRegistrationGradientNormalSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
#include "tags.h"
static int istringnmatch(char *str1, const char *str2, int n);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

static double tscale = 5.0 ;
static int read_median = 0 ;
static int read_grad = 0 ;
static int cropy0 = -1 ; 
static int cropy1 = -1 ;
static FILE *logfp = NULL ;
static int debug = 0, gdiagno = -1;
static int ndilates = 2 ;
static char *patch_fname = NULL ;
static char *label_fname = NULL ;
static char *vol_fname=NULL;
static char *surf_fname = NULL ;
static char *pial_fname = NULL ;
static int pial_only = 0 ;   
static char *regfile=NULL;
static char *outregfile=NULL;
static const char *interpmethod = "trilinear";
static int   interpcode = 0;
static int   sinchw;
static double   max_trans = 200 ;
static double   max_rot = 20 ;  // degrees
static int max_skip = 32 ;
static int min_skip = 8 ;

static double max_sigma = 2 ;
static double min_sigma = .5 ;

static int do_global_search = 1 ;
static int ncalls = 0 ;

static MRI *mri_reg, *mri_grad = NULL ;

MATRIX *R0;

char *SUBJECTS_DIR=NULL;
char *subject = NULL;

//static float ipr, bpr, intensity;
//static int float2int ;
static int nargs;

static int niter = 0 ;
static int write_iter = 1 ;

static char *SegRegCostFile = NULL;

static int SynthSeed = -1;
static int AddNoise = 0;
static double NoiseStd;

static int UseASeg = 0;
static int DoCrop = 0;
static int DoProfile = 0;
static int mask_size = 20 ;
static int apply_median_filter = 0 ;
#define NMAX 100
static int ntx=0, nty=0, ntz=0, nax=0, nay=0, naz=0;
static double txlist[NMAX],tylist[NMAX],tzlist[NMAX];
static double axlist[NMAX],aylist[NMAX],azlist[NMAX];

static double find_optimal_translations(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *R0, 
                                        double max_trans, int skip,
                                        double (*similarity_func)
                                        (MRI_SURFACE *mris, MRI *mri_reg, 
                                         MRI *mri_mask, MATRIX *m, int skip, 
                                         double scale, int diag)) ;
static double find_optimal_rotations(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *R0, 
                                     double max_rot, 
                                     int skip,
                                     double (*similarity_func)
                                     (MRI_SURFACE *mris, MRI *mri_reg, 
                                      MRI *mri_mask, MATRIX *m, int skip, 
                                      double scale, int diag)) ;
static double find_optimal_rigid_alignment(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *R0, 
                                           double min_trans, double max_trans, double min_rot, double max_rot, int skip,
                                           double (*similarity_func)
                                           (MRI_SURFACE *mris, MRI *mri_reg, 
                                            MRI *mri_mask, MATRIX *m, int skip,
                                            double scale, int diag)) ;

static double (*similarity_func)
     (MRI_SURFACE *mris, MRI *mri_reg, 
      MRI *mri_mask, MATRIX *m, int skip, double scale, 
      int diag) = mrisRegistrationGradientNormalSimilarity ;

/*---------------------------------------------------------------*/
int
main(int argc, char **argv) 
{
  char          *saved_pial_fname ;
  MRI_SURFACE   *mris ;
  int           skip, i, msec ;
  MRI           *mri_kernel, *mri_smooth, *mri_mask, *mri_mag ;
  double        sigma, last_similarity, similarity ;
#if 0
  MRI           *mri_dist;
  HISTOGRAM     *h ;
#endif
  MATRIX        *m_save = NULL ;
  Timer then ;

  nargs = handleVersionOption(argc, argv, "mris_register_to_volume");
  if(nargs && argc - nargs == 1) exit (0);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  then.reset() ;
  if(argc == 0) 
    usage_exit();

  parse_commandline(argc, argv);
  if(gdiagno > -1) 
    Gdiag_no = gdiagno;
  check_options();
  dump_options(stdout);

  printf("Loading surface %s\n", surf_fname);
  mris = MRISread(surf_fname) ;
  if (mris == NULL)
    ErrorExit(Gerror, "") ;
  MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
  if (pial_fname)
  {
    if (MRISreadVertexPositions(mris, pial_fname) != NO_ERROR)
      ErrorExit(Gerror, "") ;
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
    MRISsaveNormals(mris, PIAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
  }

  if (label_fname)
  {
    LABEL *area ;
    
    area = LabelRead("", label_fname) ;
    if (area == NULL)
      ErrorExit(ERROR_NOFILE, "%s: could not read label %s", Progname, label_fname) ;
    LabelRipRestOfSurface(area, mris) ;
    LabelFree(&area) ;
    MRISdilateRipped(mris, ndilates) ;
    printf("using surface label with %d vertices\n", MRISvalidVertices(mris)) ;
  }

  if (patch_fname)
  {
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISreadPatchNoRemove(mris, patch_fname) ;
    MRISdilateRipped(mris, ndilates) ;
    printf("using surface patch with %d vertices\n", MRISvalidVertices(mris)) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRIScomputeMetricProperties(mris) ;
  }
  MRISresetNeighborhoodSize(mris, 2) ;

  if (1)
  {
    if (read_median == 0) // otherwise will be read below
    {
      MRI *mri_tmp ;
      int y0, y1 ;
      printf("Loading mov %s\n", vol_fname);
      mri_reg = MRIread(vol_fname);
      if (mri_reg == NULL) 
        ErrorExit(ERROR_NOFILE, "%s: could not read movable volume from %s\n", vol_fname) ;
      
      y0 = mask_size ; 
      y1 = (mri_reg->height-2*mask_size) ;
      if (cropy0 >= 0 && cropy0 > y0)
        y0 = cropy0 ;
      if (cropy1 >= 0 && cropy1 < y1)
        y1 = cropy1 ;
      mri_tmp = MRIextract(mri_reg, NULL, mask_size, y0, mask_size,
                           (mri_reg->width-2*mask_size),
                           y1,
                           (mri_reg->depth-2*mask_size)) ;
      MRIfree(&mri_reg) ;
      mri_reg = mri_tmp ;
    }
  }
  else
    mri_reg = MRIread("reduced.mgz") ;
  if (outregfile)
  {
    char fname[STRLEN] ;
    FileNameOnly(outregfile, fname) ;
    FileNameRemoveExtension(fname, fname) ;
    strcat(fname, ".log") ;
    logfp = fopen(fname, "w") ;
  }
  if (apply_median_filter) /* -median option ... modify 
                              mri_reg using the filter */
  {
    MRI        *mri_tmp ;
    char       fname_median[STRLEN] ;

    FileNameRemoveExtension(vol_fname, fname_median) ;
    strcat(fname_median, ".median.mgz") ;
    if (read_median == 0)
    {
      printf("applying median filter to input volume...\n") ;
      mri_tmp = MRIcopy(mri_reg, NULL) ;  // copy stuff outside of range
      MRImedian(mri_reg, mri_tmp, 3, NULL) ;
      printf("writing median filtered volume to %s\n", fname_median) ;
      MRIwrite(mri_tmp, fname_median) ;
      MRIfree(&mri_reg) ;
      mri_reg = mri_tmp ;
    }
    else
    {
      printf("reading median filtered input volume %s...\n", fname_median) ;
      mri_reg = MRIread(fname_median) ;
      if (mri_reg == NULL)
        ErrorExit(ERROR_NOFILE, "%s: could not read median volume from %s\n", fname_median) ;
    }
  }
  mri_mask = MRIclone(mri_reg, NULL) ;
  MRIsetValues(mri_mask, 1) ;
  MRIeraseBorderPlanes(mri_mask, 1) ;

  if (regfile)
  {
    printf("reading in previously computed registration from %s\n", regfile);
    R0 = regio_read_surfacexform_from_register_dat(regfile, mris, mri_reg, &subject);
    if (R0 == NULL)
      exit(Gerror) ;
  }

  for (skip = max_skip ; skip >= min_skip;  skip /= 2)
  {
    for (sigma = max_sigma ; sigma >= min_sigma ; sigma /= 2)
    {
      char       fname_mag[STRLEN], fname_grad[STRLEN], fname_only[STRLEN] ;

      printf("---------------- skip = %d, sigma = %2.1f ---------------------\n", skip, sigma);
      printf("computing gradient at scale %2.3f...\n",sigma) ;
      if (logfp)
      {
        fprintf(logfp, "---------------- skip = %d, sigma = %2.1f ---------------------\n", skip, sigma);
        fprintf(logfp, "computing gradient at scale %2.3f...\n",sigma) ;
        fflush(logfp) ;
      }
      FileNameRemoveExtension(vol_fname, fname_only) ;
      sprintf(fname_grad, "%s.grad.sigma%2.3f.mgz", fname_only, sigma) ;
      sprintf(fname_mag, "%s.mag.sigma%2.3f.mgz", fname_only, sigma) ;
      if (read_grad == 0)  // compute gradient of smoothed image
      {
        mri_kernel = MRIgaussian1d(sigma/mri_reg->xsize, -1) ;
        mri_smooth = MRIconvolveGaussian(mri_reg, NULL, mri_kernel) ;
        mri_mag = MRIcloneDifferentType(mri_smooth, MRI_FLOAT) ;
        mri_grad = MRIsobel(mri_smooth, mri_grad, mri_mag) ;

        printf("writing mag (%s) and grad (%s) volumes\n", fname_mag, fname_grad) ;
        MRIwrite(mri_mag, fname_mag) ;
        MRIwrite(mri_grad, fname_grad) ;
        MRIeraseBorderPlanes(mri_grad, 1) ;
        if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
        {
          MRIwrite(mri_reg, "r.mgz") ;
        }
        MRIfree(&mri_kernel) ; MRIfree(&mri_smooth) ; 
      }
      else
      {
        printf("reading mag (%s) and grad (%s) volumes\n", fname_mag, fname_grad) ;
        mri_grad = MRIread(fname_grad) ;
        mri_mag = MRIread(fname_mag) ;
        if (mri_grad == NULL || mri_mag == NULL)
          ErrorExit(Gerror, "") ;
      }
      if (0)
      {
        MRI *mri_edge, *mri_reduced, *mri_tmp, *mri_reduced_grad  ;
        int  nreductions = 0;

        mri_reduced = MRIcopy(mri_reg, NULL) ;
        while (mri_reduced->xsize < MIN_RES)
        {
          char fname[STRLEN] ;

          mri_tmp = mri_reduced ;
          mri_reduced = MRIreduce(mri_tmp, NULL) ;
          MRIfree(&mri_tmp) ; 
          nreductions++ ;
          sprintf(fname, "red%d.mgz", nreductions) ;
          MRIwrite(mri_reduced, fname) ;
        }

        mri_kernel = MRIgaussian1d(sigma/mri_reduced->xsize, -1) ;
        mri_smooth = MRIconvolveGaussian(mri_reduced, NULL, mri_kernel) ;
        MRIwrite(mri_reduced, "red.mgz") ;
        MRIwrite(mri_smooth, "rsmooth.mgz") ;
        mri_smooth = MRIconvolveGaussian(mri_reduced, NULL, mri_kernel) ;
        mri_reduced_grad = MRIsobel(mri_smooth, NULL, NULL) ;
        MRIwrite(mri_reduced_grad, "rgrad.mgz") ;
        mri_edge = mri_surface_edges(mri_reduced, mri_reduced_grad, 1.0, 2.0, 1.0, NULL) ;
        MRIwrite(mri_edge, "edge.mgz") ;
        MRIfree(&mri_smooth) ;
        mri_smooth = MRIconvolveGaussian(mri_edge, NULL, mri_kernel) ;
        MRIfree(&mri_edge) ; mri_edge = mri_smooth ;
        MRIwrite(mri_edge, "sedge.mgz") ;
        MRIfree(&mri_kernel) ;
        while (nreductions > 0)
        {
          char fname[STRLEN] ;

          mri_tmp = MRIupsample2(mri_edge, NULL) ;
          MRIfree(&mri_edge) ;
          mri_edge = mri_tmp ;
          sprintf(fname, "up%d.mgz", nreductions) ;
          MRIwrite(mri_edge, fname) ;
          nreductions-- ;
        }
        MRIwrite(mri_edge, "eup.mgz") ;
        {
          MRI       *mri_distance ;
          HISTOGRAM *h, *hcdf ;
          int       b ;
          float     thresh ;

          h = MRIhistogram(mri_edge, 100) ;
          HISTOplot(h, "h.plt") ;
          HISTOclearZeroBin(h) ;
          hcdf = HISTOmakeCDF(h, NULL) ;
          HISTOplot(hcdf, "cdf.plt") ;
          b = HISTOfindBinWithCount(hcdf, 0.8) ;
          thresh = h->bins[b] ;
          MRIbinarize(mri_edge, mri_edge, thresh, 0, 1) ;
          mri_distance = MRIdistanceTransform(mri_edge, NULL, 1, 1000, DTRANS_MODE_UNSIGNED, NULL) ;
          MRIscalarMul(mri_distance, mri_distance, -1) ;
          MRIwrite(mri_distance, "dist.mgz") ;
          mri_edge->outside_val = -1000;
          MRIfree(&mri_edge) ; mri_edge = mri_distance ;
          HISTOfree(&h) ; HISTOfree(&hcdf) ;
        }

        find_optimal_translations(mris, mri_edge,mri_mask,R0, max_trans, skip, 
                                  mrisRegistrationOverlapSimilarity) ;
        if (0)
          find_optimal_rotations(mris, mri_edge, mri_mask, R0, max_rot, skip, 
                                 mrisRegistrationOverlapSimilarity) ;
        find_optimal_rigid_alignment(mris, mri_edge, mri_mask, R0, 0, max_trans/tscale,
                                     .1, max_rot, skip, 
                                     mrisRegistrationOverlapSimilarity);
        MRIfree(&mri_edge) ; MRIfree(&mri_reduced) ;
      }
      if (1) // threshold gradient image
      {
        int       b ;
        double    thresh ;
        HISTOGRAM *h, *hcdf ;

        h = MRIhistogram(mri_mag, 100) ;
        HISTOplot(h, "h.plt") ;
        HISTOclearZeroBin(h) ;
        hcdf = HISTOmakeCDF(h, NULL) ;
        HISTOplot(hcdf, "cdf.plt") ;
        thresh = .5 ;
        b = HISTOfindBinWithCount(hcdf, thresh) ; // only keep top 20% of edges
        thresh = h->bins[b] ;
        MRIbinarize(mri_mag, mri_mag, thresh, 0, 1) ;
        MRImask(mri_grad, mri_mag, mri_grad, 0,0) ;
        MRIwrite(mri_grad, "mgrad.mgz") ;
        MRIwrite(mri_grad, "mmag.mgz") ;
        HISTOfree(&h) ; HISTOfree(&hcdf) ;
      }

      if (do_global_search != 0)
      {
        find_optimal_translations(mris, mri_reg, mri_mask,R0, max_trans, skip, similarity_func) ;
#if 0
        find_optimal_rotations(mris, mri_reg, mri_mask, R0, max_rot, skip, similarity_func) ;
#endif
        find_optimal_rigid_alignment(mris, mri_reg, mri_mask, R0, 2*mri_reg->xsize, max_trans/tscale, 
                                     1, max_rot, skip, similarity_func);
        find_optimal_rigid_alignment(mris, mri_reg, mri_mask, R0, mri_reg->xsize/2, max_trans/10,
                                     .1, max_rot/4, skip, similarity_func);  // redo to finer scale
      }
      //      powell_minimize_rigid(mris, mri_reg, mri_mask, R0, skip,similarity_func);
      if (pial_only)
        similarity = (*similarity_func)(mris, mri_reg, mri_mask,R0,0, 1, 0) ;
      else
        similarity =  
          mrisRegistrationCNRSimilarity(mris, mri_reg, mri_mask,R0,0, 1, 0) ;

      i = 0 ;
      do
      {
        printf("-------- loop %d: starting similarity %2.1f ----------\n",
               ++i, similarity) ;
        m_save = MatrixCopy(R0, m_save) ;
        powell_minimize_rigid(mris, mri_reg, mri_mask, R0, 0, similarity_func);
        if (pial_fname && !pial_only)
        {
          saved_pial_fname = pial_fname ;
          pial_fname = NULL ;
          printf("!!!!!! disabling pial surface !!!!!!!!!!!!\n") ;
          powell_minimize_rigid(mris, mri_reg, mri_mask, R0, 0, similarity_func);
          pial_fname = saved_pial_fname ;
        }
        if (similarity_func != mrisRegistrationCNRSimilarity && !pial_only)
          powell_minimize_rigid(mris, mri_reg, mri_mask, R0, 0, 
                                mrisRegistrationCNRSimilarity);
        last_similarity = similarity ;
        if (pial_only)
          similarity = (*similarity_func)(mris, mri_reg, mri_mask,R0,0, 1, 0) ;
        else
          similarity =  
            mrisRegistrationCNRSimilarity(mris, mri_reg, mri_mask,R0,0, 1, 0) ;
        if (similarity < last_similarity)
        {
          printf("similarity decreased - restoring saved registration\n") ;
          MatrixCopy(m_save, R0) ;
        }
      } while (similarity > last_similarity) ;

      printf("saving current registration to %s\n", outregfile) ;
      regio_write_surfacexform_to_register_dat(R0, outregfile, mris, mri_reg, subject, FLT2INT_ROUND) ;
      write_lta(R0, outregfile, mris, mri_reg) ;
      if (FZERO(min_sigma) && FZERO(sigma))
        break ;
    }
    if (skip == 0)
      break ;
  }

  printf("Reg at min cost was \n");
  MatrixPrint(stdout,R0);
  printf("\n");
  
  if(outregfile){
    printf("Writing optimal reg to %s \n",outregfile);
    write_lta(R0, outregfile, mris, mri_reg) ;
    regio_write_surfacexform_to_register_dat(R0, outregfile, mris, mri_reg, subject, FLT2INT_ROUND) ;
  }

  printf("\n");
  printf("mris_register_to_volume done\n");
  msec = then.milliseconds() ;
  fprintf(stderr,
          "registration took %2.1f minutes\n", (float)msec/(60*1000.0f));

  return(0);
}


/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  int nv,n;
  double vmin, vmax, vdelta;
  

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option,      "--help"))     print_help() ;
    else if (!strcasecmp(option, "--version"))  print_version() ;
    else if (!strcasecmp(option, "--debug"))    debug = 1;
    else if (!strcasecmp(option, "--aseg"))     UseASeg = 1;
    else if (!strcasecmp(option, "--no-crop"))  DoCrop = 0;
    else if (!strcasecmp(option, "--crop"))     DoCrop = 1;
    else if (!strcasecmp(option, "--rm"))       read_median = 1 ;
    else if (!strcasecmp(option, "--rg"))       read_grad = 1 ;
    else if (!strcasecmp(option, "--noglobal")) do_global_search = 0;
    else if (!strcasecmp(option, "--profile"))  DoProfile = 1;
    else if (!strcasecmp(option, "--median"))   apply_median_filter = 1;
    else if (istringnmatch(option, "--reg",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--out-reg",0)) {
      if (nargc < 1) argnerr(option,1);
      outregfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--tscale",0)) {
      if (nargc < 1) argnerr(option,1);
      tscale = atof(pargv[0]);
      nargsused = 1;
    } else if (istringnmatch(option, "--mov",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      vol_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--surf",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      surf_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--pial",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      pial_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--pial_only",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      pial_fname = pargv[0];
      pial_only = 1 ;
      nargsused = 1;
    } else if (istringnmatch(option, "--dilate",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      ndilates = atoi(pargv[0]) ;
      nargsused = 1;
    } else if (istringnmatch(option, "--s",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      subject = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--patch",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      patch_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--label",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      label_fname = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--noise",0)) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&NoiseStd);
      AddNoise = 1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--seed")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SynthSeed);
      nargsused = 1;
    } else if (!strcasecmp(option, "--skip")) {
      if (nargc < 2) CMDargNErr(option,2);
      min_skip = atoi(pargv[0]) ;
      max_skip = atoi(pargv[1]) ;
      nargsused = 2;
    } else if (!strcasecmp(option, "--sigma")) {
      if (nargc < 2) CMDargNErr(option,2);
      min_sigma = atof(pargv[0]) ;
      max_sigma = atof(pargv[1]) ;
      nargsused = 2;
    } else if (!strcasecmp(option, "--dist")) {
      nargsused = 0;
      similarity_func = mrisRegistrationDistanceSimilarity;
    } else if (!strcasecmp(option, "--cropy")) {
      cropy0 = atoi(pargv[0]) ;
      cropy1 = atoi(pargv[1]) ;
      nargsused = 2;
    } else if (!strcasecmp(option, "--CNR")) {
      similarity_func = mrisRegistrationCNRSimilarity ;
      nargsused = 0;
    } else if (!strcasecmp(option, "--DOT")) {
      similarity_func = mrisRegistrationGradientNormalSimilarity ;
      nargsused = 0;
    } else if (!strcasecmp(option, "--GRADIENT")) {
      similarity_func = mrisRegistrationGradientSimilarity ;
      nargsused = 0;
    } else if (!strcasecmp(option, "--w")) {
      write_iter = atoi(pargv[0]) ;
      nargsused = 1;
    } else if (!strcasecmp(option, "--max_rot")) {
      if (nargc < 1) CMDargNErr(option,1);
      max_rot = atof(pargv[0]) ;
      nargsused = 1;
    } else if (!strcasecmp(option, "--max_trans")) {
      if (nargc < 1) CMDargNErr(option,1);
      max_trans = atof(pargv[0]) ;
      nargsused = 1;
    } else if (!strcasecmp(option, "--border")) {
      if (nargc < 1) CMDargNErr(option,1);
      mask_size = atoi(pargv[0]) ;
      nargsused = 1;
    } else if (istringnmatch(option, "--tx-mmd",0)) {
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("ntx = %d  %lf %lf %lf\n",nv,vmin,vmax,vdelta);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
	txlist[ntx] = vmin + vdelta*n;
	ntx++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--ty-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("nty = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
	tylist[nty] = vmin + vdelta*n;
	nty++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--tz-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("ntz = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
	tzlist[ntz] = vmin + vdelta*n;
	ntz++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--ax-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("nax = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
	axlist[nax] = vmin + vdelta*n;
	nax++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--ay-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("nay = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
	aylist[nay] = vmin + vdelta*n;
	nay++;
      }
      nargsused = 3;
    } else if (istringnmatch(option, "--az-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%lf",&vmin);
      sscanf(pargv[1],"%lf",&vmax);
      sscanf(pargv[2],"%lf",&vdelta);
      nv = (int)round((vmax-vmin)/vdelta) + 1;
      printf("naz = %d\n",nv);
      if(nv <= 0) exit(1);
      for(n=0; n < nv; n++){
        azlist[naz] = vmin + vdelta*n;
        naz++;
      }
      nargsused = 3;
    } else if ( !strcmp(option, "--gdiagno") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&gdiagno);
      nargsused = 1;
    } else if (istringnmatch(option, "--cost",0)) {
      if (nargc < 1) 
        argnerr(option,1);
      SegRegCostFile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--interp",8)) {
      if (nargc < 1) argnerr(option,1);
      interpmethod = pargv[0];
      nargsused = 1;
      if (!strcmp(interpmethod,"sinc") && CMDnthIsArg(nargc, pargv, 1)) {
        sscanf(pargv[1],"%d",&sinchw);
        nargsused ++;
      }
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if (singledash(option))
        fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
printf("\n");
printf("mris_register_to_volume\n");
printf("  --surf surface\n");
printf("  --pial pial surface name\n");
printf("  --pial_only pial surface name\n");
printf("  --reg regfile\n");
printf("  --noglobal\n");
printf("  --median\n");
printf("  --mri_reg fvol\n");
printf("\n");
printf("  --tx-mmd txmin txmax txdelta : translation (mm) in x\n");
printf("  --ty-mmd tymin tymax tydelta : translation (mm) in y\n");
printf("  --tz-mmd tzmin tzmax tzdelta : translation (mm) in z\n");
printf("  --ax-mmd axmin axmax axdelta : rotation (deg) about x\n");
printf("  --ay-mmd aymin aymax aydelta : rotation (deg) about y\n");
printf("  --az-mmd azmin azmax azdelta : rotation (deg) about z\n");
printf("\n");
printf("  --cost costfile\n");
printf("\n");
printf("  --interp interptype : interpolation trilinear or nearest (def is trilin)\n");
printf("\n");
printf("  --noise stddev : add noise with stddev to input for testing sensitivity\n");
printf("  --seed randseed : for use with --noise\n");
printf("  --skip min max  : # of vertices to skip (starting at max and reducing)\n");
printf("  --sigma min max  : size of blurring kernels to use (starting at max and reducing)\n");
printf("  --CNR           : use CNR-based similarity function\n");
printf("  --max_rot angle : max angle (degrees) to search over\n");
printf("  --max_trans dist :max translation (mm) to search over\n");
printf("  --border border : size of the border region to ignore\n");
printf("  --s subject     : specify name of subject (for register.dat file)\n");
printf("  --dilate ndil   : dilate ripflags ndil times (only with --patch)\n");
printf("  --patch <patch>   : load patch <patch> and limit calculations to it\n") ;
printf("  --label <label>   : load label <label> and limit calculations to it\n") ;
printf("\n");
printf("  --out-reg outreg : reg at lowest cost (updated continuously)\n");
printf("\n");

}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n%s\n\n",getVersion().c_str());
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if (SUBJECTS_DIR==NULL) {
    printf("ERROR: SUBJECTS_DIR undefined.\n");
    exit(1);
  }
  if (vol_fname == NULL) {
    printf("ERROR: No mov volume supplied.\n");
    exit(1);
  }
  if (surf_fname == NULL) {
    printf("ERROR: No surface supplied.\n");
    exit(1);
  }
  if (regfile == NULL) {
    printf("using identity as initial registration\n") ;
    R0 = MatrixIdentity(4, NULL) ;
  }
#if 0
  if (SegRegCostFile == NULL) {
    printf("ERROR: need --cost.\n");
    exit(1);
  }
#endif

  interpcode = MRIinterpCode(interpmethod);
  if (interpcode < 0) {
    printf("ERROR: interpolation method %s unrecognized\n",interpmethod);
    printf("       legal values are nearest, trilin, and sinc\n");
    exit(1);
  }

  if(ntx == 0) {ntx=1; txlist[0] = 0;}
  if(nty == 0) {nty=1; tylist[0] = 0;}
  if(ntz == 0) {ntz=1; tzlist[0] = 0;}
  if(nax == 0) {nax=1; axlist[0] = 0;}
  if(nay == 0) {nay=1; aylist[0] = 0;}
  if(naz == 0) {naz=1; azlist[0] = 0;}

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) 
{
  //  int n;
  fprintf(fp,"movvol %s\n",vol_fname);
  fprintf(fp,"surface %s\n",surf_fname);
  fprintf(fp,"regfile %s\n",regfile);
  if(outregfile) fprintf(fp,"outregfile %s\n",outregfile);
  fprintf(fp,"interp  %s (%d)\n",interpmethod,interpcode);
  if(interpcode == SAMPLE_SINC) fprintf(fp,"sinc hw  %d\n",sinchw);
  fprintf(fp,"Crop      %d\n",DoCrop);
  fprintf(fp,"Profile   %d\n",DoProfile);
  fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);
  fprintf(fp,"ntx %d\n",ntx);
#if 0
  if(ntx > 0){
    fprintf(fp," tx values\n");
    for(n=0; n < ntx; n++) printf("    %2d %g\n",n+1,txlist[n]);
  }
  fprintf(fp,"nty %d\n",nty);
  if(nty > 0){
    fprintf(fp," ty values\n");
    for(n=0; n < nty; n++) printf("    %2d %g\n",n+1,tylist[n]);
  }
  fprintf(fp,"ntz %d\n",ntz);
  if(ntz > 0){
    fprintf(fp," tz values\n");
    for(n=0; n < ntz; n++) printf("    %2d %g\n",n+1,tzlist[n]);
  }
  fprintf(fp,"nax %d\n",nax);
  if(nax > 0){
    fprintf(fp," ax values\n");
    for(n=0; n < nax; n++) printf("    %2d %g\n",n+1,axlist[n]);
  }
  fprintf(fp,"nay %d\n",nay);
  if(nay > 0){
    fprintf(fp," ay values\n");
    for(n=0; n < nay; n++) printf("    %2d %g\n",n+1,aylist[n]);
  }
  fprintf(fp,"naz %d\n",naz);
  if(naz > 0){
    fprintf(fp," az values\n");
    for(n=0; n < naz; n++) printf("    %2d %g\n",n+1,azlist[n]);
  }
#endif
  return;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);
  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
/*------------------------------------------------------------
  istringnmatch() - compare the first n characters of two strings,
  return a 1 if they match (ignoring case), a zero otherwise. If
  n=0, then do a full comparison.
  ------------------------------------------------------------*/
static int istringnmatch(char *str1, const char *str2, int n) {
  if (n > 0  && ! strncasecmp(str1,str2,n)) return(1);
  if (n <= 0 && ! strcasecmp(str1,str2)) return(1);
  return(0);
}

/*-------------------------------------------------------*/
#if 0
double *GetCosts(MRI *mri_reg, MRI *seg, MATRIX *R0, MATRIX *R, 
		 double *p, double *costs)
{
  double angles[0];
  MATRIX *Mrot=NULL, *Mtrans=NULL, *invR=NULL,*vox2vox = NULL;
  MATRIX *Tin, *invTin, *Ttemp;
  extern MRI *out;

  if(R==NULL){
    printf("ERROR: GetCosts(): R cannot be NULL\n");
    return(NULL);
  }

  Tin      = MRIxfmCRS2XYZtkreg(mri_reg);
  invTin   = MatrixInverse(Tin,NULL);
  Ttemp    = MRIxfmCRS2XYZtkreg(seg);

  Mtrans = MatrixIdentity(4,NULL);
  Mtrans->rptr[1][4] = p[0];
  Mtrans->rptr[2][4] = p[1];
  Mtrans->rptr[3][4] = p[2];

  angles[0] = p[3]*(M_PI/180);
  angles[1] = p[4]*(M_PI/180);
  angles[2] = p[5]*(M_PI/180);
  Mrot = MRIangles2RotMat(angles);

  // R = Mtrans*Mrot*R0
  R = MatrixMultiply(Mrot,R0,R);
  R = MatrixMultiply(Mtrans,R,R);
  invR = MatrixInverse(R,invR);

  // vox2vox = invTin*R*Ttemp
  vox2vox = MatrixMultiply(invTin,R,vox2vox);
  MatrixMultiply(vox2vox,Ttemp,vox2vox);
  
  // resample
  MRIvol2Vol(mri_reg,out,vox2vox,interpcode,sinchw);
  
  // compute costs
  costs = SegRegCost(regseg,out,costs);
  
  MatrixFree(&Mrot);
  MatrixFree(&Mtrans);
  MatrixFree(&vox2vox);
  MatrixFree(&Tin);
  MatrixFree(&invTin);
  MatrixFree(&Ttemp);
  MatrixFree(&invR);

  return(costs);
}

/*---------------------------------------------------------------------*/
int Min1D(MRI *mri_reg, MRI_SURFACE *mris, MATRIX *R, double *p, 
	  char *costfile, double *costs)
{
  double q, q0, pp[6], c, copt=0, qopt=0, costsopt[8];
  int nthp, nth, n, hit;
  MATRIX *R0, *Rtmp;
  FILE *fp;

  if(R==NULL) exit(1);
  if(p==NULL) exit(1);
  if(costs==NULL) exit(1);

  for(nthp = 0; nthp < 6; nthp++) 
    pp[nthp] = p[nthp];
  R0 = MatrixCopy(R,NULL);
  Rtmp = MatrixAlloc(4,4,MATRIX_REAL);

  GetCosts(mri_reg, seg, R0, Rtmp, pp, costs);
  copt = costs[7];

  nth = 0;
  for(nthp = 0; nthp < 6; nthp++){
    qopt = 0;
    hit = 0;
    q0 = pp[nthp];

    for(q = -2; q <= 2; q += .2){
      nth ++;
      pp[nthp] = q;

      GetCosts(mri_reg, seg, R0, Rtmp, pp, costs);

      if(costfile != NULL){
	// write costs to file
	fp = fopen(costfile,"a");
	fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
	fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
	fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
	fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
	fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
	fprintf(fp,"\n");
	fclose(fp);
      }
      
      fp = stdout;
      fprintf(fp,"%5d ",nth);
      fprintf(fp,"%7.3lf %7.3lf %7.3lf ",pp[0],pp[1],pp[2]);
      fprintf(fp,"%6.3lf %6.3lf %6.3lf ",pp[3],pp[4],pp[5]);
      fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
      fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
      fprintf(fp,"%8.4lf %8.4lf   %8.4lf ",costs[6],costs[7],copt); // t, cost=1/t
      printf("\n");
      fflush(stdout);

      c = costs[7];
      if(c < copt){
	copt = c;
	qopt = q;
	MatrixCopy(Rtmp,R);
	for(n=0; n<8; n++) costsopt[n] = costs[n];
	hit = 1;
      }

    }
    if(hit) pp[nthp] = qopt;
    else    pp[nthp] = q0;
  } // loop over params

  for(nthp = 0; nthp < 6; nthp++) p[nthp] = pp[nthp];
  for(n=0; n<8; n++) costs[n] = costsopt[n];

  return(nth);
}
#endif

#define SAMPLE_DIST_MM 0.5

static double
mrisRegistrationCNRSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag)
{
  double    similarity, grad, grad_var, grad_mean, xv, yv, zv, nx, ny, nz, 
            mag, gm_mean, wm_mean, sample_dist, vertex_similarity,
            xw, yw, zw, xg, yg, zg ;
  int       vno, num, num_in_fov = 0 ;
  static VECTOR *v1 = NULL, *v2 = NULL ;
  double    gm_var, wm_var, contrast, noise, total_contrast, total_noise ;
  int       n, num_nbrs ;
  MRI       *mri = NULL ;

  if (diag)
    mri = MRIclone(mri_reg, NULL) ;

  if (v1 == NULL)
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = 1.0 ;
  }
  skip++ ;
  MRISclearMarks(mris) ;
  MRISsetVals(mris, 0) ;
  MRIScopyValToVal2(mris) ;
  MRISsetCurvature(mris, 0.0) ;
  grad_var = grad_mean = 0 ;
  sample_dist = SAMPLE_DIST_MM/mri_reg->xsize ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    VERTEX * const v = &mris->vertices[vno];
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;
    num++ ;
    V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
    v2 = MatrixMultiply(m, v1, v2) ;
    MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2), 
                                &xv, &yv, &zv) ;
    if (MRIindexNotInVolume(mri_reg, xv, yv, zv) || 
        MRIgetVoxVal(mri_mask, xv, yv, zv, 0) == 0)
      grad = 0 ;
    else  // in the volume
    {
      num_in_fov++ ;
      V3_X(v1) = v->whitex+v->nx ; V3_Y(v1) = v->whitey+v->ny ; 
      V3_Z(v1) = v->whitez+v->nz ;
      v2 = MatrixMultiply(m, v1, v2) ;
      MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2), 
                                  &nx, &ny, &nz) ;
      nx -= xv ; ny -= yv ; nz -= zv ;  // vertex normal in voxel coords
      mag = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (FZERO(mag))
        mag = 1.0 ;
      nx /= mag;  ny /= mag ; nz /= mag ;
      MRIsampleVolumeDerivativeScale(mri_reg, xv, yv, zv, nx, ny, nz, &grad, .5/mri_reg->xsize) ;
      xw = xv-sample_dist*nx ; yw = yv-sample_dist*ny ; zw = zv-sample_dist*nz;
      xg = xv+sample_dist*nx ; yg = yv+sample_dist*ny ; zg = zv+sample_dist*nz;
      if ((MRIindexNotInVolume(mri_reg, xw, yw, zw) == 0) &&
          (MRIindexNotInVolume(mri_reg, xg, yg, zg) == 0))
      {
        MRIsampleVolume(mri_reg, xw, yw, zw, &wm_mean) ;
        MRIsampleVolume(mri_reg, xg, yg, zg, &gm_mean) ;
        v->val = wm_mean ;
        v->val2 = gm_mean ;
        v->marked = 1 ;
      }
    }
    if (!FZERO(grad))
      grad /= fabs(grad) ;
    grad_mean += grad ;
    grad_var += (grad*grad) ;
  }

  if (num == 0)
    num = 1 ;
  grad_mean /= num ;
  grad_var = grad_var / num - grad_mean*grad_mean ;
  if (FZERO(grad_var))
    grad_var = 1.0 ;

  total_contrast = total_noise = 0.0 ; similarity = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno += skip)
  {
    VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[vno];
    VERTEX                * const v  = &mris->vertices         [vno];
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->curv = 0 ;
    if (v->ripflag || v->marked == 0)
      continue ;
    gm_mean = v->val2 ; gm_var = v->val2*v->val2 ;
    vertex_similarity = 0 ;
    wm_mean = v->val ;  wm_var = v->val*v->val ;
    for (n = 0, num_nbrs = 1 ; n < vt->vtotal ; n++)
    {
      VERTEX const * const vn = &mris->vertices[vt->v[n]] ;
      if (vn->marked == 0 || vn->ripflag)
        continue ;
      num_nbrs++ ;
      gm_mean += vn->val2 ;         wm_mean += vn->val ;
      gm_var += vn->val2*vn->val2 ; wm_var += vn->val*vn->val ;
    }
    if (num_nbrs > 1) // for variance estimate
    {
      gm_mean /= num_nbrs ;
      wm_mean /= num_nbrs ;
      gm_var = gm_var/num_nbrs - gm_mean*gm_mean ;
      wm_var = wm_var/num_nbrs - wm_mean*wm_mean ;
      noise = gm_var + wm_var ;
      if (FZERO(noise) == 0)
      {
        contrast = wm_mean - gm_mean ;
        if (wm_mean < gm_mean)  // enforce direction of contrast
        {
          total_contrast += contrast ;
          total_noise += noise ;
          vertex_similarity = contrast*contrast / noise ;
        }
        v->curv = contrast*contrast  ;  // for diagnostics
      }
      else
        DiagBreak() ;
    }
    if (vertex_similarity > 10000)
      DiagBreak() ;
    similarity += vertex_similarity ;
    if (mri)
    {
      float val = vertex_similarity > 0 ? vertex_similarity : -5 ;
      V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
      v2 = MatrixMultiply(m, v1, v2) ;
      MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2), 
                                  &xv, &yv, &zv) ;
      MRIfillRegion(mri, xv, yv, zv, val < 0 ? -.02 : .02, nint(.25/mri->xsize)) ;
      MRIsetVoxVal(mri, xv, yv, zv, 0, val) ;
    }
  }

  if (diag)
  {
    printf("num_in_fov = %d\n", num_in_fov) ;
    if (!(++ncalls%write_iter))
    {
      char name[STRLEN], fname[STRLEN] ;

      FileNameOnly(outregfile, name) ;
      FileNameRemoveExtension(name, name) ;

      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        VERTEX const * const v = &mris->vertices[vno] ;
        if (v->ripflag == 0)
          continue ;
        V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
        v2 = MatrixMultiply(m, v1, v2) ;
        MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2),V3_Z(v2),
                                    &xv, &yv, &zv) ;
        if (MRIindexNotInVolume(mri, xv, yv, zv) == 0)
          MRIsetVoxVal(mri, xv, yv, zv, 0, -.01) ;
      }

      sprintf(fname, "%s_hvol%03d.mgz", name, ncalls) ;
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        printf("writing hit vol to %s\n", fname) ;
        MRIwrite(mri, fname) ;
      }
    }
    MRIfree(&mri) ;
  }
  return(similarity) ;
}

static double
mrisRegistrationGradientSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag)
{
  double    similarity, grad, grad_var=0.0, grad_mean, xv, yv, zv, nx, ny, nz, 
            mag, gm_mean, wm_mean, sample_dist ;
  int       vno, num, num_in_fov = 0 ;
  VERTEX    *v ;
  static VECTOR *v1 = NULL, *v2 = NULL ;

  if (v1 == NULL)
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = 1.0 ;
  }
  skip++ ;
#if 0
  MRISclearMarks(mris) ;
  MRISsetVals(mris, 0) ;
  MRIScopyValToVal2(mris) ;
  MRISsetCurvature(mris, 0.0) ;
#endif
  grad_var = grad_mean = 0 ;
  sample_dist = SAMPLE_DIST_MM/mri_reg->xsize ;

  for (num = vno = 0 ; vno < mris->nvertices ; vno += skip)
  {
    v = &mris->vertices[vno] ;
    v->marked = 0 ;
    if (v->ripflag)
      continue ;
    num++ ;
    V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
    v2 = MatrixMultiply(m, v1, v2) ;
    MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2), &xv, &yv, &zv) ;
    if (MRIindexNotInVolume(mri_reg, xv, yv, zv) || MRIgetVoxVal(mri_mask, xv, yv, zv, 0) == 0)
      grad = 0 ;
    else  // in the volume
    {
      num_in_fov++ ;
      V3_X(v1) = v->whitex+v->nx ; V3_Y(v1) = v->whitey+v->ny ; V3_Z(v1) = v->whitez+v->nz ;
      v2 = MatrixMultiply(m, v1, v2) ;
      MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2), &nx, &ny, &nz) ;
      nx -= xv ; ny -= yv ; nz -= zv ;  // vertex normal in voxel coords
      mag = sqrt(nx*nx + ny*ny + nz*nz) ;
      if (FZERO(mag))
        mag = 1.0 ;
      nx /= mag;  ny /= mag ; nz /= mag ;
      MRIsampleVolumeDerivativeScale(mri_reg, xv, yv, zv, nx, ny, nz, 
                                     &grad, .5/mri_reg->xsize) ;
      MRIsampleVolume(mri_reg, xv-sample_dist*nx, yv-sample_dist*ny, 
                      zv-sample_dist*nz, &wm_mean) ;
      MRIsampleVolume(mri_reg, xv+sample_dist*nx, yv+sample_dist*ny, 
                      zv+sample_dist*nz, &gm_mean) ;
      v->val = wm_mean ;
      v->val2 = gm_mean ;
      v->marked = 1 ;
    }
    if (!FZERO(grad))
      grad /= fabs(grad) ;
    if (grad > 0)
      grad_mean += grad ;
    grad_var += (grad*grad) ;
  }

  if (diag)
    printf("num_in_fov = %d\n", num_in_fov) ;
  similarity = fabs(grad_mean) ;

  if (num == 0)
    num = 1;
  grad_mean /= num ;
  grad_var = grad_var / num - grad_mean*grad_mean ;
  if (FZERO(grad_var))
    grad_var = 1.0 ;
  //  similarity = grad_mean * grad_mean / grad_var ;
  return(similarity) ;  
}
static double
mrisRegistrationGradientNormalSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag)
{
  double    similarity, xv, yv, zv, nx, ny, nz, mag, dx, dy, dz, white_dot, 
            pial_dot, theta ;
  int       vno, num, num_in_fov = 0 ;
  VERTEX    *v ;
  static VECTOR *v1 = NULL, *v2 = NULL ;
  MRI       *mri = NULL;
  char fname[STRLEN] ;

  if (diag)
    mri = MRIclone(mri_reg, NULL) ;

  if (v1 == NULL)
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = 1.0 ;
  }
  for (similarity=0.0, num = vno = 0 ; vno < mris->nvertices; vno += (skip+1))
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    v->marked = 0 ;
    if (v->ripflag)
      continue ;
    num++ ;
    if (pial_only == 0)
    {
      V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
      v2 = MatrixMultiply(m, v1, v2) ;
      MRISsurfaceRASToVoxelCached(mris, mri_reg, 
                                  V3_X(v2), V3_Y(v2), V3_Z(v2), &xv, &yv, &zv);
      if (MRIindexNotInVolume(mri_reg, xv, yv, zv) || 
          MRIgetVoxVal(mri_mask, xv, yv, zv, 0) == 0)
        white_dot = 0 ;
      else  // in the volume
      {
        num_in_fov++ ;
        V3_X(v1) = v->whitex+v->nx ; V3_Y(v1) = v->whitey+v->ny ; 
        V3_Z(v1) = v->whitez+v->nz ;
        v2 = MatrixMultiply(m, v1, v2) ;
        MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2),
                                    &nx, &ny, &nz) ;
        nx -= xv ; ny -= yv ; nz -= zv ;  // vertex normal in voxel coords
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        if (FZERO(mag))
          mag = 1.0 ;
        nx /= mag;  ny /= mag ; nz /= mag ;
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 0, &dx) ;
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 1, &dy) ;
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 2, &dz) ;
        mag = sqrt(dx*dx + dy*dy + dz*dz) ; 
        if (FZERO(mag))
          mag = 1.0 ;
        dx /= mag;  dy /= mag ; dz /= mag ;
        
        v->marked = 1 ;
        white_dot = nx*dx + ny*dy + nz*dz ;
        theta = acos(white_dot) ;
#define MAX_THETA (RADIANS(60))
#if 1
        if (white_dot < 0)
          white_dot = 0 ;
#else
        if (fabs(theta) > MAX_THETA)
          white_dot = 0 ;
        else
          DiagBreak() ;
#endif
        if (mri)
        {
          float val ;
          
          if (FZERO(white_dot))
            val = -.5 ;
          else
            val = white_dot ;
          MRIfillRegion(mri, nint(xv), nint(yv), nint(zv), val < 0 ? -.02 : 0.02, nint(.25/mri->xsize)) ;
          MRIsetVoxVal(mri, nint(xv), nint(yv), nint(zv), 0, val) ;
        }
      }
    }
    else
      white_dot = 0 ;

    if (pial_fname)  // only if pial loaded and white was reasonable
    {
      V3_X(v1) = v->pialx ; V3_Y(v1) = v->pialy ; V3_Z(v1) = v->pialz ;
      v2 = MatrixMultiply(m, v1, v2) ;
      MRISsurfaceRASToVoxelCached(mris, mri_reg, 
                                  V3_X(v2), V3_Y(v2), V3_Z(v2), &xv, &yv, &zv);
      if (MRIindexNotInVolume(mri_reg, xv, yv, zv) || 
          MRIgetVoxVal(mri_mask, xv, yv, zv, 0) == 0)
        pial_dot = 0 ;
      else  // in the volume
      {
        num_in_fov++ ;
        V3_X(v1) = v->pialx+v->pnx ; V3_Y(v1) = v->pialy+v->pny ; 
        V3_Z(v1) = v->pialz+v->pnz ;
        v2 = MatrixMultiply(m, v1, v2) ;
        MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2),
                                    &nx, &ny, &nz) ;
        nx -= xv ; ny -= yv ; nz -= zv ;  // vertex normal in voxel coords
        mag = sqrt(nx*nx + ny*ny + nz*nz) ;
        if (FZERO(mag))
          mag = 1.0 ;
        nx /= mag;  ny /= mag ; nz /= mag ;
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 0, &dx) ;
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 1, &dy) ;
        MRIsampleVolumeFrame(mri_grad, xv, yv, zv, 2, &dz) ;
        mag = sqrt(dx*dx + dy*dy + dz*dz) ; 
        if (FZERO(mag))
          mag = 1.0 ;
        dx /= mag;  dy /= mag ; dz /= mag ;

        v->marked = 1 ;
        pial_dot = -1*(nx*dx + ny*dy + nz*dz) ; // pial contrast should be inverse
        if (pial_dot < 0)
          pial_dot = 0 ;
        else
          DiagBreak() ;
        if (mri)
          MRIsetVoxVal(mri, nint(xv), nint(yv), nint(zv), 0, pial_dot) ;
      }
      v->val = white_dot ;
      //      if (!FZERO(white_dot) && !FZERO(pial_dot))
      similarity += (white_dot+pial_dot)  ;
    }
    else
      similarity += white_dot ;  // just add white if no pial
    
  }

  if (diag)
  {
    //    printf("num_in_fov = %d\n", num_in_fov) ;
    if (!(++ncalls%write_iter))
    {
      char name[STRLEN] ;
      FileNameOnly(outregfile, name) ;
      FileNameRemoveExtension(name, name) ;

      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        v = &mris->vertices[vno] ;
        if (v->ripflag == 0)
          continue ;
        V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
        v2 = MatrixMultiply(m, v1, v2) ;
        MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2),V3_Z(v2),
                                    &xv, &yv, &zv) ;
        if (MRIindexNotInVolume(mri, xv, yv, zv) == 0)
          MRIsetVoxVal(mri, xv, yv, zv, 0, -.01) ;
      }
      if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
      {
        sprintf(fname, "%s_hvol%03d.mgz", name, ncalls) ;
        printf("writing hit vol to %s\n", fname) ;
        MRIwrite(mri, fname) ;
      }
    }
    MRIfree(&mri) ;
  }
    

  if (num == 0)
    num = 1;
  return(similarity) ;  
}

static double
mrisRegistrationDistanceSimilarity(MRI_SURFACE *mris, MRI *mri_dist, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag)
{
  double    similarity, xv, yv, zv, dist ;
  int       vno, num, num_in_fov = 0 ;
  VERTEX    *v ;
  static VECTOR *v1 = NULL, *v2 = NULL ;

  if (v1 == NULL)
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = 1.0 ;
  }
  skip++ ;

  for (similarity = 0.0, num = vno = 0 ; vno < mris->nvertices ; vno += skip)
  {
    v = &mris->vertices[vno] ;
    v->marked = 0 ;
    if (v->ripflag)
      continue ;
    num++ ;
    V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
    v2 = MatrixMultiply(m, v1, v2) ;
    MRISsurfaceRASToVoxelCached(mris, mri_reg, V3_X(v2), V3_Y(v2), V3_Z(v2), &xv, &yv, &zv) ;
    if (MRIindexNotInVolume(mri_dist, xv, yv, zv) || 
        MRIgetVoxVal(mri_mask, xv, yv, zv, 0) == 0)
      dist = mri_dist->outside_val ;
    else  // in the volume
    {
      num_in_fov++ ;
      MRIsampleVolume(mri_dist, xv, yv, zv, &dist) ;
      v->marked = 1 ;
    }
    if (dist < 0)
      dist = 0 ;
    similarity += -dist ;
  }

  if (diag)
    printf("num_in_fov = %d\n", num_in_fov) ;

  if (num == 0)
    num = 1;
  similarity /= (double)num ;
  return(similarity) ;  
}

static double
find_optimal_translations(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *M_min, double max_trans, 
                          int skip,
                          double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                    MRI *mri_mask, MATRIX *m, int skip, double scale, int diag))
{
  double   tx, ty, tz, scale, max_similarity, similarity,
           xtrans, ytrans, ztrans ;
  MATRIX   *m_trans, *m = NULL ;
  int      good_step = 0 ;

  if (niter == 0 && write_iter > 0)
    write_snapshot(mris, M_min, outregfile, niter) ;
  m_trans = MatrixIdentity(4, NULL) ;

  max_similarity = similarity = 
    (*similarity_func)(mris, mri_reg, mri_mask, M_min, skip, skip/8.0, 0) ;
#ifdef WHALF
#undef WHALF
#endif
#define WHALF 11
  scale = max_trans/WHALF ; 
  do
  {
    for (tx = -WHALF ; tx <= WHALF ; tx++)
      for (ty = -WHALF ; ty <= WHALF ; ty++)
        for (tz = -WHALF ; tz <= WHALF ; tz++)
        {
          xtrans = scale*tx ;
          ytrans = scale*ty ;
          ztrans = scale*tz ;
          *MATRIX_RELT(m_trans, 1, 4) = xtrans ;
          *MATRIX_RELT(m_trans, 2, 4) = ytrans ;
          *MATRIX_RELT(m_trans, 3, 4) = ztrans ;
          m = MatrixMultiply(m_trans, M_min, m) ;
          similarity = (*similarity_func)(mris, mri_reg,  mri_mask, m, skip, skip/8.0, 0) ;
          if (similarity > max_similarity)
          {
            MatrixCopy(m, M_min) ;
            max_similarity = similarity ;
            //            MatrixPrint(Gstdout, M_min) ;
            printf("%03d: new max %2.2f found at T = (%2.1f, %2.1f, %2.1f)\n",
                   niter+1, max_similarity, xtrans, ytrans, ztrans) ;
            if ((++niter % write_iter) == 0)
              write_snapshot(mris, M_min, outregfile, niter) ;
            if (logfp)
            {
              fprintf(logfp, "%03d: new max %2.2f found at T = (%2.1f, %2.1f, %2.1f)\n",
                      niter, max_similarity, xtrans, ytrans, ztrans) ;
              fflush(logfp) ;
            }
            
            similarity = (*similarity_func)(mris, mri_reg,  mri_mask, m, skip,skip/8.0, 1) ;
            good_step = 1 ;
          }
        }
    if (good_step == 0)
    {
      scale /= 2 ;
      if (scale >= mri_reg->xsize/2)
        printf("reducing scale to %2.4f mm\n", scale) ;
    }
    good_step = 0 ;
  } while (scale >= mri_reg->xsize/2);

  MatrixFree(&m) ; MatrixFree(&m_trans) ;
  return(max_similarity) ;
}

#ifdef WHALF
#undef WHALF
#endif
#define WHALF 5
static double
find_optimal_rigid_alignment(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *M_min, 
                             double min_trans, double max_trans, 
                             double min_rot, double max_rot, int skip,
                             double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                       MRI *mri_mask, MATRIX *m, int skip, double scale, int diag))

{
  double   tx, ty, tz, rot_scale, trans_scale, max_similarity, similarity, 
           xtrans, ytrans, ztrans, rx, ry, rz, whalf = WHALF ;
  MATRIX   *M_trans, *M_test = NULL, *M_rot, *M_tmp, *M_ctr, *M_ctr_inv  ;
  int      good_step = 0 ;
  double   angles[3] ;

  if (min_trans <= 0)
    min_trans = mri_reg->xsize/2 ;

  printf("computing optimal rigid alignment in R=[%2.1f, %2.1f], T=[%2.2f, %2.2f]\n", -max_rot, max_rot, -max_trans, max_trans) ;
  M_trans = MatrixIdentity(4, NULL) ;
  M_ctr = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(M_ctr, 1,4) = mris->xctr ;
  *MATRIX_RELT(M_ctr, 2,4) = mris->yctr ;
  *MATRIX_RELT(M_ctr, 3,4) = mris->zctr ;
  M_ctr_inv = MatrixInverse(M_ctr, NULL) ;
  M_tmp = MatrixCopy(M_trans, NULL) ;
  M_test = MatrixCopy(M_trans, NULL) ;
  max_similarity = similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_min, skip, skip/8.0, 0) ;
  trans_scale = max_trans/whalf ; 
  rot_scale = max_rot/whalf ;
  do
  {
    for (rot_scale = max_rot/whalf ; 
         rot_scale >= min_rot ;
         rot_scale /= 2)
      for (rx = -whalf ; rx <= whalf ; rx++)
        for (ry = -whalf ; ry <= whalf ; ry++)
          for (rz = -whalf ; rz <= whalf ; rz++)
            for (tx = -whalf ; tx <= whalf ; tx++)
              for (ty = -whalf ; ty <= whalf ; ty++)
                for (tz = -whalf ; tz <= whalf ; tz++)
                {
                  angles[0] = RADIANS(rx*rot_scale) ;
                  angles[1] = RADIANS(ry*rot_scale) ;
                  angles[2] = RADIANS(rz*rot_scale) ;
                  M_rot = MRIangles2RotMat(angles);
                  xtrans = trans_scale*tx ; ytrans = trans_scale*ty ; ztrans = trans_scale*tz ;
                  *MATRIX_RELT(M_trans, 1, 4) = xtrans ;
                  *MATRIX_RELT(M_trans, 2, 4) = ytrans ;
                  *MATRIX_RELT(M_trans, 3, 4) = ztrans ;

                  // Mtest = M_trans * M_ctr * M_rot * inv(M_ctr) * M_min
                  MatrixMultiply(M_ctr_inv, M_min, M_tmp) ;
                  MatrixMultiply(M_rot, M_tmp, M_test) ;
                  MatrixMultiply(M_ctr, M_test, M_tmp) ;
                  MatrixMultiply(M_trans, M_tmp, M_test) ;
                  MatrixFree(&M_rot) ;
                  similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_test, skip, skip/8.0, 0) ;
                  if (similarity > max_similarity)
                  {
                    MatrixCopy(M_test, M_min) ;
                    max_similarity = similarity ;
                    printf("%03d: new max %2.2f found at T=(%2.3f, %2.3f, %2.3f), "
                           "R=(%2.3f, %2.3f, %2.3f)\n",
                           niter+1, max_similarity, 
                           xtrans, ytrans, ztrans,
                           DEGREES(angles[0]), 
                           DEGREES(angles[1]), DEGREES(angles[2]));
                    if (logfp)
                    {
                      fprintf(logfp, "%03d: new max %2.2f found at T=(%2.3f, %2.3f, %2.3f), "
                              "R=(%2.3f, %2.3f, %2.3f)\n",
                              niter+1, max_similarity, 
                              xtrans, ytrans, ztrans,
                              DEGREES(angles[0]), 
                              DEGREES(angles[1]), DEGREES(angles[2]));
                      fflush(logfp) ;
                    }

                    //                    MatrixPrint(Gstdout, M_min) ;
                    if ((++niter % write_iter) == 0)
                      write_snapshot(mris, M_min, outregfile, niter) ;
                    similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_min,skip,skip/8.0,1);
                    good_step = 1 ;
                  }
                }
    if (good_step == 0)
    {
      trans_scale /= 2 ;
      if (trans_scale >= min_trans)
        printf("reducing scale to %2.4f mm, %2.4f deg\n", 
               trans_scale, max_rot/whalf);
    }
    good_step = 0 ;
  } while (trans_scale >= min_trans) ;

  MatrixFree(&M_test) ; MatrixFree(&M_trans) ; MatrixFree(&M_ctr) ; MatrixFree(&M_ctr_inv) ;
  MatrixFree(&M_tmp) ;
  return(max_similarity) ;
}

static double
find_optimal_rotations(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, 
                       MATRIX *M_min, double max_rot, int skip,
                       double (*similarity_func)(MRI_SURFACE *mris,MRI *mri_reg,
                                                 MRI *mri_mask, MATRIX *m, 
                                                 int skip, double scale, int diag))

{
  double   rot_scale, max_similarity, similarity, rx, ry,rz;
  MATRIX   *M_test = NULL, *M_rot, *M_tmp, *M_ctr, *M_ctr_inv  ;
  int      good_step = 0 ;
  double   angles[3] ;

  M_ctr = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(M_ctr, 1,4) = mris->xctr ;
  *MATRIX_RELT(M_ctr, 2,4) = mris->yctr ;
  *MATRIX_RELT(M_ctr, 3,4) = mris->zctr ;
  M_ctr_inv = MatrixInverse(M_ctr, NULL) ;
  M_tmp = MatrixCopy(M_ctr, NULL) ;
  M_test = MatrixCopy(M_ctr, NULL) ;
  max_similarity = similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_min, skip, skip/8.0, 0) ;
#ifdef WHALF
#undef WHALF
#endif
#define WHALF 11
  rot_scale = max_rot/WHALF ;
  do
  {
    for (rx = -WHALF ; rx <= WHALF ; rx++)
      for (ry = -WHALF ; ry <= WHALF ; ry++)
        for (rz = -WHALF ; rz <= WHALF ; rz++)
        {
          angles[0] = RADIANS(rx*rot_scale) ;
          angles[1] = RADIANS(ry*rot_scale) ;
          angles[2] = RADIANS(rz*rot_scale) ;
          M_rot = MRIangles2RotMat(angles);

          // Mtest = M_ctr * M_rot * inv(M_ctr) * M_min
          MatrixMultiply(M_ctr_inv, M_min, M_test) ;
          MatrixMultiply(M_rot, M_test, M_tmp) ;
          MatrixMultiply(M_ctr, M_tmp, M_test) ;
          MatrixFree(&M_rot) ;
          similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_test, skip, skip/8.0, 0) ;
          if (similarity > max_similarity)
          {
            MatrixCopy(M_test, M_min) ;
            max_similarity = similarity ;
            if (logfp)
            {
              fprintf(logfp, "%03d: new max %2.2f found at R=(%2.3f, %2.3f, %2.3f)\n",
                      niter+1, max_similarity, DEGREES(angles[0]), 
                      DEGREES(angles[1]), DEGREES(angles[2])) ;
              fflush(logfp) ;
            }

            printf("%03d: new max %2.2f found at R=(%2.3f, %2.3f, %2.3f)\n",
                   niter+1, max_similarity, DEGREES(angles[0]), 
                   DEGREES(angles[1]), DEGREES(angles[2])) ;
            //            MatrixPrint(Gstdout, M_min) ;
            if ((++niter % write_iter) == 0)
              write_snapshot(mris, M_min, outregfile, niter) ;
            similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_min, skip, skip/8.0, 1) ;
            good_step = 1 ;
          }
        }
    if (good_step == 0)
    {
      rot_scale /= 2 ;
      if (rot_scale >= .01)
        printf("reducing scale to %2.4f deg\n", rot_scale) ;
    }
    good_step = 0 ;
  } while (rot_scale >= .01);

  MatrixFree(&M_test) ; MatrixFree(&M_ctr) ; MatrixFree(&M_ctr_inv) ; MatrixFree(&M_tmp) ;
  return(max_similarity) ;
}

static int
write_snapshot(MRI_SURFACE *mris, MATRIX *m, char *name, int n)
{
  static VECTOR *v1, *v2 ;
  static char fname[STRLEN], fname_only[STRLEN], *ext ;
  VERTEX *v ;
  int    vno ;

  FileNameOnly(name, fname_only) ;
  ext = strrchr(fname_only, '.') ;
  if (ext && *(ext-1) != 'h' && (*(ext-2) != 'r' || *(ext-2) != 'l'))
    *ext = 0 ;
  MRISstoreRipFlags(mris);
  MRISunrip(mris) ;
  if (v1 == NULL)
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = 1.0 ;
  }

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;

  MRISfreeDistsButNotOrig(mris);
    // MRISsetXYZ will invalidate all of these,
    // so make sure they are recomputed before being used again!
  
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
    v2 = MatrixMultiply(m, v1, v2) ;
    MRISsetXYZ(mris,vno, V3_X(v2),V3_Y(v2),V3_Z(v2)) ;
  }

  sprintf(fname, "%s%03d.dat", fname_only,n) ;
  regio_write_surfacexform_to_register_dat(m, fname, mris, mri_reg, subject,
                                           FLT2INT_ROUND) ;
  write_lta(m, fname, mris, mri_reg) ;
  sprintf(fname, "%s%03d", fname_only,n) ;
  printf("writing snapshot to %s\n", fname) ;
  MRISwrite(mris, fname) ;
  MRISrestoreRipFlags(mris);
  sprintf(fname, "%s%03d.patch", fname_only,n) ;
  MRISwritePatch(mris, fname) ;
  MRISunrip(mris) ;
  
  if (pial_fname)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      V3_X(v1) = v->pialx ; V3_Y(v1) = v->pialy ; V3_Z(v1) = v->pialz ;
      v2 = MatrixMultiply(m, v1, v2) ;
      MRISsetXYZ(mris,vno, V3_X(v2),V3_Y(v2),V3_Z(v2)) ;
    }
    
    sprintf(fname, "%s_pial%03d", fname_only,n) ;
    printf("writing snapshot to %s\n", fname) ;
    MRISwrite(mris, fname) ;
    sprintf(fname, "%s%03d.patch", fname_only,n) ;
    MRISwritePatch(mris, fname) ;
  }

#if 0
  sprintf(fname, "%s%03d.mgz", fname_only,n) ;
  printf("writing overlay to %s\n", fname) ;
  MRISwriteCurvature(mris, fname) ;
#endif

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreRipFlags(mris);
  return(NO_ERROR) ;
}


#define NPARMS_RIGID (6)
static MRI_SURFACE *Gmris ;
static MRI *Gmri_reg, *Gmri_mask ;
static int Gskip ;
static double (*Gsimilarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                  MRI *mri_mask, MATRIX *m, int skip, 
                                  double scale, int diag) ;

static int
powell_minimize_rigid(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *mat, int skip,
                      double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                MRI *mri_mask, MATRIX *m, int skip, double scale, int diag))
{
  float *p, **xi, fret, fstart ;
  int    r, c, iter, diag, old_write ;
  double xr, yr, zr, xt, yt, zt;
  const char *sname ;

  if (similarity_func == mrisRegistrationCNRSimilarity)
    sname = "CNR" ;
  else if (similarity_func == mrisRegistrationGradientNormalSimilarity)
    sname = "DOT" ;
  else
    sname = "UNK" ;
  old_write = write_iter ; write_iter = 1 ;

  // extract rigid body parameters from matrix
  MatrixToRigidParameters(mat, &xr, &yr, &zr, &xt, &yt, &zt) ;
  printf("initial rigid parms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
         DEGREES(xr), DEGREES(yr), DEGREES(zr), xt, yt, zt) ;
  if (logfp)
  {
    fprintf(logfp, "initial rigid parms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
            DEGREES(xr), DEGREES(yr), DEGREES(zr), xt, yt, zt) ;
    fflush(logfp) ;
  }
  if ((niter == 0) && (write_iter >= 0))
    write_snapshot(mris, mat, outregfile, niter) ;
  p = vector(1, NPARMS_RIGID) ;
  xi = matrix(1, NPARMS_RIGID, 1, NPARMS_RIGID) ;
  p[1] = xr ;
  p[2] = yr ;
  p[3] = zr ;
  p[4] = xt ;
  p[5] = yt ;
  p[6] = zt ;

  Gmris = mris ;
  Gmri_reg = mri_reg ;
  Gmri_mask = mri_mask ;
  Gskip = skip ;
  for (r = 1 ; r <= NPARMS_RIGID ; r++) {
    for (c = 1 ; c <= NPARMS_RIGID ; c++) {
      xi[r][c] = r == c ? 1 : 0 ;
    }
  }

  diag = Gdiag ;
  Gdiag |= DIAG_VERBOSE ;
  Gsimilarity_func = similarity_func ;
  OpenPowell(p, xi, NPARMS_RIGID, TOL, &iter, &fret, compute_powell_rigid_sse);
  MatrixFromRigidParameters(mat, p[1],p[2],p[3],p[4],p[5],p[6]) ;
  compute_powell_rigid_sse(p) ;
  if (logfp)
  {
    fprintf(logfp, "%3.3d: best alignment after %s powell: "
            "%2.5f (%d steps)\n\tparms = R:(%2.4f, %2.4f, %2.4f) + T:(%2.2f, %2.2f, %2.2f)\n",
            niter+1,sname, fret, iter,
            DEGREES(p[1]), DEGREES(p[2]), DEGREES(p[3]), p[4], p[5], p[6]) ;
    fflush(logfp) ;
  }

  printf("%3.3d: best alignment after %s powell: "
         "%2.5f (%d steps)\n\tparms = R:(%2.4f, %2.4f, %2.4f) + T:(%2.2f, %2.2f, %2.2f)\n",
         niter+1,sname, fret, iter,
         DEGREES(p[1]), DEGREES(p[2]), DEGREES(p[3]), p[4], p[5], p[6]) ;
  if ((++niter % write_iter) == 0)
    write_snapshot(mris, mat, outregfile, niter) ;
  Gdiag = diag ;
  do {
    for (r = 1 ; r <= NPARMS_RIGID ; r++) {
      for (c = 1 ; c <= NPARMS_RIGID ; c++) {
        xi[r][c] = r == c ? 1 : 0 ;
      }
    }

    fstart = fret ;
    OpenPowell(p, xi, NPARMS_RIGID, TOL, &iter, &fret, compute_powell_rigid_sse);
    MatrixFromRigidParameters(mat, p[1],p[2],p[3],p[4],p[5],p[6]) ;
#if 0
    *MATRIX_RELT(mat, 4, 1) = 0.0 ;
    *MATRIX_RELT(mat, 4, 2) = 0.0 ;
    *MATRIX_RELT(mat, 4, 3) = 0.0 ;
    *MATRIX_RELT(mat, 4, 4) = 1.0 ;
#endif
    if (logfp)
    {
      fprintf(logfp, "%3.3d: best alignment after %s powell: "
              "%2.5f (%d steps)\n\tparms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
              niter+1,sname, fret, iter,
           DEGREES(p[1]), DEGREES(p[2]), DEGREES(p[3]), p[4], p[5], p[6]) ;
      fflush(logfp) ;
    }
    printf("%3.3d: best alignment after %s powell: "
           "%2.5f (%d steps)\n\tparms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
           niter+1,sname, fret, iter,
           DEGREES(p[1]), DEGREES(p[2]), DEGREES(p[3]), p[4], p[5], p[6]) ;
    if ((++niter % write_iter) == 0)
      write_snapshot(mris, mat, outregfile, niter) ;
  } while (fret < fstart) ;

  free_matrix(xi, 1, NPARMS_RIGID, 1, NPARMS_RIGID) ;
  free_vector(p, 1, NPARMS_RIGID) ;
  write_iter = old_write ;
  return(NO_ERROR) ;
}

static float
compute_powell_rigid_sse(float *p) 
{
  static MATRIX *mat = NULL ;
  double  similarity ;

  if (mat == NULL)
    mat = MatrixAlloc(4, 4, MATRIX_REAL) ;

  MatrixFromRigidParameters(mat, p[1],p[2],p[3],p[4],p[5],p[6]) ;
  similarity = (*Gsimilarity_func)(Gmris, Gmri_reg, Gmri_mask, mat, Gskip, 
                                   Gskip/8, 0) ;
  return(-(float)similarity) ;
}
#if 0
static int
write_register_dat(MATRIX *B, char *fname, MRI_SURFACE *mris, MRI *mri, char *subject)
{
  MATRIX *Ta, *Sa, *invTa, *A, *R, *S, *invS, *T, *m1, *m2 ;
  MRI *mri_surf = MRIallocHeader(mris->vg.width, mris->vg.height, 
                                 mris->vg.depth, MRI_UCHAR) ;

  MRIcopyVolGeomToMRI(mri_surf, &mris->vg) ;

  T = MRIxfmCRS2XYZtkreg(mri) ;
  S = MRIgetVoxelToRasXform(mri) ;
  invS = MatrixInverse(S, NULL) ;
  Ta = MRIxfmCRS2XYZtkreg(mri_surf);
  Sa = MRIgetVoxelToRasXform(mri_surf);
  invTa = MatrixInverse(Ta,NULL);
  A  = MatrixMultiply(Sa,invTa, NULL);
  
  m1 = MatrixMultiply(A, B, NULL) ;
  m2 = MatrixMultiply(invS, m1, NULL) ;
  R = MatrixMultiply(T, m2, NULL) ;
  regio_write_register(fname,subject,mri_reg->xsize,
                       mri_reg->zsize,1,R,FLT2INT_ROUND);
  write_lta(R, fname, mris, mri_reg) ;
  MatrixFree(&A) ; MatrixFree(&Ta) ; MatrixFree(&Sa) ; MatrixFree(&invTa) ;
  MatrixFree(&R) ; MatrixFree(&m1) ; MatrixFree(&m2) ; MatrixFree(&S) ; MatrixFree(&invS);
  MatrixFree(&T) ;
  MRIfree(&mri_surf) ;
  return(NO_ERROR) ;
}



#endif
static MRI *
mri_surface_edges(MRI *mri, MRI *mri_grad, float normal_dist, float tangent_dist, int dir, 
                  MRI *mri_edge)
{
  int     x, y, z, nout, nin ;
  double  dx, dy, dz, d, norm, mean_in, mean_out, var_in, var_out, step_size, 
          x1, y1, z1,val, t1x, t1y, t1z, t2x, t2y, t2z, d1, d2  ;

  if (mri_edge == NULL)
  {
    mri_edge = MRIalloc(mri->width, mri->height, mri->depth, MRI_FLOAT) ;
    MRIcopyHeader(mri, mri_edge) ;
  }

  step_size = 0.5 ;         // in voxels
  normal_dist /= mri->xsize ;  // convert to voxels from mm
  for (x = 0 ; x < mri->width ; x++)
  {
    for (y = 0 ; y < mri->height ; y++)
    {
      for (z = 0 ; z < mri->depth ; z++)
      {
        if (x == Gx && y == Gy && z == Gz)
          DiagBreak() ;
        dx = MRIgetVoxVal(mri_grad, x, y, z, 0) ;
        dy = MRIgetVoxVal(mri_grad, x, y, z, 1) ;
        dz = MRIgetVoxVal(mri_grad, x, y, z, 2) ;
        norm = sqrt(dx*dx + dy*dy + dz*dz) ;
        if (FZERO(norm))
          continue ;
        dx /= norm ; dy /= norm ; dz /= norm ;
        if (!FEQUAL(dx, dy))  // find non-colinear vector swap x and y
        { 
          t2x = dy ; t2y = dx ; t2z = dz ;
        }
        else if (!FEQUAL(dy, dz))  // swap y and z
        {
          t2x = dx ; t2y = dz ; t2z = dy ;
        }
        else  // swap x and z
        {
          t2x = dz ; t2y = dy ; t2z = dx ;
        }

        // cross product for 1st tangent
        t1x=dy*t2z-dz*t2y;  t1y=dz*t2x-dx*t2z; t1z=dx*t2y-dy*t2x;
        norm = sqrt(t1x*t1x + t1y*t1y + t1z*t1z ) ;
        t1x /= norm ; t1y /= norm; ; t1z /= norm ;

        // 2nd tangent
        t2x=dy*t1z-dz*t1y;  t2y=dz*t1x-dx*t1z; t2z=dx*t1y-dy*t1x;
        norm = sqrt(t2x*t2x + t2y*t2y + t2z*t2z ) ;
        t2x /= norm ; t2y /= norm; ; t2z /= norm ;

        mean_in = mean_out = var_in = var_out = 0.0 ;
        nout = nin = 0 ;
        for (d1 = -tangent_dist/mri->xsize  ; d1 <= tangent_dist/mri->xsize ; d1++)
          for (d2 = -tangent_dist/mri->xsize  ; d2 <= tangent_dist/mri->xsize ; d2++)
            for (d = step_size ; d <= normal_dist ; d += step_size)
            {
              x1 = x + d*dx + d1*t1x + d2*t2x ; 
              y1 = y + d*dy + d1*t1y + d2*t2y; 
              z1 = z + d*dz + d1*t1z + d2*t2z ; // outside
              if (MRIindexNotInVolume(mri, x1, y1, z1) != 1)
              {
                MRIsampleVolume(mri, x1, y1, z1, &val) ;
                mean_out += val ; var_out += (val*val) ;
                nout++ ;
              }

              x1 = x - d*dx + d1*t1x + d2*t2x ; 
              y1 = y - d*dy + d1*t1y + d2*t2y; 
              z1 = z - d*dz + d1*t1z + d2*t2z ; // outside
              if (MRIindexNotInVolume(mri, x1, y1, z1) != 1)
              {
                MRIsampleVolume(mri, x1, y1, z1, &val) ;
                mean_in += val ; var_in += (val*val) ;
                nin++ ;
              }
            }
        if (nin < 3 || nout < 3)
          continue ;
        mean_out /= nout ; mean_in /= nin ;
        var_out = var_out/nout - mean_out*mean_out ;
        var_in = var_in/nin - mean_in*mean_in ;
        val = dir*(mean_out-mean_in) / (pow((mean_out+mean_in)/2, 0.333)*sqrt(var_out+var_in)) ;
        if (val < 0)
          val = 0 ;
        if (val > 100)
          DiagBreak() ;
        MRIsetVoxVal(mri_edge, x, y, z, 0, val) ;
      }
    }
  }

  return(mri_edge) ;
}

static double
mrisRegistrationOverlapSimilarity(MRI_SURFACE *mris, MRI *mri_reg, 
                                  MRI *mri_mask, MATRIX *m, int skip, 
                                  double scale, int diag)
{
  double    similarity, xv, yv, zv, white_val ;
  int       vno, num, num_in_fov = 0 ;
  VERTEX    *v ;
  static VECTOR *v1 = NULL, *v2 = NULL ;

  if (v1 == NULL)
  {
    v1 = VectorAlloc(4, MATRIX_REAL) ;
    VECTOR_ELT(v1, 4) = 1.0 ;
  }
  skip++ ;
  for (similarity = 0.0, num = vno = 0 ; vno < mris->nvertices ; vno += skip)
  {
    v = &mris->vertices[vno] ;
    v->marked = 0 ;
    if (v->ripflag)
      continue ;
    num++ ;
    V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
    v2 = MatrixMultiply(m, v1, v2) ;
    MRISsurfaceRASToVoxelCached(mris, mri_reg, 
                                V3_X(v2), V3_Y(v2), V3_Z(v2), &xv, &yv, &zv);
    if (MRIindexNotInVolume(mri_reg, xv, yv, zv) || 
        MRIgetVoxVal(mri_mask, xv, yv, zv, 0) == 0)
      white_val = mri_reg->outside_val ;
    else  // in the volume
    {
      num_in_fov++ ;
      white_val = MRIgetVoxVal(mri_reg, xv, yv, zv, 0) ;
      v->marked = 1 ;
    }
    similarity += white_val ;
  }

  if (diag)
    printf("num_in_fov = %d\n", num_in_fov) ;

  if (num == 0)
    num = 1;
  return(similarity) ;  
}
static int
write_lta(MATRIX *m, char *fname_from_caller, MRI_SURFACE *mris, MRI *mri_reg)
{
  char *dot, fname[STRLEN] ;
  LTA  *lta ;

  strcpy(fname, fname_from_caller) ;
  dot = strrchr(fname, '.') ;
  strcpy(dot+1, "lta") ;

  lta = LTAalloc(1, NULL) ;

  MatrixCopy(m, lta->xforms[0].m_L) ;

  LTAwriteEx(lta, fname) ;
  return(NO_ERROR) ;
}

