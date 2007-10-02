/**
 * @file  mris_register_to_volume.c
 * @brief program for computing/optimizing registration of a surface to a volume
 *        
 *
 * Program to compute a rigid alignment between a surface and a volume by maximizing the gradient
 * magnitude across the gray/white boundary, divided by it's variance
 */
/*
 * Original Author: Greg Grev
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2007/10/02 18:42:27 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */


/*
BEGINUSAGE --------------------------------------------------------------

 mri_register_to_volume
  --reg regfile
  --mov fvol
  --surf surface   : surface to read in
  --pial pial surface name   : surface to read in

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

#ifdef X
#undef X
#endif

double *GetCosts(MRI *mri_reg, MRI *seg, MATRIX *R0, MATRIX *R, double *p, double *costs);
int Min1D(MRI *mri_reg, MRI_SURFACE *mris, MATRIX *R, double *p, char *costfile, double *costs);

// For some reason, this does not seemed to be defined in math.h
double round(double x);
static int write_snapshot(MRI_SURFACE *mris, MATRIX *R0, char *fname, int n) ;

static float compute_powell_rigid_sse(float *p) ;
static int powell_minimize_rigid(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *mat, int skip,
                                 double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                           MRI *mri_mask, MATRIX *m, int skip, double scale, int diag));

static int  parse_commandline(int argc, char **argv);
static double mrisRegistrationCNRSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag);
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
static int istringnmatch(char *str1, char *str2, int n);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mris_register_to_volume.c,v 1.2 2007/10/02 18:42:27 fischl Exp $";
char *Progname = NULL;

static int debug = 0, gdiagno = -1;
static int ndilates = 2 ;
static char *patch_fname = NULL ;
static char *vol_fname=NULL;
static char *surf_fname = NULL ;
static char *pial_fname = NULL ;
static char *regfile=NULL;
static char *outregfile=NULL;
static char *interpmethod = "trilinear";
static int   interpcode = 0;
static int   sinchw;
static double   max_trans = 200 ;
static double   max_rot = 20 ;  // degrees
static int max_skip = 32 ;
static int min_skip = 8 ;

static double max_sigma = 2 ;
static double min_sigma = .5 ;


static MRI *mri_reg, *mri_grad = NULL ;

MATRIX *R0;

char *SUBJECTS_DIR=NULL;
char *subject = "unknown";

static float ipr, bpr, intensity;
static int float2int, nargs;

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
                                           double max_trans, double max_rot, int skip,
                                           double (*similarity_func)
                                           (MRI_SURFACE *mris, MRI *mri_reg, 
                                            MRI *mri_mask, MATRIX *m, int skip,
                                            double scale, int diag)) ;

static double (*similarity_func)
     (MRI_SURFACE *mris, MRI *mri_reg, 
      MRI *mri_mask, MATRIX *m, int skip, double scale, 
      int diag) = mrisRegistrationGradientSimilarity ;

/*---------------------------------------------------------------*/
int
main(int argc, char **argv) 
{
  char          cmdline[CMD_LINE_LEN] ;
  MRI_SURFACE   *mris ;
  int           skip ;
  MRI           *mri_kernel, *mri_smooth, *mri_mask ;
  double        sigma ;
#if 0
  MRI           *mri_dist;
  HISTOGRAM     *h ;
#endif

  make_cmd_version_string(argc, argv,
                          "$Id: mris_register_to_volume.c,v 1.2 2007/10/02 18:42:27 fischl Exp $",
                          "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option(argc, argv,
                                "$Id: mris_register_to_volume.c,v 1.2 2007/10/02 18:42:27 fischl Exp $",
                                "$Name:  $");
  if(nargs && argc - nargs == 1) exit (0);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

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
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
  }
    

  if (patch_fname)
  {
    MRISsaveVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISreadPatchNoRemove(mris, patch_fname) ;
    MRISdilateRipped(mris, ndilates) ;
    printf("using surface patch with %d vertices\n", MRISvalidVertices(mris)) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  }
  MRISresetNeighborhoodSize(mris, 2) ;
  
  printf("Loading mov %s\n", vol_fname);
  mri_reg = MRIread(vol_fname);
  if (mri_reg == NULL) 
    exit(1);
  {
    MRI *mri_tmp ;
    mri_tmp = MRIextract(mri_reg, NULL, mask_size, mask_size, mask_size,
                         (mri_reg->width-2*mask_size),
                         (mri_reg->height-2*mask_size),
                         (mri_reg->depth-2*mask_size)) ;
    MRIfree(&mri_reg) ;
    mri_reg = mri_tmp ;
  }

  mri_mask = MRIclone(mri_reg, NULL) ;
  MRIsetValues(mri_mask, 1) ;
  MRIeraseBorderPlanes(mri_mask, 1) ;

#if 0
  h = MRIhistogram(mri_mag, 100) ;
  HISTOplot(h, "h.plt") ;
  HISTOfree(&h) ;
  MRIbinarize(mri_mag, mri_mag, 1.0, 0, 1.0) ;
  MRIwrite(mri_mag, "b.mgz") ;
  mri_dist = 
    MRIdistanceTransform(mri_mag, NULL, 1, 
                         10,
                         DTRANS_MODE_SIGNED) ;
  MRIwrite(mri_dist, "d.mgz") ;
#endif
                                  
  // Cropping reduces the size of the target volume down to the 
  // voxels that really matter. This can greatly increase the speed
  // (factor of 2-3). Requires that the registration matrix be
  // recomputed for the smaller volume.
#if 0
  if(DoCrop){
    printf("Cropping\n");

    // Prepare to adjust the input reg matrix
    //    Ttemp    = MRIxfmCRS2XYZtkreg(regseg); // Vox-to-tkRAS Matrices
    //    Stemp    = MRIxfmCRS2XYZ(regseg,0); // Vox-to-ScannerRAS Matrices
    invStemp = MatrixInverse(Stemp,NULL);

    MRIboundingBox(regseg, 0.5, &box);
    printf("BBbox start: %d %d %d, delta = %d %d %d\n",
           box.x,box.y,box.z,box.dx,box.dy,box.dz);
    mritmp = MRIcrop(regseg, box.x, box.y, box.z, box.x+box.dx, box.y+box.dy, box.z+box.dz);
    MRIfree(&regseg);
    regseg = mritmp;

    Tcrop  = MRIxfmCRS2XYZtkreg(regseg); // Vox-to-tkRAS Matrices
    invTcrop = MatrixInverse(Tcrop,NULL);
    Scrop  = MRIxfmCRS2XYZ(regseg,0); // Vox-to-ScannerRAS Matrices

    // Now adjust input reg
    // Rc = R*T*inv(S)*Sc*inv(Tc)
    R0 = MatrixMultiply(R0,Ttemp,R0);
    R0 = MatrixMultiply(R0,invStemp,R0);
    R0 = MatrixMultiply(R0,Scrop,R0);
    R0 = MatrixMultiply(R0,invTcrop,R0);
    
    MatrixFree(&Ttemp);
    MatrixFree(&Stemp);
    MatrixFree(&invStemp);
    MatrixFree(&Tcrop);
    MatrixFree(&invTcrop);
    MatrixFree(&Scrop);

    //MRIwrite(regseg,"regsegcrop.mgz");
    //regio_write_register("crop.reg",subject,mri_reg->xsize,
    //		 mri_reg->zsize,1,R0,FLT2INT_ROUND);
    //exit(1);
  }

  if(AddNoise){
    // Seed the random number generator just in case
    if (SynthSeed < 0) SynthSeed = PDFtodSeed();
    srand48(SynthSeed);
    printf("Adding noise, Seed = %d, stddev = %lf\n",SynthSeed,NoiseStd);
    noise = MRIrandn(mri_reg->width,mri_reg->height,mri_reg->depth,mri_reg->nframes,
                     0,NoiseStd,NULL);
    mri_reg = MRIadd(mri_reg, noise, mri_reg);
  }


  // Vox-to-tkRAS Matrices
  Tin      = MRIxfmCRS2XYZtkreg(mri_reg);
  invTin   = MatrixInverse(Tin,NULL);
  Ttemp    = MRIxfmCRS2XYZtkreg(regseg);
  invTtemp = MatrixInverse(Ttemp,NULL);
  R = MatrixCopy(R0,NULL);
#endif

#if 0
  write_snapshot(mris, R0, outregfile, niter) ;
  skip = 64 ;
  similarity_func = mrisRegistrationDistanceSimilarity;
  mri_reg = mri_dist ;
  powell_minimize_rigid(mris, mri_dist, mri_mask, R0, skip, mrisRegistrationDistanceSimilarity);
#endif
  for (skip = max_skip ; skip >= min_skip;  skip /= 2)
  {
    for (sigma = max_sigma ; sigma >= min_sigma ; sigma /= 2)
    {
      printf("---------------- skip = %d, sigma = %2.1f ---------------------\n", skip, sigma);
      printf("computing gradient at scale %2.3f...\n",sigma) ;
      mri_kernel = MRIgaussian1d(sigma/mri_reg->xsize, -1) ;
      mri_smooth = MRIconvolveGaussian(mri_reg, NULL, mri_kernel) ;
      MRIfree(&mri_kernel) ;
      mri_grad = MRIsobel(mri_smooth, mri_grad, NULL) ;
      MRIeraseBorderPlanes(mri_grad, 1) ;
      MRIfree(&mri_smooth) ; 
      find_optimal_translations(mris, mri_reg, mri_mask,R0, max_trans, skip, similarity_func) ;
      find_optimal_rotations(mris, mri_reg, mri_mask, R0, max_rot, skip, similarity_func) ;
      find_optimal_rigid_alignment(mris, mri_reg, mri_mask, R0, max_trans/25, 
                                   max_rot/2, skip, similarity_func);
      powell_minimize_rigid(mris, mri_reg, mri_mask, R0, skip, similarity_func);
      if (similarity_func != mrisRegistrationCNRSimilarity)
        powell_minimize_rigid(mris, mri_reg, mri_mask, R0, skip, mrisRegistrationCNRSimilarity);
      printf("saving current registration to %s\n", outregfile) ;
      regio_write_register(outregfile,subject,mri_reg->xsize,
                           mri_reg->zsize,1,R0,FLT2INT_ROUND);
    }
    if (skip == 0)
      break ;
  }

#if 0
  if(1)
  {
    // 1D minimization
    TimerStart(&mytimer) ;
    nth = 0;
    for(n=0; n < 3; n++){
      printf("n = %d --------------------------------------\n",n);
      R = MatrixCopy(R0,NULL);
      nth += Min1D(mri_reg, mris, R, p, SegRegCostFile, costs);
      printf("\n");
    }
    secCostTime = TimerStop(&mytimer)/1000.0 ;
    
    printf("Min cost was %lf\n",costs[7]);
    printf("Number of iterations %5d in %lf sec\n",nth,secCostTime);
    printf("Parameters at optimum\n");
    printf("%7.3lf %7.3lf %7.3lf %6.3lf %6.3lf %6.3lf \n",
           p[0],p[1],p[2],p[3],p[4],p[5]);
    printf("Costs at optimum\n");
    printf("%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
    printf("%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
    printf("%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
    printf("\n");
    
    printf("Reg at min cost was \n");
    MatrixPrint(stdout,R);
    printf("\n");
    
    if(outregfile){
      printf("Writing optimal reg to %s \n",outregfile);
      regio_write_register(outregfile,subject,mri_reg->xsize,
                           mri_reg->zsize,1,R,FLT2INT_ROUND);
    }
    exit(0);
  }

  ntot = ntx*nty*ntz*nax*nay*naz;
  printf("Staring loop over %d\n",ntot);
  nth = 0;
  mincost = 10e10;
  for(nthtx = 0; nthtx < ntx; nthtx++){
    tx = txlist[nthtx];
    for(nthty = 0; nthty < nty; nthty++){
      ty = tylist[nthty];
      for(nthtz = 0; nthtz < ntz; nthtz++){
        tz = tzlist[nthtz];
        for(nthax = 0; nthax < nax; nthax++){
          ax = axlist[nthax];
          for(nthay = 0; nthay < nay; nthay++){
            ay = aylist[nthay];
            for(nthaz = 0; nthaz < naz; nthaz++){
              az = azlist[nthaz];
              nth ++;

              p[0] = tx;
              p[1] = ty;
              p[2] = tz;
              p[3] = ax;
              p[4] = ay;
              p[5] = az;
              TimerStart(&mytimer) ;
              GetCosts(mri_reg, regseg, R0, R, p, costs);
              secCostTime = TimerStop(&mytimer)/1000.0 ;

              // write costs to file
              fp = fopen(SegRegCostFile,"a");
              fprintf(fp,"%7.3lf %7.3lf %7.3lf ",tx,ty,tz);
              fprintf(fp,"%6.3lf %6.3lf %6.3lf ",ax,ay,az);
              fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
              fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
              fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
              fprintf(fp,"\n");
              fclose(fp);

              fp = stdout;
              fprintf(fp,"%5d ",nth);
              fprintf(fp,"%7.3lf %7.3lf %7.3lf ",tx,ty,tz);
              fprintf(fp,"%6.3lf %6.3lf %6.3lf ",ax,ay,az);
              fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
              fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
              fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
              if(DoProfile) fprintf(fp,"%4.2lf ",secCostTime);
              printf("\n");
              fflush(stdout);

              if(mincost > costs[7]){
                mincost = costs[7];
                Rmin = MatrixCopy(R,Rmin);
                if(outregfile){
                  regio_write_register(outregfile,subject,mri_reg->xsize,
                                       mri_reg->zsize,1,Rmin,FLT2INT_ROUND);
                }
              }

              // clean up
            }
          }
        }
      }
    }
  }
  printf("min cost was %lf\n",mincost);
#endif

  printf("Reg at min cost was \n");
  MatrixPrint(stdout,R0);
  printf("\n");
  
  if(outregfile){
    printf("Writing optimal reg to %s \n",outregfile);
    regio_write_register(outregfile,subject,mri_reg->xsize,
                         mri_reg->zsize,1,R0,FLT2INT_ROUND);
  }

  printf("\n");
  printf("mris_register_to_volume done\n");

  return(0);
}


/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  int err,nv,n;
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
    else if (!strcasecmp(option, "--profile"))  DoProfile = 1;
    else if (istringnmatch(option, "--reg",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      err = regio_read_register(regfile, &subject, &ipr, &bpr,
                                &intensity, &R0, &float2int);
      if (err) exit(1);
      nargsused = 1;
    } else if (istringnmatch(option, "--out-reg",0)) {
      if (nargc < 1) argnerr(option,1);
      outregfile = pargv[0];
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
    } else if (!strcasecmp(option, "--CNR")) {
      similarity_func = mrisRegistrationCNRSimilarity ;
      nargsused = 0;
    } else if (!strcasecmp(option, "--DOT")) {
      similarity_func = mrisRegistrationGradientNormalSimilarity ;
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
printf("  --reg regfile\n");
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
printf("\n");
printf("  --out-reg outreg : reg at lowest cost (updated continuously)\n");
printf("\n");

}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n%s\n\n",vcid);
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
  printf("%s\n", vcid) ;
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
static int istringnmatch(char *str1, char *str2, int n) {
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
            mag, gm_mean, wm_mean, sample_dist ;
  int       vno, num, num_in_fov = 0 ;
  VERTEX    *v ;
  static VECTOR *v1 = NULL, *v2 = NULL ;
  VERTEX    *vn ;
  double    gm_var, wm_var, contrast, noise, total_contrast, total_noise ;
  int       n, num_nbrs ;

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
    v = &mris->vertices[vno] ;
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
      MRIsampleVolumeDerivativeScale(mri_reg, xv, yv, zv, nx, ny, nz, &grad, .5/mri_reg->xsize) ;
      MRIsampleVolume(mri_reg, xv-sample_dist*nx, yv-sample_dist*ny, zv-sample_dist*nz, &wm_mean) ;
      MRIsampleVolume(mri_reg, xv+sample_dist*nx, yv+sample_dist*ny, zv+sample_dist*nz, &gm_mean) ;
      v->val = wm_mean ;
      v->val2 = gm_mean ;
      v->marked = 1 ;
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
    v = &mris->vertices[vno] ;
    v->curv = 0 ;
    if (v->ripflag || v->marked == 0)
      continue ;
    gm_mean = v->val2 ; gm_var = v->val2*v->val2 ;
    wm_mean = v->val ;  wm_var = v->val*v->val ;
    for (n = 0, num_nbrs = 1 ; n < v->vtotal ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
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
          similarity += contrast*contrast / noise ;
        }
        v->curv = contrast*contrast  ;  // for diagnostics
      }
      else
        DiagBreak() ;
    }
  }

  if (diag)
    printf("num_in_fov = %d\n", num_in_fov) ;
  return(similarity) ;
}

static double
mrisRegistrationGradientSimilarity(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *m, int skip, double scale, int diag)
{
  double    similarity, grad, grad_var, grad_mean, xv, yv, zv, nx, ny, nz, 
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
  double    similarity, dot, xv, yv, zv, nx, ny, nz, mag, dx, dy, dz ;
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
      dot = 0 ;
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
      dot = nx*dx + ny*dy + nz*dz ;
      if (dot < 0)
        dot = 0 ;
    }
    v->val = dot ;
    similarity += dot ;

    if (pial_fname)
    {
      V3_X(v1) = v->pialx ; V3_Y(v1) = v->pialy ; V3_Z(v1) = v->pialz ;
      v2 = MatrixMultiply(m, v1, v2) ;
      MRISsurfaceRASToVoxelCached(mris, mri_reg, 
                                  V3_X(v2), V3_Y(v2), V3_Z(v2), &xv, &yv, &zv);
      if (MRIindexNotInVolume(mri_reg, xv, yv, zv) || 
          MRIgetVoxVal(mri_mask, xv, yv, zv, 0) == 0)
        dot = 0 ;
      else  // in the volume
      {
        num_in_fov++ ;
        V3_X(v1) = v->pialx+v->nx ; V3_Y(v1) = v->pialy+v->ny ; 
        V3_Z(v1) = v->pialz+v->nz ;
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
        dot = -1*(nx*dx + ny*dy + nz*dz) ; // pial contrast should be inverse
        if (dot < 0)
          dot = 0 ;
      }
      v->val = dot ;
      similarity += dot ;
    }
  }

  if (diag)
    printf("num_in_fov = %d\n", num_in_fov) ;

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
find_optimal_translations(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *R0, double max_trans, 
                          int skip,
                          double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                    MRI *mri_mask, MATRIX *m, int skip, double scale, int diag))
{
  double   tx, ty, tz, scale, max_similarity, similarity, xtrans, ytrans, ztrans ;
  MATRIX   *m_trans, *m = NULL ;
  int      good_step = 0 ;

  if (niter == 0 && write_iter > 0)
    write_snapshot(mris, R0, outregfile, niter++) ;
  m_trans = MatrixIdentity(4, NULL) ;

  max_similarity = similarity = 
    (*similarity_func)(mris, mri_reg, mri_mask, R0, skip, skip/8.0, 1) ;
#ifdef WHALF
#undef WHALF
#endif
#define WHALF 9
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
          m = MatrixMultiply(m_trans, R0, m) ;
          similarity = (*similarity_func)(mris, mri_reg,  mri_mask, m, skip, skip/8.0, 0) ;
          if (similarity > max_similarity)
          {
            MatrixCopy(m, R0) ;
            max_similarity = similarity ;
            //            MatrixPrint(Gstdout, R0) ;
            if ((++niter % write_iter) == 0)
              write_snapshot(mris, R0, outregfile, niter) ;
            printf("%03d: new max %2.2f found at T = (%2.1f, %2.1f, %2.1f)\n",
                   niter, max_similarity, xtrans, ytrans, ztrans) ;
            similarity = (*similarity_func)(mris, mri_reg,  mri_mask, m, skip,skip/8.0, 1) ;
            good_step = 1 ;
          }
        }
    if (good_step == 0)
    {
      scale /= 2 ;
      printf("reducing scale to %2.4f mm\n", scale) ;
    }
    good_step = 0 ;
  } while (scale >= mri_reg->xsize/2);

  MatrixFree(&m) ; MatrixFree(&m_trans) ;
  return(max_similarity) ;
}

static double
find_optimal_rigid_alignment(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *M_min, double max_trans, 
                             double max_rot, int skip,
                             double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                       MRI *mri_mask, MATRIX *m, int skip, double scale, int diag))

{
  double   tx, ty, tz, rot_scale, trans_scale, max_similarity, similarity, xtrans, ytrans, ztrans,
    rx, ry, rz ;
  MATRIX   *M_trans, *M_test = NULL, *M_rot, *M_tmp, *M_ctr, *M_ctr_inv  ;
  int      good_step = 0 ;
  double   angles[3] ;

  M_trans = MatrixIdentity(4, NULL) ;
  M_ctr = MatrixIdentity(4, NULL) ;
  *MATRIX_RELT(M_ctr, 1,4) = mris->xctr ;
  *MATRIX_RELT(M_ctr, 2,4) = mris->yctr ;
  *MATRIX_RELT(M_ctr, 3,4) = mris->zctr ;
  M_ctr_inv = MatrixInverse(M_ctr, NULL) ;
  M_tmp = MatrixCopy(M_trans, NULL) ;
  M_test = MatrixCopy(M_trans, NULL) ;
  max_similarity = similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_min, skip, skip/8.0, 0) ;
#ifdef WHALF
#undef WHALF
#endif
#define WHALF 3
  trans_scale = max_trans/WHALF ; 
  rot_scale = max_rot/WHALF ;
  do
  {
    for (rot_scale = max_rot/WHALF ; 
         rot_scale > .1 ;
         rot_scale /= 2)
      for (rx = -WHALF ; rx <= WHALF ; rx++)
        for (ry = -WHALF ; ry <= WHALF ; ry++)
          for (rz = -WHALF ; rz <= WHALF ; rz++)
            for (tx = -WHALF ; tx <= WHALF ; tx++)
              for (ty = -WHALF ; ty <= WHALF ; ty++)
                for (tz = -WHALF ; tz <= WHALF ; tz++)
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
                    if ((++niter % write_iter) == 0)
                      write_snapshot(mris, M_min, outregfile, niter) ;
                    printf("%03d: new max %2.2f found at T=(%2.3f, %2.3f, %2.3f), "
                           "R=(%2.3f, %2.3f, %2.3f)\n",
                           niter, max_similarity, 
                           xtrans, ytrans, ztrans,
                           DEGREES(angles[0]), 
                           DEGREES(angles[1]), DEGREES(angles[2]));
                      //                    MatrixPrint(Gstdout, M_min) ;
                    similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_min, skip, skip/8.0, 1) ;
                    good_step = 1 ;
                  }
                }
    if (good_step == 0)
    {
      trans_scale /= 2 ;
      printf("reducing scale to %2.4f mm, %2.4f deg\n", trans_scale, max_rot/WHALF);
    }
    good_step = 0 ;
  } while (trans_scale >= mri_reg->xsize/2) ;

  MatrixFree(&M_test) ; MatrixFree(&M_trans) ; MatrixFree(&M_ctr) ; MatrixFree(&M_ctr_inv) ;
  MatrixFree(&M_tmp) ;
  return(max_similarity) ;
}

static double
find_optimal_rotations(MRI_SURFACE *mris, MRI *mri_reg, MRI *mri_mask, MATRIX *M_min, 
                       double max_rot, int skip,
                       double (*similarity_func)(MRI_SURFACE *mris, MRI *mri_reg, 
                                                 MRI *mri_mask, MATRIX *m, 
                                                 int skip, double scale, int diag))

{
  double   rot_scale, max_similarity, similarity, rx, ry, rz ;
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
#define WHALF 8
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
            if ((++niter % write_iter) == 0)
              write_snapshot(mris, M_min, outregfile, niter) ;
            printf("%03d: new max %2.2f found at R=(%2.3f, %2.3f, %2.3f)\n",
                   niter, max_similarity, DEGREES(angles[0]), 
                   DEGREES(angles[1]), DEGREES(angles[2])) ;
            //            MatrixPrint(Gstdout, M_min) ;
            similarity = (*similarity_func)(mris, mri_reg,  mri_mask, M_min, skip, skip/8.0, 1) ;
            good_step = 1 ;
          }
        }
    if (good_step == 0)
    {
      rot_scale /= 2 ;
      printf("reducing scale to %2.4f deg\n", rot_scale) ;
    }
    good_step = 0 ;
  } while (rot_scale >= .001);

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
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    V3_X(v1) = v->whitex ; V3_Y(v1) = v->whitey ; V3_Z(v1) = v->whitez ;
    v2 = MatrixMultiply(m, v1, v2) ;
    v->x = V3_X(v2) ; v->y = V3_Y(v2) ; v->z = V3_Z(v2) ;
  }

  sprintf(fname, "%s%03d", fname_only,n) ;
  printf("writing snapshot to %s\n", fname) ;
  MRISwrite(mris, fname) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    V3_X(v1) = v->pialx ; V3_Y(v1) = v->pialy ; V3_Z(v1) = v->pialz ;
    v2 = MatrixMultiply(m, v1, v2) ;
    v->x = V3_X(v2) ; v->y = V3_Y(v2) ; v->z = V3_Z(v2) ;
  }

  sprintf(fname, "%s_pial%03d", fname_only,n) ;
  printf("writing snapshot to %s\n", fname) ;
  MRISwrite(mris, fname) ;

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

  old_write = write_iter ; write_iter = 1 ;

  // extract rigid body parameters from matrix
  MatrixToRigidParameters(mat, &xr, &yr, &zr, &xt, &yt, &zt) ;
  printf("initial rigid body parameters = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
         DEGREES(xr), DEGREES(yr), DEGREES(zr), xt, yt, zt) ;
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
  printf("%3.3d: best alignment at after powell: "
         "%2.5f (%d steps)\n\tparms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
         niter,fret, iter,
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
    printf("%3.3d: best alignment at after powell: "
           "%2.5f (%d steps)\n\tparms = (%2.4f, %2.4f, %2.4f) + (%2.2f, %2.2f, %2.2f)\n",
           niter,fret, iter,
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
                                   Gskip/8.0, 0) ;
  return(-(float)similarity) ;
}
