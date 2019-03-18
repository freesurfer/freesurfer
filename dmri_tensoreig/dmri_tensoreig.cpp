/**
 * @file  dti_tensoreig.c
 * @brief calculates eigensystem and fa
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author$
 *    $Date$
 *    $Revision$
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

#include <float.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/utsname.h>
#include <unistd.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <errno.h>

#include "macros.h"
#include "const.h"
#include "machine.h"
#include "fio.h"
#include "utils.h"
#include "mri.h"
#include "mri2.h"
#include "gcamorph.h"
#include "minc.h"
#include "analyze.h"
#include "mri_identify.h"
#include "error.h"
#include "diag.h"
#include "version.h"
#include "mghendian.h"
#include "fio.h"
#include "cmdargs.h"
#include "fsenv.h"
#include "DICOMRead.h"
#include "dti.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
static int fastestdim45(MRI *invol, int n1, int n2);
static void MRIavg4(MRI *invol, int n4, int n5);
static void MRIavg5(MRI *invol, int n4, int n5);
static float avgsnr(MRI *invol, MRI *mask, int frame0, int nframe);

static char vcid[] = "$Id$";

const char *Progname ;
FILE *fpin;
char *InFile = NULL;
char *OutDir = NULL;
float bValue = 0.0;
int nAcq = 1;
int nDir = 6;
char *GradFile = NULL;
char *MaskFile = NULL;
char *OutFmt = "nii";
int IsTensorInput = 0;
int debug = 0;
struct utsname uts;
char *cmdline, cwd[2000];
int DiffMode;
char tmpstr[2000];
FSENV *fsenv;

/***-------------------------------------------------------****/
int main(int argc, char *argv[]) {
  int nargs, ng, ig, i, nx, ny, nz, nf, navg, ix, iy, iz, id, nRow, nCol;
  float smax, ssum, norm, mean;
  float *grads = NULL;
  float *gp = NULL;
  float *eval = NULL;
  MATRIX *B = NULL, *Bpseudo = NULL, *dwi = NULL, *tensor = NULL, 
    *Tensor = NULL, *Evec = NULL;
  MRI *invol = NULL, *mask = NULL, *lowb = NULL, *avgdwi = NULL,
  *tenstack = NULL, *eigval = NULL, 
    *eigvec1 = NULL, *eigvec2 = NULL, *eigvec3 = NULL, 
    *trace = NULL, *fa = NULL;
  FILE *fp = NULL;
  char outfile[1024];
  const double minexp = exp(-10^35);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name$");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();

  fsenv = FSENVgetenv();

  /* parse the command line */
  parse_commandline(argc, argv);
  check_options();

  /* Read input volume */
  invol = MRIread(InFile);
  if (invol == NULL) {
    printf("ERROR: Could not open %s\n", InFile);
    exit(1);
  }
  
  nx = invol->width;
  ny = invol->height;
  nz = invol->depth;
  nf = invol->nframes;
  printf("INFO: Image size %d x %d x %d x %d\n", nx, ny, nz, nf);
  
  if( IsTensorInput == 1 ) {
    
    printf("INFO: Calculating eigensystem and fa from input tensors only.");
    
    // allocate memory for the eigensystem and the fa
    eigval = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigval);

    eigvec1 = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigvec1);
  
    eigvec2 = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigvec2);
  
    eigvec3 = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigvec3);
    
    fa = MRIalloc(nx, ny, nz, MRI_FLOAT);
    MRIcopyHeader(invol, fa);
    
    // our symmetric tensor
    Tensor = MatrixAlloc(3, 3, MATRIX_REAL);
    eval = (float *) malloc(3*sizeof(float));
    Evec = MatrixAlloc(3, 3, MATRIX_REAL);

    // go through all the tensors and calculate the eigensystem and fa
    for (iz=0; iz<nz; iz++) {
      for (iy=0; iy<ny; iy++) {
        for (ix=0; ix<nx; ix++) {

          // copy the tensor over
          id = 0;
          for( nRow=1; nRow<=3; nRow++ ) {
            for( nCol=1; nCol<=3; nCol++ ) {
              Tensor->rptr[nRow][nCol] = MRIgetVoxVal(invol, ix, iy, iz, id);
              id++;
            }
          }
                                
          // calculate eigensystem
          MatrixEigenSystem(Tensor, eval, Evec);

          // copy eigensystem
          for( nRow=1; nRow<=3; nRow++ ) {
            MRIsetVoxVal(eigval, ix, iy, iz, 0, eval[0]);
            MRIsetVoxVal(eigvec1, ix, iy, iz, nRow-1, Evec->rptr[nRow][1]);
            MRIsetVoxVal(eigvec2, ix, iy, iz, nRow-1, Evec->rptr[nRow][2]);
            MRIsetVoxVal(eigvec3, ix, iy, iz, nRow-1, Evec->rptr[nRow][3]);
          }
          
          // calculate fa
          mean = ( eval[0] + eval[1] + eval[2] ) / 3;
          norm = eval[0]*eval[0] + eval[1]*eval[1] + eval[2]*eval[2];
          ssum = 0;
          for (i=0; i<3; i++) {
            float tmp = eval[i] - mean;
            ssum += tmp*tmp;
          }
          MRIsetVoxVal(fa, ix, iy, iz, 0, sqrt(ssum/norm));
          
        }
      }
    }
            
  } else {  
    printf("INFO: Reading gradient vectors.");
    /* Read gradient vectors */
    ng = 3*(nDir+1);
    grads = (float *) malloc(ng*sizeof(float));
    memset(grads, 0, ng*sizeof(float));
    fp = fopen(GradFile, "r");
    if (fp == NULL) {
      printf("ERROR: Could not open %s\n", GradFile);
      exit(1);
    }
    ig = 0;
    gp = grads+3;
    while (fscanf(fp, "%f", gp++) != EOF)
      ig++;
    fclose(fp);
    if (ig != ng-3) {
      printf("ERROR: Expected %d values in %s, found %d\n", 
             ng-3, GradFile, ig);
      exit(1);
    }
    if (debug) {
      gp = grads;
      for (ig=nDir+1; ig>0; ig--) {
        printf("INFO: Gradient vector %g %g %g\n", *gp, *(gp+1), *(gp+2));
        gp += 3;
      }
    }
  
    /* Normalize gradients to unit maximum length */
    smax = 0;
    gp = grads+3;
    for (ig=nDir; ig>0; ig--) {
      ssum = (*gp) * (*gp);
      gp++;
      ssum += (*gp) * (*gp);
      gp++;
      ssum += (*gp) * (*gp);
      gp++;
      if (ssum>smax) smax=ssum;
    }
    smax = sqrt(smax);
    gp = grads+3;
    for (ig=ng-3; ig>0; ig--)
      *gp++ /= smax;
  
    /* Calculate B matrix */
    B = MatrixAlloc(nDir+1, 7, MATRIX_REAL);
    memset(B->rptr[1], 0, 6*sizeof(float));
    B->rptr[1][7] = 1;
    gp = grads+3;
    for (ig=2; ig<=nDir+1; ig++) {
      ssum = bValue * (*gp);
      B->rptr[ig][1] = ssum * (*gp);  /* b g_1^2 */
      B->rptr[ig][2] = 2 * ssum * (*(gp+1)); /* 2 b g_1 g_2 */
      B->rptr[ig][3] = 2 * ssum * (*(gp+2)); /* 2 b g_1 g_3 */
      gp++;
      ssum = bValue * (*gp);
      B->rptr[ig][4] = ssum * (*gp);  /* b g_2^2 */
      B->rptr[ig][5] = 2 * ssum * (*(gp+1)); /* 2 b g_2 g_3 */
      gp++;
      B->rptr[ig][6] = bValue * (*gp) * (*gp); /* b g_3^2 */
      gp++;
      B->rptr[ig][7] = 1;
    }
    if (debug) {
      printf("B =\n");
      for (ig=1; ig<=nDir+1; ig++) {
        for (i=1; i<=7; i++)
          printf("%f ", B->rptr[ig][i]);
        printf("\n");
      }
    }
  
    /* Calculate pseudoinverse of B */
    Bpseudo = MatrixAlloc(7, nDir+1, MATRIX_REAL);
    MatrixPseudoInverse(B, Bpseudo);
    if (debug) {
      printf("Bpseudo =\n");
      for (i=1; i<=7; i++) {
        for (ig=1; ig<=nDir+1; ig++)
          printf("%f ", Bpseudo->rptr[i][ig]);
        printf("\n");
      }
    }
    
    /* Average volumes to get one volume per diffusion direction */
    navg = nf / (nDir+1);
    if (navg > 1) {
      printf("INFO: Averaging %d volumes per diffusion direction\n", navg);
      if ( fastestdim45(invol, nDir+1, navg) == 1 ) {
        printf("INFO: Detected highest image dimensions %d x %d (dir x acq)\n",
               nDir+1, navg);
        MRIavg5(invol, nDir+1, navg);
      } else {
        printf("INFO: Detected highest image dimensions %d x %d (acq x dir)\n",
               navg, nDir+1);
        MRIavg4(invol, nDir+1, navg);
      }
    }
  
    /* Read binary mask */
    if (MaskFile) {
      mask = MRIread(MaskFile);
      if (mask == NULL) {
        printf("ERROR: Could not open %s\n", MaskFile);
        exit(1);
      }
      if ( (nx != mask->width) || 
           (ny != mask->height) || 
           (nz != mask->depth) ) {
        printf("ERROR: Mask size %d x %d x %d does not agree "
               "with image size\n",
               mask->width, mask->height, mask->depth);
        exit(1);
      }
    } else {
      mask = MRIalloc(nx, ny, nz, MRI_INT);
      for (iz=0; iz<nz; iz++)
        for (iy=0; iy<ny; iy++)
          for (ix=0; ix<nx; ix++) {
            int mflag = 1;
            for (id=0; id<nDir+1; id++)
              if ( MRIgetVoxVal(invol, ix, iy, iz, id) < minexp )
                mflag = 0;
            MRIsetVoxVal(mask, ix, iy, iz, 0, mflag);
          }
    }
    printf("INFO: Mask contains %g%% of total voxels\n",
           MRIvoxelMean(mask, nx/2, ny/2, nz/2, MAX(MAX(nx,ny),nz), 0) * 100);
    printf("INFO: SNR of low-b image is %g\n",
           avgsnr(invol, mask, 0, 1));
    printf("INFO: Average SNR of all DW images is %g\n",
           avgsnr(invol, mask, 0, nDir+1));
  
    /* Estimate tensors and diffusion measures */
    lowb = MRIalloc(nx, ny, nz, MRI_FLOAT);
    MRIcopyHeader(invol, lowb);
  
    avgdwi = MRIalloc(nx, ny, nz, MRI_FLOAT);
    MRIcopyHeader(invol, avgdwi);
  
    tenstack = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 9);
    MRIcopyHeader(invol, tenstack);
  
    eigval = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigval);
  
    eigvec1 = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigvec1);
  
    eigvec2 = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigvec2);
  
    eigvec3 = MRIallocSequence(nx, ny, nz, MRI_FLOAT, 3);
    MRIcopyHeader(invol, eigvec3);
  
    trace = MRIalloc(nx, ny, nz, MRI_FLOAT);
    MRIcopyHeader(invol, trace);
  
    fa = MRIalloc(nx, ny, nz, MRI_FLOAT);
    MRIcopyHeader(invol, fa);
  
    dwi = MatrixAlloc(nDir+1, 1, MATRIX_REAL);
    
    // these are horribly named
    tensor = MatrixAlloc(7, 1, MATRIX_REAL);
    Tensor = MatrixAlloc(3, 3, MATRIX_REAL);
    eval = (float *) malloc(3*sizeof(float));
    Evec = MatrixAlloc(3, 3, MATRIX_REAL);
    
    // go through each voxel and calculate tensors, trace, fa, etc.
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
        for (ix=0; ix<nx; ix++)
          if ( MRIgetVoxVal(mask, ix, iy, iz, 0) ) {
  
            /* Average and mask DWI voxel values
               (just to save to disk, not used in tensor estimation) */
            MRIsetVoxVal(lowb, ix, iy, iz, 0, 
                         MRIgetVoxVal(invol, ix, iy, iz, 0));
            mean = 0;
            for (id=1; id<nDir+1; id++)
              mean += MRIgetVoxVal(invol, ix, iy, iz, id);
            MRIsetVoxVal(avgdwi, ix, iy, iz, 0, mean/nDir);
  
            /* Estimate tensor */
            for (id=0; id<nDir+1; id++)
              dwi->rptr[id+1][1] = -log( MRIgetVoxVal(invol, ix, iy, iz, id) );
  
            MatrixMultiply(Bpseudo, dwi, tensor);
  
            // create a symmetric tensor
            Tensor->rptr[1][1] = tensor->rptr[1][1];
            Tensor->rptr[1][2] = Tensor->rptr[2][1] = tensor->rptr[2][1];
            Tensor->rptr[1][3] = Tensor->rptr[3][1] = tensor->rptr[3][1];
            Tensor->rptr[2][2] = tensor->rptr[4][1];
            Tensor->rptr[2][3] = Tensor->rptr[3][2] = tensor->rptr[5][1];
            Tensor->rptr[3][3] = tensor->rptr[6][1];
  
            MRIsetVoxVal(tenstack, ix, iy, iz, 0, Tensor->rptr[1][1]);
            MRIsetVoxVal(tenstack, ix, iy, iz, 1, Tensor->rptr[2][1]);
            MRIsetVoxVal(tenstack, ix, iy, iz, 2, Tensor->rptr[3][1]);
  
            MRIsetVoxVal(tenstack, ix, iy, iz, 3, Tensor->rptr[1][2]);
            MRIsetVoxVal(tenstack, ix, iy, iz, 4, Tensor->rptr[2][2]);
            MRIsetVoxVal(tenstack, ix, iy, iz, 5, Tensor->rptr[3][2]);
  
            MRIsetVoxVal(tenstack, ix, iy, iz, 6, Tensor->rptr[1][3]);
            MRIsetVoxVal(tenstack, ix, iy, iz, 7, Tensor->rptr[2][3]);
            MRIsetVoxVal(tenstack, ix, iy, iz, 8, Tensor->rptr[3][3]);
  
            /* Do eigen-decomposition */
            MatrixEigenSystem(Tensor, eval, Evec);
  
            // copy eigensystem
            for( nRow=1; nRow<=3; nRow++ ) {
              MRIsetVoxVal(eigval, ix, iy, iz, 0, eval[0]);
              MRIsetVoxVal(eigvec1, ix, iy, iz, nRow-1, Evec->rptr[nRow][1]);
              MRIsetVoxVal(eigvec2, ix, iy, iz, nRow-1, Evec->rptr[nRow][2]);
              MRIsetVoxVal(eigvec3, ix, iy, iz, nRow-1, Evec->rptr[nRow][3]);
            }
  
            /* Calculate diffusion measures */
            mean = eval[0] + eval[1] + eval[2];
            MRIsetVoxVal(trace, ix, iy, iz, 0, mean);
            norm = eval[0]*eval[0] + eval[1]*eval[1] + eval[2]*eval[2];
            mean /= 3;
            ssum = 0;
            for (i=0; i<3; i++) {
              float tmp = eval[i] - mean;
              ssum += tmp*tmp;
            }
            MRIsetVoxVal(fa, ix, iy, iz, 0, sqrt(ssum/norm));
          } else {
            MRIsetVoxVal(lowb, ix, iy, iz, 0, 0);
            MRIsetVoxVal(avgdwi, ix, iy, iz, 0, 0);
            for (i=0; i<9; i++) {
              MRIsetVoxVal(tenstack, ix, iy, iz, i, 0);
            }
            for (i=0; i<3; i++) {
              MRIsetVoxVal(eigval, ix, iy, iz, i, 0);
              MRIsetVoxVal(eigvec1, ix, iy, iz, i, 0);
              MRIsetVoxVal(eigvec2, ix, iy, iz, i, 0);
              MRIsetVoxVal(eigvec3, ix, iy, iz, i, 0);
            }
            MRIsetVoxVal(trace, ix, iy, iz, 0, 0);
            MRIsetVoxVal(fa, ix, iy, iz, 0, 0);
          }
  
    /* Write output files */
    sprintf(outfile, "%s/lowb.%s", OutDir, OutFmt);
    MRIwrite(lowb, outfile);
  
    sprintf(outfile, "%s/avgdwi.%s", OutDir, OutFmt);
    MRIwrite(avgdwi, outfile);
  
    sprintf(outfile, "%s/dtensor.%s", OutDir, OutFmt);
    MRIwrite(tenstack, outfile);
  
    sprintf(outfile, "%s/trace.%s", OutDir, OutFmt);
    MRIwrite(trace, outfile);
    
    /* Write mask, if it was created by this program */
    if (!MaskFile) {
      sprintf(outfile, "%s/mask.%s", OutDir, OutFmt);
      MRIwrite(mask, outfile);
    }
  } // end else IsTensorInput
  
  // these always get saved out
  sprintf(outfile, "%s/fa.%s", OutDir, OutFmt);
  MRIwrite(fa, outfile);

  sprintf(outfile, "%s/eigval.%s", OutDir, OutFmt);
  MRIwrite(eigval, outfile);

  sprintf(outfile, "%s/eigvec1.%s", OutDir, OutFmt);
  MRIwrite(eigvec1, outfile);

  sprintf(outfile, "%s/eigvec2.%s", OutDir, OutFmt);
  MRIwrite(eigvec2, outfile);

  sprintf(outfile, "%s/eigvec3.%s", OutDir, OutFmt);
  MRIwrite(eigvec3, outfile);  

  /* Write log file */
  sprintf(outfile, "%s/dti_tensoreig.log", OutDir);
  fp = fopen(outfile, "w");
  dump_options(fp);
  fclose(fp);

  /* Free allocated memory */
  if( grads != NULL ) {
    free(grads);
  }
  
  if( eval != NULL ) {
    free(eval);
  }
  
  MRIfree(&invol);
  MRIfree(&mask);
  MRIfree(&lowb);
  MRIfree(&avgdwi);
  MRIfree(&tenstack);
  MRIfree(&eigval);
  MRIfree(&eigvec1);
  MRIfree(&eigvec2);
  MRIfree(&eigvec3);
  MRIfree(&trace);
  MRIfree(&fa);
  MatrixFree(&B);
  MatrixFree(&Bpseudo);
  MatrixFree(&dwi);
  MatrixFree(&tensor);
  MatrixFree(&Tensor);
  MatrixFree(&Evec);
  exit(0);

} /* end main() */

/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused, err;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))         print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;

    else if (!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      InFile = pargv[0];
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--o")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutDir = pargv[0];
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--sdcm") || 
               !strcasecmp(option, "--infodump")) {
      // Pass it a siemens dicom or an info dump file from which parameters
      // should be extracted.
      if (nargc < 1) CMDargNErr(option,1);
      err = DTIparamsFromSiemensAscii(pargv[0], &bValue, &nDir, &DiffMode);
      if (err) exit(1);
      switch (DiffMode) {
      case 3:
        sprintf(tmpstr,"%s/diffusion/graddir/6-cube-mghnew.txt",
                fsenv->FREESURFER_HOME);
        GradFile = strcpyalloc(tmpstr);
        break;
      default:
        printf("WARNING: diffusion mode %d unrecoginzed\n",DiffMode);
      }
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--b")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%f", &bValue);
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--ndir")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &nDir);
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--nacq")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &nAcq);
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--g")) {
      if (nargc < 1) CMDargNErr(option,1);
      GradFile = pargv[0];
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--m")) {
      if (nargc < 1) CMDargNErr(option,1);
      MaskFile = pargv[0];
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--ofmt")) {
      if (nargc < 1) CMDargNErr(option,1);
      OutFmt = pargv[0];
      nargc --;
      pargv ++;
    } else if (!strcasecmp(option, "--tensor")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0], "%d", &IsTensorInput);
      nargc --;
      pargv ++;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s <options> \n",Progname) ;
  printf("\n");
  printf("   --i   file : input file \n");
  printf("   --o    dir : output directory \n");
  printf("   --ofmt fmt : format (extension) of output files\n");
  printf("   --b    num : b value \n");
  printf("   --nacq num : number of T2 weightings \n");
  printf("   --ndir num : number of diffusion directions \n");
  printf("   --g   file : gradient file \n");
  printf("   --m   file : mask file\n");
  printf("   --tensor   : create the output from the input tensors,"
         "rather than creating\nthem from the diffusion weighted images \n");
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf(
    "\n"
    "Proceed with caution.\n"
  );


  exit(1) ;
}

/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}

/* --------------------------------------------- */
static void check_options(void) {
  if (InFile==NULL) {
    printf("ERROR: Need to specify input file\n");
    exit(1);
  }
  if (OutDir==NULL) {
    printf("ERROR: Need to specify output directory\n");
    exit(1);
  } else {
    int err = mkdir(OutDir, 0777);
    if (err != 0 && errno != EEXIST) {
      printf("ERROR: Could not create output directory %s\n", OutDir);
      perror(NULL);
      exit(1);
    }
  }
  if (bValue < 0 && IsTensorInput != 1) {
    printf("ERROR: Negative b value %g\n", bValue);
    exit(1);
  }
  return;
}

/* ------------------------------------------------------ */
static void usage_exit(void) {
  print_usage() ;
  exit(1) ;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp, "\n");
  fprintf(fp, "%s\n", vcid);
  fprintf(fp, "Cwd %s\n", cwd);
  fprintf(fp, "Cmdline %s\n", cmdline);
  fprintf(fp, "Sysname %s\n", uts.sysname);
  fprintf(fp, "Hostname %s\n", uts.nodename);
  fprintf(fp, "Machine %s\n", uts.machine);
  fprintf(fp, "User %s\n", VERuser());
  fprintf(fp, "Input file %s\n", InFile);
  fprintf(fp, "Output directory %s\n", OutDir);
  fprintf(fp, "Output format %s\n", OutFmt);
  fprintf(fp, "Mask %s\n", MaskFile);
  fprintf(fp, "b value %f\n", bValue);
  fprintf(fp, "Number of T2 weightings %d\n", nAcq);
  fprintf(fp, "Number of diffusion directions %d\n", nDir);
  fprintf(fp, "Gradients %s\n", GradFile);
  fprintf(fp, "Tensor Input %d\n", IsTensorInput);

  return;
}

/***-------------------------------------------------------****/
static int fastestdim45(MRI *invol, int n1, int n2) {
  const int nx = invol->width;
  const int ny = invol->height;
  const int nz = invol->depth;
  const int nstep = 10;
  const int stepx = nx/nstep;
  const int stepy = ny/nstep;
  const int stepz = nz/nstep;
  int ix, iy, iz, i1, i2;
  float mean1, var1, avgvar1, mean2, var2, avgvar2, tmp;

  avgvar1 = avgvar2 = 0.0;
  for (iz=0; iz<nz; iz+=stepz)
    for (iy=0; iy<ny; iy+=stepy)
      for (ix=0; ix<nx; ix+=stepx)
        for (i2=0; i2<n2; i2++) {

          mean1 = mean2 = 0.0;
          for (i1=0; i1<n1; i1++) {
            mean1 += MRIgetVoxVal(invol, ix, iy, iz, i1+i2*n1);
            mean2 += MRIgetVoxVal(invol, ix, iy, iz, i2+i1*n2);
          }
          mean1 /= n1;
          mean2 /= n1;

          var1 = var2 = 0.0;
          for (i1=0; i1<n1; i1++) {
            tmp = MRIgetVoxVal(invol, ix, iy, iz, i1+i2*n1) - mean1;
            var1 += tmp * tmp;
            tmp = MRIgetVoxVal(invol, ix, iy, iz, i2+i1*n2) - mean2;
            var2 += tmp * tmp;
          }
          var1 /= (n1-1);
          var2 /= (n1-1);

          avgvar1 += var1;
          avgvar2 += var2;
        }

  if (avgvar2 > avgvar1)  return 1;
  else     return 2;
}

/***-------------------------------------------------------****/
static void MRIavg4(MRI *invol, int n4, int n5) {
  const int nx = invol->width;
  const int ny = invol->height;
  const int nz = invol->depth;
  int ix, iy, iz, i5, iframe, endframe;
  float mean;

  for (iz=0; iz<nz; iz++)
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++)
        for (i5=0; i5<n5; i5++) {
          mean = 0.0;
          endframe = (i5+1)*n4;
          for (iframe=i5*n4; iframe<endframe; iframe++)
            mean += MRIgetVoxVal(invol, ix, iy, iz, iframe) / n4;
          MRIsetVoxVal(invol, ix, iy, iz, i5, mean);
        }

  return;
}

/***-------------------------------------------------------****/
static void MRIavg5(MRI *invol, int n4, int n5) {
  const int nx = invol->width;
  const int ny = invol->height;
  const int nz = invol->depth;
  const int nframe = invol->nframes;
  int ix, iy, iz, i4, iframe;
  float mean;

  for (iz=0; iz<nz; iz++)
    for (iy=0; iy<ny; iy++)
      for (ix=0; ix<nx; ix++)
        for (i4=0; i4<n4; i4++) {
          mean = 0.0;
          iframe = i4;
          for (iframe=i4; iframe<nframe; iframe+=n4)
            mean += MRIgetVoxVal(invol, ix, iy, iz, iframe) / n5;
          MRIsetVoxVal(invol, ix, iy, iz, i4, mean);
        }

  return;
}

/***-------------------------------------------------------****/
static float avgsnr(MRI *invol, MRI *mask, int frame0, int nframe) {
  const int nx = invol->width;
  const int ny = invol->height;
  const int nz = invol->depth;
  const int frame1 = frame0 + nframe;
  int ix, iy, iz, iframe, nmask;
  float meanfg, meanbg, varbg, snr;

  snr = 0.0;
  for (iframe=frame0; iframe<frame1; iframe++) {

    nmask = 0;
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
        for (ix=0; ix<nx; ix++)
          nmask += MRIgetVoxVal(mask, ix, iy, iz, 0);

    meanfg = 0.0;
    for (iz=0; iz<nz; iz++)
      for (iy=0; iy<ny; iy++)
        for (ix=0; ix<nx; ix++)
          if ( MRIgetVoxVal(mask, ix, iy, iz, 0) )
            meanfg += MRIgetVoxVal(invol, ix, iy, iz, iframe) / nmask;

    meanbg = 0.0;
    for (iz=0; iz<5; iz++)
      for (iy=0; iy<5; iy++)
        for (ix=0; ix<5; ix++)
          meanbg += MRIgetVoxVal(invol, ix, iy, iz, iframe) / 125;

    varbg = 0.0;
    for (iz=0; iz<5; iz++)
      for (iy=0; iy<5; iy++)
        for (ix=0; ix<5; ix++) {
          float tmp = MRIgetVoxVal(invol, ix, iy, iz, iframe) - meanbg;
          varbg += tmp*tmp / 124;
        }

    snr += meanfg / sqrt(varbg) / nframe;
  }

  return(snr);
}
