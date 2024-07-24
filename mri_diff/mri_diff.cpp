/**
 * @brief Determines whether two volumes differ.
 *
 * The basic usage is something like:
 *
 *   mri_diff vol1 vol2
 *
 * It then prints to the terminal whether they differ or not.
 *
 * Volumes can differ in six ways:
 * 1. Dimension,               return status = 101
 * 2. Resolutions,             return status = 102
 * 3. Acquisition Parameters,  return status = 103
 * 4. Geometry,                return status = 104
 * 5. Precision,               return status = 105
 * 6. Pixel Data,              return status = 106
 * 7. m3z 3D Morph,            return status = 201
 * 8. Color Table,             return status = 202
 */
/*
 * Original Author: Doug Greve
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


// mri_voldiff v1 v2
//
/*
  BEGINHELP

  Determines whether two volumes differ. See below for what
  'differ' means.

  The basic usage is something like:

  mri_diff vol1 vol2

  It then prints to the terminal whether they differ or not.

  NOTE: stuff might get printed to the terminal regardless of whether
  the volumes are different or not, so you cannot just use the
  presence or absence of terminal output to determine whether they
  are different as you can with unix diff.

  There are three ways to determine whether they are different:
  1. Look for 'volumes differ' in the terminal output
  2. Program exits with status > 100
  3. Create a log file

  To create a log file, add --log yourlogfile any where in the
  command line, eg:

  mri_diff vol1 vol2 --log yourlogfile

  If yourlogfile exists, it will be immediately deleted when the
  program starts. If a difference is detected, yourlogfile will
  be created, and information about the difference will be
  printed to it. If there is no difference, yourlogfile will
  not exist when the program finishes.

  Volumes can differ in six ways:
  1. Dimension,               return status = 101
  2. Resolutions,             return status = 102
  3. Acquisition Parameters,  return status = 103
  4. Geometry,                return status = 104
  5. Precision,               return status = 105
  6. Pixel Data,              return status = 106
  7. m3z 3D Morph,            return status = 201
  8. Color Table,             return status = 202

  Dimension is number of rows, cols, slices, and frames.
  Resolution is voxel size.
  Acqusition parameters are: flip angle, TR, TE, and TI.
  Geometry checks the vox2ras matrices for equality.
  Precision is int, float, short, etc.

  By default, all of these are checked, but some can be turned off
  with certain command-line flags:

  --notallow-res  : turns off resolution checking
  --notallow-acq  : turns off acquistion parameter checking
  --notallow-geo  : turns off geometry checking
  --notallow-prec : turns off precision checking
  --notallow-pix  : turns off pixel checking

  In addition, the minimum difference in pixel value required
  to be considered different can be controlled with --thresh.
  Eg, if two volumes differ by .00001, and you do not consider
  this to be significant, then --thresh .00002 will prevent
  that difference from being considered different. The default
  threshold is 0.

  --diff diffvol

  Saves difference image to diffvol.

  QUALITY ASSURANCE

  mri_diff can be used to check that two volumes where acquired in the
  same way with the --qa flag. This turns on Res, Acq, and Prec, and
  turns off Geo and Pix. Instead of checking geometry, it checks for the
  basic orientation of the volumes (eg, RAS, LPI, etc). The idea here is
  that the pixel data and exact geometry may be different, but other
  things should be the same.

  EXIT STATUS

  0   Volumes are not different and there were no errors
  1   Errors encountered. Volumes may or may not be different
  101 Volumes differ in dimension
  102 Volumes differ in resolution
  103 Volumes differ in acquisition parameters
  104 Volumes differ in geometry
  105 Volumes differ in precision
  106 Volumes differ in pixel data
  107 Volumes differ in orientation
  201 m3z 3D Morph differ
  202 Volumes differ in Color Table

  ENDHELP
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "cma.h"

#include "mri_identify.h"
#include "gcamorph.h"

static int diff_colortable(const COLOR_TABLE *ctab1, const COLOR_TABLE *ctab2);
static void diff_mgh_morph(const char *file1, const char *file2);
static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int verbose=0;
int checkoptsonly=0;
struct utsname uts;

char *InVol1File=NULL;
char *InVol2File=NULL;
char *subject, *hemi, *SUBJECTS_DIR;
double pixthresh=0, resthresh=0, geothresh=0, rgbthresh = 0;
int diffcountthresh = 0;
char *DiffFile=NULL;
int DiffAbs=0, AbsDiff=1,DiffPct=0;
char *AvgDiffFile=NULL;

MRI *InVol1=NULL, *InVol2=NULL, *DiffVol=NULL, *DiffLabelVol=NULL;
char *DiffVolFile=NULL;
char *DiffLabelVolFile=NULL;

int CheckColortable = 0;
int CheckImag = 0;
int CheckResolution=1;
int CheckAcqParams=1;
int CheckPixVals=1;
int CheckGeo=1;
int CheckOrientation=1;
int CheckPrecision=1;
int CheckType=1;
int SegDiff = -2;
char *SegDiffFile=NULL;
MATRIX *vox2ras1,*vox2ras2;
char Orient1[4], Orient2[4];

int ExitOnDiff = 1;
int ExitStatus = 0;
int DoRSS = 0; // Compute sqrt of sum squares
int PrintSSD = 0; // Print sum of squared differences over all voxel
int PrintRMS = 0; // Print root mean squared differences over all voxel

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs, r, c, s, f;
  double val1, val2;
  FILE *fp=NULL;

  nargs = handleVersionOption(argc, argv, "mri_diff");
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
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);

  if (debug) dump_options(stdout);

  // njs: commented-out so that mri_diff is useful w/o SUBJECTS_DIR declared
  //  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  //  if(SUBJECTS_DIR == NULL){
  //    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
  //    exit(1);
  //  }

  printf("%s\n", cmdline);
  if (DiffFile) {
    if (fio_FileExistsReadable(DiffFile)) unlink(DiffFile);
    if (fio_FileExistsReadable(DiffFile)) {
      printf("ERROR: could not delete %s\n",DiffFile);
      exit(1);
    }
  }

  int fileType1 = mri_identify(InVol1File);
  if (fileType1 == MGH_MORPH)
  {
    int fileType2 = mri_identify(InVol2File);
    if (fileType1 != fileType2)
    {
      printf("Input files %s and %s are different types\n", InVol1File, InVol2File); 
      exit(1);
    }
    
    diff_mgh_morph(InVol1File, InVol2File);
    exit(ExitStatus);
  }

  InVol1 = MRIread(InVol1File);
  if (InVol1==NULL) exit(1);

  InVol2 = MRIread(InVol2File);
  if (InVol2==NULL) exit(1);

  /*- Check dimension ---------------------------------*/
  if (InVol1->width   != InVol2->width  ||
      InVol1->height  != InVol2->height ||
      InVol1->depth   != InVol2->depth  ||
      InVol1->nframes != InVol2->nframes) {
    printf("Volumes differ in dimension\n");
    printf("v1dim %d %d %d %d\n",
           InVol1->width,InVol1->height,InVol1->depth,InVol1->nframes);
    printf("v2dim %d %d %d %d\n",
           InVol2->width,InVol2->height,InVol2->depth,InVol2->nframes);
    if (DiffFile) {
      fp = fopen(DiffFile,"w");
      if (fp==NULL) {
        printf("ERROR: could not open %s\n",DiffFile);
        exit(1);
      }
      dump_options(fp);
      fprintf(fp,"v1dim %d %d %d %d\n",
              InVol1->width,InVol1->height,InVol1->depth,InVol1->nframes);
      fprintf(fp,"v2dim %d %d %d %d\n",
              InVol2->width,InVol2->height,InVol2->depth,InVol2->nframes);
      fprintf(fp,"Volumes differ in dimension\n");
      fclose(fp);
    }
    // Must exit here
    exit(101);
  }

  //---------------------------------------------------
  if (CheckResolution) {
    if (fabs(InVol1->xsize - InVol2->xsize) > resthresh  ||
        fabs(InVol1->ysize - InVol2->ysize) > resthresh  ||
        fabs(InVol1->zsize - InVol2->zsize) > resthresh) {
      printf("Volumes differ in resolution\n");
      printf("v1res %f %f %f\n",
             InVol1->xsize,InVol1->ysize,InVol1->zsize);
      printf("v2res %f %f %f\n",
             InVol2->xsize,InVol2->ysize,InVol2->zsize);
      printf("diff %f %f %f\n",
             InVol1->xsize-InVol2->xsize,
             InVol1->ysize-InVol2->ysize,
             InVol1->zsize-InVol2->zsize);
      if (DiffFile) {
        fp = fopen(DiffFile,"w");
        if (fp==NULL) {
          printf("ERROR: could not open %s\n",DiffFile);
          exit(1);
        }
        dump_options(fp);
        fprintf(fp,"v1res %f %f %f\n",
                InVol1->xsize,InVol1->ysize,InVol1->zsize);
        fprintf(fp,"v2res %f %f %f\n",
                InVol2->xsize,InVol2->ysize,InVol2->zsize);
        fprintf(fp,"Volumes differ in resolution\n");
        fclose(fp);
      }
      ExitStatus = 102;
      if (ExitOnDiff) exit(ExitStatus);
    }
  }

  //---------------------------------------------------
  if (CheckAcqParams) {
    if (InVol1->flip_angle   != InVol2->flip_angle ||
        InVol1->tr   != InVol2->tr ||
        InVol1->te   != InVol2->te ||
        InVol1->ti   != InVol2->ti) {
      printf("Volumes differ in acquisition parameters\n");
      printf("v1acq fa=%f tr=%f te=%f ti=%f\n",
             InVol1->flip_angle,InVol1->tr,InVol1->te,InVol1->ti);
      printf("v2acq fa=%f tr=%f te=%f ti=%f\n",
             InVol2->flip_angle,InVol2->tr,InVol2->te,InVol2->ti);
      if (DiffFile) {
        fp = fopen(DiffFile,"w");
        if (fp==NULL) {
          printf("ERROR: could not open %s\n",DiffFile);
          exit(1);
        }
        dump_options(fp);
        fprintf(fp,"v1acq fa=%f tr=%f te=%f ti=%f\n",
                InVol1->flip_angle,InVol1->tr,InVol1->te,InVol1->ti);
        fprintf(fp,"v2acq fa=%f tr=%f te=%f ti=%f\n",
                InVol2->flip_angle,InVol2->tr,InVol2->te,InVol2->ti);
        fprintf(fp,"Volumes differ in acquisition parameters\n");
        fclose(fp);
      }
      ExitStatus = 103;
      if (ExitOnDiff) exit(103);
    }
  }

  //------------------------------------------------------
  if (CheckGeo) {
    vox2ras1 = MRIxfmCRS2XYZ(InVol1,0);
    vox2ras2 = MRIxfmCRS2XYZ(InVol2,0);
    for (r=1; r<=4; r++) {
      for (c=1; c<=4; c++) {
        val1 = vox2ras1->rptr[r][c];
        val2 = vox2ras2->rptr[r][c];
        double diff = fabs(val1-val2);
        if (diff > geothresh) {
          printf("Volumes differ in geometry row=%d col=%d diff=%lf (%g)\n",
                 r,c,diff,diff);
          if (DiffFile) {
            fp = fopen(DiffFile,"w");
            if (fp==NULL) {
              printf("ERROR: could not open %s\n",DiffFile);
              exit(1);
            }
            dump_options(fp);
            fprintf(fp,"v1 vox2ras ----------\n");
            MatrixPrint(fp,vox2ras1);
            fprintf(fp,"v2 vox2ras ---------\n");
            MatrixPrint(fp,vox2ras2);

            fprintf(fp,"Volumes differ in geometry vox2ras r=%d c=%d %lf\n",
                    r,c,diff);
            fclose(fp);
          }
          ExitStatus = 104;
          if (ExitOnDiff) exit(104);
        }
      } // c
    } // r
  } // checkgeo

  //-------------------------------------------------------
  if (CheckPrecision) {
    if (InVol1->type != InVol2->type) {
      printf("Volumes differ in precision %d %d\n",
             InVol1->type,InVol2->type);
      if (DiffFile) {
        fp = fopen(DiffFile,"w");
        if (fp==NULL) {
          printf("ERROR: could not open %s\n",DiffFile);
          exit(1);
        }
        dump_options(fp);
        fprintf(fp,"Volumes differ in precision %d %d\n",
                InVol1->type,InVol2->type);

        fclose(fp);
      }
      ExitStatus = 105;
      if (ExitOnDiff) exit(105);
    }
  }

  //------------------------------------------------------
  // Compare pixel values
  if (CheckPixVals) {
    if(DiffVolFile) {
      if(!DoRSS) f = InVol1->nframes;
      else       f = 1;
      DiffVol = MRIallocSequence(InVol1->width,InVol1->height,
                                 InVol1->depth,MRI_FLOAT,f);
      MRIcopyHeader(InVol1,DiffVol);
    }
    if (DiffLabelVolFile) {
      DiffLabelVol = MRIallocSequence(InVol1->width,InVol1->height,
                                      InVol1->depth,MRI_FLOAT,InVol1->nframes);
      MRIcopyHeader(InVol1,DiffLabelVol);
    }

    int mode = MRI_COMPLEX_REAL;
    if (InVol1->type == MRI_FLOAT_COMPLEX && CheckImag) {
      mode = MRI_COMPLEX_IMAG;
      printf("Compare imaginary units ...\n");
    }
    
    double SumDiff=0.0;
    c=r=s=f=0;
    int notzero=0;
    val1 = MRIgetVoxVal(InVol1,c,r,s,f, mode);
    val2 = MRIgetVoxVal(InVol2,c,r,s,f, mode);
    double maxdiff = fabs(val1-val2);
    int cmax=0, rmax=0, smax=0, fmax=0;
    double SumSqDiff=0.0; // over all voxel
    int ndiff=0;
    for (c=0; c < InVol1->width; c++) {
      for (r=0; r < InVol1->height; r++) {
        for (s=0; s < InVol1->depth; s++) {
          double SumSqErr = 0.0; //over all frames
          for (f=0; f < InVol1->nframes; f++) {
            val1 = MRIgetVoxVal(InVol1,c,r,s,f, mode);
            val2 = MRIgetVoxVal(InVol2,c,r,s,f, mode);
	    if (val1 != 0 && val2 != 0) notzero++;
            double diff = val1-val2;
	    if(fabs(diff) > pixthresh) ndiff++;
            if(fabs(diff) > pixthresh && verbose) {
              printf("diff %6d %12.8f at %d %d %d %d  %g %g\n",ndiff,diff,c,r,s,f,val1,val2);
            }
            SumSqDiff += (diff*diff);
            SumSqErr  += (diff*diff);
            if(AbsDiff)   diff = fabs(diff);
            if(DiffAbs)   diff = fabs(fabs(val1)-fabs(val2));
            if(DiffPct)   diff = 100*(val1-val2)/((val1+val2)/2.0);
            if(DiffVolFile && !DoRSS)  MRIsetVoxVal(DiffVol,c,r,s,f,diff, mode);
            if(AvgDiffFile && !DoRSS)  SumDiff += diff;
            if(DiffLabelVolFile) {
              if (diff==0) MRIsetVoxVal(DiffLabelVol,c,r,s,f,val1, mode);
              else {
                MRIsetVoxVal(DiffLabelVol,c,r,s,f,SUSPICIOUS, mode);
              }
            }
            if (maxdiff < diff) {
              maxdiff = diff;
              cmax = c;
              rmax = r;
              smax = s;
              fmax = f;
            }
          }
          if(DiffVolFile && DoRSS) MRIsetVoxVal(DiffVol,c,r,s,0,sqrt(SumSqErr), mode);
          if(AvgDiffFile && DoRSS) SumDiff += sqrt(SumSqErr);
        }
      }
    }
    if(debug) printf("maxdiff %f at %d %d %d %d\n",
		     maxdiff,cmax,rmax,smax,fmax);
		     
    if(PrintSSD) printf("%f sum of squared differences\n",SumSqDiff);
    if(PrintRMS) printf("%f root mean squared differences\n",sqrt(SumSqDiff/notzero));
    if(CheckPixVals) printf("diffcount %d\n",ndiff);

    if(DiffVolFile) MRIwrite(DiffVol,DiffVolFile);      
    if(DiffLabelVolFile) MRIwrite(DiffLabelVol,DiffLabelVolFile);
    if(AvgDiffFile){
      int navg = InVol1->width * InVol1->height * InVol1->depth;
      if(! DoRSS) navg *= InVol1->nframes;
      double AvgDiff = SumDiff/navg;
      if(debug) printf("AvgStats %d %lf %lf\n",navg,SumDiff,AvgDiff);
      fp = fopen(AvgDiffFile,"w");
      fprintf(fp,"%lf\n",AvgDiff);
      fclose(fp);
    }

    if(maxdiff > pixthresh) {
      printf("Volumes differ in pixel data\n");
      printf("maxdiff %12.8f at %d %d %d %d\n",
             maxdiff,cmax,rmax,smax,fmax);
      if (DiffFile) {
        fp = fopen(DiffFile,"w");
        if (fp==NULL) {
          printf("ERROR: could not open %s\n",DiffFile);
          exit(1);
        }
        dump_options(fp);
        fprintf(fp,"maxdiff %12.8f at %d %d %d %d\n",maxdiff,cmax,rmax,smax,fmax);
        fprintf(fp,"Volumes differ in pixel value\n");
	if(CheckPixVals) fprintf(fp,"diffcount %d\n",ndiff);
        fclose(fp);
      }
      if(ndiff > diffcountthresh) {
	ExitStatus = 106;
	if (ExitOnDiff) exit(106);
      }
      else {
	printf("But diff count does not exceed diff count thresh %d\n",diffcountthresh);
      }
    }
  } // end of Compare pixel values

  //----------------------------------------------------------
  // Do a diff on a specific label in aseg.mgz
  if (SegDiff > -1) {
    MRI* SegDiffVol = MRIallocSequence(InVol1->width,InVol1->height,
                                      InVol1->depth,MRI_INT,InVol1->nframes);
    MRIcopyHeader(InVol1,SegDiffVol);
    int same = 0;
    int n1 = 0;
    int n2 = 0;
    int l1 = 0;
    int l2 = 0;
    printf("\nComputing difference for label %d ...\n",SegDiff);
    for (c=0; c < InVol1->width; c++) {
      for (r=0; r < InVol1->height; r++) {
        for (s=0; s < InVol1->depth; s++) {
          for (f=0; f < InVol1->nframes; f++) {
            val1 = MRIgetVoxVal(InVol1,c,r,s,f);
            val2 = MRIgetVoxVal(InVol2,c,r,s,f);
            if ((int)val1 == SegDiff &&(int)val2 == SegDiff )
            {
              MRIsetVoxVal(SegDiffVol,c,r,s,f,3);
              same++;
              l1++;
              l2++;
            }
            else if ((int)val1 != SegDiff &&(int)val2 != SegDiff )
              MRIsetVoxVal(SegDiffVol,c,r,s,f,0);
            else if ((int)val1 == SegDiff &&(int)val2 != SegDiff )
            {
              MRIsetVoxVal(SegDiffVol,c,r,s,f,1);
              n1++;
              l1++;
            }
            else if ((int)val1 != SegDiff &&(int)val2 == SegDiff )
            {
              MRIsetVoxVal(SegDiffVol,c,r,s,f,2);
              n2++;
              l2++;
            }
            else {//cannot happen 
              printf("ERROR: this error is not possible\n");
              exit(1);
            }
          }
        }
      }
    }
    if(SegDiffFile) MRIwrite(SegDiffVol,SegDiffFile);      
    MRIfree(&SegDiffVol);
    printf("\nDiff on label %d\n",SegDiff);
    printf("%d  identical\n",same);
    printf("%d  only in first\n",n1);
    printf("%d  only in second\n",n2);
    printf("%d - %d = %d difference (second-first)\n",l2,l1, l2-l1);
  }
  // Do a diff on all labels:
  // print label from vol1 if voxel label differes from vol2
  else if(SegDiff == -1) {
    MRI* SegDiffVol = MRIallocSequence(InVol1->width,InVol1->height,
                                      InVol1->depth,MRI_INT,InVol1->nframes);
    MRIcopyHeader(InVol1,SegDiffVol);
    printf("\nComputing difference for all labels ...\n");
    int same = 0;
    int different = 0;
    for (c=0; c < InVol1->width; c++) {
      for (r=0; r < InVol1->height; r++) {
        for (s=0; s < InVol1->depth; s++) {
          for (f=0; f < InVol1->nframes; f++) {
            val1 = MRIgetVoxVal(InVol1,c,r,s,f);
            val2 = MRIgetVoxVal(InVol2,c,r,s,f);
						if (val1 != val2) 
						{
						  MRIsetVoxVal(SegDiffVol,c,r,s,f,val1);
							different++;
						}
						else 
						{
						  MRIsetVoxVal(SegDiffVol,c,r,s,f,0);
							same++;
						}
          }
        }
      }
    }
    if(SegDiffFile) MRIwrite(SegDiffVol,SegDiffFile);      
    MRIfree(&SegDiffVol);
    printf("\nDiff:\n");
    printf("%d  identical\n",same);
    printf("%d  different\n",different);
  }
  
  
  //----------------------------------------------------------
  if (CheckOrientation) {
    MRIdircosToOrientationString(InVol1,Orient1);
    MRIdircosToOrientationString(InVol2,Orient2);
    if (strcmp(Orient1,Orient2)) {
      printf("Volumes differ in orientation %s %s\n",
             Orient1,Orient2);
      if (DiffFile) {
        fp = fopen(DiffFile,"w");
        if (fp==NULL) {
          printf("ERROR: could not open %s\n",DiffFile);
          exit(1);
        }
        dump_options(fp);
        fprintf(fp,"Volumes differ in orientation %s %s\n",
                Orient1,Orient2);
        fclose(fp);
      }
      ExitStatus = 107;
      if (ExitOnDiff) exit(107);
    }
  }

  if (CheckColortable) {
    printf("\nDiff embedded colortables ...\n");
    int nctabdiff_count = diff_colortable(InVol1->ct, InVol2->ct);
    if (nctabdiff_count > 0)
    {
      ExitStatus = 202;
      printf("****** Found %d differences between color tables\n", nctabdiff_count);
    }
    else
      printf("****** No difference found!\n");
  }

  exit(ExitStatus);
  return(ExitStatus);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--diag-debug")) Gdiag |= DIAG_INFO;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--verbose")) verbose = 1;
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;
    else if (!strcasecmp(option, "--no-exit-on-diff")) ExitOnDiff = 0;
    else if (!strcasecmp(option, "--no-exit-on-error")) ExitOnDiff = 0;

    else if (!strcasecmp(option, "--notallow-res"))  CheckResolution = 0;
    else if (!strcasecmp(option, "--notallow-acq"))  CheckAcqParams = 0;
    else if (!strcasecmp(option, "--notallow-geo"))  CheckGeo = 0;
    else if (!strcasecmp(option, "--notallow-prec")) CheckPrecision = 0;
    else if (!strcasecmp(option, "--notallow-type")) CheckType = 0;
    else if (!strcasecmp(option, "--notallow-pix"))  CheckPixVals = 0;
    else if (!strcasecmp(option, "--notallow-ori"))  CheckOrientation = 0;
    else if (!strcasecmp(option, "--no-absdiff"))  AbsDiff = 0;
    else if (!strcasecmp(option, "--absdiff"))     AbsDiff = 1; //default
    else if (!strcasecmp(option, "--diffabs"))  DiffAbs = 1;
    else if (!strcasecmp(option, "--diffpct"))  DiffPct = 1;
    else if (!strcasecmp(option, "--rss"))  DoRSS = 1;
    else if (!strcasecmp(option, "--ssd"))  PrintSSD = 1;
    else if (!strcasecmp(option, "--rms"))  {PrintRMS = 1;CheckPixVals=1;}
    else if (!strcasecmp(option, "--count"))  {CheckPixVals=1;}
    else if (!strcasecmp(option, "--diff-ctab")) CheckColortable = 1;
    else if (!strcasecmp(option, "--diff-imag")) CheckImag = 1;
    else if (!strcasecmp(option, "--qa")) {
      CheckPixVals = 0;
      CheckGeo     = 0;
    } else if (!strcasecmp(option, "--pix-only") || !strcasecmp(option, "--po") ||
	       !strcasecmp(option, "--pixonly")) {
      CheckPixVals = 1;
      CheckResolution = 0;
      CheckAcqParams = 0;
      CheckGeo     = 0;
      CheckPrecision = 0;
      CheckOrientation = 0;
    } else if (!strcasecmp(option, "--v1")) {
      if (nargc < 1) CMDargNErr(option,1);
      InVol1File = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--v2")) {
      if (nargc < 1) CMDargNErr(option,1);
      InVol2File = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&pixthresh);
      nargsused = 1;
    } else if (!strcasecmp(option, "--count-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&diffcountthresh);
      nargsused = 1;
    } else if (!strcasecmp(option, "--res-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&resthresh);
      CheckResolution=1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--geo-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&geothresh);
      CheckResolution=1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--rgb-thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&rgbthresh);
      nargsused = 1;
    } else if (!strcasecmp(option, "--diff-file") ||
               !strcasecmp(option, "--log")) {
      if (nargc < 1) CMDargNErr(option,1);
      DiffFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--diff")) {
      if (nargc < 1) CMDargNErr(option,1);
      DiffVolFile = pargv[0];
      CheckPixVals=1;
      nargsused = 1;
    } else if (!strcasecmp(option, "--avg-diff")) {
      if (nargc < 1) CMDargNErr(option,1);
      AvgDiffFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--diff_label_suspicious")) {
      if (nargc < 1) CMDargNErr(option,1);
      DiffLabelVolFile = pargv[0];
      nargsused = 1;
    } else if (!strcasecmp(option, "--merge-edits")) {
      // newauto oldauto manedits mergedvol
      if(nargc < 4) CMDargNErr(option,4);
      MRI *newauto, *oldauto, *manedit, *mergedvol;
      newauto = MRIread(pargv[0]);
      if(newauto==NULL) exit(1);
      oldauto = MRIread(pargv[1]);
      if(oldauto==NULL) exit(1);
      manedit = MRIread(pargv[2]);
      if(manedit==NULL) exit(1);
      mergedvol = MRIapplyEdits(newauto, oldauto, manedit, NULL);
      if(mergedvol==NULL) exit(1);
      int err = MRIwrite(mergedvol,pargv[3]);
      exit(err);
      nargsused = 1;
    } else if (!strcasecmp(option, "--segdiff")) {
      if (nargc < 2) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&SegDiff);
      SegDiffFile = pargv[1];
      CheckPixVals     = 0;
      CheckResolution  = 1;
      CheckAcqParams   = 1;
      CheckGeo         = 0;
      CheckPrecision   = 1;
      CheckOrientation = 1;
      nargsused = 2;
    } else {
      if (InVol1File == NULL)      InVol1File = option;
      else if (InVol2File == NULL) InVol2File = option;
      else {
        fprintf(stderr,"ERROR: Option %s unknown\n",option);
        if (CMDsingleDash(option))
          fprintf(stderr,"       Did you really mean -%s ?\n",option);
        exit(-1);
      }
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {
  if (InVol1File==NULL) {
    printf("ERROR: need to spec volume 1\n");
    exit(1);
  }
  if (InVol2File==NULL) {
    printf("ERROR: need to spec volume 2\n");
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
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"%s\n",Progname);
  fprintf(fp,"FREESURFER_HOME %s\n",getenv("FREESURFER_HOME"));
  fprintf(fp,"cwd       %s\n",cwd);
  fprintf(fp,"cmdline   %s\n",cmdline);
  fprintf(fp,"timestamp %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname   %s\n",uts.sysname);
  fprintf(fp,"hostname  %s\n",uts.nodename);
  fprintf(fp,"machine   %s\n",uts.machine);
  fprintf(fp,"user      %s\n",VERuser());
  fprintf(fp,"v1        %s\n",InVol1File);
  fprintf(fp,"v2        %s\n",InVol2File);
  fprintf(fp,"pixthresh %lf\n",pixthresh);
  fprintf(fp,"resthresh %lf\n",resthresh);
  fprintf(fp,"geothresh %lf\n",geothresh);
  fprintf(fp,"checkres  %d\n",CheckResolution);
  fprintf(fp,"checkacq  %d\n",CheckAcqParams);
  fprintf(fp,"checkpix  %d\n",CheckPixVals);
  fprintf(fp,"checkgeo  %d\n",CheckGeo);
  fprintf(fp,"checkori  %d\n",CheckOrientation);
  fprintf(fp,"checkprec %d\n",CheckPrecision);
  fprintf(fp,"absdiff   %d\n",AbsDiff);
  fprintf(fp,"diffabs   %d\n",DiffAbs);
  fprintf(fp,"diffpct   %d\n",DiffPct);
  fprintf(fp,"rss       %d\n",DoRSS);
  fprintf(fp,"ssd       %d\n",PrintSSD);
  fprintf(fp,"rms       %d\n",PrintRMS);
  fprintf(fp,"logfile   %s\n",DiffFile);
  return;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s <options> vol1file vol2file <options> \n",Progname) ;
  printf("\n");
  //printf("   --v1 volfile1 : first  input volume \n");
  //printf("   --v2 volfile2 : second input volume \n");
  printf("\n");
  printf("   --notallow-res  : do not check for resolution diffs\n");
  printf("   --notallow-acq  : do not check for acq param diffs\n");
  printf("   --notallow-geo  : do not check for geometry diffs\n");
  printf("   --notallow-prec : do not check for precision diffs\n");
  printf("   --notallow-pix  : do not check for pixel diffs\n");
  printf("   --notallow-ori  : do not check for orientation diffs\n");
  printf("   --notallow-type  : do not check for file type diffs\n");
  printf("   --no-exit-on-diff : do not exit on diff "
         "(runs thru everything)\n");
  printf("\n");
  printf("   --qa         : check res, acq, precision, "
         "and orientation only\n");
  printf("   --pix-only   : only check pixel data\n");
  printf("   --absdiff    : take abs of diff (default)\n");
  printf("   --no-absdiff : do not take abs of diff\n");
  printf("   --diffabs    : take abs before computing diff\n");
  printf("   --diffpct    : 100*(v1-v2)/((v1+v2)/2)\n");
  printf("   --rss        : save sqrt sum squares with --diff\n");
  printf("   --ssd        : print sum squared differences over all voxel\n");
  printf("   --rms        : print root mean squared diff. over all non-zero voxel\n");
  printf("   --count      : print number of differing voxels\n");
  printf("\n");
  printf("   --thresh thresh     : pix diffs must be greater than this \n");
  printf("   --res-thresh thresh : resolution diffs must be greater than this \n");
  printf("   --geo-thresh thresh : geometry diffs must be greater than this \n");
  printf("   --count-thresh nvox : there must be at least this many voxels that are diff\n");
  printf("   --log DiffFile : store diff info in this file. \n");
  printf("   --diff DiffVol : save difference image. \n");
  printf("   --diff_label_suspicious DiffVol : differing voxels replaced\n");
  printf("                                     with label SUSPICIOUS\n");
  printf("                                     (for comparing aseg.mgz's)\n\n");
  printf("   --segdiff labelIDX DiffASEG : diff on voxels with labelIDX\n");
  printf("                                 output image: 0 not in both,\n");
  printf("                                 1 only in 1st, 2 only in 2nd\n");
  printf("                                 3 in both (for aseg.mgz)\n");
  printf("                                 if labelIDX==-1 diff on all labels:\n");
  printf("                                  show labels from vol1 at voxels\n");
  printf("                                  that differ in vol2\n\n");
  printf("   --merge-edits newauto oldauto manedits mergedvol\n");
  printf("   --avg-diff avgdiff.txt : save average difference \n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --verbose   print out info on all differences found\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf("\n");
  printf("Determines whether two volumes differ. See below for what \n");
  printf("'differ' means. \n");
  printf("\n");
  printf("The basic usage is something like:\n");
  printf("\n");
  printf("mri_diff vol1 vol2\n");
  printf("\n");
  printf("It then prints to the terminal whether they differ or not.\n");
  printf("\n");
  printf("NOTE: stuff might get printed to the terminal regardless "
         "of whether\n");
  printf("the volumes are different or not, so you cannot just use the \n");
  printf("presence or absence of terminal output to determine whether they\n");
  printf("are different as you can with unix diff.\n");
  printf("\n");
  printf("There are three ways to determine whether they are different:\n");
  printf("1. Look for 'volumes differ' in the terminal output\n");
  printf("2. Program exits with status > 100\n");
  printf("3. Create a log file\n");
  printf("\n");
  printf("To create a log file, add --log yourlogfile any where in the\n");
  printf("command line, eg:\n");
  printf("\n");
  printf("mri_diff vol1 vol2 --log yourlogfile\n");
  printf("\n");
  printf("If yourlogfile exists, it will be immediately deleted when the \n");
  printf("program starts. If a difference is detected, yourlogfile will\n");
  printf("be created, and information about the difference will be\n");
  printf("printed to it. If there is no difference, yourlogfile will\n");
  printf("not exist when the program finishes.\n");
  printf("\n");
  printf("Volumes can differ in six ways:\n");
  printf("  1. Dimension,               return status = 101\n");
  printf("  2. Resolutions,             return status = 102\n");
  printf("  3. Acquisition Parameters,  return status = 103\n");
  printf("  4. Geometry,                return status = 104\n");
  printf("  5. Precision,               return status = 105\n");
  printf("  6. Pixel Data,              return status = 106\n");
  printf("  7. m3z 3D Morph,            return status = 201\n");
  printf("  8. Color Talbe,             return status = 202\n");
  printf("\n");
  printf("Dimension is number of rows, cols, slices, and frames.\n");
  printf("Resolution is voxel size.\n");
  printf("Acqusition parameters are: flip angle, TR, TE, and TI.\n");
  printf("Geometry checks the vox2ras matrices for equality.\n");
  printf("Precision is int, float, short, etc.\n");
  printf("\n");
  printf("By default, all of these are checked, but some can be turned off\n");
  printf("with certain command-line flags:\n");
  printf("\n");
  printf("--notallow-res  : turns off resolution checking\n");
  printf("--notallow-acq  : turns off acquistion parameter checking\n");
  printf("--notallow-geo  : turns off geometry checking\n");
  printf("--notallow-prec : turns off precision checking\n");
  printf("--notallow-pix  : turns off pixel checking\n");
  printf("\n");
  printf("In addition, the minimum difference in pixel value required\n");
  printf("to be considered different can be controlled with --thresh.\n");
  printf("Eg, if two volumes differ by .00001, and you do not consider\n");
  printf("this to be significant, then --thresh .00002 will prevent\n");
  printf("that difference from being considered different. The default\n");
  printf("threshold is 0.\n");
  printf("\n");
  printf("--diff diffvol\n");
  printf("\n");
  printf("Saves difference image to diffvol.\n");
  printf("\n");
  printf("QUALITY ASSURANCE\n");
  printf("\n");
  printf("mri_diff can be used to check that two volumes where "
         "acquired in the\n");
  printf("same way with the --qa flag. This turns on Res, Acq, and "
         "Prec, and\n");
  printf("turns off Geo and Pix. Instead of checking geometry, it checks "
         "for the\n");
  printf("basic orientation of the volumes (eg, RAS, LPI, etc). The idea "
         "here is\n");
  printf("that the pixel data and exact geometry may be different, "
         "but other\n");
  printf("things should be the same.\n");
  printf("\n");
  printf("EXIT STATUS\n");
  printf("\n");
  printf("0   Volumes are not different and there were no errors\n");
  printf("1   Errors encounted. Volumes may or may not be different\n");
  printf("101 Volumes differ in dimension\n");
  printf("102 Volumes differ in resolution\n");
  printf("103 Volumes differ in acquisition parameters\n");
  printf("104 Volumes differ in geometry\n");
  printf("105 Volumes differ in precision\n");
  printf("106 Volumes differ in pixel data\n");
  printf("107 Volumes differ in orientation\n");
  printf("201 m3z 3D Morph differ\n");
  printf("202 Volumes differ in Color Table\n");
  printf("\n");

  exit(1) ;
}

static void diff_mgh_morph(const char *file1, const char *file2)
{
  printf("Diff MGH_MORPH %s and %s using threshold = %.10f\n", file1, file2, pixthresh);

  GCAM *gcam1 = GCAMread(file1);
  GCAM *gcam2 = GCAMread(file2);

  if (gcam1->width  != gcam2->width  || 
      gcam1->height != gcam2->height ||
      gcam1->depth  != gcam2->depth)
    printf("Dimensions differ: (%d x %d x %d) vs (%d x %d x %d)\n",
           gcam1->width, gcam1->height, gcam1->depth, gcam2->width, gcam2->height, gcam2->depth);

  if (gcam1->spacing != gcam2->spacing)
    printf("Spacing differ: %d vs %d\n", gcam1->spacing, gcam2->spacing);

  if (gcam1->exp_k != gcam2->exp_k)
    printf("exp_k differ: %.2f vs %.2f\n", gcam1->exp_k, gcam2->exp_k);

  double maxdiffx = 0.0, maxdiffy = 0.0, maxdiffz = 0.0;
  int ndifforigx = 0, ndiffx = 0, ndiffxn = 0, ndiffinvalid = 0, ndifflabel = 0;

  int width  = (gcam1->width  <= gcam2->width)  ? gcam1->width  : gcam2->width;
  int height = (gcam1->height <= gcam2->height) ? gcam1->height : gcam2->height;
  int depth  = (gcam1->depth  <= gcam2->depth)  ? gcam1->depth  : gcam2->depth;

  GCA_MORPH_NODE *gcamn1, *gcamn2;
  for (int x = 0; x < width; x++)
  {
    for (int y = 0; y < height; y++)
    {
      for (int z = 0; z < depth; z++)
      {
        gcamn1 = &gcam1->nodes[x][y][z];
        gcamn2 = &gcam2->nodes[x][y][z];

        if (gcamn1->origx != gcamn2->origx ||
            gcamn1->origy != gcamn2->origy ||
            gcamn1->origz != gcamn2->origz)
	{
          ndifforigx++;
           
          if (verbose)
            printf("(origx, origy, origz) differ at (%03d,%03d,%03d): (%.6f %.6f %.6f) vs (%.6f %.6f %.6f)\n", 
                   x, y, z, gcamn1->origx, gcamn1->origy, gcamn1->origz, gcamn2->origx, gcamn2->origy, gcamn2->origz);
        }

	double diffx = fabs(gcamn1->x - gcamn2->x);
	double diffy = fabs(gcamn1->y - gcamn2->y);
	double diffz = fabs(gcamn1->z - gcamn2->z);
        if (diffx > pixthresh ||
            diffy > pixthresh ||
            diffz > pixthresh)
	{
          ndiffx++;

	  if (diffx > maxdiffx)
	  {
            maxdiffx = diffx;
	  }
	  if (diffy > maxdiffy)
	  {
            maxdiffy = diffy ;
	  }
	  if (diffz > maxdiffz)
	  {
            maxdiffz = diffz;
	  }
	  
          if (verbose)
	    printf("(x, y, z) differ at (%03d,%03d,%03d): (%.6f %.6f %.6f) vs (%.6f %.6f %.6f)\n", 
                   x, y, z, gcamn1->x, gcamn1->y, gcamn1->z, gcamn2->x, gcamn2->y, gcamn2->z);
        }

        if (gcamn1->xn != gcamn2->xn ||
            gcamn1->yn != gcamn2->yn ||
            gcamn1->zn != gcamn2->zn)
	{
          ndiffxn++;

          if (verbose)
            printf("(xn, yn, zn) differ at (%03d,%03d,%03d): (%d %d %d) vs (%d %d %d)\n", 
                   x, y, z, gcamn1->xn, gcamn1->yn, gcamn1->zn, gcamn2->xn, gcamn2->yn, gcamn2->zn);
        }

        if (gcamn1->invalid != gcamn2->invalid)
	{
          ndiffinvalid++;

          if (verbose)
            printf("Invalid differ at (%03d,%03d,%03d): %c vs %c\n", 
                   x, y, z, gcamn1->invalid, gcamn2->invalid);
	}

        if (gcamn1->label != gcamn2->label)
	{
          ndifflabel++;

          if (verbose)
            printf("Labels differ at (%03d,%03d,%03d): %d vs %d\n", 
                   x, y, z, gcamn1->label, gcamn2->label);
        }
      }
    }
  }

  if (ndifforigx || ndiffx || ndiffxn || ndiffinvalid || ndifflabel)
  {
    ExitStatus = 201;
    
    if (ndifforigx)
      printf("(origx, origy, origz) diff counts = %d\n", ndifforigx);
    if (ndiffx)
    {
      printf("(x, y, z)             diff counts = %d\n", ndiffx);
      printf("                      max diffx = %.6f, max diffy = %.6f, max diffz = %.6f\n", maxdiffx, maxdiffy, maxdiffz); 
    }
    if (ndiffxn)
      printf("(xn, yn, zn)          diff counts = %d\n", ndiffxn);
    if (ndiffinvalid)
      printf("invalid flag          diff counts = %d\n", ndiffinvalid);
    if (ndifflabel)
      printf("label flag            diff counts = %d\n", ndifflabel);

    if (!verbose)
      printf("use --verbose to print the details\n");
  }

  // TAG_GCAMORPH_LABELS
  if (gcam1->status != gcam2->status)
  {
    ExitStatus = 201;
    printf("Status differ:  %d vs %d\n", gcam1->status, gcam2->status);
  }

  // TAG_GCAMORPH_GEOM
  int same = vg_isEqual(&gcam1->image, &gcam2->image);
  if (!same)
  {
    ExitStatus = 201;
    
    printf("image VOL_GEOMs differ:\n");
    printf("%s geometry:\n", file1);
    writeVolGeom(stdout, &gcam1->image);
    printf("%s geometry:\n", file2);
    writeVolGeom(stdout, &gcam2->image);
  }
  
  same = vg_isEqual(&gcam1->atlas, &gcam2->atlas);
  if (!same)
  {
    ExitStatus = 201;
    
    printf("atlas VOL_GEOMs differ:\n");
    printf("%s geometry:\n", file1);
    writeVolGeom(stdout, &gcam1->atlas);
    printf("%s geometry:\n", file2);
    writeVolGeom(stdout, &gcam2->atlas);
  }

  // TAG_GCAMORPH_TYPE
  if (gcam1->type != gcam2->type)
  {
    ExitStatus = 201;
    printf("Types differ: %d vs %d\n", gcam1->type, gcam2->type);
  }

  // TAG_MGH_XFORM
  if (gcam1->m_affine == NULL && gcam2->m_affine == NULL)
    return;

  if (gcam1->m_affine != NULL && gcam2->m_affine != NULL)
  {
    for (int r = 1; r <= 4; r++) {
      for (int c = 1; c <= 4; c++) {
        double val1 = gcam1->m_affine->rptr[r][c];
        double val2 = gcam2->m_affine->rptr[r][c];
        double diff = fabs(val1-val2);
        if (diff > geothresh) 
        {
          printf("m_affine matrix differ in geometry row=%d col=%d diff=%lf (%g)\n",
                 r, c, diff, diff);
          
          printf("v1 m_affine (det = %.3f) ----------\n", gcam1->det);
          MatrixPrint(stdout, gcam1->m_affine);
          printf("v2 m_affine (det = %.3f) ---------\n", gcam2->det);
          MatrixPrint(stdout, gcam2->m_affine);

          ExitStatus = 201;
          if (ExitOnDiff) exit(ExitStatus);
        }
      } // c
    } // r
  }
  else if (gcam1->m_affine == NULL)
  {
    ExitStatus = 201;
    printf("%s doesn't contain m_affine matrix\n", file1);
  }
  else if (gcam2->m_affine == NULL)
  {
    ExitStatus = 201;
    printf("%s doesn't contain m_affine matrix\n", file2);
  }
}


static int diff_colortable(const COLOR_TABLE *ctab1, const COLOR_TABLE *ctab2)
{
  if (ctab1 == NULL || ctab2 == NULL)
  {
    printf("[INFO] No colortables to diff\n");
    return 0;
  }

  
  int nctabdiffs = 0;
  int todiff = ctab1->nentries;
  if (ctab1->nentries != ctab2->nentries)
  {
    todiff = (ctab1->nentries < ctab2->nentries) ? ctab1->nentries : ctab2->nentries;
    printf("[WARN] colortab number of entries (%d and %d) do not match. First %d entries will be compared.\n",
           ctab1->nentries, ctab2->nentries, todiff);
  }

  for (int n = 0; n < todiff; n++)
  {
    COLOR_TABLE_ENTRY *cte1 = ctab1->entries[n];
    COLOR_TABLE_ENTRY *cte2 = ctab2->entries[n];
    if ((cte1 == NULL && cte2 != NULL) || (cte1 != NULL && cte2 == NULL))
    {
      nctabdiffs++;
      if (verbose)
        printf("[INFO] differ at entry %5d, colortable #%s has no entry\n", n, (cte1 == NULL) ? "1" : "2");
    }
    else if (cte1 != NULL && cte2 != NULL)
    {
      if (strcmp(cte1->name, cte2->name) != 0)
      {
        nctabdiffs++;
        if (verbose)
          printf("[INFO] differ at entry %5d, label name (%s) vs (%s)\n", n, cte1->name, cte2->name);
      }
      else if (cte1->ri != cte2->ri ||
               cte1->gi != cte2->gi ||
               cte1->bi != cte2->bi ||
               cte1->ai != cte2->ai)
      {
        nctabdiffs++;
        if (verbose)
	{
          printf("[INFO] differ at entry %5d, RGBi (%d, %d, %d, %d)\n", n, cte1->ri, cte1->gi, cte1->bi, cte1->ai);
          printf("                                   (%d, %d, %d, %d)\n",  cte2->ri, cte2->gi, cte2->bi, cte2->ai);
	}
      }
#if 0 
      else if (fabs(cte1->rf - cte2->rf) > rgbthresh ||
               fabs(cte1->gf - cte2->gf) > rgbthresh ||
               fabs(cte1->bf - cte2->bf) > rgbthresh ||
               fabs(cte1->af - cte2->af) > rgbthresh)
      {
        nctabdiffs++;
        if (verbose)
	{
          printf("[INFO] differ at entry %5d, RGBf (%.8f, %.8f, %.8f, %.8f)\n", n, cte1->rf, cte1->gf, cte1->bf, cte1->af);
	  printf("                                   (%.8f, %.8f, %.8f, %.8f)\n",  cte2->rf, cte2->gf, cte2->bf, cte2->af);
	}
      }
#endif      
    }
  }

  return nctabdiffs;
}
