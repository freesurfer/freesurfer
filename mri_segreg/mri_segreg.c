/**
 * @file  mri_segreg.c
 * @brief program for computing/optimizing cost function of segmentation-based registration
 *        
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: Greg Grev
 * CVS Revision Info:
 *    $Author: greve $
 *    $Date: 2007/09/14 21:05:56 $
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

mri_segreg
  --reg regfile
  --f fvol
  --tx-mmd txmin txmax txdelta
  --tx 
  --tx-mmd txmin txmax txdelta
  --tx 
  --tx-mmd txmin txmax txdelta
  --tx 


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

mri_vol2vol mri_convert, tkregister2


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

#ifdef X
#undef X
#endif

// For some reason, this does not seemed to be defined in math.h
double round(double x);

static int  parse_commandline(int argc, char **argv);
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

static char vcid[] = "$Id: mri_segreg.c,v 1.2 2007/09/14 21:05:56 greve Exp $";
char *Progname = NULL;

int debug = 0, gdiagno = -1;

char *movvolfile=NULL;
char *regfile=NULL;
char *interpmethod = "trilinear";
int   interpcode = 0;
int   sinchw;

MRI *mov, *out;

MATRIX *R, *invR;
MATRIX *vox2vox, *vox2ras;
MATRIX *Tin, *invTin, *Sin, *invSin;
MATRIX *Ttemp, *invTtemp, *Stemp, *invStemp;

char *SUBJECTS_DIR=NULL;
char *subject = NULL;

float ipr, bpr, intensity;
int float2int,err, nargs;

char tmpstr[2000];

MATRIX *Mrot = NULL;
MATRIX *Mtrans = NULL;

char *SegRegCostFile = NULL;
char  *fspec;
MRI *regseg;

#define NMAX 100
int ntx=0, nty=0, ntz=0, nax=0, nay=0, naz=0;
double txlist[NMAX],tylist[NMAX],tzlist[NMAX];
double axlist[NMAX],aylist[NMAX],azlist[NMAX];

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  char cmdline[CMD_LINE_LEN] ;
  double costs[8];
  FILE *fp;

  make_cmd_version_string(argc, argv,
                          "$Id: mri_segreg.c,v 1.2 2007/09/14 21:05:56 greve Exp $",
                          "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option(argc, argv,
                                "$Id: mri_segreg.c,v 1.2 2007/09/14 21:05:56 greve Exp $",
                                "$Name:  $");
  if(nargs && argc - nargs == 1) exit (0);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  if(gdiagno > -1) Gdiag_no = gdiagno;
  check_options();
  dump_options(stdout);

  sprintf(tmpstr,"%s/%s/mri/regseg",SUBJECTS_DIR,subject);
  fspec = IDnameFromStem(tmpstr);
  regseg = MRIread(fspec);
  if(regseg == NULL) exit(1);
  free(fspec);

  mov = MRIread(movvolfile);
  if (mov == NULL) exit(1);

  // Allocate the output
  out = MRIcloneBySpace(regseg,-1,mov->nframes);
  out = MRIallocSequence(regseg->width, regseg->height, regseg->depth, MRI_FLOAT, 1);
  MRIcopyHeader(regseg,out);




  // Vox-to-tkRAS Matrices
  Tin      = MRIxfmCRS2XYZtkreg(mov);
  invTin   = MatrixInverse(Tin,NULL);
  Ttemp    = MRIxfmCRS2XYZtkreg(regseg);
  invTtemp = MatrixInverse(Ttemp,NULL);

  // Vox-to-ScannerRAS Matrices
  Sin      = MRIxfmCRS2XYZ(mov,0);
  invSin   = MatrixInverse(Sin,NULL);
  Stemp    = MRIxfmCRS2XYZ(regseg,0);
  invStemp = MatrixInverse(Stemp,NULL);

  // Only gets here if resampling
  // vox2vox converts a template vox to input vox
  // vox2vox = invTin * R * Ttemp
  vox2vox = MatrixMultiply(invTin,R,NULL);
  MatrixMultiply(vox2vox,Ttemp,vox2vox);

  printf("\n");
  printf("Vox2Vox Matrix is:\n");
  MatrixPrint(stdout,vox2vox);
  printf("\n");

  printf("Resampling\n");
  MRIvol2Vol(mov,out,vox2vox,interpcode,sinchw);


  SegRegCost(regseg,out,costs);
  fp = fopen(SegRegCostFile,"a");
  fclose(fp);

  printf("\n");
  printf("mri_segreg done\n");

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
    else if (istringnmatch(option, "--mov",0)) {
      if (nargc < 1) argnerr(option,1);
      movvolfile = pargv[0];
      nargsused = 1;
    } else if (istringnmatch(option, "--reg",0)) {
      if (nargc < 1) argnerr(option,1);
      regfile = pargv[0];
      err = regio_read_register(regfile, &subject, &ipr, &bpr,
                                &intensity, &R, &float2int);
      if (err) exit(1);
      invR = MatrixInverse(R,NULL);
      nargsused = 1;
    } else if (istringnmatch(option, "--tx-mmd",0)) {
      if (nargc < 3) argnerr(option,3);
      scanf(argv[0],"%lf",&vmin);
      scanf(argv[1],"%lf",&vmax);
      scanf(argv[2],"%lf",&vdelta);
      nv = (int)floor((vmax-vmin)/vdelta);
      if(nv == 0) exit(1);
      for(n=0; n < nv; n++){
	txlist[ntx] = vmin + vdelta*n;
	ntx++;
      }
      nargsused = 3;
    } else if ( !strcmp(option, "--gdiagno") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&gdiagno);
      nargsused = 1;
    } else if (istringnmatch(option, "--cost",0)) {
      if (nargc < 1) argnerr(option,1);
      SegRegCostFile = pargv[0];
      nargsused = 1;
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
  if (movvolfile == NULL) {
    printf("ERROR: No mov volume supplied.\n");
    exit(1);
  }
  if (regfile == NULL) {
    printf("ERROR: need --reg.\n");
    exit(1);
  }

  interpcode = MRIinterpCode(interpmethod);
  if (interpcode < 0) {
    printf("ERROR: interpolation method %s unrecognized\n",interpmethod);
    printf("       legal values are nearest, trilin, and sinc\n");
    exit(1);
  }

  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"movvol %s\n",movvolfile);
  fprintf(fp,"regfile %s\n",regfile);
  fprintf(fp,"interp  %s (%d)\n",interpmethod,interpcode);
  if(interpcode == SAMPLE_SINC) fprintf(fp,"sinc hw  %d\n",sinchw);
  fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);
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
