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
 *    $Date: 2007/09/14 22:31:48 $
 *    $Revision: 1.3 $
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
  --mov fvol
  --cost costfile
  --tx-mmd txmin txmax txdelta
  --ty-mmd tymin tymax tydelta
  --tz-mmd tzmin tzmax tzdelta
  --ax-mmd axmin axmax axdelta
  --ay-mmd aymin aymax aydelta
  --az-mmd azmin azmax azdelta

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

static char vcid[] = "$Id: mri_segreg.c,v 1.3 2007/09/14 22:31:48 greve Exp $";
char *Progname = NULL;

int debug = 0, gdiagno = -1;

char *movvolfile=NULL;
char *regfile=NULL;
char *interpmethod = "trilinear";
int   interpcode = 0;
int   sinchw;

MRI *mov, *out;

MATRIX *R0;

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
  double costs[8], angles[3];
  double tx, ty, tz, ax, ay, az; 
  int nth,nthtx, nthty, nthtz, nthax, nthay, nthaz; 
  FILE *fp;
  MATRIX *Tin, *invTin, *Sin, *invSin;
  MATRIX *Ttemp, *invTtemp, *Stemp, *invStemp;
  MATRIX *R=NULL, *invR=NULL, *vox2vox=NULL;

  make_cmd_version_string(argc, argv,
                          "$Id: mri_segreg.c,v 1.3 2007/09/14 22:31:48 greve Exp $",
                          "$Name:  $", cmdline);

  /* rkt: check for and handle version tag */
  nargs = handle_version_option(argc, argv,
                                "$Id: mri_segreg.c,v 1.3 2007/09/14 22:31:48 greve Exp $",
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

  printf("Loading regseg\n");
  sprintf(tmpstr,"%s/%s/mri/regseg",SUBJECTS_DIR,subject);
  fspec = IDnameFromStem(tmpstr);
  regseg = MRIread(fspec);
  if(regseg == NULL) exit(1);
  free(fspec);

  printf("Loading mov\n");
  mov = MRIread(movvolfile);
  if (mov == NULL) exit(1);

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
  
  // Allocate the output
  out = MRIcloneBySpace(regseg,-1,mov->nframes);
  out = MRIallocSequence(regseg->width, regseg->height, regseg->depth, MRI_FLOAT, 1);
  MRIcopyHeader(regseg,out);

  printf("Staring loop\n");
  nth = 0;
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

	      angles[0] = ax*(M_PI/180);
	      angles[1] = ay*(M_PI/180);
	      angles[2] = az*(M_PI/180);
	      Mrot = MRIangles2RotMat(angles);

	      Mtrans = MatrixIdentity(4,NULL);
	      Mtrans->rptr[1][4] = tx;
	      Mtrans->rptr[2][4] = ty;
	      Mtrans->rptr[3][4] = tz;

	      // R = Mtrans*Mrot*R0
	      R = MatrixMultiply(Mrot,R0,R);
	      R = MatrixMultiply(Mtrans,R,R);
	      invR = MatrixInverse(R,invR);

	      // vox2vox = invTin*R*Ttemp
	      vox2vox = MatrixMultiply(invTin,R,vox2vox);
	      MatrixMultiply(vox2vox,Ttemp,vox2vox);
	      
	      // resample
	      MRIvol2Vol(mov,out,vox2vox,interpcode,sinchw);
	      
	      // compute costs
	      SegRegCost(regseg,out,costs);

	      // write costs to file
	      fp = fopen(SegRegCostFile,"a");
	      fprintf(fp,"%7.3lf %7.3lf %7.3lf ",tx,ty,tz);
	      fprintf(fp,"%5.1lf %5.1lf %5.1lf ",ax,ay,az);
	      fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
	      fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
	      fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
	      fprintf(fp,"\n");
	      fclose(fp);

	      fp = stdout;
	      fprintf(fp,"%5d ",nth);
	      fprintf(fp,"%7.3lf %7.3lf %7.3lf ",tx,ty,tz);
	      fprintf(fp,"%5.1lf %5.1lf %5.1lf ",ax,ay,az);
	      fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[0],costs[1],costs[2]); // WM  n mean std
	      fprintf(fp,"%7d %10.4lf %8.4lf ",(int)costs[3],costs[4],costs[5]); // CTX n mean std
	      fprintf(fp,"%8.4lf %8.4lf ",costs[6],costs[7]); // t, cost=1/t
	      fprintf(fp,"\n");
	      fflush(stdout);

	      // clean up
	      MatrixFree(&Mrot);
	      MatrixFree(&Mtrans);
	    }
	  }
	}
      }
    }
  }

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
                                &intensity, &R0, &float2int);
      if (err) exit(1);
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
printf("\n");
printf("\n");
printf("mri_segreg\n");
printf("  --reg regfile\n");
printf("  --mov fvol\n");
printf("  --cost costfile\n");
printf("  --tx-mmd txmin txmax txdelta\n");
printf("  --ty-mmd tymin tymax tydelta\n");
printf("  --tz-mmd tzmin tzmax tzdelta\n");
printf("  --ax-mmd axmin axmax axdelta\n");
printf("  --ay-mmd aymin aymax aydelta\n");
printf("  --az-mmd azmin azmax azdelta\n");
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
  if (movvolfile == NULL) {
    printf("ERROR: No mov volume supplied.\n");
    exit(1);
  }
  if (regfile == NULL) {
    printf("ERROR: need --reg.\n");
    exit(1);
  }
  if (SegRegCostFile == NULL) {
    printf("ERROR: need --cost.\n");
    exit(1);
  }

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
  int n;
  fprintf(fp,"movvol %s\n",movvolfile);
  fprintf(fp,"regfile %s\n",regfile);
  fprintf(fp,"interp  %s (%d)\n",interpmethod,interpcode);
  if(interpcode == SAMPLE_SINC) fprintf(fp,"sinc hw  %d\n",sinchw);
  fprintf(fp,"Gdiag_no  %d\n",Gdiag_no);
  fprintf(fp,"ntx %d\n",ntx);
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
