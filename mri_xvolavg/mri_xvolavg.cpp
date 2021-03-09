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


/*----------------------------------------------------------
  Name: mri_xvolavg
  Author: Douglas Greve
  Purpose: averages multiple volumes together into a single
  volume. The volumes can be 4D.
  ----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri_identify.h"
#include "version.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static int  singledash(char *flag);

int main(int argc, char *argv[]) ;

const char *Progname = NULL;

char *defaulttypestring;
int  defaulttype     = MRI_VOLUME_TYPE_UNKNOWN;

char *volid[500];
int   nvols = 0;
char *voltypestring = NULL;
int   voltype       = MRI_VOLUME_TYPE_UNKNOWN;

char *outvolid      = NULL;
char *outtypestring = NULL;
int   outtype       = MRI_VOLUME_TYPE_UNKNOWN;

float val1, val2;

int debug = 0;

MRI *InVol, *OutVol, *OutVoltmp;

FILE *fp;

char tmpstr[2000];

/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
/*------------------------------------------------------------------*/
int main(int argc, char **argv) {
  int r,c,s,f;
  int nrows=0, ncols=0, nslices=0, nframes=0;
  int nthvol;
  int nargs;

  nargs = handleVersionOption(argc, argv, "mri_xvolavg");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();


  for (nthvol = 0; nthvol < nvols; nthvol++) {
    printf("%d/%d ------------------\n",nthvol+1,nvols);

    /* Load the input Volume */
    printf("INFO: Loading vol %s\n",volid[nthvol]);
    InVol = MRIreadType(volid[nthvol],voltype);
    if (InVol == NULL) {
      printf("ERROR: could not read %s as type %d\n",volid[nthvol],voltype);
      exit(1);
    }
    printf("INFO; Done loading volume\n");

    if (InVol->type != MRI_FLOAT) {
      printf("INFO: changing type to float\n");
      InVol = MRISeqchangeType(InVol,MRI_FLOAT,0,0,0);
    }

    if (nthvol == 0) {
      ncols   = InVol->width;
      nrows   = InVol->height;
      nslices = InVol->depth;
      nframes = InVol->nframes;
      printf("Input dimension: %d %d %d %d\n",ncols,nrows,nslices,nframes);
      //OutVol  = MRIallocSequence(ncols,nrows,nslices, MRI_FLOAT,nframes);
      OutVol  = MRIclone(InVol,NULL);
      if (OutVol == NULL) {
        printf("ERROR: could not alloc output volume\n");
        exit(1);
      }
    } else {
      /* Check that the two have the same dimension */
      if (InVol->width != ncols) {
        printf("Widths do not match: %d %d\n",InVol->width,ncols);
        exit(1);
      }
      if (InVol->height != nrows) {
        printf("Heights do not match: %d %d\n",InVol->height,nrows);
        exit(1);
      }
      if (InVol->depth != nslices) {
        printf("Depths do not match: %d %d\n",InVol->depth,nslices);
        exit(1);
      }
      if (InVol->nframes != nframes) {
        printf("Frames do not match: %d %d\n",InVol->nframes,nframes);
        exit(1);
      }
    }

    printf("Accumulating\n");
    for (c=0; c < ncols;   c++) {
      for (r=0; r < nrows;   r++) {
        for (s=0; s < nslices; s++) {
          for (f=0; f < nframes; f++) {
            MRIFseq_vox(OutVol,c,r,s,f) += MRIFseq_vox(InVol,c,r,s,f);
          }
        }
      }
    }

    MRIfree(&InVol);
  } /* end loop over input volumes */

  printf("Averaging\n");
  for (c=0; c < ncols;   c++) {
    for (r=0; r < nrows;   r++) {
      for (s=0; s < nslices; s++) {
        for (f=0; f < nframes; f++) {
          MRIFseq_vox(OutVol,c,r,s,f) /= nvols;
        }
      }
    }
  }

  if (outtype == MRI_CORONAL_SLICE_DIRECTORY) {
    printf("INFO: changing type to UCHAR\n");
    OutVoltmp = MRISeqchangeType(OutVol, MRI_UCHAR, 0, 255,1);
    if (OutVoltmp==NULL) {
      printf("ERROR: error changing type\n");
      exit(1);
    }
    MRIfree(&OutVol);
    OutVol = OutVoltmp;
  }

  printf("Saving to %s\n",outvolid);
  MRIwriteType(OutVol,outvolid,outtype);

  MRIfree(&OutVol);

  printf("Done\n");

  return(0);
  exit(0);
}
/* --------------------------------------------- */
/* ---------End Main --------------------------- */
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

    else if (!strcasecmp(option, "--version")) print_version() ;

    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--nodebug")) debug = 0;

    else if (!strcmp(option, "--vol")) {
      if (nargc < 1) argnerr(option,1);
      volid[nvols] = pargv[0];
      nvols++;
      nargsused = 1;
    } else if (!strcmp(option, "--vol_type")) {
      if (nargc < 1) argnerr(option,1);
      voltypestring = pargv[0];
      voltype = string_to_type(voltypestring);
      nargsused = 1;
    } else if (!strcmp(option, "--out")) {
      if (nargc < 1) argnerr(option,1);
      outvolid = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--out_type")) {
      if (nargc < 1) argnerr(option,1);
      outtypestring = pargv[0];
      outtype = string_to_type(outtypestring);
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
  exit(-1) ;
}
/* --------------------------------------------- */
static void print_usage(void) {
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --vol  path to input volume, ... \n");
  printf("   --vol_type  input volume format \n");
  printf("   --out       path to output volume\n");
  printf("   --out_type  output format \n");
  printf("\n");
  printf(" Other Options\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf(
    "This program will average multiple volumes together (including 4D volumes).\n"
    "\n"
    "OPTIONS\n"
    "\n"
    "  --vol path to input volume (repeat for each input volume)\n"
    "  --vol_type format type of all input volumes \n"
    "\n"
    "  --out path to output volume \n"
    "  --out_type format type of out volume (default is that of input)\n"
    "\n"
    "  --version : print version and exit.\n"
    "\n"
    "SPECIFYING THE INPUT/OUTPUT PATH and TYPE\n"
    "\n"
    "mri_xvolavg accepts all input types as mri_convert (see mri_convert --help \n"
    "for more information).\n"
    "\n"
    "NOTES\n"
    "\n"
    "There is another program called mri_average which has similar functionality.\n"
    "mri_average can align the input volumes and force them to conform to the COR\n"
    "format. However, it cannot process volumes with multiple frames.\n"
    "\n"

  ) ;
  exit(-1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(-1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void) {

  if (nvols < 2) {
    printf("ERROR: number of input volumes = %d, must be >= 2\n",nvols);
    exit(1);
  }

  if (outvolid == NULL) {
    printf("ERROR: no output volume supplied\n");
    exit(-1);
  }

  if (voltype == MRI_VOLUME_TYPE_UNKNOWN) {
    if (defaulttype == MRI_VOLUME_TYPE_UNKNOWN)
      voltype = mri_identify(volid[0]);
    else
      voltype = defaulttype;
  }
  if (voltype == MRI_VOLUME_TYPE_UNKNOWN) {
    fprintf(stderr,"ERROR: could not determine type of %s\n",volid[0]);
    exit(-1);
  }

  if (outtype == MRI_VOLUME_TYPE_UNKNOWN) {
    if (defaulttype == MRI_VOLUME_TYPE_UNKNOWN)
      outtype = voltype;
    else
      outtype = defaulttype;
  }

  return;
}

/*---------------------------------------------------------------*/
static int singledash(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}

