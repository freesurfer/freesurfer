/**
 * @file  mri_cor2label.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:06 $
 *    $Revision: 1.7 $
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


/*----------------------------------------------------------
  Name: mri_cor2label.c
  $Id: mri_cor2label.c,v 1.7 2006/12/29 02:09:06 nicks Exp $
  Author: Douglas Greve
  Purpose: Converts values in a COR file to a label.
  -----------------------------------------------------------*/
#include <stdio.h>
#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "label.h"
#include "version.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);

static char vcid[] = "$Id: mri_cor2label.c,v 1.7 2006/12/29 02:09:06 nicks Exp $";
char *Progname ;
int main(int argc, char *argv[]) ;

static char *cordir;
static char *labelfile;
static char *volfile;
static int  labelid;
MRI *COR;
LABEL *lb;

int xi,yi,zi, c, nlabel;
float x,y,z;
char *subject_name = NULL;
int doit;
int synthlabel = 0;
int verbose = 0;
float xsum, ysum, zsum;


/*----------------------------------------------------*/
int main(int argc, char **argv) {
  FILE *fp;
  int nargs;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_cor2label.c,v 1.7 2006/12/29 02:09:06 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  cordir    = NULL;
  labelfile = NULL;
  labelid   = 257;

  parse_commandline(argc, argv);
  check_options();

  /* allocate a label as big as the COR itself */
  lb = LabelAlloc(256*256*256,subject_name,labelfile);

  fprintf(stderr,"Loading COR %s\n",cordir);
  COR = MRIread(cordir);

  fprintf(stderr,"Scanning the volume\n");
  nlabel = 0;
  xsum = 0.0;
  ysum = 0.0;
  zsum = 0.0;
  for (xi=0; xi < COR->width; xi++) {
    for (zi=0; zi < COR->height; zi++) {
      for (yi=0; yi < COR->depth; yi++) {
        c = (int) MRIvox(COR,xi,zi,yi);
        /* The call to MRIvox is with arguments (COR,xi,zi,yi) instead
           of (COR,xi,yi,zi) becase
           MRIvox(mri,x,y,z) = (((BUFTYPE *)mri->slices[z][y])[x])
           but from mriio.c, L959ish,      mri->slices[y][z])[x])*/

        doit = 0;
        if (synthlabel &&
            xi > 120 && xi < 130 &&
            zi > 128 && zi < 138 &&
            yi >  40 && yi <  50) doit = 1;

        if (c == labelid && !synthlabel) doit = 1;

        if (doit) {
          x = -xi + 128.0;
          y =  yi - 128.0;
          z = -zi + 128.0;
          if (verbose)
            printf("%5d   %3d %3d %3d   %6.2f %6.2f %6.2f \n",
                   nlabel,xi,yi,zi,x,y,z);
          lb->lv[nlabel].x = x;
          lb->lv[nlabel].y = y;
          lb->lv[nlabel].z = z;
          nlabel ++;
          xsum += x;
          ysum += y;
          zsum += z;
        }
      }
    }
  }
  lb->n_points = nlabel;
  printf("Found %d label voxels\n",nlabel);

  if (volfile != NULL) {
    printf("Writing volume to  %s\n",volfile);
    fp = fopen(volfile,"w");
    if (fp == NULL) {
      printf("ERROR: could not open  %s\n",volfile);
      exit(1);
    }
    fprintf(fp,"%d\n",nlabel);
    fclose(fp);
  }

  if (labelfile == NULL) exit(0);


  if (nlabel == 0) {
    printf("ERROR: found no voxels matching id %d \n",labelid);
    exit(1);
  }
  printf("Writing label file %s\n",labelfile);
  LabelWrite(lb,labelfile);

  fprintf(stderr,"Centroid: %6.2f  %6.2f  %6.2f \n",
          xsum/nlabel,ysum/nlabel,zsum/nlabel);

  fprintf(stderr,"mri_cor2label completed SUCCESSFULLY\n");

  return(0);
  exit(0);

}
/* --------------------------------------------- */
/* --------------------------------------------- */
/* --------------------------------------------- */

static int parse_commandline(int argc, char **argv) {
  int  nargc = argc, nargs = 0 ;
  char **pargv = argv, *option ;

  while (nargc > 0) {

    option = pargv[0];
    if (!stricmp(option, "--help")) print_help() ;
    else if (!stricmp(option, "--version")) print_version() ;

    /* ---- COR directory ------------ */
    else if (!strcmp(option, "--c")) {
      if (nargc < 2) argnerr(option,1);
      cordir = pargv[1];
      nargs = 2;
    }

    /* ---- label file ---------- */
    else if (!strcmp(option, "--l")) {
      if (nargc < 2) argnerr(option,1);
      labelfile = pargv[1];
      nargs = 2;
    }

    /* ---- count file ---------- */
    else if (!strcmp(option, "--v")) {
      if (nargc < 2) argnerr(option,1);
      volfile = pargv[1];
      nargs = 2;
    }

    /* ---- label id ---------- */
    else if (!strcmp(option, "--id")) {
      if (nargc < 2) argnerr(option,1);
      sscanf(pargv[1],"%d",&labelid);
      if (labelid < 0 || labelid > 255) {
        fprintf(stderr,"ERROR: id=%d, must be 0<id<255\n",labelid);
        exit(1);
      }
      nargs = 2;
    }

    /* ---- synthesize the label ---------- */
    else if (!strcmp(option, "--synthlabel")) {
      synthlabel = 1;
      nargs = 1;
    }

    /* ---- verbose ---------- */
    else if (!strcmp(option, "--verbose")) {
      verbose = 1;
      nargs = 1;
    } else {
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      exit(-1);
    }
    nargc -= nargs;
    pargv += nargs;
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
  fprintf(stdout, "\n");
  fprintf(stdout, "USAGE: %s \n",Progname);
  fprintf(stdout,"   --c  cordir  directory to find COR volume \n");
  fprintf(stdout,"   --id labelid 0-255 value in COR volume\n");
  fprintf(stdout,"   --l  labelfile name of output file\n");
  fprintf(stdout,"   --v  volfile : write label volume in file\n");
  fprintf(stdout,"   --help print out help information\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf(
    "\n"
    "Converts values in a COR volume to a label. The program searches the\n"
    "COR volume in directory 'cordir' for values equal to 'labelid'. The\n"
    "xyz values for each point are then computed assuming 1 mm^3 voxels and\n"
    "that xyz=0 at the center of the volume. The xyz values are then stored\n"
    "in 'labelfile' in the label file format; the vertex values are set to\n"
    "zero as is the statistic value.  While this program can be used with\n"
    "any COR volume, it was designed to convert parcellation volumes, which\n"
    "happen to be stored in COR format.  See tkmedit for more information\n"
    "on parcellations. The labelid must be within the range of 0 to 255.\n"
    "The label volume in mm^3 can be written to the argument of --v.\n"
    "\n"
    "Bugs:\n"
    "\n"
    "  If the name of the label does not include a forward slash (ie, '/')\n"
    "  then the program will attempt to put the label files in\n"
    "  $SUBJECTS_DIR/subject/label.  So, if you want the labels to go into\n"
    "  the current directory, make sure to put a './' in front of the label.\n"
    "\n"
    "Example:\n"
    "\n"
    "mri_cor2label --c /space/final/frontier/myparcellation\n"
    "              --id 57 --l /home/brainmapper/spattemp/57.label\n"
    "\n"
    "This will load the COR volume found in the directory\n"
    "/space/final/frontier/myparcellation and then search the volume for\n"
    "values equaling 57.  The results are stored in 57.label in the directory \n"
    "/home/brainmapper/spattemp, which must exist prior to execution.\n"
    "\n"
    "Hidden options:\n"
    "  --synth synthesizes a label (ignores values in COR volume)\n"
    "\n"
    "See also: tkmedit, mri_label2label\n\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  fprintf(stderr, "%s\n", vcid) ;
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
/* --------------------------------------------- */
static void check_options(void) {
  if (cordir == NULL) {
    fprintf(stderr,"ERROR: must supply a COR directory\n");
    exit(1);
  }
  if (labelfile == NULL && volfile == NULL) {
    fprintf(stderr,"ERROR: must be supply a label or volume file\n");
    exit(1);
  }
  if (labelid == 257) {
    fprintf(stderr,"ERROR: must supply a label id\n");
    exit(1);
  }

  return;
}
