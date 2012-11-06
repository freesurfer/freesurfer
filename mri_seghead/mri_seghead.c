/**
 * @file  mri_seghead.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: mreuter $
 *    $Date: 2012/11/06 23:11:52 $
 *    $Revision: 1.6 $
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


/*
  Name:    mri_seghead.c
  Author:  Douglas N. Greve
  email:   analysis-bugs@nmr.mgh.harvard.edu
  Date:    2/27/02
  Purpose: Segment the head.
  $Id: mri_seghead.c,v 1.6 2012/11/06 23:11:52 mreuter Exp $
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <unistd.h>
#include <string.h>

#include "error.h"
#include "diag.h"
#include "proto.h"

#include "fio.h"
#include "mri.h"
#include "MRIio_old.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  isflag(char *flag);
static int  singledash(char *flag);
static int  stringmatch(char *str1, char *str2);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_seghead.c,v 1.6 2012/11/06 23:11:52 mreuter Exp $";
char *Progname = NULL;
int debug = 0;
char *subjid;

char *involid;
char *outvolid;
unsigned char fillval = 255;
int thresh1  = -1;
int thresh2  = -1;
int nhitsmin = 2;

int fillslices = 1;
int fillrows = 1;
int fillcols = 1;

int dilate = 0;

MRI *invol;
MRI *HitMap;
MRI *outvol, *outvol2;
MRI *cvol, *rvol, *svol;
MRI *mrgvol;
MRI *invol_orig;

char tmpstr[1000];
char *hvoldat = NULL;
FILE *fp;

/*---------------------------------------------------------------*/
int main(int argc, char **argv) {
  int c,r,s,n;
  int smin, smax, rmin, rmax, cmin, cmax, nhits;

  /* This is to shut up the compiler */
  isflag(" ");

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();

  dump_options(stdout);

  /* Make sure output directory is writable */
  if (! fio_DirIsWritable(outvolid,1) &&
      fio_DirIsWritable(outvolid,0)) {
    printf("ERROR: could not access output %s\n",outvolid);
    exit(1);
  }

  printf("Loading input volume\n");
  fflush(stdout);
  invol_orig = MRIread(involid);
  if (invol_orig == NULL) {
    printf("ERROR: could not read %s\n",involid);
    exit(1);
  }

  invol = MRIchangeType(invol_orig,MRI_UCHAR,0,255,1);
  if (invol == NULL) {
    printf("ERROR: bvolumeWrite: MRIchangeType\n");
    return(1);
  }

  cvol = MRIclone(invol,NULL);
  rvol = MRIclone(invol,NULL);
  svol = MRIclone(invol,NULL);
  mrgvol = MRIclone(invol,NULL);
  outvol = MRIclone(invol,NULL);

  if (fillcols) {
    printf("Filling Columns\n");
    fflush(stdout);
    for (r=0; r < invol->height; r++) {
      for (s=0; s < invol->depth; s++) {

        c = 0;
        nhits = 0;
        while (c < invol->width && nhits < nhitsmin) {
          if (MRIseq_vox(invol,c,r,s,0) > thresh1) nhits ++;
          else nhits = 0;
          c++;
        }
        if (nhits < nhitsmin) continue;
        cmin = c;

        c = invol->width - 1;
        nhits = 0;
        while (c > cmin && nhits < nhitsmin) {
          if (MRIseq_vox(invol,c,r,s,0) > thresh1) nhits ++;
          else nhits = 0;
          c--;
        }
        cmax = c;

        if (cmin >= dilate) cmin-=dilate;
        else cmin = 0;
        if (cmax <= invol->width-1-dilate) cmax+=dilate;
        else cmax = invol->width-1; 

        for (c = cmin; c <= cmax; c++) MRIseq_vox(cvol,c,r,s,0) = fillval;
        //printf("%3d %3d  %3d %3d\n",r,s,cmin,cmax);
      }
    }
  }

  if (fillrows) {
    printf("Filling Rows\n");
    fflush(stdout);
    for (c=0; c < invol->width; c++) {
      for (s=0; s < invol->depth; s++) {

        r = 0;
        nhits = 0;
        while (r < invol->height && nhits < nhitsmin) {
          if (MRIseq_vox(invol,c,r,s,0) > thresh1) nhits ++;
          else nhits = 0;
          r++;
        }
        if (nhits < nhitsmin) continue;
        rmin = r;

        r = invol->height - 1;
        nhits = 0;
        while (r > rmin && nhits < nhitsmin) {
          if (MRIseq_vox(invol,c,r,s,0) > thresh1) nhits ++;
          else nhits = 0;
          r--;
        }
        rmax = r;

        if (rmin >= dilate) rmin-=dilate;
        else rmin = 0;
        if (rmax <= invol->height-1-dilate) rmax+=dilate;
        else rmax = invol->height-1; 
        
        for (r = rmin; r <= rmax; r++) MRIseq_vox(rvol,c,r,s,0) = fillval;
        //printf("%3d %3d  %3d %3d\n",r,s,rmin,rmax);
      }
    }
  }

  if (fillslices) {
    printf("Filling Slices\n");
    fflush(stdout);
    for (c=0; c < invol->width; c++) {
      for (r=0; r < invol->height; r++) {

        s = 0;
        nhits = 0;
        while (s < invol->depth && nhits < nhitsmin) {
          if (MRIseq_vox(invol,c,r,s,0) > thresh1) nhits ++;
          else nhits = 0;
          s++;
        }
        if (nhits < nhitsmin) continue;
        smin = s;

        s = invol->depth - 1;
        nhits = 0;
        while (s > smin && nhits < nhitsmin) {
          if (MRIseq_vox(invol,c,r,s,0) > thresh1) nhits ++;
          else nhits = 0;
          s--;
        }
        smax = s;

        if (smin >= dilate) smin-=dilate;
        else smin = 0;
        if (smax <= invol->depth-1-dilate) smax+=dilate;
        else smax = invol->depth-1; 
        
        for (s = smin; s <= smax; s++) MRIseq_vox(svol,c,r,s,0) = fillval;
        //printf("%3d %3d  %3d %3d\n",r,s,rmin,rmax);
      }
    }
  }

  printf("Merging and Inverting\n");
  for (c=0; c < invol->width; c++) {
    for (r=0; r < invol->height; r++) {
      for (s=0; s < invol->depth; s++) {
        if (MRIseq_vox(cvol,c,r,s,0) &&
            MRIseq_vox(rvol,c,r,s,0) &&
            MRIseq_vox(svol,c,r,s,0))
          MRIseq_vox(mrgvol,c,r,s,0) = 0;
        else
          MRIseq_vox(mrgvol,c,r,s,0) = 255-MRIseq_vox(invol,c,r,s,0);
      }
    }
  }
  
  printf("Growing\n");
  fflush(stdout);
  outvol = MRIfill(mrgvol, NULL, (int)(mrgvol->width/2.0),
                   (int)(mrgvol->height/2.0), (int)(mrgvol->depth/2.),
                   thresh2, fillval, -1);

  printf("Counting\n");
  n = 0;
  double backnoise = 0.0;
  int backcount = 0;
  double backmax = 0.0;
  double val = 0.0;
  for (c=0; c < invol->width; c++) {
    for (r=0; r < invol->height; r++) {
      for (s=0; s < invol->depth; s++) {
        if (MRIseq_vox(outvol,c,r,s,0)) n++;
        else
        {
          val = MRIgetVoxVal(invol_orig,c,r,s,0);
          backnoise += val;
          if (backmax < val) backmax = val;
          backcount++;
        }
      }
    }
  }
  backnoise= backnoise/backcount;
  printf("N Head Voxels = %d\n",n);
  printf("N Back Voxels = %d\n",backcount);
  printf("Avg. Back Intensity = %f\n",backnoise);
  printf("Max. Back Intensity = %f\n",backmax);
  
  if(hvoldat){
    fp = fopen(hvoldat,"w");
    fprintf(fp,"%lf\n",n*invol->xsize*invol->ysize*invol->zsize);
    fclose(fp);
  }

  if(outvolid){
    printf("Writing output\n");
    fflush(stdout);
    MRIwrite(outvol,outvolid);
  }

  printf("Done\n");

  return(0);
}
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  int a;

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
    else if (!strcasecmp(option, "--fillrows")) fillrows = 1;
    else if (!strcasecmp(option, "--fillcols")) fillcols = 1;
    else if (!strcasecmp(option, "--fillslices")) fillslices = 1;

    else if ( stringmatch(option, "--invol") ||
              stringmatch(option, "--i") ) {
      if (nargc < 1) argnerr(option,1);
      involid = pargv[0];
      nargsused = 1;
    } else if (stringmatch(option, "--outvol") ||
               stringmatch(option, "--o") ) {
      if (nargc < 1) argnerr(option,1);
      outvolid = pargv[0];
      nargsused = 1;
    } else if ( stringmatch(option, "--subject") ||
                stringmatch(option, "--s") ) {
      if (nargc < 1) argnerr(option,1);
      subjid = pargv[0];
      nargsused = 1;
    } else if (stringmatch(option, "--fill")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&a);
      fillval = (unsigned char) a;
      nargsused = 1;
    } else if (stringmatch(option, "--thresh1")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&thresh1);
      nargsused = 1;
    } else if (stringmatch(option, "--thresh2")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&thresh2);
      nargsused = 1;
    } else if (stringmatch(option, "--thresh")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&thresh1);
      thresh2 = thresh1;
      nargsused = 1;
    } else if (stringmatch(option, "--nhitsmin")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&nhitsmin);
      nargsused = 1;
    } else if (stringmatch(option, "--hvoldat")) {
      if (nargc < 1) argnerr(option,1);
      hvoldat = pargv[0];
      nargsused = 1;
    } else if (stringmatch(option, "--dilate")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&dilate);
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --invol     input volume id (eg, T1)\n");
  printf("   --outvol    input volume id\n");
  printf("   --fill      fill value  (255)\n");
  printf("   --thresh1   threshold value  (eg, 20)\n");
  printf("   --thresh2   threshold value  (eg, 20)\n");
  printf("   --thresh    single threshold value for 1 and 2 \n");
  printf("   --nhitsmin  min number of consecutive hits (2) \n");
  printf("   --hvoldat   file : write head volume (mm3) to an ascii file \n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf("\n%s\n\n",vcid);

  printf(
    "This program binarizes the input volume such that all the voxels in the \n"
    "head are set to 255 (or whatever is passed with --fill). The result is \n"
    "stored in the output volume passed by --outvol. \n"
    " \n"
    "The program first creates a binarized mass just a few mm below the skin.\n"
    "This mass is then grown out using a connected components algorithm so \n"
    "that most of the skin details are retained. \n"
    " \n"
    "The initial mass is created in the following way. First, for every row \n"
    "and slice, the column is searched from both ends for the 'skin'. The \n"
    "skin is defined as the first consecutive nhitsmin voxels over thresh1.\n"
    "Once the skin is found coming from both directions, everything in between\n"
    "is binarized to the fill value. This process is repeated for the rows \n"
    "and slices. The initial mass is created by ANDing all three voxels.\n"
    " \n"
    "After the initial mass is defined, the original volume is modified so\n"
    "that all the voxels in the mass are set to 255. This has the effect of \n"
    "filling in all of the subepidermal features that would normally be below \n"
    "threshold. A seed point is chosen at the center of the volume. The final \n"
    "binarized volume is computed as all the voxels above thresh2 connected to  \n"
    "the seed point. \n"
    " \n"
    "--thresh threshold will set thresh1 and thresh2 to threshold. \n"
    " \n"
    "Typical values for the parameters are: thresh=20 and nhitsmin=2. \n"
    " \n"
  );

  exit(1) ;
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
/* --------------------------------------------- */
static void check_options(void) {
  if (involid == NULL) {
    printf("An input volume path must be supplied\n");
    exit(1);
  }
/*  if(outvolid == NULL && hvoldat == NULL) {
    printf("An output volume path must be supplied\n");
    exit(1);
  }*/
  if (thresh1 <= 0) {
    printf("Must specify a thresh1 > 0\n");
    exit(1);
  }
  if (thresh2 <= 0) {
    printf("Must specify a thresh2 > 0\n");
    exit(1);
  }
  if (nhitsmin <= 0) {
    printf("Must specify a nhitsmin > 0\n");
    exit(1);
  }
  return;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"input volume:  %s\n",involid);
  fprintf(fp,"output volume: %s\n",outvolid);
  fprintf(fp,"threshold1:    %d\n",thresh1);
  fprintf(fp,"threshold2:    %d\n",thresh2);
  fprintf(fp,"nhitsmin:      %d\n",nhitsmin);
  fprintf(fp,"fill value:    %d\n",(int)fillval);
  if (dilate > 0) fprintf(fp,"dilate    :    %d\n",dilate);
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
/*---------------------------------------------------------------*/
static int isflag(char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);
  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*------------------------------------------------------------*/
static int stringmatch(char *str1, char *str2) {
  if (! strcmp(str1,str2)) return(1);
  return(0);
}

