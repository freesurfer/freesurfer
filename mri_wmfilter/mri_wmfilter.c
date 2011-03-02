/**
 * @file  mri_wmfilter.c
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:26 $
 *    $Revision: 1.18 $
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


#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "macros.h"
#include "const.h"
#include "MRIio_old.h"
/*#include "typedefs.h"*/  /* not needed without matrix stuff */
#include "matrix.h"
#include "mri.h"
#include "error.h"
#include "proto.h"
#include "version.h"

static char vcid[] = "$Id: mri_wmfilter.c,v 1.18 2011/03/02 00:04:26 nicks Exp $";

/*-------------------------------------------------------------------
                                CONSTANTS
-------------------------------------------------------------------*/

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#define MAXCOR 500
#define MAXLEN 100
#define DIR_FILE "ic1.tri"

/*-------------------------------------------------------------------
                                GLOBAL DATA
-------------------------------------------------------------------*/

int wx0=100,wy0=100;
float white_hilim = 125;  /* new smaller range goes with improved */
float white_lolim = 95;   /* intensity normalization */
float gray_hilim = 100;
float xmin,xmax;
float ymin,ymax;
float zmin,zmax;
float st,ps,fov,xx0,xx1,yy0,yy1,zz0,zz1;
float ctrx,ctry,ctrz;

#if 0
FLOATTYPE **mat,**mati,*vec,*vec2;
#endif
float xcor[MAXCOR],ycor[MAXCOR],zcor[MAXCOR];
int nver,ncor;

char *Progname ;

static int central_plane = 0 ;
static int grayscale_plane = 0 ;
static char *output_name = "wm" ;
static float slope = 0.00f ;  /* 0 mimic original behavior */
static float lslope = 0.0f ;  /* slope for planar laplacian threshold mod */

/*-------------------------------------------------------------------
                             STATIC PROTOTYPES
-------------------------------------------------------------------*/

int main(int argc,char *argv[]) ;
static MRI *plane_filter(MRI *mri_src, MRI *mri_dst, int niter) ;
static int get_option(int argc, char *argv[]) ;
static void print_version(void) ;
static void print_help(void) ;

/*-------------------------------------------------------------------
                                FUNCTIONS
-------------------------------------------------------------------*/

int
main(int argc,char *argv[]) {
  int   i,j,option=1, nargs;
  float x,y,z;
  FILE  *fptr;
  char  fname[STRLEN],mfname[STRLEN],pfname[STRLEN],dfname[STRLEN];
  char  fpref[STRLEN],pname[STRLEN];
  char  *data_dir,*mri_dir;
  MRI   *mri_src, *mri_dst ;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, "$Id: mri_wmfilter.c,v 1.18 2011/03/02 00:04:26 nicks Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  for ( ; argc > 1 && (*argv[1] == '-') ; argc--, argv++) {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 2)
    print_help() ;   /* will exit */

#if 0
  mat = matrix(4,4);
  mati = matrix(4,4);
  vec = vector(4);
  vec2 = vector(4);
#endif

  if (argc<2) {
    exit(0);
  }

  sprintf(pname,"%s",argv[1]);
  if (argc>2)
    sscanf(argv[2],"%d",&option);

  data_dir = getenv("SUBJECTS_DIR");
  if (data_dir==NULL) {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(0);
  }
  mri_dir = getenv("FREESURFER_HOME");
  if (mri_dir==NULL) {
    printf("environment variable FREESURFER_HOME undefined (use setenv)\n");
    exit(0);
  }

  sprintf(fpref,"%s/%s",data_dir,pname);
  sprintf(mfname,"%s/mri/brain",fpref);
  sprintf(pfname,"%s/mri/%s", fpref,output_name);
  sprintf(dfname,"wmfilter.dat");
  fprintf(stderr, "reading input volume from %s...\n", mfname) ;
  mri_src = MRIread(mfname) ;
  if (!mri_src)
    ErrorExit(ERROR_NOFILE, "%s: could not read input volume %s\n",
              Progname, mfname) ;

  sprintf(fname,"%s",dfname);
  fptr = fopen(fname,"r");
  if (fptr) {
    fscanf(fptr,"%*s %f",&white_hilim);
    fscanf(fptr,"%*s %f",&white_lolim);
    fscanf(fptr,"%*s %f",&gray_hilim);
    fclose(fptr);
  } else
    printf("File %s not found - using defaults.\n",fname);
  printf("white_hilim = %2.1f, white_lolim = %2.1f, gray_hilim = %2.1f\n",
         white_hilim,white_lolim,gray_hilim);
  fflush(stdout) ;

  sprintf(fname,"%s/lib/bem/%s",mri_dir,DIR_FILE);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {
    printf("File %s not found.\n",fname);
    exit(0);
  }
  fscanf(fptr,"%d",&nver);
  for (i=0,ncor=0;i<nver;i++) {
    fscanf(fptr,"%*d %f %f %f",&x,&y,&z);
    for (j=0;j<ncor;j++)
      if ((x==-xcor[j]) && (y==-ycor[j]) && (z==-zcor[j])) goto L1;
    xcor[ncor] = x;
    ycor[ncor] = y;
    zcor[ncor] = z;
    ncor++;
L1:
    ;
  }
  printf("%d unique orientations\n",ncor);
  fflush(stdout) ;
  mri_dst = plane_filter(mri_src, NULL, option);
  fprintf(stderr, "writing output to %s...\n", pfname) ;
  MRIwrite(mri_dst, pfname) ;
  exit(0) ;
  return(0) ;
}
#define DEFAULT_FRAC  0.6

static MRI *
plane_filter(MRI *mri_src, MRI *mri_dst, int niter) {
  int    i,j,k,di,dj,dk,m,n,u,ws2=2,maxi,mini, wsize, width, height, depth ;
  float  numvox,numnz,numz, probably_white ;
  float  f,f2,a,b,c,s;
  double sum2,sum,var,avg,tvar,maxvar,minvar;
  double sum2v[MAXLEN],sumv[MAXLEN],avgv[MAXLEN],varv[MAXLEN],nv[MAXLEN],
  tlaplacian, laplacian ;
  float  cfrac ;
  MRI    *mri_tmp ;

  width = mri_src->width ;
  height = mri_src->height ;
  depth = mri_src->depth ;
  if (!mri_dst)
    mri_dst = MRIclone(mri_src, NULL) ;

  mri_tmp = MRIcopy(mri_src, NULL) ;

  cfrac = DEFAULT_FRAC ;
  probably_white = (white_lolim + gray_hilim) / 2.0 ;
  /*
    printf("plane_filter(%d)\n",niter);
  */
  wsize = 2*ws2+1 ;

  /*
    eliminate all voxels that are lower than white_lolim or higher
    than white_hilim, and place the results in im[][][] and fill[]
    (im2 is  still original image).
  */
  for (k=0;k<depth-1;k++)
    for (i=0;i<height-1;i++)
      for (j=0;j<width-1;j++) {
        MRIvox(mri_tmp, j, i, k) = MRIvox(mri_src,j,i,k);
        if (MRIvox(mri_src,j,i,k)>white_hilim ||
            MRIvox(mri_src,j,i,k)< probably_white)
          MRIvox(mri_src,j,i,k) = 0;
        MRIvox(mri_dst,j,i,k) = MRIvox(mri_src,j,i,k);
      }

  for (k=ws2;k<depth-1-ws2;k++) {
    if (!((k+1)%10))
      printf("processing slice %d\n",k+1);
    fflush(stdout) ;
    for (i=ws2;i<height-1-ws2;i++)
      for (j=ws2;j<width-1-ws2;j++) {
        if (k == 80 && i == 140 && j == 145) {
          int kk = 0 ;
          kk++ ;
        }

        /* if the voxel is not in the ambiguous intensity range,
           just set the output and continue
        */
        if ((MRIvox(mri_tmp, j, i, k) < white_lolim) ||
            (MRIvox(mri_tmp, j, i, k) > white_hilim)) {
          MRIvox(mri_dst,j,i,k) = 0;
          continue ;
        }
        if (MRIvox(mri_tmp, j, i, k) > gray_hilim) {
          MRIvox(mri_dst,j,i,k) = MRIvox(mri_tmp, j, i, k) ;
          continue ;
        }

        /* count number of voxels on and off in 5x5x5 cube */
        numvox = numnz = numz = 0;
        for (dk = -ws2;dk<=ws2;dk++)
          for (di = -ws2;di<=ws2;di++)
            for (dj = -ws2;dj<=ws2;dj++) {
              f = MRIvox(mri_src,j+dj,i+di,k+dk);
              s = dk*c+di*b+dj*a;
              numvox++;
              if (f!=0) {
                numnz++;
              } else
                numz++;
            }
        /* if voxel was off but majority in volume are on, or
           if voxel was on  but majority in volume are off, take
           a look at the plane of least variance for classification
        */
        if (
          ((MRIvox(mri_src,j,i,k)==0 &&
            numnz>=(DEFAULT_FRAC)*(wsize*wsize)) &&
           (MRIvox(mri_tmp, j, i, k) >= white_lolim))
          ||
          (MRIvox(mri_src,j,i,k)!=0 &&
           (numz>=(DEFAULT_FRAC)*(wsize*wsize)) &&
           (MRIvox(mri_tmp, j, i, k) <= gray_hilim))) {
          tlaplacian = laplacian = 0.0;
          maxvar = -1000000;
          minvar = 1000000;
          maxi = mini = -1;
          for (m=0;m<ncor;m++)    /* for each orientation */
          {
            /* (a,b,c) is normal (orientation) vector */
            a = xcor[m];
            b = ycor[m];
            c = zcor[m];
            sum = sum2 = n = 0;
            for (u=0;u<wsize;u++)
              sumv[u] = sum2v[u] = nv[u] = 0;
            for (dk = -ws2;dk<=ws2;dk++)
              for (di = -ws2;di<=ws2;di++)
                for (dj = -ws2;dj<=ws2;dj++) {
                  u = ws2+floor(dk*c+di*b+dj*a+0.5);
                  u = (u<0) ? 0 : (u>=wsize) ? wsize-1 : u ;
                  if (!grayscale_plane)
                    f = MRIvox(mri_src,j+dj,i+di,k+dk);   /* segmented image */
                  else
                    f = MRIvox(mri_dst,j+dj,i+di,k+dk);  /* gray-scale image */
                  sum2v[u] += f*f ;
                  sumv[u] += f;
                  nv[u] += 1;
                  n += 1;
                  sum2 += f*f;
                  sum += f;

                }
            avg = sum/n;
            var = sum2/n-avg*avg;  /* total mean and variance */
            tvar = 0;
            if (central_plane)  /* only consider variance in central plane */
            {
              u = ws2 ;
              avgv[u] = sumv[u]/nv[u];
              varv[u] = sum2v[u]/nv[u]-avgv[u]*avgv[u];
              tvar = varv[u];
              tlaplacian = 0.0f ; /* only doing central plane - no laplacian */
            } else  /* variance in all planes in stack */
            {
              for (u=0;u<wsize;u++) {
                avgv[u] = sumv[u]/nv[u];
                varv[u] = sum2v[u]/nv[u]-avgv[u]*avgv[u];
                tvar += varv[u];
                tvar /= (wsize);
              }
              /* compute planar 'laplacian' */
              for (tlaplacian = 0.0, u=0;u<wsize;u++) {
                if (u == ws2)
                  continue ;
                if (avgv[u] > avgv[ws2])
                  tlaplacian += 1.0f ;  /* darker than surround - concave */
                else
                  tlaplacian -= 1.0f ;  /* brighter than surround - convex */
              }
            }
            if (tvar>maxvar) {
              maxvar=tvar;
              maxi=m;
            }
            if (tvar<minvar) {
              minvar=tvar;
              mini=m;
              laplacian = tlaplacian ;
            }
          }

          /* (a,b,c) now orientation vector for plane of least variance */
          a = xcor[mini];
          b = ycor[mini];
          c = zcor[mini];
          /*
                  printf(" -> %d,%d: (%f,%f) %f, %f, %f\n",mini,maxi,minvar,maxvar,a,b,c);
          */
          /*
            for (m=0;m<3;m++)
            for (n=0;n<3;n++)
            mat[m][n] = 0;
            for (m=0;m<3;m++)
            vec[m] = 0;
            for (dk = -ws2;dk<=ws2;dk++)
            for (di = -ws2;di<=ws2;di++)
            for (dj = -ws2;dj<=ws2;dj++)
            {
            f = MRIvox(mri_dst,j+dj,i+di,k+dk);
            s = dk*c+di*b+dj*a;
            s2 = s*s;
            mat[0][0] += s2*s2;
            mat[0][1] = (mat[1][0] += s2*s);
            mat[0][2] = (mat[2][0] += s2);
            mat[1][1] += s*s;
            mat[1][2] = (mat[2][1] += s);
            mat[2][2] += 1;
             vec[0] += f*s2;
             vec[1] += f*s;
             vec[2] += f;
             }
             inverse(mat,mati,3);
             vector_multiply(mati,vec,vec2,3,3);
             f = MRIvox(mri_src,j,i,k);
             if (vec2[0]>2.0 && f<=gray_hilim) f = 0;
          f = (f>255)?255:(f<0)?0:f;
          MRIvox(mri_dst,j,i,k) = f;
          */
          /* count # of on voxels and # of off voxels in plane */
          numvox = numnz = numz = sum = sum2 = 0;
          for (dk = -ws2;dk<=ws2;dk++)
            for (di = -ws2;di<=ws2;di++)
              for (dj = -ws2;dj<=ws2;dj++) {
                f = MRIvox(mri_src,j+dj,i+di,k+dk);
                f2 = MRIvox(mri_dst,j+dj,i+di,k+dk);
                s = dk*c+di*b+dj*a;
                if (fabs(s)<=0.5)  /* in central plane */
                {
                  numvox++;
                  sum2 += f2;
                  if (f!=0) {
                    numnz++;
                    sum += f;
                  } else
                    numz++;
                }
              }
          if (numnz!=0) sum /= numnz;
          if (numvox!=0) sum2 /= numvox;
          f = MRIvox(mri_src,j,i,k);
          f2 = MRIvox(mri_tmp, j, i, k);
          /*
            printf("%d %d %d %d %f\n",
            k,i,j,(int)f,(int)numvox,(int)numnz,(int)numz,numz/numvox);
          */
          /*
             a positive laplacian means that the central plane is darker
             than the surrounding planes in the same orientation. This should
             bias the voxel to being classified as nonwhite. Conversely, a
             negative laplacian means that the central plane is brighter than
             the surrounding planes and that the voxel should be more likely
             to be white.
          */
          cfrac = DEFAULT_FRAC ;
          if (f != 0)  /* compute fraction needed to change it to nonwhite */
            cfrac += (float)(f2 - probably_white)*slope - laplacian*lslope ;
          else
            cfrac += (float)(probably_white - f2)*slope + laplacian*lslope ;
          if (f!=0 && numz/numvox>cfrac)
            f=0 ;   /* change preliminary classification white --> nonwhite */
          else if (f==0 && numnz/numvox>cfrac)
            f=f2 ;  /* change preliminary classification nonwhite --> white */
          MRIvox(mri_dst,j,i,k) = f ;
        }
      }
  }
  for (k=0;k<depth-1;k++)
    for (i=0;i<height-1;i++)
      for (j=0;j<width-1;j++) {
        MRIvox(mri_src,j,i,k) = MRIvox(mri_dst,j,i,k);
        if (MRIvox(mri_src,j,i,k)>white_hilim ||
            MRIvox(mri_src,j,i,k)<white_lolim)
          MRIvox(mri_src,j,i,k) = 0;
      }

  MRIfree(&mri_tmp) ;
  return(mri_dst) ;
}

/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
print_version(void) {
  fprintf(stderr, "%s version %s\n", vcid, Progname) ;
  exit(1) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
print_help(void) {
  fprintf(stderr, "Usage: %s subject_name [option ...]\n\n", Progname);
  fprintf(stderr,
          "this program will read an input volume from\n"
          "$SUBJECTS_DIR/<subject_name>/brain and set all voxels which\n"
          "are determined to be other than white matter to 0, writing\n"
          "the output volume to $SUBJECTS_DIR/<subject_name>/wm.\n") ;

  exit(1) ;
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
  if (!strcmp(option, "-help"))
    print_help() ;
  else if (!strcmp(option, "-version"))
    print_version() ;
  else if (!strcmp(option, "grayscale")) {
    grayscale_plane = 1 ;
    fprintf(stderr, "using grayscale data to compute polv\n") ;
  } else if (!strcmp(option, "slope")) {
    slope = atof(argv[2])/100.0f ;
    fprintf(stderr, "modifying voting fraction by %2.1f%%/intensity\n",
            slope*100.0f) ;
    nargs = 1 ;
  } else if (!strcmp(option, "lslope")) {
    lslope = atof(argv[2])/100.0f ;
    fprintf(stderr,"modifying voting fraction by %2.1f%% * binary laplacian\n",
            lslope*100.0f) ;
    nargs = 1 ;
  } else if (!strcmp(option, "central")) {
    central_plane = 1 ;
    fprintf(stderr, "only using planes through origin\n") ;
  } else if (!strcmp(option, "wlo")) {
    white_lolim = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using white lolim %2.1f\n", white_lolim) ;
  } else if (!strcmp(option, "ghi")) {
    gray_hilim = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using gray hilim %2.1f\n", gray_hilim) ;
  } else if (!strcmp(option, "whi")) {
    white_hilim = atoi(argv[2]) ;
    nargs = 1 ;
    fprintf(stderr, "using white hilim %2.1f\n", white_hilim) ;
  } else switch (toupper(*option)) {
    case 'O':
      output_name = argv[2] ;
      fprintf(stderr, "outputting results to %s\n", output_name) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'U':
      print_help() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}
