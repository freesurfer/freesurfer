#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>

#include "macros.h"
#include "MRIio.h"
/*#include "typedefs.h"*/  /* not needed without matrix stuff */
#include "matrix.h"
#include "proto.h"

static char vcid[] = "$Id: mri_wmfilter.c,v 1.8 1999/07/26 16:55:56 fischl Exp $";

/*-------------------------------------------------------------------
                                CONSTANTS
-------------------------------------------------------------------*/

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#define IMGSIZE 256
#define MAXIM 512
#define MAXCOR 500
#define MAXLEN 100
#define DIR_FILE "ic1.tri"

/*-------------------------------------------------------------------
                                GLOBAL DATA
-------------------------------------------------------------------*/

int xnum=256,ynum=256;
unsigned long bufsize;
unsigned char **im[MAXIM];  /* image matrix  */
unsigned char **im2[MAXIM];
unsigned char **fill[MAXIM]; 
unsigned char *buf;  /* scratch memory  */
int imnr0,imnr1,numimg;
int wx0=100,wy0=100;
float white_hilim = 130;
float white_lolim = 90;
float gray_hilim = 100;
float threshold = 30;
float xmin,xmax;
float ymin,ymax;
float zmin,zmax;
float st,ps,fov,xx0,xx1,yy0,yy1,zz0,zz1;
float ctrx,ctry,ctrz;
float dfrac;
#if 0
FLOATTYPE **mat,**mati,*vec,*vec2;
#endif
float xcor[MAXCOR],ycor[MAXCOR],zcor[MAXCOR];
int nver,ncor;

char *Progname ;

/*-------------------------------------------------------------------
                             STATIC PROTOTYPES
-------------------------------------------------------------------*/

int main(int argc,char *argv[]) ;
static void read_image_info(char *fpref) ;
static void plane_filter(int niter) ;
static void read_images(char *fpref) ;
static void write_images(char *fpref) ;
static int get_option(int argc, char *argv[]) ;
static void print_version(void) ;
static void print_help(void) ;

/*-------------------------------------------------------------------
                                FUNCTIONS
-------------------------------------------------------------------*/

int
main(int argc,char *argv[])
{
  int   i,j,option=1, nargs;
  float x,y,z;
  FILE  *fptr;
  char  fname[200],mfname[200],pfname[200],dfname[200];
  char  fpref[200],pname[200];
  char  *data_dir,*mri_dir;

  Progname = argv[0] ;
  for ( ; argc > 1 && (*argv[1] == '-') ; argc--, argv++)
  {
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
  
  if (argc<2)
  {
    exit(0);
  }
  
  sprintf(pname,"%s",argv[1]);
  if (argc>2)
    sscanf(argv[2],"%d",&option);
  
  data_dir = getenv("SUBJECTS_DIR");
  if (data_dir==NULL)
  {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(0);
  }
  mri_dir = getenv("MRI_DIR");
  if (mri_dir==NULL)
  {
    printf("environment variable MRI_DIR undefined (use setenv)\n");
    exit(0);
  }
  
  sprintf(fpref,"%s/%s",data_dir,pname);
  sprintf(mfname,"%s/mri/brain/COR-",fpref);
  sprintf(pfname,"%s/mri/wm/COR-",fpref);
  sprintf(dfname,"wmfilter.dat");
  read_image_info(mfname);
  read_images(mfname);
  
  sprintf(fname,"%s",dfname);
  fptr = fopen(fname,"r");
  if (fptr) 
  {
    fscanf(fptr,"%*s %f",&white_hilim);
    fscanf(fptr,"%*s %f",&white_lolim);
    fscanf(fptr,"%*s %f",&gray_hilim);
    fclose(fptr);
  }
  else
    printf("File %s not found - using defaults.\n",fname);
  printf("white_hilim = %f, white_lolim = %f, gray_hilim = %f\n",
         white_hilim,white_lolim,gray_hilim);
  fflush(stdout) ;

  sprintf(fname,"%s/lib/bem/%s",mri_dir,DIR_FILE);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {printf("File %s not found.\n",fname);exit(0);}
  fscanf(fptr,"%d",&nver);
  for (i=0,ncor=0;i<nver;i++)
  {
    fscanf(fptr,"%*d %f %f %f",&x,&y,&z);
    for (j=0;j<ncor;j++)
      if ((x==-xcor[j]) && (y==-ycor[j]) && (z==-zcor[j])) goto L1;
    xcor[ncor] = x;
    ycor[ncor] = y;
    zcor[ncor] = z;
    ncor++;
  L1:;
  }
  printf("%d unique orientations\n",ncor);
  fflush(stdout) ;
  plane_filter(option);
  write_images(pfname);
  exit(0) ;
  return(0) ;
}

static void
read_image_info(char *fpref)
{
  FILE *fptr;
  char fname[100];

  sprintf(fname,"%s.info",fpref);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {printf("File %s not found.\n",fname);exit(1);}
  fscanf(fptr,"%*s %d",&imnr0);
  fscanf(fptr,"%*s %d",&imnr1);
  fscanf(fptr,"%*s %*d");
  fscanf(fptr,"%*s %d",&xnum);
  fscanf(fptr,"%*s %d",&ynum);
  fscanf(fptr,"%*s %f",&fov);
  fscanf(fptr,"%*s %f",&ps);
  fscanf(fptr,"%*s %f",&st);
  fscanf(fptr,"%*s %*f"); /* locatn */
  fscanf(fptr,"%*s %f",&xx0); /* strtx */
  fscanf(fptr,"%*s %f",&xx1); /* endx */
  fscanf(fptr,"%*s %f",&yy0); /* strty */
  fscanf(fptr,"%*s %f",&yy1); /* endy */
  fscanf(fptr,"%*s %f",&zz0); /* strtz */
  fscanf(fptr,"%*s %f",&zz1); /* endz */
  fov *= 1000;
  ps *= 1000;
  st *= 1000;
  xx0 *= 1000;
  xx1 *= 1000;
  yy0 *= 1000;
  yy1 *= 1000;
  zz0 *= 1000;
  zz1 *= 1000;
  fclose(fptr);
  numimg = imnr1-imnr0+1;
  ctrx = (xx0+xx1)/2.0;
  ctry = (yy0+yy1)/2.0;
  ctrz = (zz0+zz1)/2.0;
}

static void
plane_filter(int niter)
{
  int    i,j,k,di,dj,dk,m,n,u,ws2=2,maxi,mini, wsize, probably_white ;
  float  numvox,numnz,numz;
  float  f,f2,a,b,c,s;
  double sum2,sum,var,avg,tvar,maxvar,minvar;
  double sum2v[MAXLEN],sumv[MAXLEN],avgv[MAXLEN],varv[MAXLEN],nv[MAXLEN];
  float  cfracz = 0.6;
  float  cfracnz = 0.6;
  
  probably_white = (white_lolim + gray_hilim) / 2 ;
/*
  printf("plane_filter(%d)\n",niter);
*/
  wsize = 2*ws2+1 ;
  for (k=0;k<numimg;k++)
  {
    fill[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    im2[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    for (i=0;i<IMGSIZE;i++)
    {
      fill[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
      im2[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
    }
  }

  /*
    eliminate all voxels that are lower than white_lolim or higher
    than white_hilim, and place the results in im[][][] and fill[]
    (im2 is  still original image).
  */
  for (k=0;k<numimg-1;k++)
  for (i=0;i<IMGSIZE-1;i++)
  for (j=0;j<IMGSIZE-1;j++)
  {
    im2[k][i][j] = im[k][i][j];
    if (im[k][i][j]>white_hilim || im[k][i][j]< probably_white)
      im[k][i][j] = 0;
    fill[k][i][j] = im[k][i][j];
  }

  for (k=ws2;k<numimg-1-ws2;k++)
/*
  if (k==106)
*/
  {
    if (k == 127)
    {
      int kk = 0 ;
      kk++ ;
    }
    printf("processing slice %d\n",k+1);
    fflush(stdout) ;
    for (i=ws2;i<IMGSIZE-1-ws2;i++)
    for (j=ws2;j<IMGSIZE-1-ws2;j++)
    {
      if (k == 127 && i == 46 && j == 114)
      {
        int kk = 0 ;
        kk++ ;
      }

      /* if the voxel is not in the ambiguous intensity range,
         just set the output and continue
      */
      if ((im2[k][i][j] < white_lolim) || (im2[k][i][j] > white_hilim))
      {
        fill[k][i][j] = 0;
        continue ;
      }
      if (im2[k][i][j] > gray_hilim)
      {
        fill[k][i][j] = im2[k][i][j] ;
        continue ;
      }

      /* count number of voxels on and off in 5x5x5 cube */
      numvox = numnz = numz = 0;
      for (dk = -ws2;dk<=ws2;dk++)
      for (di = -ws2;di<=ws2;di++)
      for (dj = -ws2;dj<=ws2;dj++) 
      {
         f = im[k+dk][i+di][j+dj];
         s = dk*c+di*b+dj*a;
         numvox++;
         if (f!=0) 
         {
           numnz++;
         } 
         else 
           numz++;
      }
      /* if voxel was off but majority in volume are on, or
         if voxel was on  but majority in volume are off, take
         a look at the plane of least variance for classification
      */
      if (
          ((im[k][i][j]==0 && 
            numnz>=(cfracnz)*(wsize*wsize)) &&
           (im2[k][i][j] >= white_lolim))
          || 
          (im[k][i][j]!=0 && 
           (numz>=(cfracz)*(wsize*wsize)) &&
           (im2[k][i][j] <= gray_hilim)))
      {
        maxvar = -1000000;
        minvar = 1000000;
        maxi = mini = -1;
        for (m=0;m<ncor;m++)
        {
          a = xcor[m];
          b = ycor[m];
          c = zcor[m];
          sum = sum2 = n = 0;
          for (u=0;u<wsize;u++)
            sumv[u] = sum2v[u] = nv[u] = 0;
          for (dk = -ws2;dk<=ws2;dk++)
            for (di = -ws2;di<=ws2;di++)
              for (dj = -ws2;dj<=ws2;dj++)
              {
                u = ws2+floor(dk*c+di*b+dj*a+0.5);
                u = (u<0)?0:(u>2*ws2+1-1)?2*ws2+1-1:u;
                
                f = im[k+dk][i+di][j+dj];
                
/*
  f = im2[k+dk][i+di][j+dj];
*/
                sum2v[u] += f*f;
                sumv[u] += f;
                nv[u] += 1;
                sum2 += f*f;
                sum += f;
                n += 1;
              }
          avg = sum/n;
          var = sum2/n-avg*avg;
          tvar = 0;
          for (u=0;u<2*ws2+1;u++)
          {
            avgv[u] = sumv[u]/nv[u];
            varv[u] = sum2v[u]/nv[u]-avgv[u]*avgv[u];
            tvar += varv[u];
          }
          tvar /= (wsize);
          if (tvar>maxvar) {maxvar=tvar;maxi=m;}
          if (tvar<minvar) {minvar=tvar;mini=m;}
        }
        
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
        f = im2[k+dk][i+di][j+dj];
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
         f = im[k][i][j];
         if (vec2[0]>2.0 && f<=gray_hilim) f = 0;
      f = (f>255)?255:(f<0)?0:f;
      fill[k][i][j] = f;
*/
        numvox = numnz = numz = sum = sum2 = 0;
        for (dk = -ws2;dk<=ws2;dk++)
          for (di = -ws2;di<=ws2;di++)
            for (dj = -ws2;dj<=ws2;dj++) 
            {
              f = im[k+dk][i+di][j+dj];
              f2 = im2[k+dk][i+di][j+dj];
              s = dk*c+di*b+dj*a;
              if (fabs(s)<=0.5) 
              {
                numvox++;
                sum2 += f2;
                if (f!=0) 
                {
                  numnz++;
                  sum += f;
                } 
                else 
                  numz++;
              }
            }
        if (numnz!=0) sum /= numnz;
        if (numvox!=0) sum2 /= numvox;
        f = im[k][i][j];
        f2 = im2[k][i][j];
        /*
          printf("%d %d %d %d %f\n",
          k,i,j,(int)f,(int)numvox,(int)numnz,(int)numz,numz/numvox);
        */
        if (f!=0 && numz/numvox>cfracz && f<=gray_hilim) 
        {
/*
        printf(" *** %d,%d,%d: %d->0 %d %d %d\n",
               k,i,j,(int)f,(int)numvox,(int)numnz,(int)numz);
*/
          f=0;
        }
        else if (f==0 && numnz/numvox>cfracnz) 
/*
      else if (f==0 && 
               (numnz/numvox>cfracnz || (vec2[0]<-1.0 && f2>=white_lolim-5))) 
*/
        {
/*
        printf(" !!! %d,%d,%d: %d(%d)->%d %d %d %d %f\n",
               k,i,j,(int)f,(int)f2,(int)sum,
               (int)numvox,(int)numnz,(int)numz,vec2[0]);
*/
          f=sum;
        }
/*
      f = sum2;
*/
        fill[k][i][j] = f2 /*f*/;
      }
    }
  }
  for (k=0;k<numimg-1;k++)
    for (i=0;i<IMGSIZE-1;i++)
      for (j=0;j<IMGSIZE-1;j++)
      {
        im[k][i][j] = fill[k][i][j];
        if (im[k][i][j]>white_hilim || im[k][i][j]<white_lolim) 
          im[k][i][j] = 0;
      }
  
}

static void
read_images(char *fpref)
{
  int i,k;                   /* loop counters */
  FILE *fptr;
  char fname[100];

  numimg = imnr1-imnr0+1;
  bufsize = ((unsigned long)xnum)*ynum;
  buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
  for (k=0;k<numimg;k++)
  {
    im[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    for (i=0;i<IMGSIZE;i++)
    {
      im[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
    }
  }
  for (k=0;k<numimg;k++)
  {
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"rb");
    if (fptr==NULL) {printf("File %s not found.\n",fname);exit(1);}
    fread(buf,sizeof(char),bufsize,fptr);
    buffer_to_image(buf,im[k],xnum,ynum);
    fclose(fptr);
/*
  printf("file %s read\n",fname);
*/
  }
  printf("images %s read\n",fpref);
}

static void
write_images(char *fpref)
{
  int imnr;
  char fname[200];
  FILE *fptr;

  for (imnr=0;imnr<numimg;imnr++)
  {
    file_name(fpref,fname,imnr+imnr0,"%03d");
    fptr = fopen(fname,"wb");
    image_to_buffer(im[imnr],buf,IMGSIZE,IMGSIZE);
    fwrite(buf,sizeof(char),bufsize,fptr);
    fclose(fptr);
/*
    printf("File %s written\n",fname);
*/
  }
  printf("images %s written\n",fpref);
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
print_version(void)
{
  fprintf(stderr, "%s version %s\n", vcid, Progname) ;
  exit(1) ;
}
/*----------------------------------------------------------------------
            Parameters:

           Description:
----------------------------------------------------------------------*/
static void
print_help(void)
{
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
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;
  
  option = argv[1] + 1 ;            /* past '-' */
  if (!strcmp(option, "-help"))
    print_help() ;
  else if (!strcmp(option, "-version"))
    print_version() ;
  else switch (toupper(*option))
  {
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
