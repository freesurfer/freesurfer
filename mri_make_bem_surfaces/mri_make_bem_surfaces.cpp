/**
 * @brief creates a skull surface for use with the MNE tools
 *
 */
/*
 * Original Author: Anders Dale and Martin Sereno
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


/*============================================================================
 Copyright (c) 1996
=============================================================================*/
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include <sys/stat.h>

#include "proto.h"
#include "diag.h"
#include "error.h"
#include "const.h"
#include "icosahedron.h"
#include "version.h"
#include "MRIio_old.h" // lcalloc

/* prototypes */
/* static void read_datfile(char *fname) ; */
static void init_surf_to_image(void) ;
static void make_filenames(char *lsubjectsdir) ;
static void shrink(int niter, int nsmoothsteps) ;
static float rtanh(float x) ;
static void compute_normals(void) ;
static void normal_face(int f,float *n) ;
static void write_geometry(char *fname) ;
static void read_geometry(char *fname) ;
static void read_images(char *fpref) ;
static void read_image_info(char *fpref) ;
static void estimate_thickness(int niter);

#if 0
static void write_datfile(char *fname) ;
static void normal_vector(float *v0, float *v1, float *v2, float *norm) ;
#endif

#define TRUE 1
#define FALSE 0

#ifdef TCL
#  define PR printf("%% ");fflush(stdout);
#else
#  define PR
#endif

#ifndef SQR
#define SQR(x) ((x)*(x))
#endif
#define MATCH(A,B)   (!strcmp(A,B))
#define MATCH_STR(S) (!strcmp(str,S))

#define IMGSIZE 256

#define MAXIM 512
#define MAXVERTICES 10000
#define MAXFACES 10000
#define LOWTHRESH 50
#define NAME_LENGTH STRLEN

typedef struct ss_vertex_type_
{
  float x,y,z;
  float nx,ny,nz;
  float xb,yb,zb;
  float nxb,nyb,nzb;
  float ox,oy,oz;
  float dx,dy,dz;
  float mx,my,mz;
  float nc,onc,snc;
  float thickness;
  int vnum;
  int v[10];
  int fnum;
  int f[10];
}
ss_vertex_type;

ss_vertex_type vertex[MAXVERTICES];
int face[MAXFACES][3];
int nvertices,nfaces;

int xnum=256,ynum=256;
unsigned long bufsize;
unsigned char **im[MAXIM];  /* image matrix  */
unsigned char **fill[MAXIM];
unsigned char *buf;  /* scratch memory  */
int imnr0,imnr1,numimg;
int wx0=100,wy0=100;

float avgnc = 0;
float sf=0.55;

double istilt = 0.0;    /* Inside stilt length (mm) */
double ostilt = 0.5;    /* Outside stilt length (mm) */
double fstrength=1.0;
double fzero=40.0;
double fsteepness=0.50;
double mmbsfmax=10000;
double brainval = 80;
double skullval = 20;
double scalpval = 70;
double airval = 10;
int minscalpthickness = 5;
int minskullthickness = 3;
int maxskullthickness = 30;

float whitezero=35.0,grayzero=25.0;
float white_lolim = 60;
float gray_hilim = 70;
float threshold = 30;
float xmin,xmax;
float ymin,ymax;
float zmin,zmax;
float st,ps,fov,xx0,xx1,yy0,yy1,zz0,zz1;

int MRIflag = FALSE;
int MRIloaded = FALSE;
int momentumflag = TRUE;
int smoothflag = FALSE;
int centerflag = FALSE;
int flattenflag = FALSE;
int intersectionflag = FALSE;

float update = 0.9;
float decay = 0.9;
float xlo,xhi,ylo,yhi,zlo,zhi,xctr,yctr,zctr;
float ctrx,ctry,ctrz;

double dfrac = 0.7 ;

int changed=TRUE;
int shrinkmode=1;

char *subjectsdir;   /* SUBJECTS_DIR */
char *srname;        /* sessiondir */
char *pname;         /* name */
char *mfname;        /* abs image stem--input images (~/mri/T1/COR-) */
char *bfname;        /* abs image stem--stripped out images (~/mri/brain/COR-)*/
char *gfname;        /* in surface file */
char *ofname;        /* out surface file */
char *o2fname;        /* out surface file */
char *dfname;        /* datfile */
char *sgfname;       /* abs: sessiondir/rgb/tktrishrink.rgb */
char *rfname;        /* script */

int doublebufferflag = TRUE;
int openglwindowflag = FALSE;
int peeledflag = FALSE;
int initsurftoimageflag = FALSE;

const char *Progname ;
int
main(int argc,char *argv[])
{
  /* FILE *fptr; */
  char fpref[STRLEN];
  char *data_dir,*mri_dir;
  /* struct stat buf; */
  int nargs;

  nargs = handleVersionOption(argc, argv, "mri_make_bem_surfaces");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = argv[0] ;
  DiagInit(NULL, NULL, NULL) ;
  ErrorInit(NULL, NULL, NULL) ;

  if (argc<2)
  {
    printf("\nUsage: %s name [mfile]\n",argv[0]);
    printf("\n");
    printf("                                               [vers: 050304]\n");
    exit(0);
  }

  data_dir = getenv("SUBJECTS_DIR");
  if (data_dir==NULL)
  {
    printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
    exit(0);
  }
  mri_dir = getenv("FREESURFER_HOME");
  if (mri_dir==NULL)
  {
    printf("environment variable FREESURFER_HOME undefined (use setenv)\n");
    exit(0);
  }

  make_filenames(data_dir);

  sprintf(pname,"%s",argv[1]);
  sprintf(fpref,"%s/%s",data_dir,pname);
  sprintf(ofname, "%s/bem", fpref);
  mkdir(ofname, 0777);

#if 1
  /* Shrinkwrap "brain" -> initial estimate of inner_skull.tri */
  sprintf(mfname,"%s/mri/brain/COR-",fpref);
  sprintf(ofname,"%s/bem/inner_skull_tmp.tri",fpref);
  sprintf(gfname,"%s/lib/bem/ic4.tri",mri_dir);

  centerflag = TRUE;
  read_geometry(gfname);
  read_image_info(mfname);
  read_images(mfname);

  fzero = 10;
  mmbsfmax = 255;
  dfrac = 0.7;
  istilt = 0;
  fsteepness=0.50;
  fstrength=1.0;
  if (centerflag && MRIloaded)
  {
    init_surf_to_image();
  }

  shrinkmode = 4;
  momentumflag = TRUE;
  MRIflag = TRUE;
  flattenflag = FALSE;
  shrink(100,0);
  flattenflag = TRUE;
  MRIflag = FALSE;
  shrink(50,0);
  MRIflag = TRUE;
  flattenflag = FALSE;
  shrink(100,0);
  momentumflag = FALSE;
  MRIflag = FALSE;
  flattenflag = TRUE;
  shrink(10,0);
  write_geometry(ofname);


  /* Refine inner_skull.tri based on PD data */
  sprintf(mfname,"%s/mri/flash5/COR-",fpref);
  sprintf(ofname,"%s/bem/inner_skull.tri",fpref);
  sprintf(gfname,"%s/bem/inner_skull_tmp.tri",fpref);

  /*    if (stat(mfname, &buf)) {
        printf("Cannot find flash5 data\n");

        exit(0);
      }
  */
  read_geometry(gfname);
  read_image_info(mfname);
  read_images(mfname);

  shrinkmode = 5;
  MRIflag = TRUE;
  flattenflag = FALSE;
  momentumflag = FALSE;
  shrink(40,5);
  write_geometry(ofname);


  /* Estimate skull thickness for outer_skull.tri */
  estimate_thickness(10);
  sprintf(ofname,"%s/bem/outer_skull.tri",fpref);
  write_geometry(ofname);
#endif

  /* Find outer_skin.tri */
  sprintf(mfname,"%s/mri/T1/COR-",fpref);
  sprintf(ofname,"%s/bem/outer_skin.tri",fpref);
  /*
      sprintf(mfname,"%s/mri/flash5/COR-",fpref);
      sprintf(gfname,"%s/bem/outer_skull.tri",fpref); centerflag = FALSE;
  */
  sprintf(gfname,"%s/lib/bem/ic4.tri",mri_dir);
  centerflag = TRUE;
  read_geometry(gfname);
  read_image_info(mfname);
  read_images(mfname);

  fzero = 20;
  mmbsfmax = 1000;
  dfrac = 1.5;
  istilt = 0;
  fsteepness=0.20;
  fstrength=1.0;
  if (centerflag && MRIloaded)
  {
    init_surf_to_image();
  }

  shrinkmode = 6;
  MRIflag = TRUE;
  flattenflag = FALSE;
  momentumflag = TRUE;
  shrink(200,10);
  momentumflag = FALSE;
  shrink(100,10);
  /*
      shrinkmode = 7;
      shrink(50,0);
  */
  write_geometry(ofname);

  exit(0);
}

static void
read_image_info(char *fpref)
{
  FILE *fptr;
  char fname[STRLEN];

  sprintf(fname,"%s.info",fpref);
  fptr = fopen(fname,"r");
  if (fptr==NULL)
  {
    printf("trishrink: File %s not found\n",fname);
    printf("Best guess at inner_skull.tri is in inner_skull_tmp.tri\n");
    exit(0);
  }
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
read_images(char *fpref)
{
  int i,j,k;                   /* loop counters */
  FILE *fptr;
  char fname[STRLEN];

  if (!MRIloaded)
  {
    numimg = imnr1-imnr0+1;
    bufsize = ((unsigned long)xnum)*ynum;
    buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
    for (k=0; k<numimg; k++)
    {
      im[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
      for (i=0; i<IMGSIZE; i++)
      {
        im[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
      }
    }
  }
  else
  {
    read_image_info(fpref);  /* should check if same!! */
    for (k=0; k<numimg; k++)
      for (i=0; i<IMGSIZE; i++)
        for (j=0; j<IMGSIZE; j++)
        {
          im[k][i][j]=0;
        }
  }

  for (k=0; k<numimg; k++)
  {
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"r");
    if (fptr==NULL)
    {
      if (!MRIloaded)
      {
        printf("mri_strip_skull: ### File %s not found\n",fname);
        exit(0);
      }
      else
      {
        printf("mri_strip_skull: ### File %s not found\n",fname);
        PR return;
      }
    }
    fread(buf,sizeof(char),bufsize,fptr);
    buffer_to_image(buf,im[k],xnum,ynum);
    fclose(fptr);
  }
  printf("mri_strip_skull: images %s read\n",fpref);
  PR
  MRIloaded = TRUE;
  return ;
}

static void
init_surf_to_image(void)
{
  int i,j,k;
  float x,y,z;

  if (initsurftoimageflag)
  {
    printf(
      "mri_strip_skull: ### init_surf_to_image failed:  already done (re-read ic?.tri)\n");
    PR return;
  }

  if (!MRIloaded)
  {
    printf(
      "mri_strip_skull: ### init_surf_to_image failed:  MRI data not loaded\n");
    PR
    return;
  }
  xlo = ylo = zlo = 10000;
  xhi = yhi = zhi = -10000;
  for (k=25; k<numimg-25; k++)
    for (i=25; i<IMGSIZE-25; i++)
      for (j=25; j<IMGSIZE-25; j++)
      {
        x = xx1-j*ps;
        z = zz1-i*ps;
        y = yy0+k*st;
        if (im[k][i][j]>LOWTHRESH)
        {
          if (x<xlo)
          {
            xlo = x;
          }
          if (y<ylo)
          {
            ylo = y;
          }
          if (z<zlo)
          {
            zlo = z;
          }
          if (x>xhi)
          {
            xhi = x;
          }
          if (y>yhi)
          {
            yhi = y;
          }
          if (z>zhi)
          {
            zhi = z;
          }
        }
      }
  ctrx = (xlo+xhi)/2.0;
  ctry = (ylo+yhi)/2.0;
  ctrz = zhi-0.8*(xhi-xlo)/2.0;
  printf("xlo=%f, xhi=%f, ylo=%f, yhi=%f, zlo=%f, zhi=%f\n",
         xlo,xhi,ylo,yhi,zlo,zhi);
  PR
  for (k=0; k<nvertices; k++)
  {
    vertex[k].x = dfrac*vertex[k].x*(xhi-xlo)/2+ctrx;
    vertex[k].y = dfrac*vertex[k].y*1.3*(xhi-xlo)/2+ctry;
    vertex[k].z = dfrac*vertex[k].z*1.0*(xhi-xlo)/2+ctrz;
  }
  initsurftoimageflag = TRUE;
}

static void
read_geometry(char *fname)
{
  int i,j,k,n,last,next,skiplast,skipnext;

  FILE *fp;

#if 1
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("mri_strip_skull: ### cannot open file %s\n",fname);
    PR return;
  }
  fscanf(fp,"%d",&nvertices);
  for (k=0; k<nvertices; k++)
  {
    fscanf(fp,"%*d %f %f %f",
           &vertex[k].x,&vertex[k].y,&vertex[k].z);
    vertex[k].mx = vertex[k].my = vertex[k].mz = 0;
    vertex[k].vnum = 0;
    vertex[k].fnum = 0;  /* marty */
    vertex[k].nc = 0;
    vertex[k].snc = 0;
  }
  fscanf(fp,"%d",&nfaces);
  for (k=0; k<nfaces; k++)
  {
    fscanf(fp,"%*d");
    for (n=0; n<3; n++)
    {
      fscanf(fp,"%d",&face[k][n]);
      face[k][n]--;
    }
  }
  fclose(fp);
#else
  nvertices = ICO4_NVERTICES ;
  for (k=0; k<nvertices; k++)
  {
    vertex[k].x = ic2562_vertices[k].x ;
    vertex[k].y = ic2562_vertices[k].y ;
    vertex[k].z = ic2562_vertices[k].z ;
    vertex[k].mx = vertex[k].my = vertex[k].mz = 0;
    vertex[k].vnum = 0;
    vertex[k].fnum = 0;  /* marty */
    vertex[k].nc = 0;
    vertex[k].snc = 0;
  }

  nfaces = ICO4_NFACES ;
  for (k=0; k<nfaces; k++)
  {
    for (n=0; n<3; n++)
    {
      face[k][n] = ic2562_faces[k].vno[n] ;
      face[k][n]--;
    }
  }
#endif

  printf("nvertices=%d, nfaces=%d\n",nvertices,nfaces);

  for (k=0; k<nfaces; k++)
  {
    for (i=0; i<3; i++)
    {
      vertex[face[k][i]].f[vertex[face[k][i]].fnum++] = k;
      last = (i>0)?i-1:2;
      next = (i<2)?i+1:0;
      skiplast = skipnext = FALSE;
      for (j=0; j<vertex[face[k][i]].vnum; j++)
      {
        if (vertex[face[k][i]].v[j]==face[k][last])
        {
          skiplast = TRUE;
        }
        if (vertex[face[k][i]].v[j]==face[k][next])
        {
          skipnext = TRUE;
        }
      }
      if (!skiplast)
      {
        vertex[face[k][i]].v[vertex[face[k][i]].vnum++]=face[k][last];
      }
      if (!skipnext)
      {
        vertex[face[k][i]].v[vertex[face[k][i]].vnum++]=face[k][next];
      }
    }
  }
  compute_normals();
  for (k=0; k<nvertices; k++)
  {
    vertex[k].xb = vertex[k].x;
    vertex[k].yb = vertex[k].y;
    vertex[k].zb = vertex[k].z;
    vertex[k].nxb = vertex[k].nx;
    vertex[k].nyb = vertex[k].ny;
    vertex[k].nzb = vertex[k].nz;
  }
  initsurftoimageflag = FALSE;
  printf("mri_strip_skull: triangle file %s read\n",fname);
  PR
}

static void
write_geometry(char *fname)
{
  FILE *fp;
  int  k,n;

  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("mri_strip_skull: ### File %s not found\n",fname);
    PR return;
  }
  fprintf(fp,"%5d\n",nvertices);
  for (k=0; k<nvertices; k++)
  {
    fprintf(fp,"%5d%10.4f%10.4f%10.4f\n",k+1,vertex[k].x,vertex[k].y,vertex[k].z);
  }
  fprintf(fp,"%5d\n",nfaces);
  for (k=0; k<nfaces; k++)
  {
    fprintf(fp,"%5d",k+1);
    for (n=0; n<3; n++)
    {
      fprintf(fp,"%5d",face[k][n]+1);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  printf("mri_strip_skull: triangle file %s written\n",fname);
  PR
}

static void
normal_face(int f,float *n)
{
  float v1[3],v2[3],d;

  v1[0] = vertex[face[f][0]].x-vertex[face[f][1]].x;
  v1[1] = vertex[face[f][0]].y-vertex[face[f][1]].y;
  v1[2] = vertex[face[f][0]].z-vertex[face[f][1]].z;
  v2[0] = vertex[face[f][2]].x-vertex[face[f][1]].x;
  v2[1] = vertex[face[f][2]].y-vertex[face[f][1]].y;
  v2[2] = vertex[face[f][2]].z-vertex[face[f][1]].z;
  n[0] = v1[1]*v2[2]-v1[2]*v2[1];
  n[1] = v1[2]*v2[0]-v1[0]*v2[2];
  n[2] = v1[0]*v2[1]-v1[1]*v2[0];
  d = sqrt(SQR(n[0])+SQR(n[1])+SQR(n[2]));
  n[0] /= d;
  n[1] /= d;
  n[2] /= d;
}

static void
compute_normals(void)
{
  int j,k;
  ss_vertex_type *v;
  float n[3],nt[3];

  for (k=0; k<nvertices; k++)
  {
    v = &vertex[k];
    n[0] = n[1] = n[2] = 0;
    for (j=0; j<v->fnum; j++)
    {
      normal_face(v->f[j],nt);
      n[0] += nt[0];
      n[1] += nt[1];
      n[2] += nt[2];
    }
    v->nx = n[0]/v->fnum;
    v->ny = n[1]/v->fnum;
    v->nz = n[2]/v->fnum;
  }
}


static float
rtanh(float x)
{
  return (x<0.0)?0.0:tanh(x);
}

static void
shrink(int niter, int nsmoothsteps)
{
  float x,y,z,sx,sy,sz,val,inval,outval,nc,force,force0,force1,force2;
  float d,dx,dy,dz,sval,sinval,soutval,snc,inmean,inmax,outmean, outmax,sum,nsum;
  float nx,ny,nz;
  ss_vertex_type *v;
  int imnr,i,j,iter,k,m,n,smoothstep;
  float sd,ad,dmax;
  int navg,an,nclip,inim,ini,inj,outim,outi,outj;
  int delpos, mindelpos;
  int ninside=20,noutside=3;
  int maxdelpos = 30;
  float meanerr, minmeanerr, meanerr0;
  float insamp[50],outsamp[50];
  float delx=1.0,dely=1.0,delz=1.0,valt,xt,yt,zt;
  int imt,it,jt,h;

  if (MRIflag && !MRIloaded)
  {
    printf("mri_strip_skull: ### MRIflag but MRI data not loaded... reset\n");
    PR
    MRIflag = FALSE;
  }

  val = inval = outval = force = 0.0f ;  /* to stop compiler warnings */

  for (iter=0; iter<niter; iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0; k<nvertices; k++)
    {
      v = &vertex[k];
      v->ox = v->x;
      v->oy = v->y;
      v->oz = v->z;
      v->onc = v->nc;
    }
    sval = sinval = soutval = snc = 0;
    navg = 0;
    for (k=0; k<nvertices; k++)
    {
      v = &vertex[k];
      x = v->ox;
      y = v->oy;
      z = v->oz;
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;
      sx=sy=sz=sd=0;
      n=0;
      for (m=0; m<v->vnum; m++)
      {
        sx += dx = vertex[v->v[m]].ox - x;
        sy += dy = vertex[v->v[m]].oy - y;
        sz += dz = vertex[v->v[m]].oz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      if (n>0)
      {
        sx = sx/n;
        sy = sy/n;
        sz = sz/n;
        sd = sd/n;
        navg++;
      }
      if (MRIflag)
      {
        imnr = (int)((y-yy0)/st+0.5-imnr0);
        i = (int)((zz1-z)/ps+0.5);
        j = (int)((xx1-x)/ps+0.5);
        if (imnr<0||imnr>=numimg||i<0||i>=IMGSIZE||j<0||j>=IMGSIZE)
        {
          val = 0;
        }
        else
        {
          val = im[imnr][i][j];
        }
        inim = (int)(imnr-istilt/st*v->ny+0.5);
        ini = (int)(i+istilt/ps*v->nz+0.5);
        inj = (int)(j+istilt/ps*v->nx+0.5);
        if (inim<0||inim>=numimg||ini<0||ini>=IMGSIZE||inj<0||inj>=IMGSIZE)
        {
          inval = 0;
        }
        else
        {
          inval = im[inim][ini][inj];
        }
        outim = (int)(imnr+ostilt/st*v->ny+0.5);
        outi = (int)(i-ostilt/ps*v->nz+0.5);
        outj = (int)(j-ostilt/ps*v->nx+0.5);
        if (outim<0||outim>=numimg||
            outi<0||outi>=IMGSIZE||outj<0||outj>=IMGSIZE)
        {
          outval = 0;
        }
        else
        {
          outval = im[outim][outi][outj];
        }
        if (shrinkmode==3)
        {
          ninside = 0;
          noutside = 40;
          nc = (v->x-v->xb)*v->nxb+(v->y-v->yb)*v->nyb+(v->z-v->zb)*v->nzb;
          v->snc = nc;
          for (h= -noutside; h<ninside; h++)
          {
            xt = x-nx*(h*delx);
            yt = y-ny*(h*dely);
            zt = z-nz*(h*delz);
            imt = (int)((yt-yy0)/st+0.5-imnr0);
            it = (int)((zz1-zt)/ps+0.5);
            jt = (int)((xx1-xt)/ps+0.5);
            if (imt<0||imt>=numimg||it<0||it>=IMGSIZE||jt<0||jt>=IMGSIZE)
            {
              valt = 0;
            }
            else
            {
              valt = im[imt][it][jt];
            }
            if (h<=0)
            {
              outsamp[-h] = valt;
            }
            if (h>=0)
            {
              insamp[h] = valt;
            }
          }
          force2 = -1;
          for (h=1; h<10; h++)
          {
            valt = outsamp[h];
            force2 *= rtanh((valt-fzero)*fsteepness);
          }
          force1 = -1;
          for (h=5; h<noutside; h++)
          {
            valt = outsamp[h];
            force1 *= rtanh((fzero-valt)*fsteepness);
          }
          force0 = tanh((istilt-v->snc)*0.1*fsteepness);
          if (v->snc<istilt/2)
          {
            force1=force2=0;
          }
          force = 2.0*force2+1.0*force1+1.0*force0;
        }
        else if (shrinkmode==4)
        {
          force = fstrength*tanh((inval-fzero)*fsteepness);
        }
        else if (shrinkmode==5)
        {
          ninside=20;
          noutside=3;
          minmeanerr = 1e10;
          meanerr0 = 1e10;
          mindelpos = -maxdelpos; /* by Dav */
          for (delpos= -maxdelpos; delpos<maxdelpos; delpos++)
          {
            for (h= -noutside; h<ninside; h++)
            {
              xt = x-nx*(h*delx+istilt-delpos);
              yt = y-ny*(h*dely+istilt-delpos);
              zt = z-nz*(h*delz+istilt-delpos);
              imt = (int)((yt-yy0)/st+0.5-imnr0);
              it = (int)((zz1-zt)/ps+0.5);
              jt = (int)((xx1-xt)/ps+0.5);
              if (imt<0||imt>=numimg||it<0||it>=IMGSIZE||jt<0||jt>=IMGSIZE)
              {
                valt = 0;
              }
              else
              {
                valt = im[imt][it][jt];
              }
              if (h<=0)
              {
                outsamp[-h] = valt;
              }
              if (h>=0)
              {
                insamp[h] = valt;
              }
            }
            sum = nsum = 0;
            for (h=1; h<ninside; h++)
            {
              sum += insamp[h];
              nsum += 1;
            }
            inmean = sum/nsum;
            sum = nsum = 0;
            for (h=1; h<ninside; h++)
            {
              sum += SQR(insamp[h]-inmean);
              nsum += 1;
            }
            fzero = inmean*skullval/brainval;
            for (h=1; h<noutside; h++)
            {
              sum += SQR(outsamp[h]-fzero);
              nsum += 1;
            }
            meanerr = sqrt(sum/nsum)/inmean;
            if (meanerr<minmeanerr)
            {
              minmeanerr = meanerr;
              mindelpos = delpos;
            }
            if (delpos==0)
            {
              meanerr0 = meanerr;
            }
          }
          force1 = 0;
          if (mindelpos<0)
          {
            force1 = -1;
          }
          else if (mindelpos>0)
          {
            force1 = 1;
          }
          force = 0.5*force1*tanh(2.0*(meanerr0-minmeanerr)/meanerr0);
        }
        else if (shrinkmode==6)
        {
          ninside = 20;
          noutside = 20;
          nc = (v->x-v->xb)*v->nxb+(v->y-v->yb)*v->nyb+(v->z-v->zb)*v->nzb;
          v->snc = nc;
          for (h= -noutside; h<ninside; h++)
          {
            xt = x-nx*(h*delx);
            yt = y-ny*(h*dely);
            zt = z-nz*(h*delz);
            imt = (int)((yt-yy0)/st+0.5-imnr0);
            it = (int)((zz1-zt)/ps+0.5);
            jt = (int)((xx1-xt)/ps+0.5);
            if (imt<0||imt>=numimg||it<0||it>=IMGSIZE||jt<0||jt>=IMGSIZE)
            {
              valt = 0;
            }
            else
            {
              valt = im[imt][it][jt];
            }
            if (h<=0)
            {
              outsamp[-h] = valt;
            }
            if (h>=0)
            {
              insamp[h] = valt;
            }
          }
          outmax = -1;
          sum = nsum = 0;
          for (h=1; h<noutside; h++)
          {
            valt = outsamp[h];
            if (valt>outmax)
            {
              outmax = valt;
            }
            sum += valt;
            nsum++;
          }
          outmean = sum/nsum;
          inmax = -1;
          sum = nsum = 0;
          for (h=1; h<ninside; h++)
          {
            valt = insamp[h];
            if (valt>inmax)
            {
              inmax = valt;
            }
            sum += valt;
            nsum++;
          }
          inmean = sum/nsum;
          /*
                    force = fstrength*(0.5*tanh((insamp[0]-fzero)*fsteepness)+
                                       0.2*tanh((outmax-fzero)*fsteepness)+0.2*tanh((outmean-fzero)*fsteepness));\
          */
          force = fstrength*(0.2*tanh((insamp[0]-fzero)*fsteepness)+rtanh((inmax-fzero)*fsteepness)
                             -1.0*tanh((fzero-outmean)*fsteepness)-0.5*tanh((fzero-outmax)*fsteepness));
        }
        else if (shrinkmode==7)
        {
          ninside=5;
          noutside=20;
          minmeanerr = 1e10;
          meanerr0 = 1e10;
          mindelpos = -maxdelpos; /* by Dav */
          for (delpos= -maxdelpos; delpos<maxdelpos; delpos++)
          {
            for (h= -noutside; h<ninside; h++)
            {
              xt = x-nx*(h*delx+istilt-delpos);
              yt = y-ny*(h*dely+istilt-delpos);
              zt = z-nz*(h*delz+istilt-delpos);
              imt = (int)((yt-yy0)/st+0.5-imnr0);
              it = (int)((zz1-zt)/ps+0.5);
              jt = (int)((xx1-xt)/ps+0.5);
              if (imt<0||imt>=numimg||it<0||it>=IMGSIZE||jt<0||jt>=IMGSIZE)
              {
                valt = 0;
              }
              else
              {
                valt = im[imt][it][jt];
              }
              if (h<=0)
              {
                outsamp[-h] = valt;
              }
              if (h>=0)
              {
                insamp[h] = valt;
              }
            }
            sum = nsum = 0;
            for (h=1; h<ninside; h++)
            {
              sum += insamp[h];
              nsum += 1;
            }
            inmean = sum/nsum;
            sum = nsum = 0;
            for (h=1; h<ninside; h++)
            {
              sum += SQR(insamp[h]-inmean);
              nsum += 1;
            }
            fzero = inmean*airval/scalpval;
            for (h=1; h<noutside; h++)
            {
              sum += SQR(outsamp[h]-fzero);
              nsum += 1;
            }
            meanerr = sqrt(sum/nsum)/inmean;
            if (meanerr<minmeanerr)
            {
              minmeanerr = meanerr;
              mindelpos = delpos;
            }
            if (delpos==0)
            {
              meanerr0 = meanerr;
            }
          }
          force1 = 0;
          if (mindelpos<0)
          {
            force1 = -1;
          }
          else if (mindelpos>0)
          {
            force1 = 1;
          }
          force = 0.5*force1*tanh(2.0*(meanerr0-minmeanerr)/meanerr0);
        }
      }
      else
      {
        force = 0;
      }
      nc = sx*nx+sy*ny+sz*nz;
      sx -= nc*nx;
      sy -= nc*ny;
      sz -= nc*nz;
      snc += nc;
      v->nc = nc;
      v->snc += nc;
      if (flattenflag)
      {
        avgnc = 0;
        for (m=0; m<v->vnum; m++)
        {
          avgnc += vertex[v->v[m]].onc;
        }
        avgnc /= v->vnum;
        force += tanh((nc-avgnc)*0.5);
      }
      else
      {
        avgnc = 0;
        force += tanh((nc-avgnc)*0.1);
      }
      if ((d=sqrt(sx*sx+sy*sy+sz*sz))>1.0)
      {
        sx /= d;
        sy /= d;
        sz /= d;
      }
      dx = sx*0.5 + v->nx*force;
      dy = sy*0.5 + v->ny*force;
      dz = sz*0.5 + v->nz*force;
      if (momentumflag)
      {
        dx = decay*v->mx+update*dx;
        dy = decay*v->my+update*dy;
        dz = decay*v->mz+update*dz;
      }
      if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
      {
        dx /= d;
        dy /= d;
        dz /= d;
      }
      if (momentumflag)
      {
        v->mx = dx;
        v->my = dy;
        v->mz = dz;
      }
      d=sqrt(dx*dx+dy*dy+dz*dz);
      if (d>dmax)
      {
        dmax = d;
      }
      ad += d;
      an ++;
      v->dx = dx;
      v->dy = dy;
      v->dz = dz;
      sval += val;
      sinval += inval;
      soutval += outval;
    }

    for (smoothstep=0; smoothstep<nsmoothsteps; smoothstep++)
    {
      for (k=0; k<nvertices; k++)
      {
        v = &vertex[k];
        v->xb = v->dx;
        v->yb = v->dy;
        v->zb = v->dz;
      }
      for (k=0; k<nvertices; k++)
      {
        v = &vertex[k];
        sx=sy=sz=0;
        nsum = 0;
        for (m=0; m<v->vnum; m++)
        {
          sx += vertex[v->v[m]].xb;
          sy += vertex[v->v[m]].yb;
          sz += vertex[v->v[m]].zb;
          nsum++;
        }
        v->dx = sx/nsum;
        v->dy = sy/nsum;
        v->dz = sz/nsum;
      }
    }

    for (k=0; k<nvertices; k++)
    {
      v = &vertex[k];
      v->x += v->dx;
      v->y += v->dy;
      v->z += v->dz;
    }

    if (MRIflag)
    {
      compute_normals();
      sval /= nvertices;
      sinval /= nvertices;
      soutval /= nvertices;
      avgnc = snc/nvertices;
      if (!(iter % 10))
        printf("%d: sval=%5.2f,sinval=%5.2f,soutval=%5.2f\n",
               iter,sval,sinval,soutval);
    }
    else
    {
      printf("%d: ad=%f, dmax=%f, nclip=%d\n",iter,ad/an,dmax,nclip);
    }
  }
  compute_normals();
  PR
}

static void
estimate_thickness(int niter)
{
  float x,y,z,sx,sy,sz,val,inval,outval,/*nc,*/force /*,force0,force1,force2*/;
  float /* d,*/ dx,dy,dz,sval,sinval,soutval,snc,inmean,sum,nsum;
  float nx,ny,nz;
  ss_vertex_type *v;
  int imnr,i,j,iter,k,m,n=0;
  float sd,ad,dmax;
  int navg,an,nclip,inim,ini,inj,outim,outi,outj;
  int delpos, mindelpos;
  int ninside=20,noutside=40;
  float meanerr, minmeanerr, fskull, fscalp; /* meanerr0, */
  float insamp[50],outsamp[50];
  float delx=1.0,dely=1.0,delz=1.0,valt,xt,yt,zt;
  int imt,it,jt,h;

  if (MRIflag && !MRIloaded)
  {
    printf("mri_strip_skull: ### MRIflag but MRI data not loaded... reset\n");
    PR
    MRIflag = FALSE;
  }

  val = inval = outval = force = 0.0f ;  /* to stop compiler warnings */

  ad = 0;
  dmax = 0;
  an = 0;
  nclip = 0;
  for (k=0; k<nvertices; k++)
  {
    v = &vertex[k];
    v->ox = v->x;
    v->oy = v->y;
    v->oz = v->z;
    v->onc = v->nc;
  }
  sval = sinval = soutval = snc = 0;
  navg = 0;
  for (k=0; k<nvertices; k++)
  {
    v = &vertex[k];
    x = v->ox;
    y = v->oy;
    z = v->oz;
    nx = v->nx;
    ny = v->ny;
    nz = v->nz;
    sx=sy=sz=sd=0;
    n=0;
    for (m=0; m<v->vnum; m++)
    {
      sx += dx = vertex[v->v[m]].ox - x;
      sy += dy = vertex[v->v[m]].oy - y;
      sz += dz = vertex[v->v[m]].oz - z;
      sd += sqrt(dx*dx+dy*dy+dz*dz);
      n++;
    }
    if (n>0)
    {
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
      sd = sd/n;
      navg++;
    }
    if (MRIflag)
    {
      imnr = (int)((y-yy0)/st+0.5-imnr0);
      i = (int)((zz1-z)/ps+0.5);
      j = (int)((xx1-x)/ps+0.5);
      if (imnr<0||imnr>=numimg||i<0||i>=IMGSIZE||j<0||j>=IMGSIZE)
      {
        val = 0;
      }
      else
      {
        val = im[imnr][i][j];
      }
      inim = (int)(imnr-istilt/st*v->ny+0.5);
      ini = (int)(i+istilt/ps*v->nz+0.5);
      inj = (int)(j+istilt/ps*v->nx+0.5);
      if (inim<0||inim>=numimg||ini<0||ini>=IMGSIZE||inj<0||inj>=IMGSIZE)
      {
        inval = 0;
      }
      else
      {
        inval = im[inim][ini][inj];
      }
      outim = (int)(imnr+ostilt/st*v->ny+0.5);
      outi = (int)(i-ostilt/ps*v->nz+0.5);
      outj = (int)(j-ostilt/ps*v->nx+0.5);
      if (outim<0||outim>=numimg||
          outi<0||outi>=IMGSIZE||outj<0||outj>=IMGSIZE)
      {
        outval = 0;
      }
      else
      {
        outval = im[outim][outi][outj];
      }
      for (h= -noutside; h<ninside; h++)
      {
        xt = x-nx*(h*delx+istilt);
        yt = y-ny*(h*dely+istilt);
        zt = z-nz*(h*delz+istilt);
        imt = (int)((yt-yy0)/st+0.5-imnr0);
        it = (int)((zz1-zt)/ps+0.5);
        jt = (int)((xx1-xt)/ps+0.5);
        if (imt<0||imt>=numimg||it<0||it>=IMGSIZE||jt<0||jt>=IMGSIZE)
        {
          valt = 0;
        }
        else
        {
          valt = im[imt][it][jt];
        }
        if (h<=0)
        {
          outsamp[-h] = valt;
        }
        if (h>=0)
        {
          insamp[h] = valt;
        }
      }
      sum = nsum = 0;
      for (h=1; h<ninside; h++)
      {
        sum += insamp[h];
        nsum += 1;
      }
      inmean = sum/nsum;
      minmeanerr = 1e10;
      mindelpos = minskullthickness;
      for (delpos= minskullthickness; delpos<maxskullthickness; delpos++)
      {
        sum = nsum = 0;
        for (h=1; h<ninside; h++)
        {
          sum += SQR(insamp[h]-inmean);
          nsum += 1;
        }
        fskull = inmean*skullval/brainval;
        fscalp = inmean*scalpval/brainval;
        for (h=1; h<delpos; h++)
        {
          sum += SQR(outsamp[h]-fskull);
          nsum += 1;
        }
        for (h=delpos+1; h<delpos+1+minscalpthickness; h++)
        {
          sum += SQR(outsamp[h]-fscalp);
          nsum += 1;
        }
        meanerr = sqrt(sum/nsum);
        if (meanerr<minmeanerr)
        {
          minmeanerr = meanerr;
          mindelpos = delpos;
        }
      }
      v->thickness = mindelpos;
    }
  }
  for (k=0; k<nvertices; k++)
  {
    v = &vertex[k];
    v->ox = v->nx*v->thickness;
    v->oy = v->ny*v->thickness;
    v->oz = v->nz*v->thickness;
  }
  for (iter=0; iter<niter; iter++)
  {
    for (k=0; k<nvertices; k++)
    {
      v = &vertex[k];
      v->xb = v->ox;
      v->yb = v->oy;
      v->zb = v->oz;
    }
    for (k=0; k<nvertices; k++)
    {
      v = &vertex[k];
      sx=sy=sz=0;
      nsum = 0;
      for (m=0; m<v->vnum; m++)
      {
        sx += vertex[v->v[m]].xb;
        sy += vertex[v->v[m]].yb;
        sz += vertex[v->v[m]].zb;
        nsum++;
      }
      v->ox = sx/nsum;
      v->oy = sy/nsum;
      v->oz = sz/nsum;
    }
  }
  for (k=0; k<nvertices; k++)
  {
    v = &vertex[k];
    v->x += v->ox;
    v->y += v->oy;
    v->z += v->oz;
  }
}
static void
make_filenames(char *lsubjectsdir)
{
  subjectsdir = (char *)malloc(NAME_LENGTH*sizeof(char));/* malloc for tcl */
  srname = (char *)malloc(NAME_LENGTH*sizeof(char));
  pname = (char *)malloc(NAME_LENGTH*sizeof(char));
  mfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  bfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  gfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  ofname = (char *)malloc(NAME_LENGTH*sizeof(char));
  o2fname = (char *)malloc(NAME_LENGTH*sizeof(char));
  dfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  sgfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  rfname = (char *)malloc(NAME_LENGTH*sizeof(char));

  strcpy(subjectsdir,lsubjectsdir);
  /* TODO: fix rest of infile/outfile/datfile mess, init here */
}

/* static void
read_datfile(char *fname)
{
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL) {printf("mri_strip_skull: ### File %s not found\n",fname);PR return;}
  fscanf(fp,"%*s %lf",&brainval);
  fscanf(fp,"%*s %lf",&skullval);
  fscanf(fp,"%*s %lf",&scalpval);
  fscanf(fp,"%*s %lf",&istilt);
  fclose(fp);
  printf("mri_strip_skull: datfile read: %s\n",fname);
  printf(" brainval %f\n",brainval);
  printf(" skullval %f\n",skullval);
  printf(" scalpval %f\n",scalpval);
  printf(" istilt %f\n",istilt);PR
}

static char *
lcalloc(size_t nmemb, size_t size)
{
  char *p;

  p = calloc(nmemb,size);
  if (p==NULL) print_error("Cannot calloc()\n");
  return p;
}


static void
file_name(char *fpref, char *fname, int num, char *form)
{
  char ext[100];

  sprintf(ext,form,num);
  strcpy(fname,fpref);
  strcat(fname,ext);
}


static void
buffer_to_image(unsigned char *buf, unsigned char **im, int ysize, int xsize)
{
  int i,j;
  unsigned long k;
  float sum;

  k=0;
  sum = 0;
  for (i=0;i<ysize;i++)
    for (j=0;j<xsize;j++)
    {
      im[i][j] = buf[k++];
      sum += im[i][j];
    }
}
*/
