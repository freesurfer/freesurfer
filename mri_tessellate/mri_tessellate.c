
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "fio.h"
#include "const.h"
#include "diag.h"
#include "proto.h"
#include "error.h"
#include "MRIio.h"

#define SQR(x) ((x)*(x))
#define IMGSIZE     256
#define NUMVALS     256
#define MAXIM       256
#define MAXFACES    1000000
#define MAXVERTICES 1000000

typedef struct face_type_
{
  int imnr,i,j,f;
  int num;
  int v[4];
} face_type;

typedef struct vertex_type_
{
  int imnr,i,j;
  int num;
  int f[9];
} vertex_type;

face_type *face;
int *face_index_table0;
int *face_index_table1;

vertex_type *vertex;
int *vertex_index_table;

int xnum=256,ynum=256;
unsigned long bufsize;
unsigned char **im[MAXIM];  /* image matrix  */
unsigned char *buf;  /* scratch memory  */
int imnr0,imnr1,numimg;

float st,ps,fov,xx0,xx1,yy0,yy1,zz0,zz1;

int imin=0;
int imax=IMGSIZE;
int jmin=0;
int jmax=255;

int value;

int main(int argc, char *argv[]) ;
static void read_images(char *fpref) ;
static void add_face(int imnr, int i, int j, int f, int prev_flag) ;
static int add_vertex(int imnr, int i, int j) ;
static int facep(int im0, int i0, int j0, int im1, int i1, int j1) ;
static void check_face(int im0, int i0, int j0, int im1, int i1,int j1, 
                       int f, int n, int v_ind, int prev_flag) ;
static void make_surface(void) ;
static void write_binary_surface2(char *fname) ;
#if 0
static void write_surface(void) ;
static void write_binary_surface(char *fname) ;
#endif

char *Progname ;

int
main(int argc, char *argv[])
{
    char ifpref[STRLEN],ofpref[STRLEN] /*,*data_dir*/;
    char *getenv();

    DiagInit(NULL, NULL, NULL) ;
    ErrorInit(NULL, NULL, NULL) ;
    Progname = argv[0] ;
    if (argc<4)
    {
      printf("Usage: %s <input volume> <label> <output surface>\n",Progname);
      exit(1);
    }

    sscanf(argv[2],"%d",&value);
#if 0
    data_dir = getenv("SUBJECTS_DIR");
    if (data_dir==NULL)
    {
      printf("environment variable SUBJECTS_DIR undefined (use setenv)\n");
      exit(1);
    }
    sprintf(ifpref,"%s/%s/mri/filled/COR-",data_dir,argv[1]);
    sprintf(ofpref,"%s/%s/surf/%s.orig",data_dir,argv[1],argv[3]);
#else
    sprintf(ifpref,"%s/COR-",argv[1]);
    sprintf(ofpref,"%s",argv[3]);
#endif

    face = (face_type *)lcalloc(MAXFACES,sizeof(face_type));
    face_index_table0 = (int *)lcalloc(6*65536,sizeof(int));
    face_index_table1 = (int *)lcalloc(6*65536,sizeof(int));

    vertex = (vertex_type *)lcalloc(MAXVERTICES,sizeof(vertex_type));
    vertex_index_table = (int *)lcalloc(8*65536,sizeof(int));

    read_images(ifpref);

    make_surface();
#if 0
    write_surface();
    write_binary_surface(ofpref);
#endif
    write_binary_surface2(ofpref);
    exit(0) ;
    return(0) ;
}

static void
read_images(char *fpref)
{
  int i,j,k;                   /* loop counters */
  FILE *fptr;
  char fname[STRLEN];

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
/*  printf("%e %e %e %e %e %e %e %e %e\n",fov,ps,st,xx0,xx1,yy0,yy1,zz0,zz1); */

  fclose(fptr);
  numimg = imnr1-imnr0+1;

/* Allocate memory */

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

  for (k=1;k<numimg-1;k++)
  {
    file_name(fpref,fname,k+imnr0,"%03d");
          fptr = fopen(fname,"rb");
          if (fptr==NULL) {printf("File %s not found.\n",fname);exit(1);}
          fread(buf,sizeof(char),bufsize,fptr);
          buffer_to_image(buf,im[k],xnum,ynum);
          fclose(fptr);
  }

  for (k=1;k<numimg-1;k++)
        for (i=0;i<IMGSIZE;i++)
        for (j=0;j<IMGSIZE;j++)
          if (i<imin||i>=imax-1||j<jmin||j>=jmax-1)
            im[k][i][j] = 0;
}

static int face_index, vertex_index;

static void
add_face(int imnr, int i, int j, int f, int prev_flag)
{
  int pack = f*65536+i*256+j;

  if (face_index >= MAXFACES-1)
    ErrorExit(ERROR_NOMEMORY, "%s: max faces %d exceeded", 
              Progname,MAXFACES) ;
  if (prev_flag)
    face_index_table0[pack] = face_index;
  else
    face_index_table1[pack] = face_index;
  face[face_index].imnr = imnr;
  face[face_index].i = i;
  face[face_index].j = j;
  face[face_index].f = f;
  face[face_index].num = 0;
  face_index++;
}

static int
add_vertex(int imnr, int i, int j)
{
  int pack = i*257+j;

  if (vertex_index >= MAXVERTICES-1)
    ErrorExit(ERROR_NOMEMORY, "%s: max vertices %d exceeded", 
              Progname,MAXVERTICES) ;
  vertex_index_table[pack] = vertex_index;
  vertex[vertex_index].imnr = imnr;
  vertex[vertex_index].i = i;
  vertex[vertex_index].j = j;
  vertex[vertex_index].num = 0;
  return vertex_index++;
}

static int 
facep(int im0, int i0, int j0, int im1, int i1, int j1)
{
  return (im0>=0&&im0<numimg&&i0>=imin&&i0<imax&&j0>=jmin&&j0<jmax&&
          im1>=0&&im1<numimg&&i1>=imin&&i1<imax&&j1>=jmin&&j1<jmax&&
          im[im0][i0][j0] != im[im1][i1][j1] &&
          (im[im0][i0][j0]==value||im[im1][i1][j1]==value));
}
static void
check_face(int im0, int i0, int j0, int im1, int i1,int j1, int f, int n,
           int v_ind, int prev_flag)
{
  int f_pack = f*65536+i0*256+j0;
  int f_ind;

  if ((im0>=0&&im0<numimg&&i0>=imin&&i0<imax&&j0>=jmin&&j0<jmax&&
       im1>=0&&im1<numimg&&i1>=imin&&i1<imax&&j1>=jmin&&j1<jmax&&
       (im[im0][i0][j0]==value) && (im[im1][i1][j1]!=value)))
  {
    if (n==0)
    {
      add_face(im0,i0,j0,f,prev_flag);
    }
    if (prev_flag)
      f_ind = face_index_table0[f_pack];
    else
      f_ind = face_index_table1[f_pack];
    face[f_ind].v[n] = v_ind;
    if (vertex[v_ind].num<9)
      vertex[v_ind].f[vertex[v_ind].num++] = f_ind;
  }
}

static void
make_surface(void)
{
  int imnr,i,j,f_pack,v_ind,f;

  face_index = 0;
  vertex_index = 0;

  for (imnr=0;imnr<=numimg;imnr++)
  {
    if ((vertex_index || face_index) && !(imnr % 10))
      printf("slice %d: %d vertices, %d faces\n",imnr,vertex_index,face_index);
    for (i=0;i<=IMGSIZE;i++)
    for (j=0;j<=IMGSIZE;j++)
    {
      if (facep(imnr,i-1,j-1,imnr-1,i-1,j-1) ||
          facep(imnr,i-1,j,imnr-1,i-1,j) ||
          facep(imnr,i,j,imnr-1,i,j) ||
          facep(imnr,i,j-1,imnr-1,i,j-1) ||
          facep(imnr-1,i,j-1,imnr-1,i-1,j-1) ||
          facep(imnr-1,i,j,imnr-1,i-1,j) ||
          facep(imnr,i,j,imnr,i-1,j) ||
          facep(imnr,i,j-1,imnr,i-1,j-1) ||
          facep(imnr-1,i-1,j,imnr-1,i-1,j-1) ||
          facep(imnr-1,i,j,imnr-1,i,j-1) ||
          facep(imnr,i,j,imnr,i,j-1) ||
          facep(imnr,i-1,j,imnr,i-1,j-1))
      {
        v_ind = add_vertex(imnr,i,j);
        check_face(imnr  ,i-1,j-1,imnr-1,i-1,j-1,0,2,v_ind,0);
        check_face(imnr  ,i-1,j  ,imnr-1,i-1,j  ,0,3,v_ind,0);
        check_face(imnr  ,i  ,j  ,imnr-1,i  ,j  ,0,0,v_ind,0);
        check_face(imnr  ,i  ,j-1,imnr-1,i  ,j-1,0,1,v_ind,0);
        check_face(imnr-1,i  ,j-1,imnr-1,i-1,j-1,2,2,v_ind,1);
        check_face(imnr-1,i  ,j  ,imnr-1,i-1,j  ,2,1,v_ind,1);
        check_face(imnr  ,i  ,j  ,imnr  ,i-1,j  ,2,0,v_ind,0);
        check_face(imnr  ,i  ,j-1,imnr  ,i-1,j-1,2,3,v_ind,0);
        check_face(imnr-1,i-1,j  ,imnr-1,i-1,j-1,4,2,v_ind,1);
        check_face(imnr-1,i  ,j  ,imnr-1,i  ,j-1,4,3,v_ind,1);
        check_face(imnr  ,i  ,j  ,imnr  ,i  ,j-1,4,0,v_ind,0);
        check_face(imnr  ,i-1,j  ,imnr  ,i-1,j-1,4,1,v_ind,0);

        check_face(imnr-1,i-1,j-1,imnr  ,i-1,j-1,1,2,v_ind,1);
        check_face(imnr-1,i-1,j  ,imnr  ,i-1,j  ,1,1,v_ind,1);
        check_face(imnr-1,i  ,j  ,imnr  ,i  ,j  ,1,0,v_ind,1);
        check_face(imnr-1,i  ,j-1,imnr  ,i  ,j-1,1,3,v_ind,1);
        check_face(imnr-1,i-1,j-1,imnr-1,i  ,j-1,3,2,v_ind,1);
        check_face(imnr-1,i-1,j  ,imnr-1,i  ,j  ,3,3,v_ind,1);
        check_face(imnr  ,i-1,j  ,imnr  ,i  ,j  ,3,0,v_ind,0);
        check_face(imnr  ,i-1,j-1,imnr  ,i  ,j-1,3,1,v_ind,0);
        check_face(imnr-1,i-1,j-1,imnr-1,i-1,j  ,5,2,v_ind,1);
        check_face(imnr-1,i  ,j-1,imnr-1,i  ,j  ,5,1,v_ind,1);
        check_face(imnr  ,i  ,j-1,imnr  ,i  ,j  ,5,0,v_ind,0);
        check_face(imnr  ,i-1,j-1,imnr  ,i-1,j  ,5,3,v_ind,0);
      }
    }
    for (i=0;i<IMGSIZE;i++)
    for (j=0;j<IMGSIZE;j++)
    for (f=0;f<6;f++)
    {
      f_pack = f*256*256+i*256+j;
      face_index_table0[f_pack] = face_index_table1[f_pack];
    }
  }
}

#if 0
static void
write_surface(void)
{
  int k,n;
  float x,y,z;

  printf("%d %d\n",vertex_index,face_index);
  for (k=0;k<vertex_index;k++)
  {
    x = xx1-(vertex[k].j-0.5)*ps;
    y = yy0+(vertex[k].imnr-0.5)*st;
    z = zz1-(vertex[k].i-0.5)*ps;
/*
    if (k<10000 && (k%100)==0)
      printf("k=%d: imnr=%d i=%d j=%d x=%5.1f y=%5.1f z=%5.1f\n",
              k,vertex[k].imnr,vertex[k].i,vertex[k].j,x,y,z);
*/
    printf("%6.2f %6.2f %6.2f %d\n",x,y,z,vertex[k].num);
    for (n=0;n<vertex[k].num;n++)
      printf("%d ",vertex[k].f[n]);
    printf("\n");
  }
  for (k=0;k<face_index;k++)
  {
    for (n=0;n<4;n++)
      printf("%d ",face[k].v[n]);
    printf("\n");
  }
}
#endif

#if 0
fwrite1(v,fp)
int v;
FILE *fp;
{
  unsigned char c = (unsigned char)v;

  fwrite(&c,1,1,fp);
}

fwrite2(v,fp)
int v;
FILE *fp;
{
  short s = (short)v;

  fwrite(&s,2,1,fp);
}

fwrite3(v,fp)
int v;
FILE *fp;
{
  unsigned int i = (unsigned int)(v<<8);

  fwrite(&i,3,1,fp);
}
#endif

#if 0
static void
write_binary_surface(char *fname)
{
  int k,n;
  float x,y,z;
  FILE *fp;

  fp = fopen(fname,"w");
  fwrite3(vertex_index,fp);
  fwrite3(face_index,fp);
  for (k=0;k<vertex_index;k++)
  {
    x = xx1-(vertex[k].j-0.5)*ps;
    y = yy0+(vertex[k].imnr-0.5)*st;
    z = zz1-(vertex[k].i-0.5)*ps;
    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
    fwrite1(vertex[k].num,fp);
    for (n=0;n<vertex[k].num;n++)
      fwrite3(vertex[k].f[n],fp);
  }
  for (k=0;k<face_index;k++)
  {
    for (n=0;n<4;n++)
      fwrite3(face[k].v[n],fp);
  }
  fclose(fp);
}
#endif

static void
write_binary_surface2(char *fname)
{
  int k,n;
  float x,y,z;
  FILE *fp;

  fp = fopen(fname,"w");
  fwrite3(-1,fp);
  fwrite3(vertex_index,fp);
  fwrite3(face_index,fp);
  for (k=0;k<vertex_index;k++)
  {
    x = xx1-(vertex[k].j-0.5)*ps;
    y = yy0+(vertex[k].imnr-0.5)*st;
    z = zz1-(vertex[k].i-0.5)*ps;
    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
  }
  for (k=0;k<face_index;k++)
  {
    for (n=0;n<4;n++)
      fwrite3(face[k].v[n],fp);
  }
  fclose(fp);
}
