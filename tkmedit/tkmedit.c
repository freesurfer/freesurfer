/*============================================================================
 Copyright (c) 1996 Martin Sereno and Anders Dale
=============================================================================*/
#define TCL
#define TKMEDIT 
/*#if defined(Linux) || defined(sun) || defined(SunOS) */
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
/*#endif */
#include <fcntl.h>
#include "error.h"
#include "utils.h"



/*------------------begin medit.c-------------------*/

/*============================================================================
 Copyright (c) 1996, 1997 Anders Dale and Martin Sereno
=============================================================================*/
#include <tcl.h>
#include <tk.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <sys/time.h>
#include "MRIio.h"
#include "volume_io.h"
#include "rgb_image.h"
#include "fio.h"
#include "mrisurf.h"

#define OPENGL
#define TCL



#ifdef TCL
#  define PR   {if(promptflag){fputs("% ", stdout);} fflush(stdout);}
int            tk_NumMainWindows = 0;
#else
#  define PR
#endif


#ifdef OPENGL
#include "proto.h"
#include "macros.h"
#include "xwindow.h"
#include <X11/keysym.h>
#  define bgnpolygon()     glBegin(GL_QUADS)  /*one list OK;GL_POLYGON slower*/
#  define bgnline()        glBegin(GL_LINES)
#  define bgnpoint()       glBegin(GL_POINTS)
#  define v2f(X)           glVertex2fv(X)
#  define endpolygon()     glEnd()
#  define endline()        glEnd()
#  define endpoint()       glEnd()
#  define swapbuffers()    tkoSwapBuffers()
#  define clear()          glClear(GL_COLOR_BUFFER_BIT)
#  define linewidth(X)     glLineWidth((float)(X))
#  define getorigin(X,Y)   *(X) = w.x; *(Y) = 1024 - w.y - w.h /*WINDOW_REC w;*/
#  define getsize(X,Y)     *(X) = w.w; *(Y) = w.h
#  define Colorindex       unsigned short
#  define color(R)         glIndexs(X)
#  define mapcolor(I,R,G,B)  \
           tkoSetOneColor((int)I,(float)R/255.0,(float)G/255.0,(float)B/255.0)
#  define rectwrite(X0,Y0,X1,Y1,P) \
           glRasterPos2i(X0,Y0); \
           glDrawPixels((X1-X0)+1,(Y1-Y0)+1,GL_RGBA,GL_UNSIGNED_BYTE,P) 
/*
           glDrawPixels((X1-X0)+1,(Y1-Y0)+1,GL_LUMINANCE,GL_UNSIGNED_BYTE,P) 
           glDrawPixels((X1-X0)+1,(Y1-Y0)+1,GL_COLOR_INDEX,GL_UNSIGNED_SHORT,P) */ 
/* best Extreme glDrawPixels: GL_ABGR_EXT,GL_UNSIGNED_BYTE */
#  define wintitle(X) \
           XStringListToTextProperty(&X,1,&tp); \
           XSetWMName(xDisplay,w.wMain,&tp)
#  define ortho2(X0,X1,Y0,Y1) \
           glMatrixMode(GL_PROJECTION); \
           glLoadIdentity(); \
           glOrtho(X0,X1,Y0,Y1,-1.0,1.0);
#  define BLACK   0  /* TODO: use tkoSetOneColor, glIndexs */
#  define GREEN   1
#  define RED     2
#  define YELLOW  3
#else
#  include <gl.h>
#  include <device.h>
#endif
#ifdef Linux
#define CURSOR_VAL   255
#else
#define CURSOR_VAL   255
#endif

/* #define SQR(x)       ((x)*(x)) */
#define MATCH(A,B)   (!strcmp(A,B))
#define MATCH_STR(S) (!strcmp(str,S))

#define NUMVALS 256
#define MAXIM 256
#define MAXPTS 10000
#define MAXPARS 10
#define MAPOFFSET 0
#define CORONAL    0
#define HORIZONTAL 1
#define SAGITTAL   2
#define POSTANT   0
#define INFSUP    1
#define RIGHTLEFT 2  /* radiol */
#define NAME_LENGTH  100
#define MAX_DIR_DEPTH  30
#define TMP_DIR          "tmp"             /* relative to subjectsdir/pname */
#define TRANSFORM_DIR    "mri/transforms"  /* ditto */
#define TALAIRACH_FNAME  "talairach.xfm"   /* relative to TRANSFORM_DIR */
#define WM_EDITED_OFF 1
#define WM_MIN_VAL    2  /* 1 is used for voxels that are edited to off */
#define TO_WHITE  0
#define TO_BLACK  1
#define DIR_FILE "ic1.tri"
#define MAXCOR 500
#define MAXLEN 100
#define MOTIF_XFUDGE   8
#define MOTIF_YFUDGE  32
#define CVIDBUF 25

int selectedpixval = 0;
int updatepixval = TRUE;
int plane = CORONAL;
int xnum=256,ynum=256;
int ptype;
float ps,st,xx0,xx1,yy0,yy1,zz0,zz1;
int zf, ozf;
float fsf;
static int xdim,ydim;
unsigned long bufsize;
unsigned char **im[MAXIM];
unsigned char **im_b[MAXIM];
unsigned char **im2[MAXIM];
unsigned char **fill[MAXIM];
unsigned char **dummy_im[MAXIM];
unsigned char **sim[6]; 
unsigned char **sim2[6]; 
int second_im_allocated = FALSE;
int dummy_im_allocated = FALSE;
int wmfilter_ims_allocated = FALSE;
int changed[MAXIM];
GLubyte *vidbuf; 
/* Colorindex *vidbuf;  */
Colorindex *cvidbuf;
unsigned char *buf;
unsigned char *binbuff;
int imnr0,imnr1,numimg;
int wx0=114,wy0=302;  /* (100,100), (117,90), (556,90) */
int ptsflag = FALSE;
int maxflag = FALSE;
int surfflag = FALSE;
int surfloaded = FALSE;
int editflag = TRUE;
int revfsflag = FALSE;
int fieldsignloaded = FALSE;   /* fscontour */
int fieldsignflag = FALSE; /* overrides curvflag */
int surflinewidth = 1;
int curvloaded = FALSE;
int curvflag = FALSE;
int editedimage = FALSE;
int drawsecondflag = FALSE;
int inplaneflag = TRUE;
int linearflag = FALSE;
int bwflag = FALSE;
int truncflag = FALSE;
int scrsaveflag = TRUE;
int openglwindowflag = FALSE;
int second_im_full = FALSE;
int circleflag = TRUE;
int promptflag = FALSE;
int followglwinflag = TRUE;
int initpositiondoneflag = FALSE;
int all3flag = FALSE;
int npts = 0;
int prad = 0;
int pradlast = 0;
int ndip = 0;
int dip_spacing = 10; /* voxels */
float tm[4][4];
int jold[MAXPTS],iold[MAXPTS],imold[MAXPTS];
float ptx[MAXPTS],pty[MAXPTS],ptz[MAXPTS];
float par[MAXPARS],dpar[MAXPARS];

static MRI_SURFACE *mris ;
double fthresh = 0.35;
double fsquash = 12.0;
double fscale = 255;
double fsthresh = 0.3;
double xtalairach = 0.0;
double ytalairach = 0.0;
double ztalairach = 0.0;
int white_lolim = 80;
int white_hilim = 140;
int gray_hilim = 100;
int flossflag = TRUE;
int spackleflag = TRUE;
int lim3=170,lim2=145,lim1=95,lim0=75;
double ffrac3=1.0,ffrac2=1.0,ffrac1=1.0,ffrac0=1.0;

int imc=0,ic=0,jc=0;
int impt = -1,ipt = -1,jpt = -1;
float x_click, y_click, z_click ;

char *subjectsdir;   /* SUBJECTS_DIR */
char *srname;        /* sessiondir--from cwd */
char *pname;         /* name */
char *imtype;        /* e.g., T1 */
char *imtype2;        /* e.g., T1 */
char *surface;       /* e.g., rh.smooth */
char *mfname;        /* abs single image stem name */
char *sfname;        /* abs surface name */
char *tfname;        /* (dir!) subjectsdir/name/tmp/ */
char *dipfname;      /* curr: $session/bem/brain3d.dip */
char *decfname;      /* curr: $session/bem/brain3d.dec */
char *hpfname;       /* headpts */
char *htfname;       /* headtrans */
char *sgfname;       /* rgb */
char *fsfname;       /* $session/fs/$hemi.fs */
char *fmfname;       /* $session/fs/$hemi.fm */
char *cfname;        /* $home/surf/hemi.curv */
char *xffname;       /* $home/name/mri/transforms/TALAIRACH_FNAME */
char *rfname;        /* script */

/* Talairach stuff */
General_transform talairach_transform ; /* the next two are from this struct */
Transform         *linear_transform = NULL ;
Transform         *inverse_linear_transform = NULL ;
int               transform_loaded = 0 ;
/*--------------------- prototypes ------------------------------*/
#ifdef Linux
extern void scale2x(int, int, unsigned char *);
#endif

#ifdef TCL
void do_one_gl_event(Tcl_Interp *interp) ;
int Medit(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[]);
#else
void main(int argc,char *argv[]) ;

#endif
void usecnap(int usec) ;
int read_second_images(char *imdir2) ;
void resize_window_intstep(int newzf) ;
void write_point(char *dir) ;
void goto_vertex(int vno) ;
void mark_file_vertices(char *fname) ;
void unmark_vertices(void) ;
void goto_point_coords(int imc1,int ic1,int jc1) ;
void set_cursor(float xpt,float ypt,float zpt) ;
void mri2pix(float xpt, float ypt, float zpt, int *jpt, int *ipt,int *impt);
void rotate_brain(float a,char c) ;
void translate_brain(float a,char c) ;
void optimize2(void) ;
void optimize(int maxiter) ;
float Error(int p,float dp) ;
int  imval(float px,float py,float pz) ;
void resize_buffers(int x, int y) ;
void pop_gl_window(void) ;
void redraw(void) ;
void set_scale(void) ;
void goto_point(char *dir) ;
void upslice(void) ;
void downslice(void) ;
void pix_to_rgb(char *fname) ;
void scrsave_to_rgb(char *fname) ;
void save_rgb(char *fname) ;
void edit_pixel(int action);
int open_window(char *name) ;
void draw_surface(void) ;
int read_binary_surface(char *fname) ;
int read_surface(char *fname) ;
int read_binary_surf(char *fname) ;
int show_vertex(void) ;
int dump_vertex(int vno) ;
int write_images(char *fpref) ;
int read_images(char *fpref) ;
void select_pixel(short sx, short sy, int printflag) ;
void draw_image(int imc, int ic, int jc) ;
void draw_second_image(int imc, int ic, int jc) ;
void drawpts(void) ;
void write_dipoles(char *fname) ;
void write_decimation(char *fname) ;
void read_hpts(char *fname) ;
void read_htrans(char *fname) ;
void write_htrans(char *fname) ;
void make_filenames(char *lsubjectsdir,char *lsrname, char *lpname, 
                    char *limdir, char *lsurface) ;
void mirror(void) ;
void read_fieldsign(char *fname) ;
void read_fsmask(char *fname) ;
void read_binary_curvature(char *fname) ;
void smooth_3d(int niter) ;
void flip_corview_xyz(char *newx, char *newy, char *newz) ;
void wmfilter_corslice(int imc) ;
void sagnorm_allslices(void) ;
void sagnorm_corslice(int imc) ;
void alloc_second_im(void) ;

char *Progname ;
/*--------------------- end prototypes ------------------------------*/

/*--------------------- twitzels hacks ------------------------------*/

/*--------------------- functional ----------------------------------*/
#ifndef FUNCTIONAL_C
#define FUNCTIONAL_C

#include <errno.h>

/* global variables */
int fvwidth,fvheight,fvslices,fvframes;

float fux, fuy, fuz;
float fvx, fvy, fvz;
float fcx, fcy, fcz;

unsigned char* scaledFVolume;
float** rawFVData;
fMRI_REG* fvRegistration;
float mri2fmritm[4][4];

int statsloaded;
int hot;
int do_overlay = 0;
int do_interpolate = 0;
int overlay_frame = 0;

float slfvps;
float slfvst;

double fslope = 1.00;
double fmid = 0.0;
double f2thresh =0.0;

char* overlay_file = NULL;
char* loaded_fvfile;

typedef union 
{
  long  l;
  float f;
  int i;
  char buf[4];
  short s[2];
} t_uni4;

void transformFV(float a, float b, float c, float* d, float* e, float* f);

void cleanFOV()
{
  int i;
  for(i=0; i<fvslices; i++) {
    free(rawFVData[i]);
  }
  
  free(rawFVData);
}

void setupFOVars()
{
}

void setupFVU(float tux, float tuy, float tuz)
{
  transformFV(tux,tuy,tuz,&fux,&fuy,&fuz);
  fux = fux-fcx;
  fuy = fuy-fcy;
  fuz = fuz-fcz;
}

void setupFVV(float tvx, float tvy, float tvz)
{
  transformFV(tvx,tvy,tvz,&fvx,&fvy,&fvz);
  fvx = fvx-fcx;
  fvy = fvy-fcy;
  fvz = fvz-fcz;
}

void setupFVC(float tcx, float tcy, float tcz)
{
  if(hot)
    printf("setup centerpoint to: %f %f %f\n",tcx,tcy,tcz);
  transformFV(tcx,tcy,tcz,&fcx,&fcy,&fcz);
  if(hot)
    printf("centerpoint is: %f %f %f\n",fcx,fcy,fcz);
}


int countFVSlices(const char* prefix)
{
  int slice_number = 0;
  char fname[255];
  FILE* fp;

  do
  {
    sprintf(fname, "%s_%3.3d.bfloat", prefix, slice_number) ;
    fp = fopen(fname, "r") ;
    if (fp)   /* this is a valid slice */
    {
      fclose(fp) ;
      slice_number++ ;
    } else
      break;
  } while (1) ;

  return slice_number;
}

/* TODO: fvframes ??? */
void readFVSliceHeader(const char* filename) 
{
  FILE* fp;
  int width, height, nframes;

  fp = fopen(filename, "r");
  
  if(!fp)
    printf("could not open %s !\n",filename);
  else {
    fscanf(fp, "%d %d %d", &width, &height, &nframes);
    fvwidth = width; fvheight = height; fvframes = nframes;
    fclose(fp);
  }
}

float swapFVFloat(float value)
{

  t_uni4 fliess;
  char tmp;
  short tmp2;

  fliess.f = value;
  tmp = fliess.buf[0];
  fliess.buf[0] = fliess.buf[1];
  fliess.buf[1] = tmp;

  tmp = fliess.buf[2];
  fliess.buf[2] = fliess.buf[3];
  fliess.buf[2] = tmp;

  tmp2 = fliess.s[0];
  fliess.s[0]= fliess.s[1];
  fliess.s[1]=tmp2;
 
  return fliess.f;
}

void transformFV(float x, float y, float z, float* x_1, float* y_1,float* z_1)
{

  *x_1 = x*mri2fmritm[0][0]+y*mri2fmritm[0][1]+z*mri2fmritm[0][2]+mri2fmritm[0][3];
  *y_1 = x*mri2fmritm[1][0]+y*mri2fmritm[1][1]+z*mri2fmritm[1][2]+mri2fmritm[1][3];
  *z_1 = x*mri2fmritm[2][0]+y*mri2fmritm[2][1]+z*mri2fmritm[2][2]+mri2fmritm[2][3];
  
  if(hot)
    printf("after mult: %f, %f, %f\n",*x_1,*y_1,*z_1);
  *x_1 = ((fvRegistration->in_plane_res*fvwidth)/2.0-(*x_1))/fvRegistration->in_plane_res;
  *z_1 = fvheight-((fvRegistration->in_plane_res*fvheight)/2.0-(*z_1))/fvRegistration->in_plane_res;
  *y_1 = (*y_1-(-fvRegistration->slice_thickness*fvslices)/2.0)/fvRegistration->slice_thickness;
}

/* TODO: swaponly on SGI's */
void readFVVolume(const char* prefixname)
{
  int fvi;
  int fvj;
  int fvsize;
  int fvl;
  char fvfname[255];
  int fvfile;
  float *fvptr;

  hot = 0;

  fvsize = fvwidth*fvheight*fvframes*sizeof(float);
  rawFVData = (float**)malloc(fvslices*sizeof(float*));

  for(fvi=0; fvi<fvslices; fvi++) {
    fvptr = (float*)malloc(fvwidth*fvheight*fvframes*sizeof(float));
    sprintf(fvfname, "%s_%3.3d.bfloat", prefixname, fvi);
    fvfile = open(fvfname, 0);
    if(fvfile) {
      fvl=read(fvfile,fvptr,fvsize);
      printf("slice %d read %d bytes\n",fvi,fvl);
      if( fvl<fvsize ) {
  printf("Scheisse: %s\n",strerror(errno));
      }
      close(fvfile);
      rawFVData[fvi]=fvptr;
#ifdef Linux      
      for(fvj=0; fvj<fvwidth*fvheight*fvframes; fvj++) {
  fvptr[fvj]=swapFVFloat(fvptr[fvj]);
      }
#endif
      
    } else {
      printf("Ooops! could not open file %s\n",fvfname);
      break;
    }
  }
  slfvps = fvRegistration->in_plane_res;
  slfvst = fvRegistration->slice_thickness;
}

void copyFVMatrix()
{
  MATRIX* tmp;

  tmp = fvRegistration->mri2fmri;

  mri2fmritm[0][0] = tmp->data[0];
  mri2fmritm[0][1] = tmp->data[1];
  mri2fmritm[0][2] = tmp->data[2];
  mri2fmritm[0][3] = tmp->data[3];
  mri2fmritm[1][0] = tmp->data[4];
  mri2fmritm[1][1] = tmp->data[5];
  mri2fmritm[1][2] = tmp->data[6];
  mri2fmritm[1][3] = tmp->data[7];
  mri2fmritm[2][0] = tmp->data[8];
  mri2fmritm[2][1] = tmp->data[9];
  mri2fmritm[2][2] = tmp->data[10];
  mri2fmritm[2][3] = tmp->data[11];
  mri2fmritm[3][0] = tmp->data[12];
  mri2fmritm[3][1] = tmp->data[13];
  mri2fmritm[3][2] = tmp->data[14];
  mri2fmritm[3][3] = tmp->data[15];  
}

void loadFV()
{
  char slicename[255];
  char* prefixname;
  fvRegistration = NULL;
  
  if(overlay_file == NULL) {
    printf("No overlay filename given ! please set funname !\n");
    return;
  }
  prefixname = overlay_file;

  fvRegistration = StatReadRegistration("register.dat");
 
  if(!fvRegistration) {
    printf("Could not load registration file !\n");
    exit(-1);
  }   
  fvslices = countFVSlices(prefixname);
  
  printf("found %d slices\n", fvslices);

  sprintf(slicename, "%s_%3.3d.hdr", prefixname,0);
  readFVSliceHeader(slicename);
  
  printf("read format: %d, %d, %d\n",fvwidth,fvheight,fvframes);
  copyFVMatrix();

  readFVVolume(prefixname);

  statsloaded = 1;
  loaded_fvfile = strdup(overlay_file);
}

/* sample Data interpolates trilinear */

#define V000(o) rawFVData[unten][vorn+links+o] 
#define V100(o) rawFVData[unten][vorn+rechts+o]
#define V010(o) rawFVData[unten][hinten+links+o]
#define V001(o) rawFVData[oben][vorn+links+o]
#define V101(o) rawFVData[oben][vorn+rechts+o]
#define V011(o) rawFVData[oben][hinten+links+o]
#define V110(o) rawFVData[unten][hinten+rechts+o]
#define V111(o) rawFVData[oben][hinten+rechts+o]

#define NO_VALUE  -30000.0f

float sampleData(float x, float y, float z,int frame)
{
  int ux,lx,uy,ly,uz,lz; /* Koordinaten der Nachbarvoxel */
  int rechts, links, vorn, hinten, oben, unten;
  float v;
  int offset;

  ux = ceilf(x); lx = floorf(x);
  uy = ceilf(y); ly = floorf(y);
  uz = ceilf((fvheight-1-z)); lz = floorf(fvheight-1-z);

  /* printf("Hampf: %d %d\n", uz, lz); */

  x = x - lx;
  y = y - ly;
  z = (fvheight-1-z) - lz;
  
  rechts = ux;
  links = lx;
  oben = uy;
  unten = ly; 
  vorn =  fvwidth*lz;
  hinten = fvwidth*uz;

  offset = frame*fvwidth*fvheight;

  if(ux >= fvwidth || uz >= fvheight || uy >= fvslices || lx < 0 || ly < 0 || lz < 0) {
     
      return NO_VALUE;
  }
    
  v= V000(offset)*(1-x)*(1-z)*(1-y) +
    V100(offset)*x*(1-z)*(1-y)+
    V010(offset)*(1-x)*z*(1-y)+
    V001(offset)*(1-x)*(1-z)*y+
    V101(offset)*x*(1-z)*y+
    V011(offset)*x*z*(1-y)+
    V111(offset)*x*z*y;

  /*printf("returning: %f\n",v);*/
  return v;
}

float lookupInParametricSpace(float u, float v, int frame)
{
  float x,y,z;

  if(frame<fvframes) {
    x = u*fux + v*fvx + fcx;
    y = u*fuy + v*fvy + fcy;
    z = u*fuz + v*fvz + fcz;
    if(x >= fvwidth || y >= fvslices || z >= fvheight || x< 0 || y < 0 || z < 0)
      return NO_VALUE;
    
    if(do_interpolate==1)
      return sampleData(x,y,z,frame); 
    else
      return rawFVData[(int)floor(y)][(int)floor(fvheight-z)*fvwidth+(int)floor(x)+frame*fvwidth*fvheight];
      
  } else {
    return NO_VALUE;
  }
}


float lookupInVoxelSpace(float x, float y, float z, int frame)
{
  float x1,y1,z1;

  if(frame<fvframes) {
    if(hot)
      printf("lookup got: %f, %f, %f\n",x,y,z);
    transformFV(x,y,z,&x1,&y1,&z1);
    if(hot)
      printf("resulting functional coord is: %f %f %f\n",x1,y1,z1);
    
    if(x1 >= fvwidth || x1 < 0 || y1 >=fvslices || y1 < 0 || z1 >= fvheight || z1 <= 0)
      return NO_VALUE ;
    
    if(do_interpolate == 1)
      return sampleData(x1,y1,z1,overlay_frame);
    else
      return rawFVData[(int)floor(y1)][(int)floor(63-z1)*fvwidth+(int)floor(x1)+frame*fvwidth*fvheight];
  }
  return(NO_VALUE) ; 
}
  
void printOutFunctionalCoordinate(float x, float y, float z)
{
  
  float x2,y2,z2;
  float x1,y1,z1;
  /*printf("got coord: %f, %f, %f\n",x,y,z);*/
  x1 = x;
  y1 = y;
  z1 = z;
  if(overlay_frame < fvframes) {
    /*printf("pretransform: %f, %f, %f\n",x1,y1,z1);
      hot = 1;*/
    transformFV(x1,y1,z1,&x2,&y2,&z2);
    hot = 0;
    /*printf("resulting functional coordinate: %f, %f, %f\n",x2,y2,z2);
      printf("resulting voxel coordinate is: (%f, %f, %f)\n",y2,63-z2,x2);*/
    printf("p-value: %f\n",
     (x2 >= fvwidth || x2 < 0 || y2 >=fvslices || y2 < 0 
      || z2 >= fvheight || z2 <= 0)?NO_VALUE:rawFVData[(int)floor(y2)]
     [(int)floor(63-z2)*fvwidth+(int)floor(x2)+overlay_frame*fvwidth*fvheight]); 
  }
}
#endif 
/*--------------------- drawing hacks -------------------------------*/
#ifndef HACK
#define HACK

/* NOTE: In fischl space slices are coded in z, in the fspace in y, so you
   have to swap z & y before !!! transforming
*/

/* static variables *wuerg* */
#define COLOR_WHEEL         0   /* complexval */
#define HEAT_SCALE          1   /* stat,positive,"" */
#define CYAN_TO_RED         2   /* stat,positive,"" */
#define BLU_GRE_RED         3   /* stat,positive,"" */
#define TWOCOND_GREEN_RED   4   /* complexval */
#define JUST_GRAY           5   /* stat,positive,"" */
#define BLUE_TO_RED_SIGNED  6   /* signed */
#define GREEN_TO_RED_SIGNED 7   /* signed */
#define RYGB_WHEEL          8   /* complexval */
#define NOT_HERE_SIGNED     9   /* signed */

float ux, uy, uz;
float vx, vy, vz;
float cx, cy, cz;

unsigned char* dhcache;
float* fcache;

int colscale=1;

/* setup the span vectors */

void initCache()
{
  int i;
  dhcache=(unsigned char*)malloc(512*512);
  fcache=(float*)malloc(512*512*sizeof(float));

  for(i=0; i<512*512; i++)
    fcache[i]=NO_VALUE;
}

void setupSpans()
{
  ux = 128;
  uy = 0;
  uz = 0;

  vx = 0;
  vy = 128;
  vz = 0;

  cx = 128;
  cy = 128;
  cz = 128;  
}

/* some default setups */

void setupCoronal(int slice)
{
  ux = 128;
  uy = 0;
  uz = 0;

  vx = 0;
  vy = 128;
  vz = 0;

  cx = 128;
  cy = 128;
  cz = slice;  

  if(statsloaded) {
    
    setupFVC(128-cx,128-(255-cz),cy-128);
    setupFVU(128-(cx+ux),128-(255-(cz+uz)),(cy+uy)-128);
    setupFVV(128-(cx+vx),128-(255-(cz+vz)),(cy+vy)-128);
  }
}

void setupSagittal(int slice)
{
  ux = 0;
  uy = 0;
  uz = 128;

  vx = 0;
  vy = 128;
  vz = 0;

  cx = slice;
  cy = 128;
  cz = 128;
  
  if(statsloaded) {
    setupFVC(128-cx,128-(255-cz),cy-128);
    setupFVU(128-(cx+ux),128-(255-(cz+uz)),(cy+uy)-128);
    setupFVV(128-(cx+vx),128-(255-(cz+vz)),(cy+vy)-128);
  }
   
}

void setupHorizontal(int slice)
{
  ux = 128;
  uy = 0;
  uz = 0;

  vx = 0;
  vy = 0;
  vz = 128;

  cx = 128;
  cy = slice;
  cz = 128;

  if(statsloaded) {
    setupFVC(128-cx,128-(255-cz),cy-128);
    setupFVU(128-(cx+ux),128-(255-(cz+uz)),(cy+uy)-128);
    setupFVV(128-(cx+vx),128-(255-(cz+vz)),(cy+vy)-128);
  }
}

/* getPlane 
   TODO: 
   * change this into a midpoint algorithm using the spatial 
     coherency of a plane
   * avoid a new setup for every scanline
   * think over integrating texturescaling for sgis
*/

void getPlane(char* fbuffer, int zf, int xoff, int yoff)
{
  float x, y, z;
  int w,h;
  float u,v;
  float step;
  float ostep;
  int ozf;
  int myoff;

  register float sux, svx, aux, avx;
  register float suy, svy, auy, avy;
  register float suz, svz, auz, avz;

  ozf = zf;
  myoff = 512*yoff+xoff;

#ifdef Linux
  if(zf==2) 
    zf = 1;
#endif

  step = 1.0/(128.0*zf);
  ostep = 1.0/(128.0*ozf);

  sux = step*ux; suy = step*uy; suz = step*uz;
  svx = step*vx; svy = step*vy; svz = step*vz;

  aux = -ux; auy = -uy; auz = -uz;
  avx = -vx; avy = -vy; avz = -vz;
  x = aux + avx + cx;
  y = auy + avy + cy; 
  z = auz + avz + cz;

  for(h=0; h<256*zf; ++h) {
    x = aux + avx + cx;
    y = auy + avy + cy; 
    z = auz + avz + cz;
    for(w=0; w<256*zf; ++w) {
      if(z<0 || z > 255 || x <0 || x>255 || y<0 || y>255) {
  dhcache[myoff+w]=0;
  continue;
      }
      /* this sucks on a sgi !!!
   floating point calculation is to slow on R5000 and
   this cast during memory access kills the remaining speed
      */
      dhcache[myoff+w] = im[(int)z][(int)(256-1-y)][(int)x];
      x += sux;
      y += suy; 
      z += suz;
    }
    avx += svx;
    avy += svy;
    avz += svz;

    aux = -ux;
    auy = -uy;
    auz = -uz;
    myoff += 512;
  }

  myoff = 512*yoff+xoff;
  u = v = -1;

  if(statsloaded && do_overlay==1){ 
    for(h=0; h <256*ozf; ++h) {
      for(w=0; w<256*ozf; ++w) {
  
  fcache[myoff+w] = lookupInParametricSpace(u,v,overlay_frame); 
  u += ostep; 
      }
      u = -1;
      v += ostep;
      myoff += 512;
    }
    
  } else if (!statsloaded && do_overlay==1) {
    loadFV();
  }
  
#ifdef Linux
  if(ozf==2)
    scale2x(512,512,dhcache);
#endif
}

void getMaximumProjection(int zf, int xoff, int yoff)
{
  int w, h;
  
  for(h=0; h<256; h++) {
    for(w=0; w<256; w++) {
      dhcache[(h+yoff)*512+w+xoff] = sim[2][255-h][w];
    }
  }
}

void setStatColor(float f, unsigned char *rp, unsigned char *gp, unsigned char *bp,
      float tmpoffset)
{

  float r,g,b;
  float ftmp,c1,c2;

  if (fabs(f)>f2thresh && fabs(f)<fmid)
  {
    ftmp = fabs(f);
    c1 = 1.0/(fmid-f2thresh);
    c2 = 1.0;

    ftmp = c2*(ftmp-f2thresh)+f2thresh;
    f = (f<0)?-ftmp:ftmp;
  }

  if (colscale==HEAT_SCALE)
  {
    if (f>=0)
    {
      r = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<f2thresh)?0:(f<fmid)?(f-f2thresh)/(fmid-f2thresh):1);
      g = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      b = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0);
    } else
    {
      f = -f;
      b = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<f2thresh)?0:(f<fmid)?(f-f2thresh)/(fmid-f2thresh):1);
      g = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      r = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0);
    }
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==BLU_GRE_RED)
  {  
    if (f>=0)
    {
      r = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<f2thresh)?0:(f<fmid)?(f-f2thresh)/(fmid-f2thresh):1);
      g = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      b = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
    } else
    {
      f = -f;
      b = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<f2thresh)?0:(f<fmid)?(f-f2thresh)/(fmid-f2thresh):1);
      g = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      r = tmpoffset*((f<f2thresh)?1:(f<fmid)?1-(f-f2thresh)/(fmid-f2thresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
    }
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==JUST_GRAY)
    {
      if (f<0) f = -f;
      r = g = b = f*255;
    }
  *rp = (unsigned char)r;
  *gp = (unsigned char)g;
  *bp = (unsigned char)b;
}

void compose(unsigned char* stbuffer, unsigned char* outbuffer)
{
  int w,h;
  int hax,hay;
  int i,j,k;
  int curs;
  int imnr;

  unsigned char red,green,blue;
  unsigned char ccolor1,ccolor2;

  if(all3flag) {
    ccolor1 = 255;
    ccolor2 = 255;
  } else {
    ccolor1 = NUMVALS-1;
    ccolor2 = 0;
  }

  hax = xdim/2;
  hay = ydim/2;
  
  if(all3flag)
    curs = 15;
  else 
    curs = 2;

  for(h=0; h<ydim; h++) {
    for(w=0; w<xdim; w++) {
      if(do_overlay == 1) { 
  if(fcache[h*512+w]!=NO_VALUE)
    setStatColor(fcache[h*512+w],&red,&green,&blue,stbuffer[h*512+w]/255.0);
  else {
    red = green = blue = stbuffer[h*512+w];
  }
      } else {
  red = green = blue = stbuffer[h*512+w];
      }
      outbuffer[4*h*xdim+4*w]=red;
      outbuffer[4*h*xdim+4*w+1]=green;
      outbuffer[4*h*xdim+4*w+2]=blue;
      outbuffer[4*h*xdim+4*w+3]=255;
    }
  }

  if(all3flag || plane==SAGITTAL) {
    for (i=ic-curs;i<=ic+curs;i++) {
      if (all3flag) 
  k = 4*(i/2*xdim/2+imc/2 + i/2*hax);
      else if(plane==SAGITTAL)        
  k = 4*(i*xdim+imc);
      outbuffer[k] = ccolor1 ; outbuffer[k+1] = ccolor2;
      outbuffer[k+2] = ccolor2;
      outbuffer[k+3]=255;
    }
    for (imnr=imc-curs;imnr<=imc+curs;imnr++) {
      if (all3flag) k = 4*(ic/2*xdim/2+imnr/2 + ic/2*hax);
      else if(plane==SAGITTAL)
  k = 4*(ic*xdim+imnr);
      outbuffer[k] = ccolor1; 
      outbuffer[k+1] = ccolor2;
      outbuffer[k+2] = ccolor2;
      outbuffer[k+3]=255;
    }
  }
  if(all3flag || plane==CORONAL) {
    for (i=ic-curs;i<=ic+curs;i++) {
      if (all3flag) k = 4*(xdim*hay + i/2*xdim/2+jc/2 + i/2*hax);
      else if(plane==CORONAL)      
  k = 4*(i*xdim+jc);
      vidbuf[k] = ccolor1; vidbuf[k+1] = vidbuf[k+2] = ccolor2; 
      vidbuf[k+3]= 255;
    }
    for (j=jc-curs;j<=jc+curs;j++) {
      if (all3flag) k = 4*(xdim*hay + ic/2*xdim/2+j/2 + ic/2*hax);
      else if(plane==CORONAL)
  k = 4*(ic*xdim+j);
      vidbuf[k] = ccolor1; vidbuf[k+1] = vidbuf[k+2] = ccolor2; 
      vidbuf[k+3]= 255;
    }
  }
  if(all3flag || plane==HORIZONTAL) {
    for (imnr=imc-curs;imnr<=imc+curs;imnr++) {
      if (all3flag) k = 4*(xdim*hay+hax + imnr/2*xdim/2+jc/2 + imnr/2*hax);
      else if(plane == HORIZONTAL)
  k = 4*(imnr*xdim+jc);
      vidbuf[k] = ccolor1; vidbuf[k+1] = vidbuf[k+2] = ccolor2;
      vidbuf[k+3]=255;
    }
    for (j=jc-curs;j<=jc+curs;j++) {
      if (all3flag) k = 4*(xdim*hay + hax + imc/2*xdim/2+j/2 + imc/2*hax);
      else if(plane == HORIZONTAL)
  k = 4*(imc*xdim+j);
      vidbuf[k] = ccolor1; vidbuf[k+1] = vidbuf[k+2] = ccolor2;
      vidbuf[k+3]=255;
    }
  }
}

void draw_image_hacked(int imc, int ic, int jc)
{
  memset(fcache,0,4*512*512);

  memset(dhcache,128,512*512);

  if(!all3flag) {
    switch(plane) {
    case CORONAL:
      setupCoronal(imc/zf);
      break;
    case SAGITTAL:
      setupSagittal(jc/zf);
      break;
    case HORIZONTAL:
      setupHorizontal(ic/zf);
      break;
    }
    getPlane(NULL,2,0,0);
  } else if(all3flag) {
    setupCoronal(imc/zf);
    getPlane(NULL,1,0,256);
    setupSagittal(jc/zf);
    getPlane(NULL,1,0,0);
    setupHorizontal(ic/zf);
    getPlane(NULL,1,256,256);
    getMaximumProjection(1,256,0);
  }
  
  compose(dhcache,vidbuf);
  rectwrite(0,0,xdim-1,ydim-1,vidbuf);

  /* rectwrite(0,0,xdim-1,ydim-1,dhcache);*/
}

#endif
/*--------------------- twitzels hack end ---------------------------*/

#ifdef TCL
  int Medit(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
#else
void
main(int argc,char *argv[])
#endif
{
    int i,j;
    int tclscriptflag=FALSE;
    char fname[NAME_LENGTH];
    char *lsubjectsdir;
    FILE *fp;
    char lpname[NAME_LENGTH], limdir[NAME_LENGTH], lsrname[NAME_LENGTH];
    char *getenv(), lsurface[NAME_LENGTH];
    /*%%%Sean %%%%*/
    /* char cwd[NAME_LENGTH]; */
    char *cwd; 
    char *word;
    char path[MAX_DIR_DEPTH][NAME_LENGTH];

    if (argc<2) {
printf("\n");
printf("Usage: %s name imagedir [surfname] [-] [-tcl script]\n",argv[0]);
printf("       %s local [-] [-tcl script]\n",argv[0]);
printf("\n");
printf("       name:  subjectdir name              (relative to $SUBJECTS_DIR)\n");
printf("   imagedir:  orig,T1,brain,wm,filled      (relative to $subject/mri)\n");
printf("   surfname:  ?h.orig,?h.smoothwm,?h.plump (relative to $subject/surf)\n");
printf("  lone dash:  disable editing\n");
#ifdef OPENGL
printf("                                         [vers: 272134--OpenGL]\n");
#else
printf("                                         [vers: 272134--GL]\n");
#endif
printf("\n");
exit(1);
    }
    zf = ozf = 2;
    fsf = (float)zf;
    surflinewidth = zf;

    if (MATCH(argv[argc-2],"-tcl")) {
      argc -= 2;   /* strip so like before */
      tclscriptflag = TRUE;
    }
    statsloaded = 0;

    lsubjectsdir = getenv("SUBJECTS_DIR");
    if (lsubjectsdir==NULL) {
      printf("medit: env var SUBJECTS_DIR undefined (use setenv)\n");exit(1);}

    strcpy(lpname,argv[1]);
    if (MATCH(lpname,"local"))  strcpy(limdir,"local");
    else if (MATCH(lpname,".")) strcpy(limdir,"local");
    else if (argc>2)            strcpy(limdir,argv[2]);
    else  { printf("medit: ### imagedir missing\n"); exit(1); }

    sprintf(fname,"%s/%s",lsubjectsdir,lpname);
    if (!MATCH(lpname,"local")) {
      if ((fp=fopen(fname,"r"))==NULL) {
        printf("medit: ### can't find subject %s\n",lpname); exit(1);}
      else fclose(fp);
    }

    if (argv[argc-1][0]=='-') {
      editflag = FALSE;
      printf("medit: ### editing disabled\n");
    }
    if ((editflag && argc>3) || (!editflag && argc>4)) {
      surfflag = TRUE;
      strcpy(lsurface,argv[3]);
    }
    else {
      strcpy(lsurface,"rh.orig");
    }
    if (argc>2 && MATCH(limdir,"local") && editflag)
      printf("medit: ### ignored bad arg after \"local\": %s\n",argv[2]);

    /* parse cwd to set session root (guess path to bem) */
    /* %%%%Sean %%%% */
    /* getwd(cwd); */
#ifdef Linux
    cwd = getcwd(NULL,0);
#else
    cwd = getenv("PWD");
#endif

    word = strtok(cwd,"/");
    strcpy(path[0],word);
    i = 1;
    while ((word = strtok(NULL,"/")) != NULL) {  /* save,count */
      strcpy(path[i],word);
      i++;
    }
    if (MATCH(path[i-1],"scripts") && 
        MATCH(path[i-2],lpname)) {
      printf("medit: in subjects \"scripts\" dir\n");
      j = i-1;
    } else if (MATCH(path[i-1],"scripts") &&
               MATCH(path[i-2],"eegmeg") &&
               MATCH(path[i-3],lpname)) {
      printf("medit: in subjects \"eegmeg/scripts\" dir\n");
      j = i-1;
    } else if (MATCH(path[i-1],"scripts")) {
      printf("medit: in \"scripts\" dir (not subjects,eegmeg)\n");
      j = i-1;
    } else if (MATCH(limdir,"local")) {  /* local even if in mri */
      printf("medit: local data set => session root is cwd\n");
      j = i;
    } else if (MATCH(path[i-2],"mri")) {
      printf("medit: in subdir of \"mri\" dir\n");
      j = i-1;
    } else {
      printf(
  "medit: not in \"scripts\" dir or \"mri\" subdir => session root is cwd\n");
      j = i;
    }
    sprintf(lsrname,"/%s",path[0]);  /* reassemble absolute */
    for(i=1;i<j;i++)
      sprintf(lsrname,"%s/%s",lsrname,path[i]);
    printf("medit: session root data dir ($session) set to:\n");
    printf("medit:     %s\n",lsrname);

    make_filenames(lsubjectsdir,lsrname,lpname,limdir,lsurface);

    if (surfflag)
      read_binary_surface(sfname);
    read_images(mfname);  /* sets xnum/dim,ynum/dim */

    vidbuf = (GLubyte*)lcalloc((size_t)xdim*ydim*4,(size_t)sizeof(GLubyte)); 
    /* vidbuf = (GLubyte*)lcalloc((size_t)xdim*ydim,(size_t)sizeof(GLubyte)); */
   /* vidbuf = (Colorindex *)lcalloc((size_t)xdim*ydim,(size_t)sizeof(Colorindex));  */
   cvidbuf = (Colorindex *)lcalloc((size_t)CVIDBUF,(size_t)sizeof(Colorindex));
   binbuff = (unsigned char *)lcalloc((size_t)3*xdim*ydim,(size_t)sizeof(char));

    imc = zf*imnr1/2;
    ic = ydim/2;
    jc = xdim/2;
   
    /* for twitzels stuff */
    statsloaded = 0;
    initCache();
    setupFOVars();

    /*loadFV("970121NH_02986_00003_00009_00001_001");*/
    /*loadFV("twsyn");*/
    /*overlay_file = "Sel_COND_6_0_diff";
      loadFV();*/
    /*loadFV("Sel_KP");*/
    /*loadFV("Sel_nomask_4_1_diff");*/
   
    /* end twitzels stuff */

    if (tclscriptflag) {
      /* called from tkmedit.c; do nothing (don't even open gl window) */
      /* wait for tcl interp to start; tkanalyse calls tcl script */
    }
    else {   /* open window for medit or non-script tkmedit */
      if (open_window(pname) < 0) {
  exit(1);
      }
      redraw();
    }

#ifndef TCL  /* tcl: omit event loop */
    if (tclscriptflag) {
      printf("medit: can't read tcl script  ...ignored (use tkmedit)\n");
    }
#  ifdef OPENGL
    printf("medit: ### non-tk medit event loop not converted for OpenGL\n");
    exit(1);
#  else
    while(1)
    {
        dev = qread(&val);  /* blocks here for next event */
        switch(dev)  {
          case ESCKEY:
            exit(0);
            break;
          case LEFTMOUSE:
            if (val == 0)  break;
            sx = getvaluator(MOUSEX);
            sy = getvaluator(MOUSEY);
            select_pixel(sx,sy,TRUE);
            redraw();
            break;
          case MIDDLEMOUSE:
            if (val == 0)  break;
            while (getbutton(MIDDLEMOUSE))
            {
              sx = getvaluator(MOUSEX);
              sy = getvaluator(MOUSEY);
              select_pixel(sx,sy,FALSE);
              for (k= -prad;k<=prad;k++)
              for (i= -prad;i<=prad;i++)
              for (j= -prad;j<=prad;j++)
              if (imc/zf+k>=0 && imc/zf+k<=numimg &&
                  (ydim-1-ic)/zf+i>=0 && (ydim-1-ic)/zf+i<ynum &&
                  jc/zf+j>=0 && jc/zf+j<xnum)
              {
                if (im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] < WM_MIN_VAL)
                  im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] = 255;
              }
              for (k= -prad;k<=prad;k++)
              if (imc/zf+k>=0 && imc/zf+k<=numimg)
              {
                changed[imc/zf+k] = TRUE;
                editedimage = imnr0+imc/zf+k;
              }
            }
            redraw();
            break;
          case RIGHTMOUSE:
            if (val == 0)  break;
            while (getbutton(RIGHTMOUSE))
            {
              sx = getvaluator(MOUSEX);
              sy = getvaluator(MOUSEY);
              select_pixel(sx,sy,FALSE);
              for (k= -prad;k<=prad;k++)
              for (i= -prad;i<=prad;i++)
              for (j= -prad;j<=prad;j++)
              if (imc/zf+k>=0 && imc/zf+k<=numimg &&
                  (ydim-1-ic)/zf+i>=0 && (ydim-1-ic)/zf+i<ynum &&
                  jc/zf+j>=0 && jc/zf+j<xnum)
              {
/*
                im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] = 0;
*/
/*
  AKL Erased pixels set to 1 instead of 0
      Leave pixels that are already 0, set to 0 when you try to erase
*/
                im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] = WM_EDITED_OFF;
              }
              for (k= -prad;k<=prad;k++)
              if (imc/zf+k>=0 && imc/zf+k<=numimg)
              {
                changed[imc/zf+k] = TRUE;
                editedimage = imnr0+imc/zf+k;
              }
            }
            redraw();
            break;
          case REDRAW:
            reshapeviewport();
            getsize(&lxdim,&lydim);
            xdim = (int)lxdim;
            ydim = (int)lydim;
            resize_buffers(xdim,ydim);
            zf = xdim/xnum;
            fsf = (float)zf;
            imc = (zf*imc)/ozf;
            ic = (zf*ic)/ozf;
            jc = (zf*jc)/ozf;
            redraw();
            break;
          case LEFTARROWKEY: case DOWNARROWKEY:
            if (val == 0)  break;
            if (plane==CORONAL)
              imc = (imc<zf)?imnr1*zf-zf+imc:imc-zf;
            else if (plane==HORIZONTAL)
              ic = (ic<zf)?ydim-zf+ic:ic-zf;
            else if (plane==SAGITTAL)
              jc = (jc<zf)?xdim-zf+jc:jc-zf;
      redraw();
            break;
          case RIGHTARROWKEY: case UPARROWKEY:
            if (val == 0)  break;
            if (plane==CORONAL)
              imc = (imc>=imnr1*zf-zf)?imc+zf-imnr1*zf:imc+zf;
            else if (plane==HORIZONTAL)
              ic = (ic>=ydim-zf)?ic+zf-ydim:ic+zf; 
            else if (plane==SAGITTAL)
              jc = (jc>=xdim-zf)?jc+zf-xdim:jc+zf;
      redraw();
            break;
          case KEYBD:
            switch((char)val)
            {
              case 'D': 
                        printf("enter dipole spacing (in voxels): ");
                        scanf("%d",&dip_spacing);
                        sprintf(dipfname,"%s/%s/bem/brain3d%d.dip",
                                subjectsdir,pname,dip_spacing);
                        sprintf(decfname,"%s/%s/bem/brain3d%d-%d.dec",
                                subjectsdir,pname,dip_spacing,dip_spacing);
                        write_dipoles(dipfname);
                        write_decimation(decfname);
                        break;
              case 't':
                        fthresh = 0.95;
                        break;
              case 'x': case 'X': plane=SAGITTAL;redraw();break;
              case 'y': case 'Y': plane=HORIZONTAL;redraw();break;
              case 'z': case 'Z': plane=CORONAL;redraw();break;
              case 'w': write_images(mfname);break;
              case 'b': prad = 0;printf("brush size = %d\n",prad*2+1);break;
              case 'B': prad=(prad+1)%5;
                        printf("brush size = %d\n",prad*2+1);
      break;
              case 'o': optimize(10);redraw();break;
              case 'O': optimize2();redraw();break;
              case 'M': maxflag=!maxflag;redraw();break;
              case 'm': maxflag=!maxflag;redraw();break;
              case 'I': mirror();redraw();break;
              case 'd': surfflag = !surfflag;break;
              case '*': fsquash *= 1.1; break;
              case '/': fsquash /= 1.1; break;
              case '!': fsquash = 5.0; break;
              case '+': fthresh += 0.05;
                        break;
              case '-': fthresh -= 0.05;
                        break;
              case '0': fthresh = 0.5;
                        break;
              case 'r': goto_point(tfname);
                        redraw();
                        break;
              case 'f': write_point(tfname);
                        break;
              case 'R': 
                        sprintf(fname,"%s","../bem/head2mri.hpts");
                        if ((fp=fopen(fname,"r"))!=NULL)
                        {
                          ptsflag = TRUE;
                          fgets(line,NAME_LENGTH,fp);
                          l = 0;
                          while (!feof(fp))
                          {
                            jold[l] = iold[l] = imold[l] = 0;
                            if (line[0]=='2') /* .show format? */
                              sscanf(line,"%*s%*s%*s%*s%*s%*s%*s%f%*s%f%*s%f",
                                     &ptx[l],&pty[l],&ptz[l]);
                            else /* .apos format */
                            {
                              sscanf(line,"%*s%*s%f%f%f",&ptx[l],&pty[l],&ptz[l]);
                              ptx[l] *= 1000; /* convert from m to mm */
                              pty[l] *= 1000;
                              ptz[l] *= 1000;
                            }
/*
                            printf("%d: %f %f %f\n",l,ptx[l],pty[l],ptz[l]);
*/
                            fgets(line,NAME_LENGTH,fp);
                            l++;
                          }
                          printf("head points file %s read\n",fname);
                          npts = l;
                          fclose(fp);
                          sprintf(fname,"%s","../bem/head2mri.trans");
                          if ((fp=fopen(fname,"r"))!=NULL)
                          {
                            for (i=0;i<4;i++)
                            for (j=0;j<4;j++)
                              fscanf(fp,"%f",&tm[i][j]);
                            printf("transformation file %s read\n",fname);
                          } else
                          {
                            for (i=0;i<4;i++)
                            for (j=0;j<4;j++)
                              tm[i][j] = (i==j);
                          }
                          for (l=0;l<MAXPARS;l++)
                          {
                            par[l] = dpar[l] = 0;
                          }
                        } else
                          printf("file %s not found\n",fname);
                  redraw();
                        break;
              case 'W':
                        sprintf(fname,"%s","../bem/head2mri.trans");
                        if ((fp=fopen(fname,"w"))==NULL)
                        {
                          printf("can't open file %s\n",fname);
                        } else
                        {
                          for (i=0;i<4;i++)
                          {
                            for (j=0;j<4;j++)
                              fprintf(fp,"%13.6e ",tm[i][j]);
                            fprintf(fp,"\n");
                          }
                          fclose(fp);
                        }
                        printf("transformation file %s written\n",fname);
                        break;
              /* rot z ccw */
              case '[': if (plane==SAGITTAL) rotate_brain(5.0,'x');
                        if (plane==HORIZONTAL) rotate_brain(-5.0,'z');
                        if (plane==CORONAL) rotate_brain(5.0,'y');
                        redraw();
                        break;
              case '{': if (plane==SAGITTAL) rotate_brain(50.0,'x');
                        if (plane==HORIZONTAL) rotate_brain(-50.0,'z');
                        if (plane==CORONAL) rotate_brain(50.0,'y');
                        redraw();
                        break;
              /* rot z cw */
              case ']': if (plane==SAGITTAL) rotate_brain(-5.0,'x');
                        if (plane==HORIZONTAL) rotate_brain(5.0,'z');
                        if (plane==CORONAL) rotate_brain(-5.0,'y');
                        redraw();
                        break;
              case '}': if (plane==SAGITTAL) rotate_brain(-50.0,'x');
                        if (plane==HORIZONTAL) rotate_brain(50.0,'z');
                        if (plane==CORONAL) rotate_brain(-50.0,'y');
                        redraw();
                        break;
              /* trans +y */
              case 'p': if (plane==SAGITTAL) translate_brain(0.5,'z');
                        if (plane==HORIZONTAL) translate_brain(0.5,'y');
                        if (plane==CORONAL) translate_brain(0.5,'z');
                        redraw();
                        break;
              case 'P': if (plane==SAGITTAL) translate_brain(2.5,'z');
                        if (plane==HORIZONTAL) translate_brain(2.5,'y');
                        if (plane==CORONAL) translate_brain(2.5,'z');
                        redraw();
                        break;
              /* trans -y */
              case '.': if (plane==SAGITTAL) translate_brain(-0.5,'z');
                        if (plane==HORIZONTAL) translate_brain(-0.5,'y');
                        if (plane==CORONAL) translate_brain(-0.5,'z');
                        redraw();
                        break;
              case '>': if (plane==SAGITTAL) translate_brain(-2.5,'z');
                        if (plane==HORIZONTAL) translate_brain(-2.5,'y');
                        if (plane==CORONAL) translate_brain(-2.5,'z');
                        redraw();
                        break;
              /* trans -x */
              case 'l': if (plane==SAGITTAL) translate_brain(-0.5,'y');
                        if (plane==HORIZONTAL) translate_brain(0.5,'x');
                        if (plane==CORONAL) translate_brain(0.5,'x');
                        redraw();
                        break;
              case 'L': if (plane==SAGITTAL) translate_brain(-2.5,'y');
                        if (plane==HORIZONTAL) translate_brain(2.5,'x');
                        if (plane==CORONAL) translate_brain(2.5,'x');
                        redraw();
                        break;
              /* trans +x */
              case ';': if (plane==SAGITTAL) translate_brain(0.5,'y');
                        if (plane==HORIZONTAL) translate_brain(-0.5,'x');
                        if (plane==CORONAL) translate_brain(-0.5,'x');
                        redraw();
                        break;
              case ':': if (plane==SAGITTAL) translate_brain(2.5,'y');
                        if (plane==HORIZONTAL) translate_brain(-2.5,'x');
                        if (plane==CORONAL) translate_brain(-2.5,'x');
                        redraw();
                        break;
            }
            break;
          }  /* dev switch */
    set_scale();
    }  /* while events */
#  endif
#endif  /* tcl: omit event loop */
return(0);
}

void
do_one_gl_event(Tcl_Interp *interp)   /* tcl */
{

#ifdef OPENGL     /* derived from event.c:DoNextEvent(),tkExec() */
  XEvent current, ahead;
  char buf[1000];
  char command[NAME_LENGTH];
  KeySym ks;
  static int ctrlkeypressed = FALSE;
  static int altkeypressed = FALSE;
  static int shiftkeypressed = FALSE;
  static int button1pressed = FALSE;
  static int button2pressed = FALSE;
  static int button3pressed = FALSE;
  short sx,sy; /* Screencoord sx,sy; */
  XWindowAttributes wat;
  Window junkwin;
  int rx, ry;

  if (!openglwindowflag) return;
  if (updatepixval) {
    Tcl_Eval(interp,"pixvaltitle 1 1 1");
    /*    Tcl_Eval(interp,"set selectedpixval $selectedpixval"); *//*touch for trace*/
    updatepixval = FALSE;
  }
  if (XPending(xDisplay)) {  /* do one if queue test */
  
    XNextEvent(xDisplay, &current);   /* blocks here if no event */

    switch (current.type) {

      case ConfigureNotify:
        XGetWindowAttributes(xDisplay, w.wMain, &wat);
        XTranslateCoordinates(xDisplay, w.wMain, wat.root,
                              -wat.border_width, -wat.border_width,
                              &rx, &ry, &junkwin);
        w.x = rx;  /* left orig         (current.xconfigure.x = 0 relative!) */
        w.y = ry;  /* top orig: +y down (current.xconfigure.y = 0 relative!) */
        w.w = current.xconfigure.width;
        w.h = current.xconfigure.height;
        resize_window_intstep(0);
        /* Tcl_Eval(interp,"set zf $zf"); */ /* touch for trace */
        if (followglwinflag && zf != 4) {
          sprintf(command,"wm geometry . +%d+%d",
                      w.x, w.y + w.h + MOTIF_YFUDGE /*+MOTIF_XFUDGE*/);
          Tcl_Eval(interp,command);
          /*Tcl_Eval(interp,"raise .");*/
        }
        break;

      case Expose:
        if (XPending(xDisplay)) {
          XPeekEvent(xDisplay, &ahead);
          if (ahead.type==Expose) break;  /* skip extras */
        }
        redraw();
        break;

      case ButtonPress:
        sx = current.xbutton.x;
        sy = current.xbutton.y;
        sx += w.x;   /* convert back to screen pos (ugh) */
        sy = 1024 - w.y - sy;
        if (current.xbutton.button == 1) {  /** left **/
          select_pixel(sx,sy,TRUE);
          Tcl_Eval(interp,"unzoomcoords; sendupdate");
          button1pressed = TRUE;
        }
        if (current.xbutton.button == 2) {  /** middle **/
          if (drawsecondflag) {
            printf("medit: N.B.: editing image in hidden buffer\n"); PR}
          select_pixel(sx,sy,FALSE);
          edit_pixel(TO_WHITE);
          button2pressed = TRUE;
        }
        if (current.xbutton.button == 3) {  /** right **/
          if (drawsecondflag) {
            printf("medit: N.B.: editing image in hidden buffer\n"); PR}
          select_pixel(sx,sy,FALSE);
          edit_pixel(TO_BLACK);
          button3pressed = TRUE;
        }
        break;

      case MotionNotify:
        sx = current.xmotion.x;
        sy = current.xmotion.y;
        sx += w.x;   /* convert back to screen pos (ugh) */
        sy = 1024 - w.y - sy;
        if (button2pressed) {select_pixel(sx,sy,FALSE); edit_pixel(TO_WHITE);}
        if (button3pressed) {select_pixel(sx,sy,FALSE); edit_pixel(TO_BLACK);}
        break;

      case ButtonRelease:
        if (current.xbutton.button == 1)  button1pressed = FALSE;
        if (current.xbutton.button == 2)  button2pressed = FALSE;
        if (current.xbutton.button == 3)  button3pressed = FALSE;
        redraw();
        break;

      case KeyPress:
  XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
        switch (ks) {

          case XK_c: case XK_0: case XK_apostrophe:
            if (altkeypressed) Tcl_Eval(interp,
              ".mri.main.right.cmp.aCOMPARE.bu invoke");
            break;
          case XK_M:
            if (altkeypressed) Tcl_Eval(interp,
              ".mri.main.left.view.butt.left.amaxflag.ck invoke");
            break;
          case XK_m:
            if (altkeypressed) Tcl_Eval(interp,
              ".mri.main.left.view.butt.left.amaxflag.ck invoke");
            break;
          case XK_d:
            if (altkeypressed) Tcl_Eval(interp,
              ".mri.main.left.view.butt.left.asurfflag.ck invoke");
            break;
          case XK_b:
            if (altkeypressed) Tcl_Eval(interp,
              "set prad 0");
            break;
          case XK_v:
            if (altkeypressed) Tcl_Eval(interp,
              "set pradtmp $prad; set prad $pradlast; set pradlast $pradtmp");
            break;
          case XK_B:
            if (altkeypressed) Tcl_Eval(interp,
              "set prad [expr ($prad+1)]");
            break;
          case XK_asterisk:
            if (altkeypressed) Tcl_Eval(interp,
              "set fsquash [expr $fsquash * 1.1]; set_scale");
            break;
          case XK_slash:
            if (altkeypressed) Tcl_Eval(interp,
              "set fsquash [expr $fsquash / 1.1]; set_scale");
            break;
          case XK_plus:
            if (altkeypressed) Tcl_Eval(interp,
              "set fthresh [expr $fthresh + 0.05]; set_scale");
            break;
          case XK_minus:
            if (altkeypressed) Tcl_Eval(interp,
              "set fthresh [expr $fthresh - 0.05]; set_scale");
            break;
          case XK_w:
            if (altkeypressed) Tcl_Eval(interp,
              ".mri.main.left.head.save.aSAVEIMG.bu invoke");
            break;
          case XK_x: case XK_X: plane=SAGITTAL;redraw();break;
          case XK_y: case XK_Y: plane=HORIZONTAL;redraw();break;
          case XK_z: case XK_Z: plane=CORONAL;redraw();break;
          case XK_r:
            if (altkeypressed) {
              /*Tcl_Eval(interp,"raise .");*/
              goto_point(tfname);
              redraw();
            }
            break;
          case XK_f:
            if (altkeypressed) write_point(tfname);
            break;
          case XK_I:
            if (altkeypressed) { mirror(); redraw();}
            break;
          case XK_R:
            if (altkeypressed) {
              read_hpts(hpfname);
              read_htrans(htfname);
              redraw();
            }
            break;
          case XK_W:
            if (altkeypressed) write_htrans(htfname);
            break;
          case XK_Up:
            Tcl_Eval(interp,"changeslice up; sendupdate");
            break;
          case XK_Down:
            Tcl_Eval(interp,"changeslice down; sendupdate");
            break;
          case XK_Right:
            Tcl_Eval(interp,"changeslice up; sendupdate");
            break;
          case XK_Left:
            Tcl_Eval(interp,"changeslice down; sendupdate");
            break;

          /* TODO: install X versions of rotate headpoints bindings */

          /* modifiers */
          case XK_Shift_L:   case XK_Shift_R:   shiftkeypressed=TRUE; break;
          case XK_Control_L: case XK_Control_R: ctrlkeypressed=TRUE;  break;
          case XK_Alt_L:     case XK_Alt_R:     altkeypressed=TRUE;   break;
        }
        break;

      case KeyRelease:   /* added this mask to xwindow.c */
  XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
        switch (ks) {
          case XK_Shift_L:   case XK_Shift_R:   shiftkeypressed=FALSE; break;
          case XK_Control_L: case XK_Control_R: ctrlkeypressed=FALSE;  break;
          case XK_Alt_L:     case XK_Alt_R:     altkeypressed=FALSE;   break;
        }
        break;

    }
    /* return GL_FALSE; */
  }

#else  /* use gl calls */
  short dev, val;
  static int ctrlkeypressed = FALSE;
  static int altkeypressed = FALSE;
  int i,j,k,l;
  int irad,jrad,krad;
  Screencoord sx,sy;
  long lxdim,lydim,xorig,yorig;
  char command[NAME_LENGTH];
  char fname[NAME_LENGTH], line[NAME_LENGTH];
  FILE *fp;


  if (!openglwindowflag) return;
  if (qtest()) {  /* do one event */
    dev = qread(&val);
    if (dev != LEFTMOUSE && dev != RIGHTMOUSE) /* hack: mouse zeros getbutton!*/
      ctrlkeypressed = getbutton(LEFTCTRLKEY) || getbutton(RIGHTCTRLKEY);
    switch(dev) {
        break;
      case REDRAW:
        reshapeviewport();
        getsize(&lxdim,&lydim);
        xdim = (int)lxdim;
        ydim = (int)lydim;
        resize_buffers(xdim,ydim);
        ozf = zf;
        zf = xdim/xnum;
        if (zf!=ozf && surflinewidth==ozf)
          surflinewidth=zf;
        fsf = (float)zf;
        imc = (zf*imc)/ozf;
        ic = (zf*ic)/ozf;
        jc = (zf*jc)/ozf;
        Tcl_Eval(interp,"set zf $zf");  /* touch for trace */
        redraw();
        getorigin(&xorig,&yorig);
        if (followglwinflag && zf != 4) {
          sprintf(command,"wm geometry . +%d+%d",
                      xorig, 1024 - yorig + MOTIF_YFUDGE /*+ MOTIF_XFUDGE*/);
          Tcl_Eval(interp,command);
          Tcl_Eval(interp,"raise .");
        }
        break;
      case LEFTARROWKEY: case DOWNARROWKEY:
        if (val == 0)  break;
        Tcl_Eval(interp,"changeslice down; sendupdate");
        break;
      case RIGHTARROWKEY: case UPARROWKEY:
        if (val == 0)  break;
        Tcl_Eval(interp,"changeslice up; sendupdate");
        break;
      case LEFTMOUSE:
        if (val == 0)  break;
        sx = getvaluator(MOUSEX);
        sy = getvaluator(MOUSEY);
        select_pixel(sx,sy,TRUE);
        Tcl_Eval(interp,"unzoomcoords; sendupdate");
        redraw();
        break;
      case MIDDLEMOUSE:
        if (val == 0)  break;
        if (drawsecondflag) {
          printf("medit: N.B.: editing image in hidden buffer\n"); PR}
        while (getbutton(MIDDLEMOUSE))
        {
          sx = getvaluator(MOUSEX);
          sy = getvaluator(MOUSEY);
          select_pixel(sx,sy,FALSE);
          jrad=irad=krad=prad;
          if (inplaneflag) { if(plane==CORONAL)    krad=0;
                             if(plane==HORIZONTAL) irad=0;
                             if(plane==SAGITTAL)   jrad=0; }
          for (k= -krad;k<=krad;k++)
          for (i= -irad;i<=irad;i++)
          for (j= -jrad;j<=jrad;j++)
          if (imc/zf+k>=0 && imc/zf+k<numimg &&
              (ydim-1-ic)/zf+i>=0 && (ydim-1-ic)/zf+i<ynum &&
              jc/zf+j>=0 && jc/zf+j<xnum)
          {
            if (circleflag && k*k+i*i+j*j>prad*prad) continue;
            if (im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] < WM_MIN_VAL)
              im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] = 255;
          }
          for (k= -prad;k<=prad;k++)
          if (imc/zf+k>=0 && imc/zf+k<numimg)
          {
            changed[imc/zf+k] = TRUE;
            editedimage = imnr0+imc/zf+k;
          }
        }
        redraw();
        break;
      case RIGHTMOUSE:
        if (val == 0)  break;
        if (drawsecondflag) {
          printf("medit: N.B.: editing image in hidden buffer\n"); PR}
        while (getbutton(RIGHTMOUSE))
        {
          sx = getvaluator(MOUSEX);
          sy = getvaluator(MOUSEY);
          select_pixel(sx,sy,FALSE);
          jrad=irad=krad=prad;
          if (inplaneflag) { if(plane==CORONAL)    krad=0;
                             if(plane==HORIZONTAL) irad=0;
                             if(plane==SAGITTAL)   jrad=0; }
          for (k= -krad;k<=krad;k++)
          for (i= -irad;i<=irad;i++)
          for (j= -jrad;j<=jrad;j++)
          if (imc/zf+k>=0 && imc/zf+k<numimg &&
              (ydim-1-ic)/zf+i>=0 && (ydim-1-ic)/zf+i<ynum &&
              jc/zf+j>=0 && jc/zf+j<xnum)
          {
            if (circleflag && k*k+i*i+j*j>prad*prad) continue;
            im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] = WM_EDITED_OFF;
          }
          for (k= -prad;k<=prad;k++)
          if (imc/zf+k>=0 && imc/zf+k<numimg)
          {
            changed[imc/zf+k] = TRUE;
            editedimage = imnr0+imc/zf+k;
          }
        }
        redraw();
        break;
      case LEFTALTKEY: case RIGHTALTKEY:
        if (val)  altkeypressed = TRUE;
        else      altkeypressed = FALSE;
        break;
      case KEYBD:
        if (altkeypressed)
        switch((char)val) {
          case 'c': case '0': case '\'':
                    Tcl_Eval(interp,".mri.main.right.cmp.aCOMPARE.bu invoke");
                    break;
          case 'M': Tcl_Eval(interp,
                           ".mri.main.left.view.butt.left.amaxflag.ck invoke");
                    break;
          case 'm': Tcl_Eval(interp,
                           ".mri.main.left.view.butt.left.amaxflag.ck invoke");
                    break;
          case 'd': Tcl_Eval(interp,
                           ".mri.main.left.view.butt.left.asurfflag.ck invoke");
                    break;
          case 'b': Tcl_Eval(interp,"set prad 0");
                    break;
          case 'B': Tcl_Eval(interp,"set prad [expr ($prad+1)]");
                    break;
          case 'v': Tcl_Eval(interp,
                  "set pradtmp $prad;set prad $pradlast;set pradlast $pradtmp");
                    break;
          case '*': Tcl_Eval(interp,
                           "set fsquash [expr $fsquash * 1.1]; set_scale");
                    break;
          case '/': Tcl_Eval(interp,
                           "set fsquash [expr $fsquash / 1.1]; set_scale");
                    break;
          case '+': Tcl_Eval(interp,
                           "set fthresh [expr $fthresh + 0.05]; set_scale");
                    break;
          case '-': Tcl_Eval(interp,
                           "set fthresh [expr $fthresh - 0.05]; set_scale");
                    break;
          case 'w': Tcl_Eval(interp,
                           ".mri.main.left.head.save.aSAVEIMG.bu invoke");
                    break;
          case 'x': case 'X': plane=SAGITTAL;redraw();break;
          case 'y': case 'Y': plane=HORIZONTAL;redraw();break;
          case 'z': case 'Z': plane=CORONAL;redraw();break;
          case 'r': /*Tcl_Eval(interp,"raise .");*/
                    goto_point(tfname);
                    redraw();
                    break;
          case 'f': write_point(tfname);
                    break;
          case 'I': mirror();redraw();break;
          case 'R': read_hpts(hpfname);
                    read_htrans(htfname);
                    redraw();
                    break;
          case 'W': write_htrans(htfname);
                    break;
          /* rot z ccw */
          case '[': if (plane==SAGITTAL) rotate_brain(5.0,'x');
                    if (plane==HORIZONTAL) rotate_brain(-5.0,'z');
                    if (plane==CORONAL) rotate_brain(5.0,'y');
                    redraw();
                    break;
          case '{': if (plane==SAGITTAL) rotate_brain(50.0,'x');
                    if (plane==HORIZONTAL) rotate_brain(-50.0,'z');
                    if (plane==CORONAL) rotate_brain(50.0,'y');
                    redraw();
                    break;
          /* rot z cw */
          case ']': if (plane==SAGITTAL) rotate_brain(-5.0,'x');
                    if (plane==HORIZONTAL) rotate_brain(5.0,'z');
                    if (plane==CORONAL) rotate_brain(-5.0,'y');
                    redraw();
                    break;
          case '}': if (plane==SAGITTAL) rotate_brain(-50.0,'x');
                    if (plane==HORIZONTAL) rotate_brain(50.0,'z');
                    if (plane==CORONAL) rotate_brain(-50.0,'y');
                    redraw();
                    break;
          /* trans +y */
          case 'p': if (plane==SAGITTAL) translate_brain(0.5,'z');
                    if (plane==HORIZONTAL) translate_brain(0.5,'y');
                    if (plane==CORONAL) translate_brain(0.5,'z');
                    redraw();
                    break;
          case 'P': if (plane==SAGITTAL) translate_brain(2.5,'z');
                    if (plane==HORIZONTAL) translate_brain(2.5,'y');
                    if (plane==CORONAL) translate_brain(2.5,'z');
                    redraw();
                    break;
          /* trans -y */
          case '.': if (plane==SAGITTAL) translate_brain(-0.5,'z');
                    if (plane==HORIZONTAL) translate_brain(-0.5,'y');
                    if (plane==CORONAL) translate_brain(-0.5,'z');
                    redraw();
                    break;
          case '>': if (plane==SAGITTAL) translate_brain(-2.5,'z');
                    if (plane==HORIZONTAL) translate_brain(-2.5,'y');
                    if (plane==CORONAL) translate_brain(-2.5,'z');
                    redraw();
                    break;
          /* trans -x */
          case 'l': if (plane==SAGITTAL) translate_brain(-0.5,'y');
                    if (plane==HORIZONTAL) translate_brain(0.5,'x');
                    if (plane==CORONAL) translate_brain(0.5,'x');
                    redraw();
                    break;
          case 'L': if (plane==SAGITTAL) translate_brain(-2.5,'y');
                    if (plane==HORIZONTAL) translate_brain(2.5,'x');
                    if (plane==CORONAL) translate_brain(2.5,'x');
                    redraw();
                    break;
          /* trans +x */
          case ';': if (plane==SAGITTAL) translate_brain(0.5,'y');
                    if (plane==HORIZONTAL) translate_brain(-0.5,'x');
                    if (plane==CORONAL) translate_brain(-0.5,'x');
                    redraw();
                    break;
          case ':': if (plane==SAGITTAL) translate_brain(2.5,'y');
                    if (plane==HORIZONTAL) translate_brain(-2.5,'x');
                    if (plane==CORONAL) translate_brain(-2.5,'x');
                    redraw();
                    break;
        } /* End switch (char)val */
        break; /* End case KEYBD */
    }
  }
#endif

}

int
open_window(char *name)

{
#ifdef OPENGL
  XSizeHints hin;

  if (openglwindowflag) {
    printf("medit: ### GL window already open: can't open second\n");
    PR return(0);
  }

  /* TKO_DEPTH because all OpenGL 4096 visuals have depth buffer!! */
  /* tkoInitDisplayMode(TKO_DOUBLE | TKO_INDEX | TKO_DEPTH);  */
  /*  tkoInitDisplayMode(TKO_DOUBLE | TKO_INDEX );  */
   tkoInitDisplayMode(TKO_RGB | TKO_SINGLE);   

  if (!initpositiondoneflag)
    tkoInitPosition(MOTIF_XFUDGE+wx0,(1024-wy0-ydim)+MOTIF_XFUDGE,xdim,ydim);
  if (!tkoInitWindow(name)) {
    printf("medit: ### tkoInitWindow(name) failed\n");exit(1);}
  hin.max_width = hin.max_height = 4*xnum + xnum/2;  /* maxsize */
  hin.min_aspect.x = hin.max_aspect.x = xdim;        /* keepaspect */
  hin.min_aspect.y = hin.max_aspect.y = ydim;
  hin.flags = PMaxSize|PAspect;
  XSetWMNormalHints(xDisplay, w.wMain, &hin);

  /* TODO: bitmap drawing is slower */
  /* (did: glDisable's,4bit,RasPos:+0.5,glOrtho-1,z=0,wintitle,clrdepth) */
  /* (did: glPixelStorei:GL_UNPACK_ALIGNMENT,glOrtho not -1) */
  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glOrtho(0.0, (double)(xnum-1), 0.0, (double)(ynum-1), -1.0, 1.0);
  /*  0.0 to x,y-1.0:  COR,SAG surf up2; HOR surf *correct* */
  /*  0.0 to x,y:      COR,SAG surf up/left; HOR surf down/left */
  /* -0.5 to x,y-0.5:  COR,SAG surf up/left; HOR surf down/left */
  /*  0.5 to x,y+0.5:  fails (??!!) */
  /*  0.0 to x-1,y:    COR,SAG surf up; HOR down */
  /*  0.0 to x,y-1:    COR,SAG surf up2/left; HOR left */
  /*  0.0 to x-1,y+1:  COR,SAG surf *correct*; HOR down */
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

#else  /* use gl calls */
  if (openglwindowflag) {
    printf("medit: ### GL window already open: can't open second\n");PR return;}
  prefposition(wx0+MOTIF_XFUDGE, wx0+MOTIF_XFUDGE+xdim-1,
               wy0-MOTIF_XFUDGE, wy0+ydim-1-MOTIF_XFUDGE);
  foreground();  /* tcl prompt */
  winopen(name);
  keepaspect(1,1);
  stepunit((long)xnum,(long)ynum);
  maxsize((long)(4*xnum),(long)(4*ynum));
  winconstraints();
  cmode();
  doublebuffer();
  gconfig();
  ortho2(0.0, (float)(xnum-1), 0.0, (float)(ynum-1));  /* orig:COR,HOR not -1 */
  zbuffer(FALSE);

  qdevice(ESCKEY);
  qdevice(REDRAW);
  qdevice(LEFTMOUSE);
  qdevice(RIGHTMOUSE);
  qdevice(MIDDLEMOUSE);
  qdevice(LEFTARROWKEY);
  qdevice(RIGHTARROWKEY);
  qdevice(UPARROWKEY);
  qdevice(DOWNARROWKEY);
  qdevice(LEFTALTKEY);
  qdevice(RIGHTALTKEY);
  qdevice(KEYBD);
#endif

  openglwindowflag = TRUE;
  return(0);
}

void
resize_window_intstep(int newzf)

{
#ifdef OPENGL
  int tzf;

  if (newzf==0)
    tzf = rint((float)w.w/(float)xnum); /* int constrain interactive resize */
  else
    tzf = newzf;
  tzf = (tzf<1)?1:(tzf>4)?4:tzf;
  if (w.w%xnum || w.h%ynum || newzf>0) {
    ozf = zf;
    zf = tzf;
    fsf = (float)zf;
    xdim = zf*xnum;
    ydim = zf*ynum;
    XResizeWindow(xDisplay, w.wMain, xdim, ydim);
    if (TKO_HAS_OVERLAY(w.type))
      XResizeWindow(xDisplay, w.wOverlay, xdim, ydim);
    glViewport(0, 0, xdim, ydim); 
    resize_buffers(xdim, ydim);
    w.w = xdim;
    w.h = ydim;

    surflinewidth=zf;
    imc = (zf*imc)/ozf;
    ic = (zf*ic)/ozf;
    jc = (zf*jc)/ozf;
  }
#else
  printf("medit: ### OpenGL only\n");
#endif
}

void 
move_window( int x, int y)
{
#ifdef OPENGL
  if (openglwindowflag) {
    XMoveWindow(xDisplay, w.wMain, x, y);
    w.x = x;
    w.y = y;
  }
  else if (!initpositiondoneflag) {
    tkoInitPosition(x,y,xdim,ydim);
    initpositiondoneflag = TRUE;
  }
  else ;
#endif
}
void
edit_pixel(int action)
{
  int irad,jrad,krad;
  int i,j,k;
  unsigned char pixval;

  jrad=irad=krad=prad;
  if (inplaneflag) { if(plane==CORONAL)    krad=0;
                     if(plane==HORIZONTAL) irad=0;
                     if(plane==SAGITTAL)   jrad=0; }
  for (k= -krad;k<=krad;k++)
  for (i= -irad;i<=irad;i++)
  for (j= -jrad;j<=jrad;j++) {
    if (imc/zf+k>=0 && imc/zf+k<numimg &&
        (ydim-1-ic)/zf+i>=0 && (ydim-1-ic)/zf+i<ynum &&
        jc/zf+j>=0 && jc/zf+j<xnum) {
      if (circleflag && k*k+i*i+j*j>prad*prad) continue;
      pixval = im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j];
      if      (action==TO_WHITE && pixval<WM_MIN_VAL) pixval = 255;
      else if (action==TO_BLACK)                      pixval = WM_EDITED_OFF;
      else ;
      im[imc/zf+k][(ydim-1-ic)/zf+i][jc/zf+j] = pixval;
    }
  }
  for (k= -prad;k<=prad;k++) {
    if (imc/zf+k>=0 && imc/zf+k<numimg) {
      changed[imc/zf+k] = TRUE;
      editedimage = imnr0+imc/zf+k;
    }
  }
}
void
save_rgb(char *fname)
{
  if (!openglwindowflag) {
    printf("medit: ### save_rgb failed: no gl window open\n");PR return; }

  if (scrsaveflag) { scrsave_to_rgb(fname); }
  else             { pix_to_rgb(fname);     }
}

void
scrsave_to_rgb(char *fname)  /* about 2X faster than pix_to_rgb */
{
  char command[2*NAME_LENGTH];
  FILE *fp;
  long xorig,xsize,yorig,ysize;
  int x0,y0,x1,y1;

  getorigin(&xorig,&yorig);
  getsize(&xsize,&ysize);

  x0 = (int)xorig;  x1 = (int)(xorig+xsize-1);
  y0 = (int)yorig;  y1 = (int)(yorig+ysize-1);
  fp = fopen(fname,"w");
  if (fp==NULL) {printf("medit: ### can't create file %s\n",fname);PR return;}
  fclose(fp);
  if (bwflag) sprintf(command,"scrsave %s %d %d %d %d -b\n",fname,x0,x1,y0,y1);
  else        sprintf(command,"scrsave %s %d %d %d %d\n",fname,x0,x1,y0,y1);
  system(command);
  printf("medit: file %s written\n",fname);PR
}

void
pix_to_rgb(char *fname)
{
#ifdef OPENGL
  GLint swapbytes, lsbfirst, rowlength;
  GLint skiprows, skippixels, alignment;
  RGB_IMAGE *image;
  int y,width,height,size;
  unsigned short *r,*g,*b;
  unsigned short  *red, *green, *blue;
  FILE *fp;

  if (bwflag) {printf("medit: ### bw pix_to_rgb failed--TODO\n");PR return;}

  width = (int)xdim;
  height = (int)ydim;
  size = width*height;

  red = (unsigned short *)calloc(size, sizeof(unsigned short));
  green = (unsigned short *)calloc(size, sizeof(unsigned short));
  blue = (unsigned short *)calloc(size, sizeof(unsigned short));

  glGetIntegerv(GL_UNPACK_SWAP_BYTES, &swapbytes);
  glGetIntegerv(GL_UNPACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength);
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows);
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

  glReadPixels(0, 0, width, height, GL_RED,GL_UNSIGNED_SHORT, (GLvoid *)red);
  glReadPixels(0, 0, width, height, GL_GREEN,GL_UNSIGNED_SHORT,(GLvoid *)green);  glReadPixels(0, 0, width, height, GL_BLUE,GL_UNSIGNED_SHORT, (GLvoid *)blue);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes);
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels);
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment);

  fp = fopen(fname,"w");
  if (fp==NULL){printf("medit: ### can't create file %s\n",fname);PR return;}
  fclose(fp);
  image = iopen(fname,"w",RLE(1), 3, width, height, 3);
  for(y = 0 ; y < height; y++) {
    r = red + y * width;
    g = green + y * width;
    b = blue + y * width;
    putrow(image, r, y, 0);
    putrow(image, g, y, 1);
    putrow(image, b, y, 2);
  }
  iclose(image);
  free(red); free(green); free(blue);

  printf("medit: file %s written\n",fname);PR
#else
  printf("medit: ### pix_to_rgb implemented only in OpenGL version\n");PR
  return;
#endif
}

void
downslice(void)
{
  if (plane==CORONAL)
    imc = (imc<zf)?imnr1*zf-zf+imc:imc-zf;
  else if (plane==HORIZONTAL)
    ic = (ic<zf)?ydim-zf+ic:ic-zf;
  else if (plane==SAGITTAL)
    jc = (jc<zf)?xdim-zf+jc:jc-zf;
}

void
upslice(void)
{
  if (plane==CORONAL)
    imc = (imc>=imnr1*zf-zf)?imc+zf-imnr1*zf:imc+zf;
  else if (plane==HORIZONTAL)
    ic = (ic>=ydim-zf)?ic+zf-ydim:ic+zf;
  else if (plane==SAGITTAL)
    jc = (jc>=xdim-zf)?jc+zf-xdim:jc+zf;
}

void
goto_point(char *dir)
{
  char fname[NAME_LENGTH];
  FILE *fp;
  float xpt,ypt,zpt;

  sprintf(fname,"%s/edit.dat",dir);
  fp=fopen(fname,"r");
  if (fp==NULL) {printf("medit: ### File %s not found\n",fname);PR return;}
  fscanf(fp,"%f %f %f",&xpt,&ypt,&zpt);
  fclose(fp);
  set_cursor(xpt,ypt,zpt);
}
void
unmark_vertices(void)
{
  int vno ;

  if (!mris)
  {
    fprintf(stderr, "no surface loaded.\n") ;
    return ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].marked = 0 ;
  redraw() ;
}
void
mark_file_vertices(char *fname)
{
  FILE  *fp ;
  char  line[200], *cp ;
  int   vno, nvertices, nargs ;
  float area ;

  if (!mris)
  {
    fprintf(stderr, "no surface loaded.\n") ;
    return ;
  }
  fp = fopen(fname, "r") ;
  if (!fp)
  {
    fprintf(stderr, "could not open file %s.\n", fname) ;
    return ;
  }

  fgetl(line, 150, fp) ;
  nargs = sscanf(line, "%d %f", &nvertices, &area) ;
  if (nargs == 2)
    fprintf(stderr, "marking %d vertices, %2.3f mm^2 surface area\n",
            nvertices, area) ;
  else if (nargs == 1)
    fprintf(stderr, "marking %d vertices\n", nvertices) ;

  while  ((cp = fgetl(line, 150, fp)) != NULL)
  {
    sscanf(cp, "%d", &vno) ;
    if (vno >= 0 && vno < mris->nvertices)
    {
      mris->vertices[vno].marked = 1 ;
    }
  }
  goto_vertex(vno) ;
  fclose(fp) ;
}
void
goto_vertex(int vno)
{
  VERTEX *v ;

  if (!mris)
  {
    fprintf(stderr, "surface not loaded.\n") ;
    return ;
  }
  v = &mris->vertices[vno] ;
  set_cursor(v->x, v->y, v->z);
  redraw() ;
}


void
goto_point_coords(int imc1, int ic1,int jc1)
{
  float xpt,ypt,zpt;

  xpt = xx1-ps*jc1/fsf;
  ypt = yy0+st*imc1/fsf;
  zpt = zz1-ps*(255.0-ic1/fsf);
  set_cursor(xpt,ypt,zpt);
}

void
write_point(char *dir)
{
  char fname[NAME_LENGTH];
  FILE *fp;
  float xpt,ypt,zpt;
  Real  x, y, z, x_tal, y_tal, z_tal ;

  sprintf(fname,"%s/edit.dat",dir);
  fp=fopen(fname,"w");
  if (fp==NULL) {printf("medit: ### can't create file %s\n",fname);PR return;}
  xpt = xx1-ps*jc/fsf;
  ypt = yy0+st*imc/fsf;
  zpt = zz1-ps*(255.0-ic/fsf);
  fprintf(fp,"%f %f %f\n",xpt,ypt,zpt);
  if (transform_loaded)
  {
    fprintf(stderr, "writing transformed point to file...\n") ;
    x = (Real)xpt ; y = (Real)ypt ; z = (Real)zpt ;
    transform_point(linear_transform, x, y, z, &x_tal, &y_tal, &z_tal) ;
    fprintf(fp, "%f %f %f\n", x_tal, y_tal, z_tal) ;
    xtalairach = x_tal;  ytalairach = y_tal;  ztalairach = z_tal; 
  }
  /*else { fprintf(stderr, "NOT writing transformed point to file...\n") ; }*/
  fclose(fp);
}

void coords_to_talairach(void)
{
  float xpt,ypt,zpt;
  Real  x, y, z, x_tal, y_tal, z_tal ;

  xpt = xx1-ps*jc/fsf;
  ypt = yy0+st*imc/fsf;
  zpt = zz1-ps*(255.0-ic/fsf);
  x = (Real)xpt ; y = (Real)ypt ; z = (Real)zpt ;
  if (transform_loaded) {
    transform_point(linear_transform, x, y, z, &x_tal, &y_tal, &z_tal) ;
    xtalairach = x_tal;  ytalairach = y_tal;  ztalairach = z_tal; 
  }
}

void talairach_to_coords(void)
{
  int nimc,nic,njc;
  float xpt,ypt,zpt;
  double dzf;
  Real  x, y, z, x_tal, y_tal, z_tal ;

  if (transform_loaded) {
    x_tal = (Real)xtalairach; y_tal = (Real)ytalairach; z_tal =(Real)ztalairach;
    transform_point(inverse_linear_transform, x_tal, y_tal, z_tal, &x, &y, &z) ;
    xpt = (float)x; ypt = (float)y; zpt = (float)z;
    if (ptype==0) /* Horizontal */
    {
      jpt = (int)((xx1-xpt)*zf/ps+0.5);
      ipt = (int)((ypt-yy0)*zf/ps+0.5);
      impt = (int)((zpt-zz0)*zf/st+0.5);
    } else if (ptype==2) /* Coronal */
    {
      jpt = (int)((xx1-xpt)*zf/ps+0.5);
      ipt = (int)((255.0-(zz1-zpt)/ps)*zf+0.5);
      impt = (int)((ypt-yy0)*zf/st+0.5);
    } else if (ptype==1) /* Sagittal */
    {
      jpt = (int)((xx1-xpt)*zf/ps+0.5);
      ipt = (int)((ypt-yy0)*zf/ps+0.5);
      impt = (int)((zpt-zz0)*zf/st+0.5);
    }
    dzf = (double)zf;
    nimc = zf * (int)(rint((double)impt/dzf));  /* round to slice */
    nic =  zf * (int)(rint((double)ipt/dzf));
    njc =  zf * (int)(rint((double)jpt/dzf));
    if (njc/zf>=0 && njc/zf<xnum && nimc/zf>=0 && nimc/zf<=numimg &&
        (ydim-1-nic)/zf>=0 && (ydim-1-nic)/zf<ynum) {
      jc = njc;
      imc = nimc;
      ic = nic;
    }
    jpt=ipt=impt = -1;
  }
}

void set_cursor(float xpt, float ypt, float zpt)
{
  double dzf;
  Real   x, y, z, x_tal, y_tal, z_tal ;
  int    secondpixelval, xi, yi, zi ;

  x = y = z = 0.0;
  if (ptype==0) /* Horizontal */
  {
    jpt = (int)((xx1-xpt)*zf/ps+0.5);
    ipt = (int)((ypt-yy0)*zf/ps+0.5);
    impt = (int)((zpt-zz0)*zf/st+0.5);
  } else if (ptype==2) /* Coronal */
  {
    jpt = (int)((xx1-xpt)*zf/ps+0.5);
    ipt = (int)((255.0-(zz1-zpt)/ps)*zf+0.5);
    impt = (int)((ypt-yy0)*zf/st+0.5);
  } else if (ptype==1) /* Sagittal */
  {
    jpt = (int)((xx1-xpt)*zf/ps+0.5);
    ipt = (int)((ypt-yy0)*zf/ps+0.5);
    impt = (int)((zpt-zz0)*zf/st+0.5);
  }

  dzf = (double)zf;
  imc = zf * (int)(rint((double)impt/dzf));  /* round to slice */
  ic =  zf * (int)(rint((double)ipt/dzf));
  jc =  zf * (int)(rint((double)jpt/dzf));
  /*jc = jpt;*/
  /*ic = ipt;*/
  /*imc = impt;*/
  jpt=ipt=impt = -1;  /* crosshair off */

  if (imc/zf>=0 && imc/zf<imnr1) 
  {
    xi = (int)(jc/zf) ; yi = (int)((ydim-1-ic)/zf) ; zi = (int)(imc/zf) ;
    selectedpixval = im[zi][yi][xi];
    printf("val=%d ",selectedpixval);
    if (second_im_allocated) 
    {
      secondpixelval = im2[zi][yi][xi];
      printf("(2nd val=%d ",secondpixelval);
    }
  }
  if (ptype==0) /* Horizontal */
  {
    printf(
    "imnr(P/A)=%d, i(I/S)=%df, j(R/L)=%d (x=%2.1f y=%2.1f z=%2.1f)\n",
       zi,yi,xi,xx1-ps*jc/fsf,yy0+ps*ic/fsf,zz0+st*imc/fsf);
    x = (Real)(xx1-ps*jc/fsf) ;
    y = (Real)(yy0+ps*ic/fsf) ;
    z = (Real)(zz0+st*imc/fsf) ;
  } else if (ptype==2) /* Coronal */
  {
    printf(
    "imnr(P/A)=%d, i(I/S)=%d, j(R/L)=%d (x=%2.1f y=%2.1f z=%2.1f)\n",
            zi,yi,xi,xx1-ps*jc/fsf,yy0+st*imc/fsf,zz1-ps*(255.0-ic/fsf));
    x = (Real)(xx1-ps*jc/fsf) ;
    y = (Real)(yy0+st*imc/fsf) ;
    z = (Real)(zz1-ps*(255.0-ic/fsf)) ;
  } else if (ptype==1) /* Sagittal */
  {
    printf(
    "imnr(P/A)=%d, i(I/S)=%d, j(R/L)=%d (x=%2.1f y=%2.1f z=%2.1f)\n",
         zi,yi,xi,xx1-ps*jc/fsf,yy0+ps*ic/fsf,zz0+st*imc/fsf);
    x = (Real)(xx1-ps*jc/fsf) ;
    y = (Real)(yy0+ps*ic/fsf) ;
    z = (Real)(zz0+st*imc/fsf) ;
  }
  if (transform_loaded)
  {
    transform_point(linear_transform, x, y, z, &x_tal, &y_tal, &z_tal) ;
    printf("TALAIRACH: (%2.1f, %2.1f, %2.1f)\n", x_tal, y_tal, z_tal) ;
    xtalairach = x_tal;  ytalairach = y_tal;  ztalairach = z_tal; 
  }

  /* twitzels stuff */
  if(statsloaded)
    printOutFunctionalCoordinate(x,y,z);
  updatepixval = TRUE;
  PR
}

void set_scale(void)
{
  Colorindex i;
  float f;
  short v;

  for (i=0;i<NUMVALS;i++)
  {
/*
    f = i/fscale+fthresh;
    f = ((f<0.0)?0.0:((f>1.0)?1.0:f));
    f = pow(f,fsquash);
    v = f*fscale+0.5;
*/
    if (linearflag)
      v = i;
    else {
      f = i/fscale;
      f = 1.0/(1.0+exp(-fsquash*(f-fthresh)));
      v = f*fscale+0.5;
    }
    mapcolor(i+MAPOFFSET,v,v,v);
  }
  mapcolor(NUMVALS+MAPOFFSET,v,0.0,0.0);
}

#define FIRST_TITLE "editable image: name of dir containing COR images"
#define SECOND_TITLE  "Second Buffer"


void redraw(void)
{
#ifdef OPENGL
  XTextProperty tp;
#else
  char title[5*NAME_LENGTH];
#endif


  if (!openglwindowflag) {
    printf("medit: ### redraw failed: no gl window open\n");PR return; }

  set_scale();
  color(BLACK);
  /* clear(); */
  if (drawsecondflag && second_im_allocated) {
    draw_second_image(imc,ic,jc);
    wintitle(imtype2);
  }
  else {
    draw_image_hacked(imc,ic,jc);
    /* draw_image(imc,ic,jc); */
    wintitle(imtype);
  }
  if (ptsflag && !all3flag) drawpts();
  if (surfflag && surfloaded)
    draw_surface();
  if (surfflag && !surfloaded) {
    printf("medit: ### no surface read\n");PR
    surfflag=FALSE;
  }
  /*  swapbuffers(); */
}

void pop_gl_window(void)
{
#ifdef OPENGL
  XRaiseWindow(xDisplay, w.wMain);
#else
  winpop();
#endif
}

void resize_buffers(int x,int y)
{
  free(vidbuf);
  free(cvidbuf);
  free(binbuff);
  vidbuf = (GLubyte*)lcalloc((size_t)xdim*ydim*4,(size_t)sizeof(GLubyte));  

  /*  vidbuf = (Colorindex *)lcalloc((size_t)x*y,(size_t)sizeof(Colorindex));  */
  cvidbuf = (Colorindex *)lcalloc((size_t)CVIDBUF,(size_t)sizeof(Colorindex));
  binbuff = (unsigned char *)lcalloc((size_t)3*x*y,(size_t)sizeof(char));
}

void
mri2pix(float xpt, float ypt, float zpt, int *jpt, int *ipt,int *impt)
{
  if (ptype==0) /* Horizontal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((ypt-yy0)/ps+0.5);
    *impt = (int)((zpt-zz0)/st+0.5);
  } else if (ptype==2) /* Coronal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((255.0-(zz1-zpt)/ps)+0.5);
    *impt = (int)((ypt-yy0)/st+0.5);
  } else if (ptype==1) /* Sagittal */
  {
    *jpt = (int)((xx1-xpt)/ps+0.5);
    *ipt = (int)((ypt-yy0)/ps+0.5);
    *impt = (int)((zpt-zz0)/st+0.5);
  }
}

int 
imval(float px,float py,float pz)
{
  float x,y,z;
  int j,i,imn;

  x = px*tm[0][0]+py*tm[0][1]+pz*tm[0][2]+tm[0][3];
  y = px*tm[1][0]+py*tm[1][1]+pz*tm[1][2]+tm[1][3];
  z = px*tm[2][0]+py*tm[2][1]+pz*tm[2][2]+tm[2][3];
  mri2pix(x,y,z,&j,&i,&imn);
  if (imn>=0&&imn<numimg&&i>=0&&i<ynum&&j>=0&&j<xnum)
    return(im[imn][ynum-1-i][j]);
  else return 0;
}

float
Error(int p,float dp)
{
  int i,num;
  float mu,error,sum;
  float mu1,mu2,sum1,sum2;

  if (p>=0)
    par[p] += dp;
  mu = mu1 = mu2 = 0;
  num = 0;
  for (i=0;i<npts;i++)
  {
    mu += imval(ptx[i],pty[i],ptz[i]);
    mu1 += imval(ptx[i]*0.9,pty[i]*0.9,ptz[i]*0.9);
    mu2 += imval(ptx[i]*1.05,pty[i]*1.05,ptz[i]*1.05);
    num ++;
  }
  mu /= num;
  mu1 /= num;
  mu2 /= num;
  sum = sum1 = sum2 = 0;
  num = 0;
  for (i=0;i<npts;i++)
  {
    error = imval(ptx[i],pty[i],ptz[i])-mu;
    sum += error*error;
    error = imval(ptx[i]*0.9,pty[i]*0.9,ptz[i]*0.9)-mu1;
    sum1 += error*error;
    error = imval(ptx[i]*1.05,pty[i]*1.05,ptz[i]*1.05)-mu2;
    sum2 += error*error;
    num ++;
  }
  sum = sqrt((sum+sum2)/num);
  if (p>=0)
    par[p] -= dp;
  return sum;
}

void
optimize(int maxiter)
{
  float lambda = 0.03;
  float epsilon = 0.1;
  float momentum = 0.8;
  int iter,p;
  float dE[3];
  float error;

  for (iter=0;iter<maxiter;iter++)
  {
    error = Error(-1,0);
    printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
           iter,par[0],par[1],par[2],error);PR
    for (p=0;p<3;p++)
    {
      dE[p] = tanh((Error(p,epsilon/2)-Error(p,-epsilon/2))/epsilon);
    }
    for (p=0;p<3;p++)
    {
      par[p] += (dpar[p] = momentum*dpar[p] - lambda*dE[p]);
    }
  }
  error = Error(-1,0);
  printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
         iter,par[0],par[1],par[2],error);PR
}
void
optimize2(void)
{
  float epsilon = 0.5;
  float p0,p1,p2,p0min,p1min,p2min;
  float error,minerror;

  p0min = p1min = p2min = 0.0f;
  error = Error(-1,0);
  minerror = error;
  printf("%5.2f %5.2f %5.2f %7.3f\n",
         par[0],par[1],par[2],error);PR
  for (p0 = -10;p0 <= 10;p0 += epsilon)
  for (p1 = -10;p1 <= 10;p1 += epsilon)
  for (p2 = -10;p2 <= 10;p2 += epsilon)
  {
    par[0] = p0;
    par[1] = p1;
    par[2] = p2;
    error = Error(-1,0);
    if (error<minerror)
    {
      printf("%5.2f %5.2f %5.2f %7.3f\n",
             par[0],par[1],par[2],error);PR
      minerror = error;
      p0min = p0;
      p1min = p1;
      p2min = p2;
    }
  }
  par[0] = p0min;
  par[1] = p1min;
  par[2] = p2min;
  error = Error(-1,0);
  printf("%5.2f %5.2f %5.2f %7.3f\n",
         par[0],par[1],par[2],error);PR
}

void
translate_brain(float a,char c)
{
  int i,j,k;
  float m1[4][4],m2[4][4];

  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    m1[i][j] = (i==j)?1.0:0.0;
  if (c=='y')
    m1[1][3] = a;
  else if (c=='x')
    m1[0][3] = a;
  else if (c=='z')
    m1[2][3] = a;
  else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  {
    m2[i][j] = 0;
    for (k=0;k<4;k++)
      m2[i][j] += tm[i][k]*m1[k][j];
  }
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = m2[i][j];
}

void
rotate_brain(float a,char c)
{
  int i,j,k;
  float m1[4][4],m2[4][4];
  float sa,ca;

  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    m1[i][j] = (i==j)?1.0:0.0;
  a = a*M_PI/1800;
  sa = sin(a);
  ca = cos(a);
  if (c=='y')
  {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  } else
  if (c=='x')
  {
    m1[1][1] = m1[2][2] = ca;
    m1[1][2] = -(m1[2][1] = sa);
  } else
  if (c=='z')
  {
    m1[0][0] = m1[1][1] = ca;
    m1[0][1] = -(m1[1][0] = sa);
  } else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
  {
    m2[i][j] = 0;
    for (k=0;k<4;k++)
      m2[i][j] += tm[i][k]*m1[k][j];
  }
  for (i=0;i<4;i++)
  for (j=0;j<4;j++)
    tm[i][j] = m2[i][j];
}

int
read_images(char *fpref)
{
  int i,j,k;
  FILE *fptr, *xfptr;
  char fname[100], char_buf[100];

  sprintf(fname,"%s.info",fpref);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {
    strcpy(fpref,"COR-");
    sprintf(fname,"%s.info",fpref);
    fptr = fopen(fname,"r");
    if (fptr==NULL) {
      printf("medit: ### File %s not found (tried local dir, too)\n",fname);
      exit(1);
    }
  } 

  fscanf(fptr,"%*s %d",&imnr0);
  fscanf(fptr,"%*s %d",&imnr1);
  fscanf(fptr,"%*s %d",&ptype);
  fscanf(fptr,"%*s %d",&xnum);
  fscanf(fptr,"%*s %d",&ynum);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %f",&ps);
  fscanf(fptr,"%*s %f",&st);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %f",&xx0); /* strtx */
  fscanf(fptr,"%*s %f",&xx1); /* endx */
  fscanf(fptr,"%*s %f",&yy0); /* strty */
  fscanf(fptr,"%*s %f",&yy1); /* endy */
  fscanf(fptr,"%*s %f",&zz0); /* strtz */
  fscanf(fptr,"%*s %f",&zz1); /* endz */
  fscanf(fptr, "%*s %*f") ;   /* tr */
  fscanf(fptr, "%*s %*f") ;   /* te */
  fscanf(fptr, "%*s %*f") ;   /* ti */
  fscanf(fptr, "%*s %s",char_buf);

  /* marty: ignore abs paths in COR-.info */
  xfptr = fopen(xffname,"r");
  if (xfptr==NULL)
    printf("medit: Talairach xform file not found (ignored)\n");
  else {
    fclose(xfptr);
    if (input_transform_file(xffname, &talairach_transform) != OK)
      printf("medit: ### File %s load failed\n",xffname);
    else {
      printf("medit: Talairach transform file loaded\n");
      printf("medit: %s\n",xffname);
      linear_transform = get_linear_transform_ptr(&talairach_transform) ;
      inverse_linear_transform =
                      get_inverse_linear_transform_ptr(&talairach_transform) ;
      transform_loaded = TRUE;
    }
  }

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
  xdim=xnum*zf;
  ydim=ynum*zf;

  bufsize = ((unsigned long)xnum)*ynum;
  buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
  for (k=0;k<numimg;k++)
  {
    im[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++)
    {
      im[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
  }
  for (k=0;k<6;k++)
  {
    sim[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++)
    {
      sim[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
  }

  for (k=0;k<numimg;k++)
  {
    changed[k] = FALSE;
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"rb");
    if(fptr==NULL){printf("medit: ### File %s not found.\n",fname);exit(1);}
    fread(buf,sizeof(char),bufsize,fptr);
    buffer_to_image(buf,im[k],xnum,ynum);
    fclose(fptr);
  }

  for (k=0;k<numimg;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
  {
    if (im[k][i][j]/2>sim[3][i][j]) sim[3][i][j] = im[k][i][j]/2;
    if (im[k][i][j]/2>sim[4][k][j]) sim[4][k][j] = im[k][i][j]/2;
    if (im[k][i][j]/2>sim[5][i][k]) sim[5][i][k] = im[k][i][j]/2;
  }
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
  for (k=0;k<3;k++)
    sim[k][i][j] = sim[k+3][i][j];
  return(0);
}


int
read_second_images(char *imdir2)
{
  int i,j,k,n;
  FILE *fptr;
  char fname[NAME_LENGTH], cmd[NAME_LENGTH], fpref[NAME_LENGTH];
  char fnamefirst[NAME_LENGTH], notilde[NAME_LENGTH];

  strcpy(imtype2, imdir2) ;


  /* TODO: replace w/setfile (add mri/<subdir>); tk.c: omit arg like others */
  if (imdir2[0] == '/')
    sprintf(fpref,"%s/COR-",imdir2);
  else if (imdir2[0] == '~') {
    for(n=1;n<=strlen(imdir2);n++)  notilde[n-1] = imdir2[n];
    sprintf(fpref,"%s/%s/%s/COR-",subjectsdir,pname,notilde);
  }
  else if (MATCH(pname,"local"))
    sprintf(fpref,"%s/%s/COR-",srname,imdir2);
  else
    sprintf(fpref,"%s/%s/mri/%s/COR-",subjectsdir,pname,imdir2);

  sprintf(fnamefirst,"%s.info",mfname);
  sprintf(fname,"%s.info",fpref);
  sprintf(cmd,"diff %s %s",fnamefirst,fname);
  if (system(cmd)!=0)
    printf("medit: ### second COR-.info file doesn't match first--ignored\n");

  if (!second_im_allocated)
    alloc_second_im();
 
  for (k=0;k<numimg;k++) {
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"rb");
    if(fptr==NULL){printf("medit: ### File %s not found.\n",fname);return(0);}
    fread(buf,sizeof(char),bufsize,fptr);
    buffer_to_image(buf,im2[k],xnum,ynum);
    fclose(fptr);
  }

  for (k=0;k<numimg;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    if (im2[k][i][j]/2>sim2[3][i][j]) sim2[3][i][j] = im2[k][i][j]/2;
    if (im2[k][i][j]/2>sim2[4][k][j]) sim2[4][k][j] = im2[k][i][j]/2;
    if (im2[k][i][j]/2>sim2[5][i][k]) sim2[5][i][k] = im2[k][i][j]/2;
  }
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
  for (k=0;k<3;k++)
    sim2[k][i][j] = sim2[k+3][i][j];

  drawsecondflag = TRUE;
  second_im_full = TRUE;
  redraw();
  return(0);
}

int
write_images(char *fpref)
{
  int k;
  FILE *fptr;
  char fname[100];

  if (!editflag) { 
    printf(
      "medit: ### can't write images: (set editflag TRUE)\n");return(0);PR }
  for (k=0;k<numimg;k++) {
    if (changed[k]) {
      changed[k] = FALSE;
      file_name(fpref,fname,k+imnr0,"%03d");
      fptr = fopen(fname,"wb");
      image_to_buffer(im[k],buf,xnum,ynum);
      fwrite(buf,sizeof(char),bufsize,fptr);
      fclose(fptr);
      printf("medit: file %s written\n",fname);PR
    }
  }
  editedimage = FALSE;
  return(0);
}
int
dump_vertex(int vno)
{
  VERTEX *v, *vn ;
  int    n ;
  float  dx, dy, dz, dist ;

  if (!mris)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "%s: surface must be loaded\n", Progname)) ;

  v = &mris->vertices[vno] ;
  fprintf(stderr, 
          "v %d: x = (%2.2f, %2.2f, %2.2f), n = (%2.2f, %2.2f, %2.2f)\n",
          vno, v->x, v->y, v->z, v->nx, v->ny, v->nz) ;
  fprintf(stderr, "nbrs:\n") ;
  for (n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    dx = vn->x - v->x ; dy = vn->y - v->y ; dz = vn->z - v->z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    fprintf(stderr, "\tvn %d, dist = %2.2f\n", v->v[n], dist) ;
  }
  return(NO_ERROR) ;
}

int
show_vertex(void)
{
  VERTEX   *v ;
  int      vno, min_vno ;
  float    dx, dy, dz, dist = 0.0, min_dist ;

  if (!mris)
    ErrorReturn(ERROR_BADPARM, 
                (ERROR_BADPARM, "%s: surface must be loaded\n", Progname)) ;

  min_vno = -1 ; min_dist = 1e9 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    dx = v->x - x_click ; dy = v->y - y_click ; dz = v->z - z_click ; 
    dist = dx*dx + dy*dy + dz*dz ;
    if (dist < min_dist)
    {
      min_dist = dist ;
      min_vno = vno ;
    }
  }
  v = &mris->vertices[min_vno] ;
  fprintf(stderr, "vno %d @ (%2.1f, %2.1f, %2.1f), %2.1f mm away\n",
          min_vno, v->x, v->y, v->z, sqrt(min_dist)) ;
  return(NO_ERROR) ;
}

void
select_pixel( short sx, short sy, int printflag)
{
  long ox,oy,lx,ly;
  Real  x, y, z, x_tal, y_tal, z_tal ;
  int   xi,yi, zi ;

  x = y = z = 0.0;
  getorigin(&ox,&oy);
  getsize(&lx,&ly);

  /* hack: x-y of click on tk win caught by qread when clicks buffered! */
  if (sx<ox || sx>ox+lx) sx = ox+lx/2;
  if (sy<oy || sy>oy+ly) sy = oy+ly/2;
/*
  printf("sx=%d, sy=%d, ox=%d, oy=%d, lx=%d, ly=%d\n",
         sx,sy,ox,oy,lx,ly);
*/

  if (all3flag) {
    if (sx-ox<xdim/2 && sy-oy>ydim/2) {  /* COR */
      ic = 2*((sy-oy)-ydim/2);
      jc = 2*(sx-ox);
    }
    if (sx-ox>xdim/2 && sy-oy>ydim/2) {  /* HOR */
      imc = 2*((sy-oy)-ydim/2);
      jc =  2*((sx-ox)-xdim/2);
    }
    if (sx-ox<xdim/2 && sy-oy<ydim/2) {  /* SAG */
      imc = 2*(sx-ox);
      ic = 2*(sy-oy);
    }
    if (sx-ox>xdim/2 && sy-oy<ydim/2) {  /* max=>center */
      ic = jc = imc = xdim/2;
    }
  }
  else {
    if (plane==CORONAL) {
      ic = (sy-oy);
      jc = (sx-ox);
    }
    else if (plane==HORIZONTAL) {
      imc = (sy-oy);
      jc = (sx-ox);
    }
    else if (plane==SAGITTAL) {
      imc = (sx-ox);
      ic = (sy-oy);
    }
  }
  if (printflag)
  {
    if (imc/zf>=0 && imc/zf<imnr1) 
    {
      xi = (int)(jc/zf) ; yi = (int)((ydim-1-ic)/zf) ; zi = (int)(imc/zf) ;
      selectedpixval = im[zi][yi][xi];
      printf("val=%d ",selectedpixval);
      if (second_im_allocated) 
      {
        selectedpixval =im2[zi][yi][xi];        
        printf("2nd val=%d ",selectedpixval); 
      }
    }
    if (ptype==0) /* Horizontal */
    {
      printf("imnr(P/A)=%d, i(I/S)=%d, j(R/L)=%d (x=%2.1f y=%2.1f z=%2.1f)\n",
         zi,yi,yi,xx1-ps*jc/fsf,yy0+ps*ic/fsf,zz0+st*imc/fsf);
      x = (Real)(xx1-ps*jc/fsf) ;
      y = (Real)(yy0+ps*ic/fsf) ;
      z = (Real)(zz0+st*imc/fsf) ;
    } else if (ptype==2) /* Coronal */
    {
      printf("imnr(P/A)=%d, i(I/S)=%d, j(R/L)=%d (x=%2.1f y=%2.1f z=%2.1f)\n",
             zi, yi, xi,xx1-ps*jc/fsf,yy0+st*imc/fsf,zz1-ps*(255.0-ic/fsf));
      x = (Real)(xx1-ps*jc/fsf) ;
      y = (Real)(yy0+st*imc/fsf) ;
      z = (Real)(zz1-ps*(255.0-ic/fsf)) ;
    } else if (ptype==1) /* Sagittal */
    {
      printf("imnr(P/A)=%d, i(I/S)=%d, j(R/L)=%d (x=%2.1f y=%2.1f z=%2.1f)\n",
             zi, yi, xi, xx1-ps*jc/fsf,yy0+ps*ic/fsf,zz0+st*imc/fsf);
      x = (Real)(xx1-ps*jc/fsf) ;
      y = (Real)(yy0+ps*ic/fsf) ;
      z = (Real)(zz0+st*imc/fsf) ;
    }
    x_click = x ; y_click = y ; z_click = z ;
    if (transform_loaded)
    {
      transform_point(linear_transform, x, y, z, &x_tal, &y_tal, &z_tal) ;
      printf("TALAIRACH: (%2.1f, %2.1f, %2.1f)\n", x_tal, y_tal, z_tal) ;
      xtalairach = x_tal;  ytalairach = y_tal;  ztalairach = z_tal; 
    }
    PR
  }
  /* twitzels stuff */
  if(statsloaded)
    printOutFunctionalCoordinate(x,y,z);
  /*
  hot=1;
  lookupInVoxelSpace(x,y,z,0);
  setupCoronal(imc);
  hot=0;*/
  updatepixval = TRUE;
}

void draw_image(int imc,int ic,int jc)
{
  int i,j,imnr,k;
  int hax,hay,curs;
  int zimc,zic,zjc;
  /*  Colorindex v; */
  GLubyte v;
  int idx_buf;
  hax = xdim/2;
  hay = ydim/2;
  if (all3flag) curs = 15;
  else          curs = 2;

  if (maxflag && !all3flag) {
    k = 0;
    if (plane==CORONAL)
      for (i=ydim-1;i>=0;i--)
      for (j=0;j<xdim;j++)  vidbuf[k++] = sim[0][i/zf][j/zf]+MAPOFFSET;
    else if (plane==HORIZONTAL)
      for (i=0;i<ydim;i++)
      for (j=0;j<xdim;j++)  vidbuf[k++] = sim[1][i/zf][j/zf]+MAPOFFSET;
    else if (plane==SAGITTAL)
      for (i=ydim-1;i>=0;i--)
      for (j=0;j<xdim;j++)  vidbuf[k++] = sim[2][i/zf][j/zf]+MAPOFFSET;
  }
  else {
    if (plane==CORONAL || all3flag) {
      k = 0;
      zimc = imc/zf;
      if (zimc>=0&&zimc<imnr1) {
        for (i=0;i<ydim;i++)
        for (j=0;j<xdim;j++) {
          if (all3flag && (i%2 || j%2)) continue;
          v = im[zimc][(ydim-i-1)/zf][j/zf]+MAPOFFSET;
          if (truncflag)
            if (v<white_lolim+MAPOFFSET || v>white_hilim+MAPOFFSET)
              v=MAPOFFSET;
          if (all3flag)
          {
            idx_buf = 4*((xdim*hay) + k + ((i/2)*hax)) ;
            vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf + 2] = v;
            vidbuf[idx_buf + 3]=255;
            k++;
          }
          else
          {
    /*
      vidbuf[k] = v;
      k++;
    */
            vidbuf[k] = vidbuf[k+1] = vidbuf[k+2]= v; 
            vidbuf[k+3]=255;
            k+=4;
          }
        }
        for (i=ic-curs;i<=ic+curs;i++) {
          if (all3flag) k = 4*(xdim*hay + i/2*xdim/2+jc/2 + i/2*hax);
          else      k = 4*(i*xdim+jc);
          vidbuf[k] = NUMVALS-1; vidbuf[k+1] = vidbuf[k+2] = 0; 
          vidbuf[k+3]=NUMVALS-1;
        }
        for (j=jc-curs;j<=jc+curs;j++) {
          if (all3flag) k = 4*(xdim*hay + ic/2*xdim/2+j/2 + ic/2*hax);
          else          k = 4*(ic*xdim+j);
          vidbuf[k] = NUMVALS-1; vidbuf[k+1] = vidbuf[k+2] = 0; 
          vidbuf[k+3]=NUMVALS-1;
        }
      }
      else {  /* outside=black */
        for (i=0;i<ydim;i++)
        for (j=0;j<xdim;j++) {
          if (all3flag && (i%2 || j%2)) continue;
          vidbuf[k++] = MAPOFFSET;
        }
      }
    }
    if (plane==HORIZONTAL || all3flag) 
    {
      k = 0;
      zic = (ydim-1-ic)/zf;
      if (zic>=0&&zic<ynum) {
        for (imnr=0;imnr<ydim;imnr++)
        for (j=0;j<xdim;j++) {
          if (all3flag && (imnr%2 || j%2)) continue;
          v = im[imnr/zf][zic][j/zf]+MAPOFFSET;
          if (truncflag)
            if (v<white_lolim+MAPOFFSET || v>white_hilim+MAPOFFSET)
              v=MAPOFFSET;
          if ((imnr==imc&&abs(j-jc)<=curs)||
              (j==jc&&abs(imnr-imc)<=curs)) 
            v=0 /*NUMVALS+MAPOFFSET*/;
          if (all3flag) 
          {
            idx_buf = 4*((xdim*hay) + hax + k + ((imnr/2)*hax)) ;
            vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf + 2] = v;
            vidbuf[idx_buf + 3]=255;
            k++;
          }
          else
          {
            vidbuf[k] = vidbuf[k+1] = vidbuf[k+2]= v; 
            vidbuf[k+3]=255;
            k+=4;
          }

        }
        for (imnr=imc-curs;imnr<=imc+curs;imnr++) {
          if (all3flag) k = 4*(xdim*hay+hax + imnr/2*xdim/2+jc/2 + imnr/2*hax);
          else          k = 4*(imnr*xdim+jc);
          vidbuf[k] = NUMVALS-1 ; vidbuf[k+1] = vidbuf[k+2] = 0;
          vidbuf[k+3]=255;
        }
        for (j=jc-curs;j<=jc+curs;j++) {
          if (all3flag) k = 4*(xdim*hay + hax + imc/2*xdim/2+j/2 + imc/2*hax);
          else          k = 4*(imc*xdim+j);
          vidbuf[k] = NUMVALS-1 ; vidbuf[k+1] = vidbuf[k+2] = 0;
          vidbuf[k+3]=255;
        }
      }
      else {
        for (imnr=0;imnr<ydim;imnr++)
        for (j=0;j<xdim;j++) {
          if (all3flag && (imnr%2 || j%2)) continue;
          vidbuf[k++] = MAPOFFSET;
        }
      }
    }
    if (plane==SAGITTAL || all3flag) {
      k = 0;
      zjc = jc/zf;
      if (zjc>=0&&zjc<xnum) {
        /*for (i=ydim-1;i>=0;i--)*//* TODO:4% faster, but 25% l.t. COR/HOR */
        for (i=0;i<ydim;i++)
        for (imnr=0;imnr<xdim;imnr++) {
          if (all3flag && (i%2 || imnr%2)) continue;
          /*v = im[imnr/zf][i/zf][zjc]+MAPOFFSET;*/
          v = im[imnr/zf][(ydim-1-i)/zf][zjc]+MAPOFFSET;
          if (truncflag)
            if (v<white_lolim+MAPOFFSET || v>white_hilim+MAPOFFSET)
              v=MAPOFFSET;
          if (all3flag) 
          {
            idx_buf = 4*(k + ((i/2)*hax));
            vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf + 2] = v;
            vidbuf[idx_buf + 3]=255;
            k++;
          }
          else
          {
            vidbuf[k] = vidbuf[k+1] = vidbuf[k+2]= v; 
            vidbuf[k+3]=255;
            k+=4;
          }
        }
        for (i=ic-curs;i<=ic+curs;i++) {
          if (all3flag) k = 4*(i/2*xdim/2+imc/2 + i/2*hax);
          else          k = 4*(i*xdim+imc);
          vidbuf[k] = NUMVALS-1 ; vidbuf[k+1] = vidbuf[k+2] = 0;
          vidbuf[k+3]=255;
        }
        for (imnr=imc-curs;imnr<=imc+curs;imnr++) {
          if (all3flag) k = 4*(ic/2*xdim/2+imnr/2 + ic/2*hax);
          else          k = 4*(ic*xdim+imnr);
          vidbuf[k] = NUMVALS-1 ; vidbuf[k+1] = vidbuf[k+2] = 0;
          vidbuf[k+3]=255;
        }
      }
      else {
        for (i=0;i<ydim;i++) {
          for (imnr=0;imnr<xdim;imnr++)
            if (all3flag && (i%2 || imnr%2)) continue;
          vidbuf[k++] = MAPOFFSET;
        }
      }
    }
    if (all3flag) {
      k = 0;
      for (i=0;i<ydim;i++)
      for (j=0;j<xdim;j++) {
        if (i%2 || j%2) continue;
        {
          idx_buf = 4*(hax + k + i/2*hax);
          vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf+2]= 
            sim[2][255-i/zf][j/zf]+MAPOFFSET;
          vidbuf[idx_buf+3] = 255;
          k++;
        }
      }
    }
    rectwrite(0,0,xdim-1,ydim-1,vidbuf);
  }
}

void
draw_second_image(int imc, int ic, int jc)
{
        int i,j,imnr,k;
        int hax,hay,curs;
  int idx_buf;
  GLubyte v;

        hax = xdim/2;
        hay = ydim/2;
        if (all3flag) curs = 15;
        else          curs = 2;

        if (maxflag && !all3flag)
        {
          k = 0;
          for (i=0;i<ydim;i++)
          for (j=0;j<xdim;j++)
          {
            if (plane==CORONAL)
              vidbuf[k] = sim2[0][255-i/zf][j/zf]+MAPOFFSET;
            else if (plane==HORIZONTAL)
              vidbuf[k] = sim2[1][i/zf][j/zf]+MAPOFFSET;
            else if (plane==SAGITTAL)
              vidbuf[k] = sim2[2][255-i/zf][j/zf]+MAPOFFSET;
            k++;
          }
        }
  else /*No max flag */        {
        if (plane==CORONAL || all3flag)
        {
          k = 0;
          for (i=0;i<ydim;i++)
          for (j=0;j<xdim;j++)
          {
            if (all3flag && (i%2 || j%2)) continue;
            if (imc/zf>=0&&imc/zf<imnr1)
              v = im2[imc/zf][(ydim-1-i)/zf][j/zf]+MAPOFFSET;
            else v=MAPOFFSET;
            if (truncflag)
              if (v<white_lolim+MAPOFFSET || v>white_hilim+MAPOFFSET)
                v=MAPOFFSET;
            if (i==ipt||j==jpt) v=255-(v-MAPOFFSET)+MAPOFFSET;
            if ((i==ic&&abs(j-jc)<=curs)||
                (j==jc&&abs(i-ic)<=curs)) v=0/*NUMVALS+MAPOFFSET*/;
            if (all3flag) 
      {
        idx_buf = 4*((xdim*hay) + k + ((i/2)*hax)) ;
        vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf + 2] = v;
        vidbuf[idx_buf + 3]=255;
        k++;
      }
            else
        {
        vidbuf[k] = vidbuf[k+1] = vidbuf[k+2]= v; 
        vidbuf[k+3]=255;
        k+=4;
        }
          }
        } 
        if (plane==HORIZONTAL || all3flag)
        {
          k = 0;
          for (imnr=0;imnr<ydim;imnr++)
          for (j=0;j<xdim;j++)
          {
            if (all3flag && (imnr%2 || j%2)) continue;
            if (imnr/zf>=0&&imnr/zf<imnr1)
              v = im2[imnr/zf][(ydim-1-ic)/zf][j/zf]+MAPOFFSET;
            else v=MAPOFFSET;
            if (truncflag)
              if (v<white_lolim+MAPOFFSET || v>white_hilim+MAPOFFSET)
                v=MAPOFFSET;
            if (imnr==impt||j==jpt) v=255-(v-MAPOFFSET)+MAPOFFSET;
            if ((imnr==imc&&abs(j-jc)<=curs)||
                (j==jc&&abs(imnr-imc)<=curs)) v=0/*NUMVALS+MAPOFFSET*/;
            if (all3flag)
      {
        idx_buf = 4*((xdim*hay) + hax + k + ((imnr/2)*hax)) ;
        vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf + 2] = v;
        vidbuf[idx_buf + 3]=255;
        k++;
      }
    else
        {
        vidbuf[k] = vidbuf[k+1] = vidbuf[k+2]= v; 
        vidbuf[k+3]=255;
        k+=4;
      }
          }
        } /*else*/
        if (plane==SAGITTAL || all3flag)
        {
          k = 0;
          for (i=0;i<ydim;i++)
          for (imnr=0;imnr<xdim;imnr++)
          {
            if (all3flag && (i%2 || imnr%2)) continue;
            if (imnr/zf>=0&&imnr/zf<imnr1)
              v = im2[imnr/zf][(ydim-1-i)/zf][jc/zf]+MAPOFFSET;
            else v=MAPOFFSET;
            if (truncflag)
              if (v<white_lolim+MAPOFFSET || v>white_hilim+MAPOFFSET)
                v=MAPOFFSET;
            if (imnr==impt||i==ipt) v=255-(v-MAPOFFSET)+MAPOFFSET;
            if ((imnr==imc&&abs(i-ic)<=curs)||
                (i==ic&&abs(imnr-imc)<=curs)) v=0/*NUMVALS+MAPOFFSET*/;
            if (all3flag) 
      {
        idx_buf = 4*(k + ((i/2)*hax));
        vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf + 2] = v;
        vidbuf[idx_buf + 3]=255;
        k++;
      }
    else
        {
        vidbuf[k] = vidbuf[k+1] = vidbuf[k+2]= v; 
        vidbuf[k+3]=255;
              k+=4;
      }
          }
        }
        if (all3flag) {
          k = 0;
          for (i=0;i<ydim;i++)
      for (j=0;j<xdim;j++) {
        if (i%2 || j%2) continue;
    {
      idx_buf = 4*(hax + k + i/2*hax);
      vidbuf[idx_buf] = vidbuf[idx_buf+1] = vidbuf[idx_buf+2]= 
        sim[2][255-i/zf][j/zf]+MAPOFFSET;
      vidbuf[idx_buf+3] = 255;
      k++;
    }
      }
  }
  rectwrite(0,0,xdim-1,ydim-1,vidbuf);
  }
}

int
read_surface(char *name)
{
  char fname[STRLEN], *cp ;

  cp = strchr(name, '/') ;
  if (!cp)  /* no path specified - put the path into it */
  {
    strcpy(surface, name) ;
    sprintf(sfname,"%s/%s/surf/%s",subjectsdir,pname,name);
    strcpy(fname, sfname) ;
  }
  else
    strcpy(fname, name) ;
  return(read_binary_surface(fname)) ;
}

int
read_binary_surface(char *fname)
{
#if 1
  if (mris)
    MRISfree(&mris) ;
  mris = MRISread(fname) ;
  if (!mris)
  {
    surfflag = FALSE;
    surfloaded = FALSE;
    curvflag = FALSE;
    curvloaded = FALSE;
    fieldsignflag = FALSE;
    fieldsignloaded = FALSE;
    return(ERROR_NOFILE) ;
  }
#else

  int k,n;                   /* loop counters */
  int ix,iy,iz;
  int version;
  int first,surfchanged=FALSE;
  int overtex_index, oface_index;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL) { printf("medit: ### File %s not found\n",fname);PR return(-1);}

  omris->nvertices = mris->nvertices;
  omris->nfaces = mris->nfaces;
  fread3(&first,fp);
  if (first == 16777215) {
    version = -1;
    printf("medit: new surface file format\n");
  }
  else {
    rewind(fp);
    version = 0;
    printf("medit: old surface file format\n");
  }
  fread3(&mris->nvertices,fp);
  fread3(&mris->nfaces,fp);

  printf("medit: vertices=%d, faces=%d\n",mris->nvertices,mris->nfaces);
  if (surfloaded) {
    if (omris->nvertices!=mris->nvertices || omris->nfaces!=mris->nfaces) {
      for (k=0;k<omris->nvertices;k++)
        free((int *)mris->vertices[k].f);
      free((vertex_type *)vertex);
      free((face_type *)face);
      surfchanged = TRUE;
    }
  }
  if (!surfloaded || surfchanged) {
    vertex = (vertex_type *)lcalloc(mris->nvertices,sizeof(vertex_type));
    face = (face_type *)lcalloc(mris->nfaces,sizeof(face_type));
  }

  for (k=0;k<mris->nvertices;k++)
  {
    fread2(&ix,fp);
    fread2(&iy,fp);
    fread2(&iz,fp);
    mris->vertices[k].x = ix/100.0;
    mris->vertices[k].y = iy/100.0;
    mris->vertices[k].z = iz/100.0;
    if (version==0)
    {
      fread1(&mris->vertices[k].num,fp);
      if (!surfloaded || surfchanged)
        mris->vertices[k].f = (int *)lcalloc(mris->vertices[k].num,sizeof(int));
      for (n=0;n<mris->vertices[k].num;n++)
        fread3(&mris->vertices[k].f[n],fp);
    } else mris->vertices[k].num = 0;
  }
  for (k=0;k<mris->nfaces;k++)
  {
    for (n=0;n<4;n++)
    {
      fread3(&mris->faces[k].v[n],fp);
      if (version<0)
        mris->vertices[mris->faces[k].v[n]].num++;
    }
  }
  fclose(fp);
  if (version<0)
  {
    for (k=0;k<mris->nvertices;k++)
    {
      if (!surfloaded || surfchanged)
        mris->vertices[k].f = (int *)lcalloc(mris->vertices[k].num,sizeof(int));
      mris->vertices[k].num = 0;
    }
    for (k=0;k<mris->nfaces;k++)
    {
      for (n=0;n<4;n++)
      {
        mris->vertices[mris->faces[k].v[n]].f[mris->vertices[mris->faces[k].v[n]].num++] = k;
      }
    }
  }
  for (k=0;k<mris->nvertices;k++) {   /* fscontour */
    mris->vertices[k].fieldsign = 0;
    mris->vertices[k].fsmask = 1;
    mris->vertices[k].curv = 0;
  }
#endif
  surfflag = TRUE;
  surfloaded = TRUE;
  curvflag = FALSE;
  curvloaded = FALSE;
  fieldsignflag = FALSE;
  fieldsignloaded = FALSE;
  return(0);
}
void
draw_surface(void)
{
  int k,n,ln,num;
  face_type *f;
  vertex_type *v,*lv;
  float h,xc,yc,zc,x,y,z,tx,ty,tz,lx,ly,lz,d,ld,vc[10][2];
  float fsf = zf;


  if (all3flag) linewidth(surflinewidth/2.0);
  else          linewidth(surflinewidth);

  mapcolor(NUMVALS+MAPOFFSET+1,  64,  64,  64);  /* undef gray */
  mapcolor(NUMVALS+MAPOFFSET+2,   0,   0, 255);  /* fs blue */
  mapcolor(NUMVALS+MAPOFFSET+3, 220, 220,   0);  /* fs yellow */

  mapcolor(NUMVALS+MAPOFFSET+4, 255,   0,   0);  /* curv red */
  mapcolor(NUMVALS+MAPOFFSET+5,   0, 255,   0);  /* curv green */

  if (plane==CORONAL || all3flag)
  {
    ortho2(0.0,xnum-1,0.0,ynum+1);
    yc = yy0+st*imc/fsf;
    for (k=0;k<mris->nfaces;k++)
    {
      f = &mris->faces[k];
      ln = VERTICES_PER_FACE-1;
      ld = mris->vertices[f->v[ln]].y-yc;
      num = 0;
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        v = &mris->vertices[f->v[n]];
        y = v->y;
        d = y-yc;
        if (d*ld<=0)
        {
          lv = &mris->vertices[f->v[ln]];
          x = v->x;
          z = v->z;
          lx = lv->x;
          ly = lv->y;
          lz = lv->z;
          h = (y-ly!=0)?(yc-ly)/(y-ly):1.0;
          tx = lx+h*(x-lx);
          tz = lz+h*(z-lz);
          vc[num][0] = (xx1-tx)/ps;
          vc[num][1] = ynum-((zz1-tz)/ps);
          if (all3flag) { 
            vc[num][0] = vc[num][0]/2.0;
            vc[num][1] = (vc[num][1] + (float)ynum + 1.0)/2.0;
          }
          if (fieldsignflag && fieldsignloaded) {  /* fscontour */
            if (v->fsmask<fsthresh)     color(NUMVALS+MAPOFFSET+1);
            else if (v->fieldsign>0.0)  color(NUMVALS+MAPOFFSET+2);
            else if (v->fieldsign<0.0)  color(NUMVALS+MAPOFFSET+3);
            else                        color(NUMVALS+MAPOFFSET+1); /*nondef*/
          }
          else if (curvflag && curvloaded) {
            if (v->curv>0.0)  color(NUMVALS+MAPOFFSET+4); /* TODO:colormap */
            else              color(NUMVALS+MAPOFFSET+5);
          }
          else /* color(YELLOW); */ /* just surface */
      glColor3f(1.0,1.0,0.0);
          num++;
        }
        ld = d;
        ln = n;
      }
      if (num>0)
      {
        bgnline();
        for (n=0;n<num;n++)
        {
          v2f(&vc[n][0]);
        }
        endline();
      }
    }
    for (k=0;k<npts;k++)
    {
      x = ptx[k]*tm[0][0]+pty[k]*tm[0][1]+ptz[k]*tm[0][2]+tm[0][3];
      y = ptx[k]*tm[1][0]+pty[k]*tm[1][1]+ptz[k]*tm[1][2]+tm[1][3];
      z = ptx[k]*tm[2][0]+pty[k]*tm[2][1]+ptz[k]*tm[2][2]+tm[2][3];
      if (fabs(y-yc)<st)
      {
        vc[0][0] = (xx1-x)/ps;
        vc[0][1] = ynum-((zz1-z)/ps);
        bgnpoint();
        color(GREEN);
        v2f(&vc[0][0]);
        endpoint();
      }
    }
  }
  if (plane==SAGITTAL || all3flag)
  {
    ortho2(0.0,xnum-1,0.0,ynum+1);
    xc = xx1-ps*jc/fsf;
    for (k=0;k<mris->nfaces;k++)
    {
      f = &mris->faces[k];
      ln = VERTICES_PER_FACE-1 ;
      ld = mris->vertices[f->v[ln]].x-xc;
      num = 0;
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        v = &mris->vertices[f->v[n]];
        x = v->x;
        d = x-xc;
        if (d*ld<=0)
        {
          lv = &mris->vertices[f->v[ln]];
          y = v->y;
          z = v->z;
          lx = lv->x;
          ly = lv->y;
          lz = lv->z;
          h = (x-lx!=0)?(xc-lx)/(x-lx):1.0;
          ty = ly+h*(y-ly);
          tz = lz+h*(z-lz);
          vc[num][0] = (ty-yy0)/st;
          vc[num][1] = ynum-((zz1-tz)/ps);
          if (all3flag) { 
            vc[num][0] = vc[num][0]/2.0;
            vc[num][1] = vc[num][1]/2.0;
          }
          if (fieldsignflag && fieldsignloaded) {  /* fscontour */
            if (v->fsmask<fsthresh)     color(NUMVALS+MAPOFFSET+1);
            else if (v->fieldsign>0.0)  color(NUMVALS+MAPOFFSET+2);
            else if (v->fieldsign<0.0)  color(NUMVALS+MAPOFFSET+3);
            else                        color(NUMVALS+MAPOFFSET+1);
          }
          else if (curvflag && curvloaded) {
            if (v->curv>0.0)  color(NUMVALS+MAPOFFSET+4);
            else              color(NUMVALS+MAPOFFSET+5);
          }
          else color(YELLOW);
          num++;
        }
        ld = d;
        ln = n;
      }
      if (num>0)
      {
        bgnline();
        for (n=0;n<num;n++)
        {
          v2f(&vc[n][0]);
        }
        endline();
      }
    }
    for (k=0;k<npts;k++)
    {
      x = ptx[k]*tm[0][0]+pty[k]*tm[0][1]+ptz[k]*tm[0][2]+tm[0][3];
      y = ptx[k]*tm[1][0]+pty[k]*tm[1][1]+ptz[k]*tm[1][2]+tm[1][3];
      z = ptx[k]*tm[2][0]+pty[k]*tm[2][1]+ptz[k]*tm[2][2]+tm[2][3];
      if (fabs(x-xc)<ps)
      {
        vc[0][0] = (y-yy0)/st;
        vc[0][1] = ynum-((zz1-z)/ps);
        bgnpoint();
        color(GREEN);
        v2f(&vc[0][0]);
        endpoint();
      }
    }
  }
  if (plane==HORIZONTAL || all3flag)
  {
    ortho2(0.0,xnum-1,0.0,ynum-1);
    zc = zz1-ps*(255.0-ic/fsf);
    for (k=0;k<mris->nfaces;k++)
    {
      f = &mris->faces[k];
      ln = VERTICES_PER_FACE-1 ;
      ld = mris->vertices[f->v[ln]].z-zc;
      num = 0;
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        v = &mris->vertices[f->v[n]];
        z = v->z;
        d = z-zc;
        if (d*ld<=0)
        {
          lv = &mris->vertices[f->v[ln]];
          x = v->x;
          y = v->y;
          lx = lv->x;
          ly = lv->y;
          lz = lv->z;
          h = (z-lz!=0)?(zc-lz)/(z-lz):1.0;
          tx = lx+h*(x-lx);
          ty = ly+h*(y-ly);
          vc[num][0] = (xx1-tx)/ps;
          vc[num][1] = ((ty-yy0)/st);
          if (all3flag) { 
            vc[num][0] = (vc[num][0] + (float)xnum - 1.0)/2.0;
            vc[num][1] = (vc[num][1] + (float)ynum - 1.0)/2.0;
          }
          if (fieldsignflag && fieldsignloaded) {  /* fscontour */
            if (v->fsmask<fsthresh)     color(NUMVALS+MAPOFFSET+1);
            else if (v->fieldsign>0.0)  color(NUMVALS+MAPOFFSET+2);
            else if (v->fieldsign<0.0)  color(NUMVALS+MAPOFFSET+3);
            else                        color(NUMVALS+MAPOFFSET+1);
          }
          else if (curvflag && curvloaded) {
            if (v->curv>0.0)  color(NUMVALS+MAPOFFSET+4);
            else              color(NUMVALS+MAPOFFSET+5);
          }
          else color(YELLOW);
          num++;
        }
        ld = d;
        ln = n;
      }
      if (num>0)
      {
        bgnline();
        for (n=0;n<num;n++)
        {
          v2f(&vc[n][0]);
        }
        endline();
      }
    }
    for (k=0;k<npts;k++)
    {
      x = ptx[k]*tm[0][0]+pty[k]*tm[0][1]+ptz[k]*tm[0][2]+tm[0][3];
      y = ptx[k]*tm[1][0]+pty[k]*tm[1][1]+ptz[k]*tm[1][2]+tm[1][3];
      z = ptx[k]*tm[2][0]+pty[k]*tm[2][1]+ptz[k]*tm[2][2]+tm[2][3];
      if (fabs(z-zc)<ps)
      {
        vc[0][0] = (xx1-x)/ps;
        vc[0][1] = ((y-yy0)/st);
        bgnpoint();
        color(GREEN);
        v2f(&vc[0][0]);
        endpoint();
      }
    }
  }

  /* TODO: cursor back on top */
  if (plane==CORONAL || all3flag) {
#ifndef OPENGL  /* OPENGL: cvidbuf pixels mangled?? jc,ic / zk?? */
    k = 0;
    for (i=ic-2;i<=ic+2;i++)
    for (j=jc-2;j<=jc+2;j++) {
      if ((i==ic&&abs(j-jc)<=2)||(j==jc&&abs(i-ic)<=2))
        cvidbuf[k] = NUMVALS+MAPOFFSET;
      else
        cvidbuf[k] = vidbuf[i*xdim+j];
      /*printf("k=%d cvidbuf[k]=%d\n",k,cvidbuf[k]);*/
      k++;
    }
    rectwrite((short)(jc-2),(short)(ic-2),
              (short)(jc+2),(short)(ic+2),cvidbuf);
#endif
#if 0  /* almost */
    color(NUMVALS+MAPOFFSET);
    linewidth(1);
    fjc = jc*255.0/256.0;
    fic = ic*257.0/255.0;
    bgnline();
      vc[0][0] = (fjc-2.0)/(float)zf;
      vc[0][1] = fic/(float)zf;
      v2f(vc[0]);
      vc[1][0] = (fjc+2.0)/(float)zf;
      vc[1][1] = fic/(float)zf;
      v2f(vc[1]);
    endline();
    bgnline();
      vc[0][0] = fjc/(float)zf;
      vc[0][1] = (fic-2.0)/(float)zf;
      v2f(vc[0]);
      vc[1][0] = fjc/(float)zf;
      vc[1][1] = (fic+2.0)/(float)zf;
      v2f(vc[1]);
    endline();
#endif
    ;
  }
  if (plane==SAGITTAL) {
    ;
  }
  if (plane==HORIZONTAL) {
    ;
  }

}


void
drawpts(void)
{
  int k;
  float xc,yc,zc,x,y,z;
  float v1[2],v2[2],v3[2],v4[3];
  float ptsize2=1;

  if (plane==CORONAL)
  {
    ortho2(0.0,xnum-1,0.0,ynum+1);
    yc = yy0+st*imc/fsf;
    for (k=npts-1;k>=0;k--)
    {
      x = ptx[k]*tm[0][0]+pty[k]*tm[0][1]+ptz[k]*tm[0][2]+tm[0][3];
      y = ptx[k]*tm[1][0]+pty[k]*tm[1][1]+ptz[k]*tm[1][2]+tm[1][3];
      z = ptx[k]*tm[2][0]+pty[k]*tm[2][1]+ptz[k]*tm[2][2]+tm[2][3];
      if (maxflag || fabs(y-yc)<st)
      {
        v1[0] = (xx1-x)/ps-ptsize2/ps;
        v1[1] = ynum-((zz1-z)/ps)-ptsize2/ps;
        v2[0] = (xx1-x)/ps+ptsize2/ps;
        v2[1] = ynum-((zz1-z)/ps)-ptsize2/ps;
        v3[0] = (xx1-x)/ps+ptsize2/ps;
        v3[1] = ynum-((zz1-z)/ps)+ptsize2/ps;
        v4[0] = (xx1-x)/ps-ptsize2/ps;
        v4[1] = ynum-((zz1-z)/ps)+ptsize2/ps;
        bgnpolygon();
        if (k<3)
          color(RED);
        else
          color(GREEN);
        v2f(v1);
        v2f(v2);
        v2f(v3);
        v2f(v4);
        endpolygon();
      }
    }
  }
  if (plane==SAGITTAL)
  {
    ortho2(0.0,xnum-1,0.0,ynum+1);
    xc = xx1-ps*jc/fsf;
    for (k=npts-1;k>=0;k--)
    {
      x = ptx[k]*tm[0][0]+pty[k]*tm[0][1]+ptz[k]*tm[0][2]+tm[0][3];
      y = ptx[k]*tm[1][0]+pty[k]*tm[1][1]+ptz[k]*tm[1][2]+tm[1][3];
      z = ptx[k]*tm[2][0]+pty[k]*tm[2][1]+ptz[k]*tm[2][2]+tm[2][3];
      if (maxflag || fabs(x-xc)<ps)
      {
        v1[0] = (y-yy0)/st-ptsize2/ps;
        v1[1] = ynum-((zz1-z)/ps)-ptsize2/ps;
        v2[0] = (y-yy0)/st+ptsize2/ps;
        v2[1] = ynum-((zz1-z)/ps)-ptsize2/ps;
        v3[0] = (y-yy0)/st+ptsize2/ps;
        v3[1] = ynum-((zz1-z)/ps)+ptsize2/ps;
        v4[0] = (y-yy0)/st-ptsize2/ps;
        v4[1] = ynum-((zz1-z)/ps)+ptsize2/ps;
        bgnpolygon();
        if (k<3)
          color(RED);
        else
          color(GREEN);
        v2f(v1);
        v2f(v2);
        v2f(v3);
        v2f(v4);
        endpolygon();
      }
    }
  }
  if (plane==HORIZONTAL)
  {
    ortho2(0.0,xnum-1,0.0,ynum-1);
    zc = zz1-ps*(255.0-ic/fsf);
    for (k=npts-1;k>=0;k--)
    {
      x = ptx[k]*tm[0][0]+pty[k]*tm[0][1]+ptz[k]*tm[0][2]+tm[0][3];
      y = ptx[k]*tm[1][0]+pty[k]*tm[1][1]+ptz[k]*tm[1][2]+tm[1][3];
      z = ptx[k]*tm[2][0]+pty[k]*tm[2][1]+ptz[k]*tm[2][2]+tm[2][3];
      if (maxflag || fabs(z-zc)<ps)
      {
        v1[0] = (xx1-x)/ps-ptsize2/ps;
        v1[1] = ((y-yy0)/st)-ptsize2/ps;
        v2[0] = (xx1-x)/ps+ptsize2/ps;
        v2[1] = ((y-yy0)/st)-ptsize2/ps;
        v3[0] = (xx1-x)/ps+ptsize2/ps;
        v3[1] = ((y-yy0)/st)+ptsize2/ps;
        v4[0] = (xx1-x)/ps-ptsize2/ps;
        v4[1] = ((y-yy0)/st)+ptsize2/ps;
        bgnpolygon();
        if (k<3)
          color(RED);
        else
          color(GREEN);
        v2f(v1);
        v2f(v2);
        v2f(v3);
        v2f(v4);
        endpolygon();
      }
    }
  }
}

void
write_dipoles(char *fname)
{
  int i,j,k,di,dj,dk;
  float x,y,z;
  FILE *fptr;
  int ***neighindex,index,nneighbors;

  ndip = 0;
/*
  ndip=2;
*/
  neighindex=(int ***)malloc((numimg/dip_spacing+1)*sizeof(int **));
  for (k=0;k<numimg;k+=dip_spacing)
  { 
    neighindex[k/dip_spacing]=(int **)malloc((ynum/dip_spacing+1)*sizeof(int *));
    for (i=0;i<ynum;i+=dip_spacing)
      neighindex[k/dip_spacing][i/dip_spacing]=(int *)malloc((xnum/dip_spacing+1)*
                                                             sizeof(int));
  }

  for (k=0;k<numimg;k+=dip_spacing)
  for (i=0;i<ynum;i+=dip_spacing)
  for (j=0;j<xnum;j+=dip_spacing)

  if (im[k][i][j]!=0)

/*
  if (sqrt(0.0+SQR(k-numimg/2)+SQR(i-ynum/2)+SQR(j-xnum/2))<=100)
*/
  {
    neighindex[k/dip_spacing][i/dip_spacing][j/dip_spacing] = ndip;
    ndip++;
  } else
   neighindex[k/dip_spacing][i/dip_spacing][j/dip_spacing] = -1;

  fptr = fopen(fname,"w");
  if (fptr==NULL) {printf("medit: ### can't create file %s\n",fname);PR return;}
  fprintf(fptr,"#!ascii\n");
  fprintf(fptr,"%d\n",ndip);
/*
  fprintf(fptr,"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
          30.0,-75.0,40.0,0.0,0.0,0.0);
  fprintf(fptr,"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
          45.0,10.0,70.0,0.0,0.0,0.0);
*/
  for (k=0;k<numimg;k+=dip_spacing)
  for (i=0;i<ynum;i+=dip_spacing)
  for (j=0;j<xnum;j+=dip_spacing)
  if (neighindex[(k)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing]>=0)
  {
    x = xx1-ps*j;
    y = yy0+st*k;
    z = zz1-ps*i;
    fprintf(fptr,"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",x,y,z,0.0,0.0,0.0);
  }
/*
  fprintf(fptr,"%d \n",0);
  fprintf(fptr,"%d \n",0);
*/
  for (k=0;k<numimg;k+=dip_spacing)
  for (i=0;i<ynum;i+=dip_spacing)
  for (j=0;j<xnum;j+=dip_spacing)
  if (neighindex[(k)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing]>=0)
  {
    nneighbors = 0;
    for (dk= -dip_spacing;dk<=dip_spacing;dk+=2*dip_spacing)
      if (dk>0)
      if (neighindex[(k+dk)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing]>=0)
        nneighbors++;
    for (di= -dip_spacing;di<=dip_spacing;di+=2*dip_spacing)
      if (di>0)
      if (neighindex[(k)/dip_spacing][(i+di)/dip_spacing][(j)/dip_spacing]>=0)
        nneighbors++;
    for (dj= -dip_spacing;dj<=dip_spacing;dj+=2*dip_spacing)
      if (dj>0)
      if (neighindex[(k)/dip_spacing][(i)/dip_spacing][(j+dj)/dip_spacing]>=0)
        nneighbors++;
    fprintf(fptr,"%d ",nneighbors);
    for (dk= -dip_spacing;dk<=dip_spacing;dk+=2*dip_spacing)
      if (dk>0)
      if ((index=neighindex[(k+dk)/dip_spacing][(i)/dip_spacing][(j)/dip_spacing])>=0)
        fprintf(fptr,"%d ",index);
    for (di= -dip_spacing;di<=dip_spacing;di+=2*dip_spacing)
      if (di>0)
      if ((index=neighindex[(k)/dip_spacing][(i+di)/dip_spacing][(j)/dip_spacing])>=0)
        fprintf(fptr,"%d ",index);
    for (dj= -dip_spacing;dj<=dip_spacing;dj+=2*dip_spacing)
      if (dj>0)
      if ((index=neighindex[(k)/dip_spacing][(i)/dip_spacing][(j+dj)/dip_spacing])>=0)
        fprintf(fptr,"%d ",index);
    fprintf(fptr,"\n");
  }
  fclose(fptr);
  printf("medit: dipole file %s written\n",fname);PR
}

void
write_decimation(char *fname)
{
  int k;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL) {printf("medit: ### can't create file %s\n",fname);PR return;}
  fprintf(fptr,"#!ascii\n");
  fprintf(fptr,"%d\n",ndip);
  for (k=0;k<ndip;k++)
  {
    fprintf(fptr,"1\n");
  }
  fclose(fptr);
  printf("medit: decimation file %s written\n",fname);PR
}

void
read_hpts(char *fname)      /* sprintf(fname,"%s","../bem/head2mri.hpts"); */
{
  FILE *fp;
  char line[NAME_LENGTH];
  int i;


  fp = fopen(fname,"r");
  if (fp==NULL) {printf("medit: ### File %s not found\n",fname);PR return;}

  ptsflag = TRUE;
  fgets(line,NAME_LENGTH,fp);
  i = 0;
  while (!feof(fp))
  {
    jold[i] = iold[i] = imold[i] = 0;
    if (line[0]=='2') /* .show format? */
      sscanf(line,"%*s%*s%*s%*s%*s%*s%*s%f%*s%f%*s%f",
             &ptx[i],&pty[i],&ptz[i]);
    else /* .apos format */
    {
      sscanf(line,"%*s%*s%f%f%f",&ptx[i],&pty[i],&ptz[i]);
      ptx[i] *= 1000; /* convert from m to mm */
      pty[i] *= 1000;
      ptz[i] *= 1000;
    }
    fgets(line,NAME_LENGTH,fp);
    i++;
  }
  printf("head points file %s read\n",fname);
  npts = i;
  fclose(fp);
}

void
read_htrans(char *fname)   /* sprintf(fname,"%s","../bem/head2mri.trans"); */
{
  int i,j,k;
  FILE *fp;
  fp = fopen(fname,"r");
  if (fp==NULL) {
    printf("medit: no existing htrans file: making default\n");PR
    for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      tm[i][j] = (i==j);
  }
  else {
    for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      fscanf(fp,"%f",&tm[i][j]);
    printf("transformation file %s read\n",fname);PR
  }
  for (k=0;k<MAXPARS;k++)
    par[k] = dpar[k] = 0;
}

void
write_htrans(char *fname)   /* sprintf(fname,"%s","../bem/head2mri.trans"); */
{
  int i,j;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL) {printf("medit: can't create file %s\n",fname);return;}
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      fprintf(fptr,"%13.6e ",tm[i][j]);
    fprintf(fptr,"\n");
  }
  fclose(fptr);
  printf("transformation file %s written\n",fname);
}

void
make_filenames(char *lsubjectsdir,char *lsrname, char *lpname, char *limdir,
               char *lsurface)
{
  char tmpsurface[NAME_LENGTH], surfhemi[NAME_LENGTH];
    
  subjectsdir = (char *)malloc(NAME_LENGTH*sizeof(char)); /* malloc for tcl */
  srname = (char *)malloc(NAME_LENGTH*sizeof(char));
  pname = (char *)malloc(NAME_LENGTH*sizeof(char));
  imtype = (char *)malloc(NAME_LENGTH*sizeof(char));
  imtype2 = (char *)malloc(NAME_LENGTH*sizeof(char));
  surface = (char *)malloc(NAME_LENGTH*sizeof(char));
  mfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  sfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  tfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  dipfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  decfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  hpfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  htfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  sgfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  fsfname = (char *)malloc(NAME_LENGTH*sizeof(char)); /* fscontour */
  fmfname = (char *)malloc(NAME_LENGTH*sizeof(char)); /* fscontour */
  cfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  xffname = (char *)malloc(NAME_LENGTH*sizeof(char));
  rfname = (char *)malloc(NAME_LENGTH*sizeof(char));

  strcpy(subjectsdir,lsubjectsdir);
  strcpy(srname,lsrname);
  strcpy(pname,lpname);
  strcpy(imtype,limdir);
  strcpy(surface,lsurface);
  sprintf(mfname,"%s/%s/mri/%s/COR-",subjectsdir,pname,imtype);
  sprintf(sfname,"%s/%s/surf/%s",subjectsdir,pname,surface);
  sprintf(tfname,"%s/%s/%s",subjectsdir,pname,TMP_DIR);
  /* next relative to session */
  sprintf(dipfname,"%s/bem/brain3d%d.dip",srname,dip_spacing);
  sprintf(decfname,"%s/bem/brain3d%d.dec",srname,dip_spacing);
  sprintf(hpfname,"%s/bem/head2mri.hpts",srname);
  sprintf(htfname,"%s/bem/head2mri.trans",srname);
  sprintf(sgfname,"/tmp/tkmedit.rgb");
  strcpy(tmpsurface,surface);  /* strtok writes in arg it parses! */
  strcpy(surfhemi,strtok(tmpsurface,"."));
  sprintf(fsfname,"%s/fs/%s.fs",srname,surfhemi);
  sprintf(fmfname,"%s/fs/%s.fm",srname,surfhemi);
  sprintf(cfname,"%s/%s/surf/%s.curv",subjectsdir,pname,surfhemi);
  sprintf(xffname,"%s/%s/%s/%s",
                     subjectsdir,pname,TRANSFORM_DIR,TALAIRACH_FNAME);
}

void
mirror(void)
{
  flip_corview_xyz("-x","y","z");
}

void
read_fieldsign(char *fname)  /* fscontour */
{
  int k,vnum;
  float f;
  FILE *fp;

  if (!surfloaded) {printf("medit: ### surface %s not loaded\n",fname);PR return;}

  fp = fopen(fname,"r");
  if (fp==NULL) {printf("medit: ### File %s not found\n",fname);PR return;}
  fread(&vnum,1,sizeof(int),fp);
  printf("medit: mris->nvertices = %d, vnum = %d\n",mris->nvertices,vnum);
  if (vnum!=mris->nvertices) {
    printf("medit: ### incompatible vertex number in file %s\n",fname);PR
    return; }
  for (k=0;k<vnum;k++)
  {
    fread(&f,1,sizeof(float),fp);
    mris->vertices[k].fieldsign = f;
  }
  fclose(fp);
  fieldsignloaded = TRUE;
  fieldsignflag = TRUE;
  surflinewidth = 3;
  PR
}

void
read_fsmask(char *fname)  /* fscontour */
{
  int k,vnum;
  float f;
  FILE *fp;

  if (!surfloaded) {printf("medit: ### surface %s not loaded\n",fname);PR return;}

  printf("medit: read_fsmask(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL) {printf("medit: ### File %s not found\n",fname);PR return;}
  fread(&vnum,1,sizeof(int),fp);
  if (vnum!=mris->nvertices) {
    printf("medit: ### incompatible vertex number in file %s\n",fname);PR
    return; }
  for (k=0;k<vnum;k++)
  {
    fread(&f,1,sizeof(float),fp);
    mris->vertices[k].fsmask = f;
  }
  fclose(fp);
  PR
}


void
read_binary_curvature(char *fname)
{
  float curvmin, curvmax, curv;
  int   k;

  MRISreadBinaryCurvature(mris, fname) ;

  curvmin= 1000000.0f ; curvmax = -curvmin;
  for (k=0;k<mris->nvertices;k++)
  {
    curv = mris->vertices[k].curv;
    if (curv>curvmax) curvmax=curv;
    if (curv<curvmin) curvmin=curv;
  }
  printf("medit: curvature read: min=%f max=%f\n",curvmin,curvmax);PR
  curvloaded = TRUE;
  curvflag = TRUE;
}

void
smooth_3d(int niter)
{
  int iter,i,j,k;
  int i2,j2,k2;
  int n;
  float sum;

  if (!dummy_im_allocated) {
    for (k=0;k<numimg;k++) {
      dummy_im[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
      for (i=0;i<ynum;i++) {
        dummy_im[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
      }
    }
    dummy_im_allocated = TRUE;
  }

  for (iter=0;iter<niter;iter++) {  /* smooth */
    for (j=1;j<xnum-1;j++) {
      printf(".");
      for (k=1;k<numimg-1;k++)
      for (i=1;i<ynum-1;i++) {
        sum = 0;
        n = 0;
        for(k2=k-1;k2<k+2;k2++)  /* 27 */
        for(i2=i-1;i2<i+2;i2++)
        for(j2=j-1;j2<j+2;j2++) {
          sum += im[k2][i2][j2];
          n++;
        }
        /*dummy_im[k][i][j] = (sum+im[k][i][j])/(float)(n+1);*/
        dummy_im[k][i][j] = (sum + 27*im[k][i][j])/(float)(n+27);
      }
      fflush(stdout);
    }
  }

  /* update view */
  for (k=0;k<numimg;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
    im[k][i][j] = dummy_im[k][i][j];
  for (k=0;k<numimg;k++)
    changed[k] = TRUE;
  editedimage = TRUE;

  for (k=0;k<6;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
    sim[k][i][j] = 0;
  for (k=0;k<numimg;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    if (im[k][i][j]/2>sim[3][i][j]) sim[3][i][j] = im[k][i][j]/2;
    if (im[k][i][j]/2>sim[4][k][j]) sim[4][k][j] = im[k][i][j]/2;
    if (im[k][i][j]/2>sim[5][i][k]) sim[5][i][k] = im[k][i][j]/2;
  }
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
  for (k=0;k<3;k++)
    sim[k][i][j] = sim[k+3][i][j];

  printf("\nmedit: finished %d smooth_3d steps\n",niter);PR
}

  /* args: [-]{x,y,z} [-]{x,y,z} [-]{x,y,z} */
void
flip_corview_xyz(char *newx, char *newy, char *newz)
{
  int fx,fy,fz;
  int lx,ly,lz;
  int xx,xy,xz,yx,yy,yz,zx,zy,zz;
  int i,j,k;
  int ni,nj,nk;

  ni = nj = nk = 0;
  fx=(newx[0]=='-')?1:0;
  fy=(newy[0]=='-')?1:0;
  fz=(newz[0]=='-')?1:0;
  lx = strlen(newx)-1;
  ly = strlen(newy)-1;
  lz = strlen(newz)-1;
  if(newx[lx]==newy[ly] || newx[lx]==newz[lz] || newy[ly]==newz[lz]) {
    printf("medit: ### degenerate flip\n");PR return; }
  xx=(newx[lx]=='x')?1:0; xy=(newx[lx]=='y')?1:0; xz=(newx[lx]=='z')?1:0;
  yx=(newy[ly]=='x')?1:0; yy=(newy[ly]=='y')?1:0; yz=(newy[ly]=='z')?1:0;
  zx=(newz[lz]=='x')?1:0; zy=(newz[lz]=='y')?1:0; zz=(newz[lz]=='z')?1:0;

  if (!dummy_im_allocated) {
    for (k=0;k<numimg;k++) {
      dummy_im[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
      for (i=0;i<ynum;i++) {
        dummy_im[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
      }
    }
    dummy_im_allocated = TRUE;
  }

  for (k=0;k<numimg;k++) {
    printf(".");
    for (i=0;i<ynum;i++)
    for (j=0;j<xnum;j++) {
      if(xx) nj=j; if(xy) nj=i; if(xz) nj=k;
      if(yx) ni=j; if(yy) ni=i; if(yz) ni=k;
      if(zx) nk=j; if(zy) nk=i; if(zz) nk=k;
      if(fx) nj=255-nj;
      if(fy) ni=255-ni;
      if(fz) nk=255-nk;
      dummy_im[k][i][j] = im[nk][ni][nj];
    }
    fflush(stdout);
  }
  printf("\n");

  /* update view */
  for (k=0;k<numimg;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
    im[k][i][j] = dummy_im[k][i][j];
  for (k=0;k<numimg;k++)
    changed[k] = TRUE;
  editedimage = TRUE;

  for (k=0;k<6;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
    sim[k][i][j] = 0;
  for (k=0;k<numimg;k++)
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++) {
    if (im[k][i][j]/2>sim[3][i][j]) sim[3][i][j] = im[k][i][j]/2;
    if (im[k][i][j]/2>sim[4][k][j]) sim[4][k][j] = im[k][i][j]/2;
    if (im[k][i][j]/2>sim[5][i][k]) sim[5][i][k] = im[k][i][j]/2;
  }
  for (i=0;i<ynum;i++)
  for (j=0;j<xnum;j++)
  for (k=0;k<3;k++)
    sim[k][i][j] = sim[k+3][i][j];

  printf("medit: imageflip done--to overwrite current COR: write_images\n");PR
  redraw();
}

void
wmfilter_corslice(int imc)
{
  char *getenv(),*mri_dir,fname[NAME_LENGTH];
  int nver,ncor;
  int i,j,k,k2,di,dj,dk,m,n,u,ws2=2,maxi,mini;
  float xcor[MAXCOR],ycor[MAXCOR],zcor[MAXCOR];
  float fwhite_hilim=140,fwhite_lolim=80,fgray_hilim=100;
  float numvox,numnz,numz;
  float f,f2,a,b,c,x,y,z,s;
  float cfracz = 0.6;
  float cfracnz = 0.6;
  double sum2,sum,var,avg,tvar,maxvar,minvar;
  double sum2v[MAXLEN],sumv[MAXLEN],avgv[MAXLEN],varv[MAXLEN],nv[MAXLEN];
  FILE *fptr;

  k2 = imc/zf;
  if (plane!=CORONAL) {
    printf("medit: ### can only wmfilter CORONAL slice\n");PR return; }
  if(k2<ws2 || k2>numimg-ws2) {
    printf("medit: ### slice too close to edge\n");PR return; }
  mri_dir = getenv("MRI_DIR");
  if (mri_dir==NULL) {
    printf("medit: ### env var MRI_DIR undefined (setenv, restart)\n");return;}
  if (!flossflag && !spackleflag) { 
    printf("medit: ### no spackle or floss  ...skipping wmfilter\n"); return; }

  sprintf(fname,"%s/lib/bem/%s",mri_dir,DIR_FILE);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {printf("medit: ### File %s not found\n",fname);PR return;}
  fscanf(fptr,"%d",&nver);
  for (i=0,ncor=0;i<nver;i++) {
    fscanf(fptr,"%*d %f %f %f",&x,&y,&z);
    for (j=0;j<ncor;j++)
      if ((x==-xcor[j]) && (y==-ycor[j]) && (z==-zcor[j])) goto L1;
    xcor[ncor] = x;
    ycor[ncor] = y;
    zcor[ncor] = z;
    ncor++;
  L1:;
  }
  fwhite_hilim = (float)white_hilim;
  fwhite_lolim = (float)white_lolim;
  fgray_hilim = (float)gray_hilim;

  if (!wmfilter_ims_allocated) {
    for (k=0;k<numimg;k++) {
      fill[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
      im_b[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
      for (i=0;i<ynum;i++) {
        fill[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
        im_b[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
      }
    }
    wmfilter_ims_allocated = TRUE;
  }

  if (!second_im_allocated)
    alloc_second_im();

  /*for (k=0;k<numimg-1;k++)*/
  for (k=k2-ws2;k<=k2+ws2;k++)
  for (i=0;i<ynum-1;i++)
  for (j=0;j<xnum-1;j++) {
    im_b[k][i][j] = im[k][i][j];
    im2[k][i][j] = im[k][i][j];
    if (im_b[k][i][j]>fwhite_hilim || im_b[k][i][j]<fwhite_lolim)
      im_b[k][i][j] = 0;
    fill[k][i][j] = im_b[k][i][j];
  }

  k = k2;
  printf("medit: %d unique orientations\n",ncor);PR
  printf("medit: white_lolim = %d   white_hilim = %d   gray_hilim = %d\n",
                 white_lolim,white_hilim,gray_hilim);PR
  printf("medit: wmfiltering coronal slice %d\n",k2);
  if (flossflag)   printf("medit:   floss => '.'\n");
  if (spackleflag) printf("medit:   spackle => '#'\n");
  printf("medit:    ...\n"); 

  for (i=ws2;i<ynum-1-ws2;i++)
  for (j=ws2;j<xnum-1-ws2;j++) {
    numvox = numnz = numz = 0;
    for (dk = -ws2;dk<=ws2;dk++)
    for (di = -ws2;di<=ws2;di++)
    for (dj = -ws2;dj<=ws2;dj++) {
      f = im_b[k+dk][i+di][j+dj];
      s = dk*c+di*b+dj*a;
      numvox++;
      if (f!=0) numnz++;
      else      numz++;
    }
    if ((im_b[k][i][j]==0 && numnz>=cfracnz*SQR(2*ws2+1)) ||
        (im_b[k][i][j]!=0 && numz>=cfracz*SQR(2*ws2+1))) {
      maxvar = -1000000;
      minvar = 1000000;
      maxi = mini = -1;
      for (m=0;m<ncor;m++) {
        a = xcor[m];
        b = ycor[m];
        c = zcor[m];
        sum = sum2 = n = 0;
        for (u=0;u<2*ws2+1;u++)
          sumv[u] = sum2v[u] = nv[u] = 0;
        for (dk = -ws2;dk<=ws2;dk++)
        for (di = -ws2;di<=ws2;di++)
        for (dj = -ws2;dj<=ws2;dj++) {
          u = ws2+floor(dk*c+di*b+dj*a+0.5);
          u = (u<0)?0:(u>2*ws2+1-1)?2*ws2+1-1:u;
          f = im_b[k+dk][i+di][j+dj];
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
        for (u=0;u<2*ws2+1;u++) {
          avgv[u] = sumv[u]/nv[u];
          varv[u] = sum2v[u]/nv[u]-avgv[u]*avgv[u];
          tvar += varv[u];
        }
        tvar /= (2*ws2+1);
        if (tvar>maxvar) {maxvar=tvar;maxi=m;}
        if (tvar<minvar) {minvar=tvar;mini=m;}
      }
      a = xcor[mini];
      b = ycor[mini];
      c = zcor[mini];

      numvox = numnz = numz = sum = sum2 = 0;
      for (dk = -ws2;dk<=ws2;dk++)
      for (di = -ws2;di<=ws2;di++)
      for (dj = -ws2;dj<=ws2;dj++) {
        f = im_b[k+dk][i+di][j+dj];
        f2 = im2[k+dk][i+di][j+dj];
        s = dk*c+di*b+dj*a;
        if (fabs(s)<=0.5) {
          numvox++;
          sum2 += f2;
          if (f!=0) { numnz++; sum += f; }
          else      { numz++; }
        }
      }
      if (numnz!=0) sum /= numnz;
      if (numvox!=0) sum2 /= numvox;
      f = im_b[k][i][j];
      f2 = im2[k][i][j];
      if (flossflag && f!=0 && numz/numvox>cfracz && f<=fgray_hilim) {
        f=0; printf("."); }
      else if (spackleflag && f==0 && numnz/numvox>cfracnz) {
        f=sum; printf("#"); }
      fill[k][i][j] = f;
    }
  }
  printf("\n");PR
  /*for (k=0;k<numimg-1;k++)*/
  k = k2;
  for (i=0;i<ynum-1;i++)
  for (j=0;j<xnum-1;j++) {
    im_b[k][i][j] = fill[k][i][j];
    if (im_b[k][i][j]>fwhite_hilim || im_b[k][i][j]<fwhite_lolim)
      im_b[k][i][j]=0;
  }

  /*for (k=0;k<numimg-1;k++)*/
  k = k2;
  for (k=k2-ws2;k<=k2+ws2;k++)
  for (i=0;i<ynum-1;i++)
  for (j=0;j<xnum-1;j++) {
    if (k==k2) im2[k][i][j] = im_b[k][i][j];
    else       im2[k][i][j] = 0;
  }
  printf("medit: wmfiltered slice put in 2nd image set (can't be saved)\n");PR
  drawsecondflag = TRUE;
  redraw();
}

void
norm_allslices(int normdir)
{
  int i,j,k;
  int x,y;
  float imf[256][256];
  float flim0,flim1,flim2,flim3;

  x = y = 0;
  if ((plane==CORONAL && normdir==POSTANT) ||
      (plane==SAGITTAL && normdir==RIGHTLEFT) ||
      (plane==HORIZONTAL && normdir==INFSUP)) {
    printf("medit: ### N.B.: norm gradient not visible in this view\n");
  }
  printf("medit: normalizing all slices...\n");
  flim0 = (float)lim0;
  flim1 = (float)lim1;
  flim2 = (float)lim2;
  flim3 = (float)lim3;
  for (k=0;k<numimg;k++) {
    printf(".");
    for (i=0;i<256;i++)
    for (j=0;j<256;j++) {
      if (normdir==POSTANT)   {x=i; y=255-k;}  /* SAG:x */
      if (normdir==INFSUP)    {x=j;     y=i;}  /* COR:y */
      if (normdir==RIGHTLEFT) {x=k; y=255-j;}  /* HOR:x */
      imf[y][x] = im[k][i][j];
      if ((255-y)<=flim0)
        imf[y][x]*=ffrac0;
      else if ((255-y)<=flim1)
        imf[y][x]*=(ffrac0+((255-y)-flim0)*(ffrac1-ffrac0)/(flim1-flim0));
      else if ((255-y)<=flim2)
        imf[y][x]*=(ffrac1+((255-y)-flim1)*(ffrac2-ffrac1)/(flim2-flim1));
      else if ((255-y)<=flim3)
        imf[y][x]*=(ffrac2+((255-y)-flim2)*(ffrac3-ffrac2)/(flim3-flim2));
      else
        imf[y][x]*=ffrac3;
      if (imf[y][x]>255) imf[y][x]=255;
      im[k][i][j] = floor(imf[y][x]+0.5);
    }
    fflush(stdout);
  }
  printf("\n");
  printf("medit: done (to undo: quit w/o SAVEIMG)\n");PR
  for (k=0;k<numimg;k++)
    changed[k] = TRUE;
  editedimage = TRUE;
  redraw();
}

void norm_slice(int imc, int ic,int jc, int normdir)
{
  int k,i,j;
  int x,y;
  int k0,k1,i0,i1,j0,j1;
  int imc0,ic0,jc0;
  float imf[256][256];
  float flim0,flim1,flim2,flim3;

  x = y = 0;
  k0 = k1 = i0 = i1 = j0 = j1 = 0;
  if ((plane==CORONAL && normdir==POSTANT) ||
      (plane==SAGITTAL && normdir==RIGHTLEFT) ||
      (plane==HORIZONTAL && normdir==INFSUP)) {
    printf("medit: ### N.B.: norm gradient not visible in this view\n");
  }

  if (!second_im_allocated)
    alloc_second_im();

  flim0 = (float)lim0;
  flim1 = (float)lim1;
  flim2 = (float)lim2;
  flim3 = (float)lim3;
  imc0 = imc/zf;
  ic0 = (ydim-1-ic)/zf;
  jc0 = jc/zf;
  if (plane==CORONAL)   {k0=imc0; k1=imc0+1; i0=0;   i1=256;   j0=0;  j1=256;}
  if (plane==SAGITTAL)  {k0=0;    k1=numimg; i0=0;   i1=256;   j0=jc0;j1=jc0+1;}
  if (plane==HORIZONTAL){k0=0;    k1=numimg; i0=ic0; i1=ic0+1; j0=0;  j1=256;}
  /*printf("k0=%d k1=%d i0=%d i1=%d j0=%d j1=%d\n",k0,k1,i0,i1,j0,j1);*/
  for (k=k0;k<k1;k++)
  for (i=i0;i<i1;i++)
  for (j=j0;j<j1;j++) {
    if (normdir==POSTANT)   {x=i; y=255-k;}  /* SAG:x */
    if (normdir==INFSUP)    {x=j;     y=i;}  /* COR:y */
    if (normdir==RIGHTLEFT) {x=k; y=255-j;}  /* HOR:x */
    imf[y][x] = im[k][i][j];
    if ((255-y)<=flim0)
      imf[y][x]*=ffrac0;
    else if ((255-y)<=flim1)
      imf[y][x]*=(ffrac0+((255-y)-flim0)*(ffrac1-ffrac0)/(flim1-flim0));
    else if ((255-y)<=flim2)
      imf[y][x]*=(ffrac1+((255-y)-flim1)*(ffrac2-ffrac1)/(flim2-flim1));
    else if ((255-y)<=flim3)
      imf[y][x]*=(ffrac2+((255-y)-flim2)*(ffrac3-ffrac2)/(flim3-flim2));
    else
      imf[y][x]*=ffrac3;
    if (imf[y][x]>255) imf[y][x]=255;
    im2[k][i][j] = floor(imf[y][x]+0.5);
  }
  printf("medit: normalized slice put in 2nd image set (can't be saved)\n");PR
  drawsecondflag = TRUE;
  redraw();
}

void
alloc_second_im(void)
{
  int i,k;

  for (k=0;k<numimg;k++) {
    im2[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++) {
      im2[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
  }
  for (k=0;k<6;k++) {
    sim2[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++) {
      sim2[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
  }
  second_im_allocated = TRUE;
}




/*----------------end medit.c -------------------------------*/

/* %%%%Sean %%% */


/* boilerplate wrap function defines for easier viewing */
#define WBEGIN (ClientData clientData,Tcl_Interp *interp,int argc,char *argv[]){
#define ERR(N,S)  if(argc!=N){interp->result=S;return TCL_ERROR;}
#define WEND   return TCL_OK;}
#define REND  (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL

/*=======================================================================*/
/* function wrappers and errors */
int                  W_open_window  WBEGIN
  ERR(1,"Wrong # args: open_window")
                       open_window(pname);  WEND

int                  W_resize_window_intstep  WBEGIN
  ERR(2,"Wrong # args: resize_window_intstep <mag=1,2,3>")
                       resize_window_intstep(atoi(argv[1]));  WEND

int                  W_move_window  WBEGIN
  ERR(3,"Wrong # args: move_window <x> <y>")
                       move_window(atoi(argv[1]),atoi(argv[2]));  WEND

int                  W_pop_gl_window WBEGIN
  ERR(1,"Wrong # args: pop_gl_window")
                       pop_gl_window();  WEND

int                  W_redraw  WBEGIN 
  ERR(1,"Wrong # args: redraw")
                       redraw();  WEND

int                  W_upslice WBEGIN 
  ERR(1,"Wrong # args: upslice")
                       upslice();  WEND

int                  W_downslice WBEGIN 
  ERR(1,"Wrong # args: downslice")
                       downslice();  WEND

int                  W_rotate_brain_x  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_x <deg>")
                       rotate_brain(atof(argv[1]),'x'); WEND

int                  W_rotate_brain_y  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_y <deg>")
                       rotate_brain(atof(argv[1]),'y'); WEND

int                  W_rotate_brain_z  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_z <deg>")
                       rotate_brain(atof(argv[1]),'z'); WEND

int                  W_translate_brain_x  WBEGIN
  ERR(2,"Wrong # args: translate_brain_x <mm>")
                       translate_brain(atof(argv[1]),'x');  WEND

int                  W_translate_brain_y  WBEGIN
  ERR(2,"Wrong # args: translate_brain_y <mm>")
                       translate_brain(atof(argv[1]),'y');  WEND

int                  W_translate_brain_z  WBEGIN
  ERR(2,"Wrong # args: translate_brain_z <mm>")
                       translate_brain(atof(argv[1]),'z');  WEND

int                  W_write_images  WBEGIN
  ERR(1,"Wrong # args: write_images")
                       write_images(mfname); WEND

int                  W_write_dipoles  WBEGIN
  ERR(1,"Wrong # args: write_dipoles")
                       write_dipoles(dipfname); WEND

int                  W_write_decimation  WBEGIN
  ERR(1,"Wrong # args: write_decimation")
                       write_decimation(decfname); WEND

int                  W_goto_vertex  WBEGIN 
  ERR(2,"Wrong # args: goto_vertex")
                       goto_vertex(atoi(argv[1]));  WEND

int                  W_mark_file_vertices  WBEGIN 
  ERR(2,"Wrong # args: mark_file_vertices")
                       mark_file_vertices(argv[1]);  WEND

int                  W_unmark_vertices  WBEGIN 
  ERR(1,"Wrong # args: unmark_vertices")
                       unmark_vertices();  WEND

int                  W_goto_point  WBEGIN
  ERR(1,"Wrong # args: goto_point")
                       goto_point(tfname); WEND

int                  W_goto_point_coords  WBEGIN
  ERR(4,"Wrong # args: goto_point_coords <inmr*zf> <i*zf> <j*zf>")
                       goto_point_coords((int)atof(argv[1]),(int)atof(argv[2]),
                                         (int)atof(argv[3])); WEND
int                  W_write_point  WBEGIN
  ERR(1,"Wrong # args: write_point")
                       write_point(tfname); WEND

int                  W_coords_to_talairach  WBEGIN
  ERR(1,"Wrong # args: coords_to_talairach")
                       coords_to_talairach(); WEND

int                  W_talairach_to_coords  WBEGIN
  ERR(1,"Wrong # args: talairach_to_coords")
                       talairach_to_coords(); WEND

int                  W_read_hpts  WBEGIN
  ERR(1,"Wrong # args: read_hpts")
                       read_hpts(hpfname); WEND

int                  W_read_htrans  WBEGIN
  ERR(1,"Wrong # args: read_htrans")
                       read_htrans(htfname); WEND

int                  W_write_htrans  WBEGIN
  ERR(1,"Wrong # args: write_htrans")
                       write_htrans(htfname); WEND

int                  W_save_rgb  WBEGIN
  ERR(1,"Wrong # args: save_rgb")
                       save_rgb(sgfname);  WEND

int                  W_set_scale  WBEGIN
  ERR(1,"Wrong # args: set_scale")
                       set_scale();  WEND

int                  W_mirror  WBEGIN
  ERR(1,"Wrong # args: mirror")
                       mirror();  WEND

int                  W_read_fieldsign  WBEGIN
  ERR(1,"Wrong # args: read_fieldsign")
                       read_fieldsign(fsfname);  WEND

int                  W_read_fsmask WBEGIN
  ERR(1,"Wrong # args: read_fsmask")
                       read_fsmask(fmfname);  WEND

int                  W_read_binary_curv  WBEGIN 
  ERR(1,"Wrong # args: read_binary_curv")
                       read_binary_curvature(cfname);  WEND

int                  W_smooth_3d  WBEGIN
  ERR(2,"Wrong # args: smooth_3d <steps>")
                       smooth_3d(atoi(argv[1]));  WEND

int                  W_flip_corview_xyz  WBEGIN
  ERR(4,"Wrong # args: flip_corview_xyz <[-]{x,y,z}> <[-]{x,y,z}> <[-]{x,y,z}>")
                       flip_corview_xyz(argv[1],argv[2],argv[3]);  WEND

int                  W_read_second_images  WBEGIN
  ERR(2,"Wrong # args: read_second_images <imtype2>     [T1,wm,filled,other]")
                       read_second_images(argv[1]);  WEND

int                  W_read_binary_surf  WBEGIN
  ERR(1,"Wrong # args: read_binary_surf")
                       read_binary_surface(sfname); WEND

int                  W_show_vertex  WBEGIN
  ERR(1,"Wrong # args: show_vertex")
                       show_vertex(); WEND

int                  W_dump_vertex  WBEGIN
  ERR(2,"Wrong # args: dump_vertex")
                       dump_vertex(atoi(argv[1])); WEND

int                  W_read_surface  WBEGIN
  ERR(2,"Wrong # args: read_surface <surface file>")
                       read_surface(argv[1]); WEND

int                  W_wmfilter_corslice  WBEGIN
  ERR(1,"Wrong # args: wmfilter_corslice")
                       wmfilter_corslice(imc); WEND

int                  W_norm_slice  WBEGIN
  ERR(2,"Wrong # args: norm_slice <dir: 0=PostAnt,1=InfSup,2=LeftRight>")
                       norm_slice(imc,ic,jc,atoi(argv[1])); WEND

int                  W_norm_allslices  WBEGIN
  ERR(2,"Wrong # args: norm_allslices <dir: 0=PostAnt,1=InfSup,2=LeftRight>")
                       norm_allslices(atoi(argv[1])); WEND

/*=======================================================================*/

/* for tcl/tk */
static void StdinProc _ANSI_ARGS_((ClientData clientData, int mask));
static void Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
static Tk_Window mainWindow;
static Tcl_Interp *interp;
static Tcl_DString command;
static int tty;

int main(argc, argv)   /* new main */
int argc;
char **argv;
{
  int code;
  int scriptok=FALSE;
  int cur_arg;
  static char *display = NULL;
  char tkmedit_tcl[NAME_LENGTH];
  char script_tcl[NAME_LENGTH];
  char *envptr, *getenv();
  FILE *fp, *fopen();

  /* get tkmedit tcl startup script location from environment */
  envptr = getenv("MRI_DIR");
  if (envptr==NULL) {
    printf("tkmedit: env var MRI_DIR undefined (use setenv)\n");
    printf("    [dir containing mri distribution]\n");
    exit(1);
  }
  sprintf(tkmedit_tcl,"%s/lib/tcl/%s",envptr,"tkmedit.tcl");
  if ((fp=fopen(tkmedit_tcl,"r"))==NULL) {
    printf("tkmedit: script %s not found\n",tkmedit_tcl);
    exit(1);
  }
  else fclose(fp);

  /* look for script: (1) cwd, (2) MRI_DIR/lib/tcl */
  strcpy(script_tcl,"newscript.tcl");   /* default if no script */
  for (cur_arg=1; cur_arg<argc; cur_arg++) {  /* just look for -tcl arg */
    if (!MATCH(argv[cur_arg],"-tcl"))
      continue;
    cur_arg++;
    if (cur_arg==argc) {
      printf("tkmedit: ### no script name given\n");
      exit(1);
    }
    if (argv[cur_arg][0]=='-') {
      printf("tkmedit: ### option is bad script name: %s\n",argv[cur_arg]);
      exit(1);
    }
    strcpy(script_tcl,argv[cur_arg]);
    fp = fopen(script_tcl,"r");  /* first look in cwd */
    if (fp==NULL) {
      sprintf(script_tcl,"%s/lib/tcl/%s",envptr,argv[cur_arg]);/* then libdir */      fp = fopen(script_tcl,"r");
      if (fp==NULL) {
        printf(
          "tkmedit: (1) File ./%s not found\n",argv[cur_arg]);
        printf(
          "         (2) File %s/lib/tcl/%s not found\n",envptr,argv[cur_arg]);
        exit(1);
      }
    }
    scriptok = TRUE;
  }

  /* start program, now as function; gl window not opened yet */
  printf("tkmedit: starting medit\n");
  Medit((ClientData) NULL, interp, argc, argv); /* event loop commented out */

  /* start tcl/tk; first make interpreter */
  interp = Tcl_CreateInterp();

  /* make main window (not displayed until event loop starts) */
  mainWindow = Tk_CreateMainWindow(interp, display, argv[0], "Tk");
  if (mainWindow == NULL) {
    fprintf(stderr, "%s\n", interp->result);
    exit(1); }

  /* set the "tcl_interactive" variable */
  tty = isatty(0);
  Tcl_SetVar(interp, "tcl_interactive", (tty) ? "1" : "0", TCL_GLOBAL_ONLY);
  if (tty) promptflag = TRUE;

  /* read tcl/tk internal startup scripts */
  if (Tcl_Init(interp) == TCL_ERROR) {
    fprintf(stderr, "Tcl_Init failed: %s\n", interp->result); }
  if (Tk_Init(interp)== TCL_ERROR) {
    fprintf(stderr, "Tk_Init failed: %s\n", interp->result); }

  /*=======================================================================*/
  /* register wrapped surfer functions with interpreter */
  Tcl_CreateCommand(interp, "open_window",        W_open_window,        REND);
 Tcl_CreateCommand(interp,"resize_window_intstep",W_resize_window_intstep,REND);
  Tcl_CreateCommand(interp, "move_window",        W_move_window,        REND);
  Tcl_CreateCommand(interp, "pop_gl_window",      W_pop_gl_window,      REND);
  Tcl_CreateCommand(interp, "redraw",             W_redraw,             REND);
  Tcl_CreateCommand(interp, "upslice",            W_upslice,            REND);
  Tcl_CreateCommand(interp, "downslice",          W_downslice,          REND);
  Tcl_CreateCommand(interp, "rotate_brain_x",     W_rotate_brain_x,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_y",     W_rotate_brain_y,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_z",     W_rotate_brain_z,     REND);
  Tcl_CreateCommand(interp, "translate_brain_x",  W_translate_brain_x,  REND);
  Tcl_CreateCommand(interp, "translate_brain_y",  W_translate_brain_y,  REND);
  Tcl_CreateCommand(interp, "translate_brain_z",  W_translate_brain_z,  REND);
  Tcl_CreateCommand(interp, "write_images",       W_write_images,       REND);
  Tcl_CreateCommand(interp, "write_dipoles",      W_write_dipoles,      REND);
  Tcl_CreateCommand(interp, "write_decimation",   W_write_decimation,   REND);
  Tcl_CreateCommand(interp, "goto_point",         W_goto_point,         REND);
  Tcl_CreateCommand(interp, "goto_point_coords",  W_goto_point_coords,  REND);
  Tcl_CreateCommand(interp, "write_point",        W_write_point,        REND);
  Tcl_CreateCommand(interp, "coords_to_talairach",W_coords_to_talairach,REND);
  Tcl_CreateCommand(interp, "talairach_to_coords",W_talairach_to_coords,REND);
  Tcl_CreateCommand(interp, "read_hpts",          W_read_hpts,          REND);
  Tcl_CreateCommand(interp, "read_htrans",        W_read_htrans,        REND);
  Tcl_CreateCommand(interp, "write_htrans",       W_write_htrans,       REND);
  Tcl_CreateCommand(interp, "save_rgb",           W_save_rgb,           REND);
  Tcl_CreateCommand(interp, "set_scale",          W_set_scale,          REND);
  Tcl_CreateCommand(interp, "mirror",             W_mirror,             REND);
  Tcl_CreateCommand(interp, "read_fieldsign",     W_read_fieldsign,     REND);
  Tcl_CreateCommand(interp, "read_fsmask",        W_read_fsmask,        REND);
  Tcl_CreateCommand(interp, "read_binary_curv",   W_read_binary_curv,   REND);
  Tcl_CreateCommand(interp, "goto_vertex",        W_goto_vertex,        REND);
  Tcl_CreateCommand(interp, "mark_file_vertices", W_mark_file_vertices, REND);
  Tcl_CreateCommand(interp, "unmark_vertices",    W_unmark_vertices,    REND);
  Tcl_CreateCommand(interp, "smooth_3d",          W_smooth_3d,          REND);
  Tcl_CreateCommand(interp, "flip_corview_xyz",   W_flip_corview_xyz,   REND);
  Tcl_CreateCommand(interp, "read_second_images", W_read_second_images, REND);
  Tcl_CreateCommand(interp, "show_vertex",        W_show_vertex,   REND);
  Tcl_CreateCommand(interp, "dump_vertex",        W_dump_vertex,   REND);
  Tcl_CreateCommand(interp, "read_binary_surf",   W_read_binary_surf,   REND);
  Tcl_CreateCommand(interp, "read_surface",       W_read_surface,   REND);
  Tcl_CreateCommand(interp, "wmfilter_corslice",  W_wmfilter_corslice,  REND);
  Tcl_CreateCommand(interp, "norm_slice",         W_norm_slice,         REND);
  Tcl_CreateCommand(interp, "norm_allslices",     W_norm_allslices,     REND);
  /*=======================================================================*/
  /***** link global BOOLEAN variables to tcl equivalents */
  Tcl_LinkVar(interp,"maxflag",(char *)&maxflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"surfflag",(char *)&surfflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"editflag",(char *)&editflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"fieldsignflag",(char *)&fieldsignflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"curvflag",(char *)&curvflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"drawsecondflag",(char *)&drawsecondflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"inplaneflag",(char *)&inplaneflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"linearflag",(char *)&linearflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"bwflag",(char *)&bwflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"truncflag",(char *)&truncflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"flossflag",(char *)&flossflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"spackleflag",(char *)&spackleflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"circleflag",(char *)&circleflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"promptflag",(char *)&promptflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"followglwinflag",(char *)&followglwinflag, 
                                                        TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"all3flag",(char *)&all3flag,TCL_LINK_BOOLEAN);
  /*=======================================================================*/
  /***** link global INT variables to tcl equivalents */
  Tcl_LinkVar(interp,"zf",(char *)&zf, TCL_LINK_INT);
  Tcl_LinkVar(interp,"xdim",(char *)&xdim, TCL_LINK_INT);
  Tcl_LinkVar(interp,"ydim",(char *)&ydim, TCL_LINK_INT);
  Tcl_LinkVar(interp,"imnr1",(char *)&imnr1, TCL_LINK_INT);
  Tcl_LinkVar(interp,"plane",(char *)&plane, TCL_LINK_INT);
  Tcl_LinkVar(interp,"imc",(char *)&imc, TCL_LINK_INT);
  Tcl_LinkVar(interp,"ic",(char *)&ic, TCL_LINK_INT);
  Tcl_LinkVar(interp,"jc",(char *)&jc, TCL_LINK_INT);
  Tcl_LinkVar(interp,"prad",(char *)&prad, TCL_LINK_INT);
  Tcl_LinkVar(interp,"pradlast",(char *)&pradlast, TCL_LINK_INT);
  Tcl_LinkVar(interp,"dip_spacing",(char *)&dip_spacing, TCL_LINK_INT);
  Tcl_LinkVar(interp,"editedimage",(char *)&editedimage, TCL_LINK_INT);
  Tcl_LinkVar(interp,"surflinewidth",(char *)&surflinewidth, TCL_LINK_INT);
  Tcl_LinkVar(interp,"white_lolim",(char *)&white_lolim, TCL_LINK_INT);
  Tcl_LinkVar(interp,"white_hilim",(char *)&white_hilim, TCL_LINK_INT);
  Tcl_LinkVar(interp,"gray_hilim",(char *)&gray_hilim, TCL_LINK_INT);
  Tcl_LinkVar(interp,"lim3",(char *)&lim3, TCL_LINK_INT);
  Tcl_LinkVar(interp,"lim2",(char *)&lim2, TCL_LINK_INT);
  Tcl_LinkVar(interp,"lim1",(char *)&lim1, TCL_LINK_INT);
  Tcl_LinkVar(interp,"lim0",(char *)&lim0, TCL_LINK_INT);
  Tcl_LinkVar(interp,"second_im_full",(char *)&second_im_full, TCL_LINK_INT);
  /*=======================================================================*/
  /***** link global DOUBLE variables to tcl equivalents (were float) */
  Tcl_LinkVar(interp,"fsquash",(char *)&fsquash, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fthresh",(char *)&fthresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fsthresh",(char *)&fsthresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fscale",(char *)&fscale, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ffrac3",(char *)&ffrac3, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ffrac2",(char *)&ffrac2, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ffrac1",(char *)&ffrac1, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ffrac0",(char *)&ffrac0, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"xtalairach",(char *)&xtalairach, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ytalairach",(char *)&ytalairach, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ztalairach",(char *)&ztalairach, TCL_LINK_DOUBLE);
  /*=======================================================================*/
  /***** link global malloced STRING vars */
  Tcl_LinkVar(interp,"home",        (char *)&subjectsdir,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"session",     (char *)&srname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subject",     (char *)&pname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"imtype",      (char *)&imtype,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"surface",     (char *)&surface,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"abs_imstem",  (char *)&mfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"sfname",      (char *)&sfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subjtmpdir",  (char *)&tfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"dip",         (char *)&dipfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"dec",         (char *)&decfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"hpts",        (char *)&hpfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"htrans",      (char *)&htfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"rgb",         (char *)&sgfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"fs",          (char *)&fsfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"fm",          (char *)&fmfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"curv",        (char *)&cfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"insurf",      (char *)&sfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"transform",   (char *)&xffname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"script",      (char *)&rfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"selectedpixval",(char *)&selectedpixval, TCL_LINK_INT);
  /*=======================================================================*/
  /***** twitzels stuff ****/
  Tcl_LinkVar(interp,"f2thresh",(char *)&f2thresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fslope",(char *)&fslope, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fmid",(char *)&fmid, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"do_overlay",(char *)&do_overlay,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"do_interpolate",(char *)&do_interpolate,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"overlay_frame",(char *)&overlay_frame, TCL_LINK_INT);
  Tcl_LinkVar(interp,"colscale",(char *)&colscale,TCL_LINK_INT);
  Tcl_LinkVar(interp,"floatstem",(char *)&overlay_file,TCL_LINK_STRING);

  strcpy(rfname,script_tcl);  /* save in global (malloc'ed in Program) */

  /* run tcl/tk startup script to set vars, make interface; no display yet */
  printf("tkmedit: interface: %s\n",tkmedit_tcl);
  code = Tcl_EvalFile(interp,tkmedit_tcl);
  if (*interp->result != 0)  printf(interp->result);

  /* if command line script exists, now run as batch job (possibly exiting) */
  if (scriptok) {    /* script may or may not open gl window */
    printf("tkmedit: run tcl script: %s\n",script_tcl);
    code = Tcl_EvalFile(interp,script_tcl);
    if (*interp->result != 0)  printf(interp->result);
  } else {
    ; /* tkmedit has already opened gl window if no script */
  }

  /* always start up command line shell too */
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  if (tty)
    Prompt(interp, 0);
  fflush(stdout);
  Tcl_DStringInit(&command);
  Tcl_ResetResult(interp);

  /*Tk_MainLoop();*/  /* standard */

  /* dual event loop (interface window made now) */
  while(tk_NumMainWindows > 0) {
    while (Tk_DoOneEvent(TK_ALL_EVENTS|TK_DONT_WAIT)) {
      /* do all the tk events; non-blocking */
    }
    do_one_gl_event(interp);
    /*sginap((long)1);*/   /* block for 10 msec */
    usecnap(10000);        /* block for 10 msec */
  }

  Tcl_Eval(interp, "exit");
  exit(0);
}

void usecnap(int usec)
{
  struct timeval delay;

  delay.tv_sec = 0;
  delay.tv_usec = (long)usec;
  select(0,NULL,NULL,NULL,&delay);
}

/*=== from TkMain.c ===================================================*/
static void StdinProc(clientData, mask)
  ClientData clientData;
  int mask;
{
#define BUFFER_SIZE 4000
  char input[BUFFER_SIZE+1];
  static int gotPartial = 0;
  char *cmd;
  int code, count;

  count = read(fileno(stdin), input, BUFFER_SIZE);
  if (count <= 0) {
    if (!gotPartial) {
      if (tty) {Tcl_Eval(interp, "exit"); exit(1);}
      else     {Tk_DeleteFileHandler(0);}
      return;
    }
    else count = 0;
  }
  cmd = Tcl_DStringAppend(&command, input, count);
  if (count != 0) {
    if ((input[count-1] != '\n') && (input[count-1] != ';')) {
      gotPartial = 1;
      goto prompt; }
    if (!Tcl_CommandComplete(cmd)) {
      gotPartial = 1;
      goto prompt; }
  }
  gotPartial = 0;
  Tk_CreateFileHandler(0, 0, StdinProc, (ClientData) 0);
  code = Tcl_RecordAndEval(interp, cmd, TCL_EVAL_GLOBAL);
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  Tcl_DStringFree(&command);
  if (*interp->result != 0)
    if ((code != TCL_OK) || (tty))
      puts(interp->result);
  prompt:
  if (tty)  Prompt(interp, gotPartial);
  Tcl_ResetResult(interp);
}

/*=== from TkMain.c ===================================================*/
static void Prompt(interp, partial)
  Tcl_Interp *interp;
  int partial;
{
  char *promptCmd;
  int code;

  promptCmd = Tcl_GetVar(interp,
      partial ? "tcl_prompt2" : "tcl_prompt1", TCL_GLOBAL_ONLY);
  if (promptCmd == NULL) {
    defaultPrompt:
    if (!partial)
      fputs("% ", stdout);
  }
  else {
    code = Tcl_Eval(interp, promptCmd);
    if (code != TCL_OK) {
      Tcl_AddErrorInfo(interp,
                    "\n    (script that generates prompt)");
      fprintf(stderr, "%s\n", interp->result);
      goto defaultPrompt;
    }
  }
  fflush(stdout);
}

