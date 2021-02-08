/**
 * @brief Tcl/Tk-based MRI volume registration utility
 *
 * See: http://surfer.nmr.mgh.harvard.edu/fswiki/FsTutorial/Talairach
 */
/*
 * Original Authors: Martin Sereno and Anders Dale, 1996; Doug Greve, 2002
 *
 * Copyright (C) 2002-2011, CorTechs Labs, Inc. (La Jolla, CA) and
 * The General Hospital Corporation (Boston, MA).
 *
 * Terms and conditions for use, reproduction, distribution and contribution
 * are found in the 'FreeSurfer/CorTechs Software License Agreement' contained
 * in the file 'license.cortechs.txt' found in the FreeSurfer distribution,
 * and here:
 *
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferCorTechsLicense
 *
 * Reporting: freesurfer@nmr.mgh.harvard.edu
 *
 */

// if NO_GUI is defined when built, then tkregister2 is built without any of
// the GUI elements, leaving only the command-line processing capability,
// such as that used by recon-all.  this build functionality is useful when
// when building on systems which don't readily have Tcl/Tk and OpenGL.
#ifndef NO_GUI
#define HAVE_TCL_TK_GL 1
#endif

#ifndef lint
#endif /* lint */

#ifdef HAVE_TCL_TK_GL
#define TCL
#include <tcl.h>
#include <tk.h>
#endif // HAVE_TCL_TK_GL
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>

#ifdef HAVE_TCL_TK_GL

#ifndef OPENGL
#define OPENGL
#endif

#ifndef RGB
#define RGB
#endif

#ifndef TCL
#define TCL
#endif

#ifdef OPENGL
#  include "xwindow.h"
#  include "macros.h"
#  include "xwindow.h"
#  include <X11/keysym.h>
#ifdef RGB
#  define swapbuffers()    pseudo_swapbuffers(); tkoSwapBuffers()
#else
#  define swapbuffers()    tkoSwapBuffers()
#endif
#  define frontbuffer(X)   glDrawBuffer(GL_FRONT)
#  define backbuffer(X)    glDrawBuffer(GL_BACK)
#  define clear()          glClear(GL_COLOR_BUFFER_BIT)
//#  define getorigin(X,Y)   *(X) = w.x; *(Y) = 1024 - w.y - w.h
#  define getorigin(X,Y)   *(X) = w.x; *(Y) = 1024 - w.y - w.h
/*WINDOW_REC w;*/
#  define getsize(X,Y)     *(X) = w.w; *(Y) = w.h
#  define Colorindex       unsigned short
#  define color(X)         glIndexs(X)
#  define mapcolor(I,R,G,B)  \
            tkoSetOneColor((int)I,(float)R/255.0,(float)G/255.0,(float)B/255.0)
#ifdef RGB
#  define rectwrite(X0,Y0,X1,Y1,P) \
            glRasterPos2i(X0,Y0); \
            glDrawPixels(X1+1,Y1+1,GL_RGB,GL_UNSIGNED_BYTE,P)
#else
#  define rectwrite(X0,Y0,X1,Y1,P) \
            glRasterPos2i(X0,Y0); \
            glDrawPixels(X1+1,Y1+1,GL_COLOR_INDEX,GL_UNSIGNED_SHORT,P)
#endif
#else
#  include <gl.h>
#  include <device.h>
#endif
#endif // HAVE_TCL_TK_GL

#include "proto.h"
#include "error.h"
#include "diag.h"
#include "utils.h"
#include "const.h"
#include "machine.h"
#include "MRIio_old.h"
#include "mri.h"
#include "mri2.h"
#include "mrisurf.h"
#include "mri_identify.h"
#include "registerio.h"
#include "version.h"
#include "fio.h"
#include "pdf.h"
#include "resample.h"
#include "pdf.h"
#include "fmriutils.h"
#include "mri_conform.h"

/* Prototypes */

void open_window(char name[]);
void resize_window_intstep();
void move_window(int x, int y);
void reload_buffers();
void blinkbuffers();
void record_swapbuffers();
void resize_buffers(int x, int y);
void read_reg(char fname[]);
void read_fslreg(char *fname);
void write_reg(char fname[]);
void write_fslreg(char *fname);
void write_xfmreg(char *fname);
void write_lta(char *fname);
void write_freeviewreg(char *fname);
void make_backup(char fname[]);
void save_rgb(char fname[]);
void scrsave_to_rgb(char fname[]);
void pix_to_rgb(char *fname);
void downslice();
void upslice();
void goto_point(char dir[]);
void write_point(char dir[]);
void rotate_brain(float a, char c);
void align_points();
void translate_brain(float a,char c);
void scale_brain(float s, char c);
void mirror_brain();
void set_cursor(float xpt, float ypt, float zpt);
void set_scale();
void redraw();
void pop_gl_window();
void mri2pix(float xpt, float ypt, float zpt, int *jpt, int *ipt, int *impt);
int imval(float px,float py,float pz);
float Error(int p,float dp);
void optimize(int maxiter);
void optimize2();
void read_images(char fpref[]);
void read_second_images(char fpref[]);
void select_pixel(short sx,short sy);
void transform(float x1,float y1,float z1,
               float *x2,float *y2,float *z2,
               float M[4][4]);
void draw_image(int imc,int ic,int jc);
void draw_image2(int imc,int ic,int jc);
void blur(float factor);
void make_filenames(char *lsubjectsdir);
void read_float_images(float ***fim,char *format,
                       int nslices,int nperslice,
                       int xdim,int ydim,short **buf);
void usecnap(int usec);
void initcolormap();
void pseudo_swapbuffers();
int crScreen2AnatInd(int c, int r, int *cor, int *hor, int *sag);
int AnatInd2AnatXYZ(int cor, int hor, int sag, float *x, float *y, float *z);
int crScreen2AnatXYZ(int c, int r, float *x, float *y, float *z);
int FuncXYZ2FuncInd(float x, float y, float z,
                    float *col, float *row, float *slice);
int draw_cross_hair(int rScreen, int cScreen);
int erase_cross_hair(int rScreen, int cScreen);

void UpdateMatrices(void);
MATRIX *ScreenCR2XYZMtx(MATRIX *T);

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int isflag(const char *flag);
static int singledash(char *flag);
static int stringmatch(const char *str1, const char *str2);
static int nth_is_arg(int nargc, char **argv, int nth);
static int checkfmt(char *fmt);
#ifdef HAVE_TCL_TK_GL
static int MRItagVol(MRI *mri, float val);
#endif // HAVE_TCL_TK_GL
static int MRIisConformant(MRI *vol);

MATRIX *Load4x4(char *fname);
char *Vox2VoxFName = NULL;
MATRIX *Vox2Vox = NULL;

#ifndef TRUE
#  define TRUE 1
#endif
#ifndef FALSE
#  define FALSE 0
#endif

#define APD1    1   /* imtypes */
#define BSHORT  2
#define SKIP    3
#define AFNI    4

#ifndef SQR
#define SQR(x)       ((x)*(x))
#endif
#define MATCH(A,B)   (!strcmp(A,B))
#define MATCH_STR(S) (!strcmp(str,S))

#define HEADERSIZE 80
#define NUMVALS 256
#define MAXIM 256
#define MAXPTS 10000
#define MAXPARS 10
#ifdef RGB
#define MAPOFFSET 0
#else
#define MAPOFFSET 1000
#endif
#define CORONAL 0
#define HORIZONTAL 1
#define SAGITTAL 2
#define NAME_LENGTH STRLEN
#define BLINK_DELAY 30
#define BLINK_TIME 20
#define MOTIF_XFUDGE   8
#define MOTIF_YFUDGE  32
#define TMP_DIR "tmp"
#define SLICE1_POS_UNDEF   9999.0
/* overlay modes */
#define TARGET   1
#define MOVEABLE 2

int WINDOW_ROWS = 512;
int WINDOW_COLS = 512;

float FOV = 256;
float mov_scale = 1;
int plane = CORONAL;
int plane_init = CORONAL;
int slice_init = -1;
int xnum=256,ynum=256;
int ptype;
float ps,st,xx0,xx1,yy0,yy1,zz0,zz1;
int zf,ozf;
float fsf;
int xdim,ydim;
unsigned long bufsize, bufsize_2;
unsigned char *buf, *buf_2;
unsigned char **im[MAXIM];
unsigned char **sim[6];
int changed[MAXIM];
#ifdef HAVE_TCL_TK_GL
GLubyte *vidbuf;
GLubyte *blinkbuft;
GLubyte *blinkbufm;
#endif // HAVE_TCL_TK_GL
unsigned char *binbuff;
int imnr0,imnr1,numimg;
int wx0=600,wy0=100;  /* (114,302) (117,90) */
int ptsflag = FALSE;
int maxflag = FALSE;
int updateflag = FALSE;
int blinkflag = FALSE;
int blinkdelay = BLINK_DELAY;
int blinktime = BLINK_TIME;
int overlay_mode = MOVEABLE;
int overlay_mode_init = MOVEABLE;
int visible_mode = 0;
int last_visible_mode = 0;
int visible_plane = 0;
int last_visible_plane = 0;
int editedmatrix = FALSE;
int maskflag = FALSE;
int scrsaveflag = TRUE;
int openglwindowflag = FALSE;
int promptflag = FALSE;
int followglwinflag = TRUE;
int initpositiondoneflag = FALSE;
int npts = 0;
int prad = 0;
float TM[4][4];
float tm[4][4];
MATRIX *RegMat=NULL, *XFM=NULL;
double ps_2,st_2,fscale_2=0.0; /* was float */
int float2int = 0;
int float2int_use = FLT2INT_ROUND;
float xx0_2,xx1_2,yy0_2,yy1_2,zz0_2,zz1_2;
int xnum_2,ynum_2,numimg_2;
int imnr0_2,imnr1_2,xdim_2=0,ydim_2=0;
float **fim_2[MAXIM];
float ptx[MAXPTS],pty[MAXPTS],ptz[MAXPTS];
float par[MAXPARS],dpar[MAXPARS];
int nslices=0,nperslice=0;

const char *Progname ;
double fthresh = 0.35;
double fsquash = 12.0;
double fscale = 255;
int imc=0,ic=0,jc=0;
float xc=0,yc=0,zc=0;
float xc_old=0,yc_old=0,zc_old=0;
int impt = -1,ipt = -1,jpt = -1;
int cScreenCur = 0, rScreenCur = 0;
int PixelSelected = 1;

char *freesurferhome; // FREESURFER_HOME
char *subjectsdir;   /* SUBJECTS_DIR */
char *srname;        /* sessiondir (funct: image) */
char *psrname;       /* parent sessiondir (funct: 970703MS) */
char *pname;         /* subject */
char *regfname;      /* register.dat */
char *xfmfname=NULL; /* something.xfm (minc reg mat) */
char *afname;        /* analyse.dat */
char *targpref;      /* abs single image structural stem name */
char *movformat;     /* abs single image epi structural stem name */
char *tfname;        /* (dir!) SUBJECTS_DIR/name/tmp/ */
char *sgfname;       /* (dir!) set: get from cwd: $session/rgb/ */
char *tkrtitle=NULL; /* window title */

char *int_regfname=NULL;  /* intermediate registration file */
char *int_vol_id=NULL;  /* intermediate volume */
MRI  *int_vol=NULL;      // MRI for intermediate volume
MATRIX *IntRegMat=NULL,*Int2MovRegMat=NULL;

int blinktop = 0;    /* says whats on top while blinking */
int invalid_buffers = 1;

int pswapnext = TARGET;
char colormap[512];

int use_draw_image2 = 1;
int interpmethod = SAMPLE_TRILINEAR;
int use_inorm = 1;
int use_colornorm = 0;
int DoSlicePrescription = 0;

char subjectid[1000];
int subjectidOverride = 0;
char *mov_vol_id = NULL;
int   mov_vol_fmt = MRI_VOLUME_TYPE_UNKNOWN;
const char *targ_vol_id;
int   targ_vol_fmt = MRI_VOLUME_TYPE_UNKNOWN;
char targ_vol_path[1000];
int  fstarg = 0;
int mkheaderreg = 0;
int mkheaderregCenter= 0;
#ifdef  HAVE_TCL_TK_GL
int noedit = 0;  // false by default, if gui
#endif // HAVE_TCL_TK_GL
#ifndef HAVE_TCL_TK_GL
int noedit = 1; // true by default, if without gui
#endif // HAVE_TCL_TK_GL

int fixtkreg = 1, fixonly = 0;
int identityreg = 0;
int LoadVol = 1;
int tagmov = 0;

MRI *mov_vol, *targ_vol,*targ_vol0, *mritmp, *mrisurf;
MRI_SURFACE *surf;

//float movimg[WINDOW_ROWS][WINDOW_COLS];
//float targimg[WINDOW_ROWS][WINDOW_COLS];
//int surfimg[WINDOW_ROWS][WINDOW_COLS];
float **movimg, **targimg;
int **surfimg;

int debug;

MATRIX *Ttarg, *Tmov, *invTtarg, *invTmov;
MATRIX *Tscreen, *Qtarg, *Qmov;

char *fslregfname;
char *fslregoutfname;
char *freeviewfname;
char *xfmoutfname=NULL;
MATRIX *invDmov, *FSLRegMat, *invFSLRegMat;
MATRIX *Mtc, *invMtc, *Vt2s,*Ttargcor, *invTtargcor;
MATRIX *Dtargcor, *invDtargcor, *Dtarg, *invDtarg;
MATRIX *vox2ras_targ=NULL,*ras2vox_targ=NULL;

int LoadSurf = 0, UseSurf=0;
const char *surfname = "white";
char surf_path[2000];
int lhsurf_only = 0,rhsurf_only = 0;
int fstal=0, fixxfm=1;
char talxfmfile[2000],talxfmdir[2000];
const char *talxfmname = "talairach.xfm";
char tmpstr[2000];

char *mov_ostr = NULL; // orientation string for mov
char *targ_ostr = NULL; // orientation string for targ
int mov_frame = 0;
char *pc;

float int_ipr, int_bpr, int_fscale;
int int_float2int,err;
char *xfmfileinfo=NULL;

char *ltafname = NULL;
TRANSFORM *FSXform = NULL;
LTA *lta = NULL;
LT  *linxfm = NULL;
char *ltaoutfname=NULL;

int checkreg = 0;

#ifndef HAVE_TCL_TK_GL
#define ClientData void*
#define Tcl_Interp void*
#endif // HAVE_TCL_TK_GL

int ZeroCRAS = 0;
MATRIX *Ctarg, *invCtarg, *Starg, *Mcras0, *invMcras0;

int DoFMovTarg = 0;

char *DetFile = NULL;
int AllocBuffs(void);
static int istringnmatch(const char *str1, const char *str2, int n);

char *seg_vol_id = NULL;
MRI *seg_vol = NULL;
COLOR_TABLE *ctab = NULL;
char *ctabfile = NULL;
int DoASeg=0, DoAParcASeg=0, DoWMParc=0;
int ShowSeg=0;
char *FREESURFER_HOME=NULL;
char *tkregister_tcl = NULL;
const char *fstaltarg = "mni305.cor.mgz";
int SurfRGB[3] = {0,255,0};
int invLTAOut=0;
double angles[3] = {0,0,0};
MATRIX *Mrot = NULL;
double xyztrans[3] = {0,0,0};
MATRIX *Mtrans = NULL;
int conformTarget = 0;

/**** ------------------ main() ------------------------------- ****/
int Register(ClientData clientData,
             Tcl_Interp *interp,
             int argc, char *argv[]) {
  int i,j,err;
  FILE *fp;
#ifdef HAVE_TCL_TK_GL
  int n,c,r,s;
  MATRIX *Vxyz=NULL, *Vcrs=NULL;
#endif // HAVE_TCL_TK_GL

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if (argc == 0) usage_exit();

  memmove(subjectid,"subject-unknown",strlen("subject-unknown"));

  subjectsdir = getenv("SUBJECTS_DIR");
  if (subjectsdir==NULL) {
    printf("ERROR: SUBJECTS_DIR undefined. Use setenv or -sd\n");
    exit(1);
  }
  freesurferhome = getenv("FREESURFER_HOME");
  if (freesurferhome==NULL) {
    printf("ERROR: FREESURFER_HOME undefined. \n");
    exit(1);
  }

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);
  printf("%s\n",getVersion().c_str());
  printf("Diagnostic Level %d\n",Gdiag_no);

  AllocBuffs();


  /* read the registration here to get subjectid */
  if (!mkheaderreg && fslregfname == NULL && !fstal &&
      int_vol_id == NULL && !identityreg && ltafname == NULL)
    read_reg(regfname);
  // Just use identity
  if (identityreg) RegMat = MatrixIdentity(4,NULL);

  if(Mrot){
    printf("Applying rotation matrix (R=M*R)\n");
    printf("Current Reg Matrix is:\n");
    MatrixPrint(stdout,RegMat);
    printf("  Angles (deg): %lf %lf %lf\n",angles[0]*180/M_PI,angles[1]*180/M_PI,angles[2]*180/M_PI);
    printf("  Angles (rad): %lf %lf %lf\n",angles[0],angles[1],angles[2]);
    printf("  Rotation matrix:\n");
    MatrixPrint(stdout,Mrot);
    RegMat = MatrixMultiply(Mrot,RegMat,RegMat);
  }
  if(Mtrans){
    printf("Applying translation matrix (R=M*R)\n");
    printf("Current Reg Matrix is:\n");
    MatrixPrint(stdout,RegMat);
    printf("  Trans (mm): %lf %lf %lf\n",xyztrans[0],xyztrans[1],xyztrans[2]);
    printf("  Translation matrix:\n");
    MatrixPrint(stdout,Mtrans);
    RegMat = MatrixMultiply(Mtrans,RegMat,RegMat);
  }

  if(DoASeg){
    sprintf(tmpstr,"%s/%s/mri/aseg.mgz",subjectsdir,subjectid);
    seg_vol_id = strcpyalloc(tmpstr);
  }
  if(DoAParcASeg){
    sprintf(tmpstr,"%s/%s/mri/aparc+aseg.mgz",subjectsdir,subjectid);
    seg_vol_id = strcpyalloc(tmpstr);
  }
  if(DoWMParc){
    sprintf(tmpstr,"%s/%s/mri/wmparc.mgz",subjectsdir,subjectid);
    seg_vol_id = strcpyalloc(tmpstr);
  }
  if(seg_vol_id){
    FREESURFER_HOME = getenv("FREESURFER_HOME");
    if(ctabfile == NULL){
      ctabfile = (char *) calloc(sizeof(char),1000);
      sprintf(ctabfile,"%s/FreeSurferColorLUT.txt",FREESURFER_HOME);
      ctab = CTABreadASCII(ctabfile);
    }
    seg_vol = MRIread(seg_vol_id);
    if(seg_vol == NULL) exit(1);
    printf(" \n");
    printf(" \n");
    printf(" \n");
    printf("Computing seg boundary  ... ");
    mritmp = MRIsegBoundary(seg_vol);
    printf(" done\n");
    printf(" \n");
    printf(" \n");
    printf(" \n");
    MRIfree(&seg_vol);
    seg_vol = mritmp;
    fscale_2 = 1.0;
  }

  if(fstal) {
    if (subjectsdir == NULL) subjectsdir = getenv("SUBJECTS_DIR");
    if (subjectsdir==NULL) {
      printf("ERROR: SUBJECTS_DIR undefined. Use setenv or --sd\n");
      exit(1);
    }
    // Load the talairach.xfm
    int req = snprintf(talxfmdir, 2000, "%s/%s/mri/transforms",subjectsdir,subjectid);
    if (req >= 2000) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    req = snprintf(talxfmfile, 2000, "%s/%s",talxfmdir,talxfmname);
    if (req >= 2000) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (!mkheaderreg) {
      if (!fio_DirIsWritable(talxfmdir,0)) {
        printf("\n");
        printf("\n");
        printf("WARNING: cannot write to %s.\n",talxfmdir);
        printf("You will not be able to save any edits.\n");
        printf("Hit Enter to continue: ");
        getc(stdin);
        printf(" ... continuing\n");
        printf("\n");
      }
      if (regio_read_mincxfm(talxfmfile, &RegMat, &xfmfileinfo)) exit(1);
      printf("%s ---------------------\n",talxfmname);
      MatrixPrint(stdout,RegMat);
      for (i=0;i<4;i++) {
        for (j=0;j<4;j++) {
          tm[i][j] = RegMat->rptr[i+1][j+1];
        }
      }
    }
    mov_vol_id = (char *) calloc(sizeof(char),2000);
    sprintf(mov_vol_id,"%s/average/%s",freesurferhome,fstaltarg);

    ps_2 = 1.0;
    st_2 = 1.0;
  }

  // Get full path to target
  if (fstarg) {
    // Relative path is specified, compute full
    if (subjectsdir == NULL) subjectsdir = getenv("SUBJECTS_DIR");
    if (subjectsdir==NULL) {
      printf("ERROR: SUBJECTS_DIR undefined. Use setenv or -sd\n");
      exit(1);
    }
    // First try .mgz
    int req = snprintf(targ_vol_path,1000,"%s/%s/mri/%s.mgz",
                   subjectsdir,subjectid,targ_vol_id);
    if (req >= 1000) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    if (! fio_FileExistsReadable(targ_vol_path)) {
      // Now try COR
      req = snprintf(targ_vol_path, 1000, "%s/%s/mri/%s",subjectsdir,subjectid,targ_vol_id);
      if (req >= 1000) {
        std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
      }
      if (! fio_FileExistsReadable(targ_vol_path)) {
        printf("ERROR: could not find %s as either mgz or COR\n",targ_vol_id);
        printf("%s\n",targ_vol_path);
        exit(1);
      }
    }
  } else {
    // Full path is specified
    memmove(targ_vol_path,targ_vol_id,strlen(targ_vol_id));
    if (! fio_FileExistsReadable(targ_vol_path)) {
      printf("ERROR: could not find %s\n",targ_vol_path);
      exit(1);
    }
  }

  /* extract and check types */
  if (targ_vol_fmt == MRI_VOLUME_TYPE_UNKNOWN) {
    targ_vol_fmt = mri_identify(targ_vol_path);
    if (targ_vol_fmt == MRI_VOLUME_TYPE_UNKNOWN) {
      printf("ERROR: cannot determine type of %s\n", targ_vol_path);
      exit(1);
    }
  }
  if (mov_vol_fmt == MRI_VOLUME_TYPE_UNKNOWN) {
    mov_vol_fmt = mri_identify(mov_vol_id);
    if (mov_vol_fmt == MRI_VOLUME_TYPE_UNKNOWN) {
      printf("ERROR: cannot determine type of %s\n", mov_vol_id);
      exit(1);
    }
  }

  /* Check that both volumes are accessible */
  mritmp = MRIreadHeader(targ_vol_path,targ_vol_fmt);
  if (mritmp == NULL) {
    printf("ERROR: could not read %s as %d\n",targ_vol_path,targ_vol_fmt);
    exit(1);
  }
  MRIfree(&mritmp);
  mritmp = MRIreadHeader(mov_vol_id,mov_vol_fmt);
  if (mritmp == NULL) {
    printf("ERROR: could not read %s as %d\n",mov_vol_id,mov_vol_fmt);
    exit(1);
  }
  MRIfree(&mritmp);

  /*------------------------------------------------------*/
  printf("INFO: loading target %s\n",targ_vol_path);
  if(LoadVol)  {
    targ_vol = MRIreadType(targ_vol_path,targ_vol_fmt);
    if(targ_vol == NULL) exit(1);
    if(targ_vol->nframes > 1){
      // extract the first frame if necessary
      mritmp = fMRIframe(targ_vol, 0, NULL);
      MRIfree(&targ_vol);
      targ_vol = mritmp;
    }
    if(conformTarget){
      // only do this if the registration was computed to a conformed target
      // this happens when running mri_em_register to a GCA. 
      printf("Conforming target\n");
      mritmp = MRIconform(targ_vol);
      MRIfree(&targ_vol);
      targ_vol = mritmp;
    }
  }
  else {
    targ_vol = MRIreadHeader(targ_vol_path,targ_vol_fmt);
    if(targ_vol == NULL) exit(1);
  }
  if(fstal && ZeroCRAS){
    printf("Zeroing CRAS of target\n");
    Starg = MRIxfmCRS2XYZ(targ_vol,0);
    targ_vol->c_r = 0;
    targ_vol->c_a = 0;
    targ_vol->c_s = 0;
    Ctarg = MRIxfmCRS2XYZ(targ_vol,0);
    invCtarg = MatrixInverse(Ctarg,NULL);
    Mcras0    = MatrixMultiply(Starg,invCtarg,NULL);
    invMcras0 = MatrixInverse(Mcras0,NULL);
    // At this point, RegMat holds tal.xfm
    RegMat = MatrixMultiply(RegMat,Mcras0,RegMat);
    printf("new xfm -----------------\n");
    MatrixPrint(stdout,RegMat);
    printf("---------------------\n");
    MatrixFree(&Ctarg);
    MatrixFree(&invCtarg);
    MatrixFree(&Starg);
  }
  if(targ_ostr) {
    printf("Setting targ orientation to %s\n",targ_ostr);
    MRIorientationStringToDircos(targ_vol, targ_ostr);
  }
  if (targ_vol->type != MRI_FLOAT && LoadVol) {
    printf("INFO: changing target type to float\n");
    mritmp = MRIchangeType(targ_vol,MRI_FLOAT,0,0,0);
    if (mritmp == NULL) {
      printf("ERROR: could change type\n");
      exit(1);
    }
    MRIfree(&targ_vol);
    targ_vol = mritmp;
  }

  /* Make target conformant if necessary */
  if (! MRIisConformant(targ_vol) ) {
    printf("INFO: target does not conform to COR format, so I'm going to\n");
    printf("reslice to COR. This will not affect the final registration.\n");
    mritmp = MRIallocSequence(256,256,256, MRI_FLOAT, 1) ;
    mritmp->xsize = 1;
    mritmp->ysize = 1;
    mritmp->zsize = 1;
    mritmp->x_r = -1;
    mritmp->x_a = 0;
    mritmp->x_s = 0;
    mritmp->y_r = 0;
    mritmp->y_a = 0;
    mritmp->y_s = -1;
    mritmp->z_r = 0;
    mritmp->z_a = 1;
    mritmp->z_s = 0;
    mritmp->c_r = targ_vol->c_r;
    mritmp->c_a = targ_vol->c_a;
    mritmp->c_s = targ_vol->c_s;

    Ttarg     = MRIxfmCRS2XYZtkreg(targ_vol);
    invTtarg  = MatrixInverse(Ttarg,NULL);
    Ttargcor  = MRIxfmCRS2XYZtkreg(mritmp);
    Dtargcor  = MRIxfmCRS2XYZ(mritmp,0);
    Dtarg     = MRIxfmCRS2XYZ(targ_vol,0);
    invDtarg  = MatrixInverse(Dtarg,0);
    Vt2s = MatrixMultiply(invDtarg,Dtargcor,NULL);

    invTtargcor  = MatrixInverse(Ttargcor,NULL);
    invDtargcor  = MatrixInverse(Dtargcor,NULL);
    Mtc = MatrixMultiply(Ttargcor,invDtargcor,NULL);
    MatrixMultiply(Mtc,Dtarg,Mtc);
    MatrixMultiply(Mtc,invTtarg,Mtc);

    if(LoadVol) {
      err = MRIvol2Vol(targ_vol, mritmp, Vt2s, SAMPLE_TRILINEAR, 0);
      if(err) exit(1);
      MRIfree(&targ_vol);
    }
    targ_vol = mritmp;

    // Keep a copy of the uncorformed header
    targ_vol0 = MRIreadHeader(targ_vol_path,targ_vol_fmt);

    MatrixFree(&Ttargcor);
    MatrixFree(&invTtargcor);
    MatrixFree(&Dtargcor);
    MatrixFree(&Dtarg);
    MatrixFree(&invDtarg);
    MatrixFree(&Ttarg);
    MatrixFree(&invTtarg);
    MatrixFree(&Vt2s);
    //MRIwrite(targ_vol,"cor.mgh");
  } else {
    Mtc = MatrixIdentity(4,NULL);
    targ_vol0 = targ_vol;
  }
  invMtc = MatrixInverse(Mtc,NULL);

  if(!fstal || !fixxfm) Ttarg = MRIxfmCRS2XYZtkreg(targ_vol);
  else                  Ttarg = MRIxfmCRS2XYZ(targ_vol,0);
  invTtarg = MatrixInverse(Ttarg,NULL);
  printf("Ttarg: --------------------\n");
  MatrixPrint(stdout,Ttarg);

  /*------------------------------------------------------*/
  printf("INFO: loading movable %s\n",mov_vol_id);
  if (LoadVol)  mov_vol = MRIreadType(mov_vol_id,mov_vol_fmt);
  else         mov_vol = MRIreadHeader(mov_vol_id,mov_vol_fmt);
  if (mov_vol == NULL) {
    printf("ERROR: could not read %s\n",mov_vol_id);
    exit(1);
  }
  if (mov_ostr) {
    printf("Setting mov orientation to %s\n",mov_ostr);
    MRIorientationStringToDircos(mov_vol, mov_ostr);
  }
  if (mov_vol->type != MRI_FLOAT  && LoadVol) {
    printf("INFO: changing move type to float\n");
    mritmp = MRISeqchangeType(mov_vol,MRI_FLOAT,0,0,0);
    if (mritmp == NULL) {
      printf("ERROR: could change type\n");
      exit(1);
    }
    MRIfree(&mov_vol);
    mov_vol = mritmp;
  }
  if(!fstal || !fixxfm) Tmov = MRIxfmCRS2XYZtkreg(mov_vol);
  else                  Tmov = MRIxfmCRS2XYZ(mov_vol,0);
  invTmov = MatrixInverse(Tmov,NULL);
  printf("Tmov: --------------------\n");
  //Tmov = MRIxfmCRS2XYZtkreg(mov_vol); // should this be here?
  MatrixPrint(stdout,Tmov);

  /*------------------------------------------------------*/
  printf("mkheaderreg = %d, float2int = %d\n",mkheaderreg,float2int);
  if(mkheaderreg) {
    /* Compute Reg from Header Info */
    printf("Computing reg from header (and possibly input matrix)\n");
    if(!Vox2Vox){
      RegMat = MRItkRegMtx(targ_vol,mov_vol,XFM);
      if(mkheaderregCenter){
	RegMat->rptr[1][4] = 0;
	RegMat->rptr[2][4] = 0;
	RegMat->rptr[3][4] = 0;
      }
    }
    else         RegMat = MRItkRegMtxFromVox2Vox(targ_vol,mov_vol,Vox2Vox);
  } else if (fslregfname != NULL) {
    /* Compute Reg from FSLReg */
    RegMat = MRIfsl2TkReg(targ_vol0, mov_vol, FSLRegMat);
    MatrixMultiply(RegMat,invMtc,RegMat);
  } else if (int_vol_id != NULL) {
    printf("Creating reg from intermediate\n");
    int_vol = MRIreadHeader(int_vol_id,MRI_VOLUME_TYPE_UNKNOWN);
    if (int_vol == NULL) {
      printf("ERROR: reading intermediate volume %s\n",int_vol_id);
      exit(1);
    }
    err = regio_read_register(int_regfname, &pc, &int_ipr, &int_bpr,
                              &int_fscale, &IntRegMat, &int_float2int);
    if (err) {
      printf("ERROR: reading intermediate registration %s\n",int_regfname);
      exit(1);
    }
    Int2MovRegMat = MRItkRegMtx(int_vol,mov_vol,NULL);
    RegMat = MatrixMultiply(Int2MovRegMat,IntRegMat,NULL);
    if(fscale_2 == 0.0) fscale_2 = int_fscale;
    MRIfree(&int_vol);
    MatrixFree(&IntRegMat);
    MatrixFree(&Int2MovRegMat);
  } else {
    /* Use the registration file that was loaded above */
    MatrixMultiply(RegMat,invMtc,RegMat);
    if (float2int == FLT2INT_TKREG) {
      if (fixtkreg) {
        printf("INFO: making tkreg matrix compatible with round\n");
        RegMat = MRIfixTkReg(mov_vol,RegMat);
        printf("---- New registration matrix --------\n");
        MatrixPrint(stdout,RegMat);
        printf("---------------------------------------\n");
        for (i=0;i<4;i++) {
          for (j=0;j<4;j++) {
            tm[i][j] = RegMat->rptr[i+1][j+1];
          }
        }
      } else {
        printf("INFO: this registration matrix was created using\n");
        printf("      the old float2int truncation method and you\n");
        printf("      have chosen not to fix it (--nofix), so\n");
        printf("      the registration may be shifted.\n");
      }
    }
  }


  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      tm[i][j] = RegMat->rptr[i+1][j+1];
    }
  }
  if (mkheaderreg || fslregfname || identityreg || ltafname) {
    printf("---- Input registration matrix (computed) --------\n");
    MatrixPrint(stdout,RegMat);
    printf("---------------------------------------\n");
    ps_2 = mov_vol->xsize;
    st_2 = mov_vol->zsize;
    if (fscale_2 == 0) fscale_2 = .15;
    float2int = FLT2INT_ROUND;
  }

  printf("---- Input registration matrix --------\n");
  MatrixPrint(stdout,RegMat);
  printf("Determinant %g\n",MatrixDeterminant(RegMat));
  if(DetFile){
    fp = fopen(DetFile,"w");
    fprintf(fp,"%lf\n",MatrixDeterminant(RegMat));
    fclose(fp);
  }

  pname = subjectid;
  printf("subject = %s\n",subjectid);
  if (noedit) {
    write_reg(regfname);
    exit(0);
  }

#ifndef HAVE_TCL_TK_GL
  printf("\nERROR: This tkregister2 was built without the GUI "
         "(Tcl/Tk/OpenGL excluded).\n");
  printf("       Functionality is limited to the following flags:\n"
         "         --targ --mov --reg --fslregout --regheader ...\n");
  exit(1);
#else

  /*------------------------------------------------------*/
  if (LoadSurf) {
    /* alloc a volume data struct for the surface */
    mrisurf = MRIallocSequence(256,256,256, MRI_UCHAR, 1) ;
    if (mrisurf == NULL) {
      printf("ERROR: could not alloc volume for surf\n");
      exit(1);
    }
    mrisurf->x_r = -1.0;
    mrisurf->x_a =  0.0;
    mrisurf->x_s =  0.0;
    mrisurf->y_r =  0.0;
    mrisurf->y_a = +1.0;
    mrisurf->y_s =  0.0;
    mrisurf->z_r =  0.0;
    mrisurf->z_a =  0.0;
    mrisurf->z_s = -1.0;
    mrisurf->c_r = +0.5;
    mrisurf->c_a = -0.5;
    mrisurf->c_s = +0.5;

    // For the surf, vox2ras must be tkreg-style
    vox2ras_targ = MRIxfmCRS2XYZtkreg(targ_vol);
    ras2vox_targ = MatrixInverse(vox2ras_targ,NULL);

    if (! rhsurf_only) {
      /* load in the left hemi */
      sprintf(surf_path,"%s/%s/surf/lh.%s",subjectsdir,subjectid,surfname);
      printf("Reading lh surface %s\n",surfname);
      surf = MRISread(surf_path) ;
      if (!surf) {
        ErrorExit(ERROR_NOFILE,"%s: could not read surface %s",
                  Progname, surfname) ;
        exit(1);
      }
      printf("Done reading surface\n");

      /* make a surface mask in the volume */
      Vxyz = MatrixAlloc(4,1,MATRIX_REAL);
      Vxyz->rptr[4][1] = 1.0;
      for (n=0; n < surf->nvertices; n++) {
        Vxyz->rptr[1][1] = surf->vertices[n].x;
        Vxyz->rptr[2][1] = surf->vertices[n].y;
        Vxyz->rptr[3][1] = surf->vertices[n].z;
        Vcrs = MatrixMultiply(ras2vox_targ,Vxyz,Vcrs);
        c = nint(Vcrs->rptr[1][1]);
        r = nint(Vcrs->rptr[2][1]);
        s = nint(Vcrs->rptr[3][1]);
        if (c < 0 || c > 255 || r < 0 || r > 255 || s < 0 || s > 255) {
          printf("ERROR: vertex %d (%g,%g,%g) is out of range (%d,%d,%d)\n",
                 n,surf->vertices[n].x,surf->vertices[n].y,surf->vertices[n].z,
                 c,r,s);
          exit(1);
        }
        MRIvox(mrisurf,c,r,s) = 1;
      }
      MRISfree(&surf);
    }

    if (! lhsurf_only) {
      /* load in the right hemi */
      sprintf(surf_path,"%s/%s/surf/rh.%s",subjectsdir,subjectid,surfname);
      printf("Reading rh surface %s\n",surfname);
      surf = MRISread(surf_path) ;
      if (!surf)
        ErrorExit(ERROR_NOFILE,"%s: could not read surface %s",
                  Progname, surfname) ;
      printf("Done reading surface\n");

      /* make a surface mask in the volume */
      Vxyz = MatrixAlloc(4,1,MATRIX_REAL);
      Vxyz->rptr[4][1] = 1.0;
      for (n=0; n < surf->nvertices; n++) {
        Vxyz->rptr[1][1] = surf->vertices[n].x;
        Vxyz->rptr[2][1] = surf->vertices[n].y;
        Vxyz->rptr[3][1] = surf->vertices[n].z;
        Vcrs = MatrixMultiply(ras2vox_targ,Vxyz,Vcrs);
        c = nint(Vcrs->rptr[1][1]);
        r = nint(Vcrs->rptr[2][1]);
        s = nint(Vcrs->rptr[3][1]);
        if (c < 0 || c > 255 || r < 0 || r > 255 || s < 0 || s > 255) {
          printf("ERROR: vertex %d (%g,%g,%g) is out of range (%d,%d,%d)\n",
                 n,surf->vertices[n].x,surf->vertices[n].y,surf->vertices[n].z,
                 c,r,s);
          exit(1);
        }
        MRIvox(mrisurf,c,r,s) = 2;
      }
      MRISfree(&surf);
    }

  }

  /*------------------------------------------------------*/
  if ( abs(ps_2 - mov_vol->xsize) > .001 ) {
    printf("WARNING: pixel size in regfile (%g) does not match that of "
           "movable volume (%g)\n",ps_2,mov_vol->xsize);
    printf("If the movable volume is a bshort or bfloat, make sure\n");
    printf("that a .bhdr file exists\n");
  }
  if ( abs(st_2 - mov_vol->zsize) > .001 ) {
    printf("WARNING: slice thickness in regfile (%f) does not match that of "
           "movable volume (%f)\n",st_2,mov_vol->zsize);
    printf("If the movable volume is a bshort or bfloat, make sure\n");
    printf("that a .bhdr file exists\n");
  }

  xdim_2 = mov_vol->width;
  ydim_2 = mov_vol->height;
  imnr1_2 = mov_vol->depth-1;

  xnum = targ_vol->width;
  ynum = targ_vol->height;
  numimg = targ_vol->depth;

  xdim = WINDOW_COLS;
  ydim = WINDOW_ROWS;
  bufsize = xdim*ydim;

  imnr0 = 0;
  imnr1 = 255;

  zf = (int)nint((float)xdim/xnum); /* zoom factor */
  ozf = zf;
  fsf = (float)zf;
  printf("Zoom Factor = %g, SQR() = %g\n",(float)zf,(float)SQR(zf));
  printf("FOV = %g\n",FOV);

  vidbuf    = (GLubyte *)lcalloc(3*bufsize*SQR(zf),sizeof(GLubyte));
  blinkbuft = (GLubyte *)lcalloc(3*bufsize*SQR(zf),sizeof(GLubyte));
  blinkbufm = (GLubyte *)lcalloc(3*bufsize*SQR(zf),sizeof(GLubyte));


  ps = targ_vol->xsize;
  st = targ_vol->zsize;
  xx0 = -128.0;
  xx1 = +128.0;
  yy0 = -128.0;
  yy1 = +128.0;
  zz0 = -128.0;
  zz1 = +128.0;

  binbuff = (unsigned char *)calloc(3*xdim*ydim,sizeof(char));

  if (tkrtitle == NULL) tkrtitle = subjectid;

  printf("Opening window %s\n",pname);
  fflush(stdout);
  open_window(tkrtitle);

  printf("Setting scale\n");
  fflush(stdout);
  set_scale();

  imc = zf*imnr1/2; // Cor
  ic = ydim/2;      // Ax
  jc = xdim/2;      // Sag
  if (slice_init > 0) {
    if (plane_init == CORONAL)    imc = zf*slice_init;
    if (plane_init == SAGITTAL)   jc  = (int)(xdim*(slice_init/256.0));
    if (plane_init == HORIZONTAL) ic  = (int)(ydim*(slice_init/256.0));
    printf("plane = %d, slice = %d\n",plane_init,slice_init);
    printf("imc = %d, jc = %d, ic = %d\n",imc,jc,ic);
  }

  if (overlay_mode_init == TARGET) {
    overlay_mode = MOVEABLE;
    redraw();
    overlay_mode = TARGET;
    redraw();
  } else {
    overlay_mode = TARGET;
    redraw();
    overlay_mode = MOVEABLE;
    redraw();
  }

  if (mov_scale != 1) {
    printf("movscale = %g\n",mov_scale);
    scale_brain(1.0/mov_scale,'x');
    scale_brain(1.0/mov_scale,'y');
    scale_brain(1.0/mov_scale,'z');
  }
  if (tagmov) {
    printf("Tagging mov volume\n");
    MRItagVol(mov_vol, 10000);
  }

  updateflag = TRUE;

  return GL_FALSE;
#endif // HAVE_TCL_TK_GL
}


/****----------*********----------****----------*********----------*/
/* ------------------------------------------------------------------ */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;
  char *fmt;
  char *errstr;

  isflag(" "); /* shuts up compiler */

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
    else if (!strcasecmp(option, "--fstarg"))    fstarg = 1;
    else if (!strcasecmp(option, "--nofstarg"))  fstarg = 0;
    else if (!strcasecmp(option, "--nofix"))     fixtkreg = 0;
    else if (!strcasecmp(option, "--fixonly"))   fixonly = 1;
    else if (!strcasecmp(option, "--inorm"))     use_inorm = 1;
    else if (!strcasecmp(option, "--no-inorm"))  use_inorm = 0;
    else if (!strcasecmp(option, "--regheader")) mkheaderreg = 1;
    else if (!strcasecmp(option, "--regheader-center")){
      // This performs the header registration but computes the shift
      // to align the volume centers
      mkheaderreg = 1;
      mkheaderregCenter = 1;
    }
    else if (!strcasecmp(option, "--identity"))  identityreg = 1;
    else if (!strcasecmp(option, "--noedit"))    noedit = 1;
    else if (!strcasecmp(option, "--zero-cras"))     ZeroCRAS = 1;
    else if (!strcasecmp(option, "--no-zero-cras"))  ZeroCRAS = 0;
    else if (!strcasecmp(option, "--fmov-targ"))  DoFMovTarg = 1;
    else if (!strcasecmp(option, "--fstal")) {
      fstal = 1;
      LoadSurf = 0;
      UseSurf = 0;
      fscale_2 = 1;
      ZeroCRAS = 1;
    } 
    else if (!strcasecmp(option, "--fstal-targ")) {
      if (nargc < 1) argnerr(option,1);
      fstaltarg = pargv[0];
      nargsused = 1;
      fstal = 1;
      LoadSurf = 0;
      UseSurf = 0;
      fscale_2 = 1;
      ZeroCRAS = 1;
    } 
    else if (!strcasecmp(option, "--fstal-avi")) {
      fstaltarg = "711-2C_as_mni_average_305.4dfp.img";
      fstal = 1;
      LoadSurf = 0;
      UseSurf = 0;
      fscale_2 = 1;
      ZeroCRAS = 1;
    } 
    else if (!strcasecmp(option, "--fixxfm"))    fixxfm = 1;
    else if (!strcasecmp(option, "--nofixxfm"))  fixxfm = 0;
    else if (!strcasecmp(option, "--tag"))    tagmov = 1;
    else if (!strcasecmp(option, "--notag"))  tagmov = 0;
    else if (!strcasecmp(option, "--mgz"))  ; // for backwards compat
    else if (stringmatch(option, "--fsl-targ")) {
      sprintf(tmpstr,"%s/etc/standard/avg152T1",getenv("FSLDIR"));
      printf("Trying %s\n",tmpstr);
      targ_vol_id = IDnameFromStem(tmpstr); // For FSL 4.0
      if(! fio_FileExistsReadable(targ_vol_id)){
	sprintf(tmpstr,"%s/data/standard/avg152T1",getenv("FSLDIR"));
	printf("Trying %s\n",tmpstr);
	targ_vol_id = IDnameFromStem(tmpstr); // For FSL 4.0
	if(targ_vol_id == NULL) exit(1);
      }
    } 
    else if (stringmatch(option, "--fsl-targ-lr")) {
      sprintf(tmpstr,"%s/etc/standard/avg152T1_LR-marked",getenv("FSLDIR"));
      printf("Trying %s\n",tmpstr);
      targ_vol_id = IDnameFromStem(tmpstr); // For FSL 4.0
      if(! fio_FileExistsReadable(targ_vol_id)){
	sprintf(tmpstr,"%s/data/standard/avg152T1_LR-marked",getenv("FSLDIR"));
	printf("Trying %s\n",tmpstr);
	targ_vol_id = IDnameFromStem(tmpstr); // For FSL 4.0
	if(targ_vol_id == NULL) exit(1);
      } 
    }
    else if (!strcasecmp(option, "--lh-only")){ lhsurf_only=1 ; LoadSurf = 1; UseSurf  = 1;}
    else if (!strcasecmp(option, "--rh-only")){ rhsurf_only=1 ; LoadSurf = 1; UseSurf  = 1;}
    else if (!strcasecmp(option, "--check-reg") ||
	     !strcasecmp(option, "--check") ||
	     !strcasecmp(option, "--junk")){
      sprintf(tmpstr,"/tmp/reg.tmp.%ld.dat",PDFtodSeed());
      regfname = strcpyalloc(tmpstr);
      checkreg = 1;
    }
    else if (!strcasecmp(option, "--2")){
      WINDOW_ROWS = 2*512;
      WINDOW_COLS = 2*512;
    }
    else if (!strcasecmp(option, "--size")){
      double size;
      sscanf(pargv[0],"%lf",&size);
      WINDOW_ROWS = nint(size*512);
      WINDOW_COLS = nint(size*512);
      nargsused = 1;
    }
    else if (stringmatch(option, "--targ")) {
      if (nargc < 1) argnerr(option,1);
      targ_vol_id = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        fmt = pargv[1];
        nargsused ++;
        targ_vol_fmt = checkfmt(fmt);
      }
    } 
    else if (!strcasecmp(option, "--conf-targ")) conformTarget = 1;
    else if (!strcmp(option, "--seg")) {
      if (nargc < 1) argnerr(option,1);
      seg_vol_id = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--aseg")) DoASeg = 1;
    else if (!strcmp(option, "--aparc+aseg")) DoAParcASeg = 1;
    else if (!strcmp(option, "--wmparc")) DoWMParc = 1;
    else if (!strcmp(option, "--targ-orientation")) {
      if (nargc < 1) argnerr(option,1);
      targ_ostr = pargv[0];
      errstr = MRIcheckOrientationString(targ_ostr);
      if (errstr) {
        printf("ERROR: with target orientation string %s\n",targ_ostr);
        printf("%s\n",errstr);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--mov")) {
      if (nargc < 1) argnerr(option,1);
      mov_vol_id = pargv[0];
      nargsused = 1;
      if (nth_is_arg(nargc, pargv, 1)) {
        fmt = pargv[1];
        nargsused ++;
        mov_vol_fmt = checkfmt(fmt);
      }
    } else if (!strcmp(option, "--int")) {
      if (nargc < 2) argnerr(option,2);
      int_vol_id = pargv[0];
      int_regfname = pargv[1];
      nargsused = 2;
    } else if (!strcmp(option, "--mov-orientation")) {
      if (nargc < 1) argnerr(option,1);
      mov_ostr = pargv[0];
      errstr = MRIcheckOrientationString(mov_ostr);
      if (errstr) {
        printf("ERROR: with mov orientation string %s\n",mov_ostr);
        printf("%s\n",errstr);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--movbright") || !strcmp(option, "--fmov")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fscale_2);
      use_inorm = 0;
      nargsused = 1;
    } else if (!strcmp(option, "--movscale")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&mov_scale);
      nargsused = 1;
    } else if (!strcmp(option, "--movframe")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&mov_frame);
      nargsused = 1;
    } else if (!strcmp(option, "--slice")) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&slice_init);
      nargsused = 1;
    } else if (!strcmp(option, "--plane")) {
      if (nargc < 1) argnerr(option,1);
      plane_init = -1000;
      if ( strcmp(pargv[0],"cor") == 0) plane_init = CORONAL;
      if ( strcmp(pargv[0],"sag") == 0) plane_init = SAGITTAL;
      if ( strcmp(pargv[0],"ax")  == 0) plane_init = HORIZONTAL;
      if ( strcmp(pargv[0],"hor") == 0) plane_init = HORIZONTAL;
      if (plane_init == -1000) {
        printf("ERROR: orientation %s unrecognized\n",pargv[0]);
        exit(1);
      }
      nargsused = 1;
    } else if (!strcmp(option, "--volview")) {
      if (nargc < 1) argnerr(option,1);
      overlay_mode_init = -1000;
      if ( strcmp(pargv[0],"mov") == 0)  overlay_mode_init = MOVEABLE;
      if ( strcmp(pargv[0],"targ") == 0) overlay_mode_init = TARGET;
      if (overlay_mode_init == -1000) {
        printf("ERROR: volview %s unrecognized\n",pargv[0]);
        exit(1);
      }
      nargsused = 1;
    } 
    else if (!strcmp(option, "--surf") || !strcmp(option, "--surfs")) {
      LoadSurf = 1;
      UseSurf  = 1;
      nargsused = 0;
      if (nth_is_arg(nargc, pargv, 0)) {
        surfname = pargv[0];
        nargsused ++;
        printf("surfname set to %s\n",surfname);
      }
    } 
    else if (!strcmp(option, "--surf-rgb")) {
      if(nargc < 3) argnerr(option,3);
      sscanf(pargv[0],"%d",&SurfRGB[0]);
      sscanf(pargv[1],"%d",&SurfRGB[1]);
      sscanf(pargv[2],"%d",&SurfRGB[2]);
      nargsused = 3;
    }
    else if (!strcmp(option, "--talxfmname")) {
      if (nargc < 1) argnerr(option,1);
      talxfmname = pargv[0];
      fstal = 1;
      LoadSurf = 0;
      UseSurf = 0;
      fscale_2 = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--xfm")) {
      if (nargc < 1) argnerr(option,1);
      xfmfname = pargv[0];
      printf("INFO: reading xfm file %s, trying as MINC xfm \n",xfmfname);
      err = regio_read_mincxfm(xfmfname, &XFM,NULL);
      if (err) exit(1);
      mkheaderreg = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--ixfm")) {
      MATRIX *m_tmp ;
      if (nargc < 1) argnerr(option,1);
      xfmfname = pargv[0];
      printf("INFO: reading xfm file %s, trying as MINC xfm \n",xfmfname);
      err = regio_read_mincxfm(xfmfname, &m_tmp,NULL);
      if (err) exit(1);
      printf("inverting registration\n") ;
      XFM = MatrixInverse(m_tmp, NULL) ; MatrixFree(&m_tmp) ;
      mkheaderreg = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--reg")) {
      if (nargc < 1) argnerr(option,1);
      regfname = pargv[0];
      nargsused = 1;
    } 
    else if(!strcmp(option, "--lta") || !strcmp(option, "--lta-inv")) {
      // Having a separate flag for --lta allows conversion
      if (nargc < 1) argnerr(option,1);
      ltafname = pargv[0];
      FSXform = TransformRead(ltafname);
      if(FSXform == NULL) exit(1);
      lta = (LTA*) FSXform->xform;
      if(lta->type != LINEAR_RAS_TO_RAS){
        printf("INFO: LTA input is not RAS to RAS...converting...\n");
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      } 
     if(lta->type != LINEAR_RAS_TO_RAS){
        printf("ERROR: LTA input is not RAS to RAS\n");
        exit(1);
      }
      linxfm = &(lta->xforms[0]);
      // Assume RAS2RAS and uses vox2ras from input volumes:
      // Note: This ignores the volume geometry in the LTA file.
      XFM = lta->xforms[0].m_L;
      printf("lta->subject %s\n",lta->subject);
      if(lta->subject[0] != 0) {
	memset(subjectid,'\0',1000);
	memmove(subjectid,lta->subject,strlen(lta->subject));
      }
      if(!strcmp(option, "--lta-inv")){
	printf("Inverting LTA\n");
	XFM = MatrixInverse(XFM,XFM);
      }
      mkheaderreg = 1;
      nargsused = 1;
    } 
    else if (istringnmatch(option, "--rot",0)) {
      if (nargc < 3) argnerr(option,3);
      // Angles are in degrees
      sscanf(pargv[0],"%lf",&angles[0]);
      sscanf(pargv[1],"%lf",&angles[1]);
      sscanf(pargv[2],"%lf",&angles[2]);
      angles[0] *= (M_PI/180);
      angles[1] *= (M_PI/180);
      angles[2] *= (M_PI/180);
      Mrot = MRIangles2RotMat(angles);
      nargsused = 3;
    } 
    else if (istringnmatch(option, "--trans",0)) {
      if (nargc < 3) argnerr(option,3);
      // Translation in mm
      sscanf(pargv[0],"%lf",&xyztrans[0]);
      sscanf(pargv[1],"%lf",&xyztrans[1]);
      sscanf(pargv[2],"%lf",&xyztrans[2]);
      Mtrans = MatrixIdentity(4,NULL);
      Mtrans->rptr[1][4] = xyztrans[0];
      Mtrans->rptr[2][4] = xyztrans[1];
      Mtrans->rptr[3][4] = xyztrans[2];
      nargsused = 3;
    } 
    else if (!strcmp(option, "--det")) {
      if (nargc < 1) argnerr(option,1);
      DetFile = pargv[0];
      nargsused = 1;
    } 
    else if (stringmatch(option, "--fsfeat")) {
      if (nargc < 1) argnerr(option,1);
      //pargv[0] is featdir
      sprintf(tmpstr,"%s/reg/freesurfer/register.dat",pargv[0]);
      if(! fio_FileExistsReadable(tmpstr)){
	printf("ERROR: cannot find %s\n",tmpstr);
	exit(1);
      }
      regfname = strcpyalloc(tmpstr);
      sprintf(tmpstr,"%s/example_func",pargv[0]);
      mov_vol_id = IDnameFromStem(tmpstr);
      nargsused = 1;
    } 
    else if (stringmatch(option, "--feat")) {
      if (nargc < 1) argnerr(option,1);
      //pargv[0] is featdir
      sprintf(tmpstr,"%s/etc/standard/avg152T1",getenv("FSLDIR"));
      if(! fio_FileExistsReadable(tmpstr))
	sprintf(tmpstr,"%s/data/standard/avg152T1",getenv("FSLDIR"));
      targ_vol_id = IDnameFromStem(tmpstr); // For FSL 4.0
      sprintf(tmpstr,"%s/example_func",pargv[0]);
      mov_vol_id = IDnameFromStem(tmpstr);
      sprintf(tmpstr,"%s/reg/example_func2standard.mat",pargv[0]);
      fslregfname = strcpyalloc(tmpstr);
      read_fslreg(fslregfname);
      fslregoutfname = fslregfname;
      sprintf(tmpstr,"/tmp/feat.exf2std.reg.%d",(int)PDFtodSeed());
      regfname = strcpyalloc(tmpstr);
      tagmov = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--fslreg") || !strcmp(option, "--fsl")) {
      if (nargc < 1) argnerr(option,1);
      fslregfname = pargv[0];
      read_fslreg(fslregfname);
      nargsused = 1;
    } 
    else if (!strcmp(option, "--gca-skull")) {
      if(nargc < 1) argnerr(option,1);
      sprintf(subjectid,"%s",pargv[0]);
      sprintf(tmpstr,"%s/%s/mri/transforms/talairach_with_skull.lta",subjectsdir,subjectid);
      FSXform = TransformRead(tmpstr);
      if(FSXform == NULL) exit(1);
      ltaoutfname = strcpyalloc(tmpstr);
      lta = (LTA*) FSXform->xform;
      if(lta->type != LINEAR_RAS_TO_RAS){
        printf("INFO: LTA input is not RAS to RAS...converting...\n");
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      }
      if(lta->type != LINEAR_RAS_TO_RAS){
        printf("ERROR: LTA input is not RAS to RAS\n");
        exit(1);
      }
      linxfm = &(lta->xforms[0]);
      // Assume RAS2RAS and uses vox2ras from input volumes:
      // Note: This ignores the volume geometry in the LTA file.
      XFM = MatrixInverse(lta->xforms[0].m_L,NULL);
      sprintf(tmpstr,"%s/%s/mri/nu.mgz",subjectsdir,subjectid);
      targ_vol_id = strcpyalloc(tmpstr);
      mov_vol_id = strcpyalloc(lta->xforms->dst.fname);
      if(! fio_FileExistsReadable(mov_vol_id)){
	sprintf(tmpstr,"%s/average/RB_all_withskull_2008-03-26.gca",
		getenv("FREESURFER_HOME"));
	free(mov_vol_id);
	mov_vol_id = strcpyalloc(tmpstr);
      }
      if(regfname == NULL){
	sprintf(tmpstr,"/tmp/tkregister2.%s.junk",subjectid);
	regfname = strcpyalloc(tmpstr);
      }
      mkheaderreg = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--gca")) {
      if(nargc < 1) argnerr(option,1);
      sprintf(subjectid,"%s",pargv[0]);
      sprintf(tmpstr,"%s/%s/mri/transforms/talairach.lta",subjectsdir,subjectid);
      FSXform = TransformRead(tmpstr);
      if(FSXform == NULL) exit(1);
      lta = (LTA*) FSXform->xform;
      if(lta->type != LINEAR_RAS_TO_RAS){
        printf("INFO: LTA input is not RAS to RAS...converting...\n");
        lta = LTAchangeType(lta, LINEAR_RAS_TO_RAS);
      }
      if(lta->type != LINEAR_RAS_TO_RAS){
        printf("ERROR: LTA input is not RAS to RAS\n");
        exit(1);
      }
      linxfm = &(lta->xforms[0]);
      // Assume RAS2RAS and uses vox2ras from input volumes:
      // Note: This ignores the volume geometry in the LTA file.
      XFM = MatrixInverse(lta->xforms[0].m_L,NULL);
      sprintf(tmpstr,"%s/%s/mri/T1.mgz",subjectsdir,subjectid);
      targ_vol_id = strcpyalloc(tmpstr);
      mov_vol_id = strcpyalloc(lta->xforms->dst.fname);
      mkheaderreg = 1;
      nargsused = 1;
    } 
    else if (!strcmp(option, "--vox2vox")){
      if (nargc < 1) argnerr(option,1);
      Vox2VoxFName = pargv[0];
      Vox2Vox = Load4x4(Vox2VoxFName);
      if(!Vox2Vox) exit(1);
      printf("Vox2Vox Matrix \n");
      MatrixPrint(stdout,Vox2Vox);
      mkheaderreg = 1;
      nargsused = 1;
    } else if (!strcmp(option, "--ltaout")) {
      if(nargc < 1) argnerr(option,1);
      ltaoutfname = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--ltaout-inv")) {
      // has no effect without --ltaout
      invLTAOut = 1;
    } 
    else if (!strcmp(option, "--fslregout")) {
      if (nargc < 1) argnerr(option,1);
      fslregoutfname = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--freeview")) {
      if (nargc < 1) argnerr(option,1);
      freeviewfname = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--xfmout")) {
      if (nargc < 1) argnerr(option,1);
      xfmoutfname = pargv[0];
      nargsused = 1;
    } else if (!strcmp(option, "--float2int")) {
      if (nargc < 1) argnerr(option,1);
      float2int_use = float2int_code(pargv[0]);
      if (float2int_use == -1) {
        printf("ERROR: float2int method %s unrecognized\n",pargv[0]);
        exit(1);
      }
      nargsused = 1;
    } else if ( !strcmp(option, "--subject") ||
                !strcmp(option, "--s") ) {
      if (nargc < 1) argnerr(option,1);
      memset(subjectid,'\0',1000);
      memmove(subjectid,pargv[0],strlen(pargv[0]));
      subjectidOverride = 1;
      nargsused = 1;
    } else if ( !strcmp(option, "--sd") ) {
      if (nargc < 1) argnerr(option,1);
      setenv("SUBJECTS_DIR", pargv[0], 1);
      nargsused = 1;
    } else if ( !strcmp(option, "--title") ) {
      if (nargc < 1) argnerr(option,1);
      tkrtitle = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--tcl") ) {
      if(nargc < 1) argnerr(option,1);
      tkregister_tcl = pargv[0];
      nargsused = 1;
    } 
    else if ( !strcmp(option, "--gdiagno") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&Gdiag_no);
      nargsused = 1;
    } else if ( !strcmp(option, "--contrast") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fsquash);
      nargsused = 1;
    } else if ( !strcmp(option, "--midpoint") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%lf",&fthresh);
      nargsused = 1;
    } else if ( !strcmp(option, "--fov") ) {
      if (nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&FOV);
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --help : usage and documentation\n");
  printf("\n");
  printf("   --mov  movable volume  <fmt> \n");
  printf("   --targ target volume <fmt>\n");
  printf("   --fstarg : target is relative to subjectid/mri\n");
  printf("   --reg  register.dat : input/output registration file\n");
  printf("   --check-reg : only check, no --reg needed\n");
  printf("   --regheader : compute regstration from headers\n");
  printf("   --regheader-center : same as --regheader but aligns volume centers\n");
  printf("   --fsl-targ : use FSLDIR/data/standard/avg152T1.nii.gz\n");
  printf("   --fsl-targ-lr : use FSLDIR/data/standard/avg152T1_LR-marked.nii.gz\n");
  printf("   --gca subject : check linear GCA registration  \n");
  printf("   --gca-skull subject : check linear 'with skull' GCA registration  \n");
  printf("   --no-zero-cras : do not zero target cras (done with --fstal)\n");
  printf("   --movbright  f : brightness of movable volume\n");
  printf("   --no-inorm  : turn off intensity normalization\n");
  printf("   --fmov fmov : set mov brightness \n");
  printf("   --fmov-targ : apply fmov brightness to the target\n");
  printf("   --plane  orient  : startup view plane <cor>, sag, ax\n");
  printf("   --slice  sliceno : startup slice number\n");
  printf("   --volview volid  : startup with targ or mov\n");
  printf("   --fov FOV  : window FOV in mm (default is 256)\n");
  printf("   --movscale scale : scale size of mov by scale\n");
  printf("   --surf surfname : display surface as an overlay \n");
  printf("   --surf-rgb R G B : set surface color (0-255) \n");
  printf("   --lh-only : only load/display left hemi \n");
  printf("   --rh-only : only load/display right hemi \n");
  printf("   --fstal : set mov to be tal and reg to be taliarach.xfm  \n");
  printf("   --talxfmname talxfmname : set mov to be tal and reg to be talxfmname  \n");
  printf("   --ixfm file : MNI-style inverse registration input matrix\n");
  printf("   --xfm file : MNI-style registration input matrix\n");
  printf("   --xfmout file : MNI-style registration output matrix\n");
  printf("   --fsl file : FSL-style registration input matrix\n");
  printf("   --fslregout file : FSL-Style registration output matrix\n");
  printf("   --freeview file : FreeView registration output matrix\n");
  printf("   --vox2vox file : vox2vox matrix in ascii\n");
  printf("   --lta ltafile : Linear Transform Array\n");
  printf("   --lta-inv ltafile : Read in LTA and invert\n");
  printf("   --ltaout ltaoutfile : Output a Linear Transform Array\n");
  printf("   --ltaout-inv : invert transform in ltaoutfile\n");
  printf("   --feat featdir : check example_func2standard registration\n");
  printf("   --fsfeat featdir : check reg/freesurfer/register.dat registration\n");
  printf("   --identity : use identity as registration matrix\n");
  printf("   --s subjectid : set subject id \n");
  printf("   --sd dir : use dir as SUBJECTS_DIR\n");
#ifdef HAVE_TCL_TK_GL  
  printf("   --noedit : do not open edit window (exit) - for conversions\n");
#endif // HAVE_TCL_TK_GL
  printf("   --nofix : don't fix old tkregister matrices\n");
  printf("   --float2int code : spec old tkregister float2int\n");
  printf("   --title title : set window title\n");
  printf("   --tag : tag mov vol near the col/row origin\n");
  printf("   --mov-orientation ostring : supply orientation string for mov\n");
  printf("   --targ-orientation ostring : supply orientation string for targ\n");
  printf("   --int intvol intreg : use registration from intermediate volume \n");
  printf("   --2 : double window size \n");
  printf("   --size scale : scale window by scale (eg, 0.5, 1.5) \n");
  printf("   --det  detfile : save determinant of reg mat here\n");
  printf("   --aseg : load aseg (hit 'd' to toggle)\n");
  printf("   --aparc+aseg : load aparc+aseg (hit 'c' to toggle)\n");
  printf("   --wmparc : load wmparc (hit 'c' to toggle)\n");
  printf("   --gdiagno n : set debug level\n");
  printf("   --trans Tx Ty Tz : translation (mm) to apply to reg matrix\n");
  printf("   --rot   Ax Ay Az : rotation angles (deg) to apply to reg matrix\n");
  printf("   --conf-targ : conform target (assumes reg computed to conf target, eg, GCA)\n");
  printf("\n");
  //printf("   --svol svol.img (structural volume)\n");
}


/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;

  printf("\n");
  std::cout << getVersion() << std::endl;

  printf(

    "SUMMARY\n"
    "\n"
    "Note: consider using tkregisterfv which is a frontend for freeview\n"
    "\n"
    "tkregister2 is a tool to assist in the manual tuning of the linear\n"
    "registration between two volumes, mainly for the purpose of\n"
    "interacting with the FreeSurfer anatomical stream.\n"
    "The output register.dat maps from target tkregRAS to mov tkregRAS\n"
    "Note that this is reverse of what you might expect\n\n"
    "The GUI displays a\n"
    "window in which the user can toggle between the two volumes "
    "to see how\n"
    "well they line up. It is also possible to display "
    "the cortical surface\n"
    "to assist in alignment. The user can edit the registration by\n"
    "translation, rotation, and stretching.  "
    "The initial registration can\n"
    "be based on: (1) a previously existing registration, "
    "(2) the spatial\n"
    "information contained in the headers of the two input volumes, "
    "or (3)\n"
    "an FSL registration matrix. The output takes the form of a "
    "FreeSurfer\n"
    "registration matrix. It is also possible to output an "
    "FSL registration\n"
    "matrix. It is possible to run tkregister2 without the "
    "manual editing\n"
    "stage (ie, it exits immediately) in order to compute "
    "and/or convert\n"
    "registration matrices. tkregister2 can also compute registration\n"
    "matrices from two volumes aligned in SPM or FSL.\n"
    "\n"
    "  \n"
    "FLAGS AND OPTIONS\n"
    "  \n"
    "  --targ target-volume-id\n"
    "  \n"
    "  This is the path to the target volume. It can read "
    "anything readable \n"
    "  by mri_convert. See also --fstarg.\n"
    "  \n"
    "  --fstarg\n"
    "  \n"
    "  This flag implicitly sets the target volume to be the "
    "T1 volume found\n"
    "  in the subjects directory under the subject name found in the \n"
    "  input registration matrix. \n"
    "  \n"
    "  --mov movable-volume-id\n"
    "  \n"
    "  This is the path to the target volume. "
    "It can read anything readable \n"
    "  by mri_convert. See also --fstal.\n"
    "  \n"
    "  --fstal\n"
    "\n"
    "  Check and edit the talairach registration that was created during\n"
    "  the FreeSurfer reconstruction. Sets the movable volume to be \n"
    "  $FREESURFER_HOME/average/mni305.cor.mgz and sets the registration file to be\n"
    "  $SUBJECTS_DIR/subjectid/transforms/talairach.xfm. User must have\n"
    "  write permission to this file. Do not specify --reg with this\n"
    "  flag. It is ok to specify --regheader with this flag. The format\n"
    "  of the anatomical is automatically detected as mgz or COR. By default,\n"
    "  the target c_ras is temporarily set to 0 to assure that the target\n"
    "  is properly centered. This is taken into account when computing \n"
    "  and writing the output xfm. To turn this off, add --no-zero-cras.\n"
    "\n"
    "  --plane <orientation>\n"
    "\n"
    "  Set the initial orientation. Valid values are cor, sag, and ax.\n"
    "  The default is cor.\n"
    "\n"
    "  --slice sliceno\n"
    "\n"
    "  Set the initial slice to view. \n"
    "\n"
    "  --fov FOV\n"
    "\n"
    "  Set the view port field-of-view. Default is 256. Note, the zoom\n"
    "  can also be controlled interactively with - and =.\n"
    "\n"
    "  --movscale scale\n"
    "\n"
    "  Adjust registration matrix to scale mov volume by scale. "
    "This has\n"
    "  the same effect as adjusting the scale manually. It also flags\n"
    "  the matrix as being edited, and so you will be asked to save it\n"
    "  even if you have not made any manual edits.\n"
    "\n"
    "  --movframe frame\n"
    "\n"
    "  Use frame from moveable volume.\n"
    "\n"
    "  --surf <surfacename>\n"
    "\n"
    "  Load the cortical surface (both hemispheres) and display as an\n"
    "  overlay on the volumes. If just --surf is supplied, "
    "then the white \n"
    "  surface is loaded, otherwise surfacename is loaded. The subject\n"
    "  name (as found in the input registration file or as specified\n"
    "  on the command-line) must be valid, and the surface must exist\n"
    "  in the subject's anatomical directory. The surface display can\n"
    "  be toggled by clicking in the window and hitting 's'.\n"
    "\n"
    "  --reg registration-file\n"
    "  \n"
    "  Path to the input or output FreeSurfer registration file. "
    "If the user\n"
    "  specifies that the initial registration be computed from "
    "the input\n"
    "  headers or from an FSL matrix, then this will be the output "
    "file.\n"
    "  Otherwise, the initial registration will be read from this file\n"
    "  (it will still be the output file). Note: a FreeSurfer "
    "registration\n"
    "  file must be specified even if you only want the "
    "FSL registration\n"
    "  file as the output. See also --fstal.\n"
    "  \n"
    "  --regheader\n"
    "  \n"
    "  Compute the initial registration from the information "
    "in the headers\n"
    "  of the input volumes. The contents of the registration file, "
    "if it\n"
    "  exists, will be ignored. This mostly makes sense when the "
    "two volumes\n"
    "  were acquired at the same time (ie, when the head motion "
    "is minimal).\n"
    "  \n"
    "  --fsl FSL-registration-file\n"
    "  \n"
    "  Use the matrix produced by the FSL routines as the initial registration.\n"
    "  It should be an ascii file with a 4x4 matrix. Note: the matrix should\n"
    "  map from the mov to the target. See also --feat and --fslregout.\n"
    "  \n"
    "  --xfm MNI-Style registration matrix\n"
    "  \n"
    "  Use the matrix produced by an MNI prgram as the "
    "  initial registration.\n"
    "  Note: the matrix should map from the mov to the target.\n"
    "  \n"
    "  --lta ltafile \n"
    "  \n"
    "  RAS-to-RAS linear transform array file."
    "  Note: the matrix should map from the mov to the target.\n"
    "  \n"
    "  --ltaout ltaoutfile \n"
    "  \n"
    "  RAS-to-RAS linear transform array file."
    "  \n"
    "  --vox2vox vox2voxfile \n"
    "  \n"
    "  Input registration is a vox2vox ascii file. Vox2Vox maps target\n"
    "  indices (c,r,s) to mov indices. Lines that begin with '#' are ignored.\n"
    "  \n"
    "  --s identity\n"
    "  \n"
    "  Use identity as input registration. Same as simply creating\n"
    "  the identity matrix in a register.dat.\n"
    "  \n"
    "  --s subjectid\n"
    "  \n"
    "  Subject identifier string that will be printed in the "
    "output registration\n"
    "  file when no the input registration file is not specified "
    "(ie, with\n"
    "  --regheader or --fsl).\n"
    "  \n"
    "  --sd subjectsdir\n"
    "  \n"
    "  Set the path to the parent directory of the FreeSurfer "
    "anatomical \n"
    "  reconstructions. If unspecified, this will default to the "
    "SUBJECTS_DIR \n"
    "  environment variable. This only has an effect with --fstarg.\n"
    "  \n"
    "  --noedit\n"
    "  \n"
    "  Do not bring up the GUI, just print out file(s) and exit. "
    "This is mainly\n"
    "  useful saving the header-computer registration matrix "
    "(or for converting\n"
    "  from FreeSurfer to FSL, or the other way around).\n"
    "  \n"
    "  --fslregout FSL-registration-file\n"
    "  \n"
    "  Compute an FSL-compatible registration matrix based on either the\n"
    "  FreeSurfer matrix or the header. This can be helpful for initializing\n"
    "  the FSL registration routines.\n"
    "  \n"
    "  --feat featdir\n"
    "  \n"
    "  View/modify the FSL FEAT registration to standard space (ie,\n"
    "  example_func2standard.mat. Manual edits to the registration\n"
    "  will change this file. There will also be a 'tag', ie, a grid of\n"
    "  grid dots near the col-row-slice origin of the example_func.\n"
    "  This might be helpful for determining if the reg is left-right\n"
    "  reversed.\n"
    "  \n"
    "  --nofix\n"
    "  \n"
    "  This is only here for debugging purposes in order to "
    "simulate the behavior\n"
    "  of the old tkregister.\n"
    "  \n"
    "  --float2int code\n"
    "  \n"
    "  This is only here for debugging purposes in order to "
    "simulate the \n"
    "  behavior of the old tkregister.\n"
    "  \n"
    "  --title title\n"
    "  \n"
    "  Set the window titles to title. Default is subjectname.\n"
    "  \n"
    "  --tag\n"
    "  \n"
    "  Creates a hatched pattern in the mov volume prior to resampling.\n"
    "  This pattern is in all slices near the col/row origin (ie, near\n"
    "  col=0,row=0). This can help to determine if there is a "
    "left-right\n"
    "  reversal. Think of this as a synthetic fiducial. Can be good in \n"
    "  combination with --mov-orientation.\n"
    "  \n"
    "  --mov-orientation ostring\n"
    "  --targ-orientation ostring\n"
    "  \n"
    "  Supply the orientation information in the form of an "
    "orientation string \n"
    "  (ostring). The ostring is three letters that roughly "
    "describe how the volume\n"
    "  is oriented. This is usually described by the direction "
    "cosine information\n"
    "  as originally derived from the dicom but might not be "
    "available in all data\n"
    "  sets. --mov-orientation will have no effect unless "
    "--regheader is specified.\n"
    "  The first  character of ostring determines the "
    "direction of increasing column.\n"
    "  The second character of ostring determines the "
    "direction of increasing row.\n"
    "  The third  character of ostring determines the "
    "direction of increasing slice.\n"
    "  Eg, if the volume is axial starting inferior and "
    "going superior the slice \n"
    "  is oriented such that nose is pointing up and the "
    "right side of the subject\n"
    "  is on the left side of the image, then this would "
    "correspond to LPS, ie,\n"
    "  as the column increases, you move to the patients left; "
    "as the row increases,\n"
    "  you move posteriorly, and as the slice increases, "
    "you move superiorly. Valid\n"
    "  letters are L, R, P, A, I, and S. There are "
    "48 valid combinations (eg, RAS\n"
    "  LPI, SRI). Some invalid ones are DPS (D is not a valid letter), "
    "RRS (can't\n"
    "  specify R twice), RAP (A and P refer to the same axis). "
    "Invalid combinations\n"
    "  are detected immediately, an error printed, and the "
    "program exits. Case-\n  insensitive.\n"
    "  \n"
    "  \n"
    "  --int volid regmat\n"
    "\n"
    "  Use registration from an intermediate volume. This can be "
    "useful when the\n"
    "  FOV of the moveable volume does not cover the entire brain. "
    "In this case, \n"
    "  you can register a full-volume COLLECTED IN THE SAME "
    "SESSION AS THE MOVEABLE\n"
    "  to the target. Then specify this volume and its "
    "registration with --int.\n"
    "  regmat will be the registration resulting from a separate "
    "invocation of \n"
    "  tkregister2 in which the intermediate volume is specified "
    "as the moveable. \n"
    "  \n"
    "  --gdiagno N\n"
    "\n"
    "  Set the diagnostic/debug level. Default is 0.\n"
    "\n"
    "FREESURFER REGISTRATION CONVENTIONS\n"
    "\n"
    "For the purposes of FreeSurfer, the registration matrix maps "
    "the XYZ\n"
    "of the anatomical reference (ie, the subjects brain as found in \n"
    "$SUBJECTS_DIR) to the XYZ of the functional volume. The anatomical\n"
    "reference is the 'target' volume (argument of --targ) and the \n"
    "functional volume is the 'movable' volume (argument of --mov).\n"
    "The XYZ of a given col, row, and slice is defined\n"
    "based on the field-of-view by the following matrix:\n"
    "\n"
    "          (-dc 0   0  dc*Nc/2)\n"
    "     T =  ( 0  0  ds -ds*Ns/2) \n"
    "          ( 0 -dr  0  dr*Nr/2)\n"
    "          ( 0  0   0     1   )\n"
    "\n"
    "where dc, dr, and ds are the voxel resolutions in the column, row,\n"
    "and slice directions, respectively, and  Nc, Nr, and Ns are the\n"
    "number of columns, rows, and slices.  Under this convention, \n"
    "\n"
    "  XYZ = T*[r c s 1]\n"
    "\n"
    "The FreeSurfer registration matrix is then defined by:\n"
    "\n"
    "   XYZmov = R*XYZtarg\n"
    "\n"
    "FREESURFER REGISTRATION FILE FORMAT\n"
    "\n"
    "The FreeSurfer registration is stored in an ASCII file with the \n"
    "following format:\n"
    "\n"
    "      subjectname\n"
    "      in-plane-res-mm\n"
    "      between-plane-res-mm\n"
    "      intensity\n"
    "      m11 m12 m13 m14\n"
    "      m21 m22 m23 m24\n"
    "      m31 m32 m33 m34\n"
    "      0   0   0   1\n"
    "\n"
    "The subject name is that as found in the SUBJECTS_DIR. "
    "The in-plane-res-mm\n"
    "is the in-plane pixels size in mm. The between-plane-res-mm is\n"
    "the distance between slices in mm. Intensity only affect the "
    "display\n"
    "of the movable in the GUI.\n"
    "\n"
    "USING THE GUI\n"
    "\n"
    "Two volumes are compared slice-by-slice by hitting the "
    "Compare button;\n"
    "this will alternatively display one volume and then the other. "
    "If held\n"
    "down, it will flash. The relative position of the movable "
    "volume can\n"
    "be changed with three sliders (1) Translate, (2) Rotate, "
    "and (3) Scale.\n"
    "Each operates in-plane; the Rotate rotates about the "
    "red cross-hair.\n"
    "The current orientation can be changed by hitting either "
    "the CORONAL,\n"
    "SAGITTAL, or HORIZONTAL buttons. The matrix can be saved "
    "by hitting the\n"
    "SAVE REG button. The intensity of the movable can be "
    "changed by editing\n"
    "the value of the 'fmov:' text box. The brightntess and "
    "contrast of the\n"
    "target can be changed by editing the 'contrast:' and "
    "'midpoint:' text\n"
    "boxes. The target can be masked off to the FOV of the movable by \n"
    "pressing the 'masktarg' radio button. If a surface has "
    "been loaded, it\n"
    "can be toggled by clicking in the display window and hitting 's'.\n"
    "\n"
    "SUMMARY OF KEYPRESS COMMANDS\n"
    "\n"
    "0 swap buffers (same as Compare)\n"
    "1 dispaly target\n"
    "2 dispaly moveable\n"
    "a increase moveable frame by 1\n"
    "b decrease moveable frame by 1\n"
    "c toggle colorization (inorm only)\n"
    "d toggle segmentation visibility\n"
    "e toggle slice prescription indicator\n"
    "i intensity normalize images\n"
    "n use nearest neighbor interpolation\n"
    "t use trilinear interpolation\n"
    "s toggle display of cortical surface\n"
    "x show sagittal view\n"
    "y show horizontal view\n"
    "z show coronal view\n"
    "- or _ zoom out\n"
    "+ or = zoom in\n"
    "p translate up\n"
    ". translate down\n"
    "l translate left\n"
    "; translate right\n"
    "g translate into screen\n"
    "h translate out of screen\n"
    "[ rotate counter-clockwise about image normal\n"
    "] rotate clockwise about image normal\n"
    "q rotate about horiztonal image axis (neg)\n"
    "w rotate about horiztonal image axis (pos)\n"
    "r rotate about vertical image axis (neg)\n"
    "f rotate about vertical image axis (pos)\n"
    "Insert increase scale horizontally\n"
    "Delete decrease scale horizontally\n"
    "Home   increase scale vertically\n"
    "End    decrease scale vertically\n"
    "\n"
    "\n"
    "USING WITH FSL and SPM \n"
    "\n"
    "In order to apply the FreeSurfer tools to data analyzed in "
    "FSL or SPM, it\n"
    "is necessary that you create a FreeSurfer-style registration "
    "matrix. This\n"
    "matrix is stored in a text file which is usually called "
    "'register.dat'.\n"
    "It is not the same as an SPM .mat file. It is not the same "
    "as an FSL .mat \n"
    "file. It is not the same as an analyze .hdr file. You can obtain a\n"
    "FreeSurfer-style registration matrix from those things but "
    "it is not the\n"
    "same as those things. If you do not obtain a FreeSurfer-style "
    "registration \n"
    "matrix, then you will not be able to use the FreeSurfer tools "
    "with your\n"
    "functional data.\n"
    "\n"
    "You can obtain a FreeSurfer-style registration matrix by "
    "directly editing the\n"
    "registration from within tkregister2, however, it will be "
    "easier if you\n"
    "use FSL or SPM to perform the registration and then use "
    "tkregister2 to\n"
    "generate the FreeSurfer-style registration matrix. To do this, "
    "first convert\n"
    "the anatomical (ie, the one in SUBJECTS_DIR/subjectname/mri/orig) \n"
    "to analyze format (use "
    "mri_convert SUBJECTS_DIR/subjectname/mri/orig cor.img).\n"
    "Convert the functional to analyze format (eg, f.img). \n"
    "\n"
    "When using the SPM registration tool, select the functional "
    "(f.img) as the object\n"
    "and the anatomical (cor.img) as the target, and select the "
    "'Coregister only' option. \n"
    "The registration will be written into the functional .mat file. "
    "Run tkregister2 \n"
    "with the --regheader option. Use the --noedit option to suppress\n"
    "the GUI. Editing the registration will not affect the registration\n"
    "as seen by SPM. An example command-line is\n"
    "\n"
    "  tkregister2 --mov f.img --s yoursubject --regheader --noedit "
    "--reg register.dat \n"
    "\n"
    "This will create/overwrite register.dat.\n"
    "\n"
    "For FSL, coregister the functional (f.img) to the anatomical "
    "(cor.img) (ie, \n"
    "use the anatomical (cor.img) as the reference). This will produce "
    "an FSL \nregistration "
    "file (an ASCII file with the 4x4 matrix). Run tkregister2 "
    "specifying\n"
    "this file as the argument to the --fsl flag. Note, it is "
    "possible\n"
    "to generate an FSL matrix from the headers, which can be useful to\n"
    "initialize the FSL registration routines. To do this, just run\n"
    "tkregister2 with the --fslregout fsl.mat and --regheader flags,\n"
    "where fsl.mat is the name of the FSL matrix file. Use the --noedit\n"
    "option to suppress the GUI. Editing the registration will not "
    "affect\n"
    "the registration as seen by FSL unless --fslregout is specfied.\n"
    "An example command-line is\n"
    "\n"
    "  tkregister2 --mov f.img --fsl fsl.mat --s yoursubject "
    "--regheader\n"
    "      --noedit --reg register.dat \n"
    "\n"
    "This will create/overwrite register.dat.\n"
    "\n"
    "If you don't include the --noedit, then the GUI will come up, "
    "and you\n"
    "can flip back and forth between the functional and the structural."
    " If\n"
    "they are not well aligned then you are going to see garbage "
    "when you\n"
    "when you use tkmedit and tksurfer to view your data. So, "
    "CHECK YOUR\n"
    "REGISTRATION!!!\n"
    "\n"
    "BUGS\n"
    "\n"
    "It used to be the case that the GUI would not work if the target\n"
    "was not 256x256x256 voxels, isotropic at 1mm. This has been fixed\n"
    "by resampling the target into this space. This does not affect the\n"
    "final registration, ie, the registration matrix maps the original\n"
    "target (regarless of size) to the movable.\n"
    "\n"
    "AUTHORS\n"
    "\n"
    "The original tkregister was written by Martin Sereno and "
    "Anders Dale\n"
    "in 1996. The original was modified by Douglas Greve in 2002.\n"
    );

  exit(1) ;
}


/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
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
  if (targ_vol_id == NULL) {
    printf("INFO: no target volume specified, assuming "
           "FreeSurfer orig volume.\n");
    targ_vol_id = "orig";
    fstarg =1;
  }
  if (mov_vol_id == NULL && !fstal) {
    printf("ERROR: no movable volume specified\n");
    exit(1);
  }
  if (regfname == NULL && !fstal) {
    printf("ERROR: no registration file specified\n");
    exit(1);
  }
  if (strlen(subjectid) == 0 && fstal) {
    printf("ERROR: must spec subjectid with --fstal\n");
    exit(1);
  }
  if (mkheaderreg && fstal) {
    printf("ERROR: cannot spec both --regheader and --fstal\n");
    exit(1);
  }
  if (xfmfname != NULL && fslregfname != NULL) {
    printf("ERROR: cannot make reg from xfm AND fslreg \n");
    exit(1);
  }
  if (mkheaderreg && fslregfname != NULL) {
    printf("ERROR: cannot make reg from header and fslreg \n");
    exit(1);
  }

  if (mkheaderreg && int_vol_id != NULL) {
    printf("ERROR: cannot make reg from header and use an "
           "intermediate volume\n");
    exit(1);
  }

  if (noedit) LoadVol = 0;

  return;
}


/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"target  volume %s\n",targ_vol_id);
  fprintf(fp,"movable volume %s\n",mov_vol_id);
  fprintf(fp,"reg file       %s\n",regfname);
  fprintf(fp,"LoadVol        %d\n",LoadVol);
  fprintf(fp,"ZeroCRAS       %d\n",ZeroCRAS);

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
static int isflag(const char *flag) {
  int len;
  len = strlen(flag);
  if (len < 2) return(0);

  if (flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}


/*------------------------------------------------------------*/
static int stringmatch(const char *str1, const char *str2) {
  if (! strcmp(str1,str2)) return(1);
  return(0);
}


/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth) {
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if (nargc <= nth) return(0);

  /* check whether the nth arg is a flag */
  if (isflag(argv[nth])) return(0);

  return(1);
}


/*------------------------------------------------------------*/
static int checkfmt(char *fmt) {
  int fmtid;
  fmtid = string_to_type(fmt);
  if (fmtid == MRI_VOLUME_TYPE_UNKNOWN) {
    printf("ERROR: format string %s unrecognized\n",fmt);
    exit(1);
  }
  return(fmtid);
}


#ifdef HAVE_TCL_TK_GL

/*-----------------------------------------------
  imc = current cor - in screen coords
  ic  = current hor - in screen coords
  jc  = current sag - in screen coords
  -----------------------------------------------*/
void draw_image2(int imc,int ic,int jc) {
  extern GLubyte *blinkbuft;
  extern GLubyte *blinkbufm;
  extern int interpmethod;
  extern int plane;
  extern int xdim, ydim, zf; /* window size */
  extern int overlay_mode;
  extern double ps_2, st_2, fscale_2;
  extern int xdim_2, ydim_2, imnr1_2; /* func Nc, Ns, Nr */
  //extern float movimg[WINDOW_ROWS][WINDOW_COLS];
  //extern float targimg[WINDOW_ROWS][WINDOW_COLS];
  //extern int surfimg[WINDOW_ROWS][WINDOW_COLS];
  extern MATRIX *Qmov, *Qtarg;
  static int firstpass = 1;
  static int PrevPlane, PrevImc, PrevIc, PrevJc, PrevOverlayMode;
  static int PrevInorm, PrevInterp, PrevMaskFlag;

  int update_needed;
  int r,c,k;
  char *planestring = NULL;
  unsigned char voxval;
  int ivoxval;
  int cor, ax, sag, toob=0;
  unsigned char* lvidbuf;
  float f=0;
  int NcFunc, NrFunc, NsFunc;
  double rVoxVal;
  double dVoxVal;
  float targimgmax, movimgmax;
  float targimgmin, movimgmin;
  float targimgrange, movimgrange;

  float fcTarg,frTarg,fsTarg;
  float fcMov,frMov,fsMov;
  int icTarg,irTarg,isTarg;
  int icMov,irMov,isMov;
  int valid;

  int segval = 0;

  if (firstpass)    update_needed = 1;
  else {
    update_needed = 0;
    if (PrevPlane   != plane)        update_needed = 1;
    //if(PrevImc     != imc)          update_needed = 2;
    if (PrevInorm   != use_inorm)    update_needed = 3;
    if (PrevInterp  != interpmethod) update_needed = 4;
    if (PrevMaskFlag    != maskflag) update_needed = 5;
    if (PrevOverlayMode != overlay_mode) update_needed = 6;
  }
  if (PixelSelected) update_needed = 0;
  else              update_needed = 1;
  PixelSelected = 0;

  /* Set the current point */
  switch (plane) {
  case SAGITTAL:
    cScreenCur = imc;
    rScreenCur = ic;
    break;
  case HORIZONTAL:
    cScreenCur = jc;
    rScreenCur = imc;
    break;
  case CORONAL:
    cScreenCur = jc;
    rScreenCur = ic;
    break;
  }

  NcFunc = xdim_2;
  NrFunc = ydim_2;
  NsFunc = imnr1_2 + 1;

  if (overlay_mode==TARGET) {
    lvidbuf = blinkbuft;
    pswapnext = TARGET;
  } else {
    lvidbuf = blinkbufm;
    pswapnext = MOVEABLE;
  }

  switch (plane) {
  case CORONAL:
    planestring = "Coronal";
    break;
  case HORIZONTAL:
    planestring = "Horizontal";
    break;
  case SAGITTAL:
    planestring = "Sagittal";
    break;
  }

  if (Gdiag_no > 0) {
    printf("%d ----------------------------------------\n",update_needed);
    crScreen2AnatInd(jc, ic, &cor, &ax, &sag);
    printf("imc = %d, ic = %d, jc = %d\n",imc,ic,jc);
    printf("cor = %d, ax = %d, sag = %d, toob = %d\n",cor,ax,sag,toob);
    printf("Func: Nc = %d, Nr = %d, Ns = %d\n",NcFunc,NrFunc,NsFunc);
    printf("Anat: Nc = %d, Nr = %d, Ns = %d\n",xnum_2,ynum_2,numimg_2);
    printf("ps_2 = %g, st_2 = %g, fscale_2 = %g\n",ps_2,st_2,fscale_2);
    printf("maskflag = %d, inorm = %d, interp=%d\n",maskflag,use_inorm,
           interpmethod);
    printf("xdim = %d, ydim = %d, zf = %d\n",xdim,ydim,zf);
    printf("mode = %d, plane = %s (%d), %d %d %d\n",
           overlay_mode, planestring, plane, imc, ic, jc);
    fflush(stdout);
  }

  UpdateMatrices();
  if (Gdiag_no > 0) {
    printf("-------- Tscreen ----------\n");
    MatrixPrint(stdout,Tscreen);
    printf("-------- invTtarg ----------\n");
    MatrixPrint(stdout,invTtarg);
    printf("-------- invTmov ----------\n");
    MatrixPrint(stdout,invTmov);
    printf("-------- Qtarg ----------\n");
    MatrixPrint(stdout,Qtarg);
    printf("-------- Qmov ----------\n");
    MatrixPrint(stdout,Qmov);
  }

  if (update_needed != 0) {
    targimgmax = -100000;
    targimgmin = 1000000000.0;
    movimgmax  = -100000;
    movimgmin  = 1000000000.0;
    for (r = 0; r < ydim; r++) {
      for (c = 0; c < xdim; c++) {
        targimg[r][c] = 0.0;
        movimg[r][c]  = 0.0;
        surfimg[r][c] = 0;

        /* crsTarg = Qtarg * [c r 0 1]' */
        fcTarg = Qtarg->rptr[1][1]*c + Qtarg->rptr[1][2]*r + Qtarg->rptr[1][4];
        frTarg = Qtarg->rptr[2][1]*c + Qtarg->rptr[2][2]*r + Qtarg->rptr[2][4];
        fsTarg = Qtarg->rptr[3][1]*c + Qtarg->rptr[3][2]*r + Qtarg->rptr[3][4];

        /* use (int) here to make the halfs consistent */
        icTarg = (int)fcTarg;
        irTarg = (int)frTarg;
        isTarg = (int)fsTarg;
        //icTarg = nint(fcTarg);
        //irTarg = nint(frTarg);
        //isTarg = nint(fsTarg);

        if (icTarg < 0 || icTarg >= targ_vol->width ||
            irTarg < 0 || irTarg >= targ_vol->height ||
            isTarg < 0 || isTarg >= targ_vol->depth)  continue;

	if(! ShowSeg ){
	  /* Could interp here, but why bother? */
	  //targimg[r][c] = MRIFseq_vox(targ_vol,icTarg,irTarg,isTarg,0);
	  /* This implements the trilinear interp - makes it so that
	     when the same volume is loaded as targ and mov, the
	     registration is perfect */
	  //MRIsampleVolume(targ_vol,fcTarg,frTarg,fsTarg,&rVoxVal);
	  MRIsampleVolumeFrame(targ_vol,fcTarg,frTarg,fsTarg,0,&dVoxVal);
	  rVoxVal = dVoxVal;
	  targimg[r][c] = rVoxVal;
	}
	else {
	  segval = MRIgetVoxVal(seg_vol,icTarg,irTarg,isTarg,0);
	  targimg[r][c] = segval;
	}

        if (UseSurf) surfimg[r][c] = MRIvox(mrisurf,icTarg,irTarg,isTarg);

        /* Keep track of the max and min */
        if (targimgmax < targimg[r][c]) targimgmax = targimg[r][c];
        if (targimgmin > targimg[r][c]) targimgmin = targimg[r][c];

        /* crsMov = Qmov * [c r 0 1]' */
        fcMov = Qmov->rptr[1][1]*c + Qmov->rptr[1][2]*r + Qmov->rptr[1][4];
        frMov = Qmov->rptr[2][1]*c + Qmov->rptr[2][2]*r + Qmov->rptr[2][4];
        fsMov = Qmov->rptr[3][1]*c + Qmov->rptr[3][2]*r + Qmov->rptr[3][4];

        icMov = nint(fcMov);
        irMov = nint(frMov);
        isMov = nint(fsMov);
        if (icMov < 0 || icMov >= mov_vol->width ||
            irMov < 0 || irMov >= mov_vol->height ||
            isMov < 0 || isMov >= mov_vol->depth ) {
          if (maskflag) targimg[r][c] = 0;
          continue;
        }

        if (Gdiag_no > 3) {
          printf("Scrn cr:  %d, %d\n",c,r);
          printf("Anat crs: %g, %g, %g\n",fcTarg,frTarg,fsTarg);
          printf("Anat crs: %d, %d, %d\n",icTarg,irTarg,isTarg);
          printf("Func crs: %g, %g, %g\n",fcMov,frMov,fsMov);
          printf("Func crs: %d, %d, %d\n",icMov,irMov,isMov);
          fflush(stdout);
        }

        switch (interpmethod) {
        case  SAMPLE_NEAREST:
          switch (float2int_use) {
          case FLT2INT_ROUND:
            f = (float) MRIFseq_vox(mov_vol,icMov,irMov,isMov,mov_frame);
            break;
          case FLT2INT_TKREG:
            /* This is provided for testing only */
            icMov = floor(fcMov);
            irMov =  ceil(frMov);
            isMov = floor(fsMov);
            if (icMov < 0 || icMov >= mov_vol->width ||
                irMov < 0 || irMov >= mov_vol->height ||
                isMov < 0 || isMov >= mov_vol->depth)
              f = 0;
            else
              f = (float) MRIFseq_vox(mov_vol,icMov,irMov,isMov,mov_frame);
            break;
          }
          //printf("Func crs: %d, %d, %d \n",cFunc,rFunc,sFunc);
          //f = (float) fscale_2*fim_2[sFunc][rFunc][cFunc];
          break;
        case  SAMPLE_TRILINEAR:
          MRIsampleVolumeFrame(mov_vol,fcMov,frMov,fsMov,mov_frame,&rVoxVal);
          f = rVoxVal;
          break;
        case  SAMPLE_SINC:
          MRIsincSampleVolume(mov_vol,fcMov,frMov,fsMov,5,&rVoxVal);
          f = rVoxVal;
          break;
        }
        movimg[r][c] = f;
        if(DoSlicePrescription){
          if(icMov == 0 || irMov == 0   || isMov == 0 ||
             icMov == mov_vol->width-1  || irMov == mov_vol->height-1 ||
             isMov == mov_vol->depth-1  || isMov%4 == 0) {
            targimg[r][c] = 255;
            movimg[r][c] = 255;
          }
        }

        if (movimgmax  < movimg[r][c])  movimgmax  = movimg[r][c];
        if (movimgmin  > movimg[r][c])  movimgmin  = movimg[r][c];

      } // image column
    } // image row

    /* Auto-control of intensity - this does not work very well*/
    targimgrange = targimgmax - targimgmin;
    if (targimgrange == 0) targimgrange = 1;
    movimgrange = movimgmax - movimgmin;
    if (movimgrange == 0) movimgrange = 1;
    for (r = 0; r < ydim; r++) {
      for (c = 0; c < xdim; c++) {
        if(use_inorm) {
	  if(! ShowSeg){
	    targimg[r][c] = 255*(targimg[r][c] - targimgmin)/targimgrange;
	    if (targimg[r][c] < 0) targimg[r][c] = 0;
	  }
          movimg[r][c]  = 255*(movimg[r][c]  - movimgmin)/movimgrange;
          if (movimg[r][c] < 0) movimg[r][c] = 0;
        } 
	else {
	  movimg[r][c] *= fscale_2;
	  if(DoFMovTarg) targimg[r][c] *= fscale_2;
	}
      }
    }

    if (Gdiag_no > 0) {
      printf("targimg: %g %g %g\n",targimgmin,targimgmax,targimgrange);
      printf("movimg:  %g %g %g\n",movimgmin,movimgmax,movimgrange);
    }
  }/*--------- end if update needed ---------------------*/

  k = 0;
  for (r = 0; r < ydim; r++) {
    for (c = 0; c < xdim; c++) {
      if (UseSurf && surfimg[r][c]) {
        lvidbuf[3*k]   = SurfRGB[0];
        lvidbuf[3*k+1] = SurfRGB[1];
        lvidbuf[3*k+2] = SurfRGB[2];
      } else {
        switch (overlay_mode) {
        case TARGET:
          f = targimg[r][c];
          break;
        case MOVEABLE:
          f = movimg[r][c];
          break;
        }
        if (f < 0)        f = 0;
        else if (f > 255) f = 255;
        voxval = (unsigned char)(nint(f));
        if(!use_inorm) {
          lvidbuf[3*k]   = colormap[voxval];
          lvidbuf[3*k+1] = colormap[voxval];
          lvidbuf[3*k+2] = colormap[voxval];
	  if(use_colornorm){
	    if(voxval > 128){
	      lvidbuf[3*k]   = 0;
	      lvidbuf[3*k+1] = colormap[voxval];
	      lvidbuf[3*k+2] = 0;
	    } else {
	      lvidbuf[3*k]   = 0;
	      lvidbuf[3*k+1] = 0;
	      lvidbuf[3*k+2] = colormap[voxval];
	    }
	  }
        } else {
          lvidbuf[3*k]   = voxval;
          lvidbuf[3*k+1] = voxval;
          lvidbuf[3*k+2] = voxval;
        }
	if(overlay_mode == TARGET && ShowSeg) {
	  ivoxval = targimg[r][c];
          CTABisEntryValid(ctab,ivoxval,&valid);
          if(!valid) {
            printf("ERROR: cannot find seg id %d (%d,%d)\n",ivoxval,r,c);
            exit(1);
          }
          lvidbuf[3*k]   = ctab->entries[ivoxval]->ri;
          lvidbuf[3*k+1] = ctab->entries[ivoxval]->gi;
          lvidbuf[3*k+2] = ctab->entries[ivoxval]->bi;
	}
      }
      k++;
    }
  }

  /* draw the cross hair */
  draw_cross_hair(rScreenCur,cScreenCur);

  swapbuffers();
  invalid_buffers=1;

  PrevPlane = plane;
  PrevImc = imc;
  PrevIc = ic;
  PrevJc = jc;
  PrevOverlayMode = overlay_mode;
  PrevInorm = use_inorm;
  PrevInterp = interpmethod;
  PrevMaskFlag = maskflag;
  firstpass = 0;
}

#endif // HAVE_TCL_TK_GL


/*------------------------------------------------------------------
  ScreenCR2XYZMtx() - creates the matrix that converts screen
  column and row to to Anatomical XYZ. CR origin at bottom left.
  Anatomical XYZ is the same as the registration space.
  ----------------------------------------------------------*/
MATRIX *ScreenCR2XYZMtx(MATRIX *T) {
  extern int plane; /* view port orientation */
  extern int imc, ic, jc; /* current cor, hor, and sag, sort of */
  extern int xdim,ydim,zf; /* Screen dim and zoom factor */
  extern int xnum,ynum,numimg; /* Anat: Nr, Nc, Ns */
  static int first = 1;
  extern float FOV;
  static MATRIX *Pxyz, *dcCol, *dcRow, *crsCur, *xyzCur;
  float deltaCol, deltaRow;
  int NcScreen, NrScreen;
  int Nsag, Nhor, Ncor;
  int sagCur, horCur, corCur,i;

  if (first) {
    first = 0;
    Pxyz = MatrixAlloc(3,1,MATRIX_REAL);
    dcCol = MatrixAlloc(3,1,MATRIX_REAL);
    dcRow = MatrixAlloc(3,1,MATRIX_REAL);
    crsCur = MatrixAlloc(4,1,MATRIX_REAL);
    crsCur->rptr[4][1] = 1;
    xyzCur = MatrixAlloc(4,1,MATRIX_REAL);
  }

  if (T==NULL) T = MatrixAlloc(4,4,MATRIX_REAL);

  NcScreen = xdim;
  NrScreen = ydim;
  Nsag = xnum;
  Nhor = ynum;
  Ncor = numimg;

  corCur = (int)((float)imc/zf);
  horCur = (int)((float)(NrScreen-ic)/zf);
  sagCur = (int)((float)jc/zf);

  crsCur->rptr[1][1] = sagCur;
  crsCur->rptr[2][1] = horCur;
  crsCur->rptr[3][1] = corCur;

  xyzCur = MatrixMultiply(Ttarg,crsCur,xyzCur);

  deltaCol = FOV/NcScreen;
  deltaRow = FOV/NrScreen;

  switch (plane) {
  case CORONAL:
    dcCol->rptr[1][1] = -1.0;
    dcCol->rptr[2][1] =  0.0;
    dcCol->rptr[3][1] =  0.0;
    dcRow->rptr[1][1] =  0.0;
    dcRow->rptr[2][1] =  0.0;
    dcRow->rptr[3][1] = +1.0;
    //Pxyz->rptr[1][1] = +128.0; /*x*/
    Pxyz->rptr[1][1] = +FOV/2; /*x*/
    Pxyz->rptr[2][1] = xyzCur->rptr[2][1]; /*y*/
    //Pxyz->rptr[3][1] = -127.0; /*z*/
    Pxyz->rptr[3][1] = -(FOV-2)/2; /*z*/
    break;
  case SAGITTAL:
    dcCol->rptr[1][1] =  0.0;
    dcCol->rptr[2][1] = +1.0;
    dcCol->rptr[3][1] =  0.0;
    dcRow->rptr[1][1] =  0.0;
    dcRow->rptr[2][1] =  0.0;
    dcRow->rptr[3][1] = +1.0;
    Pxyz->rptr[1][1] = xyzCur->rptr[1][1]; /*x*/
    //Pxyz->rptr[2][1] = -128.0; /*y*/
    Pxyz->rptr[2][1] = -FOV/2; /*y*/
    //Pxyz->rptr[3][1] = -127.0; /*z*/
    Pxyz->rptr[3][1] = -(FOV-2)/2; /*z*/
    break;
  case HORIZONTAL:
    dcCol->rptr[1][1] = -1.0;
    dcCol->rptr[2][1] =  0.0;
    dcCol->rptr[3][1] =  0.0;
    dcRow->rptr[1][1] =  0.0;
    dcRow->rptr[2][1] = +1.0;
    dcRow->rptr[3][1] =  0.0;
    //Pxyz->rptr[1][1] = +128.0;
    Pxyz->rptr[1][1] = +FOV/2;
    //Pxyz->rptr[2][1] = -127.0;
    Pxyz->rptr[2][1] = -(FOV-2)/2;
    Pxyz->rptr[3][1] = xyzCur->rptr[3][1]; /*z*/
    break;
  }

  for (i=1;i<=3;i++) {
    T->rptr[i][1] = (dcCol->rptr[i][1])*deltaCol;
    T->rptr[i][2] = (dcRow->rptr[i][1])*deltaRow;
    T->rptr[i][3] = 0;
    T->rptr[i][4] = Pxyz->rptr[i][1];
    T->rptr[4][i] = 0.0;
  }
  T->rptr[4][4] = 1;

  return(T);
}


/*----------------------------------------------------------
  UpdateMatrices() - updates relevant transform matrices. They
  may need to be updated in response to changes in the view
  (ie, changes to Tscreen) or changes in the registration
  (ie, RegMat or tm).
  ----------------------------------------------------------*/
void UpdateMatrices(void) {
  extern MATRIX *Tscreen, *RegMat, *invTtarg, *invTmov;
  extern MATRIX *Qtarg,*Qmov;
  extern float tm[4][4];
  int r,c;

  if (RegMat==NULL) RegMat = MatrixAlloc(4,4,MATRIX_REAL);
  for (r=0;r<4;r++) {
    for (c=0;c<4;c++) {
      RegMat->rptr[r+1][c+1] = tm[r][c];
    }
  }

  /* Tscreen converts from Screen CR to Targ XYZ */
  Tscreen = ScreenCR2XYZMtx(Tscreen);

  /* Qtarg converts from Screen CR to Targ CRS */
  /* Qtarg = inv(Ttarg)*Tscreen */
  Qtarg = MatrixMultiply(invTtarg,Tscreen,Qtarg);

  /* Qmov converts from Screen CR to Mov CRS */
  /* Qmov = inv(Tmov)*R*Tscreen */
  Qmov  = MatrixMultiply(invTmov,RegMat,Qmov);
  Qmov  = MatrixMultiply(Qmov,Tscreen,Qmov);
}


/*------------------------------------------------------------------
  int crScreenInd2AnatInd - CR origin at bottom left.
  col -> Sag, row -> hor, slice -> Cor. Returns 1 if out-of-bounds,
  otherwise returns 0.
  ----------------------------------------------------------*/
int crScreen2AnatInd(int c, int r, int *cor, int *hor, int *sag) {
  extern MATRIX *Qtarg;
  extern int xnum,ynum,numimg; /* Anat: Nr, Nc, Ns */
  int Nsag, Nhor, Ncor;
  float fcTarg, frTarg, fsTarg;

  UpdateMatrices();

  fcTarg = Qtarg->rptr[1][1]*c + Qtarg->rptr[1][2]*r + Qtarg->rptr[1][4];
  frTarg = Qtarg->rptr[2][1]*c + Qtarg->rptr[2][2]*r + Qtarg->rptr[2][4];
  fsTarg = Qtarg->rptr[3][1]*c + Qtarg->rptr[3][2]*r + Qtarg->rptr[3][4];

  *cor = (int)fsTarg;
  *hor = (int)frTarg;
  *sag = (int)fcTarg;

  Nsag = xnum;
  Nhor = ynum;
  Ncor = numimg;

  if (*sag < 0)     return(1);
  if (*sag >= Nsag) return(1);
  if (*cor < 0)     return(1);
  if (*cor >= Ncor) return(1);
  if (*hor < 0)     return(1);
  if (*hor >= Nhor) return(1);

  return(0);
}


/*------------------------------------------------------------------
  int crScreenInd2AnatXYZ - CR origin at bottom left. Better than
  AnatCRS2AnatXYZ because screen has more resolution.
  ----------------------------------------------------------*/
int crScreen2AnatXYZ(int c, int r, float *x, float *y, float *z) {
  extern MATRIX *Tscreen;

  UpdateMatrices();
  *x = Tscreen->rptr[1][1]*c + Tscreen->rptr[1][2]*r + Tscreen->rptr[1][4];
  *y = Tscreen->rptr[2][1]*c + Tscreen->rptr[2][2]*r + Tscreen->rptr[2][4];
  *z = Tscreen->rptr[3][1]*c + Tscreen->rptr[3][2]*r + Tscreen->rptr[3][4];

  return(0); /*------------------------------------------*/
}


/*-----------------------------------------------*/
int AnatInd2AnatXYZ(int cor, int hor, int sag,
                    float *x, float *y, float *z) {
  extern MATRIX *Ttarg;
  static MATRIX *crs, *xyz;

  if (crs==NULL) {
    crs = MatrixAlloc(4,1,MATRIX_REAL);
    crs->rptr[4][1] = 1;
  }

  crs->rptr[1][1] = sag;
  crs->rptr[2][1] = hor;
  crs->rptr[3][1] = cor;

  //UpdateMatrices(); /* no need to update here */
  xyz = MatrixMultiply(Ttarg,crs,xyz);
  *x = xyz->rptr[1][1];
  *y = xyz->rptr[2][1];
  *z = xyz->rptr[3][1];

  return(0);
}


/*-----------------------------------------------*/
int FuncXYZ2FuncInd(float x, float y, float z,
                    float *col, float *row, float *slice) {
  extern MATRIX *invTmov;
  extern int xdim_2, ydim_2, imnr1_2; /* func Nc, Ns, Nr */
  int Nc, Nr, Ns;
  static MATRIX *xyz, *crs;

  if (xyz == NULL) {
    xyz = MatrixAlloc(4,1,MATRIX_REAL);
    xyz->rptr[4][1] = 1.0;
  }

  //UpdateMatrices(); /* no need to update here */

  xyz->rptr[1][1] = x;
  xyz->rptr[2][1] = y;
  xyz->rptr[3][1] = z;
  crs = MatrixMultiply(invTmov,xyz,crs);
  *col   = crs->rptr[1][1];
  *row   = crs->rptr[2][1];
  *slice = crs->rptr[3][1];

  Nc = xdim_2;
  Nr = ydim_2;
  Ns = imnr1_2 + 1;

  if (nint(*col)   <    0) return(1);
  if (nint(*col)   >=  Nc) return(1);
  if (nint(*row)   <    0) return(1);
  if (nint(*row)   >=  Nr) return(1);
  if (nint(*slice) <    0) return(1);
  if (nint(*slice) >=  Ns) return(1);

  return(0);
}


#ifdef HAVE_TCL_TK_GL

/*-----------------------------------------------------------*/
int draw_cross_hair(int rScreen, int cScreen) {
  int r,c,k;
  unsigned char *lvidbuf;

  //printf("Drawing crosshair at row = %d, col = %d\n",rScreen,cScreen);

  if(overlay_mode==TARGET) lvidbuf = blinkbuft;
  else                     lvidbuf = blinkbufm;

  if (lvidbuf == NULL) return(1);

  for (r = rScreen - 4; r < rScreen + 4; r++) {
    if (r < 0 || r >= WINDOW_ROWS) continue;
    k =  ((r*WINDOW_COLS)+cScreen);
    lvidbuf[3*k]   = 255;
    lvidbuf[3*k+1] = 0;
    lvidbuf[3*k+2] = 0;
    //printf("r = %d, c = %d, k = %d\n",r,cScreen,k);
  }
  for (c = cScreen - 4; c < cScreen + 4; c++) {
    if (c < 0 || c >= WINDOW_COLS) continue;
    k = ((rScreen*WINDOW_COLS)+c);
    lvidbuf[3*k]   = 255;
    lvidbuf[3*k+1] = 0;
    lvidbuf[3*k+2] = 0;
    //printf("r = %d, c = %d, k = %d\n",rScreen,c,k);
  }

  return(0);
}


/*-----------------------------------------------------------*/
int erase_cross_hair(int rScreen, int cScreen) {
  int r,c,k;
  unsigned char *lvidbuf, voxval;
  float f=0;

  if (overlay_mode==TARGET) lvidbuf = blinkbuft;
  else                     lvidbuf = blinkbufm;

  if (lvidbuf == NULL) return(1);

  for (r = rScreen - 4; r < rScreen + 4; r++) {
    if (r < 0 || r >= WINDOW_ROWS) continue;
    k = 3 * ((r*WINDOW_COLS)+cScreen);
    switch (overlay_mode) {
    case TARGET:
      f = targimg[r][cScreen];
      break;
    case MOVEABLE:
      f = movimg[r][cScreen];
      break;
    }
    if (f < 0)        f = 0;
    else if (f > 255) f = 255;
    voxval = (unsigned char)f;
    lvidbuf[k]   = colormap[voxval];
    lvidbuf[k+1] = colormap[voxval];
    lvidbuf[k+2] = colormap[voxval];
  }
  for (c = cScreen - 4; c < cScreen + 4; c++) {
    if (c < 0 || c >= WINDOW_COLS) continue;
    k = 3 * ((rScreen*WINDOW_COLS)+c);
    switch (overlay_mode) {
    case TARGET:
      f = targimg[rScreen][c];
      break;
    case MOVEABLE:
      f = movimg[rScreen][c];
      break;
    }
    if (f < 0)        f = 0;
    else if (f > 255) f = 255;
    voxval = (unsigned char)f;
    lvidbuf[k]   = colormap[voxval];
    lvidbuf[k+1] = colormap[voxval];
    lvidbuf[k+2] = colormap[voxval];
  }

  //printf("Erase: done\n");

  return(0);
}


/*-----------------------------------------------------------*/
/* sx,sy are with respect to the desktop */
void select_pixel(short sx, short sy) {
  extern int cScreenCur, rScreenCur;
  long ox,oy,lx,ly;
  int cor, hor, sag;
  float xAnat, yAnat, zAnat;
  float xFunc, yFunc, zFunc;
  float cfFunc, rfFunc, sfFunc;
  int cFunc, rFunc, sFunc;
  int soob, foob;
  int kScreen, cScreen, rScreen;
  unsigned char *lvidbuf;

  getorigin(&ox,&oy);
  getsize(&lx,&ly);
  //printf("select_pix: sx = %d, sy = %d, ox = %d, oy = %d, lx=%d, ly = %d\n",
  // sx,sy,ox,oy,lx,ly);

  if (overlay_mode==TARGET) lvidbuf = blinkbuft;
  else                     lvidbuf = blinkbufm;
  //printf("sx=%d, sy=%d, ox=%ld, oy=%ld, lx=%ld, ly=%ld, dx = %ld, dy = %ld\n",
  // sx,sy,ox,oy,lx,ly,sx-ox,sy-oy);

  /* hack: x-y of click on tk win caught by qread when clicks buffered! */
  if (sx<ox || sx>ox+lx) sx = ox+lx/2;
  if (sy<oy || sy>oy+ly) sy = oy+ly/2;

  if (plane==CORONAL) {
    ic = (sy-oy);
    jc = (sx-ox);
  } else
    if (plane==HORIZONTAL) {
      imc = (sy-oy);
      jc = (sx-ox);
    } else
      if (plane==SAGITTAL) {
        imc = (sx-ox);
        ic = (sy-oy);
      }
  xc = xx1-ps*jc/fsf;
  yc = yy0+st*imc/fsf;
  //zc = zz1-ps*(255.0-ic/fsf);
  zc = zz1-ps*(((WINDOW_ROWS/2.0)-1)-ic/fsf);
  //printf("select_pixel: xc=%g, yc=%g,  zc=%g\n",xc,yc,zc);
  //printf("select_pixel: jc=%d, imc=%d, ic=%d\n",jc,imc,ic);
  //printf("xx1 = %g, ps = %g, fsf = %g, st=%g, yy0 = %g, zz1 = %g\n",
  //     xx1,ps,fsf,st,yy0,zz1);

  erase_cross_hair(rScreenCur,cScreenCur);
  cScreen = sx-ox;
  rScreen = sy-oy;
  kScreen = (rScreen * WINDOW_COLS) + cScreen;
  cScreenCur = cScreen;
  rScreenCur = rScreen;

  soob = crScreen2AnatInd(cScreen, rScreen, &cor, &hor, &sag);
  AnatInd2AnatXYZ(cor, hor, sag, &xAnat, &yAnat, &zAnat);
  transform(xAnat,yAnat,zAnat,&xFunc,&yFunc,&zFunc,tm);
  foob = FuncXYZ2FuncInd(xFunc,yFunc,zFunc,&cfFunc,&rfFunc,&sfFunc);
  cFunc = nint(cfFunc);
  rFunc = nint(rfFunc);
  sFunc = nint(sfFunc);

  printf("------------------------------------------------------------\n");
  printf("  Screen:  %3d %3d (%d,%d,%d), inorm = %d, mov_frame = %d \n",
         cScreen,rScreen,lvidbuf[3*kScreen],
         lvidbuf[3*kScreen+1], lvidbuf[3*kScreen+2],
         use_inorm,mov_frame);
  printf("  Anat:    (%3d %3d %3d)   (%6.1f %6.1f %6.1f)  ",
         sag,hor,cor,xAnat,yAnat,zAnat);
  if (!soob) printf("%8.4f  %6.1f\n",MRIFseq_vox(targ_vol,sag,hor,cor,0),
                    targimg[rScreen][cScreen]);
  else      printf("OutOfBounds\n");

  printf("  Func:    (%3d %3d %3d)   (%6.1f %6.1f %6.1f)  ",
         cFunc,rFunc,sFunc,xFunc,yFunc,zFunc);
  if (!foob) printf("%8.4f  %6.1f\n",MRIFseq_vox(mov_vol,cFunc,rFunc,sFunc,0),
                    movimg[rScreen][cScreen]);
  else      printf("OutOfBounds\n");
  printf("--------------------------------------------------------------\n");

  draw_cross_hair(rScreenCur,cScreenCur);
  invalid_buffers=1;

}


/*-----------------------------------------------------------*/
int do_one_gl_event(Tcl_Interp *interp)   /* tcl */
{
  XEvent current, ahead;
  char buf[1000], c=0;
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
  float d,r;

  blinkbuffers();

  if (updateflag) {
    set_scale();
    redraw();
    updateflag = FALSE;
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

      resize_window_intstep();
      Tcl_Eval(interp,"set zf $zf");  /* touch for trace */
      if (followglwinflag && zf != 4) {
        /* below */
        /*sprintf(command,"wm geometry . +%d+%d",w.x,w.y+w.h+MOTIF_YFUDGE);*/
        /* above */
        sprintf(command,"putaboveglwin %d %d", w.x, w.y-MOTIF_YFUDGE);
        Tcl_Eval(interp,command);
        /* Tcl_Eval(interp,"raise ."); */
      }
      if (visible_mode!=overlay_mode) /* suppress mode change */
        overlay_mode = visible_mode;
      glViewport(0, 0, xdim, ydim);
      /* updateflag = TRUE; */
      break;

    case Expose:
      if (XPending(xDisplay)) {
        XPeekEvent(xDisplay, &ahead);
        if (ahead.type==Expose) break;  /* skip extras */
      }
      if (visible_mode!=overlay_mode)
        overlay_mode = visible_mode;
      updateflag = TRUE;
      break;

    case ButtonPress:
      sx = current.xbutton.x;
      sy = current.xbutton.y;
      sx += w.x;   /* convert back to screen pos (ugh) */
      sy = 1024 - w.y - sy;
      if (current.xbutton.button == 1) {  /** left **/
        select_pixel(sx,sy);
        Tcl_Eval(interp,"unzoomcoords $plane");
        button1pressed = TRUE;
        PixelSelected = 1;
      }
      if (current.xbutton.button == 2) {  /** middle **/
        button2pressed = TRUE;
      }
      if (current.xbutton.button == 3) {  /** right **/
        button3pressed = TRUE;
      }
      if (visible_mode!=overlay_mode)
        overlay_mode = visible_mode;
      updateflag = TRUE;
      break;

    case ButtonRelease:
      if (current.xbutton.button == 1)  button1pressed = FALSE;
      if (current.xbutton.button == 2)  button2pressed = FALSE;
      if (current.xbutton.button == 3)  button3pressed = FALSE;
      break;

    case KeyPress:
      XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
      //printf("button press %c\n",buf[0]);
      switch (ks) {

        /* numbers */
      case XK_0:
        record_swapbuffers();
        updateflag = TRUE;
        break;
      case 'd':
        if(altkeypressed) {
          if (overlay_mode == TARGET) overlay_mode = MOVEABLE;
          else if (overlay_mode == MOVEABLE) overlay_mode = TARGET;
          updateflag = TRUE;
        }
	else{
	  if(seg_vol != NULL){
	    ShowSeg = ~ShowSeg;
	    updateflag = TRUE;
	  }
	}
        break;

      case XK_1:
        overlay_mode = TARGET;
        updateflag = TRUE;
        break;
      case XK_2:
        overlay_mode = MOVEABLE;
        updateflag = TRUE;
        break;

        /* upper case */
      case XK_A:
        if (altkeypressed)  ;
        break;

        /* lower case */
      case XK_c:
	use_colornorm = !use_colornorm;
	SurfRGB[0] = 255;
	SurfRGB[1] = 0;
	SurfRGB[2] = 0;
        updateflag = TRUE;
        break;
      case XK_i:
        use_inorm = !use_inorm;
        updateflag = TRUE;
        break;
      case XK_n:
        interpmethod = SAMPLE_NEAREST;
        updateflag = TRUE;
        break;
      case XK_t:
        interpmethod = SAMPLE_TRILINEAR;
        updateflag = TRUE;
        break;

      case 'e':
        DoSlicePrescription = !DoSlicePrescription;
        updateflag = TRUE;
        break;

        /* translate in-plane up or down */
      case XK_p:
      case '.':
        if (plane == SAGITTAL)   c = 'z';
        if (plane == CORONAL)    c = 'z';
        if (plane == HORIZONTAL) c = 'y';
        if (ks == 'p') translate_brain(-.5,c);
        if (ks == '.') translate_brain(+.5,c);
        updateflag = TRUE;
        break;

        /* translate left or right */
      case XK_l:
      case ';':
        d = 0.5;
        if (plane == SAGITTAL)   {
          c = 'y';
          d = -0.5;
        }
        if (plane == CORONAL)    c = 'x';
        if (plane == HORIZONTAL) c = 'x';
        if (ks == ';') translate_brain(+d,c);
        if (ks == 'l') translate_brain(-d,c);
        updateflag = TRUE;
        break;

        /* translate thru-plane in or out */
      case 'g':
      case 'h':
        if (plane == SAGITTAL)   c = 'x';
        if (plane == CORONAL)    c = 'y';
        if (plane == HORIZONTAL) c = 'z';
        if (ks == 'g') translate_brain(-.5,c);
        if (ks == 'h') translate_brain(+.5,c);
        updateflag = TRUE;
        break;

        /* rotation about normal to image plane */
      case '[':
      case ']':
        r = +2.0;
        if (plane == SAGITTAL)   c = 'x';
        if (plane == CORONAL)    c = 'y';
        if (plane == HORIZONTAL) {
          c = 'z';
          r = -2.0;
        }
        if (ks == ']') rotate_brain(+r,c); /* in tenths of deg */
        if (ks == '[') rotate_brain(-r,c); /* in tenths of deg */
        updateflag = TRUE;
        break;

        /* rotation about horizontal image axis */
      case 'q': /* top comes out, bottom goes in */
      case 'w': /* top goes in, comes out */
        //printf("ks = %c\n",ks);
        r = +10.0;
        if (plane == SAGITTAL)   c = 'y';
        if (plane == CORONAL)    c = 'x';
        if (plane == HORIZONTAL) {
          c = 'x';
          r = -10.0;
        }
        if (ks == 'q') rotate_brain(+r,c); /* in tenths of deg */
        if (ks == 'w') rotate_brain(-r,c); /* in tenths of deg */
        updateflag = TRUE;
        break;

        /* rotation about vertical image axis */
      case 'r': /* left goes in, right comes out*/
      case 'f':
        //printf("ks = %c\n",ks);
        r = +10.0;
        if (plane == SAGITTAL)   c = 'z';
        if (plane == CORONAL)    c = 'z';
        if (plane == HORIZONTAL) {
          c = 'y';
          r = -10.0;
        }
        if (ks == 'r') rotate_brain(+r,c); /* in tenths of deg */
        if (ks == 'f') rotate_brain(-r,c); /* in tenths of deg */
        updateflag = TRUE;
        break;

      case 'a': /* advance frame */
        mov_frame ++;
        if (mov_frame >= mov_vol->nframes) mov_frame = 0;
        //printf("mov_frame = %d\n",mov_frame);
        updateflag = TRUE;
        break;

      case 'b': /* recede frame */
        mov_frame --;
        if (mov_frame < 0) mov_frame = mov_vol->nframes-1;
        //printf("mov_frame = %d\n",mov_frame);
        updateflag = TRUE;
        break;

        /* scale horizontally */
      case XK_Insert:
      case XK_Delete:
        if (plane == SAGITTAL)   c = 'y';
        if (plane == CORONAL)    c = 'x';
        if (plane == HORIZONTAL) c = 'x';
        if (ks == XK_Insert) scale_brain(+0.995,c);
        if (ks == XK_Delete) scale_brain(+1.005,c);
        updateflag = TRUE;
        break;

        /* scale vertically */
      case XK_Home:
      case XK_End:
        if (plane == SAGITTAL)   c = 'z';
        if (plane == CORONAL)    c = 'z';
        if (plane == HORIZONTAL) c = 'y';
        if (ks == XK_Home) scale_brain(+0.995,c);
        if (ks == XK_End) scale_brain(+1.005,c);
        updateflag = TRUE;
        break;

      case XK_x:
        plane=SAGITTAL;
        updateflag = TRUE;
        break;
      case XK_y:
        plane=CORONAL;
        updateflag = TRUE;
        break;
      case XK_z:
        plane=HORIZONTAL;
        updateflag = TRUE;
        break;
      case XK_s:
        if (LoadSurf) UseSurf = !UseSurf;
        updateflag = TRUE;
        break;

      case '+':
      case '=':
        FOV = FOV/2;
        updateflag = TRUE;
        break;
      case '-':
      case '_':
        FOV = 2*FOV;
        updateflag = TRUE;
        break;


        /* others */
      case XK_Up:
        upslice();
        updateflag = TRUE;
        //Tcl_Eval(interp,
        //         "set fscale_2 [expr $fscale_2 * 1.5]; set updateflag TRUE");
        break;
      case XK_Down:
        downslice();
        updateflag = TRUE;
        //Tcl_Eval(interp,
        //       "set fscale_2 [expr $fscale_2 / 1.5]; set updateflag TRUE");
        break;
      case XK_Right:
        Tcl_Eval(interp,"changeslice up $plane; set updateflag TRUE");
        break;
      case XK_Left:
        Tcl_Eval(interp,"changeslice down $plane; set updateflag TRUE");
        break;

        /* modifiers */
      case XK_Shift_L:
      case XK_Shift_R:
        shiftkeypressed=TRUE;
        break;
      case XK_Control_L:
      case XK_Control_R:
        ctrlkeypressed=TRUE;
        break;
      case XK_Alt_L:
      case XK_Alt_R:
        altkeypressed=TRUE;
        break;
      }
      break;

    case KeyRelease:   /* added this mask to xwindow.c */
      XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
      switch (ks) {
      case XK_Shift_L:
      case XK_Shift_R:
        shiftkeypressed=FALSE;
        break;
      case XK_Control_L:
      case XK_Control_R:
        ctrlkeypressed=FALSE;
        break;
      case XK_Alt_L:
      case XK_Alt_R:
        altkeypressed=FALSE;
        break;
      case XK_0:
        break;
      }
      break;
    }
    return GL_FALSE;
  }
  return GL_FALSE;
}


/*-----------------------------------------------------*/
void  open_window(char *name) {
  XSizeHints hin;

  if (openglwindowflag) {
    printf("medit: ### GL window already open: can't open second\n");
    return;
  }

  /* TKO_DEPTH because all OpenGL 4096 visuals have depth buffer!! */

#ifdef RGB
  tkoInitDisplayMode(TKO_DOUBLE | TKO_RGB | TKO_DEPTH);
#else
  tkoInitDisplayMode(TKO_DOUBLE | TKO_INDEX | TKO_DEPTH);
#endif

  if (!initpositiondoneflag)
    tkoInitPosition(MOTIF_XFUDGE+wx0,(1024-wy0-ydim)+MOTIF_XFUDGE,xdim,ydim);

  if (!tkoInitWindow(name)) {
    printf("register: ### tkoInitWindow(name) failed\n");
    exit(1);
  }

  printf("Opening %s, xnum = %d, xdim = %d\n",name,xnum,xdim);

  //hin.max_width = hin.max_height = 3*xnum + xnum/2;  /* maxsize */
  hin.max_width = hin.max_height = xdim;
  hin.min_aspect.x = hin.max_aspect.x = xdim;        /* keepaspect */
  hin.min_aspect.y = hin.max_aspect.y = ydim;
  hin.flags = PMaxSize|PAspect;
  XSetWMNormalHints(xDisplay, w.wMain, &hin);

  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glOrtho(0.0, (double)(xnum-1), 0.0, (double)(ynum-1), -1.0, 1.0);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();

  openglwindowflag = TRUE;
}


/*-----------------------------------------------------*/
void  resize_window_intstep() {
#ifdef OPENGL
  int tzf;

  tzf = rint((float)w.w/(float)xnum);
  tzf = (tzf<1)?1:(tzf>3)?3:tzf;
  if (w.w%xnum || w.h%ynum) {
    ozf = zf;
    zf = tzf;
    fsf = (float)zf;
    xdim = zf*xnum;
    ydim = zf*ynum;
    XResizeWindow(xDisplay, w.wMain, xdim, ydim);
    if (TKO_HAS_OVERLAY(w.type))
      XResizeWindow(xDisplay, w.wOverlay, xdim, ydim);
    resize_buffers(xdim, ydim);
    w.w = xdim;
    w.h = ydim;

    imc = (zf*imc)/ozf;
    ic = (zf*ic)/ozf;
    jc = (zf*jc)/ozf;
  }
#endif
}


/*-----------------------------------------------------*/
void move_window(int x, int y) {
#ifdef OPENGL
  if (openglwindowflag) {
    XMoveWindow(xDisplay, w.wMain, x, y);
    w.x = x;
    w.y = y;
  } else if (!initpositiondoneflag) {
    tkoInitPosition(x,y,xdim,ydim);
    initpositiondoneflag = TRUE;
  } else ;
#endif
}


/*-----------------------------------------------------*/
void blinkbuffers() {
  if (blinkflag) {
#ifdef RGB
    /* if(blinkdelay==1) { */
    /* printf("reaload for delay 1\n"); */
    /* reload_buffers(); */
    /* printf("saving frontbuffer");
       glReadPixels(0,0,512,512,GL_RGB,GL_UNSIGNED_BYTE,blinkbuf); */
    /* } */
#endif
    if (blinkdelay<0) {
      record_swapbuffers();
      /*sginap(blinktime);*/
      /*usecnap(blinktime*10000);*/
    } else
      blinkdelay--;
  }
}


/*-----------------------------------------------------*/
void record_swapbuffers()    /* called by compare button */
{
  int swaptmp;

#ifdef RGB
  /* if(blinkflag) {
     if(invalid_buffers == 1) {
     reload_buffers();
     } else {
     if(blinktop==1) {
     rectwrite(0,0,xdim-1,ydim-1,blinkbuft);
     } else {
     rectwrite(0,0,xdim-1,ydim-1,blinkbufm);
     }
     }
     blinktop = (blinktop == 0) ? 1 : 0;
     } */
#endif
  swapbuffers();
  swaptmp = visible_plane;
  visible_plane = last_visible_plane;
  last_visible_plane = swaptmp;

  swaptmp = visible_mode;
  visible_mode = last_visible_mode;
  last_visible_mode = swaptmp;
}


/*-----------------------------------------------------*/
void resize_buffers(int x,int y) {
  free(vidbuf);
  free(binbuff);
  free(blinkbuft);
  free(blinkbufm);

#ifdef RGB
  vidbuf = (GLubyte *)lcalloc(3*(size_t)(x*y),(size_t)sizeof(GLubyte));
  blinkbuft = (GLubyte *)lcalloc(3*(size_t)(x*y),(size_t)sizeof(GLubyte));
  blinkbufm = (GLubyte *)lcalloc(3*(size_t)(x*y),(size_t)sizeof(GLubyte));
#else
  vidbuf = (Colorindex *)lcalloc((size_t)(x*y),(size_t)sizeof(Colorindex));
#endif
  binbuff = (unsigned char *)lcalloc((size_t)(3*x*y),(size_t)sizeof(char));
}


#endif // HAVE_TCL_TK_GL


/*-----------------------------------------------------*/
void  read_reg(char *fname) {
  extern char subjectid[1000];
  extern double ps_2, st_2, fscale_2;
  extern float tm[4][4];
  extern int float2int;
  extern MATRIX *RegMat;

  int i,j,err;
  char *tmpstr;
  float ipr, bpr, fscale;

  err = regio_read_register(fname, &tmpstr, &ipr, &bpr,
                            &fscale, &RegMat, &float2int);
  if (err) {
    printf("ERROR: reading %s\n",fname);
    exit(1);
  }

  ps_2 = ipr;
  st_2 = bpr;
  if (fscale_2 == 0) fscale_2 = fscale;
  if(subjectidOverride == 0) strcpy(subjectid,tmpstr);

  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      tm[i][j] = RegMat->rptr[i+1][j+1];
    }
  }

  printf("---- Input registration matrix --------\n");
  MatrixPrint(stdout,RegMat);
  printf("float2int = %d\n",float2int);
  printf("---------------------------------------\n");

  free(tmpstr);
}


/*-----------------------------------------------------*/
void  read_fslreg(char *fname) {
  extern MATRIX *FSLRegMat;
  FILE *fp;
  int i,j,n;

  fp = fopen(fname,"r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n",fname);
    exit(1);
  }

  if (FSLRegMat==NULL) FSLRegMat = MatrixAlloc(4,4,MATRIX_REAL);

  for (i=0;i<4;i++) {
    for (j=0;j<4;j++) {
      n = fscanf(fp,"%f",&(FSLRegMat->rptr[i+1][j+1]));
      if (n != 1) {
        printf("ERROR: reading %s, row %d, col %d\n",fname,i,j);
        exit(1);
      }
    }
  }

  printf("---- FSL registration matrix --------\n");
  MatrixPrint(stdout,FSLRegMat);
  printf("---------------------------------------\n");

  return;
}


/*-----------------------------------------------------*/
void write_reg(char *fname) {
  extern char *fslregoutfname, *subjectsdir, *pname;
  extern char *freeviewfname;
  extern int fstal;
  extern char talxfmfile[2000];
  extern MATRIX *RegMat, *Mtc;
  static MATRIX *RegMatTmp=NULL;
  //  int i,j;
  FILE *fp;
  char touchfile[1000];

  editedmatrix = FALSE;

  RegMatTmp = MatrixMultiply(RegMat,Mtc,RegMatTmp);

  printf("RegMat ---------------------------\n");
  MatrixPrint(stdout,RegMatTmp);

  if(fname != NULL){
    if(fscale_2 == 0.0) fscale_2 = .1;
    make_backup(fname);
    {
      LTA *lta = LTAalloc(1, NULL) ;
      strcpy(lta->subject, pname) ;
      lta->fscale = fscale_2 ;
      lta->xforms[0].m_L = MatrixCopy(RegMatTmp, NULL) ;
      getVolGeom(mov_vol,   &lta->xforms[0].src);
      getVolGeom(targ_vol0, &lta->xforms[0].dst);
      lta->type = REGISTER_DAT;
      //lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);
      if(LTAwrite(lta, fname) != NO_ERROR)
        printf("register: ### can't create file %s\n",fname);
      LTAfree(&lta) ;
    }
  }    

  if(fslregoutfname != NULL) write_fslreg(fslregoutfname);
  if(freeviewfname != NULL) write_freeviewreg(freeviewfname);
  if(ltaoutfname != NULL) write_lta(ltaoutfname);

  if(xfmoutfname != NULL) write_xfmreg(xfmoutfname);
  else if(fstal) {
    if(ZeroCRAS){
      printf("UnZeroing CRAS for fstal output xfm\n");
      RegMatTmp = MatrixMultiply(RegMatTmp,invMcras0,RegMatTmp);
    }
    make_backup(talxfmfile);
    regio_write_mincxfm(talxfmfile,RegMatTmp,xfmfileinfo);
    sprintf(touchfile,"%s/%s/touch",subjectsdir,pname);
    if (!fio_IsDirectory(touchfile)) mkdir(touchfile,0777);
    sprintf(touchfile,"%s/%s/touch/talairach.tkregister2.touch",
            subjectsdir,pname);
    fp = fopen(touchfile,"w");
    fprintf(fp,"talairach registration %s edited by tkregister2\n",talxfmfile);
    fclose(fp);
  }

  return;
}


/*-----------------------------------------------------*/
void write_fslreg(char *fname) {
  extern MRI *mov_vol, *targ_vol0;
  extern MATRIX *RegMat, *Mtc;
  static MATRIX *RegMatTmp=NULL;
  int i,j;
  FILE *fp;
  MATRIX *Mfsl;

  RegMatTmp = MatrixMultiply(RegMat,Mtc,RegMatTmp);

  Mfsl = MRItkreg2FSL(targ_vol0, mov_vol, RegMatTmp);

  fp = fopen(fname,"w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s for writing\n",fslregoutfname);
    return;
  }
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      fprintf(fp,"%13.8f ",Mfsl->rptr[i+1][j+1]);
    fprintf(fp,"\n");
  }
  fclose(fp);

  return;
}


/*-----------------------------------------------------*/
void write_freeviewreg(char *fname) {
  extern MRI *mov_vol, *targ_vol0;
  extern MATRIX *RegMat, *Mtc;
  static MATRIX *RegMatTmp=NULL;
  int i,j;
  FILE *fp;
  MATRIX *Mfv=NULL,*Ttarg, *Tmov, *Starg, *Smov,*InvTmov,*InvStarg;

  RegMatTmp = MatrixMultiply(RegMat,Mtc,RegMatTmp);

  Tmov  = MRIxfmCRS2XYZtkreg(mov_vol);
  Ttarg = MRIxfmCRS2XYZtkreg(targ_vol0);
  Smov  = MRIxfmCRS2XYZ(mov_vol,0);
  Starg = MRIxfmCRS2XYZ(targ_vol0,0);

  InvStarg = MatrixInverse(Starg,NULL);
  InvTmov = MatrixInverse(Tmov,NULL);

  Mfv = MatrixMultiply(Ttarg,InvStarg,NULL);
  Mfv = MatrixMultiply(Mfv,Smov,Mfv);
  Mfv = MatrixMultiply(Mfv,InvTmov,Mfv);
  Mfv = MatrixMultiply(Mfv,RegMatTmp,Mfv);

  printf("FreeView Matrix ---------------------------\n");
  MatrixPrint(stdout,Mfv);

  fp = fopen(fname,"w");
  if (fp==NULL) {
    printf("register: ### can't create file %s\n",fname);
    return;
  }
  fprintf(fp,"%s\n",pname);
  fprintf(fp,"%f\n",ps_2);
  fprintf(fp,"%f\n",st_2);
  if(fscale_2 == 0.0) fscale_2 = .1;
  fprintf(fp,"%f\n",fscale_2);
  for (i=0;i<4;i++) {
    for (j=0;j<4;j++)
      //fprintf(fp,"%e ",tm[i][j]);
      fprintf(fp,"%e ",Mfv->rptr[i+1][j+1]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"round\n");
  printf("register: file %s written\n",fname);
  fclose(fp);

  MatrixFree(&Mfv);
  MatrixFree(&Ttarg);
  MatrixFree(&Tmov);
  MatrixFree(&InvTmov);
  MatrixFree(&Starg);
  MatrixFree(&InvStarg);
  MatrixFree(&Smov);
  MatrixFree(&RegMatTmp);

  return;
}


/*-----------------------------------------------------*/
void write_xfmreg(char *fname) {
  extern MRI *mov_vol, *targ_vol0;
  extern MATRIX *RegMat;
  int i,j;
  FILE *fp;
  MATRIX *Mxfm=NULL, *RegMatTmp=NULL;

  RegMatTmp = MatrixMultiply(RegMat,Mtc,NULL);
  Mxfm = MRItkReg2Native(targ_vol0, mov_vol, RegMatTmp);

  fp = fopen(fname,"w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s for writing\n",fname);
    return;
  }
  fprintf(fp,"MNI Transform File\n");
  fprintf(fp,"%% tkregister2\n");
  fprintf(fp,"\n");
  fprintf(fp,"Transform_Type = Linear;\n");
  fprintf(fp,"Linear_Transform =\n");
  for (i=0;i<3;i++) {
    for (j=0;j<4;j++)
      fprintf(fp,"%13.8f ",Mxfm->rptr[i+1][j+1]);
    if (i != 2) fprintf(fp,"\n");
    else       fprintf(fp,";\n");
  }
  fclose(fp);

  MatrixFree(&Mxfm);
  MatrixFree(&RegMatTmp);

  return;
}


/*-----------------------------------------------------*/
void write_lta(char *fname) {
  extern MRI *mov_vol, *targ_vol0;
  extern MATRIX *RegMat, *Mtc;
  extern char *pname;
  extern int invLTAOut;
  LTA              *lta ;
  MATRIX *RegMatTmp=NULL;

  RegMatTmp = MatrixMultiply(RegMat,Mtc,NULL);
  lta = LTAalloc(1, NULL) ;
  strcpy(lta->subject, pname) ;
  lta->fscale = fscale_2 ;
  if(! invLTAOut) {
    /* The definitions of mov=src and ref=dst are consistent with
       LTAchangeType() and ltaReadRegisterDat(). This is an unfortunate
       definition because the registration matrix actually does from
       ref to mov. But this was an error introduced a long time ago
       and the rest of the code base has built up around it. */
    lta->xforms[0].m_L = MatrixCopy(RegMatTmp, NULL) ;
    getVolGeom(mov_vol,   &lta->xforms[0].src);
    getVolGeom(targ_vol0, &lta->xforms[0].dst);
  } 
  else {
    // Note: cannot just run LTAfillInverse()
    lta->xforms[0].m_L = MatrixInverse(RegMatTmp, NULL) ;
    getVolGeom(mov_vol,   &lta->xforms[0].dst);
    getVolGeom(targ_vol0, &lta->xforms[0].src);
  }
  lta->type = REGISTER_DAT;
  lta = LTAchangeType(lta, LINEAR_VOX_TO_VOX);

  if (LTAwrite(lta, fname) != NO_ERROR)
    printf("register: ### can't create file %s\n",fname);
  LTAfree(&lta) ;
  return;
}


/*-----------------------------------------------------*/
void make_backup(char *fname) {
  char command[2*NAME_LENGTH];
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp!=NULL) {
    sprintf(command,"cp %s %s~",fname,fname);
    system(command);
    fclose(fp);
  }
}


#ifdef HAVE_TCL_TK_GL

/*-----------------------------------------------------*/
void save_rgb(char *fname) {
  if (scrsaveflag) {
    scrsave_to_rgb(fname);
  } else             {
    pix_to_rgb(fname);
  }
}
void pix_to_rgb(char *fname) {}


/*-----------------------------------------------------*/
void scrsave_to_rgb(char *fname)  /* about 2X faster than pix_to_rgb */
{
  char command[2*NAME_LENGTH];
  FILE *fp;
  long xorig,xsize,yorig,ysize;
  int x0,y0,x1,y1;

  getorigin(&xorig,&yorig);
  getsize(&xsize,&ysize);

  x0 = (int)xorig;
  x1 = (int)(xorig+xsize-1);
  y0 = (int)yorig;
  y1 = (int)(yorig+ysize-1);
  fp = fopen(fname,"w");
  if (fp==NULL) {
    printf("register: ### can't create file %s\n",fname);
    return;
  }
  fclose(fp);
  sprintf(command,"scrsave %s %d %d %d %d\n",fname,x0,x1,y0,y1);
  system(command);
  printf("register: file %s written\n",fname);
}

#endif // HAVE_TCL_TK_GL


/*-----------------------------------------------------*/
void downslice() {
  if (plane==CORONAL)
    imc = (imc<zf)?imnr1*zf-zf+imc:imc-zf;
  else if (plane==HORIZONTAL)
    ic = (ic<zf)?ydim-zf+ic:ic-zf;
  else if (plane==SAGITTAL)
    jc = (jc<zf)?xdim-zf+jc:jc-zf;
}


/*-----------------------------------------------------*/
void upslice() {
  if (plane==CORONAL)
    imc = (imc>=imnr1*zf-zf)?imc+zf-imnr1*zf:imc+zf;
  else if (plane==HORIZONTAL)
    ic = (ic>=ydim-zf)?ic+zf-ydim:ic+zf;
  else if (plane==SAGITTAL)
    jc = (jc>=xdim-zf)?jc+zf-xdim:jc+zf;
}


/*-----------------------------------------------------*/
void goto_point(char *dir) {
  char fname[NAME_LENGTH];
  FILE *fp;
  float xpt,ypt,zpt;

  sprintf(fname,"%s/edit.dat",dir);
  fp=fopen(fname,"r");
  if (fp==NULL) {
    printf("register: ### File %s not found\n",fname);
    return;
  }
  fscanf(fp,"%f %f %f",&xpt,&ypt,&zpt);
  fclose(fp);
  set_cursor(xpt,ypt,zpt);
}


/*-----------------------------------------------------*/
void write_point(char *dir) {
  char fname[NAME_LENGTH];
  FILE *fp;
  float xpt,ypt,zpt;

  sprintf(fname,"%s/edit.dat",dir);
  fp=fopen(fname,"w");
  if (fp==NULL) {
    printf("register: ### can't create file %s\n",fname);
    return;
  }
  xpt = xx1-ps*jc/fsf;
  ypt = yy0+st*imc/fsf;
  zpt = zz1-ps*(255.0-ic/fsf);
  fprintf(fp,"%f %f %f\n",xpt,ypt,zpt);
  fclose(fp);
}


/*-----------------------------------------------------*/
void rotate_brain(float a,char c) {
  int i,j,k;
  float m1[4][4],m2[4][4];
  float sa,ca;

  if (debug) {
    printf("rotating: a = %g deg, c = %c\n",a/10.0,c);
    printf("xc=%g, yc=%g, zc=%g\n",xc,yc,zc);
  }

  if (c=='x') {
    translate_brain(yc,'y');
    translate_brain(zc,'z');
  } else
    if (c=='z') {
      translate_brain(xc,'x');
      translate_brain(yc,'y');
    } else
      if (c=='y') {
        translate_brain(xc,'x');
        translate_brain(zc,'z');
      }

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  a = a*M_PI/1800;
  sa = sin(a);
  ca = cos(a);
  if (c=='y') {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  } else
    if (c=='x') {
      m1[1][1] = m1[2][2] = ca;
      m1[1][2] = -(m1[2][1] = sa);
    } else
      if (c=='z') {
        m1[0][0] = m1[1][1] = ca;
        m1[0][1] = -(m1[1][0] = sa);
      } else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
    for (j=0;j<4;j++) {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += tm[i][k]*m1[k][j];
    }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      tm[i][j] = m2[i][j];

  if (c=='x') {
    translate_brain(-yc,'y');
    translate_brain(-zc,'z');
  } else
    if (c=='z') {
      translate_brain(-xc,'x');
      translate_brain(-yc,'y');
    } else
      if (c=='y') {
        translate_brain(-xc,'x');
        translate_brain(-zc,'z');
      }
  editedmatrix = TRUE;
}

void align_points() {
  if (plane==SAGITTAL) {
    translate_brain(-(yc-yc_old),'y');
    translate_brain(-(zc-zc_old),'z');
  }
  if (plane==HORIZONTAL) {
    translate_brain(-(xc-xc_old),'x');
    translate_brain(-(yc-yc_old),'y');
  }
  if (plane==CORONAL) {
    translate_brain(-(xc-xc_old),'x');
    translate_brain(-(zc-zc_old),'z');
  }
}


/*-----------------------------------------------------*/
void translate_brain(float a, char c) {
  int i,j,k;
  float m1[4][4],m2[4][4];

  if (debug) printf("translating: a = %g, c = %c\n",a,c);

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
    for (j=0;j<4;j++) {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += tm[i][k]*m1[k][j];
    }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      tm[i][j] = m2[i][j];
  editedmatrix = TRUE;
}


/*-----------------------------------------------------*/
void scale_brain(float s, char c) {
  int i,j,k;
  float m1[4][4],m2[4][4];

  if (debug) printf("scaling: s = %g, c = %c\n",s,c);

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  if (c=='x')
    m1[0][0] = s;
  else if (c=='y')
    m1[1][1] = s;
  else if (c=='z')
    m1[2][2] = s;
  else printf("Illegal axis %c\n",c);
  for (i=0;i<4;i++)
    for (j=0;j<4;j++) {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += tm[i][k]*m1[k][j];
    }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      tm[i][j] = m2[i][j];
  editedmatrix = TRUE;
}


/*-----------------------------------------------------*/
void mirror_brain() {
  int j;

  for (j=0;j<4;j++)
    tm[0][j] = -tm[0][j];
  editedmatrix = TRUE;
}


/*-----------------------------------------------------*/
void set_cursor(float xpt,float ypt,float zpt) {
  double dzf;

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
  jpt=ipt=impt = -1;
}


#ifdef HAVE_TCL_TK_GL

/*-----------------------------------------------------*/
void set_scale() 
{
  Colorindex i;
  float f;
  short v;

  for (i=0;i<NUMVALS;i++) {
    /*
      f = i/fscale+fthresh;
      f = ((f<0.0)?0.0:((f>1.0)?1.0:f));
      f = pow(f,fsquash);
      v = f*fscale+0.5;
    */
    f = i/fscale;
    f = 1.0/(1.0+exp(-fsquash*(f-fthresh)));
    v = f*fscale+0.5;
#ifdef RGB
    colormap[i] = (unsigned char)v;
#else
    mapcolor(i+MAPOFFSET,v,v,v);
#endif
  }
  mapcolor(NUMVALS+MAPOFFSET,v,0.0,0.0);
}


/*-----------------------------------------------------*/
void redraw() {
  color(0);
  clear();
  draw_image2(imc,ic,jc);

  last_visible_mode = visible_mode;
  visible_mode = overlay_mode;

  last_visible_plane = visible_plane;
  visible_plane = plane;

}


/*-----------------------------------------------------*/
void pop_gl_window() {
#ifdef OPENGL
  XRaiseWindow(xDisplay, w.wMain);
#else
  winpop();
#endif
}

#endif // HAVE_TCL_TK_GL



/*-----------------------------------------------------*/
void mri2pix(float xpt,float ypt,float zpt,int *jpt,int *ipt,int *impt) {
  if (ptype==0) /* Horizontal */
  {
    *jpt = nint(((xx1-xpt)/ps+0.5));
    *ipt = nint((ypt-yy0)/ps+0.5);
    *impt = nint((zpt-zz0)/st+0.5);
  } else if (ptype==2) /* Coronal */
  {
    *jpt = nint((xx1-xpt)/ps+0.5);
    *ipt = nint((255.0-(zz1-zpt)/ps)+0.5);
    *impt = nint((ypt-yy0)/st+0.5);
  } else if (ptype==1) /* Sagittal */
  {
    *jpt = nint((xx1-xpt)/ps+0.5);
    *ipt = nint((ypt-yy0)/ps+0.5);
    *impt = nint((zpt-zz0)/st+0.5);
  }
}


/*-----------------------------------------------------*/
int imval(float px,float py,float pz) {
  float x,y,z;
  int j=0,i=0,imn=0;

  x = px*TM[0][0]+py*TM[0][1]+pz*TM[0][2]+TM[0][3]+par[0];
  y = px*TM[1][0]+py*TM[1][1]+pz*TM[1][2]+TM[1][3]+par[1];
  z = px*TM[2][0]+py*TM[2][1]+pz*TM[2][2]+TM[2][3]+par[2];
  mri2pix(x,y,z,&j,&i,&imn);
  if (imn>=0&&imn<numimg&&i>=0&&i<ynum&&j>=0&&j<xnum)
    return(im[imn][ynum-1-i][j]);
  else return 0;
}


/*-----------------------------------------------------*/
float Error(int p,float dp) {
  int i,num;
  float mu,error,sum;
  float mu1,mu2,sum1,sum2;

  if (p>=0)
    par[p] += dp;
  mu = mu1 = mu2 = 0;
  num = 0;
  for (i=0;i<npts;i++) {
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
  for (i=0;i<npts;i++) {
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


/*-----------------------------------------------------*/
void optimize(int maxiter) {
  float lambda = 0.03;
  float epsilon = 0.1;
  float momentum = 0.8;
  int iter,p;
  float dE[3];
  float error;

  for (iter=0;iter<maxiter;iter++) {
    error = Error(-1,0);
    printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
           iter,par[0],par[1],par[2],error);
    for (p=0;p<3;p++) {
      dE[p] = tanh((Error(p,epsilon/2)-Error(p,-epsilon/2))/epsilon);
    }
    for (p=0;p<3;p++) {
      par[p] += (dpar[p] = momentum*dpar[p] - lambda*dE[p]);
    }
  }
  error = Error(-1,0);
  printf("%d: %5.2f %5.2f %5.2f %7.3f\n",
         iter,par[0],par[1],par[2],error);
}


/*-----------------------------------------------------*/
void optimize2() {
  float epsilon = 0.5;
  /* int p; */
  float p0,p1,p2,p0min,p1min,p2min;
  float error,minerror;

  error = Error(-1,0);
  minerror = error;
  p0min = p1min = p2min = 0 ;
  printf("%5.2f %5.2f %5.2f %7.3f\n",
         par[0],par[1],par[2],error);
  for (p0 = -10;p0 <= 10;p0 += epsilon)
    for (p1 = -10;p1 <= 10;p1 += epsilon)
      for (p2 = -10;p2 <= 10;p2 += epsilon) {
        par[0] = p0;
        par[1] = p1;
        par[2] = p2;
        error = Error(-1,0);
        if (error<minerror) {
          printf("%5.2f %5.2f %5.2f %7.3f\n",
                 par[0],par[1],par[2],error);
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
         par[0],par[1],par[2],error);
}


/*-----------------------------------------------------*/
void read_images(char *fpref) {
  int i,j,k;                   /* loop counters */
  FILE *fptr;
  char fname[NAME_LENGTH];

  sprintf(fname,"%s.info",fpref);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {
    printf("register: ### File %s not found\n",fname);
    exit(1);
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
  /*
    printf("imnr0=%d,imnr1=%d,numimg=%d,xdim=%d,ydim=%d\n",
    imnr0,imnr1,numimg,xdim,ydim);
  */

  /* Allocate memory */

  bufsize = ((unsigned long)xnum)*ynum;
  buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
  for (k=0;k<numimg;k++) {
    im[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++) {
      im[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
  }
  for (k=0;k<6;k++) {
    sim[k] = (unsigned char **)lcalloc(ynum,sizeof(char *));
    for (i=0;i<ynum;i++) {
      sim[k][i] = (unsigned char *)lcalloc(xnum,sizeof(char));
    }
  }

  for (k=0;k<numimg;k++) {
    changed[k] = FALSE;
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"rb");
    if (fptr==NULL) {
      printf("register: ### File %s not found\n",fname);
      exit(1);
    }
    fread(buf,sizeof(char),bufsize,fptr);
    buffer_to_image(buf,im[k],xnum,ynum);
    fclose(fptr);
    /*
      printf("file %s read in\n",fname);
    */
  }
  printf("register: done reading target COR images\n");

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
}


void read_second_images(char *fpref) {
  unsigned long n;
  int i,j,k,xdim_2b,ydim_2b;
  double ps_2b,st_2b;
  FILE *fptr;
  char fname[NAME_LENGTH];

  xdim_2b=xdim_2;
  ydim_2b=ydim_2;
  ps_2b=ps_2;
  st_2b=st_2;

  sprintf(fname,"%s.info",fpref);
  fptr = fopen(fname,"r");
  if (fptr==NULL) {
    printf("register: ### File %s not found\n",fname);
    exit(1);
  }
  fscanf(fptr,"%*s %d",&imnr0_2);
  fscanf(fptr,"%*s %d",&imnr1_2);
  fscanf(fptr,"%*s %d",&ptype);
  fscanf(fptr,"%*s %d",&xnum_2);
  fscanf(fptr,"%*s %d",&ynum_2);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %lf",&ps_2);
  fscanf(fptr,"%*s %lf",&st_2);
  fscanf(fptr,"%*s %*f");
  fscanf(fptr,"%*s %f",&xx0_2); /* strtx */
  fscanf(fptr,"%*s %f",&xx1_2); /* endx */
  fscanf(fptr,"%*s %f",&yy0_2); /* strty */
  fscanf(fptr,"%*s %f",&yy1_2); /* endy */
  fscanf(fptr,"%*s %f",&zz0_2); /* strtz */
  fscanf(fptr,"%*s %f",&zz1_2); /* endz */
  ps_2 *= 1000;
  st_2 *= 1000;
  xx0_2 *= 1000;
  xx1_2 *= 1000;
  yy0_2 *= 1000;
  yy1_2 *= 1000;
  zz0_2 *= 1000;
  zz1_2 *= 1000;
  fclose(fptr);
  numimg_2 = imnr1_2-imnr0_2+1;
  xdim_2 = xnum_2;   /* no zoom here */
  ydim_2 = ynum_2;

  if (nslices!=numimg_2) {
    printf("register ### numslices mismatch--analyse.dat ignored \n");
    printf("  nslices: analyse.dat = %d   COR-.info = %d\n",
           nslices,numimg_2);
  }
  if (fabs(ps_2b-ps_2)>0.00001 || fabs(st_2b-st_2)>0.00001) {
    printf("register ### slicethick,pixsize mismatch--register.dat ignored\n");
    printf("  thick:   register.dat = %f   COR-.info = %f\n",st_2b,st_2);
    printf("  pixsize: register.dat = %f   COR-.info = %f\n",ps_2b,ps_2);
  }
  if (xdim_2b!=xdim_2 || ydim_2b!=ydim_2) {
    printf("register ### xdim,ydim mismatch--analyse.dat ignored\n");
    printf("  xdim:  analyse.dat = %d   COR-.info = %d\n",xdim_2b,xdim_2);
    printf("  ydim:  analyse.dat = %d   COR-.info = %d\n",ydim_2b,ydim_2);
  }

  bufsize_2 = ((unsigned long)xnum_2)*ynum_2;
  buf_2 = (unsigned char *)lcalloc(bufsize_2,sizeof(char));
  for (k=0;k<numimg_2;k++) {
    fim_2[k] = (float **)lcalloc(ynum_2,sizeof(float *));
    for (i=0;i<ynum_2;i++)
      fim_2[k][i] = (float *)lcalloc(xnum_2,sizeof(float));
  }

  for (k=0;k<numimg_2;k++) {
    file_name(fpref,fname,k+imnr0_2,"%03d");
    fptr = fopen(fname,"rb");
    if (fptr==NULL) {
      printf("register: ### File %s not found\n",fname);
      exit(1);
    }
    fread(buf_2,sizeof(float),bufsize_2,fptr);
    fclose(fptr);
    n=0;
    for (i=0;i<ynum_2;i++)
      for (j=0;j<xnum_2;j++) {
        fim_2[k][i][j] = buf_2[n++];  /* 64Meg <- 16Meg! */
      }
  }
  printf("register: done reading to-be-registered functional images\n");
}


void transform(float x1,float y1,float z1,
               float *x2,float *y2,float *z2,
               float M[4][4]) {
  *x2 = x1*M[0][0]+y1*M[0][1]+z1*M[0][2]+M[0][3];
  *y2 = x1*M[1][0]+y1*M[1][1]+z1*M[1][2]+M[1][3];
  *z2 = x1*M[2][0]+y1*M[2][1]+z1*M[2][2]+M[2][3];
}


#ifdef HAVE_TCL_TK_GL
/*------------------------------------------------------*/
void blur(float factor)  /* test hack */
{
  FILE *fp;
  long xorig,xsize,yorig,ysize;
  int x0,y0,x1,y1;
  int i,j,k;
  /* short r; */
  char command[2*NAME_LENGTH];

  getorigin(&xorig,&yorig);
  getsize(&xsize,&ysize);
  x0 = (int)xorig;
  x1 = (int)(xorig+xsize-1);
  y0 = (int)yorig;
  y1 = (int)(yorig+ysize-1);

  sprintf(command,"scrsave /tmp/tmp1.rgb %d %d %d %d\n",x0,x1,y0,y1);
  system(command);
  sprintf(command,"blur /tmp/tmp1.rgb /tmp/tmp2.rgb %d\n",
          (int)((float)xsize/factor));
  system(command);
  sprintf(command,"tobin /tmp/tmp2.rgb /tmp/tmp2.bin\n");
  system(command);

  fp = fopen("/tmp/tmp2.bin","r");
  fread(binbuff,3,xdim*ydim,fp);
  k = 0;
  for (i=0;i<ydim;i++)
    for (j=0;j<xdim;j++) {
      /*vidbuf[k] = binbuff[(i*xdim+j)*3+0]+MAPOFFSET;*/  /* x9 */

#ifdef RGB
      vidbuf[3*k] =
        binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
      vidbuf[3*k+1] =
        binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
      vidbuf[3*k+2] =
        binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
#else
      vidbuf[k] = binbuff[i*xdim+j]+MAPOFFSET;  /* stretched dark to light */
#endif
      k++;
    }
  backbuffer(FALSE);
  frontbuffer(TRUE);
  rectwrite(0,0,xdim-1,ydim-1,vidbuf);
  backbuffer(TRUE);
  frontbuffer(FALSE);
  sprintf(command,"rm -f /tmp/tmp1.rgb /tmp/tmp2.rgb /tmp/tmp2.bin\n");
  system(command);
}
#endif // HAVE_TCL_TK_GL


void make_filenames(char *lsubjectsdir) {
  subjectsdir = (char *)malloc(NAME_LENGTH*sizeof(char)); /* malloc for tcl */
  srname = (char *)malloc(NAME_LENGTH*sizeof(char));
  psrname = (char *)malloc(NAME_LENGTH*sizeof(char));
  pname = (char *)malloc(NAME_LENGTH*sizeof(char));
  regfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  afname = (char *)malloc(NAME_LENGTH*sizeof(char));
  targpref = (char *)malloc(NAME_LENGTH*sizeof(char));
  movformat = (char *)malloc(NAME_LENGTH*sizeof(char));
  tfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  sgfname = (char *)malloc(NAME_LENGTH*sizeof(char));

  strcpy(subjectsdir,lsubjectsdir);
  strcpy(movformat,"");
  /* TODO: init others here */
}


void read_float_images
(
  float ***fim,
  char *format,
  int nslices,
  int nperslice,
  int xdim,
  int ydim,
  short **buf) /* was MRIio.c */
{
#ifdef Linux
  int scount ;
#endif
  int i,j,k,n,t,imnr,imtype,bufsize;
  char fname[NAME_LENGTH],*ext;
  float max,min,sum,sum2; /* ,avg,stdev, */
  float f;
  FILE *fp = NULL;
  long offset;

  ext = &format[strlen(format)-1];
  while (ext>format && ext[0]!='.') ext--;
  if (!strcmp(ext,".im")) {
    imtype = APD1;
    printf("register: input file type: im (x,y: <pref><z*tcnt+t>.im)\n");
  } else if (!strcmp(ext,".bshort")) {
    imtype = BSHORT;
    printf("register: input file type: bshort (x,y,t: <prefix><z>.bshort)\n");
  } else if (!strcmp(ext,".skip")) {
    imtype = SKIP;
    printf("register: input file type: skip (x,y: <prefix>.<z>.<t>.skip)\n");
  } else if (!strcmp(ext,".BRIK")) {
    imtype = AFNI;
    printf("register: input file type: AFNI (x,y,z,t: <prefix>.BRIK)\n");
  } else if (!strcmp(ext,".brik")) {
    imtype = AFNI;
    printf("register: input file type: AFNI (x,y,z,t: <prefix>.brik)\n");
  } else {
    printf("register: ### file format %s (%s) not supported\n",format,ext);
    exit(1);
  }

  for (k=0;k<nslices;k++) {
    fim[k] = (float **)calloc(ydim,sizeof(float *));
    for (i=0;i<ydim;i++) {
      fim[k][i] = (float *)calloc(xdim,sizeof(float));
    }
  }
  bufsize = xdim*ydim;
  if ((*buf)==NULL)
    *buf = (short *)calloc(bufsize,sizeof(float));

  if (imtype==AFNI) { /* one file per scan (x,y,z,t)  */
    sprintf(fname,"%s",format);
    fp = fopen(fname,"r");
    if (fp==NULL) {
      printf("register: ### File %s not found\n",fname);
      return;
    }
  }
  t = 0;
  for (k=0;k<nslices;k++) {
    if (imtype==APD1) { /* x,y file skip header */
      imnr = k*nperslice + t;
      sprintf(fname,format,imnr);
      fp = fopen(fname,"r");
      if (fp==NULL) {
        printf("register: ### File %s not found\n",fname);
        return;
      }
      fread(*buf,sizeof(char),HEADERSIZE,fp);
    }
    if (imtype==BSHORT) { /* one file per slice (x,y,t) */
      imnr = k;
      sprintf(fname,format,imnr);
      fp = fopen(fname,"r");
      if (fp==NULL) {
        printf("register: ### File %s not found\n",fname);
        return;
      }
    }
    if (imtype==SKIP) { /* x,y file */
      sprintf(fname,format,k,t);   /* assumes zero-based */
      fp = fopen(fname,"r");
      if (fp==NULL) {
        printf("register: ### File %s not found\n",fname);
        return;
      }
    }
    if (imtype==AFNI) { /* brick: skip t*zcnt + z */
      offset = (t*nslices + k) * bufsize * sizeof(short);
      fseek(fp,offset,SEEK_SET);
    }
    if (k==0) printf("First slice read: %s\n",fname);
    fflush(stdout);
    fread(*buf,sizeof(short),bufsize,fp);
#ifdef Linux
    for (scount=0; scount<bufsize; scount++) {
      (*buf)[scount]=swapShort((*buf)[scount]);
    }
#endif
    if (imtype==APD1 || imtype==BSHORT || imtype==SKIP) fclose(fp);
    max = -1.0e10;
    min = 1.0e10;
    sum = sum2 = 0;
    n = 0;
    for (i=0;i<ydim;i++)
      for (j=0;j<xdim;j++) {
        f = fim[k][i][j] = (*buf)[n++];
        if (f<min) min = f;
        if (f>max) max = f;
        sum += f;
        sum2 += f*f;
      }
    sum /= xdim*ydim;
    sum2 = sqrt(sum2/(xdim*ydim)-sum*sum);
    /*printf("File %s read\n",fname);*/
    /*printf("min=%f,max=%f,avg=%f,stdev=%f\n",min,max,sum,sum2);*/
  }
  if (imtype==AFNI) fclose(fp);
  printf("Last slice read:  %s\n",fname);
  fflush(stdout);
}

/* boilerplate wrap function defines for easier viewing */
#define WBEGIN (ClientData clientData,Tcl_Interp *interp,\
int argc,char *argv[]){
#define ERR(N,S)  if(argc!=N){interp->result=S;return TCL_ERROR;}
#define WEND   return TCL_OK;}
#define REND  (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL


#ifdef HAVE_TCL_TK_GL
/*=======================================================================*/
/* function wrappers and errors */
int                  W_move_window  WBEGIN
ERR(3,"Wrong # args: move_window <x> <y>")
  move_window(atoi(argv[1]),atoi(argv[2]));
  WEND

  int                  W_pop_gl_window  WBEGIN
  ERR(1,"Wrong # args: pop_gl_window")
  pop_gl_window();
  WEND

  int                  W_redraw  WBEGIN
  ERR(1,"Wrong # args: redraw")
  redraw();
  WEND

  int                  W_upslice WBEGIN
  ERR(1,"Wrong # args: upslice")
  upslice();
  WEND

  int                  W_downslice WBEGIN
  ERR(1,"Wrong # args: downslice")
  downslice();
  WEND

  int                  W_rotate_brain_x  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_x <deg>")
  rotate_brain(atof(argv[1]),'x');
  WEND

  int                  W_rotate_brain_y  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_y <deg>")
  rotate_brain(atof(argv[1]),'y');
  WEND

  int                  W_rotate_brain_z  WBEGIN
  ERR(2,"Wrong # args: rotate_brain_z <deg>")
  rotate_brain(atof(argv[1]),'z');
  WEND

  int                  W_translate_brain_x  WBEGIN
  ERR(2,"Wrong # args: translate_brain_x <mm>")
  translate_brain(atof(argv[1]),'x');
  WEND

  int                  W_translate_brain_y  WBEGIN
  ERR(2,"Wrong # args: translate_brain_y <mm>")
  translate_brain(atof(argv[1]),'y');
  WEND

  int                  W_translate_brain_z  WBEGIN
  ERR(2,"Wrong # args: translate_brain_z <mm>")
  translate_brain(atof(argv[1]),'z');
  WEND

  int                  W_scale_brain_x  WBEGIN
  ERR(2,"Wrong # args: scale_brain_x <mm>")
  scale_brain(atof(argv[1]),'x');
  WEND

  int                  W_scale_brain_y  WBEGIN
  ERR(2,"Wrong # args: scale_brain_y <mm>")
  scale_brain(atof(argv[1]),'y');
  WEND

  int                  W_scale_brain_z  WBEGIN
  ERR(2,"Wrong # args: scale_brain_z <mm>")
  scale_brain(atof(argv[1]),'z');
  WEND

  int                  W_goto_point  WBEGIN
  ERR(1,"Wrong # args: goto_point")
  goto_point(tfname);
  WEND

  int                  W_write_point  WBEGIN
  ERR(1,"Wrong # args: write_point")
  write_point(tfname);
  WEND

  int                  W_save_rgb  WBEGIN
  ERR(1,"Wrong # args: save_rgb")
  save_rgb(sgfname);
  WEND

  int                  W_record_swapbuffers  WBEGIN
  ERR(1,"Wrong # args: record_swapbuffers")
  record_swapbuffers();
  WEND

  int                  W_set_scale  WBEGIN
  ERR(1,"Wrong # args: set_scale")
  set_scale();
  WEND

  int                  W_blur  WBEGIN
  ERR(2,"Wrong # args: blur <factor>")
  blur(atof(argv[1]));
  WEND

  int                  W_read_reg  WBEGIN
  ERR(1,"Wrong # args: read_reg")
  read_reg(regfname);
  WEND

  int                  W_write_reg  WBEGIN
  ERR(1,"Wrong # args: write_reg")
  write_reg(regfname);
  WEND

  int                  W_align_points  WBEGIN
  ERR(1,"Wrong # args: align_points")
  align_points();
  WEND

  int                  W_mirror_brain  WBEGIN
  ERR(1,"Wrong # args: mirror_brain")
  mirror_brain();
  WEND
/*===============================================================*/

/* for tcl/tk */
  static void StdinProc _ANSI_ARGS_((ClientData clientData, int mask));
static void Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
static Tcl_Interp *interp;
static Tcl_DString command;
static int tty;

#else
static void* interp;
#endif // HAVE_TCL_TK_GL


int main(int argc, char **argv)   /* new main */
{
  int nargs;

  nargs = handleVersionOption(argc, argv, "tkregister2");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

#ifdef HAVE_TCL_TK_GL
  initcolormap();
#endif // HAVE_TCL_TK_GL

  /* start program, now as function; gl window not opened yet */
  //printf("tkregister: starting register\n");
  Register((ClientData) NULL, &interp, argc, argv);/* event loop commented out*/

#ifdef HAVE_TCL_TK_GL

  FILE *fp ;

  /* get tkregister tcl startup script location from environment */
  char *envptr = getenv("FREESURFER_HOME");
  if (envptr==NULL) {
    printf("tkregister: env var FREESURFER_HOME undefined (use setenv)\n");
    printf("    [dir containing mri distribution]\n");
    exit(1);
  }
  if(tkregister_tcl == NULL){
    sprintf(tmpstr,"%s/tktools/%s",envptr,"tkregister2.tcl");
    tkregister_tcl = strcpyalloc(tmpstr);
  }
  printf("tkregister_tcl %s\n",tkregister_tcl);
  if ((fp=fopen(tkregister_tcl,"r"))==NULL) {
    printf("tkregister2: startup script %s not found\n",tkregister_tcl);
    exit(1);
  } else fclose(fp);

  /* start tcl/tk; first make interpreter */
  interp = Tcl_CreateInterp();

  /* make main window (not displayed until event loop starts) */

  /* set the "tcl_interactive" variable */
  tty = isatty(0);
  Tcl_SetVar(interp, "tcl_interactive", (tty) ? "1" : "0", TCL_GLOBAL_ONLY);
  if (tty) promptflag = TRUE;

  /* read tcl/tk internal startup scripts */
  if (Tcl_Init(interp) == TCL_ERROR) {
    fprintf(stderr, "Tcl_Init failed: %s\n", interp->result);
  }
  if (Tk_Init(interp)== TCL_ERROR) {
    fprintf(stderr, "Tk_Init failed: %s\n", interp->result);
  }

  /*=======================================================================*/
  /* register wrapped surfer functions with interpreter */
  Tcl_CreateCommand(interp, "redraw",
                    (Tcl_CmdProc *)   W_redraw,             REND);
  Tcl_CreateCommand(interp, "move_window",
                    (Tcl_CmdProc *)   W_move_window,        REND);
  Tcl_CreateCommand(interp, "pop_gl_window",
                    (Tcl_CmdProc *)   W_pop_gl_window,      REND);
  Tcl_CreateCommand(interp, "upslice",
                    (Tcl_CmdProc *)   W_upslice,            REND);
  Tcl_CreateCommand(interp, "downslice",
                    (Tcl_CmdProc *)   W_downslice,          REND);
  Tcl_CreateCommand(interp, "rotate_brain_x",
                    (Tcl_CmdProc *)   W_rotate_brain_x,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_y",
                    (Tcl_CmdProc *)   W_rotate_brain_y,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_z",
                    (Tcl_CmdProc *)   W_rotate_brain_z,     REND);
  Tcl_CreateCommand(interp, "translate_brain_x",
                    (Tcl_CmdProc *)  W_translate_brain_x,  REND);
  Tcl_CreateCommand(interp, "translate_brain_y",
                    (Tcl_CmdProc *)  W_translate_brain_y,  REND);
  Tcl_CreateCommand(interp, "translate_brain_z",
                    (Tcl_CmdProc *)  W_translate_brain_z,  REND);
  Tcl_CreateCommand(interp, "scale_brain_x",
                    (Tcl_CmdProc *)  W_scale_brain_x,      REND);
  Tcl_CreateCommand(interp, "scale_brain_y",
                    (Tcl_CmdProc *)  W_scale_brain_y,      REND);
  Tcl_CreateCommand(interp, "scale_brain_z",
                    (Tcl_CmdProc *)  W_scale_brain_z,      REND);
  Tcl_CreateCommand(interp, "goto_point",
                    (Tcl_CmdProc *)  W_goto_point,         REND);
  Tcl_CreateCommand(interp, "write_point",
                    (Tcl_CmdProc *)  W_write_point,        REND);
  Tcl_CreateCommand(interp, "save_rgb",
                    (Tcl_CmdProc *)  W_save_rgb,           REND);
  Tcl_CreateCommand(interp, "record_swapbuffers",
                    (Tcl_CmdProc *) W_record_swapbuffers, REND);
  Tcl_CreateCommand(interp, "set_scale",
                    (Tcl_CmdProc *)  W_set_scale,          REND);
  Tcl_CreateCommand(interp, "blur",
                    (Tcl_CmdProc *)  W_blur,               REND);
  Tcl_CreateCommand(interp, "read_reg",
                    (Tcl_CmdProc *)  W_read_reg,           REND);
  Tcl_CreateCommand(interp, "write_reg",
                    (Tcl_CmdProc *)  W_write_reg,          REND);
  Tcl_CreateCommand(interp, "align_points",
                    (Tcl_CmdProc *)  W_align_points,       REND);
  Tcl_CreateCommand(interp, "mirror_brain",
                    (Tcl_CmdProc *)  W_mirror_brain,       REND);
  /*=======================================================================*/
  /***** link global BOOLEAN variables to tcl equivalents */
  Tcl_LinkVar(interp,"maxflag",(char *)&maxflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"updateflag",(char *)&updateflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"blinkflag",(char *)&blinkflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"editedmatrix",(char *)&editedmatrix, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"maskflag",(char *)&maskflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"promptflag",(char *)&promptflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"followglwinflag",(char *)&followglwinflag,
              TCL_LINK_BOOLEAN);
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
  Tcl_LinkVar(interp,"blinktop",(char *)&blinktop, TCL_LINK_INT);
  Tcl_LinkVar(interp,"overlay_mode",(char *)&overlay_mode, TCL_LINK_INT);
  Tcl_LinkVar(interp,"visible_mode",(char *)&visible_mode, TCL_LINK_INT);
  Tcl_LinkVar(interp,"last_visible_mode",(char *)&last_visible_mode,
              TCL_LINK_INT);
  Tcl_LinkVar(interp,"visible_plane",(char *)&visible_plane, TCL_LINK_INT);
  Tcl_LinkVar(interp,"last_visible_plane",(char *)&last_visible_plane,
              TCL_LINK_INT);
  Tcl_LinkVar(interp,"blinkdelay",(char *)&blinkdelay, TCL_LINK_INT);
  Tcl_LinkVar(interp,"blinktime",(char *)&blinktime, TCL_LINK_INT);
  Tcl_LinkVar(interp,"mov_frame",(char *)&mov_frame, TCL_LINK_INT);
  /*=======================================================================*/
  /***** link global DOUBLE variables to tcl equivalents (were float) */
  Tcl_LinkVar(interp,"fsquash",(char *)&fsquash, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fthresh",(char *)&fthresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fscale_2",(char *)&fscale_2, TCL_LINK_DOUBLE);
  /*=======================================================================*/
  /***** link global malloced STRING vars */
  Tcl_LinkVar(interp,"home",        (char *)&subjectsdir,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"session",     (char *)&srname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subject",     (char *)&pname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"tkrtitle",    (char *)&tkrtitle,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"registerdat", (char *)&regfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"analysedat",  (char *)&afname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subjtmpdir",  (char *)&tfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"rgb",         (char *)&sgfname,TCL_LINK_STRING);
  /*=======================================================================*/

  /* run tcl/tk startup script to set vars, make interface; no display yet */
  printf("tkregister2: interface: %s\n",tkregister_tcl);
  Tcl_EvalFile(interp,tkregister_tcl);
  if (*interp->result != 0)  printf("%s",interp->result);
  plane = plane_init;

  /* always start up command line shell too */
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
  if (tty)
    Prompt(interp, 0);
  fflush(stdout);
  Tcl_DStringInit(&command);
  Tcl_ResetResult(interp);

  /*Tk_MainLoop();*/  /* standard */

  /* dual event loop (interface window made now) */
  while (Tk_GetNumMainWindows() > 0) {
    while (Tk_DoOneEvent(TK_ALL_EVENTS|TK_DONT_WAIT)) {
      /* do all the tk events; non-blocking */
    }
    do_one_gl_event(interp);
    /*sginap((long)1);*/   /* block for 10 msec */
    usecnap(10000);     /* block for 10 msec */
  }

  Tcl_Eval(interp, "exit");
#endif // HAVE_TCL_TK_GL
  exit(0);
}

void usecnap(int usec) {
  struct timeval delay;

  delay.tv_sec = 0;
  delay.tv_usec = (long)usec;
  select(0,NULL,NULL,NULL,&delay);
}

#ifdef HAVE_TCL_TK_GL
void initcolormap() 
{
  int i;
  for (i=0; i<256; i++)    colormap[i]=i;
  for (i=256; i<512; i++)  colormap[i]=255;
}

void pseudo_swapbuffers() {
  if (pswapnext==MOVEABLE) {
    rectwrite(0,0,xdim-1,ydim-1,blinkbufm);
    pswapnext=TARGET;
  } else {
    rectwrite(0,0,xdim-1,ydim-1,blinkbuft);
    pswapnext=MOVEABLE;
  }
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
      if (tty) {
        Tcl_Eval(interp, "exit");
        exit(1);
      } else     {
        Tk_DeleteFileHandler(0);
      }
      return;
    } else count = 0;
  }
  cmd = Tcl_DStringAppend(&command, input, count);
  if (count != 0) {
    if ((input[count-1] != '\n') && (input[count-1] != ';')) {
      gotPartial = 1;
      goto prompt;
    }
    if (!Tcl_CommandComplete(cmd)) {
      gotPartial = 1;
      goto prompt;
    }
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

  promptCmd =
    (char*)Tcl_GetVar
    (interp,partial ? "tcl_prompt2" : "tcl_prompt1",TCL_GLOBAL_ONLY);
  if (promptCmd == NULL) {
  defaultPrompt:
    if (!partial)
      fputs("% ", stdout);
  } else {
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

/*----------------------------------------------*/
static int MRItagVol(MRI *mri, float val) {
  int r,c,s;
  float min, max;

  MRIlimits(mri,&min,&max);

  for (c=0; c < 10; c+=2)
    for (r=0; r < 10; r+=2)
      for (s=0; s < mri->depth; s++)
        MRIsetVoxVal(mri,c,r,s,0,.9*max);
  return(0);
}

#endif // HAVE_TCL_TK_GL


/*----------------------------------------------------
  MRIisConformant() - checks whether the input volume
  conforms to the COR format, ie, coronally sliced,
  256^3 which voxel size 1mm^3.
  ----------------------------------------------------*/
static int MRIisConformant(MRI *vol) {
  // Voxel size should be 1mm^3
  if (fabs(vol->xsize - 1) > .001) return(0);
  if (fabs(vol->ysize - 1) > .001) return(0);
  if (fabs(vol->zsize - 1) > .001) return(0);

  // Column Direction Cosine should be -1 0 0
  if (fabs(vol->x_r + 1) > .001) return(0);
  if (fabs(vol->x_a) > .001) return(0);
  if (fabs(vol->x_s) > .001) return(0);

  // Row Direction Cosine should be 0 0 -1
  if (fabs(vol->y_r) > .001) return(0);
  if (fabs(vol->y_a) > .001) return(0);
  if (fabs(vol->y_s + 1) > .001) return(0);

  // Slice Direction Cosine should be 0 1 0
  if (fabs(vol->z_r) > .001) return(0);
  if (fabs(vol->z_a - 1) > .001) return(0);
  if (fabs(vol->z_s) > .001) return(0);

  // Dimension should be 256^3
  if (vol->width  != 256) return(0);
  if (vol->depth  != 256) return(0);
  if (vol->height != 256) return(0);

  return(1);
}


/*!
  \fn MATRIX *Load4x4(char *fname)
  \brief Reads in a 4x4 matrix. fname is an ascii file. Lines that begin with
  # are considered comments. Blank lines are ok.
*/
MATRIX *Load4x4(char *fname)
{
  FILE *fp;
  char s[1000];
  int r,c,m;
  MATRIX *mat;
  
  fp = fopen(fname,"r");
  if(!fp) {
    printf("ERROR: cannot open %s for reading\n",fname);
    return(NULL);
  }

  mat = MatrixAlloc(4,4,MATRIX_REAL);
  r = 1;
  c = 1;
  while(1){
    m = fgetc(fp);
    if(m == EOF){
      printf("ERROR: reading %s: EOF at r=%d c=%d \n",fname,r,c);
      fclose(fp);
      return(NULL);
    }
    if(m == '#') {
      // # is a comment, skip entire line
      fgets(s,1000,fp);
      continue; 
    }
    if(m == '\n' || m == '\r') continue;

    ungetc(m,fp); // put it back

    m = fscanf(fp,"%f",&(mat->rptr[r][c]));
    if(m != 1){
      printf("ERROR: reading %s: item r=%d c=%d \n",fname,r,c);
      fclose(fp);
      return(NULL);
    }
    c++;
    if(c == 5){
      c = 1;
      r++;
      if(r == 5) break;
    }
  }
  fclose(fp);
  return(mat);
}

/*-----------------------------------------------------------*/
int AllocBuffs(void)
{
  extern float **movimg, **targimg;
  extern int **surfimg;
  int n;

  movimg = (float **) calloc(WINDOW_ROWS,sizeof(float*));
  targimg = (float **) calloc(WINDOW_ROWS,sizeof(float*));
  surfimg = (int **) calloc(WINDOW_ROWS,sizeof(int*));

  for(n=0; n < WINDOW_ROWS; n++){
    movimg[n] = (float *) calloc(WINDOW_COLS,sizeof(float));
    targimg[n] = (float *) calloc(WINDOW_COLS,sizeof(float));
    surfimg[n] = (int *) calloc(WINDOW_COLS,sizeof(int));
  }
  return(0);
}
static int istringnmatch(const char *str1, const char *str2, int n) {
  if (n > 0  && ! strncasecmp(str1,str2,n)) return(1);
  if (n <= 0 && ! strcasecmp(str1,str2)) return(1);
  return(0);
}
