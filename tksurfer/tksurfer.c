 /*
 * @file  tksurfer.c
 * @brief Tcl/Tk-based cortical surface viewer
 *
 * TkSurfer displays surface data and allows the user to navigate through
 * that data and view it from different orientations. TkSurfer also displays
 * other data types such as functional data and curvature as overlays onto
 * this surface data.
 * See: http://surfer.nmr.mgh.harvard.edu/fswiki/TkSurferGuide
 */
/*
 * Original Author: Martin Sereno and Anders Dale, 1996
 * CVS Revision Info:
 *    $Author: fischl $
 *    $Date: 2016/12/11 14:33:47 $
 *    $Revision: 1.364 $
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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif /* HAVE_CONFIG_H */
#undef VERSION

#include <sys/types.h>
#include <sys/stat.h>
#define TCL
#define TKSURFER
#define TCL8
#include "xDebug.h"
#include "proto.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "const.h"
#include "fio.h"
#include "error.h"
#include "tritri.h"
#include "diag.h"
#include "rgb_image.h"
#include "const.h"
#include "version.h"
#include "fsgdf_wrap.h"
#include "tiffio.h"
#include "gcsa.h"
#include "mri2.h"
#include "path.h"
#include "fsenv.h"
#include "mrishash_internals.h"

#define bzero(b,len) (memset((b), '\0', (len)), (void) 0)

//////////////////////////////////////////////////////

static void resize_brain(float surface_area) ;
static void set_value_label_name(char *label_name, int field) ;
static void drawcb(void) ;
static void read_imag_vals(char *fname) ;
static void read_soltimecourse(char *fname) ;
static void sol_plot(int timept, int plot_type) ;

static void remove_triangle_links(void) ;
static int draw_curvature_line(void) ;

static void move_window(int x,int y) ;
int MRIStransformBrain(MRI_SURFACE *mris,
                       float exx, float exy, float exz,
                       float eyx, float eyy, float eyz,
                       float ezx, float ezy, float ezz) ;
static void save_rgbfile(char *fname, int width, int height,
                         unsigned short *red,
                         unsigned short *green, unsigned short *blue) ;
static void grabPixels(unsigned int width, unsigned int height,
                       unsigned short *red, unsigned short *green,
                       unsigned short *blue) ;

static int  is_val_file(char *fname) ;
void show_flat_regions(char *surf_name, double thresh) ;
void set_area_thresh(float area) ;
void val_to_mark(void) ;
void transform_brain(void) ;
void curv_to_val(void) ;
int read_curv_to_val(char *fname) ;
int read_parcellation(char *parc_fname, char *lut_fname) ;
int read_and_smooth_parcellation(char *parc_fname, char *lut_fname,
                                 int siter, int miter) ;
void read_stds(int cond_no) ;
void val_to_curv(void) ;
void val_to_stat(void) ;
void stat_to_val(void) ;

static void label_from_stats(int field) ;
static void label_to_stat(int which_overlay) ;
static void label_set_stats(float val ) ;
static void f_to_t(void) ;
static void t_to_p(int dof) ;
static void f_to_p(int numer_dof, int denom_dof) ;
int mask_label(char *label_name) ;

static int zero_mean = 0 ;
MRI_SURFACE *mris = NULL, *mris2 = NULL ;
static char *sdir = NULL ;
static char *sphere_reg ;
static char *sphere_reg_contra ;
MRI *mrismask = NULL;
double mrismaskthresh = -1;
char *mrismaskfile = NULL;

static char *mask_cortex_name = NULL ;
static GCSA *Ggcsa = NULL ;
#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_VERSION_MAGIC_NUMBER  16777215

/*#include "surfer.c"*/


/* start of surfer.c */


#if defined(TCL) && defined(TKSURFER)
#include <tcl.h>
#include <tk.h>
#include <tix.h>
/* begin rkt */
#include <blt.h>
#include <signal.h>
/* end rkt */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef Windows_NT
int stricmp(const char* const str1, const char* const  str2)
{
  return strcasecmp(str1, str2);
}
#endif // Windows_NT
#include <unistd.h>
#include <GL/gl.h>
#include "glut.h"
#include "typedefs.h"
#include "mgh_matrix.h"
#include "label.h"
#include "fio.h"
#include "MRIio_old.h"
#include "rgb_image.h"
#include "transform.h"
#include "proto.h"
#include "diag.h"
#include "macros.h"
#include "utils.h"
#include "machine.h"

/* used to ifdef Linux'd, but should always be true */
#ifndef OPENGL
#define OPENGL
#endif
#define TCL

const char *Progname = "surfer" ;

#if defined(Linux) || defined(SunOS) || defined(sun)
#define GL_ABGR_EXT                         0x8000

#endif

#ifdef TCL
#  define PR   {if(promptflag){fputs("% ", stdout);} fflush(stdout);}
#else
#  define PR
#endif

#ifdef OPENGL
#include <X11/keysym.h>
#include "xwindow.h"
#  define RGBcolor(R,G,B)  glColor3ub((GLubyte)(R),(GLubyte)(G),(GLubyte)(B))
#  define translate(X,Y,Z) glTranslatef(X,Y,Z)
/*one list OK;GL_POLYGON slower*/
#if VERTICES_PER_FACE == 4
#  define bgnpolygon()     glBegin(GL_QUADS)
#else
#  define bgnpolygon()     glBegin(GL_TRIANGLES)
#  define bgnquadrangle()     glBegin(GL_QUADS)
#endif
#  define bgnline()        glBegin(GL_LINES)
#  define bgnpoint()       glBegin(GL_POINTS)
#  define v3f(X)           glVertex3fv(X)
#  define n3f(X)           glNormal3fv(X)
#  define endpolygon()     glEnd()
#  define endline()        glEnd()
#  define endpoint()       glEnd()
#  define swapbuffers()    tkoSwapBuffers()
#  define lrectread(X0,Y0,X1,Y1,P) \
            glReadPixels(X0,Y0,X1-X0+1,Y1-Y0+1,GL_ABGR_EXT,GL_UNSIGNED_BYTE,P)
#  define czclear(A,B)     glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
#  define pushmatrix()     glPushMatrix()
#  define popmatrix()      glPopMatrix()
#  define linewidth(X)     glLineWidth((float)(X))
#  define getmatrix(X)     glGetFloatv(GL_MODELVIEW_MATRIX,(float *)(X))
#  define loadmatrix(X)    glLoadMatrixf((float *)(X))
#  define getorigin(X,Y)   *(X) = w.x; *(Y) = 1024 - w.y - w.h
/*WINDOW_REC w;*/
#  define getsize(X,Y)     *(X) = w.w; *(Y) = w.h
#else
#  include <gl.h>
#  include <device.h>
#endif

/* begin rkt */
// Do the following only for RedHat Enterprise Linux
// It seems that the later version of Tix uses ITcl and ITk.
// stupid tix people who cannot handle version.   I had to use gcc version.
// you cannot include itk.h either(producing so many unknowns) either.
#if NEEDS_ITCL_ITK
#ifndef Itcl_Init
int Itcl_Init(Tcl_Interp* interp);
#endif
#ifndef Itk_Init
int Itk_Init(Tcl_Interp* interp);
#endif
#endif

/* make these decls if the headers are screwed up */
#ifndef Blt_Init
int Blt_Init ( Tcl_Interp* interp );
#endif
#ifndef Blt_SafeInit
int Blt_SafeInit ( Tcl_Interp* interp );
#endif

#ifndef Tix_SafeInit
int Tix_SafeInit ( Tcl_Interp* interp );
#endif

/* end rkt */

#ifndef TCL
#  ifndef OPENGL
#    include <gl.h>     /* non-TCL event loop not moved=>X for OpenGL */
#    include <device.h> /* ditto */
#  endif
#endif

#ifndef SQR
#define SQR(x)       ((x)*(x))
#endif

#define MATCH(A,B)   (!strcmp(A,B))
#define MATCH_STR(S) (!strcmp(str,S))
#ifndef FALSE
#  define FALSE  0
#endif
#ifndef TRUE
#  define TRUE   1
#endif

#define IMGSIZE      256
#define NUMVALS      256
#define MAXIM        256
#define MAXMARKED    500000
#define NLABELS      256
#define CMFBINS       30
#define GRAYBINS      15
#define GRAYINCR     0.25   /* mm */
#define MAXVAL       60.0
#define MOTIF_XFUDGE   8
#define MOTIF_YFUDGE  32
#define NAME_LENGTH   STRLEN
#define MAX_DIR_DEPTH  30
#define SCALE_UP_MOUSE   2.0
#define SCALE_UP_LG      1.25
#define SCALE_UP_SM      1.05
#define BLINK_DELAY 30
#define BLINK_TIME 20

#define MESH_LINE_PIX_WIDTH    1
#define CURSOR_LINE_PIX_WIDTH  2
#define SCALEBAR_WIDTH        6
#define SCALEBAR_BRIGHT     128
#define SCALEBAR_MM          10
#define SCALEBAR_XPOS       0.9
#define SCALEBAR_YPOS       0.9
#define COLSCALEBAR_XPOS    0.925
#define COLSCALEBAR_YPOS   -0.95
#define COLSCALEBAR_WIDTH   0.05
#define COLSCALEBAR_HEIGHT  0.5

#define TMP_DIR         "tmp"              /* relative to subjectsdir/pname */
#define LABEL_DIR       "label"            /* ditto */
#define IMAGE_DIR       "rgb"              /* ditto */
#define BRAINSURF_DIR   "surf"             /* ditto */
#define BEM_DIR         "bem"              /* ditto */
#define DEVOLT1_DIR     "mri/T1"           /* ditto */
#define FILLDIR_STEM    "morph/fill"       /* ditto */
#define CURVDIR_STEM    "morph/curv"       /* ditto */
#define EEGMEG_DIR      "eegmeg"           /* ditto */
#define TRANSFORM_DIR   "mri/transforms"   /* ditto */
#define FS_DIR          "fs"               /* relative to session dir */
#define TALAIRACH_FNAME "talairach.xfm"    /* relative to TRANSFORM_DIR */

/* surfcolor */
#define NONE                 0
#define CURVATURE_OR_SULCUS  1
#define AREAL_DISTORTION     2
#define SHEAR_DISTORTION     3
/* colscales for: set_{complexval_,"",stat_,positive_,signed_}color */
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
/* set_color modes */
#define GREEN_RED_CURV     0   /* default: red,green,gray curvature */
#define REAL_VAL           1   /* specified further by some colscales */
#define FIELDSIGN_POS      4   /* blue */
#define FIELDSIGN_NEG      5   /* yellow */
#define BORDER             6   /* yellow */
#define MARKED             7   /* white */


/* project directions in MRI coords (flatten) */
#define RIGHT       1    /* x */
#define LEFT       -1
#define ANTERIOR    2    /* y */
#define POSTERIOR  -2
#define SUPERIOR    3    /* z */
#define INFERIOR   -3
/* draw_spokes options */
#define RADIUS    0
#define THETA     1
/* for curv/fill image flags  */
#define CURVIM_FILLED        0x01
#define CURVIM_DEFINED       0x02
#define CURVIM_DEFINED_OLD   0x04
#define CURVIM_FIXEDVAL      0x08
/* truncated max min for normalize curv (scales byte image) */
#define CURVIM_NORM_MIN  -2.0
#define CURVIM_NORM_MAX   2.0
/* filltypes */
#define RIPFILL     0
#define STATFILL    1
#define CURVFILL    2

/* vrml modes */
#define VRML_POINTS             0
#define VRML_QUADS_GRAY         1
#define VRML_QUADS_GRAY_SMOOTH  2
#define VRML_QUADS_CURV         3
#define VRML_QUADS_CURV_SMOOTH  4

/* begin rkt */
/* value field labels for swap_vertex_fields(), also used in tcl script */
#define FIELD_CURV      0
#define FIELD_CURVBAK   1
#define FIELD_VAL       2
#define FIELD_IMAG_VAL  3
#define FIELD_VAL2      4
#define FIELD_VALBAK    5
#define FIELD_VAL2BAK   6
#define FIELD_STAT      7

/* label sets for update_labels, also used in tcl script */
#define LABELSET_CURSOR     0
#define LABELSET_MOUSEOVER  1

/* label fields for update_labels, also used in tcl script */
enum
{
  LABEL_VERTEXINDEX = 0,
  LABEL_DISTANCE,
  LABEL_COORDS_RAS,
  LABEL_COORDS_MNITAL,
  LABEL_COORDS_TAL,
  LABEL_COORDS_INDEX,
  LABEL_COORDS_NORMAL,
  LABEL_COORDS_SPHERE_XYZ,
  LABEL_COORDS_SPHERE_RT,
  LABEL_CURVATURE,
  LABEL_FIELDSIGN,
  LABEL_VAL,
  LABEL_VAL2,
  LABEL_VALBAK,
  LABEL_VAL2BAK,
  LABEL_VALSTAT,
  LABEL_IMAGVAL,
  LABEL_MEAN,
  LABEL_MEANIMAG,
  LABEL_STDERROR,
  LABEL_AMPLITUDE,
  LABEL_ANGLE,
  LABEL_DEGREE,
  LABEL_LABEL,
  LABEL_ANNOTATION,
  LABEL_MRIVALUE,
  LABEL_PARCELLATION_NAME
};

/* extra error codes */
#define ERROR_UNDO_ACTION -998      /* invalid action struct */
#define ERROR_UNDO_LIST_STATE  -999 /* inapproriate operation for undo
list state */
#define ERROR_ARRAY       -1000     /* problem with xGrowableArray */
#define ERROR_FUNC        -1001     /* other error with func volume */
#define ERROR_NOT_INITED  -1002     /* something that should have been
previously inited, etc, hasn't been */

/* end rkt */

int file_type = MRIS_BINARY_QUADRANGLE_FILE ;

#if 0
FLOATTYPE **E,**B,**D,**We,**Wb,**Wd;

FLOATTYPE **A,**At;
FLOATTYPE **C;
FLOATTYPE *R;
FLOATTYPE **M,**Mi,**Mb;
FLOATTYPE **W;
FLOATTYPE **_Y,**Yi;
FLOATTYPE **Cov;
FLOATTYPE *dipval,*dipval2,*eegval,*eegval2;
FLOATTYPE *megval,*megval2,*depval,*depval2;
FLOATTYPE *dipvalbak,*eegvalbak,*megvalbak,*depvalbak;
FLOATTYPE *eegx,*eegy,*eegz,*megx,*megy,*megz,*depx,*depy,*depz;
FLOATTYPE *tmp1,*tmp2;

int *dipindex;

int ndip,neeg,nmeg,ndep,nsens;

float dipscale,eegscale,megscale,depscale,sensscale;
float dfthresh,dfsquash;
#endif

/* global names for tcl links */
char val_dir[STRLEN] ; /* directory that value (.w) file was read from */
char *subjectsdir; /* $subjectsdir: from getenv */
char *srname;      /* $session(dir!) abs:#/951206MS/image,#/MARTY0928/08192*/
char *pname;       /* name: $home = $subjectsdir/$name */
char *tkstitle=NULL;/* image window title, pname by default, change with -title*/
char *stem;        /* hemisphere (head: e.g., rh from rh.wmsooth) */
char *ext;         /* surf (suffix: e.g., wmsmooth from rh.wmsmooth) */
char *fpref;       /* $home/surf/hemi. (for assembly) */
char *ifname;      /* $home/surf/hemi.surf: read */
char *if2name;     /* $home/surf/hemi.surf: read (curv target) */
char *ofname;      /* set: $home/surf/hemi.surf: write */
char *cfname;      /* $home/surf/hemi.curv */
char *cf2name;     /* $home/surf/hemi.curv--target surface */
char *sfname;      /* $home/surf/hemi.sub */
char *dfname;      /* $home/surf/hemi.dip */
char *kfname;      /* $home/surf/hemi.sulc */
char *mfname;      /* $home/mri/T1/COR- */
char *vfname;      /* set: $home/surf/hemi.val */
char *fsfname;     /* $session/fs/hemi.fs */
char *fmfname;     /* $session/fs/hemi.fm */
char *pname2;      /* name2 */
char *stem2;       /* hemi2. */
char *ext2;        /* surf2 */
char *tf2name;     /* $home/morph/{fill}hemi2.surf2/COR- */
char *afname;      /* $home/surf/hemi.area */
char *pfname;      /* $home/surf/hemi.patch */
char *tfname;      /* (dir!) ~/name/tmp/ */
char *gfname;      /* (dir!) ~/name/rgb/ */
char *sgfname;     /* (dir!) set: get from cwd: $session/rgb/ */
char *agfname;     /* $home/rgb/tksurfer.rgb */
char *fifname;     /* prefix: $home/morph/{fill}.hemi.surf/COR- */
char *cif2name;    /* prefix: $targhome/morph/{curv}.hemi.surf/COR- */
char *rfname;      /* set: scriptfile (undef at start!) */
char *nfname;      /* annotation.rgb */
char *orfname;     /* $home/surf/hemi.orig */
char *paint_fname; /* $home/surf/hemi.smoothwm */
char *elfname;     /* $home/surf/hemi.1000a2ell */
char *lfname;      /* $home/label/hemi-name.label */
char *vrfname;     /* $home/name/surf/rh.smoothwm.wrl */
char *xffname;     /* $home/name/mri/transforms/TALAIRACH_FNAME */
char *orig_suffix = "orig" ;
char *white_suffix = "white" ;
char *sphere_reg_suffix = "sphere.reg" ;
char *sphere_reg_contra_suffix = "sphere.left_right" ;

FILE *fpvalfile;              /* mult frames */
int openvalfileflag = FALSE;
int twocond_flag = FALSE ;    /* are we in two-condition mode */
int disc_flag = FALSE ;
int cond0 = 0 ;
int cond1 = 1 ;
int openedcmpfilenum = 0;
int cmpfilenamedframe = -1;
float cmpfilenamedfirstlat;
FILE *fpcmpfilenamed;        /* for open cmp movie file */

unsigned long *framebuff, *framebuff2, *framebuff3;
unsigned char *binbuff;
int framenum = 0;
long frame_xdim = 600;
long frame_ydim = 600;
int ilat = 0;  /* latency */

int sub_num ;
#if 0
int vertex_index, face_index;
FACE *face;
VERTEX *vertex; int vertex2_index,face2_index;
face2_type *face2;
VERTEX *vertex2;
float total_area;
#endif

int xnum=256,ynum=256;
unsigned long bufsize;
unsigned char **im[MAXIM];
unsigned char **fill[MAXIM];
unsigned char *buf = NULL;  /* scratch */
int imnr0,imnr1,numimg;
int wx0=650,wy0=416;

unsigned char **curvim[MAXIM];
unsigned char **ocurvim[MAXIM];
unsigned char **curvimflags[MAXIM];
#if 0
float curvmin2,curvmax2;
float curvmin,curvmax;
#endif
float avgcurvmin,avgcurvmax;
double icstrength=1.0;   /* icurv force */
int curvim_allocated=FALSE;
int curvim_averaged=0;

static float sf=0.55;      /* initial scale factor */
static float zf=1.0;       /* current scale */
static int reassign = 0 ; // reassign label vertices when reading

double scalebar_xpos = SCALEBAR_XPOS;
double scalebar_ypos = SCALEBAR_YPOS;
int scalebar_bright = SCALEBAR_BRIGHT;
double colscalebar_width = COLSCALEBAR_WIDTH;
double colscalebar_height = COLSCALEBAR_HEIGHT;
double colscalebar_xpos = COLSCALEBAR_XPOS;
double colscalebar_ypos = COLSCALEBAR_YPOS;

int long_config_overlay = TRUE;

double cthk = 1.0;  /* cortical thickness (mm) */
float ostilt = 1.0; /* outside stilt length (mm) */
int mingm = 50;     /* outer edge gray matter */

double mstrength = 0.125;
double mmid = 45.0;       /* was fzero=35.0 */
double mslope = 0.05;     /* was fsteepness */
double whitemid = 45.0;   /* was whitezero=35.0 */
double graymid = 30.0;    /* was grayzero=25.0; */
double fslope = 0.50;     /* was fsquash */
double fmid = 4.0;        /* was fzero */
double foffset = 0 ;
double fcurv = 0.00;
double cslope = 1.00;     /* was fsquash */
double cmid = 0.0;
double cmax = 100000.0f ;
double cmin = -100000.0f ;

float cup = 0.1;   /* dist cursor above surface */
float mup = 0.05;  /* dist edge above surface */
float pup = 0.07;  /* dist point above surface */

int mesh_linewidth = MESH_LINE_PIX_WIDTH ;
int meshr = 255;   /* mesh color */
int meshg = 255;
int meshb = 255;

float xmin,xmax;
float ymin,ymax;
float zmin,zmax;
float st,ps,fov,xx0,xx1,yy0,yy1,zz0,zz1;
float dipscale = 1.0;

static int selection = -1;
static int mouseover_vno = -1;
int nmarked = 0;
int *marked;

int autoflag = FALSE;
int autoscaleflag = FALSE;
int MRIflag = TRUE;
int MRIloaded = TRUE;
int electrodeflag = FALSE;
int explodeflag = FALSE;
int expandflag = FALSE;
int momentumflag = FALSE;
int sulcflag = FALSE;
int avgflag = FALSE;
int stiffnessflag = FALSE;
int areaflag = FALSE;
int flag2d = FALSE;
int shrinkfillflag = FALSE;
int ncthreshflag = FALSE;
int complexvalflag = FALSE;
int fieldsignflag = FALSE;
int wholeflag = FALSE;
int overlayflag = FALSE;
int verticesflag = FALSE;
int revfsflag = FALSE;
int revphaseflag = FALSE;
int invphaseflag = FALSE;
int rectphaseflag = FALSE;
int truncphaseflag = FALSE;
int scalebarflag = FALSE;
int colscalebarflag = FALSE;
int colscalebartextflag = TRUE;
int colscalebartickflag = TRUE;
char *colscalebar_label[4] = {NULL, NULL, NULL, NULL} ;
int colscalebar_font_size = 1; /* 1 is small, 2 is med, 3 is big */
int colscalebarvertflag = TRUE;
int colscalebaruselabelsflag = FALSE;
int linkvertexmode = 0;
int numvertices = 0;
int surfaceflag = TRUE;
int pointsflag = FALSE;
int statflag = FALSE; /* vertex (fMRI) stats read in ? */
int isocontourflag = FALSE;
int phasecontourflag = FALSE;
int curvimflag = FALSE;
int curvimloaded = FALSE;
int curvloaded = FALSE;
int secondsurfaceloaded = FALSE;
int origsurfloaded = FALSE ;
int annotationloaded = FALSE ;
int canonsurfloaded = FALSE ;
int canonsurffailed = FALSE ;
int sphericalsurfloaded = FALSE ;
int white_surf_loaded = FALSE ;
int inflated_surf_loaded = FALSE ;
int pial_surf_loaded = FALSE ;
int secondcurvloaded = FALSE;
int doublebufferflag = FALSE;
int openglwindowflag = FALSE;
int blinkflag = FALSE;
int renderoffscreen = FALSE;
int renderoffscreen1 = FALSE;
int blackcursorflag = FALSE;
int bigcursorflag = FALSE;
int vrml2flag = FALSE;
int showorigcoordsflag = FALSE;
int scrsaveflag = FALSE;
int phasecontourmodflag = FALSE;
int promptflag = FALSE;
int followglwinflag = TRUE;
int initpositiondoneflag = FALSE;
/* begin rkt */
int curvflag = TRUE; /* draw curv if loaded */
int mouseoverflag = TRUE; /* show mouseover information */
int redrawlockflag = FALSE; /* redraw on window uncover events */
int simpledrawmodeflag = TRUE; /* draw based on scalar values */
int forcegraycurvatureflag = FALSE; /* always draw grayscale curvature */
int drawcursorflag = TRUE; /* draw the cyan cursor */
int ignorezeroesinhistogramflag = FALSE; /* if true, don't count 0s in
overlay histogram */
int labels_before_overlay_flag = FALSE; /* draw labels under overlay or not */

Tcl_Interp *g_interp = NULL;

int curwindowleft = 0; /* keep track of window position, updated on move */
int curwindowbottom = 0;
int dontloadspherereg = FALSE; /* if true, don't try loading sphere.reg */
int scriptok=FALSE; /* global flag for signifying to parse a script */
char script_tcl[NAME_LENGTH]; /* name of the script to run */

/* end rkt */

int blinkdelay = BLINK_DELAY;
int blinktime = BLINK_TIME;

float contour_spacing[3];
double phasecontour_min = 0.3;
double phasecontour_max = 0.35;
int phasecontour_bright = 255;

int drawmask = 1;  /* same as next 6 */
int surfcolor = CURVATURE_OR_SULCUS;
int shearvecflag = FALSE;
int normvecflag = FALSE;
int movevecflag = FALSE;
int project = NONE;  /* flatten */

int computed_shear_flag = FALSE;

double update = 0.9;
double decay = 0.9;
float cweight = 0.33;
float aweight = 0.33;
float dweight = 0.33;
float epsilon = 0.000001;

double dip_spacing = 1.0;

int nrip = 0;
double stressthresh = 1e10;
float maxstress = 1e10;
float avgstress = 0;
int senstype = 0;
double fthresh = 2.0;
double fthreshmax = 5.0;
double tksfmax;
double wt = 0.5;
double wa = 0.5;
double ws = 0.5;
double wn = 0.5;
double wc = 0.0;
double wsh = 0.5;
double wbn = 0.5;
double ncthresh = 0.02;

int fixed_border = FALSE;
float fillscale = 1.5;
float curvscale;
double angle_offset = 0.0;
double angle_cycles = 1.0;
int smooth_cycles = 5;
int colscale = HEAT_SCALE;   /* changed by BRF - why was it COLOR_WHEEL??? */
double blufact = 1.0;
double cvfact = 1.5;
double fadef = 0.7;


float normal1[3] =
  {
    1.0,0.0,0.0
  };
float normal2[3] =
  {
    0.0,1.0,0.0
  };
float xpos[3] =
  {
    1.0,0.0,0.0
  };
float xneg[3] =
  {
    -1.0,0.0,0.0
  };
float ypos[3] =
  {
    0.0,1.0,0.0
  };
float yneg[3] =
  {
    0.0,-1.0,0.0
  };
float zpos[3] =
  {
    0.0,0.0,1.0
  };
float zneg[3] =
  {
    0.0,0.0,-1.0
  };
float v1[3],v2[3],v3[3],v4[3];
float v1[3],v2[3],v3[3],v4[3];

float gradavg[NLABELS][CMFBINS],gradsum[NLABELS][CMFBINS];
float gradnum[NLABELS][CMFBINS],gradtot[NLABELS];
float val1avg[NLABELS][CMFBINS],val2avg[NLABELS][CMFBINS];
float angleavg[NLABELS][CMFBINS];
float valnum[NLABELS][CMFBINS];
float gradx_avg[NLABELS],grady_avg[NLABELS],valtot[NLABELS],radmin[NLABELS];

double dipavg,dipvar,logaratavg,logaratvar,logshearavg,logshearvar;
double dipavg2;

#if 0
float xlo,xhi,ylo,yhi,zlo,zhi,xctr,yctr,zctr;
float xctr2,yctr2,zctr2;
#endif

#define PLOT_XDIM 200
#define PLOT_YDIM 200
#define PLOT_XFUDGE 7
#define PLOT_YFUDGE 30
#define MAX_RECFILES 20
#define MAX_NPLOTLIST 100

int *sol_dipindex,sol_neeg=0,sol_nmeg=0,sol_nchan=0,sol_ndipfiles=0,sol_ndip;
FLOATTYPE **sol_W,**sol_A,**sol_M,**sol_Mi;
FLOATTYPE **sol_Data[MAX_RECFILES],**sol_NoiseCovariance;
FLOATTYPE *sol_pval,*sol_prior,*sol_sensval,*sol_sensvec1,*sol_sensvec2;
FLOATTYPE *sol_lat,**sol_dipcmp_val[MAX_RECFILES],*sol_dipfact;
int sol_plotlist[MAX_NPLOTLIST],sol_ndec=0,sol_rectpval=0;
int sol_ntime,sol_ptime,sol_nnz,sol_nperdip,sol_plot_type=1,sol_proj_type=0;
int sol_allocated=FALSE,sol_nrec=0,sol_nplotlist=0,sol_trendflag=TRUE;
double sol_sample_period,sol_epoch_begin_lat,sol_snr_rms=10;
double sol_baseline_period=100,sol_baseline_end=0,sol_lat0=0,sol_lat1=100000;
double sol_loflim=2.0,sol_hiflim=10.0;
double sol_pthresh=0.0,sol_pslope=0.0,sol_maxrat=10.0;
int vertex_nplotlist=0,vertex_plotlist[MAX_NPLOTLIST];
int LeftRightRev = 0;

#define LIGHT0_BR  0.4 /* was 0.2 */
#define LIGHT1_BR  0.0
#define LIGHT2_BR  0.6 /* was 0.3 */
#define LIGHT3_BR  0.2 /* was 0.1 */
#define OFFSET 0.35   /* was 0.15, then was .25 (changed to .4 by BRF) */
#define BACKGROUND 0x00000000
double offset = OFFSET;
double light0_br,light1_br,light2_br,light3_br;

float idmat[4][4] =
  {
    {
      1.0,0.0,0.0,0.0
    }
    ,  /* Matrix idmat = */
    {0.0,1.0,0.0,0.0},
    {0.0,0.0,1.0,0.0},
    {0.0,0.0,0.0,1.0}
  };

/* accumulate really_ tranforms here */
float reallymat[4][4] =
  {
    {
      1.0,0.0,0.0,0.0
    }
    ,   /* Matrix reallymat = */
    {0.0,1.0,0.0,0.0},
    {0.0,0.0,1.0,0.0},
    {0.0,0.0,0.0,1.0}
  };

/* Talairach stuff */
LINEAR_TRANSFORM_ARRAY  *lta  = NULL;
int               transform_loaded = 0 ;

/* parcellation stuff */
static char parc_red[256] ;
static char parc_green[256] ;
static char parc_blue[256] ;
static char *parc_names[256] ;
static int parc_flag = 0 ;

static double dlat = 0 ;
static int surface_compiled = -1 ;
static int use_display_lists = 0 ;
static int FS_Brain_List = 1;
static int vertex_array_dirty = 0;
static int use_vertex_arrays = 1;

static int color_scale_changed = TRUE;
static char tmpstr[STRLEN];

/*---------------------- PROTOTYPES (should be static) ----------------*/
#ifdef TCL
int do_one_gl_event(Tcl_Interp *interp) ;
#endif
void compute_CMF(void) ;
void fix_nonzero_vals(void) ;
void smooth_momentum(int niter) ;
void smooth_logarat(int niter) ;
void smooth_shear(int niter) ;
void smooth_boundary_normals(int niter) ;
void scaledist(float sf) ;
float rtanh(float x) ;
void shrink(int niter) ;
void taubin_smooth(MRI_SURFACE *mris, int niter, double mu, double lambda) ;
void curv_shrink_to_fill(int niter) ;
void shrink_to_fill(int niter) ;
void transform(float *xptr, float *yptr, float *zptr, float nx, float ny,
               float nz, float d) ;
void really_translate_brain(float x, float y, float z) ;
void really_scale_brain(float x, float y, float z) ;
void really_rotate_brain(float a, char axis) ;
void align_sphere(MRI_SURFACE *mris) ;
void really_align_brain(void) ;  /* trans cent first -> cent targ */
void really_center_brain(void) ;
void really_center_second_brain(void) ;
void print_real_transform_matrix(void) ;
void write_really_matrix(char *dir) ;
void read_really_matrix(char *dir) ;
void flatten(char *dir) ;
void area_shrink(int niter) ; /* consider area */
void shrink2d(int niter) ;
void sphere_shrink(int niter, float rad) ;
void ellipsoid_project(float a, float b, float c) ;
void ellipsoid_morph(int niter, float a, float b, float c) ;
void ellipsoid_shrink(int niter, float a, float b, float c) ;
void ellipsoid_shrink_bug(int niter, float rad, float len) ;
void compute_curvature(void) ;
void clear_curvature(void) ;
void normalize_area(void) ;
void normalize_curvature(int which_norm) ;
void normalize_surface(void) ;
void load_brain_coords(float x, float y, float z, float v[]) ;
int outside(float x,float y, float z) ;
void draw_surface(void) ;
void draw_ellipsoid_latlong(float a, float b, float c) ;
void draw_second_surface(void) ;
void draw_scalebar(void) ;
void draw_colscalebar(void) ;
void set_stat_color(float f, float *rp, float *gp, float *bp, float tmpoffset);
void set_positive_color(float f, float *rp, float *gp, float *bp,
                        float tmpoffset) ;
void set_signed_color(float f, float *rp, float *gp,float *bp,float tmpoffset);
void set_color(float val, float curv, int mode) ;
void set_complexval_color(float x, float y, float stat, float curv) ;
void draw_spokes(int option) ;
void set_vertex_color(float r, float th, int option) ;
void set_color_wheel(float a, float a_offset, float a_cycles, int mode,
                     int logmode, float fscale) ;

int dngheat(float f, float *r, float *g, float *b);
int dngcolorwheel(float f, float *r, float *g, float *b);
int UseNewOverlay = 0;
static void fill_color_array2(MRI_SURFACE *mris, float *colors);
int LabelColor(int vno, float* r, float* g, float* b);

void restore_ripflags(int mode) ;
void dilate_ripped(void) ;
void rip_unmarked_vertices(void) ;
void floodfill_marked_patch(int filltype) ;
/* begin rkt */
/* Replacements for the floodfill_marked_patch stuff. */
int rip_all_vertices_except_contiguous_upripped ();
int mark_contiguous_vertices_with_similar_curvature();
int mark_contiguous_vertices_over_thresh();
/* end rkt */
void clear_ripflags(void) ;
void clear_vals(void) ;
void cut_marked_vertices(int closedcurveflag) ;
void cut_plane(void) ;
void cut_vertex(void) ;
void cut_line(int closedcurveflag) ;
void plot_curv(int closedcurveflag) ;
void draw_fundus(int bdry_index) ;
void plot_marked(char *fname) ;
void put_retinotopy_stats_in_vals(void) ;
void draw_vector(char *fname) ;
void clear_vertex_marks(void) ;
void clear_all_vertex_marks(void) ;
/* begin rkt */
void find_closest_marked_vertex (int screen_x, int screen_y,
                                 int* closest_index, int* closest_vno );

/* end rkt */
void mark_vertex(int vindex, int onoroff) ;
void mark_translated_vertex(int vindex, int onoroff, char *sphere_surf) ;
void mark_face(int fno) ;
void mark_annotation(int selection) ;
void mark_faces(int vno) ;
void prompt_for_parameters(void) ;
void prompt_for_drawmask(void) ;
void fill_triangle(float x0,float y0,float z0,float x1,float y1,float z1,
                   float x2,float y2,float z2) ;
void fill_surface(void) ;
void fill_second_surface(void) ;
void resize_buffers(long frame_xdim, long frame_ydim) ;
void open_window(char *name) ;
void swap_buffers(void);
void to_single_buffer(void) ;
void to_double_buffer(void) ;
void do_lighting_model(float lite0,float lite1,float lite2,float lite3,
                       float newoffset) ;
void make_filenames(char *lsubjectsdir,char *lsrname,char *lpname,char *lstem,
                    char *lext) ;
void print_help_tksurfer(void) ;
int Surfer(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[]);
int main(int argc, char *argv[]) ;
double dipval(int cond_num, int nzdip_num, int ilat) ;
void read_ncov(char *fname) ;
void regularize_matrix(FLOATTYPE **m, int n, float fact) ;
void read_iop(char *fname, int dipfilenum) ;
void read_rec(char *fname) ;
void filter_recs(void) ;
void filter_rec(int np) ;
void remove_trend(FLOATTYPE *v, int n) ;
void compute_timecourses(void) ;
void compute_timecourse(int num) ;
void plot_all_time_courses(void) ;
void read_plot_list(char *fname) ;
void read_vertex_list(char *fname) ;
void plot_nearest_nonzero(int vnum) ;
void find_nearest_nonzero(int vnum, int *indxptr) ;
void find_nearest_fixed(int vnum, int *vnumptr) ;
void plot_sol_time_course(int imin, int plot_xdim, int plot_ydim,
                          int plot_xloc, int plot_yloc, int plotnum) ;
void load_vals_from_sol(float tmid, float dt, int num) ;
void load_var_from_sol(int num) ;
void variance_ratio(void) ;
void load_pvals(void) ;
void compute_pval_fwd(float pthresh) ;
void compute_select_fwd(float maxdist) ;
void compute_select_crosstalk(void) ;
void compute_all_crosstalk(int weighttype) ;
void compute_select_pointspread(void) ;
void compute_all_pointspread(void) ;
void recompute_select_inverse(void) ;
void compute_pval_inv(void) ;
void normalize_vals(void);
void scale_vals(float scale) ;
void normalize_time_courses(int normtype) ;
void normalize_inverse(void) ;
void setsize_window(int pix);
void resize_window(int pix) ;
void save_rgb(char *fname) ;
void save_rgb_num(char *dir) ;
void save_rgb_named_orig(char *dir, char *name) ;
void scrsave_to_rgb(char *fname) ;
void pix_to_rgb(char *fname) ;
void read_annotated_image(char *fpref, int frame_xdim, int frame_ydim) ;
void read_annotations(char *fname) ;
int diff(int a,int b) ;
void save_rgb_cmp_frame(char *dir, int ilat) ;
void open_rgb_cmp_named(char *dir, char *name) ;
void save_rgb_cmp_frame_named(float lat) ;
void do_rgb_cmp_frame(long xsize,long ysize, FILE *fp) ;
void close_rgb_cmp_named(void) ;
void redraw(void) ;
void twocond(int c0, int c1) ;
void redraw_second(void) ;
void blinkbuffers(void) ;
void redraw_overlay(void) ;
void draw_cursor(int vindex,int onoroff ) ;
void draw_all_cursor(void) ;
void draw_all_vertex_cursor(void) ;
void clear_all_vertex_cursor(void) ;
void select_vertex(short sx,short sy) ;
void select_vertex_by_vno(int vno) ;
void find_vertex_at_screen_point(short sx,short sy,int* ovno, float* od) ;
void invert_vertex(int vno) ;
void invert_face(int fno) ;
void orient_sphere(void) ;
void dump_vertex(int vno) ;
void dump_faces(int vno) ;
static void load_gcsa(char *fname) ;
void left_click(short sx,short sy) ;
void sample_annotated_image(void) ;
void restore_zero_position(void) ;
void restore_initial_position(void) ;
void make_lateral_view(char *stem) ;
void make_medial_view(char *stem) ;
void make_lateral_view_second(char *stem) ;
void write_val_histogram(float min, float max, int nbins) ;
void write_view_matrix(char *dir) ;
void read_view_matrix(char *dir) ;
void translate_brain(float x, float y, float z) ;
void scale_brain(float s) ;
void rotate_brain(float a, char c) ;
void read_image_info(char *fpref) ;
void read_talairach(char *fname) ;
void read_images(char *fpref) ;
void alloc_curv_images(void) ;
void read_curv_images(char *fpref) ;
void curv_to_curvim(void) ;
void second_surface_curv_to_curvim(void) ;
void swap_curv(void);
void curvim_to_surface(void) ;
void curvim_to_second_surface(void) ;
void smooth_curvim(int window) ;
void add_subject_to_average_curvim(char *name, char *morphsubdir)   ;
void smooth_curvim_sparse(int niter) ;
unsigned char floattobyte(float f, float min, float max) ;
float bytetofloat(unsigned char c, float min, float max) ;
void write_images(unsigned char ***mat,char *fpref) ;
int  read_binary_surface(char *fname) ;
/* begin rkt */
int vset_read_vertex_set(int set, char* fname ) ;
/* end rkt */
void read_positions(char *name) ;
void save_surf(void) ;
void restore_surf(void) ;
void read_second_binary_surface(char *fname) ;
void read_second_binary_curvature(char *fname) ;
void normalize_second_binary_curvature(void) ;
void show_surf(char *surf_name) ;
int  read_orig_vertex_coordinates(char *fname) ;
int  read_inflated_vertex_coordinates(void) ;
int  read_white_vertex_coordinates(void) ;
int read_pial_vertex_coordinates(void) ;
int read_canon_vertex_coordinates(char *fname) ;
void send_spherical_point(char *subject_name,char *canon_name,
                          char *orig_fname);
void send_contralateral_point(char *canon_name, char *orig_name) ;
void send_to_subject(char *subject_name) ;
void send_to_other_hemi(void) ;
static void resend_to_subject(void) ;
void invert_surface(void) ;
void read_ellipsoid_vertex_coordinates(char *fname,float a,float b,float c) ;
void find_orig_vertex_coordinates(int vindex) ;
void select_talairach_point(int *vindex,float x_tal,float y_tal,float z_tal);
void select_orig_vertex_coordinates(int *vindex) ;
void print_nearest_vertex_to_talairach_point(float x_tal, float y_tal, float z_tal);
void read_curvim_at_vertex(int vindex) ;
int  write_binary_surface(char *fname) ;
void write_binary_patch(char *fname) ;
void read_binary_patch(char *fname) ;
void write_labeled_vertices(char *fname) ;
void read_and_color_labeled_vertices(int r, int g, int b) ;
void read_labeled_vertices(char *fname) ;
void write_binary_curvature(char *fname) ;
void write_binary_areas(char *fname) ;
void write_binary_values(char *fname) ;
/* begin rkt */
void read_binary_values(char *fname) ;
void read_binary_values_frame(char *fname) ;
/* end rkt */
void swap_stat_val(void) ;
void swap_val_val2(void) ;
void shift_values(void) ;
void swap_values(void) ;
void read_binary_decimation(char *fname) ;
void write_binary_decimation(char *fname) ;
void write_decimation(char *fname) ;
void read_binary_dipoles(char *fname) ;
void write_binary_dipoles(char *fname) ;
void write_dipoles(char *fname) ;
void write_subsample(char *fname) ;
void write_binary_subsample(char *fname) ;
void read_binary_subsample(char *fname) ;
void read_binary_curvature(char *fname) ;
void normalize_binary_curvature(void) ;
void read_binary_areas(char *fname) ;
void read_fieldsign(char *fname);
void write_fieldsign(char *fname) ;
void read_fsmask(char *fname) ;
void write_fsmask(char *fname) ;
void write_vrml(char *fname,int mode) ;
void rip_faces(void) ;
void normalize(float v[3]) ;
void normal_face(int fac,int n,float *norm);
float triangle_area(int fac,int n) ;
void find_neighbors(void) ;
void compute_normals(void) ;
void compute_shear(void) ;
void double_angle(float x,float y,float *dx,float *dy) ;
void halve_angle(float x,float y,float *dx,float *dy);
void compute_boundary_normals(void) ;
void subsample_dist(int spacing) ;
void subsample_orient(float spacing) ;
void smooth_curv(int niter) ;
void smooth_val_sparse(int niter) ;
void smooth_val(int niter) ;
void smooth_fs(int niter) ;
void fatten_border(int niter) ;
void compute_angles(void) ;
float circsubtract(float a,float b) ;
void compute_fieldsign(void) ;
void compute_cortical_thickness(void) ;
static void read_disc(char *subject_name) ;
static void deconvolve_weights(char *weight_fname, char *scale_fname) ;
/* begin rkt */
static int init_vertex_arrays(MRI_SURFACE *mris) ;

void flip_normals (char *axes);

void swap_vertex_fields(int a, int b); /* params are FIELD_*  */
static void print_vertex_data(int vno, FILE *fp, float dmin) ;
static void send_current_labels() ;
static void update_labels(int label_set, int vno, float dmin) ; /* LABELSET_ */

/* set the rip value of a face or vertex. */
void set_face_rip(int fno, int rip, int undoable);
void set_vertex_rip(int vno, int rip, int undoable);

/* adds another marked vertex to close the loop. */
void close_marked_vertices ();

/* draws a crosshair at a vertex in the current color. */
void draw_vertex_hilite (int vno);

/* draws all marked verts with white hilites, using draw_vertex_hilite(). */
void draw_marked_vertices ();
/* end rkt */

void LoadMRISMask(void);

/* krish -- if linktimepoint mode is enabled */
void link_timepoint_ROI(int vno);

/* external prototypes */
void buffer_to_image(unsigned char *buf,unsigned char**im,int ysize,int xsize);
void image_to_buffer(unsigned char **im,unsigned char*buf,int ysize,int xsize);
void file_name(char *fpref,char *fname,int num,char *form) ;
void inverse(FLOATTYPE **a,FLOATTYPE **y, int n) ;
/* void vector_multiply(FLOATTYPE **a, FLOATTYPE *b,
   FLOATTYPE *c, int n, int m); */
void bpfilter(FLOATTYPE **data, int nchan, int nsamp,float lo,float hi);


/* begin rkt */

/* -------------------------------------------------- the window and events */

#ifdef USE_XGLUT_WINDOW

#include "xGLutWindow.h"

xGLutWindowRef gWindow = NULL;

/* This was experimental code to make the main window be a GLut window
   instead of a straight XWindow with Motif event handling. It wasn't
   compatible with all systems. */
void wndw_create (int x, int y, int width, int height);
void wndw_set_title (char* title);
void wndw_handle_event (void* data, xGWin_tEventRef event);

#endif

/* ------------------------------------------------------------------------ */

/* -------------------------------------------------- ctrl-c cancel support */

int cncl_listening = 0;
int cncl_canceled = 0;

/* This code lets you trap a ctrl-c during a long operation in case
   the user wants to cancel it. Call cncl_start_listening() at the
   beginning of your operation, and cncl_user_canceled() to check if
   the user has hit ctrl-c. When you're done, call
   cncl_stop_listening(). cncl_initialize() installs
   cncl_handle_sigint() as the handler function for ctrl-c, and it
   will set cncl_canceled to 1 if it receives that signal. */
void cncl_initialize ();
void cncl_start_listening ();
void cncl_stop_listening ();
int cncl_user_canceled ();

void cncl_handle_sigint (int signal);

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------- menu set support */

/* menu sets */
#define MENUSET_VSET_INFLATED_LOADED   0
#define MENUSET_VSET_WHITE_LOADED      1
#define MENUSET_VSET_PIAL_LOADED       2
#define MENUSET_VSET_ORIGINAL_LOADED   3
#define MENUSET_TIMECOURSE_LOADED      4
#define MENUSET_OVERLAY_LOADED         5
#define MENUSET_CURVATURE_LOADED       6
#define MENUSET_LABEL_LOADED           7
#define MENUSET_FIELDSIGN_LOADED       8
#define MENUSET_FIELDMASK_LOADED       9

int enable_menu_set (int set, int enable);

/* ------------------------------------------------------------------------ */

/* -------------------------------------------- multiple vertex set support */

/*
 * this code lets the user load in alternate vertex sets and switch between
 * them. it doesn't use the vertices in the MRIS strucutre. each set is
 * loaded in with MRISread but then the vertices are copied into external
 * storage, leaaving the MRIS unchanged (except for the tmp vertices swap).
 * the main vertices are also copied into external storage. when the user
 * wants to display a different set, the vertices are copied from external
 * storage into the MRIS and the normals are recomputed.
 */

/* vertex sets */
#define VSET_MAIN      0
#define VSET_INFLATED  1
#define VSET_WHITE     2
#define VSET_PIAL      3
#define VSET_ORIGINAL  4
#define NUM_VERTEX_SETS  5

/* the current set */
static int vset_current_set = VSET_MAIN ;

/* storage for other sets, a node and a list. the list is allocated
   dynamically, getting the number of verts from the MRIS. */
typedef struct
{
  float x,y,z;
}
VSET_VERTEX;

VSET_VERTEX* vset_vertex_list[NUM_VERTEX_SETS];

/* call at startup. sets the storage pointers to null. */
int vset_initialize ();

/* reads in a vertex set from a file. actually copies the MRIS main verts
   into tmp space in the MRIS, loads the surface, copies the verts into
   external storage, and restores the MRIS verts from the tmp storage. */
int vset_read_vertex_set(int set, char* fname);

/* copies and restores main verts from MRIS into and out of
   external storage. */
int vset_load_surface_vertices(int set);
int vset_save_surface_vertices(int set) ;

/* changes the current set, calling vset_load_surface_vertices and recomputes
   the normals. marks the vertex array dirty for the next redraw. */
int vset_set_current_set(int set);

/* ------------------------------------------------------------------------ */

/* -------------------------------------------------- coordinate conversion */

MATRIX* conv_mnital_to_tal_m_ltz = NULL;
MATRIX* conv_mnital_to_tal_m_gtz = NULL;
MATRIX* conv_tal_to_mnital_m_ltz = NULL;
MATRIX* conv_tal_to_mnital_m_gtz = NULL;
MATRIX* conv_tmp1_m = NULL;
MATRIX* conv_tmp2_m = NULL;
MATRIX* surfaceRAStoRAS = NULL;

MRI* orig_mri_header = NULL;

/* Handles conversion to and from Talairach coordinates. */
int conv_initialize ();
int conv_ras_to_mnital (float rasx, float rasy, float rasz,
                        float* talx, float* taly, float* talz);
int conv_ras_to_tal    (float rasx, float rasy, float rasz,
                        float* talx, float* taly, float* talz);
int conv_mnital_to_ras (float talx, float taly, float talz,
                        float* rasx, float* rasy, float* rasz);
int conv_tal_to_ras    (float talx, float taly, float talz,
                        float* rasx, float* rasy, float* rasz);

/* ------------------------------------------------------------------------ */

/* ----------------------------------------------------------- undo support */

#include "xGrowableArray.h"

/*
 * an undo action (UNDO_ACTION) represents a single logical action
 * that can be undone in user langugage, such as a cut that effected
 * many verices and faces. an action has a type and a list of action
 * nodes. the type represents that type of action that can be undone
 * (i.e. UNDO_CUT), and has an associated node type
 * (i.e. UNDO_CUT_NODE).
 *
 * to start an undoable action, call undo_begin_action to allocate the
 * action. then pass pointers to the action node to undo_copy_action_node,
 * in which it will be copied to the list. when done adding nodes, call
 * undo_finish_action.
 *
 * use the undo_get_action_string to get a human readable description
 * of the action to be undo. call undo_do_first_action to undo the
 * last action. multiple undos can be done, up to NUM_UNDOS. undoing
 * an action removes the action from the undo list - it's not
 * 'redoable'.
 */

/* types of undo actions */
#define UNDO_INVALID      0
#define UNDO_NONE         1
#define UNDO_CUT          2
#define NUM_UNDO_ACTIONS  3

char *undo_action_strings[NUM_UNDO_ACTIONS] =
  {
    "INVALID UNDO ACTION",
    "Nothing to Undo",
    "Undo Cut"
  };

/* represents an undoable action, a list of action nodes of specific types */
typedef struct
{
  int undo_type; /* UNDO_* */
  xGrowableArrayRef node_list;
}
UNDO_ACTION;

/* storage of undo actions. when the list is open, undo_list[0] is the list
   being created and added to. when closed, undo_list[0] is the first
   action that will be undone. */
#define UNDO_LIST_POS_FIRST 0
#define UNDO_LIST_POS_LAST  3
#define NUM_UNDO_LISTS      4
static UNDO_ACTION* undo_list[NUM_UNDO_LISTS];

/* initialize anything that needs to get initialized */
int undo_initialize ();

/* gets the size of the node associated with this action type. returns -1
   if the action type is invalid. */
int undo_get_action_node_size(int action_type); /* UNDO_* */

/* return the string that describes the action available to be undone.
   returns a null string if the action type is invalid. */
char* undo_get_action_string(int action_type); /* UNDO_* */

/* creates a new undo list and inserts it in the undo_list array. can possibly
   remove and delete an action if the list is full. add_undoable_action_node
   can be called to add individual action nodes. returns en error if the
   action is already begun. */
int undo_begin_action(int action_type);

/* finishes the list creation. */
int undo_finish_action();

/* sends a message to tcl setting the undo menu item. */
int undo_send_first_action_name();

/* undo_begin_action sets to open, undo_finish_action sets to
   closed. shouldn't be set directly. */
#define UNDO_LIST_STATE_CLOSED    0
#define UNDO_LIST_STATE_OPEN      1
static int undo_list_state=UNDO_LIST_STATE_CLOSED;

/* copies the node tothe list being created. returns an error if the list
   is closed. */
int undo_copy_action_node(void* node);

/* undoes the first action if the list is closed. returns an error if the
   list is open. */
int undo_do_first_action();

/* UNDO_CUT */
/* type of cut node */
#define UNDO_CUT_VERTEX 0
#define UNDO_CUT_FACE   1

typedef struct
{
  int cut_type;   /* UNDO_CUT_* */
  int index;      /* face or vertex */
  int rip_value;  /* the value that will be restored */
}
UNDO_CUT_NODE;

/* creates an undo cut node add and copies it to the current open list. */
int undo_new_action_cut(int cut_type, int index, int rip_value);

/* undoes an action node, i.e. sets the referenced face's or vertex's
   rip value to the value in the node. */
int undo_do_action_cut(UNDO_ACTION* action);

/* ---------------------------------------------------------------------- */

/* -------------------------------------------- functional volume support */

#include "mriFunctionalDataAccess.h"
#include "xUtilities.h"

/* we keep a separate list of  */
typedef struct
{
  float x, y, z;
  int vno;
}
FUNC_SELECTED_VOXEL;

static mriFunctionalDataRef func_timecourse = NULL;
static mriFunctionalDataRef func_timecourse_offset = NULL;
static int func_use_timecourse_offset = FALSE;
static int func_sub_prestim_avg = FALSE;
static xGrowableArrayRef func_selected_ras = NULL;
static tBoolean func_is_scalar_volume = TRUE;

#define knLengthOfGraphDataItem              18 // for "100.1 1000.12345 "
#define knLengthOfGraphDataHeader            20 // for cmd name + cond + {}
#define knMaxCommandLength                   50

/* tcl vars */
static double func_time_resolution = 0;
static int func_num_prestim_points = 0;
static int func_num_conditions = 0;
static int func_num_timepoints = 0;

/* Controls whether we are graphing the current selected vertex or the
   average of the marked vertices or the average of the label. */
#define FUNC_GRAPH_AVG_MODE_SELECTED 0
#define FUNC_GRAPH_AVG_MODE_MARKED 1
#define FUNC_GRAPH_AVG_MODE_LABEL 2
static int func_graph_avg_mode = FUNC_GRAPH_AVG_MODE_SELECTED;

int func_initialize ();

int func_load_timecourse (char* fname, FunD_tRegistrationType reg_type,
                          char* registration);
int func_load_timecourse_offset (char* fname, FunD_tRegistrationType reg_type,
                                 char* registration);

/* Adds voxels to the functional vertex selection. The first adds only
   the selected vertex, the second adds all marked vertices, and the
   third adds all the verts in the current label. */
int func_select_selected_vertex ();
int func_select_marked_vertices ();
int func_select_label ();

/* Clear the selection. */
int func_clear_selection();

/* Add a voxel to the selection. */
int func_select_voxel (int vno, float x, float y, float z);

/* Graph the currently selected vertices. If > 1, will graph the
   average. */
int func_graph_timecourse_selection ();

/* Calcs the avg values and deviations of the selected voxels. */
int func_calc_avg_timecourse_values (int condition, int* num_good_voxels,
                                     float values[], float deviations[] );

/* Normalizes the values in the time course volume. */
int func_normalize ();

/* Calcs the correlation of the specified vertex and the time
   course through all time points and conditions. It writes the
   results to an overlay layer. */
int func_calc_correlation_and_write_to_overlay (int vno, int field);

/* Prints summary data about the currently selected verts to a
   file. */
int func_print_timecourse_selection (char* filename);

int func_convert_error (FunD_tErr error);

/* ---------------------------------------------------------------------- */

/* --------------------------------------------------- scalar value mgmnt */


#if 0
#define SCLV_FIELDSIGN    6
#define SCLV_FSMASK       7
#endif

static char *sclv_field_names [NUM_SCALAR_VALUES] =
  {
    "val", "val2", "valbak", "val2bak", "valstat", "imagval",
    "mean", "meanimag", "std_error"
  };

typedef struct
{
  int is_functional_volume; /* use func_volume */
  tBoolean is_scalar_volume;     /* use vno,0,0 as func index */
  int cur_timepoint;
  int cur_condition;
  int num_timepoints; /* always 1 for .w files */
  int num_conditions; /* always 1 for .w files */
  mriFunctionalDataRef func_volume;
  float fthresh;
  float fmid;
  float fslope;
  float foffset;
  float min_value;
  float max_value;

  int num_freq_bins;
  int ***frequencies; /* the frequency of values in num_freq_bins for
                             each time point and condition i.e.
                             frequency[cond][tp][bin] */
  int num_zeroes_in_zero_bin;   /* This contains the number of 0s in the 0 */
  int zero_bin_index;           /* bin. We do this separately so we
                                       can turn off the 0 displa  easily. */

}
SCLV_FIELD_INFO;

static int sclv_current_field = SCLV_VAL;
static SCLV_FIELD_INFO sclv_field_info[NUM_SCALAR_VALUES];
static double sclv_value_min = 0;
static double sclv_value_max = 0;
static int sclv_num_timepoints = 0;
static int sclv_num_conditions = 0;
static int sclv_cur_timepoint = 0;
static int sclv_cur_condition = 0;

static double sclv_overlay_alpha = 1.0;
static int sclv_opaque = 1; /* If on, min->mid will be a solid color */

#define sclv_set_value(v,i,n) \
 switch((i)) { \
     case SCLV_VAL: (v)->val = (n); break; \
     case SCLV_VAL2: (v)->val2 = (n); break; \
     case SCLV_VALBAK: (v)->valbak = (n); break; \
     case SCLV_VAL2BAK: (v)->val2bak = (n); break; \
     case SCLV_VALSTAT: (v)->stat = (n); break; \
     case SCLV_IMAG_VAL: (v)->imag_val = (n); break; \
     case SCLV_MEAN: (v)->mean = (n); break; \
     case SCLV_MEAN_IMAG: (v)->mean_imag = (n); break; \
     case SCLV_STD_ERROR: (v)->std_error = (n); break; \
 }

#define sclv_get_value(v,i,n) \
 switch(i) { \
     case SCLV_VAL: (*n) = (v)->val; break; \
     case SCLV_VAL2: (*n) = (v)->val2; break; \
     case SCLV_VALBAK: (*n) = (v)->valbak; break; \
     case SCLV_VAL2BAK: (*n) = (v)->val2bak; break; \
     case SCLV_VALSTAT: (*n) = (v)->stat; break; \
     case SCLV_IMAG_VAL: (*n) = (v)->imag_val; break; \
     case SCLV_MEAN: (*n) = (v)->mean ; break; \
     case SCLV_MEAN_IMAG: (*n) = (v)->mean_imag; break; \
     case SCLV_STD_ERROR: (*n) = (v)->std_error; break; \
 }

int sclv_initialize ();

int sclv_unload_field (int field);

int sclv_read_from_dotw (char* fname, int field);
int sclv_read_from_dotw_frame (char* fname, int field);

int sclv_read_from_volume (char* fname, FunD_tRegistrationType reg_type,
                           char* registration, int field);

int sclv_read_from_annotcorr (char* fname, int field);

/* Creates a new overlay and fills it with the stat values from a
   label. */
int sclv_new_from_label (int field, int label);

/* Initializes an empty field. Pass a name for the field or NULL. */
int sclv_new_empty (int field, char* name);

/* writes .w files only */
int sclv_write_dotw (char* fname, int field);

/* make a pass through all the conditions and timepoints and calculate
   frequencies for everything. */
#define SCLV_NUM_FREQUENCY_BINS 250
int sclv_calc_frequencies (int field);

/* generic version of smooth_val, smooth_fs, etc */
int sclv_smooth (int niter, int field);

/* sets the overlay alpha value. */
int sclv_set_overlay_alpha (double alpha);

/* Sets the current field that will be drawn to the surface. */
int sclv_set_current_field (int field);

int sclv_send_current_field_info ();

/* Sets the timepoint of an overlay with a volume. this will actually
   paint the timepoint's data onto the proper sclv field in the
   volume.*/
int sclv_set_timepoint_of_field (int field, int timepoint, int condition);

/* Fills out an array of float nvertices long with the proper
   functional values for the given field and timepoint. */
int sclv_get_values_for_field_and_timepoint (int field, int timepoint,
    int condition, float* values);

/* this stuff is kind of experimental, it works but it's not very
   useful as the percentages have to be very precise. */
int sclv_set_current_threshold_from_percentile (float thresh, float mid,
    float max);
int sclv_set_threshold_from_percentile         (int field, float thresh,
    float mid, float max);
int sclv_get_value_for_percentile              (int field, float percentile,
    float* value);

/* sets the threshold using FDR. */
int sclv_set_threshold_using_fdr (int field, float rate, int only_marked);

/* copies field settings from one field to another */
int sclv_copy_view_settings_from_current_field (int field);
int sclv_copy_all_view_settings_from_current_field ();
int sclv_copy_view_settings_from_field (int field, int fromfield);

/* swaps values in two fields */
int sclv_swap_fields (int fielda, int fieldb);

/* sends a tcl command with the information necessary to build a
   histogram. */
int sclv_send_histogram ( int field );

/* gets a single 0-1 color value for a scalar value. used by tcl to
   color a histogram. */
int sclv_get_normalized_color_for_value (int field, float value,
    float *outRed,
    float *outGreen,
    float *outBlue);

/* Applies a functional value's heat scale color to the given color,
   overlaying them in the given opacity. */
int sclv_apply_color_for_value (float fval, float opacity,
                                GLubyte* r, GLubyte* g, GLubyte* b );

int sclv_load_label_value_file (char* fname, int field);

/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------ multiple labels */

int labl_debug = 0;

typedef enum {
  LABL_STYLE_FILLED = 0,
  LABL_STYLE_OUTLINE,
  LABL_NUM_STYLES
} LABL_DRAW_STYLE;

/* a label without an assigned structure. */
#define LABL_TYPE_FREE -1

/* the default color for free labels. */
#define LABL_DEFAULT_COLOR_R 100
#define LABL_DEFAULT_COLOR_G 100
#define LABL_DEFAULT_COLOR_B 200
#define LABL_DEFAULT_COLOR_A 255

typedef struct
{
  LABEL* label;
  int structure;    /* the structure assigned to this label, or
                           LABL_TYPE_FREE for a free color label */
  int r, g, b, a;   /* from color of structure in lookup table or if
                           LABL_TYPE_FREE, assigned by the user */
  int visible;
  char name[NAME_LENGTH];  /* name of this label (not neccessarily the
                                  same as the structure name, if
                                  assigned) */
  float min_x, max_x;    /* the bounding cube of this label. */
  float min_y, max_y;
  float min_z, max_z;

  int* border_vno;      /* a list of border vnos */
  int  num_border_vnos;

  int cached;   /* whether or not this label is in the cache */
}
LABL_LABEL;

/* a fixed array of labels. */
#define LABL_MAX_LABELS 30000
LABL_LABEL labl_labels[LABL_MAX_LABELS];
int labl_num_labels;

/* the color lookup table. */
COLOR_TABLE* labl_ctab = NULL;
int labl_num_structures;

/* the currently selected label */
#define LABL_NONE_SELECTED -1
int labl_selected_label = LABL_NONE_SELECTED;

/* a global count of labels crated for naming purposes. */
int labl_num_labels_created;

/* style in which to draw the labels. */
LABL_DRAW_STYLE labl_draw_style = LABL_STYLE_FILLED;

/* the selection outline color */
int labl_outline_color[3] = {255, 255, 0}; /* yellow default */

/* whether or not to draw labels. */
int labl_draw_flag = 1;

/* the fudge for the bounding cube for each label. */
#define LABL_FUDGE 4.0

/* name of the color table. */
char* labl_color_table_name;

/* cache of labels at vnos */
int   labl_cache_updated;
int*  labl_num_labels_at_cache_vno;
int** labl_cache;   // label_cache[vno][0..labl_num_labels_at_cache_vno[vno]]

/* Whether or not clicking selects labels. */
int labl_select_flag = 1;

/* initialize everything. */
int labl_initialize ();

/* makes and updates the label cache */
int labl_update_cache (int force_rebuild_all);

/* loads a new color table. if any labels had existing indicies that
   are now out of bounds, they are set to 0. sends an update to the
   tcl list with the structure names. */
int labl_load_color_table (char* fname);

/* sends the entries from either the surface's private color table or
   the external color table we loaded to tcl. */
int labl_send_color_table_info ();

/* NOTE: read_labeled_vertices now routed here. */
int labl_load (char* fname);
int labl_save (int index, char* fname);
int labl_save_all (char* prefix);

/* finds the border of a label and sets the border flags appropriately. */
int labl_find_and_set_all_borders ();
int labl_find_and_set_border (int index);

/* returns whether a vno is a border in a label. */
int labl_vno_is_border (int index, int vno);

/* reads an annotation file and makes multiple labels, one for each
   annotation. tries to associate a structure index by matching the
   color with a structure in the current lookup table. otherwise marks
   it 0. NOTE: read_annotations now routed here. */
int labl_import_annotation (char* fname);

/* saves all labels in one annotation file. if any labels have shared
   voxels, the last label in the list will overwrite any previous
   ones. */
int labl_export_annotation (char* fname);

/* makes the marked vertices a new label. */
int labl_new_from_marked_vertices ( int *new_index_out);

/* add marked verts to an exisiting label */
int labl_add_marked_vertices_to_label (int index);

/* remove marked verts from an exisiting label */
int labl_remove_marked_vertices_from_label (int index);

/* marks all the vertices in a label. */
int labl_mark_vertices (int index);

/* marks thresholded vertices in a label. */
int labl_mark_threshold_vertices (int index, float threshold);

/* selects a label, drawing it with apn outline and hilighting it in
   the label list tcl window. */
int labl_select (int index);

/* sets the label index name from the structure name in the color table */
int labl_set_name_from_table (int index);

/* changes information about a label. */
int labl_set_info (int index, char* name, int structure, int visible,
                   int r, int g, int b);
/* changes the color of a label. only valid if it is a free label. */
int labl_set_color (int index, int r, int g, int b);
int labl_set_alpha (int index, int a);

/* sends a label's information to tcl */
int labl_send_info (int index);

/* makes a new LABL_LABEL entry with this label data. note that this
   function uses the label you pass in, so don't delete it. passes
   back the index of the label. sends an update to the tcl list with
   the new name. */
int labl_add (LABEL* label, int* new_index);

/* call when a label's internal coordinates have changed. recalcs the
   extent, finds the borders, all that stuff. */
int labl_changed (int new_index, int vertices_were_removed);

/* Call when the vertex set has changed, as all label coords have
   changed. */
int labl_all_changed ();

/* removes and deletes the label and bumps down all other labels in
   the list. */
int labl_remove (int index);

/* Thresholds a label according to a user entered scalar value */
int labl_threshold (int index, float threshold);

/* removes and deletes all labels. */
int labl_remove_all ();

/* checks if this vno is in a label. passes back the label index. */
int labl_find_label_by_vno (int vno, int min_label,
                            int* index_array, int array_size, int* num_found);

/* figures out if a click is in a label. if so, selects it. */
int labl_select_label_by_vno (int vno);

/* if this vno is in a label, changes the color of this vertex
   accordingly. */
int labl_apply_color_to_vertex (int vno, GLubyte* r, GLubyte* g, GLubyte* b );

/* morphology operations. */
int labl_erode (int index);
int labl_dilate (int index);
int labl_fill_holes (int index);

/* prints the label list. */
int labl_print_list ();

/* prints the color table. */
int labl_print_table ();

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------- paths */

int path_debug = 0;

typedef struct
{
  int num_vertices;      /* number of vertices in path. */
  int* vertices;         /* array of vertices in path. */
  float min_x, max_x;    /* the bounding cube of this path. */
  float min_y, max_y;
  float min_z, max_z;

}
PATH_PATH;

/* array of paths. */
#define PATH_MAX_PATHS 200
PATH_PATH path_paths[PATH_MAX_PATHS];
int path_num_paths = 0;

/* currently selected path. */
#define PATH_NONE_SELECTED -1
int path_selected_path = PATH_NONE_SELECTED;

static PATH_PATH* path_copy(PATH_PATH *bsrc, PATH_PATH *bdst) ;
static int        path_set_marks(PATH_PATH *b, MRI_SURFACE *mris,int mark) ;

/* Find length of path. */
static double     path_sse(PATH_PATH *b, MRI_SURFACE *mris,
                           double target_curv, double l_curv, double l_len) ;
static double     path_length(PATH_PATH *b, MRI_SURFACE *mris) ;

#define PATH_FUDGE 4.0

/* Array of flags for each mris vertex, whether or not a surface
   vertex is in a path. */
char* path_is_path = NULL;

/* Initialize all the data structures for the path stuff. */
int path_initialize ();

/* Creates a new path between the currently marked verices. */
int path_new_path_from_marked_vertices ();

/* Removes the currently selected path, deleting it. */
int path_remove_selected_path ();

/* Creates and removes paths from the path_path array. */
int path_add (int num_vertices, int* vertices, int* new_index);
int path_remove (int index);

/* Update the bounding cube of all paths, necesssary when the vertex
   set swithces as the path now lies on different coords. */
int path_update_bounding_boxes ();

/* Mark the path as selected. If this changes the selected path, cause
   a redraw event. */
int path_select (int index);

/* Mark the vertices in the path. */
int path_mark_selected_path ();
int path_mark (int index);

/* Updates the path_is_path array after a path has been added or
   removed. */
int path_update_surface_paths ();

/* Returns true if this vertex is on a fill path. */
char path_is_vertex_on_path (int vno);

/* Figures out if a click is on or near a path. if so, selects it. */
#define PATH_DISTANCE_TO_SELECT 25
int path_select_path_by_vno (int vno);

/* If this vno is a path, changes the color of this vertex
   accordingly. */
int path_apply_color_to_vertex (int vno, GLubyte* r, GLubyte* g, GLubyte* b );

/* Writes and reads a path files. */
int path_save (char* fname);
int path_load (char* fname);

/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------- floodfill mark */

#define FILL_NO_ACTION_JUST_MARK 0
#define FILL_ACTION_NEW_LABEL    1
#define FILL_ACTION_ADD_LABEL    2
#define FILL_ACTION_REMOVE_LABEL 3
#define NUM_FILL_ACTIONS         4

typedef struct
{

  char dont_cross_path;
  char dont_cross_label;
  char dont_cross_cmid;
  char dont_cross_fthresh;
  char dont_fill_unlabeled;

  char use_multiple_seeds; /* If set, use all marked[] vnos. */

  int action;
  int argument;

  int new_label_index; /* set if the action is NEW_LABEL */

}
FILL_PARAMETERS;

/* Performs a flood fill from the parameter vertex, using the
   FILL_PARAMETERS to determine what and how the fill should do. The
   options are hopefully self-explanatory. */
int fill_flood_from_seed (int vno, FILL_PARAMETERS* params);

/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------- vertex editing */

#define VEDIT_NO_ACTION           0
#define VEDIT_ACTION_NEW_LABEL    1
#define VEDIT_ACTION_ADD_LABEL    2
#define VEDIT_ACTION_REMOVE_LABEL 3
#define VEDIT_ACTION_CLEAR_LABELS 4
#define NUM_VEDIT_ACTIONS         5

int edit_vertex_at_cursor ( int action, int argument );
int edit_vertex ( int vno, int action, int argument );

/* ---------------------------------------------------------------------- */

/* --------------------------------------------------------- path finding */

/* Takes a list of vnos and a message to display to the shell, then
   finds a list of vnos that form a connected path between the input
   vnos. Returns them in path. Will not return more than
   max_path. Returns the length of the path in path_length. */
int find_path ( int* vert_vno, int num_vno, char* message, int max_path_length,
                int* path, int* path_length );

/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------------------ tcl */

char** cached_tcl_commands = NULL;
int num_cached_tcl_commands = 0;
int max_num_cached_tcl_commands = 0;
void send_tcl_command (char* cmd);
void send_cached_tcl_commands ();

/* ---------------------------------------------------------------------- */

/* -------------------------------------------------------------- caption */

/* The caption is a string that is drawn on the top of the window when
   enabled. It can contain a mix of static text and substituted
   values, such as "My brain" or "My brain: Vertex !V" where !V is
   substituted with the vertex number that was clicked. Use the
   function cptn_set_format_string() to set the string, and the codes
   listed below. There is also a Tk interface with a list of codes. */

/* Flag to control whether it's drawn or not. */
int cptn_draw_flag = FALSE;

/* Location of the caption rectangle and string. This puts it in the
   top of the window. */
#define CPTN_FONT GLUT_BITMAP_8_BY_13
#define CPTN_LOC_RECTANGLE_LEFT -0.97
#define CPTN_LOC_RECTANGLE_RIGHT 0.97
#define CPTN_LOC_RECTANGLE_TOP 0.97
#define CPTN_LOC_RECTANGLE_BOTTOM 0.9
#define CPTN_LOC_CAPTION_X -0.95
#define CPTN_LOC_CAPTION_Y 0.92

/* The format string uses the codes listed below to stand in for
   values that will be calcualted when update_labels() is run. The
   value string is the working result of that, and is drawn to the
   screen. */
#define CPTN_STRING_LEN 2000
char* cptn_format_string = NULL;
char* cptn_value_string = NULL;

#define CPTN_CODE_VERTEXINDEX "!V" /* vertex index */
#define CPTN_CODE_DISTANCE "!D" /* distance */
#define CPTN_CODE_COORDS_RAS "!R" /* RAS coords */
#define CPTN_CODE_COORDS_MNITAL "!M" /* mni tal coords */
#define CPTN_CODE_COORDS_TAL "!T" /* tal coords */
#define CPTN_CODE_COORDS_INDEX "!I" /* MRI index */
#define CPTN_CODE_COORDS_NORMAL "!N" /* normal */
#define CPTN_CODE_COORDS_SPHERE_XYZ "!sxyz" /* spherical XYZ */
#define CPTN_CODE_COORDS_SPHERE_RT "!srt" /* spherical RT */
#define CPTN_CODE_CURVATURE "!C" /* curvature */
#define CPTN_CODE_FIELDSIGN "!F" /* fieldsign */
#define CPTN_CODE_FIELD_PREFIX "!o"
#define CPTN_CODE_FIELD0 "!o1" /* overlay layer 1 */
#define CPTN_CODE_FIELD1 "!o2" /* overlay layer 2 */
#define CPTN_CODE_FIELD2 "!o3" /* overlay layer 3 */
#define CPTN_CODE_FIELD3 "!o4" /* overlay layer 4 */
#define CPTN_CODE_FIELD4 "!o5" /* overlay layer 5 */
#define CPTN_CODE_FIELD5 "!o6" /* overlay layer 6 */
#define CPTN_CODE_FIELD6 "!o7" /* overlay layer 7 */
#define CPTN_CODE_FIELD7 "!o8" /* overlay layer 8 */
#define CPTN_CODE_FIELD8 "!o9" /* overlay layer 9 */
#define CPTN_CODE_CONDITION  "!cond" /* condition */
#define CPTN_CODE_TIME_POINT "!tp" /* time point */
#define CPTN_CODE_AMPLITUDE "!amp" /* amplitude */
#define CPTN_CODE_ANGLE "!ang" /* angle */
#define CPTN_CODE_DEGREE "!deg" /* degree */
#define CPTN_CODE_LABEL "!L" /* label */
#define CPTN_CODE_ANNOTATION "!A" /* annotation */
#define CPTN_CODE_MRIVALUE "!mriv" /* MRI value */
#define CPTN_CODE_PARCELLATION_NAME "!P" /* parcellation */

/* Initializes strings. */
int cptn_initialize ();

/* Call to set the format string for the caption, using the codes
   listed above. */
int cptn_set_format_string ( char* in );

/* Call to clear the value string. This is done because the actual
   caption is built up in the update_labels() function as the values
   are calculated. This will set the caption string to the format
   string.*/
int cptn_clear_value_string ();

/* Fill out a value for a code using sprintf format string and
   arguments. */
int cptn_sprintf_for_code (const char* code, const char* format, ... );

/* Draw the current caption to the window. */
int cptn_draw ();

/* Substitute a code in the caption string with an actual value. */
int cptn_substitute_code_with_value (char* caption, int caption_size,
                                     const char* code, char* value);

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------- debug */

/* Some drawing tools for drawing certain faces and verts. */
#define DEBUG_DRAWING_TOOLS 0
#if DEBUG_DRAWING_TOOLS
int* ddt_hilited_vnos = NULL;
int* ddt_hilited_fnos = NULL;
int ddt_initialize();
int ddt_clear();
int ddt_hilite_vertex (int vno, int type);
int ddt_get_hilite_vertex_color (int vno, GLubyte* r, GLubyte* g, GLubyte* b);
/* Note that we get color by vertex index, but set it by face index. */
int ddt_hilite_face (int fno, int type);
int ddt_get_hilite_face_color (int vno, GLubyte* r, GLubyte* g, GLubyte* b);
#else
#define ddt_initialize()
#define ddt_clear()
#define ddt_hilite_vertex(vno,type)
#define ddt_get_hilite_vertex_color(vno,r,g,b)
#define ddt_hilite_face(fno,type)
#define ddt_get_hilite_face_color(vno,r,g,b)
#endif

int link_tool_and_image_windows_flag = 1;

/* ---------------------------------------------------------------------- */

int save_tiff (char* fname);

/* Generates the edit.dat filename to use based on whether the subject
   name has been defined, etc. Same code exists in tkmedit to ensure
   the same names are generated. */
void copy_edit_dat_file_name (char* fname, int len);
static   char *gd_fnames[100] ;

/* end rkt */

/* ---------------------------- main() --------------------------------*/
#ifdef TCL
int Surfer(ClientData clientData, Tcl_Interp *interp, int argc, char *argv[])
#else
int  main(int argc,char *argv[])
#endif
{
  int tclscriptflag;
  int i;
  int j;
  long last_frame_xdim;
  FILE *fplog, *pplog;
  char str[NAME_LENGTH] = "";
  char logbuf[NAME_LENGTH];
  char lsubjectsdir[NAME_LENGTH];
  char lsrname[NAME_LENGTH];
  char lpname[NAME_LENGTH];
  char lstem[NAME_LENGTH];
  char lext[NAME_LENGTH];
  char *envptr;
  char cwd[100*NAME_LENGTH];
  char *word;
  char path[MAX_DIR_DEPTH][NAME_LENGTH];
  int  nargs = 0;
  char *functional_fnames[100], *patch_name = NULL, *overlayannotfname = NULL;
  int  noverlays = 0 ;
  /* begin rkt */
  FunD_tRegistrationType overlay_reg_type = FunD_tRegistration_Identity;
  char overlay_reg[NAME_LENGTH];

  int load_timecourse = FALSE;
  FunD_tRegistrationType timecourse_reg_type = FunD_tRegistration_Identity;
  char timecourse_fname[NAME_LENGTH];
  char timecourse_reg[NAME_LENGTH];

  int load_timecourse_offset = FALSE;
  FunD_tRegistrationType timecourse_offset_reg_type =
    FunD_tRegistration_Identity;
  char timecourse_offset_fname[NAME_LENGTH];
  char timecourse_offset_reg[NAME_LENGTH];

  int load_label = FALSE;
  int load_annotation = FALSE;
  int load_curv = FALSE;
  int load_sulc = FALSE;
  int overlayannotflag = FALSE; 
  char annotation_fname[NAME_LENGTH] = "";
  char *label_fnames[100] ;
  int  nlabels=0 ;

  int setlateralview = FALSE;
  int setmedialview = FALSE;

  int take_snapshot = FALSE;
  double pthresh = -1, pmid=0, pmax=0 ;
  char snap_fname[NAME_LENGTH] = "";

  int load_colortable = FALSE;
  char colortable_fname[NAME_LENGTH] = "";

  char tcl_cmd[1024];
  int err;

  InitDebugging("tksurfer") ;
  EnableDebuggingOutput ;

  // Read in env defaults
  if(getenv("FS_TKFTHRESH")) sscanf(getenv("FS_TKFTHRESH"),"%lf",&fthresh);
  if(getenv("FS_TKFMAX"))
  {
    sscanf(getenv("FS_TKFMAX"),"%lf",&tksfmax);
    fmid = (tksfmax+fthresh)/2.0;
    fslope = 1.0/(tksfmax-fthresh);
  }
  if(getenv("FS_TKS_OV_WIDE") && !strcmp(getenv("FS_TKS_OV_WIDE"),"1")) long_config_overlay = FALSE;

  /* This is a pre-pass of our arguments and just looks for -option
     options and copies their information for later use. It actually
     memmoves the argv contents so that the old arg parsing code can
     work as it did historically. In other words, this is a dirty hack
     that should be rewritten. */
  for (i = 0 ; i < argc ; i++)
  {
    /*      fprintf(stderr, "argv[%d] = %s\n", i, argv[i]);*/
    if (!stricmp(argv[i], "-o") || !stricmp(argv[i], "-ov") ||
        !stricmp(argv[i], "-overlay"))
    {
      nargs = 2 ;
      functional_fnames[noverlays] = argv[i+1] ;
      if (!fio_FileExistsReadable(functional_fnames[noverlays]))
      {
        printf("ERROR: cannot find %s\n",functional_fnames[noverlays]);
        exit(1);
      }
      gd_fnames[noverlays] = NULL ;
      load_curv = TRUE;
      noverlays++ ;
      forcegraycurvatureflag = TRUE;
      labl_draw_style = LABL_STYLE_OUTLINE;
    }
    else if (!stricmp(argv[i], "-fsgd") || !stricmp(argv[i], "-gdf") ||
        !stricmp(argv[i], "-gd") || !stricmp(argv[i], "-fsgdf"))
    {
      nargs = 3 ;
      gd_fnames[noverlays] = argv[i+1] ;
      functional_fnames[noverlays] = argv[i+2] ;
      if (!fio_FileExistsReadable(functional_fnames[noverlays]))
      {
        printf("ERROR: cannot find %s\n",functional_fnames[noverlays]);
        exit(1);
      }
      if (!fio_FileExistsReadable(gd_fnames[noverlays]))
      {
        printf("ERROR: cannot find %s\n",gd_fnames[noverlays]);
        exit(1);
      }
      load_curv = TRUE;
      noverlays++ ;
      forcegraycurvatureflag = TRUE;
      labl_draw_style = LABL_STYLE_OUTLINE;
    }
    else if (!stricmp(argv[i], "-overlay-annot-corrmatrix"))
    {
      nargs = 2;
      overlayannotfname = argv[i+1] ;
      if (!fio_FileExistsReadable(overlayannotfname))
      {
	printf("ERROR: cannot find %s\n",overlayannotfname);
	exit(1);
      }
      load_curv = TRUE;
      forcegraycurvatureflag = TRUE;
      overlayannotflag = TRUE;
      labl_draw_style = LABL_STYLE_OUTLINE;
    }
    else if (!stricmp(argv[i], "-overlay-reg") ||
             !stricmp(argv[i], "-orf") || !stricmp(argv[i], "-ovreg"))
    {
      nargs = 2 ;
      strncpy (overlay_reg, argv[i+1], sizeof(overlay_reg) );
      overlay_reg_type = FunD_tRegistration_File;
      if (!fio_FileExistsReadable(overlay_reg))
      {
        printf("ERROR: cannot find %s\n",overlay_reg);
        exit(1);
      }
    }
    else if (!stricmp(argv[i], "-reg"))
    {
      nargs = 2 ;
      strncpy (overlay_reg, argv[i+1], sizeof(overlay_reg) );
      overlay_reg_type = FunD_tRegistration_File;
      if (!fio_FileExistsReadable(overlay_reg)){
        printf("ERROR: cannot find %s\n",overlay_reg);
        exit(1);
      }
      strncpy (timecourse_reg, argv[i+1], sizeof(timecourse_reg));
      timecourse_reg_type = FunD_tRegistration_File;
    }
    else if (!stricmp(argv[i], "-overlay-reg-find"))
    {
      nargs = 1 ;
      overlay_reg_type = FunD_tRegistration_Find;
    }
    else if (!stricmp(argv[i], "-overlay-reg-identity"))
    {
      nargs = 1 ;
      overlay_reg_type = FunD_tRegistration_Identity;
    }
    else if (!stricmp(argv[i], "-mask_cortex"))
    {
      nargs = 2 ;
      mask_cortex_name = argv[i+1] ;
      printf("masking cortex label %s\n", mask_cortex_name) ;
    }
    else if ( !stricmp(argv[i], "-mni152reg" ) ){
      sprintf(overlay_reg,"%s/average/mni152.register.dat",
	      getenv("FREESURFER_HOME"));
      sprintf(timecourse_reg,"%s/average/mni152.register.dat",
	      getenv("FREESURFER_HOME"));
      overlay_reg_type = FunD_tRegistration_File;
      timecourse_reg_type = FunD_tRegistration_File;
      nargs = 1;
    } 
    else if (!stricmp(argv[i], "-zm"))
    {
      zero_mean = 1 ;
      printf("zero meaning input overlays\n") ;
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-fslope"))
    {
      nargs = 2 ;
      fslope = atof(argv[i+1]) ;
      fprintf(stderr, "setting fslope to %2.4f\n", fslope) ;
    }
    else if (!stricmp(argv[i], "-colscalebarflag"))
    {
      nargs = 2 ;
      colscalebarflag = atoi(argv[i+1]) ;
      fprintf(stderr, "setting colscalebarflag to %d\n", colscalebarflag) ;
    }
    else if (!stricmp(argv[i], "-colscaletext"))
    {
      nargs = 2 ;
      colscalebartextflag = atoi(argv[i+1]) ;
      fprintf(stderr, "setting colscalebartextflag to %d\n",
              colscalebartextflag) ;
    }
    else if (!stricmp(argv[i], "-colscalebarvertflag"))
    {
      nargs = 2 ;
      colscalebarvertflag = atoi(argv[i+1]) ;
      fprintf(stderr, "setting colscalebarvertflag to %d\n",
              colscalebarvertflag) ;
    }
     else if (!stricmp(argv[i], "-linkvertexmode"))
    {
      nargs = 2 ;
      linkvertexmode = atoi(argv[i+1]) ;
      fprintf(stderr, "setting linkvertexmode to %d\n",
              linkvertexmode) ;
    }
    else if (!stricmp(argv[i], "-scalebarflag"))
    {
      nargs = 2 ;
      scalebarflag = atoi(argv[i+1]) ;
      fprintf(stderr, "setting scalebarflag to %d\n", scalebarflag) ;
    }
    else if (!stricmp(argv[i], "-truncphaseflag"))
    {
      nargs = 2 ;
      truncphaseflag = atoi(argv[i+1]) ;
      fprintf(stderr, "setting truncphaseflag to %d\n", truncphaseflag) ;
    }
    else if (!stricmp(argv[i], "-revphaseflag"))
    {
      nargs = 2 ;
      revphaseflag = atoi(argv[i+1]) ;
      fprintf(stderr, "setting revphaseflag to %d\n", revphaseflag) ;
    }
    else if (!stricmp(argv[i], "-invphaseflag"))
    {
      nargs = 2 ;
      invphaseflag = atoi(argv[i+1]) ;
      fprintf(stderr, "setting invphaseflag to %d\n", invphaseflag) ;
    }
    else if (!stricmp(argv[i], "-fthresh"))
    {
      nargs = 2 ;
      fthresh = atof(argv[i+1]) ;
      fprintf(stderr, "setting fthresh to %2.4f\n", fthresh) ;
    }
    else if (!stricmp(argv[i], "-fminmax"))
    {
      nargs = 3 ;
      fthresh = atof(argv[i+1]) ;
      tksfmax = atof(argv[i+2]) ;
      fmid = (tksfmax+fthresh)/2.0;
      fslope = 1.0/(tksfmax-fthresh);
      printf("thresholds min=%g max=%g slope=%g mid=%g\n", fthresh,tksfmax,fslope,fmid) ;
    }
    else if (!stricmp(argv[i], "-patch"))
    {
      nargs = 2 ;
      patch_name = argv[i+1] ;
      fprintf(stderr, "displaying patch %s...\n", patch_name) ;
    }
    else if (!stricmp(argv[i], "-mask"))
    {
      nargs = 2 ;
      mrismaskfile = argv[i+1] ;
      fprintf(stderr, "mrismaskfile %s...\n", mrismaskfile) ;
      if (!fio_FileExistsReadable(mrismaskfile))
      {
        printf("ERROR: cannot find %s\n",mrismaskfile);
        exit(1);
      }
    }
    else if (!stricmp(argv[i], "-mask-thresh"))
    {
      nargs = 2 ;
      mrismaskthresh = atof(argv[i+1]) ;
      fprintf(stderr, "setting mrismaskthresh to %2.4f\n", mrismaskthresh) ;
    }
    else if (!stricmp(argv[i], "-fmid"))
    {
      nargs = 2 ;
      fmid = atof(argv[i+1]) ;
      fprintf(stderr, "setting fmid to %2.4f\n", fmid) ;
    }
    else if (!stricmp(argv[i], "-foffset"))
    {
      nargs = 2 ;
      foffset = atof(argv[i+1]) ;
      fprintf(stderr, "setting foffset to %2.4f\n", foffset) ;
    }
    else if (!stricmp(argv[i], "-offset"))
    {
      nargs = 2 ;
      offset = atof(argv[i+1]) ;
      fprintf(stderr, "setting offset to %2.4f\n", offset) ;
    }
    else if (!stricmp(argv[i], "-long-config-overlay"))
    {
      long_config_overlay = TRUE;
      nargs = 1 ;
      fprintf(stderr, "long view of configure overlay dialog enabled\n") ;
    }
    else if (!stricmp(argv[i], "-wide-config-overlay"))
    {
      long_config_overlay = FALSE;
      nargs = 1 ;
      fprintf(stderr, "wide view of configure overlay dialog enabled\n") ;
    }
    else if (!stricmp(argv[i], "-sdir"))
    {
      nargs = 2 ;
      sdir = argv[i+1] ;
      fprintf(stderr, "using SUBJECTS_DIR %s\n", sdir) ;
      FSENVsetSUBJECTS_DIR(sdir);
    }
    else if (!stricmp(argv[i], "-reassign"))
    {
      reassign = 1 ;
      nargs = 1 ;
      fprintf(stderr, "reassigning label vertex #s\n") ;
    }
    else if (!stricmp(argv[i], "-orig"))
    {
      nargs = 2 ;
      orig_suffix = argv[i+1] ;
      fprintf(stderr, "using orig suffix %s\n", orig_suffix) ;
    }
    else if (!stricmp(argv[i], "-sphere"))
    {
      nargs = 2 ;
      sphere_reg_suffix = argv[i+1] ;
      fprintf(stderr, "using sphere_reg suffix %s\n", sphere_reg_suffix) ;
    }
    else if (!stricmp(argv[i], "-white"))
    {
      nargs = 2 ;
      white_suffix = argv[i+1] ;
      fprintf(stderr, "using white suffix %s\n", white_suffix) ;
    }
    else if (!stricmp(argv[i], "-delink"))
    {
      link_tool_and_image_windows_flag = 0;
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-ovnew"))
    {
      UseNewOverlay = 1;
      nargs = 1 ;
    }
    /* begin rkt */
    else if (!stricmp(argv[i], "-timecourse") || !stricmp(argv[i], "-t"))
    {
      nargs = 2;
      strncpy (timecourse_fname, argv[i+1], sizeof(timecourse_fname));
      load_timecourse = TRUE;
      if (!fio_FileExistsReadable(timecourse_fname))
      {
        printf("ERROR: cannot find %s\n",timecourse_fname);
        exit(1);
      }
    }
    else if (!stricmp(argv[i], "-timecourse-reg") ||
             !stricmp(argv[i], "-treg"))
    {
      nargs = 2;
      strncpy (timecourse_reg, argv[i+1], sizeof(timecourse_reg));
      timecourse_reg_type = FunD_tRegistration_File;
      if (!fio_FileExistsReadable(timecourse_reg))
      {
        printf("ERROR: cannot find %s\n",timecourse_reg);
        exit(1);
      }
    }
    else if (!stricmp(argv[i], "-timecourse-reg-find"))
    {
      nargs = 1 ;
      timecourse_reg_type = FunD_tRegistration_Find;
    }
    else if (!stricmp(argv[i], "-timecourse-reg-identity"))
    {
      nargs = 1 ;
      timecourse_reg_type = FunD_tRegistration_Identity;
    }
    else if (!stricmp(argv[i], "-timecourse-offset"))
    {
      nargs = 2;
      strncpy (timecourse_offset_fname, argv[i+1],
               sizeof(timecourse_offset_fname));
      load_timecourse_offset = TRUE;
    }
    else if (!stricmp(argv[i], "-timecourse-offset-reg-file"))
    {
      nargs = 2;
      strncpy (timecourse_offset_reg, argv[i+1], 
               sizeof(timecourse_offset_reg));
      timecourse_offset_reg_type = FunD_tRegistration_File;
    }
    else if (!stricmp(argv[i], "-timecourse-offset-reg-find"))
    {
      nargs = 1 ;
      timecourse_offset_reg_type = FunD_tRegistration_Find;
    }
    else if (!stricmp(argv[i], "-timecourse-offset-reg-identity"))
    {
      nargs = 1 ;
      timecourse_offset_reg_type = FunD_tRegistration_Identity;
    }
    else if (!stricmp(argv[i], "-annotation") ||
             !stricmp(argv[i], "-annot"))
    {
      nargs = 2 ;
      strncpy (annotation_fname, argv[i+1], sizeof(annotation_fname));
      load_annotation = TRUE;
    }
    else if (!stricmp(argv[i], "-aparc"))
    {
      strncpy (annotation_fname, "aparc.annot", sizeof(annotation_fname));
      load_annotation = TRUE;
      labl_draw_style = LABL_STYLE_OUTLINE;
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-a2009s"))
    {
      strncpy (annotation_fname, "aparc.a2009s.annot", sizeof(annotation_fname));
      load_annotation = TRUE;
      labl_draw_style = LABL_STYLE_OUTLINE;
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-gldebug"))
    {
      // let glutInit parse this flag
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-snap") ||!stricmp(argv[i], "-snapshot"))
    {
      strcpy (snap_fname, argv[i+1]) ;
      take_snapshot = TRUE;
      nargs = 2 ;
    }
    else if (!stricmp(argv[i], "-pthresh"))
    {
      pthresh = atof(argv[i+1]) ;
      pmid = atof(argv[i+2]) ;
      pmax = atof(argv[i+3]) ;
      printf("setting percentile thresholds (%2.2f, %2.2f, %2.2f)\n",
             pthresh, pmid, pmax) ;
      nargs = 4 ;
    }
    else if (!stricmp(argv[i], "-lateral"))
    {
      setlateralview = TRUE;
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-medial"))
    {
      setmedialview = TRUE;
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-lrrev"))
    {
      LeftRightRev = TRUE;
      nargs = 1 ;
    }
    else if (!stricmp(argv[i], "-title"))
    {
      nargs = 2 ;
      tkstitle = argv[i+1];
    }
    else if (!stricmp(argv[i], "-curv"))
    {
      nargs = 1 ;
      load_curv = TRUE;
    }
    else if (!stricmp(argv[i], "-sulc"))
    {
      nargs = 1 ;
      load_sulc = TRUE;
    }
    else if (!stricmp(argv[i], "-gray"))
    {
      nargs = 1 ;
      load_curv = TRUE;
      forcegraycurvatureflag = TRUE;
    }
    else if (!stricmp(argv[i], "-labels-under"))
    {
      nargs = 1 ;
      labels_before_overlay_flag = TRUE;
    }
    else if (!stricmp(argv[i], "-label"))
    {
      nargs = 2 ;
      load_label = TRUE ;
      label_fnames[nlabels] = argv[i+1] ;
#if 0
      if (!fio_FileExistsReadable(label_fnames[nlabels]))
      {
        printf("ERROR: cannot find %s\n",label_fnames[nlabels]);
        exit(1);
      }
#endif
      load_curv = TRUE;
      nlabels++ ;
    }
    else if (!stricmp(argv[i], "-label-outline"))
    {
      nargs = 1 ;
      labl_draw_style = LABL_STYLE_OUTLINE;
    }
    else if (!stricmp(argv[i], "-colortable") ||
             !stricmp(argv[i], "-ctab"))
    {
      nargs = 2 ;
      strncpy (colortable_fname, argv[i+1], sizeof(colortable_fname));
      load_colortable = TRUE;
    }
    else if (!stricmp(argv[i], "-tcl"))
    {
      nargs = 2;
      strncpy (script_tcl, argv[i+1], sizeof(script_tcl));
      scriptok = TRUE;
    }
    /* end rkt */
    else
    {
      nargs = 0 ;
      if (i >= 4)
      {
        printf("WARNING: flag %s unrecognized\n",argv[i]);
        if (getenv("TK_EXIT_ON_CMD_ERROR")!=NULL)
        {
          printf("  ... and exiting because of it.\n");
          exit(1);
        }
      }
    }
    if (nargs > 0)
    {
      if (argc-(nargs+i) > 0)
        memmove(argv+i, argv+i+nargs, (argc-(nargs+i))*sizeof(char*)) ;
      argc -= nargs ;
      i-- ;
    }
  }

  /* args */
  if (argc==2)
    strcpy(str,argv[1]);
  if (MATCH_STR("-help") || MATCH_STR("-h"))
  {
    strcpy(str,argv[0]);
    if (MATCH_STR("tksurfer")) print_help_tksurfer();
    exit(1);
  }
  if (argc<4 || argc>6)
  {
    strcpy(str,argv[0]);
    if (MATCH_STR("surfer"))
    {
      printf("\n");
      printf("Usage: %s [-]name hemi surf\n",argv[0]);
      printf("       %s -help, -h [after startup: h,?]\n",argv[0]);
      printf("\n");
    }
    if (MATCH_STR("tksurfer"))
    {
      printf("\n");
      printf("Usage: %s [-]name hemi surf [-tcl script]\n",argv[0]);
      printf("       %s -help, -h [after startup: help]\n",argv[0]);
      printf("\n");
      printf("   name:  subjname (dash prefix => don't read MRI images)\n");
      printf("   hemi:  lh,rh\n");
      printf("   surf:  orig,smoothwm,plump,1000,1000a\n");
      printf("\n");
    }
    /* begin rkt */
    print_help_tksurfer();
    /* end rkt */
    exit(1);
  }
  if (argv[1][0]=='-')
  {
    MRIflag = FALSE;
    MRIloaded = FALSE;
    sprintf(lpname,"%s",argv[1]+1);
  }
  else
  {
#if 0
    MRIflag = TRUE;
    MRIloaded = TRUE;
#else
    MRIflag = FALSE;    /* don't load MRI volume - BRF */
    MRIloaded = FALSE;
#endif
    sprintf(lpname,"%s",argv[1]);
  }
  sprintf(lstem,"%s",argv[2]);
  sprintf(lext,"%s",argv[3]);
  tclscriptflag = FALSE;

  printf("subject is %s\n",lpname);
  printf("hemi    is %s\n",lstem);
  printf("surface is %s\n",lext);

  /* rkt: commented this part out. i'm not sure how it was accepting
     any other command line options with it active. */
#if 0
  if (argc>=5)
  {
    strcpy(str,argv[4]);
    if (MATCH_STR("-tcl"))
      tclscriptflag = TRUE;
    else
    {
      option = atoi(argv[4]);
      if (option==5)
      {
        printf("surfer: ### option 5 defunct--instead use:\n\n");
        printf("    tksurfer %s %s %s -tcl inflate.tcl\n\n",
               argv[1],argv[2],argv[3]);
      }
      else if (option==12)
      {
        printf("surfer: ### option 12 defunct--instead use:\n\n");
        printf("    tksurfer %s %s %s -tcl flatten.tcl\n\n",
               argv[1],argv[2],argv[3]);
      }
      else
      {
        printf("surfer: defunct option\n");
      }
      exit(1);
    }
  }
#endif

  /* subjects dir from environment */
  if (sdir)
    envptr = sdir ;
  else
    envptr = getenv("SUBJECTS_DIR");
  if (envptr==NULL)
  {
    printf("surfer: env var SUBJECTS_DIR undefined (use setenv)\n");
    exit(1);
  }
  strcpy(lsubjectsdir,envptr);
  printf("surfer: current subjects dir: %s\n",lsubjectsdir);

  /* OK, now we're going to do something silly and check whether the subject
   actually exists in SUBJECTS_DIR. If it does not, then we're going to do
   something REALLY silly and print out a coherent msg to the terminal.*/
  sprintf(tmpstr,"%s/%s",lsubjectsdir,lpname);
  if(!fio_IsDirectory(tmpstr)){
    printf("Checking %s\n",lpname);
    if(fio_FileExistsReadable(lpname)){
      /* Check whether subjectname is an actual file. Try opening it and reading
	 a string if it is */
      FILE *fp00;
      printf("Getting subjectname from %s\n",lpname);
      fp00 = fopen(lpname,"r");
      fscanf(fp00,"%s",lpname);
      fclose(fp00);
      printf("subject is %s\n",lpname);
      sprintf(tmpstr,"%s/%s",lsubjectsdir,lpname);
      if(!fio_IsDirectory(tmpstr)){
	printf("ERROR: cannot find subject %s in %s\n",lpname,lsubjectsdir);
	exit(1);
      }
    }
    else{
      printf("ERROR: cannot find subject %s in %s\n",lpname,lsubjectsdir);
      exit(1);
    }
  }

  /* guess datadir from cwd path */
  getcwd(cwd,100*NAME_LENGTH);
  word = strtok(cwd,"/");
  strcpy(path[0],word);
  i = 1;
  while ((word = strtok(NULL,"/")) != NULL)
  {  /* save,count */
    strcpy(path[i],word);
    i++;
  }
  if (MATCH(path[i-1],"scripts") && MATCH(path[i-2],"image"))
  {
    printf("surfer: in new (APD2) format \"scripts\" dir\n");
    j = i-1;
  }
  else if (MATCH(path[i-1],"scripts") && strspn(path[i-2],"0123456789")==5)
  {
    printf("surfer: in old (APD1) format \"scripts\" dir\n");
    j = i-1;
  }
  else if (MATCH(path[i-1],"scripts") && MATCH(path[i-2],lpname))
  {
    printf("surfer: in subjects \"scripts\" dir\n");
    j = i-1;
  }
  else if (MATCH(path[i-1],"scripts"))
  {
    printf("surfer: in \"scripts\" dir (not APD1,APD2,subjects format)\n");
    j = i-1;
  }
  else if (strstr(path[i-1],"scripts")!=NULL)
  {
    printf("surfer: in dir with \"scripts\" in name\n");
    j = i-1;
  }
  else
  {
    printf("surfer: not in \"scripts\" dir ==> using cwd for session root\n");
    j = i;
  }
  sprintf(lsrname,"/%s",path[0]);  /* reassemble absolute */
  for (i=1;i<j;i++)
  {
    // dont do this: sprintf(lsrname,"%s/%s",lsrname,path[i]);
    // instead, do this:
    strcat(lsrname,"/");
    strcat(lsrname,path[i]);
  }
  printf("surfer: session root data dir ($session) set to:\n");
  printf("surfer:     %s\n",lsrname);

  /* logfile */
  fplog = fopen("surfer.log","a");
  if (fplog==NULL)
  {
    printf("surfer: can't create file surfer.log in cwd\n");
    sprintf(str,"/tmp/surfer.log.%s",getenv("USER"));
    fplog = fopen(str,"a");
    if (fplog==NULL)
    {
      printf("surfer: can't create surfer.log--written to stdout\n");
      fplog = stdout;  /* give up and write log to standard out */
    }
    else
    {
      printf("surfer: surfer.log created in /tmp\n");
      strcpy(lsrname,"/tmp");
      printf("surfer: session root data dir ($session) reset to:\n");
      printf("surfer:     %s\n",lsrname);
    }
  }
  fprintf(fplog,"\n\n############################\n");
  pplog = popen("date","r");    /* date */
  fgets(logbuf,NAME_LENGTH,pplog);
  fprintf(fplog,"%s",logbuf);
  pclose(pplog);
  pplog = popen("pwd","r");     /* pwd */
  fgets(logbuf,NAME_LENGTH,pplog);
  fprintf(fplog,"%s",logbuf);
  pclose(pplog);
  for (i=0;i<argc;i++)           /* argv */
    fprintf(fplog,"%s ",argv[i]);
  fprintf(fplog,"\n############################\n");
  fflush(fplog);

  framebuff = (unsigned long *)calloc(frame_xdim*frame_ydim,sizeof(long));
  framebuff2 = (unsigned long *)calloc(frame_xdim*frame_ydim,sizeof(long));
  framebuff3 = (unsigned long *)calloc(frame_xdim*frame_ydim,sizeof(long));
  binbuff = (unsigned char *)calloc(3*frame_xdim*frame_ydim,sizeof(char));
  last_frame_xdim = frame_xdim;
  sf=0.55;

  printf("checking for nofix files in '%s'\n", lext) ;
  if ((strcmp(lext, "orig.nofix") == 0) || (strcmp(lext, "inflated.nofix") == 0))
  {
    printf("nofix surface detected - using nofix for orig and white\n") ;
    white_suffix = "orig.nofix" ;
    orig_suffix = "orig.nofix" ;
  }
  make_filenames(lsubjectsdir,lsrname,lpname,lstem,lext);

  if (MATCH(cwd,srname))
  {   /* if elsewhere, write rgbs in cwd */
    strcpy(gfname,srname);
    strcpy(sgfname,srname);
    sprintf(agfname,"%s/%s",srname,"surfer.rgb");
    printf("surfer: all rgb files written to %s\n",srname);
  }

  read_image_info(mfname);

  xmin = xx0;
  xmax = xx1;
  ymin = yy0;
  ymax = yy1;
  zmin = zz0;
  zmax = zz1;

  read_talairach(xffname);

  if (MRIflag) read_images(mfname);

  if (read_binary_surface(ifname) != NO_ERROR)
    ErrorExit(Gerror, "%s: could not read surface file %s.",Progname,ifname);

  if (read_orig_vertex_coordinates(orfname) == NO_ERROR)
  {
    char fname[STRLEN] ;
    sprintf(fname, "%s.%s", fpref, sphere_reg_suffix) ;
    if (FileExists(fname))
      read_canon_vertex_coordinates(fname) ;
  }

  if (noverlays > 0)
  {  /* -o specified on command line */
    int i ;
    char *functional_fname ;
    for (i = 0 ; i < noverlays ; i++)
    {
      functional_fname = functional_fnames[i] ;
      if (strlen(functional_fname) > 2 &&
          strcmp (&functional_fname[strlen(functional_fname)-2], ".w") == 0)
      {
        read_binary_values(functional_fname) ;
      }
      else
      {
        sclv_read_from_volume(functional_fname, overlay_reg_type,
                              overlay_reg, i);
      }
      if (gd_fnames[i])
      {
	sprintf(tcl_cmd, "GDF_Load %s\n", gd_fnames[i]) ;
	send_tcl_command (tcl_cmd);
//	gdfs[i] = gdfRead(gd_fnames[i], 1) ;
      }
    }
    overlayflag = TRUE ;
    colscale = HEAT_SCALE ;
  }

  if (load_curv)
  {
    read_binary_curvature(cfname) ;
    val_to_stat() ;
  }
  if (load_sulc)
  {
    read_binary_curvature(kfname) ;
    val_to_stat() ;
  }

#if 0
  read_binary_areas(afname);
#endif

  /* begin rkt */
  if (load_timecourse)
  {
    err = func_load_timecourse (timecourse_fname,
                                timecourse_reg_type, timecourse_reg);

    if ( ERROR_NONE == err )
    {
      /* send the number of conditions */
      sprintf (tcl_cmd, "Graph_SetNumConditions %d", func_num_conditions);
      send_tcl_command (tcl_cmd);

      /* show the graph window */
      send_tcl_command ("Graph_ShowWindow");
    }
  }

  if (load_timecourse_offset)
  {

    err = func_load_timecourse_offset (timecourse_offset_fname,
                                       timecourse_offset_reg_type,
                                       timecourse_offset_reg);

    if ( ERROR_NONE == err )
    {
      /* turn on the offset options */
      send_tcl_command ("Graph_ShowOffsetOptions 1");
    }
  }

  if (load_colortable)
  {
    labl_load_color_table (colortable_fname);
  }

  if (load_annotation)
  {
    err = labl_import_annotation (annotation_fname);
    if (err && getenv("TK_EXIT_ON_CMD_ERROR")!=NULL) exit(1);
  }
  
  if (overlayannotfname)
  {
    sclv_read_from_annotcorr(overlayannotfname, SCLV_VAL);
  }

  /* If we didn't load an annotation or color table filename, load the
     default color table. */
  if (!load_colortable && !load_annotation)
  {
#if 0
    char *freesurfer_home_envptr = getenv( "FREESURFER_HOME" );
    if ( NULL != freesurfer_home_envptr )
    {
      sprintf (colortable_fname, "%s/surface_labels.txt",
               freesurfer_home_envptr);
      fprintf( stderr, "Loading %s\n", colortable_fname );
      labl_load_color_table( colortable_fname );
    }
#else
    // there is no 'default' colortable, since there are multiple parc schemes
    //fprintf( stderr, "\n\nWARNING: No colortable found!\n\n");
#endif
  }

  /* end rkt */

  if (mask_cortex_name)
    mask_label(mask_cortex_name) ;
  if (tclscriptflag)
  {
    /* tksurfer tcl script */
    /* called from tksurfer.c; do nothing (don't even open gl window) */
    /* wait for tcl interp to start; tksurfer calls tcl script */
  }
  else
  {

    /* open window for surfer or non-script tksurfer (a few envs) */
//#ifndef Linux
    if ((envptr=getenv("doublebufferflag"))!=NULL)
    { /*tmp:TODO OGL toggle*/
      if (MATCH("1",envptr))     doublebufferflag = TRUE;
      if (MATCH("TRUE",envptr))  doublebufferflag = TRUE;
    }
    if ((envptr=getenv("renderoffscreen"))!=NULL)
    {
      if (MATCH("1",envptr))     renderoffscreen = TRUE;
      if (MATCH("TRUE",envptr))  renderoffscreen = TRUE;
    }
//#endif
    if (tkstitle == NULL) tkstitle = pname;
    open_window(tkstitle);
    if (stem[0]=='r'&&stem[1]=='h')
      rotate_brain(-90.0,'y');
    else
      rotate_brain(90.0,'y');
    redraw();
  }

  if (patch_name)
  {
    strcpy(pfname, patch_name) ;
    read_binary_patch(patch_name) ;
    restore_zero_position() ;
    rotate_brain(-90.0, 'x') ;
    redraw() ;
  }

  if (setlateralview)
  {
    printf("setting view to lateral...\n") ;
    make_lateral_view(stem);
    redraw() ;
  }

  if (setmedialview)
  {
    printf("setting view to medial...\n") ;
    make_medial_view(stem);
    redraw() ;
  }

  
  printf("setting percentile thresholds (%2.2f, %2.2f, %2.2f)\n",
         pthresh, pmid, pmax) ;
  if (pthresh >= 0)
  {
    sclv_set_current_threshold_from_percentile(pthresh, pmid, pmax) ;
    redraw() ;
  }

  if (take_snapshot)
  {
    printf("saving snapshot to %s and exiting\n", snap_fname) ;
    save_tiff(snap_fname) ;
    exit(0) ;
  }

  if (load_label)
    {
      int  i ;
      for (i = 0 ; i < nlabels ; i++)
	{
	  printf("loading label %s\n", label_fnames[i]) ;
	  err = labl_load (label_fnames[i]);
	  if (err && getenv("TK_EXIT_ON_CMD_ERROR")!=NULL) exit(1);
	}
    }

  sclv_send_current_field_info ();  // sends fmid/fthresh to gui
  return(0) ;
}

#ifdef TCL
int
do_one_gl_event(Tcl_Interp *interp)   /* tcl */
{

#ifdef OPENGL     /* derived from event.c:DoNextEvent(),tkExec() */
  XEvent current, ahead;
  char buf[1000],*tclvar;
  char command[NAME_LENGTH];
  KeySym ks;
  static int ctrlkeypressed = FALSE;
  static int altkeypressed = FALSE;
  static int shiftkeypressed = FALSE;
  static int button1pressed = FALSE;
  static int button2pressed = FALSE;
  static int button3pressed = FALSE;
  short sx,sy;
  long ox,oy,lx,ly;
  float wx, wy;
  XWindowAttributes wat;
  Window junkwin;
  int rx, ry;
  /* begin rkt */
  int mvno;
  float md;
  int vno;
  /* end rkt */
  int tpvno;
  float tpd;

  if (!openglwindowflag) return(0);

  if (XPending(xDisplay))
  {  /* do one if queue test */

    XNextEvent(xDisplay, &current);   /* blocks here if no event */

    switch (current.type)
    {

    case ConfigureNotify:
      // Generated when window moves
      XGetWindowAttributes(xDisplay, w.wMain, &wat);
      XTranslateCoordinates(xDisplay, w.wMain, wat.root,
                            -wat.border_width, -wat.border_width,
                            &rx, &ry, &junkwin);
      w.x = rx;  /* left orig         (current.xconfigure.x = 0 relative!) */
      w.y = ry;  /* top orig: +y down (current.xconfigure.y = 0 relative!) */
      w.w = current.xconfigure.width;
      w.h = current.xconfigure.height;
      resize_window(0);
      tclvar =
        (char*)Tcl_GetVar(interp,(char*)"tksurferinterface",TCL_GLOBAL_ONLY);
      /* begin rkt */

      if (link_tool_and_image_windows_flag)
      {
        /* link tool window with image window */
        sprintf(command,"MoveToolWindow %d %d",
                w.x, w.y + w.h + MOTIF_YFUDGE /*+MOTIF_XFUDGE*/);
        send_tcl_command (command);
      }

      /* update window position */
      curwindowleft = w.x;
      curwindowbottom = w.y + w.h;

      /* if our redraw lock flag is on, redraw the window. check to
         see if there are any expose or configure events ahead of us,
         and if so, don't redraw now. this saves us from having
         multiple redraw flashes. */
      if (redrawlockflag)
      {
        if (XPending(xDisplay))
        {
          XPeekEvent(xDisplay, &ahead);
          if (ahead.type==Expose || ahead.type==ConfigureNotify)
            break;
        }
        send_tcl_command ("UpdateAndRedraw");
      }
      /* end rkt */
      break;

    case Expose:
      /* begin rkt */
      /* if our redraw lock flag is on, redraw the window. check to
         see if there are any expose or configure events ahead of us,
         and if so, don't redraw now. this saves us from having
         multiple redraw flashes. */
      if (redrawlockflag)
      {
        if (XPending(xDisplay))
        {
          XPeekEvent(xDisplay, &ahead);
          if (ahead.type==Expose || ahead.type==ConfigureNotify)
            break;
        }
        send_tcl_command ("UpdateAndRedraw");
      }
#if 0
      if (XPending(xDisplay))
      {
        XPeekEvent(xDisplay, &ahead);
        if (ahead.type==Expose || ahead.type==ConfigureNotify) break;
      }
      /*send_tcl_command ("redrawbutton");*/ /* still multiple */
#endif
      /* end rkt */
      break;

      /* begin rkt */
    case MotionNotify:
      if (mouseoverflag)
      {
        if (XPending(xDisplay))
        {
          XPeekEvent(xDisplay, &ahead);
          if (ahead.type==MotionNotify)
            break;
        }
        sx = current.xmotion.x;
        sy = current.xmotion.y;
        sx += w.x;   /* convert back to screen pos (ugh) */
        sy = 1024 - w.y - sy;
        find_vertex_at_screen_point( sx, sy, &mvno, &md );
        if ( mvno >= 0 )
          mouseover_vno = mvno;
        update_labels(LABELSET_MOUSEOVER, mouseover_vno, md);
      }
      break;
      /* end rkt */

    case ButtonPress:
      sx = current.xbutton.x;
      sy = current.xbutton.y;
      sx += w.x;   /* convert back to screen pos (ugh) */
      sy = 1024 - w.y - sy;
      if (current.xbutton.button == 1)
      {  /** left **/
        button1pressed = TRUE;
        if (ctrlkeypressed)
        { /* scale around click */
          getorigin(&ox,&oy);
          getsize(&lx,&ly);
          wx = sf*(sx-ox-lx/2.0)*2.0*fov/lx;
          wy = sf*(sy-oy-ly/2.0)*2.0*fov/ly;
          /* begin rkt */
#if 0
          sprintf(command,"set xtrans %d; set ytrans %d",(int)-wx,(int)-wy);
          send_tcl_command (command);
          sprintf(command,"set scalepercent %d",(int)SCALE_UP_MOUSE*100);
          send_tcl_command (command);
          send_tcl_command ("redrawbutton");
#endif
          translate_brain (-wx, -wy, 0);
          scale_brain (SCALE_UP_MOUSE);
          redraw();
          /* end rkt */
        }
        else if (shiftkeypressed)
        {  /* curvim */
          select_vertex(sx,sy);
          read_curvim_at_vertex(selection);
          draw_cursor(selection,TRUE);
        }
        else
        {
	  /* if the vertex index is linked to the timepoint in the config overlay dialog */
	  if (linkvertexmode == 1)
	  {
            find_vertex_at_screen_point(sx, sy, &tpvno, &tpd);
	    sclv_cur_timepoint = tpvno;
            sclv_set_timepoint_of_field(sclv_current_field, sclv_cur_timepoint, sclv_cur_condition);
            send_tcl_command ("UpdateAndRedraw");
	  }

	  /* if the vertex index is linked to the avg of the ROI in the config overlay dialog */
	  if (linkvertexmode == 2 || linkvertexmode == 3)
          {
            find_vertex_at_screen_point(sx, sy, &tpvno, &tpd);
            link_timepoint_ROI(tpvno);
	  }
          /* begin rkt */
          if (selection>=0)
          {
            draw_vertex_hilite(selection);
            draw_cursor(selection,FALSE);

          }
          select_vertex(sx,sy);
          if (selection>=0)
          {
            mark_vertex(selection,TRUE);
            draw_cursor(selection,TRUE);
          }
          /* draw all the marked verts. */
          draw_marked_vertices ();
          /* end rkt */
          if (showorigcoordsflag)
            find_orig_vertex_coordinates(selection);
        }
      }
      if (current.xbutton.button == 2)
      {  /** middle **/
        button2pressed = TRUE;

        /* begin rkt */
        find_closest_marked_vertex (sx, sy, NULL, &vno);
        if (vno>=0 && vno < mris->nvertices)
        {
          fprintf (stderr, "Unmarking %d\n", vno);
          mark_vertex (vno, FALSE);
          draw_cursor (vno, FALSE);
        }
        /* end rkt */
#if 0
        select_vertex(sx,sy);
        if (selection>=0)
          mark_vertex(selection,FALSE);
#endif

      }
      if (current.xbutton.button == 3)
      {  /** right **/
        button3pressed = TRUE;
        if (ctrlkeypressed)
        {
          getorigin(&ox,&oy);
          getsize(&lx,&ly);
          wx = sf*(sx-ox-lx/2.0)*2.0*fov/lx;
          wy = sf*(sy-oy-ly/2.0)*2.0*fov/ly;
          /* begin rkt */
#if 0
          sprintf(command,"set xtrans %d; set ytrans %d",(int)-wx,(int)-wy);
          send_tcl_command (command);
          sprintf(command,"set scalepercent %d",(int)(1.0/SCALE_UP_MOUSE*100));
          send_tcl_command (command);
          send_tcl_command ("redrawbutton");
#endif
          translate_brain (-wx, -wy, 0);
          scale_brain (1.0/SCALE_UP_MOUSE);
          redraw();
          /* end rkt */
        }
        else
        {
          clear_all_vertex_marks();
          /* begin rkt */
          /* deselect label. */
          labl_select(-1);
          redraw();
          /* end rkt */
        }
      }

      break;

    case ButtonRelease:
      if (current.xbutton.button == 1)  button1pressed = FALSE;
      if (current.xbutton.button == 2)  button2pressed = FALSE;
      if (current.xbutton.button == 3)  button3pressed = FALSE;
      break;


    case KeyPress:
      XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
      switch (ks)
      {

        /* numbers */
      case XK_0:
        if (altkeypressed)  ;
        break;

        /* upper case */
      case XK_A:
        if (altkeypressed)  ;
        break;

        /* lower case */
        /* begin rkt */
#if 0
      case XK_r:
        if (altkeypressed) send_tcl_command ("redrawbutton");
        break;
#endif
      case XK_r:
        if (altkeypressed)
          send_tcl_command ("UpdateAndRedraw");
        break;

      case XK_q:
        if (ctrlkeypressed)
          exit (0);
        break;
        /* end rkt */

        /* others */
        /* begin rkt */
#if 0
      case XK_Up:
        if (altkeypressed) send_tcl_command ("set xrot [expr $xrot+18.0]");
        break;
      case XK_Down:
        if (altkeypressed) send_tcl_command ("set xrot [expr $xrot-18.0]");
        break;
      case XK_Right:
        if (altkeypressed) send_tcl_command ("set yrot [expr $yrot-18.0]");
        break;
      case XK_Left:
        if (altkeypressed) send_tcl_command ("set yrot [expr $yrot+18.0]");
        break;
#endif
      case XK_Up:
        if (altkeypressed)
          send_tcl_command ("set gNextTransform(rotate,x) "
                            "[expr $gNextTransform(rotate,x)+18.0]");
        break;
      case XK_Down:
        if (altkeypressed)
          send_tcl_command ("set gNextTransform(rotate,x) "
                            "[expr $gNextTransform(rotate,x)-18.0]");
        break;
      case XK_Right:
        if (altkeypressed)
          send_tcl_command ("set gNextTransform(rotate,y) "
                            "[expr $gNextTransform(rotate,y)-18.0]");
        break;
      case XK_Left:
        if (altkeypressed)
          send_tcl_command ("set gNextTransform(rotate,y) "
                            "[expr $gNextTransform(rotate,y)+18.0]");
        break;
      case XK_Home:
        sprintf(command,"MoveToolWindow %d %d",
                curwindowleft, curwindowbottom + MOTIF_YFUDGE);
        send_tcl_command (command);
        break;
        /* end rkt */

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

    case KeyRelease:   /* added this mask to owindow.c */
      XLookupString(&current.xkey, buf, sizeof(buf), &ks, 0);
      switch (ks)
      {
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
      }
      break;

    }
    return GL_FALSE;
  }

#else  /* use gl calls */
  // This stuff should be deletable 10/5/07
  short dev, val;
  static int ctrlkeypressed = FALSE;
  static int altkeypressed = FALSE;
  static int shiftkeypressed = FALSE;
  short sx,sy; /* Screencoord sx,sy; */
  float wx,wy; /* Coord wx,wy; */
  long ox,oy,lx,ly;
  char fname[200],fpref[200];
  int i;
  FILE *fp;
  float a,b,c;

  float tmid,dt,f1,f2,f3;
  int num;

  blinkbuffers();

  if (qtest())
  {  /* do one event */
    dev = qread(&val);
    if (dev != LEFTMOUSE &&
        dev != RIGHTMOUSE) /* hack: mouse zeros getbutton! */
      ctrlkeypressed = getbutton(LEFTCTRLKEY) || getbutton(RIGHTCTRLKEY);
    switch (dev)
    {
    case REDRAW:
      resize_window(0);   /* reshape flag */
      break;
    case LEFTMOUSE:
      if (val == 0)  break;
      sx = getvaluator(MOUSEX);
      sy = getvaluator(MOUSEY);
      if (ctrlkeypressed)
      { /* scale around click */
        getorigin(&ox,&oy);
        getsize(&lx,&ly);
        wx = sf*(sx-ox-lx/2.0)*2.0*fov/lx;
        wy = sf*(sy-oy-ly/2.0)*2.0*fov/ly;
        translate_brain(-wx,-wy,0.0);
        scale_brain(SCALE_UP_MOUSE);
        redraw();
      }
      else if (shiftkeypressed)
      {
        select_vertex(sx,sy);
        read_curvim_at_vertex(selection);
        draw_cursor(selection,TRUE);
      }
      else
      {
        if (selection>=0)
          draw_cursor(selection,FALSE);
        select_vertex(sx,sy);
        if (selection>=0)
          mark_vertex(selection,TRUE);
        if (selection>=0)
          draw_cursor(selection,TRUE);
      }
      break;
    case MIDDLEMOUSE:
      if (val == 0)  break;
      sx = getvaluator(MOUSEX);
      sy = getvaluator(MOUSEY);
      select_vertex(sx,sy);
      if (selection>=0)
        mark_vertex(selection,FALSE);
      break;
    case RIGHTMOUSE:
      if (val == 0)  break;
      sx = getvaluator(MOUSEX);
      sy = getvaluator(MOUSEY);
      if (ctrlkeypressed)
      {
        getorigin(&ox,&oy);
        getsize(&lx,&ly);
        wx = sf*(sx-ox-lx/2.0)*2.0*fov/lx;
        wy = sf*(sy-oy-ly/2.0)*2.0*fov/ly;
        translate_brain(-wx,-wy,0.0);
        scale_brain(1.0/SCALE_UP_MOUSE);
        redraw();
      }
      else
      {
        clear_vertex_marks();
      }
      break;
    case UPARROWKEY:
      if (!altkeypressed) break;
      if (val == 0)  break;
      rotate_brain(18.0,'x');
      break;
    case DOWNARROWKEY:
      if (!altkeypressed) break;
      if (val == 0)  break;
      rotate_brain(-18.0,'x');
      break;
    case LEFTARROWKEY:
      if (!altkeypressed) break;
      if (val == 0)  break;
      rotate_brain(18.0,'y');
      break;
    case RIGHTARROWKEY:
      if (!altkeypressed) break;
      if (val == 0)  break;
      rotate_brain(-18.0,'y');
      break;
    case PAD2:
      if (!altkeypressed) break;
      if (val == 0)  break;
      translate_brain(0.0,-10.0,0.0);
      break;
    case PAD8:
      if (!altkeypressed) break;
      if (val == 0)  break;
      translate_brain(0.0,10.0,0.0);
      break;
    case PAD4:
      if (!altkeypressed) break;
      if (val == 0)  break;
      translate_brain(-10.0,0.0,0.0);
      break;
    case PAD6:
      if (!altkeypressed) break;
      if (val == 0)  break;
      translate_brain(10.0,0.0,0.0);
      break;
    case DELKEY:
      if (!altkeypressed) break;
      if (val == 0)  break;
      rotate_brain(18.0,'z');
      break;
    case PAGEDOWNKEY:
      if (!altkeypressed) break;
      if (val == 0)  break;
      rotate_brain(-18.0,'z');
      break;
    case LEFTALTKEY:
    case RIGHTALTKEY:
      if (val)  altkeypressed = TRUE;
      else      altkeypressed = FALSE;
      break;
    case LEFTSHIFTKEY:
    case RIGHTSHIFTKEY:
      if (val)  shiftkeypressed = TRUE;
      else      shiftkeypressed = FALSE;
      break;
    case KEYBD:
      if (altkeypressed)
        switch ((char)val)
        {
        case 'f':
          find_orig_vertex_coordinates(selection);
          break;
        case 'F':
          select_orig_vertex_coordinates(&selection);
          break;
        case 'o':
          read_orig_vertex_coordinates(orfname);
          break;
        case 'O':
          printf("enter a, b, c: ");
          scanf("%f %f %f",&a,&b,&c);
          read_ellipsoid_vertex_coordinates(elfname,a,b,c);
          break;
        case 'i':
          floodfill_marked_patch(RIPFILL);
          printf("floodfill_marked_patch(%d)\n",RIPFILL);
          break;
        case 'I':
          floodfill_marked_patch(STATFILL);
          printf("floodfill_marked_patch(%d)\n",STATFILL);
          break;
        case 'u':
          restore_ripflags(1);
          printf("restore_ripflags(1)\n");
          break;
        case 'j':
          cut_line(FALSE);
          printf("cut_line\n");
          break;
        case 'J':
          cut_plane();
          printf("cut_plane\n");
          break;
        case 'L':
          flatten(tfname);
          flag2d=TRUE;
          printf("flatten\n");
          break;
        case 'p':
          cut_vertex();
          break;
        case 'P':
          cut_line(TRUE);
          break;
        case 'r':
          redraw();
          PR
          break;
        case 'R':
          printf("enter .iop-name "
                 "and dipole-file-number: ");
          scanf("%s %d",fname,&num);
          read_iop(fname,num);
          printf("enter number of rec files: ");
          scanf("%d",&num);
          sol_nrec = 0;
          for (i=0;i<num;i++)
          {
            printf("enter .rec-name: ");
            scanf("%s",fname);
            read_rec(fname);
            filter_rec(sol_nrec-1);
          }
          compute_timecourses();
          PR
          break;
        case 'T':
          if (selection>=0)
            draw_cursor(selection,FALSE);
          find_nearest_fixed(selection,&selection);
          if (selection>=0)
            draw_cursor(selection,TRUE);
          PR
          break;
        case 't':
          plot_nearest_nonzero(selection);
          break;
        case 'l':
          printf("enter .vlist file name: ");
          scanf("%s",fname);
          read_plot_list(fname);
          PR
          break;
          /*
          case 'I': printf("enter .pri file name: ");
          scanf("%s",fname);
          read_binary_values(fname);
          PR
          break;
          */
        case 'D':
          printf("enter .dec file name: ");
          scanf("%s",fpref);
          sprintf(fname,"%s",fpref);
          read_binary_decimation(fname);
          printf("enter .dip file name: ");
          scanf("%s",fpref);
          sprintf(fname,"%s",fpref);
          read_binary_dipoles(fname);
          PR
          break;
        case 'd':
          printf("enter tmid dt num: ");
          scanf("%f %f %d",&tmid,&dt,&num);
          load_vals_from_sol(tmid,dt,num);
          break;
        case 'C':
          if (sol_nrec>1)
          {
            printf("enter num: ");
            scanf("%d",&num);
          }
          else
            num = 0;
          load_var_from_sol(num);
          PR
          break;
        case 'c':
          printf("enter plot type: ");
          scanf("%d",&sol_plot_type);
          PR
          break;
        case 'S':
          printf("enter type: ");
          scanf("%d",&num);
          normalize_time_courses(num);
          break;
        case 's':
          printf("enter pthresh (%f) pslope "
                 "(%f) maxrat (%f): ",
                 sol_pthresh,sol_pslope,sol_maxrat);
          scanf("%f %f %f",&f1,&f2,&f3);
          sol_pthresh = f1;
          sol_pslope = f2;
          sol_maxrat = f3;
          PR
          break;
        case 'x':
          compute_select_crosstalk();
          PR
          break;
        case 'X':
          printf("enter type (0=unweighted, "
                 "1=distance-weighted): ");
          scanf("%d",&num);
          compute_all_crosstalk(num);
          PR
          break;
        case 'g':
          compute_select_pointspread();
          PR
          break;
        case 'G':
          compute_all_pointspread();
          PR
          break;
        case 'H':
          if (selection>=0)
            draw_cursor(selection,FALSE);
          printf("enter vertex index (%d): ",selection);
          scanf("%d",&selection);
          if (selection>=0)
            draw_cursor(selection,TRUE);
          PR
          break;
        case 'h':
          printf("enter min max nbins: ");
          scanf("%f %f %d",&f1,&f2,&num);
          write_val_histogram(f1,f2,num);
          PR
          break;
        }
    }
  }
#endif
  return(0) ;
}
#endif


void
read_ncov(char *fname)
{
  int i,j,tnmeg,tneeg;
  float f;
  FILE *fptr;
  float regconst=0.01;

  printf("read_ncov(%s)\n",fname);
  fptr = fopen(fname,"r");
  if (fptr==NULL)
  {
    printf("can't find file %s\n",fname);
    exit(1);
  }
  fscanf(fptr,"%*s");
  fscanf(fptr,"%d %d",&tnmeg,&tneeg);
  if ((sol_nmeg!=0 && tnmeg!=sol_nmeg) || (sol_nmeg!=0 && tnmeg!=sol_nmeg))
  {
    printf("incompatible nmeg or meeg in ncov file\n");
    exit(1);
  }
  sol_NoiseCovariance = MGH_matrix(sol_nchan,sol_nchan);
  for (i=0;i<tneeg+tnmeg;i++)
    for (j=0;j<tneeg+tnmeg;j++)
    {
      fscanf(fptr,"%f",&f);
      if (sol_nmeg==0)
      {
        if (i<tneeg && j<tneeg)
          sol_NoiseCovariance[i][j] = f;
      }
      else
        if (sol_neeg==0)
        {
          if (i>=tneeg && j>=tneeg)
            sol_NoiseCovariance[i-tneeg][j-tneeg] = f;
        }
        else
          sol_NoiseCovariance[i][j] = f;
    }
  /*
    for (i=0;i<sol_nchan;i++)
    printf("%d: %9.4e\n",i,sol_NoiseCovariance[i][i]);
  */
  regularize_matrix(sol_NoiseCovariance,sol_nchan,regconst);
  fclose(fptr);
}

void
regularize_matrix(FLOATTYPE **m, int n, float fact)
{
  int i;                   /* loop counters */
  float sum;

  sum = 0;
  for (i=0;i<n;i++)
    sum += m[i][i];
  sum /= n;
  /*
    printf("avg. diag. = %g\n",sum);
  */
  if (sum==0) sum=1;
  if (sum!=0)
    for (i=0;i<n;i++)
      m[i][i] += sum*fact;
}

void
read_iop(char *fname, int dipfilenum)
{
  int i,j,k,jc,d;
  FILE *fp;
  char c,str[200];
  float f;
  int iop_vernum = 0;

  printf("read_iop(%s,%d)\n",fname,dipfilenum);

  sol_pthresh = 1000;
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### can't find file %s\n",fname);
    PR return;
  }
  c = fgetc(fp);
  if (c=='#')
  {
    fscanf(fp,"%s",str);
    if (MATCH_STR("version"))
      fscanf(fp,"%d",&iop_vernum);
    printf("iop version = %d\n",iop_vernum);
    fscanf(fp,"%d %d %d %d",&sol_neeg,&sol_nmeg,&sol_nperdip,&sol_ndipfiles);
    sol_nchan = sol_neeg+sol_nmeg;
    for (i=1;i<=dipfilenum;i++)
    {
      fscanf(fp,"%d %d",&sol_ndip,&sol_nnz);
      if (i==dipfilenum)
      {
        if (sol_ndip != sol_ndec)
        {
          printf(".dec and .iop file mismatch (%d!=%d)\n",
                 sol_ndip,sol_ndec);
          fclose(fp);
          return;
        }
        sol_W = MGH_matrix(sol_nnz*sol_nperdip,sol_nchan);
        if (iop_vernum==1)
          sol_A = MGH_matrix(sol_nchan,sol_nnz*sol_nperdip);

        sol_M = MGH_matrix(sol_nchan,sol_nchan);
        /* temp. space for xtalk */
        sol_Mi = MGH_matrix(sol_nchan,sol_nchan);
        sol_sensvec1 = MGH_vector(sol_nchan);
        sol_sensvec2 = MGH_vector(sol_nchan);

        sol_sensval = MGH_vector(sol_nchan);
        sol_dipindex = MGH_ivector(sol_nnz);
        sol_pval = MGH_vector(sol_nnz);
        sol_prior = MGH_vector(sol_nnz);
        sol_dipfact = MGH_vector(sol_nnz);
      }
      for (j=0;j<sol_nnz;j++)
      {
        if (i==dipfilenum)
        {
          fscanf(fp,"%d",&d);
          sol_dipindex[j] = d;
        }
        else
          fscanf(fp,"%*d");
      }
      for (j=0;j<sol_nnz;j++)
      {
        if (i==dipfilenum)
        {
          fscanf(fp,"%f",&f);
          sol_pval[j] = f;
          f = fabs(f);
          if (f<sol_pthresh) sol_pthresh = f;
        }
        else
          fscanf(fp,"%*f");
      }
      for (j=0;j<sol_nnz;j++)
      {
        if (i==dipfilenum)
        {
          fscanf(fp,"%f",&f);
          sol_prior[j] = f;
          mris->vertices[sol_dipindex[j]].val = f;
        }
        else
          fscanf(fp,"%*f");
      }
      for (j=0;j<sol_nnz;j++)
        for (jc=0;jc<sol_nperdip;jc++)
          for (k=0;k<sol_nchan;k++)
          {
            if (i==dipfilenum)
            {
              fscanf(fp,"%f",&f);
              sol_W[j*sol_nperdip+jc][k] = f;
            }
            else
              fscanf(fp,"%*f");
          }
      if (iop_vernum==1)
        for (j=0;j<sol_nnz;j++)
          for (jc=0;jc<sol_nperdip;jc++)
            for (k=0;k<sol_nchan;k++)
            {
              if (i==dipfilenum)
              {
                fscanf(fp,"%f",&f);
                sol_A[k][j*sol_nperdip+jc] = f;
              }
              else
                fscanf(fp,"%*f");
            }
    }
  }
  else
  {
    printf("Can't read binary .iop files\n");
  }
  fclose(fp);
}

void
read_rec(char *fname)
{
  int i,j,tntime,tptime,tnmeg,tneeg,tnchan;
  float f;
  FILE *fp;

  printf("read_rec(%s)\n",fname);

  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("can't find file %s\n",fname);
    exit(1);
  }
  if (sol_nmeg+sol_nmeg<=0)
  {
    printf("iop-file must be read before rec-file\n");
    return;
  }
  fscanf(fp,"%*s");
  fscanf(fp,"%d %d %d",&tntime,&tnmeg,&tneeg);
  tnchan = tnmeg+tneeg;
  if (sol_ntime>0 && sol_ntime!=tntime)
  {
    printf("ntime does not match tntime (%d != %d)\n",sol_ntime,tntime);
    exit(1);
  }
  tptime = (int)(2*pow(2.0,ceil(log((float)tntime)/log(2.0))));
  printf("tntime=%d\n",tntime);
  printf("tptime=%d\n",tptime);
  sol_ntime = tntime;
  sol_ptime = tptime;
  if (sol_nmeg>0 && sol_nmeg!=tnmeg)
  {
    printf("nmeg does not match tnmeg (%d != %d)\n",sol_nmeg,tnmeg);
    exit(1);
  }
  if (sol_neeg>0 && sol_neeg!=tneeg)
  {
    printf("neeg does not match tnmeg (%d != %d)\n",sol_neeg,tneeg);
    exit(1);
  }
  if (sol_nrec==0)
    sol_lat = MGH_vector(sol_ptime);
  sol_Data[sol_nrec] = MGH_matrix(tnchan,sol_ptime);
  for (j=0;j<tntime;j++)
  {
    fscanf(fp,"%f",&f);
    sol_lat[j] = f;
    for (i=0;i<tneeg;i++)
    {
      fscanf(fp,"%f",&f);
      if (sol_neeg>0)
        sol_Data[sol_nrec][i][j] = f;
    }
    for (i=0;i<tnmeg;i++)
    {
      fscanf(fp,"%f",&f);
      if (sol_nmeg>0)
        sol_Data[sol_nrec][i+sol_neeg][j] = f;
    }
  }
  fclose(fp);
  sol_sample_period = sol_lat[1]-sol_lat[0];
  sol_epoch_begin_lat = sol_lat[0];
  sol_dipcmp_val[sol_nrec] = MGH_matrix(sol_nnz*sol_nperdip,sol_ntime);
  sol_nrec++;
}
void
filter_recs(void)
{
  int i;

  for (i=0;i<sol_nrec;i++)
    filter_rec(i);
}

void
filter_rec(int np)
{
  int i,j,k,n,i0,i1;
  float sum;

  if (sol_trendflag)
    for (i=0;i<sol_nchan;i++)
      remove_trend(sol_Data[np][i],sol_ntime);
  for (i=0;i<sol_nchan;i++)
    for (j=sol_ntime;j<sol_ptime;j++)
      sol_Data[np][i][j] =
        sol_Data[np][i][sol_ntime-1]*
        (1-(j-(sol_ntime-1))/((float)sol_ptime-(sol_ntime-1)))+
        sol_Data[np][i][0]*(j-(sol_ntime-1))/((float)sol_ptime-(sol_ntime-1));

  bpfilter(sol_Data[np],
           sol_nchan,
           sol_ptime,
           sol_loflim*sol_ptime*sol_sample_period/1000.0,
           sol_hiflim*sol_ptime*sol_sample_period/1000.0);

  i1 = floor((sol_baseline_end-sol_epoch_begin_lat)/sol_sample_period+0.5);
  i0 = floor(i1-sol_baseline_period/sol_sample_period+0.5);
  if (i0<0) i0 = 0;
  for (k=0;k<sol_nchan;k++)
  {
    for (j=i0,n=0,sum=0;j<=i1;j++,n++)
      sum += sol_Data[np][k][j];
    sum /= n;
    for (j=0;j<sol_ptime;j++)
      sol_Data[np][k][j] -= sum;
  }
}

void
remove_trend(FLOATTYPE *v, int n)
{
  int k;
  double a,b,c1,c2,c11,c12,c21,c22,f;

  c1 = c2 = c11 = c12 = c21 = c22 = 0;
  for (k=0;k<n;k++)
  {
    f = v[k];
    c1 += f;
    c2 += k*f;
    c11 += k;
    c12 += 1;
    c21 += k*k;
    c22 += k;
  }
  a = (c1*c22-c2*c12)/(c11*c22-c12*c21);
  b = (c2*c11-c1*c21)/(c11*c22-c12*c21);
  for (k=0;k<n;k++)
    v[k] -= a*k+b;
  /*
    printf("remove trend: a=%f, b=%f\n",a,b);
  */
}

void
compute_timecourses(void)
{
  int i;

  for (i=0;i<sol_nrec;i++)
    compute_timecourse(i);
}

void
compute_timecourse(int num)
{
  int i,j,k,jc;
  double sum;

  printf("compute_timecourse(%d)\n",num);
  for (i=0;i<sol_ntime;i++)
  {
    for (j=0;j<sol_nnz;j++)
      for (jc=0;jc<sol_nperdip;jc++)
      {
        sum = 0;
        for (k=0;k<sol_nchan;k++)
          sum += sol_W[j*sol_nperdip+jc][k]*sol_Data[num][k][i];
        sol_dipcmp_val[num][j*sol_nperdip+jc][i] = sum;
      }
  }
}

void
plot_all_time_courses(void)
{
  int j,xpos,ypos;

  printf("plot_all_time_courses()\n");
  for (j=0;j<sol_nplotlist;j++)
  {
    xpos = PLOT_XFUDGE+(j%6)*(2*PLOT_XFUDGE+PLOT_XDIM);
    ypos = PLOT_YFUDGE+(j/6)*(PLOT_YFUDGE+PLOT_XFUDGE+PLOT_YDIM);
    plot_sol_time_course(sol_plotlist[j],PLOT_XDIM,PLOT_YDIM,xpos,ypos,j);
  }
}

void
read_plot_list(char *fname)
{
  /*
    AKL This differs from megsurfer.c
    The plot list has a label associated with each vertex number
  */

  int j,vnum;
  FILE *fp;

  printf("read_plot_list(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### can't find file %s\n",fname);
    PR return;
  }
  fscanf(fp,"%d",&sol_nplotlist);
  for (j=0;j<sol_nplotlist;j++)
  {
    fscanf(fp,"%d",&vnum);
    find_nearest_nonzero(vnum,&sol_plotlist[j]);
  }
  fclose(fp);
}

void
read_vertex_list(char *fname)
{
  int j;
  FILE *fp;

  printf("read_vertex_list(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### can't find file %s\n",fname);
    PR return;
  }
  fscanf(fp,"%d",&vertex_nplotlist);
  printf("  vertex_nplotlist = %d\n  ",vertex_nplotlist);
  for (j=0;j<vertex_nplotlist;j++)
  {
    fscanf(fp,"%d",&vertex_plotlist[j]);
    printf("%d ",vertex_plotlist[j]);
  }
  fclose(fp);
  printf("\n");
}

void
plot_nearest_nonzero(int vnum)
{
  int imin;

  find_nearest_nonzero(vnum,&imin);
  plot_sol_time_course(imin,0,0,0,0,0);
}

void
find_nearest_nonzero(int vnum, int *indxptr)
{
  int i,imin;
  float d,mindist;

  imin = 0 ;
  printf("find_nearest_nonzero(%d,*)\n",vnum);
  mindist = 1000000;
  for (i=0;i<sol_nnz;i++)
  {
    d = SQR(mris->vertices[sol_dipindex[i]].x-mris->vertices[vnum].x)+
        SQR(mris->vertices[sol_dipindex[i]].y-mris->vertices[vnum].y)+
        SQR(mris->vertices[sol_dipindex[i]].z-mris->vertices[vnum].z);
    if (d<mindist)
    {
      mindist=d;
      imin=i;
    }
  }
  printf("%d => %d (%d) (mindist = %f)\n",
         vnum,sol_dipindex[imin],imin,sqrt(mindist));
  *indxptr = imin;
}

void
find_nearest_fixed(int vnum, int *vnumptr)
{
  int i,imin= -1;
  float d,mindist;

  printf("find_nearest_fixed(%d,*)\n",vnum);
  mindist = 1000000;
  for (i=0;i<mris->nvertices;i++)
    if (mris->vertices[i].fixedval)
    {
      d = SQR(mris->vertices[i].x-mris->vertices[vnum].x)+
          SQR(mris->vertices[i].y-mris->vertices[vnum].y)+
          SQR(mris->vertices[i].z-mris->vertices[vnum].z);
      if (d<mindist)
      {
        mindist=d;
        imin=i;
      }
    }
  *vnumptr = imin;
  printf("selection = %d\n",imin);
}

void
plot_sol_time_course(int imin, int plot_xdim, int plot_ydim, int plot_xloc,
                     int plot_yloc, int plotnum)
{
  int i,k;
  float val;
  FILE *fp;
  char fname[STRLEN],command[STRLEN],cmd[STRLEN];

  printf("plot_sol_time_course(%d,%d,%d,%d,%d)\n",
         imin,plot_xdim,plot_ydim,plot_xloc,plot_yloc);
  sprintf(command,"xplot");
  for (i=0;i<sol_nrec;i++)
  {
    sprintf(fname,"/tmp/tmp-%d-%d.xplot",plotnum,i);
    sprintf(cmd,"%s %s",command,fname);
    sprintf(command,"%s",cmd);
    fp = fopen(fname,"w");
    if (fp==NULL)
    {
      printf("surfer: ### can't create file %s\n",fname);
      PR return;
    }

    fprintf(fp,"/color %d\n",i+1);

    fprintf(fp,"/plotname %d\n",i);
    fprintf(fp,"/ytitle %d\n",plotnum);
    if (plot_xdim>0)
      fprintf(fp,"/geometry %dx%d+%d+%d\n",
              plot_xdim,plot_ydim,plot_xloc,plot_yloc);
    for (k=0;k<sol_ntime;k++)
    {
      val = dipval(i,imin,k);
      fprintf(fp,"%f %f\n",sol_lat[k],val);
    }
    fclose(fp);
  }
  sprintf(cmd,"%s &\n",command);
  sprintf(command,"%s",cmd);
  system(command);
}

void
load_vals_from_sol(float tmid, float dt, int num)
{
  int j,k,ilat0,ilat1;
  double avgval,nsum,val,f;

  printf("load_vals_from_sol(%f,%f,%d)\n",tmid,dt,num);
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  ilat0 = floor((tmid-dt/2-sol_epoch_begin_lat)/sol_sample_period+0.5);
  ilat1 = floor((tmid+dt/2-sol_epoch_begin_lat)/sol_sample_period+0.5);
  printf("ilat0=%d, ilat1=%d\n",ilat0,ilat1);
  nsum = ilat1-ilat0+1;
  for (k=0;k<sol_nnz;k++)
  {
    avgval = 0;
    for (j=ilat0;j<=ilat1;j++)
    {
      val = dipval(num,k,j);
      avgval += val;
    }
    f = fabs(sol_pval[k]);
    if (sol_pslope!=0)
      f = tanh(sol_pslope*(f-sol_pthresh));
    else
      f = 1;
    /*
      if (f<0) f = 0;
    */
    mris->vertices[sol_dipindex[k]].val = f*avgval/nsum;
  }
}

void
load_var_from_sol(int num)
{
  int j,k,ilat2;
  double val,maxval,f;

  printf("load_var_from_sol(%d)\n",num);
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  ilat2 = floor((0-sol_epoch_begin_lat)/sol_sample_period+0.5);
  for (k=0;k<sol_nnz;k++)
  {
    maxval = 0;
    for (j=ilat2;j<sol_ntime;j++)
    {
      val = fabs(dipval(num,k,j));
      if (val>maxval) maxval = val;
    }
    f = fabs(sol_pval[k]);
    if (sol_pslope!=0)
      f = tanh(sol_pslope*(f-sol_pthresh));
    else
      f = 1;
    /*
      if (f<0) f = 0;
    */
    mris->vertices[sol_dipindex[k]].val = f*maxval;
  }
}

void
variance_ratio(void)
{
  int j,k,ilat2;
  double sum0,sum1,nsum,val;

  printf("variance_ratio()\n");
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  ilat2 = floor((0-sol_epoch_begin_lat)/sol_sample_period+0.5);
  for (k=0;k<sol_nnz;k++)
  {
    sum0 = sum1 = 0;
    nsum = 0;
    for (j=ilat2;j<sol_ntime;j++)
    {
      val = dipval(0,k,j);
      sum0 += SQR(val);
      val = dipval(1,k,j);
      sum1 += SQR(val);
      nsum++;
    }
    if (sum0==0 || sum1==0) sum0 = sum1 = 1;
    mris->vertices[sol_dipindex[k]].val = 10*log10(sum0/sum1);
  }
}

void
load_pvals(void)
{
  int k;
  double f;

  printf("load_pvals()\n");
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  for (k=0;k<sol_nnz;k++)
  {
    f = sol_pval[k];
    if (sol_rectpval)
      f = fabs(f);
    mris->vertices[sol_dipindex[k]].val = f;
  }
}

void
compute_pval_fwd(float pthresh)
{
#if 0
  int j,k;
  float f;

  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  for (j=0;j<sol_nchan;j++)
    sol_sensval[j] = 0;
  for (k=0;k<sol_nnz;k++)
  {
    f = sol_pval[k];
    if ((pthresh>=0 && f<pthresh) || (pthresh<0 && -f<-pthresh)) f = 0;
    if (f<0) f = -f;
    mris->vertices[sol_dipindex[k]].val = f;
    if (f>0)
    {
      if (sol_nperdip==1)
        for (j=0;j<sol_nchan;j++)
          sol_sensval[j] += f*sol_A[j][sol_nperdip*k];
      else
        for (j=0;j<sol_nchan;j++)
          sol_sensval[j] +=
            f*sol_A[j][sol_nperdip*k+0]*
            mris->vertices[sol_dipindex[k]].dipnx +
            f*sol_A[j][sol_nperdip*k+1]*
            mris->vertices[sol_dipindex[k]].dipny +
            f*sol_A[j][sol_nperdip*k+2]*
            mris->vertices[sol_dipindex[k]].dipnz;
    }
  }
#endif
}

void
compute_select_fwd(float maxdist)
{
#if 0
  int j,k;
  float dist;

  for (j=0;j<sol_nchan;j++)
    sol_sensval[j] = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  find_nearest_fixed(selection,&selection);
  for (k=0;k<sol_nnz;k++)
  {
    dist = sqrt(SQR(mris->vertices[sol_dipindex[k]].\
                    x-mris->vertices[selection].x)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    y-mris->vertices[selection].y)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    z-mris->vertices[selection].z));
    if (dist<=maxdist)
    {
      mris->vertices[sol_dipindex[k]].val = 1;
      if (sol_nperdip==1)
        for (j=0;j<sol_nchan;j++)
          sol_sensval[j] += sol_A[j][sol_nperdip*k];
      else
        for (j=0;j<sol_nchan;j++)
          sol_sensval[j] +=
            sol_A[j][sol_nperdip*k+0]*
            mris->vertices[sol_dipindex[k]].dipnx +
            sol_A[j][sol_nperdip*k+1]*
            mris->vertices[sol_dipindex[k]].dipny +
            sol_A[j][sol_nperdip*k+2]*
            mris->vertices[sol_dipindex[k]].dipnz;
    }
  }
#endif
}

void
compute_select_crosstalk(void)
{
#if 0
  int j,k,l,lc;
  float dist,maxdist=1e-10;
  double sum,xtalk,xtalk0,sumwtdist,sumwt,hvol;

  for (j=0;j<sol_nchan;j++)
    sol_sensval[j] = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  find_nearest_fixed(selection,&selection);
  sumwtdist = sumwt = hvol = 0;
  for (k=0;k<sol_nnz;k++)
  {
    dist = sqrt(SQR(mris->vertices[sol_dipindex[k]].\
                    x-mris->vertices[selection].x)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    y-mris->vertices[selection].y)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    z-mris->vertices[selection].z));
    if (dist<=maxdist)
    {
      for (l=0;l<sol_nnz;l++)
      {
        xtalk = 0;
        for (lc=0;lc<sol_nperdip;lc++)
        {
          sum = 0;
          for (j=0;j<sol_nchan;j++)
            sum +=
              sol_W[k*sol_nperdip+lc][j]*sol_A[j][sol_nperdip*l+lc];
          xtalk += sum*sum;
        }
        mris->vertices[sol_dipindex[l]].val = xtalk;
      }
      xtalk0 = mris->vertices[sol_dipindex[k]].val;
      for (l=0;l<sol_nnz;l++)
      {
        xtalk = mris->vertices[sol_dipindex[l]].val /= xtalk0;
        /* normalize by auto-crosstalk */
        dist = sqrt(SQR(mris->vertices[sol_dipindex[l]].\
                        dipx-mris->vertices[sol_dipindex[k]].dipx)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipy-mris->vertices[sol_dipindex[k]].dipy)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipz-mris->vertices[sol_dipindex[k]].dipz));
        sumwtdist += xtalk*dist;
        sumwt += xtalk;
        if (xtalk>0.5) hvol++;
      }
    }
  }
  printf("mean half-max area: %f, xtalk-weighted dist: %f,  mean xtalk: %f)\n",
         hvol/sol_nnz*mris->total_area,sumwtdist/sumwt,sqrt(sumwt/sol_nnz));
#endif
}

#if 0
void
compute_all_crosstalk(void)
{
  int i,j,k,l,lc;
  float dist,maxdist=1e-10;
  double sum,xtalk,xtalk0,sumwtdist,sumwt,hvol,sumhvol,sumxtalk,nsum;

  printf("compute_all_crosstalk()\n");

  for (j=0;j<sol_nchan;j++)
    sol_sensval[j] = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  sumhvol = nsum = 0;
  for (k=0;k<sol_nnz;k++)
  {
    printf("\rcountdown %4d  ",sol_nnz-k);
    fflush(stdout);
    sumwtdist = sumwt = hvol = 0;
    for (l=0;l<sol_nnz;l++)
    {
      xtalk = 0;
      for (lc=0;lc<sol_nperdip;lc++)
      {
        sum = 0;
        for (j=0;j<sol_nchan;j++)
          sum += sol_W[k*sol_nperdip+lc][j]*sol_A[j][sol_nperdip*l+lc];
        xtalk += sum*sum;
      }
      mris->vertices[sol_dipindex[l]].val2 = xtalk;
    }
    xtalk0 = mris->vertices[sol_dipindex[k]].val2;
    for (l=0;l<sol_nnz;l++)
    {
      xtalk = mris->vertices[sol_dipindex[l]].val2 /= xtalk0;
      /* normalize by auto-cross talk */
      sumwt += xtalk;
      if (xtalk>0.5) hvol++;
    }
    mris->vertices[sol_dipindex[k]].val = hvol/sol_nnz*mris->total_area;
    sumxtalk += sumwt/sol_nnz;
    sumhvol += hvol;
    nsum++;
  }
  printf("\n");
  sumhvol /= nsum;
  sumxtalk /= nsum;
  printf("mean half-max area: %f, mean xtalk: %f)\n",
         sumhvol/sol_nnz*mris->total_area,sqrt(sumxtalk));
}
#endif

void
compute_all_crosstalk(int weighttype)
{
#if 0
  int j,k,l,lc;
  float dist;
  double sum,xtalk,xtalk0,sumwtxtalk,sumwt,sumxtalk;

  printf("compute_all_crosstalk()\n");

  for (j=0;j<sol_nchan;j++)
    sol_sensval[j] = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  for (k=0;k<sol_nnz;k++)
  {
    printf("\rcountdown %4d  ",sol_nnz-k);
    fflush(stdout);
    sumwtxtalk = sumxtalk = sumwt = 0;
    for (l=0;l<sol_nnz;l++)
    {
      xtalk = 0;
      for (lc=0;lc<sol_nperdip;lc++)
      {
        sum = 0;
        for (j=0;j<sol_nchan;j++)
          sum += sol_W[k*sol_nperdip+lc][j]*sol_A[j][sol_nperdip*l+lc];
        xtalk += sum*sum;
      }
      mris->vertices[sol_dipindex[l]].val2 = xtalk;
    }
    if (weighttype==0)
    {
      xtalk0 = mris->vertices[sol_dipindex[k]].val2;
      for (l=0;l<sol_nnz;l++)
      {
        xtalk = mris->vertices[sol_dipindex[l]].val2 /= xtalk0;
        /* normalize by auto-crosstalk */
        wt = 1;
        sumwt += wt;
        sumxtalk += xtalk;
        sumwtxtalk += wt*xtalk;
      }
    }
    else if (weighttype==1)
    {
      xtalk0 = mris->vertices[sol_dipindex[k]].val2;
      for (l=0;l<sol_nnz;l++)
      {
        dist = sqrt(SQR(mris->vertices[sol_dipindex[l]].\
                        dipx-mris->vertices[sol_dipindex[k]].dipx)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipy-mris->vertices[sol_dipindex[k]].dipy)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipz-mris->vertices[sol_dipindex[k]].dipz));
        xtalk = mris->vertices[sol_dipindex[l]].val2 /= xtalk0;
        /* normalize by auto-crosstalk */
        wt = dist;
        sumwt += wt;
        sumxtalk += xtalk;
        sumwtxtalk += wt*xtalk;
      }
    }
    else if (weighttype==2)
    {
      xtalk0 = mris->vertices[sol_dipindex[k]].val2;
      for (l=0;l<sol_nnz;l++)
      {
        dist = sqrt(SQR(mris->vertices[sol_dipindex[l]].\
                        dipx-mris->vertices[sol_dipindex[k]].dipx)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipy-mris->vertices[sol_dipindex[k]].dipy)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipz-mris->vertices[sol_dipindex[k]].dipz));
        xtalk = mris->vertices[sol_dipindex[l]].val2 /= xtalk0;
        /* normalize by auto-crosstalk */
        wt = dist;
        sumwt += xtalk;
        sumxtalk += xtalk;
        sumwtxtalk += wt*xtalk;
      }
    }
    mris->vertices[sol_dipindex[k]].val = sumwtxtalk/sumwt;
  }
  printf("\n");
#endif
}

void
compute_select_pointspread(void)
{
#if 0
  int j,k,l,lc;
  float dist,maxdist=1e-10;
  double sum,xtalk,xtalk0,sumwtdist,sumwt,hvol;

  for (j=0;j<sol_nchan;j++)
    sol_sensval[j] = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  find_nearest_fixed(selection,&selection);
  sumwtdist = sumwt = hvol = 0;
  for (k=0;k<sol_nnz;k++)
  {
    dist = sqrt(SQR(mris->vertices[sol_dipindex[k]].\
                    x-mris->vertices[selection].x)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    y-mris->vertices[selection].y)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    z-mris->vertices[selection].z));
    if (dist<=maxdist)
    {
      for (l=0;l<sol_nnz;l++)
      {
        xtalk = 0;
        for (lc=0;lc<sol_nperdip;lc++)
        {
          sum = 0;
          for (j=0;j<sol_nchan;j++)
            sum +=
              sol_W[l*sol_nperdip+lc][j]*sol_A[j][sol_nperdip*k+lc];
          xtalk += sum*sum;
        }
        mris->vertices[sol_dipindex[l]].val = xtalk;
      }
      xtalk0 = mris->vertices[sol_dipindex[k]].val;
      for (l=0;l<sol_nnz;l++)
      {
        xtalk = mris->vertices[sol_dipindex[l]].val;
        dist = sqrt(SQR(mris->vertices[sol_dipindex[l]].\
                        dipx-mris->vertices[sol_dipindex[k]].dipx)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipy-mris->vertices[sol_dipindex[k]].dipy)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipz-mris->vertices[sol_dipindex[k]].dipz));
        sumwtdist += xtalk*dist*dist;
        xtalk /= xtalk0; /* normalize by auto-crosstalk */
        sumwt += xtalk;
        if (xtalk>0.5) hvol++;
        mris->vertices[sol_dipindex[l]].val = xtalk;
      }
    }
  }
  printf("mean half-max area: %f, mean pointspread: "
         "%f, pointspread-weighted dist2: %f)\n",
         hvol/sol_nnz*mris->total_area,sqrt(sumwt/sol_nnz),sumwtdist/sumwt);
#endif
}

void
compute_all_pointspread(void)
{
#if 0
  int i,j,k,l,lc;
  float dist2;
  double sum,xtalk,xtalk0,sumwtdist,sumwt,hvol,sumhvol,sumxtalk,totwtdist,nsum;

  printf("compute_all_crosstalk()\n");
  /*
    for (k=0;k<sol_nnz;k++)
    for (kc=0;kc<sol_nperdip;kc++)
    {
    sum = 0;
    for (i=0;i<sol_nchan;i++)
    for (l=0;l<sol_nnz;l++)
    for (lc=0;lc<sol_nperdip;lc++)
    sum += sol_W[k*sol_nperdip+kc][i]*sol_A[i][l*sol_nperdip+lc];
    printf("%d(%d): area = %f\n",k,kc,sum);
    }
  */
  printf("AW\n");
  for (i=0;i<10;i++)
  {
    for (j=0;j<10;j++)
    {
      sum = 0;
      for (l=0;l<sol_nnz;l++)
        for (lc=0;lc<sol_nperdip;lc++)
          sum += sol_A[i][l*sol_nperdip+lc]*sol_W[l*sol_nperdip+lc][j];
      printf("%9.1e",sum);
    }
    printf("\n");
  }
  printf("A\n");
  for (i=0;i<10;i++)
  {
    for (j=0;j<10;j++)
    {
      printf("%9.1e",sol_A[i][j]);
    }
    printf("\n");
  }
  printf("W\n");
  for (i=0;i<10;i++)
  {
    for (j=0;j<10;j++)
    {
      printf("%9.1e",sol_W[i][j]);
    }
    printf("\n");
  }
  return;

  for (j=0;j<sol_nchan;j++)
    sol_sensval[j] = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  sumhvol = nsum = totwtdist = 0;
  for (k=0;k<sol_nnz;k++)
  {
    printf("\rcountdown %4d  ",sol_nnz-k);
    fflush(stdout);
    sumwtdist = sumwt = hvol = 0;
    for (l=0;l<sol_nnz;l++)
    {
      xtalk = 0;
      for (lc=0;lc<sol_nperdip;lc++)
      {
        sum = 0;
        for (j=0;j<sol_nchan;j++)
          sum += sol_W[l*sol_nperdip+lc][j]*sol_A[j][sol_nperdip*k+lc];
        xtalk += sum*sum;
      }
      mris->vertices[sol_dipindex[l]].val2 = xtalk;
    }
    xtalk0 = mris->vertices[sol_dipindex[k]].val2;
    for (l=0;l<sol_nnz;l++)
    {
      dist2 =
        SQR(mris->vertices[sol_dipindex[l]].\
            dipx-mris->vertices[sol_dipindex[k]].dipx)+
        SQR(mris->vertices[sol_dipindex[l]].\
            dipy-mris->vertices[sol_dipindex[k]].dipy)+
        SQR(mris->vertices[sol_dipindex[l]].\
            dipz-mris->vertices[sol_dipindex[k]].dipz);
      xtalk = mris->vertices[sol_dipindex[l]].val2;
      sumwtdist += dist2*xtalk;
      xtalk /= xtalk0; /* normalize by auto-crosstalk */
      sumwt += xtalk;
      if (xtalk>0.5) hvol++;
    }
    mris->vertices[sol_dipindex[k]].val = hvol/sol_nnz*mris->total_area;
    sumxtalk += sumwt/sol_nnz;
    totwtdist += sumwtdist/sol_nnz;
    sumhvol += hvol;
    nsum++;
  }
  printf("\n");
  sumhvol /= nsum;
  sumxtalk /= nsum;
  totwtdist /= nsum;
  printf("mean half-max area: %f, mean pointspread: "
         "%f, mean wt. dist^2: %f)\n",
         sumhvol/sol_nnz*mris->total_area,sqrt(sumxtalk),totwtdist);
#endif
}

void
recompute_select_inverse(void)
{
#if 0
  int i,j,k,l,lc;
  float dist,maxdist=1e-10;
  double sum;
  /*
    if (sol_nperdip!=1) {printf("sol_nperdip must be 1!\n");return;}
  */
  find_nearest_fixed(selection,&selection);
  for (i=0;i<sol_nchan;i++)
    for (j=0;j<sol_nchan;j++)
      sol_M[i][j] = 0;
  for (k=0;k<sol_nnz;k++)
  {
    dist = sqrt(SQR(mris->vertices[sol_dipindex[k]].\
                    x-mris->vertices[selection].x)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    y-mris->vertices[selection].y)+
                SQR(mris->vertices[sol_dipindex[k]].\
                    z-mris->vertices[selection].z));
    if (dist<=maxdist)
    {
      for (l=0;l<sol_nnz;l++)
      {
        dist = sqrt(SQR(mris->vertices[sol_dipindex[l]].\
                        dipx-mris->vertices[sol_dipindex[k]].dipx)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipy-mris->vertices[sol_dipindex[k]].dipy)+
                    SQR(mris->vertices[sol_dipindex[l]].\
                        dipz-mris->vertices[sol_dipindex[k]].dipz));
        for (i=0;i<sol_nchan;i++)
          for (j=0;j<sol_nchan;j++)
            for (lc=0;lc<sol_nperdip;lc++)
              sol_M[i][j] +=
                sol_A[i][sol_nperdip*l+lc]*
                sol_A[j][sol_nperdip*l+lc]*sqrt(dist);
      }
      regularize_matrix(sol_M,sol_nchan,0.0001);
      inverse(sol_M,sol_Mi,sol_nchan);
      for (lc=0;lc<sol_nperdip;lc++)
      {
        for (i=0;i<sol_nchan;i++)
          sol_sensvec1[i] = sol_A[i][sol_nperdip*k+lc];
        vector_multiply(sol_Mi,
                        sol_sensvec1,
                        sol_sensvec2,
                        sol_nchan,
                        sol_nchan);
        sum = 0;
        for (i=0;i<sol_nchan;i++)
          sum += sol_sensvec1[i]*sol_sensvec2[i];
        if (sum!=0) sum = 1/sum;
        for (i=0;i<sol_nchan;i++)
          sol_W[sol_nperdip*k+lc][i] = sol_sensvec2[i]*sum;
      }
    }
  }
#endif
}

void
compute_pval_inv(void)
{
#if 0
  int j,k,jc;
  double sum,sum2;

  printf("compute_pval_inv()\n");
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].val = 0;
    if (!mris->vertices[k].fixedval)
      mris->vertices[k].undefval = TRUE;
  }
  if (sol_plot_type==0)
  {
    for (j=0;j<sol_nnz;j++)
    {
      sum2 = 0;
      for (jc=0;jc<sol_nperdip;jc++)
      {
        sum = 0;
        for (k=0;k<sol_nchan;k++)
          sum += sol_W[j*sol_nperdip+jc][k]*sol_sensval[k];
        sum2 += sum*sum;
      }
      mris->vertices[sol_dipindex[j]].val = sqrt(sum2);
    }
  }
  if (sol_plot_type==1 || sol_plot_type==2)
  {
    for (j=0;j<sol_nnz;j++)
    {
      sum2 = 0;
      for (jc=0;jc<sol_nperdip;jc++)
      {
        sum = 0;
        for (k=0;k<sol_nchan;k++)
          sum += sol_W[j*sol_nperdip+jc][k]*sol_sensval[k];
        if (sol_nperdip==1)
          sum2 += sum;
        else
        {
          if (jc==0)
            sum2 += sum*mris->vertices[sol_dipindex[j]].dipnx;
          else if (jc==1)
            sum2 += sum*mris->vertices[sol_dipindex[j]].dipny;
          else if (jc==2)
            sum2 += sum*mris->vertices[sol_dipindex[j]].dipnz;
        }
      }
      mris->vertices[sol_dipindex[j]].val = sum2;
    }
  }
  normalize_vals();
#endif
}

void
normalize_vals(void)
{
  int k;
  float val,maxval=0;

  for (k=0;k<mris->nvertices;k++)
  {
    val = fabs(mris->vertices[k].val);
    if (val>maxval) maxval = val;
  }
  if (maxval!=0)
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].val /= maxval;
}

void
normalize_time_courses(int normtype)
{
  double val,sum,sum2,nsum,noise_rms,maxval;
  int i,j,k,ic,ilat0,ilat1,ilat2;
  FLOATTYPE tmpvec1[2048];

  sum = maxval = 0.0f ;
  printf("normalize_time_courses(%d)\n",normtype);
  if (normtype==0)
  {
    for (i=0;i<sol_nnz;i++)
      sol_dipfact[i] = 1.0;
  }
  else if (normtype==1) /* normalize for baseline noise covariance */
  {
    ilat1 = floor((sol_baseline_end-sol_epoch_begin_lat)/
                  sol_sample_period+0.5);
    ilat0 = floor(ilat1-sol_baseline_period/sol_sample_period+0.5);
    ilat2 = floor((0-sol_epoch_begin_lat)/sol_sample_period+0.5);
    for (k=0;k<sol_nnz;k++)
    {
      sol_dipfact[k] = 1;
      sum2 = 0;
      nsum = 0;
      for (i=0;i<sol_nrec;i++)
      {
        for (j=ilat0;j<=ilat1;j++)
        {
          val = dipval(i,k,j);
          sum2 += SQR(val);
          val = fabs(val);
          if (val>maxval) maxval=val;
          nsum++;
        }
      }
      noise_rms = sqrt(sum2/nsum);
      maxval = 0;
      for (i=0;i<sol_nrec;i++)
      {
        for (j=ilat2;j<sol_ntime;j++)
        {
          val = dipval(i,k,j);
          val = fabs(val);
          if (val>maxval) maxval=val;
        }
      }
      sol_dipfact[k] = 1.0/(noise_rms);
    }
  } else if (normtype==2) /* normalize for white noise w.w = 1 */
  {
    for (i=0;i<sol_nnz;i++)
    {
      sum2 = 0;
      for (ic=0;ic<sol_nperdip;ic++)
      {
        for (j=0;j<sol_nchan;j++)
        {
          sum = 0;
          for (k=0;k<sol_nchan;k++)
            sum += SQR(sol_W[i*sol_nperdip+ic][k]);
        }
        sum2 += sum;
      }
      noise_rms = sqrt(sum2);
      sol_dipfact[i] = 1.0/noise_rms;
    }
  } else if (normtype==3) /* normalize by read noise covariance matrix */
  {
    for (i=0;i<sol_nnz;i++)
    {
      sum2 = 0;
      for (ic=0;ic<sol_nperdip;ic++)
      {
        for (j=0;j<sol_nchan;j++)
        {
          sum = 0;
          for (k=0;k<sol_nchan;k++)
            sum +=
              sol_W[i*sol_nperdip+ic][k]*sol_NoiseCovariance[j][k];
          tmpvec1[j] = sum;
        }
        sum = 0;
        for (j=0;j<sol_nchan;j++)
          sum += tmpvec1[j]*sol_W[i*sol_nperdip+ic][j];
        sum2 += sum;
      }
      noise_rms = sqrt(sum2);
      sol_dipfact[i] = 1.0/noise_rms;
      /*
        printf("sol_dipfact[%d] = %e\n",i,sol_dipfact[i]);
      */
    }
  } else if (normtype==4) /* "Unity Gain  wa = 1 */
  {
    for (i=0;i<sol_nnz;i++)
    {
      for (ic=0;ic<sol_nperdip;ic++)
      {
        sum = 0;
        for (k=0;k<sol_nchan;k++)
          sum +=
            sol_W[i*sol_nperdip+ic][k]*sol_A[k][i*sol_nperdip+ic];
        if (sum!=0)
          sum = 1/sum;
        for (k=0;k<sol_nchan;k++)
          sol_W[i*sol_nperdip+ic][k] *= sum;
      }
      sol_dipfact[i] = 1.0;
    }
  }
  for (i=0;i<sol_nnz;i++)
  {
    for (ic=0;ic<sol_nperdip;ic++)
      for (k=0;k<sol_nchan;k++)
        sol_W[i*sol_nperdip+ic][k] *= sol_dipfact[i];
  }
}

void
normalize_inverse(void)
{
  int j,k,jc;
  double sum;

  printf("normalize_inverse()\n");
  for (j=0;j<sol_nnz;j++)
  {
    sum = 0;
    for (jc=0;jc<sol_nperdip;jc++)
      for (k=0;k<sol_nchan;k++)
        sum += SQR(sol_W[j*sol_nperdip+jc][k]);
    sum = sqrt(sum);
    for (jc=0;jc<sol_nperdip;jc++)
      for (k=0;k<sol_nchan;k++)
        sol_W[j*sol_nperdip+jc][k] /= sum;
  }
}

void
setsize_window(int pix)
{
  if (openglwindowflag)
  {
    printf("surfer: ### setsize_window failed: gl window already open\n");
    PR return;
  }

  frame_xdim = (long)pix;
  frame_ydim = (long)pix;
}

void
resize_window(int pix)
{
  if (!openglwindowflag)
  {
    printf("surfer: ### resize_window failed: no gl window open\n");
    PR return;
  }

#ifdef OPENGL
  if (renderoffscreen1)
  {
    printf("surfer: ### resize_window failed: can't resize offscreen win\n");
    PR printf("surfer: ### use setsize_window <pix> before open_window\n");
    PR return;
  }

#ifndef USE_XGLUT_WINDOW
  if (pix>0)
  {  /* command line (not mouse) resize */
    XResizeWindow(xDisplay, w.wMain, pix, pix);
    if (TKO_HAS_OVERLAY(w.type))
      XResizeWindow(xDisplay, w.wOverlay, pix, pix);
    w.w = w.h = pix;
  }
#endif /* USE_XGLUT_WINDOW */

  frame_xdim = w.w;
  frame_ydim = w.h;
  glViewport(0, 0, frame_xdim, frame_ydim);
  resize_buffers(frame_xdim, frame_ydim);

#else
  if (pix>0)
  {   /* tcl: zeropix -> flag to reshape after mouse resize */
    prefposition(0,(short)pix,0,(short)pix);
    winconstraints();
    reshapeviewport();
    getsize(&frame_xdim,&frame_ydim);
    keepaspect(1,1);
    winconstraints();  /* call again to keep resizable */
    resize_buffers(frame_xdim,frame_ydim);
  }
  else
  {  /* tcl: zeropix flag->reshape w/mouse resize (REDRAW event) */
    reshapeviewport();
    getsize(&frame_xdim,&frame_ydim);
    resize_buffers(frame_xdim,frame_ydim);
  }
#endif
}

void
save_rgb(char *fname)
{
#if 1
  unsigned short *red, *blue, *green ;
  int             width, height, size ;

  width = (int)frame_xdim;
  height = (int)frame_ydim;
  size = width*height;
  red = (unsigned short *)calloc(size, sizeof(unsigned short)) ;
  green = (unsigned short *)calloc(size, sizeof(unsigned short)) ;
  blue = (unsigned short *)calloc(size, sizeof(unsigned short)) ;
  grabPixels(frame_xdim, frame_ydim, red, green, blue) ;
  save_rgbfile(fname, frame_xdim, frame_ydim, red, green, blue) ;
  free(blue) ;
  free(red) ;
  free(green) ;
#else
  if (!openglwindowflag)
  {
    printf("surfer: ### save_rgb failed: no gl window open\n");
    PR return;
  }

  if (renderoffscreen1)
    pix_to_rgb(fname);
  else
  {
    if (scrsaveflag)
    {
      scrsave_to_rgb(fname);
    }
    else
    {
      pix_to_rgb(fname);
    }
  }
#endif
}

static void
move_window(int x,int y)
{
#ifdef OPENGL
  if (openglwindowflag)
  {
    XMoveWindow(xDisplay, w.wMain, x, y);
    w.x = x;
    w.y = y;
  }
  else if (!initpositiondoneflag)
  {
    tkoInitPosition(x,y,frame_xdim,frame_ydim);
    initpositiondoneflag = TRUE;
  }
  else ;
#endif
}

void
save_rgb_num(char *dir)
{
  char fname[NAME_LENGTH];
  FILE *fp;

  if (!openglwindowflag)
  {
    printf("surfer: ### save_rgb_num failed: no gl window open\n");
    PR return;
  }

  fp = stdin;
  framenum = 0;
  while (fp!=NULL && framenum<999)
  {
    framenum++;
    sprintf(fname,"%s/im-%03d.rgb",dir,framenum);
    fp = fopen(fname,"r");
    if (fp!=NULL) fclose(fp);
  }

  if (renderoffscreen1)
    pix_to_rgb(fname);
  else
  {
    if (scrsaveflag)
    {
      scrsave_to_rgb(fname);
    }
    else
    {
      pix_to_rgb(fname);
    }
  }
}

void
save_rgb_named_orig(char *dir, char *name)
{
  char fname[NAME_LENGTH];

  if (!openglwindowflag)
  {
    printf("surfer: ### save_rgb_named_orig failed: no gl window open\n");
    PR return;
  }

  sprintf(fname,"%s/%s",dir,name);

  if (renderoffscreen1)
    pix_to_rgb(fname);
  else
  {
    if (scrsaveflag)
    {
      scrsave_to_rgb(fname);
    }
    else
    {
      pix_to_rgb(fname);
    }
  }
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
  x0 = (int)xorig;
  x1 = (int)(xorig+xsize-1);
  y0 = (int)yorig;
  y1 = (int)(yorig+ysize-1);
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fclose(fp);
  sprintf(command,"scrsave %s %d %d %d %d\n",fname,x0,x1,y0,y1);
  system(command);
  printf("surfer: file %s written\n",fname);
  PR
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

  width = (int)frame_xdim;
  height = (int)frame_ydim;
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
  glReadPixels(0, 0, width, height,
               GL_GREEN,GL_UNSIGNED_SHORT,(GLvoid *)green);
  glReadPixels(0, 0, width, height, GL_BLUE,GL_UNSIGNED_SHORT, (GLvoid *)blue);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes);
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst);
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength);
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows);
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels);
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment);

  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fclose(fp);
#ifdef IRIX
  image = iopen(fname,"w",RLE(1), 3, width, height, 3);
#else
  image = iopen(fname,"w",UNCOMPRESSED(1), 3, width, height, 3);
#endif
  if (!image)
    return ;
  for (y = 0 ; y < height; y++)
  {
    r = red + y * width;
    g = green + y * width;
    b = blue + y * width;
    putrow(image, r, y, 0);
    putrow(image, g, y, 1);
    putrow(image, b, y, 2);
  }
  iclose(image);
  free(red);
  free(green);
  free(blue);

  printf("surfer: file %s written\n",fname);
  PR
#else
  printf("surfer: ### pix_to_rgb implemented only in OpenGL version\n");
  PR return;
#endif
}

void
read_annotations(char *fname)
{
  /* begin rkt */
  /* now handled as importing into multiple labels. */
#if 1
  labl_import_annotation (fname);
#else
  if (MRISreadAnnotation(mris, fname) != NO_ERROR)
    return ;
  annotationloaded = TRUE ;
  surface_compiled = 0 ;
#endif
  /* end rkt */
}
void
read_annotated_image(char *fpref, int frame_xdim, int frame_ydim)
{
  char fname[NAME_LENGTH],command[NAME_LENGTH*2];
  FILE *fp;
  int i,j,k;

  sprintf(command,"tobin %s.rgb %s.bin\n",fpref,fpref);
  system(command);
  sprintf(fname,"%s.bin",fpref);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  k = fread(binbuff,3,frame_xdim*frame_ydim,fp);
  if (k != frame_xdim*frame_ydim)
  {
    printf("surfer: ### %s.rgb does not match current window size\n",fpref);
    PR return;
  }
  fclose(fp);
  sprintf(command,"rm %s.bin\n",fpref);
  system(command);
  for (i=0;i<frame_ydim;i++)
    for (j=0;j<frame_xdim;j++)
    {
      /*
        framebuff[i*frame_xdim+j] = (binbuff[(i*frame_xdim+j)*3+2]<<16) |
        (binbuff[(i*frame_xdim+j)*3+1]<<8) |
        (binbuff[(i*frame_xdim+j)*3+0]);
      */
      framebuff[i*frame_xdim+j] =
        (binbuff[(i*frame_xdim+j)+2*(frame_ydim*frame_xdim)]<<16) |
        (binbuff[(i*frame_xdim+j)+1*(frame_ydim*frame_xdim)]<<8) |
        (binbuff[(i*frame_xdim+j)+0]);
    }
  for (i=0;i<frame_ydim;i++)
    for (j=0;j<frame_xdim;j++)
    {
      binbuff[(i*frame_xdim+j)*3+2] =
        (framebuff[i*frame_xdim+j]&0xff0000)>>16;
      binbuff[(i*frame_xdim+j)*3+1] = (framebuff[i*frame_xdim+j]&0xff00)>>8;
      binbuff[(i*frame_xdim+j)*3+0] = (framebuff[i*frame_xdim+j]&0xff);
    }
  sprintf(fname,"%s.bin~",fpref);
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fwrite(binbuff,3,frame_xdim*frame_ydim,fp);
  fclose(fp);
  sprintf(command,"frombin %s.bin~ %s.rgb~ %d %d 3 -i\n",
          fpref,fpref,frame_xdim,frame_ydim);
  system(command);
  sprintf(command,"rm %s.bin~\n",fpref);
  system(command);
}

int diff(int a,int b)
{
  int d;
  unsigned char a1,a2,a3,b1,b2,b3,d1,d2,d3;

  a1 = ((a>>16)&0xff);
  a2 = ((a>>8)&0xff);
  a3 = (a&0xff);
  b1 = ((b>>16)&0xff);
  b2 = ((b>>8)&0xff);
  b3 = (b&0xff);
  d1 = a1-b1;
  d2 = a2-b2;
  d3 = a3-b3;
  d = (((int)d1)<<16)|(((int)d2)<<8)|((int)d3);
  return d;
}

void
save_rgb_cmp_frame(char *dir, int ilat)    /* open/close version */
{
  int numframes,file_xdim,file_ydim;
  float flat;
  static float firstflat;
  static long beginlastframe;
  long xsize,ysize;
  short x0,y0,x1,y1;
  char fname[NAME_LENGTH];
  FILE *fp;

  if (renderoffscreen1)
  {
    printf("surfer: ### save_rgb_cmp_frame failed: TODO: offscreen\n");
    PR return;
  }

  getsize(&xsize,&ysize);
  x0 = 0;
  y0 = 0;
  x1 = xsize-1;
  y1 = ysize-1;
  flat = (float)ilat;

  if (!openglwindowflag)
  {
    printf("surfer: ### save_rgb_cmp_frame failed: no gl window open\n");
    PR return;
  }

  /* open next avail new file, else append (this session), get numframes */
  numframes = -1;  /* last==first: n frames written = n-1 images */
  if (!openedcmpfilenum)
  {
    fp = stdin;
    while (fp!=NULL && openedcmpfilenum<999)
    {
      openedcmpfilenum++;
      sprintf(fname,"%s/movie-%03d.cmp",dir,openedcmpfilenum);
      fp = fopen(fname,"r");
      if (fp!=NULL) fclose(fp);
    }
    fp = fopen(fname,"w");
    if (fp==NULL)
    {
      printf("surfer: ### can't create file %s\n",fname);
      PR openedcmpfilenum = FALSE;
      return;
    }
    printf("surfer: opening new compressed movie:\n");
    printf("surfer:   %s\n",fname);
    numframes = 0;
    fwrite2((int)xsize,fp);
    fwrite2((int)ysize,fp);
    fwrite2(numframes,fp);
    lrectread(x0,y0,x1,y1,framebuff3);  /* save for later */
    firstflat = flat;
  }
  else
  { /* get frames written, reopen to append */
    sprintf(fname,"%s/movie-%03d.cmp",dir,openedcmpfilenum);
    fp = fopen(fname,"r");
    if (fp==NULL)
    {
      printf("surfer: ### File %s not found\n",fname);
      PR openedcmpfilenum = FALSE;
      return;
    }
    fread2(&file_xdim,fp); /* ignored */
    fread2(&file_ydim,fp); /* ignored */
    fread2(&numframes,fp);
    fclose(fp);
    fp = fopen(fname,"r+"); /* "r+": overwrite ("a": no backseek;"w": trunc) */
    if (fp==NULL)
    {
      printf("surfer: ### can't create file %s\n",fname);
      PR return;
    }
  }

  /* write current frame over previous last (=first) frame */
  if (numframes > 0)
    fseek(fp,beginlastframe,SEEK_SET);
#if 0
  fwrite(&flat,1,sizeof(float),fp);
#else
  fwriteFloat(flat, fp) ;
#endif
  lrectread(x0,y0,x1,y1,framebuff);
  do_rgb_cmp_frame(xsize,ysize,fp);
  lrectread(x0,y0,x1,y1,framebuff2);
  numframes++;

  /* copy saved first frame (framebuff3) to last */
  beginlastframe = ftell(fp);
#if 0
  fwrite(&firstflat,1,sizeof(float),fp);
#else
  fwriteFloat(firstflat, fp) ;
#endif
  memmove(framebuff,framebuff3,frame_xdim*frame_ydim*sizeof(long));
  do_rgb_cmp_frame(xsize,ysize,fp);
  /* don't save first=last frame to framebuff2 (because backup/overwrite) */
  fclose(fp);

  /* update num frames written */
  fp = fopen(fname,"r+");
  fseek(fp,4L,SEEK_SET);
  fwrite2(numframes,fp);
  fseek(fp,0L,SEEK_END);
  fclose(fp);
  printf("surfer: %s: frame %d done\n",fname,numframes+1);
  PR
}

/* open file return if already open */
void
open_rgb_cmp_named(char *dir, char *name)
{
  char str[NAME_LENGTH];

  if (renderoffscreen1)
  {
    printf("surfer: ### open_rgb_cmp_named failed: TODO: offscreen\n");
    PR return;
  }
  if (!openglwindowflag)
  {
    printf("surfer: ### open_rgb_cmp_named failed: no gl window open\n");
    PR return;
  }
  if (cmpfilenamedframe != -1)
  {
    printf("surfer: ### named compressed rgb movie already open\n");
    PR return;
  }
  else
  {
    sprintf(str,"%s/%s",dir,name);
    fpcmpfilenamed = fopen(str,"w");
    if (fpcmpfilenamed==NULL)
    {
      printf("surfer: ### can't create file %s\n",str);
      PR return;
    }
    cmpfilenamedframe = 0;
    printf("surfer: opening new compressed movie:\n");
    printf("surfer:   %s\n",str);
    PR
  }
  resize_buffers(frame_xdim, frame_ydim);
}

void
save_rgb_cmp_frame_named(float lat)
{
  long xsize,ysize;
  short x0,y0,x1,y1;

  if (renderoffscreen1)
  {
    printf("surfer: ### save_rgb_cmp_frame_named failed: TODO: offscreen\n");
    PR return;
  }

  getsize(&xsize,&ysize);
  x0 = 0;
  y0 = 0;
  x1 = xsize-1;
  y1 = ysize-1;

  if (cmpfilenamedframe == -1)
  {
    printf("surfer: ### can't write frame: cmp movie file not opened yet\n");
    PR return;
  }
  if (cmpfilenamedframe == 0)
  {  /* save 1st framebuff for end; write header */
    lrectread(x0,y0,x1,y1,framebuff3);
    fwrite2((int)xsize,fpcmpfilenamed);
    fwrite2((int)ysize,fpcmpfilenamed);
    fwrite2(cmpfilenamedframe,fpcmpfilenamed);  /* overwritten at end */
    cmpfilenamedfirstlat = lat;
  }

  /* write compressed frame */
#if 0
  fwrite(&lat,1,sizeof(float),fpcmpfilenamed);
#else
  fwriteFloat(lat, fpcmpfilenamed) ;
#endif
  lrectread(x0,y0,x1,y1,framebuff);
  do_rgb_cmp_frame(xsize,ysize,fpcmpfilenamed);
  lrectread(x0,y0,x1,y1,framebuff2);
  cmpfilenamedframe++;
  printf("surfer: cmp movie frame %d (lat=%3.3f) written\n",
         cmpfilenamedframe,lat);
  PR
}

void
do_rgb_cmp_frame(long xsize,long ysize, FILE *fp)
{
  int i,lo,hi,pos,change;
  long32 fb1, fb2 ;

  lo = 0;
  hi = 0;
  pos = 0;
  while (lo<xsize*ysize)
  {
    fb1 = orderLong32Bytes(framebuff[hi])&0xffffff ;
    fb2 = orderLong32Bytes(framebuff2[hi])&0xffffff ;
    while ((hi<xsize*ysize)&&
           (hi-lo<32767)&&
           (fb1==fb2))
    {
      hi++;
      if (hi < xsize*ysize)
      {
        fb1 = orderLong32Bytes(framebuff[hi])&0xffffff ;
        fb2 = orderLong32Bytes(framebuff2[hi])&0xffffff ;
      }
    }
    if (hi>lo)
    {
      i = ((hi-lo)|0x8000);
      fwrite2(i,fp);
      pos += 2;
    }
    else
    {
      fb1 = orderLong32Bytes(framebuff[lo])&0xffffff ;
      fb2 = orderLong32Bytes(framebuff2[lo])&0xffffff ;
      change = diff(fb1,fb2);
      fb1 = orderLong32Bytes(framebuff[hi])&0xffffff ;
      fb2 = orderLong32Bytes(framebuff2[hi])&0xffffff ;
      while ((hi<xsize*ysize)&&
             (hi-lo<32767)&&
             (diff(fb1,fb2)==change))
      {
        hi++;
        if (hi < xsize*ysize)
        {
          fb1 = orderLong32Bytes(framebuff[hi])&0xffffff ;
          fb2 = orderLong32Bytes(framebuff2[hi])&0xffffff ;
        }
      }
      i = hi-lo;
      fwrite2(i,fp);
      fwrite3(change,fp);
      pos += 5;
    }
    lo = hi;
  }
}

void
close_rgb_cmp_named(void)      /* load first; save; write total frames */
{
  long xsize,ysize;
  short x0,y0,x1,y1;

  if (cmpfilenamedframe == -1)
  {
    printf("surfer: ### can't close_rgb_cmp_named: file not opened yet\n");
    PR return;
  }

  getsize(&xsize,&ysize);
  x0 = 0;
  y0 = 0;
  x1 = xsize-1;
  y1 = ysize-1;

  /* write last compressed frame (same as first) */
#if 0
  fwrite(&cmpfilenamedfirstlat,1,sizeof(float),fpcmpfilenamed);
#else
  fwriteFloat(cmpfilenamedfirstlat,fpcmpfilenamed) ;
#endif
  memmove(framebuff,framebuff3,frame_xdim*frame_ydim*sizeof(long));
  do_rgb_cmp_frame(xsize,ysize,fpcmpfilenamed);
  fseek(fpcmpfilenamed,4L,SEEK_SET);
  fwrite2(cmpfilenamedframe,fpcmpfilenamed);
  fseek(fpcmpfilenamed,0L,SEEK_END);  /* else file truncated! */
  fclose(fpcmpfilenamed);
  printf("surfer: closed compressed movie file\n");
  PR cmpfilenamedframe = -1;
}

void
redraw(void)
{
  int i,navg;

  if (!openglwindowflag)
  {
    printf("surfer: ### redraw failed: no gl window open\n");
    PR return;
  }

  if (overlayflag)
  {
    redraw_overlay();
    return ;
  }

  czclear(BACKGROUND,getgconfig(GC_ZMAX));

  dipscale = 0;
  dipavg = dipvar = logaratavg = logaratvar = logshearavg = logshearvar = 0;
  navg = 0;
  for (i=0;i<mris->nvertices;i++)
    if (!mris->vertices[i].ripflag)
    {
      if (fabs(mris->vertices[i].curv)>dipscale)
        dipscale=fabs(mris->vertices[i].curv);
      dipavg += mris->vertices[i].curv;
      dipvar += SQR(mris->vertices[i].curv);
#if 0
      logaratavg += mris->vertices[i].logarat;
      logaratvar += SQR(mris->vertices[i].logarat);
      logshearavg += mris->vertices[i].logshear;
      logshearvar += SQR(mris->vertices[i].logshear);
#endif
      navg++;
    }
  dipavg /= navg;
  dipvar = sqrt(dipvar/navg - dipavg*dipavg);
  logaratavg /= navg;
  logaratvar = sqrt(logaratvar/navg - logaratavg*logaratavg);
  logshearavg /= navg;
  logshearvar = sqrt(logshearvar/navg - logshearavg*logshearavg);

  draw_surface();

  if (selection>=0) draw_cursor(selection,TRUE);

  /* begin rkt */
  draw_marked_vertices ();
  /* end rkt */

  if (doublebufferflag) swapbuffers();

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
    printf("surfer: dipscale=%f, dipavg=%f, dipvar=%f\n",
           dipscale,dipavg,dipvar);
  dipscale = (dipscale!=0)?1/dipscale:1.0;
#if 0
  if (areaflag)
  {
    printf("surfer: logaratavg=%f, logaratvar=%f\n",logaratavg,logaratvar);
    PR
  }
  if (flag2d)
  {
    printf("surfer: logshearavg=%f, logshearvar=%f\n",logshearavg,logshearvar);
    PR
  }
#endif
}

void
redraw_second(void)
{
  int i,navg;

  if (!openglwindowflag)
  {
    printf("surfer: ### redraw_second failed: no gl window open\n");
    PR return;
  }
  if (!secondsurfaceloaded)
  {
    printf("surfer: ### redraw_second failed: no second surface read\n");
    PR return;
  }

  czclear(BACKGROUND,getgconfig(GC_ZMAX));

  dipavg2 = 0;
  navg = 0;
  for (i=0;i<mris2->nvertices;i++)
    if (!mris2->vertices[i].ripflag)
    {
      dipavg2 += mris2->vertices[i].curv;
      navg++;
    }
  dipavg2 /= navg;

  draw_second_surface();
  if (doublebufferflag) swapbuffers();
}
void
blinkbuffers(void)
{
  if (blinkflag)
  {
    if (blinkdelay<0)
    {
      if (doublebufferflag) swapbuffers();
#ifdef Irix
      sginap(blinktime);
#else
      sleep(blinktime);
#endif
    }
    else
      blinkdelay--;
  }
}

void
redraw_overlay(void)
{
  int i;
  float curvavg;

  czclear(BACKGROUND,getgconfig(GC_ZMAX));

  curvavg = 0;
  for (i=0;i<mris->nvertices;i++)
    curvavg += mris->vertices[i].curv;
  curvavg /= mris->nvertices;
  if (avgflag)
  {
    for (i=0;i<mris->nvertices;i++)
      mris->vertices[i].curv -= curvavg;
  }
  if (autoscaleflag)
  {
    dipscale = 0;
    for (i=0;i<mris->nvertices;i++)
      if (fabs(mris->vertices[i].val)>dipscale)
        dipscale=fabs(mris->vertices[i].val);
    printf("surfer: dipscale=%f\n",dipscale);
    PR dipscale = (dipscale!=0)?1/dipscale:1.0;
  }
  curvscale = 0;
  for (i=0;i<mris->nvertices;i++)
    if (fabs(mris->vertices[i].curv)>curvscale)
      curvscale=fabs(mris->vertices[i].curv);
  curvscale = (curvscale!=0)?1/curvscale:1.0;

  if (selection>=0) draw_cursor(selection,TRUE);

  draw_surface();

  if (doublebufferflag) swapbuffers();

  /*printf("curvscale=%f, curvavg=%f\n",curvscale,curvavg);PR*/
}


void
draw_cursor(int vindex,int onoroff)
{
  /* begin rkt */
# if 1
  /* don't draw if our flag isn't on. */
  if (drawcursorflag)
  {
    if (onoroff)
    {
      RGBcolor (0, 255, 255);
      draw_vertex_hilite (vindex);
    }
    else
    {
      set_color (0, 0, GREEN_RED_CURV);
      draw_vertex_hilite (vindex);
    }
  }

  return;

#else

  int i,k,n;
  FACE *f;
  VERTEX *v,*vselect;

  /* offscreen render opengl bug: RGBcolor(white) => surrounding faces color */

  if ((vindex > mris->nvertices)||(vindex<0))
  {
    if (vindex != -1)
      printf ("surfer: ### vertex index %d out of bounds\n",vindex);
    return;
  }
  vselect = &mris->vertices[vindex];
  if (onoroff==FALSE)
    set_color(0.0,0.0,GREEN_RED_CURV);
  else
  {
    if (blackcursorflag)  RGBcolor(0,0,0);
    else                  RGBcolor(0,255,255);
  }

  for (i=0;i<vselect->vnum;i++)
  {
    v = &mris->vertices[vselect->v[i]];
    linewidth(CURSOR_LINE_PIX_WIDTH);
    bgnline();
    load_brain_coords(v->nx,v->ny,v->nz,v1);
    n3f(v1);
    load_brain_coords(v->x+cup*v->nx,v->y+cup*v->ny,v->z+cup*v->nz,v1);
    v3f(v1);
    load_brain_coords(vselect->nx,vselect->ny,vselect->nz,v1);
    n3f(v1);
    load_brain_coords(vselect->x+cup*v->nx,
                      vselect->y+cup*v->ny,
                      vselect->z+cup*v->nz,v1);
    v3f(v1);
    endline();
  }

  if (bigcursorflag)
  {
    for (k=0;k<vselect->num;k++)
    {
      bgnquadrangle();
      f = &mris->faces[vselect->f[k]];
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        v = &mris->vertices[f->v[n]];
        load_brain_coords(v->nx,v->ny,v->nz,v1);
        n3f(v1);
        load_brain_coords(v->x+cup*v->nx,v->y+cup*v->ny,v->z+cup*v->nz,v1);
        v3f(v1);
      }
      endpolygon();
    }
  }
  glFlush() ;
#endif
  /* end rkt */
}

void
draw_all_cursor(void)
{
  int j;

  for (j=0;j<sol_nplotlist;j++)
  {
    draw_cursor(sol_dipindex[sol_plotlist[j]],TRUE);
  }
} /*end draw_all_cursor*/

void
draw_all_vertex_cursor(void)
{
  int j;

  for (j=0;j<vertex_nplotlist;j++)
  {
    draw_cursor(vertex_plotlist[j],TRUE);
  }
} /*end draw_all_vertex_cursor*/

void
clear_all_vertex_cursor(void)
{
  int j;

  for (j=0;j<vertex_nplotlist;j++)
  {
    draw_cursor(vertex_plotlist[j],FALSE);
  }
} /*end clear_all_vertex_cursor*/


void
invert_vertex(int vno)
{
  VERTEX *v ;
  int    n ;

  v = &mris->vertices[vno] ;
  v->nx *= -1.0f ;
  v->ny *= -1.0f ;
  v->nz *= -1.0f ;
  for (n = 0 ; n < v->num ; n++)
  {
    int fno = v->f[n];
    FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
    setFaceNorm(mris, fno, -fNorm->nx, -fNorm->ny, -fNorm->nz);
  }
}

void
invert_face(int fno)
{
  if (fno < 0 || fno >= mris->nfaces)
    return ;
  FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
  setFaceNorm(mris, fno, -fNorm->nx, -fNorm->ny, -fNorm->nz);
}

void
orient_sphere(void)
{
  int    vno, n ;
  VERTEX *v ;

  mris->status = MRIS_SPHERE ;
  MRIScomputeMetricProperties(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->stat = v->curv ;
    v->curv = 0.0 ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    for (n = 0 ; n < v->num ; n++)
      if (mris->faces[v->f[n]].area < 0)
        v->curv = -1 ;
  }
}

static void
load_gcsa(char *fname)
{
  Ggcsa = GCSAread(fname) ;
  if (Ggcsa == NULL)
    ErrorPrintf(ERROR_NOFILE, "load_gcsa(%s): could not read file",
                fname) ;
}
void
dump_faces(int vno)
{
  VERTEX const * const v = &mris->vertices[vno] ;

  int n ;
  for (n = 0 ; n < v->num ; n++)
  {
    FACE               const * f     = &mris->faces[     v->f[n]];
    FaceNormCacheEntry const * fNorm = getFaceNorm(mris, v->f[n]);
    fprintf(stderr,
            "face %6d [%d,%d,%d], n = (%2.1f,%2.1f,%2.1f), area = %2.1f\n",
            v->f[n], f->v[0],f->v[1],f->v[2],fNorm->nx, fNorm->ny, fNorm->nz, f->area) ;
  }
}

void
dump_vertex(int vno)
{
  VERTEX *v ;

  v = &mris->vertices[vno] ;
  fprintf(stderr, "v %d, n = (%2.2f, %2.2f, %2.2f), area = %2.2f\n",
          vno, v->nx, v->ny, v->nz, v->area) ;
  print_vertex_data(vno, stdout, 0.0) ;
}

/* Screencoord sx,sy; */
void
select_vertex(short sx,short sy)
{

  /* begin rkt */
  int vno;
  float d;
  char command[NAME_LENGTH];

  /* sets d to the distance of the vertex found */
  find_vertex_at_screen_point(sx, sy, &vno, &d);
  if (vno>=0)
  {
    selection = vno;
    print_vertex_data(selection, stdout, d) ;
  }

  /* select the label at this vertex, if there is one. */
  if (labl_draw_flag && labl_select_flag)
  {
    labl_select_label_by_vno (vno);
  }

  /* if we have functional data... */
  if (func_timecourse)
  {
    func_clear_selection ();

    /* Depending on our avg mode, select the proper voxels and graph
       them. */
    switch (func_graph_avg_mode)
    {
    case FUNC_GRAPH_AVG_MODE_SELECTED:
      func_select_selected_vertex();
      break;
    case FUNC_GRAPH_AVG_MODE_MARKED:
      func_select_marked_vertices();
      break;
    case FUNC_GRAPH_AVG_MODE_LABEL:
      func_select_label();
      break;
    }

    func_graph_timecourse_selection();
  }

  /* select the path at this vertex, if there is one. */
  path_select_path_by_vno (vno);

  /* let the tcl side of things respond. */
  sprintf(command,"SelectVertex %d", vno);
  send_tcl_command (command);

  /* finally, update the labels. */
  if (vno>=0)
  {
    update_labels(LABELSET_CURSOR, vno, d);
  }

  /* Redraw our caption if we need to. */
  if (cptn_draw_flag)
    cptn_draw();

  /* end rkt */
}

/* begin rkt */
void
select_vertex_by_vno (int vno)
{
  VERTEX* v;

  if (vno < 0 || vno >= mris->nvertices)
  {
    return;
  }

  selection = vno;
  print_vertex_data(selection, stdout, 0) ;

  /* if we have functional data... */
  if (func_timecourse)
  {
    v = &(mris->vertices[selection]);

    /* select only this voxel and graph it */
    func_clear_selection ();
    func_select_voxel (selection, v->origx, v->origy, v->origz);
    func_graph_timecourse_selection ();
  }
  
  /* if the vertex index is linked to the timepoint in the config overlay dialog (linkvertexmode is 1 )*/
  if (linkvertexmode ==1)
  {
    sclv_cur_timepoint = vno;
    sclv_set_timepoint_of_field(sclv_current_field, sclv_cur_timepoint, sclv_cur_condition);
    send_tcl_command("UpdateAndRedraw");
  }
  
  /* if the vertex index is linked to the avg (or normalized avg) of ROI in the config overlay dialog (linkvertexmode is 2 or 3)*/
  if (linkvertexmode == 2 || linkvertexmode==3)
  {
    link_timepoint_ROI(vno);
  }

  /* select the label at this vertex, if there is one. */
  if (labl_select_flag)
    labl_select_label_by_vno (vno);

  /* select the path at this vertex, if there is one. */
  path_select_path_by_vno (vno);

  /* finally, update the labels. */
  if (vno>=0)
  {
    update_labels(LABELSET_CURSOR, vno, 0);
  }

  /* Redraw our caption if we need to. */
  if (cptn_draw_flag)
    cptn_draw();
}

void
find_vertex_at_screen_point (short sx, short sy, int* ovno, float* od)
{
  float m[4][4];                /* Surface -> screen matrix */
  int i, j;                     /* Matrix index counters */
  float mf;                     /* Temp matrix value */
  long ox, oy;                  /* Origin of window in global coords */
  long lx, ly;                  /* Size Of/ window */
  float wx, wy;                 /* x/y of click point in screen space */
  float p1[3], p2[3];           /* Segment from screen to back of screen */
  float sImin;                  /* Distance of plane found at closest vert. */
  int fmin;                     /* Closest face. */
  float xmin[3];                /* Closest face intersection point. */
  int fno;                      /* Face counter */
  FACE *f;                      /* Current face */
  VERTEX *v0;                   /* First vertex on face, vert on face plane */
  float plane[3];               /* Point on face plane */
  float n[3];                   /* Normal of face plane */
  float uu[3], ww[3];           /* Vectors of segment and seg->plane */
  float D, N;                   /* Dot values */
  float sI;                     /* Intersection factor */
  float x[3];                   /* Intersection point */
  int vno;                      /* Vertex counter */
  float vs[3];                  /* Vertex to test, in screen space */
  float min[3], max[3];         /* Bounding rect of the face */
  int imin;                     /* vno of closest vertex found */
  float dmin;                   /* Distance of closest vertex found */
  VERTEX *v;                    /* Current vertex when testng int with plane */
  float dx, dy, dz, d;          /* Distance of currect vert to intersection */

  xmin[0] = 0;
  xmin[1] = 0;
  xmin[2] = 0;

  if (Gdiag)
    ddt_clear();

  /* Get the view transformation matrix */
  getmatrix(m);
#ifdef OPENGL
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      mf = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = mf;
    }
#endif

  /* This normalizes our screen point to -137->137 in both directions
     (legacy code, blah) */
  getorigin(&ox,&oy);
  getsize(&lx,&ly);
  wx = (sx-ox-lx/2.0)*2*fov*sf/lx;
  wy = (sy-oy-ly/2.0)*2*fov*sf/ly;

  /* From this, make endpoints of a segment going through the view
     plane. */
  p1[0] = wx;
  p1[1] = wy;
  p1[2] = -10000;
  p2[0] = wx;
  p2[1] = wy;
  p2[2] = 10000;

  /* For each face... */
  sImin = 10000;
  fmin = -1;
  for (fno = 0; fno < mris->nfaces; fno++)
  {
    f = &mris->faces[fno];
    FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
    
    v0 = &mris->vertices[f->v[0]];

    /* We want to get a plane representing this face. Take the first
       vertex in the face as the point on the plane and the plane
       nomral. Multiply by the view transformation matrix. For the
       normal, don't use the translation element of the matrix. */
    plane[0] =   -m[0][0]*v0->x + m[1][0]*v0->z + m[2][0]*v0->y + m[3][0];
    plane[1] =   -m[0][1]*v0->x + m[1][1]*v0->z + m[2][1]*v0->y + m[3][1];
    plane[2] = -(-m[0][2]*v0->x + m[1][2]*v0->z + m[2][2]*v0->y + m[3][2]);

    n[0] =   -m[0][0]*fNorm->nx + m[1][0]*fNorm->nz + m[2][0]*fNorm->ny;
    n[1] =   -m[0][1]*fNorm->nx + m[1][1]*fNorm->nz + m[2][1]*fNorm->ny;
    n[2] = -(-m[0][2]*fNorm->nx + m[1][2]*fNorm->nz + m[2][2]*fNorm->ny);

    /* Make sure the normal's z < 0, so the face is facing us. */
    /* RKT: Removed this test as it caused vertices 'inside' folds
       on the non-inflated surface to be not selected, even if in
       the inflated view. */
    /*       if (n[2] > 0) */
    /*         continue; */

    /* Do an initial distance test in the xy to weed out unnecessary
       intersection tests. */
    dx = plane[0] - wx;
    dy = plane[1] - wy;
    d = sqrt (dx*dx + dy*dy);
#if 0
    if (d > zf*5.0)
      continue;
#endif

    /* Intersect our segment with our plane. */
    /* uu = p2 - p1 */
    uu[0] = p2[0] - p1[0];
    uu[1] = p2[1] - p1[1];
    uu[2] = p2[2] - p1[2];

    /* ww = p1 - plane */
    ww[0] = p1[0] - plane[0];
    ww[1] = p1[1] - plane[1];
    ww[2] = p1[2] - plane[2];

    /* D = n dot uu */
    D = n[0]*uu[0] + n[1]*uu[1] + n[2]*uu[2];

    /* N = - (n dot ww) */
    N = - (n[0]*ww[0] + n[1]*ww[1] + n[2]*ww[2]);

    /* If intersection... */
    if (!(fabs(D) < 0.0001 || fabs(N) < 0.0001))
    {
      /* If the intersection is on the segment... */
      sI = N / D;
      if (sI >= 0.0 && sI <= 1.0)
      {

        /* Get the intersection point */
        /* x = p1 + sI*uu */
        x[0] = p1[0] + sI*uu[0];
        x[1] = p1[1] + sI*uu[1];
        x[2] = p1[2] + sI*uu[2];

        /* Get the bounds of the face. */
        min[0] = min[1] = min[2] = 9999;
        max[0] = max[1] = max[2] = -9999;
        for (vno = 0; vno < VERTICES_PER_FACE; vno++)
        {
          v = &mris->vertices[f->v[vno]];

          if (v->ripflag)
            continue;

          vs[0] =
            -m[0][0]*v->x + m[1][0]*v->z + m[2][0]*v->y + m[3][0];
          vs[1] =
            -m[0][1]*v->x + m[1][1]*v->z + m[2][1]*v->y + m[3][1];

          min[0] = MIN(vs[0],min[0]);
          min[1] = MIN(vs[1],min[1]);

          max[0] = MAX(vs[0],max[0]);
          max[1] = MAX(vs[1],max[1]);
        }

        /* If the intersection is within the bounds, it's a
           hit. NOTE This is a rough estimate, but good enough
           in most cases. */
        if (x[0] >= min[0] && x[0] <= max[0] &&
            x[1] >= min[1] && x[1] <= max[1])
        {
          if (Gdiag && DIAG_VERBOSE_ON)
          {
            ddt_hilite_face (fno, 1);
            fprintf (stderr,"Hit fno %d sI %f\n"
                     "\tbounds %f %f, %f %f\n"
                     "\tint %f %f %f\n",
                     fno, sI,
                     min[0], max[0], min[1], max[1],
                     x[0], x[1], x[2]);
          }

          /* Save the closest face. */
          if (sI < sImin &&
              !f->ripflag)
          {
            sImin = sI;
            fmin = fno;
            xmin[0] = x[0];
            xmin[1] = x[1];
            xmin[2] = x[2];
          }
        }
      }
    }
  }

  /* If we couldn't find a face here, bail. */
  if (-1 ==fmin)
  {
    *ovno = -1;
    *od = -1;

    if (Gdiag & DIAG_VERBOSE_ON)
      fprintf (stderr,"No face found\n");
    return;
  }

  /* Get the closest face and find the vertex closest to the
     intersection point, because we're hitting vertices, not arbitrary
     points. */
  dmin = 1000;
  imin = -1;
  f = &mris->faces[fmin];
  for (vno = 0; vno < VERTICES_PER_FACE; vno++)
  {
    v = &mris->vertices[f->v[vno]];

    vs[0] =   -m[0][0]*v->x + m[1][0]*v->z + m[2][0]*v->y + m[3][0];
    vs[1] =   -m[0][1]*v->x + m[1][1]*v->z + m[2][1]*v->y + m[3][1];
    vs[2] = -(-m[0][2]*v->x + m[1][2]*v->z + m[2][2]*v->y + m[3][2]);

    dx = xmin[0] - vs[0];
    dy = xmin[1] - vs[1];
    dz = xmin[2] - vs[2];
    d = sqrt (dx*dx + dy*dy + dz*dz);

    ddt_hilite_vertex (f->v[vno], 1);

    /* If this is the closest vertex, remember
       it. But not if it's ripped. */
    if (d < dmin &&
        !v->ripflag)
    {
      dmin = d;
      imin = f->v[vno];

      if (Gdiag && DIAG_VERBOSE_ON)
      {
        fprintf (stderr,"\t** Found close vno %d d %f\n"
                 "\t   vs %f %f %f\n"
                 "\t   dx %f dy %f dz %f\n",
                 f->v[vno], d,
                 vs[0], vs[1], vs[2],
                 dx, dy, dz );
        ddt_hilite_vertex (f->v[vno], 2);
      }

    }
  }

  *ovno = imin;
  *od = dmin;

  if (Gdiag && DIAG_VERBOSE_ON)
    fprintf (stderr,"Got vno %d d %f\n", imin, dmin);
}
/* end rkt */

void
left_click(short sx,short sy)
{
#ifdef OPENGL
  sx += w.x;
  sy = 1024 - w.y - sy;
#endif

  if (selection>=0)
    draw_cursor(selection,FALSE);
  select_vertex(sx,sy);
  if (selection>=0)
    mark_vertex(selection,TRUE);
  if (selection>=0)
    draw_cursor(selection,TRUE);
}

void
sample_annotated_image(void)
{
  int i,j;
  int sx,sy;
  int c1,c2,c3;
  float cx,cy,cz,f;
  long ox,oy,lx,ly;
  float m[4][4]; /* Matrix m; */
  VERTEX *v;

  getmatrix(m);
#ifdef OPENGL
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      f = m[i][j];
      m[i][j] = m[j][i];
      m[j][i] = f;
    }
#endif
  getorigin(&ox,&oy);
  getsize(&lx,&ly);
  for (i=0;i<mris->nvertices;i++)
    if (!mris->vertices[i].ripflag)
    {
      v = &mris->vertices[i];
      cx = -m[0][0]*v->x+m[1][0]*v->z+m[2][0]*v->y+m[3][0];
      cy = -m[0][1]*v->x+m[1][1]*v->z+m[2][1]*v->y+m[3][1];
      cz = -(-m[0][2]*v->x+m[1][2]*v->z+m[2][2]*v->y+m[3][2]);
      sx = cx*lx/(2*fov*sf)+lx/2.0;
      sy = cy*ly/(2*fov*sf)+ly/2.0;
      /*
        printf("%d: x=%f y=%f z=%f, sx=%d, sy=%d\n",i,v->x,v->y,v->z,sx,sy);
      */
      if (sy>=0 && sy<frame_ydim && sx>=0 && sx<frame_xdim)
      {
        c1 = framebuff[sy*frame_xdim+sx]&0xff;
        c2 = (framebuff[sy*frame_xdim+sx]>>8)&0xff;
        c3 = (framebuff[sy*frame_xdim+sx]>>16)&0xff;
        if (c1!=0 && c1==c2 && c2==c3)
        {
          v->annotation = framebuff[sy*frame_xdim+sx]&0xff;
          /*
            printf("mris->vertices[%d].annotation=%06x\n",
            i,v->annotation);
          */
        }
      }
    }
}

void
restore_zero_position(void)
{
  if (!openglwindowflag)
  {
    printf("surfer: ### restore_zero_position failed: no gl window open\n");
    PR;
    return;
  }

  loadmatrix(idmat);
  zf = 1.0;
}

void
restore_initial_position(void)
{
  if (!openglwindowflag)
  {
    printf("surfer: ### restore_initial_position failed: no gl window open\n");
    PR;
    return;
  }

  loadmatrix(idmat);
  translate(mris->xctr,-mris->zctr,-mris->yctr);
  zf = 1.0;
}

void
make_lateral_view(char *stem)
{
  if (!openglwindowflag)
  {
    printf("surfer: ### redraw failed: no gl window open\n");
    PR
    return;
  }

  loadmatrix(idmat);
  translate(mris->xctr,-mris->zctr,-mris->yctr);
  if (stem[0]=='r'&&stem[1]=='h')
    rotate_brain(-90.0,'y');
  else
    rotate_brain(90.0,'y');
  zf = 1.0;
}

void
make_medial_view(char *stem)
{
  if (!openglwindowflag)
  {
    printf("surfer: ### redraw failed: no gl window open\n");
    PR
    return;
  }

  loadmatrix(idmat);
  translate(mris->xctr,-mris->zctr,-mris->yctr);
  if (stem[0]=='r'&&stem[1]=='h')
    rotate_brain(90.0,'y');
  else
    rotate_brain(-90.0,'y');
  zf = 1.0;
}

void
make_lateral_view_second(char *stem)
{
  if (!openglwindowflag)
  {
    printf("surfer: ### make_lateral_view_second failed: no gl window open\n");
    PR return;
  }
  if (!secondsurfaceloaded)
  {
    printf("surfer: ### make_lateral_view_second failed: "
           "no second surface read\n");
    PR return;
  }

  loadmatrix(idmat);
  translate(mris2->xctr,-mris2->zctr,-mris2->yctr);
  if (stem[0]=='r'&&stem[1]=='h')
    rotate_brain(-90.0,'y');
  else
    rotate_brain(90.0,'y');
  zf = 1.0;
}

void
write_val_histogram(float min, float max, int nbins)
{
  int i,num,index;
  FILE *fp;
  char fname[200];
  float sum,hist[10000],chist[10000];

  for (i=0;i<nbins;i++)
    hist[i] = chist[i] = 0;
  num = 0;
  for (i=0;i<mris->nvertices;i++)
    if (!mris->vertices[i].ripflag)
    {
      index = floor(nbins*(mris->vertices[i].val-min)/(max-min)+0.5);
      if (index>=0 && index<nbins)
        hist[index]++;
      /*
        if (index<0) index = 0;
        if (index>nbins-1) index = nbins-1;
      */
      num++;
    }
  sum = 0;
  for (i=0;i<nbins;i++)
  {
    hist[i] /= num;
    sum += hist[i];
    chist[i] = sum;
  }
  sprintf(fname,"hist.tmp");
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  for (i=0;i<nbins;i++)
    fprintf(fp,"%f %f\n",min+(i+0.5)*(max-min)/nbins,hist[i]);
  fclose(fp);
  sprintf(fname,"chist.tmp");
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  for (i=0;i<nbins;i++)
    fprintf(fp,"%f %f\n",min+(i+0.5)*(max-min)/nbins,chist[i]);
  fclose(fp);
}

void
print_view_matrix()
{
  int      i, j ;
  float m[4][4];

  getmatrix(m);
  printf("----- view matrix\n");
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
    {
      printf("%f ",m[j][i]);
    }
    printf("\n");
  }
  printf("-----------------\n");
}

void
write_view_matrix(char *dir)
{
  int i,j;
  float m[4][4]; /* Matrix m; */
  char fname[NAME_LENGTH];
  FILE *fp;

  if (!openglwindowflag)
  {
    printf("surfer: ### write_view_matrix failed: no gl window open\n");
    PR return;
  }

  sprintf(fname,"%s/surfer.mat",dir);
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  getmatrix(m);
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
    {
      fprintf(fp,"%13.3e ",m[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  printf("surfer: file %s written\n",fname);
  PR
}

void
read_view_matrix(char *dir)
{
  int i;
  float m[4][4]; /* Matrix m; */
  float a,b,c,d;
  char line[NAME_LENGTH];
  char fname[NAME_LENGTH];
  FILE *fp;

  if (!openglwindowflag)
  {
    printf("surfer: ### read_view_matrix failed: no gl window open\n");
    PR return;
  }

  sprintf(fname,"%s/surfer.mat",dir);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }

  i = 0;
  while (fgets(line,NAME_LENGTH,fp) != NULL)
  {
    if (sscanf(line,"%f %f %f %f",&a,&b,&c,&d) == 4)
    {
      m[i][0] = a;
      m[i][1] = b;
      m[i][2] = c;
      m[i][3] = d;
      i++;
    }
    else
    {
      printf("surfer: ### couldn't parse this line in matrix file:  %s",line);
      printf("surfer: ###   ...read_view_matrix() failed\n");
      PR return;
    }
  }
  loadmatrix(m);
}

void
translate_brain(float x, float y, float z)
{
  int i,j,k;
  float m[4][4], m1[4][4], m2[4][4]; /* Matrix m,m1,m2; */

  if (!openglwindowflag)
  {
    printf("surfer: ### translate_brain failed: no gl window open\n");
    PR return;
  }

  getmatrix(m);
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  m1[3][0] = x;
  m1[3][1] = y;
  m1[3][2] = z;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += m[i][k]*m1[k][j];
    }
  loadmatrix(m2);
}

void
scale_brain(float s)
{
  int i,j,k;
  float m[4][4], m1[4][4], m2[4][4]; /* Matrix m,m1,m2; */

  if (!openglwindowflag)
  {
    printf("surfer: ### scale_brain failed: no gl window open\n");
    PR return;
  }

  zf *= s;

  getmatrix(m);
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  m1[0][0] = s;
  m1[1][1] = s;
  m1[2][2] = s;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += m[i][k]*m1[k][j];
    }
  loadmatrix(m2);
}

void
rotate_brain(float a, char c)
{
  int i,j,k;
  float m[4][4], m1[4][4], m2[4][4]; /* Matrix m,m1,m2; */
  float sa,ca;

  if (!openglwindowflag)
  {
    printf("surfer: ### rotate_brain failed: no gl window open\n");
    PR return;
  }

  getmatrix(m);
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      m1[i][j] = (i==j)?1.0:0.0;
  a = a*M_PI/180;
  sa = sin(a);
  ca = cos(a);
  if (c=='y')
  {
    m1[0][0] = m1[2][2] = ca;
    m1[2][0] = -(m1[0][2] = sa);
  }
  else if (c=='x')
  {
    m1[1][1] = m1[2][2] = ca;
    m1[1][2] = -(m1[2][1] = sa);
  }
  else if (c=='z')
  {
    m1[0][0] = m1[1][1] = ca;
    m1[1][0] = -(m1[0][1] = sa);
  }
  else
  {
    printf("surfer: ### Illegal axis %c\n",c);
    return;
    PR
  }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      m2[i][j] = 0;
      for (k=0;k<4;k++)
        m2[i][j] += m[i][k]*m1[k][j];
    }
  loadmatrix(m2);
}

void read_image_info(char *fpref)
{
  char fname[NAME_LENGTH];
  MRI* mri_header;
  FILE* fTest;

  printf("Reading image info (%s/%s)\n",subjectsdir, pname);

  mri_header = NULL;
  sprintf (fname, "%s.info", fpref);
  fTest = fopen (fname, "r");
  if (NULL != fTest)
  {
    fclose (fTest);
    printf("Reading %s\n",fname);
    mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
  }

  if (NULL == mri_header){
    // Start with orig since it is most likely to be there
    sprintf (fname, "%s/%s/mri/orig.mgz", subjectsdir, pname);
    fTest = fopen (fname, "r");
    if (NULL != fTest)
    {
      fclose (fTest);
      printf("Reading %s\n",fname);
      mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
    }
  }
  if (NULL == mri_header){
    sprintf (fname, "%s/%s/mri/orig.mgh", subjectsdir, pname);
    fTest = fopen (fname, "r");
    if (NULL != fTest)
    {
      fclose (fTest);
      printf("Reading %s\n",fname);
      mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
    }
  }
  if (NULL == mri_header) {
    sprintf (fname, "%s/%s/mri/T1.mgh", subjectsdir, pname);
    fTest = fopen (fname, "r");
    if (NULL != fTest)
    {
      fclose (fTest);
      printf("Reading %s\n",fname);
      mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
    }
  }
  if (NULL == mri_header){
    sprintf (fname, "%s/%s/mri/T1.mgz", subjectsdir, pname);
    fTest = fopen (fname, "r");
    if (NULL != fTest)
    {
      fclose (fTest);
      printf("Reading %s\n",fname);
      mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
    }
  }
  if (NULL == mri_header){
    printf ("ERROR: could not read header info from T1 or orig in %s/%s/mri\n",
            subjectsdir, pname);
    exit(1);
  }
  if (mri_header){
    printf ("surfer: Reading header info from %s\n", fname);
    imnr0 = mri_header->imnr0;
    imnr1 = mri_header->imnr1;
    xnum = mri_header->width;
    ynum = mri_header->height;
    ps = mri_header->xsize;
    st = mri_header->ysize;
    xx0 = mri_header->xstart;
    xx1 = mri_header->xend;
    yy0 = mri_header->ystart;
    yy1 = mri_header->yend;
    zz0 = mri_header->zstart;
    zz1 = mri_header->zend;
    fov = MRIfovCol( mri_header );
    MRIfree (&mri_header);
  }

  /* RKT: Check for fov == 0, which is incorrect. If it is, set it to
     0.256, which is a reasonable default.  */
  /* Krish: For older brains, 256 seems to work than 0.256 
     Also, the fovs of standard brains like bert and fsaverage are 256 */
  if (fabs(fov) < 0.00001)
  {
    printf("surfer: WARNING: fov was ~0, setting to 256\n");
    fov = 256;
  }

  numimg = imnr1-imnr0+1;
}

void
read_talairach(char *fname)    /* marty: ignore abs paths in COR-.info */
{
  lta = LTAreadEx(fname) ;
  if (lta==NULL)
    printf("surfer: Talairach xform file not found (ignored)\n");
  else
  {
    transform_loaded = TRUE;
    if ( lta->type == LINEAR_VOX_TO_VOX )
    {
      lta->xforms[0].m_L = DevolveXFM(pname, lta->xforms[0].m_L, fname);
    }
  }
}

void
read_images(char *fpref)
{
  int i,k;
  FILE *fptr;
  char fname[NAME_LENGTH];

  numimg = imnr1-imnr0+1;
  bufsize = ((unsigned long)xnum)*ynum;
  if (buf==NULL) buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
  for (k=0;k<numimg;k++)
  {
    im[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    for (i=0;i<IMGSIZE;i++)
    {
      im[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
    }
  }
  printf("surfer: allocated image buffer (16 Meg)\n");
  PR for (k=0;k<numimg;k++)
  {
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"r");
    if (fptr==NULL)
    {
      printf("surfer: ### File %s not found\n",fname);
      exit(1);
    }
    fread(buf,sizeof(char),bufsize,fptr);
    buffer_to_image(buf,im[k],xnum,ynum);
    fclose(fptr);
    printf("read image %d (%s)\n",k,fname);
  }
}

void
alloc_curv_images(void)
{
  int i,k;

  for (k=0;k<numimg;k++)
  {
    curvim[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    ocurvim[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    curvimflags[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    for (i=0;i<IMGSIZE;i++)
    {
      curvim[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
      ocurvim[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
      curvimflags[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
    }
  }
  curvim_allocated = TRUE;
  printf("surfer: allocated curv,ocurv,flags images (48 Meg)\n");
  PR
}

void
read_curv_images(char *fpref)/* assumes norm'ed curvim:{CURVIM_NORM_MIN,MAX} */
{
  int i,j,k;
  char fname[NAME_LENGTH];
  FILE *fptr;

  if (!curvim_allocated) alloc_curv_images();

  for (k=0;k<numimg;k++)
  {
    file_name(fpref,fname,k+imnr0,"%03d");
    fptr = fopen(fname,"r");
    if (fptr==NULL)
    {
      printf("surfer: ### File %s not found\n",fname);
      PR return;
    }
    if (buf==NULL)
      buf = (unsigned char *)lcalloc(IMGSIZE*IMGSIZE,sizeof(char));
    fread(buf,sizeof(char),IMGSIZE*IMGSIZE,fptr);
    buffer_to_image(buf,curvim[k],xnum,ynum);
    fclose(fptr);
    printf("surfer: read curv image %d (%s)\n",k,fname);
  }
  /* next two used by:  ellipsoid_shrink, second_surface_curv_to_curvim */
  /*   smooth_curvim_sparse, read_curvim_at_vertex, curv_shrink_to_fill */
  mris2->min_curv = CURVIM_NORM_MIN;
  mris2->max_curv = CURVIM_NORM_MAX;

  /* hack: byte0 => UNDEFINED; should save flags to read w/images!! */
  /* no FIXEDVAL allows smooth; should write smooth_curvim() */
  for (k=0;k<numimg;k++)
    for (i=0;i<IMGSIZE;i++)
      for (j=0;j<IMGSIZE;j++)
      {
        if (curvim[k][i][j]!=0)
          curvimflags[k][i][j] |= CURVIM_DEFINED;
      }
  printf("surfer: non-zero bytes DEFINED (not FIXEDVAL) in curv image\n");

  curvimflag = TRUE;
  curvimloaded = TRUE;
}

void
curv_to_curvim(void)
{
  VERTEX *v;
  int imnr,i,j,k,vdef,pdef;
  float x,y,z;

  if (!curvim_allocated)
  {
    alloc_curv_images();
  }
  else
  {
    for (k=0;k<numimg;k++)
      for (i=0;i<IMGSIZE;i++)
        for (j=0;j<IMGSIZE;j++)
          curvim[k][i][j] = ocurvim[k][i][j] = curvimflags[k][i][j] = 0;
    printf("surfer: cleared curv,ocurv,flags images");
    PR
  }

  vdef = pdef = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;
    imnr = (int)((y-yy0)/st+0.5);
    i = (int)((zz1-z)/ps+0.5);
    j = (int)((xx1-x)/ps+0.5);
    imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
    i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
    j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
    curvim[imnr][i][j] = floattobyte(v->curv,CURVIM_NORM_MIN,CURVIM_NORM_MAX);
    if (!(curvimflags[imnr][i][j] & CURVIM_DEFINED))
      pdef++;
    curvimflags[imnr][i][j] = CURVIM_DEFINED | CURVIM_FIXEDVAL;
    vdef++;
  }
  printf("surfer: curv image made--%d pix set\n",vdef);
  PR printf("surfer:                  %d unique pix\n",pdef);
  PR curvimflag = TRUE;
  curvimloaded = TRUE;
}

void
second_surface_curv_to_curvim(void)
{
  VERTEX *v;
  int imnr,i,j,k,vdef,pdef;
  float x,y,z;

  if (!secondsurfaceloaded)
  {
    printf("surfer: ### second surface not loaded!\n");
    PR return;
  }
  if (!secondcurvloaded)
  {
    printf("surfer: ### second curv not loaded!\n");
    PR return;
  }

  if (!curvim_allocated)
  {
    alloc_curv_images();
  }
  else
  {
    for (k=0;k<numimg;k++)
      for (i=0;i<IMGSIZE;i++)
        for (j=0;j<IMGSIZE;j++)
          curvim[k][i][j] = ocurvim[k][i][j] = curvimflags[k][i][j] = 0;
    printf("surfer: cleared curv,ocurv,flags images");
    PR
  }

  vdef = pdef = 0;
  for (k=0;k<mris2->nvertices;k++)
  {
    v = &mris2->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;
    imnr = (int)((y-yy0)/st+0.5);
    i = (int)((zz1-z)/ps+0.5);
    j = (int)((xx1-x)/ps+0.5);
    imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
    i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
    j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
    curvim[imnr][i][j] = floattobyte(v->curv,mris2->min_curv,mris2->max_curv);
    if (!(curvimflags[imnr][i][j] & CURVIM_DEFINED))
      pdef++;
    curvimflags[imnr][i][j] = CURVIM_DEFINED | CURVIM_FIXEDVAL;
    vdef++;
  }
  printf("surfer: curv image made from second surface--%d pix set\n",vdef);
  PR printf("surfer:                                      %d unique pix\n",
            pdef);
  PR curvimflag = TRUE;
  curvimloaded = TRUE;
}

void
swap_curv(void)
{
  VERTEX *v;
  int k;
  float tmp;

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    tmp = v->curv;
    v->curv = v->curvbak;
    v->curvbak = tmp;
  }
}

/* begin rkt */
void
swap_vertex_fields(int typea, int typeb)
{
  VERTEX *v;
  int k;
  float *a;
  float tmp;

  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);

  /* for every vertex, swap the values specified in the parameters */
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    /* point us in the direction of the first field */
    a = NULL;
    switch (typea)
    {
    case FIELD_CURV:
      a = &(v->curv);
      break;
    case FIELD_CURVBAK:
      a = &(v->curvbak);
      break;
    case FIELD_VAL:
      a = &(v->val);
      break;
    case FIELD_IMAG_VAL:
      a = &(v->imag_val);
      break;
    case FIELD_VAL2:
      a = &(v->val2);
      break;
    case FIELD_VALBAK:
      a = &(v->valbak);
      break;
    case FIELD_VAL2BAK:
      a = &(v->val2bak);
      break;
    case FIELD_STAT:
      a = &(v->stat);
      break;
    default:
      printf("### surfer: can't switch surface fields, invalid type %d\n",
             typea);
    }
    /* save the value */
    tmp = *a;
    /* get the second field. set the value of the first field to the second
       field and set the second to the saved value. */
    switch (typeb)
    {
    case FIELD_CURV:
      *a = v->curv;
      v->curv = tmp;
      break;
    case FIELD_CURVBAK:
      *a = v->curvbak;
      v->curvbak = tmp;
      break;
    case FIELD_VAL:
      *a = v->val;
      v->val = tmp;
      break;
    case FIELD_IMAG_VAL:
      *a = v->imag_val;
      v->imag_val = tmp;
      break;
    case FIELD_VAL2:
      *a = v->val2;
      v->val2 = tmp;
      break;
    case FIELD_VALBAK:
      *a = v->valbak;
      v->valbak = tmp;
      break;
    case FIELD_VAL2BAK:
      *a = v->val2bak;
      v->val2bak = tmp;
      break;
    case FIELD_STAT:
      *a = v->stat;
      v->stat = tmp;
      break;
    default:
      printf("### surfer: can't switch surface fields, invalid type %d\n",
             typeb);
    }
  }

  /* TODO: update_labels for last selected vertex */
}
/* end rkt */

void
curvim_to_surface(void)   /* assumes norm'ed curvim:{CURVIM_NORM_MIN,MAX} */
{
  VERTEX *v;
  int imnr,i,j,k;
  float x,y,z;

  if (!curvimloaded)
  {
    printf("surfer: ### curvim not loaded!\n");
    PR return;
  }

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;
    imnr = (int)((y-yy0)/st+0.5);
    i = (int)((zz1-z)/ps+0.5);
    j = (int)((xx1-x)/ps+0.5);
    imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
    i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
    j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;

    /* hack: byte 0 => UNDEFINED; should save flags w/images!! */
    /* no FIXEDVAL allows smooth; should write smooth_curvim() */
    if (curvim[imnr][i][j]==0)
    {
      v->curv = 0.0;
    }
    else
    {
      v->curv= bytetofloat(curvim[imnr][i][j],CURVIM_NORM_MIN,CURVIM_NORM_MAX);
      curvimflags[imnr][i][j] |= CURVIM_DEFINED; /* Why? AMD */
    }
  }
}

/* assumes norm'ed curvim:{CURVIM_NORM_MIN,MAX} */
void
curvim_to_second_surface(void)
{
  VERTEX *v;
  int imnr,i,j,k;
  float x,y,z;

  if (!curvimloaded)
  {
    printf("surfer: ### curvim not loaded!\n");
    PR return;
  }
  if (!secondsurfaceloaded)
  {
    printf("surfer: ### second surface not loaded!\n");
    PR return;
  }

  for (k=0;k<mris2->nvertices;k++)
  {
    v = &mris2->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;
    imnr = (int)((y-yy0)/st+0.5);
    i = (int)((zz1-z)/ps+0.5);
    j = (int)((xx1-x)/ps+0.5);
    imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
    i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
    j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;

    /* hack: byte 0 => UNDEFINED; should save flags w/images!! */
    /* no FIXEDVAL allows smooth; should write smooth_curvim() */
    if (curvim[imnr][i][j]==0)
    {
      v->curv = 0.0;
    }
    else
    {
      v->curv = bytetofloat(curvim[imnr][i][j],
                            CURVIM_NORM_MIN,CURVIM_NORM_MAX);
      curvimflags[imnr][i][j] |= CURVIM_DEFINED;
    }
  }
  secondcurvloaded = TRUE;
}

void
smooth_curvim(int window)
{
  int i,j,k,di,dj,dk;
  double avgcurv,numcurv;

  printf("smooth_curvim(%d)\n",window);

  for (k=0;k<numimg;k++)
    for (i=0;i<IMGSIZE;i++)
      for (j=0;j<IMGSIZE;j++)
        ocurvim[k][i][j] = curvim[k][i][j];

  for (i=0;i<IMGSIZE;i++)
  {
    printf(".");
    fflush(stdout);
    for (j=0;j<IMGSIZE;j++)
    {
      for (k=0;k<numimg;k++)
        if (curvimflags[k][i][j] & CURVIM_DEFINED)
        {
          avgcurv = numcurv = 0;
          for (dk= -window;dk<=window;dk++)
            if ((k+dk>=0)&&(k+dk<numimg))
              if (curvimflags[k+dk][i][j] & CURVIM_DEFINED)
              {
                avgcurv +=
                  bytetofloat(ocurvim[k+dk][i][j],
                              CURVIM_NORM_MIN,CURVIM_NORM_MAX);
                numcurv++;
              }
          if (numcurv>0)
            avgcurv /= numcurv;
          curvim[k][i][j] = floattobyte(avgcurv,
                                        CURVIM_NORM_MIN,CURVIM_NORM_MAX);
        }
    }
  }
  printf("\n");

  for (k=0;k<numimg;k++)
    for (i=0;i<IMGSIZE;i++)
      for (j=0;j<IMGSIZE;j++)
        ocurvim[k][i][j] = curvim[k][i][j];

  for (k=0;k<numimg;k++)
  {
    printf(".");
    fflush(stdout);
    for (j=0;j<IMGSIZE;j++)
    {
      for (i=0;i<IMGSIZE;i++)
        if (curvimflags[k][i][j] & CURVIM_DEFINED)
        {
          avgcurv = numcurv = 0;
          for (di= -window;di<=window;di++)
            if ((i+di>=0)&&(i+di<IMGSIZE))
              if (curvimflags[k][i+di][j] & CURVIM_DEFINED)
              {
                avgcurv +=
                  bytetofloat(ocurvim[k][i+di][j],
                              CURVIM_NORM_MIN,CURVIM_NORM_MAX);
                numcurv++;
              }
          if (numcurv>0)
            avgcurv /= numcurv;
          curvim[k][i][j] =
            floattobyte(avgcurv,CURVIM_NORM_MIN,CURVIM_NORM_MAX);
        }
    }
  }
  printf("\n");

  for (k=0;k<numimg;k++)
    for (i=0;i<IMGSIZE;i++)
      for (j=0;j<IMGSIZE;j++)
        ocurvim[k][i][j] = curvim[k][i][j];

  for (k=0;k<numimg;k++)
  {
    printf(".");
    fflush(stdout);
    for (i=0;i<IMGSIZE;i++)
    {
      for (j=0;j<IMGSIZE;j++)
        if (curvimflags[k][i][j] & CURVIM_DEFINED)
        {
          avgcurv = numcurv = 0;
          for (dj= -window;dj<=window;dj++)
            if ((j+dj>=0)&&(j+dj<IMGSIZE))
              if (curvimflags[k][i][j+dj] & CURVIM_DEFINED)
              {
                avgcurv +=
                  bytetofloat(ocurvim[k][i][j+dj],
                              CURVIM_NORM_MIN,CURVIM_NORM_MAX);
                numcurv++;
              }
          if (numcurv>0)
            avgcurv /= numcurv;
          curvim[k][i][j] =
            floattobyte(avgcurv,CURVIM_NORM_MIN,CURVIM_NORM_MAX);
        }
    }
  }
  printf("\n");
}

#if 0
smooth_curvim(window)
int window;
{
  int i,j,k,di,dj,dk;
  double curv,avgcurv,numcurv;

  printf("smooth_curvim(%d)\n",window);
  for (k=0;k<numimg;k++)
    for (i=0;i<IMGSIZE;i++)
      for (j=0;j<IMGSIZE;j++)
        ocurvim[k][i][j] = curvim[k][i][j];

  for (k=0;k<numimg;k++)
  {
    printf(".");
    fflush(stdout);
    for (i=0;i<IMGSIZE;i++)
      for (j=0;j<IMGSIZE;j++)
        if (curvimflags[k][i][j] & CURVIM_DEFINED)
        {
          avgcurv = numcurv = 0;
          for (dk= -window;dk<=window;dk++)
            for (di= -window;di<=window;di++)
              for (dj= -window;dj<=window;dj++)
                if ((k+dk>=0)&&(k+dk<numimg)&&
                    (i+di>=0)&&(i+di<IMGSIZE)&&
                    (j+dj>=0)&&(j+dj<IMGSIZE))
                  if (curvimflags[k+dk][i+di][j+dj] & CURVIM_DEFINED)
                  {
                    avgcurv +=
                      bytetofloat(ocurvim[k+dk][i+di][j+dj],
                                  mris2->min_curv,mris2->max_curv);
                    numcurv++;
                  }
          avgcurv /= numcurv;
          curvim[k][i][j] =
            floattobyte(avgcurv,mris2->min_curv,mris2->max_curv);
        }
  }
  printf("\n");
}
#endif

/* assumes norm'd curvim */
void
add_subject_to_average_curvim(char *name, char *morphsubdir)
{
  int i,j,k;
  float curv,avgcurv;
  char fname[NAME_LENGTH],fpref[NAME_LENGTH];
  FILE* test;

  /* check if images there first */
  /*sprintf(fpref,"%s/%s/%s.%s.%s/COR-",
    subjectsdir,name,CURVDIR_STEM,stem,ext);*/
  sprintf(fname,"%s/%s/morph/%s/COR-001",subjectsdir,name,morphsubdir);
  if ((test = fopen(fname,"r"))==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  fclose (test);

  /* avgcurv->old */
  if (!curvim_allocated)
    alloc_curv_images();
  else
  {
    for (k=0;k<numimg;k++)
      for (i=0;i<IMGSIZE;i++)
        for (j=0;j<IMGSIZE;j++)
        {
          ocurvim[k][i][j] = curvim[k][i][j];
        }
  }
  /*sprintf(fpref,"%s/%s/%s.%s.%s/COR-",
    subjectsdir,name,CURVDIR_STEM,stem,ext);*/
  sprintf(fpref,"%s/%s/morph/%s/COR-",subjectsdir,name,morphsubdir);
  read_curv_images(fpref);

  if (curvim_averaged)
  {
    for (k=0;k<numimg;k++)
      for (i=0;i<IMGSIZE;i++)
        for (j=0;j<IMGSIZE;j++)
        {
          avgcurv =
            bytetofloat(ocurvim[k][i][j],CURVIM_NORM_MIN,CURVIM_NORM_MAX);
          curv = bytetofloat(curvim[k][i][j],CURVIM_NORM_MIN,CURVIM_NORM_MAX);
          avgcurv = (avgcurv*curvim_averaged + curv)/(curvim_averaged+1);
          curvim[k][i][j] =
            floattobyte(avgcurv,CURVIM_NORM_MIN,CURVIM_NORM_MAX);
        }
  }
  else
  {
    curvimflag = TRUE;
    curvimloaded = TRUE;
  }
  curvim_averaged++;
  printf("surfer: %d subjects in average curvim\n",curvim_averaged);
}

void
smooth_curvim_sparse(int niter)
{
  int iter,i,j,k,n;
  int i2,j2,k2;
  int ndef;
  float sum;

  printf("surfer: smooth_curvim_sparse:\n");
  PR for (iter=0;iter<niter;iter++)
  {
    ndef = 0;
    for (k=0;k<numimg;k++)
      for (i=0;i<IMGSIZE;i++)
        for (j=0;j<IMGSIZE;j++)
        {
          ocurvim[k][i][j] = curvim[k][i][j];
          if (curvimflags[k][i][j] & CURVIM_DEFINED)
          {
            curvimflags[k][i][j] |= CURVIM_DEFINED_OLD;
            ndef++;
          }
        }
    printf("surfer:  iter = %d  defined curv pix = %d\n",iter,ndef);
    PR for (k=1;k<numimg-1;k++)
    {
      for (i=1;i<IMGSIZE-1;i++)
      {
        for (j=1;j<IMGSIZE-1;j++)
        {
          if (!(curvimflags[k][i][j] & CURVIM_FIXEDVAL))
          {
            sum = 0;
            n = 0;
#if 0
            if (curvimflags[k][i][j] & CURVIM_DEFINED_OLD)
            {  /* center + 6 */
              sum +=
                bytetofloat(ocurvim[k2][i][j],mris2->min_curv,mris2->max_curv);
              n++;
            }
            for (k2=k-1;k2<k+2;k2+=2)
            {  /* 6 */
              if (curvimflags[k2][i][j] & CURVIM_DEFINED_OLD)
              {
                sum +=
                  bytetofloat(ocurvim[k2][i][j],
                              mris2->min_curv,mris2->max_curv);
                n++;
              }
            }
            for (i2=i-1;i2<i+2;i2+=2)
            {
              if (curvimflags[k][i2][j] & CURVIM_DEFINED_OLD)
              {
                sum += bytetofloat(ocurvim[k][i2][j],
                                   mris2->min_curv,mris2->max_curv);
                n++;
              }
            }
            for (j2=j-1;j2<j+2;j2+=2)
            {
              if (curvimflags[k][i][j2] & CURVIM_DEFINED_OLD)
              {
                sum += bytetofloat(ocurvim[k][i][j2],
                                   mris2->min_curv,mris2->max_curv);
                n++;
              }
            }
#endif
            for (k2=k-1;k2<k+2;k2++)  /* 27 */
              for (i2=i-1;i2<i+2;i2++)
                for (j2=j-1;j2<j+2;j2++)
                {
                  if (curvimflags[k2][i2][j2] & CURVIM_DEFINED_OLD)
                  {
                    sum += bytetofloat(ocurvim[k2][i2][j2],
                                       mris2->min_curv,mris2->max_curv);
                    n++;
                  }
                }
            if (n>0)
            {
              curvim[k][i][j] =
                floattobyte(sum/n,mris2->min_curv,mris2->max_curv);
              curvimflags[k][i][j] |= CURVIM_DEFINED;
            }
          }
        }
      }
    }
  }
}

unsigned char
floattobyte(float f, float min, float max)
{
  f = (f>max)?max:(f<min)?min:f;  /* needed? */
  return (unsigned char)((f-min)*(255/(max-min)));
}

float
bytetofloat(unsigned char c, float min, float max)
{
  /*if ((int)c>255) {
    printf("surfer: bad input byte %d for current 255\n",(int)c);
    return HUGE_VAL;
    }*/
  return (float)c/(255/(max-min)) + min;
}

/* tksurfer.c: write_{curv,fill}_images */
void
write_images(unsigned char ***mat,char *fpref)
{
  int k;
  FILE *fptr;
  char fname[NAME_LENGTH];

  bufsize = ((unsigned long)xnum)*ynum;
  if (buf==NULL) buf = (unsigned char *)lcalloc(bufsize,sizeof(char));
  for (k=0;k<numimg;k++)
  {
    file_name(fpref,fname,k+1,"%03d");
    fptr = fopen(fname,"w");
    if (fptr==NULL)
    {
      printf("surfer: ### can't create file %s\n",fname);
      PR return;
    }
    image_to_buffer(mat[k],buf,xnum,ynum);
    fwrite(buf,sizeof(char),bufsize,fptr);
    fclose(fptr);
    printf("surfer: file %s written\n",fname);
    PR
  }
}

void
save_surf(void)
{
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
}
void
restore_surf(void)
{
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeNormals(mris) ;
}
void
read_positions(char *name)
{
  if (MRISreadVertexPositions(mris, name) == NO_ERROR)
  {
    mris->status = MRIS_SURFACE ;
    MRIScomputeMetricProperties(mris) ;
  }
  surface_compiled = 0 ;
}
/*------------------------------------------------------- */
int read_binary_surface(char *fname)
{
  if (mris)  MRISfree(&mris) ;
  mris = MRISread(ifname) ;
  if (!mris)  return(Gerror) ;

  if (LeftRightRev)
  {
    printf("Applying Left-Right reversal\n");
    MRISreverse(mris, REVERSE_X, 1) ;
  }

  marked = (int *)calloc(mris->nvertices, sizeof(int)) ;
  if (!marked) ErrorExit(ERROR_BADPARM, "%s: could not allocate %d vertex list array",
                           Progname, mris->nvertices) ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  flag2d = FALSE ;
  if (!mris)
    return(ERROR_NOFILE) ;
  MRISclearCurvature(mris) ;
  printf("surfer: vertices=%d, faces=%d\n",mris->nvertices,mris->nfaces);
  numvertices = mris->nvertices;
  surface_compiled = 0 ;
  vertex_array_dirty = 1; /* rkt - added this */

  /* begin rkt */
  /* since we just changed the surface and the number of verts, a lot
     of stuff that has initialized itself based on the number of verts
     is no longer valid. so we call the init functions to re-init a
     lot of this stuff. */
  init_vertex_arrays(mris);
  vset_initialize();
  func_initialize();
  sclv_initialize();
  conv_initialize();
  labl_initialize();
  path_initialize();
  ddt_initialize();
  /* end rkt */

  return(NO_ERROR) ;
}


void
read_second_binary_surface(char *fname)   /* inlined hack */
{
  mris2 = MRISread(fname) ;
  PR
}

void
read_second_binary_curvature(char *fname)
{
  MRISreadCurvatureFile(mris2, fname) ;
  secondcurvloaded = TRUE;
  surface_compiled = 0 ;
  printf("surfer: second curvature read: min=%f max=%f\n",
         mris2->min_curv,mris2->max_curv);
}


void
normalize_second_binary_curvature(void)
{
  int k;
  float curv,min,max;
  float sum,avg,sum_sq,sd,n;

  min = max = 0.0f ;
  if (!secondcurvloaded)
  {
    printf("surfer: ### second curv not loaded!\n");
    PR return;
  }

  sum = 0;
  for (k=0;k<mris2->nvertices;k++)
    sum += mris2->vertices[k].curv;
  avg = sum/mris2->nvertices;

  n = (float)mris2->nvertices;
  sum = sum_sq = 0.0;
  for (k=0;k<mris2->nvertices;k++)
  {
    mris2->vertices[k].curv -= avg;
    curv = mris2->vertices[k].curv;
    sum += curv;
    sum_sq += curv*curv;
  }
  sd = sqrt((n*sum_sq - sum*sum)/(n*(n-1.0)));

  for (k=0;k<mris2->nvertices;k++)
  {
    curv = (mris2->vertices[k].curv)/sd;
    if (k==0) min=max=curv;
    if (curv>max) max=curv;
    if (curv<min) min=curv;
    if (curv<CURVIM_NORM_MIN) curv = CURVIM_NORM_MIN;
    if (curv>CURVIM_NORM_MAX) curv = CURVIM_NORM_MAX;
    mris2->vertices[k].curv = curv;
  }
  mris2->min_curv = CURVIM_NORM_MIN;
  mris2->max_curv = CURVIM_NORM_MAX;
  printf("surfer: second curvature normalized: avg=%f sd=%f\n",avg,sd);
  printf("surfer: min=%f max=%f trunc (%f,%f)\n",
         min,max,mris2->min_curv,mris2->max_curv);
  PR
}

static void
read_imag_vals(char *fname)
{
  MRISreadImagValues(mris, fname) ;
}

static char *last_subject_name = NULL ;

static void
resend_to_subject(void)
{
  if (!last_subject_name)
  {
    printf("must send_to_subject to specify subject name first.\n") ;
    return ;
  }

  send_to_subject(last_subject_name) ;
}

void
send_to_subject(char *subject_name)
{
  char canon_name[STRLEN], orig_name[STRLEN] ;

  sprintf(canon_name, "%s.%s", stem, sphere_reg) ;
  sprintf(orig_name, "%s.white", stem) ;
#if 0
  send_spherical_point(subject_name, canon_name, orig_name) ;
#else
  send_spherical_point(subject_name, canon_name,
                       FileNameOnly(orfname, orig_name)) ;
#endif
}
void
send_spherical_point(char *subject_name, char *canon_name, char *orig_name)
{
  float           x, y, z, dx, dy, dz, dist, min_dist ;
  SMALL_VERTEX    *sv ;
  VERTEX          *v ;
  char            fname[STRLEN] ;
  SMALL_SURFACE   *mriss ;
  int             vno, min_vno ;
  FILE            *fp ;

  if (selection < 0)
  {
    printf("must select a vertex.\n") ;
    return ;
  }
  if (canonsurfloaded == FALSE)
  {
    if (DIAG_VERBOSE_ON)
      printf("reading canonical vertex positions from %s...\n", canon_name) ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    if (MRISreadVertexPositions(mris, canon_name) != NO_ERROR)
    {
      canonsurffailed = TRUE ;
      return ;
    }
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    canonsurfloaded = TRUE ;
  }

  sprintf(fname, "%s/%s/surf/%s", subjectsdir, subject_name, canon_name) ;
  if (DIAG_VERBOSE_ON)
    printf("reading spherical coordinates for subject %s from %s...\n",
           subject_name, fname) ;

  mriss = MRISSread(fname) ;
  if (!mriss)
  {
    fprintf(stderr, "### could not open surface file %s\n", fname) ;
    return ;
  }

  v = &mris->vertices[selection] ;
  x = v->cx ;
  y = v->cy ;
  z = v->cz ;
  min_dist = 1000000.0f ;
  min_vno = -1 ;
  for (vno = 0 ; vno < mriss->nvertices ; vno++)
  {
    sv = &mriss->vertices[vno] ;
    dx = sv->x - x ;
    dy = sv->y - y ;
    dz = sv->z - z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (dist < min_dist)
    {
      min_dist = dist ;
      min_vno  = vno ;
    }
  }
  MRISSfree(&mriss) ;

  printf("closest vertex %d is %2.1f mm distant.\n", min_vno, min_dist) ;

  sprintf(fname, "%s/%s/surf/%s", subjectsdir, subject_name, orig_name) ;
  if (DIAG_VERBOSE_ON)
    printf("reading original coordinates for subject %s from %s...\n",
           subject_name, fname) ;

  mriss = MRISSread(fname) ;
  if (!mriss)
  {
    fprintf(stderr, "### could not open surface file %s\n", fname) ;
    return ;
  }

  sv = &mriss->vertices[min_vno] ;

#if 0
  copy_edit_dat_file_name (fname, sizeof(fname));
#else
  sprintf(fname, "%s/%s/tmp/edit.dat", subjectsdir, subject_name) ;
#endif
  if (DIAG_VERBOSE_ON)
    printf("writing coordinates to file %s\n", fname) ;
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  if (DIAG_VERBOSE_ON)
    printf("vertex %d coordinates:\n", min_vno);
#if 0
  if (transform_loaded)
    printf("TALAIRACH (%2.1f %2.1f %2.1f)\n",x_tal,y_tal,z_tal);
#endif
  if (DIAG_VERBOSE_ON)
    printf("ORIGINAL  (%2.1f %2.1f %2.1f)\n",x,y,z);

  /* Write the point to the file. */
  fprintf(fp,"%f %f %f\n",sv->x,sv->y,sv->z);

#if 0
  fprintf(fp,"%f %f %f\n",x_tal,y_tal,z_tal);
#endif
  fclose(fp);
}
void
send_to_other_hemi(void) 
{
  send_contralateral_point(sphere_reg_contra, white_suffix) ;
}

void
send_contralateral_point(char *canon_name, char *orig_name)
{
  float           x, y, z, dx, dy, dz, dist, min_dist ;
  SMALL_VERTEX    *sv ;
  VERTEX          *v ;
  char            fname[STRLEN], *ohemi ;
  SMALL_SURFACE   *mriss ;
  int             vno, min_vno ;
  FILE            *fp ;

  if (strcmp(stem, "lh") == 0)
    ohemi = "rh" ;
  else
    ohemi = "lh" ;

  if (selection < 0)
  {
    printf("must select a vertex.\n") ;
    return ;
  }

  // read in left/right canonical vertices and save current and current canonical
  sprintf(fname, "%s/%s/surf/%s.%s", subjectsdir, pname, stem, canon_name) ;
  printf("reading canonical vertex positions from %s...\n", fname) ;
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISsaveVertexPositions(mris, TMP2_VERTICES) ;
  if (MRISreadVertexPositions(mris, fname) != NO_ERROR)
    return ;
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;

  sprintf(fname, "%s/%s/surf/%s.%s", subjectsdir, pname, ohemi, canon_name) ;
  printf("reading contralateral canonical vertex positions from %s...\n", fname) ;
  mriss = MRISSread(fname) ;
  if (!mriss)
  {
    fprintf(stderr, "### could not open surface file %s\n", fname) ;
    return ;
  }

  v = &mris->vertices[selection] ;
  x = v->cx ; y = v->cy ; z = v->cz ;
  min_dist = 1000000.0f ;
  min_vno = -1 ;
  for (vno = 0 ; vno < mriss->nvertices ; vno++)
  {
    sv = &mriss->vertices[vno] ;
    dx = sv->x - x ; dy = sv->y - y ; dz = sv->z - z ;
    dist = sqrt(dx*dx + dy*dy + dz*dz) ;
    if (dist < min_dist)
    {
      min_dist = dist ;
      min_vno  = vno ;
    }
  }
  MRISSfree(&mriss) ;

  printf("closest vertex %d is %2.1f mm distant.\n", min_vno, min_dist) ;

  sprintf(fname, "%s/%s/surf/%s.%s", subjectsdir, pname, ohemi, orig_name) ;
  printf("reading original coordinates from %s...\n", fname) ;

  mriss = MRISSread(fname) ;
  if (!mriss)
  {
    fprintf(stderr, "### could not open surface file %s\n", fname) ;
    return ;
  }

  sv = &mriss->vertices[min_vno] ;

#if 0
  copy_edit_dat_file_name (fname, sizeof(fname));
#else
  sprintf(fname, "%s/%s/tmp/edit.dat", subjectsdir, pname) ;
#endif
  if (DIAG_VERBOSE_ON)
    printf("writing coordinates to file %s\n", fname) ;
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  if (DIAG_VERBOSE_ON)
    printf("vertex %d coordinates:\n", min_vno);
#if 0
  if (transform_loaded)
    printf("TALAIRACH (%2.1f %2.1f %2.1f)\n",x_tal,y_tal,z_tal);
#endif
  if (DIAG_VERBOSE_ON)
    printf("ORIGINAL  (%2.1f %2.1f %2.1f)\n",x,y,z);

  /* Write the point to the file. */
  fprintf(fp,"%f %f %f\n",sv->x,sv->y,sv->z);

#if 0
  fprintf(fp,"%f %f %f\n",x_tal,y_tal,z_tal);
#endif
  fclose(fp);
  MRISrestoreVertexPositions(mris, TMP2_VERTICES) ;  // restore old canonical
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;  // restore old current
}
/*--------------------------------------------------------- */
int read_white_vertex_coordinates(void)
{
  int n;
  fprintf(stderr, "reading white matter vertex locations...\n") ;
#if 1
  if (MRISreadWhiteCoordinates(mris, white_suffix) == NO_ERROR)
  {
    if (LeftRightRev)
    {
      printf("Applying Left-Right reversal\n");
      for (n=0; n < mris->nvertices; n++) mris->vertices[n].whitex *= -1.0;
      // Don't reverse faces
    }
    white_surf_loaded = TRUE ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris, WHITE_VERTICES) ;
    vset_save_surface_vertices( VSET_WHITE );
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    enable_menu_set( MENUSET_VSET_WHITE_LOADED, 1 );
  }
  else
    return(Gerror) ;
#else
  if (MRISreadOriginalProperties(mris, white_suffix) != NO_ERROR)
    return(Gerror) ;
  surface_compiled = 0 ;
  white_surf_loaded = 1 ;
#endif
  return(NO_ERROR) ;
}

int read_inflated_vertex_coordinates(void)
{
  int n;
  fprintf(stderr, "reading inflated vertex locations...\n") ;
  if (MRISreadInflatedCoordinates(mris, "inflated") != NO_ERROR)
    return(Gerror) ;
  if (LeftRightRev)
  {
    printf("Applying Left-Right reversal\n");
    for (n=0; n < mris->nvertices; n++) mris->vertices[n].infx *= -1.0;
    // Don't reverse faces
  }
  surface_compiled = 0 ;
  inflated_surf_loaded = 1 ;
  return(NO_ERROR) ;
}

int read_pial_vertex_coordinates(void)
{
  fprintf(stderr, "reading pial surface vertex locations...\n") ;
  /*  MRISsaveVertexPositions(mris, TMP_VERTICES) ;*/
  if (MRISreadVertexPositions(mris, "pial") != NO_ERROR)
  {
    fprintf(stderr, "could not read canonical surface from 'pial'\n") ;
    return(Gerror) ;
  }
  MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
  /*  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;*/
  pial_surf_loaded = 1 ;
  surface_compiled = 0 ;
  return(NO_ERROR) ;
}

void show_surf(char *surf_name)
{
  if (!stricmp(surf_name, "white"))
  {
    if (!white_surf_loaded)
      read_white_vertex_coordinates() ;
    MRISrestoreVertexPositions(mris, ORIG_VERTICES) ;
  }
  else if (!stricmp(surf_name, "pial") || !stricmp(surf_name, "folded"))
  {
    if (!pial_surf_loaded)
      read_pial_vertex_coordinates() ;
    MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  }
  else if (!stricmp(surf_name, "inflated"))
  {
    if (!inflated_surf_loaded)
      read_inflated_vertex_coordinates() ;
    MRISrestoreVertexPositions(mris, INFLATED_VERTICES) ;
  }
  else
  {
    fprintf(stderr, "unknown surface %s\n", surf_name) ;
    return ;
  }
  MRIScomputeNormals(mris) ;
  surface_compiled = 0 ;
  vertex_array_dirty = 1 ;
  redraw() ;
}

int read_orig_vertex_coordinates(char *fname)
{
  int n;
  if (MRISreadOriginalProperties(mris, fname) == NO_ERROR)
  {
    if (LeftRightRev)
    {
      printf("Applying Left-Right reversal\n");
      for (n=0; n < mris->nvertices; n++) mris->vertices[n].origx *= -1.0;
      // Don't reverse faces
    }
    origsurfloaded = TRUE ;
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    vset_save_surface_vertices( VSET_ORIGINAL );
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    enable_menu_set( MENUSET_VSET_ORIGINAL_LOADED, 1 );
  }
  else
    return(Gerror) ;
  return(NO_ERROR) ;
}
int
read_canon_vertex_coordinates(char *fname)
{
  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  if (MRISreadVertexPositions(mris, fname) == NO_ERROR)
  {
    MRISsaveVertexPositions(mris, CANONICAL_VERTICES) ;
    canonsurfloaded = TRUE ;
  }
  else
    canonsurffailed = TRUE ;

  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  if (canonsurfloaded)
    return(NO_ERROR);
  else
  {
    canonsurffailed = TRUE ;
    return(Gerror) ;
  }
}

void
read_ellipsoid_vertex_coordinates(char *fname,float a,float b,float c)
{
#if 0
  int k,n,num,dummy;
  float x,y,z,ctrx,ctry,ctrz,phi,theta;
  int ix,iy,iz;
  int version;
  FILE *fp;
  int first;

  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }

  /* marty */
  fread3(&first,fp);
  if (first == QUAD_FILE_MAGIC_NUMBER)
  {
    version = -1;
    printf("surfer: new surface file format\n");
  }
  else
  {
    rewind(fp);
    version = 0;
    printf("surfer: old surface file format\n");
  }
  fread3(&mris->nvertices,fp);
  fread3(&mris->nfaces,fp);
  printf("surfer: read_ellipsoid_vertex_coordinates(%s,%f,%f,%f)\n",
         fname,a,b,c);
  printf("surfer: vertices=%d, faces=%d\n",mris->nvertices,mris->nfaces);
  ctrx = ctry = ctrz = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    fread2(&ix,fp);
    fread2(&iy,fp);
    fread2(&iz,fp);
    x = ix/100.0;
    y = iy/100.0;
    z = iz/100.0;
    ctrx += x;
    ctry += y;
    ctrz += z;
  }
  ctrx /= mris->nvertices;
  ctry /= mris->nvertices;
  ctrz /= mris->nvertices;

  rewind(fp);
  fread3(&first,fp);
  if (first == QUAD_FILE_MAGIC_NUMBER)
  {
    version = -1;
  }
  else
  {
    rewind(fp);
    version = 0;
  }
  fread3(&mris->nvertices,fp);
  fread3(&mris->nfaces,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    fread2(&ix,fp);
    fread2(&iy,fp);
    fread2(&iz,fp);
    x = ix/100.0-ctrx;
    y = iy/100.0-ctry;
    z = iz/100.0-ctrz;
    phi = asin((y/b)/sqrt(SQR(x/a)+SQR(y/b)+SQR(z/c)));
    theta = atan2(x/a,z/c);
    mris->vertices[k].val = phi*180/M_PI;
    mris->vertices[k].val2 = theta*180/M_PI;
    mris->vertices[k].coords[0] = phi*180/M_PI;
    mris->vertices[k].coords[1] = 2+theta*180/M_PI;
    mris->vertices[k].coords[2] = 0;
    if (version==0)
    {
      fread1(&num,fp);
      for (n=0;n<num;n++)
        fread3(&dummy,fp);
    }
  }
  fclose(fp);
  isocontourflag = TRUE;
  contour_spacing[0] = contour_spacing[1] = contour_spacing[2] = 10;
#endif
}

static float e0[3] =
  {
    0.08, -0.73, 0.67
  } ;
static float e1[3] =
  {
    0.57, 0.58, 0.57
  } ;
static float e2[3] =
  {
    -0.82, 0.34, 0.47
  } ;

void
find_orig_vertex_coordinates(int vindex)
{
  float x,y,z;
  char  fname[NAME_LENGTH];
  float x_tal, y_tal, z_tal ;
  FILE  *fp ;
  int   error = 0 ;

  if (vindex < 0 || vindex >= mris->nvertices)
  {
    fprintf(stderr, "no vertex selected.\n") ;
    return ;
  }

  x_tal = y_tal = z_tal = 0.0 ;

  if (white_surf_loaded == FALSE){
    if (origsurfloaded == FALSE)
      {
	printf("surfer: reading original coordinates from\n");
	printf("surfer:   %s\n",orfname);
      }
    /* read coordinates from .orig file and put them in the .tx fields */
    if (origsurfloaded == FALSE &&
	read_orig_vertex_coordinates(orfname) != NO_ERROR)
      {
	error = 1 ;
	printf("surfer: wrong number of vertices/faces in file %s\n",orfname);
	PR printf("surfer: writing current coordinate (not orig) to file\n");
	PR x = mris->vertices[vindex].x ;
	y = mris->vertices[vindex].y ;
	z = mris->vertices[vindex].z ;
      }
    else  /* read file successfully */
      {
	x = mris->vertices[vindex].origx ;
	y = mris->vertices[vindex].origy ;
	z = mris->vertices[vindex].origz ;
      }
  } else {
    printf("Reading  coordinates from %s\n",white_suffix);
    x = mris->vertices[vindex].whitex ;
    y = mris->vertices[vindex].whitey ;
    z = mris->vertices[vindex].whitez ;
  }

  if (transform_loaded)
    conv_ras_to_tal(x, y, z, &x_tal, &y_tal, &z_tal) ;
  copy_edit_dat_file_name (fname, sizeof(fname));
  printf("writing coordinates to file %s\n", fname) ;
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  printf("vertex %d coordinates %5.2f %5.2f %5.2f:\n", 
	 vindex,x,y,z);
  if(!error){
    if (transform_loaded)
      printf("TALAIRACH (%2.1f %2.1f %2.1f)\n",x_tal,y_tal,z_tal);
    //printf("ORIGINAL  (%2.1f %2.1f %2.1f)\n",x,y,z);
  }
  else
    printf("CURRENT   (%2.1f %2.1f %2.1f)\n",x,y,z);
  fprintf(fp,"%f %f %f\n",x,y,z);
  fprintf(fp,"%f %f %f\n",x_tal,y_tal,z_tal);
  fclose(fp);


  if (canonsurfloaded == FALSE)
  {
    sprintf(fname, "%s.%s", fpref, sphere_reg_suffix) ;
    if (FileExists(fname))
    {
      printf("surfer: reading canonical coordinates from\n");
      printf("surfer:   %s\n",fname);
    }
  }

  if (canonsurfloaded == TRUE ||
      (FileExists(fname) &&
       read_canon_vertex_coordinates(fname) == NO_ERROR))
  {
    float sx, sy, sz, r, d, phi, theta ;

    x = mris->vertices[vindex].cx ;
    y = mris->vertices[vindex].cy ;
    z = mris->vertices[vindex].cz ;
    sx = x*e0[0] + y*e0[1] + z*e0[2] ;
    sy = x*e1[0] + y*e1[1] + z*e1[2] ;
    sz = x*e2[0] + y*e2[1] + z*e2[2] ;
    x = sx ;
    y = sy ;
    z = sz ;

    r = sqrt(x*x + y*y + z*z) ;
    d = r*r-z*z ;
    if (d < 0.0)
      d = 0.0 ;

    phi = atan2(sqrt(d), z) ;
    theta = atan2(y/r, x/r) ;
#if 0
#if 0
    if (theta < 0.0f)
      theta = 2 * M_PI + theta ;  /* make it 0 --> 2*PI */
#else
    if (theta > M_PI)
      theta -= 2 * M_PI ;  /* make it -PI --> PI */
#endif
#endif

    printf("SPHERICAL (%2.1f %2.1f %2.1f): (%2.1f, %2.1f)\n",sx,sy,sz,
           DEGREES(phi), DEGREES(theta));
  }
  PR
}

void copy_edit_dat_file_name (char* fname, int len)
{
  static int warned_local = 0;
  char* local_file  = NULL;
  int found = 0;
  FILE* test_file = NULL;
  char file_name[NAME_LENGTH] = "";

  /* First check if the local edit.dat file exists. If not, use the
     normal one. */
  found = FALSE;
  local_file = getenv ("FS_SAVE_GOTO_POINT");
  if (NULL != local_file)
  {

    sprintf (file_name, "%s-%s", local_file, pname);

    test_file = fopen (file_name, "a");
    if (test_file)
    {
      found = TRUE;
      fclose (test_file);

      if (!warned_local)
      {
        printf ("tksurfer: Using local edit.dat file %s\n", file_name);
        warned_local = TRUE;
      }
    }
  }

  if (!found)
  {

    /* Make the normal file name. */
    sprintf (file_name, "%s/edit.dat", tfname);
  }

  /* Return the file name. */
  strncpy (fname, file_name, len);
}


void
select_talairach_point(int *vindex,float x_tal,float y_tal,float z_tal)
{
  float x,y,z;
  char fname[NAME_LENGTH];
  FILE *fp;

  if (!transform_loaded)
  {
    printf("surfer: ### select_talairach_point failed: transform not "
           "loaded\n");
    PR return;
  }

  conv_tal_to_ras(x_tal, y_tal, z_tal, &x, &y, &z) ;

  copy_edit_dat_file_name (fname, sizeof(fname));
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fprintf(fp,"%f %f %f\n",x,y,z);
  fprintf(fp,"%f %f %f\n",x_tal,y_tal,z_tal);
  fclose(fp);

  select_orig_vertex_coordinates(vindex);
}

void
select_orig_vertex_coordinates(int *vindex)
{
  int   k;
  float d=0;
  float x,y,z,px,py,pz,mind;
  char  fname[NAME_LENGTH];
  FILE  *fp;

  copy_edit_dat_file_name (fname, sizeof(fname));
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  fscanf(fp,"%f %f %f\n",&px,&py,&pz);
  fclose(fp);

  if (selection>=0)
    draw_cursor(selection,FALSE);

  if (origsurfloaded == FALSE)
  {
    printf("surfer: reading coordinates from %s\n",orfname);
    read_orig_vertex_coordinates(orfname) ;
  }

  mind = 1e10;
  if (origsurfloaded == FALSE)
  {
    printf("surfer: ### wrong number of vertices/faces in file %s\n",fname);
    PR;
    printf("        using current position instead of orig...\n");
    PR;
    for (k = 0 ; k < mris->nvertices ; k++)
    {
      x=mris->vertices[k].x;
      y = mris->vertices[k].y;
      z = mris->vertices[k].z;
      d = SQR(x-px)+SQR(y-py)+SQR(z-pz);
      if (d<mind)
      {
        mind=d;
        *vindex=k;
      }
    }
  }
  else   /* find closest original vertex */
  {
    for (k=0;k<mris->nvertices;k++)
    {
      x = mris->vertices[k].origx ;
      y = mris->vertices[k].origy ;
      z = mris->vertices[k].origz ;
      d = SQR(x-px)+SQR(y-py)+SQR(z-pz);
      if (d<mind)
      {
        mind=d;
        *vindex=k ;
      }
    }
  }
  printf("surfer: vertex %d: dist = %f\n",*vindex,sqrt(mind));
  PR;

  /* begin rkt */
  if (selection>=0)
    draw_cursor(selection,TRUE);
  if (*vindex>=0)
    update_labels(LABELSET_CURSOR, *vindex, sqrt(mind));
  /* end rkt */

  print_vertex_data(*vindex, stdout, sqrt(mind)) ;
}

void
print_nearest_vertex_to_talairach_point(float x_tal, float y_tal, float z_tal)
{
  float x,y,z;
  int   k, vindex=-1;
  float d=0;
  float px=0.0,py=0.0,pz=0.0,mind;

  if (!transform_loaded)
  {
    printf("surfer: ### select_talairach_point failed: transform not "
           "loaded\n");
    PR return;
  }

  conv_tal_to_ras(x_tal, y_tal, z_tal, &px, &py, &pz) ;

  mind = 1e10;
  if (origsurfloaded == FALSE)
  {
    printf("surfer: ### print_nearest_vertex_to_talairach_point failed: orig surf is not loaded\n");
    PR return;
  }
  else
  {
    for (k=0;k<mris->nvertices;k++)
    {
      x = mris->vertices[k].origx ;
      y = mris->vertices[k].origy ;
      z = mris->vertices[k].origz ;
      d = SQR(x-px)+SQR(y-py)+SQR(z-pz);
      if (d<mind)
      {
        mind=d;
        vindex=k ;
      }
    }
  }

  printf("surfer: vertex %d\n", vindex);
  PR;
}




/* print curv,icurv,icurvnei,=>cfact */
void
read_curvim_at_vertex(int vindex)
{
  VERTEX *v;
  int imnr,i,j,m,n;
  int delcurvdefined;
  float x,y,z,sx,sy,sz,sd;
  float dx,dy,dz;
  float xnei,ynei,znei;
  float curv,icurv,icurvnei,cfact;
  float icrange,crange;

  icurv = icurvnei = 0.0f ;
  if (!curvloaded)
  {
    printf("surfer: ### curv not loaded!\n");
    PR return;
  }
  if (!curvimloaded)
  {
    printf("surfer: ### curvim not loaded!\n");
    PR return;
  }

  icrange = mris2->max_curv-mris2->min_curv;
  crange = mris->max_curv-mris->min_curv;

  v = &mris->vertices[vindex];
  x = v->x;
  y = v->y;
  z = v->z;
  sx=sy=sz=sd=0;
  n=0;
  delcurvdefined = TRUE;
  imnr = (int)((y-yy0)/st+0.5);
  i = (int)((zz1-z)/ps+0.5);
  j = (int)((xx1-x)/ps+0.5);
  imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
  i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
  j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
  if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
    icurv = bytetofloat(curvim[imnr][i][j],mris2->min_curv,mris2->max_curv);
  else
    delcurvdefined = FALSE;
  curv = v->curv;
  for (m=0;m<v->vnum;m++)
  {
    xnei = mris->vertices[v->v[m]].x;
    ynei = mris->vertices[v->v[m]].y;
    znei = mris->vertices[v->v[m]].z;
    imnr = (int)((ynei-yy0)/st+0.5);
    i = (int)((zz1-znei)/ps+0.5);
    j = (int)((xx1-xnei)/ps+0.5);
    imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
    i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
    j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
    if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
      icurvnei = bytetofloat(curvim[imnr][i][j],
                             mris2->min_curv,mris2->max_curv);
    else
      delcurvdefined = FALSE;
    /*curvnei = mris->vertices[v->v[m]].curv;*/   /* two dels?? */
    /* del im w/only sign of im/loc diff */
    cfact = 1.0;
    if (delcurvdefined)
      cfact += (icurvnei-icurv)/icrange
               * copysign(icstrength,mris2->min_curv+curv*(icrange/crange) -
                          icurv);
    sx += dx = (xnei - x)*cfact;
    sy += dy = (ynei - y)*cfact;
    sz += dz = (znei - z)*cfact;
    sd += sqrt(dx*dx+dy*dy+dz*dz);
    n++;
    if (delcurvdefined)
    {
      printf("surfer: ### nei vertex number: %d\n",n);
      PR printf("surfer: curv: %f\n",curv);
      PR printf("surfer: icurv: %f\n",icurv);
      PR printf("surfer: icurvnei: %f\n",icurvnei);
      PR printf("surfer: cfact: %f\n",cfact);
      PR
    }
  }
  if (!delcurvdefined)
  {
    printf("surfer: del curv not defined somewhere around this vertex\n");
    PR
  }
}

#if 0
find_nearest_vertices(rad)
/* marty: get peaked distribution of nearest pts */
float rad;
{
  FILE *fp;
  float xpt,ypt,zpt;
  float sqrad,sqdist;
  int k;

  fp=fopen(tfname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",tfname);
    PR return;
  }
  fscanf(fp,"%f %f %f",&xpt,&ypt,&zpt);
  fclose(fp);
  sqrad = SQR(rad);

  /* read vertex coords from orig surface file */

  for (k=0;k<mris->nvertices;k++)
  {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
    sqdist = SQR(x-xpt) + SQR(y-ypt) + SQR(z-zpt);
    if (dist < sqrad)
    {
      mris->vertices[k].val = ;
      mris->vertices[k].val2 = ;
    }
  }
}
#endif

int
write_binary_surface(char *fname)
{
#if 1
  MRISwrite(mris, fname) ;
#else
  int k, type;
  float x,y,z;
  FILE *fp;

  type = mrisFileNameType(fname) ;
  if (type == MRIS_ASCII_QUADRANGLE_FILE)
    return(-1/*MRISwriteAscii(mris, fname)*/) ;
  else if (type == MRIS_GEO_TRIANGLE_FILE)
    return(-1/*MRISwriteGeo(mris, fname)*/) ;
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return(-1);
  }
  else
  {
    return(MRISwriteTriangularSurface(fname)) ;
  }
  fwrite3(-1,fp);
  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces/2,fp);  /* # of quadrangles */
  for (k=0;k<mris->nvertices;k++)
  {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
    fwrite2((int)(x*100),fp);
    fwrite2((int)(y*100),fp);
    fwrite2((int)(z*100),fp);
  }
  for (k=0;k<mris->nfaces;k+=2)
  {
    fwrite3(mris->faces[k].v[0],fp);
    fwrite3(mris->faces[k].v[1],fp);
    fwrite3(mris->faces[k+1].v[0],fp);
    fwrite3(mris->faces[k].v[2],fp);
  }
  fclose(fp);
  printf("surfer: file %s written\n",fname);
  PR
#endif
  return(0) ;
}

void
write_binary_patch(char *fname)
{
#if 1
  MRISwritePatch(mris, fname) ;
#else
  int k,i,npts;
  float x,y,z;
  FILE *fp;

  npts = 0;
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag) npts++;
  printf("npts=%d\n",npts);
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fwriteInt(npts,fp);
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      i = (mris->vertices[k].border)?-(k+1):k+1;
      fwriteInt(i,fp);
      x = mris->vertices[k].x;
      y = mris->vertices[k].y;
      z = mris->vertices[k].z;
      fwrite2((int)(x*100),fp);
      fwrite2((int)(y*100),fp);
      fwrite2((int)(z*100),fp);
      /*
      printf("k=%d, i=%d\n",k,i);
      */
    }
  fclose(fp);
#endif
}

void
read_binary_patch(char *fname)
{
  if (mris->patch)
  {
#if 1
    MRISunrip(mris) ;
#else
    char mris_fname[500] ;

    strcpy(mris_fname, mris->fname) ;
    MRISfree(&mris) ;
    mris = MRISread(mris_fname) ;
#endif
  }
  MRISreadPatchNoRemove(mris, fname) ;
  {
    int vno ;
    VERTEX *v ;
    flag2d = TRUE;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      if (v->ripflag)
        continue ;
      if (!FZERO(v->z))
      {
        printf("surface not flattened - disabling 2d code...\n");
        flag2d = FALSE ;
        break ;
      }
    }
  }
  surface_compiled = 0 ;
  vertex_array_dirty = 1 ;
}

void
read_and_color_labeled_vertices(int r, int g, int b)
{
  meshr = r ;
  meshg = g ;
  meshb = b ;
  read_labeled_vertices(lfname) ;
  surface_compiled = 0 ;
}
/* begin rkt */

/* read_labeled_vertices now just calls labl_load to support multiple
   labels. */

#if 1

void
read_labeled_vertices(char *fname)
{
  labl_load (fname);
}

#else

void
read_labeled_vertices(char *fname)
{
  if (area)
    LabelFree(&area) ;
  area = LabelRead(pname,fname) ;
  if (!area)
    return ;
  if (!origsurfloaded)
    read_orig_vertex_coordinates(orfname) ;
  if (reassign)
    LabelUnassign(area) ;
  LabelFillUnassignedVertices(mris, area, WHITE_VERTICES);
  LabelMark(area, mris) ;
  {
    int n, nonzero = 0 ;

    for (n = 0 ; n < area->n_points ; n++)
      if (!FZERO(area->lv[n].stat))
        nonzero++ ;
    printf("********************      %d nonzero vertices found ********************\n", nonzero) ;
    if (nonzero== 0)
    {
      printf("label stat field identically zero - setting to 1\n") ;
      for (n = 0 ; n < area->n_points ; n++)
        area->lv[n].stat = 1 ;
    }
  }

  surface_compiled = 0 ;
  redraw() ;
}
#endif

/* end rkt */

/* begin rkt */

/* replaced write_labeled_vertices here so that what used to be done
   in LabelFromMarkedSurfaces is now done here, so we can fill the
   area->lv[n].stat field with the current overlay value instead of
   v->stat. */
void
write_labeled_vertices(char *fname)
{
  int    vno, npoints, n ;
  VERTEX *v ;
  LABEL* area;

  fprintf(stderr, "generating label from marked vertices...\n") ;

  /* count the number of marked vertices. */
  for (npoints = vno = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].marked)
      npoints++ ;

  /* if we didn't get any, return. */
  if (!npoints)
  {
    fprintf(stderr, "no marked vertices...\n") ;
    return;
  }

  /* allocate a label. */
  area = LabelAlloc(npoints, NULL, NULL) ;

  /* Copy the subject name. */
  strncpy( area->subject_name, pname, 100 );

  /* for every vertex, if it's marked, save its vertex coords,
     index, and fill the value of the current overlay. */
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (!v->marked)
      continue ;
    area->lv[n].x = v->x ;
    area->lv[n].y = v->y ;
    area->lv[n].z = v->z ;
    area->lv[n].vno = vno ;
    sclv_get_value (v, sclv_current_field, &(area->lv[n].stat) );
    n++ ;
  }
  area->n_points = npoints ;

  fprintf(stderr, "writing %d labeled vertices to %s.\n",
          area ? area->n_points : -1, fname) ;
  LabelToWhite(area, mris) ;
  LabelWrite(area, fname) ;
  LabelFree(&area) ;
}

#if 0
void
write_labeled_vertices(char *fname)
{
#if 1
  if (area)
  {
    if (origsurfloaded == FALSE)
    {
      fprintf(stderr, "reading original vertex locations...\n") ;
      MRISreadOriginalProperties(mris, NULL) ;
      origsurfloaded = TRUE ;
    }
    fprintf(stderr, "writing %d labeled vertices to %s.\n",
            area ? area->n_points : -1, fname) ;
    LabelToWhite(area, mris) ;
    LabelWrite(area, fname) ;
  }
  else
  {
    fprintf(stderr, "generating label from marked vertices...\n") ;
    area = LabelFromMarkedSurface(mris) ;
    if (!area)
    {
      fprintf(stderr, "no marked vertices...\n") ;
      return ;
    }
    fprintf(stderr, "writing %d labeled vertices to %s.\n",
            area ? area->n_points : -1, fname) ;
    LabelToWhite(area, mris) ;
    LabelWrite(area, fname) ;
    LabelFree(&area) ;
  }
#else
  int k,npts,asciiflag=TRUE;
  FILE *fp;

  read_orig_vertex_coordinates(orfname);
  npts = 0;
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag) npts++;
  printf("npts=%d\n",npts);
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  if (asciiflag)
  {
    fprintf(fp,"#!ascii\n");
    fprintf(fp,"%d\n",npts);
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
        fprintf(fp,"%d %f %f %f %f\n",
                k,mris->vertices[k].coords[0],mris->vertices[k].coords[1],
                mris->vertices[k].coords[2],mris->vertices[k].stat);
  }
  else
  {
    fputc('\0',fp);
    fwrite2(npts,fp);
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        fwriteInt(k,fp);
        fwriteFloat(mris->vertices[k].coords[0],fp);
        fwriteFloat(mris->vertices[k].coords[1],fp);
        fwriteFloat(mris->vertices[k].coords[2],fp);
        fwriteFloat(mris->vertices[k].stat,fp);
      }
  }
  fclose(fp);
  printf("surfer: file %s written\n",fname);
  PR
#endif
}
#endif
/* end rkt */


void
write_binary_curvature(char *fname)
{
#if 1
  MRISwriteCurvature(mris, fname) ;
#else
  int k;
  FILE *fp;

  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    fwrite2((int)(mris->vertices[k].curv*100),fp);
  }
  fclose(fp);
#endif
  printf("surfer: file %s written\n",fname);
  PR
}

void
write_binary_areas(char *fname)
{
#if 1
  MRISwriteArea(mris, fname) ;
#else
  int k;
  FILE *fp;

  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fwrite3(mris->nvertices,fp);
  fwrite3(mris->nfaces,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    fwriteFloat((mris->vertices[k].area),fp);
    mris->vertices[k].origarea = mris->vertices[k].area;
  }
  fclose(fp);
#endif
  printf("surfer: file %s written\n",fname);
  PR
}

/* 1/29/96: from paint.c */
void
write_binary_values(char *fname)
{
#if 1
  MRISwriteValues(mris, fname) ;
#else
  int k,num;
  float f;
  FILE *fp;
  double sum=0,sum2=0,max= -1000,min=1000;

  fp = fopen(fname,"wb");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  for (k=0,num=0;k<mris->nvertices;k++) if (mris->vertices[k].val!=0) num++;
  printf("num = %d\n",num);
  fwrite2(0,fp);
  fwrite3(num,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    if (mris->vertices[k].val!=0)
    {
      fwrite3(k,fp);
      f = mris->vertices[k].val;
      fwriteFloat(f,fp);
      sum += f;
      sum2 += f*f;
      if (f>max) max=f;
      if (f<min) min=f;
    }
  }
  fclose(fp);
  sum /= num;
  sum2 = sqrt(sum2/num-sum*sum);
  printf("avg = %f, stdev = %f, min = %f, max = %f\n",sum,sum2,min,max);
#endif
  printf("surfer: file %s written\n",fname);
  PR

}
void
read_stds(int cond_no)
{
  char  fname[STRLEN] ;
  int   dof, vno ;
  FILE  *fp  ;
  VERTEX *v ;
  double std ;

  sprintf(fname, "sigvar%d-%s.w", cond_no, stem) ;
  surface_compiled = 0 ;
  MRISreadValues(mris, fname) ;
  sprintf(fname, "sigavg%d.dof", cond_no) ;

  fp = fopen(fname, "r") ;
  if (!fp)
  {
    fprintf(stderr, "##tksurfer: could not open dof file %s\n", fname) ;
    return ;
  }
  fscanf(fp, "%d", &dof) ;
  fclose(fp) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    std = sqrt(v->val * (float)dof) ;
    v->val = std ;
  }
}

#include "mriFunctionalDataAccess.h"
#include "mri_transform.h"
#include "xVoxel.h"

/* begin rkt */
/* this is a dirty hack. i don't want to mess with read_binary_values too
 * too much. this version with the old signature read values into the v->val
 * field. it calls the new version passing SCLV_VAL. the new version is just
 * the old vesion with two differences. first, when reading a bfloat, it
 * switches on the field and sets the correct one. second, before calling
 * MRISreadValues, if the field is not SCLV_VAL, it saves all the v->vals,
 * and after, moves all the values in v->vals to the correct field, and
 * restores the values of v->val. something similar is done for
 * read_binary_values_frame.
 */

void
read_binary_values(char *fname)
{
  sclv_read_from_dotw (fname, SCLV_VAL);
}

int
sclv_read_from_dotw(char *fname, int field)  /* marty: openclose */
{
  int                   vno ;
  int                   error_code;
  VERTEX                *v ;
  float*                saved_vals = NULL;
  char                  val_name[STRLEN];
  char                  cmd[STRLEN];
  float                 min, max, mean;

  /* unload this field if it already exists */
  sclv_unload_field (field);

  /* save all the v->val values */
  if (field != SCLV_VAL)
  {
    saved_vals = (float*) calloc (mris->nvertices, sizeof(float));
    if (saved_vals==NULL)
      ErrorReturn(ERROR_NOMEMORY,
                  (ERROR_NOMEMORY,
                   "sclv_read_from_dotw: calloc with %d elmnts failed\n",
                   mris->nvertices));

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      saved_vals[vno] = v->val;
    }
  }

  /* save the directory for later */
  FileNamePath (fname, val_dir );

  /* read the file. if not found, bail. This sets all the v->val
     values to the new surface values. */
  error_code = MRISreadValues(mris, fname) ;
  if (error_code != NO_ERROR)
  {
    printf ("surfer: ### File %s could not be opened.\n", fname);
    return (error_code);
  }

  // remove the mean if specified by the user
  if (zero_mean)  
  {
    mean = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      mean += v->val ;
    }
    mean /= mris->nvertices ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->val -= mean ;
    }
  }
  /* look for the min and max values. */
  min = 1000000;
  max = -1000000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;

    /* look for the min or max value */
    if (v->val < min)
      min = v->val;
    if (v->val > max)
      max = v->val;
    mean += v->val ;
  }

  /* move the v->vals into the proper field and restore the saved
     values. */
  if (field != SCLV_VAL && saved_vals != NULL)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;

      /* look for the min or max value */
      if (v->val < min)
        min = v->val;
      if (v->val > max)
        max = v->val;

      /* set the field value from v->val */
      sclv_set_value (v, field, v->val);

      /* restore the value */
      v->val = saved_vals[vno];
    }
    free (saved_vals);
  }

  /* save the range */
  sclv_field_info[field].min_value = min;
  sclv_field_info[field].max_value = max;

  /* dummy info for time point and conditions, since .w files only have
     one plane of info */
  sclv_field_info[field].cur_timepoint = 0;
  sclv_field_info[field].cur_condition = 0;
  sclv_field_info[field].num_timepoints = 1;
  sclv_field_info[field].num_conditions = 1;

  /* calc the frquencies */
  sclv_calc_frequencies (field);

  /* request a redraw. turn on the overlay flag and select this value set */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (field);

  /* set the field name to the name of the file loaded */
  FileNameOnly (fname, val_name);
  sprintf (cmd, "UpdateValueLabelName %d \"%s\"", field, val_name);
  send_tcl_command (cmd);
  sprintf (cmd, "ShowValueLabel %d 1", field);
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup view");
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);

  /* enable the menu items */
  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);

  return (ERROR_NONE);
}

int
sclv_read_from_volume (char* fname, FunD_tRegistrationType reg_type,
                       char* registration, int field)
{
  FunD_tErr volume_error;
  mriFunctionalDataRef volume;
  char cmd[STRLEN];
  char val_name[STRLEN];
  Volm_tErr volm_err = Volm_tErr_NoErr;
  mriVolumeRef volm = NULL;
  int good = 0;

  /* unload this field if it already exists */
  sclv_unload_field (field);

  if (FunD_tRegistration_None == reg_type)
  {
    printf ("surfer: ERROR: Must specify registration type for overlay.\n"
            "Use -overlay-reg <file>, -overlay-reg-find, "
            "or -overlay-reg-identity.\n");
    return ERROR_BADPARM;
  }

  /* Unless they selected the none-needed registration method, we need
     to load up a volume for them to use as the base for the
     transform. So we use the orig header, which should already have
     been loaded, to get one. */
  if (FunD_tRegistration_NoneNeeded != reg_type)
  {
    good = 0;
    if (NULL != orig_mri_header)
    {
      volm_err = Volm_New (&volm);
      if (Volm_tErr_NoErr == volm_err)
      {
        volm_err = Volm_ImportData (volm, orig_mri_header->fname);
        if (Volm_tErr_NoErr == volm_err)
        {
          good = 1;
        }
      }
    }

    if (!good)
    {
      if (NULL != volm)
        Volm_Delete (&volm);
      printf ("surfer: ERROR: You specified registration type identity,\n"
              "but tksurfer cannot find an anatomical volume with which\n"
              "to calculate the identity transform. Please try another\n"
              "registration method.\n");
      return ERROR_BADPARM;
    }
  }

  if ((FunD_tRegistration_File == reg_type) && 
      (registration) &&
      (strlen(registration)==0))
  {
    printf ("surfer: sclv_read_from_volume,  ERROR: "
            "missing registration filename\n");
    return ERROR_BADPARM;
  }

  /* create volume. */
  volume_error = FunD_New (&volume,
                           fname,
                           reg_type,
                           registration,
                           mris->nvertices, /* Try to be scalar */
                           volm,
                           mris->hemisphere == LEFT_HEMISPHERE);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    printf("surfer: couldn't load %s.\n "
           "If you were trying to load a functional volume, make sure\n"
           "you selected the right registration method, and if necessary,\n"
           "that the registration file exists. \n"
           "If you were trying to load a volume-encoded value file,\n"
           "make sure it has the same number of values as this surface\n"
           "does vertices (%d).\n", fname, mris->nvertices);

    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_volume: error in FunD_New\n"));
  }

  /* Notify the volume we're in tkreg space, so our transform is
     correct. */
  if (zero_mean)
  {
    double mean = MRImeanFrame(volume->mpData, -1);
    printf("removing mean %2.2f from input volume\n", mean) ;
    MRIzeroMean(volume->mpData, volume->mpData) ;
  }
  volume_error = FunD_ClientSpaceIsTkRegRAS (volume);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    if (NULL != volume)
      FunD_Delete (&volume);
    sclv_unload_field (field);
    printf("surfer: couldn't load %s\n",fname);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_volume: error in "
                 "FunD_ClientSpaceIsTkRegRAS\n"));
  }

  // see if it is a correlation volume
  if (volume->mpData->width*volume->mpData->height*volume->mpData->depth == 
      2*mris->nvertices)
  {
    MRI *mri_tmp ;
    volume->mbScalar = 1 ;
    if (mris->hemisphere == LEFT_HEMISPHERE)
      mri_tmp = MRIextractInto(volume->mpData, NULL, 0, 0, 0, 
			       volume->mpData->width/2, volume->mpData->height, volume->mpData->depth,
			       0,0,0) ;
    else
      mri_tmp = MRIextractInto(volume->mpData, NULL, volume->mpData->width/2, 0, 0, 
			       volume->mpData->width/2, volume->mpData->height, volume->mpData->depth,
			       0,0,0) ;

    linkvertexmode = 1 ;
    fprintf(stderr, "setting linkvertexmode to %d\n", linkvertexmode) ;
    MRIfree(&volume->mpData) ; volume->mpData = mri_tmp ;
  }

  /* See if it's scalar */
  FunD_IsScalar (volume, &sclv_field_info[field].is_scalar_volume);

  if (sclv_field_info[field].is_scalar_volume)
  {
    printf ("surfer: Interpreting overlay volume %s "
            "as encoded scalar volume.\n", fname);
    if ((volume->mpData->nframes == mris->nvertices) ||
	(volume->mpData->nframes == 2*mris->nvertices))
    {
      char cmd[STRLEN] ;
      enable_menu_set (MENUSET_OVERLAY_LOADED, 1);
      enable_menu_set(MENUSET_TIMECOURSE_LOADED, 1) ;
      sprintf (cmd, "UpdateLinkedVarGroup overlay");
      send_tcl_command (cmd);
      linkvertexmode = 1 ;
      fprintf(stderr, "setting linkvertexmode to %d\n", linkvertexmode) ;
      func_load_timecourse(fname, FunD_tRegistration_Identity, NULL) ;
    }
  }
  else
  {
    printf ("surfer: Interpreting overlay volume %s "
            "as registered functional volume.\n", fname);
  }

  /* save the volume and mark this field as binary */
  sclv_field_info[field].is_functional_volume = TRUE;
  sclv_field_info[field].func_volume = volume;

  /* get the range information */
  volume_error =
    FunD_GetValueRange (sclv_field_info[field].func_volume,
                        &sclv_field_info[field].min_value,
                        &sclv_field_info[field].max_value);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    if (NULL != volume)
      FunD_Delete (&volume);
    sclv_unload_field (field);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_volume: error in FunD_GetValueRange\n"));
  }
  
  sclv_field_info[field].func_volume->mNumTimePoints= volume->mpData->nframes;
  volume_error =
    FunD_GetNumConditions (sclv_field_info[field].func_volume,
                           &sclv_field_info[field].num_conditions);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    if (NULL != volume)
      FunD_Delete (&volume);
    sclv_unload_field (field);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_volume: error in FunD_GetNumConditions\n"));
  }
  volume_error =
    FunD_GetNumTimePoints (sclv_field_info[field].func_volume,
                           &sclv_field_info[field].num_timepoints);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    if (NULL != volume)
      FunD_Delete (&volume);
    sclv_unload_field (field);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_volume: error in FunD_GetNumTimePoints\n"));
  }

  /* paint the first condition/timepoint in this field */
  sclv_set_timepoint_of_field (field, 0, 0);
  if (sclv_field_info[field].num_timepoints > 1)
  {
    enable_menu_set(MENUSET_TIMECOURSE_LOADED, 1) ;
    enable_menu_set(MENUSET_OVERLAY_LOADED, 1) ;
  }

  /* calc the frquencies */
  sclv_calc_frequencies (field);

  /* turn on the overlay flag and select this value set */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (field);

  /* set the field name to the name of the stem loaded */
  FileNameOnly (fname, val_name);
  sprintf (cmd, "UpdateValueLabelName %d \"%s\"", field, val_name);
  send_tcl_command (cmd);
  sprintf (cmd, "ShowValueLabel %d 1", field);
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup view");
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);

  /* enable the menu items */
  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);

  if (NULL != volm)
    Volm_Delete (&volm);

  return (ERROR_NONE);
}

int sclv_read_from_annotcorr (char* fname, int field)
{
  FunD_tErr volume_error;
  FunD_tErr eResult = FunD_tErr_NoError;
  FunD_tRegistrationType reg_type = FunD_tRegistration_NoneNeeded;
  char cmd[STRLEN];
  char val_name[STRLEN];
  int li_x, li_y, vnox, vnoy, annotformat=0;
  int lvnox, lvnoy;
  float val;
  MRI *mriannotoverlay, *mrinverts;
  LABEL *labelx, *labely;
  mriFunctionalDataRef this;
  mriVolumeRef volm = NULL;
  
  /* unload this field if it already exists */
  sclv_unload_field (field);
  
  this = (mriFunctionalDataRef) malloc(sizeof(mriFunctionalData));
  this->mSignature = FunD_kSignature;
  
  strcpy( this->msFileName, fname );
  /* Init values. */
  this->mSampleType = FunD_tSampleType_Nearest;
  this->mConvMethod = FunD_tConversionMethod_FFF;
  this->mNumTimePoints = -1;
  this->mNumConditions = -1;
  this->mbNullConditionPresent = FALSE;
  this->mpData = NULL;
  this->mMinValue = 10000;
  this->mMaxValue = -10000;
  this->mTimeResolution = 0;
  this->mNumPreStimTimePoints = 0;
  this->mIdxToIdxTransform = NULL;
  this->mOriginalIdxToIdxTransform = NULL;
  this->mbErrorDataPresent = FALSE;
  this->mCovMtx = NULL;
  this->mpResampledData = NULL;
  this->mClientXMin = 0;
  this->mClientYMin = 0;
  this->mClientZMin = 0;
  this->mClientXMax = -1;
  this->mClientYMax = -1;
  this->mClientZMax = -1;
  this->mbHaveClientBounds = FALSE;
  this->mbScalar = FALSE;

  this->mFrequencies = NULL;
  this->mNumBins = 0;
  
  /* Init transform objects */
  DebugNote( ("Creating idx to idx transform") );
  Trns_New( &(this->mIdxToIdxTransform) );

  /* Load the data. And create a nvertices x nvertices timepoint volume */
  mriannotoverlay = MRIread(fname);
  //printf("labl_num_labels : %d\n", labl_num_labels);
  //printf("mriannotoverlay->width : %d\n", mriannotoverlay->width);
  //printf("mriannotoverlay->height : %d\n", mriannotoverlay->height);
  /* the annotoverlay volume has to have num_labels x 1 x 1 x num_labels size 
   * or num_labels x num_labels x 1 x 1 size */
  if (mriannotoverlay->width > 1 && mriannotoverlay->height == 1 &&
      mriannotoverlay->depth == 1 && mriannotoverlay->nframes > 1 )
  {
    annotformat = 1;
    if ( mriannotoverlay->width != labl_num_labels || mriannotoverlay->nframes != labl_num_labels)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_read_from_annotcorr: every non-zero dimension of the annot overlay file must be the number of labels"));
  }
  else if (mriannotoverlay->width > 1 && mriannotoverlay->height > 1 &&
      mriannotoverlay->depth == 1 && mriannotoverlay->nframes == 1 )
  {
    annotformat = 2;
    if ( mriannotoverlay->width != labl_num_labels || mriannotoverlay->height != labl_num_labels)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_read_from_annotcorr: every non-zero dimension of the annot overlay file must be the number of labels"));
  }
  else
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_read_from_annotcorr: annot correlation matrix wrong format"));
  
  mrinverts  = MRIallocSequence( mris->nvertices, 1, 1, 
		                 mriannotoverlay->type,
				 mris->nvertices);
  MRIcopyHeader(mriannotoverlay, mrinverts);
  /* create a nverts x 1 x 1 x nverts volume from the annotation volume */
  for ( li_x = 0; li_x < labl_num_labels ; li_x++)
  {
    labelx = labl_labels[li_x].label;
    for (lvnox = 0; lvnox < labelx->n_points; lvnox++)
    {
      vnox = labelx->lv[lvnox].vno;
      if ( vnox < 0 || vnox >= mris->nvertices )
        continue;
      for ( li_y = 0; li_y < labl_num_labels ; li_y++)
      {
	if ( annotformat == 1)      
          val = MRIgetVoxVal(mriannotoverlay, li_x, 0, 0, li_y);
	else val = MRIgetVoxVal(mriannotoverlay, li_x, li_y, 0, 0);
        labely = labl_labels[li_y].label;
        for (lvnoy = 0; lvnoy < labely->n_points; lvnoy++)
        {
          vnoy = labely->lv[lvnoy].vno;
	  MRIsetVoxVal(mrinverts, vnox, 0, 0, vnoy, val);
        }
      }
    }
  }

  this->mpData = mrinverts;

  /* Try to parse a stem header if we have a stem specified. */
  DebugNote( ("Trying stem header") );
  eResult = FunD_FindAndParseStemHeader_( this );
  if ( FunD_tErr_NoError != eResult )
  {

    DebugNote( ("Guess meta information") );
    eResult = FunD_GuessMetaInformation_( this );
  }
  
  /* Try to reshape if we can. */
  DebugNote( ("Trying reshape with %d values", mris->nvertices) );
  FunD_ReshapeIfScalar_( this, mris->nvertices, NULL , mris->hemisphere == LEFT_HEMISPHERE);

  
  /* If we're not scalar by now, parse the registration file */
  if ( !this->mbScalar )
  {
    DebugNote( ("Parsing registration file") );
    eResult = FunD_ParseRegistrationAndInitMatricies_( this,
              reg_type,
              volm );
  }

  /* Get the value range. */
  DebugNote( ("Getting value range") );
  MRIvalRange( this->mpData, &this->mMinValue, &this->mMaxValue );

  volume_error = FunD_ClientSpaceIsTkRegRAS (this);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != this)
      FunD_Delete (&this);
    sclv_unload_field (field);
    printf("surfer: couldn't load %s\n",fname);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_volume: error in "
                 "FunD_ClientSpaceIsTkRegRAS\n"));
  }
  
  /* set the scalar volume field */
  sclv_field_info[field].is_scalar_volume = TRUE;
  
  printf ("surfer: Interpreting overlay volume %s "
            "as an encoded annotation correlation volume.\n", fname);
  
  /* save the volume and mark this field as binary */
  sclv_field_info[field].is_functional_volume = TRUE;
  sclv_field_info[field].func_volume = this;
  
  /* get the range information */
  volume_error =
    FunD_GetValueRange (sclv_field_info[field].func_volume,
                        &sclv_field_info[field].min_value,
                        &sclv_field_info[field].max_value);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != this)
      FunD_Delete (&this);
    sclv_unload_field (field);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_annotcorr: error in FunD_GetValueRange\n"));
  }
  volume_error =
    FunD_GetNumConditions (sclv_field_info[field].func_volume,
                           &sclv_field_info[field].num_conditions);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != this)
      FunD_Delete (&this);
    sclv_unload_field (field);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_annotcorr: error in FunD_GetNumConditions\n"));
  }
  volume_error =
    FunD_GetNumTimePoints (sclv_field_info[field].func_volume,
                           &sclv_field_info[field].num_timepoints);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != this)
      FunD_Delete (&this);
    sclv_unload_field (field);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "sclv_read_from_annotcorr: error in FunD_GetNumTimePoints\n"));
  }
  
  /* paint the first condition/timepoint in this field */
  sclv_set_timepoint_of_field (field, 0, 0);
  if (sclv_field_info[field].num_timepoints > 1)
  {
    enable_menu_set(MENUSET_TIMECOURSE_LOADED, 1) ;
    enable_menu_set(MENUSET_OVERLAY_LOADED, 1) ;
  }

  /* calc the frquencies */
  sclv_calc_frequencies (field);

  /* turn on the overlay flag and select this value set */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  colscale = HEAT_SCALE ;
  linkvertexmode = 1;
  sclv_set_current_field (field);

  /* set the field name to the name of the stem loaded */
  FileNameOnly (fname, val_name);
  sprintf (cmd, "UpdateValueLabelName %d \"%s\"", field, val_name);
  send_tcl_command (cmd);
  sprintf (cmd, "ShowValueLabel %d 1", field);
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup view");
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);

  /* enable the menu items */
  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);

  MRIfree(&mriannotoverlay);
  return (ERROR_NONE);
}

int sclv_new_from_label (int field, int label)
{

  VERTEX *v;
  LABEL  *area;
  int    n;
  float  min, max, f;
  char   cmd[STRLEN];

  if (field < 0 || field > NUM_SCALAR_VALUES)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_new_from_label: field was out of bounds: %d)",
                 field));

  if (label < 0 || label > labl_num_labels)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_new_from_label: label was out of bounds: %d)",
                 field));

  /* Unload this field if it already exists */
  sclv_unload_field (field);

  /* Get the label structure and go through the points. For each one,
     set our overlay value. Also check for min/max values. */
  area = labl_labels[label].label ;
  min = max = 0.0 ;
  for (n = 0  ; n < mris->nvertices ; n++)
  {
    v = &mris->vertices[n] ;
    sclv_set_value(v, field, 0) ;
  }
  for (n = 0  ; n < area->n_points ; n++)
  {
    if (area->lv[n].vno > 0 && area->lv[n].vno < mris->nvertices)
    {
      f = area->lv[n].stat ;
      v = &mris->vertices[area->lv[n].vno] ;

      /* Set the value. */
      sclv_set_value(v, field, f) ;

      /* Check for min/max. */
      if (n == 0)
        min = max = f ;
      if (f > max)
        max = f ;
      if (f < min)
        min = f ;
    }
  }

  /* Set up some initial stuff here. */
  sclv_field_info[field].min_value = min;
  sclv_field_info[field].max_value = max + epsilon;
  sclv_field_info[field].cur_timepoint = 0;
  sclv_field_info[field].cur_condition = 0;
  sclv_field_info[field].num_timepoints = 1;
  sclv_field_info[field].num_conditions = 1;

  /* calc the frquencies */
  sclv_calc_frequencies (field);

  /* request a redraw. turn on the overlay flag and select this value set */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (field);

  /* set the field name to the name of the file loaded */
  sprintf (cmd, "UpdateValueLabelName %d \"%s\"",
           field, labl_labels[labl_selected_label].name);
  send_tcl_command (cmd);
  sprintf (cmd, "ShowValueLabel %d 1", field);
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup view");
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);

  /* Enable menu items. */
  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);

  return (NO_ERROR);
}

int
sclv_new_empty (int field, char* name)
{
  int    vno;
  VERTEX *v;
  char   cmd[STRLEN];

  if (field < 0 || field > NUM_SCALAR_VALUES)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_new_empty: field was out of bounds: %d)",
                 field));


  /* Unload this field if it already exists */
  sclv_unload_field (field);

  /* Init the field to 0s. */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    sclv_set_value (v, field, 0);
  }

  /* Set up some initial stuff here. */
  sclv_field_info[field].min_value = 0;
  sclv_field_info[field].max_value = 0;
  sclv_field_info[field].cur_timepoint = 0;
  sclv_field_info[field].cur_condition = 0;
  sclv_field_info[field].num_timepoints = 1;
  sclv_field_info[field].num_conditions = 1;

  /* calc the frquencies */
  sclv_calc_frequencies (field);

  /* request a redraw. turn on the overlay flag and select this value set */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (field);

  /* set the field name to the name of the file loaded */
  if (NULL != name)
    sprintf (cmd, "UpdateValueLabelName %d \"%s\"", field, name);
  else
    sprintf (cmd, "UpdateValueLabelName %d \"New Layer\"", field);
  send_tcl_command (cmd);
  sprintf (cmd, "ShowValueLabel %d 1", field);
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup view");
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);

  /* Enable menu items. */
  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);

  return (NO_ERROR);
}

void
read_binary_values_frame(char *fname)
{
  sclv_read_from_dotw_frame (fname, SCLV_VAL);
}

int
sclv_read_from_dotw_frame(char *fname,
                          int field)  /* marty: open/leave open */
{
  int i,k,num;
  float f;
  float lat;
  VERTEX* v;

  surface_compiled = 0 ;

  if (!openvalfileflag)
  {
    fpvalfile = fopen(fname,"r");
    if (fpvalfile==NULL)
    {
      printf("surfer: ### File %s not found\n",fname);
      PR
      return(ERROR_NOFILE);
    }
    else
      openvalfileflag = TRUE;
  }

  fread2(&ilat,fpvalfile);
  lat = ilat/10.0;
  printf("surfer: latency = %6.1f\n",lat);

  for (k=0;k<mris->nvertices;k++)
  {
    v = &(mris->vertices[k]);
    sclv_set_value (v, field, 0);
#if 0
    mris->vertices[k].val=0;
#endif
  }

  fread3(&num,fpvalfile);
  printf("surfer: num=%d\n",num);
  for (i=0;i<num;i++)
  {
    fread3(&k,fpvalfile);
    f = freadFloat(fpvalfile);
    if (k>=mris->nvertices||k<0)
      printf("surfer: vertex index out of range: %d f=%f\n",k,f);
    else if (mris->vertices[k].d!=0)
      printf("surfer: subsample and data file mismatch\n");
    else
    {
      v = &(mris->vertices[k]);
      sclv_set_value (v, field, f);
#if 0
      mris->vertices[k].val = f;
#endif
      mris->vertices[k].d=0;
    }
  }
  PR

  return (ERROR_NONE);
}

int
sclv_write_dotw(char *fname, int field)
{
  float *saved_vals = NULL;
  VERTEX* v = NULL;
  int vno;

  /* first save all the current val values */
  if (field != SCLV_VAL)
  {
    saved_vals = (float*) calloc (mris->nvertices, sizeof(float));
    if (saved_vals==NULL)
      ErrorReturn(ERROR_NOMEMORY,
                  (ERROR_NOMEMORY,
                   "sclv_write_dotw: calloc with %d elmnts failed\n",
                   mris->nvertices));

    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;

      /* save the value. */
      saved_vals[vno] = v->val;

      /* write the val from the requested field into val field */
      sclv_get_value (v, field, &(v->val));
    }
  }

  /* write the values */
  MRISwriteValues(mris, fname) ;

  /* restore the saved values */
  if (field != SCLV_VAL && saved_vals != NULL)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      v->val = saved_vals[vno];
    }
    free (saved_vals);
  }

  printf("surfer: file %s written\n",fname);
  PR

  return (ERROR_NONE);
}


void
swap_stat_val(void)
{
  sclv_swap_fields (SCLV_VAL, SCLV_VALSTAT);

#if 0
  int k;
  float tmp;
  /* begin rkt */
  char cmd[STRLEN] ;
  /* end rkt */
  statflag = TRUE;
  for (k=0;k<mris->nvertices;k++)
  {
    tmp = mris->vertices[k].stat;
    mris->vertices[k].stat = mris->vertices[k].val;
    mris->vertices[k].val = tmp;
  }

  /* begin rkt */
  /* swap the names of the val and stat labels */
  sprintf (cmd, "SwapValueLabelNames %d %d", SCLV_VAL, SCLV_VALSTAT);
      load_curv = TRUE;
      forcegraycurvatureflag = TRUE;
  send_tcl_command (cmd);
  /* end rkt */
#endif
}

void
swap_val_val2(void)
{
  sclv_swap_fields (SCLV_VAL, SCLV_VAL2);

#if 0
  int k;
  float tmp;

  for (k=0;k<mris->nvertices;k++)
  {
    tmp = mris->vertices[k].val2;
    mris->vertices[k].val2 = mris->vertices[k].val;
    mris->vertices[k].val = tmp ;
  }
  {
    char cmd[STRLEN] ;
    sprintf (cmd, "SwapValueLabelNames %d %d", SCLV_VAL, SCLV_VAL2);
    send_tcl_command (cmd);
  }
#endif
}

void
shift_values(void)   /* push one: val -> val2 */
{
  int k;

  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].val2 = mris->vertices[k].val;
  {
    char cmd[STRLEN] ;
    sprintf (cmd, "SwapValueLabelNames %d %d", SCLV_VAL, SCLV_VAL2);
    send_tcl_command (cmd);
  }
}

void
swap_values(void)   /* swap complex: val,valbak <-> val2,val2bak */
{
  sclv_swap_fields (SCLV_VAL, SCLV_VALBAK);
  sclv_swap_fields (SCLV_VAL2, SCLV_VAL2BAK);

#if 0
  int k;
  float h;

  for (k=0;k<mris->nvertices;k++)
  {
    h = mris->vertices[k].valbak;
    mris->vertices[k].valbak = mris->vertices[k].val;
    mris->vertices[k].val = h;
    h = mris->vertices[k].val2bak;
    mris->vertices[k].val2bak = mris->vertices[k].val2;
    mris->vertices[k].val2 = h;
  }
  {
    char cmd[STRLEN] ;
    sprintf (cmd, "SwapValueLabelNames %d %d", SCLV_VAL, SCLV_VALBAK);
    send_tcl_command (cmd);
    sprintf (cmd, "SwapValueLabelNames %d %d", SCLV_VAL2BAK, SCLV_VAL2);
    send_tcl_command (cmd);
  }
#endif
}

void
read_binary_decimation(char *fname)
{
  int vno ;
  surface_compiled = 0 ;
  sol_ndec = MRISreadDecimation(mris, fname) ;
  MRISclearMarks(mris) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].fixedval == TRUE)
      mris->vertices[vno].marked = 1 ;
  printf("surfer: decimation file %s read\n",fname);
  PR
}

void
write_binary_decimation(char *fname)
{
  MRISwriteDecimation(mris, fname) ;
  printf("surfer: decimation file %s written\n",fname);
  PR
}

void
write_decimation(char *fname)
{
  int k;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fprintf(fptr,"#!ascii\n");
  fprintf(fptr,"%d\n",mris->nvertices);
  for (k=0;k<mris->nvertices;k++)
  {
    if (mris->vertices[k].d==0)
      fprintf(fptr,"1\n");
    else
      fprintf(fptr,"0\n");
  }
  fclose(fptr);
  printf("surfer: decimation file %s written\n",fname);
  PR
}

void
read_binary_dipoles(char *fname)
{
#if 0
  int k,d;
  char c;
  FILE *fptr;

  fptr = fopen(fname,"r");
  if (fptr==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  c = fgetc(fptr);
  if (c=='#')
  {
    fscanf(fptr,"%*s");
    fscanf(fptr,"%d",&d);
    if (d!=mris->nvertices)
    {
      printf("surfer: ### diople file mismatch %s\n",fname);
      PR return;
    }
    for (k=0;k<mris->nvertices;k++)
    {
      fscanf(fptr,"%f %f %f %f %f %f",
             &mris->vertices[k].dipx,
             &mris->vertices[k].dipy,
             &mris->vertices[k].dipz,
             &mris->vertices[k].dipnx,
             &mris->vertices[k].dipny,
             &mris->vertices[k].dipnz);
    }
  }
  else
  {
    d = freadInt(fptr);
    if (d!=mris->nvertices)
    {
      printf("surfer: ### dipole file mismatch %s\n",fname);
      PR return;
    }
    for (k=0;k<mris->nvertices;k++)
    {
      mris->vertices[k].dipx = freadFloat(fptr);
      mris->vertices[k].dipy = freadFloat(fptr);
      mris->vertices[k].dipz = freadFloat(fptr);
      mris->vertices[k].dipnx = freadFloat(fptr);
      mris->vertices[k].dipny = freadFloat(fptr);
      mris->vertices[k].dipnz = freadFloat(fptr);
    }
  }
  fclose(fptr);
  printf("surfer: dipole file %s read\n",fname);
  PR
#endif
}

void
write_binary_dipoles(char *fname)
{
  int i,k;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fputc('\0',fptr);
  fwriteInt(mris->nvertices,fptr);
  for (k=0;k<mris->nvertices;k++)
  {
    fwriteFloat(mris->vertices[k].x,fptr);
    fwriteFloat(mris->vertices[k].y,fptr);
    fwriteFloat(mris->vertices[k].z,fptr);
    fwriteFloat(mris->vertices[k].nx,fptr);
    fwriteFloat(mris->vertices[k].ny,fptr);
    fwriteFloat(mris->vertices[k].nz,fptr);
  }
  for (k=0;k<mris->nvertices;k++)
  {
    fwriteInt(mris->vertices[k].vnum,fptr);
    for (i=0;i<mris->vertices[k].vnum;i++)
      fwriteInt(mris->vertices[k].v[i],fptr);
  }
  fclose(fptr);
  printf("surfer: dipole file %s written\n",fname);
  PR
}

void
write_dipoles(char *fname)
{
  int i,k;
  float x,y,z;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fprintf(fptr,"#!ascii\n");
  fprintf(fptr,"%d\n",mris->nvertices);
  for (k=0;k<mris->nvertices;k++)
  {
    x = mris->vertices[k].x;
    y = mris->vertices[k].y;
    z = mris->vertices[k].z;
    fprintf(fptr,"%6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
            x,y,z,
            mris->vertices[k].nx,
            mris->vertices[k].ny,
            mris->vertices[k].nz);
  }
  for (k=0;k<mris->nvertices;k++)
  {
    fprintf(fptr,"%d ",mris->vertices[k].vnum);
    for (i=0;i<mris->vertices[k].vnum;i++)
      fprintf(fptr,"%d ",mris->vertices[k].v[i]);
    fprintf(fptr,"\n");
  }
  fclose(fptr);
  printf("surfer: dipole file %s written\n",fname);
  PR
}

void
write_subsample(char *fname)
{
  int k;
  float x,y,z;
  FILE *fptr;

  fptr = fopen(fname,"w");
  if (fptr==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fprintf(fptr,"%d\n",sub_num);
  for (k=0;k<mris->nvertices;k++)
  {
    if (mris->vertices[k].d==0)
    {
      x = mris->vertices[k].x;
      y = mris->vertices[k].y;
      z = mris->vertices[k].z;
      fprintf(fptr,"%d %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
              k,
              x,y,z,
              mris->vertices[k].nx,
              mris->vertices[k].ny,
              mris->vertices[k].nz);
    }
  }
  fclose(fptr);
}

void
write_binary_subsample(char *fname)
{
  int k;
  float x,y,z,nx,ny,nz;
  FILE *fp;

  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fwrite3(sub_num,fp);
  for (k=0;k<mris->nvertices;k++)
  {
    if (mris->vertices[k].d==0)
    {
      x = mris->vertices[k].x;
      y = mris->vertices[k].y;
      z = mris->vertices[k].z;
      nx = mris->vertices[k].nx;
      ny = mris->vertices[k].ny;
      nz = mris->vertices[k].nz;
      fwrite3(k,fp);
      fwrite2((int)(x*100),fp);
      fwrite2((int)(y*100),fp);
      fwrite2((int)(z*100),fp);
      fwrite2((int)(nx*10000),fp);
      fwrite2((int)(ny*10000),fp);
      fwrite2((int)(nz*10000),fp);
    }
  }
  fclose(fp);
  printf("surfer: file %s written\n",fname);
  PR
}

void
read_binary_subsample(char *fname)
{
  int i,k;
  int ix,iy,iz,inx,iny,inz;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  fread3(&sub_num,fp);
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].d = 1;
  for (i=0;i<sub_num;i++)
  {
    fread3(&k,fp);
    fread2(&ix,fp);
    fread2(&iy,fp);
    fread2(&iz,fp);
    fread2(&inx,fp);
    fread2(&iny,fp);
    fread2(&inz,fp);
    mris->vertices[k].d = 0;
  }
  fclose(fp);
}

void
read_binary_curvature(char *fname)
{
#if 1
  int    vno, n ;
  VERTEX *v ;
  int    error;

  surface_compiled = 0 ;
  error = MRISreadCurvatureFile(mris, fname) ;
  if (NO_ERROR != error)
  {
    printf ("surfer: error reading curvature file\n");
    return;
  }
  for (dipavg = 0.0, n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    dipavg += v->curv ;
    n++ ;
  }
  dipavg /= (float)n ;
#else
  int k,i,vnum,fnum;
  float curv;
  FILE *fp;

  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  fread3(&vnum,fp);
  fread3(&fnum,fp);
  if (vnum!=mris->nvertices)
  {
    printf("surfer: ### incompatible vertex number in file %s\n",fname);
    PR return;
  }
  for (k=0;k<vnum;k++)
  {
    fread2(&i,fp);
    curv = i/100.0;
    if (k==0) mris->min_curv=mris->max_curv=curv;
    if (curv>mris->max_curv) mris->max_curv=curv;
    if (curv<mris->min_curv) mris->min_curv=curv;
    mris->vertices[k].curv = curv;
  }
  fclose(fp);
#endif
  printf("surfer: curvature read: min=%f max=%f\n",
         mris->min_curv,mris->max_curv);
  PR curvloaded = TRUE;
  /* begin rkt */
  /* save the curv min and max. */
  cmin = mris->min_curv;
  cmax = mris->max_curv;

  /* enable our menu options. */
  enable_menu_set (MENUSET_CURVATURE_LOADED, curvloaded);
  send_tcl_command ("ShowLabel kLabel_Curvature 1");

  /* turn on the curvflag */
  curvflag = 1;

  /* end rkt */
}

void
normalize_binary_curvature(void)
{
  int k;
  float curv,min,max;
  float sum,avg,sum_sq,sd,n;

  min = max = 0.0f ;
  if (!curvloaded)
  {
    printf("surfer: ### curv not loaded!\n");
    PR return;
  }

  sum = 0;
  for (k=0;k<mris->nvertices;k++)
    sum += mris->vertices[k].curv;
  avg = sum/mris->nvertices;

  n = (float)mris->nvertices;
  sum = sum_sq = 0.0;
  for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].curv -= avg;
    curv = mris->vertices[k].curv;
    sum += curv;
    sum_sq += curv*curv;
  }
  sd = sqrt((n*sum_sq - sum*sum)/(n*(n-1.0)));

  for (k=0;k<mris->nvertices;k++)
  {
    curv = (mris->vertices[k].curv)/sd;
    if (k==0) min=max=curv;
    if (curv<min) min=curv;
    if (curv>max) max=curv;
    if (curv<CURVIM_NORM_MIN) curv = CURVIM_NORM_MIN;
    if (curv>CURVIM_NORM_MAX) curv = CURVIM_NORM_MAX;
    mris->vertices[k].curv = curv;
  }
  mris->min_curv = CURVIM_NORM_MIN;
  mris->max_curv = CURVIM_NORM_MAX;
  printf("surfer: curvature normalized: avg=%f sd=%f\n",avg,sd);
  printf("surfer: min=%f max=%f trunc to (%f,%f)\n",
         min,max,mris->min_curv,mris->max_curv);
  PR
}

void
read_binary_areas(char *fname)
{
#if 1
  surface_compiled = 0 ;
  MRISreadBinaryAreas(mris, fname) ;
#else
  int k,vnum,fnum;
  float f;
  FILE *fp;

  mris->total_area = 0;
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: no area file %s\n",fname);
    PR return;
  }
  fread3(&vnum,fp);
  fread3(&fnum,fp);
  if (vnum!=mris->nvertices)
  {
    printf("surfer: ### incompatible vertex number in file %s\n",fname);
    PR printf("   ...file ignored\n");
    PR return;
  }
  for (k=0;k<vnum;k++)
  {
    f = freadFloat(fp);
    mris->vertices[k].origarea = f;
    mris->total_area += f;
  }
  fclose(fp);
  mris->total_area /= 2;
  /* hack to correct overest. of area in compute_normals */
#endif
}

void
read_fieldsign(char *fname)
{
  int k,vnum;
  float f;
  FILE *fp;

  printf("surfer: read_fieldsign(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  vnum = freadInt(fp);
  printf("surfer: mris->nvertices = %d, vnum = %d\n",mris->nvertices,vnum);
  if (vnum!=mris->nvertices)
  {
    printf("surfer: ### incompatible vertex number in file %s\n",fname);
    printf("   ...file read anyway\n");
  }
  for (k=0;k<mris->nvertices;k++)
  {
    f = freadFloat(fp);
    mris->vertices[k].fieldsign = f;
  }
  fclose(fp);
  fieldsignflag = TRUE;

  enable_menu_set (MENUSET_FIELDSIGN_LOADED, 1);
  send_tcl_command ("ShowLabel kLabel_Fieldsign 1");
}

void write_fieldsign(char *fname)
{
  int k,vnum;
  float f;
  FILE *fp;
  MRI *mri;
  char tmpstr[2000];

  vnum = mris->nvertices;
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fwriteInt(vnum,fp);
  printf("surfer: mris->nvertices = %d, vnum = %d\n",mris->nvertices,vnum);
  for (k=0;k<vnum;k++)
  {
    f = mris->vertices[k].fieldsign;
    fwriteFloat(f,fp);
  }
  fclose(fp);
  printf("surfer: file %s written\n",fname);

  // This writes it out as an mgz file
  sprintf(tmpstr,"%s.mgz",fname);
  mri = MRIcopyMRIS(NULL, mris, 0, "fieldsign");
  MRIwrite(mri,tmpstr);
  MRIfree(&mri);

  PR
}

void read_fsmask(char *fname)
{
  int k,vnum;
  float f;
  FILE *fp;

  printf("surfer: read_fsmask(%s)\n",fname);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }
  vnum = freadInt(fp);
  if (vnum!=mris->nvertices)
    printf("surfer: warning: incompatible vertex number in file %s\n",fname);
  for (k=0;k<mris->nvertices;k++)
  {
    f = freadFloat(fp);
    mris->vertices[k].fsmask = f;
  }
  fclose(fp);
  PR

  enable_menu_set (MENUSET_FIELDMASK_LOADED, 1);
}

void write_fsmask(char *fname)
{
  int k,vnum;
  float f;
  FILE *fp;
  MRI *mri;
  char tmpstr[2000];


  vnum = mris->nvertices;
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  fwriteInt(vnum,fp);
  for (k=0;k<vnum;k++)
  {
    f = mris->vertices[k].fsmask;
    fwriteFloat(f,fp);
  }
  fclose(fp);
  printf("surfer: file %s written\n",fname);

  // This writes it out as an mgz file
  sprintf(tmpstr,"%s.mgz",fname);
  mri = MRIcopyMRIS(NULL, mris, 0, "fsmask");
  MRIwrite(mri,tmpstr);
  MRIfree(&mri);

  PR
}

void
write_vrml(char *fname,int mode)
{
  FILE *fp;
  int k,n;

  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }

  if (mode==VRML_POINTS)
  {
    /* color here */
    if (vrml2flag)
    {
      fprintf(fp,"#VRML V2.0 utf8\n");
      fprintf(fp,"Collision { collide FALSE\n");
      fprintf(fp,"children [ Group { children [ Shape {\n");
      fprintf(fp,"appearance Appearance { "
              "material DEF _DefMat Material { } }\n");
      fprintf(fp,"geometry PointSet { coord Coordinate {\n");
    }
    else
    {
      fprintf(fp,"#VRML V1.0 ascii\n");
      fprintf(fp,"Coordinate3 {\n");
    }

    fprintf(fp,"\tpoint [\n");
    for (k=0;k<mris->nvertices;k++)
    {
      fprintf(fp,"%3.2f %3.2f %3.2f",
              mris->vertices[k].x,
              mris->vertices[k].y,
              mris->vertices[k].z);
      if (k!=mris->nvertices-1) fprintf(fp,",\n");
    }
    fprintf(fp,"\n\t]\n");

    if (vrml2flag)
    {
      fprintf(fp,"} } } ] } ] }\n");
    }
    else
    {
      fprintf(fp,"}\n");
      fprintf(fp,"PointSet {\n");
      fprintf(fp,"\tstartIndex 0\n");
      fprintf(fp,"\tnumPoints -1\n");
      fprintf(fp,"}\n");
    }
  }
  if (mode==VRML_QUADS_GRAY || mode==VRML_QUADS_GRAY_SMOOTH)
  {
    fprintf(fp,"#VRML V1.0 ascii\n");
    fprintf(fp,"Coordinate3 {\n");
    fprintf(fp,"\tpoint [\n");
    for (k=0;k<mris->nvertices;k++)
    {
      fprintf(fp,"%3.2f %3.2f %3.2f",
              mris->vertices[k].x,
              mris->vertices[k].y,
              mris->vertices[k].z);
      if (k!=mris->nvertices-1) fprintf(fp,",\n");
    }
    fprintf(fp,"\n\t]\n");
    fprintf(fp,"}\n");

    if (mode==VRML_QUADS_GRAY_SMOOTH)
    {
      fprintf(fp,"ShapeHints {\n");
      /*fprintf(fp,"\tvertexOrdering COUNTERCLOCKWISE\n");*/
      /*fprintf(fp,"\tshapeType SOLID\n");*/  /* no backface speedup */
      fprintf(fp,"\tcreaseAngle 0.8\n");
      fprintf(fp,"}\n");
    }

    fprintf(fp,"IndexedFaceSet {\n");
    fprintf(fp,"\tcoordIndex [\n");
    for (k=0;k<mris->nfaces;k++)
    {
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        fprintf(fp,"%d",mris->faces[k].v[n]);
        if (n<3)
        {
          fprintf(fp,", ");
        }
        else
        {
          if (k!=mris->nfaces-1) fprintf(fp,", -1,\n");
        }
      }
    }
    fprintf(fp,"\n\t]\n");
    fprintf(fp,"}\n");
  }
  if (mode==VRML_QUADS_CURV)
  {
    if (!curvloaded) read_binary_curvature(cfname);
#if 0
    /*** TODO ***/
    fprintf(fp,"#VRML V1.0 ascii\n");
    fprintf(fp,"Coordinate3 {\n");
    fprintf(fp," point [\n");
    for (k=0;k<mris->nvertices;k++)
    {
      fprintf(fp,"%3.2f %3.2f %3.2f",
              mris->vertices[k].x,
              mris->vertices[k].y,
              mris->vertices[k].z);
      if (k!=mris->nvertices-1) fprintf(fp,",\n");
    }
    fprintf(fp," ]\n");
    fprintf(fp,"}\n");

    fprintf(fp,"Material {\n");
    fprintf(fp," diffuseColor [\n");
    fprintf(fp," ]\n");
    fprintf(fp,"}\n");

    fprintf(fp,"MaterialBinding {\n");
    fprintf(fp," value PER_VERTEX\n");
    /* kill overlap: PER_MRIS->NVERTICESED */
    fprintf(fp,"}\n");

    fprintf(fp,"IndexedFaceSet {\n");
    fprintf(fp," coordIndex [\n");
    for (k=0;k<mris->nfaces;k++)
    {
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        fprintf(fp,"%d",mris->faces[k].v[n]);
        if (n<3)
        {
          fprintf(fp,", ");
        }
        else
        {
          if (k!=mris->nfaces-1) fprintf(fp,", -1,\n");
        }
      }
    }
    fprintf(fp," ]\n");
    fprintf(fp,"}\n");
    fprintf(fp," materialIndex [\n");
    /* ... */
    fprintf(fp," ]\n");
    fprintf(fp,"}\n");
#endif
  }
  fclose(fp);
  if (vrml2flag) printf("surfer: vrml 2.0 file %s written",fname);
  else           printf("surfer: vrml 1.0 file %s written",fname);
  PR
}

void
rip_faces(void)
{
  int n,k;
  FACE *f;

  for (k=0;k<mris->nfaces;k++)
    mris->faces[k].ripflag = FALSE;
  for (k=0;k<mris->nfaces;k++)
  {
    f = &mris->faces[k];
    for (n=0;n<VERTICES_PER_FACE;n++)
      if (mris->vertices[f->v[n]].ripflag)
        f->ripflag = TRUE;
  }
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].border = FALSE;
  for (k=0;k<mris->nfaces;k++)
    if (mris->faces[k].ripflag)
    {
      f = &mris->faces[k];
      for (n=0;n<VERTICES_PER_FACE;n++)
        mris->vertices[f->v[n]].border = TRUE;
    }
}

void
normalize(float v[3])
{
  float d;

  d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if (d>0)
  {
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
  }
}

void
normal_face(int fac,int n,float *norm)
{
  int n0,n1;
  FACE *f;
  float v0[3],v1[3];

  n0 = (n==0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n==VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
  normalize(v0);
  normalize(v1);
  norm[0] = -v1[1]*v0[2] + v0[1]*v1[2];
  norm[1] = v1[0]*v0[2] - v0[0]*v1[2];
  norm[2] = -v1[0]*v0[1] + v0[0]*v1[1];
  /*
    printf("[%5.2f,%5.2f,%5.2f] x [%5.2f,%5.2f,%5.2f] = [%5.2f,%5.2f,%5.2f]\n",
    v0[0],v0[1],v0[2],v1[0],v1[1],v1[2],norm[0],norm[1],norm[2]);
  */
}

float triangle_area(int fac,int n)
{
  int n0,n1;
  FACE *f;
  float v0[3],v1[3],d1,d2,d3;

  n0 = (n==0)                   ? VERTICES_PER_FACE-1 : n-1;
  n1 = (n==VERTICES_PER_FACE-1) ? 0                   : n+1;
  f = &mris->faces[fac];
  v0[0] = mris->vertices[f->v[n]].x - mris->vertices[f->v[n0]].x;
  v0[1] = mris->vertices[f->v[n]].y - mris->vertices[f->v[n0]].y;
  v0[2] = mris->vertices[f->v[n]].z - mris->vertices[f->v[n0]].z;
  v1[0] = mris->vertices[f->v[n1]].x - mris->vertices[f->v[n]].x;
  v1[1] = mris->vertices[f->v[n1]].y - mris->vertices[f->v[n]].y;
  v1[2] = mris->vertices[f->v[n1]].z - mris->vertices[f->v[n]].z;
  d1 = -v1[1]*v0[2] + v0[1]*v1[2];
  d2 = v1[0]*v0[2] - v0[0]*v1[2];
  d3 = -v1[0]*v0[1] + v0[0]*v1[1];
  return sqrt(d1*d1+d2*d2+d3*d3)/2;
}

#define MAX_NEIGHBORS 300  /* ridiculously large */


void
find_neighbors(void)
{
  int n0,n1,i,k,m,n;
  FACE *f;
  VERTEX *v;

#if 1
  mrisCompleteTopology(mris);

#else
  //
  // This code is a duplicate of mrisFindNeighbors2 in mri_mc.c and several other places
  // and that operation is now being done in one place - mrisCompleteTopology
  //
  int vtmp[MAX_NEIGHBORS];

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->vnum = 0;
    for (m=0;m<v->num;m++)
    {
      n = v->n[m];
      f = &mris->faces[v->f[m]];
      n0 = (n==0)                   ? VERTICES_PER_FACE-1 : n-1;
      n1 = (n==VERTICES_PER_FACE-1) ? 0                   : n+1;
      for (i=0;i<v->vnum && vtmp[i]!=f->v[n0];i++);
      if (i==v->vnum)
        vtmp[(int)v->vnum++] = f->v[n0];
      for (i=0;i<v->vnum && vtmp[i]!=f->v[n1];i++);
      if (i==v->vnum)
        vtmp[(int)v->vnum++] = f->v[n1];
    }
    mris->vertices[k].v = (int *)lcalloc(mris->vertices[k].vnum,sizeof(int));
    for (i=0;i<v->vnum;i++)
    {
      v->v[i] = vtmp[i];
    }
    /*
      if (v->num != v->vnum)
      printf("%d: num=%d vnum=%d\n",k,v->num,v->vnum);
    */
  }
#endif
  for (k=0;k<mris->nfaces;k++)
  {
    f = &mris->faces[k];
    for (m=0;m<VERTICES_PER_FACE;m++)
    {
      v = &mris->vertices[f->v[m]];
      for (i=0;i<v->num && k!=v->f[i];i++);
      if (i==v->num)
        printf("mris->faces[%d].v[%d] = %d\n",k,m,f->v[m]);
    }
  }
}

void
compute_normals(void)  /* no triangle area in msurfer, no explodeflag here */
{
  int k,n;
  VERTEX *v;
  FACE *f;
  float norm[3],snorm[3];

  for (k=0;k<mris->nfaces;k++)
    if (mris->faces[k].ripflag)
    {
      f = &mris->faces[k];
      for (n=0;n<VERTICES_PER_FACE;n++)
        mris->vertices[f->v[n]].border = TRUE;
    }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      snorm[0]=snorm[1]=snorm[2]=0;
      v->area = 0;
      for (n=0;n<v->num;n++)
        if (!mris->faces[v->f[n]].ripflag)
        {
          normal_face(v->f[n],v->n[n],norm);
          snorm[0] += norm[0];
          snorm[1] += norm[1];
          snorm[2] += norm[2];
          v->area += triangle_area(v->f[n],v->n[n]);
          /* Note: overest. area by 2! */
        }
      normalize(snorm);

      if (v->origarea<0)
        v->origarea = v->area;

      v->nx = snorm[0];
      v->ny = snorm[1];
      v->nz = snorm[2];
    }
}

void
flip_normals (char *axes)
{
  int flip_x, flip_y, flip_z;
  int    vno ;
  VERTEX *v ;

  printf ("surfer: flipping normals for ");
  flip_x = flip_y = flip_z = 0;
  if (NULL != strstr (axes, "x"))
  {
    flip_x = 1;
    printf ("x ");
  }
  if (NULL != strstr (axes, "y"))
  {
    flip_y = 1;
    printf ("y ");
  }
  if (NULL != strstr (axes, "z"))
  {
    flip_z = 1;
    printf ("z ");
  }
  if (!flip_x && !flip_y && !flip_z)
  {
    printf ("no axes. Please specify x y or z or a combination, i.e. xz\n");
    return;
  }
  printf ("\n");

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (flip_x) v->nx = -v->nx;
    if (flip_y) v->ny = -v->ny;
    if (flip_z) v->nz = -v->nz;
  }
  vertex_array_dirty = 1;
  color_scale_changed = 1;
}

void
compute_shear(void)
{
#if 0
  int k,n,sumn;
  float x0,y0,x1,y1,x2,y2,x3,y3;
  float x01,y01,x12,y12,x23,y23,x30,y30,x02,y02,x13,y13;
  float d01,d12,d23,d30,d02,d13,d,r,r1,r2,shearx,sheary;
  float sx,sy,dsx,dsy,dshearx,dsheary;
  FACE *f;
  VERTEX *v;

  for (k=0;k<mris->nfaces;k++)
    if (!mris->faces[k].ripflag)
    {
      f = &mris->faces[k];
      x0 = mris->vertices[f->v[0]].x;
      y0 = mris->vertices[f->v[0]].y;
      x1 = mris->vertices[f->v[1]].x;
      y1 = mris->vertices[f->v[1]].y;
      x2 = mris->vertices[f->v[2]].x;
      y2 = mris->vertices[f->v[2]].y;
#if VERTICES_PER_FACE == 4
      x3 = mris->vertices[f->v[3]].x;
      y3 = mris->vertices[f->v[3]].y;
#else
      /* don't know what to do here!! */
      x3 = mris->vertices[f->v[0]].x;
      y3 = mris->vertices[f->v[0]].y;
#endif
      x01 = x1-x0;
      y01 = y1-y0;
      x12 = x2-x1;
      y12 = y2-y1;
      x23 = x3-x2;
      y23 = y3-y2;
      x30 = x0-x3;
      y30 = y0-y3;
      x02 = x2-x0;
      y02 = y2-y0;
      x13 = x3-x1;
      y13 = y3-y1;
      d01 = sqrt(x01*x01+y01*y01);
      d12 = sqrt(x12*x12+y12*y12);
      d23 = sqrt(x23*x23+y23*y23);
      d30 = sqrt(x30*x30+y30*y30);
      d02 = sqrt(x02*x02+y02*y02);
      d13 = sqrt(x13*x13+y13*y13);
      r1 = (((d01+d23)*(d12+d30))>0)?log((d01+d23)/(d12+d30)):0;
      r2 = ((d02*d13)>0)?log(d02/d13):0;
      if (fabs(r1)>fabs(r2))
      {
        r = fabs(r1);
        if (r1>0)
        {
          shearx = (x01-x23-y12+y30);
          sheary = (y01-y23+x12-x30);
        }
        else
        {
          shearx = (x12-x30+y01-y23);
          sheary = (y12-y30-x01+x23);
        }
      }
      else
      {
        r = fabs(r2);
        if (r2>0)
        {
          shearx = (x02-y13);
          sheary = (y02+x13);
        }
        else
        {
          shearx = (x13+y02);
          sheary = (y13-x02);
        }
      }
      d = sqrt(shearx*shearx+sheary*sheary);
      /*
        printf("k=%d: r1=%f, r2=%f, r=%f,
        shear=(%f,%f)\n",k,r1,r2,r,shearx,sheary);
      */
      shearx = (d>0)?r*shearx/d:0;
      sheary = (d>0)?r*sheary/d:0;
      f->logshear = r;
      f->shearx = shearx;
      f->sheary = sheary;
    }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      dshearx = dsheary = sumn = 0;
      for (n=0;n<v->num;n++)
        if (!mris->faces[v->f[n]].ripflag)
        {
          sx = mris->faces[v->f[n]].shearx;
          sy = mris->faces[v->f[n]].sheary;
          double_angle(sx,sy,&dsx,&dsy);
          dshearx += dsx;
          dsheary += dsy;
          sumn++;
        }
      dshearx = (sumn>0)?dshearx/sumn:0;
      dsheary = (sumn>0)?dsheary/sumn:0;
      halve_angle(dshearx,dsheary,&shearx,&sheary);
      mris->vertices[k].shearx = shearx;
      mris->vertices[k].sheary = sheary;
      r = sqrt(shearx*shearx+sheary*sheary);
      mris->vertices[k].logshear = r;
      /*
        printf("k=%d: logshear=%f shear=(%f,%f)\n",k,r,shearx,sheary);
      */
    }
#endif
}
void
double_angle(float x,float y,float *dx,float *dy)
{
  float a,d;

  d = sqrt(x*x+y*y);
  a = atan2(y,x);
  *dx = d*cos(2*a);
  *dy = d*sin(2*a);
  /*
    printf("double_angle(%5.2f,%5.2f,%5.2f,%5.2f):d=%5.2f,a=%6.1f,2*a=%6.1f\n",
    x,y,*dx,*dy,d,a*180/M_PI,2*a*180/M_PI);
  */
}

void
halve_angle(float x,float y,float *dx,float *dy)
{
  float a,d;

  d = sqrt(x*x+y*y);
  a = atan2(y,x);
  *dx = d*cos(a/2);
  *dy = d*sin(a/2);
  /*
    printf("halve_angle(%5.2f,%5.2f,%5.2f,%5.2f):d=%5.2f,a=%6.1f,a/2=%6.1f\n",
    x,y,*dx,*dy,d,a*180/M_PI,a/2*180/M_PI);
  */
}

void
compute_boundary_normals(void)
{
#if 0
  int k,m,n;
  VERTEX *v;
  float sumx,sumy,r,nx,ny,f;

  for (k=0;k<mris->nvertices;k++)
    if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
    {
      v = &mris->vertices[k];
      n = 0;
      sumx = 0;
      sumy = 0;
      for (m=0;m<v->vnum;m++)
        if (!mris->vertices[v->v[m]].ripflag)
        {
          sumx += v->x-mris->vertices[v->v[m]].x;
          sumy += v->y-mris->vertices[v->v[m]].y;
          n++;
        }
      v->bnx = (n>0)?sumx/n:0;
      v->bny = (n>0)?sumy/n:0;
    }
  for (k=0;k<mris->nvertices;k++)
    if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
    {
      v = &mris->vertices[k];
      n = 0;
      sumx = 0;
      sumy = 0;
      for (m=0;m<v->vnum;m++)
        if ((!mris->vertices[v->v[m]].ripflag)&&
            mris->vertices[v->v[m]].border)
        {
          nx = -(v->y-mris->vertices[v->v[m]].y);
          ny = v->x-mris->vertices[v->v[m]].x;
          f = nx*v->bnx+ny*v->bny;
          /*
            f = (f<0)?-1.0:(f>0)?1.0:0.0;
          */
          sumx += f*nx;
          sumy += f*ny;
          n++;
        }
      v->bnx = (n>0)?sumx/n:0;
      v->bny = (n>0)?sumy/n:0;
    }
  for (k=0;k<mris->nvertices;k++)
    if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
    {
      r = sqrt(SQR(mris->vertices[k].bnx)+SQR(mris->vertices[k].bny));
      if (r>0)
      {
        mris->vertices[k].bnx /= r;
        mris->vertices[k].bny /= r;
      }
    }
#endif
}

void
subsample_dist(int spacing)
{
  int k ;

  sub_num = MRISsubsampleDist(mris, spacing) ;
  MRISclearMarks(mris) ;
  for (k = 0 ; k < mris->nvertices ; k++)
    if (mris->vertices[k].d == 0)
      mris->vertices[k].marked = 1 ;

  printf("surfer: surface distance subsampled to %d vertices\n",sub_num);
  PR
}

/* 1-adjdot */
void
subsample_orient(float spacing)
{
  int k,m,n;
  VERTEX *v,*v2,*v3;
  float d;

  sub_num = 0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->d = 10000;
    v->val = 0;
  }
  for (k=0;k<mris->nvertices;k++)
    if (mris->vertices[k].d > 0)
    {
      v = &mris->vertices[k];
      for (m=0;m<v->vnum;m++)
      {
        v2 = &mris->vertices[v->v[m]];
        if (v2->d == 0)
        {
          d = 1.0 - (v->nx*v2->nx+v->ny*v2->ny+v->nz*v2->nz);
          if (v->d > d) v->d = d;
        }
        for (n=0;n<v2->vnum;n++)
        {
          v3 = &mris->vertices[v2->v[n]];
          if (v3->d == 0)
          {
            d = 1.0 - (v->nx*v3->nx+v->ny*v3->ny+v->nz*v3->nz);
            if (v->d > d) v->d = d;
          }
        }
      }
      if (v->d>=spacing)
      {
        v->d = 0;
        v->val = 1;
        v->fixedval = TRUE;
        sub_num++;
      }
    }
  for (k=0;k<mris->nvertices;k++)
    if (mris->vertices[k].d == 0)
    {
      v = &mris->vertices[k];
      n = 0;
      for (m=0;m<v->vnum;m++)
      {
        if (mris->vertices[v->v[m]].d > 0)
          n++;
      }
      if (n == 0)
      {
        v->d = 1;
        v->val = 0;
        sub_num--;
      }
    }

  printf("surfer: surface orientation subsampled to %d vertices\n",sub_num);
  PR
}

void
smooth_curv(int niter)
{
  int iter,k,m,n;
  VERTEX *v;
  float sum;

  for (iter=0;iter<niter;iter++)
  {
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      sum = v->curv;
      n = 1;
      for (m=0;m<v->vnum;m++)
        if (!mris->vertices[v->v[m]].ripflag)
        {
          if (mris->vertices[v->v[m]].d <= v->d)
          {
            sum += mris->vertices[v->v[m]].curv;
            n++;
          }
        }
      if (n>0) v->curv = sum/n;
    }
  }
}

#if 1
void
smooth_val_sparse(int niter)
{
  int    iter,k,m,n;
  VERTEX *v;
  float  sum;

  printf("surfer: smooth_val_sparse(%d)\n",niter);
  for (iter=0;iter<niter;iter++)
  {
    printf(".");
    fflush(stdout);
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        mris->vertices[k].tdx = mris->vertices[k].val;
        mris->vertices[k].old_undefval = mris->vertices[k].undefval;
      }
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
        if (!mris->vertices[k].fixedval)
        {
          v = &mris->vertices[k];
          if (k == 3)
            DiagBreak() ;
          sum=0;
          n = 0;
          if (!v->old_undefval)
          {
            sum = v->val;
            n = 1;
          }
          for (m=0;m<v->vnum;m++)
          {
            /*        if (mris->vertices[v->v[m]].dist < v->dist) */
            /*        if (mris->vertices[v->v[m]].dist <= v->dist) */
            if (!mris->vertices[v->v[m]].old_undefval)
            {
              sum += mris->vertices[v->v[m]].tdx;
              n++;
            }
          }
          if (n>0)
          {
            v->val = sum/n;
            v->undefval = FALSE;
          }
        }
  }
  printf("\n");
  PR
}
#else
void
smooth_val_sparse(int niter)
{
#if 0
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->marked = v->fixedval ;
  }
  MRISsoapBubbleVals(mris, niter) ;
  MRISclearMarks(mris) ;

#else
int iter,k,m,n;
VERTEX *v;
float sum;

printf("surfer: smooth_val_sparse(%d)\n",niter);
for (iter=0;iter<niter;iter++)
{
  printf(".");
  fflush(stdout);
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      mris->vertices[k].tdx = mris->vertices[k].val;
      mris->vertices[k].old_undefval = mris->vertices[k].undefval;
    }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
      if (!mris->vertices[k].fixedval)
      {
        v = &mris->vertices[k];
        sum=0;
        n = 0;
        if (!v->old_undefval)
        {
          sum = v->val;
          n = 1;
        }
        for (m=0;m<v->vnum;m++)
        {
          /*        if (mris->vertices[v->v[m]].d < v->d) */
          /*        if (mris->vertices[v->v[m]].d <= v->d) */
          if (!mris->vertices[v->v[m]].old_undefval)
          {
            sum += mris->vertices[v->v[m]].tdx;
            n++;
          }
        }
        if (n>0)
        {
          v->val = sum/n;
          v->undefval = FALSE;
        }
      }
}
printf("\n");
PR
#endif
}
#endif

void
smooth_val(int niter)
{
  int iter,k,m,n;
  VERTEX *v;
  float sum;

  surface_compiled = 0 ;
  printf("surfer: smooth_val(%d)\n",niter);
  for (iter=0;iter<niter;iter++)
  {
    printf(".");
    fflush(stdout);
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].tdx = mris->vertices[k].val;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      sum=v->tdx;
      if (k == Gdiag_no)
        DiagBreak() ;
      n = 1;
      for (m=0;m<v->vnum;m++)
      {
        sum += mris->vertices[v->v[m]].tdx;
        n++;
      }
      if (!finite(sum))
        DiagBreak() ;
      if (n>0)
        v->val = sum/n;
      if (!finite(v->val))
        DiagBreak() ;
    }
  }
  printf("\n");
  PR
}

void
smooth_fs(int niter)
{
  int iter,k,m,n;
  VERTEX *v;
  float sum;

  printf("surfer: smooth_fs(%d)\n",niter);
  for (iter=0;iter<niter;iter++)
  {
    printf(".");
    fflush(stdout);
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].tdx = mris->vertices[k].fieldsign;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      sum=v->tdx;
      n = 1;
      for (m=0;m<v->vnum;m++)
      {
        sum += mris->vertices[v->v[m]].tdx;
        n++;
      }
      if (n>0) v->fieldsign = sum/n;
    }
  }
  printf("\n");
  PR
}

/* begin rkt */

int
sclv_smooth(int niter, int field)
{
  int iter,k,m,n;
  VERTEX *v;
  float sum, average;

  surface_compiled = 0 ;
  printf("surfer: sclv_smooth(%d,%s)\n",niter,sclv_field_names[field]);
  for (iter=0;iter<niter;iter++)
  {
    printf(".");
    fflush(stdout);
    for (k=0;k<mris->nvertices;k++)
      sclv_get_value( (&(mris->vertices[k])),
                      field, &(mris->vertices[k].tdx));
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      sum=v->tdx;
      if (k == Gdiag_no)
        DiagBreak() ;
      n = 1;
      for (m=0;m<v->vnum;m++)
      {
        sum += mris->vertices[v->v[m]].tdx;
        n++;
      }
      average = 0;
      if ( n != 0 )
        average = sum / (float)n;
      if (!finite(sum))
        DiagBreak() ;
      if (!finite(average))
        DiagBreak() ;
      sclv_set_value (v, field, average);
    }
  }
  printf("\n");
  PR;

  /* values have changed, need to recalc frequencies */
  sclv_calc_frequencies (field);

  return (ERROR_NONE);
}

/* end rkt */

void
fatten_border(int niter)
{
  int iter,k,m,n;
  VERTEX *v;

  printf("surfer: fatten_border(%d)\n",niter);
  n = 0;
  for (iter=0;iter<niter;iter++)
  {
    for (k=0;k<mris->nvertices;k++)
      mris->vertices[k].d = mris->vertices[k].border;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      for (m=0;m<v->vnum;m++)
        if (mris->vertices[v->v[m]].d!=0)
        {
          v->border = TRUE;
          n++;
        }
    }
  }
  printf("surfer: %d border points added\n",n);
  PR
}

void
compute_angles(void)
{
  int k;
  float val,valbak,val2,val2bak;

  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      val = atan2(mris->vertices[k].val2,mris->vertices[k].val);
      val2 = sqrt(SQR(mris->vertices[k].val2)+SQR(mris->vertices[k].val));
      valbak = atan2(mris->vertices[k].val2bak,mris->vertices[k].valbak);
      val2bak = sqrt(SQR(mris->vertices[k].val2bak)+
                     SQR(mris->vertices[k].valbak));
      mris->vertices[k].val = val;
      mris->vertices[k].val2 = val2;
      mris->vertices[k].valbak = valbak;
      mris->vertices[k].val2bak = val2bak;
    }
}

float circsubtract(float a,float b)
{
  float h = a-b;
  if (h<-M_PI) h = h+2*M_PI;
  else if (h>M_PI) h = h-2*M_PI;
  return h;
}

void
compute_fieldsign(void)
{
  int k,m,n;
  VERTEX *v;
  float dv1,dv2,dx,dy,dv1dx,dv1dy,dv2dx,dv2dy;
  float m11,m12,m13,m22,m23,z1,z2,z3,z1b,z2b,z3b,denom;

  MRISremoveTriangleLinks(mris) ;
  printf("surfer: compute_fieldsign()\n");
  PR for (k=0;k<mris->nvertices;k++)
  {
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      dv1dx = dv1dy = dv2dx = dv2dy = 0;
      m11 = m12 = m13 = m22 = m23 = z1 = z2 = z3 = z1b = z2b = z3b = 0;
      n = 0;
      for (m=0;m<v->vnum;m++)
        if (!mris->vertices[v->v[m]].ripflag)
        {
          dv1 = circsubtract(v->val,mris->vertices[v->v[m]].val);
          dv2 = circsubtract(v->valbak,mris->vertices[v->v[m]].valbak);
          dx = v->x-mris->vertices[v->v[m]].x;
          dy = v->y-mris->vertices[v->v[m]].y;
          m11 += dx*dx;
          m12 += dx*dy;
          m13 += dx;
          m22 += dy*dy;
          m23 += dy;
          z1 += dx*dv1;
          z2 += dy*dv1;
          z3 += dv1;
          z1b += dx*dv2;
          z2b += dy*dv2;
          z3b += dv2;
        }
      dv1dx = (m22*z1-m23*m23*z1-m12*z2+m13*m23*z2-m13*m22*z3+m12*m23*z3);
      dv2dx = (m22*z1b-m23*m23*z1b-
               m12*z2b+m13*m23*z2b-
               m13*m22*z3b+m12*m23*z3b);
      dv1dy = (-m12*z1+m13*m23*z1+m11*z2-m13*m13*z2+m12*m13*z3-m11*m23*z3);
      dv2dy = (-m12*z1b+m13*m23*z1b+
               m11*z2b-m13*m13*z2b+
               m12*m13*z3b-m11*m23*z3b);
      denom = -m12*m12+m11*m22-m13*m13*m22+2*m12*m13*m23-m11*m23*m23;
      if (denom!=0)
        v->fieldsign = (dv1dx*dv2dy-dv2dx*dv1dy)/(denom*denom);
      else
        v->fieldsign = 0;

      v->fieldsign =  ((v->fieldsign<0)?-1:(v->fieldsign>0)?1:0);
      if (revfsflag)
        v->fieldsign = -(v->fieldsign);

      v->fsmask = sqrt(v->val2*v->val2bak);  /* geom mean of r,th power */
    }
    fieldsignflag = TRUE;
  }
}

void
compute_cortical_thickness(void)   /* marty */
{
  VERTEX *v;
  int imnr,i,j;
  int h,k,m,n;
  int outi,outj,outim;
  int gmval;
  float sum;
  float gmvalhist[GRAYBINS];
  float gmdist;
  float x,y,z;

  if (!MRIflag || !MRIloaded)
  {
    printf("surfer: ### need MRI data to compute cortical thickness\n");
    PR return;
  }

  compute_normals();

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->d = 1;   /* for smooth_val_sparse */
    v->val = 0.0;
    x = v->x;
    y = v->y;
    z = v->z;
    for (h=0;h<GRAYBINS;h++)
    {
      gmdist = (float)h*GRAYINCR;

      imnr = (int)((y-yy0)/st+0.5-imnr0);
      i = (int)((zz1-z)/ps+0.5);
      j = (int)((xx1-x)/ps+0.5);
      imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
      i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
      j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
      outim = (int)(imnr+gmdist/st*v->ny+0.5);
      outi = (int)(i-gmdist/ps*v->nz+0.5);
      outj = (int)(j-gmdist/ps*v->nx+0.5);
      if (outim>=0&&outim<numimg) gmval = im[outim][outi][outj];
      else                        gmval = 0;

      sum = gmval;
      n = 1;

      for (m=0;m<v->vnum;m++)
      {
        x = mris->vertices[v->v[m]].x;
        y = mris->vertices[v->v[m]].y;
        z = mris->vertices[v->v[m]].z;
        imnr = (int)((y-yy0)/st+0.5-imnr0);
        i = (int)((zz1-z)/ps+0.5);
        j = (int)((xx1-x)/ps+0.5);
        imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
        i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
        j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
        outim = (int)(imnr+gmdist/st*v->ny+0.5);
        outi = (int)(i-gmdist/ps*v->nz+0.5);
        outj = (int)(j-gmdist/ps*v->nx+0.5);
        if (outim>=0&&outim<numimg) gmval = im[outim][outi][outj];
        else                        gmval = 0;

        sum += gmval;
        n++;
      }
      if (n>0) gmvalhist[h] = sum/n;
      else     gmvalhist[h] = 0.0;
      /*** printf("avggmval[%d]=%f ",h,gmvalhist[h]); */
    }
    for (h=0;h<GRAYBINS;h++)
    {
      gmdist = (float)h*GRAYINCR;
      if (gmvalhist[h]<=mingm)
      {
        gmdist -= 0.5 * (GRAYBINS-1)*GRAYINCR;
        /* green thin, red thick */
        v->val = gmdist;  /* no gm dropoff: val=0.0 -  */
        break;
      }
    }
    /*** printf("gmdist=%f\n\n",gmdist); */
  }
}

void
compute_CMF(void)  /* TODO: update to use new labels! */
{
  int i,j,k,m,n,label,bin;
  VERTEX *v;
  float dv1,dv2,dx,dy,dv1dx,dv1dy,dv2dx,dv2dy;
  float m11,m12,m13,m22,m23,z1,z2,z3,z1b,z2b,z3b,denom;
  float sum,a,val1,val2,rad,radmax = 0.0f;
  char fname[200];
  FILE *fp;
  float maxrad=10.0,minrad=0.05*maxrad,afact;

  afact = pow(maxrad/minrad,1.0/(maxrad-minrad));
  for (i=0;i<NLABELS;i++)
  {
    valtot[i] = 0;
    radmin[i] = 1000000;
    gradx_avg[i] = 0;
    grady_avg[i] = 0;
    for (j=0;j<CMFBINS;j++)
    {
      val1avg[i][j] = 0;
      val2avg[i][j] = 0;
      valnum[i][j] = 0;
    }
  }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag && !mris->vertices[k].border)
    {
      v = &mris->vertices[k];
      a = v->val/(2*M_PI);
      a += 0.5;
      a += angle_offset;
      a = fmod(a,1.0);
      v->val = a*2*M_PI;
    }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag && !mris->vertices[k].border)
      if (mris->vertices[k].annotation!=0)
      {
        v = &mris->vertices[k];
        dv1dx = dv1dy = dv2dx = dv2dy = 0;
        m11 = m12 = m13 = m22 = m23 = z1 = z2 = z3 = z1b = z2b = z3b = 0;
        n = 0;
        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].ripflag &&
              !mris->vertices[v->v[m]].border)
          {
            dv1 = circsubtract(v->val,mris->vertices[v->v[m]].val);
            dv2 = circsubtract(v->valbak,mris->vertices[v->v[m]].valbak);
            dx = v->x-mris->vertices[v->v[m]].x;
            dy = v->y-mris->vertices[v->v[m]].y;
            m11 += dx*dx;
            m12 += dx*dy;
            m13 += dx;
            m22 += dy*dy;
            m23 += dy;
            z1 += dx*dv1;
            z2 += dy*dv1;
            z3 += dv1;
            z1b += dx*dv2;
            z2b += dy*dv2;
            z3b += dv2;
          }
        dv1dx = (m22*z1-m23*m23*z1-m12*z2+m13*m23*z2-m13*m22*z3+m12*m23*z3);
        dv2dx =
          (m22*z1b-m23*m23*z1b-m12*z2b+m13*m23*z2b-m13*m22*z3b+m12*m23*z3b);
        dv1dy = (-m12*z1+m13*m23*z1+m11*z2-m13*m13*z2+m12*m13*z3-m11*m23*z3);
        dv2dy =
          (-m12*z1b+m13*m23*z1b+m11*z2b-m13*m13*z2b+m12*m13*z3b-m11*m23*z3b);
        denom = -m12*m12+m11*m22-m13*m13*m22+2*m12*m13*m23-m11*m23*m23;
        if (denom!=0)
        {
          dv1dx /= denom;
          dv2dx /= denom;
          dv1dy /= denom;
          dv2dy /= denom;
        }
        label = mris->vertices[k].annotation;
        gradx_avg[label] += dv1dx;
        grady_avg[label] += dv1dy;
      }
  for (i=0;i<NLABELS;i++)
  {
    a = sqrt(SQR(gradx_avg[i])+SQR(grady_avg[i]));
    if (a!=0)
    {
      gradx_avg[i] /= a;
      grady_avg[i] /= a;
    }
  }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag && !mris->vertices[k].border)
      if (mris->vertices[k].annotation!=0)
      {
        v = &mris->vertices[k];
        label = mris->vertices[k].annotation;
        rad = (v->x*gradx_avg[label]+v->y*grady_avg[label]);
        if (rad<radmin[label]) radmin[label] = rad;
      }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag && !mris->vertices[k].border)
      if (mris->vertices[k].annotation!=0)
      {
        v = &mris->vertices[k];
        label = mris->vertices[k].annotation;
        rad = (v->x*gradx_avg[label]+v->y*grady_avg[label])-radmin[label];
        bin = floor(rad/MAXVAL*CMFBINS);
        if (bin<0)
        {
          printf("surfer: --- bin=%d too low\n",bin);
          bin=0;
        }
        if (bin>CMFBINS-1)
        {
          printf("surfer: --- bin=%d too high\n",bin);
          bin=CMFBINS-1;
        }
        val1 = cos(v->val)*v->val2;
        val2 = sin(v->val)*v->val2;
        val1avg[label][bin] += val1;
        val2avg[label][bin] += val2;
        valnum[label][bin] += 1;
        valtot[label] += 1;
      }
  for (i=0;i<NLABELS;i++)
    if (valtot[i]!=0)
    {
      printf("surfer: label=%d\n",i);
      sum = 0;
      for (j=0;j<CMFBINS;j++)
      {
        if (valnum[i][j]!=0)
        {
          a = atan2(val2avg[i][j],val1avg[i][j]);
          if (a<0) a += 2*M_PI;
          a = 2*M_PI-a;
          a = a*maxrad/(2*M_PI);
          a = minrad*pow(afact,a-minrad);
          angleavg[i][j] = a;
        }
      }
      sprintf(fname,"cmftmp-%d.xplot",i);
      fp = fopen(fname,"w");
      for (j=0;j<CMFBINS;j++)
        if (valnum[i][j]!=0)
          radmax = j*MAXVAL/CMFBINS;
      for (j=0;j<CMFBINS;j++)
        if (valnum[i][j]!=0)
        {
          fprintf(fp,"%f %f\n",radmax-j*MAXVAL/CMFBINS,angleavg[i][j]);
        }
      fclose(fp);
    }
  PR
}

void
smooth_momentum(int niter)
{
#if 0
  int iter,k,m,n;
  VERTEX *v;
  float sumx,sumy,sumz;

  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      v->smx = v->odx;
      v->smy = v->ody;
      v->smz = v->odz;
    }
  for (iter=0;iter<niter;iter++)
  {
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        v = &mris->vertices[k];
        v->osmx = v->smx;
        v->osmy = v->smy;
        v->osmz = v->smz;
      }
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        v = &mris->vertices[k];
        sumx = v->osmx;
        sumy = v->osmy;
        sumz = v->osmz;
        n = 1;
        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].ripflag)
          {
            sumx += mris->vertices[v->v[m]].osmx;
            sumy += mris->vertices[v->v[m]].osmy;
            sumz += mris->vertices[v->v[m]].osmz;
            n++;
          }
        if (n>0) v->smx = sumx/n;
        if (n>0) v->smy = sumy/n;
        if (n>0) v->smz = sumz/n;
      }
  }
#endif
}

void
smooth_logarat(int niter)
{
#if 0
  int iter,k,m,n;
  VERTEX *v;
  float sum,arat;

  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      arat = mris->vertices[k].area/mris->vertices[k].origarea;
      mris->vertices[k].logarat = (arat>0)?log(arat):0;
      /*
        printf("%d: area=%f, origarea=%f, arat=%f, logarat=%f\n",
        k,mris->vertices[k].area,mris->vertices[k].origarea,arat,
        mris->vertices[k].logarat);
      */
    }
  for (iter=0;iter<niter;iter++)
  {
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
        mris->vertices[k].ologarat = mris->vertices[k].logarat;
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        v = &mris->vertices[k];
        sum = v->ologarat;
        n = 1;
        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].ripflag)
          {
            sum += mris->vertices[v->v[m]].ologarat;
            n++;
          }
        if (n>0) v->logarat = sum/n;
      }
  }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      mris->vertices[k].sqrtarat = exp(0.5*mris->vertices[k].logarat);
    }
#endif
}

void
smooth_shear(int niter)
{
#if 0
  int iter,k,m,n;
  VERTEX *v;
  float sumx,sumy,r;

  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      double_angle(mris->vertices[k].shearx,mris->vertices[k].sheary,
                   &mris->vertices[k].shearx,&mris->vertices[k].sheary);
    }
  for (iter=0;iter<niter;iter++)
  {
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        mris->vertices[k].oshearx = mris->vertices[k].shearx;
        mris->vertices[k].osheary = mris->vertices[k].sheary;
      }
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        v = &mris->vertices[k];
        sumx = v->oshearx;
        sumy = v->osheary;
        n = 1;
        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].ripflag)
          {
            sumx += mris->vertices[v->v[m]].oshearx;
            sumy += mris->vertices[v->v[m]].osheary;
            n++;
          }
        v->shearx = (n>0)?sumx/n:0;
        v->sheary = (n>0)?sumy/n:0;
      }
  }
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      halve_angle(mris->vertices[k].shearx,mris->vertices[k].sheary,
                  &mris->vertices[k].shearx,&mris->vertices[k].sheary);
      r = sqrt(SQR(mris->vertices[k].shearx)+SQR(mris->vertices[k].sheary));
      mris->vertices[k].logshear = r;
    }
#endif
}

void
smooth_boundary_normals(int niter)
{
#if 0
  int iter,k,m,n;
  VERTEX *v;
  float sumx,sumy,r;

  for (iter=0;iter<niter;iter++)
  {
    for (k=0;k<mris->nvertices;k++)
      if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
      {
        mris->vertices[k].obnx = mris->vertices[k].bnx;
        mris->vertices[k].obny = mris->vertices[k].bny;
      }
    for (k=0;k<mris->nvertices;k++)
      if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
      {
        v = &mris->vertices[k];
        n = 1;
        sumx = v->obnx;
        sumy = v->obny;
        for (m=0;m<v->vnum;m++)
          if ((!mris->vertices[v->v[m]].ripflag)&&
              mris->vertices[v->v[m]].border)
          {
            sumx += mris->vertices[v->v[m]].obnx;
            sumy += mris->vertices[v->v[m]].obny;
            n++;
          }
        v->bnx = (n>0)?sumx/n:0;
        v->bny = (n>0)?sumy/n:0;
      }
  }
  for (k=0;k<mris->nvertices;k++)
    if ((!mris->vertices[k].ripflag)&&mris->vertices[k].border)
    {
      r = sqrt(SQR(mris->vertices[k].bnx)+SQR(mris->vertices[k].bny));
      if (r>0)
      {
        mris->vertices[k].bnx /= r;
        mris->vertices[k].bny /= r;
      }
    }
#endif
}

void
scaledist(float sf)
{
  int k;
  VERTEX *v;

  printf("surfer: scaledist(%f)\n",sf);
  PR for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->x = (v->x+mris->xctr)*sf-mris->xctr;
    v->y = (v->y-mris->yctr)*sf+mris->yctr;
    v->z = (v->z-mris->zctr)*sf+mris->zctr;
  }
}

float rtanh(float x)
{
  return (x<0.0)?0.0:tanh(x);
}
static void
Laplacian(MRI_SURFACE *mris, int vno, double *plx, double *ply, double *plz)
{
  VERTEX *v, *vn ;
  double lx, ly, lz, w ;
  int    n ;

  v = &mris->vertices[vno] ;
  if (v->vnum == 0)
    return ;
  w = 1.0 / (double)v->vnum ;
  for (lx = ly = lz = 0.0, n = 0 ; n < v->vnum ; n++)
  {
    vn = &mris->vertices[v->v[n]] ;
    lx += w * (v->x - vn->x) ;
    ly += w * (v->y - vn->y) ;
    lz += w * (v->z - vn->z) ;
  }
  *plx = lx ; *ply = ly ; *plz = lz ;
  return ;
}
#define K_bp   0.1
#define Lambda .63
#define Mu     -.67236

void
taubin_smooth(MRI_SURFACE *mris, int niter, double mu, double lambda)
{
  int    iter, vno ;
  double dx, dy, dz ;
  VERTEX *v ;

  dx = dy = dz = 0 ;
  for (iter = 0 ; iter < niter ; iter++)
  {
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      Laplacian(mris, vno, &dx, &dy, &dz) ;
      if (EVEN(iter))
      {
        v->tx = v->x + lambda * dx ;
        v->ty = v->y + lambda * dy ;
        v->tz = v->z + lambda * dz ;
      }
      else
      {
        v->tx = v->x + mu * dx ;
        v->ty = v->y + mu * dy ;
        v->tz = v->z + mu * dz ;
      }
    }
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  }


  MRIScomputeMetricProperties(mris) ;
  vset_save_surface_vertices(vset_current_set) ;
  vset_set_current_set(vset_current_set) ;
  vertex_array_dirty = 1 ;
  redraw() ;
}


void
shrink(int niter)   /* shrinktypes: nei,mri,curv,exp,area,2d,fill,sphe,ell */
{
  float x,y,z,sx,sy,sz,val,inval,outval,nc,force;
  float d,dx,dy,dz,sval,sinval,soutval,snc,anc, orig_area, scale ;
  float nx,ny,nz;
  VERTEX *v;
  int imnr,i,j,iter,k,m,n;
  float stress,sd,ad,dmax;
  int navg,an,nclip,inim,ini,inj,outim,outi,outj;

  orig_area = mris->total_area ;
  val = inval = outval = 0.0f ;
  if (shrinkfillflag)
  {
    shrink_to_fill(niter);
    return;
  }
  if (flag2d)
  {
    shrink2d(niter);
    return;
  }
  if (areaflag)
  {
    area_shrink(niter);
    return;
  }
  if (MRIflag && !MRIloaded)
    printf("surfer: MRIflag on but no MRI data loaded... MRIflag ignored\n");

  for (iter=0;iter<niter;iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
      /*      v->onc = v->nc;*/
      v->oripflag = v->ripflag;
    }
    sval = sinval = soutval = snc = 0;
    maxstress = avgstress = 0;
    navg = 0;
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].oripflag)
      {
        v = &mris->vertices[k];
        x = v->tx;
        y = v->ty;
        z = v->tz;
        nx = v->nx;
        ny = v->ny;
        nz = v->nz;
        sx=sy=sz=sd=0;
        n=0;
        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].oripflag)
          {
            sx += dx = mris->vertices[v->v[m]].tx - x;
            sy += dy = mris->vertices[v->v[m]].ty - y;
            sz += dz = mris->vertices[v->v[m]].tz - z;
            sd += sqrt(dx*dx+dy*dy+dz*dz);
            n++;
          }
        if (n>0)
        {
          sx = sx/n;
          sy = sy/n;
          sz = sz/n;
          sd = sd/n;
          stress = sd;
          if (explodeflag && stress>=stressthresh)
          {
            nrip++;
            v->ripflag = TRUE;
          }
          if (stress>maxstress)
            maxstress = stress;
          avgstress += stress;
          navg++;
#if 0
          v->stress = stress;
#endif
        }
        if (MRIflag && MRIloaded)
        {
          imnr = (int)((y-yy0)/st+0.5);
          i = (int)((zz1-z)/ps+0.5);
          j = (int)((xx1-x)/ps+0.5);
          imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
          i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
          j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
          val = im[imnr][i][j];
          inim = (int)(imnr-cthk/st*v->ny+0.5);
          ini = (int)(i+cthk/ps*v->nz+0.5);
          inj = (int)(j+cthk/ps*v->nx+0.5);
          if (inim>=0&&inim<numimg)
            inval = im[inim][ini][inj];
          else
            inval = 0;
          outim = (int)(imnr+ostilt/st*v->ny+0.5);
          outi = (int)(i-ostilt/ps*v->nz+0.5);
          outj = (int)(j-ostilt/ps*v->nx+0.5);
          if (outim>=0&&outim<numimg)
            outval = im[outim][outi][outj];
          else
            outval = 0;
          /*      force = mstrength*tanh((val-mmid)*mslope); */
          /*      force = mstrength*tanh((inval-mmid)*mslope); */
          /*
            force = mstrength*(tanh((inval-mmid)*mslope)
            -rtanh((outval-mmid)*mslope));
          */

          force = (mstrength*tanh((inval-whitemid)*mslope)+
                   mstrength*tanh((val-graymid)*mslope))/2;

        }
        else
          if (expandflag)
          {
            force = 0.1;
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
#if 0
        if (stiffnessflag)
        {
          anc = 0;
          for (m=0;m<v->vnum;m++)
            anc += mris->vertices[v->v[m]].onc;
          anc /= v->vnum;
          force += tanh((nc-anc)*0.2);
        }
        else
        {
          anc = 0;
          force += tanh((nc-anc)*0.5);
        }
#else
        anc = 0;
        force += tanh((nc-anc)*0.5);
#endif
        dx = sx*0.5 + v->nx*force;
        dy = sy*0.5 + v->ny*force;
        dz = sz*0.5 + v->nz*force;
        if (momentumflag)
        {
          dx = decay*v->odx+update*dx;
          dy = decay*v->ody+update*dy;
          dz = decay*v->odz+update*dz;
          if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
          {
            nclip ++;
            dx /= d;
            dy /= d;
            dz /= d;
          }
          v->odx = dx;
          v->ody = dy;
          v->odz = dz;
        }
        if (sulcflag)
          v->curv += dx*v->nx + dy*v->ny + dz*v->nz;
        d=sqrt(dx*dx+dy*dy+dz*dz);
        if (d>dmax) dmax = d;
        ad += d;
        an ++;
        v->x += dx;
        v->y += dy;
        v->z += dz;
        sval += val;
        sinval += inval;
        soutval += outval;
      }
    snc /= mris->nvertices;
    if (MRIflag && MRIloaded)
    {
      compute_normals();
      sval /= mris->nvertices;
      sinval /= mris->nvertices;
      soutval /= mris->nvertices;
      printf("surfer: %d: sval=%5.2f,sinval=%5.2f,soutval=%5.2f,"
             "snc=%5.2f\n",
             iter,sval,sinval,soutval,snc);
    }
    else
    {
      printf("surfer: %d: ad=%f, dmax=%f, snc=%f\n",iter,ad/an,dmax,snc);
      if (expandflag || explodeflag)
        compute_normals();
    }
    if (navg>0)
      avgstress /= navg;
    if (explodeflag)
      printf("surfer: %d: max = %f, avg = %f, threshold = %f, nrip = %d\n",
             iter,maxstress,avgstress,stressthresh,nrip);
    PR
  }
  if (!(MRIflag || expandflag || explodeflag))
  {
    compute_normals();
  }

  MRIScomputeMetricProperties(mris) ;
  scale = sqrt(orig_area / (mris->total_area+mris->neg_area)) ;
  MRISscaleBrain(mris, mris, scale) ;
  MRIScomputeMetricProperties(mris) ;
  vset_save_surface_vertices(vset_current_set) ;
  vset_set_current_set(vset_current_set) ;
}

void
curv_shrink_to_fill(int niter)
{
  VERTEX *v;
  int imnr,i,j,iter,k,m,n;
  int an,nclip,delcurvdefined;
  float x,y,z,sx,sy,sz,nc,force;
  float d,dx,dy,dz,snc;
  float nx,ny,nz;
  float val,sval;
  float sd,ad,dmax;
  float xnei,ynei,znei;
  float curv,icurv,icurvnei,cfact;
  float icrange,crange;

  icurv = curv = icurvnei = 0.0f ;
  if (!curvloaded)
  {
    printf("surfer: ### curv not loaded!\n");
    PR return;
  }
  if (!curvimloaded)
  {
    printf("surfer: ### curvim not loaded!\n");
    PR return;
  }
  if (!curvimflag)
  {
    printf("surfer: ### curvimflag not set!\n");
    PR return;
  }

  icrange = mris2->max_curv-mris2->min_curv;
  crange = mris->max_curv-mris->min_curv;
  for (iter=0;iter<niter;iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
      /*      v->onc = v->nc;*/
      v->oripflag = v->ripflag;
    }
    sval = snc = 0;
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].oripflag)
      {
        v = &mris->vertices[k];
        x = v->tx;
        y = v->ty;
        z = v->tz;
        nx = v->nx;
        ny = v->ny;
        nz = v->nz;
        sx=sy=sz=sd=0;
        n=0;

        /**** nei curve force */
        delcurvdefined = TRUE;
        imnr = (int)((y-yy0)/st+0.5);
        i = (int)((zz1-z)/ps+0.5);
        j = (int)((xx1-x)/ps+0.5);
        if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
          icurv =
            bytetofloat(curvim[imnr][i][j],mris2->min_curv,mris2->max_curv);
        else
          delcurvdefined = FALSE;
        curv = v->curv;
        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].oripflag)
          {
            xnei = mris->vertices[v->v[m]].tx;
            ynei = mris->vertices[v->v[m]].ty;
            znei = mris->vertices[v->v[m]].tz;
            imnr = (int)((ynei-yy0)/st+0.5);
            i = (int)((zz1-znei)/ps+0.5);
            j = (int)((xx1-xnei)/ps+0.5);
            if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
              icurvnei =
                bytetofloat(curvim[imnr][i][j],
                            mris2->min_curv,mris2->max_curv);
            else
              delcurvdefined = FALSE;
            cfact = 1.0;
            if (delcurvdefined)
              cfact += (icurvnei-icurv)/icrange
                       * copysign(icstrength,
                                  mris2->min_curv+curv*(icrange/crange) -
                                  icurv);
            sx += dx = (xnei - x)*cfact;
            sy += dy = (ynei - y)*cfact;
            sz += dz = (znei - z)*cfact;
            sd += sqrt(dx*dx+dy*dy+dz*dz);
            n++;
          }
        if (n>0)
        {
          sx = sx/n;
          sy = sy/n;
          sz = sz/n;
          sd = sd/n;
        }

        /**** norm shrink_to_fill force */
        imnr = (int)(y+numimg/2.0+0.5);
        i = (int)(IMGSIZE/2.0-z+0.5);
        j = (int)(IMGSIZE/2.0-x+0.5);
        imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
        i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
        j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
        val = (curvimflags[imnr][i][j] & CURVIM_FILLED)?128:0;
        force = mstrength*tanh((val-mmid)*mslope);
        nc = sx*nx+sy*ny+sz*nz;
        snc += nc;
        v->nc = nc;
        force += tanh(nc*0.5);

        dx = sx*0.5 + v->nx*force;
        dy = sy*0.5 + v->ny*force;
        dz = sz*0.5 + v->nz*force;
        if (momentumflag)
        {
          dx = decay*v->odx+update*dx;
          dy = decay*v->ody+update*dy;
          dz = decay*v->odz+update*dz;
          if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
          {
            nclip ++;
            dx /= d;
            dy /= d;
            dz /= d;
          }
          v->odx = dx;
          v->ody = dy;
          v->odz = dz;
        }
        d=sqrt(dx*dx+dy*dy+dz*dz);
        if (d>dmax) dmax = d;
        ad += d;
        an ++;
        v->x += dx;
        v->y += dy;
        v->z += dz;
        sval += val;
      }
    compute_normals();  /* only inside in shrink_to_fill */
    snc /= mris->nvertices;
    printf("surfer: %d: sval=%5.2f,snc=%5.2f\n",iter,sval,snc);
    PR
  }
}

void
shrink_to_fill(int niter)
{
  float x,y,z,sx,sy,sz,val,inval,outval,nc,force;
  float d,dx,dy,dz,sval,sinval,soutval,snc;
  float nx,ny,nz;
  VERTEX *v;
  int imnr,i,j,iter,k,m,n;
  float sd,ad,dmax;
  int navg,an,nclip,inim,ini,inj,outim,outi,outj;

  outval = inval = val = 0.0f ;
  for (iter=0;iter<niter;iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
      /*      v->onc = v->nc;*/
      v->oripflag = v->ripflag;
    }
    sval = sinval = soutval = snc = 0;
    navg = 0;
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].oripflag)
      {
        v = &mris->vertices[k];
        x = v->tx;
        y = v->ty;
        z = v->tz;
        nx = v->nx;
        ny = v->ny;
        nz = v->nz;
        sx=sy=sz=sd=0;
        n=0;
        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].oripflag)
          {
            sx += dx = mris->vertices[v->v[m]].tx - x;
            sy += dy = mris->vertices[v->v[m]].ty - y;
            sz += dz = mris->vertices[v->v[m]].tz - z;
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
        if (MRIflag && MRIloaded)
        {
          imnr = (int)(y+numimg/2.0+0.5);
          i = (int)(IMGSIZE/2.0-z+0.5);
          j = (int)(IMGSIZE/2.0-x+0.5);
          imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
          i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
          j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
          val = im[imnr][i][j];
          inim = (int)(imnr-cthk*v->ny+0.5);
          ini = (int)(i+cthk*v->nz+0.5);
          inj = (int)(j+cthk*v->nx+0.5);
          if (inim>=0&&inim<numimg)
            inval = im[inim][ini][inj];
          else
            inval = 0;
          outim = (int)(imnr+ostilt*v->ny+0.5);
          outi = (int)(i-ostilt*v->nz+0.5);
          outj = (int)(j-ostilt*v->nx+0.5);
          if (outim>=0&&outim<numimg)
            outval = im[outim][outi][outj];
          else
            outval = 0;
          force = mstrength*tanh((val-mmid)*mslope);
        }
        else force = 0;
        nc = sx*nx+sy*ny+sz*nz;
        /*
          sx -= nc*nx;
          sy -= nc*ny;
          sz -= nc*nz;
        */
        snc += nc;
        v->nc = nc;
        /*
          if (stiffnessflag)
          {
          anc = 0;
          for (m=0;m<v->vnum;m++)
          anc += mris->vertices[v->v[m]].onc;
          anc /= v->vnum;
          force += tanh((nc-anc)*0.2);
          }
          else
          {
          anc = 0;
          force += tanh((nc-anc)*0.5);
          }
        */
        dx = sx*0.5 + v->nx*force;
        dy = sy*0.5 + v->ny*force;
        dz = sz*0.5 + v->nz*force;
        if (momentumflag)
        {
          dx = decay*v->odx+update*dx;
          dy = decay*v->ody+update*dy;
          dz = decay*v->odz+update*dz;
          if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
          {
            nclip ++;
            dx /= d;
            dy /= d;
            dz /= d;
          }
          v->odx = dx;
          v->ody = dy;
          v->odz = dz;
        }
        d=sqrt(dx*dx+dy*dy+dz*dz);
        if (d>dmax) dmax = d;
        ad += d;
        an ++;
        v->x += dx;
        v->y += dy;
        v->z += dz;
        sval += val;
        sinval += inval;
        soutval += outval;
      }
    snc /= mris->nvertices;
    compute_normals();
    if (MRIflag && MRIloaded)
    {
      sval /= mris->nvertices;
      sinval /= mris->nvertices;
      soutval /= mris->nvertices;
      printf("surfer: %d: sval=%5.2f,sinval=%5.2f,"
             "soutval=%5.2f,snc=%5.2f\n",
             iter,sval,sinval,soutval,snc);
    }
    else
    {
      printf("surfer: %d: ad=%f, dmax=%f, snc=%f\n",iter,ad/an,dmax,snc);
    }
    PR
  }
}

void
transform(float *xptr, float *yptr, float *zptr, float nx, float ny, float nz,
          float d)  /* 2 vects ortho to summed normal */
{
  float x = *xptr, y = *yptr, z = *zptr;

  *zptr = nx*x + ny*y + nz*z;
  *yptr = -ny/d*x + nx/d*y;
  *xptr = nx*nz/d*x + ny*nz/d*y - d*z;
  /*
    printf("transform {%f,%f,%f} -> {%f,%f,%f}\n",
    x,y,z,*xptr,*yptr,*zptr);
  */
}

void
really_translate_brain(float x, float y, float z)
{
  VERTEX *v;
  int i,j,k;
  float curr[4][4], accum[4][4]; /* Matrix curr,accum; */

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      curr[i][j] = idmat[i][j];
  curr[3][0] = x;
  curr[3][1] = y;
  curr[3][2] = z;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      accum[i][j] = 0;
      for (k=0;k<4;k++)
        accum[i][j] += curr[i][k]*reallymat[k][j];
    }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      reallymat[i][j] = accum[i][j];
  print_real_transform_matrix();

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->x += x;
    v->y += y;
    v->z += z;
  }
}

void
really_scale_brain(float x, float y, float z)
{
  VERTEX *v;
  int i,j,k;
  float curr[4][4], accum[4][4]; /* Matrix curr,accum; */

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      curr[i][j] = idmat[i][j];
  curr[0][0] = x;
  curr[1][1] = y;
  curr[2][2] = z;
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      accum[i][j] = 0;
      for (k=0;k<4;k++)
        accum[i][j] += curr[i][k]*reallymat[k][j];
    }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      reallymat[i][j] = accum[i][j];
  print_real_transform_matrix();

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->x *= x;
    v->y *= y;
    v->z *= z;
  }
}

int
MRIStransformBrain(MRI_SURFACE *mris,
                   float exx, float exy, float exz,
                   float eyx, float eyy, float eyz,
                   float ezx, float ezy, float ezz)
{
  int     vno ;
  VERTEX  *v ;
  float   x, y, z ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    v->x = x*exx + y*exy + z*exz ;
    v->y = x*eyx + y*eyy + z*eyz ;
    v->z = x*ezx + y*ezy + z*ezz ;
  }

  return(NO_ERROR) ;
}
void
align_sphere(MRI_SURFACE *mris)
{
  int     vno1, vno2 ;
  VERTEX  *v1, *v2 ;
  double  ex[3], ey[3], ez[3], len, p2[3], dot, tmp[3] ;

  /*
    two points must be selected on the spherical surface. The first
    will specify the z-axis (the 1st pole of the sphere), and the second
    the x-axis (the orientation of the vertical meridian.
  */
  if (nmarked != 2)
    fprintf(stderr, "must pick origin and alignment points (previous two)\n");
  else
  {
    vno1 = marked[nmarked-2] ;
    vno2 = marked[nmarked-1] ;
    v1 = &mris->vertices[vno1] ;
    v2 = &mris->vertices[vno2] ;
    ez[0] = v1->x ;
    ez[1] = v1->y ;
    ez[2] = v1->z ;
    len = VLEN(ez) ;
    SCALAR_MUL(ez, 1.0f/len, ez) ;
    p2[0] = v2->x ;
    p2[1] = v2->y ;
    p2[2] = v2->z ;
    dot = DOT(ez, p2) ;
    SCALAR_MUL(tmp, -dot, ez) ;   /* take ez component out of ex */
    ADD(ex, p2, tmp) ;
    len = VLEN(ex) ;
    SCALAR_MUL(ex, 1.0f/len, ex) ;
    CROSS(ey, ez, ex) ;
    len = VLEN(ey) ;
    SCALAR_MUL(ey, 1.0f/len, ey) ;
    dot = DOT(ez, ex) ;
    dot = DOT(ez, ey) ;
    dot = DOT(ey, ex) ;
    MRIStransformBrain(mris,
                       ex[0],ex[1],ex[2],
                       ey[0],ey[1],ey[2],
                       ez[0],ez[1],ez[2]);
    MRIScomputeMetricProperties(mris) ;
    redraw() ;
    printf("e0: (%2.2f, %2.2f, %2.2f), "
           "e1: (%2.2f, %2.2f, %2.2f), "
           "e2: (%2.2f, %2.2f, %2.2f)\n",
           ex[0], ex[1], ex[2],
           ey[0], ey[1], ey[2],
           ez[0], ez[1], ez[2]) ;
    PR
  }
}

void
really_rotate_brain(float a, char axis)
{
  VERTEX *v;
  float x,y,z;
  float cosa,sina;
  int i,j,k;
  float curr[4][4], accum[4][4]; /* Matrix curr,accum; */

  a = a*M_PI/180;
  cosa = cos(a);
  sina = sin(a);

  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      curr[i][j] = idmat[i][j];
  if (axis=='y')
  {
    curr[0][0] = curr[2][2] = cosa;
    curr[2][0] = -(curr[0][2] = sina);
  }
  else if (axis=='x')
  {
    curr[1][1] = curr[2][2] = cosa;
    curr[1][2] = -(curr[2][1] = sina);
  }
  else if (axis=='z')
  {
    curr[0][0] = curr[1][1] = cosa;
    curr[1][0] = -(curr[0][1] = sina);
  }
  else
  {
    printf("surfer: ### Illegal axis %c\n",axis);
    return;
    PR
  }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
    {
      accum[i][j] = 0;
      for (k=0;k<4;k++)
        accum[i][j] += curr[i][k]*reallymat[k][j];
    }
  for (i=0;i<4;i++)
    for (j=0;j<4;j++)
      reallymat[i][j] = accum[i][j];
  print_real_transform_matrix();

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;
    if (axis=='x')
    {
      v->y = y*cosa - z*sina;
      v->z = y*sina + z*cosa;
    }
    else if (axis=='y')
    {
      v->x = x*cosa + z*sina;
      v->z = -x*sina + z*cosa;
    }
    else if (axis=='z')
    {
      v->x = x*cosa - y*sina;
      v->y = x*sina + y*cosa;
    }
    else
    {
      printf("surfer: ### Illegal axis %c\n",axis);
      return;
      PR
    }
  }
}

void
really_align_brain(void)  /* trans cent first -> cent targ */
{
  VERTEX *v;
  VERTEX *v2;
  int k;
  float cx,cy,cz;
  float tcx,tcy,tcz;

  if (!secondsurfaceloaded)
  {
    printf("surfer: ### really_align brain failed: no second surface read\n");
    PR return;
  }

  cx = cy = cz = 0.0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    cx += v->x;
    cy += v->y;
    cz += v->z;
  }
  cx /= mris->nvertices;
  cy /= mris->nvertices;
  cz /= mris->nvertices;

  tcx = tcy = tcz = 0.0;
  for (k=0;k<mris2->nvertices;k++)
  {
    v2 = &mris2->vertices[k];
    tcx += v2->x;
    tcy += v2->y;
    tcz += v2->z;
  }
  tcx /= mris2->nvertices;
  tcy /= mris2->nvertices;
  tcz /= mris2->nvertices;

  really_translate_brain(tcx-cx,tcy-cy,tcz-cz);
  printf("surfer: center of first brain aligned to target brain\n");
  PR
}

void
really_center_brain(void)
{
  VERTEX *v;
  int k;
  float cx,cy,cz;

  cx = cy = cz = 0.0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    cx += v->x;
    cy += v->y;
    cz += v->z;
  }
  cx /= mris->nvertices;
  cy /= mris->nvertices;
  cz /= mris->nvertices;

  mris->xctr = mris->zctr = mris->yctr = 0.0;
  really_translate_brain(-cx,-cy,-cz);
  printf("surfer: avg vertex (%f, %f, %f) moved to (0,0,0)\n",cx,cy,cz);
  PR
}

void
really_center_second_brain(void)
{
  VERTEX *v2;
  int k;
  float cx,cy,cz;

  if (!secondsurfaceloaded)
  {
    printf("surfer: ### really_center_second_brain failed: "
           "no second surface read\n");
    PR return;
  }

  cx = cy = cz = 0.0;
  for (k=0;k<mris2->nvertices;k++)
  {
    v2 = &mris2->vertices[k];
    cx += v2->x;
    cy += v2->y;
    cz += v2->z;
  }
  cx /= mris2->nvertices;
  cy /= mris2->nvertices;
  cz /= mris2->nvertices;

  mris2->xctr = mris2->zctr = mris2->yctr = 0.0;

  for (k=0;k<mris2->nvertices;k++)
  {  /* really translate */
    v2 = &mris2->vertices[k];
    v2->x += -cx;
    v2->y += -cy;
    v2->z += -cz;
  }
  printf("surfer: avg vertex (%f, %f, %f) surf2 moved to (0,0,0)\n",
         cx,cy,cz);
  PR
}

void
print_real_transform_matrix(void)
{
  int i,j;

  printf("surfer: accumulated real transforms\n");
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
    {
      printf(" %13.3e ",reallymat[i][j]);
    }
    printf("\n");
  }
  PR
}

void
write_really_matrix(char *dir)
{
  int i,j;
  char fname[NAME_LENGTH];
  FILE *fp;

  sprintf(fname,"%s/really.mat",dir);
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  for (i=0;i<4;i++)
  {
    for (j=0;j<4;j++)
    {
      fprintf(fp,"%13.3e ",reallymat[i][j]);
    }
    fprintf(fp,"\n");
  }
  fclose(fp);
  printf("surfer: file %s written\n",fname);
  PR
}

void
read_really_matrix(char *dir)
{
  VERTEX *v;
  int i,j,k;
  float a,b,c,d;
  float or[4],r[4];
  char line[NAME_LENGTH];
  char fname[NAME_LENGTH];
  FILE *fp;

  sprintf(fname,"%s/really.mat",dir);
  fp = fopen(fname,"r");
  if (fp==NULL)
  {
    printf("surfer: ### File %s not found\n",fname);
    PR return;
  }

  i = 0;
  while (fgets(line,NAME_LENGTH,fp) != NULL)
  {
    if (sscanf(line,"%f %f %f %f",&a,&b,&c,&d) == 4)
    {
      reallymat[i][0] = a;
      reallymat[i][1] = b;
      reallymat[i][2] = c;
      reallymat[i][3] = d;
      i++;
    }
    else
    {
      printf("surfer: ### couldn't parse this line in matrix file:  %s",line);
      printf("surfer: ###   ...read_really_matrix() failed\n");
      PR return;
    }
  }
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    or[0] = v->x;
    or[1] = v->y;
    or[2] = v->z;
    or[3] = 1.0;
    r[0] = r[1] = r[2] = r[3] = 0.0;
    for (i=0;i<4;i++)
      for (j=0;j<4;j++)
      {
        r[i] += reallymat[i][j]*or[j];
      }
    v->x = r[0];
    v->y = r[1];
    v->z = r[2];
  }
  print_real_transform_matrix();
  printf("surfer: file %s read,applied\n",fname);
  PR
}

void
flatten(char *dir)
{
  float x,y,z,d,d1,d2;
  float nx,ny,nz;
  VERTEX *v;
  int k,an;
  FILE *fp;
  char fname[NAME_LENGTH];

  x = y = z = nx = ny = nz = 0;
  an = 0;
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      x += v->x;
      y += v->y;
      z += v->z;
      nx += v->nx;
      ny += v->ny;
      nz += v->nz;
      an++;
    }
  x /= an;
  y /= an;
  z /= an;
  printf("surfer: avg p = {%f,%f,%f}\n",x,y,z);
  printf("surfer: sum n = {%f,%f,%f}\n",nx,ny,nz);
  /* or override with direct front,back */
  if (project==POSTERIOR)
  {
    nx = nz = 0.0;
    ny = -1.0;
  }
  if (project==ANTERIOR)
  {
    nx = nz = 0.0;
    ny = 1.0;
  }
  d = sqrt(nx*nx+ny*ny+nz*nz);
  nx /= d;
  ny /= d;
  nz /= d;
  d = sqrt(nx*nx+ny*ny);
  printf("surfer: norm n = {%f,%f,%f}\n",nx,ny,nz);
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      v->x -= x;
      v->y -= y;
      v->z -= z;
      d1 = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
      transform(&v->x,&v->y,&v->z,nx,ny,nz,d);
      d2 = sqrt(v->x*v->x+v->y*v->y+v->z*v->z);
      if (fabs(d1-d2)>0.0001) printf("surfer: d1=%f, d2=%f\n",d1,d2);
      transform(&v->nx,&v->ny,&v->nz,nx,ny,nz,d);
    }

  /* print transform matrix in tmp dir */
  sprintf(fname,"%s/surfer.mat",dir);
  fp = fopen(fname,"w");
  if (fp==NULL)
  {
    printf("surfer: ### can't create file %s\n",fname);
    PR return;
  }
  else
  {
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",  nx*nz/d,  -nx,  ny/d,   0.0);
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",    d,       nz,   0.0,   0.0);
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",  -ny*nz/d,  ny,   nx/d,  0.0);
    fprintf(fp,"%13.3e %13.3e %13.3e %13.3e\n",   0.0,      0.0,  0.0,   1.0);
    fclose(fp);
    printf("surfer: file %s written\n",fname);
  }

  transform(&nx,&ny,&nz,nx,ny,nz,d);
  printf("surfer: transformed n = {%f,%f,%f}\n",nx,ny,nz);
  PR for (k=0;k<mris->nvertices;k++)
  {
    mris->vertices[k].z = 0;
  }
  compute_normals();
}

void
area_shrink(int niter)  /* consider area */
{
#if 0
  float x,y,z,sx,sy,sz,nc,force,nforce;
  float d,dx,dy,dz,sval,sinval,soutval,snc;
  float nx,ny,nz;
  VERTEX *v,*v1;
  int iter,k,m,n;
  float sd,ad,dmax,dd;
  int navg,an,nclip;
  double area,logarat,avgarea,vararea,minrat,maxrat;
  float ax,ay,az,tx,ty,tz;
  int nneg;

  for (iter=0;iter<niter;iter++)
  {
    smooth_logarat(10);
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
      /*      v->onc = v->nc;*/
      v->oripflag = v->ripflag;
    }
    sval = sinval = soutval = snc = 0;
    navg = 0;
    nneg = 0;
    avgarea = vararea = 0;
    maxrat = -100;
    minrat = 100;
    dd = 0;
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].oripflag)
      {
        v = &mris->vertices[k];
        x = v->tx;
        y = v->ty;
        z = v->tz;
        nx = v->nx;
        ny = v->ny;
        nz = v->nz;
        area = v->area;
        logarat = v->logarat;
        sx=sy=sz=sd=0;
        ax=ay=az=0;
        tx=ty=tz=0;
        n=0;
        avgarea += logarat;
        vararea += logarat*logarat;
        if (logarat>maxrat) maxrat = logarat;
        if (logarat<minrat) minrat = logarat;
        if (area/v->origarea<0) nneg++;

        /*v->curv = logarat;*/   /* marty: out 10/9/97 */

        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].oripflag)
          {
            v1 = &mris->vertices[v->v[m]];
            sx += dx = v1->tx - x;
            sy += dy = v1->ty - y;
            sz += dz = v1->tz - z;

            d = sqrt(dx*dx+dy*dy+dz*dz);
            if (d>0)
            {
              ax += (v1->tx - x)/d*(v1->logarat-logarat);
              ay += (v1->ty - y)/d*(v1->logarat-logarat);
              az += (v1->tz - z)/d*(v1->logarat-logarat);
            }
            if (v->border)
            {
              force = -logarat;
              d = sqrt(x*x+y*y+z*z);
              tx += x/d*force;
              ty += y/d*force;
              tz += z/d*force;
            }
            sd += d = sqrt(dx*dx+dy*dy+dz*dz);
            dd += SQR(d-1.0);
            n++;
          }
        if (n>0)
        {
          sx = sx/n;
          sy = sy/n;
          sz = sz/n;
          tx = tx/n;
          ty = ty/n;
          tz = tz/n;
          ax = ax/n;
          ay = ay/n;
          az = az/n;
          sd = sd/n;
        }
        if ((d=sqrt(ax*ax+ay*ay+az*az))>1.0)
        {
          ax /= d;
          ay /= d;
          az /= d;
        }
        force = 0;
        nc = sx*nx+sy*ny+sz*nz;
        sx -= nc*nx;
        sy -= nc*ny;
        sz -= nc*nz;
        snc += nc;
        v->nc = nc;
        nforce = 0;

        v->val = nc;

        if (logarat<0)
        {
          nforce = -logarat;
        }
        if (ncthreshflag)
        {
          if (nc<0)
            nc = (nc<-ncthresh)?nc+ncthresh:0;
          else
            nc = (nc>ncthresh)?nc-ncthresh:0;
        }
        dx = tx*wt + ax*wa + sx*ws + nx*nc*wn + nx*nforce*wc;
        dy = ty*wt + ay*wa + sy*ws + ny*nc*wn + ny*nforce*wc;
        dz = tz*wt + az*wa + sz*ws + nz*nc*wn + nz*nforce*wc;
        if (momentumflag)
        {
          dx = decay*v->odx+update*dx;
          dy = decay*v->ody+update*dy;
          dz = decay*v->odz+update*dz;
        }
        d=sqrt(dx*dx+dy*dy+dz*dz);
        /*
          if (!(d<1))
          printf("surfer: k=%d: d=%f, (%f,%f,%f),
          s=(%f,%f,%f), a=(%f,%f,%f)\n",
          k,d,dx,dy,dz,sd,sy,sz,ax,ay,az);
        */
        if (d>1.0)
        {
          nclip ++;
          dx /= d;
          dy /= d;
          dz /= d;
        }
        v->odx = dx;
        v->ody = dy;
        v->odz = dz;
        d=sqrt(dx*dx+dy*dy+dz*dz);
        if (d>dmax) dmax = d;
        ad += d;
        an ++;
        if (!(v->border&&fixed_border))
        {
          v->x += dx;
          v->y += dy;
          v->z += dz;
        }
      }
    ad /= an;
    snc /= an;
    avgarea /= an;
    vararea = sqrt(vararea/an - avgarea*avgarea);
    dd = sqrt(dd/an);

    compute_normals();

    printf("surfer: %d: ad=%f, dmax=%f, snc=%f dd=%f\n",iter,ad,dmax,snc,dd);
    printf("surfer:    avg=%f, var=%f, min=%f, max=%f, nneg=%d\n",
           avgarea,vararea,minrat,maxrat,nneg);
    PR
  }
#endif
}

void
shrink2d(int niter)
{
#if 0
  float x,y,z,sx,sy,sz,nc,force,nforce;
  float d,dx,dy,dz,sval,sinval,soutval,snc;
  float nx,ny,nz;
  VERTEX *v,*v1;
  int iter,k,m,n;
  float sd,ad,dmax,dd;
  int navg,an,nclip;
  double area,logarat,avgarea,vararea,minrat,maxrat;
  double shearx,sheary,logshear,avgshear,varshear,minshear,maxshear;
  float ax,ay,az,shx,shy,shz,tx,ty,tz,tnx,tny,tnz,f1,f2;
  int nneg;

  if (!computed_shear_flag)
  {
    smooth_logarat(10);
    compute_shear();
    smooth_shear(10);
    compute_boundary_normals();
    smooth_boundary_normals(10);
    computed_shear_flag = TRUE;
  }

  for (iter=0;iter<niter;iter++)
  {
    /*
      smooth_momentum(10);
    */
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
      /*      v->onc = v->nc;*/
      v->oripflag = v->ripflag;
    }
    sval = sinval = soutval = snc = 0;
    navg = 0;
    nneg = 0;
    avgarea = vararea = 0;
    maxrat = -100;
    minrat = 100;
    avgshear = varshear = 0;
    maxshear = -100;
    minshear = 100;
    dd = 0;
    for (k=0;k<mris->nvertices;k++)
      if (!mris->vertices[k].ripflag)
      {
        v = &mris->vertices[k];
        x = v->tx;
        y = v->ty;
        z = v->tz;
        nx = v->nx;
        ny = v->ny;
        nz = v->nz;
        area = v->area;
        logarat = v->logarat;
        logshear = v->logshear;
        shearx = v->shearx;
        sheary = v->sheary;
        sx=sy=sz=sd=0;
        ax=ay=az=0;
        shx=shy=shz=0;
        tnx=tny=tnz=0;
        tx=ty=tz=0;
        n=0;
        avgarea += logarat;
        vararea += logarat*logarat;
        if (logarat>maxrat) maxrat = logarat;
        if (logarat<minrat) minrat = logarat;

        avgshear += logshear;
        varshear += logshear*logshear;
        if (logshear>maxshear) maxshear = logshear;
        if (logshear<minshear) minshear = logshear;
        /*
          v->curv = logarat;
          v->curv = logshear;
        */

        for (m=0;m<v->vnum;m++)
          if (!mris->vertices[v->v[m]].oripflag)
          {
            v1 = &mris->vertices[v->v[m]];
            sx += dx = v1->tx - x;
            sy += dy = v1->ty - y;
            sz += dz = v1->tz - z;

            d = sqrt(dx*dx+dy*dy+dz*dz);
            if (!v->border)
            {
              if (d>0)
              {
                ax += (v1->tx - x)/d*(v1->logarat-logarat);
                ay += (v1->ty - y)/d*(v1->logarat-logarat);
                az += (v1->tz - z)/d*(v1->logarat-logarat);
              }
              f1 = fabs((v1->tx-x)*v1->shearx+(v1->ty-y)*v1->sheary)-
                   fabs(-(v1->tx-x)*v1->sheary+(v1->ty-y)*v1->shearx);
              f2 = fabs((v1->tx-x)*shearx+(v1->ty-y)*sheary)-
                   fabs(-(v1->tx-x)*sheary+(v1->ty-y)*shearx);
              if (d>0)
              {
                shx += (v1->tx-x)/d*(f1-f2);
                shy += (v1->ty-y)/d*(f1-f2);
                shz += 0;
              }
            }
            if (v->border)
            {
              force = -logarat;
              if (v->nz<0)
              {
                force = 0.0;
                nneg++;
              }
              tx += v->bnx*force;
              ty += v->bny*force;
              tz += 0;
              f1 = -(fabs((v->bnx)*shearx+(v->bny)*sheary)-
                     fabs(-(v->bnx)*sheary+(v->bny)*shearx));
              if (v->nz<0)
              {
                f1 = 0.0;
              }
              tnx += v->bnx*f1;
              tny += v->bny*f1;
              tnz += 0;
            }
            sd += d = sqrt(dx*dx+dy*dy+dz*dz);
            dd += SQR(d-1.0);
            n++;
          }
        if (n>0)
        {
          sx = sx/n;
          sy = sy/n;
          sz = sz/n;
          tx = tx/n;
          ty = ty/n;
          tz = tz/n;
          tnx = tnx/n;
          tny = tny/n;
          tnz = tnz/n;
          ax = ax/n;
          ay = ay/n;
          az = az/n;
          shx = shx/n;
          shy = shy/n;
          shz = shz/n;
          sd = sd/n;
        }
        if ((d=sqrt(ax*ax+ay*ay+az*az))>1.0)
        {
          ax /= d;
          ay /= d;
          az /= d;
        }
        if ((d=sqrt(shx*shx+shy*shy+shz*shz))>1.0)
        {
          shx /= d;
          shy /= d;
          shz /= d;
        }
        force = 0;
        nc = sx*nx+sy*ny+sz*nz;
        sx -= nc*nx;
        sy -= nc*ny;
        sz -= nc*nz;
        snc += nc;
        v->nc = nc;
        nforce = 0;
        /*
          v->val = nc;
        */
        /* mgh: omit !, then comment out */
        /*
          if (!v->border)
          {
          bnc = sx*v->bnx+sy*v->bny;
          sx -= bnc*nx;
          sy -= bnc*ny;
          }
        */

        if (logarat<0)
        {
          nforce = -logarat;
        }
        if (nc<0)
          nc = (nc<-ncthresh)?nc+ncthresh:0;
        else
          nc = (nc>ncthresh)?nc-ncthresh:0;
        dx = tx*wt + ax*wa + sx*ws + shx*wsh + tnx*wbn;
        dy = ty*wt + ay*wa + sy*ws + shy*wsh + tny*wbn;
        dz = tz*wt + az*wa + sz*ws + shz*wsh + tnz*wbn;
        if (momentumflag)
        {
          dx = decay*v->odx+update*dx;
          dy = decay*v->ody+update*dy;
          dz = decay*v->odz+update*dz;
        }
        d=sqrt(dx*dx+dy*dy+dz*dz);
        if (d>1.0)
        {
          nclip ++;
          dx /= d;
          dy /= d;
          dz /= d;
        }
        v->odx = dx;
        v->ody = dy;
        v->odz = dz;
        /*
          v->odx = 0.8*dx+0.5*v->smx;
          v->ody = 0.8*dy+0.5*v->smy;
          v->odz = 0.8*dz+0.5*v->smz;
        */
        d=sqrt(v->odx*v->odx+v->ody*v->ody+v->odz*v->odz);
        if (d>1.0)
        {
          v->odx /= d;
          v->ody /= d;
          v->odz /= d;
        }
        d=sqrt(dx*dx+dy*dy+dz*dz);
        if (d>dmax) dmax = d;

        if (!((d>=0)||(d<=0)))
          printf("surfer: d=%f, a:(%f,%f,%f), "
                 "sh:(%f,%f,%f), s:(%f,%f,%f)\n",
                 d,ax,ay,az,shx,shy,shz,sx,sy,sz);

        ad += d;
        an ++;
        if (!(v->border&&fixed_border))
        {
          v->x += dx;
          v->y += dy;
          v->z += dz;
        }
      }
    ad /= an;
    snc /= an;
    avgarea /= an;
    vararea = sqrt(vararea/an - avgarea*avgarea);
    avgshear /= an;
    varshear = sqrt(varshear/an - avgshear*avgshear);
    dd = sqrt(dd/an);

    compute_normals();
    smooth_logarat(10);
    compute_shear();
    smooth_shear(10);
    compute_boundary_normals();
    smooth_boundary_normals(10);

    printf("surfer: %d: ad=%f, dmax=%f, dd=%f nneg=%d\n",
           iter,ad,dmax,dd,nneg);
    printf("surfer:     area:  avg=%f, var=%f, min=%f, max=%f\n",
           avgarea,vararea,minrat,maxrat);
    printf("surfer:     shear: avg=%f, var=%f, min=%f, max=%f\n",
           avgshear,varshear,minshear,maxshear);
    PR
  }
#endif
}

void
sphere_shrink(int niter, float rad)
{
  float x,y,z,sx,sy,sz;
  float d,dx,dy,dz;
  float xc,yc,zc,r,dr;
  VERTEX *v;
  int iter,k,m,n;
  float sd,ad,dmax;
  int an,nclip;

  for (iter=0;iter<niter;iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
    }
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      x = v->tx;
      y = v->ty;
      z = v->tz;
      sx=sy=sz=sd=0;
      n=0;
      for (m=0;m<v->vnum;m++)
      {
        sx += dx = mris->vertices[v->v[m]].tx - x;
        sy += dy = mris->vertices[v->v[m]].ty - y;
        sz += dz = mris->vertices[v->v[m]].tz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      if (n>0)
      {
        sx = sx/n;
        sy = sy/n;
        sz = sz/n;
        sd = sd/n;
      }
      xc = x+mris->xctr;
      yc = y-mris->yctr;
      zc = z-mris->zctr;
      r = sqrt(xc*xc+yc*yc+zc*zc);
      if (r==0) r=0.00001;
      dr = (rad-r)/rad;
      dx = sx*0.5 + xc/r*dr;
      dy = sy*0.5 + yc/r*dr;
      dz = sz*0.5 + zc/r*dr;
      if (momentumflag)
      {
        dx = decay*v->odx+update*dx;
        dy = decay*v->ody+update*dy;
        dz = decay*v->odz+update*dz;
        if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
        {
          nclip ++;
          dx /= d;
          dy /= d;
          dz /= d;
        }
        v->odx = dx;
        v->ody = dy;
        v->odz = dz;
      }
      d=sqrt(dx*dx+dy*dy+dz*dz);
      if (d>dmax) dmax = d;
      ad += d;
      an ++;
      v->x += dx;
      v->y += dy;
      v->z += dz;
    }
    printf("surfer: %d: ad=%f, dmax=%f, nclip=%d\n",
           iter,ad/an,dmax,nclip);
    PR
  }
  compute_normals();
}

/* a=rh/lh, b=ant/post,  c=sup/inf */
void
ellipsoid_project(float a, float b, float c)
{
  VERTEX *v;
  int k;
  float x,y,z,x2,y2,z2,dx,dy,dz,a2,b2,c2,a4,b4,c4,a6,b6,c6;
  float f,g,h,d,dist,avgdist=0;

  printf("ellipsoid_project(%f,%f,%f)\n",a,b,c);

  a2 = a*a;
  b2 = b*b;
  c2 = c*c;
  a4 = a2*a2;
  b4 = b2*b2;
  c4 = c2*c2;
  a6 = a2*a4;
  b6 = b2*b4;
  c6 = c2*c4;

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    /*
      printf("%6d: before: %6.2f\n",k,SQR(v->x/a)+SQR(v->y/b)+SQR(v->z/c));
    */
    x = v->x;
    y = v->y;
    z = v->z;
    x2 = x*x;
    y2 = y*y;
    z2 = z*z;
    f = x2/a6+y2/b6+z2/c6;
    g = 2*(x2/a4+y2/b4+z2/c4);
    h = x2/a2+y2/b2+z2/c2-1;
    d = (-g+sqrt(g*g-4*f*h))/(2*f);
    dx = d*x/a2;
    dy = d*y/b2;
    dz = d*z/c2;
    v->x = x+dx;
    v->y = y+dy;
    v->z = z+dz;
    dist = sqrt(dx*dx+dy*dy+dz*dz);
    avgdist += dist;
    /*
      printf("%6d: after: %6.2f\n",k,SQR(v->x/a)+SQR(v->y/b)+SQR(v->z/c));
    */
  }
  /*
    printf("ellipsoid_project: avgdist = %f\n",avgdist/mris->nvertices);
  */
  compute_normals();
}

/* a=rh/lh, b=ant/post,  c=sup/inf */
void
ellipsoid_morph(int niter, float a, float b, float c)
{
  VERTEX *v;
  int imnr,i,j,iter,k,m,n,dk,di,dj;
  int an,nclip,delcurvdefined;
  float x,y,z,sx,sy,sz,val;
  float xnei,ynei,znei,gradx,grady,gradz;
  float d,dx,dy,dz;
  float sd,ad,dmax;
  float curv,icurv = 0.0f;
  float icrange,crange;
  double gradnormsum,rmscurv,rmsicurv,rmscurverr,curvfact,rmsinum;

  if (curvloaded && curvimloaded)
  {
    icrange = mris2->max_curv-mris2->min_curv;
    crange = mris->max_curv-mris->min_curv;
  }

  printf("ellipsoid_morph(%d,%f,%f,%f)\n",niter,a,b,c);

  ellipsoid_project(a,b,c);
  rmscurv = rmsicurv = rmsinum = 0;
  if (curvimflag && curvimloaded)
    for (imnr=0;imnr<numimg;imnr++)
      for (i=0;i<IMGSIZE;i++)
        for (j=0;j<IMGSIZE;j++)
          if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
          {
            icurv = bytetofloat(curvim[imnr][i][j],
                                mris2->min_curv,mris2->max_curv);
            rmsicurv += SQR(icurv);
            rmsinum++;
          }
  rmsicurv = sqrt(rmsicurv/rmsinum);
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    curv = v->curv;
    rmscurv += SQR(curv);
  }
  rmscurv = sqrt(rmscurv/mris->nvertices);
  curvfact = rmsicurv/rmscurv;
  printf("rmscurv = %f, rmsicurv = %f, curvfact = %f\n",
         rmscurv,rmsicurv,curvfact);
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->curv *= curvfact;
  }
  for (iter=0;iter<niter;iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
    }
    gradnormsum = rmsicurv = rmscurverr = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      x = v->tx;
      y = v->ty;
      z = v->tz;
      sx=sy=sz=sd=0;
      n=0;
      if (curvimflag && curvimloaded)
      {
        delcurvdefined = TRUE;
        imnr = (int)((y-yy0)/st+0.5);
        i = (int)((zz1-z)/ps+0.5);
        j = (int)((xx1-x)/ps+0.5);
        imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
        i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
        j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
        if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
          icurv = bytetofloat(curvim[imnr][i][j],
                              mris2->min_curv,mris2->max_curv);
        else
          delcurvdefined = FALSE;
        curv = v->curv;
        rmsicurv += SQR(icurv);
        rmscurverr += SQR(curv-icurv);
        gradx = grady = gradz = 0;
        for (dk=-1;dk<=1;dk++)
          for (di=-1;di<=1;di++)
            for (dj=-1;dj<=1;dj++)
            {
              if (!(curvimflags[imnr+dk][i+di][j+dj] & CURVIM_DEFINED))
                delcurvdefined = FALSE;
              val = bytetofloat(curvim[imnr+dk][i+di][j+dj],
                                mris2->min_curv,mris2->max_curv);
              gradx += -dj*val;
              grady +=  dk*val;
              gradz += -di*val;
            }
        if (!delcurvdefined)
          printf("undefined gradient at vertex %d (%f,%f,%f)\n",k,x,y,z);
        gradx /= 3*3*2;
        grady /= 3*3*2;
        gradz /= 3*3*2;
        sx += icstrength*gradx*(curv-icurv);
        sy += icstrength*grady*(curv-icurv);
        sz += icstrength*gradz*(curv-icurv);
        gradnormsum += sqrt(gradx*gradx+grady*grady+gradz*gradz);
      }
      for (m=0;m<v->vnum;m++)
      {
        xnei = mris->vertices[v->v[m]].tx;
        ynei = mris->vertices[v->v[m]].ty;
        znei = mris->vertices[v->v[m]].tz;
        sx += dx = xnei - x;
        sy += dy = ynei - y;
        sz += dz = znei - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      if (n>0)
      {
        sx = sx/n;
        sy = sy/n;
        sz = sz/n;
        sd = sd/n;
      }
      dx = sx;
      dy = sy;
      dz = sz;
      if (momentumflag)
      {
        dx = decay*v->odx+update*dx;
        dy = decay*v->ody+update*dy;
        dz = decay*v->odz+update*dz;
        if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
        {
          nclip ++;
          dx /= d;
          dy /= d;
          dz /= d;
        }
        v->odx = dx;
        v->ody = dy;
        v->odz = dz;
      }
      d=sqrt(dx*dx+dy*dy+dz*dz);
      if (d>dmax) dmax = d;
      ad += d;
      an ++;
      v->x += dx;
      v->y += dy;
      v->z += dz;
    }
    printf("surfer: %d: ad=%f, dmax=%f, nclip=%d\n",
           iter,ad/an,dmax,nclip);
    PR ellipsoid_project(a,b,c);
    printf("average gradient = %f, curvature error = %e\n",
           gradnormsum/mris->nvertices,sqrt(rmscurverr/rmsicurv));
  }
  compute_normals();
}

/* a=rh/lh, b=ant/post,  c=sup/inf */
void
ellipsoid_shrink(int niter, float a, float b, float c)
{
  VERTEX *v;
  int imnr,i,j,iter,k,m,n;
  int an,nclip,delcurvdefined = 0;
  float x,y,z,sx,sy,sz;
  float xnei,ynei,znei;
  float d,dx,dy,dz;
  float xc,yc,zc,r,dr;
  float sd,ad,dmax;
  float acx,acy,acz;
  float sqd,rad;
  float curv,icurv,icurvnei,cfact;
  float icrange,crange;

  curv = icurv = icurvnei = 0.0f ;
  if (curvloaded && curvimloaded)
  {
    icrange = mris2->max_curv-mris2->min_curv;
    crange = mris->max_curv-mris->min_curv;
  }

  acx = acy = acz = 0.0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    acx += v->x;
    acy += v->y;
    acz += v->z;
  }
  acx /= mris->nvertices;
  acy /= mris->nvertices;
  acz /= mris->nvertices;

  for (iter=0;iter<niter;iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
    }
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      x = v->tx;
      y = v->ty;
      z = v->tz;
      sx=sy=sz=sd=0;
      n=0;
      if (curvimflag && curvimloaded)
      {
        delcurvdefined = TRUE;
        imnr = (int)((y-yy0)/st+0.5);
        i = (int)((zz1-z)/ps+0.5);
        j = (int)((xx1-x)/ps+0.5);
        imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
        i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
        j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
        if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
          icurv = bytetofloat(curvim[imnr][i][j],
                              mris2->min_curv,mris2->max_curv);
        else
          delcurvdefined = FALSE;
        curv = v->curv;
      }
      for (m=0;m<v->vnum;m++)
      {
        xnei = mris->vertices[v->v[m]].tx;
        ynei = mris->vertices[v->v[m]].ty;
        znei = mris->vertices[v->v[m]].tz;
        if (curvimflag && curvimloaded)
        {
          imnr = (int)((ynei-yy0)/st+0.5);
          i = (int)((zz1-znei)/ps+0.5);
          j = (int)((xx1-xnei)/ps+0.5);
          imnr = (imnr<0)?0:(imnr>=numimg)?numimg-1:imnr;
          i = (i<0)?0:(i>=IMGSIZE)?IMGSIZE-1:i;
          j = (j<0)?0:(j>=IMGSIZE)?IMGSIZE-1:j;
          if (curvimflags[imnr][i][j] & CURVIM_DEFINED)
            icurvnei = bytetofloat(curvim[imnr][i][j],
                                   mris2->min_curv,mris2->max_curv);
          else
            delcurvdefined = FALSE;
          /*curvnei = mris->vertices[v->v[m]].curv;*/   /* two dels?? */
          cfact = 1.0;
          if (delcurvdefined)
            cfact += icstrength*(icurvnei-icurv)*(curv-icurv);
          /* del im w/only sign of im/loc diff */
          /*cfact += (icurvnei-icurv) * copysign(icstrength,curv-icurv);*/
          /*cfact += (icurvnei-icurv)/icrange
           * copysign(icstrength,mris2->min_curv+curv*
           (icrange/crange) - icurv);*/
          sx += dx = (xnei - x)*cfact;
          sy += dy = (ynei - y)*cfact;
          sz += dz = (znei - z)*cfact;
        }
        else
        {
          sx += dx = xnei - x;
          sy += dy = ynei - y;
          sz += dz = znei - z;
        }
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      if (n>0)
      {
        sx = sx/n;
        sy = sy/n;
        sz = sz/n;
        sd = sd/n;
      }
      xc = x-acx;
      yc = y-acy;
      zc = z-acz;
      r = sqrt(xc*xc+yc*yc+zc*zc);
      if (r==0) r=0.00001;
      sqd = SQR(xc/a) + SQR(yc/b) + SQR(zc/c);
      rad = sqrt(SQR(xc)/sqd + SQR(yc)/sqd + SQR(zc)/sqd);  /* ellipsoid */
      dr = (rad-r)/rad;
      dx = sx*0.5 + xc/r*dr;      /* radial (tangential comp on ellipsoid) */
      dy = sy*0.5 + yc/r*dr;
      dz = sz*0.5 + zc/r*dr;
#if 0
      if (dr > 0.01 || dr < -0.01)
      {
        dx = sx*0.5 + xc/r*dr;   /* radial approach */
        dy = sy*0.5 + yc/r*dr;
        dz = sz*0.5 + zc/r*dr;
      }
      else
      {
        dx = sx*0.5 + dr*v->nx;  /* normal there (unstable if folded) */
        dy = sy*0.5 + dr*v->ny;
        dz = sz*0.5 + dr*v->nz;
      }
#endif
      if (momentumflag)
      {
        dx = decay*v->odx+update*dx;
        dy = decay*v->ody+update*dy;
        dz = decay*v->odz+update*dz;
        if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
        {
          nclip ++;
          dx /= d;
          dy /= d;
          dz /= d;
        }
        v->odx = dx;
        v->ody = dy;
        v->odz = dz;
      }
      d=sqrt(dx*dx+dy*dy+dz*dz);
      if (d>dmax) dmax = d;
      ad += d;
      an ++;
      v->x += dx;
      v->y += dy;
      v->z += dz;
    }
    printf("surfer: %d: ad=%f, dmax=%f, nclip=%d\n",iter,ad/an,dmax,nclip);
    PR
  }
  compute_normals();
}

/* 50 ellipsoid_shrink(2,100,150); */
void
ellipsoid_shrink_bug(int niter, float rad, float len)
{
  float x,y,z,sx,sy,sz;
  float d,dx,dy,dz;
  float xc,yc,zc,r,dr;
  float ex,ey,ez;
  float acx,acy,acz;
  VERTEX *v;
  int iter,k,m,n;
  float sd,ad,dmax;
  int an,nclip;

  ex = 0.0;
  ey = (len-rad)/2.0;
  ez = 0.0;

  acx = acy = acz = 0.0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    acx += v->x;
    acy += v->y;
    acz += v->z;
  }
  acx /= (float)mris->nvertices;
  acy /= (float)mris->nvertices;
  acz /= (float)mris->nvertices;

  for (iter=0;iter<niter;iter++)
  {
    ad = 0;
    dmax = 0;
    an = 0;
    nclip = 0;
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
    }
    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k];
      x = v->tx;
      y = v->ty;
      z = v->tz;
      sx=sy=sz=sd=0;
      n=0;
      for (m=0;m<v->vnum;m++)
      {
        sx += dx = mris->vertices[v->v[m]].tx - x;
        sy += dy = mris->vertices[v->v[m]].ty - y;
        sz += dz = mris->vertices[v->v[m]].tz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      if (n>0)
      {
        sx = sx/n;
        sy = sy/n;
        sz = sz/n;
        sd = sd/n;
      }
      /*
        xc = x+mris->xctr;
        yc = y-mris->yctr;
        zc = z-mris->zctr;
      */
      xc = x-acx + (x-acx)*ex;  /* sphere + sphere dot ellipse */
      yc = y-acy + (y-acy)*ey;
      zc = z-acz + (z-acz)*ez;
      r = sqrt(xc*xc+yc*yc+zc*zc);
      if (r==0) r=0.00001;
      dr = (rad-r)/rad;
      dx = sx*0.5 + xc/r*dr;
      dy = sy*0.5 + yc/r*dr;
      dz = sz*0.5 + zc/r*dr;
      if (momentumflag)
      {
        dx = decay*v->odx+update*dx;
        dy = decay*v->ody+update*dy;
        dz = decay*v->odz+update*dz;
        if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0)
        {
          nclip ++;
          dx /= d;
          dy /= d;
          dz /= d;
        }
        v->odx = dx;
        v->ody = dy;
        v->odz = dz;
      }
      d=sqrt(dx*dx+dy*dy+dz*dz);
      if (d>dmax) dmax = d;
      ad += d;
      an ++;
      v->x += dx;
      v->y += dy;
      v->z += dz;
    }
    printf("surfer: %d: ad=%f, dmax=%f, nclip=%d\n",iter,ad/an,dmax,nclip);
    PR
  }
  compute_normals();
}

void
compute_curvature(void)
{
  float x,y,z,dx,dy,dz,r;
  VERTEX *v;
  int k,m;
  float min_curv, max_curv;

  min_curv = 99999;
  max_curv = -99999;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    x = v->x;
    y = v->y;
    z = v->z;
    v->curv = 0;
    for (m=0;m<v->vnum;m++)
    {
      dx = mris->vertices[v->v[m]].x-x;
      dy = mris->vertices[v->v[m]].y-y;
      dz = mris->vertices[v->v[m]].z-z;
      r = sqrt(dx*dx+dy*dy+dz*dz);
      if (r>0)
        v->curv += (dx*v->nx+dy*v->ny+dz*v->nz)/r;
    }

    if (v->curv < min_curv) min_curv = v->curv;
    if (v->curv > max_curv) max_curv = v->curv;
  }

  /* begin rkt */

  /* save the curv min and max. */
  mris->min_curv = cmin = min_curv;
  mris->max_curv = cmax = max_curv;

  curvloaded = TRUE;
  enable_menu_set (MENUSET_CURVATURE_LOADED, 1);
  send_tcl_command ("ShowLabel kLabel_Curvature 1");

  /* turn on the curvflag */
  curvflag = 1;

  /* end rkt */
}

void
clear_curvature(void)
{
  VERTEX *v;
  int k;

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->curv = 0;
  }
}

void
clear_vals(void)
{
  VERTEX *v;
  int k;

  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    v->val = 0;
  }
}

void
normalize_curvature(int which_norm)
{
  MRISnormalizeCurvature(mris, which_norm) ;
}
void
normalize_area(void)
{
  VERTEX *v;
  int k;
  float a,oa,f;

  if (MRIflag && MRIloaded)
  {
    printf("surfer: ### normalize_area not done w/MRIflag: "
           "misaligns surface\n");
    PR return;
  }

  oa = a = 0;
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      a += v->area;
      oa += v->origarea;
    }
  f = sqrt(oa/a);
  printf("surfer: oa=%f sqmm, a=%f sqmm, f=%f\n",oa/2,a/2,f);
  PR for (k=0;k<mris->nvertices;k++)
  {
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      v->x *= f;
      v->y *= f;
      v->z *= f;
    }
  }
  compute_normals();
}

void
normalize_surface(void)
{
  VERTEX *v;
  int i,j,k;
  float x,y,z;
  float minx,maxx,miny,maxy,minz,maxz;
  float minx2,maxx2,miny2,maxy2,minz2,maxz2;

  minx2 = miny2 = minz2 = 1000;
  maxx2 = maxy2 = maxz2 = -1000;
  for (k=1;k<numimg-1;k++)
    for (i=1;i<IMGSIZE-1;i++)
      for (j=1;j<IMGSIZE-1;j++)
      {
        x = IMGSIZE/2.0-j;
        z = IMGSIZE/2.0-i;
        y = k-numimg/2.0;
        if (im[k][i][j]!=0)
        {
          if (x>maxx2) maxx2 = x;
          if (x<minx2) minx2 = x;
          if (y>maxy2) maxy2 = y;
          if (y<miny2) miny2 = y;
          if (z>maxz2) maxz2 = z;
          if (z<minz2) minz2 = z;
        }
      }
  minx = miny = minz = 1000;
  maxx = maxy = maxz = -1000;
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      if (v->x>maxx) maxx = v->x;
      if (v->x<minx) minx = v->x;
      if (v->y>maxy) maxy = v->y;
      if (v->y<miny) miny = v->y;
      if (v->z>maxz) maxz = v->z;
      if (v->z<minz) minz = v->z;
    }
  printf("surfer: minx2=%f, maxx2=%f\n",minx2,maxx2);
  printf("surfer: miny2=%f, maxy2=%f\n",miny2,maxy2);
  printf("surfer: minz2=%f, maxz2=%f\n",minz2,maxz2);
  printf("surfer: minx=%f, maxx=%f\n",minx,maxx);
  printf("surfer: miny=%f, maxy=%f\n",miny,maxy);
  printf("surfer: minz=%f, maxz=%f\n",minz,maxz);
  PR
  for (k=0;k<mris->nvertices;k++)
    if (!mris->vertices[k].ripflag)
    {
      v = &mris->vertices[k];
      v->x = minx2+(v->x-minx)*(maxx2-minx2)/(maxx-minx);
      v->y = miny2+(v->y-miny)*(maxy2-miny2)/(maxy-miny);
      v->z = minz2+(v->z-minz)*(maxz2-minz2)/(maxz-minz);
    }
  compute_normals();
}

void
load_brain_coords(float x, float y, float z, float v[])
{
  v[0] = -x;
  v[1] = z;
  v[2] = y;
}

int outside(float x,float y, float z)
{
  return (x<xmin||x>xmax
          ||y<ymin||y>ymax
          ||-z<zmin||-z>zmax);
}

static GLfloat *vertices = NULL ;
static GLfloat *normals = NULL ;

static float *colors=0;
static float *mesh_colors=0;

static void fill_color_array(MRI_SURFACE *mris, float *colors);
static GLuint *faces ;
static int init_vertex_arrays(MRI_SURFACE *mris) ;

static int get_color_vals(float val,
                          float curv,
                          int mode,
                          GLubyte *r, GLubyte *g, GLubyte *b) ;
static void get_complexval_color_vals
(
  float x, float y,
  float stat,
  float curv, GLubyte *r, GLubyte *g, GLubyte *b
) ;
static int fill_vertex_arrays(MRI_SURFACE *mris) ;
static int
init_vertex_arrays(MRI_SURFACE *mris)
{
#ifndef IRIX
  colors = (float *)calloc(3*VERTICES_PER_FACE*mris->nfaces, sizeof(float));
  if (!colors)
    ErrorExit(ERROR_NOMEMORY, "init_vertex_arrays: calloc failed") ;

  /* begin rkt */
  if (NULL != vertices)
    free (vertices);
  if (NULL != normals)
    free (normals);
  if (NULL != faces)
    free (faces);
  /* end rkt */

  vertices = (GLfloat *)calloc(3*mris->nvertices, sizeof(GLfloat)) ;
  normals = (GLfloat *)calloc(3*mris->nvertices, sizeof(GLfloat)) ;
  faces = (GLuint *)calloc(VERTICES_PER_FACE*mris->nfaces,
                           sizeof(unsigned int)) ;
  if (!vertices || !faces || !colors || !normals)
    ErrorExit(ERROR_NOMEMORY, "init_vertex_arrays: calloc failed") ;

  /* begin rkt */
  /* commented this stuff out */
#if 0
  fill_vertex_arrays(mris) ;

  /* glEnableClientState ( GL_NORMAL_ARRAY ); */

  glEnableClientState ( GL_VERTEX_ARRAY );
  glEnableClientState ( GL_COLOR_ARRAY );
  glEnableClientState ( GL_NORMAL_ARRAY );
  glVertexPointer(3, GL_FLOAT, 0, vertices) ;
  glNormalPointer(GL_FLOAT, 0, normals) ;
  glColorPointer(3, GL_FLOAT, 0, colors);
#endif
  /* end rkt */

#endif
  return(NO_ERROR) ;
}

static int init_mesh_colors(MRI_SURFACE *mris)
{
  int i;
  mesh_colors = (float *)calloc(3*VERTICES_PER_FACE*mris->nfaces,
                                sizeof(float));
  if (!mesh_colors)
    ErrorExit(ERROR_NOMEMORY, "init_mesh_colors: calloc failed") ;

  for (i=0; i<mris->nvertices; i++)
  {
    mesh_colors[3*i] = (float)meshr/255.0;
    mesh_colors[3*i+1] = (float)meshg/255.0;
    mesh_colors[3*i+2] = (float)meshb/255.0;
  }
  return(NO_ERROR) ;
}

static int fill_mesh_colors()
{
  int i;
  GLubyte r, g, b;
  for (i=0; i<mris->nvertices; i++)
  {
    r = meshr;
    g = meshg;
    b = meshb;
    ddt_get_hilite_face_color (i, &r, &g, &b);
    mesh_colors[3*i] = (float)r/255.0;
    mesh_colors[3*i+1] = (float)g/255.0;
    mesh_colors[3*i+2] = (float)b/255.0;
  }
  return(NO_ERROR) ;
}
static int
fill_vertex_arrays(MRI_SURFACE *mris)
{
  int    vno, fno, i, n ;
  VERTEX *v ;
  FACE   *f ;

  for (i = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      i+=3;
      continue ;
    }
    vertices[i++] = -v->x ;
    vertices[i++] = v->z ;
    vertices[i++] = v->y ;
  }

  for (i = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
    {
      i+=3;
      continue ;
    }
    normals[i++] = -v->nx ;
    normals[i++] = v->nz ;
    normals[i++] = v->ny ;
  }

  for (i =  0; i<mris->nfaces*VERTICES_PER_FACE; i++)
    faces[i]=0;

  for (i = fno = 0 ; fno < mris->nfaces ; fno++)
  {
    f = &mris->faces[fno] ;
    if (f->ripflag)
      continue ;
    for (n = 0 ; n < VERTICES_PER_FACE ; n++)
      faces[i++] = f->v[n] ;
  }
  return(NO_ERROR) ;
}

void
draw_surface(void)  /* marty: combined three versions */
{

  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glDeleteLists(FS_Brain_List,1);
  if (!colors) init_vertex_arrays(mris) ;
  if (vertex_array_dirty==1)
  {
    fill_vertex_arrays(mris);
    glEnableClientState ( GL_VERTEX_ARRAY );
    glEnableClientState ( GL_COLOR_ARRAY );
    glEnableClientState ( GL_NORMAL_ARRAY );
    glVertexPointer(3, GL_FLOAT, 0, vertices) ;
    glNormalPointer(GL_FLOAT, 0, normals) ;
    glColorPointer(3, GL_FLOAT, 0, colors);
    vertex_array_dirty = 0;
  }
  if (color_scale_changed)
  {
    if(UseNewOverlay)
      fill_color_array2(mris, colors) ; // draws the overlay
    else
      fill_color_array(mris, colors) ; // draws the overlay
    color_scale_changed = TRUE;
    glColorPointer  ( 3, GL_FLOAT, 0, colors );
  }
  glPolygonOffset(1.0,1.0);
  glEnable(GL_POLYGON_OFFSET_FILL);

  /* Draw the object*/
  if (surfaceflag)
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    if (VERTICES_PER_FACE==3)
    {
      glDrawElements ( GL_TRIANGLES, 3*mris->nfaces, GL_UNSIGNED_INT,faces );
    }
    else
    {
      glDrawElements ( GL_QUADS, 4*mris->nfaces, GL_UNSIGNED_INT, faces );
    }
  }
  glDisable(GL_POLYGON_OFFSET_FILL);
  if (pointsflag)
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
    if (!mesh_colors)
      init_mesh_colors(mris);
    fill_mesh_colors();
    glColorPointer( 3, GL_FLOAT, 0, mesh_colors);
    if (VERTICES_PER_FACE==3)
    {
      glDrawElements ( GL_TRIANGLES, 3*mris->nfaces, GL_UNSIGNED_INT,faces );
    }
    else
    {
      glDrawElements ( GL_QUADS, 4*mris->nfaces, GL_UNSIGNED_INT, faces );
    }
  }
  if (verticesflag)
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    if (!mesh_colors)
      init_mesh_colors(mris);
    fill_mesh_colors();
    glColorPointer( 3, GL_FLOAT, 0, mesh_colors);
    if (VERTICES_PER_FACE==3)
    {
      glDrawElements( GL_TRIANGLES, 3*mris->nfaces, GL_UNSIGNED_INT, faces );
    }
    else
    {
      glDrawElements ( GL_QUADS, 4*mris->nfaces, GL_UNSIGNED_INT, faces );
    }
  }

  glFlush() ;
  if (scalebarflag)    draw_scalebar();
  if (colscalebarflag) draw_colscalebar();
  if (cptn_draw_flag)  cptn_draw();
  glFlush() ;
}

/*--------------------------------------------------------------*/
static void fill_color_array2(MRI_SURFACE *mris, float *colors)
{
  extern MRI *mrismask;
  int n;
  VERTEX *v;
  float  r, g, b, f;
  float min, mid, max;
  int maskout=0;

  min = (float)(fthresh);
  mid = (float)(fmid);
  max = (0.5 / (float)fslope) + (float)fmid;

  LoadMRISMask();

  // Go through each vertex
  for (n=0;n<mris->nvertices;n++){
    v = &mris->vertices[n];

    maskout = 0;
    if(mrismask) {
      if(fabs(MRIgetVoxVal(mrismask,n,0,0,0)) < min) maskout = 1;
      else if(fabs(v->val) < 1e-6) maskout = 1; // dont display 0s
    }
    else{
      if(fabs(v->val) < min) maskout = 1;
    }

    if(maskout){
      if(v->curv > 0) {r = 0.75; g = 0.75; b = 0.75;}
      else            {r = 0.25; g = 0.25; b = 0.25;}
    }
    else {
      if(mrismask) f = v->val/max;
      else         f = (v->val-min)/(max-min);
      dngheat(f,&r,&g,&b);
    }
    LabelColor(n, &r,&g,&b);

    colors[3*n]   = r;
    colors[3*n+1] = g;
    colors[3*n+2] = b;
  } // end loop over vertices

  return;
}
/*-------------------------------------------------------------------*/
int LabelColor(int vno, float* r, float* g, float* b)
{
  int label_index_array[LABL_MAX_LABELS];
  int num_labels_found, found_label_index;
  int label_index;

  /* if our display flag is off, do nothing. */
  if( !labl_draw_flag ) return(0);

  /* try and find a label. if found... */
  labl_find_label_by_vno (vno, 0, label_index_array,
                          LABL_MAX_LABELS, &num_labels_found);
  if (num_labels_found > 0){
    for (found_label_index = 0; found_label_index < num_labels_found;
         found_label_index++) {
      label_index = label_index_array[found_label_index];

      /* if this is the selected label and this is a border of width 1
         or 2, make it our outline color. */
      if (labl_selected_label == label_index &&
          labl_vno_is_border(labl_selected_label, vno)){
        *r = labl_outline_color[0];
        *g = labl_outline_color[1];
        *b = labl_outline_color[2];
      }
      /* else if this label is visible... */
      else if (labl_labels[label_index].visible) {
        /* color it in the given drawing style. */
        switch (labl_draw_style) {
        case LABL_STYLE_FILLED:
          *r = labl_labels[label_index].r/255.0;
          *g = labl_labels[label_index].g/255.0;
          *b = labl_labels[label_index].b/255.0;
          break;
        case LABL_STYLE_OUTLINE:
          /* if this is a border of width 1, color it the color of the
             label. */
          if (labl_vno_is_border (label_index, vno)){
            *r = labl_labels[label_index].r;
            *g = labl_labels[label_index].g;
            *b = labl_labels[label_index].b;
          }
          break;
        default:
          ;
        }
      }
    }
  }

  return(1);
}



/*!
  \fn static void fill_color_array(MRI_SURFACE *mris, float *colors)
  \brief The purpose is to set the colors vector. This is a
  3*nvertices vector of RGB triples (ie, colors[0] = r0, colors[1] =
  g0, colors[2] = b0, colors[3] = r1, etc, where rN is the color of
  vertex N. The RGB values are between 0 and 255. Simple stats
  overlays when overlayflag=1 and complexvalflag=0. The value used for
  stats is based on the sclv_current_field. The value must be >
  fthresh or < -fthresh, otherwise it gets the color of the underlying
  surface. Things that can affect surface color are: background (curv:
  red/green or gray), annotation, path, stat, debug hilite, and
  something about the surface normal.
*/
static void fill_color_array(MRI_SURFACE *mris, float *colors)
{
  int n;
  float curv;
  VERTEX *v;
#if VERTICES_PER_FACE == 4
  float intdiv0,intdiv1,intdiv2,intdiv3,frac0,frac1,frac2,frac3,nx,ny,nz;
  float cross_x[4],cross_y[4],cross_z[4];
  int crossnum;
#endif
  GLubyte  r, g, b ;
  /* begin rkt */
  GLubyte r_base, b_base, g_base;
  GLubyte r_overlay, b_overlay, g_overlay;
  float val, val2, m[4][4], nz;
  int maskout=0, mode;
  extern MRI *mrismask;
  extern double mrismaskthresh;
  float min, mid, max;

  min = (float)(fthresh);
  mid = (float)(fmid);
  max = (0.5 / (float)fslope) + (float)fmid;

  /* So we can get vertex normals in camera space and color backfaces
     red. */
  getmatrix(m);

  LoadMRISMask();

  // modes:
  //  REAL_VAL -> colscale = HEAT_SCALE, CYAN_TO_RED, BLU_GRE_RED, JUST_GRAY
  //  GREEN_RED_CURV
  //  FIELDSIGN_POS
  //  FIELDSIGN_NEG
  //  BORDER?
  //  MARKED?

  // Go through each vertex
  for (n=0;n<mris->nvertices;n++)
  {
    v = &mris->vertices[n];

    r_base = g_base = 0 ; 
    b_base = 255;

    if (mrismask)
    {
      if (fabs(MRIgetVoxVal(mrismask,n,0,0,0)) < mrismaskthresh) maskout = 1;
      else maskout = 0;
    }

    /* begin rkt */
    if (simpledrawmodeflag)
    {
      /* if surfcolor (curvature) is on and there is an overlay,
         or if the force grayscale curvature flag is on, get a
         grayscale value. if just surfcolor and curvflag are on,
         get a red/green color based on the curvature. else just
         use the solid background color. */
      if (surfcolor && curvflag && forcegraycurvatureflag)
      {
        /* grayscale curvature */
        mode = REAL_VAL;
        val2 = v->curv;
      }
      else if (surfcolor && curvflag)
      {
        /* red green curvature */
        mode = GREEN_RED_CURV;
        val2 = v->curv;
      }
      else
      {
        /* solid background color */
        mode = REAL_VAL;
        val2 = 0;
      }

      val = 0;
      if (fieldsignflag)
      {
        val = v->fsmask * v->fieldsign;
        if (v->fieldsign > 0)
        {
          if (!revphaseflag)  mode = FIELDSIGN_POS;
          else                mode = FIELDSIGN_NEG;
        }
        else
        {
          val = -val;
          if (!revphaseflag)  mode = FIELDSIGN_NEG;
          else                mode = FIELDSIGN_POS;
        }
      }

      /* This gets the RGB of the background only.  The RGB may be
         overwritten or blended later.  Note: val is NOT ignored, so
         we need to pass in the current overlay value. */
      if (overlayflag) sclv_get_value(v, sclv_current_field, &val);
      if (maskout) val = 0;
      get_color_vals(val, val2, mode, &r_base, &g_base, &b_base);

      /* save the base color for later comparison, but set our
         final rgb values to it for now. */
      r = r_base;
      g = g_base;
      b = b_base;

      /* This replaces the RGB with that of the label. The color may
         be overwritten below if a stat is above threshold thus putting
         the stat map ABOVE the label. */
      if (labels_before_overlay_flag)
        /* get any label color for this vertex. this will not apply
           any color if there is no label. */
        labl_apply_color_to_vertex (n, &r, &g, &b ); // n = vertex no

      /* if overlay flag is on... */
      if(overlayflag) {
        if(complexvalflag) {
          /* if complexvalflag is on, we have to do this
             special drawing thing. this is for compatibility
             with the two-cond stuff. assumes that val, val2,
             and stat are loaded. */
          if (surfcolor)
            get_complexval_color_vals(v->val,v->val2,v->stat,v->curv,
                                      &r_overlay,&g_overlay,&b_overlay);
          else
            get_complexval_color_vals(v->val,v->val2,v->stat,0,
                                      &r_overlay,&g_overlay,&b_overlay);

          r = r_overlay;
          g = g_overlay;
          b = b_overlay;
        }
        else
        {  // not complex
          /* get a color based on the currently selected field
             if it is above fthresh. */
          sclv_get_value(v, sclv_current_field, &val);
          if (!maskout)
          {
            /* This will blend the functional color into the
               input color. rgb are currently the background color*/
            sclv_apply_color_for_value(val, sclv_overlay_alpha,&r, &g, &b);
          }
        }
      }

      // This replaces the RGB with that of the label. This will
      // overwrite the color of a stat set above thus putting
      // the stat map BELOW the label.
      if (!labels_before_overlay_flag)
        /* get any label color for this vertex. this will not apply
           any color if there is no label. */
        labl_apply_color_to_vertex (n, &r, &g, &b );

      /* let the path code color this vertex, if it wants to. */
      path_apply_color_to_vertex (n, &r, &g, &b);

      /* Apply debug hilite color */
      ddt_get_hilite_vertex_color (n, &r, &g, &b);

      /* Transform the normal by the viewing transform and see if
         our normal is pointing away from is in camera space. If
         so, color this poly red. */
      nz = -(-m[0][2]*v->nx + m[1][2]*v->nz + m[2][2]*v->ny);
      if ( nz > 0 )
      {
        r = 255;
      }
      /* end rkt */

    }
    else
    { // not simpledrawmodeflag
      // probably will never get here. dng
      /**** msurfer: single val data on gray curvature */
      if (overlayflag && !complexvalflag)
      {
        if (v->annotation)
        {
          r = g = b = v->annotation ;
        }
        else if (fieldsignflag)
        {
          if (v->fieldsign>0.0)
          {
            if (revphaseflag)
              get_color_vals(v->fsmask*v->fieldsign,
                             v->curv,FIELDSIGN_NEG,
                             &r, &g, &b);
            else
              get_color_vals(v->fsmask*v->fieldsign,
                             v->curv,FIELDSIGN_POS,
                             &r, &g, &b);
          }
          else
          {
            if (revphaseflag)
              get_color_vals(-v->fsmask*v->fieldsign,
                             v->curv,FIELDSIGN_POS,
                             &r, &g, &b);
            else
              get_color_vals(-v->fsmask*v->fieldsign,
                             v->curv,FIELDSIGN_NEG,
                             &r, &g, &b);
          }
        }
        else if (surfcolor)
          get_color_vals(v->val,v->curv,REAL_VAL, &r, &g, &b);
        else
          get_color_vals(v->val,0.0,REAL_VAL, &r, &g, &b);
      }

      /**** msurfer: complex val data on gray curvature */
      else if (overlayflag && complexvalflag)
      {
        if (surfcolor)
          get_complexval_color_vals(v->val,v->val2,v->stat,v->curv,
                                    &r, &g, &b);
        else
          get_complexval_color_vals(v->val,v->val2,
                                    v->stat,0.0, &r, &g, &b);
      }

      /**** nsurfer: curvature, etc. red/green curv, 2d */
      else
      {
        if (v->annotation)
        {
          /* int  r, g, b ; */
          r = v->annotation & 0xff ;
          g = (v->annotation >> 8) & 0xff ;
          b = (v->annotation >> 16) & 0xff ;
        }
        else
        {
          if (surfcolor==CURVATURE_OR_SULCUS)
            curv = (avgflag)?v->curv-dipavg:v->curv;
          else
            curv = 0.0;

          if (surfcolor)
          {
            get_color_vals(0.0,curv,GREEN_RED_CURV, &r, &g, &b);
          }
          else
            get_color_vals(0.0,0.0,GREEN_RED_CURV, &r, &g, &b);

          if (v->border)
            get_color_vals(0.0,0.0,BORDER, &r, &g, &b);
        }
      }
    } // end not simpledrawmodeflag

    if (v->marked)
      get_color_vals(0.0,0.0,MARKED+v->marked-1, &r, &g, &b);



    colors[3*n]   = ((float)r)/255.0;
    colors[3*n+1] = ((float)g)/255.0;
    colors[3*n+2] = ((float)b)/255.0;
  } // end loop over vertices

}

/*!
  \fn int get_color_vals()
  \brief Appears to set the RGB based only on the curv. The "val"
  arg is immediately offset with foffset (set on the gui). 
  Colors of the overlay are controlled with sclv_apply_color_for_value().
*/
static int get_color_vals(float val, float curv, int mode,
                          GLubyte *pr, GLubyte *pg, GLubyte *pb)
{
  short r,g,b;
  float f,fr,fg,fb,tmpoffset;

  val -= foffset ;

  /* Refresh the curv min and max. This is kind of hacky since we
     could just reference the value in mris, but for some reason we
     have a separate variable for cmin and cmax, and who am I to
     judge. */
  cmin = mris->min_curv;
  cmax = mris->max_curv;

  r = g = b = 0 ;
  /* rkt: changed curv<0 to curv<cmid so that the grayscale curv
     display would use cmid */
  if (curv<cmid)  tmpoffset = cvfact*offset;
  else         tmpoffset = offset;

  if (mode==GREEN_RED_CURV)
  {
    if (curv < cmin)
    {
      b = 255 * (offset/blufact + 0.95*(1-offset/blufact)); /* yellow */
      r = g = 0 ;

    }
    else if (curv > cmax)
    {
      r = g = 255 * (offset/blufact + 0.95*(1-offset/blufact));
      /* yellow */
      b = 0 ;
    }
    else
    {
      f = tanh(cslope*(curv-cmid));
      if (f>0)
      {
        r = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
        g = 255 * (offset/blufact*(1 - fabs(f)));
      }
      else
      {
        r = 255 * (offset/blufact*(1 - fabs(f)));
        g = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
      }
      b = 255 * (offset*blufact*(1 - fabs(f)));
    }
  }

  if (mode==REAL_VAL)   /* single val positive or signed */
  {
    if (colscale==HEAT_SCALE)  /* stat */
    {
      set_stat_color(val,&fr,&fg,&fb,tmpoffset);
      r=fr;
      g=fg;
      b=fb;
    } else  /* positive */
      if (colscale==CYAN_TO_RED ||
          colscale==BLU_GRE_RED ||
          colscale==JUST_GRAY)
      {
        if (val<fthresh)
        {
          r = g = 255 * (tmpoffset/blufact);
          b =     255 * (tmpoffset*blufact);
        }
        else
        {
          if (fslope!=0)
            f = (tanh(fslope*fmid)+
                 tanh(fslope*(val-fmid)))/(2-tanh(fslope*fmid));
          else
            f = (val<0)?0:((val>1)?1:val);
          set_positive_color(f,&fr,&fg,&fb,tmpoffset);
          r=fr;
          g=fg;
          b=fb;
        }
      }
      else /* signed */
      {
        if (fabs(val)<fthresh)
        {
          r = g = 255 * (tmpoffset/blufact);
          b =     255 * (tmpoffset*blufact);
        }
        else
        {
          if (fslope!=0)
          {
            if (fmid==0)
              f = tanh(fslope*(val));
            else
            {
              if (val<0)
                f = -(tanh(fslope*fmid) + tanh(fslope*(-val-fmid)))/
                    (2-tanh(fslope*fmid));
              else
                f = (tanh(fslope*fmid) + tanh(fslope*( val-fmid)))/
                    (2-tanh(fslope*fmid));
            }
          }
          else
            f = (val<-1)?-1:((val>1)?1:val);
          if (revphaseflag)
            f = -f;
          if (truncphaseflag)
          {
            if (f<0.0) f = 0.0;
          }
          set_signed_color(f,&fr,&fg,&fb,tmpoffset);
          r=fr;
          g=fg;
          b=fb;
        }
      }
  }

  if (mode==FIELDSIGN_POS || mode==FIELDSIGN_NEG)
  {
    if (val<fthresh)
    {
      r = g = 255 * (tmpoffset/blufact);
      b =     255 * (tmpoffset*blufact);
    }
    else
    {
      f = (1.0 + tanh(fslope*(val-fmid)))/2.0;
      if (mode==FIELDSIGN_POS)
      {
        b = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
        r = g = 255* (tmpoffset*(1 - fabs(f)));
      }
      else
      {
        b = 255 * (tmpoffset*(1 - fabs(f)));
        r = g = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
      }
    }
  }

  if (mode==BORDER)  /* AMD 5/27/95 */
  {
    r = 255;
    g = 255;
    b = 0;
  }

  if (mode==MARKED)
  {
#if 0
    r = 255;
    g = 255;
    b = 255;
#else
    r = meshr;
    g = meshg;
    b = meshb;
#endif
  }

  if (mode > MARKED)
  {
    if (EVEN(mode-MARKED))
    {
      r = 255 ;
      g = 255 ;
      b = 0 ;
    }
    else
    {
      r = 0 ;
      g = 0 ;
      b = 255 ;
    }
  }

  r = (r<0)?0:(r>255)?255:r;
  g = (g<0)?0:(g>255)?255:g;
  b = (b<0)?0:(b>255)?255:b;
  *pr = (unsigned char)r ;
  *pg = (unsigned char)g ;
  *pb = (unsigned char)b ;
  return(NO_ERROR) ;
}

void
get_complexval_color_vals(float x, float y, float stat, float curv,
                          GLubyte *pr, GLubyte *pg, GLubyte *pb)
{
  short sr,sg,sb;
  float f,a,r,g,b;
  float tmpoffset,fscale;
  float a_cycles, oa = 0.0f;

  stat -= foffset ;

  if (statflag || sclv_current_field == SCLV_VALSTAT)
    f = stat;
  else
    f = sqrt(x*x+y*y);

  if (curv<0.0) tmpoffset = cvfact*offset;
  else          tmpoffset = offset;

  if (fabs(f)<fthresh)  /* trunc */
  {
    r = g = 255 * (tmpoffset/blufact);
    b =     255 * (tmpoffset*blufact);
  } else  /* use complex (or stat vals which ignore complex!!) */
  {
    if (!statflag && sclv_current_field != SCLV_VALSTAT)
    {
      if (fslope!=0)
        f = (1.0 + tanh(fslope*(f-fmid)))/2.0;
      else
        f = (f<0)?0:((f>1)?1:f);

      if (truncphaseflag)
      {
        a = atan2(y,x)/(2*M_PI);
        if (revphaseflag)
          a = -a;
        if (invphaseflag)
          a += 0.5;
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        if (a>0.5)
          f = 0;
      }
    }

    fscale = f;

    if (colscale==HEAT_SCALE || colscale==CYAN_TO_RED ||
        colscale==BLU_GRE_RED || colscale==JUST_GRAY)
    {
      if (statflag || sclv_current_field == SCLV_VALSTAT)
        set_stat_color(f,&r,&g,&b,tmpoffset);
      else
        set_positive_color(f,&r,&g,&b,tmpoffset);
    }
    else if (colscale==TWOCOND_GREEN_RED)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a -= angle_offset;       /* pos num cancels delay */
      a = fmod(a,1.0);         /* -1 to 1 */
      if (a<0) a += 1;         /* make positive */
      r = g = b = 0;
      f = sin(a*2*M_PI);
      if (f>0.0)
        r = 1;
      else
        g = 1;
      f = fabs(f)*fscale;
      r = 255 * (tmpoffset/blufact*(1-f)+f*r);
      g = 255 * (tmpoffset/blufact*(1-f)+f*g);
      b = 255 * (tmpoffset*blufact*(1-f));
    }
    else if (colscale==COLOR_WHEEL)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a_cycles = angle_cycles;
      oa = a;

      if (fmod(angle_cycles,1.0)==0.0)
        /* integral cycles (eccentricity) */
      {
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a -= 0.333333;           /* center on blue (1/3)*/
        a = a_cycles*a;          /* allow multiple */
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a = 3*a;
        r = g = b = 0;
        if (a<1.0)
        {
          r = 1-a;
          b = a;
        }
        else if (a<2.0)
        {
          b = 1-(a-1);
          g = a-1;
        }
        else
        {
          r = a-2;
          g = 1-(a-2);
        }
      }
      else /* non-integral cycles (polar angle) */
      {
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a -= 0.5;          /* center on blue (1/2) */
        a = a_cycles*a;
        r = g = b = 0;
        if (a<-0.33333)
        {
          r = 1;
        }
        else if (a<0.0)
        {
          r = 1-(a-(-0.33333))/(0.0-(-0.33333));
          b = (a-(-0.33333))/(0.0-(-0.33333));
        }
        else if (a<0.33333)
        {
          b = 1-(a)/(0.33333);
          g = (a)/(0.33333);
        }
        else
        {
          g = 1;
        }

        if (a>fadef*a_cycles/2)
        {
          f = 1-(a-fadef*a_cycles/2)/(a_cycles/2-fadef*a_cycles/2);
          r = (tmpoffset*(1-f)+f*r);
          g = (tmpoffset*(1-f)+f*g);
          b = (tmpoffset*(1-f)+f*b);
        }
        if (a<-fadef*a_cycles/2)
        {
          f = (a-(-a_cycles/2))/(a_cycles/2-fadef*a_cycles/2);
          r = (tmpoffset*(1-f)+f*r);
          g = (tmpoffset*(1-f)+f*g);
          b = (tmpoffset*(1-f)+f*b);
        }

      } /* end non-integral */
      r = (tmpoffset*(1-fscale)+fscale*r)*255;
      b = (tmpoffset*(1-fscale)+fscale*b)*255;
      g = (tmpoffset*(1-fscale)+fscale*g)*255;
    }  /* end color wheel */

    else if (colscale==RYGB_WHEEL)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a_cycles = angle_cycles;
      oa = a;

      a -= angle_offset;
      a = fmod(a,1.0);
      if (a<0) a += 1;
      a -= 0.25;           /* center on blue (1/4)*/
      a = a_cycles*a;          /* allow multiple */
      a = fmod(a,1.0);
      if (a<0) a += 1;
      a = 4*a;
      r = g = b = 0;
      if (a<1.0)
      {
        r = 1.0;
        g = a;
      }
      else if (a<2.0)
      {
        r = 1-(a-1);
        g = 1.0;
      }
      else if (a<3.0)
      {
        g = 1-(a-2);
        b = a-2;
      }
      else
      {
        r = a-3;
        b = 1-(a-3);
      }
      r = (tmpoffset*(1-fscale)+fscale*r)*255;
      b = (tmpoffset*(1-fscale)+fscale*b)*255;
      g = (tmpoffset*(1-fscale)+fscale*g)*255;
    }  /* end RYGB wheel */

    if (phasecontourflag)
    {
      if (phasecontour_min < phasecontour_max)
      {
        if (oa>phasecontour_min&&oa<phasecontour_max)
        {
          if (phasecontourmodflag)
            r = g = b = (tmpoffset*(1-fscale)+fscale*1.0)*255;
          else
            r = g = b = phasecontour_bright;
        }
      }
      else
      { /* wrap */
        if (oa>phasecontour_min||oa<phasecontour_max)
        {
          if (phasecontourmodflag)
            r = g = b = (tmpoffset*(1-fscale)+fscale*1.0)*255;
          else
            r = g = b = phasecontour_bright;
        }
      }
    }
  }
  sr = (r<0)?0:(r>255)?255:r;
  sg = (g<0)?0:(g>255)?255:g;
  sb = (b<0)?0:(b>255)?255:b;
  *pr = sr ;
  *pg = sg ;
  *pb = sb ;
}

const float DEG2RAD = 3.14159/180;

void
drawCircle(float radius)
{
  int i;
  float degInRad;

  glBegin(GL_LINE_LOOP);
  for (i=0; i < 360; i++)
  {
    degInRad = i*DEG2RAD;
    glVertex2f(cos(degInRad)*radius,sin(degInRad)*radius);
  }
  glEnd();
}

void
draw_ellipsoid_latlong(float a, float b, float c)   /* 50.0,140.0,80.0 */
{
#ifdef OPENGL
  VERTEX *v;
  int k;
  float cx,cy,cz;
  float lati,longi;

  a *= 1.01;  /* brain x,  50, lhrh     */
  b *= 1.01;  /* brain y, 140  ant/post */
  c *= 1.01;  /* brain z,  80, sup/inf  */
  cx = cy = cz = 0.0;
  for (k=0; k<mris->nvertices; k++)
  {
    v = &mris->vertices[k];
    cx += v->x;
    cy += v->y;
    cz += v->z;
  }
  cx /= mris->nvertices;
  cy /= mris->nvertices;
  cz /= mris->nvertices;
  printf("surfer: ellipsoid center: x=%f y=%f z=%f\n",cx,cy,cz);

  glPushMatrix();

  /* center:cx->lh/rh, cz->sup/inf, cy->ant/post */
  glTranslatef (-cx, -cz, cy);

  for (lati= -90.0; lati<90.0; lati+=10.0)
  {

    glPushMatrix();

    glTranslatef (0.0, 0.0, b*sin(M_PI/180.0*lati));
    glScalef (a/c, 1.0, 1.0);

    if (lati==0.0)
    {
      glLineWidth (6.0);
      glColor3ub (255, 255, 0);
    }
    else
    {
      glLineWidth (2.0);
      glColor3ub (217, 170, 0);
    }
    drawCircle (c*cos(M_PI/180.0*lati));
    glPopMatrix();
  }

  glColor3ub (190, 145, 255);
  glLineWidth (2.0);

  for (longi=0.0; longi<180.0; longi+=10.0)
  {
    glPushMatrix();
    glRotatef (90.0, 0, 1, 0 ); /* Rotate around y axis */
    glScalef (1.0, c/b, a/b);
    glRotatef (longi, 1, 0, 0); /* Rotate around the x axis */
    drawCircle (b);
    glPopMatrix();
  }

  glPopMatrix();


#else
  VERTEX *v;
  int k;
  float x,y,z;
  float cx,cy,cz;
  float lati,longi;

  a *= 1.01;  /* brain x,  50, lhrh     */
  b *= 1.01;  /* brain y, 140  ant/post */
  c *= 1.01;  /* brain z,  80, sup/inf  */
  cx = cy = cz = 0.0;
  for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    cx += v->x;
    cy += v->y;
    cz += v->z;
  }
  cx /= mris->nvertices;
  cy /= mris->nvertices;
  cz /= mris->nvertices;
  printf("surfer: ellipsoid center: x=%f y=%f z=%f\n",cx,cy,cz);

  blendfunction(BF_SA, BF_MSA);  /* no antialias despite these 4 lines */
  linesmooth(SML_ON);
  setlinestyle(0xFFFF);
  lsrepeat(1);
  pushmatrix();    /* circ() draws in brain transverse plane */
  translate(-cx,cz,cy);  /* center:cx->lh/rh, cz->sup/inf, cy->ant/post */
  for (lati= -90.0; lati<90.0; lati+=10.0)
  {
    pushmatrix();
    translate(0.0,0.0,b*sin(M_PI/180.0*lati));
    scale(a/c,1.0,1.0);
    if (lati==0.0)
    {
      linewidth(6);
      RGBcolor(255,255,0);
    }
    else
    {
      linewidth(2);
      RGBcolor(217,170,0);
    }
    circ(0.0,0.0,c*cos(M_PI/180.0*lati));
    popmatrix();
  }
  RGBcolor(190,145,255);
  linewidth(2);
  for (longi=0.0; longi<180.0; longi+=10.0)
  {
    pushmatrix();
    rot(90.0,'y');
    scale(1.0,c/b,a/b);
    rot(longi,'x');
    circ(0.0,0.0,b);
    popmatrix();
  }
  popmatrix();
  linesmooth(SML_OFF);
#endif
}

void
draw_second_surface(void)   /* for blink doublebuffer--FIX: no 2d */
{
  int k,n;
  FACE *f;
  VERTEX *v;
  float curv;

  if (flag2d)
  {
    printf("surfer: ### can't yet blink flat\n");
    PR return;
  }

  for (k=0;k<mris2->nfaces;k++)
    if (!mris2->faces[k].ripflag)
    {
      f = &mris2->faces[k];
      bgnpolygon();
      for (n=0;n<VERTICES_PER_FACE;n++)
        if (!flag2d || mris2->vertices[f->v[n]].nz>0)
        {
          v = &mris2->vertices[f->v[n]];
          if (surfcolor==CURVATURE_OR_SULCUS)
            curv = (avgflag)?v->curv-dipavg2:v->curv;
          else
            curv = 0.0;
          if (surfcolor)
            set_color(0.0,curv,GREEN_RED_CURV);
          else
            set_color(0.0,0.0,GREEN_RED_CURV);
          if (v->border)
            set_color(0.0,0.0,BORDER);
          load_brain_coords(v->nx,v->ny,v->nz,v1);
          n3f(v1);
          load_brain_coords(v->x,v->y,v->z,v1);
          v3f(v1);
        }
      endpolygon();
    }
}

void
draw_scalebar(void)
{
  float old_width;
  float v[3], tmpzf;

  pushmatrix();
  tmpzf = zf;  /* push zf */
  glGetFloatv (GL_LINE_WIDTH, &old_width);
  linewidth(SCALEBAR_WIDTH);
  RGBcolor(scalebar_bright,scalebar_bright,scalebar_bright);
  restore_zero_position();  /* zf => 1.0 */
  scale_brain(tmpzf);
  bgnline();
  v[0] = fov*sf*scalebar_xpos/tmpzf;
  v[1] = -fov*sf*scalebar_ypos/tmpzf;
  v[2] = fov*sf*9.99/tmpzf;
  v3f(v);
  v[0] -= SCALEBAR_MM;

  /* rkt: scale for surface area */
  if (mris->group_avg_surface_area > 0)
  {
    v[0] *= sqrt (mris->group_avg_surface_area / mris->total_area);
  }

  v3f(v);
  endline();
  popmatrix();
  zf = tmpzf;
  linewidth (old_width);
}

void
draw_colscalebar(void)
{
  int i, j;
  float v[3], tmpzf, stat, maxval;
  int NSEGMENTS = 10000 ;
  void *glut_font;
  float func_per_segment;
  char label[256];
  int cur_char;
  float bar_min_value;
  float bar_max_value;
  float bar_range;
  int label_at_segment[4] ;

  maxval = fmid+0.5/fslope;
  pushmatrix();
  tmpzf = zf;  /* push zf */
  restore_zero_position();  /* zf => 1.0 */
  v[0] = v[1] = 0;
  v[2] = 1.0;
  n3f(v);

  /* This is an array of segment indices at which we will draw
     labels. Init them to -1 for now. */
  for (i=0; i<4; i++)
  {
    label_at_segment[i] = -1;
  }

  /* Find the min and max value for the bar depending on our display
     flags. */
  if (truncphaseflag)
  {
    bar_min_value = 0;
    bar_max_value = maxval;
  }
  else
  {
    bar_min_value = -maxval;
    bar_max_value = maxval;
  }

  /* The full range of the bar. */
  bar_range = bar_max_value - bar_min_value;

  /* Find where we should draw our labels. */
  func_per_segment = bar_range / (float)NSEGMENTS;
  if (truncphaseflag)
  {
    label_at_segment[0] = 0;
    label_at_segment[1] = (fthresh / func_per_segment);
    label_at_segment[2] = NSEGMENTS-1;
  }
  else
  {
    label_at_segment[0] = 0;
    label_at_segment[1] = (NSEGMENTS/2) - (-fthresh / func_per_segment);
    label_at_segment[2] = (NSEGMENTS/2) + (-fthresh / func_per_segment);
    label_at_segment[3] = NSEGMENTS-1;
  }

  /* Find our font. */
  switch (colscalebar_font_size)
  {
  case 1:
    glut_font = ((void*)GLUT_BITMAP_8_BY_13);
    break;
  case 2:
    glut_font = ((void*)GLUT_BITMAP_9_BY_15);
    break;
  case 3:
    glut_font = ((void*)GLUT_BITMAP_TIMES_ROMAN_24);
    break;
  default:
    glut_font = ((void*)GLUT_BITMAP_8_BY_13);
    break;
  }

  /* For each segment... */
  for (i=0;i<NSEGMENTS;i++)
  {
    /*
      stat = fthresh+i*(maxval-fthresh)/(NSEGMENTS-1.0);
    */
    stat = bar_min_value + (float)i*bar_range/(float)(NSEGMENTS-1) +foffset;
    if (statflag || sclv_current_field == SCLV_VALSTAT)
      set_complexval_color(0.0,0.0,stat,0.0);
    else
    {
      if (complexvalflag)
        set_complexval_color(stat,0.0,0.0,0.0);
      else
        set_color(stat,0.0,REAL_VAL);
    }
    bgnquadrangle() ;
    if (colscalebarvertflag)
    {
      v[0] = fov*sf*colscalebar_xpos;
      v[1] = fov*sf*(colscalebar_ypos+colscalebar_height*(i/(NSEGMENTS-1.0)));
      v[2] = fov*sf*9.99;
      v3f(v);
      v[0] = fov*sf*(colscalebar_xpos+colscalebar_width);
      v3f(v);
      v[1] = fov*sf*(colscalebar_ypos+
                     colscalebar_height*
                     ((i+1)/(NSEGMENTS-1.0)));
      v3f(v);
      v[0] = fov*sf*colscalebar_xpos;
      v3f(v);
    }
    else
    {
      v[0] = fov*sf*(colscalebar_xpos+colscalebar_width*(i/(NSEGMENTS-1.0)));
      v[1] = fov*sf*colscalebar_ypos;
      v[2] = fov*sf*9.99;
      v3f(v);
      v[0] = fov*sf*(colscalebar_xpos+
                     colscalebar_width*
                     ((i+1)/(NSEGMENTS-1.0)));
      v3f(v);
      v[1] = fov*sf*(colscalebar_ypos+colscalebar_height);
      v3f(v);
      v[0] = fov*sf*(colscalebar_xpos+colscalebar_width*(i/(NSEGMENTS-1.0)));
      v3f(v);
    }
    endpolygon();

    /* Check our list of segments at which to draw labels, and see if
       this is one of them. */
    if (colscalebartextflag || colscalebartickflag)
      for (j=0; j <4; j++)
        if (label_at_segment[j] == i)
        {

          if (colscalebartickflag)
          {
            /* Draw an extra little line to our label. */
            glBegin (GL_LINES);
            glVertex3f (v[0], v[1], v[2]);
            if (colscalebarvertflag)
              glVertex3f (v[0]-2, v[1], v[2]);
            else
              glVertex3f (v[0], v[1]+2, v[2]);
            glEnd ();
          }

          if (colscalebartextflag)
          {
            /* Create the label string. Otherwise generate a format
               string from the decimal size of our value and print the
               value to the string. */
            if (colscalebaruselabelsflag && colscalebar_label[j])
            {
              strcpy (label, colscalebar_label[j]);
            }
            else
            {
              /* Force label to have 3 decimal places*/
              sprintf (label, "%2.3f", stat);
            }

            /* Figure out a good label position based. Here,
               strlen(label)*3.1 + strlen(label)*colscalebar_font_size*0.6 is
               a good rough estimate as to the width of the string. */
            glDisable(GL_LIGHTING); // necessary to get color setting to work:
            // http://www.devmaster.net/forums/showthread.php?p=55473
            glColor3f (1.0, 1.0, 1.0);// white
            if (colscalebarvertflag)
            {
              glRasterPos3i (v[0] - (((float)strlen(label)*3.1) +
                                     ((float)strlen(label)*
                                      (float)colscalebar_font_size*0.6))
                             - 2,
                             v[1], v[2]);
            }
            else
            {
              glRasterPos3i (v[0] - (((float)strlen(label)*3.1) +
                                     ((float)strlen(label)*
                                      (float)colscalebar_font_size*0.6))
                             / 2,
                             v[1] + 5, v[2]);
            }

            /* Draw the string. */
            for (cur_char = 0; cur_char < strlen(label); cur_char++)
            {
              glutBitmapCharacter (glut_font, label[cur_char]);
            }
            glEnable(GL_LIGHTING); // restore (was disabled to set color
          }
        }
  }

  popmatrix();
  zf = tmpzf;
}

void
set_stat_color(float f, float *rp, float *gp, float *bp, float tmpoffset)
{
  float r,g,b;
  float ftmp,c1,c2;
  float min, mid, max;
  float or, ob, og;

  r = g = b = 0.0f ;
  if (invphaseflag)
    f = -f;
  if (truncphaseflag && f<0)
    f = 0;
  if (rectphaseflag)
    f = fabs(f);

  /* rkt: same way values are calc'd in tkmedit. The main difference
     is that max is 0.5/slope + mid instead of 1/slope + mid, to make
     the linear version work better. */
  min = (float)(fthresh);
  mid = (float)(fmid);
  max = (0.5 / (float)fslope) + (float)fmid;

  if (fabs(f)>fthresh && fabs(f)<fmid)
  {
    ftmp = fabs(f);
    c1 = 1.0/(fmid-fthresh);
    if (fcurv!=1.0)
      c2 = (fmid-fthresh-fcurv*c1*SQR(fmid-fthresh))/
           ((1-fcurv)*(fmid-fthresh));
    else
      c2 = 0;
    ftmp = fcurv*c1*SQR(ftmp-fthresh)+c2*(1-fcurv)*(ftmp-fthresh)+fthresh;
    f = (f<0)?-ftmp:ftmp;
  }

  if (colscale==HEAT_SCALE)
  {
    if (f>=0)
    {
      /* rkt: changed this to make it match up with the tkmedit
         method. */
      if (sclv_opaque)
      {
        /* If opaque, don't use blending at all. Min->mid is all
           red, and mid->max gets yellower. */
        r = ((f<min) ? tmpoffset : 1.0);
        g = ((f<min) ? tmpoffset : (f<mid) ?
             0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
        b = ((f<min) ? tmpoffset : 0);
      }
      else
      {
        or = tmpoffset *
             ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
        og = tmpoffset *
             ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
        ob = tmpoffset *
             ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
        r = or +
            ((f<min) ? 0.0 : (f<mid) ? (f-min)/(mid-min) : 1.0);
        g = og +
            ((f<mid) ? 0.0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
        b = ob;
      }
    }
    else
    {
      f = -f;
      if (sclv_opaque)
      {
        b = ((f<min) ? tmpoffset : 1.0);
        g = ((f<min) ? tmpoffset : (f<mid) ?
             0.0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
        r = ((f<min) ? tmpoffset : 0);
      }
      else
      {
        or = tmpoffset *
             ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
        og = tmpoffset *
             ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
        ob = tmpoffset *
             ( (f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0 );
        b = ob +
            ((f<min) ? 0.0 : (f<mid) ? (f-min)/(mid-min) : 1.0);
        g = og +
            ((f<mid) ? 0.0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
        r = or;
      }
    }
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==BLU_GRE_RED)
  {
    /*
      if (f<0) f = -f;
      b = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
      ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):
      (f<fmid+0.25/fslope)?1-4*(f-fmid)*fslope:
      (f<fmid+0.75/fslope)?0:
      (f<fmid+1.00/fslope)?4*(f-(fmid+0.75/fslope))*fslope:1);
      g = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
      ((f<fmid)?0:(f<fmid+0.25/fslope)?4*(f-fmid)*fslope:
      (f<fmid+0.50/fslope)?1-4*(f-(fmid+0.25/fslope))*fslope:
      (f<fmid+0.75/fslope)?4*(f-(fmid+0.50/fslope))*fslope:1);
      r = tmpoffset*((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
      ((f<fmid+0.25/fslope)?0:(f<fmid+0.50/fslope)?4*(f-(fmid+0.25/fslope))
      *fslope:1);
      */
    if (f>=0)
    {
      r = tmpoffset*
          ((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):1);
      g = tmpoffset*
          ((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      b = tmpoffset*
          ((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
    }
    else
    {
      f = -f;
      b = tmpoffset*
          ((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fthresh)?0:(f<fmid)?(f-fthresh)/(fmid-fthresh):1);
      g = tmpoffset*
          ((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
          ((f<fmid)?0:(f<fmid+1.00/fslope)?1*(f-fmid)*fslope:1);
      r = tmpoffset*
          ((f<fthresh)?1:(f<fmid)?1-(f-fthresh)/(fmid-fthresh):0) +
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
  *rp = r;
  *gp = g;
  *bp = b;
}

void
set_positive_color(float f, float *rp, float *gp, float *bp, float tmpoffset)
{
  float r,g,b;

  f -= foffset ;
  r = g = b = 0 ;
  f = fabs(f);
  if (colscale==HEAT_SCALE)
  {
    r = tmpoffset + f*(1-tmpoffset);
    g = tmpoffset + ((f>0.333)?(f-0.333)/(1.0-0.333):0)*(1-tmpoffset);
    b = tmpoffset + ((f>0.666)?(f-0.666)/(1.0-0.666):0)*(1-tmpoffset);
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==CYAN_TO_RED)
  {
    b = g = tmpoffset + (1-tmpoffset)*((f<0.3)?f:0.3*(1-(f-0.3)/(1-0.3)));
    r = tmpoffset + (1-tmpoffset)*((f>0.3)?(f-0.3)/(1-0.3):0);
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==BLU_GRE_RED)
  {
    b = tmpoffset + (1-tmpoffset) *
        ((f<0.25)?f:((f<0.50)?(0.25)*(1-(f-0.25)/(1-0.25)):0));
    g = tmpoffset + (1-tmpoffset) *
        ((f<0.25)?0:((f<0.50)?2*(f-0.25):2*(0.50-0.25)*(1-(f-0.50)/(1-0.50))));
    r = tmpoffset + (1-tmpoffset)*((f<0.50)?0:(f-0.50)/(1-0.50));
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==JUST_GRAY)
  {
    r = g = b = f*255;
  }
  *rp = r;
  *gp = g;
  *bp = b;
}

void
set_signed_color(float f, float *rp, float *gp, float *bp, float tmpoffset)
{
  float r,g,b;

  if (colscale==BLUE_TO_RED_SIGNED)
  {
    b = tmpoffset + (1-tmpoffset) *
        ((f<0)?-2*f:(f>0.75)?((0.5)*(f-0.75)/(1-0.75)):0);
    r = tmpoffset + (1-tmpoffset) *
        ((f>0)?2*f:(f<-0.75)?((0.5)*(-f-0.75)/(1-0.75)):0);
    g = tmpoffset + (1-tmpoffset) *
        ((f<-0.5)?(-f-0.5)/(1.0-0.5):
         (f>0.50)?(f-0.5)/(1.0-0.5):0);
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else
    if (colscale==GREEN_TO_RED_SIGNED)
    {
      b = tmpoffset + (1-tmpoffset) *
          ((f<-0.5)?(1-(-f-0.5)/(1-0.5)):(f<0.0)?(-f/0.5):
           (f>0.5)?(1-(f-0.5)/(1-0.5)):(f/0.5));
      r = tmpoffset + (1-tmpoffset) * ((f>0)?2*f:0);
      g = tmpoffset + (1-tmpoffset) * ((f<0)?-2*f:0);
      r = r*255;
      g = g*255;
      b = b*255;
    }
    else  /*  */
    {
      if (f>0)
      {
        r = 255 *
            (tmpoffset/blufact + 0.95*(1-tmpoffset/blufact)*fabs(f));
        g = 255 * (tmpoffset/blufact*(1 - fabs(f)));
      }
      else
      {
        r = 255 * (tmpoffset/blufact*(1 - fabs(f)));
        g = 255 *
            (tmpoffset/blufact + 0.95*(1-tmpoffset/blufact)*fabs(f));
      }
      b = 255 * (tmpoffset*blufact*(1 - fabs(f)));
    }
  *rp = r;
  *gp = g;
  *bp = b;
}

void
set_color(float val, float curv, int mode)
{
  short r,g,b;
  float f,fr,fg,fb,tmpoffset;

  val -= foffset ;

  r = g = b = 0 ;
  if (curv<0)  tmpoffset = cvfact*offset;
  else         tmpoffset = offset;

  if (mode==GREEN_RED_CURV)
  {
    if (curv < cmin)
    {
      b = 255 * (offset/blufact + 0.95*(1-offset/blufact)); /* yellow */
      r = g = 0 ;

    }
    else if (curv > cmax)
    {
      r = g = 255 *
              (offset/blufact + 0.95*(1-offset/blufact)); /* yellow */
      b = 0 ;
    }
    else
    {
      f = tanh(cslope*(curv-cmid));
      if (f>0)
      {
        r = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
        g = 255 * (offset/blufact*(1 - fabs(f)));
      }
      else
      {
        r = 255 * (offset/blufact*(1 - fabs(f)));
        g = 255 * (offset/blufact + 0.95*(1-offset/blufact)*fabs(f));
      }
      b = 255 * (offset*blufact*(1 - fabs(f)));
    }
  }

  if (mode==REAL_VAL)   /* single val positive or signed */
  {
    if (colscale==HEAT_SCALE)  /* stat */
    {
      set_stat_color(val,&fr,&fg,&fb,tmpoffset);
      r=fr;
      g=fg;
      b=fb;
    } else  /* positive */
      if (colscale==CYAN_TO_RED ||
          colscale==BLU_GRE_RED ||
          colscale==JUST_GRAY)
      {
        if (val<fthresh)
        {
          r = g = 255 * (tmpoffset/blufact);
          b =     255 * (tmpoffset*blufact);
        }
        else
        {
          if (fslope!=0)
            f = (tanh(fslope*fmid)+
                 tanh(fslope*(val-fmid)))/(2-tanh(fslope*fmid));
          else
            f = (val<0)?0:((val>1)?1:val);
          set_positive_color(f,&fr,&fg,&fb,tmpoffset);
          r=fr;
          g=fg;
          b=fb;
        }
      }
      else /* signed */
      {
        if (fabs(val)<fthresh)
        {
          r = g = 255 * (tmpoffset/blufact);
          b =     255 * (tmpoffset*blufact);
        }
        else
        {
          if (fslope!=0)
          {
            if (fmid==0)
              f = tanh(fslope*(val));
            else
            {
              if (val<0)
                f = -(tanh(fslope*fmid) + tanh(fslope*(-val-fmid)))/
                    (2-tanh(fslope*fmid));
              else
                f = (tanh(fslope*fmid) + tanh(fslope*( val-fmid)))/
                    (2-tanh(fslope*fmid));
            }
          }
          else
            f = (val<-1)?-1:((val>1)?1:val);
          if (revphaseflag)
            f = -f;
          if (truncphaseflag)
          {
            if (f<0.0) f = 0.0;
          }
          set_signed_color(f,&fr,&fg,&fb,tmpoffset);
          r=fr;
          g=fg;
          b=fb;
        }
      }
  }

  if (mode==FIELDSIGN_POS || mode==FIELDSIGN_NEG)
  {
    if (val<fthresh)
    {
      r = g = 255 * (tmpoffset/blufact);
      b =     255 * (tmpoffset*blufact);
    }
    else
    {
      f = (1.0 + tanh(fslope*(val-fmid)))/2.0;
      if (mode==FIELDSIGN_POS)
      {
        b = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
        r = g = 255* (tmpoffset*(1 - fabs(f)));
      }
      else
      {
        b = 255 * (tmpoffset*(1 - fabs(f)));
        r = g = 255 * (tmpoffset + 0.95*(1-tmpoffset)*fabs(f));
      }
    }
  }

  if (mode==BORDER)  /* AMD 5/27/95 */
  {
    r = 255;
    g = 255;
    b = 0;
  }

  if (mode==MARKED)
  {
#if 0
    r = 255;
    g = 255;
    b = 255;
#else
    r = meshr;
    g = meshg;
    b = meshb;
#endif
  }
  if (mode > MARKED)
  {
    if (EVEN(mode-MARKED))
    {
      r = 255 ;
      g = 255 ;
      b = 0 ;
    }
    else
    {
      r = 0 ;
      g = 0 ;
      b = 255 ;
    }
  }

  r = (r<0)?0:(r>255)?255:r;
  g = (g<0)?0:(g>255)?255:g;
  b = (b<0)?0:(b>255)?255:b;
  RGBcolor(r,g,b);
}

void
set_complexval_color(float x, float y, float stat, float curv)
{
  short sr,sg,sb;
  float f,a,r,g,b;
  float tmpoffset,fscale;
  float a_cycles, oa = 0.0f;

  if (statflag || sclv_current_field == SCLV_VALSTAT)
    f = stat;
  else
    f = sqrt(x*x+y*y);

  f -= foffset ;

  if (curv<0.0) tmpoffset = cvfact*offset;
  else          tmpoffset = offset;

  if (fabs(f)<fthresh)  /* trunc */
  {
    r = g = 255 * (tmpoffset/blufact);
    b =     255 * (tmpoffset*blufact);
  } else  /* use complex (or stat vals which ignore complex!!) */
  {
    if (!statflag && (sclv_current_field != SCLV_VALSTAT))
    {
      if (fslope!=0)
        f = (1.0 + tanh(fslope*(f-fmid)))/2.0;
      else
        f = (f<0)?0:((f>1)?1:f);

      if (truncphaseflag)
      {
        a = atan2(y,x)/(2*M_PI);
        if (revphaseflag)
          a = -a;
        if (invphaseflag)
          a += 0.5;
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        if (a>0.5)
          f = 0;
      }
    }

    fscale = f;

    if (colscale==HEAT_SCALE || colscale==CYAN_TO_RED ||
        colscale==BLU_GRE_RED || colscale==JUST_GRAY)
    {
      if (statflag || sclv_current_field == SCLV_VALSTAT)
        set_stat_color(f,&r,&g,&b,tmpoffset);
      else
        set_positive_color(f,&r,&g,&b,tmpoffset);
    }
    else if (colscale==TWOCOND_GREEN_RED)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a -= angle_offset;       /* pos num cancels delay */
      a = fmod(a,1.0);         /* -1 to 1 */
      if (a<0) a += 1;         /* make positive */
      r = g = b = 0;
      f = sin(a*2*M_PI);
      if (f>0.0)
        r = 1;
      else
        g = 1;
      f = fabs(f)*fscale;
      r = 255 * (tmpoffset/blufact*(1-f)+f*r);
      g = 255 * (tmpoffset/blufact*(1-f)+f*g);
      b = 255 * (tmpoffset*blufact*(1-f));
    }
    else if (colscale==COLOR_WHEEL)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a_cycles = angle_cycles;
      oa = a;

      if (fmod(angle_cycles,1.0)==0.0)
        /* integral cycles (eccentricity) */
      {
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a -= 0.333333;           /* center on blue (1/3)*/
        a = a_cycles*a;          /* allow multiple */
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a = 3*a;
        r = g = b = 0;
        if (a<1.0)
        {
          r = 1-a;
          b = a;
        }
        else if (a<2.0)
        {
          b = 1-(a-1);
          g = a-1;
        }
        else
        {
          r = a-2;
          g = 1-(a-2);
        }
      }
      else /* non-integral cycles (polar angle) */
      {
        a -= angle_offset;
        a = fmod(a,1.0);
        if (a<0) a += 1;
        a -= 0.5;          /* center on blue (1/2) */
        a = a_cycles*a;
        r = g = b = 0;
        if (a<-0.33333)
        {
          r = 1;
        }
        else if (a<0.0)
        {
          r = 1-(a-(-0.33333))/(0.0-(-0.33333));
          b = (a-(-0.33333))/(0.0-(-0.33333));
        }
        else if (a<0.33333)
        {
          b = 1-(a)/(0.33333);
          g = (a)/(0.33333);
        }
        else
        {
          g = 1;
        }

        if (a>fadef*a_cycles/2)
        {
          f = 1-(a-fadef*a_cycles/2)/(a_cycles/2-fadef*a_cycles/2);
          r = (tmpoffset*(1-f)+f*r);
          g = (tmpoffset*(1-f)+f*g);
          b = (tmpoffset*(1-f)+f*b);
        }
        if (a<-fadef*a_cycles/2)
        {
          f = (a-(-a_cycles/2))/(a_cycles/2-fadef*a_cycles/2);
          r = (tmpoffset*(1-f)+f*r);
          g = (tmpoffset*(1-f)+f*g);
          b = (tmpoffset*(1-f)+f*b);
        }

      } /* end non-integral */
      r = (tmpoffset*(1-fscale)+fscale*r)*255;
      b = (tmpoffset*(1-fscale)+fscale*b)*255;
      g = (tmpoffset*(1-fscale)+fscale*g)*255;
    }  /* end color wheel */

    else if (colscale==RYGB_WHEEL)
    {
      a = atan2(y,x)/(2*M_PI);
      if (revphaseflag)
        a = -a;
      if (invphaseflag)
        a += 0.5;
      a_cycles = angle_cycles;
      oa = a;

      a -= angle_offset;
      a = fmod(a,1.0);
      if (a<0) a += 1;
      a -= 0.25;           /* center on blue (1/4)*/
      a = a_cycles*a;          /* allow multiple */
      a = fmod(a,1.0);
      if (a<0) a += 1;
      a = 4*a;
      r = g = b = 0;
      if (a<1.0)
      {
        r = 1.0;
        g = a;
      }
      else if (a<2.0)
      {
        r = 1-(a-1);
        g = 1.0;
      }
      else if (a<3.0)
      {
        g = 1-(a-2);
        b = a-2;
      }
      else
      {
        r = a-3;
        b = 1-(a-3);
      }
      r = (tmpoffset*(1-fscale)+fscale*r)*255;
      b = (tmpoffset*(1-fscale)+fscale*b)*255;
      g = (tmpoffset*(1-fscale)+fscale*g)*255;
    }  /* end RYGB wheel */

    if (phasecontourflag)
    {
      if (phasecontour_min < phasecontour_max)
      {
        if (oa>phasecontour_min&&oa<phasecontour_max)
        {
          if (phasecontourmodflag)
            r = g = b = (tmpoffset*(1-fscale)+fscale*1.0)*255;
          else
            r = g = b = phasecontour_bright;
        }
      }
      else
      { /* wrap */
        if (oa>phasecontour_min||oa<phasecontour_max)
        {
          if (phasecontourmodflag)
            r = g = b = (tmpoffset*(1-fscale)+fscale*1.0)*255;
          else
            r = g = b = phasecontour_bright;
        }
      }
    }
  }
  sr = (r<0)?0:(r>255)?255:r;
  sg = (g<0)?0:(g>255)?255:g;
  sb = (b<0)?0:(b>255)?255:b;
  RGBcolor(sr,sg,sb);
}

void
draw_spokes(int option)
{
  float /*r0=0.1*/r0=0.07,r1=1.0,th0=0,th1=2*M_PI;
  int nr=100,nth=100;
  float scf=50.0;
  int thi,ri;
  double cosa,sina,cosb,sinb,ra,rb,tha,thb,a;
  float v0[3],v1[3],v2[3],v3[3],v4[3];

  RGBcolor(0,0,0);
  czclear(BACKGROUND,getgconfig(GC_ZMAX));
  v0[2] = v1[2] = v2[2] = v3[2] = 0;
  v4[0] = v4[1] = 0;
  v4[2] = 1;
  a = pow(r1/r0,1.0/(r1-r0)); /* marty */
  for (thi=0;thi<nth-1;thi++)
  {
    tha = th0+thi*(th1-th0)/(nth-1);
    thb = th0+(thi+1)*(th1-th0)/(nth-1);
    cosa = cos(tha);
    sina = sin(tha);
    cosb = cos(thb);
    sinb = sin(thb);
    for (ri=1;ri<nr+1;ri++)
    {
      ra = r0+(ri-1)*(r1-r0)/(nr);
      rb = r0+(ri)*(r1-r0)/(nr);
      ra = (ra<r0)?r0:(ra>r1)?r1:ra;
      rb = (rb<r0)?r0:(rb>r1)?r1:rb;
      if (option==RADIUS) /* marty: just change size */
      {
        ra = r0*pow(a,ra-r0);
        rb = r0*pow(a,rb-r0);
      }
      v0[0] = cosa*ra*scf;
      v0[1] = sina*ra*scf;
      v1[0] = cosa*rb*scf;
      v1[1] = sina*rb*scf;
      v2[0] = cosb*rb*scf;
      v2[1] = sinb*rb*scf;
      v3[0] = cosb*ra*scf;
      v3[1] = sinb*ra*scf;
      n3f(v4);
      bgnquadrangle();
      set_vertex_color((ri-1.0)/(nr),(thi-0.0)/(nth-1.0),option);
      v3f(v0);
      set_vertex_color((ri-0.0)/(nr),(thi-0.0)/(nth-1.0),option);
      v3f(v1);
      set_vertex_color((ri-0.0)/(nr),(thi+1.0)/(nth-1.0),option);
      v3f(v2);
      set_vertex_color((ri-1.0)/(nr),(thi+1.0)/(nth-1.0),option);
      v3f(v3);
      endpolygon();
    }
  }
}

void
set_vertex_color(float r, float th, int option)
{
  if (option==RADIUS) set_color_wheel(r,0.0,angle_cycles,0,TRUE,1.0);
  if (option==THETA)  set_color_wheel(th,0.0,angle_cycles,0,FALSE,1.0);
}

/* msurferCMF.c (msurferMARTY.c hack) */
void
set_color_wheel(float a, float a_offset, float a_cycles, int mode, int logmode,
                float fscale)
{
  float r,g,b,f,fadef=0.7;
  short sr,sg,sb;
  float maxrad=1.0,minrad=0.05,afact,offset=0.0f;

  if (mode==-1) offset=cvfact*offset;
  else          offset=offset;
  afact = pow(maxrad/minrad,1.0/(maxrad-minrad));
  if (fmod(a_cycles,1.0)==0)
  {
    a += 0.5;
    a += a_offset;
    a = fmod(a,1.0);
    if (logmode)
      a += 0.5;  /* marty: color center of half-duty ring */
    /*a = (log(a/minrad)/log(afact)+minrad-minrad)/(1-minrad); as below */
    /*a = minrad*pow(afact,a-minrad); */
    a = fmod(a,1.0);  /* marty: OK if not logmode?? */
    a = 3*a;
    r = g = b = 0;
    if (a<1.0)
    {
      r = 1-a;
      b = a;
    }
    else
      if (a<2.0)
      {
        b = 1-(a-1);
        g = a-1;
      }
      else
      {
        r = a-2;
        g = 1-(a-2);
      }
  }
  else  /* with angle offset */
  {
    a += 0.5;
    a += a_offset;
    a = fmod(a,1.0);
    a -= 0.5;
    a = a_cycles*a;
    /*a = (a<-1)?-1:(a>1)?1:a;*/
    r = g = b = 0;
    if (a<-0.5)
    {
      r = 1;
    }
    else
      if (a<0.0)
      {
        r = 1-(a-(-0.5))/(0.0-(-0.5));
        b = (a-(-0.5))/(0.0-(-0.5));
      }
      else
        if (a<0.5)
        {
          b = 1-(a)/(0.5);
          g = (a)/(0.5);
        }
        else
        {
          g = 1;
        }
    if (a>fadef*a_cycles/2)
    {
      f = 1-(a-fadef*a_cycles/2)/(a_cycles/2-fadef*a_cycles/2);
      /*
        r *= f;
        g *= f;
        b *= f;
      */
      r = (offset*(1-f)+f*r);
      g = (offset*(1-f)+f*g);
      b = (offset*(1-f)+f*b);
    }
    if (a<-fadef*a_cycles/2)
    {
      f = (a-(-a_cycles/2))/(a_cycles/2-fadef*a_cycles/2);
      /*r *= f;g *= f;b *= f;*/
      r = (offset*(1-f)+f*r);
      g = (offset*(1-f)+f*g);
      b = (offset*(1-f)+f*b);
    }
  }
  /*r = fscale*r*255;b = fscale*b*255;g = fscale*g*255;*/
  r = (offset*(1-fscale)+fscale*r)*255;
  b = (offset*(1-fscale)+fscale*b)*255;
  g = (offset*(1-fscale)+fscale*g)*255;
  sr = (r<0)?0:(r>255)?255:r;
  sg = (g<0)?0:(g>255)?255:g;
  sb = (b<0)?0:(b>255)?255:b;
  RGBcolor(sr,sg,sb);
}

void
restore_ripflags(int mode)
{
  int k,h;

  printf("restore_ripflags(%d)\n",mode);
  if (mode==1)
    for (k=0;k<mris->nvertices;k++)
    {
      h = mris->vertices[k].oripflag;
      mris->vertices[k].oripflag = mris->vertices[k].ripflag;
      mris->vertices[k].ripflag = h;
    }
  else if (mode==2)
    for (k=0;k<mris->nvertices;k++)
    {
      mris->vertices[k].oripflag = mris->vertices[k].ripflag;
      mris->vertices[k].ripflag = mris->vertices[k].origripflag;
    }
  rip_faces();
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].marked = FALSE;
  clear_vertex_marks();
  vertex_array_dirty = 1;
}

void
dilate_ripped(void)
{
  int    vno, n, nripped ;
  VERTEX *v, *vn ;

  MRISclearD(mris) ;
  nripped = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag == 0)
      continue ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->ripflag == 0)
        vn->d = 1.0 ;
    }
    if (v->border && v->d < 0.5)
    {
      v->d = 1 ;
    }
  }
  undo_begin_action (UNDO_CUT);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->d > 0.5)
    {
      v->d = 0 ;
      set_vertex_rip(vno, TRUE, TRUE) ;
      nripped++ ;
    }
  }
  rip_faces() ;
  vertex_array_dirty = 1;
  undo_finish_action ();
  /* might need to reuild vertex positions */
  vset_set_current_set(vset_current_set) ;
  printf("%d vertices ripped\n", nripped) ;
}

void 
rip_unmarked_vertices(void)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked == 0)
      v->ripflag = 1 ;
  }
  rip_faces() ;
  surface_compiled = 0;
  vertex_array_dirty = 1;

  vset_set_current_set (vset_current_set);
}

void
floodfill_marked_patch(int filltype)
{
  /* Compatibility function. */
  switch (filltype)
  {
  case RIPFILL:
    rip_all_vertices_except_contiguous_upripped ();
    break;
  case CURVFILL:
    mark_contiguous_vertices_with_similar_curvature ();
    break;
  case STATFILL:
    mark_contiguous_vertices_over_thresh ();
    break;
  }

#if 0
  FILL_PARAMETERS params;
  VERTEX *v, *vn;
  int i,k,m,filled,totfilled=0,kmin,kmax,kstep;
  /* begin rkt */
  float value = 0;
  /* end rkt */

  if (nmarked==0) return;

  params.dont_cross_path = 0;
  params.dont_cross_label = 0;
  params.dont_cross_cmid = 0;
  params.dont_cross_fthresh = 0;
  params.dont_fill_unlabeled = 0;
  params.use_multiple_seeds = 0;
  params.action = FILL_NO_ACTION_JUST_MARK;
  fill_flood_from_seed (marked[0], &params);

  if (filltype == CURVFILL)
  {
    area = LabelAlloc(mris->nvertices, pname, lfname) ;
    strncpy( area->subject_name, pname, 100 );
    fprintf(stderr, "subject_name=%s, lfname=%s\n", pname,lfname) ;
    LabelCurvFill(area, marked, nmarked, mris->nvertices, mris) ;
    LabelMark(area, mris) ;
    /*    redraw() ;*/
    return ;
  }
  else if (filltype == RIPFILL)
  {

    /* begin rkt */

    undo_begin_action (UNDO_CUT);

    for (k=0;k<mris->nvertices;k++)
    {
      v = &mris->vertices[k] ;
      if (!v->marked)
        set_vertex_rip (k, TRUE, TRUE);
      else
        set_vertex_rip (k, FALSE, TRUE);
    }
    rip_faces() ;

    /* begin rkt */
    undo_finish_action ();

    clear_vertex_marks();

    /* might need to reuild vertex positions */
    vset_set_current_set(vset_current_set) ;

    /* end rkt */

    redraw() ;
    return ;
  }


  filled = nmarked;
  i = 0;
  while (filled>0)
  {
    filled = 0;
    if ((i%2)==0)
    {
      kmin = 0;
      kmax = mris->nvertices-1;
      kstep = 1;
    }
    else
    {
      kmin = mris->nvertices-1;
      kmax = 0;
      kstep = -1;
    }
    for (k=kmin;k!=kmax;k+=kstep)
    {
      if (mris->vertices[k].marked)   /* look for an unmarked neighbor */
      {
        v = &mris->vertices[k];
        for (m=0;m<v->vnum;m++)
        {
          vn = &mris->vertices[v->v[m]] ;
          /* begin rkt */
#if 0
          if (!vn->marked&&
              !vn->border&&
              (filltype!=STATFILL || fabs(vn->stat)>=fthresh))
#else
          sclv_get_value(vn,sclv_current_field,&value);
          if (!vn->marked&&
              !vn->border&&
              (filltype!=STATFILL || fabs(value)>=fthresh))
#endif
            /* end rkt */
          {
            vn->marked = TRUE;
            filled++;
            totfilled++;
          }
        }
      }
    }
    printf("surfer: %d: filled = %d, total = %d (of %d)\n",
           i,filled,totfilled,mris->nvertices);
    draw_surface();
    i++;
  }
  area = LabelFromMarkedSurface(mris) ;
  PR
#endif
}

int
mark_contiguous_vertices_over_thresh ()
{
  FILL_PARAMETERS params;
  int error;

  if (nmarked==0)
  {
    printf ("surfer: need a seed point, please mark a vertex\n");
    return (ERROR_BADPARM);
  }

  params.dont_cross_path     = 0;
  params.dont_cross_label    = 0;
  params.dont_cross_cmid     = 0;
  params.dont_cross_fthresh  = 1;
  params.dont_fill_unlabeled = 0;
  params.use_multiple_seeds  = 0;
  params.action = FILL_NO_ACTION_JUST_MARK;

  error = fill_flood_from_seed (marked[0], &params);
  if (NO_ERROR != error)
    printf ("surfer: Error marking.\n");

  surface_compiled = 0 ;
  vertex_array_dirty = 1;

  return (NO_ERROR);
}

int
mark_contiguous_vertices_with_similar_curvature ()
{
  FILL_PARAMETERS params;
  int error;

  if (nmarked==0)
  {
    printf ("surfer: need a seed point, please mark a vertex\n");
    return (ERROR_BADPARM);
  }

  params.dont_cross_path     = 0;
  params.dont_cross_label    = 0;
  params.dont_cross_cmid     = 1;
  params.dont_cross_fthresh  = 0;
  params.dont_fill_unlabeled = 0;
  params.use_multiple_seeds  = 0;
  params.action = FILL_NO_ACTION_JUST_MARK;

  error = fill_flood_from_seed (marked[0], &params);
  if (NO_ERROR != error)
    printf ("surfer: Error marking.\n");

  surface_compiled = 0 ;
  vertex_array_dirty = 1;

  return (NO_ERROR);
}

int
rip_all_vertices_except_contiguous_upripped ()
{
  FILL_PARAMETERS params;
  int error;
  int vno;
  int count;

  if (nmarked==0)
  {
    printf ("surfer: need a seed point, please mark a vertex\n");
    return (ERROR_BADPARM);
  }

  params.dont_cross_path     = 0;
  params.dont_cross_label    = 0;
  params.dont_cross_cmid     = 0;
  params.dont_cross_fthresh  = 0;
  params.dont_fill_unlabeled = 0;
  params.use_multiple_seeds  = 0;
  params.action = FILL_NO_ACTION_JUST_MARK;

  error = fill_flood_from_seed (marked[0], &params);
  if (NO_ERROR != error)
    printf ("surfer: Error ripping.\n");

  undo_begin_action (UNDO_CUT);

  count = 0;
  for (vno = 0; vno < mris->nvertices; vno++)
  {
    if (!mris->vertices[vno].marked)
    {
      count++;
      set_vertex_rip (vno, TRUE, TRUE);
    }
  }

  undo_finish_action ();

  clear_all_vertex_marks ();

  printf ("surfer: ripped %d vertices\n", count);

  rip_faces();
  surface_compiled = 0;
  vertex_array_dirty = 1;

  vset_set_current_set (vset_current_set);

  return (NO_ERROR);
}

void
clear_ripflags(void)
{
  int k;

  /* begin rkt */
  undo_begin_action (UNDO_CUT);

  for (k=0;k<mris->nvertices;k++)
    /* mris->vertices[k].ripflag = FALSE;*/
    set_vertex_rip (k, FALSE, TRUE);

  undo_finish_action();
  /* end rkt */

  rip_faces();
}

void
cut_marked_vertices(int closedcurveflag)
{
  int i,j,k,m,mmin;
  float d,d1,d2,dx,dy,dz,x1,y1,z1,x2,y2,z2,dmin;
  VERTEX *v1,*v2;

  /* begin rkt */
  printf ("******************** cut_marked_vertices called\n");
  undo_begin_action (UNDO_CUT);
  /* end rkt */

  for (i=0;i<((closedcurveflag)?nmarked:nmarked-1);i++)
  {
    j = (i==nmarked-1)?0:i+1;
    v2 = &mris->vertices[marked[j]];
    x2 = v2->x;
    y2 = v2->y;
    z2 = v2->z;
    k = marked[i];
    while (k!=marked[j])
    {
      v1 = &mris->vertices[k];
      v1->border = 1;
      /* begin rkt */
      /*v1->ripflag = 1;*/
      set_vertex_rip (k, TRUE, TRUE);
      /* end rkt */
      x1 = v1->x;
      y1 = v1->y;
      z1 = v1->z;
      dmin = 100000;
      mmin = 0;
      for (m=0;m<v1->vnum;m++)
        if (!mris->vertices[v1->v[m]].border||v1->v[m]==marked[0])
        {
          dx = mris->vertices[v1->v[m]].x - x1;
          dy = mris->vertices[v1->v[m]].y - y1;
          dz = mris->vertices[v1->v[m]].z - z1;
          d1 = sqrt(dx*dx+dy*dy+dz*dz);
          dx = mris->vertices[v1->v[m]].x - x2;
          dy = mris->vertices[v1->v[m]].y - y2;
          dz = mris->vertices[v1->v[m]].z - z2;
          d2 = sqrt(dx*dx+dy*dy+dz*dz);
          if ((d=d1+d2)<dmin)
          {
            dmin = d;
            mmin = m;
          }
          printf("surfer:       k=%d, m=%d, v1->v[m]=%d, "
                 "d1=%f, d2=%f, d=%f\n",
                 k,m,v1->v[m],d1,d2,d);
        }
      if (mmin==100000)
      {
        printf("surfer: failed (infinite loop)\n");
        PR return;
      }
      k = v1->v[mmin];
      printf("surfer: best: mmin=%d, k=%d, dmin=%f\n",mmin,k,dmin);
    }
    /* begin rkt */
    /*v2->border = 1;*/
    /*v2->ripflag = 1;*/
    set_vertex_rip (marked[j], TRUE, TRUE);
    /* end rkt */
  }
  PR
  rip_faces();

  /* begin rkt */
  undo_begin_action (UNDO_CUT);
  /* end rkt */

  clear_vertex_marks();
}

void
cut_plane(void)
{
  int k;
  float dx0,dy0,dz0,dx1,dy1,dz1,nx,ny,nz,poffset,sign;
  VERTEX *v,*v0,*v1,*v2,*v3;

  /* begin rkt */
  undo_begin_action (UNDO_CUT);
  /* end rkt */

  if (nmarked!=4)
  {
    printf("surfer: needs 4 marked vertices\n");
    PR return;
  }
  /* begin rkt */
#if 0
  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].oripflag = mris->vertices[k].ripflag;
#endif
  /* end rkt */
  v0 = &mris->vertices[marked[0]];
  v1 = &mris->vertices[marked[1]];
  v2 = &mris->vertices[marked[2]];
  v3 = &mris->vertices[marked[3]];
  dx0 = v1->x-v0->x;
  dy0 = v1->y-v0->y;
  dz0 = v1->z-v0->z;
  dx1 = v2->x-v1->x;
  dy1 = v2->y-v1->y;
  dz1 = v2->z-v1->z;
  nx = -dy1*dz0 + dy0*dz1;
  ny = dx1*dz0 - dx0*dz1;
  nz = -dx1*dy0 + dx0*dy1;
  poffset = nx*v0->x+ny*v0->y+nz*v0->z;
  sign = (nx*v3->x+ny*v3->y+nz*v3->z)-poffset;
  printf("surfer: poffset = %f, sign = %f, n = {%f,%f,%f}\n",
         poffset,sign,nx,ny,nz);
  PR for (k=0;k<mris->nvertices;k++)
  {
    v = &mris->vertices[k];
    if ((((nx*v->x+ny*v->y+nz*v->z)-poffset)/sign)<0)
      /* begin rkt */
      /*v->ripflag = TRUE;*/
      set_vertex_rip (k, TRUE, TRUE);
    /* end rkt */
  }
  rip_faces();
  clear_vertex_marks();
  vertex_array_dirty = 1;

  /* begin rkt */
  undo_finish_action();
  /* end rkt */
}

void
cut_vertex(void)
{
  /* begin rkt */
# if 0
  int k;

  undo_begin_action (UNDO_CUT);

  for (k=0;k<mris->nvertices;k++)
    mris->vertices[k].oripflag = mris->vertices[k].ripflag;
  mris->vertices[marked[nmarked-1]].ripflag = TRUE;
#endif
  set_vertex_rip (marked[nmarked-1], TRUE, TRUE);
  rip_faces();
  clear_vertex_marks();

  undo_finish_action();
  /* end rkt */
}

void
cut_line(int closedcurveflag)
{
  int vno;
  int* path;
  int path_length;

  if (nmarked<2)
  {
    printf("surfer: needs at least 2 marked vertices\n");
    PR return;
  }

  if (closedcurveflag)
    close_marked_vertices();

  path = (int*) calloc (mris->nvertices, sizeof(int));

  find_path (marked, nmarked, "making cut", mris->nvertices,
             path, &path_length);

  undo_begin_action (UNDO_CUT);

  for (vno = 0; vno < path_length; vno++)
  {
    set_vertex_rip (path[vno], TRUE, TRUE);
  }

  rip_faces();
  clear_vertex_marks();
  vertex_array_dirty = 1;

  undo_finish_action();
}

/* begin rkt */

void set_face_rip(int fno, int rip, int undoable)
{
  FACE *f;

  f = &mris->faces[fno];

  /* add the action to undo list if the values are different */
  if (undoable && rip!=f->ripflag)
    undo_new_action_cut (UNDO_CUT_FACE, fno, f->ripflag);

  f->ripflag = rip;
}

void set_vertex_rip(int vno, int rip, int undoable)
{
  VERTEX *v;

  v = &mris->vertices[vno];

  /* add the action to undo list if the values are different */
  if (undoable && rip!=v->ripflag)
    undo_new_action_cut (UNDO_CUT_VERTEX, vno, v->ripflag);

  v->ripflag = rip;

}

void close_marked_vertices ()
{
  /* just add a copy of the first vertex to the end of the list. */
  marked[nmarked] = marked[0];
  nmarked++;
}

void draw_vertex_hilite (int vno)
{
  VERTEX* v;
  VERTEX* vn;
  int neighbor_index;

  if (vno < 0 || vno >= mris->nvertices)
    return;

  v = &(mris->vertices[vno]);

  for (neighbor_index = 0; neighbor_index < v->vnum; neighbor_index++)
  {
    vn = &(mris->vertices[v->v[neighbor_index]]);

    glLineWidth (CURSOR_LINE_PIX_WIDTH);

    glBegin (GL_LINES);

    glNormal3f (-(vn->nx), vn->nz, vn->ny);

    glVertex3f (-(vn->x + cup * vn->nx),
                vn->z + cup * vn->nz,
                vn->y + cup * vn->ny);

    glNormal3f (-(v->nx), v->nz, v->ny);

    glVertex3f (-(v->x + cup * vn->nx),
                v->z + cup * vn->nz,
                v->y + cup * vn->ny);
    glEnd ();
  }

  glFlush ();
}

void draw_marked_vertices ()
{
  int marked_index;

  /* draw all our marked verts in white. */
  RGBcolor (255, 255, 255);
  for (marked_index = 0; marked_index < nmarked-1; marked_index++)
    draw_vertex_hilite (marked[marked_index]);
}

/* end rkt */

void
draw_vector(char *fname)
{
  float   dx, dy, len ;
  float   v[3] ;
  FILE    *fp ;
  static int ncalls = 0 ;

  fp = fopen(fname, "r") ;

  if (!fp)
  {
    printf("could not open file %s.\n", fname) ;
    return ;
  }

  fscanf(fp, "%f %f", &dx, &dy) ;
  len = sqrt(dx*dx + dy*dy) ;
  dx /= len ;
  dy /= len ;

  if ((ncalls++ % 2) == 0)
    RGBcolor(0, 255, 255) ;
  else
    RGBcolor(255, 255, 0) ;
  linewidth(CURSOR_LINE_PIX_WIDTH);
  bgnline() ;
  load_brain_coords(0,0,1,v);
  n3f(v);

  load_brain_coords(0,0,0,v);
  v3f(v);
  load_brain_coords(dx*10,dy*10,0,v);
  v3f(v);

  endline() ;
  fclose(fp) ;
}

void
put_retinotopy_stats_in_vals(void)
{
  int      vno ;
  VERTEX   *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = hypot(v->val, v->val2) ;
  }
}

void
plot_marked(char *fname)
{
  FILE   *fp ;
  int    vno, total ;
  VERTEX *v ;

  fprintf(stderr, "generating data file %s with entries:\n", fname) ;
  fprintf(stderr, "vno x y z distance curv val val2 stat amp "
          "deg normalized_deg radians\n") ;
  fp = fopen(fname, "w") ;
  if (!fp)
  {
    fprintf(stderr, "### could not open file %s\n", fname) ;
    return ;
  }

  for (vno = total = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag || !v->marked)
      continue ;
    fprintf(fp, "%d %f %f %f %f %f %f %f %f %f %f %f %f\n",
            vno, v->x, v->y, v->z, 0.0, v->curv, v->val, v->val2, v->stat,
            hypot(v->val,v->val2), (float)(atan2(v->val2,v->val)*180/M_PI),
            (float)(atan2(v->val2,v->val)/(2*M_PI)),
            (float)(atan2(v->val2,v->val))) ;

  }
  fclose(fp) ;
}

static int find_best_path(MRI_SURFACE *mris, int start_vno, int end_vno,
                          int *path_indices,
                          double target_curv, double l_curv,
                          double l_len) ;

void
draw_fundus(int bdry_index)
{
  PATH_PATH  *borig, *bnew, *bold ;
  double         last_sse, sse, curv_error, initial_length, total_length,
  l_curv, l_len, min_curv_error, target_curv ;
  int            i, vno, n, min_n, *indices/*, new_index, old_index*/, \
  v_in_path, start_vno, end_vno, niter ;
  VERTEX         *v = NULL, *v_last, *vn ;

  MRISclearMarks(mris) ;
  target_curv = mris->max_curv ;
  if (bdry_index == PATH_NONE_SELECTED)
  {
    printf("### no boundary selected\n") ;
    return ;
  }
  borig = &path_paths[bdry_index] ;
  v_last = NULL ;
  l_curv = 1.0 ;
  l_len = 0.01 ;
  initial_length = path_length(borig, mris) ;

  bnew = (PATH_PATH *)calloc(1, sizeof(PATH_PATH)) ;
  bnew->vertices = (int *)calloc(mris->nvertices, sizeof(int)) ;
  bold = path_copy(borig, NULL) ;
  indices = (int *)calloc(mris->nvertices, sizeof(int)) ;
  if (indices == NULL || bnew->vertices == NULL)
    ErrorExit(ERROR_NOMEMORY, "%s: could not allocate index array", Progname) ;
  path_copy(bold, bnew) ;
  niter = 0 ;
  do
  {
    /*    path_set_marks(bold, mris, 1) ;*/

    /* compute sse of current line */
    last_sse = path_sse(bold, mris, target_curv, l_curv, l_len) ;

    bnew->num_vertices = 0 ;

    /* first move endpoints in direction of max curvature */
    path_set_marks(bold, mris, 1) ;
    start_vno = bold->vertices[0] ;
    min_n = -1 ;
    min_curv_error = v->curv-target_curv ;
    min_curv_error *= min_curv_error ;
    v = &mris->vertices[start_vno] ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked)
        continue ;  /* already in path */
      curv_error = vn->curv - target_curv ;
      curv_error *= curv_error ;
      if (min_n < 0 || curv_error < min_curv_error)
      {
        min_curv_error = curv_error ;
        min_n = n ;
      }
    }
    if (min_n >= 0)
    {
      printf("modifying starting vertex to be %d (was %d)...\n",
             v->v[min_n],start_vno) ;
      start_vno = v->v[min_n] ;
    }
    else   /* try searching nbrs of nbr of endpoint */
    {
      min_n = -1 ;
      v = &mris->vertices[bold->vertices[0]] ;
      min_curv_error = v->curv-target_curv ;
      min_curv_error *= min_curv_error ;
      v = &mris->vertices[bold->vertices[1]] ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked)
          continue ;  /* already in path */
        curv_error = vn->curv - target_curv ;
        curv_error *= curv_error ;
        if (min_n < 0 || curv_error < min_curv_error)
        {
          min_curv_error = curv_error ;
          min_n = n ;
        }
      }
      if (min_n >= 0)
      {
        printf("modifying starting vertex to be %d (was %d)...\n",
               v->v[min_n],start_vno) ;
        start_vno = v->v[min_n] ;
      }
    }

    end_vno = bold->vertices[bold->num_vertices-1];
    v = &mris->vertices[end_vno] ;
    min_n = -1 ;
    min_curv_error = v->curv-target_curv ;
    min_curv_error *= min_curv_error ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vn = &mris->vertices[v->v[n]] ;
      if (vn->marked)
        continue ;  /* already in path */
      curv_error = vn->curv - target_curv ;
      curv_error *= curv_error ;
      if (min_n < 0 || curv_error < min_curv_error)
      {
        min_curv_error = curv_error ;
        min_n = n ;
      }
    }
    if (min_n >= 0)
    {
      printf("modifying end point to be %d (was %d)...\n",
             v->v[min_n], end_vno) ;
      end_vno = v->v[min_n] ;
    }
    else   /* try searching nbrs of nbr of endpoint */
    {
      min_n = -1 ;
      v = &mris->vertices[bold->vertices[bold->num_vertices-1]] ;
      min_curv_error = v->curv-target_curv ;
      min_curv_error *= min_curv_error ;
      v = &mris->vertices[bold->vertices[bold->num_vertices-2]] ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vn = &mris->vertices[v->v[n]] ;
        if (vn->marked)
          continue ;  /* already in path */
        curv_error = vn->curv - target_curv ;
        curv_error *= curv_error ;
        if (min_n < 0 || curv_error < min_curv_error)
        {
          min_curv_error = curv_error ;
          min_n = n ;
        }
      }
      if (min_n >= 0)
      {
        printf("modifying ending vertex to be %d (was %d)...\n",
               v->v[min_n],end_vno) ;
        end_vno = v->v[min_n] ;
      }
    }


    /* consider removing this vertex and putting a new path in */
    v_in_path = find_best_path(mris, start_vno, end_vno, indices,
                               target_curv, l_curv, l_len) ;
    /* indices[0] = startvno and indices[v_in_path-1] = endvno */
    for (i = 0 ; i < v_in_path ; i++)
    {
      vno = indices[i] ;
      if (vno == 0 || vno == Gdiag_no)
        DiagBreak() ;
      mris->vertices[vno].marked = 1 ;  /* mark it as part of path now */
      bnew->vertices[i] = vno ;
    }
    bnew->num_vertices = v_in_path ;

    sse = path_sse(bnew, mris, target_curv, l_curv, l_len) ;
    total_length = path_length(bnew, mris) ;
    if (niter++ > 20 && last_sse < sse)
      break ;
    path_set_marks(bnew, mris, 1) ;
    path_set_marks(bold, mris, 0) ;
    path_copy(bnew, bold) ;


    /*    if (sse < last_sse)*/
    path_copy(bnew, bold) ;

    printf("sse = %2.4f (%2.4f, %2.3f%%), %d vertices in path\n",
           sse, last_sse, 100*(last_sse-sse)/last_sse,bnew->num_vertices) ;
    if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
      redraw() ;
  }
  while (!FZERO(sse-last_sse)) ;
  free(indices) ;

}
static double
path_length(PATH_PATH *b, MRI_SURFACE *mris)
{
  double total_length, length ;
  int    i, vno ;
  VERTEX *v, *v_last = NULL ;

  /* compute length of current line */
  for (total_length = 0.0, i = 0 ; i < b->num_vertices ; i++)
  {
    vno = b->vertices[i] ;
    v = &mris->vertices[vno] ;
    if (i > 0)
    {
      length = sqrt(SQR(v->x-v_last->x) +
                    SQR(v->y-v_last->y) + SQR(v->z-v_last->z)) ;
      total_length += length ;
    }
    v_last = v ;
  }
  return(total_length) ;
}

static double
path_sse(PATH_PATH *b, MRI_SURFACE *mris, double target_curv,
         double l_curv, double l_len)
{
  double sse, curv_error, length ;
  int    i, vno ;
  VERTEX *v, *v_last = NULL ;

  /* compute length of current line */
  for (sse = 0.0, i = 0 ; i < b->num_vertices ; i++)
  {
    vno = b->vertices[i] ;
    v = &mris->vertices[vno] ;
    curv_error = (target_curv - v->curv) ;
    sse += l_curv * (curv_error*curv_error) ;
    if (i > 0)
    {
      length = sqrt(SQR(v->x-v_last->x) +
                    SQR(v->y-v_last->y) + SQR(v->z-v_last->z)) ;
      sse += l_len * length ;
    }
    v_last = v ;
  }
  return(sse) ;
}

static int
find_best_path(MRI_SURFACE *mris, int start_vno, int end_vno,
               int *path_indices,
               double target_curv, double l_curv, double l_len)
{
  int    *nbrs, *new_nbrs, num_new_nbrs, num_nbrs, n, done = 0, i, vno,
      min_n ;
  double sse, min_d, d ;
  VERTEX *v, *vn ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
    mris->vertices[vno].d = 0.0 ;

  nbrs = (int *)calloc(mris->nvertices, sizeof(int)) ;
  new_nbrs = (int *)calloc(mris->nvertices, sizeof(int)) ;
  v = &mris->vertices[end_vno] ;
  done = FALSE ;
  for (num_nbrs = n = 0 ; n < v->vnum ; n++)
  {
    vno = v->v[n] ;
    vn = &mris->vertices[vno] ;
    if (vn->marked && vno != start_vno)
      continue ;
    vn->d =
      l_len * sqrt(SQR(v->x-vn->x) + SQR(v->y-vn->y) + SQR(v->z-vn->z)) +
      l_curv * ((vn->curv-target_curv)*(vn->curv-target_curv));
    if (vno == start_vno)
      done = TRUE ;
    nbrs[num_nbrs++] = vno ;
  }

  do  /* keep expanding front until start_vno is found */
  {
    for (num_new_nbrs = i = 0 ; i < num_nbrs ; i++)
    {
      v = &mris->vertices[nbrs[i]] ;
      for (n = 0 ; n < v->vnum ; n++)
      {
        vno = v->v[n] ;
        vn = &mris->vertices[vno] ;
        if (vn->marked && (vno != start_vno))  /* already visited */
          continue ;
        d = v->d +
            l_len * sqrt(SQR(v->x-vn->x) +
                         SQR(v->y-vn->y) +
                         SQR(v->z-vn->z)) +
            l_curv * ((vn->curv-target_curv)*(vn->curv-target_curv));
        if (!FZERO(vn->d))  /* already added */
        {
          if (d < vn->d)
            vn->d = d ;
          continue ;    /* don't add it to nbr list twice */
        }
        if (vno == start_vno)
        {
          done = TRUE ;
          break ;
        }
        vn->d = d ;
        new_nbrs[num_new_nbrs++] = vno ;
      }
    }

    /* copy new nbr list into old one and do it again */
    if (done == FALSE)
      for (num_nbrs = 0 ; num_nbrs < num_new_nbrs ; num_nbrs++)
        nbrs[num_nbrs] = new_nbrs[num_nbrs] ;
  }
  while (done == FALSE) ;

  /* now search from start_vno to find shortest total path */
  done = FALSE ;
  sse = 0.0 ;
  path_indices[0] = vno = start_vno ;
  num_nbrs = 1 ;
  v = &mris->vertices[vno] ;
  do
  {
    min_n = -1 ;
    min_d = 10000000 ;
    for (n = 0 ; n < v->vnum ; n++)
    {
      vno = v->v[n] ;
      vn = &mris->vertices[vno] ;
      if (vno == end_vno)
        /* special case - only examine end vertex with d==0 */
      {
        done = TRUE ;
        min_n = n ;
        min_d = vn->d ;
        break ;
      }
      if (FZERO(vn->d))   /* wasn't visited */
        continue ;
      if (vn->d < min_d || min_n < 0)
      {
        min_d = vn->d ;
        min_n = n ;
      }
    }
    vno = v->v[min_n] ;
    path_indices[num_nbrs++] = vno ;
    v = &mris->vertices[vno] ;
    sse += v->d ;
  }
  while (done == FALSE) ;
  return(num_nbrs) ;
}

void
plot_curv(int closedcurveflag)
{
  float x0, y0, z0, dx, dy, dz, dist;
  VERTEX *v;
  FILE *fp;
  int vno;
  int* path;
  int path_length;


#define LINE_FNAME "surfer_curv.dat"

  if (closedcurveflag)
    close_marked_vertices();

  /* Find the path between the marked verts. */
  path = (int*) calloc (mris->nvertices, sizeof(int));

  find_path (marked, nmarked, "plotting curv", mris->nvertices,
             path, &path_length);

  /* Make sure we got a path. */
  if (path_length < 2)
  {
    printf ("surfer: needs at least 2 marked vertices\n");
    free (path);
    return;
  }

  /* Open the dest file. */
  fprintf(stderr, "generating data file %s with entries:\n", LINE_FNAME) ;
  fprintf(stderr, "vno x y z distance curv val val2 stat amp "
          "deg normalized_deg radians\n") ;
  if (nmarked<2)
  {
    printf("surfer: needs at least 2 marked vertices\n");
    PR return;
  }
  fp = fopen(LINE_FNAME, "w") ;

  /* Get the first location */
  x0 = mris->vertices[path[0]].x ;
  y0 = mris->vertices[path[0]].y ;
  z0 = mris->vertices[path[0]].z ;

  /* For each vert in the path, write some info about it. */
  for (vno = 0; vno < path_length; vno++)
  {
    v = &mris->vertices[path[vno]];

    /* Get the info from the first point. */
    dx = v->x - x0;
    dy = v->y - y0;
    dz = v->z - z0;
    dist = sqrt(dx*dx + dy*dy + dz*dz);

    /* Draw the vertices in blue. */
    v->marked = 2;

    /* Print info about this vertex to the file. */
    fprintf(fp, "%d %f %f %f %f %f %f %f %f %f %f %f %f\n",
            path[vno], v->x, v->y, v->z,
            dist, v->curv, v->val, v->val2, v->stat,
            hypot(v->val,v->val2),
            (float)(atan2(v->val2,v->val)*180/M_PI),
            (float)(atan2(v->val2,v->val)/(2*M_PI)),
            (float)(atan2(v->val2,v->val))) ;

  }

  /* Close data file. */
  fclose (fp);

  /* Clear the marked verts. */
  clear_vertex_marks();

  /* redraw screen with blue path. */
  redraw ();

  free (path);
}

void
clear_vertex_marks(void)
{
  int i;

  for (i=0;i<nmarked;i++)
    mris->vertices[marked[i]].marked = FALSE;
  nmarked = 0;
}
void
clear_all_vertex_marks(void)
{
  int i;

  for (i=0;i<mris->nvertices;i++)
    mris->vertices[i].marked = FALSE;
  nmarked = 0;
}

/* begin rkt */
void
find_closest_marked_vertex (int screen_x, int screen_y,
                            int* closest_index, int* closest_vno )
{
  int vno;
  float d;
  VERTEX* v;
  float ras[3];
  int test_marked;
  int test_vno;
  VERTEX* test_v;
  float dx, dy, dz;
  float dmin;
  int closest_marked;

  find_vertex_at_screen_point (screen_x, screen_y, &vno, &d);
  if (vno < 0)
  {
    if (closest_index) *closest_index = -1;
    if (closest_vno) *closest_vno = -1;
    return;
  }

  v = &(mris->vertices[vno]);
  ras[0] = v->x;
  ras[1] = v->y;
  ras[2] = v->z;

  dmin = 1000;
  closest_marked = -1;
  for (test_marked = 0; test_marked < nmarked; test_marked++)
  {
    test_vno = marked[test_marked];
    test_v = &(mris->vertices[test_vno]);
    dx = test_v->x - ras[0];
    dy = test_v->y - ras[1];
    dz = test_v->z - ras[2];
    d = sqrt (dx*dx + dy*dy + dz*dz);
    if (d < dmin)
    {
      dmin = d;
      closest_marked = test_marked;
    }
  }

  if (closest_index) *closest_index = closest_marked;
  if (closest_vno) *closest_vno = marked[closest_marked];
}
/* end rkt */

void
mark_translated_vertex(int vindex, int onoroff, char *surf_fname)
{
  MRI_SURFACE *mris2 ;
  float       x0, y0, z0, dist, min_dist, x, y, z ;
  int         vno, min_vno=0;
  char        surf_name[STRLEN] ;
  VERTEX      *v ;

  mris2 = MRISread(surf_fname) ;
  if (!mris2)
    ErrorExit(ERROR_NOFILE, "%s: could not load translation surface %s",
              Progname, surf_fname) ;
  FileNameOnly(surf_fname, surf_name) ;
  read_canon_vertex_coordinates(surf_name)  ;
  x0 = mris2->vertices[vindex].x ;
  y0 = mris2->vertices[vindex].y ;
  z0 = mris2->vertices[vindex].z ;
  min_dist = 1e10 ;
  min_vno = -1 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    x = v->cx ;
    y = v->cy ;
    z = v->cz ;
    dist = sqrt(SQR(x-x0) + SQR(y-y0) + SQR(z-z0)) ;
    if (dist < min_dist)
    {
      min_dist = dist ;
      min_vno = vno ;
    }
  }
  printf("translating to vertex %d (min dist = %2.1f)\n", min_vno, min_dist) ;
  mark_vertex(min_vno, onoroff) ;
}

void
mark_vertex(int vindex, int onoroff)
{
  int i,j;
  VERTEX *v ;

  mris->vertices[vindex].marked = onoroff;

  /* Find the index in marked[] of vindex. */
  for (i=0; i < nmarked && marked[i] != vindex; i++);

  v = &mris->vertices[vindex] ;
  if ((onoroff==FALSE)&&(i<nmarked))
  {
    v->marked = 0 ;
    nmarked--;
    for (j=i;j<nmarked;j++) marked[j] = marked[j+1];
    printf("surfer: vertex %d unmarked\n",vindex);
  }
  if ((onoroff==TRUE)&&(i==nmarked))
  {
    if (nmarked==mris->nvertices-1)
      printf("surfer: too many marked vertices\n");
    else
    {
      marked[nmarked] = vindex;
      printf("surfer: vertex %d marked (curv=%f, stat=%f)\n",
             vindex,mris->vertices[vindex].curv,
             mris->vertices[vindex].stat);
      printf("x = (%2.1f, %2.1f, %2.1f), n = (%2.1f, %2.1f, %2.1f).\n",
             v->x, v->y, v->z, v->nx, v->ny, v->nz) ;
      nmarked++;
    }
  }
  else if (onoroff > 1)
  {
    v = &mris->vertices[vindex] ;
    v->marked = onoroff ;
  }
}

void
mark_annotation(int vno_annot)
{
  int    vno, annotation ;
  VERTEX *v ;

  if (vno_annot < 0)
  {
    fprintf(stderr, "no vertex currently selected...\n") ;
    return ;
  }

  v = &mris->vertices[vno_annot] ;
  annotation = v->annotation ;
  fprintf(stderr, "marking annotation %d...\n", annotation) ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (v->annotation == annotation)
      v->marked = 1  ;
  }
}
void
mark_faces(int vno)
{
  int    n ;
  VERTEX *v;

  if (vno < 0 || vno >= mris->nvertices)
    return ;
  v = &mris->vertices[vno] ;
  for (n = 0 ; n < v->num ; n++)
    mark_face(v->f[n]) ;
}
void
mark_face(int fno)
{
  int    i ;
  FACE   *f;
  VERTEX *v;

  /* offscreen render opengl bug: RGBcolor(white) => surrounding faces color */

  if ((fno >= mris->nfaces)||(fno<0))
  {
    if (fno != -1)
      printf ("surfer: ### face index %d out of bounds\n",fno);
    return;
  }

#if 0
  if (onoroff==FALSE)
    set_color(0.0,0.0,GREEN_RED_CURV);
  else
#endif
    // (onoroff==TRUE)
  {
    if (blackcursorflag)  RGBcolor(0,0,0);
    else                  RGBcolor(0,255,255);
  }

  f = &mris->faces[fno] ;
  FaceNormCacheEntry const * fNorm = getFaceNorm(mris, fno);
  RGBcolor(0,255,255);
  bgnpolygon();
  for (i = 0 ; i< VERTICES_PER_FACE ; i++)
  {
    v = &mris->vertices[f->v[i]];
    load_brain_coords(v->nx,v->ny,v->nz,v1);
    n3f(v1);
    load_brain_coords(v->x+mup*fNorm->nx,v->y+mup*fNorm->ny,v->z+mup*fNorm->nz,v1);
    v3f(v1);
  }
#if 0
  linewidth(CURSOR_LINE_PIX_WIDTH);
  load_brain_coords(fNorm->nx,fNorm->ny,fNorm->nz,v1);
  n3f(v1);
#endif
  endpolygon();

#if 0
  if (bigcursorflag)
  {
    for (k=0;k<vselect->num;k++)
    {
      bgnquadrangle();
      f = &mris->faces[vselect->f[k]];
      for (n=0;n<VERTICES_PER_FACE;n++)
      {
        v = &mris->vertices[f->v[n]];
        load_brain_coords(v->nx,v->ny,v->nz,v1);
        n3f(v1);
        load_brain_coords(v->x+cup*v->nx,
                          v->y+cup*v->ny,
                          v->z+cup*v->nz,v1);
        v3f(v1);
      }
      endpolygon();
    }
  }
#endif
}
void
prompt_for_parameters(void)
{
  float  f1, f2, f3, f4, f5, f6, f7 ;
  printf("surfer: enter wt,wa,ws,wn,wc,wsh,wbn ");
  printf("surfer: (%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f)\n",
         wt,wa,ws,wn,wc,wsh,wbn);
  scanf("%f %f %f %f %f %f %f",&f1,&f2,&f3,&f4,&f5,&f6,&f7);
  wt = f1 ;
  wa = f2 ;
  ws = f3 ;
  wn = f4 ;
  wc = f5 ;
  wsh = f6 ;
  wbn = f7 ;
  PR
}

void
prompt_for_drawmask(void)
{
  printf("surfer: enter mask (curv,area,shear,shearvec,"
         "normvec,movevec) = %d: ",
         drawmask);
  scanf("%d",&drawmask);
  if ((drawmask&1)!=0)   surfcolor = CURVATURE_OR_SULCUS;
  if ((drawmask&2)!=0)   surfcolor = AREAL_DISTORTION;
  if ((drawmask&4)!=0)   surfcolor = SHEAR_DISTORTION;
  if ((drawmask&8)!=0)   shearvecflag = TRUE;
  if ((drawmask&16)!=0)  normvecflag = TRUE;
  if ((drawmask&32)!=0)  movevecflag = TRUE;
  PR
}

void
fill_triangle(float x0,float y0,float z0,float x1,float y1,float z1,
              float x2,float y2,float z2)
{
  int i,j,imnr;
  float d0,d1,d2,dmax,u,v;
  float px,py,pz,px0,py0,pz0,px1,py1,pz1;
  int numu,numv;

  d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
  d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
  d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
  dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
  numu = ceil(2*d0);
  numv = ceil(2*dmax);
  if (numu==0) numu=1;
  if (numv==0) numv=1;
  for (v=0;v<=numv;v++)
  {
    px0 = x0 + (x2-x0)*v/numv;
    py0 = y0 + (y2-y0)*v/numv;
    pz0 = z0 + (z2-z0)*v/numv;
    px1 = x1 + (x2-x1)*v/numv;
    py1 = y1 + (y2-y1)*v/numv;
    pz1 = z1 + (z2-z1)*v/numv;
    for (u=0;u<=numu;u++)
    {
      px = px0 + (px1-px0)*u/numu;
      py = py0 + (py1-py0)*u/numu;
      pz = pz0 + (pz1-pz0)*u/numu;
      if (curvim_allocated)
      {
        imnr = (int)((py-mris2->yctr)/fillscale+numimg/2.0+0.5);
        i = (int)(IMGSIZE/2.0-(pz-mris2->zctr)/fillscale+0.5);
        j = (int)(IMGSIZE/2.0-(px-mris2->xctr)/fillscale+0.5);
        curvim[imnr][i][j] = 255;
      }
      else
      {
        imnr = (int)((py-mris->yctr)/fillscale+numimg/2.0+0.5);
        i = (int)(IMGSIZE/2.0-(pz-mris->zctr)/fillscale+0.5);
        j = (int)(IMGSIZE/2.0-(px-mris->xctr)/fillscale+0.5);
        fill[imnr][i][j] = 255;
      }
    }
  }
}

void
fill_surface(void)
{
  int i,j,k;
  int totalfilled,newfilled;

  numimg = IMGSIZE;
  for (k=0;k<numimg;k++)
  {
    fill[k] = (unsigned char **)lcalloc(IMGSIZE,sizeof(char *));
    for (i=0;i<IMGSIZE;i++)
    {
      fill[k][i] = (unsigned char *)lcalloc(IMGSIZE,sizeof(char));
    }
  }

  for (k=0;k<mris->nfaces;k++)
  {
    fill_triangle(mris->vertices[mris->faces[k].v[0]].x,
                  mris->vertices[mris->faces[k].v[0]].y,
                  mris->vertices[mris->faces[k].v[0]].z,
                  mris->vertices[mris->faces[k].v[1]].x,
                  mris->vertices[mris->faces[k].v[1]].y,
                  mris->vertices[mris->faces[k].v[1]].z,
                  mris->vertices[mris->faces[k].v[2]].x,
                  mris->vertices[mris->faces[k].v[2]].y,
                  mris->vertices[mris->faces[k].v[2]].z);
#if VERTICES_PER_FACE == 4
    fill_triangle(mris->vertices[mris->faces[k].v[0]].x,
                  mris->vertices[mris->faces[k].v[0]].y,
                  mris->vertices[mris->faces[k].v[0]].z,
                  mris->vertices[mris->faces[k].v[2]].x,
                  mris->vertices[mris->faces[k].v[2]].y,
                  mris->vertices[mris->faces[k].v[2]].z,
                  mris->vertices[mris->faces[k].v[3]].x,
                  mris->vertices[mris->faces[k].v[3]].y,
                  mris->vertices[mris->faces[k].v[3]].z);
#endif
  }

  fill[1][1][1] = 64;
  totalfilled = newfilled = 1;
  while (newfilled>0)
  {
    newfilled = 0;
    for (k=1;k<numimg-1;k++)
      for (i=1;i<IMGSIZE-1;i++)
        for (j=1;j<IMGSIZE-1;j++)
          if (fill[k][i][j]==0)
            if (fill[k-1][i][j]==64||
                fill[k][i-1][j]==64||
                fill[k][i][j-1]==64)
            {
              fill[k][i][j] = 64;
              newfilled++;
            }
    for (k=numimg-2;k>=1;k--)
      for (i=IMGSIZE-2;i>=1;i--)
        for (j=IMGSIZE-2;j>=1;j--)
          if (fill[k][i][j]==0)
            if (fill[k+1][i][j]==64||
                fill[k][i+1][j]==64||
                fill[k][i][j+1]==64)
            {
              fill[k][i][j] = 64;
              newfilled++;
            }
    totalfilled += newfilled;
    printf("surfer: filled (fill): new=%d, total=%d of %d\n",
           newfilled,totalfilled,(numimg-2)*(IMGSIZE-2)*(IMGSIZE-2));
  }
  for (k=1;k<numimg-1;k++)
    for (i=1;i<IMGSIZE-1;i++)
      for (j=1;j<IMGSIZE-1;j++)
        fill[k][i][j] = (fill[k][i][j]==0)?128:0;
  PR
}

void
fill_second_surface(void)
{
  int i,j,k;
  int totalfilled,newfilled;

  if (!curvim_allocated)
  {
    printf("surfer: ### fill_second_surface failed: first make curv images\n");
    PR return;
  }

  for (k=0;k<mris2->nfaces;k++)
  {
    fill_triangle(mris2->vertices[mris2->faces[k].v[0]].x,
                  mris2->vertices[mris2->faces[k].v[0]].y,
                  mris2->vertices[mris2->faces[k].v[0]].z,
                  mris2->vertices[mris2->faces[k].v[1]].x,
                  mris2->vertices[mris2->faces[k].v[1]].y,
                  mris2->vertices[mris2->faces[k].v[1]].z,
                  mris2->vertices[mris2->faces[k].v[2]].x,
                  mris2->vertices[mris2->faces[k].v[2]].y,
                  mris2->vertices[mris2->faces[k].v[2]].z);
#if VERTICES_PER_FACE == 4
    fill_triangle(mris2->vertices[mris2->faces[k].v[0]].x,
                  mris2->vertices[mris2->faces[k].v[0]].y,
                  mris2->vertices[mris2->faces[k].v[0]].z,
                  mris2->vertices[mris2->faces[k].v[2]].x,
                  mris2->vertices[mris2->faces[k].v[2]].y,
                  mris2->vertices[mris2->faces[k].v[2]].z,
                  mris2->vertices[mris2->faces[k].v[3]].x,
                  mris2->vertices[mris2->faces[k].v[3]].y,
                  mris2->vertices[mris2->faces[k].v[3]].z);
#endif
  }

  ocurvim[1][1][1] = 64;   /* ocurvim for tmp: don't overwrite loaded curvim */
  totalfilled = newfilled = 1;
  while (newfilled>0)
  {
    newfilled = 0;
    for (k=1;k<numimg-1;k++)
      for (i=1;i<IMGSIZE-1;i++)
        for (j=1;j<IMGSIZE-1;j++)
          if (ocurvim[k][i][j]==0)
            if (ocurvim[k-1][i][j]==64||
                ocurvim[k][i-1][j]==64||
                ocurvim[k][i][j-1]==64)
            {
              ocurvim[k][i][j] = 64;
              newfilled++;
            }
    for (k=numimg-2;k>=1;k--)
      for (i=IMGSIZE-2;i>=1;i--)
        for (j=IMGSIZE-2;j>=1;j--)
          if (ocurvim[k][i][j]==0)
            if (ocurvim[k+1][i][j]==64||
                ocurvim[k][i+1][j]==64||
                ocurvim[k][i][j+1]==64)
            {
              ocurvim[k][i][j] = 64;
              newfilled++;
            }
    totalfilled += newfilled;
    printf("surfer: filled (curvim/second): new=%d, total=%d of %d\n",
           newfilled,totalfilled,(numimg-2)*(IMGSIZE-2)*(IMGSIZE-2));
  }
  for (k=1;k<numimg-1;k++)
    for (i=1;i<IMGSIZE-1;i++)
      for (j=1;j<IMGSIZE-1;j++)
        if (ocurvim[k][i][j]==0)
          curvimflags[k][i][j] |= CURVIM_FILLED;
  PR
}

void
resize_buffers(long frame_xdim, long frame_ydim)
{
  free(framebuff);
  free(framebuff2);
  free(framebuff3);
  free(binbuff);
  framebuff = (unsigned long *)calloc(frame_xdim*frame_ydim,sizeof(long));
  framebuff2 = (unsigned long *)calloc(frame_xdim*frame_ydim,sizeof(long));
  framebuff3 = (unsigned long *)calloc(frame_xdim*frame_ydim,sizeof(long));
  binbuff = (unsigned char *)calloc(3*frame_xdim*frame_ydim,sizeof(char));
}

void
open_window(char *name)
{
#ifdef OPENGL
#ifndef USE_XGLUT_WINDOW
  XSizeHints hin;
  int success;
#endif /* USE_XGLUT_WINDOW */

  if (openglwindowflag)
  {
    printf("surfer: ### GL window already open: can't open second\n");
    PR return;
  }
  xmin=ymin=zmin= -(xmax=ymax=zmax=fov);

#ifdef USE_XGLUT_WINDOW

  wndw_create (MOTIF_XFUDGE + 5, MOTIF_YFUDGE + 5, frame_xdim, frame_ydim);
  wndw_set_title (name);

#else /* USE_XGLUT_WINDOW */

  if (doublebufferflag)
  {
    tkoInitDisplayMode(TKO_DOUBLE | TKO_RGB | TKO_DEPTH);
    printf("surfer: double buffered window\n");
  }
  else
  {
    tkoInitDisplayMode(TKO_SINGLE | TKO_RGB | TKO_DEPTH);
    printf("surfer: single buffered window\n");
  }

  if (renderoffscreen)
  {
    renderoffscreen1 = renderoffscreen; /* val at winopen time safe from tcl */
    tkoInitPosition(MOTIF_XFUDGE+5,MOTIF_YFUDGE+5,frame_xdim,frame_ydim);
    tkoInitPixmap(frame_xdim,frame_ydim);
    printf("surfer: ### rendering to offscreen window ###\n");
    PR
  }
  else
  {
    if (!initpositiondoneflag)
      tkoInitPosition(MOTIF_XFUDGE+wx0,(1024-wy0-frame_ydim)+MOTIF_XFUDGE,
                      frame_xdim,frame_ydim);
    /* begin rkt */
    /* if we don't get a display, try again with the other kind of
    buffer. */
    printf("surfer: tkoInitWindow(%s)\n",name);
    success = tkoInitWindow(name);
    if (!success)
    {
      if (doublebufferflag)
      {
        tkoInitDisplayMode(TKO_SINGLE | TKO_RGB | TKO_DEPTH);
        printf("surfer: failed, trying single buffered window\n");
      }
      else
      {
        tkoInitDisplayMode(TKO_DOUBLE | TKO_RGB | TKO_DEPTH);
        printf("surfer: failed, trying double buffered window\n");
      }
      success = tkoInitWindow(name);
    }
    if (!success)
    {
      printf("surfer: failed, no suitable display found\n");
      exit(1);
    }
    hin.max_width = hin.max_height = 1024;             /* maxsize */
    hin.min_aspect.x = hin.max_aspect.x = frame_xdim;  /* keepaspect */
    hin.min_aspect.y = hin.max_aspect.y = frame_ydim;
    hin.flags = PMaxSize|PAspect;
    XSetWMNormalHints(xDisplay, w.wMain, &hin);
  }
#endif /* USE_XGLUT_WINDOW */

  glLoadIdentity();
  glMatrixMode(GL_PROJECTION);
  glOrtho(-fov*sf, fov*sf, -fov*sf, fov*sf, -10.0*fov*sf, 10.0*fov*sf);
  glClearColor(0.0, 0.0, 0.0, 1.0);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  do_lighting_model(-1.0,-1.0,-1.0,-1.0,-1.0);
  glDepthFunc(GL_LEQUAL);
  glEnable(GL_DEPTH_TEST);
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
  glTranslatef(mris->xctr,-mris->zctr,-mris->yctr);
  /* MRI coords (vs. viewer:translate_brain_ */


#else  /* use gl calls */
  if (openglwindowflag)
  {
    printf("surfer: ### gl window already open: can't open second\n");
    PR return;
  }
  if (renderoffscreen)
  {
    printf("surfer: ### render offscreen request ignored "
           "(OpenGL version only)\n");
    renderoffscreen1 = FALSE;
  }
  xmin=ymin=zmin= -(xmax=ymax=zmax=fov);

  prefposition(0,frame_xdim-1,0,frame_ydim-1);
  foreground();  /* tcl prompt */
  winopen(name);
  keepaspect(1,1);     /* resize_window assumes square */
  maxsize(1024,1024);
  winconstraints();    /* makes re-sizeable */
  drawmode(NORMALDRAW);
  RGBmode();
  if (doublebufferflag) doublebuffer();
  else                  singlebuffer();
  gconfig();
  zbuffer(TRUE);
  mmode(MVIEWING);
  ortho(-fov*sf, fov*sf, -fov*sf, fov*sf, -10.0*fov*sf, 10.0*fov*sf);
  do_lighting_model(-1.0,-1.0,-1.0,-1.0,-1.0);
  qdevice(ESCKEY);
  qdevice(REDRAW);
  qdevice(LEFTMOUSE);
  qdevice(RIGHTMOUSE);
  qdevice(MIDDLEMOUSE);
  qdevice(UPARROWKEY);
  qdevice(DOWNARROWKEY);
  qdevice(LEFTARROWKEY);
  qdevice(RIGHTARROWKEY);
  qdevice(DELKEY);
  qdevice(PAGEDOWNKEY);
  qdevice(PAD2);
  qdevice(PAD4);
  qdevice(PAD6);
  qdevice(PAD8);
  qdevice(LEFTSHIFTKEY);
  qdevice(RIGHTSHIFTKEY);
  qdevice(LEFTCTRLKEY);
  qdevice(RIGHTCTRLKEY);
  qdevice(LEFTALTKEY);
  qdevice(RIGHTALTKEY);
  qdevice(KEYBD);
  translate(mris->xctr,-mris->zctr,-mris->yctr);
#endif

  openglwindowflag = TRUE;
}

void
swap_buffers(void)  /* so tcl sees OpenGL define */
{
  swapbuffers();
}

void
to_single_buffer(void)
{
  if (!doublebufferflag) return;
#ifdef OPENGL
  /* TODO: to toggle single/double in X: init second, glXChooseVisual */
  printf("surfer: ### toggle single/double not implemented in OPENGL\n");
  PR printf("surfer: ### tcl: set doublebufferflag FALSE before "
            "open_window\n");
  PR printf("surfer:     csh: setenv doublebufferflag 0; "
            "restart tksurfer\n");
  PR
#else
  if (!openglwindowflag)
  {
    printf("surfer: ### gl window not open: use open_window\n");
    PR return;
  }
  doublebufferflag = FALSE;
  singlebuffer();
  gconfig();
#endif
}

void
to_double_buffer(void)
{
  if (doublebufferflag) return;
#ifdef OPENGL
  /* TODO: to toggle single/double in X: init second, glXChooseVisual */
  printf("surfer: ### toggle single/double not implemented in OPENGL\n");
  PR printf("surfer: ### tcl: set doublebufferflag TRUE before "
            "open_window\n");
  PR printf("surfer:     csh: setenv doublebufferflag 1; "
            "restart tksurfer\n");
  PR
#else
  if (!openglwindowflag)
  {
    printf("surfer: ### gl window not open: use open_window\n");
    PR
    return;
  }
  doublebufferflag = TRUE;
  doublebuffer();
  gconfig();
#endif
}

void
do_lighting_model(float lite0,float lite1,float lite2,float lite3,
                  float newoffset)
{
#ifdef OPENGL
  /* init static arrays for non-default vars */
  static GLfloat mat0_ambient[] =
    {
      0.0, 0.0, 0.0, 1.0
    };
  static GLfloat mat0_diffuse[] =
    {
      OFFSET, OFFSET, OFFSET, 1.0
    };
  static GLfloat mat0_emission[] =
    {
      0.0, 0.0, 0.0, 1.0
    };
  static GLfloat light0_diffuse[] =
    {
      LIGHT0_BR, LIGHT0_BR, LIGHT0_BR, 1.0
    };
  static GLfloat light1_diffuse[] =
    {
      LIGHT1_BR, LIGHT1_BR, LIGHT1_BR, 1.0
    };
  static GLfloat light2_diffuse[] =
    {
      LIGHT2_BR, LIGHT2_BR, LIGHT2_BR, 1.0
    };
  static GLfloat light3_diffuse[] =
    {
      LIGHT3_BR, LIGHT3_BR, LIGHT3_BR, 1.0
    };
  static GLfloat light0_position[] =
    {
      0.0, 0.0, 1.0, 0.0
    };
  /* w=0:directional*/
  static GLfloat light1_position[] =
    {
      0.0, 0.0,-1.0, 0.0
    };
  static GLfloat light2_position[] =
    {
      0.6, 0.6, 1.6, 0.0
    }
    ;/* 0.6->1.6 */
  static GLfloat light3_position[] =
    {
      -1.0, 0.0, 0.0, 0.0
    };
  static GLfloat lmodel_ambient[] =
    {
      0.0, 0.0, 0.0, 0.0
    };
  /* cancel 0.2 amb */

  /* neg lite => no change (func args override interface light?_br update) */
  if (lite0 < 0.0)     lite0 = light0_diffuse[0];
  if (lite1 < 0.0)     lite1 = light1_diffuse[0];
  if (lite2 < 0.0)     lite2 = light2_diffuse[0];
  if (lite3 < 0.0)     lite3 = light3_diffuse[0];
  if (newoffset < 0.0) newoffset = offset;
  else                 offset = newoffset;
  light0_br = lite0;   /* update global double copies for tcl */
  light1_br = lite1;
  light2_br = lite2;
  light3_br = lite3;

  /* material: change DIFFUSE,EMISSION (purpler: EMISSIONg=0.05*newoffset) */
  /* NOTE: on top of (!) direct (set_color,etc) effects of offset on color */
  mat0_diffuse[0] = mat0_diffuse[1] = mat0_diffuse[2] = newoffset;
  mat0_emission[0] = mat0_emission[1] = mat0_emission[2] = 0.1*newoffset;

  /* lights: change DIFFUSE */
  light0_diffuse[0] = light0_diffuse[1] = light0_diffuse[2] = lite0;
  light1_diffuse[0] = light1_diffuse[1] = light1_diffuse[2] = lite1;
  light2_diffuse[0] = light2_diffuse[1] = light2_diffuse[2] = lite2;
  light3_diffuse[0] = light3_diffuse[1] = light3_diffuse[2] = lite3;

  glMatrixMode(GL_MODELVIEW);
  pushmatrix();
  glLoadIdentity();

  glMaterialfv(GL_FRONT, GL_AMBIENT, mat0_ambient);
  glMaterialfv(GL_FRONT, GL_DIFFUSE, mat0_diffuse);
  glMaterialfv(GL_FRONT, GL_EMISSION, mat0_emission);

  glLightfv(GL_LIGHT0, GL_DIFFUSE, light0_diffuse);
  glLightfv(GL_LIGHT1, GL_DIFFUSE, light1_diffuse);
  glLightfv(GL_LIGHT2, GL_DIFFUSE, light2_diffuse);
  glLightfv(GL_LIGHT3, GL_DIFFUSE, light3_diffuse);

  /* might need to move */
  glLightfv(GL_LIGHT0, GL_POSITION, light0_position);
  glLightfv(GL_LIGHT1, GL_POSITION, light1_position);
  glLightfv(GL_LIGHT2, GL_POSITION, light2_position);
  glLightfv(GL_LIGHT3, GL_POSITION, light3_position);

  glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
  glEnable(GL_NORMALIZE);
  glEnable(GL_COLOR_MATERIAL); /* def: follow AMBIENT_AND_DIFFUSE */
  glEnable(GL_LIGHTING);
  if (light0_br > 0.0) glEnable(GL_LIGHT0);
  else glDisable(GL_LIGHT0);
  if (light1_br > 0.0) glEnable(GL_LIGHT1);
  else glDisable(GL_LIGHT1);
  if (light2_br > 0.0) glEnable(GL_LIGHT2);
  else glDisable(GL_LIGHT2);
  if (light3_br > 0.0) glEnable(GL_LIGHT3);
  else glDisable(GL_LIGHT3);
  popmatrix();

#else  /* use gl calls */
  /* init float arrays */
  static float material0[] =
    {
      ALPHA,1.0,
      AMBIENT,0.0,0.0,0.0,
      DIFFUSE,OFFSET,OFFSET,OFFSET,
      EMISSION,0.0,0.0,0.0,
      SHININESS,0.0,
      LMNULL
    };
  static float light0[] =
    {
      AMBIENT,0.0,0.0,0.0,
      LCOLOR,LIGHT0_BR,LIGHT0_BR,LIGHT0_BR,
      POSITION,0.0,0.0,1.0,0.0,
      LMNULL
    };
  static float light1[] =
    {
      AMBIENT,0.0,0.0,0.0,
      LCOLOR,LIGHT1_BR,LIGHT1_BR,LIGHT1_BR,
      POSITION,0.0,0.0,-1.0,0.0,
      LMNULL
    };
  static float light2[] =
    {
      AMBIENT,0.0,0.0,0.0,
      LCOLOR,LIGHT2_BR,LIGHT2_BR,LIGHT2_BR,
      POSITION,0.6,0.6,0.6,0.0,
      LMNULL
    };
  static float light3[] =
    {
      AMBIENT,0.0,0.0,0.0,
      LCOLOR,LIGHT3_BR,LIGHT3_BR,LIGHT3_BR,
      POSITION,-1.0,0.0,0.0,0.0,
      LMNULL
    };

  if (lite0 < 0.0)     lite0 = light0[5];
  if (lite1 < 0.0)     lite1 = light1[5];
  if (lite2 < 0.0)     lite2 = light2[5];
  if (lite3 < 0.0)     lite3 = light3[5];
  if (newoffset < 0.0) newoffset = offset;
  else                 offset = newoffset;

  /* material: change DIFFUSE,EMISSION (purpler:material0[12]=0.05*newoffset)*/
  material0[7] = material0[8] = material0[9] = newoffset;
  material0[11] = material0[12] = material0[13] = 0.1*newoffset;
  /* material0[12] = 0.05*newoffset; *//* purpler */

  /* lights: change LCOLOR */
  light0[5] = light0[6] = light0[7] = lite0;
  light1[5] = light1[6] = light1[7] = lite1;
  light2[5] = light2[6] = light2[7] = lite2;
  light3[5] = light3[6] = light3[7] = lite3;

  lmdef(DEFMATERIAL,1,17,material0);  /* was: 100,17 */
  lmdef(DEFLIGHT,1,14,light0);
  lmdef(DEFLIGHT,2,14,light1);
  lmdef(DEFLIGHT,3,14,light2);
  lmdef(DEFLIGHT,4,14,light3);
  lmdef(DEFLMODEL,1,0,NULL);  /* default? */

  pushmatrix();
  loadmatrix(idmat);
  lmbind(LIGHT0,1);
  lmbind(LIGHT1,2);
  lmbind(LIGHT2,3);
  lmbind(LIGHT3,4);
  lmbind(LMODEL,1);
  lmbind(MATERIAL,1);  /* was: 100 */
  popmatrix();

  lmcolor(LMC_DIFFUSE);
#endif

  if (Gdiag & DIAG_SHOW && DIAG_VERBOSE_ON)
  {
    printf("surfer: offset=%3.2f\n",newoffset);
    printf("surfer: light0=%3.2f, light1=%3.2f, light2=%3.2f, "
           "light3=%3.2f\n",
           lite0,lite1,lite2,lite3);
  }
}

void
invert_surface(void)
{
  int     vno ;
  VERTEX  *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->nx *= -1.0 ;
    v->ny *= -1.0 ;
    v->nz *= -1.0 ;
  }
}

void
make_filenames(char *lsubjectsdir,char *lsrname,char *lpname,char *lstem,
               char *lext)
{
  /* malloc for tcl */
  sphere_reg_contra = (char *)malloc(NAME_LENGTH*sizeof(char));
  sphere_reg = (char *)malloc(NAME_LENGTH*sizeof(char));
  subjectsdir = (char *)malloc(NAME_LENGTH*sizeof(char));
  srname = (char *)malloc(NAME_LENGTH*sizeof(char));
  pname = (char *)malloc(NAME_LENGTH*sizeof(char));
  stem = (char *)malloc(NAME_LENGTH*sizeof(char));
  ext = (char *)malloc(NAME_LENGTH*sizeof(char));
  fpref = (char *)malloc(NAME_LENGTH*sizeof(char));
  ifname = (char *)malloc(NAME_LENGTH*sizeof(char));
  if2name = (char *)malloc(NAME_LENGTH*sizeof(char));
  ofname = (char *)malloc(NAME_LENGTH*sizeof(char));
  cfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  cf2name = (char *)malloc(NAME_LENGTH*sizeof(char));
  sfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  dfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  kfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  mfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  vfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  fsfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  fmfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  afname = (char *)malloc(NAME_LENGTH*sizeof(char));
  pfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  tfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  gfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  sgfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  agfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  fifname = (char *)malloc(NAME_LENGTH*sizeof(char));
  cif2name = (char *)malloc(NAME_LENGTH*sizeof(char));
  orfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  paint_fname = (char *)malloc(NAME_LENGTH*sizeof(char));
  elfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  lfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  vrfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  xffname = (char *)malloc(NAME_LENGTH*sizeof(char));
  /* following not set below */
  nfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  nfname[0]=0;
  rfname = (char *)malloc(NAME_LENGTH*sizeof(char));
  rfname[0]=0;
  pname2 = (char *)malloc(NAME_LENGTH*sizeof(char));
  pname2[0]=0;
  stem2 = (char *)malloc(NAME_LENGTH*sizeof(char));
  stem2[0]=0;
  ext2 = (char *)malloc(NAME_LENGTH*sizeof(char));
  ext2[0]=0;
  tf2name = (char *)malloc(NAME_LENGTH*sizeof(char));
  tf2name[0]=0;

  /* make default names */
  strcpy(sphere_reg_contra, sphere_reg_contra_suffix) ;
  strcpy(sphere_reg, sphere_reg_suffix) ;
  strcpy(subjectsdir,lsubjectsdir);
  strcpy(srname,lsrname);
  strcpy(pname,lpname);
  strcpy(stem,lstem);
  strcpy(ext,lext);
  sprintf(fpref,"%s/%s/%s/%s",subjectsdir,pname,BRAINSURF_DIR,stem);
  sprintf(ifname,"%s.%s",fpref,ext);
  sprintf(if2name,"%s.%s",fpref,ext);
  sprintf(ofname,"%s.surftmp",fpref);  /* [return to] default */
  sprintf(cfname,"%s.curv",fpref);
  sprintf(cf2name,"%s.curv",fpref);
  sprintf(sfname,"%s/%s/%s-%d.dec",srname,BEM_DIR,stem,(int)dip_spacing);
  sprintf(dfname,"%s/%s/%s.dip",srname,BEM_DIR,stem);
  sprintf(kfname,"%s.sulc",fpref);
  sprintf(mfname,"%s/%s/%s/COR-",subjectsdir,pname,DEVOLT1_DIR);
  sprintf(vfname,"%s.val",fpref);
  sprintf(vfname,"%s/3/%s.val",srname,stem);
  sprintf(fsfname,"%s/%s/%s.fs",srname,FS_DIR,stem);
  sprintf(fmfname,"%s/%s/%s.fm",srname,FS_DIR,stem);
  sprintf(afname,"%s.area",fpref);
  sprintf(pfname,"%s.patch",fpref);
  sprintf(tfname,"%s/%s/%s",subjectsdir,pname,TMP_DIR);
  sprintf(gfname,"%s/%s/%s",subjectsdir,pname,IMAGE_DIR);
  sprintf(sgfname,"%s/%s",srname,IMAGE_DIR);
  sprintf(agfname,"%s/%s/%s",srname,IMAGE_DIR,"tksurfer.rgb");
  sprintf(vrfname,"%s/%s/%s/%s.%s.wrl",
          subjectsdir,pname,BRAINSURF_DIR,stem,ext);
  sprintf(xffname,"%s/%s/%s/%s",
          subjectsdir,pname,TRANSFORM_DIR,TALAIRACH_FNAME);
  /* ~/morph/curv.rh.1000a2adj: curv (or fill), need hemi */
  sprintf(fifname,"%s/%s/%s.%s.%s/COR-",
          subjectsdir,pname,FILLDIR_STEM,stem,ext);
  sprintf(cif2name,"%s/%s/%s.%s.%s/COR-",
          subjectsdir,pname,CURVDIR_STEM,stem,ext);

  if (getenv("USE_WHITE") == NULL)
    sprintf(orfname,"%s.%s",fpref, orig_suffix);
  else
    sprintf(orfname,"%s.white",fpref);
  sprintf(paint_fname,"%s.smoothwm",fpref);
  sprintf(elfname,"%s.1000a2ell",fpref);
#if 1
  sprintf(lfname,"%s/%s/%s/%s-area.label",subjectsdir,pname,LABEL_DIR,stem);
#else
  sprintf(lfname,"%s-area.label",stem);
#endif
}

void
print_help_tksurfer(void)
{
  printf("usage: tksurfer subject hemisphere surface [options]\n");
  printf("\n");
  printf("subject    : a subject directory in the SUBJECTS_DIR path\n");
  printf("hemipshere : rh or lh\n");
  printf("surface    : a surface file name relative to the "
         "subject's surf directory\n");
  printf("\n");

  printf("Options\n");
  printf("\n");
  printf("-title title : set window title to title. Def is subject name.\n");

  printf("-reassign      : resample labels onto surface (set vnos=-1)\n");
  printf("-sdir <path>   : sets the subjects directory path\n");
  printf("-orig <suffix> : sets the orig suffix string\n");
  printf("-sphere <suffix>:sets the sphere.reg suffix string\n");
  printf("\n");

  printf("-patch <filename> : load a patch\n");
  printf("\n");

  printf("-tcl <filename> : run a script\n");
  printf("\n");

  printf("-annotation <filename> : load an annotation\n");
  printf("-aparc : set annotation to aparc.annot and use outline mode\n");
  printf("-a2009s : set annotation to aparc.a2009s.annot and use outline mode\n");
  printf("-snap <filename> : save initial view to <filename> and exit\n");
  printf("-colortable <filename> : load a color table file\n");
  printf("-labels-under : display labels under any overlay\n");
  printf("-label-outline : draw labels as outlines\n");
  printf("\n");
  printf("-curv : automatically load ?h.curv\n");
  printf("-gray : automatically load ?h.curv and make it gray\n");
  printf("\n");

  printf("-overlay          <filename> : load an overlay volume\n");
  printf("-overlay-reg      <filename> : use a file for the overlay registration\n");
  printf("-overlay-reg-find            : look in the data directory for a register.dat\n");
  printf("-fsgd <gd fname> <overlay fname> : load an FSGD file and an overlay volume\n");
  printf("                             : file\n");
  printf("-overlay-reg-identity        : calculate an identity transform for registration\n");
  printf("-mni152reg : for use with average subject when overlay is mni152\n");
  printf("-zm                          : remove mean from overlays\n") ;
  printf("-fminmax min max             : set the overlay threshold min (fthresh) and max\n");
  printf("-fslope <value>              : set the overlay threshold slope value\n");
  printf("-fmid <value>                : set the overlay threshold midpoint value\n");
  printf("-fthresh <value>             : set the overlay threshold minimum value\n");
  printf("-foffset <value>             : set the overlay threshold offset value\n");
  printf("-wide-config-overlay         : enable the wide mode of Configure Overlay Display dialog\n");
  printf("             can also 'setenv FS_TKS_OV_WIDE 1'\n");
  printf("-long-config-overlay         : enable the vertical mode of Configure Overlay Display dialog (default)\n");
  printf("-colscalebarflag <1|0>       : display color scale bar\n");
  printf("-colscaletext <1|0>          : display text in color scale bar\n");
  printf("-truncphaseflag <1|0>        : truncate the overlay display\n");
  printf("-revphaseflag <1|0>          : reverse the overlay display\n");
  printf("-invphaseflag <1|0>          : invert the overlay display\n");
  printf("\n");

  printf("-timecourse          <filename>        : load an timecourse volume\n");
  printf("-timecourse-reg      <filename>        : use a file for the timecourse\n");
  printf("                                       : registration\n");
  printf("-timecourse-reg-find                   : look in the data directory for a\n");
  printf("                                       : register.dat file\n");
  printf("-timecourse-reg-identity               : calculate an identity transform for\n");
  printf("                                       : regisrtation\n");
  printf("-timecourse-offset          <filename> : load an timecourse offset volume\n");
  printf("-timecourse-offset-reg-file <filename> : use a file for the timecourse offset\n");
  printf("                                       : registration\n");
  printf("-timecourse-offset-reg-find            : look in the data directory for a\n");
  printf("                                       : register.dat file\n");
  printf("-timecourse-offset-reg-identity        : calculate an identity transform for\n");
  printf("                                       : registration\n");
  printf("\n");
  printf("-delink : do not move tool window with image window\n");

  printf("-scalebarflag <1|0> : display the scale bar\n");
  printf("\n");
}

#endif /* defined(TCL) && defined(TKSURFER) */


/* end of surfer.c and start of tksurfer.c */



#if 1 /*defined(Linux) || defined(SunOS) || defined(sun)*/
#include <sys/time.h>
#include <sys/types.h>
#include <unistd.h>
#endif
#if 1 /* defined(Linux) || defined(SunOS) || defined(sun)*/
int tk_NumMainWindows = 1 ;
#endif

/* boilerplate wrap function defines for easier viewing */
#define WBEGIN (ClientData clientData,\
Tcl_Interp *interp,int argc,char *argv[]){
/* rkt begin - changed to use Tcl_SetResult and TCL_VOLATILE */
#define ERR(N,S)  if(argc!=N){Tcl_SetResult(interp,S,TCL_VOLATILE);\
return TCL_ERROR;}
// #define ERR(N,S)  if(argc!=N){interp->result=S;return TCL_ERROR;}
/* rkt end */
#define WEND   return TCL_OK;}
#define REND  (ClientData) NULL, (Tcl_CmdDeleteProc *) NULL
#define PARM (ClientData clientData,Tcl_Interp *interp,int argc,char *argv[])

/* Prototypes */
int W_swap_buffers PARM;
int W_to_single_buffer PARM;
int W_to_double_buffer PARM;
int W_open_window PARM;
int W_help  PARM;
int W_redraw  PARM;
int W_dump_vertex  PARM;
int W_show_flat_regions  PARM;
int W_transform_brain  PARM;
int W_resize_brain  PARM;
int W_set_area_thresh  PARM;
int W_val_to_mark  PARM;
int W_val_to_curv  PARM;
int W_val_to_stat  PARM;
int W_set_vals  PARM;
int W_stat_to_val  PARM;

int W_read_imag_vals  PARM;
int W_read_soltimecourse  PARM;
int W_sol_plot  PARM;

int W_remove_triangle_links  PARM;

int W_taubin_smooth  PARM;

int W_label_to_stat  PARM;
int W_label_from_stats  PARM;
int W_set_stats  PARM;
int W_f_to_t  PARM;
int W_t_to_p  PARM;
int W_f_to_p  PARM;

int W_scale_vals  PARM;
int W_curv_to_val  PARM;
int W_read_curv_to_val  PARM;
int W_read_parcellation  PARM;
int W_read_and_smooth_parcellation  PARM;
int W_read_disc  PARM;
int W_deconvolve_weights  PARM;
int W_read_stds  PARM;
int W_mark_label  PARM;
int W_orient_sphere  PARM;
int W_dump_faces  PARM;
int W_load_gcsa  PARM;
int W_mark_annotation  PARM;
int W_mark_face  PARM;
int W_invert_face  PARM;
int W_mark_faces  PARM;
int W_invert_vertex  PARM;
int W_redraw_second  PARM;
int W_shrink  PARM;
int W_area_shrink  PARM;
int W_sphere_shrink  PARM;
int W_ellipsoid_project  PARM;
int W_ellipsoid_morph  PARM;
int W_ellipsoid_shrink  PARM;
int W_ellipsoid_shrink_bug  PARM;
int W_curv_shrink_to_fill  PARM;
int W_smooth_curvim  PARM;
int W_smooth_curv  PARM;
int W_smooth_val  PARM;
int W_smooth_val_sparse  PARM;
int W_smooth_curvim_sparse  PARM;
int W_smooth_fs  PARM;
int W_add_subject_to_average_curvim  PARM;
int W_read_curv_images PARM;
int W_read_second_binary_surf  PARM;
int W_read_second_binary_curv  PARM;
int W_normalize_second_binary_curv  PARM;
int W_curv_to_curvim PARM;
int W_second_surface_curv_to_curvim PARM;
int W_curvim_to_second_surface PARM;
int W_swap_curv PARM;
int W_curvim_to_surface PARM;
int W_read_binary_curv  PARM;
int W_read_binary_sulc  PARM;
int W_read_binary_values  PARM;
int W_read_binary_values_frame  PARM;
int W_read_annotated_image  PARM;
int W_read_annotateds  PARM;
int W_read_binary_patch  PARM;
int W_read_fieldsign  PARM;
int W_read_fsmask PARM;
int W_write_binary_areas  PARM;
int W_write_binary_surface  PARM;
int W_write_binary_curv  PARM;
int W_write_binary_sulc  PARM;
int W_write_binary_values  PARM;
int W_write_binary_patch  PARM;
int W_write_labeled_vertices  PARM;
int W_read_labeled_vertices  PARM;
int W_read_and_color_labeled_vertices  PARM;
int W_write_fieldsign  PARM;
int W_write_fsmask PARM;
int W_write_vrml PARM;
int W_write_binary_dipoles PARM;
int W_write_binary_decimation PARM;
int W_write_dipoles PARM;
int W_write_decimation PARM;
int W_write_curv_images PARM;
int W_write_fill_images PARM;
int W_fill_second_surface  PARM;
int W_subsample_dist PARM;
int W_subsample_orient PARM;
int W_write_subsample PARM;
int W_write_binary_subsample PARM;
int W_compute_curvature  PARM;
int W_compute_CMF  PARM;
int W_compute_cortical_thickness PARM;
int W_clear_curvature  PARM;
int W_clear_vals  PARM;
int W_clear_ripflags  PARM;
int W_restore_ripflags  PARM;
int W_floodfill_marked_patch  PARM;
int W_rip_unmarked_vertices  PARM;
int W_dilate_ripped  PARM;
int W_twocond  PARM;
int W_cut_line  PARM;
int W_plot_curv  PARM;
int W_draw_fundus  PARM;
int W_plot_marked  PARM;
int W_put_retinotopy_stats_in_vals  PARM;
int W_draw_vector  PARM;
int W_cut_plane  PARM;
int W_flatten  PARM;
int W_normalize_binary_curv  PARM;
int W_normalize_curvature  PARM;
int W_normalize_area  PARM;
int W_shift_values  PARM;
int W_swap_values  PARM;
int W_swap_stat_val  PARM;
int W_swap_val_val2  PARM;
int W_compute_angles  PARM;
int W_compute_fieldsign  PARM;
int W_draw_radius PARM;
int W_draw_theta PARM;
int W_save_rgb  PARM;
int W_save_rgb_num  PARM;
int W_save_rgb_named_orig  PARM;
int W_save_rgb_cmp_frame  PARM;
int W_open_rgb_cmp_named  PARM;
int W_save_rgb_cmp_frame_named  PARM;
int W_close_rgb_cmp_named  PARM;
int W_rotate_brain_x  PARM;
int W_rotate_brain_y  PARM;
int W_rotate_brain_z  PARM;
int W_translate_brain_x  PARM;
int W_translate_brain_y  PARM;
int W_translate_brain_z  PARM;
int W_scale_brain  PARM;
int W_resize_window  PARM;
int W_setsize_window  PARM;
int W_do_lighting_model  PARM;
int W_restore_zero_position  PARM;
int W_restore_initial_position  PARM;
int W_make_lateral_view  PARM;
int W_make_lateral_view_second  PARM;
int W_read_view_matrix  PARM;
int W_write_view_matrix  PARM;
int W_read_really_matrix  PARM;
int W_write_really_matrix  PARM;
int W_really_translate_brain  PARM;
int W_really_scale_brain  PARM;
int W_align_sphere  PARM;
int W_really_rotate_brain_x  PARM;
int W_really_rotate_brain_y  PARM;
int W_really_rotate_brain_z  PARM;
int W_really_center_brain  PARM;
int W_really_center_second_brain  PARM;
int W_really_align_brain  PARM;
int W_read_binary_decimation  PARM;
int W_read_binary_dipoles  PARM;
int W_load_vals_from_sol  PARM;
int W_load_var_from_sol  PARM;
int W_read_iop  PARM;
int W_read_rec  PARM;
int W_read_ncov  PARM;
int W_filter_recs  PARM;
int W_normalize_time_courses  PARM;
int W_compute_timecourses  PARM;
int W_normalize_inverse  PARM;
int W_compute_pval_fwd  PARM;
int W_compute_select_fwd  PARM;
int W_compute_pval_inv  PARM;
int W_find_orig_vertex_coordinates  PARM;
int W_select_orig_vertex_coordinates  PARM;
int W_select_talairach_point  PARM;
int W_read_orig_vertex_coordinates  PARM;

int W_save_surf  PARM;
int W_restore_surf  PARM;
int W_read_pial_vertex_coordinates  PARM;
int W_read_white_vertex_coordinates  PARM;

int W_read_canon_vertex_coordinates  PARM;
int W_send_spherical_point  PARM;
int W_send_contralateral_point  PARM;
int W_send_to_subject  PARM;
int W_send_to_subject  PARM;
int W_resend_to_subject  PARM;
int W_drawcb  PARM;
int W_read_ellipsoid_vertex_coordinates  PARM;
int W_invert_surface  PARM;
int W_fix_nonzero_vals  PARM;
int W_draw_ellipsoid_latlong  PARM;
int W_left_click PARM;
int W_plot_all_time_courses PARM;
int W_read_plot_list  PARM;
int W_read_vertex_list  PARM;
int W_draw_cursor  PARM;
int W_draw_marked_vertices  PARM;
int W_mark_vertex  PARM;
int W_mark_translated_vertex  PARM;
int W_draw_all_cursor  PARM;
int W_draw_all_vertex_cursor  PARM;
int W_clear_all_vertex_cursor  PARM;
/* begin rkt */
int W_send_current_labels PARM;
int W_select_vertex_by_vno PARM;
int W_sclv_read_from_dotw  PARM;
int W_sclv_read_from_dotw_frame  PARM;
int W_sclv_read_from_volume  PARM;
int W_sclv_write_dotw PARM;
int W_sclv_smooth  PARM;
int W_sclv_set_overlay_alpha  PARM;
int W_sclv_unload_field  PARM;
int W_sclv_set_current_field  PARM;
int W_sclv_set_current_timepoint  PARM;
int W_sclv_copy_view_settings_from_current_field PARM;
int W_sclv_copy_all_view_settings_from_current_field PARM;
int W_sclv_copy_view_settings_from_field PARM;
int W_sclv_set_current_threshold_from_percentile PARM;
int W_sclv_set_current_threshold_using_fdr PARM;
int W_sclv_send_histogram PARM;
int W_sclv_send_current_field_info PARM;
int W_sclv_get_normalized_color_for_value PARM;
int W_clear_vertex_marks PARM;
int W_swap_vertex_fields PARM;
int W_undo_last_action PARM;
int W_get_marked_vnos PARM;
int W_get_selected_path_vnos PARM;
int W_save_tiff PARM;
int W_flip_normals PARM;
int W_mark_contiguous_vertices_over_thresh PARM;
int W_mark_contiguous_vertices_with_similar_curvature PARM;
int W_rip_all_vertices_except_contiguous_upripped PARM;
int W_func_calc_correlation_and_write_to_overlay PARM;
int W_func_normalize PARM;
int W_cptn_set_format_string PARM;
/* end rkt */

#define TkCreateMainWindow Tk_CreateMainWindow


/*=======================================================================*/
/* function wrappers and errors */
int                  W_swap_buffers
WBEGIN
ERR(1,"Wrong # args: swap_buffers")
swap_buffers();
WEND

int                  W_to_single_buffer  WBEGIN
ERR(1,"Wrong # args: to_single_buffer")
to_single_buffer();
WEND

int                  W_to_double_buffer  WBEGIN
ERR(1,"Wrong # args: to_double_buffer")
to_double_buffer();
WEND

int                  W_open_window  WBEGIN
ERR(1,"Wrong # args: open_window")
open_window(pname);
WEND

int                  W_help  WBEGIN
ERR(1,"Wrong # args: help")
print_help_tksurfer();
WEND

int                  W_redraw  WBEGIN
ERR(1,"Wrong # args: redraw")
redraw();
WEND

int                  W_redraw_second  WBEGIN
ERR(1,"Wrong # args: redraw_second")
redraw_second();
WEND

int                  W_shrink  WBEGIN
ERR(2,"Wrong # args: shrink <steps>")
shrink(atoi(argv[1]));
WEND

int                  W_area_shrink  WBEGIN
ERR(2,"Wrong # args: area_shrink <steps>")
area_shrink(atoi(argv[1]));
WEND

int                  W_sphere_shrink  WBEGIN
ERR(3,"Wrong # args: sphere_shrink <steps> <mm rad>")
sphere_shrink(atoi(argv[1]),atof(argv[2]));
WEND

int                 W_ellipsoid_project  WBEGIN
ERR(4,"Wrong # args: ellipsoid_project <mm:rh/lh> <ant/post> <sup/inf>")
ellipsoid_project(atof(argv[1]),atof(argv[2]),
                  atof(argv[3]));
WEND

int                 W_ellipsoid_morph  WBEGIN
ERR(5,"Wrong # args: ellipsoid_morph <steps> "
    "<mm:rh/lh> <ant/post> <sup/inf>")
ellipsoid_morph(atoi(argv[1]),atof(argv[2]),
                atof(argv[3]),atof(argv[4]));
WEND

int                 W_ellipsoid_shrink  WBEGIN
ERR(5,"Wrong # args: ellipsoid_shrink <steps> "
    "<mm:rh/lh> <ant/post> <sup/inf>")
ellipsoid_shrink(atoi(argv[1]),atof(argv[2]),
                 atof(argv[3]),atof(argv[4]));
WEND
int                  W_ellipsoid_shrink_bug  WBEGIN
ERR(4,"Wrong # args: ellipsoid_shrink_bug <steps> <mm rad> <mm len>")
ellipsoid_shrink_bug(atoi(argv[1]),atof(argv[2]),
                     atof(argv[3]));
WEND

int                  W_curv_shrink_to_fill  WBEGIN
ERR(2,"Wrong # args: curv_shrink_to_fill <steps>")
curv_shrink_to_fill(atoi(argv[1]));
WEND

int                  W_smooth_curvim  WBEGIN
ERR(2,"Wrong # args: smooth_curvim <steps>")
smooth_curvim(atoi(argv[1]));
WEND

int                  W_smooth_curv  WBEGIN
ERR(2,"Wrong # args: smooth_curv <steps>")
smooth_curv(atoi(argv[1]));
WEND

int                  W_smooth_val  WBEGIN
ERR(2,"Wrong # args: smooth_val <steps>")
smooth_val(atoi(argv[1]));
WEND

int                  W_smooth_val_sparse  WBEGIN
ERR(2,"Wrong # args: smooth_val_sparse <steps>")
smooth_val_sparse(atoi(argv[1]));
WEND

int                  W_smooth_curvim_sparse  WBEGIN
ERR(2,"Wrong # args: smooth_curvim_sparse <steps>")
smooth_curvim_sparse(atoi(argv[1]));
WEND

int                  W_smooth_fs  WBEGIN
ERR(2,"Wrong # args: smooth_fs <steps>")
smooth_fs(atoi(argv[1]));
WEND

int                  W_add_subject_to_average_curvim  WBEGIN
ERR(3,"Wrong # args: add_subject_to_average_curvim "
    "<subjname> <morphsubdir>")
add_subject_to_average_curvim(argv[1],argv[2]);
WEND

int                  W_read_curv_images WBEGIN
ERR(1,"Wrong # args: read_curv_images")
read_curv_images(cif2name);
WEND

int                  W_read_stds WBEGIN
ERR(2,"Wrong # args: read_stds")
read_stds(atoi(argv[1]));
WEND

int                  W_read_second_binary_surf  WBEGIN
ERR(1,"Wrong # args: read_second_binary_surf")
read_second_binary_surface(if2name);
WEND

int                  W_read_second_binary_curv  WBEGIN
ERR(1,"Wrong # args: read_second_binary_curv")
read_second_binary_curvature(cf2name);
WEND

int                  W_normalize_second_binary_curv  WBEGIN
ERR(1,"Wrong # args: normalize_second_binary_curv")
normalize_second_binary_curvature();
WEND

int                  W_curv_to_curvim WBEGIN
ERR(1,"Wrong # args: curv_to_curvim")
curv_to_curvim();
WEND

int                  W_second_surface_curv_to_curvim WBEGIN
ERR(1,"Wrong # args: second_surface_curv_to_curvim")
second_surface_curv_to_curvim();
WEND

int                  W_curvim_to_second_surface WBEGIN
ERR(1,"Wrong # args: curvim_to_second_surface")
curvim_to_second_surface();
WEND

int                  W_swap_curv WBEGIN
ERR(1,"Wrong # args: swap_curv")
swap_curv();
WEND

int                  W_curvim_to_surface WBEGIN
ERR(1,"Wrong # args: curvim_to_surface")
curvim_to_surface();
WEND

int                  W_read_binary_surf  WBEGIN
ERR(1,"Wrong # args: read_binary_surf")
read_binary_surface(ifname);
WEND

int                  W_read_surf  WBEGIN
ERR(2,"Wrong # args: read_surf")
read_positions(argv[1]);
WEND

int                  W_save_surf  WBEGIN
ERR(1,"Wrong # args: save_surf")
save_surf();
WEND

int                  W_restore_surf  WBEGIN
ERR(1,"Wrong # args: restore_surf")
restore_surf();
WEND

int                  W_show_surf  WBEGIN
ERR(2,"Wrong # args: show_surf")
show_surf(argv[1]);
WEND

int                  W_read_binary_curv  WBEGIN
ERR(1,"Wrong # args: read_binary_curv")
read_binary_curvature(cfname);
WEND

int                  W_read_binary_sulc  WBEGIN
ERR(1,"Wrong # args: read_binary_sulc")
read_binary_curvature(kfname);
WEND

int                  W_read_binary_values  WBEGIN
ERR(1,"Wrong # args: read_binary_values")
read_binary_values(vfname);
WEND

int                  W_read_binary_values_frame  WBEGIN
ERR(1,"Wrong # args: read_binary_values_frame")
read_binary_values_frame(vfname);
WEND

int                  W_read_annotations  WBEGIN
ERR(2,"Wrong # args: read_annotations")
read_annotations(argv[1]);
WEND

int                  W_read_annotated_image  WBEGIN
ERR(1,"Wrong # args: read_annotated_image")
read_annotated_image(nfname,frame_xdim,frame_ydim);
WEND

int                  W_read_binary_patch  WBEGIN
ERR(1,"Wrong # args: read_binary_patch")
read_binary_patch(pfname);
WEND

int                  W_read_fieldsign  WBEGIN
ERR(1,"Wrong # args: read_fieldsign")
read_fieldsign(fsfname);
WEND

int                  W_read_fsmask WBEGIN
ERR(1,"Wrong # args: read_fsmask")
read_fsmask(fmfname);
WEND

int                  W_write_binary_areas  WBEGIN
ERR(1,"Wrong # args: write_binary_areas")
write_binary_areas(afname);
WEND

int                  W_write_binary_surface  WBEGIN
ERR(1,"Wrong # args: write_binary_surface")
write_binary_surface(ofname);
WEND

int                  W_write_binary_curv  WBEGIN
ERR(1,"Wrong # args: write_binary_curv")
write_binary_curvature(cfname);
WEND

int                  W_write_binary_sulc  WBEGIN
ERR(1,"Wrong # args: write_binary_sulc")
write_binary_curvature(kfname);
WEND

int                  W_write_binary_values  WBEGIN
ERR(1,"Wrong # args: write_binary_values")
write_binary_values(vfname);
WEND

int                  W_write_binary_patch  WBEGIN
ERR(1,"Wrong # args: write_binary_patch")
write_binary_patch(pfname);
WEND

int                  W_write_labeled_vertices  WBEGIN
ERR(1,"Wrong # args: write_labeled_vertices")
write_labeled_vertices(lfname);
WEND

int                  W_read_labeled_vertices  WBEGIN
ERR(1,"Wrong # args: read_labeled_vertices")
read_labeled_vertices(lfname);
WEND

int                  W_read_and_color_labeled_vertices  WBEGIN
ERR(4,"Wrong # args: read_and_color_labeled_vertices")
read_and_color_labeled_vertices(atoi(argv[1]),
                                atoi(argv[2]),
                                atoi(argv[3]));
WEND

int                  W_write_fieldsign  WBEGIN
ERR(1,"Wrong # args: write_fieldsign")
write_fieldsign(fsfname);
WEND

int                  W_write_fsmask WBEGIN
ERR(1,"Wrong # args: write_fsmask")
write_fsmask(fmfname);
WEND

int                  W_write_vrml WBEGIN
ERR(2,"Wrong # args: write_vrml <mode:0=pts,1=gray,2=curv>")
write_vrml(vrfname,atoi(argv[1]));
WEND

int                  W_write_binary_dipoles WBEGIN
ERR(1,"Wrong # args: write_binary_dipoles")
write_binary_dipoles(dfname);
WEND

int                  W_write_binary_decimation WBEGIN
ERR(1,"Wrong # args: write_binary_decimation")
write_binary_decimation(sfname);
WEND

int                  W_write_dipoles WBEGIN
ERR(1,"Wrong # args: write_dipoles")
write_dipoles(dfname);
WEND

int                  W_write_decimation WBEGIN
ERR(1,"Wrong # args: write_decimation")
write_decimation(sfname);
WEND

int                  W_write_curv_images WBEGIN
ERR(1,"Wrong # args: write_curv_images")
write_images(curvim,cif2name);
WEND

int                  W_write_fill_images WBEGIN
ERR(1,"Wrong # args: write_fill_images")
write_images(fill,fifname);
WEND

int                  W_fill_second_surface  WBEGIN
ERR(1,"Wrong # args: fill_second_surface")
fill_second_surface();
WEND

int                  W_subsample_dist WBEGIN
ERR(2,"Wrong # args: subsample_dist <mm>")
subsample_dist(atoi(argv[1]));
WEND

int                  W_subsample_orient WBEGIN
ERR(2,"Wrong # args: subsample_orient <orientthresh>")
subsample_orient(atof(argv[1]));
WEND

int                  W_write_subsample WBEGIN
ERR(1,"Wrong # args: write_subsample")
write_subsample(sfname);
WEND

int                  W_write_binary_subsample WBEGIN
ERR(1,"Wrong # args: write_binary_subsample")
write_binary_subsample(sfname);
WEND

int                  W_compute_curvature  WBEGIN
ERR(1,"Wrong # args: compute_curvature")
compute_curvature();
WEND

int                  W_compute_CMF  WBEGIN
ERR(1,"Wrong # args: compute_CMF")
compute_CMF();
WEND

int                  W_compute_cortical_thickness WBEGIN
ERR(1,"Wrong # args: compute_cortical_thickness")
compute_cortical_thickness();
WEND

int                  W_clear_curvature  WBEGIN
ERR(1,"Wrong # args: clear_curvature")
clear_curvature();
WEND

int                  W_clear_vals  WBEGIN
ERR(1,"Wrong # args: clear_vals")
clear_vals();
WEND

int                  W_clear_ripflags  WBEGIN
ERR(1,"Wrong # args: clear_ripflags")
clear_ripflags();
WEND

int                  W_restore_ripflags  WBEGIN
ERR(2,"Wrong # args: restore_ripflags <mode>")
restore_ripflags(atoi(argv[1]));
WEND

int                  W_floodfill_marked_patch  WBEGIN
ERR(2,"Wrong # args: floodfill_marked_patch "
    "<0=cutborder,1=fthreshborder,2=curvfill>")
floodfill_marked_patch(atoi(argv[1]));
WEND

int                  W_rip_unmarked_vertices  WBEGIN
ERR(1,"Wrong # args: rip_unmarked_vertices ")
rip_unmarked_vertices();
WEND

int                  W_dilate_ripped  WBEGIN
ERR(1,"Wrong # args: dilate_ripped")
dilate_ripped();
WEND

int                  W_draw_curvature_line  WBEGIN
ERR(1,"Wrong # args: draw_curvature_line")
draw_curvature_line();
WEND

int                  W_twocond  WBEGIN
ERR(3,"Wrong # args: twocond <cond #0> <cond #1>")
twocond(atoi(argv[1]), atoi(argv[2]));
WEND

int                  W_cut_line  WBEGIN
ERR(2,"Wrong # args: cut_line <0=open,1=closed>")
cut_line(atoi(argv[1]));
WEND

int                  W_plot_marked  WBEGIN
ERR(2,"Wrong # args: plot_marked <file name>")
plot_marked(argv[1]);
WEND

int                  W_plot_curv  WBEGIN
ERR(2,"Wrong # args: plot_curv <0=open,1=closed>")
plot_curv(atoi(argv[1]));
WEND

int                  W_draw_fundus  WBEGIN
ERR(2,"Wrong # args: draw_fundus")
draw_fundus(path_selected_path);
WEND

int                  W_put_retinotopy_stats_in_vals  WBEGIN
ERR(1,"Wrong # args: put_retinotopy_stats_in_vals")
put_retinotopy_stats_in_vals();
WEND

int                  W_draw_vector  WBEGIN
ERR(2,"Wrong # args: draw_vector <file name>")
draw_vector(argv[1]);
WEND

int                  W_cut_plane  WBEGIN
ERR(1,"Wrong # args: cut_plane     "
    "[select 4 pnts: 3=>plane,1=>what to save]")
cut_plane();
WEND

int                  W_flatten  WBEGIN
ERR(1,"Wrong # args: flatten")
flatten(tfname);
WEND

int                  W_normalize_binary_curv  WBEGIN
ERR(1,"Wrong # args: normalize_binary_curv")
normalize_binary_curvature();
WEND

int                  W_normalize_area  WBEGIN
ERR(1,"Wrong # args: normalize_area")
normalize_area();
WEND

int                  W_normalize_curvature  WBEGIN
ERR(2,"Wrong # args: normalize_curvature")
normalize_curvature(atoi(argv[1]));
WEND

int                  W_shift_values  WBEGIN
ERR(1,"Wrong # args: shift_values")
shift_values();
WEND

int                  W_swap_values  WBEGIN
ERR(1,"Wrong # args: swap_values")
swap_values();
WEND

int                  W_swap_stat_val  WBEGIN
ERR(1,"Wrong # args: swap_stat_val")
swap_stat_val();
WEND

int                  W_swap_val_val2  WBEGIN
ERR(1,"Wrong # args: swap_val_val2")
swap_val_val2();
WEND

int                  W_compute_angles  WBEGIN
ERR(1,"Wrong # args: compute_angles")
compute_angles();
WEND

int                  W_compute_fieldsign  WBEGIN
ERR(1,"Wrong # args: compute_fieldsign")
compute_fieldsign();
WEND

int                  W_draw_radius WBEGIN
ERR(1,"Wrong # args: draw_radius")
draw_spokes(RADIUS);
WEND

int                  W_draw_theta WBEGIN
ERR(1,"Wrong # args: draw_theta")
draw_spokes(THETA);
WEND

int                  W_save_rgb  WBEGIN
ERR(1,"Wrong # args: save_rgb")
save_rgb(agfname);
WEND

int                  W_save_rgb_num  WBEGIN
ERR(1,"Wrong # args: save_rgb_num")
save_rgb_num(gfname);
WEND

int                  W_save_rgb_named_orig  WBEGIN
ERR(2,"Wrong # args: save_rgb_named_orig <relfilename>")
save_rgb_named_orig(sgfname,argv[1]);
WEND

int                  W_save_rgb_cmp_frame  WBEGIN
ERR(1,"Wrong # args: save_rgb_cmp_frame")
save_rgb_cmp_frame(gfname,ilat);
WEND

int                  W_open_rgb_cmp_named  WBEGIN
ERR(2,"Wrong # args: open_rgb_cmp_named <filename>")
open_rgb_cmp_named(sgfname,argv[1]);
WEND

int                  W_save_rgb_cmp_frame_named  WBEGIN
ERR(2,"Wrong # args: save_rgb_cmp_frame_named <framelatency>")
save_rgb_cmp_frame_named(atof(argv[1]));
WEND

int                  W_close_rgb_cmp_named  WBEGIN
ERR(1,"Wrong # args: close_rgb_cmp_named")
close_rgb_cmp_named();
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
translate_brain(atof(argv[1]),0.0,0.0);
WEND

int                  W_translate_brain_y  WBEGIN
ERR(2,"Wrong # args: translate_brain_y <mm>")
translate_brain(0.0,atof(argv[1]),0.0);
WEND

int                  W_translate_brain_z  WBEGIN
ERR(2,"Wrong # args: translate_brain_z <mm>")
translate_brain(0.0,0.0,atof(argv[1]));
WEND

int                  W_resize_brain  WBEGIN
ERR(2,"Wrong # args: resize_brain <mm>")
resize_brain(atof(argv[1]));
WEND

int                  W_resize_window  WBEGIN
ERR(2,"Wrong # args: resize_window <pix>")
resize_window(atoi(argv[1]));
WEND

int                  W_setsize_window  WBEGIN
ERR(2,"Wrong # args: setsize_window <pix>")
setsize_window(atoi(argv[1]));
WEND

int                  W_move_window  WBEGIN
ERR(3,"Wrong # args: move_window <x> <y>")
move_window(atoi(argv[1]),atoi(argv[2]));
WEND

int                  W_do_lighting_model  WBEGIN
ERR(6,"Wrong # args: do_lighting_model "
    "<lit0> <lit1> <lit2> <lit3> <offset>")
do_lighting_model(atof(argv[1]),atof(argv[2]),
                  atof(argv[3]),atof(argv[4]),atof(argv[5]));
WEND

int                  W_restore_zero_position  WBEGIN
ERR(1,"Wrong # args: restore_zero_position")
restore_zero_position();
WEND

int                  W_restore_initial_position  WBEGIN
ERR(1,"Wrong # args: restore_initial_position")
restore_initial_position();
WEND

int                  W_make_lateral_view  WBEGIN
ERR(1,"Wrong # args: make_lateral_view")
make_lateral_view(stem);
WEND

int                  W_make_lateral_view_second  WBEGIN
ERR(1,"Wrong # args: make_lateral_view_second")
make_lateral_view_second(stem);
WEND

int                  W_read_view_matrix  WBEGIN
ERR(1,"Wrong # args: read_view_matrix")
read_view_matrix(tfname);
WEND

int                  W_write_view_matrix  WBEGIN
ERR(1,"Wrong # args: write_view_matrix")
write_view_matrix(tfname);
WEND

int                  W_read_really_matrix  WBEGIN
ERR(1,"Wrong # args: read_really_matrix")
read_really_matrix(tfname);
WEND

int                  W_write_really_matrix  WBEGIN
ERR(1,"Wrong # args: write_really_matrix")
write_really_matrix(tfname);
WEND

int                  W_really_translate_brain  WBEGIN
ERR(4,"Wrong # args: really_translate_brain "
    "<mm:rh/lh> <ant/post> <sup/inf>")
really_translate_brain(atof(argv[1]),atof(argv[2]),
                       atof(argv[3]));
WEND
int                  W_really_scale_brain  WBEGIN
ERR(4,"Wrong # args: really_scale_brain <rh/lh> <ant/post> <sup/inf>")
really_scale_brain(atof(argv[1]),atof(argv[2]),
                   atof(argv[3]));
WEND
int                  W_really_rotate_brain_x  WBEGIN
ERR(2,"Wrong # args: really_rotate_brain_x <deg>    [N.B.: lh/rh axis]")
really_rotate_brain(atof(argv[1]),'x');
WEND

int                  W_align_sphere  WBEGIN
ERR(1,"Wrong # args: align_sphere <deg>    [N.B.: lh/rh axis]")
align_sphere(mris);
WEND

int                  W_really_rotate_brain_y  WBEGIN
ERR(2,"Wrong # args: really_rotate_brain_y <deg>    "
    "[N.B.: ant/post axis]")
really_rotate_brain(atof(argv[1]),'y');
WEND

int                  W_really_rotate_brain_z  WBEGIN
ERR(2,"Wrong # args: really_rotate_brain_z <deg>    [N.B.: sup/inf axis]")
really_rotate_brain(atof(argv[1]),'z');
WEND

int                  W_really_center_brain  WBEGIN
ERR(1,"Wrong # args: really_center_brain")
really_center_brain();
WEND

int                  W_really_center_second_brain  WBEGIN
ERR(1,"Wrong # args: really_center_second_brain")
really_center_second_brain();
WEND

int                  W_really_align_brain  WBEGIN
ERR(1,"Wrong # args: really_align_brain")
really_align_brain();
WEND

int                  W_read_binary_decimation  WBEGIN
ERR(1,"Wrong # args: read_binary_decimation")
read_binary_decimation(sfname);
WEND

int                  W_read_binary_dipoles  WBEGIN
ERR(1,"Wrong # args: read_binary_dipoles")
read_binary_dipoles(dfname);
WEND

int                  W_load_vals_from_sol  WBEGIN
ERR(4,"Wrong # args: load_vals_from_sol <tmid> <dt> <cond_num>")
load_vals_from_sol(atof(argv[1]),atof(argv[2]),
                   atoi(argv[3]));
WEND

int                  W_load_var_from_sol  WBEGIN
ERR(2,"Wrong # args: load_var_from_sol <cond_num>")
load_var_from_sol(atoi(argv[1]));
WEND

int                  W_read_iop  WBEGIN
ERR(3,"Wrong # args: read_iop <iop_name> <hemi_num>")
read_iop(argv[1],atoi(argv[2]));
WEND

int                  W_read_rec  WBEGIN
ERR(2,"Wrong # args: read_rec <rec_file>")
read_rec(argv[1]);
WEND

int                  W_read_ncov  WBEGIN
ERR(2,"Wrong # args: read_ncov <ncov_file>")
read_ncov(argv[1]);
WEND

int                  W_filter_recs  WBEGIN
ERR(1,"Wrong # args: filter_recs")
filter_recs();
WEND

int                  W_normalize_time_courses  WBEGIN
ERR(2,"Wrong # args: normalize_time_courses <type>")
normalize_time_courses(atoi(argv[1]));
WEND

int                  W_compute_timecourses  WBEGIN
ERR(1,"Wrong # args: compute_timecourses")
compute_timecourses();
WEND

int                  W_normalize_inverse  WBEGIN
ERR(1,"Wrong # args: normalize_inverse")
normalize_inverse();
WEND

int                  W_compute_pval_fwd  WBEGIN
ERR(2,"Wrong # args: compute_pval_fwd <pvalthresh>")
compute_pval_fwd(atof(argv[1]));
WEND

int                  W_compute_select_fwd  WBEGIN
ERR(2,"Wrong # args: compute_select_fwd <maxdist (mm)>")
compute_select_fwd(atof(argv[1]));
WEND

int                  W_compute_pval_inv  WBEGIN
ERR(1,"Wrong # args: compute_pval_inv")
compute_pval_inv();
WEND

int                  W_find_orig_vertex_coordinates  WBEGIN
ERR(1,"Wrong # args: find_orig_vertex_coordinates")
find_orig_vertex_coordinates(selection);
WEND

int                  W_select_orig_vertex_coordinates  WBEGIN
ERR(1,"Wrong # args: select_orig_vertex_coordinates")
select_orig_vertex_coordinates(&selection);
WEND

int                  W_select_talairach_point  WBEGIN
ERR(4,"Wrong # args: select_talairach_point <xtal> <ytal> <ztal>")
select_talairach_point(&selection,
                       atof(argv[1]),atof(argv[2]),atof(argv[3]));
WEND

int                  W_print_nearest_vertex_to_talairach_point  WBEGIN
ERR(4,"Wrong # args: print_nearest_vertex_to_talairach_point <xtal> <ytal> <ztal>")
print_nearest_vertex_to_talairach_point(
                       atof(argv[1]),atof(argv[2]),atof(argv[3]));
WEND

int                  W_read_orig_vertex_coordinates  WBEGIN
ERR(1,"Wrong # args: read_orig_vertex_coordinates")
read_orig_vertex_coordinates(orfname);
WEND

int                  W_read_white_vertex_coordinates  WBEGIN
ERR(1,"Wrong # args: read_white_vertex_coordinates")
read_white_vertex_coordinates();
WEND

int                  W_read_pial_vertex_coordinates  WBEGIN
ERR(1,"Wrong # args: read_pial_vertex_coordinates")
read_pial_vertex_coordinates();
WEND

int                  W_read_canon_vertex_coordinates  WBEGIN
ERR(2,"Wrong # args: read_canon_vertex_coordinates")
read_canon_vertex_coordinates(argv[1]);
WEND

int                  W_send_spherical_point  WBEGIN
ERR(4,"Wrong # args: send_spherical_point <subject> <sphere> <orig>")
send_spherical_point(argv[1], argv[2], argv[3]);
WEND
int                  W_send_contralateral_point  WBEGIN

ERR(3,"Wrong # args: send_contralateral_point <left/right surface> <orig surface>")
send_contralateral_point(argv[1], argv[2]);
WEND

int                  W_send_to_subject  WBEGIN
ERR(2,"Wrong # args: send_to_subject <subject name>")
send_to_subject(argv[1]);
WEND

int                  W_send_to_other_hemi  WBEGIN
ERR(1,"Wrong # args: send_to_other_hemi")
send_to_other_hemi();
WEND

int                  W_resend_to_subject  WBEGIN
ERR(1,"Wrong # args: resend_to_subject")
resend_to_subject();
WEND

int                  W_drawcb  WBEGIN
ERR(1,"Wrong # args: drawcb")
drawcb();
WEND

int                  W_read_ellipsoid_vertex_coordinates  WBEGIN
ERR(4,"Wrong # args: read_ellipsoid_vertex_coordinates <a> <b> <c>")
read_ellipsoid_vertex_coordinates(elfname,
                                  atof(argv[1]),
                                  atof(argv[2]),
                                  atof(argv[3]));
WEND

int                  W_invert_surface  WBEGIN
ERR(1,"Wrong # args: invert_surface ")
invert_surface();
WEND

int                  W_fix_nonzero_vals  WBEGIN
ERR(1,"Wrong # args: fix_nonzero_vals ")
fix_nonzero_vals();
WEND

int                  W_invert_vertex  WBEGIN
ERR(2,"Wrong # args: invert_vertex ")
invert_vertex(atoi(argv[1]));
WEND

int                  W_invert_face  WBEGIN
ERR(2,"Wrong # args: invert_face ")
invert_face(atoi(argv[1]));
WEND

int                  W_mark_annotation  WBEGIN
ERR(1,"Wrong # args: mark_annotation ")
mark_annotation(selection);
WEND

int                  W_mark_faces  WBEGIN
ERR(2,"Wrong # args: mark_faces ")
mark_faces(atoi(argv[1]));
WEND

int                  W_mark_face  WBEGIN
ERR(2,"Wrong # args: mark_face ")
mark_face(atoi(argv[1]));
WEND

int                  W_dump_vertex  WBEGIN
ERR(2,"Wrong # args: dump_vertex ")
dump_vertex(atoi(argv[1]));
WEND

int                  W_val_to_mark  WBEGIN
ERR(1,"Wrong # args: val_to_mark ")
val_to_mark();
WEND

int                  W_set_area_thresh  WBEGIN
ERR(2,"Wrong # args: set_area_thresh ")
set_area_thresh(atof(argv[1]));
WEND


int                  W_scale_brain  WBEGIN
ERR(2,"Wrong # args: scale_brain ")
scale_brain(atof(argv[1]));
WEND

int                  W_transform_brain  WBEGIN
ERR(1,"Wrong # args: transform_brain ")
transform_brain();
WEND

int                  W_show_flat_regions  WBEGIN
ERR(3,"Wrong # args: show_flat_regions ")
show_flat_regions(argv[1], atof(argv[2]));
WEND

int                  W_val_to_curv  WBEGIN
ERR(1,"Wrong # args: val_to_curv ")
val_to_curv();
WEND

int                  W_val_to_stat  WBEGIN
ERR(1,"Wrong # args: val_to_stat ")
val_to_stat();
WEND

int                  W_set_vals  WBEGIN
ERR(2,"Wrong # args: set_vals ")
MRISsetVals(mris,atof(argv[2]));
WEND

int                  W_stat_to_val  WBEGIN
ERR(1,"Wrong # args: stat_to_val ")
stat_to_val();
WEND


int                  W_f_to_p  WBEGIN
ERR(3,"Wrong # args: f_to_p(numer_dof, denom_dof) ")
f_to_p(atoi(argv[1]), atoi(argv[2]));
WEND

int                  W_t_to_p  WBEGIN
ERR(2,"Wrong # args: t_to_p(dof) ")
t_to_p(atoi(argv[1]));
WEND

int                  W_f_to_t  WBEGIN
ERR(1,"Wrong # args: f_to_t ")
f_to_t();
WEND

int                  W_label_to_stat  WBEGIN
ERR(2,"Wrong # args: label_to_stat ")
label_to_stat(atoi(argv[1]));
WEND

int                  W_taubin_smooth  WBEGIN
ERR(4,"Wrong # args: taubin_smooth <niter> <lambda> <mu>  where  mu < -lambda < 0 ")
  taubin_smooth(mris, atoi(argv[1]), atof(argv[2]), atof(argv[3]));
WEND

int                  W_label_from_stats  WBEGIN
ERR(2,"Wrong # args: label_from_stats ")
label_from_stats(atoi(argv[1]));
WEND

int                  W_label_set_stats  WBEGIN
ERR(2,"Wrong # args: label_set_stats ")
label_set_stats(atof(argv[1]));
WEND

int                  W_remove_triangle_links  WBEGIN
ERR(1,"Wrong # args: remove_triangle_links ")
remove_triangle_links();
WEND

int                  W_read_soltimecourse  WBEGIN
ERR(2,"Wrong # args: read_soltimecourse ")
read_soltimecourse(argv[1]);
WEND

int                  W_read_imag_vals  WBEGIN
ERR(2,"Wrong # args: read_imag_vals ")
read_imag_vals(argv[1]);
WEND

int                  W_sol_plot  WBEGIN
ERR(3,"Wrong # args: sol_plot ")
sol_plot(atoi(argv[1]), atoi(argv[2]));
WEND

int                  W_read_curv_to_val  WBEGIN
ERR(2,"Wrong # args: read_curv_to_val ")
read_curv_to_val(argv[1]);
WEND

int                  W_read_and_smooth_parcellation  WBEGIN
ERR(5,"Wrong # args: read_and_smooth_parcellation ")
read_and_smooth_parcellation(argv[1],argv[2],
                             atoi(argv[3]),
                             atoi(argv[4]));
WEND

int                  W_read_parcellation  WBEGIN
ERR(3,"Wrong # args: read_parcellation ")
read_parcellation(argv[1],argv[2]);
WEND

int                  W_read_disc  WBEGIN
ERR(2,"Wrong # args: read_disc ")
read_disc(argv[1]);
WEND

int                  W_deconvolve_weights  WBEGIN
ERR(3,"Wrong # args: deconvolve_weights(weight_fname, scale_fname) ")
deconvolve_weights(argv[1], argv[2]);
WEND

int                  W_curv_to_val  WBEGIN
ERR(1,"Wrong # args: curv_to_val ")
curv_to_val();
WEND

int                  W_scale_vals  WBEGIN
ERR(2,"Wrong # args: scale_vals ")
     scale_vals(atof(argv[1]));
WEND

int                  W_mask_label  WBEGIN
ERR(2,"Wrong # args: mask_label ")
mask_label(argv[1]);
WEND

int                  W_dump_faces  WBEGIN
ERR(2,"Wrong # args: dump_faces ")
dump_faces(atoi(argv[1]));
WEND

int                  W_load_gcsa  WBEGIN
ERR(2,"Wrong # args: load_gcsa ")
load_gcsa(argv[1]);
WEND

int                  W_orient_sphere  WBEGIN
ERR(1,"Wrong # args: orient_sphere ")
orient_sphere();
WEND

int                  W_draw_ellipsoid_latlong  WBEGIN
ERR(4,"Wrong # args: draw_ellipsoid_latlong <a> <b> <c>")
draw_ellipsoid_latlong(atof(argv[1]),atof(argv[2]),
                       atof(argv[3]));
WEND
int                  W_left_click WBEGIN
ERR(3,"Wrong # args: left_click <x> <y>      [relative to window 0,0]")
left_click(
  (short)atoi(argv[1]),(short)atoi(argv[2]));

WEND

int                  W_plot_all_time_courses WBEGIN
ERR(1,"Wrong # args: plot_all_time_courses")
plot_all_time_courses();
WEND

int                  W_read_plot_list  WBEGIN
ERR(2,"Wrong # args: read_plot_list <vlist_file> ")
read_plot_list(argv[1]);
WEND

int                  W_read_vertex_list  WBEGIN
ERR(2,"Wrong # args: read_vertex_list <vlist_file> ")
read_vertex_list(argv[1]);
WEND

int                  W_draw_cursor  WBEGIN
ERR(3,"Wrong # args: draw_cursor <vertex_number> <on_or_off>")
draw_cursor(atoi(argv[1]),atoi(argv[2]));
WEND

int                  W_draw_marked_vertices  WBEGIN
ERR(1,"Wrong # args: draw_marked_vertices")
draw_marked_vertices();
WEND

int                  W_mark_vertex  WBEGIN
ERR(3,"Wrong # args: mark_vertex <vertex_number> <on_or_off>")
mark_vertex(atoi(argv[1]),atoi(argv[2]));
WEND

int                  W_mark_translated_vertex  WBEGIN
ERR(4,"Wrong # args: mark_translated_vertex <vertex_number> <on_or_off> <reg surf>")
mark_translated_vertex(atoi(argv[1]),atoi(argv[2]), argv[3]);
WEND

int                  W_draw_all_cursor  WBEGIN
ERR(1,"Wrong # args: draw_all_cursor")
draw_all_cursor();
WEND

int                  W_draw_all_vertex_cursor  WBEGIN
ERR(1,"Wrong # args: draw_all_vertex_cursor")
draw_all_vertex_cursor();
WEND

int                  W_clear_all_vertex_cursor  WBEGIN
ERR(1,"Wrong # args: clear_all_vertex_cursor")
clear_all_vertex_cursor();
WEND
/* begin rkt */
int W_send_current_labels WBEGIN
ERR(1, "Wrong # args: send_current_labels")
send_current_labels ();
WEND

int W_select_vertex_by_vno WBEGIN
ERR(2, "Wrong # args: select_vertex_by_vno vno")
select_vertex_by_vno (atoi(argv[1]));
WEND

int W_swap_vertex_fields WBEGIN
ERR(3,"Wrong # args: swap_vertex_fields <typea> <typeb>")
swap_vertex_fields(atoi(argv[1]),atoi(argv[2]));
WEND
int W_clear_vertex_marks WBEGIN
ERR(1,"Wrong # args: clear_vertex_marks")
clear_vertex_marks();
WEND
int W_clear_all_vertex_marks WBEGIN
ERR(1,"Wrong # args: clear_all_vertex_marks")
clear_all_vertex_marks();
WEND

int W_close_marked_vertices WBEGIN
ERR(1,"Wrong # args: close_marked_vertices")
close_marked_vertices();
WEND

int W_undo_last_action WBEGIN
ERR(1,"Wrong # args: undo_last_action")
undo_do_first_action();
WEND

int                  W_sclv_read_from_dotw  WBEGIN
ERR(2,"Wrong # args: sclv_read_from_dotw field")
sclv_read_from_dotw(vfname,atoi(argv[1]));
WEND

int                  W_sclv_read_from_dotw_frame  WBEGIN
ERR(2,"Wrong # args: sclv_read_from_dotw_frame field")
sclv_read_from_dotw_frame(vfname,atoi(argv[1]));
WEND

int                  W_sclv_write_dotw  WBEGIN
ERR(2,"Wrong # args: sclv_write_dotw field")
sclv_write_dotw(vfname,atoi(argv[1]));
WEND

int                  W_read_surface_vertex_set  WBEGIN
ERR(3,"Wrong # args: read_surface_vertex_set field file")
vset_read_vertex_set(atoi(argv[1]), argv[2]);
WEND

int                  W_set_current_vertex_set  WBEGIN
ERR(2,"Wrong # args: set_current_vertex_set field")
vset_set_current_set(atoi(argv[1]));
WEND

int                  W_sclv_load_label_value_file  WBEGIN
ERR(3,"Wrong # args: sclv_load_label_value_file fname field")
sclv_load_label_value_file(argv[1],atoi(argv[2]));
WEND


int W_func_load_timecourse (ClientData clientData,Tcl_Interp *interp,
                            int argc,char *argv[])
{
  FunD_tRegistrationType reg_type = FunD_tRegistration_None;

  if (argc!=3 && argc!=4)
  {
    Tcl_SetResult(interp,"Wrong # args: func_load_timecourse volumeFileName"
                  "registrationType [registrationFileName]",TCL_VOLATILE);
    return TCL_ERROR;
  }

  reg_type = (FunD_tRegistrationType) atoi (argv[2]);

  if (argc==3)
    func_load_timecourse (argv[1], reg_type, NULL);
  /* even if we have 3 args, tcl could have passed us a blank string
     for the 4th. if it's blank, call the load function with a null
     registration, otherwise it will think it's a valid file name. */
  if (argc==4)
  {
    if (strcmp(argv[3],"")==0)
      func_load_timecourse (argv[1], reg_type, NULL);
    else
      func_load_timecourse (argv[1], reg_type, argv[2]);
  }
  return TCL_OK;
}

int W_func_load_timecourse_offset (ClientData clientData,Tcl_Interp *interp,
                                   int argc,char *argv[])
{
  FunD_tRegistrationType reg_type = FunD_tRegistration_None;

  if (argc!=3 && argc!=4)
  {
    Tcl_SetResult(interp,"Wrong # args: func_load_timecourse_offset fname"
                  "registrationType [registration]",TCL_VOLATILE);
    return TCL_ERROR;
  }

  reg_type = (FunD_tRegistrationType) atoi (argv[2]);

  if (argc==3)
    func_load_timecourse_offset (argv[1], reg_type, NULL);
  /* even if we have 3 args, tcl could have passed us a blank string
     for the 4th. if it's blank, call the load function with a null
     registration, otherwise it will think it's a valid file name. */
  if (argc==4)
  {
    if (strcmp(argv[3],"")==0)
      func_load_timecourse_offset (argv[1], reg_type, NULL);
    else
      func_load_timecourse_offset (argv[1], reg_type, argv[2]);
  }

  return TCL_OK;
}

int W_sclv_read_from_volume (ClientData clientData,Tcl_Interp *interp,
                             int argc,char *argv[])
{
  FunD_tRegistrationType reg_type = FunD_tRegistration_None;

  if (argc!=4&&argc!=5)
  {
    Tcl_SetResult(interp,"Wrong # args: sclv_read_from_volume field fname "
                  "registrationType [registration]",TCL_VOLATILE);
    return TCL_ERROR;
  }

  reg_type = (FunD_tRegistrationType) atoi (argv[3]);

  if (argc==4)
    sclv_read_from_volume (argv[2], reg_type, NULL, atoi(argv[1]));
  /* even if we have 4 args, tcl could have passed us a blank string
     for the 4th. if it's blank, call the read function with a null
     registration, otherwise it will think it's a valid file name. */
  if (argc==5)
  {
    if (strcmp(argv[4],"")==0)
      sclv_read_from_volume (argv[2], reg_type, NULL, atoi(argv[1]));
    else
      sclv_read_from_volume (argv[2], reg_type, argv[4], atoi(argv[1]));
  }
  return TCL_OK;
}
int W_sclv_smooth  WBEGIN
ERR(3,"Wrong # args: sclv_smooth steps field")
sclv_smooth(atoi(argv[1]), atoi(argv[2]));
WEND
int W_sclv_set_overlay_alpha  WBEGIN
ERR(2,"Wrong # args: sclv_set_overlay_alpha alpha")
sclv_set_overlay_alpha(atof(argv[1]));
WEND
int W_sclv_set_current_field  WBEGIN
ERR(2,"Wrong # args: sclv_set_current_field field")
sclv_set_current_field(atoi(argv[1]));
WEND
int W_sclv_unload_field  WBEGIN
ERR(2,"Wrong # args: sclv_unload_field field")
sclv_unload_field(atoi(argv[1]));
WEND
int W_sclv_set_current_timepoint  WBEGIN
ERR(3,"Wrong # args: sclv_set_current_timepoint timepoint condition")
sclv_set_timepoint_of_field(sclv_current_field,
                            atoi(argv[1]),atoi(argv[2]));
WEND
int W_sclv_copy_view_settings_from_field  WBEGIN
ERR(3,"Wrong # args: sclv_copy_view_settings_from_field field fromfield")
sclv_copy_view_settings_from_field(atoi(argv[1]),atoi(argv[2]));
WEND
int W_sclv_copy_view_settings_from_current_field  WBEGIN
ERR(2,"Wrong # args: sclv_copy_view_settings_from_current_field field")
sclv_copy_view_settings_from_current_field(atoi(argv[1]));
WEND
int W_sclv_copy_all_view_settings_from_current_field  WBEGIN
ERR(1,"Wrong # args: sclv_copy_all_view_settings_from_current_field")
sclv_copy_all_view_settings_from_current_field();
WEND
int W_sclv_set_current_threshold_from_percentile  WBEGIN
ERR(4,"Wrong # args: sclv_set_current_threshold_from_percentile "
    "thresh mid slope")
sclv_set_current_threshold_from_percentile(
  atof(argv[1]),
  atof(argv[2]),
  atof(argv[3]));
WEND
int W_sclv_set_current_threshold_using_fdr  WBEGIN
ERR(3,"Wrong # args: sclv_set_current_threshold_using_fdr rate marked")
sclv_set_threshold_using_fdr(sclv_current_field,atof(
                               argv[1]),
                             atoi(argv[2]));
WEND
int W_sclv_send_histogram  WBEGIN
ERR(2,"Wrong # args: sclv_send_histogram field")
sclv_send_histogram(atoi(argv[1]));
WEND
int W_sclv_send_current_field_info  WBEGIN
ERR(1,"Wrong # args: sclv_send_current_field_info")
sclv_send_current_field_info();
WEND


int W_sclv_get_normalized_color_for_value (ClientData clientData,
    Tcl_Interp *interp,
    int argc,char *argv[])
{
  float r, g, b;
  char stringResult[256];

  if (argc!=2)
  {
    Tcl_SetResult(interp,"Wrong # args: sclv_get_normalized_color_for_value "
                  "value",TCL_VOLATILE);
    return TCL_ERROR;
  }

  /* get the rgb values here and print them to a string */
  sclv_get_normalized_color_for_value (sclv_current_field,
                                       atof(argv[1]), &r, &g, &b);
  sprintf (stringResult, "%f %f %f", r, g, b);

  /* set tcl result, volatile so tcl will make a copy of it. */
  Tcl_SetResult( interp, stringResult, TCL_VOLATILE );

  return TCL_OK;
}


int W_func_select_marked_vertices  WBEGIN
ERR(1,"Wrong # args: func_select_marked_vertices")
func_select_marked_vertices();
WEND
int W_func_select_selected_vertex  WBEGIN
ERR(1,"Wrong # args: func_select_selected_vertex")
func_select_selected_vertex();
WEND
int W_func_select_label  WBEGIN
ERR(1,"Wrong # args: func_select_label")
func_select_label();
WEND
int W_func_clear_selection  WBEGIN
ERR(1,"Wrong # args: func_clear_selection")
func_clear_selection();
WEND
int W_func_graph_timecourse_selection  WBEGIN
ERR(1,"Wrong # args: func_graph_timecourse_selection")
func_graph_timecourse_selection();
WEND
int W_func_print_timecourse_selection  WBEGIN
ERR(2,"Wrong # args: func_print_timecourse_selection filename")
func_print_timecourse_selection(argv[1]);
WEND


int W_labl_load_color_table WBEGIN
ERR(2,"Wrong # args: labl_load_color_table filename")
labl_load_color_table (argv[1]);
WEND
int W_labl_load WBEGIN
ERR(2,"Wrong # args: labl_load filename")
labl_load (argv[1]);
WEND
int W_labl_save WBEGIN
ERR(3,"Wrong # args: labl_save index filename")
labl_save (atoi(argv[1]), argv[2]);
WEND
int W_labl_save_all WBEGIN
ERR(1,"Wrong # args: labl_save_all prefix")
labl_save_all (argv[1]);
WEND
int W_labl_import_annotation WBEGIN
ERR(2,"Wrong # args: labl_import_annotation filename")
labl_import_annotation (argv[1]);
WEND
int W_labl_export_annotation WBEGIN
ERR(2,"Wrong # args: labl_export_annotation filename")
labl_export_annotation (argv[1]);
WEND
int W_labl_new_from_marked_vertices WBEGIN
ERR(1,"Wrong # args: labl_new_from_marked_vertices")
labl_new_from_marked_vertices (NULL);
WEND
int W_labl_mark_vertices WBEGIN
ERR(2,"Wrong # args: labl_mark_vertices index")
labl_mark_vertices (atoi(argv[1]));
WEND
int W_labl_select WBEGIN
ERR(2,"Wrong # args: labl_select index")
labl_select (atoi(argv[1]));
WEND
int W_labl_set_name_from_table WBEGIN
ERR(2,"Wrong # args: labl_set_name_from_table index")
labl_set_name_from_table (atoi(argv[1]));
WEND
int W_labl_set_info WBEGIN
ERR(8,"Wrong # args: labl_set_info index name structure visibility")
labl_set_info (atoi(argv[1]),
               argv[2], atoi(argv[3]),
               atoi(argv[4]), atoi(argv[5]),
               atoi(argv[6]), atoi(argv[7]) );
WEND
int W_labl_set_color WBEGIN
ERR(5,"Wrong # args: labl_set_color index r g b")
labl_set_color (atoi(argv[1]), atoi(argv[2]),
                atoi(argv[3]), atoi(argv[4]) );
WEND
int W_labl_remove WBEGIN
ERR(2,"Wrong # args: labl_remove index")
labl_remove (atoi(argv[1]));
WEND
int W_labl_threshold WBEGIN
ERR(3,"Wrong # args: labl_threshold index threshold")
labl_threshold (atoi(argv[1]), atof(argv[2]));
WEND
int W_labl_remove_all WBEGIN
ERR(1,"Wrong # args: labl_remove_all")
labl_remove_all ();
WEND
int W_labl_select_label_by_vno WBEGIN
ERR(2,"Wrong # args: labl_select_label_by_vno vno")
labl_select_label_by_vno (atoi(argv[1]));
WEND
int W_labl_erode WBEGIN
ERR(2,"Wrong # args: labl_erode index")
labl_erode (atoi(argv[1]));
WEND
int W_labl_dilate WBEGIN
ERR(2,"Wrong # args: labl_dilate index")
labl_dilate (atoi(argv[1]));
WEND
int W_labl_fill_holes WBEGIN
ERR(2,"Wrong # args: labl_fill_holes index")
labl_fill_holes (atoi(argv[1]));
WEND
int W_labl_print_list WBEGIN
ERR(1,"Wrong # args: labl_print_list")
labl_print_list();
WEND
int W_labl_print_table WBEGIN
ERR(1,"Wrong # args: labl_print_table")
labl_print_table();
WEND
int W_path_select WBEGIN
ERR(2,"Wrong # args: path_select")
path_select( atoi(argv[1]));
WEND
int W_path_new_path_from_marked_vertices WBEGIN
ERR(1,"Wrong # args: path_new_path_from_marked_vertices")
path_new_path_from_marked_vertices();
WEND
int W_path_remove_selected_path WBEGIN
ERR(1,"Wrong # args: path_remove_selected_path")
path_remove_selected_path();
WEND
int W_path_mark_selected_path WBEGIN
ERR(1,"Wrong # args: path_mark_selected_path")
path_mark_selected_path();
WEND
int W_path_save WBEGIN
ERR(2,"Wrong # args: path_save fname")
path_save(argv[1]);
WEND
int W_path_load WBEGIN
ERR(2,"Wrong # args: path_load fname")
path_load(argv[1]);
WEND
int W_edit_vertex_at_cursor WBEGIN
ERR(3,"Wrong # args: edit_vertex_at_cursor action argument")
edit_vertex_at_cursor(atoi(argv[1]), atoi(argv[2]));
WEND

int W_fill_flood_from_cursor (ClientData clientData,Tcl_Interp *interp,
                              int argc,char *argv[])
{
  FILL_PARAMETERS params;
  char first_fill;
  int seed;
  int seeds[mris->nvertices];
  int nseeds;
  int n;

  if (argc != 9)
  {
    Tcl_SetResult(interp,"Wrong # args: fill_flood_from_cursor "
                  "dont_cross_path dont_cross_label dont_fill_unlabeled "
                  "dont_cross_cmid dont_cross_fthresh  "
                  "use_multiple_seeds action argument",
                  TCL_VOLATILE);
    return TCL_ERROR;
  }

  params.dont_cross_path     = atoi(argv[1]);
  params.dont_cross_label    = atoi(argv[2]);
  params.dont_fill_unlabeled = atoi(argv[3]);
  params.dont_cross_cmid     = atoi(argv[4]);
  params.dont_cross_fthresh  = atoi(argv[5]);
  params.use_multiple_seeds  = atoi(argv[6]);
  params.action              = atoi(argv[7]);
  params.argument            = atoi(argv[8]);

  if (params.use_multiple_seeds)
  {
    /* Since making a label will mess with the marked array, we need to
       copy the marked verts to use them as seeds. */
    for (n = 0; n < nmarked; n++)
      seeds[n] = marked[n];
    nseeds = nmarked;

    /* Clear the marked verts so any extraneous ones don't end up in the
       label. */
    clear_all_vertex_marks();

    /* Fill for each seed point. */
    first_fill = TRUE;
    for (n = 0; n < nseeds; n++)
    {
      fill_flood_from_seed (seeds[n], &params);

      /* If this was the first fill and this was a NEW_LABEL
         action, take the new label index that was pasesd back in
         the params block and make it the new argument of a
         ADD_LABEL action, so that the next label is added to the
         one we just made. */
      if (first_fill && FILL_ACTION_NEW_LABEL == params.action )
      {
        first_fill = FALSE;
        params.action = FILL_ACTION_ADD_LABEL;
        params.argument = params.new_label_index;
      }
    }
  }
  else
  {
    /* Use the selection as the seed. */
    seed = selection;

    /* Clear the marked verts so any extraneous ones don't end up in the
       label. */
    clear_all_vertex_marks();

    /* Fill using the first seed. */
    fill_flood_from_seed (seed, &params);
  }

  return TCL_OK;
}

int W_get_marked_vnos ( ClientData clientData, Tcl_Interp *interp,
                        int argc, char *argv[] )
{
  Tcl_Obj *list;
  int vno;
  VERTEX* v = NULL;

  if (argc != 1)
  {
    Tcl_SetResult(interp, "Wrong # args: get_marked_vnos", TCL_VOLATILE);
    return TCL_ERROR;
  }

  list = Tcl_NewListObj(0,NULL);

  for (vno = 0; vno < mris->nvertices; vno++)
  {
    v = &mris->vertices[vno];
    if (v->marked)
    {
      Tcl_ListObjAppendElement(interp,list,Tcl_NewIntObj(vno));
    }
  }

  Tcl_SetObjResult(interp,list);

  return TCL_OK;
}

int W_get_selected_path_vnos ( ClientData clientData, Tcl_Interp *interp,
                               int argc, char *argv[] )
{
  Tcl_Obj *list;
  int path_vno, vno;

  if (argc != 1)
  {
    Tcl_SetResult(interp, "Wrong # args: get_selected_path_vnos",
                  TCL_VOLATILE);
    return TCL_ERROR;
  }

  if (path_selected_path < 0 || path_selected_path >= path_num_paths)
  {
    Tcl_SetResult(interp, "No path selected.", TCL_VOLATILE);
    return TCL_ERROR;
  }

  list = Tcl_NewListObj(0,NULL);

  for (path_vno = 0;
       path_vno < path_paths[path_selected_path].num_vertices;
       path_vno++)
  {
    vno = path_paths[path_selected_path].vertices[path_vno];
    Tcl_ListObjAppendElement(interp,list,Tcl_NewIntObj(vno));
  }

  Tcl_SetObjResult(interp,list);

  return TCL_OK;
}

int W_save_tiff WBEGIN
ERR(2,"Wrong # args: save_tiff filename")
save_tiff(argv[1]);
WEND

int W_flip_normals WBEGIN
ERR(2,"Wrong # args: flip_normals axes")
flip_normals(argv[1]);
WEND

int W_mark_contiguous_vertices_over_thresh WBEGIN
ERR(1,"Wrong # args: mark_contiguous_vertices_over_thresh")
mark_contiguous_vertices_over_thresh();
WEND

int W_mark_contiguous_vertices_with_similar_curvature WBEGIN
ERR(1,"Wrong # args: mark_contiguous_vertices_with_similar_curvature")
mark_contiguous_vertices_with_similar_curvature();
WEND

int W_rip_all_vertices_except_contiguous_upripped WBEGIN
ERR(1,"Wrong # args: rip_all_vertices_except_contiguous_upripped")
rip_all_vertices_except_contiguous_upripped ();
WEND

int W_func_calc_correlation_and_write_to_overlay WBEGIN
ERR(2,"Wrong # args: func_calc_correlation_and_write_to_overlay field")
func_calc_correlation_and_write_to_overlay(selection,atoi(argv[1]));
WEND

int W_func_normalize WBEGIN
ERR(1,"Wrong # args: func_normalize")
func_normalize();
WEND

int W_cptn_set_format_string WBEGIN
ERR(2,"Wrong # args: cptn_set_format_string string")
cptn_set_format_string(argv[1]);
WEND



/* end rkt */
/*===================================================================*/

/* licensing */

#ifdef USE_LICENSE
#ifndef IRIX
extern char *crypt(char *, char *) ;
#endif
void checkLicense(char* dirname)
{
  FILE* lfile;
  char* email;
  char* magic;
  char* key;
  char* gkey;
  char* lfilename;

  lfilename = (char*)malloc(512);
  email = (char*)malloc(512);
  magic = (char*)malloc(512);
  key = (char*)malloc(512);
  gkey = (char*)malloc(1024);

  sprintf(lfilename,"%s/.license",dirname);

  lfile = fopen(lfilename,"r");
  if (lfile)
  {
    fscanf(lfile,"%s\n",email);
    fscanf(lfile,"%s\n",magic);
    fscanf(lfile,"%s\n",key);

    sprintf(gkey,"%s.%s",email,magic);
    if (strcmp(key,crypt(gkey,"*C*O*R*T*E*C*H*S*0*1*2*3*"))!=0)
    {
      printf("No valid license key !\n");
      exit(-1);
    }
  }
  else
  {
    printf("License file not found !\n");
    exit(-1);
  }
  free(email);
  free(magic);
  free(key);
  free(gkey);
  free(lfilename);
  return;
}

#endif

/* for tcl/tk */
#ifndef Windows_NT
static void StdinProc _ANSI_ARGS_((ClientData clientData, int mask));
#endif
static void Prompt _ANSI_ARGS_((Tcl_Interp *interp, int partial));
#ifndef TCL8
static Tk_Window mainWindow;
#endif
static Tcl_Interp *interp;
static Tcl_DString command;
static int tty;

int main(int argc, char *argv[])   /* new main */
{
  int code;
  int aliasflag=FALSE;
#ifndef TCL8
  static char *display = NULL;
#endif
  char tksurfer_tcl[NAME_LENGTH];
  char str[NAME_LENGTH];
  char alias_tcl[NAME_LENGTH];
  char *envptr;
  FILE *fp;
#ifndef USE_XGLUT_WINDOW
  struct timeval tv;
#endif /* USE_XGLUT_WINDOW */
  /* begin rkt */
  int found_script = FALSE;
  char* tksurfer_scripts_dir = NULL;
  int nargs;
  char tcl_cmd[STRLEN] = "";
  /* end rkt */

  /* rkt: check for and handle version tag */
  nargs =
    handle_version_option
    (argc, argv,
     "$Id: tksurfer.c,v 1.364 2016/12/11 14:33:47 fischl Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
    exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  //#ifdef USE_XGLUT_WINDOW
  //NJS note: glut needs to initialized, even if USE_XGLUT_WINDOW is not
  //defined, because if freeglut is used, then glutinit needs be to called
  //in order for certain tksurfer functions to work, like the 'Show Color
  //Scale Bar' button.
  /* init glut */
  DebugNote( ("Initializing glut") );
  glutInit( &argc, argv );
  glutInitDisplayMode( GLUT_RGBA | GLUT_DOUBLE );
  //#endif

  /* begin rkt */
  undo_initialize();
  vset_initialize();
  func_initialize();
  sclv_initialize();
  conv_initialize();
  labl_initialize();
  cncl_initialize();
  path_initialize();
  cptn_initialize();
  /* end rkt */

  /* get tksurfer tcl startup script location from environment */
  envptr = getenv("FREESURFER_HOME");
  if (envptr==NULL)
  {
    printf("tksurfer: env var FREESURFER_HOME undefined (use setenv)\n");
    printf("    [dir containing mri distribution]\n");
    exit(0);
  }

#ifdef USE_LICENSE
  checkLicense(envptr);
#endif

  /* begin rkt */

  /* here is the priority of the tksurfer.tcl files:
     1) tksurfer.new.tcl in TKSURFER_SCRIPTS_DIR
     2) tksurfer.tcl in TKSURFER_SCRIPTS_DIR
     3) tksurfer.new.tcl in local dir
     4) tksurfer.tcl in local dir
     5) tksurfer.new.tcl in FREESURFER_HOME/tktools
     6) tksurfer.tcl in FREESURFER_HOME/tktools
  */

  found_script = FALSE;
  tksurfer_scripts_dir = getenv ("TKSURFER_SCRIPTS_DIR");

  if (!found_script && tksurfer_scripts_dir)
  {
    sprintf (tksurfer_tcl, "%s/tksurfer.new.tcl", tksurfer_scripts_dir);
    if ((fp=fopen(tksurfer_tcl,"r"))!=NULL)
    {
      fclose(fp);
      found_script = TRUE;
    }
  }
  if (!found_script && tksurfer_scripts_dir)
  {
    sprintf (tksurfer_tcl, "%s/tksurfer.tcl", tksurfer_scripts_dir);
    if ((fp=fopen(tksurfer_tcl,"r"))!=NULL)
    {
      fclose(fp);
      found_script = TRUE;
    }
  }
  if (!found_script)
  {
    strcpy (tksurfer_tcl, "tksurfer.new.tcl");
    if ((fp=fopen(tksurfer_tcl,"r"))!=NULL)
    {
      fclose(fp);
      found_script = TRUE;
    }
  }
  if (!found_script)
  {
    strcpy (tksurfer_tcl, "tksurfer.tcl");
    if ((fp=fopen(tksurfer_tcl,"r"))!=NULL)
    {
      fclose(fp);
      found_script = TRUE;
    }
  }
  if (!found_script && envptr)
  {
    sprintf (tksurfer_tcl, "%s/tktools/tksurfer.new.tcl", envptr);
    if ((fp=fopen(tksurfer_tcl,"r"))!=NULL)
    {
      fclose(fp);
      found_script = TRUE;
    }
  }
  if (!found_script && envptr)
  {
    sprintf (tksurfer_tcl, "%s/tktools/tksurfer.tcl", envptr);
    if ((fp=fopen(tksurfer_tcl,"r"))!=NULL)
    {
      fclose(fp);
      found_script = TRUE;
    }
  }
  if (!found_script)
  {
    printf ("surfer: cannot find tksurfer.tcl script\n");
    exit (1);
  }

  /* end rkt */

  /* look for script: (1) cwd,
     (2) FREESURFER_HOME/lib/tcl, (3) [same]/alias.tcl */
  sprintf(script_tcl,"%s/lib/tcl/twocond-views.tcl",envptr);  /* default */
  if (argc==6)
  {
    strcpy(str,argv[4]);
    if (MATCH_STR("-tcl")) /* if command line script, run as batch job */
    {
      char *cp ;

      if (argc!=6)
      {
        printf("Usage: tksurfer [-]name hemi surf [-tcl script]\n");
        exit(0);
      }

      strcpy(script_tcl,argv[5]);
      fp = fopen(script_tcl,"r");  /* first, look in cwd */
      if (fp==NULL)                /* then FREESURFER_HOME/lib/tcl dir */
      {
        sprintf(script_tcl,"%s/lib/tcl/%s",envptr,argv[5]);
        fp = fopen(script_tcl,"r");
        if (fp == NULL)            /* then FREESURFER_HOME/tktools dir */
        {
          sprintf(script_tcl,"%s/tktools/%s",envptr,argv[5]);
          fp = fopen(script_tcl,"r");
          if (fp == NULL)            /* then TKSURFER_TCL_SCRIPTS dir */
          {
            cp = getenv("TKSURFER_TCL_SCRIPTS") ;
            if (cp)
              /* see if script is in users has own scripts directory */
            {
              sprintf(script_tcl,"%s/%s",cp,argv[5]);
              fp = fopen(script_tcl,"r");
            } else
              fp = NULL ;
            if (fp==NULL)              /* then aliases */
            {
              aliasflag = TRUE;
              sprintf(script_tcl,"%s/lib/tcl/alias.tcl",envptr);
              fp = fopen(script_tcl,"r");
              if (fp==NULL)            /* couldn't find it anywhere */
              {
                printf("tksurfer: (1) File ./%s not found\n",
                       argv[5]);
                printf("          (2) File %s/lib/tcl/%s "
                       "not found\n",
                       envptr,argv[5]);
                printf("          (3) File %s/lib/tcl/alias.tcl "
                       "not found\n",
                       envptr);
                exit(0);
              }
            }
          }
        }
      }
      scriptok = TRUE;
    } else
    {
      ; /* ignore 6 arg command lines without -tcl option */
    }
  }

  /* start surfer, now as function; gl window not opened yet */
  /*  printf("tksurfer: starting surfer\n");*/
  Surfer((ClientData) NULL, interp, argc, argv); /* event loop commented out */

  /* start tcl/tk; first make interpreter */
  interp = Tcl_CreateInterp();
  /* begin rkt */
  g_interp = interp;

  /* end rkt */

  /* make main window (not displayed until event loop starts) */
  /*mainWindow = TkCreateMainWindow(interp, display, argv[0], "Tk");
    if (mainWindow == NULL) {
    fprintf(stderr, "%s\n", interp->result);
    exit(1); }
  */
  /* set the "tcl_interactive" variable */
  tty = isatty(0);
  Tcl_SetVar(interp, "tcl_interactive", (tty) ? "1" : "0", TCL_GLOBAL_ONLY);
  if (tty) promptflag = TRUE;  /* no-CR blocks pipe log read */

  /* read tcl/tk internal startup scripts */
  if (Tcl_Init(interp) == TCL_ERROR)
  {
    fprintf(stderr, "Tcl_Init failed: %s\n", interp->result);
  }
  if (Tk_Init(interp) == TCL_ERROR)
  {
    fprintf(stderr, "Tk_Init failed: %s\n", interp->result);
  }

  // later version of Tix needs these.
  // Unfortunately Tix does not define Major and Minor
  // to distinguish what it is (stupidity).  I had to use gnu ;-(....
  // Do the following only for RedHat Enterprise Linux only

#if NEEDS_ITCL_ITK
  if (Itcl_Init(interp) == TCL_ERROR)
  {
    fprintf(stderr, "Itcl_Init failed: %s\n", interp->result);
  }
  if (Itk_Init(interp) == TCL_ERROR)
  {
    fprintf(stderr, "Itk_Init failed: %s\n", interp->result);
  }
#endif

  /* begin rkt */
  if (Tix_Init(interp) == TCL_ERROR)
  {
    fprintf(stderr, "Tix_Init failed: %s\n", interp->result);
  }
  if (Blt_Init(interp) == TCL_ERROR)
  {
    fprintf(stderr, "Blt_Init failed: %s\n", interp->result);
  }

  /* Initialize our Fsgdf functions. This is in fsgdf_wrap.c */
  if (Fsgdf_Init(interp) == TCL_ERROR)
  {
    fprintf(stderr, "Fsgdf_Init failed: %s\n", interp->result);
  }

  Tcl_StaticPackage( interp, "BLT", Blt_Init, Blt_SafeInit );
#ifdef Windows_NT
#define Tix_SafeInit NULL
#endif
  Tcl_StaticPackage( interp, "Tix", Tix_Init, Tix_SafeInit );
  /* end rkt */

  /*=======================================================================*/
  /* register wrapped surfer functions with interpreter */
  Tcl_CreateCommand(interp, "swap_buffers",
                    (Tcl_CmdProc*) W_swap_buffers,       REND);
  Tcl_CreateCommand(interp, "to_single_buffer",
                    (Tcl_CmdProc*) W_to_single_buffer,   REND);
  Tcl_CreateCommand(interp, "to_double_buffer",
                    (Tcl_CmdProc*) W_to_double_buffer,   REND);
  Tcl_CreateCommand(interp, "open_window",
                    (Tcl_CmdProc*) W_open_window,        REND);
  Tcl_CreateCommand(interp, "help",
                    (Tcl_CmdProc*) W_help,               REND);
  Tcl_CreateCommand(interp, "redraw",
                    (Tcl_CmdProc*) W_redraw,             REND);
  Tcl_CreateCommand(interp, "redraw_second",
                    (Tcl_CmdProc*) W_redraw_second,      REND);
  Tcl_CreateCommand(interp, "shrink",
                    (Tcl_CmdProc*) W_shrink,             REND);
  Tcl_CreateCommand(interp, "area_shrink",
                    (Tcl_CmdProc*) W_area_shrink,        REND);
  Tcl_CreateCommand(interp, "sphere_shrink",
                    (Tcl_CmdProc*) W_sphere_shrink,      REND);
  Tcl_CreateCommand(interp, "ellipsoid_project",
                    (Tcl_CmdProc*) W_ellipsoid_project,  REND);
  Tcl_CreateCommand(interp, "ellipsoid_morph",
                    (Tcl_CmdProc*) W_ellipsoid_morph,     REND);
  Tcl_CreateCommand(interp, "ellipsoid_shrink",
                    (Tcl_CmdProc*) W_ellipsoid_shrink,     REND);
  Tcl_CreateCommand(interp, "ellipsoid_shrink_bug",
                    (Tcl_CmdProc*) W_ellipsoid_shrink_bug,REND);
  Tcl_CreateCommand(interp, "curv_shrink_to_fill",
                    (Tcl_CmdProc*) W_curv_shrink_to_fill,REND);
  Tcl_CreateCommand(interp, "smooth_curvim",
                    (Tcl_CmdProc*) W_smooth_curvim,      REND);
  Tcl_CreateCommand(interp, "smooth_curv",
                    (Tcl_CmdProc*) W_smooth_curv,        REND);
  Tcl_CreateCommand(interp, "smooth_val",
                    (Tcl_CmdProc*) W_smooth_val,         REND);
  Tcl_CreateCommand(interp, "smooth_val_sparse",
                    (Tcl_CmdProc*) W_smooth_val_sparse,  REND);
  Tcl_CreateCommand(interp, "smooth_curvim_sparse",
                    (Tcl_CmdProc*) W_smooth_curvim_sparse,    REND);
  Tcl_CreateCommand(interp, "smooth_fs",
                    (Tcl_CmdProc*) W_smooth_fs,          REND);
  Tcl_CreateCommand(interp, "add_subject_to_average_curvim",
                    (Tcl_CmdProc*) W_add_subject_to_average_curvim,REND);
  Tcl_CreateCommand(interp, "read_curv_images",
                    (Tcl_CmdProc*) W_read_curv_images,   REND);
  Tcl_CreateCommand(interp, "read_stds",
                    (Tcl_CmdProc*) W_read_stds,          REND);
  Tcl_CreateCommand(interp, "read_second_binary_surf",
                    (Tcl_CmdProc*) W_read_second_binary_surf,          REND);
  Tcl_CreateCommand(interp, "read_second_binary_curv",
                    (Tcl_CmdProc*) W_read_second_binary_curv,          REND);
  Tcl_CreateCommand(interp, "normalize_second_binary_curv",
                    (Tcl_CmdProc*) W_normalize_second_binary_curv,     REND);
  Tcl_CreateCommand(interp, "curv_to_curvim",
                    (Tcl_CmdProc*) W_curv_to_curvim,                   REND);
  Tcl_CreateCommand(interp, "second_surface_curv_to_curvim",
                    (Tcl_CmdProc*) W_second_surface_curv_to_curvim,    REND);
  Tcl_CreateCommand(interp, "curvim_to_second_surface",
                    (Tcl_CmdProc*) W_curvim_to_second_surface,         REND);
  Tcl_CreateCommand(interp, "swap_curv",
                    (Tcl_CmdProc*) W_swap_curv,                        REND);
  Tcl_CreateCommand(interp, "curvim_to_surface",
                    (Tcl_CmdProc*) W_curvim_to_surface,                REND);
  Tcl_CreateCommand(interp, "read_binary_surf",
                    (Tcl_CmdProc*) W_read_binary_surf,   REND);
  Tcl_CreateCommand(interp, "read_surf",
                    (Tcl_CmdProc*) W_read_surf,   REND);
  Tcl_CreateCommand(interp, "save_surf",
                    (Tcl_CmdProc*) W_save_surf,   REND);
  Tcl_CreateCommand(interp, "store_surf",
                    (Tcl_CmdProc*) W_save_surf,   REND);
  Tcl_CreateCommand(interp, "restore_surf",
                    (Tcl_CmdProc*) W_restore_surf,   REND);
  Tcl_CreateCommand(interp, "surf",
                    (Tcl_CmdProc*) W_show_surf,   REND);
  Tcl_CreateCommand(interp, "read_binary_curv",
                    (Tcl_CmdProc*) W_read_binary_curv,   REND);
  Tcl_CreateCommand(interp, "read_binary_sulc",
                    (Tcl_CmdProc*) W_read_binary_sulc,   REND);
  Tcl_CreateCommand(interp, "read_binary_values",
                    (Tcl_CmdProc*) W_read_binary_values, REND);
  Tcl_CreateCommand(interp, "read_binary_values_frame",
                    (Tcl_CmdProc*) W_read_binary_values_frame,        REND);
  Tcl_CreateCommand(interp, "read_annotated_image",
                    (Tcl_CmdProc*) W_read_annotated_image,REND);
  Tcl_CreateCommand(interp, "read_annotations",
                    (Tcl_CmdProc*) W_read_annotations,REND);
  Tcl_CreateCommand(interp, "read_binary_patch",
                    (Tcl_CmdProc*) W_read_binary_patch,  REND);
  Tcl_CreateCommand(interp, "read_fieldsign",
                    (Tcl_CmdProc*) W_read_fieldsign,     REND);
  Tcl_CreateCommand(interp, "read_fsmask",
                    (Tcl_CmdProc*) W_read_fsmask,        REND);
  Tcl_CreateCommand(interp, "write_binary_areas",
                    (Tcl_CmdProc*) W_write_binary_areas, REND);
  Tcl_CreateCommand(interp, "write_binary_surface",
                    (Tcl_CmdProc*) W_write_binary_surface,REND);
  Tcl_CreateCommand(interp, "write_binary_curv",
                    (Tcl_CmdProc*) W_write_binary_curv,REND);
  Tcl_CreateCommand(interp, "write_binary_sulc",
                    (Tcl_CmdProc*) W_write_binary_sulc,REND);
  Tcl_CreateCommand(interp, "write_binary_values",
                    (Tcl_CmdProc*) W_write_binary_values,REND);
  Tcl_CreateCommand(interp, "write_binary_patch",
                    (Tcl_CmdProc*) W_write_binary_patch, REND);
  Tcl_CreateCommand(interp, "write_labeled_vertices",
                    (Tcl_CmdProc*) W_write_labeled_vertices,    REND);
  Tcl_CreateCommand(interp, "read_labeled_vertices",
                    (Tcl_CmdProc*) W_read_labeled_vertices,     REND);
  Tcl_CreateCommand(interp, "read_and_color_labeled_vertices",
                    (Tcl_CmdProc*) W_read_and_color_labeled_vertices, REND);
  Tcl_CreateCommand(interp, "write_fieldsign",
                    (Tcl_CmdProc*) W_write_fieldsign,    REND);
  Tcl_CreateCommand(interp, "write_fsmask",
                    (Tcl_CmdProc*) W_write_fsmask,       REND);
  Tcl_CreateCommand(interp, "write_vrml",
                    (Tcl_CmdProc*) W_write_vrml,         REND);
  Tcl_CreateCommand(interp, "write_binary_dipoles",
                    (Tcl_CmdProc*) W_write_binary_dipoles,REND);
  Tcl_CreateCommand(interp, "write_binary_decimation",
                    (Tcl_CmdProc*) W_write_binary_decimation,     REND);
  Tcl_CreateCommand(interp, "write_dipoles",
                    (Tcl_CmdProc*) W_write_dipoles,      REND);
  Tcl_CreateCommand(interp, "write_decimation",
                    (Tcl_CmdProc*) W_write_decimation,   REND);
  Tcl_CreateCommand(interp, "write_curv_images",
                    (Tcl_CmdProc*) W_write_curv_images,  REND);
  Tcl_CreateCommand(interp, "write_fill_images",
                    (Tcl_CmdProc*) W_write_fill_images,  REND);
  Tcl_CreateCommand(interp, "fill_second_surface",
                    (Tcl_CmdProc*) W_fill_second_surface,REND);
  Tcl_CreateCommand(interp, "subsample_dist",
                    (Tcl_CmdProc*) W_subsample_dist,     REND);
  Tcl_CreateCommand(interp, "subsample_orient",
                    (Tcl_CmdProc*) W_subsample_orient,   REND);
  Tcl_CreateCommand(interp, "write_subsample",
                    (Tcl_CmdProc*) W_write_subsample,    REND);
  Tcl_CreateCommand(interp, "compute_curvature",
                    (Tcl_CmdProc*) W_compute_curvature,  REND);
  Tcl_CreateCommand(interp, "compute_CMF",
                    (Tcl_CmdProc*) W_compute_CMF,        REND);
  Tcl_CreateCommand(interp, "compute_cortical_thickness",
                    (Tcl_CmdProc*) W_compute_cortical_thickness,     REND);
  Tcl_CreateCommand(interp, "clear_curvature",
                    (Tcl_CmdProc*) W_clear_curvature,    REND);
  Tcl_CreateCommand(interp, "clear_vals",
                    (Tcl_CmdProc*) W_clear_vals,    REND);
  Tcl_CreateCommand(interp, "clear_ripflags",
                    (Tcl_CmdProc*) W_clear_ripflags,     REND);
  Tcl_CreateCommand(interp, "restore_ripflags",
                    (Tcl_CmdProc*) W_restore_ripflags,   REND);
  Tcl_CreateCommand(interp, "floodfill_marked_patch",
                    (Tcl_CmdProc*) W_floodfill_marked_patch,    REND);
  Tcl_CreateCommand(interp, "rip_unmarked_vertices",
                    (Tcl_CmdProc*) W_rip_unmarked_vertices,    REND);
  Tcl_CreateCommand(interp, "dilate_ripped",
                    (Tcl_CmdProc*) W_dilate_ripped,    REND);
  Tcl_CreateCommand(interp, "twocond",
                    (Tcl_CmdProc*) W_twocond,                    REND);
  Tcl_CreateCommand(interp, "cut_line",
                    (Tcl_CmdProc*) W_cut_line,           REND);
  Tcl_CreateCommand(interp, "plot_curv",
                    (Tcl_CmdProc*) W_plot_curv,          REND);
  Tcl_CreateCommand(interp, "draw_fundus",
                    (Tcl_CmdProc*) W_draw_fundus,        REND);
  Tcl_CreateCommand(interp, "plot_marked",
                    (Tcl_CmdProc*) W_plot_marked,        REND);
  Tcl_CreateCommand(interp, "put_retinotopy_stats_in_vals",
                    (Tcl_CmdProc*) W_put_retinotopy_stats_in_vals,    REND);
  Tcl_CreateCommand(interp, "draw_vector",
                    (Tcl_CmdProc*) W_draw_vector,        REND);
  Tcl_CreateCommand(interp, "cut_plane",
                    (Tcl_CmdProc*) W_cut_plane,          REND);
  Tcl_CreateCommand(interp, "flatten",
                    (Tcl_CmdProc*) W_flatten,            REND);
  Tcl_CreateCommand(interp, "normalize_binary_curv",
                    (Tcl_CmdProc*) W_normalize_binary_curv,     REND);
  Tcl_CreateCommand(interp, "normalize_area",
                    (Tcl_CmdProc*) W_normalize_area,     REND);
  Tcl_CreateCommand(interp, "normalize_curvature",
                    (Tcl_CmdProc*) W_normalize_curvature,REND);
  Tcl_CreateCommand(interp, "shift_values",
                    (Tcl_CmdProc*) W_shift_values,       REND);
  Tcl_CreateCommand(interp, "swap_values",
                    (Tcl_CmdProc*) W_swap_values,        REND);
  Tcl_CreateCommand(interp, "swap_stat_val",
                    (Tcl_CmdProc*) W_swap_stat_val,      REND);
  Tcl_CreateCommand(interp, "swap_val_val2",
                    (Tcl_CmdProc*) W_swap_val_val2,      REND);
  Tcl_CreateCommand(interp, "compute_angles",
                    (Tcl_CmdProc*) W_compute_angles,     REND);
  Tcl_CreateCommand(interp, "compute_fieldsign",
                    (Tcl_CmdProc*) W_compute_fieldsign,  REND);
  Tcl_CreateCommand(interp, "draw_radius",
                    (Tcl_CmdProc*) W_draw_radius,        REND);
  Tcl_CreateCommand(interp, "draw_theta",
                    (Tcl_CmdProc*) W_draw_theta,         REND);
  Tcl_CreateCommand(interp, "save_rgb",
                    (Tcl_CmdProc*) W_save_rgb,           REND);
  Tcl_CreateCommand(interp, "save_rgb_named_orig",
                    (Tcl_CmdProc*) W_save_rgb_named_orig,REND);
  Tcl_CreateCommand(interp, "save_rgb_cmp_frame",
                    (Tcl_CmdProc*) W_save_rgb_cmp_frame, REND);
  Tcl_CreateCommand(interp, "open_rgb_cmp_named",
                    (Tcl_CmdProc*) W_open_rgb_cmp_named, REND);
  Tcl_CreateCommand(interp, "save_rgb_cmp_frame_named",
                    (Tcl_CmdProc*) W_save_rgb_cmp_frame_named,     REND);
  Tcl_CreateCommand(interp, "close_rgb_cmp_named",
                    (Tcl_CmdProc*) W_close_rgb_cmp_named,REND);
  Tcl_CreateCommand(interp, "rotate_brain_x",
                    (Tcl_CmdProc*) W_rotate_brain_x,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_y",
                    (Tcl_CmdProc*) W_rotate_brain_y,     REND);
  Tcl_CreateCommand(interp, "rotate_brain_z",
                    (Tcl_CmdProc*) W_rotate_brain_z,     REND);
  Tcl_CreateCommand(interp, "translate_brain_x",
                    (Tcl_CmdProc*) W_translate_brain_x,  REND);
  Tcl_CreateCommand(interp, "translate_brain_y",
                    (Tcl_CmdProc*) W_translate_brain_y,  REND);
  Tcl_CreateCommand(interp, "translate_brain_z",
                    (Tcl_CmdProc*) W_translate_brain_z,  REND);
  Tcl_CreateCommand(interp, "scale_brain",
                    (Tcl_CmdProc*) W_scale_brain,        REND);
  Tcl_CreateCommand(interp, "resize_window",
                    (Tcl_CmdProc*) W_resize_window,      REND);
  Tcl_CreateCommand(interp, "setsize_window",
                    (Tcl_CmdProc*) W_setsize_window,     REND);
  Tcl_CreateCommand(interp, "move_window",
                    (Tcl_CmdProc*) W_move_window,        REND);
  Tcl_CreateCommand(interp, "do_lighting_model",
                    (Tcl_CmdProc*) W_do_lighting_model,  REND);
  Tcl_CreateCommand(interp, "restore_zero_position",
                    (Tcl_CmdProc*) W_restore_zero_position,  REND);
  Tcl_CreateCommand(interp, "restore_initial_position",
                    (Tcl_CmdProc*) W_restore_initial_position,   REND);
  Tcl_CreateCommand(interp, "make_lateral_view",
                    (Tcl_CmdProc*) W_make_lateral_view,  REND);
  Tcl_CreateCommand(interp, "make_lateral_view_second",
                    (Tcl_CmdProc*) W_make_lateral_view_second,      REND);
  Tcl_CreateCommand(interp, "read_view_matrix",
                    (Tcl_CmdProc*) W_read_view_matrix,   REND);
  Tcl_CreateCommand(interp, "write_view_matrix",
                    (Tcl_CmdProc*) W_write_view_matrix,  REND);
  Tcl_CreateCommand(interp, "read_really_matrix",
                    (Tcl_CmdProc*) W_read_really_matrix, REND);
  Tcl_CreateCommand(interp, "write_really_matrix",
                    (Tcl_CmdProc*) W_write_really_matrix,REND);
  Tcl_CreateCommand(interp, "really_translate_brain",
                    (Tcl_CmdProc*) W_really_translate_brain,      REND);
  Tcl_CreateCommand(interp, "really_scale_brain",
                    (Tcl_CmdProc*) W_really_scale_brain, REND);
  Tcl_CreateCommand(interp, "align_sphere",
                    (Tcl_CmdProc*) W_align_sphere,         REND);
  Tcl_CreateCommand(interp, "really_rotate_brain_x",
                    (Tcl_CmdProc*) W_really_rotate_brain_x,    REND);
  Tcl_CreateCommand(interp, "really_rotate_brain_y",
                    (Tcl_CmdProc*) W_really_rotate_brain_y,     REND);
  Tcl_CreateCommand(interp, "really_rotate_brain_z",
                    (Tcl_CmdProc*) W_really_rotate_brain_z,    REND);
  Tcl_CreateCommand(interp, "really_center_brain",
                    (Tcl_CmdProc*) W_really_center_brain, REND);
  Tcl_CreateCommand(interp, "really_center_second_brain",
                    (Tcl_CmdProc*) W_really_center_second_brain, REND);
  Tcl_CreateCommand(interp, "really_align_brain",
                    (Tcl_CmdProc*) W_really_align_brain, REND);
  Tcl_CreateCommand(interp, "read_binary_decimation",
                    (Tcl_CmdProc*) W_read_binary_decimation,    REND);
  Tcl_CreateCommand(interp, "read_binary_dipoles",
                    (Tcl_CmdProc*) W_read_binary_dipoles,     REND);
  Tcl_CreateCommand(interp, "load_vals_from_sol",
                    (Tcl_CmdProc*) W_load_vals_from_sol,   REND);
  Tcl_CreateCommand(interp, "load_var_from_sol",
                    (Tcl_CmdProc*) W_load_var_from_sol,     REND);
  Tcl_CreateCommand(interp, "compute_timecourses",
                    (Tcl_CmdProc*) W_compute_timecourses, REND);
  Tcl_CreateCommand(interp, "filter_recs",
                    (Tcl_CmdProc*) W_filter_recs,                REND);
  Tcl_CreateCommand(interp, "read_rec",
                    (Tcl_CmdProc*) W_read_rec,                      REND);
  Tcl_CreateCommand(interp, "read_iop",
                    (Tcl_CmdProc*) W_read_iop,                      REND);
  Tcl_CreateCommand(interp, "read_ncov",
                    (Tcl_CmdProc*) W_read_ncov,                    REND);
  Tcl_CreateCommand(interp, "normalize_time_courses",
                    (Tcl_CmdProc*) W_normalize_time_courses,   REND);
  Tcl_CreateCommand(interp, "compute_timecourses",
                    (Tcl_CmdProc*) W_compute_timecourses,  REND);
  Tcl_CreateCommand(interp, "compute_pval_fwd",
                    (Tcl_CmdProc*) W_compute_pval_fwd,      REND);
  Tcl_CreateCommand(interp, "compute_select_fwd",
                    (Tcl_CmdProc*) W_compute_select_fwd,  REND);
  Tcl_CreateCommand(interp, "compute_pval_inv",
                    (Tcl_CmdProc*) W_compute_pval_inv,      REND);
  Tcl_CreateCommand(interp, "normalize_inverse",
                    (Tcl_CmdProc*) W_normalize_inverse,    REND);
  Tcl_CreateCommand(interp, "find_orig_vertex_coordinates",
                    (Tcl_CmdProc*) W_find_orig_vertex_coordinates, REND);
  Tcl_CreateCommand(interp, "select_orig_vertex_coordinates",
                    (Tcl_CmdProc*) W_select_orig_vertex_coordinates, REND);
  Tcl_CreateCommand(interp, "select_talairach_point",
                    (Tcl_CmdProc*) W_select_talairach_point, REND);
  Tcl_CreateCommand(interp, "print_nearest_vertex_to_talairach_point",
                    (Tcl_CmdProc*) W_print_nearest_vertex_to_talairach_point, REND);
  Tcl_CreateCommand(interp, "read_white_vertex_coordinates",
                    (Tcl_CmdProc*) W_read_white_vertex_coordinates, REND);
  Tcl_CreateCommand(interp, "read_pial_vertex_coordinates",
                    (Tcl_CmdProc*) W_read_pial_vertex_coordinates, REND);
  Tcl_CreateCommand(interp, "read_orig_vertex_coordinates",
                    (Tcl_CmdProc*) W_read_orig_vertex_coordinates, REND);
  Tcl_CreateCommand(interp, "read_canon_vertex_coordinates",
                    (Tcl_CmdProc*) W_read_canon_vertex_coordinates, REND);
  Tcl_CreateCommand(interp, "send_spherical_point",
                    (Tcl_CmdProc*) W_send_spherical_point, REND);
  Tcl_CreateCommand(interp, "send_contralateral_point",
                    (Tcl_CmdProc*) W_send_contralateral_point, REND);
  Tcl_CreateCommand(interp, "send_to_subject",
                    (Tcl_CmdProc*) W_send_to_subject, REND);
  Tcl_CreateCommand(interp, "send_to_other_hemi",
                    (Tcl_CmdProc*) W_send_to_other_hemi, REND);
  Tcl_CreateCommand(interp, "resend_to_subject",
                    (Tcl_CmdProc*) W_resend_to_subject, REND);
  Tcl_CreateCommand(interp, "drawcb",
                    (Tcl_CmdProc*) W_drawcb, REND);
  Tcl_CreateCommand(interp, "read_ellipsoid_vertex_coordinates",
                    (Tcl_CmdProc*) W_read_ellipsoid_vertex_coordinates, REND);
  Tcl_CreateCommand(interp, "invert_surface",
                    (Tcl_CmdProc*) W_invert_surface, REND);
  Tcl_CreateCommand(interp, "fix_nonzero_vals",
                    (Tcl_CmdProc*) W_fix_nonzero_vals, REND);
  Tcl_CreateCommand(interp, "invert_vertex",
                    (Tcl_CmdProc*) W_invert_vertex, REND);
  Tcl_CreateCommand(interp, "invert_face",
                    (Tcl_CmdProc*) W_invert_face, REND);
  Tcl_CreateCommand(interp, "mark_annotation",
                    (Tcl_CmdProc*) W_mark_annotation, REND);
  Tcl_CreateCommand(interp, "mark_faces",
                    (Tcl_CmdProc*) W_mark_faces, REND);
  Tcl_CreateCommand(interp, "mark_face",
                    (Tcl_CmdProc*) W_mark_face, REND);
  Tcl_CreateCommand(interp, "dump_vertex",
                    (Tcl_CmdProc*) W_dump_vertex, REND);
  Tcl_CreateCommand(interp, "val_to_mark",
                    (Tcl_CmdProc*) W_val_to_mark, REND);
  Tcl_CreateCommand(interp, "set_area_thresh",
                    (Tcl_CmdProc*) W_set_area_thresh, REND);
  Tcl_CreateCommand(interp, "resize_brain",
                    (Tcl_CmdProc*) W_resize_brain, REND);
  Tcl_CreateCommand(interp, "transform_brain",
                    (Tcl_CmdProc*) W_transform_brain, REND);
  Tcl_CreateCommand(interp, "show_flat_regions",
                    (Tcl_CmdProc*) W_show_flat_regions, REND);
  Tcl_CreateCommand(interp, "val_to_stat",
                    (Tcl_CmdProc*) W_val_to_stat, REND);
  Tcl_CreateCommand(interp, "set_vals",
                    (Tcl_CmdProc*) W_set_vals, REND);
  Tcl_CreateCommand(interp, "stat_to_val",
                    (Tcl_CmdProc*) W_stat_to_val, REND);
  Tcl_CreateCommand(interp, "scale_vals",
                    (Tcl_CmdProc*) W_scale_vals, REND);
  Tcl_CreateCommand(interp, "read_soltimecourse",
                    (Tcl_CmdProc*) W_read_soltimecourse, REND);
  Tcl_CreateCommand(interp, "read_imag_vals",
                    (Tcl_CmdProc*) W_read_imag_vals, REND);
  Tcl_CreateCommand(interp, "sol_plot",
                    (Tcl_CmdProc*) W_sol_plot, REND);
  Tcl_CreateCommand(interp, "remove_triangle_links",
                    (Tcl_CmdProc*) W_remove_triangle_links, REND);
  Tcl_CreateCommand(interp, "f_to_t",
                    (Tcl_CmdProc*) W_f_to_t, REND);
  Tcl_CreateCommand(interp, "label_to_stat",
                    (Tcl_CmdProc*) W_label_to_stat, REND);
  Tcl_CreateCommand(interp, "taubin_smooth",
                    (Tcl_CmdProc*) W_taubin_smooth, REND);

  Tcl_CreateCommand(interp, "label_from_stats",
                    (Tcl_CmdProc*) W_label_from_stats, REND);
  Tcl_CreateCommand(interp, "label_set_stats",
                    (Tcl_CmdProc*) W_label_set_stats, REND);

  Tcl_CreateCommand(interp, "t_to_p",
                    (Tcl_CmdProc*) W_t_to_p, REND);
  Tcl_CreateCommand(interp, "f_to_p",
                    (Tcl_CmdProc*) W_f_to_p, REND);
  Tcl_CreateCommand(interp, "val_to_curv",
                    (Tcl_CmdProc*) W_val_to_curv, REND);
  Tcl_CreateCommand(interp, "curv_to_val",
                    (Tcl_CmdProc*) W_curv_to_val, REND);
  Tcl_CreateCommand(interp, "read_curv_to_val",
                    (Tcl_CmdProc*) W_read_curv_to_val, REND);
  Tcl_CreateCommand(interp, "read_and_smooth_parcellation",
                    (Tcl_CmdProc*) W_read_and_smooth_parcellation, REND);
  Tcl_CreateCommand(interp, "read_parcellation",
                    (Tcl_CmdProc*) W_read_parcellation, REND);
  Tcl_CreateCommand(interp, "deconvolve_weights",
                    (Tcl_CmdProc*) W_deconvolve_weights, REND);
  Tcl_CreateCommand(interp, "read_disc",
                    (Tcl_CmdProc*) W_read_disc, REND);
  Tcl_CreateCommand(interp, "mask_label",
                    (Tcl_CmdProc*) W_mask_label, REND);
  Tcl_CreateCommand(interp, "orient_sphere",
                    (Tcl_CmdProc*) W_orient_sphere, REND);
  Tcl_CreateCommand(interp, "dump_faces",
                    (Tcl_CmdProc*) W_dump_faces, REND);
  Tcl_CreateCommand(interp, "load_gcsa",
                    (Tcl_CmdProc*) W_load_gcsa, REND);
  Tcl_CreateCommand(interp, "draw_ellipsoid_latlong",
                    (Tcl_CmdProc*) W_draw_ellipsoid_latlong, REND);
  Tcl_CreateCommand(interp, "left_click",
                    (Tcl_CmdProc*) W_left_click, REND);
  Tcl_CreateCommand(interp, "plot_all_time_courses",
                    (Tcl_CmdProc*) W_plot_all_time_courses, REND);
  Tcl_CreateCommand(interp, "read_plot_list",
                    (Tcl_CmdProc*) W_read_plot_list, REND);
  Tcl_CreateCommand(interp, "read_vertex_list",
                    (Tcl_CmdProc*) W_read_vertex_list, REND);
  Tcl_CreateCommand(interp, "draw_cursor",
                    (Tcl_CmdProc*) W_draw_cursor, REND);
  Tcl_CreateCommand(interp, "draw_marked_vertices",
                    (Tcl_CmdProc*) W_draw_marked_vertices, REND);
  Tcl_CreateCommand(interp, "mark_vertex",
                    (Tcl_CmdProc*) W_mark_vertex, REND);
  Tcl_CreateCommand(interp, "mark_translated_vertex",
                    (Tcl_CmdProc*) W_mark_translated_vertex, REND);
  Tcl_CreateCommand(interp, "draw_all_cursor",
                    (Tcl_CmdProc*) W_draw_all_cursor, REND);
  Tcl_CreateCommand(interp, "draw_all_vertex_cursor",
                    (Tcl_CmdProc*) W_draw_all_vertex_cursor, REND);
  Tcl_CreateCommand(interp, "clear_all_vertex_cursor",
                    (Tcl_CmdProc*) W_clear_all_vertex_cursor, REND);
  /* begin rkt */
  Tcl_CreateCommand(interp, "send_current_labels",
                    (Tcl_CmdProc*) W_send_current_labels, REND);

  Tcl_CreateCommand(interp, "select_vertex_by_vno",
                    (Tcl_CmdProc*) W_select_vertex_by_vno, REND);

  Tcl_CreateCommand(interp, "swap_vertex_fields",
                    (Tcl_CmdProc*) W_swap_vertex_fields, REND);

  Tcl_CreateCommand(interp, "clear_vertex_marks",
                    (Tcl_CmdProc*) W_clear_vertex_marks, REND);
  Tcl_CreateCommand(interp, "clear_all_vertex_marks",
                    (Tcl_CmdProc*) W_clear_all_vertex_marks, REND);

  Tcl_CreateCommand(interp, "close_marked_vertices",
                    (Tcl_CmdProc*) W_close_marked_vertices, REND);

  Tcl_CreateCommand(interp, "undo_last_action",
                    (Tcl_CmdProc*) W_undo_last_action, REND);

  Tcl_CreateCommand(interp, "sclv_read_from_dotw",
                    (Tcl_CmdProc*) W_sclv_read_from_dotw, REND);
  Tcl_CreateCommand(interp, "sclv_read_binary_values",
                    (Tcl_CmdProc*) W_sclv_read_from_dotw, REND);
  Tcl_CreateCommand(interp, "sclv_read_from_dotw_frame",
                    (Tcl_CmdProc*) W_sclv_read_from_dotw_frame, REND);
  Tcl_CreateCommand(interp, "sclv_read_binary_values_frame",
                    (Tcl_CmdProc*) W_sclv_read_from_dotw_frame, REND);
  Tcl_CreateCommand(interp, "sclv_read_from_volume",
                    (Tcl_CmdProc*) W_sclv_read_from_volume, REND);
  Tcl_CreateCommand(interp, "sclv_read_bfile_values",
                    (Tcl_CmdProc*) W_sclv_read_from_volume, REND);
  Tcl_CreateCommand(interp, "sclv_write_dotw",
                    (Tcl_CmdProc*) W_sclv_write_dotw, REND);
  Tcl_CreateCommand(interp, "sclv_load_label_value_file",
                    (Tcl_CmdProc*) W_sclv_load_label_value_file, REND);
  Tcl_CreateCommand(interp, "sclv_smooth",
                    (Tcl_CmdProc*) W_sclv_smooth, REND);
  Tcl_CreateCommand(interp, "sclv_set_overlay_alpha",
                    (Tcl_CmdProc*) W_sclv_set_overlay_alpha, REND);
  Tcl_CreateCommand(interp, "sclv_set_current_field",
                    (Tcl_CmdProc*) W_sclv_set_current_field, REND);
  Tcl_CreateCommand(interp, "sclv_unload_field",
                    (Tcl_CmdProc*) W_sclv_unload_field, REND);
  Tcl_CreateCommand(interp, "sclv_set_current_timepoint",
                    (Tcl_CmdProc*) W_sclv_set_current_timepoint, REND);
  Tcl_CreateCommand(interp, "sclv_copy_view_settings_from_current_field",
                    (Tcl_CmdProc*) W_sclv_copy_view_settings_from_current_field, REND);
  Tcl_CreateCommand(interp, "sclv_copy_all_view_settings_from_current_field",
                    (Tcl_CmdProc*) W_sclv_copy_all_view_settings_from_current_field, REND);
  Tcl_CreateCommand(interp, "sclv_copy_view_settings_from_field",
                    (Tcl_CmdProc*) W_sclv_copy_view_settings_from_field, REND);
  Tcl_CreateCommand(interp, "sclv_set_current_threshold_from_percentile",
                    (Tcl_CmdProc*) W_sclv_set_current_threshold_from_percentile, REND);
  Tcl_CreateCommand(interp, "sclv_set_current_threshold_using_fdr",
                    (Tcl_CmdProc*) W_sclv_set_current_threshold_using_fdr, REND);
  Tcl_CreateCommand(interp, "sclv_send_histogram",
                    (Tcl_CmdProc*) W_sclv_send_histogram, REND);
  Tcl_CreateCommand(interp, "sclv_send_current_field_info",
                    (Tcl_CmdProc*) W_sclv_send_current_field_info, REND);
  Tcl_CreateCommand(interp, "sclv_get_normalized_color_for_value",
                    (Tcl_CmdProc*) W_sclv_get_normalized_color_for_value, REND);

  Tcl_CreateCommand(interp, "read_surface_vertex_set",
                    (Tcl_CmdProc*) W_read_surface_vertex_set, REND);
  Tcl_CreateCommand(interp, "set_current_vertex_set",
                    (Tcl_CmdProc*) W_set_current_vertex_set, REND);

  Tcl_CreateCommand(interp, "func_load_timecourse",
                    (Tcl_CmdProc*) W_func_load_timecourse, REND);
  Tcl_CreateCommand(interp, "func_load_timecourse_offset",
                    (Tcl_CmdProc*) W_func_load_timecourse_offset, REND);
  Tcl_CreateCommand(interp, "func_select_selected_vertex",
                    (Tcl_CmdProc*) W_func_select_selected_vertex, REND);
  Tcl_CreateCommand(interp, "func_select_marked_vertices",
                    (Tcl_CmdProc*) W_func_select_marked_vertices, REND);
  Tcl_CreateCommand(interp, "func_select_label",
                    (Tcl_CmdProc*) W_func_select_label, REND);
  Tcl_CreateCommand(interp, "func_clear_selection",
                    (Tcl_CmdProc*) W_func_clear_selection, REND);
  Tcl_CreateCommand(interp, "func_graph_timecourse_selection",
                    (Tcl_CmdProc*) W_func_graph_timecourse_selection, REND);
  Tcl_CreateCommand(interp, "func_print_timecourse_selection",
                    (Tcl_CmdProc*) W_func_print_timecourse_selection, REND);
  Tcl_CreateCommand(interp, "func_calc_correlation_and_write_to_overlay",
                    (Tcl_CmdProc*) W_func_calc_correlation_and_write_to_overlay, REND);
  Tcl_CreateCommand(interp, "func_normalize",
                    (Tcl_CmdProc*) W_func_normalize, REND);

  Tcl_CreateCommand(interp, "labl_load_color_table",
                    (Tcl_CmdProc*) W_labl_load_color_table, REND);
  Tcl_CreateCommand(interp, "labl_load",
                    (Tcl_CmdProc*) W_labl_load, REND);
  Tcl_CreateCommand(interp, "labl_save",
                    (Tcl_CmdProc*) W_labl_save, REND);
  Tcl_CreateCommand(interp, "labl_save_all",
                    (Tcl_CmdProc*) W_labl_save_all, REND);
  Tcl_CreateCommand(interp, "labl_import_annotation",
                    (Tcl_CmdProc*) W_labl_import_annotation, REND);
  Tcl_CreateCommand(interp, "labl_export_annotation",
                    (Tcl_CmdProc*) W_labl_export_annotation, REND);
  Tcl_CreateCommand(interp, "labl_new_from_marked_vertices",
                    (Tcl_CmdProc*) W_labl_new_from_marked_vertices, REND);
  Tcl_CreateCommand(interp, "labl_mark_vertices",
                    (Tcl_CmdProc*) W_labl_mark_vertices, REND);
  Tcl_CreateCommand(interp, "labl_select",
                    (Tcl_CmdProc*) W_labl_select, REND);
  Tcl_CreateCommand(interp, "labl_set_name_from_table",
                    (Tcl_CmdProc*) W_labl_set_name_from_table, REND);
  Tcl_CreateCommand(interp, "labl_set_info",
                    (Tcl_CmdProc*) W_labl_set_info, REND);
  Tcl_CreateCommand(interp, "labl_set_color",
                    (Tcl_CmdProc*) W_labl_set_color, REND);
  Tcl_CreateCommand(interp, "labl_remove",
                    (Tcl_CmdProc*) W_labl_remove, REND);
  Tcl_CreateCommand(interp, "labl_threshold",
                    (Tcl_CmdProc*) W_labl_threshold, REND);
  Tcl_CreateCommand(interp, "labl_remove_all",
                    (Tcl_CmdProc*) W_labl_remove_all, REND);
  Tcl_CreateCommand(interp, "labl_erode",
                    (Tcl_CmdProc*) W_labl_erode, REND);
  Tcl_CreateCommand(interp, "labl_dilate",
                    (Tcl_CmdProc*) W_labl_dilate, REND);
  Tcl_CreateCommand(interp, "labl_fill_holes",
                    (Tcl_CmdProc*) W_labl_fill_holes, REND);
  Tcl_CreateCommand(interp, "labl_select_label_by_vno",
                    (Tcl_CmdProc*) W_labl_select_label_by_vno, REND);
  Tcl_CreateCommand(interp, "labl_print_list",
                    (Tcl_CmdProc*) W_labl_print_list, REND);
  Tcl_CreateCommand(interp, "labl_print_table",
                    (Tcl_CmdProc*) W_labl_print_table, REND);

  Tcl_CreateCommand(interp, "path_select",
                    (Tcl_CmdProc*) W_path_select, REND);
  Tcl_CreateCommand(interp, "path_new_path_from_marked_vertices",
                    (Tcl_CmdProc*) W_path_new_path_from_marked_vertices, REND);
  Tcl_CreateCommand(interp, "path_remove_selected_path",
                    (Tcl_CmdProc*) W_path_remove_selected_path, REND);
  Tcl_CreateCommand(interp, "path_mark_selected_path",
                    (Tcl_CmdProc*) W_path_mark_selected_path, REND);
  Tcl_CreateCommand(interp, "path_save",
                    (Tcl_CmdProc*) W_path_save, REND);
  Tcl_CreateCommand(interp, "path_load",
                    (Tcl_CmdProc*) W_path_load, REND);

  Tcl_CreateCommand(interp, "fill_flood_from_cursor",
                    (Tcl_CmdProc*) W_fill_flood_from_cursor, REND);

  Tcl_CreateCommand(interp, "edit_vertex_at_cursor",
                    (Tcl_CmdProc*) W_edit_vertex_at_cursor, REND);

  Tcl_CreateCommand(interp, "draw_curvature_line",
                    (Tcl_CmdProc*) W_draw_curvature_line, REND);

  Tcl_CreateCommand(interp, "get_marked_vnos",
                    (Tcl_CmdProc*) W_get_marked_vnos, REND);

  Tcl_CreateCommand(interp, "get_selected_path_vnos",
                    (Tcl_CmdProc*) W_get_selected_path_vnos, REND);

  Tcl_CreateCommand(interp, "save_tiff",
                    (Tcl_CmdProc*) W_save_tiff, REND);

  Tcl_CreateCommand(interp, "flip_normals",
                    (Tcl_CmdProc*) W_flip_normals, REND);

  Tcl_CreateCommand(interp, "mark_contiguous_vertices_over_thresh",
                    (Tcl_CmdProc*) W_mark_contiguous_vertices_over_thresh, REND);

  Tcl_CreateCommand(interp, "mark_contiguous_vertices_with_similar_curvature",
                    (Tcl_CmdProc*) W_mark_contiguous_vertices_with_similar_curvature, REND);

  Tcl_CreateCommand(interp, "rip_all_vertices_except_contiguous_upripped",
                    (Tcl_CmdProc*) W_rip_all_vertices_except_contiguous_upripped, REND);

  Tcl_CreateCommand(interp, "cptn_set_format_string",
                    (Tcl_CmdProc*) W_cptn_set_format_string, REND);

  /* end rkt */
  /*=======================================================================*/
  /***** link global surfer BOOLEAN variables to tcl equivalents */
  Tcl_LinkVar(interp,"use_vertex_arrays",
              (char *)&use_vertex_arrays, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"use_display_lists",
              (char *)&use_display_lists, TCL_LINK_BOOLEAN);

  Tcl_LinkVar(interp,"dlat",(char *)&dlat, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"momentumflag",(char *)&momentumflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"expandflag",(char *)&expandflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"MRIflag",(char *)&MRIflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"sulcflag",(char *)&sulcflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"avgflag",(char *)&avgflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"areaflag",(char *)&areaflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"complexvalflag",
              (char *)&complexvalflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"statflag",(char *)&statflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"overlayflag",(char *)&overlayflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"ncthreshflag",(char *)&ncthreshflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"verticesflag",(char *)&verticesflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"shearvecflag",(char *)&shearvecflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"normvecflag",(char *)&normvecflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"movevecflag",(char *)&movevecflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"autoscaleflag",(char *)&autoscaleflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"revfsflag",(char *)&revfsflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"revphaseflag",(char *)&revphaseflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"invphaseflag",(char *)&invphaseflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"rectphaseflag",(char *)&rectphaseflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"truncphaseflag",
              (char *)&truncphaseflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"fieldsignflag",(char *)&fieldsignflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"flag2d",(char *)&flag2d, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"scalebarflag",(char *)&scalebarflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"colscalebarflag",(char *)&colscalebarflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"colscalebartextflag",(char *)&colscalebartextflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"colscalebartickflag",(char *)&colscalebartickflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"colscalebaruselabelsflag",(char *)&colscalebaruselabelsflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"colscalebarvertflag",(char *)&colscalebarvertflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"pointsflag",(char *)&pointsflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"surfaceflag",(char *)&surfaceflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"phasecontourflag",(char *)&phasecontourflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"curvimflag",(char *)&curvimflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"doublebufferflag",(char *)&doublebufferflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"surface_compiled",(char *)&surface_compiled,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"openglwindowflag",(char *)&openglwindowflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"blinkflag",(char *)&blinkflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"renderoffscreen",(char *)&renderoffscreen,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"transform_loaded",(char *)&transform_loaded,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"blackcursorflag",(char *)&blackcursorflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"bigcursorflag",(char *)&bigcursorflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"vrml2flag",(char *)&vrml2flag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"showorigcoordsflag",(char *)&showorigcoordsflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"phasecontourmodflag",(char *)&phasecontourmodflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"promptflag",(char *)&promptflag,TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"followglwinflag",(char *)&followglwinflag,
              TCL_LINK_BOOLEAN);
  /* begin rkt */
  Tcl_LinkVar(interp,"curvflag",(char *)&curvflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"mouseoverflag",(char *)&mouseoverflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"simpledrawmodeflag",
              (char *)&simpledrawmodeflag, TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"redrawlockflag",(char *)&redrawlockflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"selectlabelflag",(char *)&labl_select_flag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"drawlabelflag",(char *)&labl_draw_flag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"forcegraycurvatureflag",(char *)&forcegraycurvatureflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"drawcursorflag",(char *)&drawcursorflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"ignorezeroesinhistogramflag",
              (char *)&ignorezeroesinhistogramflag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"fopaqueflag",(char *)&sclv_opaque,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"labels_before_overlay_flag",
              (char *)&labels_before_overlay_flag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"cptn_draw_flag",(char *)&cptn_draw_flag,
              TCL_LINK_BOOLEAN);
  Tcl_LinkVar(interp,"long_config_overlay",(char *)&long_config_overlay,
              TCL_LINK_BOOLEAN);
  /* end rkt */
  /*=======================================================================*/
  /***** link global surfer INT variables to tcl equivalents */
  Tcl_LinkVar(interp,"nperdip",(char *)&sol_nperdip, TCL_LINK_INT);
  Tcl_LinkVar(interp,"numvertices",(char *)&numvertices, TCL_LINK_INT);
  Tcl_LinkVar(interp,"scrsaveflag",(char *)&scrsaveflag, TCL_LINK_INT);
  Tcl_LinkVar(interp,"surfcolor",(char *)&surfcolor, TCL_LINK_INT);
  Tcl_LinkVar(interp,"mingm",(char *)&mingm, TCL_LINK_INT);
  Tcl_LinkVar(interp,"colscale",(char *)&colscale, TCL_LINK_INT);
  Tcl_LinkVar(interp,"ilat",(char *)&ilat, TCL_LINK_INT);
  Tcl_LinkVar(interp,"mesh_linewidth",(char *)&mesh_linewidth, TCL_LINK_INT);
  Tcl_LinkVar(interp,"meshr",(char *)&meshr, TCL_LINK_INT);
  Tcl_LinkVar(interp,"meshg",(char *)&meshg, TCL_LINK_INT);
  Tcl_LinkVar(interp,"meshb",(char *)&meshb, TCL_LINK_INT);
  Tcl_LinkVar(interp,"scalebar_bright",(char *)&scalebar_bright, TCL_LINK_INT);
  Tcl_LinkVar(interp,"project",(char *)&project, TCL_LINK_INT);
  Tcl_LinkVar(interp,"sol_plot_type",(char *)&sol_plot_type, TCL_LINK_INT);
  Tcl_LinkVar(interp,"phasecontour_bright",(char *)&phasecontour_bright,
              TCL_LINK_INT);
  Tcl_LinkVar(interp,"blinkdelay",(char *)&blinkdelay, TCL_LINK_INT);
  Tcl_LinkVar(interp,"blinktime",(char *)&blinktime, TCL_LINK_INT);
  Tcl_LinkVar(interp,"select",(char *)&selection, TCL_LINK_INT);
  Tcl_LinkVar(interp,"linkvertexmode",(char *)&linkvertexmode,
              TCL_LINK_INT);

  /* begin rkt */
  Tcl_LinkVar(interp,"vertexset",(char *)&vset_current_set, TCL_LINK_INT);
  Tcl_LinkVar(interp,"numprestimpoints",
              (char *)&func_num_prestim_points, TCL_LINK_INT);
  Tcl_LinkVar(interp,"currentvaluefield",(char *)&sclv_current_field,
              TCL_LINK_INT);
  Tcl_LinkVar(interp,"ftimepoint",(char *)&sclv_cur_timepoint, TCL_LINK_INT);
  Tcl_LinkVar(interp,"fcondition",(char *)&sclv_cur_condition, TCL_LINK_INT);
  Tcl_LinkVar(interp,"fnumtimepoints",
              (char *)&sclv_num_timepoints, TCL_LINK_INT);
  Tcl_LinkVar(interp,"fnumconditions",
              (char *)&sclv_num_conditions, TCL_LINK_INT);
  Tcl_LinkVar(interp,"labelstyle",(char *)&labl_draw_style, TCL_LINK_INT);
  Tcl_LinkVar(interp,"labeloutlinered",
              (char *)&labl_outline_color[0], TCL_LINK_INT);
  Tcl_LinkVar(interp,"labeloutlinegreen",
              (char *)&labl_outline_color[1], TCL_LINK_INT);
  Tcl_LinkVar(interp,"labeloutlineblue",
              (char *)&labl_outline_color[2], TCL_LINK_INT);
  Tcl_LinkVar(interp,"func_graph_avg_mode",(char *)&func_graph_avg_mode,
              TCL_LINK_INT);
  Tcl_LinkVar(interp,"colscalebar_font_size",(char *)&colscalebar_font_size,
              TCL_LINK_INT);
  /* end rkt */
  /*=======================================================================*/
  /***** link global surfer DOUBLE variables to tcl equivalents (were float) */
  Tcl_LinkVar(interp,"decay",(char *)&decay, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fthreshmax",(char *)&fthreshmax, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fthresh",(char *)&fthresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fslope",(char *)&fslope, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fcurv",(char *)&fcurv, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"foffset",(char *)&foffset, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fmid",(char *)&fmid, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"cslope",(char *)&cslope, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"cmid",(char *)&cmid, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"cmax",(char *)&cmax, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"cmin",(char *)&cmin, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"mslope",(char *)&mslope, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"mmid",(char *)&mmid, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"whitemid",(char *)&whitemid, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"graymid",(char *)&graymid, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"mstrength",(char *)&mstrength, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"stressthresh",(char *)&stressthresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"icstrength",(char *)&icstrength, TCL_LINK_DOUBLE);
  /*Tcl_LinkVar(interp,"dipscale",(char *)&dipscale, TCL_LINK_DOUBLE);*/
  Tcl_LinkVar(interp,"angle_cycles",(char *)&angle_cycles, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"angle_offset",(char *)&angle_offset, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"wt",(char *)&wt, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"wa",(char *)&wa, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ws",(char *)&ws, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"wn",(char *)&wn, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"wc",(char *)&wc, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"wsh",(char *)&wsh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"wbn",(char *)&wbn, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"update",(char *)&update, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"decay",(char *)&decay, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"cthk",(char *)&cthk, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"light0",(char *)&light0_br, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"light1",(char *)&light1_br, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"light2",(char *)&light2_br, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"light3",(char *)&light3_br, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"offset",(char *)&offset, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"blufact",(char *)&blufact, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"cvfact",(char *)&cvfact, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fadef",(char *)&fadef, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"ncthresh",(char *)&ncthresh, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"dipavg",(char *)&dipavg, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"scalebar_xpos",(char *)&scalebar_xpos, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"scalebar_ypos",(char *)&scalebar_ypos, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"colscalebar_height",(char *)&colscalebar_height,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"colscalebar_width",(char *)&colscalebar_width,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"colscalebar_xpos",(char *)&colscalebar_xpos,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"colscalebar_ypos",(char *)&colscalebar_ypos,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_lat0",(char *)&sol_lat0,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_lat1",(char *)&sol_lat1,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_pthresh",(char *)&sol_pthresh,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_pslope",(char *)&sol_pslope,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_maxrat",(char *)&sol_maxrat,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_baseline_period",(char *)&sol_baseline_period,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_baseline_end",(char *)&sol_baseline_end,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_loflim",(char *)&sol_loflim,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_hiflim",(char *)&sol_hiflim,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"sol_snr_rms",(char *)&sol_snr_rms,TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"phasecontour_min",(char *)&phasecontour_min,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"phasecontour_max",(char *)&phasecontour_max,
              TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"dip_spacing",(char *)&dip_spacing,TCL_LINK_DOUBLE);

  /* begin rkt */
  Tcl_LinkVar(interp,"timeresolution",
              (char *)&func_time_resolution, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"falpha",(char *)&sclv_overlay_alpha, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fmin",(char *)&sclv_value_min, TCL_LINK_DOUBLE);
  Tcl_LinkVar(interp,"fmax",(char *)&sclv_value_max, TCL_LINK_DOUBLE);
  /* end rkt */

  /*=======================================================================*/
  /***** link global malloc'ed STRING vars */
  Tcl_LinkVar(interp,"spherereg_contra",        (char *)&sphere_reg_contra,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"spherereg",        (char *)&sphere_reg,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"home",        (char *)&subjectsdir,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subject",     (char *)&pname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"hemi",        (char *)&stem,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"ext",         (char *)&ext,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"curv",        (char *)&cfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"sulc",        (char *)&kfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"patch",       (char *)&pfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"label",       (char *)&lfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"fs",          (char *)&fsfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"fm",          (char *)&fmfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"dip",         (char *)&dfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"dec",         (char *)&sfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"num_rgbdir",  (char *)&gfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"named_rgbdir",(char *)&sgfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"rgb",         (char *)&agfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"insurf",      (char *)&ifname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"outsurf",     (char *)&ofname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"targsurf",    (char *)&if2name,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"targcurv",    (char *)&cf2name,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"targcurvim",  (char *)&cif2name,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"val",         (char *)&vfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"annot",       (char *)&nfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"area",        (char *)&afname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"session",     (char *)&srname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"script",      (char *)&rfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"origcoords",  (char *)&orfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"paintcoords",  (char *)&paint_fname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"ellcoords",   (char *)&elfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"vrmlsurf",    (char *)&vrfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"subjtmpdir",  (char *)&tfname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"transform",   (char *)&xffname,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"colscalebar_label1",(char *)&colscalebar_label[0],
              TCL_LINK_STRING);
  Tcl_LinkVar(interp,"colscalebar_label2",(char *)&colscalebar_label[1],
              TCL_LINK_STRING);
  Tcl_LinkVar(interp,"colscalebar_label3",(char *)&colscalebar_label[2],
              TCL_LINK_STRING);
  Tcl_LinkVar(interp,"colscalebar_label4",(char *)&colscalebar_label[3],
              TCL_LINK_STRING);

  /* begin rkt */
  Tcl_LinkVar(interp,"colortablename",
              (char *)&labl_color_table_name,TCL_LINK_STRING);
  Tcl_LinkVar(interp,"captionformat",
              (char *)&cptn_format_string,TCL_LINK_STRING);
  /* end rkt */

  /*=======================================================================*/

  strcpy(rfname,script_tcl);  /* save in global (malloc'ed in Surfer) */
  if (aliasflag)
  {
    sprintf(alias_tcl,"%s/lib/tcl/%s",envptr,argv[5]);
    Tcl_SetVar(interp,"aliasedscript",alias_tcl,0);
  }

  /* run tcl/tk startup script to set vars, make interface; no display yet */
  printf("surfer: using interface %s\n",tksurfer_tcl);
  code = Tcl_EvalFile(g_interp, tksurfer_tcl);
  if (code != TCL_OK)
    printf("Error sourcing %s:\n\t%s\n", tksurfer_tcl,
           Tcl_GetStringResult(interp));

  /* begin rkt */

  /* disable certain menu sets */
  enable_menu_set (MENUSET_VSET_INFLATED_LOADED, 0);
  enable_menu_set (MENUSET_VSET_PIAL_LOADED, 0);
  if (NULL == func_timecourse)
    enable_menu_set (MENUSET_TIMECOURSE_LOADED, 0);
  if (!overlayflag)
    enable_menu_set (MENUSET_OVERLAY_LOADED, 0);
  enable_menu_set (MENUSET_CURVATURE_LOADED, 0);
  enable_menu_set (MENUSET_LABEL_LOADED, 0);
  enable_menu_set (MENUSET_FIELDSIGN_LOADED, 0);
  enable_menu_set (MENUSET_FIELDMASK_LOADED, 0);

  /* Now that we've sourced the interface file, send our cached
     commands. */
  send_cached_tcl_commands ();

  /* Try to read the white coords. These are used to load labels. */
  read_white_vertex_coordinates();

  /* end rkt */

  /* Update all the linked var sets with our values here in case they
     are different from the defaults set in tksurfer.tcl */
  sprintf (tcl_cmd, "UpdateLinkedVarGroup all");
  send_tcl_command (tcl_cmd);

  /* if command line script exists, now run as batch job (possibly exiting) */
  if (scriptok)
  {    /* script may or may not open gl window */
    printf("tksurfer: run tcl script: %s\n",script_tcl);
    code = Tcl_EvalFile (g_interp, script_tcl);
    if (code != TCL_OK)  printf("%s",Tcl_GetStringResult(interp));
  }
  else
  {
    ; /* surfer has already opened gl window if no script */
  }

  /* always start up command line shell too (if script doesn't exit) */
#ifndef Windows_NT
  Tk_CreateFileHandler(0, TK_READABLE, StdinProc, (ClientData) 0);
#endif
  if (tty) Prompt(interp, 0);
  fflush(stdout);
  Tcl_DStringInit(&command);
  Tcl_ResetResult(interp);

  /*Tk_MainLoop();*/  /* standard */

  if (statflag)   /* uggh rkt and brf */
  {
    char cmd[STRLEN] ;

    statflag = 0 ;  /* I'm really, really sorry */
    sclv_set_current_field(SCLV_VALSTAT) ;
    sprintf (cmd, "ShowValueLabel %d 1", SCLV_VALSTAT);
    send_tcl_command (cmd);
    sprintf (cmd, "UpdateLinkedVarGroup view");
    send_tcl_command (cmd);
  }

#ifdef USE_XGLUT_WINDOW
  glutMainLoop(); /* never returns */
#else /* USE_XGLUT_WINDOW */
  /* dual event loop (interface window made now) */
  while (tk_NumMainWindows > 0)
  {
    while (Tk_DoOneEvent(TK_ALL_EVENTS|TK_DONT_WAIT))
    {
      /* do all the tk events; non-blocking */
    }
    do_one_gl_event(interp);

#ifndef sgi
    tv.tv_sec = 0;
    tv.tv_usec = 10000;
    select(0, NULL, NULL, NULL, &tv);
#else
  sginap((long)1);   /* block for 10 msec */
#endif
  }
#endif  /* USE_XGLUT_WINDOW */

  send_tcl_command("exit");
  exit(0);
  return(0) ;   /* for ansi */
}

#ifndef Windows_NT
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
  if (count <= 0)
  {
    if (!gotPartial)
    {
      if (tty)
      {
        send_tcl_command("exit");
        exit(1);
      }
      else
      {
        Tk_DeleteFileHandler(0);
      }
      return;
    }
    else count = 0;
  }
  cmd = Tcl_DStringAppend(&command, input, count);
  if (count != 0)
  {
    if ((input[count-1] != '\n') && (input[count-1] != ';'))
    {
      gotPartial = 1;
      goto prompt;
    }
    if (!Tcl_CommandComplete(cmd))
    {
      gotPartial = 1;
      goto prompt;
    }
  }
  gotPartial = 0;
  Tk_CreateFileHandler(0, 0, StdinProc, (ClientData) 0);
  code = Tcl_RecordAndEval (interp, cmd, TCL_EVAL_GLOBAL);
  Tk_CreateFileHandler (0, TK_READABLE, StdinProc, (ClientData) 0);
  Tcl_DStringFree(&command);
  if (*interp->result != 0)
    if ((code != TCL_OK) || (tty))
      puts(interp->result);
prompt:
  if (tty)  Prompt(interp, gotPartial);
  Tcl_ResetResult(interp);
}
#endif // Windows_NT

/*=== from TkMain.c ===================================================*/
static void Prompt(interp, partial)
Tcl_Interp *interp;
int partial;
{
  char *promptCmd;
  int code;

  promptCmd =
    (char*)Tcl_GetVar(interp,
                      partial ? (char*)"tcl_prompt2" : (char*)"tcl_prompt1",
                      TCL_GLOBAL_ONLY);
  if (promptCmd == NULL)
  {
defaultPrompt:
    if (!partial)
      fputs("% ", stdout);
  }
  else
  {
    code = Tcl_Eval (g_interp, promptCmd);
    if (code != TCL_OK)
    {
      Tcl_AddErrorInfo(interp,
                       "\n    (script that generates prompt)");
      fprintf(stderr, "%s\n", Tcl_GetStringResult(interp));
      goto defaultPrompt;
    }
  }
  fflush(stdout);
}

static void
grabPixels(unsigned int width, unsigned int height, unsigned short *red,
           unsigned short *green, unsigned short *blue)
{
  GLint    swapbytes, lsbfirst, rowlength ;
  GLint    skiprows, skippixels, alignment ;

  glGetIntegerv(GL_UNPACK_SWAP_BYTES, &swapbytes) ;
  glGetIntegerv(GL_UNPACK_LSB_FIRST, &lsbfirst) ;
  glGetIntegerv(GL_UNPACK_ROW_LENGTH, &rowlength) ;
  glGetIntegerv(GL_UNPACK_SKIP_ROWS, &skiprows) ;
  glGetIntegerv(GL_UNPACK_SKIP_PIXELS, &skippixels) ;
  glGetIntegerv(GL_UNPACK_ALIGNMENT, &alignment) ;

  glPixelStorei(GL_UNPACK_SWAP_BYTES, GL_FALSE) ;
  glPixelStorei(GL_UNPACK_ROW_LENGTH, 0) ;
  glPixelStorei(GL_UNPACK_SKIP_ROWS, 0) ;
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, 0) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1) ;

  glReadPixels(0, 0, width, height,GL_RED,GL_UNSIGNED_SHORT, (GLvoid *)red);
  glReadPixels(0, 0,width,height,GL_GREEN,GL_UNSIGNED_SHORT,(GLvoid *)green);
  glReadPixels(0, 0, width, height,GL_BLUE,GL_UNSIGNED_SHORT,(GLvoid *)blue);

  glPixelStorei(GL_UNPACK_SWAP_BYTES, swapbytes) ;
  glPixelStorei(GL_UNPACK_LSB_FIRST, lsbfirst) ;
  glPixelStorei(GL_UNPACK_ROW_LENGTH, rowlength) ;
  glPixelStorei(GL_UNPACK_SKIP_ROWS, skiprows) ;
  glPixelStorei(GL_UNPACK_SKIP_PIXELS, skippixels) ;
  glPixelStorei(GL_UNPACK_ALIGNMENT, alignment) ;
}
static void
save_rgbfile(char *fname, int width, int height, unsigned short *red,
             unsigned short *green, unsigned short *blue)
{
  RGB_IMAGE  *image ;
  int    y ;
  unsigned short *r, *g, *b ;

#ifdef IRIX
  image = iopen(fname,"w",RLE(1), 3, width, height, 3);
#else
  image = iopen(fname,"w",UNCOMPRESSED(1), 3, width, height, 3);
#endif
  if (!image)
    return ;
  for (y = 0 ; y < height; y++)
  {
    r = red + y * width ;
    g = green + y * width ;
    b = blue + y * width ;

    /* fill rbuf, gbuf, and bbuf with pixel values */
    putrow(image, r, y, 0);    /* red row */
    putrow(image, g, y, 1);    /* green row */
    putrow(image, b, y, 2);    /* blue row */
  }
  iclose(image);
}

void
fix_nonzero_vals(void)
{
  int    vno ;
  VERTEX *v ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (!FZERO(v->val))
      v->fixedval = TRUE ;
    else
      v->fixedval = FALSE ;
  }
}

void
curv_to_val(void)
{
  int    vno ;
  VERTEX *v ;
  char   label[1024];
  double min, max;

  /* Initialize our field. */
  strncpy (label, "From curv", sizeof(label));
  sclv_new_empty (SCLV_VAL, label);

  /* Copy the data in. */
  min = 1000000000;
  max = -min;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    sclv_set_value (v, SCLV_VAL, v->curv);

    /* Update the min and max. */
    if (v->curv < min)
      min = v->curv;
    if (v->curv > max)
      max = v->curv;
  }

  /* Set the min and max. */
  sclv_field_info[SCLV_VAL].min_value = min;
  sclv_field_info[SCLV_VAL].max_value = max;

  /* Calc the frquencies */
  sclv_calc_frequencies (SCLV_VAL);

  /* Request a redraw. Turn on the overlay flag and select this value
     set. */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (SCLV_VAL);
}
void
scale_vals(float scale)
{
  int    vno ;
  VERTEX *v ;
  float  min, max ;

  min = 1e10 ; max = -min ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    sclv_get_value (v, SCLV_VAL, &(v->val));
    v->val *= scale ;
    sclv_set_value (v, SCLV_VAL, v->val);

    /* Update the min and max. */
    if (v->val < min)
      min = v->val;
    if (v->val > max)
      max = v->val;
  }

  /* Set the min and max. */
  sclv_field_info[SCLV_VAL].min_value = min;
  sclv_field_info[SCLV_VAL].max_value = max;

  /* Calc the frquencies */
  sclv_calc_frequencies (SCLV_VAL);

  /* Request a redraw. Turn on the overlay flag and select this value
     set. */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (SCLV_VAL);
}
int
read_parcellation(char *parc_name, char *lut_name)
{
  return(read_and_smooth_parcellation(parc_name, lut_name, 25, 25)) ;
}
int
read_and_smooth_parcellation(char *parc_name, char *lut_name,
                             int siter, int miter)
{
  char   *cp, fname[STRLEN], path[STRLEN], name[STRLEN] ;
  MRI    *mri ;
  int    vno, index, rd, gn, bl, xv, yv, zv ;
  FILE   *fp ;
  char   r[256], g[256], b[256], line[STRLEN] ;
  VERTEX *v ;
  double x, y, z ;

  MRISclearMarks(mris) ;
  cp = strchr(parc_name, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/../mri/%s", path, parc_name) ;
  } else
    strcpy(fname, parc_name) ;  /* path specified explcitly */
  fprintf(stderr, "reading parcellation from %s...\n", fname) ;
  mri = MRIread(fname) ;
  if (!mri)
  {
    fprintf(stderr, "### could not read parcellation file %s\n", fname) ;
    return(Gerror) ;
  }

  cp = strchr(lut_name, '/') ;
  if (!cp)                 /* no path - use same one as mris was read from */
  {
    FileNamePath(mris->fname, path) ;
    sprintf(fname, "%s/../label/%s", path, lut_name) ;
  } else
    strcpy(fname, lut_name) ;  /* path specified explcitly */

  fprintf(stderr, "reading color lut from %s...\n", fname) ;
  fp = fopen(fname, "r") ;
  if (!fp)
  {
    fprintf(stderr, "### could not read parcellation lut from %s\n", fname) ;
    return(Gerror) ;
  }


  while ((cp = fgetl(line, 199, fp)) != NULL)
  {
    sscanf(cp, "%d %s %d %d %d %*s\n", &index, name, &rd, &gn, &bl) ;
    r[index] = (char)rd ;
    g[index] = (char)gn ;
    b[index] = (char)bl ;
    parc_names[index] = (char *)calloc(strlen(name)+1, sizeof(char)) ;
    strcpy(parc_names[index], name) ;
    parc_red[index] = (char)rd ;
    parc_green[index] = (char)gn ;
    parc_blue[index] = (char)bl ;
    parc_flag = 1 ;
  }
  fclose(fp) ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  read_orig_vertex_coordinates("pial") ;
  MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
  fprintf(stderr, "painting parcellation onto surface...\n") ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    MRISvertexToVoxel(mris,v, mri, &x, &y, &z) ;
    xv = nint(x) ;
    yv = nint(y) ;
    zv = nint(z) ;
    xv = MAX(0,xv) ;
    yv = MAX(0,yv) ;
    zv = MAX(0,zv) ;
    xv = MIN(mri->width-1, xv) ;
    yv = MIN(mri->height-1, yv) ;
    zv = MIN(mri->depth-1, zv) ;
    index = MRIvox(mri, xv, yv, zv) ;
    if (index < 1)
      continue ;
    v->marked = 1 ;
    v->val = (float)index ;
    rd = r[index] ;
    gn = g[index] ;
    bl = b[index] ;
    v->annotation = rd + (gn << 8) + (bl << 16) ;
  }

  MRISsoapBubbleVals(mris, siter) ;
  fprintf(stderr, "applying mode filter...\n") ;
  MRISmodeFilterVals(mris, miter) ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    index = nint(v->val) ;
    rd = r[index] ;
    gn = g[index] ;
    bl = b[index] ;
    v->annotation = rd + (gn << 8) + (bl << 16) ;
  }
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;

  PR ;
  annotationloaded = TRUE ;
  MRIfree(&mri) ;
  MRISclearMarks(mris) ;
  return(NO_ERROR) ;
}
int
read_curv_to_val(char *fname)
{
  int    vno ;
  VERTEX *v ;

  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->tx = v->curv ;
  }
  surface_compiled = 0 ;
  if (MRISreadCurvatureFile(mris, fname) != NO_ERROR)
    return(Gerror) ;

  printf("surfer: values read: min=%f max=%f\n",
         mris->min_curv,mris->max_curv);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = v->curv ;
    v->curv = v->tx ;
  }
  return(NO_ERROR) ;
}

void
val_to_curv(void)
{
  int    vno ;
  VERTEX *v ;

  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->curv = v->val ;
  }
}

void
val_to_stat(void)
{
  int    vno ;
  VERTEX *v ;

  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->stat = v->val ;
  }
}

void
stat_to_val(void)
{
  int    vno ;
  VERTEX *v ;

  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    v->val = v->stat ;
  }
}

void
show_flat_regions(char *surf_name, double thresh)
{
  int    vno, nfound = 0 ;
  VERTEX *v ;
  float  mx ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  if (MRISreadVertexPositions(mris, surf_name) != NO_ERROR)
  {
    ErrorPrintf(ERROR_NOFILE,
                "show_flat_regions(%s): could not read surface",
                surf_name) ;
    return ;
  }

  MRIScomputeMetricProperties(mris) ;  /* compute normals */

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    mx = MAX(MAX(fabs(v->nx), fabs(v->ny)), fabs(v->nz)) ;
    if (mx > thresh)
    {
      nfound++ ;
      v->val = mx ;
    }
    else
      v->val = 0 ;
  }
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  MRIScomputeMetricProperties(mris) ;
  printf("%d vertices marked...\n", nfound) ;
  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);
}

static void
resize_brain(float surface_area)
{
  float scale ;

  MRIScomputeMetricProperties(mris) ;
  scale = sqrt(surface_area / mris->total_area) ;
  MRISscaleBrain(mris, mris, scale) ;
  MRIScomputeMetricProperties(mris) ;
  vset_save_surface_vertices(VSET_MAIN) ;
  vset_set_current_set(vset_current_set) ;
  redraw() ;
}

void
val_to_mark(void)
{
  int    vno ;
  VERTEX *v ;
  static int mark = 2 ;

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    if (fabs(v->val) > fthresh && v->marked < 2)
      v->marked = mark ;
  }
  mark++ ;
}

void
twocond(int c0, int c1)
{
  char  fname[STRLEN], val_name[STRLEN], *cp ;
  int   dof ;
  FILE  *fp ;

  cond0 = c0 ;
  cond1 = c1 ;

  /* put the averages in val and val2, and the variances in valbak, val2bak */
  twocond_flag = 1 ;
  sprintf(fname, "%s/sigavg%d-%s.w", val_dir, c1, stem) ;
  fprintf(stderr, "reading cond %d means from %s\n", c1, fname) ;
  surface_compiled = 0 ;
  MRISreadValues(mris, fname) ;
  MRIScopyValToVal2(mris) ;

  sprintf(fname, "%s/sigavg%d.dof", val_dir, c1) ;
  fp = fopen(fname, "r") ;
  if (!fp)
  {
    printf("### surfer: could not open dof file %s\n", fname) ;
    return ;
  }
  fscanf(fp, "%d", &dof) ;
  fclose(fp) ;
  sprintf(fname, "%s/sigvar%d-%s.w", val_dir, c1, stem) ;
  fprintf(stderr, "reading cond %d variances from %s and scaling by %d\n",
          c1, fname, dof) ;
  MRISreadValues(mris, fname) ;
  /* turn squared standard errors into variances */
  MRISmulVal(mris, (float)dof) ;
  MRISsqrtVal(mris) ;
  MRIScopyValToVal2Bak(mris) ;

  sprintf(fname, "%s/sigavg%d.dof", val_dir, c0) ;
  fp = fopen(fname, "r") ;
  if (!fp)
  {
    printf("### surfer: could not open dof file %s\n", fname) ;
    return ;
  }
  fscanf(fp, "%d", &dof) ;
  fclose(fp) ;
  sprintf(fname, "%s/sigvar%d-%s.w", val_dir, c0, stem) ;
  fprintf(stderr, "reading cond %d variances from %s and scaling by %d\n",
          c0, fname,dof) ;
  MRISreadValues(mris, fname) ;
  /* turn squared standard errors into variances */
  MRISmulVal(mris, (float)dof) ;
  MRISsqrtVal(mris) ;
  MRIScopyValToValBak(mris) ;

  sprintf(fname, "%s/sigavg%d-%s.w", val_dir, c0, stem) ;
  fprintf(stderr, "reading condition %d means from %s\n",
          c0, fname) ;
  MRISreadValues(mris, fname) ;
  sprintf(val_name, "GROUP%d", c0) ;
  cp = getenv(val_name) ;
  if (cp)
    sprintf(val_name, "%s mean", cp) ;
  else
    sprintf(val_name, "group %d mean", c0) ;
  set_value_label_name(val_name, SCLV_VAL) ;
  if (cp)
    sprintf(val_name, "%s std", cp) ;
  else
    sprintf(val_name, "group %d std", c0) ;
  set_value_label_name(val_name, SCLV_VALBAK) ;
  sprintf(val_name, "GROUP%d", c1) ;
  cp = getenv(val_name) ;
  if (cp)
    sprintf(val_name, "%s mean", cp) ;
  else
    sprintf(val_name, "group %d mean", c1) ;
  set_value_label_name(val_name, SCLV_VAL2) ;
  if (cp)
    sprintf(val_name, "%s std", cp) ;
  else
    sprintf(val_name, "group %d std", c1) ;
  set_value_label_name(val_name, SCLV_VAL2BAK) ;
}

int
mask_label(char *label_name)
{
  LABEL  *area ;
  int     vno ;
  VERTEX  *v ;

  surface_compiled = 0 ;
  MRISclearMarks(mris) ;
  area = LabelRead(pname, label_name) ;
  if (!area)
  {
    fprintf(stderr, "unable to read label file %s\n", label_name) ;
    return(NO_ERROR) ;
  }

  LabelMarkSurface(area, mris) ;  /* mark all points in label */

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->marked || v->ripflag)
      continue ;
    v->stat = v->val = v->imag_val = v->val2 = v->valbak = v->val2bak = 0.0 ;
  }
  LabelFree(&area) ;
  MRISclearMarks(mris) ;
  redraw() ;
  return(NO_ERROR) ;
}

#include "mrishash.h"
static void
deconvolve_weights(char *weight_fname, char *scale_fname)
{
  MHT     *mht ;
  MHBT    *bucket ;
  MHB     *bin ;
  int     vno, n, i ;
  VERTEX  *v, *vn ;
  double  angle, circumference, norm, sigma_sq, wt, mean, var, mn, mx,
  sphere_area ;
  float   radius, dist, sigma, x0, y0, z0, dscale, dx, dy, dz ;
  VECTOR  *v1, *v2 ;


  surface_compiled = 0 ;
  if (is_val_file(weight_fname))
  {
    if (MRISreadValues(mris,weight_fname) != NO_ERROR)
      return ;
  }
  else if (read_curv_to_val(weight_fname) != NO_ERROR)
    return ;

  fprintf(stderr, "loading spherical coordinate system...\n") ;
  if (read_canon_vertex_coordinates("sphere") != NO_ERROR)
    return ;

  MRISsaveVertexPositions(mris, TMP_VERTICES) ;
  MRISrestoreVertexPositions(mris, CANONICAL_VERTICES) ;
  radius = MRISaverageRadius(mris) ;
  sphere_area = M_PI * radius * radius * 4.0 ;
  dscale = sqrt(mris->total_area / sphere_area) ;
  circumference = M_PI * 2.0 * radius ;
  v1 = VectorAlloc(3, MATRIX_REAL) ;
  v2 = VectorAlloc(3, MATRIX_REAL) ;
  MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
  fprintf(stderr, "average radius = %2.1f\n", radius) ;

  MRIScopyValToVal2(mris) ;   /* put weights in val2 field */
  if (MRISreadValues(mris, scale_fname) != NO_ERROR)
  {
    fprintf(stderr,
            "### surfer: could not open scale file %s\n",scale_fname);
    PR ;
    return ;
  }
  MRIScopyValToVal2Bak(mris) ;   /* put spatial scale in val2bak field */

#define MAX_DIST  (1.0*sqrt(500))

  fprintf(stderr, "building spatial LUT...\n") ;
  mht = MHTcreateVertexTable_Resolution(mris, CANONICAL_VERTICES, MAX_DIST) ;

  fprintf(stderr, "deconvolving weights...\n") ;
  MRISsetVals(mris, 0.0) ;  /* clear all values */

  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->tdx = 0 ;  /* used for normalization below */
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    if (!(vno % (mris->nvertices/50)))
      fprintf(stderr, ".") ;
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag || (fabs(v->val2) < fthresh))
      continue ;

    x0 = v->cx ;
    y0 = v->cy ;
    z0 = v->cz ;
    sigma = sqrt(v->val2bak) ;
    if (FZERO(sigma))
      norm = 1.0 ;
    else
      norm = 1.0 / (sigma * sqrt(2.0*M_PI)) ;
    sigma_sq = sigma*sigma ;
    MRISclearMarks(mris) ;

    VECTOR_LOAD(v1, x0, y0, z0) ;  /* radius vector */

    for (dx = -MAX_DIST ; dx <= MAX_DIST ; dx += MAX_DIST)
      for (dy = -MAX_DIST ; dy <= MAX_DIST ; dy += MAX_DIST)
        for (dz = -MAX_DIST ; dz <= MAX_DIST ; dz += MAX_DIST)
        {
          x0 = v->cx + dx ;
          y0 = v->cy + dy ;
          z0 = v->cz + dz ;

          bucket = MHTacqBucket(mht, x0, y0, z0) ;
          if (!bucket)
            continue ;
          for (bin = bucket->bins, i = 0 ; i < bucket->nused ; i++, bin++)
          {
            n = bin->fno ;
            vn = &mris->vertices[n] ;
            if (vn->marked)
              continue ;
            if (n == Gdiag_no)
              DiagBreak() ;
            vn->marked = 1 ;  /* only process it once */
            VECTOR_LOAD(v2, vn->cx, vn->cy, vn->cz) ;
            /* radius vector */
            angle = fabs(Vector3Angle(v1, v2)) ;
            dist = dscale * circumference * angle / (2.0 * M_PI) ;
            if (dist > 3*sigma)
              continue ;
            if (n == Gdiag_no)
              DiagBreak() ;
            if (FZERO(sigma_sq)) /* make it a unit delta function */
            {
              if (FZERO(dist))
                wt = norm ;
              else
                wt = 0.0 ;
            } else
              wt = norm*exp(-0.5 * dist*dist / sigma_sq) ;
#if 0
            if (wt < fthresh/10)
              continue ;
#endif
            if (n == Gdiag_no)
              DiagBreak() ;
            vn->val += v->val2*wt ;
            vn->tdx += wt ;            /* for normalization later */
            if (!finite(vn->val))
              DiagBreak() ;
          }
          MHTrelBucket(&bucket);
        }
  }
  fprintf(stderr, "\n") ;

  mean = var = mn = mx = 0.0f ;
  for (i = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (!FZERO(v->tdx))
    {
      i++ ;
      v->val /= v->tdx ;
      v->tdx = 0 ;
      if (v->val < mn)
        mn = v->val ;
      if (v->val > mx)
        mx = v->val ;
      mean += v->val ;
      var += v->val * v->val ;
    }
  }

  if (i)
  {
    mean /= (double)i ;
    var = var / i - mean*mean ;
    fprintf(stderr, "%d non-zero vals, [%2.2f --> %2.2f], %2.2f += %2.2f\n",
            i, mn, mx, mean, sqrt(var)) ;
  }
  overlayflag = TRUE ;
  colscale = HEAT_SCALE ;
  val_to_stat() ;
  complexvalflag = FALSE ;

  MRISclearMarks(mris) ;
  VectorFree(&v1) ;
  VectorFree(&v2) ;
  MHTfree(&mht) ;
  redraw() ;
}

static void
read_disc(char *subject_name)
{
  char   fname[300] ;
  int    vno ;
  VERTEX *v ;
  double proj ;

  surface_compiled = 0 ;
  sprintf(fname, "./%s.offset_%s", stem, subject_name) ;
  if (read_curv_to_val(fname) != NO_ERROR)
  {
    fprintf(stderr, "### surfer: could not open %s\n", fname) ;
    return ;
  }
  shift_values() ;
  sprintf(fname, "./%s.mdiff", stem) ;
  if (read_curv_to_val(fname) != NO_ERROR)
  {
    fprintf(stderr, "### surfer: could not open %s\n", fname) ;
    return ;
  }

  swap_values() ;  /* mdiff ->valbak, offset ->val2bak */

  sprintf(fname, "control_thickness/%s.thickness_%s", stem, subject_name) ;
  if (read_curv_to_val(fname) != NO_ERROR)
  {
    fprintf(stderr, "### surfer: could not open %s\n", fname) ;
    return ;
  }
  shift_values() ;    /* thickness -> val2 */
  sprintf(fname, "./%s.disc_%s", stem, subject_name) ; /* disc -> val */
  if (read_curv_to_val(fname) != NO_ERROR)
  {
    fprintf(stderr, "### surfer: could not open %s\n", fname) ;
    return ;
  }

  MRIScopyValuesToImagValues(mris) ;   /* disc -> imag_val */

  /*
    valbak       - difference between means
    val2bak      - offset (mean of two groups)
    val2         - thickness of individual
    imag_val     - discriminant weight for subject
    val          - projection of mean difference onto discriminant
  */

  disc_flag = 1 ;
  mask_label("lh-cortex") ;
  for (proj = 0.0, vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    v->val = (v->val2 - v->val2bak) * v->imag_val ;
    proj += (double)v->val ;
  }
  printf("projection onto discriminant: %2.3f\n", proj) ;
  PR
}

/* begin rkt */

static void
print_vertex_data(int vno, FILE *fp, float dmin)
{
  VERTEX   *v ;
  int      i, j, imnr ;

  v = &mris->vertices[vno] ;
  if(twocond_flag){
    if (!FZERO(v->imag_val))
      fprintf(fp, "cond %d: %2.3f +- %2.3f, cond %d: %2.3f +- %2.3f, "
              "scale=%2.0f\n",
              cond0, v->val, v->valbak, cond1, v->val2, v->val2bak,
              v->imag_val);
    else
      fprintf(fp, "cond %d: %2.3f +- %2.3f, cond %d: %2.3f +- %2.3f\n",
              cond0, v->val, v->valbak, cond1, v->val2, v->val2bak);
    PR;
  }
  else if (disc_flag){
    fprintf(fp, "v %d:\n\tdisc:\t\t%2.3f\n\tmdiff:"
            "\t\t%2.3f\n\tthickness:\t%2.3f\n\toffset:\t\t%2.3f\n"
            "\tdiff:\t\t%2.3f\n\tproj:\t\t%2.3f\n",
            vno, v->imag_val, v->valbak, v->val2, v->val2bak,
            v->val2-v->val2bak, v->val);
    PR;
  }
  else
  {
    fprintf(fp, "-----------------------------------\n");
    fprintf(fp, "selected vertex %d out of %d\n",vno,mris->nvertices);
    fprintf(fp, "current  %6.2f %6.2f %6.2f\n",v->x,v->y,v->z);
    fprintf(fp, "orig     %6.2f %6.2f %6.2f\n",v->origx,v->origy,v->origz);
    fprintf(fp, "pial     %6.2f %6.2f %6.2f\n",v->pialx,v->pialy,v->pialz);
    fprintf(fp, "white    %6.2f %6.2f %6.2f\n",v->whitex,v->whitey,v->whitez);
    fprintf(fp, "inflated %6.2f %6.2f %6.2f\n",v->infx,v->infy,v->infz);
    fprintf(fp, "flat     %6.2f %6.2f %6.2f\n",v->fx,v->fy,v->fz);
    fprintf(fp, "normals  %6.2f %6.2f %6.2f\n",v->nx,v->ny,v->nz);
    fprintf(fp, "nneighbors  %d\n",v->vtotal);
    fprintf(fp, "ripflag     %d\n",v->ripflag);
    fprintf(fp, "surfer: dmin=%3.4f, vno=%d, x=%3.4f, y=%3.4f, z=%3.4f\n",
            dmin,vno,v->x, v->y, v->z);
    fprintf(fp, "surfer: curv=%f, fs=%f\n",v->curv,v->fieldsign);
    fprintf(fp, "surfer: val=%f, val2=%f\n",v->val,v->val2);
    fprintf(fp, "surfer: amp=%f, angle=%f deg (%f)\n",hypot(v->val,v->val2),
            (float)(atan2(v->val2,v->val)*180/M_PI),
            (float)(atan2(v->val2,v->val)/(2*M_PI)));
  }
  if (annotationloaded)
  {
    int  r, g, b ;
    r = v->annotation & 0x00ff ;
    g = (v->annotation >> 8) & 0x00ff ;
    b = (v->annotation >> 16) & 0x00ff ;
    fprintf(fp, "annot = %d (%d, %d, %d)\n", v->annotation, r, g, b) ;
  }
  if (MRIflag && MRIloaded)
  {
    imnr = (int)((v->y-yy0)/st+0.5);
    i = (int)((zz1-v->z)/ps+0.5);
    j = (int)((xx1-v->x)/ps+0.5);
    fprintf(fp, "surfer: mrival=%d imnr=%d, i=%d, j=%d\n",
            im[imnr][i][j],imnr,i,j);
    PR;
  }
  if (parc_flag && v->val > 0 && parc_names[(int)nint(v->val)])
    fprintf(stderr, "parcellation name: %s\n", parc_names[(int)nint(v->val)]) ;
  if (Ggcsa != NULL)
  {
    VERTEX *v_classifier, *v_prior ;
    int     vno_classifier, vno_prior ;
    GCSA_NODE *gcsan ;
    CP_NODE *cpn ;
    float   x, y, z ;

    v = &mris->vertices[vno] ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    v->x = v->cx ;
    v->y = v->cy ;
    v->z = v->cz ;
    v_prior = GCSAsourceToPriorVertex(Ggcsa, v) ;
    v_classifier = GCSAsourceToClassifierVertex(Ggcsa, v_prior) ;
    vno_classifier = v_classifier - Ggcsa->mris_classifiers->vertices ;
    vno_prior = v_prior - Ggcsa->mris_priors->vertices ;
    gcsan = &Ggcsa->gc_nodes[vno_classifier] ;
    cpn = &Ggcsa->cp_nodes[vno_prior] ;
    dump_gcsan(gcsan, cpn, stdout, 0) ;
    v->x = x ;
    v->y = y ;
    v->z = z ;
  }
}

static void
send_current_labels()
{
  /* We've cached the current vnos of the cursor and mouseover, so
     just send those. */
  update_labels (LABELSET_MOUSEOVER, mouseover_vno, 0);
  update_labels (LABELSET_CURSOR, selection, 0);
}

static void
update_labels(int label_set, int vno, float dmin)
{
  int err;
  char command[NAME_LENGTH];
  VERTEX *v;
  int i, j, imnr;
  int r, g, b;
  float x_mni, y_mni, z_mni;
  float x_tal, y_tal, z_tal;
  char fname[NAME_LENGTH];
  float x, y, z, sx, sy, sz, rr, dd, phi, theta;
  int field;
  int label_index_array[LABL_MAX_LABELS];
  int num_labels_found;
  float value=0;

  /* make sure we an interpreter to send commands to */
  if (g_interp==NULL)
    return;

  /* Make sure we have a valid vno. */
  if ( vno < 0 || vno >= mris->nvertices )
    return;

  /* If this is our cursor label set, we'll be updating captions
     too. */
  if (0 == label_set)
    cptn_clear_value_string();

  /* get the vertex */
  v = &mris->vertices[vno];

  /* send each label value. If this is the cursor label set, update
     the caption, too. */
  sprintf(command, "UpdateLabel %d %d %d", label_set, LABEL_VERTEXINDEX, vno);
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_VERTEXINDEX, "%d", vno);

  sprintf(command, "UpdateLabel %d %d %f", label_set, LABEL_DISTANCE, dmin);
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_DISTANCE, "%f", dmin);

  sprintf(command,"UpdateLabel %d %d \"(%.2f  %.2f  %.2f)\"",
          label_set, LABEL_COORDS_RAS, v->origx, v->origy, v->origz);
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_COORDS_RAS, "(%.2f, %.2f, %.2f)",
                           v->origx, v->origy, v->origz);

  if (MRIflag && MRIloaded)
  {
    imnr = (int)((v->y-yy0)/st+0.5);
    i = (int)((zz1-v->z)/ps+0.5);
    j = (int)((xx1-v->x)/ps+0.5);
    sprintf(command, "UpdateLabel %d %d \"(%d  %d  %d)\"",
            label_set, LABEL_COORDS_INDEX, imnr, i, j);
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_COORDS_INDEX, "(%d, %d, %d)",
                             imnr, i, j);

    sprintf(command, "UpdateLabel %d %d %d",
            label_set, LABEL_MRIVALUE, im[imnr][i][j]);
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_MRIVALUE, "%d", im[imnr][i][j]);
  }
  /* if we have a tal transform, compute the tal. */
  if (transform_loaded)
  {
    conv_ras_to_tal( v->origx, v->origy, v->origz, &x_tal, &y_tal, &z_tal );
    sprintf(command, "UpdateLabel %d %d \"(%.2f  %.2f  %.2f)\"",
            label_set, LABEL_COORDS_TAL, x_tal, y_tal, z_tal );
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_COORDS_TAL, "(%.2f, %.2f, %.2f)",
                             x_tal, y_tal, z_tal);

    conv_ras_to_mnital( v->origx, v->origy, v->origz,
                        &x_mni, &y_mni, &z_mni );
    sprintf(command, "UpdateLabel %d %d \"(%.2f  %.2f  %.2f)\"",
            label_set, LABEL_COORDS_MNITAL, x_mni, y_mni, z_mni );
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_COORDS_MNITAL, "(%.2f, %.2f, %.2f)",
                             x_mni, y_mni, z_mni);
  }

  sprintf(command, "UpdateLabel %d %d \"(%.2f  %.2f  %.2f)\"",
          label_set, LABEL_COORDS_NORMAL, v->nx, v->ny, v->nz);
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_COORDS_NORMAL, "(%.2f, %.2f, %.2f)",
                           v->nx, v->ny, v->nz);

  /* if a canon surface isn't loaded, make a name and see if it exists. */
  if (canonsurfloaded == FALSE && canonsurffailed == FALSE)
  {
    sprintf(fname, "%s.%s", fpref, sphere_reg_suffix) ;
    if (FileExists(fname))
    {
      printf("surfer: reading canonical coordinates from\n");
      printf("surfer:   %s\n",fname);
    }
    else
    {
      /* set our don't load flag here so we don't
         try and load it every time */
      dontloadspherereg = TRUE;
    }
  }
  if (canonsurffailed)
    dontloadspherereg = TRUE;

  /* if the canon surface is loaded _or_ if the filename made above exists
     and we can read in the vertex coords and our don't load flag isn't
     on... */
  if (canonsurfloaded == TRUE ||
      (dontloadspherereg == FALSE &&
       FileExists(fname) &&
       read_canon_vertex_coordinates(fname) == NO_ERROR) )
  {
    x = mris->vertices[vno].cx ;
    y = mris->vertices[vno].cy ;
    z = mris->vertices[vno].cz ;
    sx = x*e0[0] + y*e0[1] + z*e0[2] ;
    sy = x*e1[0] + y*e1[1] + z*e1[2] ;
    sz = x*e2[0] + y*e2[1] + z*e2[2] ;
    x = sx ;
    y = sy ;
    z = sz ;

    rr = sqrt(x*x + y*y + z*z) ;
    dd = rr*rr-z*z ;
    if (dd < 0.0)
      dd = 0.0 ;

    phi = atan2(sqrt(dd), z) ;
    theta = atan2(y/rr, x/rr) ;

    sprintf(command, "UpdateLabel %d %d \"(%.2f  %.2f  %.2f)\"",
            label_set, LABEL_COORDS_SPHERE_XYZ, sx, sy, sz );
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_COORDS_SPHERE_XYZ, "(%.2f, %.2f, %.2f)",
                             sx, sy, sz);

    sprintf(command, "UpdateLabel %d %d \"(%2.1f  %2.1f)\"",
            label_set, LABEL_COORDS_SPHERE_RT,
            DEGREES(phi), DEGREES(theta) );
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_COORDS_SPHERE_RT, "(%2.1f, %2.1f)",
                             DEGREES(phi), DEGREES(theta));
  }

  sprintf(command, "UpdateLabel %d %d %f",
          label_set, LABEL_CURVATURE, v->curv);
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_CURVATURE, "%f", v->curv);

  sprintf(command, "UpdateLabel %d %d %f",
          label_set, LABEL_FIELDSIGN, v->fieldsign);
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_FIELDSIGN, "%f", v->fieldsign);

  /* overlay labels. draw the current one in stars. */
  for (field = 0; field < NUM_SCALAR_VALUES; field++ )
  {
    sclv_get_value(v,field,&value);
    if (field == sclv_current_field)
    {
      sprintf(command, "UpdateLabel %d %d \"** %f **\"", label_set,
              LABEL_VAL + field, value);
    }
    else
    {
      sprintf(command, "UpdateLabel %d %d %f", label_set,
              LABEL_VAL + field, value);
    }
    send_tcl_command(command);

    if (LABELSET_CURSOR == label_set)
    {
      /* Build the string for the code here be affixing the field
         number + 1 to the prefix, e.g. field 0 -> !o1 */
      sprintf (command, "%s%d", CPTN_CODE_FIELD_PREFIX, field+1);
      cptn_sprintf_for_code (command, "%f", value);
    }
  }

  sprintf(command, "UpdateLabel %d %d %f", label_set, LABEL_AMPLITUDE,
          hypot(v->val,v->val2));
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_AMPLITUDE, "%f", hypot(v->val,v->val2));

  sprintf(command, "UpdateLabel %d %d %f", label_set, LABEL_ANGLE,
          (float)(atan2(v->val2,v->val)*180.0/M_PI));
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_ANGLE, "%f",
                           (float)(atan2(v->val2,v->val)*180.0/M_PI));

  sprintf(command, "UpdateLabel %d %d %f", label_set, LABEL_DEGREE,
          (float)(atan2(v->val2,v->val)/(2*M_PI)));
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_DEGREE, "%f",
                           (float)(atan2(v->val2,v->val)/(2*M_PI)));

  /* send label update. */
  err = labl_find_label_by_vno( vno, 0, label_index_array,
                                2, &num_labels_found );
  if (err == ERROR_NONE && num_labels_found > 0)
  {
    if (num_labels_found > 1)
    {
      sprintf(command, "UpdateLabel %d %d \"%s, %d others\"",
              label_set, LABEL_LABEL,
              labl_labels[label_index_array[0]].name,
              num_labels_found-1);
      if (LABELSET_CURSOR == label_set)
        cptn_sprintf_for_code (CPTN_CODE_LABEL, "%s, %d others",
                               labl_labels[label_index_array[0]].name,
                               num_labels_found-1 );

    }
    else
    {
      sprintf(command, "UpdateLabel %d %d \"%s\"", label_set, LABEL_LABEL,
              labl_labels[label_index_array[0]].name );
      cptn_sprintf_for_code (CPTN_CODE_LABEL, "%s",
                             labl_labels[label_index_array[0]].name);
    }
    send_tcl_command(command);
  }
  else
  {
    sprintf(command, "UpdateLabel %d %d \"None.\"", label_set, LABEL_LABEL );
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_LABEL, "None");

  }

  //  if (annotationloaded)
  //    {
  r = v->annotation & 0x00ff ;
  g = (v->annotation >> 8) & 0x00ff ;
  b = (v->annotation >> 16) & 0x00ff ;
  sprintf(command, "UpdateLabel %d %d \"%d (%d  %d  %d)\"",
          label_set, LABEL_ANNOTATION, v->annotation, r, g, b);
  send_tcl_command(command);
  if (LABELSET_CURSOR == label_set)
    cptn_sprintf_for_code (CPTN_CODE_ANNOTATION, "%d (%d, %d, %d)",
                           v->annotation, r, g, b);
  //    }

  if (parc_flag && v->val > 0 && parc_names[(int)nint(v->val)])
  {
    sprintf(command, "UpdateLabel %d %d \"%s\"",
            label_set, LABEL_PARCELLATION_NAME,
            parc_names[(int)nint(v->val)]);
    send_tcl_command(command);
    if (LABELSET_CURSOR == label_set)
      cptn_sprintf_for_code (CPTN_CODE_PARCELLATION_NAME, "%s",
                             parc_names[(int)nint(v->val)] );
  }

  /* Although we don't update the time point and condition in the
     label section, we kind of need to update the caption with their
     values here, so do that. */
  if ( sclv_current_field != -1 )
  {
    cptn_sprintf_for_code (CPTN_CODE_TIME_POINT, "%d",
                           sclv_field_info[sclv_current_field].cur_timepoint );

    cptn_sprintf_for_code (CPTN_CODE_CONDITION, "%d",
                           sclv_field_info[sclv_current_field].cur_condition );
  }
}

/* ------------------------------------------------------------- the window */

#ifdef USE_XGLUT_WINDOW

void wndw_create (int x, int y, int width, int height)
{
  xGWin_New (&gWindow, width, height, "tksurfer");
  xGWin_SetEventHandlerFunc (gWindow, wndw_handle_event, gWindow);
  xGWin_ActivateIdleEvents (gWindow);
  xGWin_ActivatePassiveMotionEvents (gWindow);

  openglwindowflag = 1;
}

void wndw_set_title (char* title)
{
  xGWin_SetWindowTitle (gWindow, title);
}

void wndw_handle_event (void* data, xGWin_tEventRef event)
{
  char command[STRLEN];
  int screen_x;
  int screen_y;
  int origin_x;
  int origin_y;
  int size_x;
  int size_y;
  float translate_x;
  float translate_y;
  float md;
  int mvno;
  int vno;
  struct timeval tv;

  switch (event->mType)
  {
  case xGWin_tEventType_KeyDown:
    switch (event->mKey)
    {
    case 'r':
      if (event->mbAltKey && g_interp)
      {
        send_tcl_command("UpdateAndRedraw");
      }
      break;
    case 'q':
      if (event->mbCtrlKey)
      {
        exit (0);
      }
      break;
    case xGWin_tKey_UpArrow:
      if (event->mbAltKey && g_interp)
      {
        send_tcl_command("set gNextTransform(rotate,x) "
                         "[expr $gNextTransform(rotate,x)+18.0]");
        send_tcl_command (command);
      }
      break;
    case xGWin_tKey_DownArrow:
      if (event->mbAltKey && g_interp)
      {
        send_tcl_command("set gNextTransform(rotate,x) "
                         "[expr $gNextTransform(rotate,x)-18.0]");
        send_tcl_command (command);
      }
      break;
    case xGWin_tKey_RightArrow:
      if (event->mbAltKey && g_interp)
      {
        send_tcl_command("set gNextTransform(rotate,y) "
                         "[expr $gNextTransform(rotate,y)-18.0]");
        send_tcl_command (command);
      }
      break;
    case xGWin_tKey_LeftArrow:
      if (event->mbAltKey && g_interp)
      {
        send_tcl_command("set gNextTransform(rotate,y) "
                         "[expr $gNextTransform(rotate,y)+18.0]");
        send_tcl_command (command);
      }
      break;
    case xGWin_tKey_Home:
      sprintf (command,"MoveToolWindow %d %d",
               glutGet (GLUT_WINDOW_X),
               glutGet (GLUT_WINDOW_Y) + w.h + MOTIF_YFUDGE);
      send_tcl_command (command);
      break;
    }
    break; /* case xGWin_tEventType_KeyDown */

  case xGWin_tEventType_MouseDown:
    screen_x = event->mWhere.mnX + w.x;
    screen_y = 1024 - w.y - event->mWhere.mnY ;
    if (1 == event->mButton)
    {
      /* Scale around click */
      if (event->mbCtrlKey)
      {
        origin_x = w.x;
        origin_y = 1024 - w.y - w.h;

        size_x = w.w;
        size_y = w.h;

        translate_x = sf *
                      (screen_x - origin_x - size_x / 2.0 ) *
                      2.0 * fov / size_x;
        translate_y = sf *
                      (screen_y - origin_y - size_y / 2.0 ) *
                      2.0 * fov / size_y;

        translate_brain (-translate_x, -translate_y, 0);
        scale_brain (SCALE_UP_MOUSE);
        redraw();
      }
      /* curvim */
      else if (event->mbShiftKey)
      {
        select_vertex (screen_x, screen_y);
        read_curvim_at_vertex (selection);
        draw_cursor (selection, TRUE);
      }
      else
      {
        /* If something is already selected, hilite it (because
           it's marked now) and deselect it. */
        if (selection>=0)
        {
          draw_vertex_hilite (selection);
          draw_cursor (selection, FALSE);
        }

        /* Select the vertex at this point. If we got one, mark
           it and draw the cursor. */
        select_vertex (screen_x, screen_y);
        if (selection >= 0)
        {
          mark_vertex (selection, TRUE);
          draw_cursor (selection, TRUE);
        }

        /* Draw all the marked verts. */
        draw_marked_vertices ();

        if (showorigcoordsflag)
          find_orig_vertex_coordinates(selection);
      }
    } /* if button 1 */

    /* Button 2, just select and mark this vertex. */
    else if (2 == event->mButton)
    {
      find_closest_marked_vertex ((int)sx, (int)sy, NULL, &vno);
      if (vno>=0)
      {
        fprintf (stderr, "Unmarking %d\n", vno);
        mark_vertex (vno, FALSE);
        draw_cursor (vno, FALSE);
      }
    }

    /* Button 3, if ctrl, zoom out, else clear all the selections. */
    else if (3 == event->mButton)
    {
      if (event->mbCtrlKey)
      {
        origin_x = w.x;
        origin_y = 1024 - w.y - w.h;

        size_x = w.w;
        size_y = w.h;

        translate_x = sf *
          (screen_x - origin_x - size_x / 2.0 ) * 2.0 * fov / size_x;
        translate_y = sf *
          (screen_y - origin_y - size_y / 2.0 ) * 2.0 * fov / size_y;

        translate_brain (-translate_x, -translate_y, 0);
        scale_brain (1.0/SCALE_UP_MOUSE);
        redraw();
      }
      else
      {
        clear_all_vertex_marks();
        labl_select(-1);
        redraw();
      }
    }
    break;

  case xGWin_tEventType_MouseUp:
    break;

  case xGWin_tEventType_MouseMoved:
    screen_x = event->mWhere.mnX + w.x;
    screen_y = 1024 - w.y - event->mWhere.mnY ;
    find_vertex_at_screen_point (screen_x, screen_y, &mvno, &md);
    if (mvno >= 0)
    {
      mouseover_vno = vno;
      update_labels (LABELSET_MOUSEOVER, mouseover_vno, md);
    }
    break;

  case xGWin_tEventType_Resize:
    /* Get the new dimensions. */
    w.x = glutGet( GLUT_WINDOW_X );
    w.y = glutGet( GLUT_WINDOW_Y );
    w.w = event->mWhere.mnX;
    w.h = event->mWhere.mnY;
    curwindowleft = w.x;
    curwindowbottom = w.y + w.h;

    /* Resize the opengl port. */
    glutReshapeWindow (w.w, w.h);

    /* Do some other legacy stuff. */
    resize_window (0);

    /* Move the tool window under us. */
    sprintf (command,"MoveToolWindow %d %d", w.x, w.y + w.h + MOTIF_YFUDGE);
    send_tcl_command (command);

    break;

  case xGWin_tEventType_Draw:
    if (redrawlockflag)
      redraw();
    break;

  case xGWin_tEventType_Idle:
    /* Call the Tk event handling function. */
    while (Tk_DoOneEvent (TK_ALL_EVENTS | TK_DONT_WAIT))
    {}

#ifndef sgi
    tv.tv_sec = 0;
    tv.tv_usec = 10000;
    select(0, NULL, NULL, NULL, &tv);
#else
    sginap((long)1);   /* block for 10 msec */
#endif

    break;

  default:
    break;

  }
}

#endif /* USE_XGLUT_WINDOW */

/* -------------------------------------------------- ctrl-c cancel support */

void cncl_initialize ()
{
  /* init the flags and register our handler. */
  cncl_listening = 0;
  cncl_canceled = 0;
  signal (SIGINT, cncl_handle_sigint);
}

void cncl_start_listening ()
{
  /* set our listening flag. */
  cncl_listening = 1;
}

void cncl_stop_listening ()
{
  /* stop listening and reset the canceled flag. */
  cncl_listening = 0;
  cncl_canceled = 0;
}

int cncl_user_canceled ()
{
  /* just return the canceled flag. */
  return (cncl_canceled);
}

void cncl_handle_sigint (int signal)
{

  /* if we're listening, set the flag, if not, exit normally. */
  if (cncl_listening)
  {
    cncl_canceled = 1;
  }
  else
  {
    printf ("Killed\n");
    fflush ( stdout );
    exit (1);
  }
}

/* ------------------------------------------------------------------------ */

/* ------------------------------------------------------- menu set support */

int
enable_menu_set (int set, int enable)
{
  char tcl_cmd[1024];

  strncpy (tcl_cmd, "tkm_SetEnableGroupStatus", sizeof(tcl_cmd));
  switch (set)
  {
  case MENUSET_VSET_INFLATED_LOADED:
    strcat (tcl_cmd, " mg_InflatedVSetLoaded");
    break;
  case MENUSET_VSET_WHITE_LOADED:
    strcat (tcl_cmd, " mg_WhiteVSetLoaded");
    break;
  case MENUSET_VSET_PIAL_LOADED:
    strcat (tcl_cmd, " mg_PialVSetLoaded");
    break;
  case MENUSET_VSET_ORIGINAL_LOADED:
    strcat (tcl_cmd, " mg_OriginalVSetLoaded");
    break;
  case MENUSET_TIMECOURSE_LOADED:
    strcat (tcl_cmd, " mg_TimeCourseLoaded");
    break;
  case MENUSET_OVERLAY_LOADED:
    strcat (tcl_cmd, " mg_OverlayLoaded");
    break;
  case MENUSET_CURVATURE_LOADED:
    strcat (tcl_cmd, " mg_CurvatureLoaded");
    break;
  case MENUSET_LABEL_LOADED:
    strcat (tcl_cmd, " mg_LabelLoaded");
    break;
  case MENUSET_FIELDSIGN_LOADED:
    strcat (tcl_cmd, " mg_FieldSignLoaded");
    break;
  case MENUSET_FIELDMASK_LOADED:
    strcat (tcl_cmd, " mg_FieldMaskLoaded");
    break;
  default:
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,
                               "enable_menu_set: bad set %d\n",set));
  }
  char enable_cmd[100];
  sprintf (enable_cmd, " %d", enable);
  strcat (tcl_cmd, enable_cmd);
  send_tcl_command(tcl_cmd);

  return(ERROR_NONE);
}

/* ------------------------------------------------------------------------ */

/* -------------------------------------------- multiple vertex set support */

int
vset_initialize()
{
  int i;
  /* set them all to null */
  for (i=0; i<NUM_VERTEX_SETS; i++)
  {
    if (NULL != vset_vertex_list[i])
      free (vset_vertex_list[i]);

    vset_vertex_list[i]=NULL;
  }

  return(NO_ERROR);
}

int vset_read_vertex_set(int set, char* fname)
{
  /* copy the current main verts into tmp */
  MRISsaveVertexPositions( mris, TMP_VERTICES );

  /* read the file */
  if ( MRISreadVertexPositions( mris, fname ) != NO_ERROR )
  {

    /* restore the vertices */
    MRISrestoreVertexPositions( mris, TMP_VERTICES );
    printf("surfer: couldn't load %s!\n", fname);
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "vset_read_vertex_set: MRISreadVertexPositions failed\n"));
  }

  printf("surfer: loaded %s into set %d\n", fname, set);

  /* save the verts into external storage */
  vset_save_surface_vertices(set);

  if(strcmp(fname,"orig") == 0)
    MRISsaveVertexPositions(mris, ORIG_VERTICES) ;
  if(strcmp(fname,"pial") == 0)
    MRISsaveVertexPositions(mris, PIAL_VERTICES) ;
  if(strcmp(fname,"white") == 0)
    MRISsaveVertexPositions(mris, WHITE_VERTICES) ;
  if(strcmp(fname,"inflated") == 0)
    MRISsaveVertexPositions(mris, INFLATED_VERTICES) ;

  /* copy the saved verts back into main space */
  MRISrestoreVertexPositions( mris, TMP_VERTICES );


  /* enable the set menu */
  switch (set)
  {
  case VSET_INFLATED:
    enable_menu_set (MENUSET_VSET_INFLATED_LOADED,1);
    break;
  case VSET_WHITE:
    enable_menu_set (MENUSET_VSET_WHITE_LOADED,1);
    break;
  case VSET_PIAL:
    enable_menu_set (MENUSET_VSET_PIAL_LOADED,1);
    break;
  case VSET_ORIGINAL:
    enable_menu_set (MENUSET_VSET_ORIGINAL_LOADED,1);
    break;
  default:
    break;
  }

  return(NO_ERROR);
}

int
vset_save_surface_vertices(int set)
{

  int vno,nvertices;
  VERTEX *v;

  if (set < 0 || set > NUM_VERTEX_SETS)
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_load_surface_vertices: invalid set %d\n",set));

  /* allocate storage if not done so already */
  if (vset_vertex_list[set]==NULL)
  {
    vset_vertex_list[set]=(VSET_VERTEX*)
                          calloc(mris->nvertices,sizeof(VSET_VERTEX));
    if (vset_vertex_list[set]==NULL)
      ErrorReturn(ERROR_NOMEMORY,
                  (ERROR_NOMEMORY,
                   "vset_save_surface_vertices: allocation of "
                   "vset_vertex_list[%d] failed (size=%d)\n",
                   set,mris->nvertices));
  }

  /* save all vertex values into storage */
  nvertices=mris->nvertices;
  for (vno=0;vno<nvertices;vno++)
  {
    v = &mris->vertices[vno];
    vset_vertex_list[set][vno].x = v->x;
    vset_vertex_list[set][vno].y = v->y;
    vset_vertex_list[set][vno].z = v->z;
  }

  return(NO_ERROR);
}

int
vset_load_surface_vertices(int set)
{

  int vno,nvertices;
  VERTEX* v;

  if (set < 0 || set > NUM_VERTEX_SETS)
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_load_surface_vertices: invalid set %d\n",set));

  /* if not allocated, no verts stored there */
  if (vset_vertex_list[set]==NULL)
  {
    printf("surfer: vertex set not loaded.\n");
    return(NO_ERROR);
  }

  /* load all vertex values from storage */
  nvertices=mris->nvertices;
  for (vno=0;vno<mris->nvertices;vno++)
  {
    v = &mris->vertices[vno];
    v->x = vset_vertex_list[set][vno].x;
    v->y = vset_vertex_list[set][vno].y;
    v->z = vset_vertex_list[set][vno].z;
  }

  return(NO_ERROR);
}

int vset_set_current_set(int set)
{
  int err;

  if (set < 0 || set > NUM_VERTEX_SETS)
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_set_current_set: invalid set %d\n",set));

  if (NULL == vset_vertex_list[set])
  {
    switch (set)
    {
    case VSET_PIAL:
      err = vset_read_vertex_set(set, "pial");
      if (err) ErrorReturn
                 (ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_set_current_set: set %d not loaded\n",set));
      break;
    case VSET_ORIGINAL:
      err = vset_read_vertex_set(set, "orig");
      if (err) ErrorReturn
                 (ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_set_current_set: set %d not loaded\n",set));
      break;
    case VSET_WHITE:
      err = vset_read_vertex_set(set, "white");
      if (err) ErrorReturn
                 (ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_set_current_set: set %d not loaded\n",set));
      break;
    case VSET_INFLATED:
      err = vset_read_vertex_set(set, "inflated");
      if (err) ErrorReturn
                 (ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_set_current_set: set %d not loaded\n",set));
      break;
    }
  }

  /* Make sure this set is loaded unless it's the main set which is
     always loaded.. */
  if (NULL == vset_vertex_list[set] && VSET_MAIN != set)
      ErrorReturn(ERROR_BADPARM,
                  (ERROR_BADPARM,
                   "vset_set_current_set: set %d not loaded\n",set));

  /* save the main set if we haven't done so yet */
  if (vset_current_set==VSET_MAIN && vset_vertex_list[VSET_MAIN]==NULL)
    vset_save_surface_vertices(VSET_MAIN);

  /* load in over the main set */
  vset_load_surface_vertices(set);

  /* save the new vertex set. recompute the normals and mark our
     vertex array dirty so it will rebuild on next redraw. */
  vset_current_set = set;
  vertex_array_dirty = 1;
  compute_normals();

  /* Find the path bounding boxes again. */
  path_update_bounding_boxes();

  /* Find the label bounding boxes again. */
  labl_all_changed();

  return(NO_ERROR);
}

/* ------------------------------------------------------------------------ */

/* -------------------------------------------------- coordinate conversion */

int conv_initialize()
{
  char fname[STRLEN] = "";
  char surf_path[STRLEN] = "";
  FILE* fTest;


  /* allocate our conversion matrices. */
  if (NULL != conv_mnital_to_tal_m_ltz)
    MatrixFree (&conv_mnital_to_tal_m_ltz);

  conv_mnital_to_tal_m_ltz = MatrixIdentity (4, NULL);
  stuff_four_by_four (conv_mnital_to_tal_m_ltz,
                      0.99,       0,     0, 0,
                      0.00,  0.9688, 0.042, 0,
                      0.00, -0.0485, 0.839, 0,
                      0.00,       0,     0, 1);

  if (NULL != conv_mnital_to_tal_m_gtz)
    MatrixFree (&conv_mnital_to_tal_m_gtz);

  conv_mnital_to_tal_m_gtz = MatrixIdentity (4, NULL);
  stuff_four_by_four (conv_mnital_to_tal_m_gtz,
                      0.99,       0,      0, 0,
                      0.00,  0.9688,  0.046, 0,
                      0.00, -0.0485, 0.9189, 0,
                      0.00,       0,      0, 1);

  /* calc the inverses. */
  conv_tal_to_mnital_m_ltz = MatrixInverse (conv_mnital_to_tal_m_ltz, NULL );
  conv_tal_to_mnital_m_gtz = MatrixInverse (conv_mnital_to_tal_m_gtz, NULL );

  /* allocate our temporary matrices. */
  if (NULL != conv_tmp1_m)
    MatrixFree (&conv_tmp1_m);
  conv_tmp1_m = MatrixAlloc (4, 1, MATRIX_REAL);

  if (NULL != conv_tmp2_m)
    MatrixFree (&conv_tmp2_m);
  conv_tmp2_m = MatrixAlloc (4, 1, MATRIX_REAL);


  /* We need to read in the COR- header from the orig volume, if
     available, and get the extract_i_to_r transform. */
  if ( NULL != mris )
  {
    FileNamePath (mris->fname, surf_path);
    sprintf (fname, "%s/../mri/orig/COR-.info", surf_path);

    fTest = fopen (fname, "r" );
    if (NULL != fTest)
    {
      fclose (fTest );
      orig_mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
    }
    if ( NULL == orig_mri_header )
    {
      sprintf (fname, "%s/../mri/orig.mgh", surf_path);
      fTest = fopen (fname, "r" );
      if (NULL != fTest)
      {
        fclose (fTest );
        orig_mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
      }
      if ( NULL == orig_mri_header )
      {
        sprintf (fname, "%s/../mri/orig.mgz", surf_path);
        fTest = fopen (fname, "r" );
        if (NULL != fTest)
        {
          fclose (fTest );
          orig_mri_header = MRIreadHeader (fname, MRI_VOLUME_TYPE_UNKNOWN);
        }
        if ( NULL == orig_mri_header )
        {
          printf ("WARNING: Could not load orig volume.\n"
                  "         Talairach coords will be incorrect.\n" );
        }
      }
    }

    if (NULL != orig_mri_header)
    {
      surfaceRAStoRAS = surfaceRASFromRAS_( orig_mri_header );
    }
  }

  return(NO_ERROR);
}


int conv_ras_to_mnital(float srasx, float srasy, float srasz,
                       float* mnix, float* mniy, float* mniz)
{
  double rasx, rasy, rasz;

  /* If we have the original MRI volume and this surface doesn't have
     the useRealRAS flag, use it to go from surface RAS coords to
     normal RAS coords. Otherwise just use the surface RAS coords. */
  if (NULL != orig_mri_header && !mris->useRealRAS)
  {
    MRIsurfaceRASToRAS (orig_mri_header, srasx, srasy, srasz,
                        &rasx, &rasy, &rasz);
  }
  else
  {
    rasx = srasx;
    rasy = srasy;
    rasz = srasz;
  }

  /* Run the talairach transformation. */
  if (transform_loaded &&
      NULL != lta &&
      lta->type == LINEAR_RAS_TO_RAS)
  {
    LTAworldToWorldEx (lta, rasx, rasy, rasz, mnix, mniy, mniz);
  }

  return(NO_ERROR);
}

int conv_ras_to_tal(float srasx, float srasy, float srasz,
                    float* talx, float* taly, float* talz)
{
  float mnix, mniy, mniz;

  conv_ras_to_mnital (srasx, srasy, srasz, &mnix, &mniy, &mniz);

  *MATRIX_RELT(conv_tmp1_m,1,1) = mnix;
  *MATRIX_RELT(conv_tmp1_m,2,1) = mniy;
  *MATRIX_RELT(conv_tmp1_m,3,1) = mniz;
  *MATRIX_RELT(conv_tmp1_m,4,1) = 1.0;

  if (mniz > 0)
  {
    MatrixMultiply (conv_mnital_to_tal_m_gtz, conv_tmp1_m, conv_tmp2_m);
  }
  else
  {
    MatrixMultiply (conv_mnital_to_tal_m_ltz, conv_tmp1_m, conv_tmp2_m);
  }

  *talx = *MATRIX_RELT(conv_tmp2_m,1,1);
  *taly = *MATRIX_RELT(conv_tmp2_m,2,1);
  *talz = *MATRIX_RELT(conv_tmp2_m,3,1);

  return(NO_ERROR);
}

int conv_mnital_to_ras(float mnix, float mniy, float mniz,
                       float* osrasx, float* osrasy, float* osrasz)
{

  float rasx, rasy, rasz;
  double srasx, srasy, srasz;

  /* Run the talairach transformation. */
  if (transform_loaded &&
      NULL != lta &&
      lta->type == LINEAR_RAS_TO_RAS)
  {
    LTAinverseWorldToWorldEx (lta, mnix, mniy, mniz, &rasx, &rasy, &rasz);
  }

  if (NULL != orig_mri_header && !mris->useRealRAS)
  {
    MRIRASToSurfaceRAS (orig_mri_header, rasx, rasy, rasz,
                        &srasx, &srasy, &srasz);
  }
  else
  {
    srasx = rasx;
    srasy = rasy;
    srasz = rasz;
  }

  *osrasx = (float) srasx;
  *osrasy = (float) srasy;
  *osrasz = (float) srasz;

  return(NO_ERROR);
}

int conv_tal_to_ras(float talx, float taly, float talz,
                    float* srasx, float* srasy, float* srasz)
{

  float mnix, mniy, mniz;


  *MATRIX_RELT(conv_tmp1_m,1,1) = talx;
  *MATRIX_RELT(conv_tmp1_m,2,1) = taly;
  *MATRIX_RELT(conv_tmp1_m,3,1) = talz;
  *MATRIX_RELT(conv_tmp1_m,4,1) = 1.0;

  if (talz > 0)
  {
    MatrixMultiply (conv_tal_to_mnital_m_gtz, conv_tmp1_m, conv_tmp2_m);
  }
  else
  {
    MatrixMultiply (conv_tal_to_mnital_m_ltz, conv_tmp1_m, conv_tmp2_m);
  }

  mnix = *MATRIX_RELT(conv_tmp2_m,1,1);
  mniy = *MATRIX_RELT(conv_tmp2_m,2,1);
  mniz = *MATRIX_RELT(conv_tmp2_m,3,1);

  conv_mnital_to_ras (mnix, mniy, mniz, srasx, srasy, srasz);

  return(NO_ERROR);
}

/* ------------------------------------------------------------------------ */

/* ----------------------------------------------------------- undo support */

int undo_initialize()
{
  int undo_index;
  for (undo_index = UNDO_LIST_POS_FIRST;
       undo_index <= UNDO_LIST_POS_LAST;
       undo_index++ )
    undo_list[undo_index] = NULL;

  undo_send_first_action_name();

  return(NO_ERROR);
}

int undo_get_action_node_size(int action_type)
{
  int size=-1;

  /* return the size of the specific action node impelementation, or
     -1 as an error condition. */
  switch (action_type)
  {
  case UNDO_CUT:
    size = sizeof(UNDO_CUT_NODE);
    break;
  default:
    size = -1;
    break;
  }

  return(size);
}

char* undo_get_action_string(int action_type)
{
  char* string=NULL;

  /* if the action type is in bounds, return a string from the
     undo_action_strings array, else return a null string. */
  if (action_type >= UNDO_INVALID && action_type < NUM_UNDO_ACTIONS)
    string=undo_action_strings[action_type];
  else
    string=NULL;

  return(string);
}

int undo_begin_action(int action_type)
{
  UNDO_ACTION* action = NULL;
  xGArr_tErr array_error = xGArr_tErr_NoErr;
  int action_size = 0;
  int undo_index = 0;

  if (undo_list_state!=UNDO_LIST_STATE_CLOSED)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "undo_begin_action: undo list state is not closed\n"));

  if (action_type<=UNDO_NONE || action_type>=NUM_UNDO_ACTIONS)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "undo_begin_action: invalid action type\n"));

  /* create a new undo action */
  action = (UNDO_ACTION*)malloc(sizeof(UNDO_ACTION));
  if (action==NULL)
    ErrorReturn(ERROR_NO_MEMORY,
                (ERROR_NO_MEMORY,
                 "undo_begin_action: malloc UNDO_ACTION failed\n"));
  action->undo_type = action_type;

  /* get the action node size and allocate the node list */
  action_size = undo_get_action_node_size(action_type);
  if (action_size<=0)
    ErrorReturn(ERROR_UNDO_ACTION,
                (ERROR_UNDO_ACTION,
                 "undo_begin_action: action_size was <= 0\n"));
  array_error = xGArr_New (&(action->node_list), action_size, 64);
  if (array_error!=xGArr_tErr_NoErr)
    ErrorReturn(ERROR_ARRAY,
                (ERROR_ARRAY,
                 "undo_begin_action: xGArr_New failed\n"));

  /* if there was a list in the last slot, delete it */
  if (undo_list[UNDO_LIST_POS_LAST]!=NULL)
  {
    array_error = xGArr_Delete (&(undo_list[UNDO_LIST_POS_LAST]->node_list));
    if (array_error!=xGArr_tErr_NoErr)
      ErrorPrintf(0,"undo_begin_action: xGArr_Delete failed\n");
    free (undo_list[UNDO_LIST_POS_LAST]);
  }

  /* move the undo lists down */
  for (undo_index=UNDO_LIST_POS_LAST;
       undo_index>UNDO_LIST_POS_FIRST;
       undo_index--)
  {
    undo_list[undo_index] = undo_list[undo_index-1];
  }

  /* save the new action in the list */
  undo_list[UNDO_LIST_POS_FIRST] = action;

  /* set our state */
  undo_list_state = UNDO_LIST_STATE_OPEN;

  return(NO_ERROR);
}

int undo_finish_action()
{
  /* set our state */
  undo_list_state = UNDO_LIST_STATE_CLOSED;

  undo_send_first_action_name ();

  return(NO_ERROR);
}

int undo_send_first_action_name()
{
  char* string = NULL;
  char command[1024] = "";

  /* get the undo action string and send it to the interface */
  if (undo_list[UNDO_LIST_POS_FIRST]!=NULL)
  {
    string = undo_get_action_string
             (undo_list[UNDO_LIST_POS_FIRST]->undo_type);
  }
  else
  {
    string = undo_get_action_string (UNDO_NONE);
  }

  if (string==NULL)
    ErrorReturn(ERROR_UNDO_ACTION,
                (ERROR_UNDO_ACTION,
                 "undo_send_first_action_name: "
                 "undo_get_action_string returned null\n"));

  /* build and send the command */
  sprintf (command, "UpdateUndoItemLabel \"%s\"", string);
  send_tcl_command (command);

  return(NO_ERROR);
}

int undo_copy_action_node(void* node)
{
  xGArr_tErr array_error = xGArr_tErr_NoErr;

  /* make sure the list is open */
  if (undo_list_state!=UNDO_LIST_STATE_OPEN)
    ErrorReturn(ERROR_UNDO_LIST_STATE,
                (ERROR_UNDO_LIST_STATE,
                 "undo_copy_action_node: undo list state was not open\n"));

  /* copy the node to the action node list */
  array_error = xGArr_Add (undo_list[UNDO_LIST_POS_FIRST]->node_list, node);
  if (array_error!=xGArr_tErr_NoErr)
    ErrorReturn(ERROR_ARRAY,
                (ERROR_ARRAY,"undo_copy_action_node: xGArr_Add failed\n"));

  return(NO_ERROR);
}

int undo_do_first_action()
{
  int error = NO_ERROR;
  int undo_index = 0;
  xGArr_tErr array_error;

  /* return if we have no list */
  if (undo_list[UNDO_LIST_POS_FIRST]==NULL)
    return(NO_ERROR);

  /* make sure the list is closed */
  if (undo_list_state!=UNDO_LIST_STATE_CLOSED)
    ErrorReturn(ERROR_UNDO_LIST_STATE,
                (ERROR_UNDO_LIST_STATE,
                 "undo_do_first_action: undo list state was not closed\n"));

  /* get the action and switch on its type. call the proper handler. if it
     returns an error, print a msg but don't bail. */
  switch (undo_list[UNDO_LIST_POS_FIRST]->undo_type)
  {
  case UNDO_CUT:
    error = undo_do_action_cut (undo_list[UNDO_LIST_POS_FIRST]);
    break;
  default:
    ErrorReturn(ERROR_UNDO_ACTION,
                (ERROR_UNDO_ACTION,
                 "undo_do_first_action: "
                 "invalid action type in non-null action\n"));
    break;
  }
  if (error!=NO_ERROR)
    ErrorPrintf(ERROR_UNDO_ACTION,
                "undo_do_first_action: undo handler returned an error\n");

  /* delete that action */
  array_error = xGArr_Delete (&(undo_list[UNDO_LIST_POS_FIRST]->node_list) );
  if (array_error!=xGArr_tErr_NoErr)
    ErrorPrintf(ERROR_ARRAY,"undo_do_first_action: xGArr_Delete failed\n");
  free (undo_list[UNDO_LIST_POS_FIRST]);

  /* remove the action from the list and bump everything up. set the
     last pos to null. */
  for (undo_index=UNDO_LIST_POS_FIRST;
       undo_index<UNDO_LIST_POS_LAST;
       undo_index++)
  {
    undo_list[undo_index] = undo_list[undo_index+1];
  }
  undo_list[UNDO_LIST_POS_LAST] = NULL;

  /* refresh the undo menu item with the new first item's name */
  undo_send_first_action_name();

  /* redraw the view */
  vertex_array_dirty = 1;
  rip_faces();
  redraw();

  return(NO_ERROR);
}

int undo_new_action_cut(int cut_type, int index, int rip_value)
{
  UNDO_CUT_NODE cut_action;

  /* make sure the list is open */
  if (undo_list_state!=UNDO_LIST_STATE_OPEN)
    ErrorReturn(ERROR_UNDO_LIST_STATE,
                (ERROR_UNDO_LIST_STATE,
                 "undo_new_action_cut: undo list state was not open\n"));

  /* fill out the fields */
  cut_action.cut_type = cut_type;
  cut_action.index = index;
  cut_action.rip_value = rip_value;

  /* copy it to the list */
  undo_copy_action_node ((void*)&cut_action);

  return(NO_ERROR);
}

int undo_do_action_cut(UNDO_ACTION* action)
{
  UNDO_CUT_NODE cut_action;
  xGArr_tErr array_error;

  if (action==NULL)
    ErrorReturn(ERROR_UNDO_ACTION,
                (ERROR_UNDO_ACTION,
                 "undo_do_action_cut: action was null\n"));

  if (action->undo_type!=UNDO_CUT)
    ErrorReturn(ERROR_UNDO_ACTION,
                (ERROR_UNDO_ACTION,
                 "undo_do_action_cut: action type was not cut\n"));

  /* go thru every node in the array... */
  xGArr_ResetIterator (action->node_list);
  array_error = xGArr_tErr_NoErr;
  while ((array_error=xGArr_NextItem(action->node_list,(void*)&cut_action))
         ==xGArr_tErr_NoErr)
  {
    /* switch on the cut action type and set the face or vertex rip flag to
       the stored value */
    switch (cut_action.cut_type)
    {
    case UNDO_CUT_VERTEX:
      set_vertex_rip (cut_action.index, cut_action.rip_value, FALSE);
      break;
    case UNDO_CUT_FACE:
      set_face_rip (cut_action.index, cut_action.rip_value, FALSE);
      break;
    default:
      ErrorPrintf(ERROR_UNDO_ACTION,
                  "undo_do_action_cut: invalid cut type\n");
    }
  }

  return(NO_ERROR);
}

/* ---------------------------------------------------------------------- */

/* -------------------------------------------- functional volume support */

int func_initialize()
{
  xGArr_tErr array_error = xGArr_tErr_NoErr;

  /* init our list */
  array_error = xGArr_New (&func_selected_ras,
                           sizeof(FUNC_SELECTED_VOXEL), 64);
  if (array_error!=xGArr_tErr_NoErr)
    ErrorReturn(ERROR_ARRAY,
                (ERROR_ARRAY,
                 "func_initialize: xGArr_New failed\n"));

  func_use_timecourse_offset = FALSE;
  func_sub_prestim_avg = FALSE;

  return(ERROR_NONE);
}

int func_load_timecourse (char* fname, FunD_tRegistrationType reg_type,
                          char* registration)
{
  FunD_tErr volume_error;
  char tcl_cmd[1024];
  float time_resolution;
  Volm_tErr volm_err = Volm_tErr_NoErr;
  mriVolumeRef volm = NULL;
  int good = 0;

  if (fname==NULL)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "func_load_timecourse: fname was null\n"));

  /* delete existing mriFunctionalDataRef or MRI. */
  if (func_timecourse!=NULL)
  {
    volume_error = FunD_Delete(&func_timecourse);
    if (volume_error!=FunD_tErr_NoError)
      ErrorPrintf(func_convert_error(volume_error),
                  "func_load_timecourse: error in FunD_Delete\n");
  }

  if (FunD_tRegistration_None == reg_type)
  {
    printf ("surfer: ERROR: Must specify registration type for "
            "time course. Use -timecourse-reg <file>, "
            "-timecourse-reg-find, or -timecourse-reg-identity.\n");
    return ERROR_BADPARM;
  }

  /* Unless they selected the none-needed registration method, we need
     to load up a volume for them to use as the base for the
     transform. So we use the orig header, which should already have
     been loaded, to get one. */
  if (FunD_tRegistration_NoneNeeded != reg_type)
  {
    good = 0;
    if (NULL != orig_mri_header)
    {
      volm_err = Volm_New (&volm);
      if (Volm_tErr_NoErr == volm_err)
      {
        volm_err = Volm_ImportData (volm, orig_mri_header->fname);
        if (Volm_tErr_NoErr == volm_err)
        {
          good = 1;
        }
      }
    }

    if (!good)
    {
      if (NULL != volm)
        Volm_Delete (&volm);
      printf ("surfer: ERROR: You specified a registration type, "
              "but tksurfer cannot find an anatomical volume with which "
              "to calculate the registration. Please make sure the "
              "orig volume is present in the subject's mri/ directory.\n");
      return ERROR_BADPARM;
    }
  }

  if ((FunD_tRegistration_File == reg_type) && 
      (registration) &&
      (strlen(registration)==0))
  {
    printf ("surfer: func_load_timecourse,  ERROR: "
            "missing registration filename\n");
    return ERROR_BADPARM;
  }

  /* create new volume */
  volume_error = FunD_New (&func_timecourse,
                           fname,
                           reg_type,
                           registration,
                           mris->nvertices,
                           volm,
                           mris->hemisphere == LEFT_HEMISPHERE);

  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    printf("### surfer: couldn't load %s\n",fname);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "func_load_timecourse: error in FunD_New\n"));
  }

  /* Notify the volume we're in tkreg space, so our transform is
     correct. */
  volume_error = FunD_ClientSpaceIsTkRegRAS (func_timecourse);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    printf("surfer: couldn't load %s\n",fname);
    ErrorReturn
      (func_convert_error(volume_error),
       (func_convert_error(volume_error),
        "func_load_timecourse: error in FunD_ClientSpaceIsTkRegRAS\n"));
  }

  printf("surfer: loaded timecourse %s\n",fname);

  /* See if it's scalar */
  FunD_IsScalar (func_timecourse, &func_is_scalar_volume);

  if (func_is_scalar_volume)
  {
    printf ("surfer: Interpreting time course volume %s "
            "as encoded scalar volume.\n", fname);
    if ((func_timecourse->mpData->nframes == mris->nvertices) ||
	(func_timecourse->mpData->nframes == 2*mris->nvertices))
    {
      char cmd[STRLEN] ;
      MRI *mri_tmp, *mri_tmp2 ;
      if (mris->hemisphere == LEFT_HEMISPHERE)
      {
	mri_tmp = MRIextractInto(func_timecourse->mpData, NULL, 0, 0, 0, 
				 func_timecourse->mpData->width/2, 
				 func_timecourse->mpData->height, 
				 func_timecourse->mpData->depth,
				 0,0,0) ;
	mri_tmp2 = MRIcopyFrames(mri_tmp, NULL, 0, mris->nvertices-1, 0) ;
      }
      else
      {
	mri_tmp = MRIextractInto(func_timecourse->mpData, NULL, 0, 0, 0, 
				 func_timecourse->mpData->width/2, 
				 func_timecourse->mpData->height, 
				 func_timecourse->mpData->depth,
				 mris->nvertices,0,0) ;
	mri_tmp2 = MRIcopyFrames(mri_tmp, NULL, mris->nvertices, 
				 2*mris->nvertices-1, 0) ;
      }
      MRIfree(&mri_tmp) ;
      MRIfree(&func_timecourse->mpData) ; func_timecourse->mpData = mri_tmp2 ;
      func_timecourse->mNumTimePoints = mri_tmp2->nframes ;
      enable_menu_set (MENUSET_OVERLAY_LOADED, 1);
      enable_menu_set(MENUSET_TIMECOURSE_LOADED, 1) ;
      sprintf (cmd, "UpdateLinkedVarGroup overlay");
      send_tcl_command (cmd);
      linkvertexmode = 1 ;
      fprintf(stderr, "setting linkvertexmode to %d\n", linkvertexmode) ;
    }
  }
  else
  {
    printf ("surfer: Interpreting time course volume %s "
            "as registered functional volume.\n", fname);
  }

  /* get the time res, num conditions, and num presitm points */
  FunD_GetNumPreStimTimePoints (func_timecourse, &func_num_prestim_points);
  FunD_GetTimeResolution (func_timecourse, &time_resolution);
  FunD_GetNumConditions (func_timecourse, &func_num_conditions);
  FunD_GetNumTimePoints (func_timecourse, &func_num_timepoints);
  func_time_resolution = (double)time_resolution;

  /* send the number of conditions */
  sprintf (tcl_cmd, "Graph_SetNumConditions %d", func_num_conditions);
  send_tcl_command (tcl_cmd);

  /* show the graph window */
  send_tcl_command ("Graph_ShowWindow");

  /* enable the related menu items */
  enable_menu_set (MENUSET_TIMECOURSE_LOADED, 1);

  if (NULL != volm)
    Volm_Delete (&volm);

  return(ERROR_NONE);
}

int func_load_timecourse_offset (char* fname, FunD_tRegistrationType reg_type,
                                 char* registration)
{
  FunD_tErr volume_error;
  Volm_tErr volm_err = Volm_tErr_NoErr;
  mriVolumeRef volm = NULL;
  int good = 0;

  if (fname==NULL)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "func_load_timecourse_offset: fname was null\n"));

  /* delete existing */
  if (func_timecourse_offset!=NULL)
  {
    volume_error = FunD_Delete(&func_timecourse_offset);
    if (volume_error!=FunD_tErr_NoError)
      ErrorPrintf(func_convert_error(volume_error),
                  "func_load_timecourse_offset: error in FunD_Delete\n");
  }

  if (FunD_tRegistration_None == reg_type)
  {
    printf ("surfer: ERROR: Must specify registration type for "
            "time course offset. Use -timecourse-offset-reg <file>, "
            "-timecourse-offset-reg-find, or "
            "-timecourse-offset-reg-identity.\n");
    return ERROR_BADPARM;
  }

  /* Unless they selected the none-needed registration method, we need
     to load up a volume for them to use as the base for the
     transform. So we use the orig header, which should already have
     been loaded, to get one. */
  if (FunD_tRegistration_NoneNeeded != reg_type)
  {
    good = 0;
    if (NULL != orig_mri_header)
    {
      volm_err = Volm_New (&volm);
      if (Volm_tErr_NoErr == volm_err)
      {
        volm_err = Volm_ImportData (volm, orig_mri_header->fname);
        if (Volm_tErr_NoErr == volm_err)
        {
          good = 1;
        }
      }
    }

    if (!good)
    {
      if (NULL != volm)
        Volm_Delete (&volm);
      printf ("surfer: ERROR: You specified a registration type, "
              "but tksurfer cannot find an anatomical volume with which "
              "to calculate the registration. Please make sure the "
              "orig volume is present in the subject's mri/ directory.\n");
      return ERROR_BADPARM;
    }
  }

  if ((FunD_tRegistration_File == reg_type) && 
      (registration) &&
      (strlen(registration)==0))
  {
    printf ("surfer: func_load_timecourse_offset,  ERROR: "
            "missing registration filename\n");
    return ERROR_BADPARM;
  }

  /* create new volume */
  volume_error = FunD_New (&func_timecourse_offset,
                           fname,
                           reg_type,
                           registration,
                           mris->nvertices, /* Try to be scalar */
                           volm,
                           mris->hemisphere == LEFT_HEMISPHERE);

  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    ErrorReturn(func_convert_error(volume_error),
                (func_convert_error(volume_error),
                 "func_load_timecourse_offset: error in FunD_New\n"));
  }

  /* Notify the volume we're in tkreg space, so our transform is
     correct. */
  volume_error = FunD_ClientSpaceIsTkRegRAS (func_timecourse_offset);
  if (volume_error!=FunD_tErr_NoError)
  {
    if (NULL != volm)
      Volm_Delete (&volm);
    printf("surfer: couldn't load %s\n",fname);
    ErrorReturn
      (func_convert_error(volume_error),
       (func_convert_error(volume_error),
        "func_load_timecourse_offset: error in FunD_ClientSpaceIsTkRegRAS\n"));
  }

  /* enable offset display */
  func_use_timecourse_offset = TRUE;

  printf("surfer: loaded timecourse offset %s\n",fname);

  /* turn on the offset options */
  send_tcl_command("Graph_ShowOffsetOptions 1");

  if (NULL != volm)
    Volm_Delete (&volm);

  return(ERROR_NONE);
}

int func_select_selected_vertex()
{
  VERTEX* v = NULL;
  char tcl_cmd[1024];

  v = &(mris->vertices[selection]);
  func_select_voxel (selection, v->origx,v->origy,v->origz);

  sprintf (tcl_cmd, "Graph_SetLabel \"VNO %d (%.2f, %.2f, %.2f)\"",
           selection, v->origx, v->origy, v->origz);
  send_tcl_command (tcl_cmd);

  return(ERROR_NONE);
}

int func_select_marked_vertices()
{
  int vno, count;
  VERTEX* v = NULL;
  char tcl_cmd[1024];

  count = 0;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &(mris->vertices[vno]);
    if (v->marked)
    {
      func_select_voxel (vno, v->origx,v->origy,v->origz);
      count++;
    }
  }

  sprintf (tcl_cmd, "Graph_SetLabel \"Averaging %d marked vertices\"", count);
  send_tcl_command (tcl_cmd);

  printf("surfer: averaging %d marked vertices\n",count);
  return(ERROR_NONE);
}

int func_select_label()
{
  int n;
  int count;
  VERTEX* v = NULL;
  LABEL* area;
  char tcl_cmd[1024];

  if (LABL_NONE_SELECTED == labl_selected_label)
  {
    sprintf (tcl_cmd, "Graph_SetLabel \"No label selected\"");
    send_tcl_command (tcl_cmd);
    printf ("surfer: No label selected.\n" );
    return (ERROR_BADPARM);
  }

  count = 0;
  area = labl_labels[labl_selected_label].label ;
  for (n = 0  ; n < area->n_points ; n++)
  {
    if (area->lv[n].vno > 0 && area->lv[n].vno < mris->nvertices)
    {
      v = &mris->vertices[area->lv[n].vno] ;
      func_select_voxel (area->lv[n].vno, v->origx, v->origy, v->origz);
      count++;
    }
  }

  sprintf (tcl_cmd, "Graph_SetLabel \"Averaging label %s\"\n",
           labl_labels[labl_selected_label].name);
  send_tcl_command (tcl_cmd);

  printf("surfer: averaging label %s\"\n",
         labl_labels[labl_selected_label].name);
  return(ERROR_NONE);
}

int func_clear_selection()
{
  xGArr_tErr array_error = xGArr_tErr_NoErr;

  if (func_selected_ras==NULL)
    ErrorReturn(ERROR_NOT_INITED,
                (ERROR_NOT_INITED,
                 "func_clear_selection: func_selected_ras is null\n"));

  /* clear our selection list */
  array_error = xGArr_Clear (func_selected_ras);
  if (array_error!=xGArr_tErr_NoErr)
    ErrorReturn(ERROR_ARRAY,
                (ERROR_ARRAY,
                 "func_clear_selection: xGArr_Clear failed\n"));

  return(ERROR_NONE);
}

int func_select_voxel (int vno, float x, float y, float z)
{
  xGArr_tErr array_error = xGArr_tErr_NoErr;
  FUNC_SELECTED_VOXEL voxel;

  if (func_selected_ras==NULL)
    ErrorReturn(ERROR_NOT_INITED,
                (ERROR_NOT_INITED,
                 "func_select_voxel: func_selected_ras is null\n"));

  /* build a voxel */
  voxel.vno = vno;
  voxel.x = x;
  voxel.y = y;
  voxel.z = z;

  /* add this voxel to the list */
  array_error = xGArr_Add (func_selected_ras, &voxel);
  if (array_error!=xGArr_tErr_NoErr)
    ErrorReturn(ERROR_ARRAY,
                (ERROR_ARRAY,
                 "func_select_voxel: xGArr_Add failed\n"));

  return(ERROR_NONE);
}

int func_graph_timecourse_selection ()
{

  int cond;
  int tp;
  int num_good_voxels;
  float* values;
  float* deviations;
  char* tcl_cmd;
  int tcl_cmd_size;
  FunD_tErr func_error = FunD_tErr_NoError;
  float second;

  /* make sure we have a volume */
  if (func_timecourse==NULL)
  {
    printf ("surfer: time course not loaded\n");
    return (ERROR_NONE);
  }

  /* make sure the graph window is open. */
  send_tcl_command ("Graph_ShowWindow");

  /* allocate storage arrays. */
  values = calloc (func_num_timepoints,sizeof(float));
  if (values==NULL)
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "func_graph_timecourse_selection: "
                 "calloc(%d,float) failed for values\n",
                 func_num_timepoints));
  deviations = calloc (func_num_timepoints,sizeof(float));
  if (deviations==NULL)
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "func_graph_timecourse_selection: "
                 "calloc(%d,float) failed for deviations\n",
                 func_num_timepoints));

  /* allocate the argument string.  */
  tcl_cmd_size = (func_num_timepoints * knLengthOfGraphDataItem) +
                 knLengthOfGraphDataHeader;
  tcl_cmd = (char*)malloc(sizeof(char)*tcl_cmd_size);
  if (tcl_cmd==NULL )
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "func_graph_timecourse_selection: "
                 "failed to alloc %d char string\n",tcl_cmd_size));

  /* clear the graph so if nothing draws, we won't have any old stuff there. */
  send_tcl_command ("Graph_ClearGraph");
  send_tcl_command ("Graph_BeginData");

  /* for each condition... */
  for (cond=0;cond<func_num_conditions;cond++)
  {

    /* get the average values for this condition. we also find out from
       this function how many of our selected voxels were actually in
       functional space. */
    func_calc_avg_timecourse_values (cond, &num_good_voxels,
                                     values, deviations);

    /* if we had any good voxels, build a list of values and send to graph */
    if (num_good_voxels>0)
    {
      /* write the cmd name, condition number and first brace */
      sprintf (tcl_cmd, "Graph_SetPointsData %d {", cond);

      /* for each time point... */
      for (tp=0; tp < func_num_timepoints; tp++)
      {

        /* convert to a second. If this is a
           mriFunctionalDataRef, let it convert for us, using TR
           if present, otherwise our second is just our time
           point. */
        if (func_timecourse)
        {
          func_error =
            FunD_ConvertTimePointToSecond
            (func_timecourse, tp, &second);
          if (func_error!=FunD_tErr_NoError)
            ErrorPrintf(func_convert_error(func_error),
                        "func_graph_timecourse_selection: "
                        "error in FunD_ConvertTimePointToSecond "
                        "tp=%d\n",tp);
        }

        /* write the second and value to the arg list */
        sprintf (tcl_cmd, "%s %1.1f %2.5f", tcl_cmd, second, values[tp]);
      }

      /* write the last brace and send the cmd. */
      sprintf (tcl_cmd, "%s}", tcl_cmd);
      send_tcl_command (tcl_cmd);

      /* send the error bars. write the cmd name. */
      sprintf (tcl_cmd, "Graph_SetErrorData %d {", cond);

      /* for each time point, write the deviation value. */
      for (tp=0; tp < func_num_timepoints; tp++)
        sprintf (tcl_cmd, "%s %2.5f", tcl_cmd, deviations[tp]);

      /* write the last brace and send the cmd. */
      sprintf (tcl_cmd, "%s}", tcl_cmd);
      send_tcl_command (tcl_cmd);
    }
  }

  send_tcl_command ("Graph_EndData");

  return(ERROR_NONE);
}

int func_print_timecourse_selection (char* fname)
{
  FILE* fp = NULL;
  int cond;
  int tp;
  int num_good_voxels;
  float* values = NULL;
  float* deviations;

  if (fname==NULL)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "func_print_timecourse_selection: file name was null\n"));

  /* make sure we have a volume */
  if (func_timecourse==NULL)
    ErrorReturn(ERROR_NOT_INITED,
                (ERROR_NOT_INITED,
                 "func_print_timecourse_selection: No timecourse volume.\n"));

  /* allocate storage arrays. */
  /* allocate storage arrays. */
  values = calloc (func_num_timepoints,sizeof(float));
  if (values==NULL)
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "func_print_timecourse_selection: "
                 "calloc(%d,float) failed for values\n",
                 func_num_timepoints));

  deviations = calloc (func_num_timepoints,sizeof(float));
  if (deviations==NULL)
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "func_print_timecourse_selection: "
                 "calloc(%d,float) failed for deviations\n",
                 func_num_timepoints));

  fp = fopen (fname, "w");
  if (fp==NULL)
    ErrorReturn(ERROR_NOFILE,
                (ERROR_NOFILE,
                 "func_print_timecourse_selection: "
                 "file %s couldn't be opened\n", fname));

  /* get the num of conditions. for each one... */
  FunD_GetNumConditions (func_timecourse,&func_num_conditions);
  for (cond=0;cond<func_num_conditions;cond++)
  {
    fprintf (fp,"Condition %d/%d\n",cond,func_num_conditions);

    /* get the average values for this condition. we also find out from
       this function how many of our selected voxels were actually in
       functional space. */
    func_calc_avg_timecourse_values (cond, &num_good_voxels,
                                     values, deviations);

    fprintf (fp,"%d voxels in range.\n",num_good_voxels);

    /* if we had any good voxels, print out a summary */
    if (num_good_voxels>0)
      for (tp=0;tp<func_num_timepoints;tp++)
        fprintf (fp, "%d\t%d\t%f\t%f\n",
                 tp, func_num_timepoints, values[tp], deviations[tp]);

  }

  fclose (fp);

  if (NULL != values)
    free (values);
  if (NULL != deviations)
    free (deviations);

  return(ERROR_NONE);
}

int func_calc_avg_timecourse_values (int condition, int* num_good_voxels,
                                     float values[], float deviations[] )
{

  int tp = 0;
  float* sums;
  tBoolean present;
  FUNC_SELECTED_VOXEL selected_voxel;
  xVoxel voxel;
  xGArr_tErr array_error = xGArr_tErr_NoErr;
  FunD_tErr func_error = FunD_tErr_NoError;
  float offset,offset_sum;

  /* allocate sums */
  sums = (float*) calloc (func_num_timepoints,sizeof(float));
  if (sums==NULL)
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "func_calc_avg_timecourse_values: "
                 "calloc(%d,float) failed\n",func_num_timepoints));

  /* no good values yet */
  (*num_good_voxels) = 0;

  /* get the values at all time points. add them
     in the sums array. keep track of how many we have. */
  offset_sum = 0;
  xGArr_ResetIterator (func_selected_ras);
  array_error = xGArr_tErr_NoErr;
  while ((array_error=xGArr_NextItem(func_selected_ras,(void*)&selected_voxel))
         ==xGArr_tErr_NoErr)
  {

    /* If it's scalar, our index is vno,0,0, otherwise use the voxel. */
    if (func_is_scalar_volume)
    {
      xVoxl_Set (&voxel, selected_voxel.vno, 0, 0);
    }
    else
    {
      xVoxl_SetFloat (&voxel,
                      selected_voxel.x,selected_voxel.y,selected_voxel.z);
    }

    /* get all values at this voxel.*/
    if (func_timecourse)
    {
      func_error = FunD_GetDataForAllTimePoints( func_timecourse, &voxel,
                   condition, values );
      /* if it was out of bounds, continue. */
      if (func_error!=FunD_tErr_NoError)
      {
        continue;
      }
    }

    /* if we are displaying offsets and we have offset data... */
    if (func_use_timecourse_offset && func_timecourse_offset)
    {
      /* get the offset at this value. only one plane in offset data. */
      func_error = FunD_GetData( func_timecourse_offset,
                                 &voxel, 0, 0, &offset );
      if (func_error==FunD_tErr_NoError )
      {
        /* divide all functional values by the offset and
           mult by 100 to get a percent */
        for (tp = 0; tp < func_num_timepoints; tp++)
          values[tp] = (values[tp]/offset)*100.0;

        /* get a sum off the offset. */
        offset_sum += offset;
      }
    }

    /* add all values to our sums array */
    for (tp=0; tp < func_num_timepoints; tp++)
      sums[tp] += values[tp];

    /* inc our count */
    (*num_good_voxels)++;
  }


  /* if we don't have any values at this point, our whole selections
     is out of range. */
  if ((*num_good_voxels)==0)
    return(ERROR_NONE);

  /* divide everything by the number of values to find the average */
  for (tp=0; tp < func_num_timepoints; tp++)
    values[tp] = sums[tp] / (float)(*num_good_voxels);

  /* if we have offset values, divide the offset sum by the number of
     values. */
  if (func_use_timecourse_offset)
    offset = offset_sum / (float)(*num_good_voxels);

  /* if there is error data present.. */
  if (func_timecourse)
  {
    FunD_IsErrorDataPresent(func_timecourse, &present);
    if (present)
    {
      /* go through the voxel list again. */
      xGArr_ResetIterator (func_selected_ras);
      array_error = xGArr_tErr_NoErr;
      while ((array_error=xGArr_NextItem(func_selected_ras,
                                         (void*)&selected_voxel))
             ==xGArr_tErr_NoErr)
      {

        /* If it's scalar, our index is vno,0,0, otherwise use
           the voxel. */
        if (func_is_scalar_volume)
        {
          xVoxl_Set (&voxel, selected_voxel.vno, 0, 0);
        }
        else
        {
          xVoxl_SetFloat (&voxel,
                          selected_voxel.x,
                          selected_voxel.y,
                          selected_voxel.z);
        }


        /* get the deviations at all time points */
        func_error =
          FunD_GetDeviationForAllTimePoints (func_timecourse, &voxel,
                                             condition, deviations);

        if (func_error!=FunD_tErr_NoError)
          ErrorPrintf(func_convert_error(func_error),
                      "func_calc_avg_timecourse_values: "
                      "error in FunD_GetDeviationForAllTimePoints\n");

        /* if we have offset values... */
        if (func_use_timecourse_offset && func_timecourse_offset!=NULL)
        {

          /* divide all deviations by the offset and mult by 100 to
             get a percent */
          for (tp=0; tp < func_num_timepoints; tp++)
            deviations[tp] = (deviations[tp]/offset)*100.0;
        }
      }
    }
  }
  else
  {
    /* fill deviations with 0s */
    for (tp=0; tp < func_num_timepoints; tp++)
      deviations[tp] = 0;
  }

  return(ERROR_NONE);
}

int func_calc_correlation_and_write_to_overlay (int selection_vno, int field)
{

  char   label[1024];
  int    vno;
  VERTEX *v;
  int    tp, cond;
  float  val1, val2;
  double cor, min, max, num, d1, d2;
  xVoxel selection_voxel, cur_voxel;
  int    nframes;
  FunD_tErr func_error = FunD_tErr_NoError;

  /* make sure we have a volume */
  if (func_timecourse==NULL)
    ErrorReturn(ERROR_NOT_INITED,
                (ERROR_NOT_INITED,
                 "func_print_timecourse_selection: No timecourse volume.\n"));

  if (selection_vno < 0 || selection_vno >= mris->nvertices)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "surfer: Please select a vertex.\n"));

  /* Initialize our field. */
  sprintf (label, "Corrltn %d", selection_vno);
  sclv_new_empty (field, label);

  /* Build the voxel for the selected point. If it's scalar, our index
     is vno,0,0, otherwise use the orig coords. */
  if (func_is_scalar_volume)
    xVoxl_Set (&selection_voxel, selection_vno, 0, 0);
  else
  {
    v = &mris->vertices[selection_vno];
    xVoxl_SetFloat (&selection_voxel, v->origx, v->origy, v->origz);
  }

  /* Calc the number of frames. */
  nframes = (func_num_conditions * func_num_timepoints);

  /* For every vertex... */
  min = 1000000000;
  max = -min;
  for (vno = 0; vno < mris->nvertices; vno++)
  {
    v = &mris->vertices[vno];

    /* If it's scalar, our index is vno,0,0, otherwise use the
       orig coords. */
    if (func_is_scalar_volume)
      xVoxl_Set (&cur_voxel, vno, 0, 0);
    else
      xVoxl_SetFloat (&cur_voxel, v->origx, v->origy, v->origz);

    num = d1 = d2 = cor = 0.0;
    for (cond=0; cond < func_num_conditions; cond++)
    {
      for (tp=0; tp < func_num_timepoints; tp++)
      {

        /* Get the value at the selection point and the current
           point. This will return an error if it's a spatial
           volume and these coords are out of bounds. */
        func_error = FunD_GetData( func_timecourse,
                                   &selection_voxel, cond, tp, &val1 );
        if (func_error!=FunD_tErr_NoError )
        {
          continue;
        }
        func_error = FunD_GetData( func_timecourse,
                                   &cur_voxel, cond, tp, &val2 );
        if (func_error!=FunD_tErr_NoError )
        {
          continue;
        }

        /* Math stuff */
        cor += (val1-val2)*(val1-val2);
        num += val1*val2;
        d1 += (val1*val1);
        d2 += (val2*val2);
      }
    }

    /* Calculation the correlation. */
    cor = 1 - (sqrt(cor/nframes));
    cor = (num/nframes) / sqrt((d1/nframes)*(d2/nframes));

    /* Set the value in the field. */
    sclv_set_value(v, field, cor);

    /* Update the min and max. */
    if (cor < min)
      min = cor;
    if (cor > max)
      max = cor;
  }

  /* Set the min and max. */
  sclv_value_min = sclv_field_info[field].min_value = min;
  sclv_value_max = sclv_field_info[field].max_value = max;

  /* Set the threshold so we go from -1 -> 1 */
  fthresh = sclv_field_info[field].fthresh = 0;
  fmid    = sclv_field_info[field].fmid = 0.5;
  fslope  = sclv_field_info[field].fslope = 1;

  /* Calc the frquencies */
  sclv_calc_frequencies (field);

  /* Request a redraw. Turn on the overlay flag and select this value
     set. */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (field);

  return (ERROR_NONE);
}

int func_normalize ()
{
  FunD_tErr func_error = FunD_tErr_NoError;

  /* make sure we have a volume */
  if (func_timecourse==NULL)
    ErrorReturn(ERROR_NOT_INITED,
                (ERROR_NOT_INITED,
                 "func_normalize: No timecourse volume.\n"));

  func_error = FunD_NormalizeOverAll( func_timecourse );
  if (func_error!=FunD_tErr_NoError)
    ErrorPrintf(func_convert_error(func_error),
                "func_normalize: error in FunD_NormalizeOverAll");

  return(ERROR_NONE);
}


int func_convert_error (FunD_tErr volume_error)
{
  int error = ERROR_NONE;
  switch (volume_error)
  {
  case FunD_tErr_PathNotFound:
  case FunD_tErr_CouldntGuessStem:
  case FunD_tErr_DataNotFound:
  case FunD_tErr_HeaderNotFound:
  case FunD_tErr_UnrecognizedHeaderFormat:
  case FunD_tErr_QuestionableHeaderFormat:
  case FunD_tErr_CouldntDetermineDataType:
    error = ERROR_NOFILE;
    break;
  case FunD_tErr_CouldntAllocateVolume:
  case FunD_tErr_CouldntAllocateStorage:
  case FunD_tErr_CouldntAllocateMatrix:
    error = ERROR_NO_MEMORY;
    break;
  default:
    error = ERROR_FUNC;
    break;
  }
  return(error);
}
/* ---------------------------------------------------------------------- */

/* --------------------------------------------------- scalar value mgmnt */

int sclv_initialize ()
{
  int field;

  /* no layers loaded, clear all the stuff. */
  for (field = 0; field < NUM_SCALAR_VALUES; field++)
  {
    sclv_field_info[field].is_functional_volume = FALSE;
    sclv_field_info[field].is_scalar_volume = FALSE;
    sclv_field_info[field].cur_timepoint = -1;
    sclv_field_info[field].cur_condition = -1;
    sclv_field_info[field].func_volume = NULL;
    sclv_field_info[field].fthresh = fthresh;
    sclv_field_info[field].fmid = fmid;
    sclv_field_info[field].foffset = foffset;
    sclv_field_info[field].fslope = fslope;
  }

  return (ERROR_NONE);
}

int sclv_unload_field (int field)
{

  FunD_tErr volume_error;

  if (field < 0 || field > NUM_SCALAR_VALUES)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_unload_field: field was out of bounds: %d)",field));

  /* if we have a binary volume, delete it */
  if (sclv_field_info[field].is_functional_volume &&
      sclv_field_info[field].func_volume != NULL )
  {
    volume_error = FunD_Delete (&(sclv_field_info[field].func_volume));
    if (volume_error!=FunD_tErr_NoError)
      ErrorReturn(func_convert_error(volume_error),
                  (func_convert_error(volume_error),
                   "sclv_unload_field: error in FunD_Delete\n"));

    sclv_field_info[field].is_functional_volume = FALSE;

    /* force a repaint next load */
    sclv_field_info[field].cur_timepoint = -1;
    sclv_field_info[field].cur_condition = -1;
  }

  return (ERROR_NONE);
}

int
sclv_calc_frequencies(int field)
{
  int timepoint, condition;
  float* values;
  int bin;
  float num_values;
  int vno;
  float valPerBin;

  if (field < 0 || field > NUM_SCALAR_VALUES)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_unload_field: field was out of bounds: %d)",field));

  sclv_field_info[field].num_freq_bins = SCLV_NUM_FREQUENCY_BINS;

  /* get the value range. */
  num_values = (sclv_field_info[field].max_value -
                sclv_field_info[field].min_value);
  valPerBin = num_values / (float)sclv_field_info[field].num_freq_bins;

  /* allocate storage for each time point and condition... */
  if (NULL != sclv_field_info[field].frequencies)
    free (sclv_field_info[field].frequencies);

  sclv_field_info[field].frequencies =
    calloc( sclv_field_info[field].num_conditions, sizeof(int**) );

  values = (float*) calloc( mris->nvertices, sizeof(float) );
  if (values == NULL)
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,
                 "sclv_calc_frequencies: couldn't allocate values storage."));

  /* Calc the zero bin and save it. */
  sclv_field_info[field].zero_bin_index =
    (float)(0 - sclv_field_info[field].min_value) / (float)valPerBin;
  sclv_field_info[field].num_zeroes_in_zero_bin = 0;

  for (condition = 0;
       condition < sclv_field_info[field].num_conditions; condition++)
  {
    sclv_field_info[field].frequencies[condition] =
      calloc( sclv_field_info[field].num_timepoints, sizeof(int*) );

    for (timepoint = 0;
         timepoint < sclv_field_info[field].num_timepoints; timepoint++)
    {

      /* allocate an array of num_freq_bins ints */
      sclv_field_info[field].frequencies[condition][timepoint] =
        calloc( sclv_field_info[field].num_freq_bins, sizeof(int) );

      /* Get the values at this timepoint and condition. */
      sclv_get_values_for_field_and_timepoint (field, timepoint, condition,
          values);

      /* for each vno, find the bin the value should go in and inc
         the count in that bin. */
      for (vno = 0 ; vno < mris->nvertices ; vno++)
      {
        if (values[vno] != 0)
        {
          bin =
            (float)(values[vno] -
                    sclv_field_info[field].min_value) \
            / (float)valPerBin;
          if (bin >= 0 && bin < sclv_field_info[field].num_freq_bins)
            sclv_field_info[field].\
            frequencies[condition][timepoint][bin]++;
        }
        else
        {
          /* This is a zero so inc the number of zeroes. */
          sclv_field_info[field].num_zeroes_in_zero_bin++;
        }
      }
    }
  }

  free (values);

  return (ERROR_NONE);
}

int sclv_set_overlay_alpha (double alpha)
{
  if (alpha < 0 || alpha > 1.0)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_set_overlay_alpha: alpha must be 0-1"));

  sclv_overlay_alpha = alpha;

  return (ERROR_NONE);
}

int sclv_set_current_field (int field)
{
  if (field < 0 || field > NUM_SCALAR_VALUES)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_set_current_field: field was out of bounds: %d)",
                 field));

  if (gd_fnames[field])
  {
    char tcl_cmd[STRLEN];
    sprintf(tcl_cmd, "GDF_Load %s\n", gd_fnames[field]) ;
    send_tcl_command (tcl_cmd);
//	gdfs[i] = gdfRead(fsgd_fnames[i], 1) ;
  }

  /* save the current threshold */
  sclv_field_info[sclv_current_field].fthresh = fthresh;
  sclv_field_info[sclv_current_field].fmid = fmid;
  ;
  sclv_field_info[sclv_current_field].foffset = foffset;
  ;
  sclv_field_info[sclv_current_field].fslope = fslope;
  sclv_field_info[sclv_current_field].cur_condition = sclv_cur_condition;
  sclv_field_info[sclv_current_field].cur_timepoint = sclv_cur_timepoint;

  /* update the field */
  sclv_current_field = field;

  /* set the shared vars. */
  fthresh = sclv_field_info[sclv_current_field].fthresh;
  fmid = sclv_field_info[sclv_current_field].fmid;
  ;
  foffset = sclv_field_info[sclv_current_field].foffset;
  ;
  fslope = sclv_field_info[sclv_current_field].fslope;
  sclv_cur_condition = sclv_field_info[sclv_current_field].cur_condition;
  sclv_cur_timepoint = sclv_field_info[sclv_current_field].cur_timepoint;
  sclv_num_timepoints = sclv_field_info[sclv_current_field].num_timepoints;
  sclv_num_conditions = sclv_field_info[sclv_current_field].num_conditions;
  sclv_value_min = sclv_field_info[sclv_current_field].min_value;
  sclv_value_max = sclv_field_info[sclv_current_field].max_value;

  /* send the info for this field */
  sclv_send_current_field_info ();

  return (ERROR_NONE);
}

int sclv_send_current_field_info ()
{

  char cmd[1024];

  /*  printf("sending info for field=%d\n\tmin=%f
      max=%f\n\tfthresh=%f fmid=%f fslope=%f\n\ttp=%d cn=%d
      ntps+%d ncns=%d\n", sclv_current_field, sclv_value_min,
      sclv_value_max, fthresh, fmid, fslope, sclv_cur_timepoint,
      sclv_cur_condition, sclv_num_timepoints, sclv_num_conditions); */

  /* Send the current histogram info here as well. */
  sclv_send_histogram (sclv_current_field);

  sprintf(cmd, "set gaLinkedVar(fslope) %f", fslope);
  sprintf(cmd, "set gaLinkedVar(fmid) %f", fmid);
  sprintf(cmd, "set gaLinkedVar(fthreshmax) %f", fthreshmax);
  sprintf(cmd, "set gaLinkedVar(fthresh) %f", fthresh);
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);
  sprintf (cmd, "OverlayLayerChanged");
  send_tcl_command (cmd);

  return (ERROR_NONE);
}

int sclv_set_timepoint_of_field (int field,
                                 int timepoint, int condition)
{
  int vno;
  VERTEX* v;
  float* values;

  if (field < 0 || field > NUM_SCALAR_VALUES)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_set_timepoint_of_field: field was out of bounds: %d)",
                 field));

  /* check if this field has a binary volume */
  if (sclv_field_info[field].is_functional_volume == FALSE)
  {
    /* commented out because tksurfer.tcl will call this function
       even when there is a .w file in this layer. */
    /*      ErrorReturn(ERROR_BADPARM,
            (ERROR_BADPARM,
            "sclv_set_timepoint_of_field: "
            "field %d doesn't have a binary volume",field)); */
    return (ERROR_NONE);
  }

  /* make sure it actually has one */
  if ((sclv_field_info[field].is_functional_volume &&
       sclv_field_info[field].func_volume == NULL))
  {
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_set_timepoint_of_field: "
                 "field %d thinks it has a binary volume "
                 "but doesn't really",field));
  }

  if (timepoint < 0 || timepoint > sclv_field_info[field].num_timepoints)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_set_timepoint_of_field: "
                 "timepoint was out of bounds: %d (max %d)",
                 timepoint, sclv_field_info[field].num_timepoints));
  if (condition < 0 || condition > sclv_field_info[field].num_conditions)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_set_timepoint_of_field: "
                 "condition was out of bounds: %d",condition));

  /* check the timepoint and condition. if they're not what we're already
     using...*/
  if (timepoint != sclv_field_info[field].cur_timepoint ||
      condition != sclv_field_info[field].cur_condition )
  {
    values = (float*) calloc (mris->nvertices, sizeof(float));
    if (values == NULL)
      ErrorReturn(ERROR_NOMEMORY,
                  (ERROR_NOMEMORY,
                   "sclv_set_timepoint_of_field: "
                   "couldn't allocate values storage."));

    /* Get the values here. */
    sclv_get_values_for_field_and_timepoint (field, timepoint, condition,
        values);

    /* For each vertex, set the value.. */
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      v = &mris->vertices[vno] ;
      sclv_set_value (v, field, values[vno]);
    }

    /* save the timepoint and condition. */
    sclv_field_info[field].cur_timepoint = timepoint;
    sclv_field_info[field].cur_condition = condition;

    /* send the info for the current field */
    if (field == sclv_current_field)
      sclv_send_current_field_info();

    free (values);
  }

  /* This is kind of a hack to update the caption in case it has time
     point and condition codes in it. Otherwise it won't get updated
     until we click a vertex. */
  update_labels (LABELSET_CURSOR, selection, 0);

  return (ERROR_NONE);
}

int sclv_get_values_for_field_and_timepoint (int field, int timepoint,
    int condition, float* values)
{
  int vno;
  VERTEX* v;
  xVoxel voxel;
  FunD_tErr volume_error;
  float func_value;

  if (field < 0 || field > NUM_SCALAR_VALUES)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_get_values_for_field_and_timepoint: "
                 "field was out of bounds: %d)",field));

  /* Check timepoint and condition */
  if (timepoint < 0 || timepoint > sclv_field_info[field].num_timepoints)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_get_values_for_field_and_timepoint: "
                 "timepoint was out of bounds: %d",timepoint));
  if (condition < 0 || condition > sclv_field_info[field].num_conditions)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_get_values_for_field_and_timepoint: "
                 "condition was out of bounds: %d",condition));

  /* Make sure we have some output. */
  if (values == NULL )
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_get_values_for_field_and_timepoint: "
                 "values was NULL"));

  DisableDebuggingOutput;

  /* for each vertex, grab a value out of the volume and stick it
     in the field */
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;

    /* skip ripped verts */
    if (v->ripflag)
    {
      values[vno] = 0;
      continue;
    }

    /* Get value from the right volume. */
    if (sclv_field_info[field].is_functional_volume &&
        sclv_field_info[field].func_volume )
    {

      /* If it's scalar, our index is vno,0,0, otherwise use the
         orig coords. */
      if (sclv_field_info[field].is_scalar_volume)
      {
        xVoxl_Set (&voxel, vno, 0, 0);
      }
      else
      {
        xVoxl_SetFloat (&voxel, v->origx, v->origy, v->origz);
      }

      /* If the voxel is valid here, use the value, else set
         to 0. */
      volume_error = FunD_GetData(sclv_field_info[field].func_volume,
                                  &voxel, condition, timepoint,
                                  &func_value);
      if (volume_error == FunD_tErr_NoError)
      {
        values[vno] = func_value;
      }
      else
      {
        values[vno] = 0;
      }
    }
    else
    {
      /* There is no volume so we can just copy the proper sclv
         field value into our output. */
      sclv_get_value (v, field, &values[vno] );
    }
  }
  EnableDebuggingOutput;

  return (ERROR_NONE);
}


int sclv_set_current_threshold_from_percentile (float thresh, float mid,
    float max)
{
  HISTOGRAM *h ;
  int       bthresh, bmid, bmax ;
  float     thresh_value, mid_value, max_value, slope ;

  h = MRISgetHistogram(mris, 1000, sclv_current_field);
  if (ignorezeroesinhistogramflag)
    HISTOclearZeroBin(h) ;
  HISTOmakeCDF(h, h) ;
  if (Gdiag & DIAG_WRITE && DIAG_VERBOSE_ON)
    HISTOplot(h, "h.plt") ;
  bthresh = HISTOfindBinWithCount(h, thresh) ;
  bmid = HISTOfindBinWithCount(h, mid) ;
  bmax = HISTOfindBinWithCount(h, max) ;
  thresh_value = h->bins[bthresh] ; mid_value = h->bins[bmid] ; max_value = h->bins[bmax] ;
  if (FEQUAL(max_value, thresh_value))
    slope = 1.0 ;
  else
    slope = .5 / (max_value - mid_value) ;
  if (thresh_value < 0)
  {
    sclv_field_info[sclv_current_field].foffset = foffset = thresh_value - (mid_value-thresh_value) ;
    printf("setting foffset to %f\n", foffset) ;
    thresh_value -= foffset ;
    mid_value -= foffset ;
    max_value -= foffset ;
    printf("HISTO vals should be %f, %f, %f (slope = .5/(mx-md) =  %2.2f)\n", 
         thresh_value, mid_value, max_value, slope) ;
  }
  else
  {

    sclv_field_info[sclv_current_field].foffset = foffset = 0;
    printf("HISTO vals should be %f, %f, %f (slope = .5/(mx-md) =  %2.2f)\n", 
           thresh_value, mid_value, max_value, slope) ;
  }

    
  HISTOfree(&h) ;
  fthresh = sclv_field_info[sclv_current_field].fthresh = thresh_value;
  fmid = sclv_field_info[sclv_current_field].fmid = mid_value;
  fslope = sclv_field_info[sclv_current_field].fslope = slope ;
  fthreshmax = sclv_field_info[sclv_current_field].max_value = max_value ;
  vertex_array_dirty = 1;
#if 0
  sclv_set_threshold_from_percentile(sclv_current_field,
                                     thresh, mid, max);
#endif

  sclv_send_current_field_info();
  redraw() ;

  return (ERROR_NONE);
}

int sclv_set_threshold_from_percentile (int field, float thresh, float mid,
                                        float max)
{
  float thresh_value=0, mid_value=0, max_value=0;

  if (field < 0 || field >= NUM_SCALAR_VALUES)
    return (ERROR_BADPARM);
  if (thresh < 0 || thresh > 100.0)
    return (ERROR_BADPARM);
  if (mid < 0 || mid > 100.0)
    return (ERROR_BADPARM);
  if (max < 0 || max > 100.0)
    return (ERROR_BADPARM);

  sclv_get_value_for_percentile (field, thresh, &thresh_value);
  sclv_get_value_for_percentile (field, mid, &mid_value);
  sclv_get_value_for_percentile (field, max, &max_value);

  printf("%2.2f --> %2.2f\n", thresh, thresh_value) ;
  printf("%2.2f --> %2.2f\n", mid, mid_value) ;
  printf("%2.2f --> %2.2f\n", max, max_value) ;

  sclv_field_info[field].fthresh = thresh_value;
  sclv_field_info[field].fmid = mid_value;
  if (max_value - mid_value < epsilon)
  {
    sclv_field_info[field].fslope = 1.0;
  }
  else
  {
    sclv_field_info[field].fslope = 1.0 / (max_value - mid_value);
  }

  vertex_array_dirty = 1;

  printf ("fthresh %.2f fmid %.2f max %.2f slope %.2f\n",
          thresh_value, mid_value, max_value, sclv_field_info[field].fslope);

  return (ERROR_NONE);
}

int sclv_get_value_for_percentile (int field, float percentile, float* value)
{
  int target_count;
  int bin, sum;

  if (field < 0 || field >= NUM_SCALAR_VALUES)
    return (ERROR_BADPARM);
  if (percentile < 0 || percentile > 100.0)
    return (ERROR_BADPARM);

  target_count = (float)mris->nvertices * (percentile / 100.0);

  sum = 0;
  bin = 0;
  while (sum < target_count && bin < sclv_field_info[field].num_freq_bins)
  {
    sum += sclv_field_info[field].\
           frequencies[sclv_field_info[field].\
                       cur_condition][sclv_field_info[field].\
                                      cur_timepoint][bin];
    bin++;
  }

  *value = sclv_field_info[field].min_value +
           ( ((sclv_field_info[field].max_value -
               sclv_field_info[field].min_value + 1) * bin) \
             / sclv_field_info[field].num_freq_bins);

  return (ERROR_NONE);
}

int sclv_set_threshold_using_fdr (int field, float rate, int only_marked)
{
  float *saved_val;
  float *saved_val2;
  float current_value=0;
  double threshold;
  int err;
  int sign;
  int *saved_undefval = NULL;
  int vno;
  int num_marked;
  VERTEX* v;

  /* if they are truncating, they only want the positive or negative
     values. if inverse is on, they want the negative values, else
     they just want the positive. */
  sign = 0;
  if (truncphaseflag)
  {
    if (invphaseflag)
      sign = -1;
    else
      sign = 1;
  }

  /* since the fdr function only works on val and overwrites val2, we
     need to back up these fields, and then write our current field
     into val. we'll restore everything later. */
  saved_val = (float*) calloc (mris->nvertices, sizeof(float));
  if (NULL == saved_val)
  {
    ErrorReturn(ERROR_NO_MEMORY,
                (ERROR_NO_MEMORY,
                 "sclv_set_threshold_using_fdr: "
                 "couldn't allocated saved_val array\n"));
  }
  saved_val2 = (float*) calloc (mris->nvertices, sizeof(float));
  if (NULL == saved_val2)
  {
    ErrorReturn(ERROR_NO_MEMORY,
                (ERROR_NO_MEMORY,
                 "sclv_set_threshold_using_fdr: "
                 "couldn't allocated saved_val2 array\n"));
  }
  for (vno = 0; vno < mris->nvertices; vno++)
  {
    /* save val and val2. */
    v = &mris->vertices[vno];
    sclv_get_value (v, SCLV_VAL, &saved_val[vno]);
    sclv_get_value (v, SCLV_VAL2, &saved_val2[vno]);

    /* get the current value and put it into val. */
    sclv_get_value (v, field, &current_value);
    sclv_set_value (v, SCLV_VAL, current_value);
  }

  /* if we're only doing marked verts, go through the surface. save
     undefval, set undefval to 1 if marked. */
  if (only_marked)
  {
    saved_undefval = (int*) calloc (mris->nvertices, sizeof(int));
    if (NULL == saved_undefval)
    {
      ErrorReturn(ERROR_NO_MEMORY,
                  (ERROR_NO_MEMORY,
                   "sclv_set_threshold_using_fdr: "
                   "couldn't allocated saved_undefval array\n"));
    }
    num_marked = 0;
    for (vno = 0; vno < mris->nvertices; vno++)
    {
      v = &mris->vertices[vno];
      saved_undefval[vno] = v->undefval;
      if (v->marked)
      {
        v->undefval = 1;
        num_marked++;
      }
      else
        v->undefval = 0;
    }
    printf ("surfer: performing FDR on %d vertices\n", num_marked);
  }

  /* call the fdr function. we pass our surface, the rate, the sign we
     got before, only_marked for masked (we set undefval before), and
     get the threshold back. */
  err = MRISfdr2vwth(mris, rate, sign, 1, only_marked, &threshold);
  fprintf (stderr, "MRISfdr2vwth(rate=%f, sign=%d, 1, only_marked=%d) = %f\n",
           rate, sign, only_marked, threshold);
  if ( err )
  {
    printf ("surfer: Error calculating threshold with FDR.\n");
    if (only_marked) free (saved_undefval);
    return (err);
  }
  printf ("surfer: MRISfdr2vwth with rate %.2f and sign %d returned "
          "threshold %f\n", rate, sign, threshold );

  /* we we're only doing marked verts, go through and restore the
     undefval values. */
  if (only_marked)
  {
    for (vno = 0; vno < mris->nvertices; vno++)
    {
      v = &mris->vertices[vno];
      v->undefval = saved_undefval[vno];
    }
    free (saved_undefval);
  }

  /* restore val and val2 */
  for (vno = 0; vno < mris->nvertices; vno++)
  {
    v = &mris->vertices[vno];
    sclv_set_value (v, SCLV_VAL, saved_val[vno]);
    sclv_set_value (v, SCLV_VAL2, saved_val2[vno]);
  }
  free (saved_val);
  free (saved_val2);

  sclv_field_info[field].fthresh = threshold;
  sclv_field_info[field].fmid = threshold + 1.5;
  sclv_field_info[field].fslope = 0.66;

  if (field == sclv_current_field)
  {
    fthresh = sclv_field_info[field].fthresh;
    fmid = sclv_field_info[field].fmid;
    fslope = sclv_field_info[field].fslope;
    vertex_array_dirty = 1;
    sclv_send_current_field_info();
  }

  return (ERROR_NONE);
}

int sclv_copy_view_settings_from_field (int field, int fromfield)
{
  if (field < 0 || field >= NUM_SCALAR_VALUES)
    return (ERROR_BADPARM);
  if (fromfield < 0 || fromfield >= NUM_SCALAR_VALUES)
    return (ERROR_BADPARM);

  /* if fromfield is the current field, update from the shared variables. */
  if (fromfield == sclv_current_field)
  {
    sclv_field_info[fromfield].fthresh = fthresh;
    sclv_field_info[fromfield].fmid = fmid;
    sclv_field_info[fromfield].foffset = foffset;
    sclv_field_info[fromfield].fslope = fslope;
  }

  sclv_field_info[field].fthresh = sclv_field_info[fromfield].fthresh;
  sclv_field_info[field].fmid = sclv_field_info[fromfield].fmid;
  sclv_field_info[field].foffset = sclv_field_info[fromfield].foffset;
  sclv_field_info[field].fslope = sclv_field_info[fromfield].fslope;

  /* if this is the current field, update the shared variables and
     send the current info. */
  if (field == sclv_current_field)
  {
    fthresh = sclv_field_info[sclv_current_field].fthresh;
    fmid = sclv_field_info[sclv_current_field].fmid;
    ;
    foffset = sclv_field_info[sclv_current_field].foffset;
    ;
    fslope = sclv_field_info[sclv_current_field].fslope;

    sclv_send_current_field_info();
  }

  return (ERROR_NONE);
}

int sclv_copy_view_settings_from_current_field (int field)
{
  sclv_copy_view_settings_from_field (field, sclv_current_field);

  return (ERROR_NONE);
}

int sclv_copy_all_view_settings_from_current_field ()
{
  int field;

  for (field = 0; field < NUM_SCALAR_VALUES; field++)
  {
    if (field != sclv_current_field)
      sclv_copy_view_settings_from_field (field, sclv_current_field);
  }

  return (ERROR_NONE);
}

int sclv_swap_fields ( int fielda, int fieldb )
{

  SCLV_FIELD_INFO swap_field;
  char cmd[STRLEN];
  int k;
  float a, b;
  a=0;
  b=0;

  if (fielda < 0 || fielda >= NUM_SCALAR_VALUES)
    return (ERROR_BADPARM);
  if (fieldb < 0 || fieldb >= NUM_SCALAR_VALUES)
    return (ERROR_BADPARM);

  /* swap the field values */
  for (k=0;k<mris->nvertices;k++)
  {
    sclv_get_value(&(mris->vertices[k]), fielda, &a );
    sclv_get_value(&(mris->vertices[k]), fieldb, &b );
    sclv_set_value(&(mris->vertices[k]), fielda, b );
    sclv_set_value(&(mris->vertices[k]), fieldb, a );
  }

  /* swap the field data */
  swap_field = sclv_field_info[fielda];
  sclv_field_info[fielda] = sclv_field_info[fieldb];
  sclv_field_info[fieldb] = swap_field;

  /* swap the field names in the interface */
  sprintf (cmd, "SwapValueLabelNames %d %d", fielda, fieldb);
  send_tcl_command (cmd);

  return (ERROR_NONE);
}

int sclv_send_histogram ( int field )
{

  float increment;
  char *tcl_cmd;
  int condition, timepoint;
  int bin;
  int count;

  /* calculate the number of values and the increment between
     each value */
  increment =
    (sclv_field_info[field].max_value - sclv_field_info[field].min_value) /
    (float)(sclv_field_info[field].num_freq_bins);

  /* start a string of the proper size; give us the length of the
     command, and then 10 characters per number, begin + end +
     increment + num values */
  // Careful -- tcl_cmd string len can easily exceed STRLEN
  tcl_cmd = (char*)calloc(21 + (sclv_field_info[field].num_freq_bins + 4) * 10,
                          sizeof(char));
  if (NULL == tcl_cmd)
  {
    return (ERROR_NO_MEMORY);
  }

  /* add the command name to the string, the min value, the max value,
     the number of values, and the increment */
  sprintf (tcl_cmd, "UpdateHistogramData %.5f %.5f %.5f %d {",
           sclv_field_info[field].min_value, sclv_field_info[field].max_value,
           increment, sclv_field_info[field].num_freq_bins);

  /* for each frequency bin, add the value. */
  condition = sclv_field_info[field].cur_condition;
  timepoint = sclv_field_info[field].cur_timepoint;
  for (bin = 0; bin < sclv_field_info[field].num_freq_bins; bin++)
  {
    count = sclv_field_info[field].frequencies[condition][timepoint][bin];

    /* If this is the zero bin, switch on
       ignorezeroesinhistogramflag to see whether or not we include
       the count in there. If we are ignoring, don't include it,
       otherwise add it in. */
    if (bin == sclv_field_info[field].zero_bin_index)
    {
      if (!ignorezeroesinhistogramflag)
        count += sclv_field_info[field].num_zeroes_in_zero_bin;
    }

    /* Add to the command string. */
    sprintf (tcl_cmd, "%s %d", tcl_cmd, count);
  }

  /* close up the command and send it off */
  strcat (tcl_cmd, "}");
  send_tcl_command (tcl_cmd);

  free (tcl_cmd);

  return (ERROR_NONE);
}

int sclv_get_normalized_color_for_value (int field, float value,
    float *outRed,
    float *outGreen,
    float *outBlue)
{
  GLubyte r, g, b;
  get_color_vals (value, 0, REAL_VAL, &r, &g, &b);
  *outRed = ((float)r / 255.0);
  *outGreen = ((float)g / 255.0);
  *outBlue = ((float)b / 255.0);
  return (ERROR_NONE);
}

int sclv_apply_color_for_value (float f, float opacity,
                                GLubyte* pr, GLubyte* pg, GLubyte* pb )
{
  float r,g,b;
  float ftmp,c1,c2;
  float min, mid, max;
  float or, ob, og;
  float br, bg, bb;
  float tmpoffset, f2, fr, fg, fb;
  // extern double fcurv; // sets curv thresh

  /* Adjust by foffset. */
  f -= foffset; // foffset default is 0 (can be changed on gui)

  r = g = b = 0.0f ;
  if (invphaseflag)           f = -f;
  if (truncphaseflag && f<0)  f = 0;
  if (rectphaseflag)          f = fabs(f);

  /* rkt: same way values are calc'd in tkmedit. The main difference
     is that max is 0.5/slope + mid instead of 1/slope + mid, to make
     the linear version work better. */
  min = (float)(fthresh);
  mid = (float)(fmid);
  max = (0.5 / (float)fslope) + (float)fmid;

  /* Calculate the background colors. */
  br = (float)*pr / 255.0;
  bg = (float)*pg / 255.0;
  bb = (float)*pb / 255.0;

  // Apparently, fcurv is always 0, which would mean this 
  // section of code does nothing (dng)
  if (fabs(f)>fthresh && fabs(f)<fmid) {
    ftmp = fabs(f);
    c1 = 1.0/(fmid-fthresh);
    if (fcurv!=1.0)
      c2 = (fmid-fthresh-fcurv*c1*SQR(fmid-fthresh))/
           ((1-fcurv)*(fmid-fthresh));
    else
      c2 = 0;
    ftmp = fcurv*c1*SQR(ftmp-fthresh)+c2*(1-fcurv)*(ftmp-fthresh)+fthresh;
    f = (f<0)?-ftmp:ftmp;
  }

  if(colscale==HEAT_SCALE) {
    if(f>=0){
      if(sclv_opaque){
        /* If opaque, don't use blending at all. Min->mid is all
           red, and mid->max gets yellower. Who decided that this
	   was a good idea? */
	// br,bg,bb are background values
	// f<min means f>0 AND f<min
        //r = ((f<min) ? br : 1.0);
        //g = ((f<min) ? bg : ( f<mid) ? 0 : (f<max) ? (f-mid)/(max-mid) : 1.0 );
        //b = ((f<min) ? bb : 0);
	// Here's my version, it actually works.
	if(f >= min){
	  ftmp = (f-min)/(max-min); // normalize
	  dngheat(ftmp, &r, &g, &b);
	} else {
	  r=br; g=bg; b=br;
	}
      }
      else{
        /* the offset is a portion of the color that is 'blended'
           into the functional color so that a func value right at
           the threshold doesn't look black, but translucent. the
           rest is a standard interpolated color scale. */
        or = br * ((f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0);
        og = bg * ((f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0);
        ob = bb * ((f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0);
        r = or + ((f<min) ? 0.0 : (f<mid) ? (f-min)/(mid-min) : 1.0);
        g = og + ((f<mid) ? 0.0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
        b = ob;
      }
    }
    else { // f < 0
      f = -f;
      if(sclv_opaque) {
        //b = ((f<min) ? bb : 1.0);
        //g = ((f<min) ? bg : (f<mid) ? 0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
        //r = ((f<min) ? br : 0);
	if(f >= min){
	  ftmp = (f-min)/(max-min); // normalize
	  dngheat(-ftmp, &r, &g, &b);
	} else {
	  r=br; g=bg; b=br;
	}
      }
      else
      {
        or = br * ((f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0);
        og = bg * ((f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0);
        ob = bb * ((f<min) ? 1.0 : (f<mid) ? 1.0 - (f-min)/(mid-min) : 0.0);
        b = ob + ((f<min) ? 0.0 : (f<mid) ? (f-min)/(mid-min) : 1.0);
        g = og + ((f<mid) ? 0.0 : (f<max) ? (f-mid)/(max-mid) : 1.0);
        r = or;
      }
    }
    r = r*255;
    g = g*255;
    b = b*255;
  }
  else if (colscale==CYAN_TO_RED ||
           colscale==BLU_GRE_RED ||
           colscale==JUST_GRAY)
  {
    tmpoffset = (float)*pr;
    if (f<fthresh)
    {
      r = g = 255 * (tmpoffset/blufact);
      b =     255 * (tmpoffset*blufact);
    }
    else
    {
      if (fslope!=0)
        f2 = (tanh(fslope*fmid)+
              tanh(fslope*(f-fmid)))/(2-tanh(fslope*fmid));
      else
        f2 = (f<0)?0:((f>1)?1:f);
      set_positive_color(f2,&fr,&fg,&fb,tmpoffset);
      r=fr;
      g=fg;
      b=fb;
    }
  }
  else if (colscale==BLUE_TO_RED_SIGNED ||
           colscale==GREEN_TO_RED_SIGNED)
  {
    tmpoffset = (float)*pr;
    if (fabs(f)>fthresh)
    {
      if (fslope!=0)
      {
        if (fmid==0)
          f2 = tanh(fslope*(f));
        else
        {
          if (f<0)
            f2 = -(tanh(fslope*fmid) + tanh(fslope*(-f-fmid)))/
                 (2-tanh(fslope*fmid));
          else
            f2 = (tanh(fslope*fmid) + tanh(fslope*( f-fmid)))/
                 (2-tanh(fslope*fmid));
        }

      }
      else
      {
        f2 = (f<-1)?-1:((f>1)?1:f);
      }
      set_signed_color(f2,&fr,&fg,&fb,tmpoffset);
      r=fr;
      g=fg;
      b=fb;
    }
  }

  /* Blend the color into the input color with the given opacity.*/
  *pr = (int)((1.0 - opacity) * (float)*pr) + (opacity * r);
  *pg = (int)((1.0 - opacity) * (float)*pg) + (opacity * g);
  *pb = (int)((1.0 - opacity) * (float)*pb) + (opacity * b);

  return (NO_ERROR);
}


int sclv_load_label_value_file (char *fname, int field)
{
  FILE* fp = NULL;
  int line_number = 0;
  char line[1024] = "";
  int num_read = 0;
  char label_name[256] = "";
  float value = 0;
  char hemisphere[256] = "";
  float min, max;
  int label_index = 0;
  LABEL* label = NULL;
  int label_vno = 0;
  int vno = 0;
  VERTEX* v = NULL;
  char val_name[1024];
  char cmd[1024];

  /* unload this field if it already exists */
  sclv_unload_field (field);

  /* Go through the file line by line. */
  fp = fopen (fname, "r");
  if (NULL == fp)
  {
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "sclv_load_label_value_file: "
                 "couldn't read file %s\n", fname));
  }

  line_number = 1;
  min = 9999;
  max = -9999;
  while (!feof(fp))
  {
    /* Read a label name, value, and hemisphere. */
    fgetl (line, 1024, fp);
    num_read = sscanf (line, "%s %f %s", label_name, &value, hemisphere);
    if (3 != num_read)
    {
      printf ("sclv_load_label_value_file: error reading line "
              "%d: %s\n", line_number, line);
      continue;
    }

    /* If this is our hemisphere... */
    if (0 == strcmp (hemisphere, stem))
    {

      /* Update min/max. */
      if (value < min)
        min = value;
      if (value > max)
        max = value;

      /* Go through the labels. If one's name matches this label name... */
      for (label_index = 0; label_index < labl_num_labels; label_index++)
      {
        if (0 == strcmp (labl_labels[label_index].name, label_name))
        {
          label = labl_labels[label_index].label;

          /* Go through the label points and set the sclv to
             this value. */
          for (label_vno = 0; label_vno < label->n_points; label_vno++)
          {
            vno = label->lv[label_vno].vno;
            if ( vno < 0 || vno >= mris->nvertices )
              continue;

            v = &mris->vertices[vno];

            sclv_set_value (v, field, value);
          }
        }
      }
    }

    line_number++;
  }

  fclose (fp);


  /* mark this field as not binary */
  sclv_field_info[field].is_functional_volume = FALSE;
  sclv_field_info[field].is_scalar_volume = FALSE;

  /* save the range */
  sclv_field_info[field].min_value = min;
  sclv_field_info[field].max_value = max;

  /* dummy info for time point and conditions, since .w files only have
     one plane of info */
  sclv_field_info[field].cur_timepoint = 0;
  sclv_field_info[field].cur_condition = 0;
  sclv_field_info[field].num_timepoints = 1;
  sclv_field_info[field].num_conditions = 1;

  /* calc the frquencies */
  sclv_calc_frequencies (field);

  /* request a redraw. turn on the overlay flag and select this value set */
  vertex_array_dirty = 1 ;
  overlayflag = TRUE;
  sclv_set_current_field (field);

  /* set the field name to the name of the file loaded */
  if (NULL != g_interp)
  {
    FileNameOnly (fname, val_name);
    sprintf (cmd, "UpdateValueLabelName %d \"%s\"", field, val_name);
    send_tcl_command (cmd);
    sprintf (cmd, "ShowValueLabel %d 1", field);
    send_tcl_command (cmd);
    sprintf (cmd, "UpdateLinkedVarGroup view");
    send_tcl_command (cmd);
    sprintf (cmd, "UpdateLinkedVarGroup overlay");
    send_tcl_command (cmd);
  }

  /* enable the menu items */
  enable_menu_set (MENUSET_OVERLAY_LOADED, 1);

  return (ERROR_NONE);
}

/* krish -- linktimepoint with ROI(avg or normalized avg ) mode */
void link_timepoint_ROI(int vno)
{
    LABEL* area;
    int n, current_vno, tmpvno;
    float *values, *average, *stddev;
    int label_index = -1;
    int num_found;

    /* try and find the label corresponding to the clicked vertex number */
    labl_find_label_by_vno (vno, labl_selected_label+1,
			  &label_index, 1, &num_found);
    /* idea is to get the list of vertices belonging to this label (area) */ 
    area = labl_labels[label_index].label ;
    /* allocate memory for values and average*/
    values = (float*) calloc( mris->nvertices, sizeof(float) );
    average = (float*) calloc( mris->nvertices, sizeof(float) );
    stddev = (float*) calloc( mris->nvertices, sizeof(float) );
    if (values == NULL || average == NULL || stddev == NULL){
      printf("link_timepoint_ROI(): couldn't allocate one or more of values/average/stddev storage .\n");
      return;
    }
    /* for every vertex in the area ( label ) , get the associated timepoint values 
     * add them all to the average array ( also calculate variance)*/
    for (n = 0  ; n < area->n_points ; n++)
    {
      if (area->lv[n].vno > 0 && area->lv[n].vno < mris->nvertices)
      {
	current_vno = area->lv[n].vno ;
	sclv_get_values_for_field_and_timepoint (sclv_current_field, current_vno,
						 sclv_cur_condition, values);
	for (tmpvno=0; tmpvno < mris->nvertices; tmpvno++) 
	  average[tmpvno] += values[tmpvno] / area->n_points;
	/* calculate the variance */
        if (linkvertexmode == 3)
          for (tmpvno=0; tmpvno < mris->nvertices; tmpvno++) 
	    stddev[tmpvno] += pow(values[tmpvno]-average[tmpvno], 2);
      }
    }
        
      
    /* in every vertex v, set average[vno] or stddev[vno] as the timepoint value */
    for (tmpvno = 0 ; tmpvno < mris->nvertices ; tmpvno++)
    { 
      VERTEX *v;
      v = &mris->vertices[tmpvno] ;
      if (linkvertexmode == 2) 
        sclv_set_value(v, sclv_current_field, average[tmpvno] );
      if ( linkvertexmode == 3 ) {
        stddev[tmpvno] = sqrt(stddev[tmpvno] / area->n_points);
	if ( fabs(stddev[tmpvno]) < epsilon ) {
	  sclv_set_value(v, sclv_current_field, average[tmpvno] / stddev[tmpvno] );
	}
	else {
	  sclv_set_value(v, sclv_current_field, average[tmpvno] );
	}
      }
    }

    free(values);
    free(average);
    free(stddev);
    send_tcl_command ("UpdateAndRedraw");
}
/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------ multiple labels */

int labl_initialize ()
{

  int label;

  /* initialize the array of labels to empty. */
  for (label = 0; label < LABL_MAX_LABELS; label++)
  {
    labl_labels[label].label = NULL;
    labl_labels[label].structure = -1;
    labl_labels[label].r = 0;
    labl_labels[label].g = 0;
    labl_labels[label].b = 0;
    labl_labels[label].visible = 0;
    labl_labels[label].border_vno = NULL;
    labl_labels[label].num_border_vnos = 0;
    labl_labels[label].cached = 0;
  }

  labl_num_labels = 0;
  labl_selected_label = LABL_NONE_SELECTED;
  labl_ctab = NULL;
  //labl_draw_style = LABL_STYLE_FILLED; // initiallized in declaration
  labl_num_labels_created = 0;
  labl_color_table_name = (char*) calloc (NAME_LENGTH, sizeof(char));
  labl_draw_flag = 1;
  labl_cache = NULL;
  labl_num_labels_at_cache_vno = NULL;
  labl_cache_updated = FALSE;

  return (ERROR_NONE);
}

int labl_update_cache (int force_rebuild_all)
{

  int label_index;
  LABEL* label;
  int label_vno;
  int vno;
  int* new;

  /* if our basic stuff is not inited, do so now. */
  if (NULL == labl_cache)
  {
    labl_cache = (int**) calloc (mris->nvertices, sizeof(int*));
  }
  if (NULL == labl_num_labels_at_cache_vno)
  {
    labl_num_labels_at_cache_vno =
      (int*) calloc (mris->nvertices, sizeof(int));
  }

  if (force_rebuild_all)
  {
    for (vno = 0; vno < mris->nvertices; vno++)
      labl_num_labels_at_cache_vno[vno] = 0;
  }

  for (label_index = 0; label_index < labl_num_labels; label_index++)
  {
    label = labl_labels[label_index].label;
    if (!labl_labels[label_index].cached || force_rebuild_all)
    {
      for (label_vno = 0; label_vno < label->n_points; label_vno++)
      {
        vno = label->lv[label_vno].vno;
        if ( vno < 0 || vno >= mris->nvertices )
          continue;

        if (NULL == labl_cache[vno])
        {
          labl_cache[vno] = (int*) calloc (1, sizeof(int));
          labl_cache[vno][0] = label_index;
          labl_num_labels_at_cache_vno[vno] = 1;
        }
        else
        {
          new = (int*)
                calloc (labl_num_labels_at_cache_vno[vno]+1, sizeof(int));
          memmove (new, labl_cache[vno],
                  labl_num_labels_at_cache_vno[vno] * sizeof(int));
          free (labl_cache[vno]);
          labl_cache[vno] = new;
          labl_cache[vno][labl_num_labels_at_cache_vno[vno]] =
            label_index;
          labl_num_labels_at_cache_vno[vno]++;
        }
      }
      labl_labels[label_index].cached = TRUE;
    }
  }

  labl_cache_updated = TRUE;

  return (NO_ERROR);
}

int labl_load_color_table (char* fname)
{
  COLOR_TABLE* ctab;
  int label_index;
  LABL_LABEL* label;
  int r, g, b, a;

  /* Attempt to read the color table. */
  ctab = CTABreadASCII (fname);
  if (NULL == ctab)
  {
    ErrorReturn(ERROR_BADFILE, (ERROR_BADFILE, "Couldn't open %s\n", fname));
  }

  /* Save it in the surface. */
  if (NULL != mris->ct)
  {
    CTABfree (&mris->ct);
  }
  mris->ct = ctab;

  /* save the name of the color table */
  strcpy (labl_color_table_name, fname);

  /* send the color table to tcl. */
  labl_send_color_table_info ();

  /* Update all currently loaded structural labels. */
  for (label_index = 0; label_index < labl_num_labels; label_index++)
  {
    label = &(labl_labels[label_index]);
    if (LABL_TYPE_FREE != label->structure)
    {
      CTABcopyName (mris->ct, label->structure,
                    label->name, sizeof(label->name));
      r = 255;
      g = b = 0;
      a = 255;
      CTABrgbaAtIndexi (mris->ct, label->structure, &r, &g, &b, &a);
      labl_set_color (label_index, r, g, b);
      labl_set_alpha (label_index, a);
      labl_send_info (label_index);
    }
  }

  redraw();

  return (ERROR_NONE);
}

int
labl_send_color_table_info ()
{
  COLOR_TABLE* ctab;
  int num_valid_entries;
  int num_entries;
  int structure;
  char structure_label[1024];
  char* structure_label_list = NULL;
  int valid;

  /* if we have our own table, get the names from there, otherwise
     use our external table. */
  if (NULL != mris->ct)
    ctab = mris->ct;
  else
    ctab = labl_ctab;

  /* Find out how many valid and total entries we have. */
  CTABgetNumberOfValidEntries (ctab, &num_valid_entries);
  CTABgetNumberOfTotalEntries (ctab, &num_entries);

  /* allocate a string long enough for the update command and all
     our labels. */
  structure_label_list = (char*)
                         malloc (256 * num_valid_entries * sizeof(char));
  if (NULL != structure_label_list)
  {
    /* build a string out of all the label names and send them to the
       tcl label list. */
    strcpy (structure_label_list, "LblLst_SetStructures {");

    /* Iterate over all the entries, but only get names for the
       valid ones. */
    for (structure = 0; structure < num_entries; structure++ )
    {
      /* If not valid, skip it. */
      CTABisEntryValid (ctab, structure, &valid);
      if (!valid)
        continue;

      CTABcopyName (ctab, structure,
                    structure_label, sizeof(structure_label));
      sprintf (structure_label_list, "%s %d %s",
               structure_label_list, structure, structure_label);
    }
    sprintf (structure_label_list, "%s }", structure_label_list);
    send_tcl_command (structure_label_list);

    free( structure_label_list );
  }
  else
  {
    fprintf (stderr, "labl_send_color_table_info: couldn't allocate "
             "string for %d structres\n", labl_num_structures );
    return (ERROR_NO_MEMORY);
  }

  return (NO_ERROR);
}

int labl_load (char* fname)
{
  LABEL* label = NULL;
  LABEL *lnew;
  int label_index;
  char name[NAME_LENGTH];
  int unassigned;

  if (white_surf_loaded == 0)
    read_white_vertex_coordinates() ;

  if (NULL == fname)
    return (ERROR_BADPARM);

  if (MRISfileNameType(fname) == MRIS_GIFTI_FILE)
  {
    printf("\n\n INFO: use File->Label->Import Annotation to load Gifti "
           "format label files!\n\n");
    return(ERROR_NO_FILE);
  }

  /* load label file. */
  label = LabelRead (pname, fname);
  if (NULL == label)
  {
    return (ERROR_NO_FILE);
  }
  {
    int n, nonzero = 0 ;

    for (n = 0 ; n < label->n_points ; n++)
      if (!FZERO(label->lv[n].stat))
        nonzero++ ;
    printf("********************      %d nonzero vertices found ********************\n", nonzero) ;
    if (nonzero== 0)
    {
      printf("label stat field identically zero - setting to 1\n") ;
      for (n = 0 ; n < label->n_points ; n++)
        label->lv[n].stat = 1 ;
    }
  }

  if (reassign)
    LabelUnassign(label) ;

  /* load the orig vertex positions if we haven't already. */
  if (!origsurfloaded)
    read_orig_vertex_coordinates(orfname) ;
  LabelToWhite (label, mris);

  /* See if the label is completely unassigned. If it is, we'll fill
     it after we assign the verts. */
  LabelIsCompletelyUnassigned (label, &unassigned);

  /* assign mris vertex numbers to unnumbered vertices based on their
     locations. */
  LabelFillUnassignedVertices (mris, label, WHITE_VERTICES);

  /* If we were unassigned before, fill it now. */
  if (unassigned)
  {
    lnew = LabelFillHoles(label, mris, WHITE_VERTICES) ;
    LabelFree(&label) ;
    label = lnew ;
  }

  /* make a new entry in the label list. */
  labl_add (label, &label_index);

  /* set the name to the tail of the filename. make this a free label,
     i.e. with no assigned strucutre. make it visible. */
  FileNameOnly (fname, name);
  labl_set_info (label_index, name, LABL_TYPE_FREE, 1,
                 LABL_DEFAULT_COLOR_R, LABL_DEFAULT_COLOR_G,
                 LABL_DEFAULT_COLOR_B );

  /* Unmark it. */
  LabelUnmark (label, mris);

  /* select this label */
  labl_select (label_index);

  surface_compiled = 0 ;


  return (ERROR_NONE);
}

int labl_save (int index, char* fname)
{
  LABEL* label = NULL;
  int n;
  int vno;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);
  if (NULL == fname)
    return (ERROR_BADPARM);

  /* Get a copy of our label. */
  printf("label at index %d with %d points\n",
         index, labl_labels[index].label->n_points) ;
  label = LabelCopy (labl_labels[index].label, NULL);
  if (NULL == label)
    return (ERROR_BADPARM);

  /* read original vertex positions if we haven't already. */
  if (!origsurfloaded)
    read_orig_vertex_coordinates(orfname) ;

  /* fill in the current overlay value in the stat field of our label
     verts. this is for the LabelWrite call. */
  for (n = 0 ; n < label->n_points ; n++)
  {
    vno = label->lv[n].vno;
    if ( vno < 0 || vno >= mris->nvertices )
      continue;

    sclv_get_value (&(mris->vertices[vno]),
                    sclv_current_field, &(label->lv[n].stat) );
  }

  /* write the label. */
  fprintf (stderr, "writing %d labeled vertices to %s.\n",
           label->n_points, fname) ;
  LabelToWhite (label, mris);
  LabelWrite (label, fname);

  /* Delete the label. */
  LabelFree (&label);

  return (ERROR_NONE);
}

int labl_save_all (char* prefix)
{
  int label;
  char fname[NAME_LENGTH];

  /* for each label we have, build a decent name based on the prefix
     we got, and save it normally. */
  for (label = 0; label < labl_num_labels; label++ )
  {
    sprintf (fname, "%s-%d", prefix, label);
    labl_save (label, fname);
  }

  return (ERROR_NONE);
}


int labl_find_and_set_all_borders ()
{
  int label;

  /* for each label we have, find the border */
  for (label = 0; label < labl_num_labels; label++ )
  {
    labl_find_and_set_border (label);
  }

  return (ERROR_NONE);
}

int labl_find_and_set_border (int index)
{
  int label_vno;
  LABEL* label;
  VERTEX* v;
  int vno;
  int neighbor_vno;
  char* border = NULL;
  int num_borders = 0;
  int label_index_array[LABL_MAX_LABELS];
  int num_labels_found, found_label_index;
  int vno_in_label, vno_in_other_label;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* make an array of border flags for just the outline. */
  border = (char*) calloc (mris->nvertices, sizeof(char));
  num_borders = 0;

  /* for each vertex in the label... */
  label = labl_labels[index].label;
  for (label_vno = 0; label_vno < label->n_points; label_vno++)
  {
    /* The label vno could be < 0 if it couldn't find a close vertex
       for this label point. */
    vno = label->lv[label_vno].vno;
    if (vno < 0 || vno >= mris->nvertices)
      continue;

    /* get the vno and look at this vertex in the mris. for each
       neighbor, if we still haven't determined that this vno is a
       border, the neighbor is not in the same label... */
    v = &(mris->vertices[vno]);
    for (neighbor_vno = 0;
         neighbor_vno < v->vnum && !border[vno];
         neighbor_vno++ )
    {
      /* See if there's a label at this neighbor vno. */
      labl_find_label_by_vno (v->v[neighbor_vno], 0, label_index_array,
                              LABL_MAX_LABELS, &num_labels_found);
      if ( num_labels_found > 0 )
      {
        /* We found at least one label, check to see if it's
           this label or somethign else. */
        vno_in_other_label = 0;
        vno_in_label = 0;
        for (found_label_index = 0;
             found_label_index < num_labels_found;
             found_label_index++)
        {
          if (label_index_array[found_label_index] != index)
            vno_in_other_label = 1;
          if (label_index_array[found_label_index] == index)
            vno_in_label = 1;
        }

        /* it's a border if this vno is in a different (or no)
           label and not in this label. mark it so in the real
           border flag array and our temp one. */
        if (vno_in_other_label && !vno_in_label &&
            !border[vno])
        {
          border[vno] = TRUE;
          num_borders++;
        }
      }
      else
      {
        /* it's also a border if the vno is in no label. */
        border[vno] = TRUE;
        num_borders++;
      }
    }

    /* Now we want to expand the border. If we determined this is a
       border, check its surrounding vertices. For each one, if it's
       not already also a border, and if it's inside the label, mark
       it in the border list. */
    if (border[vno])
    {
      for (neighbor_vno = 0;
           neighbor_vno < v->vnum && !border[v->v[neighbor_vno]];
           neighbor_vno++ )
      {
        /* If it's inside the label...*/
        labl_find_label_by_vno (v->v[neighbor_vno], 0, label_index_array,
                                LABL_MAX_LABELS, &num_labels_found);
        if ( num_labels_found > 0 )
        {
          for (found_label_index = 0;
               found_label_index < num_labels_found;
               found_label_index++)
          {
            if (label_index_array[found_label_index] == index &&
                !border[v->v[neighbor_vno]])
            {
              /* Also make it a border. */
              border[v->v[neighbor_vno]] = TRUE;
              num_borders++;
            }
          }
        }
      }
    }
  }

  /* Allocate the array on the label and go through the border array,
     writing index numbers to the label array. */
  if (NULL != labl_labels[index].border_vno)
    free (labl_labels[index].border_vno);

  labl_labels[index].border_vno = (int*) calloc (num_borders, sizeof(int));

  labl_labels[index].num_border_vnos = 0;
  for (label_vno = 0; label_vno < label->n_points; label_vno++)
  {
    vno = label->lv[label_vno].vno;
    if (vno < 0 || vno >= mris->nvertices)
      continue;

    v = &(mris->vertices[vno]);
    if (border[vno])
    {
      labl_labels[index].border_vno[labl_labels[index].num_border_vnos] =
        vno;
      labl_labels[index].num_border_vnos++;

      /* By unmarking it now, we make sure that a duplicate label
         point will not also be added. */
      border[vno] = FALSE;
    }
  }

  free (border);

  return (ERROR_NONE);
}

int labl_vno_is_border (int index, int vno)
{
  int border_index;

  if (index < 0 || index >= labl_num_labels)
    return 0;

  if (NULL == labl_labels[index].border_vno)
    return 0;

  for (border_index = 0;
       border_index < labl_labels[index].num_border_vnos; border_index++)
  {
    if (labl_labels[index].border_vno[border_index] == vno)
      return 1;
  }

  return 0;
}

int labl_import_annotation (char *fname)
{
  int mris_err;
  int ctab_err;
  COLOR_TABLE* ctab;
  int annotation_vno;
  int vno;
  unsigned int annotation, max_annot;
  int num_verts_in_annotation;
  LABEL* label = NULL;
  int label_vno;
  VERTEX* v = NULL;
  int new_index;
  char name[NAME_LENGTH];
  int r, g, b;
  int structure;
  unsigned int* done;
  int num_labels;

  /* init our done array. */
  r = g = b = 255;
  MRISRGBToAnnot(r,g,b,max_annot);
  done = (unsigned int*) calloc( max_annot, sizeof(int) );
  if ( NULL == done )
  {
    printf( "calloc of size %d failed\n", max_annot );
    return (ERROR_NO_MEMORY);
  }
  num_labels = 0;

  /* read the annotation. */
  for (vno = 0; vno < mris->nvertices; vno++)
    mris->vertices[vno].annotation = 0;
  mris_err = MRISreadAnnotation (mris, fname);
  if (mris_err)
  {
    printf("\n");
    printf("ERROR: could not load %s\n",fname);
    printf("\n");
    return (ERROR_NO_FILE);
  }

  /* Check if we got an embedded color table. */
  if (mris->ct)
  {
    printf ("Found embedded color table in annotation.\n");
  }
  else
  {
    printf ("No embedded color table found in annotation.\n");
  }

  int vertices_wo_annotation = 0; // count vertices that dont have annotation

  /* check all annotations... */
  for (annotation_vno = 0; annotation_vno < mris->nvertices; annotation_vno++)
  {
    /* get the annotation. if there is one... */
    annotation = mris->vertices[annotation_vno].annotation;

    if (annotation)
    {
      /* get the rgb colors. */
      if (annotation > max_annot)
      {
        printf("Warning: vertex %d with annotation %x - out of range!\n",
               annotation_vno, annotation) ;
        annotation &= max_annot ;
      }
      MRISAnnotToRGB( annotation, r, g, b );

      /* if we haven't imported this label yet... */
      if ( !done[annotation] )
      {

        /* mark it imported. */
        done[annotation] = 1;
        num_labels++;

        /* find out how many verts have this annotation value. */
        num_verts_in_annotation = 0;
        for (vno = 0; vno < mris->nvertices; vno++)
          if (mris->vertices[vno].annotation == annotation)
            num_verts_in_annotation++;

        /* make a new label, and go through again, setting the label
           values. */
        label = LabelAlloc(num_verts_in_annotation, NULL, NULL);
        if (NULL != label)
        {
          strncpy( label->subject_name, pname, 100 );
          label->n_points = num_verts_in_annotation;
          label_vno = 0;
          for (vno = 0; vno < mris->nvertices; vno++)
            if (mris->vertices[vno].annotation == annotation)
            {
              v = &mris->vertices[vno];
              label->lv[label_vno].x = v->x;
              label->lv[label_vno].y = v->y;
              label->lv[label_vno].z = v->z;
              label->lv[label_vno].vno = vno;
              label_vno++;
            }

          /* add the label to our list. */
          labl_add (label, &new_index);

          /* now we need to set the information about the
             label from a color table. older parcellation
             files use the external color table, but newer
             ones have their own. so we'll check the ct
             member; if it's null, use the external, otherwise
             use the color info specified in the mris. */
          if (mris->ct)
            ctab = mris->ct;
          else
            ctab = labl_ctab;

          /* If not found, structure will be -1. */
          CTABfindRGBi (ctab, r, g, b, &structure);

          /* make a name for it. if we got a color from the
             color table, get the label, else use the color. */
          if (structure != -1)
          {
            ctab_err = CTABcopyName (ctab, structure,
                                     name, sizeof(name) );
            if (NO_ERROR != ctab_err)
              sprintf (name, "Parcellation %d, %d, %d", r, g, b);
          }
          else
          {
            sprintf (name, "Parcellation %d, %d, %d", r, g, b);
          }

          /* set its other data. set the color; if we found a
             structure index from the LUT, it will use that,
             otherwise it will color it as a free label with
             the given colors (which really has the same
             effect, just doesn't give it a valid structure
             index. */
          labl_set_info (new_index, name, structure, 1, r, g, b);
        }
      }
    }
    else
    {
      vertices_wo_annotation++;  // vertex does not have annotation
    }
  }

  if (vertices_wo_annotation)
  {
    printf("%d vertices did not have an annotation!\n",vertices_wo_annotation);
  }

  /* any labels imported? */
  if (num_labels > 0)
  {

    /* if we have our own color table, now is the time to send it to the
       tcl side of things. */
    if (mris->ct)
      labl_send_color_table_info ();

    free (done);

    /* show the label label in the interface. */
    send_tcl_command ("ShowLabel kLabel_Label 1");
    labl_draw_flag = 1;
    send_tcl_command ("UpdateLinkedVarGroup label");

  }
  else
  {
    printf ("surfer: WARNING: no labels imported; annotation was empty\n" );
  }

  return(ERROR_NONE);
}

int labl_export_annotation (char *fname)
{
  int vno;
  int label_index;
  int color;
  LABL_LABEL* label;
  int label_vno;

  if (NULL == fname)
    return (ERROR_BADPARM);

  for (vno = 0; vno < mris->nvertices; vno++)
    mris->vertices[vno].annotation = 0;

  /* for each label.. */
  for (label_index = 0; label_index < labl_num_labels; label_index++)
  {

    label = &(labl_labels[label_index]);

    /* if this is not a free label... */
    if ( LABL_TYPE_FREE != label->structure)
    {

      /* make the composed color int for this label. */
      MRISRGBToAnnot (label->r, label->g, label->b, color);

      /* for every vertex in the label... */
      for (label_vno = 0; label_vno < label->label->n_points; label_vno++)
      {
        if ( label->label->lv[label_vno].vno < 0 ||
             label->label->lv[label_vno].vno >= mris->nvertices )
          continue;

        /* set the annotation value. */
        mris->vertices[label->label->lv[label_vno].vno].annotation =
          color;
      }

      if (labl_debug)
      {
        printf( "saved label %d with %d vertices, color %d %d %d "
                "anot value %d\n", label_index, label->label->n_points,
                label->r, label->g, label->b, color );
      }
    }
  }

  /* write out the annotation. */
  MRISwriteAnnotation (mris, fname);

  return (ERROR_NONE);
}

int labl_new_from_marked_vertices (int *new_index_out)
{
  LABEL* label = NULL;
  int num_marked_verts;
  int vno;
  int label_vno;
  VERTEX* v = NULL;
  int new_index;
  char tcl_command[NAME_LENGTH + 50];
  float val = 0;

  /* count the number of marked vertices. */
  num_marked_verts = 0;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].marked)
      num_marked_verts++;

  /* if we didn't get any, return. */
  if (0 == num_marked_verts)
  {
    fprintf (stderr, "no marked vertices...\n");
    return (ERROR_NONE);
  }

  /* allocate a label. */
  label = LabelAlloc(num_marked_verts, NULL, NULL);
  strncpy( label->subject_name, pname, 100 );

  /* for every vertex, if it's marked, save its vertex coords,
     index. don't fill the value of the current overlay, as that
     should only be done when the label is actually written to
     file. */
  label_vno = 0;
  for (vno = 0; vno < mris->nvertices; vno++)
  {
    v = &mris->vertices[vno];
    if (v->marked)
    {
      sclv_get_value(v, sclv_current_field, &val) ;
      label->lv[label_vno].x = v->x;
      label->lv[label_vno].y = v->y;
      label->lv[label_vno].z = v->z;
      label->lv[label_vno].stat = val ;
      label->lv[label_vno].vno = vno;
      label_vno++;
    }
  }
  label->n_points = num_marked_verts;

  /* convert to original positions. */
  if (!origsurfloaded)
    read_orig_vertex_coordinates(orfname);
  LabelToWhite (label, mris);

  /* add this label to our list. */
  labl_add (label, &new_index);

  /* select this label */
  labl_select (new_index);

  /* clear the marked verts. */
  clear_all_vertex_marks ();

  surface_compiled = 0 ;

  /* return the new index if they want it. */
  if (NULL != new_index_out)
    *new_index_out = new_index;

  /* if the fill dlog is open, this will update it. */
  sprintf (tcl_command, "LabelsChanged");
  send_tcl_command (tcl_command);

  return (ERROR_NONE);
}

int labl_add_marked_vertices_to_label (int index)
{
  int* vnos_to_add;
  int* found_labels;
  int in_label;
  int num_found;
  int found_label;
  LABEL* curlabel = NULL;
  LABEL* newlabel = NULL;
  int num_marked_verts;
  int num_label_verts;
  int num_new_verts;
  int curlabel_vno;
  int newlabel_vno;
  int vno;
  VERTEX* v = NULL;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* Look at all the marked verts and find the ones that aren't
     already in this label. Make a list of them. */
  vnos_to_add = (int*) calloc (mris->nvertices, sizeof(int));
  if (NULL == vnos_to_add)
  {
    return (ERROR_NO_MEMORY);
  }
  found_labels = (int*) calloc (LABL_MAX_LABELS, sizeof(int));
  if (NULL == found_labels)
  {
    free (vnos_to_add);
    return (ERROR_NO_MEMORY);
  }
  num_marked_verts = 0;
  for (vno = 0; vno < mris->nvertices; vno++)
  {
    v = &mris->vertices[vno];
    if (v->marked)
    {
      labl_find_label_by_vno (vno, 0, found_labels,
                              sizeof(found_labels), &num_found);

      in_label = FALSE;
      if (num_found > 0)
        for (found_label = 0; found_label < num_found; found_label++)
          if (found_labels[found_label] == index)
            in_label = TRUE;

      if (!in_label)
        vnos_to_add[num_marked_verts++] = vno;

    }
  }
  free (found_labels);

  /* if we didn't get any, return. */
  if (0 == num_marked_verts)
  {
    free (vnos_to_add);
    printf ("surfer: no marked vertices (that aren't already in the label)");
    return (ERROR_NONE);
  }

  /* get the number of points in the label */
  curlabel = labl_labels[index].label;
  num_label_verts = curlabel->n_points;

  /* create a new label with the combined number of points */
  num_new_verts = num_marked_verts + num_label_verts;
  newlabel = LabelAlloc (num_new_verts, NULL, NULL);
  if (NULL == newlabel)
  {
    free (vnos_to_add);
    return (ERROR_NO_MEMORY);
  }

  strncpy( newlabel->subject_name, pname, 100 );

  /* add the label's points (note we're going through both the
     curlabel and newlabel vnos in this loop */
  newlabel_vno = 0;
  for (curlabel_vno = 0; curlabel_vno < curlabel->n_points; curlabel_vno++)
  {
    newlabel->lv[newlabel_vno].x = curlabel->lv[curlabel_vno].x;
    newlabel->lv[newlabel_vno].y = curlabel->lv[curlabel_vno].y;
    newlabel->lv[newlabel_vno].z = curlabel->lv[curlabel_vno].z;
    newlabel->lv[newlabel_vno].vno = curlabel->lv[curlabel_vno].vno;
    newlabel_vno++;
  }

  /* add the marked verts */
  for (vno = 0; vno < num_marked_verts; vno++)
  {
    v = &mris->vertices[vnos_to_add[vno]];
    newlabel->lv[newlabel_vno].x = v->x;
    newlabel->lv[newlabel_vno].y = v->y;
    newlabel->lv[newlabel_vno].z = v->z;
    newlabel->lv[newlabel_vno].vno = vnos_to_add[vno];
    newlabel_vno++;
  }
  newlabel->n_points = num_new_verts;

  /* delete the old label */
  LabelFree (&labl_labels[index].label);

  /* point to the newlabel in the labl structure */
  labl_labels[index].label = newlabel;

  /* update this label. */
  labl_changed (index, FALSE);

  free (vnos_to_add);

  return (ERROR_NONE);
}

int labl_remove_marked_vertices_from_label (int index)
{

  LABEL* curlabel = NULL;
  LABEL* newlabel = NULL;
  int num_marked_verts;
  int num_label_verts;
  int num_new_verts;
  int curlabel_vno;
  int newlabel_vno;
  int vno;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* count the number of marked vertices. */
  num_marked_verts = 0;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
    if (mris->vertices[vno].marked)
      num_marked_verts++;

  /* if we didn't get any, return. */
  if (0 == num_marked_verts)
  {
    printf ("surfer: no marked vertices\n");
    return (ERROR_NONE);
  }

  /* get the number of points in the label */
  curlabel = labl_labels[index].label;
  num_label_verts = curlabel->n_points;

  /* calc the number of points in the new label. start with the same
     num as the label. for each point in the label, see if it's marked
     on the surface. if so, dec the count. */
  num_new_verts = num_label_verts;
  for (curlabel_vno = 0; curlabel_vno < curlabel->n_points; curlabel_vno++)
  {
    vno = curlabel->lv[curlabel_vno].vno;
    if ( vno < 0 || vno >= mris->nvertices )
      continue;

    if (mris->vertices[vno].marked)
    {
      num_new_verts--;
    }
  }

  /* if the count is the same as the original label, there were no
     intersecting marked and label points. */
  if (num_new_verts == num_label_verts)
  {
    printf ("surfer: no intersection with label\n");
    return (ERROR_NONE);
  }

  /* create a new label with the difference of points */
  newlabel = LabelAlloc (num_new_verts, NULL, NULL);
  if (NULL == newlabel)
  {
    return (ERROR_NO_MEMORY);
  }

  strncpy (newlabel->subject_name, pname, 100);

  /* for each of the labels verts, add it to the new label if the
     corresponding surface vert is not marked. */
  newlabel_vno = 0;
  for (curlabel_vno = 0; curlabel_vno < curlabel->n_points; curlabel_vno++)
  {
    vno = curlabel->lv[curlabel_vno].vno;
    if ( vno < 0 || vno >= mris->nvertices )
      continue;

    if (!mris->vertices[vno].marked && newlabel_vno < num_new_verts)
    {
      newlabel->lv[newlabel_vno].x = curlabel->lv[curlabel_vno].x;
      newlabel->lv[newlabel_vno].y = curlabel->lv[curlabel_vno].y;
      newlabel->lv[newlabel_vno].z = curlabel->lv[curlabel_vno].z;
      newlabel->lv[newlabel_vno].vno = curlabel->lv[curlabel_vno].vno;
      newlabel_vno++;
    }
  }
  newlabel->n_points = num_new_verts;

  /* delete the old label */
  LabelFree (&labl_labels[index].label);

  /* point to the newlabel in the labl structure */
  labl_labels[index].label = newlabel;

  /* recalc stuff */
  labl_changed (index, TRUE);

  return (ERROR_NONE);
}

int labl_mark_vertices (int index)
{
  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* mark the verts. */
  LabelMark (labl_labels[index].label, mris);

  return (ERROR_NONE);
}

int labl_select (int index)
{
  char tcl_command[NAME_LENGTH + 50];
  int old_selected;
  /* mark this label as selected. */
  old_selected = labl_selected_label;
  labl_selected_label = index;

  /* if something changed... */
  if (old_selected != labl_selected_label)
  {
    /* if something was selected, send the select cpmmand and the
       update command to the tcl label list with this label's
       information. */
    if ((index >= 0 && index < labl_num_labels) && g_interp)
    {
      sprintf (tcl_command, "LblLst_SelectLabel %d", index);
      send_tcl_command (tcl_command);

      labl_send_info (index);
    }

    /* redraw. */
    redraw ();
  }

  return (ERROR_NONE);
}

int labl_set_name_from_table (int index)
{
  char name[NAME_LENGTH];
  LABL_LABEL* label;
  COLOR_TABLE* ctab;
  int ctab_err;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  if (NULL == labl_ctab)
    return (ERROR_NONE);

  label = &(labl_labels[index]);

  /* if the surface has a color table, use that to get the name,
     otherwise use the external file. */
  if (mris->ct)
    ctab = mris->ct;
  else
    ctab = labl_ctab;

  ctab_err = CTABcopyName (ctab, label->structure, name, sizeof(name));
  if (NO_ERROR == ctab_err)
    labl_set_info (index, name, label->structure, label->visible,
                   label->r, label->g, label->b);

  return (ERROR_NONE);
}

int labl_set_info (int index, char* name, int structure, int visible,
                   int ir, int ig, int ib)
{
  int ctab_err;
  COLOR_TABLE* ctab;
  int r, g, b, a;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);
  if (NULL == name)
    return (ERROR_BADPARM);
  if (0 != visible && 1 != visible)
    return (ERROR_BADPARM);

  /* set the info in this label. */
  strncpy (labl_labels[index].name, name, NAME_LENGTH);
  labl_labels[index].structure = structure;
  labl_labels[index].visible = visible;

  /* if we have a table (in mris or external), and the structure is
     not free, get the color from the table. otherwise, use the color
     that they gave us. */
  if ((labl_ctab || mris->ct) &&
      LABL_TYPE_FREE != labl_labels[index].structure)
  {
    if (mris->ct)
      ctab = mris->ct;
    else
      ctab = labl_ctab;

    ctab_err = CTABrgbaAtIndexi (ctab, structure, &r, &g, &b, &a);
    if (NO_ERROR == ctab_err)
    {
      labl_set_color (index, r, g, b);
      labl_set_alpha (index, a);
    }
  }
  else
  {
    labl_set_color (index, ir, ig, ib);
    labl_set_alpha (index, 255);
  }

  /* send the label info to tcl */
  labl_send_info (index);

  return (ERROR_NONE);
}

int labl_set_color (int index, int r, int g, int b)
{
  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);
  if (r < 0 || r > 255 || g < 0 || g > 255 || b < 0 || b > 255)
    return (ERROR_BADPARM);

  /* set the color */
  labl_labels[index].r = r;
  labl_labels[index].g = g;
  labl_labels[index].b = b;

  /* send the label info to tcl */
  labl_send_info (index);

  return (ERROR_NONE);
}

int labl_set_alpha (int index, int alpha)
{
  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);
  if (alpha < 0 || alpha > 255)
    return (ERROR_BADPARM);

  /* set the color */
  labl_labels[index].a = alpha;

  /* send the label info to tcl */
  labl_send_info (index);

  return (ERROR_NONE);
}

int labl_send_info (int index)
{
  char tcl_command[NAME_LENGTH + 50];

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* send the update command to the tcl label list with this label's
     information. */
  sprintf (tcl_command, "LblLst_UpdateInfo %d \"%s\" %d %d %d %d %d",
           index, labl_labels[index].name,
           labl_labels[index].structure, labl_labels[index].visible,
           labl_labels[index].r, labl_labels[index].g,
           labl_labels[index].b );
  send_tcl_command (tcl_command);

  /* if the fill dlog is open, this will update it. */
  sprintf (tcl_command, "LabelsChanged");
  send_tcl_command (tcl_command);

  return (ERROR_NONE);
}

int labl_add (LABEL* label, int* new_index)
{
  int index;
  char tcl_command[NAME_LENGTH + 50];

  if (NULL == label)
    return (ERROR_BADPARM);

  /* make sure we can add one. */
  if (labl_num_labels >= LABL_MAX_LABELS)
  {
    printf("max labels exceeded!\n") ;
    return (ERROR_NO_MEMORY);
  }

  /* add a new node at the end of the list. */
  index = labl_num_labels;
  labl_num_labels++;

  /* copy the label in and init the other data. make it a free label
     by default, with the default color. */
  labl_labels[index].label = label;
  labl_labels[index].structure = LABL_TYPE_FREE;
  labl_labels[index].r = LABL_DEFAULT_COLOR_R;
  labl_labels[index].g = LABL_DEFAULT_COLOR_G;
  labl_labels[index].b = LABL_DEFAULT_COLOR_B;
  labl_labels[index].a = LABL_DEFAULT_COLOR_A;
  labl_labels[index].visible = 1;
  labl_labels[index].border_vno = NULL;
  labl_labels[index].num_border_vnos = 0;

  /* set the name using a global counter and increment the counter. */
  labl_num_labels_created++;
  sprintf (labl_labels[index].name, "Label %d", labl_num_labels_created);

  /* calc the extent and stuff. */
  labl_changed (index, FALSE);

  /* notify tcl of the new label. */
  sprintf (tcl_command, "LblLst_AddLabel \"%s\"", labl_labels[index].name);
  send_tcl_command (tcl_command);

  /* if the fill dlog is open, this will update it. */
  sprintf (tcl_command, "LabelsChanged");
  send_tcl_command (tcl_command);

  /* enable our label menu items. */
  enable_menu_set (MENUSET_LABEL_LOADED, 1);

  /* if they want the new index, return it. */
  if (NULL != new_index)
    *new_index = index;

  return (ERROR_NONE);
}

int labl_changed (int index, int vertices_were_removed)
{
  LABEL* label;
  float min_x, max_x;
  float min_y, max_y;
  float min_z, max_z;
  int label_vno;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  label = labl_labels[index].label;

  /* go through the label and find the bounds. */
  min_x = min_y = min_z = 500;
  max_x = max_y = max_z = -500;
  for (label_vno = 0; label_vno < label->n_points; label_vno++)
  {
    if (label->lv[label_vno].x < min_x)
      min_x = label->lv[label_vno].x;
    if (label->lv[label_vno].y < min_y)
      min_y = label->lv[label_vno].y;
    if (label->lv[label_vno].z < min_z)
      min_z = label->lv[label_vno].z;
    if (label->lv[label_vno].x > max_x)
      max_x = label->lv[label_vno].x;
    if (label->lv[label_vno].y > max_y)
      max_y = label->lv[label_vno].y;
    if (label->lv[label_vno].z > max_z)
      max_z = label->lv[label_vno].z;
  }

  /* set the bounds in the label, modifying it by the fudge value. */
  labl_labels[index].min_x = min_x - LABL_FUDGE;
  labl_labels[index].min_y = min_y - LABL_FUDGE;
  labl_labels[index].min_z = min_z - LABL_FUDGE;
  labl_labels[index].max_x = max_x + LABL_FUDGE;
  labl_labels[index].max_y = max_y + LABL_FUDGE;
  labl_labels[index].max_z = max_z + LABL_FUDGE;

  /* this label needs to be recached */
  labl_labels[index].cached = FALSE;
  labl_cache_updated = FALSE;

  /* if vertices were removed, we need to force rebuild the entire
     cache. */
  if (vertices_were_removed)
    labl_update_cache (TRUE);

  /* find the border. */
  labl_find_and_set_border (index);

  return (ERROR_NONE);
}

int labl_all_changed ()
{

  int label;

  for ( label = 0; label < labl_num_labels; label++ )
  {
    labl_changed( label, 0 );
  }

  return (ERROR_NONE);
}

int labl_remove (int index)
{
  int next;
  char tcl_command[NAME_LENGTH];

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* free the label here. */
  LabelFree (&labl_labels[index].label);

  /* Free the border storage. */
  free (labl_labels[index].border_vno);

  /* for every label above it, copy it one slot down. */
  next = index + 1;
  while (next < labl_num_labels)
  {
    labl_labels[next-1] = labl_labels[next];
    next++;
  }

  /* decrement the number of labels. */
  labl_num_labels--;

  /* rebuild the cache index. */
  labl_update_cache (TRUE);

  /* notify tcl that this label is gone. */
  sprintf (tcl_command, "LblLst_RemoveLabel %d", index );
  send_tcl_command (tcl_command);

  /* if the fill dlog is open, this will update it. */
  sprintf (tcl_command, "LabelsChanged");
  send_tcl_command (tcl_command);

  /* if this was our selected label, select nothing. */
  if (labl_selected_label == index)
    labl_select (-1);

  /* if our selection was above label selected, decrement it. */
  else if (labl_selected_label > index)
    labl_selected_label--;

  /* find all the borders again. */
  labl_find_and_set_all_borders ();

  return (ERROR_NONE);
}

int labl_mark_threshold_vertices (int index, float threshold)
{
  int    n, vno ;
  VERTEX* v;
  LABEL* area = labl_labels[index].label ;

  for (n = 0 ; n < area->n_points ; n++)
  {
    /* for the nth vertex in the label */
    vno = area->lv[n].vno ;
    if (vno < 0 || vno >= mris->nvertices)
      return (ERROR_BADPARM) ;

    /* mark the vertex if it's less than threshold */
    if ( area->lv[n].stat < threshold )
    {
      v = &mris->vertices[vno] ;
      v->marked = 1 ;
    }
  }
  return (ERROR_NONE) ;
}

int labl_threshold (int index, float threshold)
{
  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* mark all the vertices below the threshold */
  labl_mark_threshold_vertices(index, threshold);

  /* remove the marked vertices from the label */
  labl_remove_marked_vertices_from_label(index);

  labl_update_cache (TRUE);
  
  labl_find_and_set_border (index);
  
  /* clear the marked vertices */
  clear_all_vertex_marks();
  
  return (ERROR_NONE);
}

int labl_remove_all ()
{

  int labels;

  /* delete all the labels. */
  labels = labl_num_labels;
  while (labels > 0)
  {
    labl_remove (0);
    labels--;
  }
  labl_num_labels_created = 0;
  labl_num_labels = 0 ;

  return (ERROR_NONE);
}

int labl_find_label_by_vno (int vno, int min_label, int* index_array,
                            int array_size, int* out_num_found)
{
  int num_to_copy, num_to_copy_from_beginning;
  int i;

  /* if they passed -1 for min_label, they really want 0. (ugh) */
  if (min_label < 0)
    min_label = 0;

  /* update the cache if necessary */
  if (!labl_cache_updated)
    labl_update_cache (TRUE);

  /* return some points. */
  num_to_copy = MIN(array_size,labl_num_labels_at_cache_vno[vno]);
  num_to_copy_from_beginning = 0;
  if (num_to_copy > 0)
  {
    if (min_label != 0)
    {
      /* find the first label index <= min_label */
      i = 0;
      while (labl_cache[vno][i] < min_label &&
             i < labl_num_labels_at_cache_vno[vno])
        i++;

      num_to_copy = labl_num_labels_at_cache_vno[vno] - i;
      if (num_to_copy > array_size) num_to_copy = array_size;

      num_to_copy_from_beginning = i;
      if (num_to_copy_from_beginning > array_size - num_to_copy)
        num_to_copy_from_beginning = array_size - num_to_copy;
      if ( num_to_copy_from_beginning < 0 )
        num_to_copy_from_beginning = 0;


      memmove (index_array, &(labl_cache[vno][i]),
              num_to_copy * sizeof(int));

      memmove (&(index_array[num_to_copy]), labl_cache[vno],
              num_to_copy_from_beginning * sizeof(int));
    }
    else
    {
      memmove (index_array, labl_cache[vno], num_to_copy * sizeof(int));
    }
  }

  *out_num_found = num_to_copy + num_to_copy_from_beginning;

  return (ERROR_NONE);
}

int labl_select_label_by_vno (int vno)
{
  int label_index = -1;
  int num_found;

  if (vno < 0 || vno >= mris->nvertices)
    return (ERROR_BADPARM);

  /* try and find a label. if found, select it. */
  labl_debug = 1;
  labl_find_label_by_vno (vno, labl_selected_label+1,
                          &label_index, 1, &num_found);
  labl_debug = 0;
  if (num_found > 0)
    labl_select (label_index);
  else
    labl_select (-1);

  return (ERROR_NONE);
}

int labl_apply_color_to_vertex (int vno, GLubyte* r, GLubyte* g, GLubyte* b )
{
  int label_index_array[LABL_MAX_LABELS];
  int num_labels_found, found_label_index;
  float br, bg, bb;
  float lr, lg, lb;
  int label_index;

  if (vno < 0 || vno >= mris->nvertices)
    return (ERROR_BADPARM);

  /* if our display flag is off, do nothing. */
  if ( !labl_draw_flag )
  {
    return (ERROR_NONE);
  }

  /* try and find a label. if found... */
  labl_find_label_by_vno (vno, 0, label_index_array,
                          LABL_MAX_LABELS, &num_labels_found);
  if (num_labels_found > 0)
  {
    for (found_label_index = 0; found_label_index < num_labels_found;
         found_label_index++)
    {
      label_index = label_index_array[found_label_index];

      /* if this is the selected label and this is a border of width 1
         or 2, make it our outline color. */
      if (labl_selected_label == label_index &&
          labl_vno_is_border(labl_selected_label, vno))
      {
        *r = labl_outline_color[0];
        *g = labl_outline_color[1];
        *b = labl_outline_color[2];
      }
      /* else if this label is visible... */
      else if (labl_labels[label_index].visible)
      {
        /* color it in the given drawing style. */
        switch (labl_draw_style)
        {
        case LABL_STYLE_FILLED:
          /* If this is filled, we're going to blend the
             background with the label color with the alpha
             level from the color table. */
          br = ((float)(*r) / 255.0) *
               (1.0 - ((float)labl_labels[label_index].a / 255.0));
          bg = ((float)(*g) / 255.0) *
               (1.0 - ((float)labl_labels[label_index].a / 255.0));
          bb = ((float)(*b) / 255.0) *
               (1.0 - ((float)labl_labels[label_index].a / 255.0));

          lr = ((float)labl_labels[label_index].r / 255.0) *
               ((float)labl_labels[label_index].a / 255.0);
          lg = ((float)labl_labels[label_index].g / 255.0) *
               ((float)labl_labels[label_index].a / 255.0);
          lb = ((float)labl_labels[label_index].b / 255.0) *
               ((float)labl_labels[label_index].a / 255.0);

          *r = (GLubyte)((br + lr) * 255.0);
          *g = (GLubyte)((bg + lg) * 255.0);
          *b = (GLubyte)((bb + lb) * 255.0);
          break;
        case LABL_STYLE_OUTLINE:
          /* if this is a border of width 1, color it the color of the
             label. */
          if (labl_vno_is_border (label_index, vno))
          {
            *r = labl_labels[label_index].r;
            *g = labl_labels[label_index].g;
            *b = labl_labels[label_index].b;
          }
          break;
        default:
          ;
        }
      }
    }
  }

  return (ERROR_NONE);
}

int labl_erode (int index)
{
  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* Perform the erode. */
  LabelErode (labl_labels[index].label, mris, 1);

  /* Label is changed and we removed points. */
  labl_changed (index, TRUE);

  return (ERROR_NONE);
}

int labl_dilate (int index)
{
  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* Dilate the label. */
  LabelDilate (labl_labels[index].label, mris, 1, CURRENT_VERTICES);

  /* Label is changed but we didn't remove points. */
  labl_changed (index, FALSE);

  return (ERROR_NONE);
}

int labl_fill_holes (int index)
{
  LABEL* filled;

  if (index < 0 || index >= labl_num_labels)
    return (ERROR_BADPARM);

  /* Fill the holes in the label. */
#if 1
  filled = LabelFillHoles (labl_labels[index].label, mris, WHITE_VERTICES);
#else
  filled = LabelFillHolesWithOrig (labl_labels[index].label, mris);
#endif
  if (NULL != filled)
  {
    /* Save the new label. */
    LabelFree (&labl_labels[index].label);
    labl_labels[index].label = filled;

    /* Label is changed. */
    labl_changed (index, TRUE);
  }

  /* This marks verts as a byproduct, so unmark them. */
  clear_all_vertex_marks();

  return (ERROR_NONE);
}

int labl_print_list ()
{
  int label_index;
  LABL_LABEL* label;

  printf ("Num labels: %d (%d) Num strucutres: %d\n",
          labl_num_labels, labl_selected_label, labl_num_structures );
  for (label_index = 0; label_index < labl_num_labels; label_index++)
  {
    label = &(labl_labels[label_index]);

    printf ("Label %d\n", label_index );
    printf ("\tName: %s, %d vertices\n",
            label->name, label->label->n_points );
    printf ("\tStructure: %d Color: %d %d %d Visible: %d\n",
            label->structure, label->r, label->g, label->b, label->visible );
    printf ("\tBounds: x %f %f y %f %f z %f %f\n",
            label->min_x, label->max_x,
            label->min_y, label->max_y,
            label->min_z, label->max_z);
  }

  return (ERROR_NONE);
}

int labl_print_table ()
{
  int structure_index;
  COLOR_TABLE* ctab;
  int r, g, b, a;
  char structure_label[1024];

  printf( "Num strucutres: %d\n", labl_num_structures );
  for (structure_index = 0; structure_index < labl_num_structures;
       structure_index++)
  {
    if (mris->ct)
      ctab = mris->ct;
    else
      ctab = labl_ctab;

    CTABrgbaAtIndexi (ctab, structure_index, &r, &g, &b, &a);
    CTABcopyName (ctab, structure_index,
                  structure_label, sizeof(structure_label));

    printf ("Structure %d: %s r %d g %d b %d a %d",
            structure_index, structure_label, r, g, b, a);
  }

  return (ERROR_NONE);
}

/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------ fill paths */
static int
path_set_marks(PATH_PATH *b, MRI_SURFACE *mris, int mark)
{
  int i, vno ;

  for (i = 0 ; i < b->num_vertices ; i++)
  {
    vno = b->vertices[i] ;
    if (vno < 0 || vno >= mris->nvertices)
      continue ;
    mris->vertices[vno].marked = mark ;
  }
  return(NO_ERROR) ;

}

PATH_PATH *
path_copy(PATH_PATH *bsrc, PATH_PATH *bdst)
{
  int   vno ;

  if (!bdst)
  {
    bdst = (PATH_PATH *)calloc(1, sizeof(PATH_PATH)) ;
    bdst->vertices = (int *)calloc(bsrc->num_vertices, sizeof(int)) ;
    if (!bdst->vertices)
      ErrorExit(ERROR_NOMEMORY, "%s: could not copy path with %d vertices",
                Progname, bdst->num_vertices) ;
  }

  bdst->num_vertices = bsrc->num_vertices ;
  bdst->min_x = bsrc->min_x ;
  bdst->max_x = bsrc->max_x ;
  bdst->min_y = bsrc->min_y ;
  bdst->max_y = bsrc->max_y ;
  bdst->min_z = bsrc->min_z ;
  bdst->max_z = bsrc->max_z ;
  for (vno = 0 ; vno < bdst->num_vertices ; vno++)
    bdst->vertices[vno] = bsrc->vertices[vno] ;
  return(bdst) ;
}

int path_initialize ()
{
  int path;

  /* Clear the path array. */
  for (path = 0; path < PATH_MAX_PATHS; path++)
  {
    path_paths[path].num_vertices = 0;
    path_paths[path].vertices = NULL;
  }

  /* No paths present. */
  path_selected_path = PATH_NONE_SELECTED;
  path_num_paths = 0;

  /* This will be itialized in path_update_surface_paths. */
  if (NULL != path_is_path)
    free (path_is_path);
  path_is_path = NULL;

  return (ERROR_NONE);
}

int path_new_path_from_marked_vertices ()
{
  int* path;
  int path_length;

  /* Initialize a new path. */
  path = (int*) calloc (mris->nvertices, sizeof(int));

  /* Find the actual vnos in this path between the currently marked
     vertices. */
  find_path (marked, nmarked, "making path", mris->nvertices,
             path, &path_length);

  /* Add the vertices we found to a new path. */
  path_add (path_length, path, NULL);

  /* Free the storage for our path finder. */
  free (path);

  return (ERROR_NONE);
}

int path_remove_selected_path ()
{
  /* Call path_remove on currently selected path. */
  if (path_selected_path >= 0 &&
      path_selected_path < path_num_paths)
    path_remove (path_selected_path);

  return (ERROR_NONE);
}

int path_add (int num_vertices, int* vertices, int* new_index)
{
  int index;
  float min_x, max_x;
  float min_y, max_y;
  float min_z, max_z;
  int path;
  VERTEX* v;

  /* Make sure we have vertices. */
  if (num_vertices <= 0)
    return (ERROR_BADPARM);
  if (NULL == vertices)
    return (ERROR_BADPARM);

  /* Make sure we can add one. */
  if (path_num_paths >= PATH_MAX_PATHS)
    return (ERROR_NO_MEMORY);

  /* Add a new one at the end of our array of paths. */
  index = path_num_paths;
  path_num_paths++;

  /* Allocate the vertices array and copy everything in. */
  path_paths[index].num_vertices = num_vertices;
  path_paths[index].vertices = (int*) calloc (num_vertices, sizeof(int));
  if (NULL == path_paths[index].vertices)
  {
    printf ("path_add: calloc failed with %d elements\n", num_vertices);
    return (ERROR_NO_MEMORY);
  }
  memmove (path_paths[index].vertices, vertices, num_vertices * sizeof(int));

  /* Go through the path and find the bounds of the vertices. */
  min_x = min_y = min_z = 500;
  max_x = max_y = max_z = -500;
  for (path = 0;
       path < path_paths[index].num_vertices;
       path++)
  {
    v = &(mris->vertices[ path_paths[index].vertices[path] ]);

    if (v->x < min_x)
      min_x = v->x;
    if (v->y < min_y)
      min_y = v->y;
    if (v->z < min_z)
      min_z = v->z;
    if (v->x > max_x)
      max_x = v->x;
    if (v->y > max_y)
      max_y = v->y;
    if (v->z > max_z)
      max_z = v->z;
  }

  /* Set the bounds in the path, modifying it by the fudge value. */
  path_paths[index].min_x = min_x - PATH_FUDGE;
  path_paths[index].min_y = min_y - PATH_FUDGE;
  path_paths[index].min_z = min_z - PATH_FUDGE;
  path_paths[index].max_x = max_x + PATH_FUDGE;
  path_paths[index].max_y = max_y + PATH_FUDGE;
  path_paths[index].max_z = max_z + PATH_FUDGE;

  /* Update the path flags. This will mark all vertices that are in
     paths in our path_is_path array. */
  path_update_surface_paths ();

  /* Return the new index if they want it. */
  if (NULL != new_index)
    *new_index = index;

  return (ERROR_NONE);
}

int path_remove (int index)
{
  int next;

  if (index < 0 || index >= path_num_paths)
    return (ERROR_BADPARM);

  /* Free the vertices array. */
  free (path_paths[index].vertices);

  /* Bump everything above this index down one. */
  next = index + 1;
  while (next < path_num_paths)
  {
    path_paths[next-1] = path_paths[next];
    next++;
  }

  /* Decrement the number of paths. */
  path_num_paths--;

  /* Update the paths flags */
  path_update_surface_paths ();

  /* If this was our selected path, select nothing. */
  if (path_selected_path == index)
    path_select (-1);

  /* If our selection was above path selected, decrement it so the
     same path is still selected. */
  else if (path_selected_path > index)
    path_selected_path--;

  return (ERROR_NONE);
}

int path_update_bounding_boxes ()
{

  VERTEX* v;
  int path, path_vno;
  float min_x, min_y, min_z;
  float max_x, max_y, max_z;

  /* For every path... */
  for (path = 0; path < path_num_paths; path++)
  {

    /* Go through the path and find the bounds of the vertices. */
    min_x = min_y = min_z = 500;
    max_x = max_y = max_z = -500;
    for (path_vno = 0;
         path_vno < path_paths[path].num_vertices;
         path_vno++ )
    {

      v = &(mris->vertices[ path_paths[path].vertices[path_vno] ]);

      if (v->x < min_x)
        min_x = v->x;
      if (v->y < min_y)
        min_y = v->y;
      if (v->z < min_z)
        min_z = v->z;
      if (v->x > max_x)
        max_x = v->x;
      if (v->y > max_y)
        max_y = v->y;
      if (v->z > max_z)
        max_z = v->z;
    }

    /* Set the bounds in the path, modifying it by the fudge value. */
    path_paths[path].min_x = min_x - PATH_FUDGE;
    path_paths[path].min_y = min_y - PATH_FUDGE;
    path_paths[path].min_z = min_z - PATH_FUDGE;
    path_paths[path].max_x = max_x + PATH_FUDGE;
    path_paths[path].max_y = max_y + PATH_FUDGE;
    path_paths[path].max_z = max_z + PATH_FUDGE;

  }

  return (ERROR_NONE);
}

int path_select (int index)
{
  int old_selected;

  /* Select this path. */
  old_selected = path_selected_path;
  path_selected_path = index;

  /* If this changes which path is selected, redraw. */
  if (old_selected != path_selected_path)
    redraw();

  return (ERROR_NONE);
}

int path_mark_selected_path ()
{

  /* Call path_mark on currently selected path. */
  path_mark (path_selected_path);

  return (ERROR_NONE);
}

int path_mark (int index)
{

  int path_vno;

  if (path_selected_path < 0 ||
      path_selected_path >= path_num_paths)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM, "path_mark: path index not valid"));

  /* Mark all vertices in the path. */
  for (path_vno = 0;
       path_vno < path_paths[index].num_vertices;
       path_vno++ )
  {
    mark_vertex (path_paths[index].vertices[path_vno], TRUE);
  }

  return (ERROR_NONE);
}

int path_update_surface_paths ()
{
  int path;
  int path_vno;

  /* Make our path flag array if we haven't already. It's just the
     size of nvertices. All are initialized to FALSE. */
  if (NULL == path_is_path)
  {
    path_is_path = (char*) calloc (mris->nvertices, sizeof(char));
    if (NULL == path_is_path)
      return (ERROR_NO_MEMORY);
  }
  else
  {
    /* Erase the current array. */
    bzero (path_is_path, mris->nvertices*sizeof(char) );
  }

  /* For every path... */
  for (path = 0; path < path_num_paths; path++)
  {
    /* For every vertex, mark that vertex in the array as a path. */
    for (path_vno = 0;
         path_vno < path_paths[path].num_vertices;
         path_vno++ )
    {
      path_is_path[path_paths[path].vertices[path_vno]] = TRUE;
    }
  }

  return (ERROR_NONE);
}

char path_is_vertex_on_path (int vno)
{
  if (vno < 0 || vno >= mris->nvertices)
    return (FALSE);
  if (NULL == path_is_path)
    return (FALSE);

  /* Return the value of the path flag here. */
  return (path_is_path[vno]);
}

int path_select_path_by_vno (int vno)
{
  int path;
  int path_vno;
  float x, y, z;
  VERTEX *u;
  VERTEX *v;
  float dist_uv;
  int smallest_distance;
  int closest_path=-1;

  if (vno < 0 || vno >= mris->nvertices)
    return (FALSE);
  if (NULL == path_is_path)
    return (FALSE);

  /* Get this vertex. */
  v = &(mris->vertices[vno]);

  /* Get the RAS coords. */
  x = v->x;
  y = v->y;
  z = v->z;

  /* For each path.. */
  smallest_distance = 1000;
  path = -1;
  for (path = 0; path < path_num_paths; path++)
  {
    /* Preflight test sees if this point is in the bounding cube of
       the path. */
    if (x >= path_paths[path].min_x &&
        x <= path_paths[path].max_x &&
        y >= path_paths[path].min_y &&
        y <= path_paths[path].max_y &&
        z >= path_paths[path].min_z &&
        z <= path_paths[path].max_z)
    {
      /* For each vertex in the path... */
      for (path_vno = 0;
           path_vno < path_paths[path].num_vertices;
           path_vno++)
      {
        /* Get the vertex. */
        u = &(mris->vertices[path_paths[path].vertices[path_vno]]);

        /* Find the distance between our input vertex and this one. */
        dist_uv = ((v->x - u->x) * (v->x - u->x)) +
                  ((v->y - u->y) * (v->y - u->y)) +
                  ((v->z - u->z) * (v->z - u->z));

        /* If any is less than the select distance and closer
           than the closest path so far, select this path. */
        if (dist_uv < PATH_DISTANCE_TO_SELECT &&
            dist_uv < smallest_distance)
        {
          smallest_distance = dist_uv;
          closest_path = path;
        }
      }
    }
  }

  /* If we didn't find anything, path will be -1 and this will
     deselect the current path. */
  path_select (closest_path);
  return (ERROR_NONE);
}

int path_apply_color_to_vertex (int vno, GLubyte* r, GLubyte* g, GLubyte* b )
{
  float x, y, z;
  char selected;
  int path_vno;

  if (vno < 0 || vno >= mris->nvertices)
    return (ERROR_BADPARM);
  if (NULL == r || NULL == g || NULL == b)
    return (ERROR_BADPARM);

  /* If this is a path.. */
  if (path_is_vertex_on_path(vno))
  {

    /* Get the location of this vertex. */
    x = mris->vertices[vno].x;
    y = mris->vertices[vno].y;
    z = mris->vertices[vno].z;

    /* See if it's in the selected bound cube, then in the selectd
       path vertex list. */
    selected = FALSE;
    if (x >= path_paths[path_selected_path].min_x &&
        x <= path_paths[path_selected_path].max_x &&
        y >= path_paths[path_selected_path].min_y &&
        y <= path_paths[path_selected_path].max_y &&
        z >= path_paths[path_selected_path].min_z &&
        z <= path_paths[path_selected_path].max_z)
      for (path_vno = 0;
           path_vno < path_paths[path_selected_path].num_vertices;
           path_vno++)
        if (vno ==
            path_paths[path_selected_path].vertices[path_vno])
          selected = TRUE;

    if (selected)
    {
      /* If it's selected, mark it yellow. */
      *r = 255;
      *g = 255;
      *b = 0;
    }
    else
    {
      /* Else just color it red. */
      *r = 255;
      *g = 0;
      *b = 0;
    }
  }

  return (ERROR_NONE);
}

int path_save (char* fname)
{
  PATH** paths;
  int    path_index;
  int    path_del_index;
  int    path_vno;
  int    vno;
  int    err = ERROR_NONE;

  /* We need to convert out vertex paths to PATH object. Make a list
     of them of the size of the number of paths we have. */
  paths = (PATH**) calloc (path_num_paths, sizeof(PATH));
  if (NULL == paths)
  {
    ErrorReturn (ERROR_NO_MEMORY,
                 (ERROR_NO_MEMORY, "Couldn't allocate %d paths",
                  path_num_paths));

  }

  for (path_index = 0; path_index < path_num_paths; path_index++)
  {

    /* Allocate a path. */
    paths[path_index] =
      PathAlloc (path_paths[path_index].num_vertices, "");
    if (NULL == paths[path_index])
    {
      for (path_del_index = 0; path_index < path_index; path_del_index++)
      {
        PathFree( &paths[path_del_index] );
      }
      free (paths);
      ErrorReturn (ERROR_NO_MEMORY,
                   (ERROR_NO_MEMORY, "Couldn't path with %d points",
                    path_paths[path_index].num_vertices));
    }

    /* Put the orig coordinates into the path we just created. */
    for (path_vno = 0;
         path_vno < path_paths[path_index].num_vertices; path_vno++)
    {
      vno = path_paths[path_index].vertices[path_vno];

      paths[path_index]->points[path_vno].x = mris->vertices[vno].origx;
      paths[path_index]->points[path_vno].y = mris->vertices[vno].origy;
      paths[path_index]->points[path_vno].z = mris->vertices[vno].origz;
      paths[path_index]->points[path_vno].vno = vno;
    }
  }

  /* Write the paths. */
  err = PathWriteMany (fname, path_num_paths, paths);
  if (err != ERROR_NONE)
  {
    printf ("ERROR: Couldn't write paths.\n");
    err = ERROR_NO_FILE;
  }

  /* Delete the stuff we allocated. */
  for (path_index = 0; path_index < path_num_paths; path_index++)
  {
    PathFree( &paths[path_index] );
  }
  free (paths);

  return err;
}

int path_load (char* fname)
{
  int    err;
  int    num_read;
  int*   vertices = NULL;
  int*   connected_verts = NULL;
  int    size_connected_path;
  PATH** paths    = NULL;
  int    path_index;
  PATH*  path     = NULL;
  LABEL* label    = NULL;
  int    pno;
  int    path_vno;

  /* Try to read the paths file. */
  err = PathReadMany (fname, &num_read, &paths);
  if (err != ERROR_NONE ||
      num_read <= 0 || paths == NULL)
  {
    ErrorReturn (ERROR_BADFILE,
                 (ERROR_BADFILE, "Couldn't read any paths."));
  }

  /* Allocate temp vno storage for our conversions. */
  vertices = (int*) calloc (mris->nvertices, sizeof(int));
  if (NULL == vertices)
  {
    ErrorReturn (ERROR_NO_MEMORY,
                 (ERROR_NO_MEMORY, "Couldn't allocate vertex storage."));
  }

  connected_verts = (int*) calloc (mris->nvertices, sizeof(int));
  if (NULL == connected_verts)
  {
    ErrorReturn (ERROR_NO_MEMORY,
                 (ERROR_NO_MEMORY, "Couldn't allocate vertex storage."));
  }

  for (path_index = 0; path_index < num_read; path_index++)
  {
    /* Get the path. */
    path = paths[path_index];

    /* We'll do our coordinate->vertex conversion using the
       LabelFillUnassignedVertices, but it only takes LABEL
       structures, so we'll fill one of those out with the
       path's coords. */
    label = LabelAlloc (path->n_points, NULL, NULL);
    if (NULL == label)
    {
      free (vertices);
      ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,
                                 "Couldn't allocate label for %d vertices\n",
                                 path->n_points));
    }
    label->n_points = path->n_points;

    /* Put the path coords in the label. */
    for (pno = 0; pno < path->n_points; pno++)
    {
      label->lv[pno].x = path->points[pno].x;
      label->lv[pno].y = path->points[pno].y;
      label->lv[pno].z = path->points[pno].z;
      label->lv[pno].vno = path->points[pno].vno;
    }

    /* This will find vertex numbers for all those points. */
    LabelFillUnassignedVertices (mris, label, WHITE_VERTICES);

    /* Copy the vertex numbers from the label into our vertices
       array. */
    for (path_vno = 0; path_vno < path->n_points; path_vno++)
    {
      vertices[path_vno] = label->lv[path_vno].vno;
    }
    LabelFree (&label);

    size_connected_path = 0;

    /* Make sure the path is connected. */
    find_path (vertices, path->n_points, "Connecting path",
               mris->nvertices, connected_verts, &size_connected_path);

    /* Make a new path from our num_vertices and vertices
       array. */
    path_add (size_connected_path, connected_verts, NULL);

    /* Free this path. */
    PathFree (&paths[path_index]);

  }

  free (vertices);
  free (connected_verts);
  free (paths);

  return (ERROR_NONE);
}

/* ---------------------------------------------------------------------- */

/* ------------------------------------------------------- floodfill mark */

int fill_flood_from_seed (int seed_vno, FILL_PARAMETERS* params)
{
  char* filled;
  int num_filled_this_iter;
  int num_filled;
  int iter;
  int min_vno, max_vno, step_vno;
  int vno;
  int this_label = 0;
  int neighbor_index;
  int neighbor_vno;
  VERTEX* v;
  VERTEX* neighbor_v;
  float fvalue = 0;
  float seed_curv = 0;
  float seed_fvalue = 0;
  int new_index;
  int label_index_array[LABL_MAX_LABELS];
  int num_labels_found, found_label_index;
  int skip;
  int count;

  if (seed_vno < 0 || seed_vno >= mris->nvertices)
    return (ERROR_BADPARM);
  if (NULL == params)
    return (ERROR_BADPARM);
  if (params->action < 0 || params->action >= NUM_FILL_ACTIONS)
    return (ERROR_BADPARM);

  /* init filled array. */
  filled = (char*) calloc (mris->nvertices, sizeof(char));

  /* start with the seed filled.*/
  filled[seed_vno] = TRUE;

  /* find seed values for some conditions. */
  if (params->dont_cross_label)
    this_label = labl_selected_label;
  if (params->dont_cross_cmid)
    seed_curv = mris->vertices[seed_vno].curv;
  if (params->dont_cross_fthresh)
    sclv_get_value (&mris->vertices[seed_vno],
                    sclv_current_field, &seed_fvalue);

  /* Start listing for cancel */
  printf ("surfer: filling (ctrl-c to cancel)");
  cncl_start_listening ();

  /* while we're still filling stuff in a pass... */
  num_filled_this_iter = 1;
  num_filled = 0;
  iter = 0;
  while (num_filled_this_iter > 0)
  {
    if (cncl_user_canceled())
    {
      goto cancel;
    }

    num_filled_this_iter = 0;

    /* switch between iterating forward and backwards. */
    if ((iter%2)==0)
    {
      min_vno = 0;
      max_vno = mris->nvertices-1;
      step_vno = 1;
    }
    else
    {
      min_vno = mris->nvertices-1;
      max_vno = 0;
      step_vno = -1;
    }

    /* for each vertex, if it's filled, check its neighbors. for the
       rules that are up-to-and-including, make the check on this
       vertex. for the rules that are up-to-and-not-including, check
       on the neighbor. */
    for (vno = min_vno; vno != max_vno; vno += step_vno)
    {
      if (filled[vno])
      {

        /* check the neighbors... */
        v = &mris->vertices[vno];

        /* if this vert is ripped, move on. */
        if (v->ripflag)
        {
          continue;
        }

        /* if we're not crossing paths, check if this is a
           path. if so, move on. */
        if (params->dont_cross_path &&
            path_is_vertex_on_path (vno))
        {
          continue;
        }

        /* if we're not crossing the cmid, see if the cmid at this
           vertex is on the other side of the cmid as the seed
           point. if so, move on. */
        if (params->dont_cross_cmid &&
            ((seed_curv <= cmid && v->curv > cmid) ||
             (seed_curv >= cmid && v->curv < cmid)))
        {
          continue;
        }

        for (neighbor_index = 0;
             neighbor_index < v->vnum;
             neighbor_index++)
        {
          neighbor_vno = v->v[neighbor_index];
          neighbor_v = &mris->vertices[neighbor_vno] ;

          /* if the neighbor is filled, move on. */
          if (filled[neighbor_vno])
            continue;

          /* if we're not crossing labels, check if the label at
             this vertex is the same as the one at the seed. if not,
             move on. */
          if (params->dont_cross_label ||
              params->dont_fill_unlabeled)
          {

            labl_find_label_by_vno (neighbor_vno, 0,
                                    label_index_array,
                                    LABL_MAX_LABELS,
                                    &num_labels_found);
            if (num_labels_found > 0 &&
                params->dont_cross_label)
            {
              skip = 0;
              for (found_label_index = 0;
                   found_label_index < num_labels_found;
                   found_label_index++)
              {
                if (label_index_array[found_label_index] !=
                    this_label)
                {
                  skip = 1;
                  break;
                }
              }
              if (skip) continue;
            }
            if (num_labels_found == 0 &&
                params->dont_fill_unlabeled)
            {
              continue;
            }
          }

          /* if we're not crossing the fthresh, make sure this
             point is above it, or, if our initial functional
             value was negative, make sure it's not above
             -fthresh. if not, move on. */
          if (params->dont_cross_fthresh)
          {
            sclv_get_value (neighbor_v, sclv_current_field, &fvalue);
            if ((fthresh != 0 &&
                 seed_fvalue > 0 &&
                 fvalue < fthresh) ||
                (fthresh != 0 &&
                 seed_fvalue < 0 &&
                 fvalue > -fthresh) ||
                (fthresh == 0 && (fvalue * seed_fvalue < 0)))
            {
              continue;
            }
          }

          /* mark this vertex as filled. */
          filled[neighbor_vno] = TRUE;
          num_filled_this_iter++;
          num_filled++;
        }
      }
    }

    iter++;

    if ((iter % 2) == 0)
    {
      printf (".");
      fflush (stdout);
    }
  }

  /* mark all filled vertices. */
  count = 0;
  for (vno = 0; vno < mris->nvertices; vno++ )
    if (filled[vno])
    {
      mris->vertices[vno].marked = TRUE;
      count++;
    }

  switch (params->action)
  {
  case FILL_ACTION_NEW_LABEL:
    /* make a new label from the marked vertices. */
    labl_new_from_marked_vertices (&new_index);
    params->new_label_index = new_index;
    break;
  case FILL_ACTION_ADD_LABEL:
    /* add the marked vertices to the label spec'd in the argument. */
    labl_add_marked_vertices_to_label (params->argument);
    break;
  case FILL_ACTION_REMOVE_LABEL:
    /* remove the marked vertices from the label spec'd in the argument. */
    labl_remove_marked_vertices_from_label (params->argument);
    break;
  default:
    break;
  }

  printf (" done, %d vertices filled\n", count);
  fflush (stdout);

  goto done;

cancel:
  printf (" canceled\n");
  fflush (stdout);

done:
  cncl_stop_listening ();

  free (filled);

  return (ERROR_NONE);
}
/* ---------------------------------------------------------------------- */

int edit_vertex_at_cursor ( int action, int argument )
{
  return edit_vertex (selection, action, argument);
}

int edit_vertex ( int vno, int action, int argument )
{
  int i;
  int* saved_marks;
  int saved_nmarked;
  int found_labels[LABL_MAX_LABELS];
  int num_found;
  int label;

  if (vno < 0 || vno > mris->nvertices)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "edit_vertex_at_cursor: vno is oob: %d)", vno));

  if (action < 0 || action >= NUM_VEDIT_ACTIONS)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "edit_vertex_at_cursor: action was invalid: %d)", action));

  /* Save the old marks. */
  saved_marks = (int*) calloc (mris->nvertices, sizeof(int));
  memmove (saved_marks, marked, sizeof(int) * mris->nvertices);
  saved_nmarked = nmarked;

  /* Clear marks. */
  clear_all_vertex_marks ();

  /* Mark this vertex. */
  mark_vertex (vno, TRUE);

  /* Work on stuff based on this one marked vertex. */
  switch (action)
  {
  case VEDIT_ACTION_NEW_LABEL:
    labl_new_from_marked_vertices (NULL);
    break;
  case VEDIT_ACTION_ADD_LABEL:
    labl_add_marked_vertices_to_label (argument);
    break;
  case VEDIT_ACTION_REMOVE_LABEL:
    labl_remove_marked_vertices_from_label (argument);
    break;
  case VEDIT_ACTION_CLEAR_LABELS:
    labl_find_label_by_vno (vno, 0, found_labels,
                            sizeof(found_labels), &num_found);
    for (label = 0; label < num_found; label++)
      labl_remove_marked_vertices_from_label (found_labels[label]);
    break;
  default:
    break;
  }

  /* Restore the marks. */
  memmove (marked, saved_marks, sizeof(int) * mris->nvertices);
  for (i = 0; i < mris->nvertices; i++)
    mris->vertices[i].marked = saved_marks[i];
  nmarked = saved_nmarked;

  return (ERROR_NONE);
}

/* ---------------------------------------------------------------------- */

int find_path ( int* vert_vno, int num_vno, char* message, int max_path_length,
                int* path, int* path_length )
{


  int cur_vert_vno;
  int src_vno;
  int dest_vno;
  int vno;
  char* check;
  float* dist;
  int* pred;
  char done;
  VERTEX* v;
  VERTEX* u;
  float closest_dist;
  int closest_vno;
  int neighbor;
  int neighbor_vno;
  float dist_uv;
  int path_vno;
  int num_path = 0;
  int num_checked;
  float vu_x, vu_y, vu_z;

  dist = (float*) calloc (mris->nvertices, sizeof(float));
  pred = (int*) calloc (mris->nvertices, sizeof(int));
  check = (char*) calloc (mris->nvertices, sizeof(char));
  num_path = 0;
  num_checked = 0;
  (*path_length) = 0;

  if (NULL != message)
  {
    fprintf (stdout, "surfer: %s (ctrl-c to cancel)", message);
  }
  cncl_start_listening ();
  for (cur_vert_vno = 0; cur_vert_vno < num_vno-1; cur_vert_vno++)
  {
    if (cncl_user_canceled())
    {
      goto cancel;
    }

    /* clear everything */
    for (vno = 0; vno < mris->nvertices; vno++)
    {
      dist[vno] = 999999;
      pred[vno] = -1;
      check[vno] = FALSE;
    }

    /* Set src and dest */
    src_vno = vert_vno[cur_vert_vno+1];
    dest_vno = vert_vno[cur_vert_vno];

    /* make sure both are in range. */
    if (src_vno < 0 || src_vno >= mris->nvertices ||
        dest_vno < 0 || dest_vno >= mris->nvertices)
      continue;

    if (src_vno == dest_vno)
      continue;

    /* pull the src vertex in. */
    dist[src_vno] = 0;
    pred[src_vno] = vno;
    check[src_vno] = TRUE;

    done = FALSE;
    while (!done)
    {
      if (cncl_user_canceled())
      {
        goto cancel;
      }

      /* find the vertex with the shortest edge. */
      closest_dist = 999999;
      closest_vno = -1;
      for (vno = 0; vno < mris->nvertices; vno++)
        if (check[vno])
          if (dist[vno] < closest_dist)
          {
            closest_dist = dist[vno];
            closest_vno = vno;
          }
      v = &(mris->vertices[closest_vno]);
      check[closest_vno] = FALSE;

      /* if this is the dest node, we're done. */
      if (closest_vno == dest_vno)
      {
        done = TRUE;
      }
      else
      {
        /* relax its neighbors. */
        for (neighbor = 0; neighbor < v->vnum; neighbor++)
        {
          neighbor_vno = v->v[neighbor];
          u = &(mris->vertices[neighbor_vno]);

          /* calc the vector from u to v. */
          vu_x = u->x - v->x;
          vu_y = u->y - v->y;
          if (flag2d)
            vu_z = 0;
          else
            vu_z = u->z - v->z;

          /* recalc the weight. */
          if (flag2d)
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
                           ((v->y - u->y) * (v->y - u->y)));
          else
            dist_uv = sqrt(((v->x - u->x) * (v->x - u->x)) +
                           ((v->y - u->y) * (v->y - u->y)) +
                           ((v->z - u->z) * (v->z - u->z)));

          /* if this is a new shortest path, update the predecessor,
             weight, and add it to the list of ones to check next. */
          if (dist_uv + dist[closest_vno] < dist[neighbor_vno])
          {
            pred[neighbor_vno] = closest_vno;
            dist[neighbor_vno] = dist_uv + dist[closest_vno];
            check[neighbor_vno] = TRUE;
          }
        }
      }
      num_checked++;
      if ((num_checked % 100) == 0)
      {
        printf (".");
        fflush (stdout);
      }
    }

    /* add the predecessors from the dest to the src to the path. */
    path_vno = dest_vno;
    path[(*path_length)++] = dest_vno;
    while (pred[path_vno] != src_vno &&
           (*path_length) < max_path_length )
    {
      path[(*path_length)++] = pred[path_vno];
      path_vno = pred[path_vno];
    }
  }
  printf (" done\n");
  fflush (stdout);

  goto done;

cancel:
  printf (" canceled\n");
  fflush (stdout);
  *path_length = 0;

done:
  cncl_stop_listening ();

  free (dist);
  free (pred);
  free (check);

  return (ERROR_NONE);
}

/* ---------------------------------------------------------------------- */

#define TCL_CACHE_INC 10
void send_tcl_command (char* cmd)
{
  char** new_cached_tcl_commands = NULL;
  char* new_cmd = NULL;
  int n;
  int tcl_err;

  /* If g_interp is NULL we'll cache the command. */
  if (NULL == g_interp)
  {

    /* If we need more space... */
    if (num_cached_tcl_commands >= max_num_cached_tcl_commands)
    {
      /* Allocate a new batch of cached command space. */
      new_cached_tcl_commands =
        (char**)calloc (max_num_cached_tcl_commands +
                        TCL_CACHE_INC,sizeof(char*));

      /* Copy in old commands if necessary. Note we just copy the
         pointers, we don't allocate new strings.*/
      if (NULL != cached_tcl_commands)
        for (n = 0; n < num_cached_tcl_commands; n++ )
        {
          new_cached_tcl_commands[n] = cached_tcl_commands[n];
        }

      /* Free the old storage, we're still pointing to the strings
         in the new array. */
      free (cached_tcl_commands);

      /* Store new array. */
      cached_tcl_commands = new_cached_tcl_commands;

      max_num_cached_tcl_commands =
        max_num_cached_tcl_commands + TCL_CACHE_INC;
    }

    /* Make a copy of the command. Don't use strdup because that's
       been replaced by tixStrDup and we might not have called
       tixInit yet. */
    new_cmd = (char*) calloc( strlen(cmd) + 1, sizeof(char) );
    if (NULL != new_cmd)
    {
      strcpy( new_cmd, cmd );
      cached_tcl_commands[num_cached_tcl_commands] = new_cmd;
      num_cached_tcl_commands++;
    }
    else
    {
      printf ("surfer: Error creating command storage of size %d\n",
              (int)strlen(cmd));
    }
  }
  else
  {
    /* Send the command. */
    tcl_err = Tcl_Eval (g_interp, cmd);

    /* Print an error message if we got one. */
    if (TCL_OK != tcl_err)
    {
      printf ("surfer: Error sending tcl command %s:\n\t%s\n",
              cmd, Tcl_GetStringResult (g_interp));
    }
  }
}

void send_cached_tcl_commands ()
{
  int n;

  if (NULL == g_interp)
  {
    printf ("surfer: Can't send cached commands, no interpreter present.\n");
  }

  /* For each command, send it. Then free the string. */
  for (n = 0; n < num_cached_tcl_commands; n++ )
  {
    send_tcl_command (cached_tcl_commands[n]);
    free (cached_tcl_commands[n]);
  }

  num_cached_tcl_commands = 0;
}

/* ---------------------------------------------------------------------- */

int
cptn_initialize ()
{
  cptn_format_string = (char*) calloc (CPTN_STRING_LEN+1, sizeof(char));
  if (NULL == cptn_format_string)
    ErrorReturn
      (ERROR_NO_MEMORY,
       (ERROR_NO_MEMORY,
        "cptn_initialize: Couldn't init format string"));

  cptn_value_string = (char*) calloc (CPTN_STRING_LEN+1, sizeof(char));
  if (NULL == cptn_value_string)
    ErrorReturn
      (ERROR_NO_MEMORY,
       (ERROR_NO_MEMORY,
        "cptn_initialize: Couldn't init value string"));

  /* Initial value. */
  strncpy (cptn_format_string, "Vertex: !V", CPTN_STRING_LEN);

  return ERROR_NONE;
}

int
cptn_set_format_string ( char* in )
{
  if (NULL == in)
    ErrorReturn
      (ERROR_BADPARM,
       (ERROR_BADPARM,
        "cptn_set_format_string: in was NULL"));

  if (NULL == cptn_format_string)
    ErrorReturn
      (ERROR_BADPARM,
       (ERROR_BADPARM,
        "cptn_set_format_string: format string was NULL"));

  strncpy (cptn_format_string, in, CPTN_STRING_LEN);

  return ERROR_NONE;
}

int cptn_clear_value_string ()
{
  if (NULL == cptn_value_string)
    ErrorReturn
      (ERROR_BADPARM,
       (ERROR_BADPARM,
        "cptn_clear_value_string: value string was NULL"));

  strncpy (cptn_value_string, cptn_format_string, CPTN_STRING_LEN);

  return ERROR_NONE;
}


int
cptn_sprintf_for_code (const char* code, const char* format, ... )
{
  va_list args;
  char value[1024];

  if (NULL == cptn_value_string)
    ErrorReturn
      (ERROR_BADPARM,
       (ERROR_BADPARM,
        "cptn_sprintf_for_code: value string was NULL"));

  /* sprintf the value into a string */
  va_start (args, format);
  vsnprintf (value, sizeof(value), format, args);
  va_end (args);

  /* Make the substitution in the caption. */
  cptn_substitute_code_with_value (cptn_value_string,
                                   CPTN_STRING_LEN, code, value);

  return ERROR_NONE;
}

int
cptn_draw ()
{
  int cur_char;

  if (NULL == cptn_value_string)
    ErrorReturn (ERROR_BADPARM,
                 (ERROR_BADPARM, "cptn_draw: value string was NULL"));

  glPushMatrix();
  glLoadIdentity();

  /* Draw a rectangle under our caption, both to clear the old caption
     and to provide a contrasting background. */
  glNormal3f (0.0, 0.0, 1.0);
  glColor3f (0.0, 0.0, 0.0);
  glBegin (GL_QUADS);
  glVertex3f (fov*sf * CPTN_LOC_RECTANGLE_LEFT,
              fov*sf * CPTN_LOC_RECTANGLE_BOTTOM, fov*sf * 10 );
  glVertex3f (fov*sf * CPTN_LOC_RECTANGLE_RIGHT,
              fov*sf * CPTN_LOC_RECTANGLE_BOTTOM, fov*sf * 10 );
  glVertex3f (fov*sf * CPTN_LOC_RECTANGLE_RIGHT,
              fov*sf * CPTN_LOC_RECTANGLE_TOP, fov*sf * 10 );
  glVertex3f (fov*sf * CPTN_LOC_RECTANGLE_LEFT,
              fov*sf * CPTN_LOC_RECTANGLE_TOP, fov*sf * 10 );
  glEnd ();

  /* Note that the color of bitmaps is taken from the raster color,
     which is only set when glRasterPos is called. So set the color
     first, then call glRasterPos. */
  glDisable(GL_LIGHTING); // necessary to get color setting to work
  glNormal3f (0.0, 0.0, 1.0);
  glColor3f (1.0, 1.0, 1.0);
  glRasterPos3f (fov*sf * CPTN_LOC_CAPTION_X,
                 fov*sf * CPTN_LOC_CAPTION_Y, fov*sf * 10 );

  /* Draw the caption. */
  for (cur_char = 0; cur_char < strlen(cptn_value_string); cur_char++)
  {
    glutBitmapCharacter (CPTN_FONT, cptn_value_string[cur_char]);
  }
  glEnable(GL_LIGHTING); // restore

  glFinish();

  glPopMatrix();

  return ERROR_NONE;
}

int
cptn_substitute_code_with_value (char* caption, int caption_size,
                                 const char* code, char* value)
{

  char tmp[1024];
  char* found;

  if (strlen(code) < 2)
    ErrorReturn(ERROR_BADPARM,
                (ERROR_BADPARM,
                 "cptn_subtitute_code_with_value: code invalid.\n"));

  /* Copy the caption into a tmp string, find the code and replace it
     with %s, then use sprintf to copy the value in there. Then copy
     it back into the caption and back into the tmp string, and try to
     find it again..*/
  strncpy (tmp, caption, sizeof(tmp));

  found = strstr (tmp, code);
  while (NULL != found)
  {
    found[0] = '%'; /* So "Vertex !V" turns into "Vertex %s" */
    found[1] = 's'; /* which we can use as a sprintf argument. */

    /* If the code is longer than 2 chars, we need to copy the
       string back for each extra char to take up the space. */
    if (strlen(code) > 2)
      strcpy (&found[2], &found[strlen(code)]);

    /* Sprintf the value in. */
    snprintf (caption, caption_size, tmp, value);

    /* Copy caption into tmp and find again. */
    strncpy (tmp, caption, sizeof(tmp));
    found = strstr (tmp, code );
  }

  return ERROR_NONE;
}

/* ---------------------------------------------------------------------- */

/* This code only exists if DEBUG_DRAWING_TOOLS is set. */
#if DEBUG_DRAWING_TOOLS

int
ddt_initialize ()
{
  if (NULL != mris)
  {
    if (NULL != ddt_hilited_vnos)
      free (ddt_hilited_vnos);
    if (NULL != ddt_hilited_fnos)
      free (ddt_hilited_fnos);

    ddt_hilited_vnos = (int*)calloc (mris->nvertices, sizeof(int));
    ddt_hilited_fnos = (int*)calloc (mris->nvertices, sizeof(int));

    if (NULL == ddt_hilited_vnos || NULL == ddt_hilited_fnos)
      ErrorExit(ERROR_NOMEMORY,
                "Couldn't allocate ddt storage for verts or fnos");
  }
  else
  {
    ErrorReturn(ERROR_NOT_INITED,
                (ERROR_NOT_INITED,
                 "ddt_initialize called without mris.\n"));
  }

  return (ERROR_NONE);
}

int
ddt_clear ()
{
  if (NULL != ddt_hilited_vnos)
    bzero (ddt_hilited_vnos, sizeof(int) * mris->nvertices);
  if (NULL != ddt_hilited_fnos)
    bzero (ddt_hilited_fnos, sizeof(int) * mris->nvertices);

  return (ERROR_NONE);
}

int ddt_hilite_vertex (int vno, int type)
{
  if (vno < 0 || vno >= mris->nvertices)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"ddt_hilite_vertex: bad vno.\n"));

  /* Mark it, but only with a higher type value. */
  if (NULL != ddt_hilited_vnos &&
      ddt_hilited_vnos[vno] < type)
    ddt_hilited_vnos[vno] = type;

  return (ERROR_NONE);
}

int ddt_get_hilite_vertex_color (int vno, GLubyte* r, GLubyte* g, GLubyte* b)
{
  if (ddt_hilited_vnos[vno])
  {
    switch (ddt_hilited_vnos[vno])
    {
    case 1:
      *r = 200;
      *g = 0;
      *b = 0;
      break;
    case 2:
      *r = 0;
      *g = 255;
      *b = 0;
      break;
    }
  }
  return (ERROR_NONE);
}

int ddt_hilite_face (int fno, int type)
{
  FACE* f;
  int fvno, vno;

  if (fno < 0 || fno >= mris->nfaces)
    ErrorReturn(ERROR_BADPARM,(ERROR_BADPARM,"ddt_hilite_face: bad fno.\n"));

  if (NULL != ddt_hilited_fnos)
  {
    f = &mris->faces[fno];
    for (fvno = 0; fvno < VERTICES_PER_FACE; fvno++)
    {
      vno = f->v[fvno];
      /* Mark it, but only with a higher type value. */
      if (ddt_hilited_fnos[vno] < type)
        ddt_hilited_fnos[vno] = type;
    }
  }

  return (ERROR_NONE);
}


int ddt_get_hilite_face_color (int vno, GLubyte* r, GLubyte* g, GLubyte* b)
{

  if (ddt_hilited_fnos[vno])
  {
    switch (ddt_hilited_fnos[vno])
    {
    case 1:
      *r = 100;
      *g = 0;
      *b = 0;
      break;
    case 2:
      *r = 0;
      *g = 100;
      *b = 0;
      break;
    }
  }
  return (ERROR_NONE);
}

#endif

/* ---------------------------------------------------------------------- */

int save_tiff (char* fname)
{
  GLint rowlength, skiprows, skippixels, alignment;
  GLboolean swapbytes, lsbfirst;
  GLubyte* pixel_data = NULL;
  GLenum gl_error;
  TIFF *tiff = NULL;
  tsize_t line_bytes;
  unsigned char* line_buffer = NULL;
  int scan_line_size;
  int strip_size;
  int row;
  int height, width;

  width = frame_xdim;
  height = frame_ydim;

  /* Allocate a buffer for pixels. */
  pixel_data = (GLubyte*) malloc (width * height * 3);
  if (NULL == pixel_data)
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,"Error allocating pixel storage."));

  /* Read from the front buffer. */
  glReadBuffer (GL_FRONT);

  /* Save our unpack attributes. */
  glGetBooleanv (GL_PACK_SWAP_BYTES, &swapbytes);
  glGetBooleanv (GL_PACK_LSB_FIRST, &lsbfirst);
  glGetIntegerv (GL_PACK_ROW_LENGTH, &rowlength);
  glGetIntegerv (GL_PACK_SKIP_ROWS, &skiprows);
  glGetIntegerv (GL_PACK_SKIP_PIXELS, &skippixels);
  glGetIntegerv (GL_PACK_ALIGNMENT, &alignment);

  /* Set them. */
  glPixelStorei (GL_PACK_SWAP_BYTES, GL_FALSE);
  glPixelStorei (GL_PACK_ROW_LENGTH, 0);
  glPixelStorei (GL_PACK_SKIP_ROWS, 0);
  glPixelStorei (GL_PACK_SKIP_PIXELS, 0);
  glPixelStorei (GL_PACK_ALIGNMENT, 1);

  /* Read RGB pixel data. */
  glReadPixels (0, 0, width, height, GL_RGB,
                GL_UNSIGNED_BYTE, (GLvoid*)pixel_data);

  /* Check error at this point. */
  gl_error = glGetError ();

  /* Restore the attributes. */
  glPixelStorei (GL_PACK_SWAP_BYTES, swapbytes);
  glPixelStorei (GL_PACK_LSB_FIRST, lsbfirst);
  glPixelStorei (GL_PACK_ROW_LENGTH, rowlength);
  glPixelStorei (GL_PACK_SKIP_ROWS, skiprows);
  glPixelStorei (GL_PACK_SKIP_PIXELS, skippixels);
  glPixelStorei (GL_PACK_ALIGNMENT, alignment);

  /* Handle error now. We can bail safely as we've restored the GL
     context. */
  if (GL_NO_ERROR != gl_error)
  {
    free (pixel_data);
    ErrorReturn(ERROR_NOMEMORY,(ERROR_NOMEMORY,"Error reading pixels."));
  }

  /* Open a TIFF. */
  tiff = TIFFOpen( fname, "w" );
  if (NULL == tiff)
  {
    free (pixel_data);
    ErrorReturn(ERROR_NOFILE,(ERROR_NOFILE,"Couldn't create file."));
  }

  /* Set the TIFF info. */
  TIFFSetField (tiff, TIFFTAG_IMAGEWIDTH, width);
  TIFFSetField (tiff, TIFFTAG_IMAGELENGTH, height);
  TIFFSetField (tiff, TIFFTAG_SAMPLESPERPIXEL, 3);
  TIFFSetField (tiff, TIFFTAG_BITSPERSAMPLE, 8);
  TIFFSetField (tiff, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
  TIFFSetField (tiff, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
  TIFFSetField (tiff, TIFFTAG_PHOTOMETRIC, PHOTOMETRIC_RGB);

  /* Calculate some sizes and allocate a line buffer. */
  line_bytes = 3 * width;
  line_buffer = NULL;
  scan_line_size = TIFFScanlineSize (tiff);
  if (scan_line_size != line_bytes)
  {
    fprintf (stderr,"surfer: scan_line_size %d, line_bytes %d\n",
             scan_line_size, (int)line_bytes);
  }

  line_buffer = (unsigned char*) _TIFFmalloc( scan_line_size  );
  if (NULL == line_buffer)
  {
    free (pixel_data);
    TIFFClose (tiff);
    ErrorReturn(ERROR_NOMEMORY,
                (ERROR_NOMEMORY,"Couldn't create tiff storage."));
  }

  /* Set the strip size to default. */
  strip_size = TIFFDefaultStripSize (tiff, width * 3);
  TIFFSetField (tiff, TIFFTAG_ROWSPERSTRIP, strip_size);

  /* Write line by line (bottom to top). */
  for (row = 0; row < height; row++)
  {
    memmove (line_buffer, &pixel_data[(height-row-1) * line_bytes],
            line_bytes);
    TIFFWriteScanline (tiff, line_buffer, row, 0);
  }

  /* Close the tiff file and free the line buffer. */
  TIFFClose (tiff);
  _TIFFfree (line_buffer);

  free (pixel_data);

  return (NO_ERROR);
}

/* end rkt */

static int
is_val_file(char *fname)
{
  char   *dot ;

  dot = strrchr(fname, '.') ;
  if (!dot)
    return(0) ;
  return (!stricmp(dot+1, "w"));
}

#include "cdflib.h"
#include "sig.h"

static void
f_to_p(int numer_dof, int denom_dof)
{
  int    vno ;
  VERTEX *v ;
  float sig ;

  surface_compiled = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    sig = sigf(v->val, numer_dof, denom_dof) ;
    if (v->val > 0)
      v->stat = -log10(sig) ;
    else
      v->stat = log10(sig) ;
  }
}

static void
t_to_p(int dof)
{
  int    vno ;
  VERTEX *v ;
  float  sig ;

  surface_compiled = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    sig = sigt(v->val, dof) ;
    if (v->val > 0)
      v->stat = -log10(sig) ;
    else
      v->stat = log10(sig) ;
  }
}

static void
f_to_t(void)
{
  int    vno ;
  VERTEX *v ;
  float  f ;

  surface_compiled = 0 ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->ripflag)
      continue ;
    f = v->val ;
    v->val = sqrt(abs(f)) ;
    if (f < 0)
      v->val = -v->val ;
  }
}

static void
label_to_stat(int field)
{
  sclv_new_from_label (field, labl_selected_label);
}

#include "stc.h"

static double get_dipole_value(int dipno, int ilat, int plot_type, STC *stc);
static STC  *stc ;


static void
read_soltimecourse(char *fname)
{

  surface_compiled = 0 ;
  stc = StcRead(fname) ;
  if (!stc)
    return ;
  if (stc->nvertices < stc->m_vals->rows)
    sol_nperdip = 3 ;
  else
    sol_nperdip = 1 ;
  printf("setting nperdip to %d\n", sol_nperdip) ;
  complexvalflag = 0 ;
  colscale = HEAT_SCALE ;
  overlayflag = TRUE ;
}

static void
sol_plot(int timept, int plot_type)
{
  int    ilat0,ilat1, nsum, ilat ;
  double avgval;
  int    vno, dip ;
  VERTEX *v ;
  double val ;

  if (!stc)
  {
    fprintf(stderr, "stc file not loaded.\n") ;
    return ;
  }

  surface_compiled = 0 ;
  if (sol_nperdip < 1)
  {
    fprintf(stderr,
            "set # of dipoles per location (nperdip) before plotting\n") ;
    return ;
  }

  ilat0 = floor((timept-dlat/2-stc->epoch_begin_lat)/stc->sample_period+0.5);
  ilat1 = floor((timept+dlat/2-stc->epoch_begin_lat)/stc->sample_period+0.5);
  printf("ilat0=%d, ilat1=%d\n",ilat0,ilat1);
  nsum = ilat1-ilat0+1;

  MRISsetVals(mris, 0.0) ;
  MRISclearMarks(mris) ;
  for (dip = 0 ; dip < stc->nvertices ; dip++)
  {
    vno = stc->vertices[dip] ;
    if (vno < 0 || vno >= mris->nvertices)
    {
      printf("vertex # %d out of range\n", vno);
      return ;
    }

    v = &mris->vertices[vno] ;
    avgval = 0.0 ;
    for (ilat = ilat0 ; ilat <= ilat1 ; ilat++)
    {
      val = get_dipole_value(dip, ilat, plot_type, stc) ;
      avgval += val ;
    }
    avgval /= nsum ;

    /* put fMRI weighting in here potentially */

    v->val = avgval ;
    v->fixedval = 1 ;
  }
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (v->fixedval || v->ripflag)
      continue ;
    v->undefval = TRUE ;
  }
}

static double
get_dipole_value(int dipno, int ilat, int plot_type, STC *stc)
{
  double   sum, val ;
  int      i ;

  val = *MATRIX_RELT(stc->m_vals, dipno+1, ilat+1) ;
  switch (plot_type)
  {
  default:
  case 0:
    if (sol_nperdip == 1)
      val = fabs(val) ;
    else
    {
      /* take sqrt of sum of squares */
      sum = 0.0 ;
      for (i = 0 ; i < sol_nperdip ; i++)
      {
        val = *MATRIX_RELT(stc->m_vals, dipno+1+i, ilat+1) ;
        sum += (val*val) ;
      }
      val = sqrt(sum) ;
    }
    break ;
  case 1:
    if (sol_nperdip == 3)
    {
#if 0
      val =
        *MATRIX_RELT(stc->m_vals, sol_nperdip*dip+0+1,ilat+1)*v->dipnx+
        *MATRIX_RELT(stc->m_vals, sol_nperdip*dip+1+1,ilat+1)*v->dipny+
        *MATRIX_RELT(stc->m_vals, sol_nperdip*dip+2+1,ilat+1)*v->dipnz;
#endif
    }
    break ;
  }
  return(val) ;
}
#if 1
double
dipval(int cond_num, int nzdip_num, int ilat)
{
#if 0
  int ic;
  double sum,val;

  val = 0 ;
  if (sol_plot_type==0)
  {
    sum = 0;
    for (ic=0;ic<sol_nperdip;ic++)
      sum += SQR(sol_dipcmp_val[cond_num][sol_nperdip*nzdip_num+ic][ilat]);
    val = sqrt(sum);
  }
  else
    if (sol_plot_type==1 || sol_plot_type==2)
    {
      if (sol_nperdip==1)
        val = sol_dipcmp_val[cond_num][sol_nperdip*nzdip_num][ilat];
      else
      {
        val = sol_dipcmp_val[cond_num][sol_nperdip*nzdip_num+0][ilat]*
              mris->vertices[sol_dipindex[nzdip_num]].dipnx +
              sol_dipcmp_val[cond_num][sol_nperdip*nzdip_num+1][ilat]*
              mris->vertices[sol_dipindex[nzdip_num]].dipny +
              sol_dipcmp_val[cond_num][sol_nperdip*nzdip_num+2][ilat]*
              mris->vertices[sol_dipindex[nzdip_num]].dipnz;
      }
      if (sol_plot_type==2) val = fabs(val);
    }
  return val;
#else
  return 0.0 ;
#endif
}
#endif
static void
remove_triangle_links(void)
{
  MRISremoveTriangleLinks(mris) ;
  surface_compiled = 0 ;
}

static void
drawcb(void)
{
  draw_colscalebar();
  glFlush() ;
}


static void
set_value_label_name(char *label_name, int field)
{
  char cmd[STRLEN] ;

  sprintf (cmd, "UpdateValueLabelName %d \"%s\"", field, label_name);
  send_tcl_command (cmd);
  sprintf (cmd, "ShowValueLabel %d 1", field);
  send_tcl_command (cmd);
  sprintf (cmd, "UpdateLinkedVarGroup view");
  send_tcl_command (cmd);
}
static int
draw_curvature_line(void)
{
  static int firsttime = 1 ;
  int    start_vno, current_vno, end_vno, n, best_n ;
  double dot, dx, dy, dz, odx, ody, odz, best_dot, len ;
  VERTEX *vn, *vend, *vstart, *vcurrent ;

  if (nmarked < 2)
  {
    fprintf(stderr, "must  origin and end points (previous two)\n");
    return(NO_ERROR) ;
  }

  if (firsttime)
  {
    MRISsaveVertexPositions(mris, TMP_VERTICES) ;
    read_inflated_vertex_coordinates() ;
    read_white_vertex_coordinates() ;
    MRISrestoreVertexPositions(mris, ORIGINAL_VERTICES) ;
    MRISaverageVertexPositions(mris, 25) ;
    MRIScomputeMetricProperties(mris) ;
    MRIScomputeSecondFundamentalForm(mris) ;
    /* compute local basis in tangent bundle */
    MRISrestoreVertexPositions(mris, TMP_VERTICES) ;
    firsttime = 0 ;
  }

  start_vno = current_vno = marked[0] ;
  vstart = &mris->vertices[start_vno] ;
  end_vno = marked[1] ;
  vend = &mris->vertices[end_vno] ;
  vend->marked = 0 ;

  odx = vend->x - vstart->x ;
  ody = vend->y - vstart->y ;
  odz = vend->z - vstart->z ;
  do
  {
    vcurrent = &mris->vertices[current_vno] ;
    vcurrent->marked = 1 ;
    odx = vend->infx - vcurrent->infx ;
    ody = vend->infy - vcurrent->infy ;
    odz = vend->infz - vcurrent->infz ;
    best_n = -1 ;
    best_dot = 0.0 ;
    for (n = 0 ; n < vcurrent->vnum ; n++)
    {
      vn = &mris->vertices[vcurrent->v[n]] ;
      if (vn->marked)
        continue ;   /* already in line */
      dx = vn->infx - vcurrent->infx ;
      dy = vn->infy - vcurrent->infy ;
      dz = vn->infz - vcurrent->infz ;
      dot = dx*odx + dy*ody + dz*odz ;
      if (dot < 0)
        continue ;
      dx = vn->x - vcurrent->x ;
      dy = vn->y - vcurrent->y ;
      dz = vn->z - vcurrent->z ;
      len = sqrt(dx*dx + dy*dy + dz*dz) ;
      if (FZERO(len))
        continue ;
      dx /= len ;
      dy /= len ;
      dz /= len ;
      dot = dx*vcurrent->e2x + dy*vcurrent->e2y + dz*vcurrent->e2z ;
      if (fabs(dot) > best_dot)
      {
        best_dot = fabs(dot) ;
        best_n = n ;
      }
    }
    if (best_n < 0)
      break ;
    draw_cursor(current_vno, FALSE) ;
    current_vno = vcurrent->v[best_n] ;
    if (current_vno != end_vno)
    {
#if 1
      marked[nmarked++] = current_vno ;
#else
      mark_vertex(current_vno, TRUE) ;
      draw_cursor(current_vno, TRUE) ;
      /*    redraw() ;*/
#endif
    }
  }
  while (current_vno != end_vno) ;

  redraw() ;
  return(NO_ERROR) ;
}

void
transform_brain(void)
{
  int         vno ;
  VERTEX      *v ;
  float        x, y, z, xt, yt, zt ;
  float       xlo, ylo, zlo, xhi, yhi, zhi ;

  if (!lta)
    return ;

  xhi=yhi=zhi= -10000;
  xlo=ylo=zlo= 10000;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    x = v->x ;
    y = v->y ;
    z = v->z ;
    LTAworldToWorld(lta, x, y, z, &xt, &yt, &zt) ;
    v->x = xt ;
    v->y = yt ;
    v->z = zt ;
    if (v->x > xhi) xhi = v->x;
    if (v->x < xlo) xlo = v->x;
    if (v->y > yhi) yhi = v->y;
    if (v->y < ylo) ylo = v->y;
    if (v->z > zhi) zhi = v->z;
    if (v->z < zlo) zlo = v->z;
  }

  mris->xlo = xlo ;
  mris->ylo = ylo ;
  mris->zlo = zlo ;
  mris->xctr = (xhi + xlo)/2 ;
  mris->yctr = (yhi + ylo)/2 ;
  mris->zctr = (zhi + zlo)/2 ;
  MRIScomputeMetricProperties(mris) ;
  vset_save_surface_vertices(VSET_MAIN) ;
  vset_set_current_set(vset_current_set) ;
  redraw() ;
}

static void
label_set_stats(float val)
{
  int   i ;
  LABEL *area ;

  if (labl_selected_label< 0)
    return;
  area = labl_labels[labl_selected_label].label ;
  for (i = 0 ; i < area->n_points ; i++)
    area->lv[i].stat = val ;
}

static void
label_from_stats(int field)
{
  int    n, vno, new_index ;
  VERTEX *v ;
  float  val =0.0;
  LABEL  *l ;
  char   tcl_command[STRLEN] ;

  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    sclv_get_value(v, field, &val) ;
    if  (!FZERO(val))
      n++ ;
  }

  l = LabelAlloc(n, NULL, NULL) ;
  strncpy( l->subject_name, pname, 100 );
  for (n = vno = 0 ; vno < mris->nvertices ; vno++)
  {
    v = &mris->vertices[vno] ;
    if (vno == Gdiag_no)
      DiagBreak() ;
    if (v->ripflag)
      continue ;

    sclv_get_value(v, field, &val) ;
    if  (!FZERO(val))
    {
      l->lv[n].vno = vno ;
      l->lv[n].x = v->x ;
      l->lv[n].y = v->y ;
      l->lv[n].z = v->z ;
      l->lv[n].stat = val ;
      n++ ;
    }
  }
  l->n_points = n ;
  LabelToWhite (l, mris);
  /* add this label to our list. */
  labl_add (l, &new_index);
  printf("new label created at index %d with %d points\n", new_index, n) ;

  /* select this label */
  labl_select (new_index);

  surface_compiled = 0 ;

  /* if the fill dlog is open, this will update it. */
  sprintf (tcl_command, "LabelsChanged");
  send_tcl_command (tcl_command);
}

/*!
  \fn void LoadMRISMask(void)
  \brief Part of a crude hack to allow masking
*/
void LoadMRISMask(void)
{
  extern MRI *mrismask;
  extern double mrismaskthresh;
  extern char *mrismaskfile;
  extern MRIS *mris;
  int reshapefactor;
  char *threshstring;
  MRI *mritmp;

  // check whether already loaded
  if (mrismask != NULL) return;

  if (mrismaskfile == NULL)
  {
    mrismaskfile = getenv("TKS_MRIS_MASK_FILE");
    if (mrismaskfile == NULL) return;
  }

  printf("Reading mris mask %s\n",mrismaskfile);
  mrismask = MRIread(mrismaskfile);
  if (mrismask == NULL) exit(1);

  if (mrismask->height != 1 || mrismask->depth != 1)
  {
    reshapefactor = mrismask->height * mrismask->depth;
    printf("Reshaping %d\n",reshapefactor);
    mritmp = mri_reshape(mrismask, reshapefactor*mrismask->width,
                         1, 1, mrismask->nframes);
    MRIfree(&mrismask);
    mrismask = mritmp;
    reshapefactor = 0; /* reset for output */
  }
  if (mrismask->width != mris->nvertices)
  {
    printf("ERROR: dimension mismatch between mask and surf\n");
    exit(1);
  }

  if (mrismaskthresh < 0)
  {
    threshstring = getenv("TKS_MRIS_MASK_THRESH");
    if (threshstring == NULL) return;
    sscanf(threshstring,"%lf",&mrismaskthresh);
    printf("mris mask thresh = %lf\n",mrismaskthresh);
  }
}

int dngheat(float f, float *r, float *g, float *b)
{
  static float x0 = 0.0;
  static float x1 = 1.0/3.0;
  static float x2 = 1.0;
  static float a1 = 0.5625;
  static float a2 = 0.4375; //1-a1
  float absf, c1=0, cg=0, c3=0;

  absf = fabs(f);

  if(absf >= x0 && absf <= x1){
    c1 = a1 + a2*absf/(x1-x0);
    cg = 0.0;
    c3 = 0.0;
  }
  if(absf > x1 && absf <= x2){
    c1 = 1;
    cg = (absf-x1)/(x2-x1);
    c3 = 0.0;
  }
  if(absf > x2){
    c1 = 1.0;
    cg = 1.0;
    c3 = 0.0;
  }

  *g = cg;
  if(f >= 0.0){
    *r = c1;
    *b = c3;
  } else {
    *r = c3;
    *b = c1;
  }
  return(0);
}

int dngcolorwheel(float f, float *r, float *g, float *b)
{
  static float x0 = 0.0;
  static float x1 = 1.0/6.0;
  static float x2 = 2.0/6.0;
  static float x3 = 3.0/6.0;
  static float x4 = 4.0/6.0;
  static float x5 = 5.0/6.0;
  static float x6 = 1.0;
  float f2, xA, xB;

  f2 = f/2; // expecting -1 <= f <= +1, so -0.5 <= f2 <= +0.5
  if(fabs(f2) > 0.5) f2 = remainder(f2,0.5); // force it
  f2 = f2 + 0.5; // recenter so that 0 <= f2 <= +1

  xA = x0;
  xB = x1;
  if(f2 >= xA && f2 < xB){
    *r = 1.0;
    *g = (f2-xA)/(xB-xA);
    *b = 0.0;
    return(0);
  }

  xA = x1;
  xB = x2;
  if(f2 >= xA && f2 < xB){
    *r = 1.0 - (f2-xA)/(xB-xA);
    *g = 1.0;
    *b = 0.0;
    return(0);
  }

  xA = x2;
  xB = x3;
  if(f2 >= xA && f2 < xB){
    *r = 0.0;
    *g = 1.0;
    *b = (f2-xA)/(xB-xA);
    return(0);
  }

  xA = x3;
  xB = x4;
  if(f2 >= xA && f2 < xB){
    *r = 0.0;
    *g = 1.0 - (f2-xA)/(xB-xA);
    *b = 1.0;
    return(0);
  }

  xA = x4;
  xB = x5;
  if(f2 >= xA && f2 < xB){
    *r = (f2-xA)/(xB-xA);
    *g = 0.0;
    *b = 1.0;
    return(0);
  }

  xA = x5;
  xB = x6;
  *r = 1.0;
  *g = 0.0;
  *b = 1.0 - (f2-xA)/(xB-xA);
  return(0);

}
#define MAX_STEPS 100
void 
set_area_thresh(float target_area)
{
  char                  cmd[STRLEN];
  int   vno ;
  float thresh_min, thresh, thresh_max, val, area, best_thresh, thresh_step, best_area ;

  val = 0 ;
  thresh_min = 1e10 ; thresh_max = -thresh_min ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    sclv_get_value( &mris->vertices[vno],sclv_current_field, &val);
    if (val < thresh_min)
      thresh_min = val ;
    if (val > thresh_max)
      thresh_max = val ;
  }

  thresh_step = (thresh_max - thresh_min) / MAX_STEPS ;
  best_thresh = thresh_min ;
  best_area = -1 ;
  for (thresh = thresh_min ; thresh <= thresh_max ; thresh += thresh_step)
  {
    area = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      sclv_get_value( &mris->vertices[vno],sclv_current_field, &val);
      if (val > thresh)
        area += mris->vertices[vno].area ;
    }
    if (best_area < 0 || fabs(area-target_area) < best_area)
    {
      best_area = area ;
      best_thresh = thresh ;
    }
  }
  fthresh = best_thresh ;
  fslope = 1/best_thresh ;
  fmid = 2*fthresh ;
  printf("best threshold at %2.5f (%2.5f), %2.2fmm\n", fthresh, fslope, best_area) ;
  sclv_field_info[sclv_current_field].fslope = fslope ;
  sclv_field_info[sclv_current_field].fthresh = fthresh;
  vertex_array_dirty = 1 ;
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);
}

void 
set_thresh(float target_area)
{
  char                  cmd[STRLEN];
  int   vno ;
  float thresh_min, thresh, thresh_max, val, area, best_thresh, thresh_step, best_area ;

  val = 0 ;
  thresh_min = 1e10 ; thresh_max = -thresh_min ;
  for (vno = 0 ; vno < mris->nvertices ; vno++)
  {
    sclv_get_value( &mris->vertices[vno],sclv_current_field, &val);
    if (val < thresh_min)
      thresh_min = val ;
    if (val > thresh_max)
      thresh_max = val ;
  }

  thresh_step = (thresh_max - thresh_min) / MAX_STEPS ;
  best_thresh = thresh_min ;
  best_area = -1 ;
  for (thresh = thresh_min ; thresh <= thresh_max ; thresh += thresh_step)
  {
    area = 0 ;
    for (vno = 0 ; vno < mris->nvertices ; vno++)
    {
      sclv_get_value( &mris->vertices[vno],sclv_current_field, &val);
      if (val > thresh)
        area += mris->vertices[vno].area ;
    }
    if (best_area < 0 || fabs(area-target_area) < best_area)
    {
      best_area = area ;
      best_thresh = thresh ;
    }
  }
  fthresh = best_thresh ;
  fslope = 1/best_thresh ;
  fmid = 2*fthresh ;
  printf("best threshold at %2.5f (%2.5f), %2.2fmm\n", fthresh, fslope, best_area) ;
  sclv_field_info[sclv_current_field].fslope = fslope ;
  sclv_field_info[sclv_current_field].fthresh = fthresh;
  vertex_array_dirty = 1 ;
  sprintf (cmd, "UpdateLinkedVarGroup overlay");
  send_tcl_command (cmd);
}

