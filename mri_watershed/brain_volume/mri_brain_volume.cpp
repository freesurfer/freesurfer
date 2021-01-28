/**
 * @brief modified mri_watershed.c
 *
 */
/*
 * Original Author: Florent Segonne & Bruce Fischl
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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <sys/types.h>
#include <sys/time.h>
#include <unistd.h>
#include <time.h>

#include "TVector.h"

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "mrisurf.h"
#include "matrix.h"
#include "proto.h"
#include "stats.h"
#include "timer.h"
#include "const.h"
#include "mrishash.h"
#include "icosahedron.h"
#include "tritri.h"
#include "timer.h"
#include "chklc.h"
#include "diag.h"

using namespace std;

#define WM_CONST 110 /* not used anymore */
#define MAX_INT 100 /*100% is a good value for the watershed algo */

#define ITER_SMOOTH 10
#define NBR_AVGS 20
#define NBR_NGBS 2

#define DEFAULT_MODE 0
#define DIST_MODE 1
#define CURV_MODE 2

#define DEBUG_CURVE 0
#define OUTPUT 0
#define OUTPUT_CURVES 0
#define OUTPUT_SURFACES 0

#define VERBOSE_MODE 0


#define MAX_MASK_VOLUMES  50

typedef struct Cell {
  unsigned char type;
  void * next;
}
Cell;

typedef struct Basincell {
  unsigned char depth;
  unsigned long size;
  unsigned long ambiguous;
}
BasinCell;

typedef struct Bound {
  unsigned char x,y,z,val;
  struct Bound *next;
}
Bound;

typedef unsigned char Coord[3];


typedef struct STRIP_PARMS {
  /*  float  fill_level;*/
  int template_deformation;
  /*to save the surfaces into the output volume*/
  int surf_dbg;
  /*to write out the brain surface into the file surfname*/
  /*the brain surface is shrank inward of h_shrk mm*/
  int brainsurf;
  /*to write out all the BEM surfaces : brain, outer and inner skull, scalp*/
  int surf;
  /*to labelize the volume into scalp, skull, csf, white and gray*/
  int label;
  /*to use the atlas validation and correction*/
  int atlas;

  /*specify T1 volume*/
  int T1;
  /*specify the center fo the brain and its radius*/
  int cx,cy,cz,rb;

  char *surfname;
  int h_shk;
  int skull_type;
  int watershed_analyze;
  int threshold_analyze;

  int seed_coord[30][4];
  int nb_seed_points/*=0*/;
  unsigned char hpf;

  int manual_params;
  int manual_CSF_MAX,manual_TRANSITION_INTENSITY,manual_GM_INTENSITY;

  double forceParam;
}
STRIP_PARMS ;


typedef struct {
  float direction[26][3];
  MRIS *mrisphere,*mris,*mris_curv,*mris_var_curv,*mris_dCOG
  ,*mris_var_dCOG;

  double xCOG,yCOG,zCOG,rad_Brain;
  double xsCOG,ysCOG,zsCOG;
  int i_global_min/*=0*/,j_global_min,k_global_min,int_global_min;
  unsigned long estimated_size/*=0*/,main_basin_size/*=0*/;
  unsigned long brain_size /*=0*/;
  unsigned long basinnumber,basinsize;

  MRI *mri_src,*mri_dst,*mri_orig;
  int width,height,depth;

  unsigned char Imax;
  int WM_INTENSITY,WM_VARIANCE,WM_HALF_MAX,WM_HALF_MIN,WM_MAX,WM_MIN;
  int CSF_intensity,CSF_HALF_MAX,CSF_MAX,CSF_MIN;
  int GM_MIN, GM_INTENSITY,TRANSITION_INTENSITY;


  unsigned long gmnumber[256];

  /*file saving - not used */
#if OUTPUT_CURVES
  FILE *fout;
#endif
#if OUTPUT_SURFACES
  FILE *fsvout,*fsfout;
#endif

  Bound *Bound1,*Bound2;

  Cell *** Basin;

  Coord** Table[256];

  unsigned char intbasin[256];
  unsigned long tabdim[256];
  unsigned long sqrdim[256];
  unsigned long count[256];

  Coord* T1Table;
  long T1nbr;

  int decision;
  float scale;

  int atlas;
  int validation;
  int verbose_mode;

}
MRI_variables;

const char *Progname;

static int type_changed = 0 ;
static int old_type ;


static void Error(const char *string);
static int get_option(int argc, char *argv[],STRIP_PARMS *parms) ;
static STRIP_PARMS* init_parms(void);
static MRI_variables* init_variables(MRI *mri_with_skull);
MRI *MRIstripSkull(MRI *mri_with_skull, MRI *mri_without_skull,
                   STRIP_PARMS *parms);
static void MRIVfree(MRI_variables *MRI_var);
/*WATERSHED FUNCTIONS*/
static int Watershed(STRIP_PARMS *parms,MRI_variables *MRI_var);
static void Allocation(MRI_variables *MRI_var);
static int calCSFIntensity(MRI_variables *MRI_var);
static int calCOGMAX(MRI_variables *MRI_var, int *x, int *y, int *z);
static int Pre_CharSorting(STRIP_PARMS *parms, MRI_variables *MRI_var);
static void analyseWM(double *tab,MRI_variables *MRI_var);
static BasinCell* CreateBasinCell(int val, unsigned long size, unsigned long ambiguous);
static int Decision(STRIP_PARMS *parms,  MRI_variables *MRI_var);
static int CharSorting(MRI_variables *MRI_var);
static int sqroot(int r);
static int Analyze(STRIP_PARMS *parms,MRI_variables *MRI_var);
static Cell* FindBasin(Cell *cell);
static int Lookat(int,int,int,unsigned char,int*,Cell**,int*,Cell* adtab[27],
                  STRIP_PARMS *parms,MRI_variables *MRI_var);
static int Test(Coord crd,STRIP_PARMS *parms,MRI_variables *MRI_var);
static Cell* TypeVoxel(Cell *cell);
static int PostAnalyze(STRIP_PARMS *parms,MRI_variables *MRI_var);
static int Merge(unsigned char i,unsigned char j,unsigned char k,int val,int *n,MRI_variables *MRI_var);
static int AddVoxel(MRI_variables *MRI_var);
static int AroundCell(unsigned char i,unsigned char j,unsigned char k,MRI_variables *MRI_var);
static int MergeRoutine(unsigned char,unsigned char,unsigned char,int,int*,MRI_variables *MRI_var);
static int FreeMem(MRI_variables *MRI_var);

/*TEMPLATE DEFORMATION FUNCTIONS*/
static void MRISfit(MRI_variables *MRI_var,
                    void (*calcforce)
                    (double &force0, double &force1, double &force,
                     double sd,
                     const TVector &Pos, const TVector &S, const TVector &N,
                     MRI_variables *MRI_var, double forceParam),
                    double forceParam);
static void calcForce1(double &force0, double &force1, double &force,
                double sd,
                const TVector &Pos, const TVector &S, const TVector &N,
                MRI_variables *MRI_var,
                double forceParam);
static void calcForce2(double &force0, double &force1, double &force,
                double sd,
                const TVector &Pos, const TVector &S, const TVector &N,
                MRI_variables *MRI_var,
                double forceParam);
static void calcForceGM(double &force0, double &force1, double &force,
                 double sd,
                 const TVector &Pos, const TVector &S, const TVector &N,
                 MRI_variables *MRI_var,
                 double forceParam);
static void calcForceMine(
  double &force0, double &force1, double &force,
  double sd,
  const TVector &Pos, const TVector &S, const TVector &N,
  MRI_variables *MRI_var,
  double forceParam);

static void  read_geometry(int type,MRI_variables *MRI_var,char *surf_fname);
static void Template_Deformation(STRIP_PARMS *parms,MRI_variables *MRI_var);
static void brain_params(MRI_variables *MRI_var);
static void init_surf_to_image(float rx, float ry, float rz,MRI_variables *MRI_var);
static void write_image(MRI_variables *MRI_var, int val);
static void init_direction(MRI_variables *MRI_var);
static void find_normal(float nx,float ny, float nz,float* n1,float *n2,  float direction[26][3]);
static unsigned long MRISpeelBrain(float h,MRI *mri_dst,MRIS *mris,unsigned char val);
static void mean(float tab[4][9],float *moy);
static void MRISsmooth_surface(MRI_SURFACE *mris,int niter);
static void MRISshrink_surface(MRIS *mris,int h);
static void MRISshrink_Outer_Skin(MRI_variables *MRI_var,MRI* mri_src);
static void label_voxels(STRIP_PARMS *parms, MRI_variables *MRI_var,MRI* mri_with_skull);
/*VALIDATION - SURFACE CORRECTION*/
static void MRISchangeCoordinates(MRI_SURFACE *mris,MRI_SURFACE *mris_orig);
/*mri->type correction*/

// declare function pointer
int  (*MyvoxelToWorld)(MRI *mri, double xv, double yv, double zv,
                       double *xw, double *yw, double *zw) ;
int  (*MyworldToVoxel)(MRI *mri, double xw, double yw, double zw,
                       double *pxv, double *pyv, double *pzv) ;

void getTimeString(char *buf) {
  time_t t;
  time(&t);
  struct tm *lTime= localtime(&t);
  sprintf(buf, "%d-%d-%d", lTime->tm_mon+1, lTime->tm_mday, lTime->tm_year+1900);
}

int main(int argc, char *argv[]) {
  char  *in_fname, *out_fname;
  int nargs;
  MRI *mri_with_skull, *mri_without_skull=NULL, *mri_mask;

  STRIP_PARMS *parms;

  Progname=argv[0];

  parms=init_parms();

  /************* Command line****************/

  fprintf(stderr,"\n");

  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
    nargs = get_option(argc, argv,parms) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc<2) {
    fprintf(stderr, "\nUsage: %s [options] input_file [output_file]", Progname);
    fprintf(stderr, "\noptional command -forceParam val: change pushout force (default 1.0)");
    fprintf(stderr,"\noptional command --version : to show the current version\n\n");

    exit(1);
  };

  if (argc ==3) {
    in_fname = argv[argc-2];
    out_fname = argv[argc-1];

    fprintf(stderr, "\n************************************************************"
            "\nInput:\t%s"
            "\nOutput:\t%s\n", in_fname,out_fname);
  } else {
    in_fname = argv[argc-1];
    out_fname = 0;

    fprintf(stderr, "\n************************************************************"
            "\nInput:\t%s", in_fname);
  }
  /*************** PROG *********************/
  char *debug = getenv("DEBUG_BRAIN");
  if (debug)
    parms->surf_dbg =1;


  /* initialisation */
  mri_with_skull = MRIread(in_fname) ;
  if (!mri_with_skull)
    Error("read failed\n");


  if (mri_with_skull->type!=MRI_UCHAR) {
    MRI *mri_tmp ;

    type_changed = 1 ;
    old_type = mri_with_skull->type ;
    printf("changing type of input volume to 8 bits/voxel...\n") ;
    mri_tmp = MRIchangeType(mri_with_skull, MRI_UCHAR, 0.0, 0.999, FALSE) ;
    MRIfree(&mri_with_skull) ;
    mri_with_skull = mri_tmp ;
  } else
    type_changed = 0 ;

  /////////////////////////////////////////////////////
  // binarize the image
  //                                      thresh low  high
  MRIbinarize(mri_with_skull, mri_with_skull, 1, 0, 128);
  if (parms->surf_dbg) {
    char timeString[256];
    getTimeString(timeString);
    char filename[256];
    int req = snprintf(filename, 256, "%s-%s.mgh", "binarized", timeString);
    if (req >= 256) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRIwrite(mri_with_skull, filename);
  }
  /* Main routine *********************/
  mri_without_skull=MRIstripSkull(mri_with_skull, mri_without_skull,parms);
  if (mri_without_skull == NULL) {
    printf("mri_brain_volume failed.\n");
    free(parms);
    return -1;
  }
  if (type_changed)  /* make output volume the same type as input volume */
  {
    mri_with_skull = MRIread(in_fname) ;
    if (!mri_with_skull)
      Error("read failed\n");
    mri_mask = mri_without_skull ;
    mri_without_skull = MRImask(mri_with_skull, mri_mask, NULL, 0, 0) ;
  } else
    mri_mask = mri_without_skull ;


  fprintf(stderr,"\n\n******************************\nSave...");

  if (out_fname)
    MRIwrite(mri_without_skull,out_fname);

  MRIfree(&mri_with_skull) ;

  fprintf(stderr,"done\n");

  MRIfree(&mri_mask) ;

  free(parms);

  return 0;
}

/*-----------------------------------------------------
        Parameters:message error

        Returns value:void

        Description: Error routine - stop the prog
------------------------------------------------------*/
static void Error(const char *string) {
  fprintf(stderr, "\nError %s\n",string) ;
  exit(1) ;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:number of options read

        Description: read the different options of the command line
------------------------------------------------------*/

static int
get_option(int argc, char *argv[],STRIP_PARMS *parms) {
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */

  option = argv[1] + 1 ;            /* past '-' */
  if (!strcmp(option, "-version")) {
    fprintf(stderr,"##version##\n");
    exit(0);
  } else if (!strcmp(option, "forceParam")) {
    parms->forceParam = atof(argv[2]);
    fprintf(stderr, "use forceParam = %.2f\n", parms->forceParam);
    nargs = 1;
  } else if (!strcmp(option, "surf_debug")) {
    parms->surf_dbg=1;
    fprintf(stderr,"Mode:          Writing out surfaces into output volume\n") ;
    nargs = 0 ;
  } else {
    printf("Mode:          unknown option %s\n", argv[1]) ;
    exit(1) ;
  }
  // always don't do watershed
  parms->watershed_analyze = 0;
  return(nargs) ;
}

/*-----------------------------------------------------
        Parameters:void

        Returns value:STRIP_PARMS*

        Description: allocate the structure containing some parameters
  necessary for the program
------------------------------------------------------*/

static STRIP_PARMS* init_parms(void) {
  STRIP_PARMS* sp=(STRIP_PARMS*)calloc(1,sizeof(STRIP_PARMS));

  /*preflooding height used in the watershed segmentation*/
  sp->hpf=25;

  /* no writing out of the different surfaces into the ouput volume
     usefull for understanding where the errors come from*/
  sp->surf_dbg=0;
  /* no brain bem surface saving*/
  sp->brainsurf=0;
  /*no bem surfaces saving*/
  sp->surf=0;
  /*no shrank bem brain surface*/
  sp->h_shk=0;
  /*level of surface: shrunk:-1, normal:0 or expanded:1 */
  sp->skull_type=0;
  /*no labelization of the tissue into the output volume*/
  sp->label=0;
  /*  mode post-watershed analyze turned on*/
  sp->watershed_analyze=1;
  /*post-watershed analyze threshold set to 100%*/
  sp->threshold_analyze=100;
  /*no seed points*/
  sp->nb_seed_points=0;
  /*mode template deformation on*/
  sp->template_deformation=1;
  /*surface outfile name NULL*/
  sp->surfname=NULL;

  /*no manual parameters entered*/
  sp->manual_params=0;
  /*no atlas analysis*/
  sp->atlas=0;

  /*no T1 volume normalization*/
  sp->T1=0;

  /*no input brain parms*/
  sp->cx=-1;
  sp->rb=-1;

  sp->forceParam = 1.0;

  return(sp);
}

/*-----------------------------------------------------
        Parameters: MRI *mri_with_skull (input image)

        Returns value:MRI_variables *MRI_Var

        Description: allocate the structure containing some important
  variables necessary for the program
------------------------------------------------------*/

static MRI_variables* init_variables(MRI *mri_with_skull) {
  MRI_variables* v=(MRI_variables*)calloc(1,sizeof(MRI_variables));

  v->mris=NULL;
  v->mris_curv=NULL;
  v->mris_var_curv=NULL;
  v->mris_dCOG=NULL;
  v->mris_var_dCOG=NULL;
  v->mrisphere=NULL;

  v->i_global_min=0;
  v->estimated_size=0;
  v->main_basin_size=0;
  v->brain_size=0;

  v->width=mri_with_skull->width;
  v->height=mri_with_skull->height;
  v->depth=mri_with_skull->depth;

  v->mri_orig=mri_with_skull;

  v->T1Table=NULL;
  v->T1nbr=0;

  return v;

}

/*-----------------------------------------------------
  FUNCTION MRIstripSkull

  Parameters:
    MRI *mri_with_skull:orig input image (orig or T1 T1-weigthed volume)
    MRI *mri_without_skull: output volume (could be NULL)
    STRIP_PARMS *parms: coudl be NULL

  Returns value:
    MRI *mri_without_skull: (if input NULL, after allocation)

    Description: strip the skull from the input image
------------------------------------------------------*/

MRI *MRIstripSkull(MRI *mri_with_skull, MRI *mri_without_skull,
                   STRIP_PARMS *parms) {
  char fname[512];
  MRI_variables *MRI_var;
  MRI *mri_tp;
  double vol_elt;

#if OUTPUT_SURFACES
  char filename[512];
#endif

  if (mri_with_skull==NULL)
    Error("\nNULL input volume !\n");
  if (mri_with_skull->type!=0)
    Error("\nThe type of the input file is not 0 : UCHAR\n");

  mri_tp=MRIclone(mri_with_skull,NULL);
  mri_tp=MRIcopy(mri_with_skull,NULL);

  if (!mri_without_skull)
    mri_without_skull=MRIclone(mri_with_skull,NULL);
  else if (mri_without_skull->type!=mri_with_skull->type)
    Error("\ndifferent types of mri structures...\n");

  MRI_var=init_variables(mri_with_skull);

  MRI_var->verbose_mode=parms->surf_dbg;

  MRI_var->mri_src=mri_tp;
  MRI_var->mri_dst=mri_without_skull;

  /*watershed process*/
  if (Watershed(parms,MRI_var)==-1) {
    free(mri_tp);
    return NULL;
  }

  // template deformation
  {
    if (parms->surf_dbg) {
      MRI_var->mri_dst=MRIclone(mri_with_skull,NULL);
      MRI_var->mri_dst=MRIcopy(mri_with_skull,NULL);
    }

    /*template process*/
    ////////////////////////////////////////////////////////////////////////////
    Template_Deformation(parms,MRI_var);

    /*in case the src volume was modified (scaling of the intensity)*/
    free(mri_tp);
    mri_tp=MRIclone(mri_with_skull,NULL);
    mri_tp=MRIcopy(mri_with_skull,NULL);
    MRI_var->mri_src = mri_tp ;

    // mri_src is modified
    // non-brain part becomes      baloon    dst         src surface  val
    // return brain_size
    MRI_var->brain_size=MRISpeelBrain(0,MRI_var->mri_src,MRI_var->mris,0);

    vol_elt=MRI_var->mri_src->xsize*MRI_var->mri_src->ysize*MRI_var->mri_src->zsize;

    fprintf(stderr, "\n\nforceParam:%s%.2f\n", "\t", parms->forceParam);
    fprintf(stderr, "Brain Size:%s%ld\tvoxels\n","\t", MRI_var->brain_size);
    fprintf(stderr, "voxel size:%s%.2f\tmm3\n", "\t", (double) vol_elt);
    printf("%.2f\n", ((double) MRI_var->brain_size)*vol_elt);

    /*save the surface of the brain*/
    if (parms->brainsurf) {

      if (parms->surf || parms->h_shk)
        MRISsmooth_surface(MRI_var->mris,5);

      if (parms->h_shk != 0)
        MRISshrink_surface(MRI_var->mris,parms->h_shk);

      /*writing out the surface*/
      sprintf(fname,"%s",parms->surfname);
      strcat(fname,"_brain_surface");
      MRI_var->mris->useRealRAS = 1;
      MRISwrite(MRI_var->mris,fname);
    }

  }

  /*find and write out the surfaces of the inner skull, scalp and outer skull*/
  if (parms->template_deformation && parms->surf) {
    //inner skull
    MRISshrink_surface(MRI_var->mris,-3);
    MRISsmooth_surface(MRI_var->mris,5);

    if (parms->surf_dbg)
      write_image(MRI_var, 128);

    /*writing out the inner skull surface*/
    sprintf(fname,"%s",parms->surfname);
    strcat(fname,"_inner_skull_surface");
    MRISwrite(MRI_var->mris,fname);


    //scalp
    MRISfree(&MRI_var->mris);
    read_geometry(1,MRI_var,NULL);

    init_surf_to_image(1.8*MRI_var->rad_Brain,1.8*MRI_var->rad_Brain,1.8*MRI_var->rad_Brain,MRI_var);

    fprintf(stderr,"\n      outer skin surface matching...");
    MRISshrink_Outer_Skin(MRI_var,mri_with_skull);

    MRISsmooth_surface(MRI_var->mris,3);

    if (parms->surf_dbg)
      write_image(MRI_var, 180);

    /*writing out the surface*/
    sprintf(fname,"%s",parms->surfname);
    strcat(fname,"_outer_skin_surface");
    MRISwrite(MRI_var->mris,fname);

    //outer skull
    MRISsmooth_surface(MRI_var->mris,3);
    MRISshrink_surface(MRI_var->mris,3);
    MRISsmooth_surface(MRI_var->mris,5);

    if (parms->surf_dbg)
      write_image(MRI_var,200);

    /*writing out the outer skull surface*/
    sprintf(fname,"%s",parms->surfname);
    strcat(fname,"_outer_skull_surface");
    MRISwrite(MRI_var->mris,fname);

    if (parms->label)
      label_voxels(parms,MRI_var,mri_with_skull);

  }


  /*save the volume with the surfaces written in it*/
  /*used to visualize the surfaces -> debuging */
  if (parms->template_deformation && parms->surf_dbg)
    mri_without_skull=MRIcopy(MRI_var->mri_dst,NULL);
  /*normal mode saving*/
  else
    mri_without_skull=MRIcopy(MRI_var->mri_src,NULL);


  MRIVfree(MRI_var);
  return mri_without_skull;
}
/******************************************************
 ******************************************************
 *****************************************************/
/*-----------------------------------------------------
        FUNCTION Watershed

        Parameters:
    STRIP_PARMS *:contains the parameters for the prog
    MRI_variables *: contains the variables

        Returns value:void

        Description: watershed segment the input volume *mri_src
------------------------------------------------------*/

static int Watershed(STRIP_PARMS *parms,MRI_variables *MRI_var) {

  if (parms->nb_seed_points)
    fprintf(stderr, "\n%d more seed points created",parms->nb_seed_points);
  fprintf(stderr,"\npreflooding height equal to %d percent", parms->hpf) ;

  fprintf(stderr,"\nSorting...");
  Allocation(MRI_var);
  if (Pre_CharSorting(parms,MRI_var)==-1)
    return -1;

  CharSorting(MRI_var);
  fprintf(stderr,"\ndone");

  fprintf(stderr,"\nAnalyze\n");
  Analyze(parms,MRI_var);
  fprintf(stderr,"\ndone");

  fprintf(stderr,"\nPostAnalyze...");
  PostAnalyze(parms,MRI_var);
  fprintf(stderr,"done\n");

  return 0;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:void

        Description: Allocation of a table of basins
------------------------------------------------------*/
static void Allocation(MRI_variables *MRI_var) {
  int k,j;

  MRI_var->Basin=(Cell ***)malloc(MRI_var->depth*sizeof(Cell **));
  if (!MRI_var->Basin)
    Error("first allocation  failed\n");

  for (k=0;k<MRI_var->depth;k++) {
    MRI_var->Basin[k]=(Cell **)malloc(MRI_var->height*sizeof(Cell*));
    if (!MRI_var->Basin[k])
      Error("second allocation  failed\n");
    for (j=0;j<MRI_var->height;j++) {
      MRI_var->Basin[k][j]=(Cell*)calloc(MRI_var->width,sizeof(Cell));
      if (!MRI_var->Basin[k][j])
        Error("third alloc failed\n");
    }
  }

  for (k=0;k<256;k++) {
    MRI_var->tabdim[k]=0;
    MRI_var->sqrdim[k]=0;
    MRI_var->count[k]=0;
    MRI_var->intbasin[k]=k;
    MRI_var->gmnumber[k]=0;
  }
}

//
// calculate CSF_intensity
//
static int calCSFIntensity(MRI_variables *MRI_var) {
  int i, j, k;
  int n;
  double intensity_percent[256];
  BUFTYPE *pb;

  /*First estimation of the CSF */
  for (k=0;k<256;k++)
    intensity_percent[k]=0;

  n=0; // counts non-zero grey voxels
  // create a histogram
  for (k=2;k<MRI_var->depth-2;k++)
    for (j=2;j<MRI_var->height-2;j++) {
      pb=&MRIvox(MRI_var->mri_src,0,j,k);
      pb+=2;
      for (i=2;i<MRI_var->width-2;i++) {
        if (*pb) // non-zeo
        {
          n++;
          intensity_percent[*pb]++;
        }
        pb++;
      }
    }

  // accumulate histogram
  for (k=1;k<256;k++) {
    intensity_percent[k]+=intensity_percent[k-1];
    /*  fprintf(stderr," %d %f ",k,intensity_percent[k]*100/n);*/
    // the max grey value to be 99.8% of the population
    // CSF_intensity to be 10% of that value.
    if (intensity_percent[k]*100<=99.8*n)
      MRI_var->CSF_intensity=k/10+1;
  }
  return 0;
}

static int calCOGMAX(MRI_variables *MRI_var, int *x, int *y, int *z) {
  int i, j, k;
  int n;
  double intensity_percent[256];
  BUFTYPE *pb;

  /*Ignore everything which is bellow CSF_intensity
    Then first estimation of the COG coord
    Find a cube in which it will search for WM */
  for (k=0;k<256;k++) {
    intensity_percent[k]=0;
  }

  n=0;
  MRI_var->xCOG = MRI_var->yCOG = MRI_var->zCOG = 0;
  for (k=2;k<MRI_var->depth-2;k++)
    for (j=2;j<MRI_var->height-2;j++) {
      pb=&MRIvox(MRI_var->mri_src,0,j,k);
      pb+=2;
      for (i=2;i<MRI_var->width-2;i++) {
        if (*pb> 0) {
          n++;
          intensity_percent[*pb]++;
          MRI_var->xCOG+=i;
          MRI_var->yCOG+=j;
          MRI_var->zCOG+=k;
        }
        pb++;
      }
    }
  if (n==0)
    Error("\n pbm in the COG calculation ");

  MRI_var->xCOG/=n;
  MRI_var->yCOG/=n;
  MRI_var->zCOG/=n;

  *x=(int)(MRI_var->xCOG+0.5);
  *y=(int)(MRI_var->yCOG+0.5);
  *z=(int)(MRI_var->zCOG+0.5);

  // calculate Imax
  MRI_var->Imax=0;
  for (k=1;k<256;k++) {
    intensity_percent[k]+=intensity_percent[k-1];
    if (intensity_percent[k]*100<=n*MAX_INT)
      MRI_var->Imax=k;
  }
  return 0;
}

static int calBrainRadius(MRI_variables *MRI_var) {
  int m=0;
  int i,j,k;
  BUFTYPE *pb;

  m=0;
  MRI_var->rad_Brain=0;
  for (k=2;k<MRI_var->depth-2;k++)
    for (j=2;j<MRI_var->height-2;j++) {
      pb=&MRIvox(MRI_var->mri_src,0,j,k);
      pb+=2;
      for (i=2;i<MRI_var->width-2;i++) {
        if ((*pb)>=MRI_var->Imax)
          *pb=MRI_var->Imax-1;
        if (*pb)      /*don't care about 0 intensity voxel*/
          MRI_var->tabdim[*pb]++;
        if (*pb>MRI_var->CSF_intensity) {
          m++;
          MRI_var->rad_Brain+=SQR(i-MRI_var->xCOG)+SQR(j-MRI_var->yCOG)+SQR(k-MRI_var->zCOG);
        }
        pb++;
      }
    }

  if (m==0)
    Error("\n pbm with the radius calculation ");

  MRI_var->rad_Brain=sqrt(MRI_var->rad_Brain/m);

  return 0;
}
/*-----------------------------------------------------
  int Pre_CharSorting

        Parameters:

        Returns value:int

        Description: Calculate some rough statistics
               Estimate the white matter parameters
               Create a global minimum in the WM
------------------------------------------------------*/
static int Pre_CharSorting(STRIP_PARMS *parms,MRI_variables *MRI_var) {
  int retVal;
  int i,j,k,n,m,u,v;
  int ig,jg,kg;
  BUFTYPE *pb,*pbc[3][3];
  unsigned long wmint=0,wmnb=0;
  unsigned long number[256];
  double intensity_percent[256];
  float ***mean_val,mean,min,max;
  float ***var_val,var;
  float ***mean_var;
  int x,y,z,r;
  int xmin,xmax,ymin,ymax,zmin,zmax,mint;
  float tmp;
  MRI *mri_tmp;

  // calculate CSF_intensity from voxel values (CSF_intensity)
  calCSFIntensity(MRI_var);

  // calculate initial estimate of COG coords and Imax (xCOG, yCOG, zCOG, and Imax)
  calCOGMAX(MRI_var, &x, &y, &z);

  // calculate intitial estimate of brain radius (rad_Brain)
  calBrainRadius(MRI_var);

  r=(int)(MRI_var->rad_Brain/2);
  if (r == 0)
    r = (int) (MRI_var->rad_Brain+1);

  fprintf(stderr,"\n      first estimation of the COG coord: x=%.1f y=%.1f z=%.1f r=%.2f",
          MRI_var->xCOG,MRI_var->yCOG,MRI_var->zCOG,MRI_var->rad_Brain);


  if (parms->cx!=-1) {
    MRI_var->xCOG=parms->cx;
    MRI_var->yCOG=parms->cy;
    MRI_var->zCOG=parms->cz;
    fprintf(stderr,"\n      modification of the brain COG: x=%d y=%d z=%d",
            (int)MRI_var->xCOG,(int)MRI_var->yCOG,(int)MRI_var->zCOG);
    x=(int)(MRI_var->xCOG+0.5);
    y=(int)(MRI_var->yCOG+0.5);
    z=(int)(MRI_var->zCOG+0.5);
  }
  if (parms->rb!=-1) {
    MRI_var->rad_Brain=parms->rb;
    r=(int)(MRI_var->rad_Brain/2);
    fprintf(stderr,"\n      modification of the brain radius to %.2f", MRI_var->rad_Brain);
  }

  MRI_var->estimated_size=(unsigned long)(4.19*pow(MRI_var->rad_Brain,3));// 4.19 = 4PI/3

  fprintf(stderr,"\n      first estimation of the main basin volume: %ld voxels",MRI_var->estimated_size);


  mint=MIN(MRI_var->width,MIN(MRI_var->height,MRI_var->width));
  if (20*r>=mint*9)
    Error("\n main radius too high");

  /*allocate the Cube memory: mean intensity, variance, mean variance */

  r*=2;
  mean_val=(float***)malloc(r*sizeof(float**));
  var_val=(float***)malloc(r*sizeof(float**));
  mean_var=(float***)malloc(r*sizeof(float**));
  for (k=0;k<r;k++) {
    mean_val[k]=(float**)malloc(r*sizeof(float*));
    var_val[k]=(float**)malloc(r*sizeof(float*));
    mean_var[k]=(float**)malloc(r*sizeof(float*));
    for (j=0;j<r;j++) {
      mean_val[k][j]=(float*)malloc(r*sizeof(float));
      var_val[k][j]=(float*)malloc(r*sizeof(float));
      mean_var[k][j]=(float*)malloc(r*sizeof(float));
    }
  }
  r/=2;


  xmin=MAX(2,x-r-1);
  xmax=MIN(MRI_var->width-2,xmin+2*r);
  xmin=xmax-2*r;
  ymin=MAX(2,y-2*r-1);
  ymax=MIN(MRI_var->height-2,ymin+2*r);
  ymin=ymax-2*r;
  zmin=MAX(2,z-r-1);
  zmax=MIN(MRI_var->depth-2,zmin+2*r);
  zmin=zmax-2*r;

  /*Calculate the mean intensity and the variance (mean = 27 voxels) */

  for (k=0;k<256;k++) {
    number[k]=0;
    intensity_percent[k]=0;
  }

  for (k=zmin;k<zmax;k++)
    for (j=ymin;j<ymax;j++) {
      for (u=0;u<3;u++)
        for (v=0;v<3;v++) {
          pbc[u][v]=&MRIvox(MRI_var->mri_src,0,j+u,k+v);
          pbc[u][v]+=xmin;
        }
      for (i=xmin;i<xmax;i++) {
        mean=0;
        for (u=0;u<3;u++)
          for (v=0;v<3;v++)
            for (n=0;n<3;n++)
              mean+=(*(pbc[u][v]+n));
        mean/=27;
        mean_val[k-zmin][j-ymin][i-xmin]=mean;

        if (mean>2*MRI_var->CSF_intensity && mean<MRI_var->Imax-1) {
          var=0;
          for (u=0;u<3;u++)
            for (v=0;v<3;v++)
              for (n=0;n<3;n++)
                var+=SQR((*(pbc[u][v]+n))-mean);

          var/=27;
          var_val[k-zmin][j-ymin][i-xmin]=var;

        } else
          var_val[k-zmin][j-ymin][i-xmin]=1000;


        for (u=0;u<3;u++)
          for (v=0;v<3;v++)
            pbc[u][v]++;
      }
    }

  /*- Find the mean variance (27 voxels)
    - And find the mean variance for each intensity
      divided by the number of voxels of the same intensity
                          -> estimation of the MRI_var->WM_INTENSITY */

  r*=2;
  min=1000;
  max=0;
  for (k=0;k<r;k++)
    for (j=0;j<r;j++)
      for (i=0;i<r;i++) {
#if GCC_VERSION > 80000
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wint-in-bool-context"
#endif
        if (!(i*j*k*(i-r+1)*(j-r+1)*(k-r+1))) {
#if GCC_VERSION > 80000
#pragma GCC diagnostic pop
#endif
          mean=1000;
        } else {
          mean=0;
          for (u=-1;u<2;u++)
            for (v=-1;v<2;v++)
              for (n=-1;n<2;n++)
                mean+=var_val[k+u][j+v][i+n];

          mean/=27;
          if (min>=mean) {
            min=mean;
            if (mean==0) {
              wmint+=(int)(mean_val[k][j][i]+0.5);
              wmnb++;
            }
          }
          if (max<mean)
            max=mean;
          number[(int)(mean_val[k][j][i]+0.5)]++;
          intensity_percent[(int)(mean_val[k][j][i]+0.5)]+=var_val[k][j][i];
        }
        mean_var[k][j][i]=mean;
      }
  if (wmnb)
    wmint/=wmnb;

  /*  fprintf(stderr," min=%f max=%f wmint=%ld",min,max,wmint);*/

  tmp=0;
  for (k=0;k<256;k++)
    if (number[k]<100)
      intensity_percent[k]=0;
    else {
      intensity_percent[k]=SQR(number[k])/intensity_percent[k];
      if (intensity_percent[k]>tmp)
        tmp=intensity_percent[k];
    }

  analyseWM(intensity_percent,MRI_var);

  if (MRI_var->WM_INTENSITY < MRI_var->CSF_intensity) {
    fprintf(stderr, "\n\n\n\n********************************************************\n");
    fprintf(stderr, "********************************************************\n");
    fprintf(stderr, "********************************************************\n");
    fprintf(stderr, "White matter intensity %d is lower than CSF intensity %d.\n",
            MRI_var->WM_INTENSITY, MRI_var->CSF_intensity);
    fprintf(stderr, "Please examine input images.  Will terminate ...\n");
    fprintf(stderr, "********************************************************\n");
    fprintf(stderr, "********************************************************\n");
    fprintf(stderr, "********************************************************\n");
    return -1;
  }

#if DEBUG_CURVE
  for (k=0;k<256;k++) {
    fprintf(stderr,"\n%d %ld %f",k,number[k],(float)(number[k])/
            MAX(1,intensity_percent[k]));
    for (n=0;n<intensity_percent[k]/tmp*50;n++)
      fprintf(stderr,".");
  }
#endif

  m=0;
  MRI_var->WM_VARIANCE=0;
  for (n=MRI_var->WM_HALF_MIN;n<=MRI_var->WM_HALF_MAX;n++) {
    m++;
    MRI_var->WM_VARIANCE += int( number[n]/intensity_percent[n]);
  }

  if (m==0)
    Error("\n Pbm with the variance calculation ");

  MRI_var->WM_VARIANCE=  int(sqrt(double(MRI_var->WM_VARIANCE/m)));

  MRI_var->WM_MIN=(MRI_var->WM_MIN+MRI_var->WM_VARIANCE)/2;

  if (MRI_var->WM_MIN<=MRI_var->WM_HALF_MIN-3*MRI_var->WM_VARIANCE/2)
    MRI_var->WM_MIN=MRI_var->WM_HALF_MIN-3*MRI_var->WM_VARIANCE/2;



  if (parms->T1) {
    MRI_var->WM_INTENSITY=WM_CONST;
    MRI_var->WM_HALF_MAX=MRI_var->WM_HALF_MIN=WM_CONST;
    MRI_var->WM_VARIANCE=5;
    MRI_var->WM_MAX=WM_CONST;
    MRI_var->WM_MIN=WM_CONST;
  }
  if ((fabs(double(MRI_var->WM_INTENSITY-WM_CONST))<=2)
      && (MRI_var->WM_VARIANCE<3)) {
    if (fabs(double(MRI_var->WM_HALF_MAX-WM_CONST))<=2) {
      MRI_var->WM_MIN=MIN(MRI_var->WM_MIN,WM_CONST);
      MRI_var->WM_HALF_MIN=MIN(MRI_var->WM_HALF_MIN,WM_CONST);
      MRI_var->WM_INTENSITY=WM_CONST;
      MRI_var->WM_HALF_MAX=WM_CONST;
      MRI_var->WM_MAX=MAX(MRI_var->WM_MIN,WM_CONST);
    }
  }

  ///////////////////////////////////////////////////////////////////
  retVal = Decision(parms,MRI_var);
  if (retVal > 0) {
    /*find the WM coord */
    tmp=max;
    for (k=1;k<r-1;k++)
      for (j=1;j<r-1;j++)
        for (i=1;i<r-1;i++)
          if (mean_val[k][j][i]==MRI_var->WM_HALF_MAX)
            if (mean_var[k][j][i]<tmp) {
              tmp=mean_var[k][j][i];
              x=i;
              y=j;
              z=k;
            }

    i=xmin+1+x;
    j=ymin+1+y;
    k=zmin+1+z;

    /*Create the global minimum*/
    if (MRI_var->i_global_min) {
      i=MRI_var->i_global_min;
      j=MRI_var->j_global_min;
      k=MRI_var->k_global_min;
    }
    MRI_var->tabdim[MRIvox(MRI_var->mri_src,i,j,k)]--;
    MRI_var->int_global_min=MRIvox(MRI_var->mri_src,i,j,k);
    MRIvox(MRI_var->mri_src,i,j,k)=MRI_var->Imax;
    MRI_var->i_global_min=i;
    MRI_var->j_global_min=j;
    MRI_var->k_global_min=k;
    MRI_var->tabdim[MRI_var->Imax]++;
    MRI_var->Basin[k][j][i].type=3;
    MRI_var->Basin[k][j][i].next=(BasinCell*)CreateBasinCell(MRI_var->Imax,1,0);


    ig=MRI_var->i_global_min;
    jg=MRI_var->j_global_min;
    kg=MRI_var->k_global_min;

    n=parms->nb_seed_points;
    while (n) {
      i=parms->seed_coord[n-1][0];
      j=parms->seed_coord[n-1][1];
      k=parms->seed_coord[n-1][2];
      parms->seed_coord[n-1][3]=MRIvox(MRI_var->mri_src,i,j,k);
      MRI_var->tabdim[parms->seed_coord[n-1][3]]--;
      MRIvox(MRI_var->mri_src,i,j,k)=MRI_var->Imax;
      MRI_var->tabdim[MRI_var->Imax]++;
      if (MRI_var->Basin[k][j][i].type!=3) {
        MRI_var->Basin[k][j][i].type=1;
        MRI_var->Basin[k][j][i].next=(Cell*)(&MRI_var->Basin[kg][jg][ig]);
        ((BasinCell*)(MRI_var->Basin[kg][jg][ig].next))->size++;
      }
      n--;
    }
    if (parms->T1) {
      /*allocate temp MRI struct and init to zero*/
      mri_tmp=MRIclone(MRI_var->mri_src,NULL);
      for (k=0;k<MRI_var->depth;k++)
        for (j=0;j<MRI_var->height;j++) {
          pb=&MRIvox(mri_tmp,0,j,k);
          for (i=0;i<MRI_var->width;i++) {
            (*pb)=0;
            pb++;
          }
        }
      /*set Voxel=WM to 1*/
      for (n=0;n<MRI_var->T1nbr;n++) {
        i=MRI_var->T1Table[n][0];
        j=MRI_var->T1Table[n][1];
        k=MRI_var->T1Table[n][2];
        MRIvox(mri_tmp,i,j,k)=1;
      }
      free(MRI_var->T1Table);
      MRI_var->T1Table=NULL;
      MRI_var->T1nbr=0;

      /*go through the whole temp struct and keep the WM inside voxels*/
      for (k=3;k<MRI_var->depth-3;k++)
        for (j=3;j<MRI_var->height-3;j++)
          for (i=3;i<MRI_var->width-3;i++) {
            r=1;
            for (u=-1;u<2;u++)
              for (v=-1;v<2;v++)
                for (m=-1;m<2;m++) {
                  r*=MRIvox(mri_tmp,i+u,j+v,k+m);
                  if (r==0)
                    break;
                }
            if (r) {
              MRI_var->tabdim[MRIvox(MRI_var->mri_src,i,j,k)]--;
              MRIvox(MRI_var->mri_src,i,j,k)=MRI_var->Imax;
              MRI_var->tabdim[MRI_var->Imax]++;
              if (MRI_var->Basin[k][j][i].type!=3) {
                MRI_var->Basin[k][j][i].type=1;
                MRI_var->Basin[k][j][i].next=(Cell*)(&MRI_var->Basin[kg][jg][ig]);
                ((BasinCell*)(MRI_var->Basin[kg][jg][ig].next))->size++;
              }
            }
          }
      MRIfree(&mri_tmp);
    }


    if (!min && wmint==110 && MRI_var->WM_MAX==110) {
      MRI_var->WM_INTENSITY=110;
      MRI_var->WM_VARIANCE=5;
      MRI_var->WM_MAX=110;
    }


    fprintf(stderr,"\n      global maximum in x=%d, y=%d, z=%d, Imax=%d",
            MRI_var->i_global_min,MRI_var->j_global_min,MRI_var->k_global_min,MRI_var->Imax);


    fprintf(stderr,"\n      CSF=%d, WM_INTENSITY=%d, WM_VARIANCE=%d",
            MRI_var->CSF_intensity,MRI_var->WM_INTENSITY,MRI_var->WM_VARIANCE);
    fprintf(stderr,"\n      WM_MIN=%d, WM_HALF_MIN=%d,WM_HALF_MAX=%d ,WM_MAX=%d ",
            MRI_var->WM_MIN,MRI_var->WM_HALF_MIN,MRI_var->WM_HALF_MAX,MRI_var->WM_MAX);

    if (MRI_var->WM_VARIANCE>20) {
      fprintf(stderr,"\n probable pbm with the variance");
      MRI_var->WM_VARIANCE=15;
    }

    if (MRI_var->WM_MIN<=2*MRI_var->CSF_intensity) {
      fprintf(stderr,"\n probable pbm with WM_MIN");
      MRI_var->WM_MIN= int(MAX(MRI_var->WM_INTENSITY/2,MRI_var->WM_INTENSITY-1.5*MRI_var->WM_VARIANCE));
      MRI_var->WM_HALF_MIN=MAX(MRI_var->WM_HALF_MIN,MRI_var->WM_INTENSITY-MRI_var->WM_VARIANCE);
    }
  }
  // end of if(Decision())

  /*free memory*/
  for (k=0;k<r;k++) {
    for (j=0;j<r;j++) {
      free(mean_val[k][j]);
      free(var_val[k][j]);
      free(mean_var[k][j]);
    }
    free(mean_val[k]);
    free(var_val[k]);
    free(mean_var[k]);
  }
  free(mean_val);
  free(var_val);
  free(mean_var);

  if (retVal == -1)
    return -1;
  else
    return 0;
}

/*-----------------------------------------------------
        Parameters:

        Returns value:void

        Description: Analyze the white matter curve
------------------------------------------------------*/
static void analyseWM(double *tab,MRI_variables *MRI_var) {
  int k,n;
  double tmp,buff;
  double mean;
  double a,b,Sxy,Sxx,Sx,Sy;

  /*for(k=0;k<256;k++)
    fprintf(MRI_var->fout, " %f",tab[k]);
    fprintf(MRI_var->fout, "\n");*/

  if (MRI_var->decision)
    for (k=(int)MIN(MRI_var->WM_MAX*MRI_var->scale,256);k<256;k++)
      tab[k]=0;

  /*tmp=0;
   for(k=0;k<256;k++)
     if (tab[k]>tmp)
       tmp=tab[k];

   for(k=0;k<256;k+=4)
   {
     fprintf(stderr,"\n%d %f",k,(float)tab[k]);
     for(n=0;n<tab[k]*50/tmp;n++)
       fprintf(stderr,".");
   }*/


  MRI_var->WM_INTENSITY=0;
  tmp=0;
  for (k=0;k<256;k++) {

    buff=tab[k];
    if (buff>tmp) {
      tmp=buff;
      MRI_var->WM_INTENSITY=k;
    }
  }

  if (MRI_var->WM_INTENSITY>240)
    Error("\nw=White Matter =Intensity too high (>240)...valid input ?");

  k=MRI_var->WM_INTENSITY;
  for (n=k;n>0;n--)
    if (tab[n]>=tab[k]/2)
      MRI_var->WM_HALF_MIN=n;
    else
      break;

  if (MRI_var->WM_HALF_MIN==1)
    Error("\n pbm with WM_HALF_MIN");

  mean=0;
  for (k=0;k<MRI_var->WM_HALF_MIN;k++)
    mean+=tab[k];
  mean/=MRI_var->WM_HALF_MIN;


  for (k=0;k<256;k++)
    tab[k]=MAX(0,tab[k]-mean);

  k=MRI_var->WM_INTENSITY;
  for (n=k;n<255;n++)
    if (tab[n]>=tab[k]/2)
      MRI_var->WM_HALF_MAX=n;
    else
      break;

  if (MRI_var->WM_HALF_MAX==254)
    Error("\n pbm with WM_HALF_MAX");

  k=MRI_var->WM_HALF_MAX;
  for (n=k;n<255;n++)
    if (tab[n]>=tab[k]/3)
      MRI_var->WM_MAX=n;
    else
      break;

  if (MRI_var->WM_MAX==254)
    Error("\n pbm with WM_MAX");


  /*least-square distance interpolation*/
  if (MRI_var->WM_INTENSITY<MRI_var->WM_MAX) {
    MRI_var->WM_MAX+=1;

    n=MRI_var->WM_MAX-MRI_var->WM_INTENSITY+1;
    Sxy = Sx = Sy = Sxx = 0;
    for (k=MRI_var->WM_INTENSITY;k<=MRI_var->WM_MAX;k++) {
      Sxy+=(float)k*tab[k];
      Sx+=k;
      Sy+=tab[k];
      Sxx+=k*k;
    }

    a=(n*Sxy-Sy*Sx)/(n*Sxx-Sx*Sx);
    b=-(a*Sx-Sy)/n;

    if (a==0)
      Error("\n Interpolation Problem in the white matter curve analysis\n");

    MRI_var->WM_MAX= -int(b/a) ;
  }


  k=MRI_var->WM_INTENSITY;
  for (n=k;n>1;n--)
    if (tab[n]>=tab[k]/2)
      MRI_var->WM_HALF_MIN=n;
    else
      break;

  if (MRI_var->WM_HALF_MIN==2)
    Error("\n pbm with WM_HALF_MIN");

  k=MRI_var->WM_HALF_MIN;
  for (n=k;n>=1;n--)
    if (tab[n]>=tab[k]/2)
      MRI_var->WM_MIN=n-1;
    else
      break;

  if (MRI_var->WM_MIN==2)
    Error("\n pbm with WM_MIN");


  /*least-square distance interpolation*/
  if (MRI_var->WM_INTENSITY>MRI_var->WM_MIN) {
    MRI_var->WM_MIN=3*MRI_var->WM_MIN/2-MRI_var->WM_INTENSITY/2;

    n=MRI_var->WM_INTENSITY-MRI_var->WM_MIN+1;
    Sxy = Sx = Sy = Sxx = 0;
    for (k=MRI_var->WM_MIN;k<=MRI_var->WM_INTENSITY;k++) {
      Sxy+=(float)k*tab[k];
      Sx+=k;
      Sy+=tab[k];
      Sxx+=k*k;
    }

    a=(n*Sxy-Sy*Sx)/(n*Sxx-Sx*Sx);
    b=-(a*Sx-Sy)/n;

    if (a==0)
      Error("\n interpolation pbm in the white matter analysis");

    MRI_var->WM_MIN=int(MAX(0,-b/a));
  }
}


/*-----------------------------------------------------
        Parameters:

        Returns value:BasinCell*

        Description: Allocate and Init a BasinCell
------------------------------------------------------*/

static BasinCell* CreateBasinCell(int val, unsigned long size, unsigned long ambiguous) {
  BasinCell *bcell;

  bcell=(BasinCell*)malloc(sizeof(BasinCell));
  bcell->depth=(unsigned char)val;
  bcell->size=size;
  bcell->ambiguous=ambiguous;

  return bcell;
}
/*-----------------------------------------------------
        Parameters:

        Returns value:int value (success or error)

        Description: Decide if the white matter parameters are valid
               If too low, scale the image intensity
         such that wm=160
------------------------------------------------------*/
static int Decision(STRIP_PARMS *parms,  MRI_variables *MRI_var) {
  int i,j,k;
  double scale;
  if (MRI_var->WM_INTENSITY+MRI_var->WM_VARIANCE>=80)
    return 1;
  else {
    scale= 160.0/((double) MRI_var->WM_MAX);
    // if scaled once, the second time gives the same scale and thus bail out
    if (fabs(scale -1.) < .001) // allow epsilon error of .1%
      return -1;

    MRI_var->decision=1;
    MRI_var->scale= (float) scale;

    for (k=0;k<MRI_var->depth;k++)
      for (j=0;j<MRI_var->height;j++)
        for (i=0;i<MRI_var->width;i++)
          MRIvox(MRI_var->mri_src,i,j,k)=(unsigned short)MIN(255,scale*MRIvox(MRI_var->mri_src,i,j,k));
    fprintf(stderr,"\nmean intensity too low !");
    fprintf(stderr,"\nModification of the image: intensity*%2.2f",scale);
    for (k=0;k<256;k++) {
      MRI_var->tabdim[k]=0;
      MRI_var->sqrdim[k]=0;
      MRI_var->count[k]=0;
      MRI_var->intbasin[k]=k;
      MRI_var->gmnumber[k]=0;
    }
    return Pre_CharSorting(parms,MRI_var);
  }
}


/*-----------------------------------------------------
        Parameters:

        Returns value:

        Description: Sorting of the voxel in an ascending order
------------------------------------------------------*/
static int CharSorting(MRI_variables *MRI_var) {
  size_t i,j,k,u,v;
  size_t l;
  ldiv_t ld;
  BUFTYPE *pb;
  unsigned char val;

  /*allocating an 256 table of Coord** in order to process the Sorting*/

  for (k=1;k< (size_t) MRI_var->Imax+1;k++) {
    l=sqroot(MRI_var->tabdim[k]);
    MRI_var->sqrdim[k]=l;
    MRI_var->Table[k]=(Coord**)malloc(l*sizeof(Coord*));
    if (!MRI_var->Table[k]) Error("Allocation first Table Echec");
    for (u=0;u < l;u++) {
      MRI_var->Table[k][u]=(Coord*)calloc(l,sizeof(Coord));
      if (!MRI_var->Table[k][u]) Error("Allocation second Table Echec");
    }
  }

  /*Sorting itself*/

  for (k=2;k< (size_t) MRI_var->depth-2;k++)
    for (j=2;j< (size_t) MRI_var->height-2;j++) {
      pb=&MRIvox(MRI_var->mri_src,0,j,k);
      pb+=2;
      for (i=2;i< (size_t) MRI_var->width-2;i++) {
        val=*pb++;
        if (val) {
          l=MRI_var->count[val]++;
          ld=ldiv(l,MRI_var->sqrdim[val]);
          u=ld.quot;
          v=ld.rem;
          MRI_var->Table[val][u][v][0]=i;
          MRI_var->Table[val][u][v][1]=j;
          MRI_var->Table[val][u][v][2]=k;
        }
      }
    }
  return 0;
}

/*sqrt routine*/
static int sqroot(int r) {
  int i;
  i=(int)sqrt(double(r));
  if (i*i<r) return (i+1);
  else return i;
}

/*******************************ANALYZE****************************/

/*routine that analyzes all the voxels sorted in an descending order*/
static int Analyze(STRIP_PARMS *parms,MRI_variables *MRI_var) {
  int pos;
  int u,v,n,d;
  unsigned long l;
  ldiv_t ld;
  double vol_elt;

  MRI_var->basinnumber=0;
  MRI_var->basinsize=0;


  n=MRI_var->sqrdim[MRI_var->Imax];
  for (u=0;u<n;u++)
    free(MRI_var->Table[MRI_var->Imax][u]);
  free(MRI_var->Table[MRI_var->Imax]);

  for (pos=MRI_var->Imax-1;pos>0;pos--) {
    l=0;
    d=MRI_var->tabdim[pos];
    n=MRI_var->sqrdim[pos];
    while (l < (unsigned long) d) {
      ld=ldiv(l,n);
      u=ld.quot;
      v=ld.rem;
      Test(MRI_var->Table[pos][u][v],parms,MRI_var);
      l++;
    }
    for (u=0;u<n;u++)
      free(MRI_var->Table[pos][u]);
    free(MRI_var->Table[pos]);

    fprintf(stderr,"\r      %3d%%... %8ld basins; main size = %8ld         ",
            (MRI_var->Imax-pos)*100/(MRI_var->Imax-1),MRI_var->basinnumber,MRI_var->basinsize);
  }


  MRI_var->main_basin_size+=((BasinCell*)MRI_var->Basin
                             [MRI_var->k_global_min]
                             [MRI_var->j_global_min]
                             [MRI_var->i_global_min].next)->size;
  for (n=0;n<parms->nb_seed_points;n++)
    MRI_var->main_basin_size+=((BasinCell*)MRI_var->Basin
                               [parms->seed_coord[n][2]]
                               [parms->seed_coord[n][1]]
                               [parms->seed_coord[n][0]].next)->size;

  vol_elt=MRI_var->mri_src->xsize*MRI_var->mri_src->ysize*MRI_var->mri_src->zsize;
  fprintf(stderr,"\n      main basin size=%8ld voxels, voxel volume =%2.3f         ",
          MRI_var->main_basin_size,(float)vol_elt);
  fprintf(stderr,"\n           = %2.0f mmm3 = %2.3f cm3"
          ,MRI_var->main_basin_size*vol_elt,(float)MRI_var->main_basin_size/1000.*vol_elt);
  MRIvox(MRI_var->mri_src,MRI_var->i_global_min,MRI_var->j_global_min,MRI_var->k_global_min)
  =MRI_var->int_global_min;
  for (n=0;n<parms->nb_seed_points;n++)
    MRIvox(MRI_var->mri_src,parms->seed_coord[n][0],parms->seed_coord[n][1],parms->seed_coord[n][2])
    =parms->seed_coord[n][3];


  return 0;
}

/*looking at a voxel, finds the corresponding basin*/
static Cell* FindBasin(Cell *cell) {
  cell=(Cell *) cell->next;
  while (cell->type==1)
    cell=(Cell *) cell->next;
  return cell;
}

/*main routine for the merging*/
static int Lookat(int i,int j,int k,unsigned char val,int *dpt,Cell* *admax,int *nb,Cell * adtab[27],STRIP_PARMS *parms, MRI_variables *MRI_var) {
  int t,n;
  unsigned char d,hp=parms->hpf;
  Cell *add=&MRI_var->Basin[k][j][i];

  /*looks if the basin has already been processing*/

  if (add->type) {
    if (add->type==1) {
      add=FindBasin(add);

      d=((BasinCell*)(add->next))->depth;

      if (d>*dpt) {
        *admax=add;
        *dpt=d;
      }
      if (100*(d-val)<hp*MRI_var->Imax) {
        t=0;
        for (n=0;n<*nb;n++)
          if (add==adtab[n])
            t=1;
        if (!t) {
          adtab[*nb]=add;
          (*nb)++;
        }
      }
    }
    else {
      d=((BasinCell*)(add->next))->depth;


      if (d>*dpt) {
        *admax=add;
        *dpt=d;
      }
      if (100*(d-val)<hp*MRI_var->Imax) {
        t=0;
        for (n=0;n<*nb;n++)
          if (add==adtab[n])
            t=1;
        if (!t) {
          adtab[*nb]=add;
          (*nb)++;
        }
      }
    }
    return 0;
  }
  else if (100*(val-MRIvox(MRI_var->mri_src,i,j,k))<hp*MRI_var->Imax)
    return 1;

  return 0;
}


/*tests a voxel, merges it or creates a new basin*/
static int Test(Coord crd,STRIP_PARMS *parms,MRI_variables *MRI_var) {
  int n,nb=0,dpt=-1;
  unsigned char val;
  int mean,var,tp=0;
  int a,b,c;

  int i=crd[0],j=crd[1],k=crd[2];

  Cell  *adtab[27],*admax=&MRI_var->Basin[k][j][i];

  val=MRIvox(MRI_var->mri_src,i,j,k);

  Lookat(i,j,k-1,val,&dpt,&admax,&nb,adtab,parms,MRI_var);
  Lookat(i,j,k+1,val,&dpt,&admax,&nb,adtab,parms,MRI_var);
  Lookat(i,j-1,k,val,&dpt,&admax,&nb,adtab,parms,MRI_var);
  Lookat(i,j+1,k,val,&dpt,&admax,&nb,adtab,parms,MRI_var);
  Lookat(i-1,j,k,val,&dpt,&admax,&nb,adtab,parms,MRI_var);
  Lookat(i+1,j,k,val,&dpt,&admax,&nb,adtab,parms,MRI_var);

  /*creates a new basin*/

  if (parms->watershed_analyze) {
    mean=0;
    var=0;
    for ( a = -1 ; a<2 ; a++)
      for ( b = -1 ; b<2 ; b++)
        for ( c = -1 ; c<2 ; c++) {
          tp=MRIvox(MRI_var->mri_src,i+a,j+b,k+c);
          mean+=tp;
          var+=SQR(tp);
        }
    mean/=27;
    var=var/27-SQR(mean);
    if (mean>=MRI_var->WM_MIN && mean<=MRI_var->WM_MAX && var<=MRI_var->WM_VARIANCE)
      tp=1;
    else tp=0;
  }

  if (dpt==-1) {
    MRI_var->Basin[k][j][i].type=2;
    if (!tp)
      MRI_var->Basin[k][j][i].next=(BasinCell*)CreateBasinCell(val,1,0);
    else
      MRI_var->Basin[k][j][i].next=(BasinCell*)CreateBasinCell(val,1,1);
    MRI_var->basinnumber++;
    return 0;
  };

  /*Merging*/

  if (admax->type==3 && val<=MRI_var->WM_MIN && val>0)
    MRI_var->gmnumber[val]++;

  MRI_var->Basin[k][j][i].type=1;
  MRI_var->Basin[k][j][i].next=(Cell*)admax;
  ((BasinCell*)(admax->next))->size++;
  if (tp && admax->type!=3)
    ((BasinCell*)(admax->next))->ambiguous++;

  if (!tp || admax->type==3) {
    for (n=0;n<nb;n++)
      if (adtab[n]!=admax) {
        adtab[n]->type=1;
        ((BasinCell*)(admax->next))->size+=((BasinCell*)(adtab[n]->next))->size;
        free((BasinCell*)(adtab[n]->next));
        adtab[n]->next=(Cell*)admax;
        MRI_var->basinnumber--;
      }
  } else {
    for (n=0;n<nb;n++)
      if (adtab[n]!=admax) {
        adtab[n]->type=1;
        ((BasinCell*)(admax->next))->size+=((BasinCell*)(adtab[n]->next))->size;
        ((BasinCell*)(admax->next))->ambiguous+=((BasinCell*)(adtab[n]->next))->ambiguous;
        free((BasinCell*)(adtab[n]->next));
        adtab[n]->next=(Cell*)admax;
        MRI_var->basinnumber--;
      }
  }


  if (((BasinCell*)(admax->next))->size>MRI_var->basinsize) {
    MRI_var->basinsize=((BasinCell*)(admax->next))->size;
  }

  return 0;
}

/**************POST_ANALYZE***********************/

/*looks if the voxel basin is the main one, ie type==3*/
static Cell* TypeVoxel(Cell *cell) {
  Cell* cell1=cell;

  while (cell1->type==1)
    cell1=(Cell *) cell1->next;

  return cell1;

}


/*Looks if the voxel is a border from the segmented brain*/
static int AroundCell(unsigned char i,unsigned char j,unsigned char k,MRI_variables *MRI_var) {
  int val=0,n=0;

  if (MRI_var->Basin[k][j][i-1].type) {
    val+=MRIvox(MRI_var->mri_src,i-1,j,k);
    n++;
  }
  if (MRI_var->Basin[k][j][i+1].type) {
    val+=MRIvox(MRI_var->mri_src,i+1,j,k);
    n++;
  }
  if (MRI_var->Basin[k][j-1][i].type) {
    val+=MRIvox(MRI_var->mri_src,i,j-1,k);
    n++;
  }
  if (MRI_var->Basin[k][j+1][i].type) {
    val+=MRIvox(MRI_var->mri_src,i,j+1,k);
    n++;
  }
  if (MRI_var->Basin[k-1][j][i].type) {
    val+=MRIvox(MRI_var->mri_src,i,j,k-1);
    n++;
  }
  if (MRI_var->Basin[k+1][j][i].type) {
    val+=MRIvox(MRI_var->mri_src,i,j,k+1);
    n++;
  }
  if (6-n)
    return ((val+MRIvox(MRI_var->mri_src,i,j,k))/(n+1));
  else
    return 0;
}


/*Merge voxels which intensity is near the intensity of border voxels*/
static int MergeRoutine(unsigned char i,unsigned char j,unsigned char k,int val,int *n,MRI_variables *MRI_var) {
  int cond=15*val;
  Bound *buff;
  Cell cell=MRI_var->Basin[k][j][i];

  if (cell.type<2)
    if (abs(MRIvox(MRI_var->mri_src,i,j,k)-val)*100<cond) {
      if (cell.type==1) {
        ((Bound*)cell.next)->val+=val;
        ((Bound*)cell.next)->val/=2;
      } else {
        cell.type=1;
        buff=(Bound*)malloc(sizeof(Bound));
        buff->x=i;
        buff->y=j;
        buff->z=k;
        buff->val=val;
        cell.next=(Bound*)buff;
        (*n)++;
        buff->next=MRI_var->Bound2;
        MRI_var->Bound2=buff;
      }
    }
  return 0;
}


static int Merge(unsigned char i,unsigned char j,unsigned char k,int val,int *n,MRI_variables *MRI_var) {

  MergeRoutine(i-1,j,k,val,n,MRI_var);
  MergeRoutine(i+1,j,k,val,n,MRI_var);
  MergeRoutine(i,j-1,k,val,n,MRI_var);
  MergeRoutine(i,j+1,k,val,n,MRI_var);
  MergeRoutine(i,j,k-1,val,n,MRI_var);
  MergeRoutine(i,j,k+1,val,n,MRI_var);

  MergeRoutine(i-1,j-1,k-1,val,n,MRI_var);
  MergeRoutine(i-1,j-1,k,val,n,MRI_var);
  MergeRoutine(i-1,j-1,k+1,val,n,MRI_var);
  MergeRoutine(i-1,j,k-1,val,n,MRI_var);
  MergeRoutine(i-1,j,k+1,val,n,MRI_var);
  MergeRoutine(i-1,j+1,k-1,val,n,MRI_var);
  MergeRoutine(i-1,j+1,k,val,n,MRI_var);
  MergeRoutine(i-1,j+1,k+1,val,n,MRI_var);


  MergeRoutine(i+1,j-1,k-1,val,n,MRI_var);
  MergeRoutine(i+1,j-1,k,val,n,MRI_var);
  MergeRoutine(i+1,j-1,k+1,val,n,MRI_var);
  MergeRoutine(i+1,j,k-1,val,n,MRI_var);
  MergeRoutine(i+1,j,k+1,val,n,MRI_var);
  MergeRoutine(i+1,j+1,k-1,val,n,MRI_var);
  MergeRoutine(i+1,j+1,k,val,n,MRI_var);
  MergeRoutine(i+1,j+1,k+1,val,n,MRI_var);

  MergeRoutine(i,j-1,k-1,val,n,MRI_var);
  MergeRoutine(i,j-1,k+1,val,n,MRI_var);
  MergeRoutine(i,j+1,k-1,val,n,MRI_var);
  MergeRoutine(i,j+1,k+1,val,n,MRI_var);

  return 0;
}

static int AddVoxel(MRI_variables *MRI_var) {
  int n=0,p=0;
  Bound *bound=MRI_var->Bound1;
  MRI_var->Bound2=NULL;
  while (bound) {
    p++;
    Merge(bound->x,bound->y,bound->z,bound->val,&n,MRI_var);
    bound=bound->next;
  }
  while (MRI_var->Bound1) {
    bound=MRI_var->Bound1;
    MRI_var->Basin[bound->z][bound->y][bound->x].type=4;
    MRI_var->Bound1=MRI_var->Bound1->next;
    free(bound);
  }
  MRI_var->Bound1=MRI_var->Bound2;

  return n;
}
#if 0
static int Mediane(int i,int j,int k,int rang) {
  int u,v,w,p,q,r;
  static unsigned char tab[27];

  p=0;
  for (u=-1;u<2;u++)
    for (v=-1;v<2;v++)
      for (w=-1;w<2;w++)
        tab[p++]=MRIvox(MRI_var->mri_src,i+u,j+v,k+w);

  for (q=26;q>1;q--)
    for (p=0;p<q;p++)
      if (tab[p]>tab[p+1]) {
        r=tab[p+1];
        tab[p+1]=tab[p];
        tab[p]=r;
      }
  MRIvox(MRI_var->mri_src,i,j,k)=tab[rang];

  return 0;
}
#endif

#if 0
static int Ambiguous(Cell* cell) {

  if (!parms->watershed_analyze)
    return 0;

  /*Add some code here if you want to take
    into account the ambiguous number*/

  return 1;
}
#endif

static int TRY(int i,int j, int k,MRI_variables *MRI_var) {
  unsigned char val;
  int a,b,c,n;

  n=0;
  for (a=-1;a<=1;a++)
    for (b=-1;b<=1;b++)
      for (c=-1;c<=1;c++)
        if (MRI_var->Basin[k+c][j+b][i+a].type==4 ||
            (MRI_var->Basin[k+c][j+b][i+a].next==MRI_var->Basin[k][j][i].next &&
             MRI_var->Basin[k+c][j+b][i+a].type==7)
            || (MRI_var->Basin[k+c][j+b][i+a].next==MRI_var->Basin[k][j][i].next &&
                MRI_var->Basin[k+c][j+b][i+a].type==8)) {
          val=MRIvox(MRI_var->mri_src,i+a,j+b,k+c);
          if (val>=MRI_var->WM_HALF_MIN  && val<=MRI_var->WM_HALF_MAX)
            n++;
        }

  if (n>=9)
    return 1;
  else
    return 0;
}

static int TRYMERGE(int i,int j, int k,MRI_variables *MRI_var) {
  int a,b,c,n;

  n=0;
  for (a=-1;a<=1;a++)
    for (b=-1;b<=1;b++)
      for (c=-1;c<=1;c++)
        if (MRI_var->Basin[k+c][j+b][i+a].next==MRI_var->Basin[k][j][i].next)
          if (MRI_var->Basin[k+c][j+b][i+a].type==7 || MRI_var->Basin[k+c][j+b][i+a].type==8)
            n++;

  if (n>=6)
    return 1;
  else
    return 0;
}


static int PostAnalyze(STRIP_PARMS *parms,MRI_variables *MRI_var) {
  int i,j,k,p,q;
  BUFTYPE *pb;
  Cell *cell,*cell1=NULL,*celltp=NULL;
  Bound* buff;
  BasinCell* bcell;
  int val;
  unsigned long n=0,added_size=0;

  /*if main_basin_size<estimated_size/5
    find a complement limited to 2*estimated size...*/


  if (MRI_var->main_basin_size<MRI_var->estimated_size/5) {
    fprintf(stderr,"\n      Main basin size probably too small..."
            "looking for complement...");
    p= int(MRI_var->rad_Brain/5.);
    n=0;
    for (i=int(MRI_var->xCOG-p);i<int(MRI_var->xCOG+p);i++)
      for (j=int(MRI_var->yCOG-2*p);j<int(MRI_var->yCOG);j++)
        for (k=int(MRI_var->zCOG-p);k<int(MRI_var->zCOG+p);k++) {
          cell=(Cell*)MRI_var->Basin[k][j][i].next;
          if (cell) {
            if (cell->type==1)
              cell1=TypeVoxel(cell);
            else
              cell1=cell;
            if (cell1->type==2 && ((BasinCell*)cell1->next)->size>n &&
                ((BasinCell*)cell1->next)->size<MRI_var->estimated_size*2) {
              n=((BasinCell*)cell1->next)->size;
              celltp=cell1;
            }
          }
        }
    if (n) {
      celltp->type=3;   /*label this basin type 3*/
      fprintf(stderr,"OK\n");
      MRI_var->main_basin_size+=((BasinCell*)celltp->next)->size;
      fprintf(stderr,"      Corrected main basin size = %3ld\n",MRI_var->main_basin_size);
    } else
      fprintf(stderr,"Could not find a correcting basin\n");
  }


  /*Post-analyze: depends if the mode watershed_analyze is on
    if yes: analyze all the type 2 basins
    if no: all the type 2 basin are labelled type 0 and freed */


  for (k=2;k<MRI_var->depth-2;k++)
    for (j=2;j<MRI_var->height-2;j++) {
      pb=&MRIvox(MRI_var->mri_src,0,j,k);
      pb+=2;
      for (i=2;i<MRI_var->width-2;i++) {
        cell=&MRI_var->Basin[k][j][i];

        if (cell->type==1)
          cell1=TypeVoxel(cell);
        else
          cell1=cell;


        /*right now: only type 0="empty", 1=node to 2 or 3,2=auxiliar basin ,3=main basin to be freed */

        switch (cell1->type) {
          /*cases brain basin*/
        case 3:
          cell->type=4;
          cell1->type=4;
          free((BasinCell*)cell1->next);
          cell->next=(unsigned char*) MRI_var->intbasin+*pb;
          break ;
        case 4:
          cell->type=4;
          cell->next=(unsigned char*) MRI_var->intbasin+*pb;
          break ;
          /*case non brain basin*/
        case 0:
          cell->type=0;
          cell->next=(unsigned char*) MRI_var->intbasin+*pb;
          break;
          /*case ambiguous basin*/
        case 2:
          if (parms->watershed_analyze) {
            if (1)      /*Ambiguous(cell1) instead of 1
                                                                                 if you want to use the ambiguous variable*/
            {
              if (cell!=cell1) {
                cell->type=5;
                cell->next=(Cell*)cell1;
              };
              cell1->type=6;
              ((BasinCell*)cell1->next)->ambiguous=0;
            } else {
              cell->type=0;
              cell1->type=0;
              free((BasinCell*)cell1->next);
              cell->next=(unsigned char*) MRI_var->intbasin+*pb;
            }
          } else /*case non watershed_analyze*/
          {
            cell->type=0;
            cell1->type=0;
            free((BasinCell*)cell1->next);
            cell->next=(unsigned char*) MRI_var->intbasin+*pb;
          }
          break;
        case 5:  /*necessary distinction between type 5 and 6
                                                           to post free the type 6 basin */
          cell->type=5;
          cell->next=(Cell*)cell1->next;
          break;
        case 6:
          if (cell!=cell1) {
            cell->type=5;
            cell->next=(Cell*)cell1;
          }
          break;
        }
        pb++;
      }
    }


  /*right now: only type 0="empty", 4= main basin but "empty"
    , 5=node to 6, 6=auxiliar basin to be freed*/


  if (parms->watershed_analyze) /*type 4, type 5 & 6, and type 0 */
  {
    n=0;
    q=0;
    do {
      p=0;
      for (k=2;k<MRI_var->depth-2;k++)
        for (j=2;j<MRI_var->height-2;j++)
          for (i=2;i<MRI_var->width-2;i++) {
            cell=&MRI_var->Basin[k][j][i];
            if (cell->type==5) {
              cell1=(Cell*)cell->next;
              bcell=(BasinCell*)cell1->next;
              if (TRY(i,j,k,MRI_var))
                cell->type=7;
            } else if (cell->type==6) {
              bcell=(BasinCell*)cell->next;
              if (TRY(i,j,k,MRI_var))
                cell->type=8;
            }
          }

      /*right now: type 0 "empty", 4 "empty", 5 node to 6 or 8,
        6 to be freed , 7 node to 6 or 8 , 8 to be freed */

      for (k=2;k<MRI_var->depth-2;k++)
        for (j=2;j<MRI_var->height-2;j++)
          for (i=2;i<MRI_var->width-2;i++) {
            cell=&MRI_var->Basin[k][j][i];
            if (cell->type==7) {
              cell1=(Cell*)cell->next;
              if (cell1->type==9)
                cell->type=9;
              else {
                /* cell 7 is pointing on type 6 or 8*/
                bcell=(BasinCell*)cell1->next;
                if (TRYMERGE(i,j,k,MRI_var)) {
                  bcell->ambiguous++;
                  if (bcell->ambiguous>=parms->threshold_analyze
                      *pow(bcell->size,1.0/3.0)/100) {
                    fprintf(stderr,"\n      ambiguous basin, merged: at least %ld ambiguous voxels; size: %ld voxels",  bcell->ambiguous, bcell->size);
                    p++;
                    n++;
                    cell->type=9;
                    cell1->type=9;
                    added_size+=bcell->size;
                    free(bcell);
                  }
                }
              }
            } else if (cell->type==8) {
              bcell=(BasinCell*)cell->next;
              if (TRYMERGE(i,j,k,MRI_var)) {
                bcell->ambiguous++;
                if (bcell->ambiguous>=parms->threshold_analyze*pow(bcell->size,1.0/3.0)/100) {
                  fprintf(stderr,"\n      ambiguous basin, merged: at least %ld ambiguous voxels; size: %ld voxels",  bcell->ambiguous, bcell->size);
                  p++;
                  n++;
                  cell->type=9;
                  added_size+=bcell->size;
                  free(bcell);
                }
              }
            }
          }
      q++;
    } while (p);

    fprintf(stderr,"\n      ***** %d basin(s) merged in %d iteration(s)"
            "\n      ***** %ld voxel(s) added to the main basin\n",(int)n,q,added_size);

  }


  for (k=2;k<MRI_var->depth-2;k++)
    for (j=2;j<MRI_var->height-2;j++)
      for (i=2;i<MRI_var->width-2;i++) {
        switch (MRI_var->Basin[k][j][i].type) {
        case 4:
          MRI_var->Basin[k][j][i].type=3;
          break;
        case 5:
          cell1=(Cell*)MRI_var->Basin[k][j][i].next;
          if (cell1->type==3 || cell1->type==9)
            MRI_var->Basin[k][j][i].type=3;
          else
            MRI_var->Basin[k][j][i].type=0;
          break;
        case 6:
          MRI_var->Basin[k][j][i].type=0;
          if (((BasinCell*)MRI_var->Basin[k][j][i].next)->ambiguous)
            fprintf(stderr,"      ambiguous basin, non merged: %ld ambiguous voxels; size: %ld voxels\n",
                    ((BasinCell*)MRI_var->Basin[k][j][i].next)->ambiguous,
                    ((BasinCell*)MRI_var->Basin[k][j][i].next)->size);
          free((BasinCell*)MRI_var->Basin[k][j][i].next);
          break;
        case 7:
          cell1=(Cell*)MRI_var->Basin[k][j][i].next;
          if (cell1->type==3 || cell1->type==9)
            MRI_var->Basin[k][j][i].type=3;
          else
            MRI_var->Basin[k][j][i].type=0;
          break;
        case 8:
          MRI_var->Basin[k][j][i].type=0;
          fprintf(stderr,"      ambiguous basin, non merged: %ld ambiguous voxels; size: %ld voxels\n",
                  ((BasinCell*)MRI_var->Basin[k][j][i].next)->ambiguous,
                  ((BasinCell*)MRI_var->Basin[k][j][i].next)->size);
          free((BasinCell*)MRI_var->Basin[k][j][i].next);
#if __GNUC__ >= 8
  [[gnu::fallthrough]];
#endif
        case 9:
          MRI_var->Basin[k][j][i].type=3;
          break;

        }
      }

  /* seems to be useless in this combined approach.
     However, could be relevant in a single watershed process.*/

  if (!parms->template_deformation) {
    MRI_var->Bound1=NULL;
    n=0;
    for (k=2;k<MRI_var->depth-2;k++)
      for (j=2;j<MRI_var->height-2;j++)
        for (i=2;i<MRI_var->width-2;i++) {
          cell=&MRI_var->Basin[k][j][i];
          if (cell->type) {
            val=AroundCell(i,j,k,MRI_var);
            if (val) {
              n++;
              cell->type=4;
              cell->next=(unsigned char*) MRI_var->intbasin+val;
              if (val>MRI_var->CSF_intensity) {
                buff=(Bound*)malloc(sizeof(Bound));
                buff->x=i;
                buff->y=j;
                buff->z=k;
                buff->val=val;
                buff->next=MRI_var->Bound1;
                MRI_var->Bound1=buff;
              }
            }
          }
        }

    for (i=0;i<3;i++)
      AddVoxel(MRI_var);
  }

  return 0;
}


/*free the allocated Basin (in the routine Allocation)*/
static int FreeMem(MRI_variables *MRI_var) {
  int k,j;

  for (k=0;k<MRI_var->depth;k++) {
    for (j=0;j<MRI_var->height;j++)
      free(MRI_var->Basin[k][j]);
    free(MRI_var->Basin[k]);
  }
  free(MRI_var->Basin);
  return 0;
}

/************************************************************************
 ***********************************************************************
 ************************************************************************/
/*-----------------------------------------------------
        FUNCTION Template_Deformation

        Parameters:
    STRIP_PARMS *:contains the parameters for the prog
    MRI_variables *: contains the variables

        Returns value:void

        Description: the different template deformations
------------------------------------------------------*/
static void Template_Deformation(STRIP_PARMS *parms,MRI_variables *MRI_var) {
  fprintf(stderr,"\n********************TEMPLATE DEFORMATION********************");

  read_geometry(0,MRI_var,NULL);
  brain_params(MRI_var);
  ////////////////////////////////////////////////////////////////////////
  if (parms->cx!=-1) {
    MRI_var->xCOG=parms->cx;
    MRI_var->yCOG=parms->cy;
    MRI_var->zCOG=parms->cz;
    fprintf(stderr,"\n      modification of the brain COG: x=%d y=%d z=%d",
            (int)MRI_var->xCOG,(int)MRI_var->yCOG,(int)MRI_var->zCOG);
  }
  if (parms->rb!=-1) {
    MRI_var->rad_Brain=parms->rb;
    fprintf(stderr,"\n      modification of the brain radius to %d",(int)MRI_var->rad_Brain);
  }
  char *realRAS = getenv("REAL_RAS");
  if (realRAS) {
    MyworldToVoxel = &MRIworldToVoxel;
    MyvoxelToWorld = &MRIvoxelToWorld;
    MRI_var->mris->useRealRAS=1;
    printf("\nINFO: surface is saved with real RAS values\n");
  } else {
    MyworldToVoxel = &MRIsurfaceRASToVoxel;
    MyvoxelToWorld = &MRIvoxelToSurfaceRAS;
    MRI_var->mris->useRealRAS=0;
    printf("\nINFO: surface is saved with conformed RAS with c_(r,a,s) = 0\n");
  }
  ////////////////////////////////////////////////////////////////////////
  init_surf_to_image(0.8*MRI_var->rad_Brain,0.8*MRI_var->rad_Brain,0.8*MRI_var->rad_Brain,MRI_var);
  //init_surf_to_image(1.2*MRI_var->rad_Brain,1.2*MRI_var->rad_Brain,1.2*MRI_var->rad_Brain,MRI_var);

  fprintf(stderr,"\n\n      smoothing...   ");
  //////////////////////////////////////////////////////////////////////////////////////
  // using smaller brain radius to expand to the surface
  MRISfit(MRI_var, calcForce1, parms->forceParam);
  char timeString[256];
  getTimeString(timeString);
  char filename[256];
  if (parms->surf_dbg) {
    int req = snprintf(filename, 256, "surface1-%s-%.2f", timeString, parms->forceParam);
    if (req >= 256) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRISwrite(MRI_var->mris, filename);
  }
  fprintf(stderr, "\n      calcForceGM... ");
  init_direction(MRI_var);
  MRISfit(MRI_var, calcForceGM, parms->forceParam);
  if (parms->surf_dbg) {
    int req = snprintf(filename, 256, "surfaceFinal-%s-%.2f", timeString, parms->forceParam);
    if (req >= 256) {
      std::cerr << __FUNCTION__ << ": Truncation on line " << __LINE__ << std::endl;
    }
    MRISwrite(MRI_var->mris, filename);
  }
  //
  FreeMem(MRI_var);  /*necessary to free the basins previously allocated*/
}

/*load a geometry from a file into mris*/
static void read_geometry(int type,MRI_variables *MRI_var,char *surf_fname) {
  char fname[100], *mri_dir;

  mri_dir = getenv("MRI_DIR");
  switch (type) {
  case 0:
    sprintf(fname,"%s/lib/bem/ic4.tri",mri_dir);
    MRI_var->mris=MRISread(fname);
    break;
  case 1:
    sprintf(fname,"%s/lib/bem/ic5.tri",mri_dir);
    MRI_var->mris=MRISread(fname);
    break;
  case 2:
    MRI_var->mris=MRISread(surf_fname);
    break;
  }
  if (!MRI_var->mris)
    ErrorExit(ERROR_NOFILE, "Could not open the file %s\n",surf_fname) ;

  MRIScomputeNormals(MRI_var->mris);
}

static void brain_params(MRI_variables *MRI_var) {
  int i,j,k,xmin,xmax,ymin,zmin,zmax;
  unsigned long n;
  // BUFTYPE *pb;
  double x,y,z;

  x=y=z=0;
  n=0;
  for (k=0;k<MRI_var->depth;k++)
    for (j=0;j<MRI_var->height;j++) {
      // pb=&MRIvox(MRI_var->mri_src,0,j,k);
      for (i=0;i<MRI_var->width;i++)
        if (MRI_var->Basin[k][j][i].type) {
          x+=i;
          y+=j;
          z+=k;
          n++;
        }
    }

  if (n==0)
    Error("\n pbm of COG calculation");

  MRI_var->xCOG=x/n;
  MRI_var->yCOG=y/n;
  MRI_var->zCOG=z/n;

  n=0;
  xmin = MRI_var->width;
  xmax = 0;
  ymin = MRI_var->height;
  zmin = MRI_var->depth;
  zmax = 0;
  for (k=0;k<MRI_var->depth;k++)
    for (j=0;j<MRI_var->height;j++) {
      // pb=&MRIvox(MRI_var->mri_src,0,j,k);
      for (i=0;i<MRI_var->width;i++)
        if (MRI_var->Basin[k][j][i].type) {
          if (xmin>i) xmin=i;
          if (xmax<i) xmax=i;
          if (ymin>j) ymin=j;
          if (zmin>k) zmin=k;
          if (zmax<k) zmax=k;
        }
    }

  xmax=int(MAX(xmax-MRI_var->xCOG,MRI_var->xCOG-xmin));
  zmax=int(MAX(zmax-MRI_var->zCOG,MRI_var->zCOG-zmax));
  ymin=int(MRI_var->yCOG-ymin);

  MRI_var->rad_Brain=MAX(xmax,MAX(zmax,ymin));

  if (MRI_var->rad_Brain<30) {
    if (MRI_var->WM_INTENSITY==110 && MRI_var->WM_VARIANCE==5 && MRI_var->WM_MAX==110)
      Error("\n Watershed Error !\n");
    else {
      fprintf(stderr,"\n      second estimation of the COG coord: x=%d,y=%d, z=%d, r=%d",(int)MRI_var->xCOG,(int)MRI_var->yCOG,(int)MRI_var->zCOG,(int)MRI_var->rad_Brain);
      Error("\n Watershed Error... Try with the T1-weighted volume\n");
    }
  }
  fprintf(stderr,"\n      second estimation of the COG coord: x=%d,y=%d, z=%d, r=%d",(int)MRI_var->xCOG,(int)MRI_var->yCOG,(int)MRI_var->zCOG,(int)MRI_var->rad_Brain);
}

static void
init_surf_to_image(float rx, float ry, float rz,MRI_variables *MRI_var) {

  double x,y,z;
  MyvoxelToWorld(MRI_var->mri_src,MRI_var->xCOG,MRI_var->yCOG,MRI_var->zCOG
                 ,&x,&y,&z);

  MRISscaleThenTranslate(MRI_var->mris, rx, ry, rz, x, y, z);
}

static void write_image(MRI_variables *MRI_var, int val) {
  int i,j,imnr,k,u,v;
  float x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax;
  float px0,py0,pz0,px1,py1,pz1,px,py,pz;
  int numu,numv;
  double tx,ty,tz;

  MRIS *mris=MRI_var->mris;

  for (k=0;k<mris->nfaces;k++) {
    x0 =mris->vertices[mris->faces[k].v[0]].x;
    y0 =mris->vertices[mris->faces[k].v[0]].y;
    z0 =mris->vertices[mris->faces[k].v[0]].z;
    x1 =mris->vertices[mris->faces[k].v[1]].x;
    y1 =mris->vertices[mris->faces[k].v[1]].y;
    z1 =mris->vertices[mris->faces[k].v[1]].z;
    x2 =mris->vertices[mris->faces[k].v[2]].x;
    y2 =mris->vertices[mris->faces[k].v[2]].y;
    z2 =mris->vertices[mris->faces[k].v[2]].z;
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;
    numu = int(ceil(2*d0));
    numv = int(ceil(2*dmax));


    for (v=0;v<=numv;v++) {
      px0 = x0 + (x2-x0)*v/numv;
      py0 = y0 + (y2-y0)*v/numv;
      pz0 = z0 + (z2-z0)*v/numv;
      px1 = x1 + (x2-x1)*v/numv;
      py1 = y1 + (y2-y1)*v/numv;
      pz1 = z1 + (z2-z1)*v/numv;
      for (u=0;u<=numu;u++) {
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;

        MyworldToVoxel(MRI_var->mri_orig,px,py,pz,&tx,&ty,&tz);

        imnr=(int)(tz+0.5);
        j=(int)(ty+0.5);
        i=(int)(tx+0.5);


        if (i>=0 && i<MRI_var->width && j>=0 && j<MRI_var->height && imnr>=0 && imnr<MRI_var->depth)
          MRIvox(MRI_var->mri_dst,i,j,imnr) = val;
      }
    }
  }
  // vertices is 80
  for (k=0;k<mris->nvertices;k++) {
    px=mris->vertices[k].x;
    py=mris->vertices[k].y;
    pz=mris->vertices[k].z;

    MyworldToVoxel(MRI_var->mri_orig,px,py,pz,&tx,&ty,&tz);

    imnr=(int)(tz+0.5);
    j=(int)(ty+0.5);
    i=(int)(tx+0.5);


    if (i>=0 && i<MRI_var->width && j>=0 && j<MRI_var->height && imnr>=0 && imnr<MRI_var->depth)
      MRIvox(MRI_var->mri_dst,i,j,imnr) = 80;
  }
}

static void init_direction(MRI_variables *MRI_var) {
  int i,j,k,p=0;
  float norm;

  for (i=-1;i<2;i++)
    for (j=-1;j<2;j++)
      for (k=-1;k<2;k++)
        if (i || j || k) {
          norm=sqrt(double(SQR(i)+SQR(j)+SQR(k)));
          MRI_var->direction[p][0]=i/norm;
          MRI_var->direction[p][1]=j/norm;
          MRI_var->direction[p++][2]=k/norm;
        }
}

/* Find 2 normals to a  vector nx, ny, nz */
static void find_normal(float nx,float ny, float nz,float* n1,float *n2,  float direction[26][3]) {
  float ps,ps_buff;
  int p,k;

  k=0;
  ps=10;
  for (p=0;p<26;p++) {
    ps_buff=direction[p][0]*nx+direction[p][1]*ny+direction[p][2]*nz;
    if (ps_buff<0)
      ps_buff=-ps_buff;
    if (ps_buff<ps) {
      ps=ps_buff;
      k=p;
    }
  }
  n1[0]=direction[k][0];
  n1[1]=direction[k][1];
  n1[2]=direction[k][2];

  n2[0]=ny*n1[2]-nz*n1[1];
  n2[1]=nz*n1[0]-nx*n1[2];
  n2[2]=nx*n1[1]-ny*n1[0];

  ps=sqrt(SQR(n2[0])+SQR(n2[1])+SQR(n2[2]));

  if (ps==0)
    Error("\n pbm in find normal ");

  n2[0]/=ps;
  n2[1]/=ps;
  n2[2]/=ps;

}


// balloon by h
static unsigned long MRISpeelBrain(float h,MRI* mri_dst,MRIS *mris,unsigned char val) {
  int i,j,k,imnr;
  float x0,y0,z0,x1,y1,z1,x2,y2,z2,d0,d1,d2,dmax,u,v;
  float px,py,pz,px0,py0,pz0,px1,py1,pz1;
  int numu,numv,totalfilled,newfilled;
  double tx,ty,tz;
  int brainsize;

  // allocate volume (initialized to be all zeros)
  int const width  = mri_dst->width;
  int const height = mri_dst->height;
  int const depth  = mri_dst->depth;

  MRI* const mri_buff = MRIalloc(width, height, depth, MRI_UCHAR) ;

  float *savedX, *savedY, *savedZ;
  MRISexportXYZ(mris, &savedX, &savedY, &savedZ);
  
  MRISblendXYZandNXYZ(mris, float(h));

  //
  for (k=0;k<mris->nfaces;k++) {
    // get three vertices of a face
    x0 =mris->vertices[mris->faces[k].v[0]].x;
    y0 =mris->vertices[mris->faces[k].v[0]].y;
    z0 =mris->vertices[mris->faces[k].v[0]].z;

    x1 =mris->vertices[mris->faces[k].v[1]].x;
    y1 =mris->vertices[mris->faces[k].v[1]].y;
    z1 =mris->vertices[mris->faces[k].v[1]].z;

    x2 =mris->vertices[mris->faces[k].v[2]].x;
    y2 =mris->vertices[mris->faces[k].v[2]].y;
    z2 =mris->vertices[mris->faces[k].v[2]].z;

    // calculate the side lengths
    d0 = sqrt(SQR(x1-x0)+SQR(y1-y0)+SQR(z1-z0));
    d1 = sqrt(SQR(x2-x1)+SQR(y2-y1)+SQR(z2-z1));
    d2 = sqrt(SQR(x0-x2)+SQR(y0-y2)+SQR(z0-z2));
    // get the max length
    dmax = (d0>=d1&&d0>=d2)?d0:(d1>=d0&&d1>=d2)?d1:d2;

    // just pick the division
    numu = (int) ceil(2*d0);
    numv = (int) ceil(2*dmax);

    for (v=0;v<=numv;v++) {
      // PX0 spans the line from X0 to X2
      px0 = x0 + (x2-x0)*v/numv;
      py0 = y0 + (y2-y0)*v/numv;
      pz0 = z0 + (z2-z0)*v/numv;

      // PX1 spans the line from X2 to X1
      px1 = x1 + (x2-x1)*v/numv;
      py1 = y1 + (y2-y1)*v/numv;
      pz1 = z1 + (z2-z1)*v/numv;

      for (u=0;u<=numu;u++) {
        // PX spans the line from PX0 to PX1
        px = px0 + (px1-px0)*u/numu;
        py = py0 + (py1-py0)*u/numu;
        pz = pz0 + (pz1-pz0)*u/numu;

        // get the value of voxel coords
        MyworldToVoxel(mri_dst,px,py,pz,&tx,&ty,&tz);

        imnr=(int)(tz+0.5);
        j=(int)(ty+0.5);
        i=(int)(tx+0.5);
        // if inside the volume, mark 255
        // i.e. that particular voxel is touching the triangle.
        if (i>=0 && i<width && j>=0 && j<height && imnr>=0 && imnr<depth)
          MRIvox(mri_buff,i,j,imnr) = 255;
      }
    }
  }

  MRIvox(mri_buff,1,1,1)= 64; // starting point
  totalfilled = newfilled = 1;

  for (int m=0; m < depth; ++m) {
    if (mri_buff->zi[m] != m)
      printf("mri_buf->zi[m] %d is different from %d\n", mri_buff->zi[m], m);
  }

  while (newfilled>0) {
    newfilled = 0;
    // going from left to right
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++)
          // if voxel == 0, then
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,mri_buff->zi[k-1])==64||MRIvox(mri_buff,i,mri_buff->yi[j-1],k)==64||
                MRIvox(mri_buff,mri_buff->xi[i-1],j,k)==64) {
              // nearby ones are 64, then 64
              MRIvox(mri_buff,i,j,k)= 64;
              newfilled++;
            }
    // going from right to left
    for (k=depth-1;k>=0;k--)
      for (j=height-1;j>=0;j--)
        for (i=width-1;i>=0;i--)
          // if voxel == 0, then
          if (MRIvox(mri_buff,i,j,k)==0)
            if (MRIvox(mri_buff,i,j,mri_buff->zi[k+1])==64||MRIvox(mri_buff,i,mri_buff->yi[j+1],k)==64||
                MRIvox(mri_buff,mri_buff->xi[i+1],j,k)==64) {
              MRIvox(mri_buff,i,j,k) = 64;
              newfilled++;
            }
    totalfilled += newfilled;
  }
  brainsize=0;
  // if zero
  if (val==0) {
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++) {
          if (MRIvox(mri_buff,i,j,k)==64)
            MRIvox(mri_dst,i,j,k) = 0;
          else                             // non 64 counts as brainsize
            brainsize++;
        }
  } else {
    for (k=0;k<depth;k++)
      for (j=0;j<height;j++)
        for (i=0;i<width;i++) {
          if (MRIvox(mri_buff,i,j,k)!=64)
            MRIvox(mri_dst,i,j,k) = val; // set to this value
          else
            brainsize++;
        }
  }

  MRISimportXYZ(mris, savedX, savedY, savedZ);
  freeAndNULL(savedX);
  freeAndNULL(savedY);
  freeAndNULL(savedZ);

  MRIScomputeNormals(mris);

  free(mri_buff);

  int tot = brainsize + totalfilled;
  int volume = width*height*depth;

  if (tot != volume)
    printf("\nError non-labeled voxel exists\n");

  return brainsize;
}

/*Initialize 26 vectors in "each direction" */

static void mean(float tab[4][9],float *moy) {
  int p;
  for (p=0;p<4;p++)
    moy[p]=(2*tab[p][4]+tab[p][1]+tab[p][3]+tab[p][5]+
            tab[p][7])/6;
}

// move vertex with the neighbor average
static void MRISsmooth_surface(MRI_SURFACE *mris,int niter) {

  int const nvertices = mris->nvertices;

  float *p0x, *p0y, *p0z;
  MRISexportXYZ(mris, &p0x, &p0y, &p0z);

  float *p1x, *p1y, *p1z;
  MRISmemalignNFloats(nvertices, &p1x, &p1y, &p1z);

  int iter;
  for (iter=0; iter < niter; iter++) {

    int k;
    for (k=0; k < nvertices; k++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      
      float x = 0, y = 0, z = 0;
      int m;
      for (m=0; m < vt->vnum; m++) {
        x += p0x[vt->v[m]];
        y += p0y[vt->v[m]];
        z += p0z[vt->v[m]];
      }
      // average
      x /= vt->vnum;
      y /= vt->vnum;
      z /= vt->vnum;
      
      // the next value is the average
      p1x[k] = (p0x[k] + x) * 0.5f;
      p1y[k] = (p0y[k] + y) * 0.5f;
      p1z[k] = (p0z[k] + z) * 0.5f;
    }

    // make the old values the new values
    float* t; 
    t = p1x; p1x = p0x; p0x = t;
    t = p1y; p1y = p0y; p0y = t;
    t = p1z; p1z = p0z; p0z = t;
  }

  MRISimportXYZ(mris, p0x, p0y, p0z);

  MRIScomputeMetricProperties(mris) ;
  
  freeAndNULL(p0x); freeAndNULL(p1x);
  freeAndNULL(p0y); freeAndNULL(p1y);
  freeAndNULL(p0z); freeAndNULL(p1z);
}


// shrink surface by h
static void MRISshrink_surface(MRIS *mris,int h) {
  MRISsaveVertexPositions(mris,TMP_VERTICES);
  
  MRISblendXYZandNXYZ(mris, -float(h));

  MRIScomputeNormals(mris);
}


static void MRIVfree(MRI_variables *MRI_var) {
  if (MRI_var->mris)
    MRISfree(&MRI_var->mris);
  if (MRI_var->mrisphere)
    MRISfree(&MRI_var->mrisphere);
  if (MRI_var->mris_curv)
    MRISfree(&MRI_var->mris_curv);
  if (MRI_var->mris_var_curv)
    MRISfree(&MRI_var->mris_var_curv);
  if (MRI_var->mris_dCOG)
    MRISfree(&MRI_var->mris_dCOG);
  if (MRI_var->mris_var_dCOG)
    MRISfree(&MRI_var->mris_var_dCOG);

  free(MRI_var);
}


/*to get the Outer Skin*/
static void MRISshrink_Outer_Skin(MRI_variables *MRI_var,MRI* mri_src) {
  float x,y,z,sx,sy,sz,sd,sxn,syn,szn,sxt,syt,szt,nc;
  float force,force1;

  float d,dx,dy,dz,nx,ny,nz;
  int iter,k,m,n;
  float samp_mean[4];
  float test_samp[4][9];
  int a,b;
  int it,jt,kt,h,niter;

  float decay=0.8,update=0.9;
  float fzero;

  MRIS *mris;
  double tx,ty,tz;


  float val;
  int int_smooth=1;

  double lm,d10m[3],d10,f1m,f2m,dm,dbuff;
  float ***dist;
  int nb_GM,nb_TR,nb_GTM;
  float cout,pcout=0,coutbuff,varbuff,mean_sd[10],mean_dist[10];
  float n1[3],n2[3];

  mris=MRI_var->mris;

  ///////////////////////////////////////////////////////////////
  // initialization
  dist = (float ***) malloc(mris->nvertices*sizeof(float**) );

  for ( it = 0; it < mris->nvertices; it++ ) {
    dist[it] = (float**) malloc( 4*sizeof(float*) );
    for ( jt = 0; jt < 4; jt++ ) {
      dist[it][jt] = (float*) calloc( 3, sizeof(float));
    }
  }

  /*should give some correct results*/
  fzero=MIN(MRI_var->WM_INTENSITY/3,MRI_var->CSF_MAX);

  for (k=0;k<mris->nvertices;k++)
    for (m=0;m<4;m++)
      for (n=0;n<3;n++)
        dist[k][m][n]=0;

  for (n=0;n<10;n++) {
    mean_sd[n]=0;
    mean_dist[n]=0;
  }

  niter =int_smooth;
  force = 0.0f ;
  pcout=0;

  for (k=0;k<mris->nvertices;k++) {
    VERTEX * const v = &mris->vertices[k];
    v->odx = 0;
    v->ody = 0;
    v->odz = 0;
  }

  /////////////////////////////////////////////////////////////////
  // iteration starts here
  for (iter=0;niter;iter++) {
    cout = lm = d10 = f1m = f2m = dm = 0;
    for (k=0;k<mris->nvertices;k++) {
      VERTEX * const v = &mris->vertices[k];
      v->tx = v->x;
      v->ty = v->y;
      v->tz = v->z;
    }

    for (k=0;k<mris->nvertices;k++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      VERTEX                * const v  = &mris->vertices         [k];
      x = v->tx;
      y = v->ty;
      z = v->tz;
      nx = v->nx;
      ny = v->ny;
      nz = v->nz;
      sx=sy=sz=sd=0;
      n=0;
      for (m=0;m<vt->vnum;m++) {
        sx += dx =mris->vertices[vt->v[m]].tx - x;
        sy += dy =mris->vertices[vt->v[m]].ty - y;
        sz += dz =mris->vertices[vt->v[m]].tz - z;
        sd += sqrt(dx*dx+dy*dy+dz*dz);
        n++;
      }
      // S is the vector points to the mean neighbor point
      sx = sx/n;
      sy = sy/n;
      sz = sz/n;
      // mean distance
      sd = sd/n;

      lm+=sd;

      nc = sx*nx+sy*ny+sz*nz;

      // normal component
      sxn = nc*nx;
      syn = nc*ny;
      szn = nc*nz;
      // tangential component
      sxt=sx-sxn;
      syt=sy-syn;
      szt=sz-szn;

      v->nc=nc;

      ///////////////////////////////////////////////////
      // force determination
      force1=0.3;

      f1m+=force1;


      /******************************/

      find_normal(nx,ny,nz,n1,n2,MRI_var->direction);
      for (h=0;h<4;h++)
        for (a=-1;a<2;a++)
          for (b=-1;b<2;b++) {
            // get the RAS value
            MyworldToVoxel(MRI_var->mri_orig,(x-nx*h+n1[0]*a+n2[0]*b),
                           (y-ny*h+n1[1]*a+n2[1]*b),
                           (z-nz*h+n1[2]*a+n2[2]*b),&tx,&ty,&tz);
            kt=(int)(tz+0.5);
            jt=(int)(ty+0.5);
            it=(int)(tx+0.5);

            if ((kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height))
              val=0;
            else
              val=MRIvox(MRI_var->mri_orig,it,jt,kt);

            test_samp[h][3*b+a+4] = val;
          }

      val=test_samp[0][4];

      if (!val)
        force=-0.25;
      else if (val<=fzero/2)
        force=-0.1;
      else {
        mean(test_samp,samp_mean);

        if (samp_mean[0]<fzero && samp_mean[1]<fzero)
          force=-0.2;
        else {
          // ??????????????
          nb_GM=0;
          nb_TR=0;
          nb_GTM=0;
          for (h=0;h<4;h++) {
            if (samp_mean[h]>=MRI_var->TRANSITION_INTENSITY)
              nb_GM++;
            if (samp_mean[h]<fzero)
              nb_TR++;
          }

          if (nb_TR>=3)
            force=-0.2;
          else if (nb_GM>=3 && samp_mean[0]>MRI_var->TRANSITION_INTENSITY)
            force=0.7;
          else if (nb_GM==2 && samp_mean[0]>MRI_var->TRANSITION_INTENSITY)
            force=0.5;
          else if (nb_TR==0)
            force=0.3;
          else {
            nb_GM=0;
            nb_TR=0;
            for (h=0;h<4;h++) {
              for (a=0;a<9;a++) {
                if (test_samp[h][a]>=MRI_var->TRANSITION_INTENSITY)
                  nb_GM++;
                else if (test_samp[h][a]<fzero)
                  nb_TR++;
                else
                  nb_GTM++;
              }
            }

            if (nb_TR>=25)
              force=-0.3;
            else if (nb_GM>=18)
              force=0.5;
            else if (nb_GM>=15)
              force=0.3;
            else {
              if (nb_GM>9 && nb_TR<9)
                force=0.5;
              else if (nb_GTM>30)
                force=0.1;
              else
                force=-0.0;
            }
          }
        }
      }

      f2m+=force;

      force1=0.5;

      // Delta = .8 x St + force1 x Sn + force X Vn
      ///////////////////////////////////////////////
      dx = sxt*0.8 + force1*sxn +v->nx*force;
      dy = syt*0.8 + force1*syn +v->ny*force;
      dz = szt*0.8 + force1*szn +v->nz*force;

      // modify Detal
      ///////////////////////////////////////////////
      dx = decay*v->odx+update*dx;
      dy = decay*v->ody+update*dy;
      dz = decay*v->odz+update*dz;

      // if too much, make Detal < 1
      if ((d=sqrt(dx*dx+dy*dy+dz*dz))>1.0) {
        dx /= d;
        dy /= d;
        dz /= d;
      }

      // cache the Delta
      v->odx = dx;
      v->ody = dy;
      v->odz = dz;


      d=sqrt(dx*dx+dy*dy+dz*dz);

      dm+=d;

      dist[k][iter%4][0]=x;
      dist[k][iter%4][1]=y;
      dist[k][iter%4][2]=z;

      d10m[0] = d10m[1] = d10m[2] = 0;

      for (n=0;n<4;n++) {
        d10m[0]+=dist[k][n][0]/4;
        d10m[1]+=dist[k][n][1]/4;
        d10m[2]+=dist[k][n][2]/4;
      }

      dbuff=0;
      for (n=0;n<4;n++)
        dbuff+=SQR(dist[k][n][0]-d10m[0])+SQR(dist[k][n][1]-d10m[1])+
               SQR(dist[k][n][2]-d10m[2]);

      d10+=dbuff/4;

      // move the position
      MRISsetXYZ(mris,k,
        v->x + dx,
        v->y + dy,
        v->z + dz);
    }

    lm /=mris->nvertices;
    f1m /=mris->nvertices;
    f2m /=mris->nvertices;
    dm /=mris->nvertices;
    d10 /=mris->nvertices;

    mean_sd[iter%10]=lm;
    mean_dist[iter%10]=d10;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_sd[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_sd[n]-coutbuff);

    cout=varbuff;

    coutbuff=0;
    for (n=0;n<10;n++)
      coutbuff+=mean_dist[n]/10;

    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_dist[n]-coutbuff);

    cout+=10*varbuff;

    coutbuff=cout;

    cout=(cout+pcout)/2;

    pcout=coutbuff;

    MRIScomputeNormals(mris);

    /*if ((niter==int_smooth) && !(iter % 5))
      fprintf(stderr,
              "%d: lm=%5.3f,f1m=%5.3f,f2m=%5.3f,dm=%5.3f,d10m=%5.3f,c=%5.3f\n"
              ,iter,lm,f1m,f2m,dm,d10,100*cout);*/

    if (niter==int_smooth) {
      if (((iter>20)&&(10000*cout<1))||(iter>200))
        niter--;
    } else
      niter--;
  }
  fprintf(stderr,"%d iterations",iter);
  MRIScomputeNormals(mris);

  /*free memory*/
  for ( it = 0; it < mris->nvertices; it++ ) {
    for ( jt = 0; jt < 4; jt++ )
      free(dist[it][jt]);
    free(dist[it]);
  }
  free(dist);
}

static void label_voxels(STRIP_PARMS *parms, MRI_variables *MRI_var,MRI *mri_with_skull) {
  char fname[512];
  int i,j,k;
  int A,B;
  unsigned long  volume_skull;
  double vol_elt;

  A=(2*MRI_var->CSF_MAX+MRI_var->TRANSITION_INTENSITY)/3;
  B=(MRI_var->GM_INTENSITY+MRI_var->WM_MIN)/2;

  fprintf(stderr,"\n      tissue label process");

  for (i=0;i<MRI_var->width;i++)
    for (j=0;j<MRI_var->height;j++)
      for (k=0;k<MRI_var->depth;k++)
        MRIvox(MRI_var->mri_src,i,j,k)=0;

  /*loading the outer skin surface*/
  strcpy(fname,parms->surfname);
  strcat(fname,"_outer_skin_surface");
  read_geometry(2,MRI_var,fname);
  MRISpeelBrain(-1,MRI_var->mri_src,MRI_var->mris,1);

  /*loading the outer outer skull surface*/
  strcpy(fname,parms->surfname);
  strcat(fname,"_outer_skull_surface");
  read_geometry(2,MRI_var,fname);
  volume_skull=MRISpeelBrain(-1,MRI_var->mri_src,MRI_var->mris,2);
  volume_skull=MRI_var->depth*MRI_var->width*MRI_var->height-volume_skull;
  vol_elt=MRI_var->mri_src->xsize*MRI_var->mri_src->ysize*MRI_var->mri_src->zsize;

  fprintf(stderr,"\n\nSkull Size = %ld voxels, voxel volume = %2.3f mm3\n"
          ,volume_skull,(float)vol_elt);
  fprintf(stderr,"           = %2.0f mmm3 = %2.3f cm3\n"
          ,volume_skull*vol_elt,(float)volume_skull/1000.*vol_elt);


  /*loading the inner skull surface*/
  strcpy(fname,parms->surfname);
  strcat(fname,"_inner_skull_surface");
  read_geometry(2,MRI_var,fname);
  MRISpeelBrain(-1,MRI_var->mri_src,MRI_var->mris,3);

  /*loading the brain surface*/
  strcpy(fname,parms->surfname);
  strcat(fname,"_brain_surface");
  read_geometry(2,MRI_var,fname);
  MRISpeelBrain(-1,MRI_var->mri_src,MRI_var->mris,4);

  for (i=0;i<MRI_var->width;i++)
    for (j=0;j<MRI_var->height;j++)
      for (k=0;k<MRI_var->depth;k++)
        if (MRIvox(MRI_var->mri_src,i,j,k)==4) {
          if (MRIvox(mri_with_skull,i,j,k)<A)
            MRIvox(MRI_var->mri_src,i,j,k)=3;
          else if (MRIvox(mri_with_skull,i,j,k)>B)
            MRIvox(MRI_var->mri_src,i,j,k)=5;
        } else if (MRIvox(MRI_var->mri_src,i,j,k)==2) {
          if (MRIvox(mri_with_skull,i,j,k)>B)
            MRIvox(MRI_var->mri_src,i,j,k)=6;
        };

}

#define CORR_THRESHOLD 5.3f

int finite(double v);

static void MRISchangeCoordinates(MRIS *mris,MRIS *mris_orig) {
  MRIScopyXYZ(mris, mris_orig);
  mris->radius=mris_orig->radius;
  mris->status=mris_orig->status;
}


// used in MRISshrink1
static void calcForce1(
  double &force0, double &force1, double &force,
  double sd,
  const TVector &Pos, const TVector &S, const TVector &N,
  MRI_variables *MRI_var, double forceParam) {
  double force2, force3;

  int ninside=15,noutside=10;

  double x = Pos.x;
  double y = Pos.y;
  double z = Pos.z;
  double nx = N.x;
  double ny = N.y;
  double nz = N.z;

  force0 = 0.8;

  ////////////////////////////////////////////////////////////////
  force1=0.7;

  ////////////////////////////////////////////////////////////////
  // voxel based force calculation
  //
  ////////////////////////////////////////////////////////////////
  // force to shrink from outside
  force2=-1; // push in
  for (int h= -noutside;h<0;h++) {
    double tx, ty, tz;
    // look at outside side voxels (h < 0)
    MyworldToVoxel(MRI_var->mri_orig,(x-nx*h),
                   (y-ny*h),(z-nz*h),&tx,&ty,&tz);
    int kt=(int)(tz+0.5);
    int jt=(int)(ty+0.5);
    int it=(int)(tx+0.5);
    // if inside the bounding box
    if (!(kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height))
      // if this voxel is inside, then don't
      //if (MRI_var->Basin[kt][jt][it].type)
      if (MRIvox(MRI_var->mri_orig, it, jt, kt) != 0)
        force2=0;
  }

  ///////////////////////////////////////////////////////////////
  // force to expand from inside
  force3 = 1; // push out
  for (int h=1;h<ninside;h++) {
    double tx, ty, tz;
    // look at inside voxes (h > 0)
    MyworldToVoxel(MRI_var->mri_orig,(x-nx*h),
                   (y-ny*h),(z-nz*h),&tx,&ty,&tz);
    int kt=(int)(tz+0.5);
    int jt=(int)(ty+0.5);
    int it=(int)(tx+0.5);
    // if outside of the bounding box, then force3 = 0
    if (kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height)
      force3 = 0;
    // if this voxel is not inside, then don't
    // else if (!MRI_var->Basin[kt][jt][it].type)
    else if (MRIvox(MRI_var->mri_orig, it, jt, kt) == 0)
      force3=0;
  }
  ///////////////////////////////////////////////////////////////
  // force = 0.2*force2+1.0*(force3-0.1);
  force = 0.2*force2+forceParam*(force3-0.1);
}

static void calcForceGM(
  double &force0, double &force1, double &force,
  double sd,
  const TVector &Pos, const TVector &S, const TVector &N,
  MRI_variables *MRI_var,
  double forceParam) {
  float samp_mean[4];
  float test_samp[4][9];
  float n1[3],n2[3];
  double val;

  double x = Pos.x;
  double y = Pos.y;
  double z = Pos.z;
  double nx = N.x;
  double ny = N.y;
  double nz = N.z;
  double tx, ty, tz;

  float fzero=MRI_var->CSF_intensity;

  int h,a,b;
  int it, jt, kt;

  ///////////////////////////////////////////////////
  // force determination
  force0 = 0.8;
  force1 = 0.8; //0.5;

  /******************************/
  // 3 x 3 x 4 region for force calculation
  find_normal(nx,ny,nz,n1,n2,MRI_var->direction);
  for (h=0;h<4;h++)
    for (a=-1;a<2;a++)
      for (b=-1;b<2;b++) {
        // get the voxel value
        MyworldToVoxel(MRI_var->mri_orig,(x-nx*h+n1[0]*a+n2[0]*b),
                       (y-ny*h+n1[1]*a+n2[1]*b),
                       (z-nz*h+n1[2]*a+n2[2]*b),&tx,&ty,&tz);
        kt=(int)(tz+0.5);
        jt=(int)(ty+0.5);
        it=(int)(tx+0.5);

        // outside the region
        if ((kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height))
          val=0;
        else
          val=MRIvox(MRI_var->mri_orig,it,jt,kt);

        test_samp[h][3*b+a+4] = val;
      }
  // test_samp[][] now has all the voxel values.
  val=test_samp[0][4];

  int nb_GM, nb_TR, nb_GTM;

  // force is the normal force
  if (!val)
    force=-0.25;
  else if (val<=fzero/2)
    force=-0.1;
  else {
    mean(test_samp,samp_mean);

    if (samp_mean[0]<fzero && samp_mean[1]<fzero)  // dark region.  pushin
      force=-0.2;
    else {
      // count various matters
      nb_GM=0;
      nb_TR=0;
      nb_GTM=0;
      for (h=0;h<4;h++) {
        if (samp_mean[h]>=MRI_var->TRANSITION_INTENSITY)
          nb_GM++;
        if (samp_mean[h]<fzero)
          nb_TR++;
      }

      if (nb_TR>=3)
        force=-0.2;
      else if (nb_GM>=3 && samp_mean[0]>MRI_var->TRANSITION_INTENSITY)  // push out
        force=2.1; // 0.7;
      else if (nb_GM==2 && samp_mean[0]>MRI_var->TRANSITION_INTENSITY)  // push out
        force=1.5; // 0.5;
      else if (nb_TR==0)
        force=0.9; // 0.3;
      else {
        nb_GM=0;
        nb_TR=0;
        for (h=0;h<4;h++) {
          for (a=0;a<9;a++) {
            if (test_samp[h][a]>=MRI_var->TRANSITION_INTENSITY)
              nb_GM++;
            else if (test_samp[h][a]<fzero)
              nb_TR++;
            else
              nb_GTM++;
          }
        }
        if (nb_TR>=25)
          force=-0.3;
        else if (nb_GM>=18)
          force=1.5; // 0.5;
        else if (nb_GM>=15)
          force=.9; // 0.3;
        else {
          if (nb_GM>9 && nb_TR<9)
            force=1.5; // 0.5;
          else if (nb_GTM>30)
            force=0.3; // 0.1;
          else
            force=-0.0;
        }
      }
    }
  }
}


static void calcForceMine(
  double &force0, double &force1, double &force,
  double sd,
  const TVector &Pos, const TVector &S, const TVector &N,
  MRI_variables *MRI_var, double forceParam) {
  int ninside=10,noutside=10;

  double x = Pos.x;
  double y = Pos.y;
  double z = Pos.z;
  double nx = N.x;
  double ny = N.y;
  double nz = N.z;

  int it, jt, kt;
  double tx, ty, tz;

  force0 = 0.8;

  ////////////////////////////////////////////////////////////////
  force1=0.7;

  ////////////////////////////////////////////////////////////////
  // voxel based force calculation
  //
  ////////////////////////////////////////////////////////////////
  // count non-zero voxels in the strip
  int countoutside=0;
  int countinside=0;
  for (int h= -noutside;h<0;h++) {

    double tx, ty, tz;
    // look at outside side voxels (h < 0)
    MyworldToVoxel(MRI_var->mri_orig,(x-nx*h),
                   (y-ny*h),(z-nz*h),&tx,&ty,&tz);
    kt=(int)(tz+0.5);
    jt=(int)(ty+0.5);
    it=(int)(tx+0.5);
    // if inside the bounding box
    if (!(kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height))
      // if this voxel is inside, then don't
      //if (MRI_var->Basin[kt][jt][it].type)
      if (MRIvox(MRI_var->mri_orig, it, jt, kt) != 0)
        countoutside++;
  }

  ///////////////////////////////////////////////////////////////
  // force to expand from inside
  for (int h=1;h<ninside;h++) {
    double tx, ty, tz;
    // look at inside voxes (h > 0)
    MyworldToVoxel(MRI_var->mri_orig,(x-nx*h),
                   (y-ny*h),(z-nz*h),&tx,&ty,&tz);
    kt=(int)(tz+0.5);
    jt=(int)(ty+0.5);
    it=(int)(tx+0.5);
    // if outside of the bounding box, then force3 = 0
    if (!(kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height))
      // if this voxel is not inside, then don't
      // else if (!MRI_var->Basin[kt][jt][it].type)
      if (MRIvox(MRI_var->mri_orig, it, jt, kt) != 0)
        countinside++;
  }
  // get the current position in voxel unit
  MyworldToVoxel(MRI_var->mri_orig,x,y,z, &tx, &ty,&tz);
  kt=(int)(tz+0.5);
  jt=(int)(ty+0.5);
  it=(int)(tx+0.5);
  // if inside the volume
  if (!(kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height)) {
    unsigned char val = MRIvox(MRI_var->mri_orig, it,jt,kt);
    // if non-zero, which way must be pushed?
    //    count number of non-zero voxels.  brighter side should be the case
    // if zero, which way should be pushed?
    //    push to the darker side
    if (val > 0) { // push to brighter side
      force = (countoutside < countinside) ? 1.0 : -1.0;
    } else { // push to the darker side
      force = (countoutside < countinside) ? -1.0 : 1.0;
    }
  } else
    force = 0;
}

// shrink2
static void calcForce2(
  double &force0, double &force1, double &force,
  double sd,
  const TVector &Pos, const TVector &S, const TVector &N,
  MRI_variables *MRI_var, double forceParam) {
  int nb_GM,nb_TR,nb_GTM;
  double r,F,E,rmin=3.33,rmax=10.;
  float n1[3], n2[3];
  int h, a, b;
  double tx, ty, tz;
  int kt, jt, it;
  int val;
  float test_samp[4][9];
  float samp_mean[4];

  double x = Pos.x;
  double y = Pos.y;
  double z = Pos.z;

  double nx = N.x;
  double ny = N.y;
  double nz = N.z;

  /////////////////////////////////////////////////////
  // surface force (tangential direction)
  /////////////////////////////////////////////////////
  force1=0;

  double nc = S*N; // sx*nx+sy*ny+sz*nz;

  E=(1/rmin+1/rmax)/2;
  F=6/(1/rmin-1/rmax);

  if (nc) {
    r= (nc>0) ? nc : -nc;
    r=SQR(sd)/(2*r);
    force1=(1+tanh(F*(1/r-E)))/2;
  } else
    Error("pbm de normale!");

  ////////////////////////////////////////////////////
  // voxel dependent force
  ///////////////////////////////////////////////////
  find_normal(nx,ny,nz,n1,n2,MRI_var->direction);
  //
  // 3x3x4 voxels to look at
  force = 0; // initialize
  for (h=0;h<4;h++)
    for (a=-1;a<2;a++)
      for (b=-1;b<2;b++) {
        MyworldToVoxel(MRI_var->mri_orig,
                       (x-nx*h+n1[0]*a+n2[0]*b),
                       (y-ny*h+n1[1]*a+n2[1]*b),
                       (z-nz*h+n1[2]*a+n2[2]*b),&tx,&ty,&tz);
        kt=(int)(tz+0.5);
        jt=(int)(ty+0.5);
        it=(int)(tx+0.5);

        // outside of the volume
        if ((kt<0||kt>=MRI_var->depth||it<0||it>=MRI_var->width||jt<0||jt>=MRI_var->height))
          val=0;
        else
          val=MRIvox(MRI_var->mri_src,it,jt,kt);

        test_samp[h][3*b+a+4] = val;
      }

  val=(int) test_samp[0][4];

  if (!val)     /*|| val>fmax)*/
    force=-0.25;
  else if (val<=MRI_var->CSF_MAX)
    force=-0.1;
  else if (val<MRI_var->TRANSITION_INTENSITY)
    force=0.0;
  else {
    mean(test_samp,samp_mean);

    if (samp_mean[1]<MRI_var->TRANSITION_INTENSITY &&
        samp_mean[2]<MRI_var->TRANSITION_INTENSITY) {
      if (samp_mean[0]*100>samp_mean[1]*90)
        if (samp_mean[1]*100>samp_mean[2]*90)
          force=-0.1;
    } else {
      nb_GM=0; // count number of grey matter
      nb_TR=0; // count number of transition voxels
      nb_GTM=0;
      for (h=0;h<4;h++) {
        if (samp_mean[h]>=MRI_var->GM_INTENSITY)
          nb_GM++;
        if (samp_mean[h]<MRI_var->TRANSITION_INTENSITY)
          nb_TR++;
      }

      if (nb_TR>=3)
        force=-0.2;
      else if (nb_GM>=3 && samp_mean[0]>MRI_var->TRANSITION_INTENSITY)
        force=0.7;
      else if (nb_GM==2 && samp_mean[0]>MRI_var->TRANSITION_INTENSITY)
        force=0.5;
      else if (nb_TR==0)
        force=0.3;
      else {
        nb_GM=0;
        nb_TR=0;
        for (h=0;h<4;h++) {
          for (a=0;a<9;a++) {
            if (test_samp[h][a]>=MRI_var->GM_INTENSITY)
              nb_GM++;
            else if (test_samp[h][a]<MRI_var->TRANSITION_INTENSITY)
              nb_TR++;
            else
              nb_GTM++;
          }
        }

        if (nb_TR>=18)
          force=-0.3;
        else if (nb_GM>=18)
          force=0.5;
        else if (nb_GM>=15)
          force=0.3;
        else {
          if (nb_GM>9 && nb_TR<9)
            force=0.5;
          else if (nb_GTM>30)
            force=0.1;
          else
            force=-0.0;
        }
      }
    }
  }

  force += tanh(nc*0.1);
}

static void MRISfit(MRI_variables *MRI_var,
                    void (*calcforce)
                    (double &force0,double &force1, double &force2,
                     double sd,
                     const TVector &Pos,
                     const TVector &S,
                     const TVector &N,
                     MRI_variables *MRI_var,
                     double forceParam),
                    double forceParam) {
  double sd, nc;
  double force,force1;
  double force0;
  double d;
  int iter,k,m,n;
  int it, jt;
  int niter;

  double decay=0.8,update=0.9;

  int int_smooth=10;

  MRIS *mris;
  //  char surf_fname[100];


  double lm,d10m[3],d10,f1m,f2m,dm,dbuff;
  float ***dist;
  float cout,cout_prec,coutbuff,varbuff,mean_sd[10],mean_dist[10];


  mris=MRI_var->mris;
  MRIScomputeNormals(mris);

  //////////////////////////////////////////////////////////////
  // initialize vars
  dist = (float ***) malloc( mris->nvertices*sizeof(float**) );

  for (it = 0; it < mris->nvertices; it++ ) {
    dist[it] = (float**) malloc( 4*sizeof(float*) );
    for ( jt = 0; jt < 4; jt++ ) {
      dist[it][jt] = (float*) calloc( 3, sizeof(float));
    }
  }

  for (k=0;k<mris->nvertices;k++)
    for (m=0;m<4;m++)
      for (n=0;n<3;n++)
        dist[k][m][n]=0;

  for (n=0;n<10;n++) {
    mean_sd[n]=0;
    mean_dist[n]=0;
  }

  niter =int_smooth;
  force = 0.0f ;

  cout_prec = 0;

  /* momentum -> 0*/
  for (k=0;k<mris->nvertices;k++) {
    VERTEX * const v = &mris->vertices[k];
    v->odx = 0;
    v->ody = 0;
    v->odz = 0;
  }

  ///////////////////////////////////////////////////////////////
  // iterations
  for (iter=0;niter;iter++) {
    lm = d10 = f1m = f2m = dm = 0;
    for (k=0;k<mris->nvertices;k++) {
      VERTEX * const v = &mris->vertices[k];
      v->tx = v->x;  // initialize t(mp)
      v->ty = v->y;
      v->tz = v->z;
    }

    for (k=0;k<mris->nvertices;k++) {
      VERTEX_TOPOLOGY const * const vt = &mris->vertices_topology[k];
      VERTEX                * const v  = &mris->vertices         [k];      // vertex position
      TVector Pos(v->tx, v->ty, v->tz);
      TVector S(0.,0.,0.);
      sd=0;
      n=0;
      // get the mean position of neighboring vertices
      // try to minimize
      TVector dX;
      for (m=0;m<vt->vnum;m++) {
        TVector Vert(mris->vertices[vt->v[m]].tx,
                     mris->vertices[vt->v[m]].ty,
                     mris->vertices[vt->v[m]].tz);
        dX = Vert - Pos;
        S += dX;
        sd += sqrt(dX*dX);
        n++;
      }
      // S=(sx,sy,sz)  points to the mean position from this particular vertex
      S /=n;
      // mean distance
      sd = sd/n;

      // cache
      lm+=sd;

      // inner product of S and N
      TVector N(v->nx, v->ny, v->nz);
      nc = S*N;

      // Normal component of S
      TVector Sn = nc*N;

      // Tangential component of S
      TVector St = S - Sn;

      v->nc=nc;

      ////////////////////////////////////////////////////////////////
      // calculate force
      calcforce(force0, force1, force, sd, Pos, S, N, MRI_var, forceParam);

      f1m+=force1;
      f2m+=force;

      ///////////////////////////////////////////////////////////////
      // keep tangential vector smaller < 1.0
      if ((d = sqrt(St*St)) > 1.0)
        St /= d;

      //////////////////////////////////////////////////////////////
      // move delta
      //  Delta = force0 * St    (+) force1 * Sn (+) force * Vn
      //          move within    smoothness      voxel dependent force
      //////////////////////////////////////////////////////////////
      // bigger force to smooth
      dX = St*force0 + Sn*force1 + N*force;
      // v->odx (last cached value)
      // decay = .8, update = .9
      TVector V(v->odx, v->ody, v->odz);
      dX = decay*V + update*dX;

      // if too big, make it small < 1.0
      if ((d = sqrt(dX*dX)) > 1.0)
        dX /= d;

      // cache the value
      v->odx = dX.x;
      v->ody = dX.y;
      v->odz = dX.z;

      // calculate the size of the movement
      d = sqrt(dX*dX);

      dm+=d;

      /////////////////////////////////////////////
      dist[k][iter%4][0]=Pos.x;
      dist[k][iter%4][1]=Pos.y;
      dist[k][iter%4][2]=Pos.z;

      d10m[0] = d10m[1] = d10m[2] = 0;

      for (n=0;n<4;n++) {
        d10m[0]+=dist[k][n][0]/4;
        d10m[1]+=dist[k][n][1]/4;
        d10m[2]+=dist[k][n][2]/4;
      }

      dbuff=0;
      for (n=0;n<4;n++)
        dbuff+=SQR(dist[k][n][0]-d10m[0])+SQR(dist[k][n][1]-d10m[1])+SQR(dist[k][n][2]-d10m[2]);

      d10+=dbuff/4;

      ////////////////////////////////////////////////////////////
      // now move vertex by (dx, dy, dz)
      MRISsetXYZ(mris,k,
        v->x + dX.x,
        v->y + dX.y,
        v->z + dX.z);
    }

    lm /=mris->nvertices;  // neighbor distance mean
    f1m /=mris->nvertices; // force1 mean
    f2m /=mris->nvertices; // force  mean
    dm /=mris->nvertices;  // movement mean
    d10 /=mris->nvertices; //

    // save the value
    mean_sd[iter%10]=lm;
    mean_dist[iter%10]=d10;


    coutbuff=0;
    // get the past 10 average sd
    for (n=0;n<10;n++)
      coutbuff+=mean_sd[n]/10;

    // variation square
    varbuff=0;
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_sd[n]-coutbuff);

    cout=varbuff;

    coutbuff=0;
    // get the past 10 average dist
    for (n=0;n<10;n++)
      coutbuff+=mean_dist[n]/10;

    varbuff=0;
    // variation square
    for (n=0;n<10;n++)
      varbuff+=SQR(mean_dist[n]-coutbuff);

    // add both sd and dist variations
    cout+=10*varbuff;

    // save
    coutbuff=cout;

    // get the previous value and the current one averaged
    cout=(cout_prec+cout)/2;

    // save
    cout_prec=coutbuff;

    MRIScomputeNormals(mris);

    /*
    if ((niter==int_smooth) && !(iter % 5))
    {
      char surf_fname[256];
      fprintf(stderr,
        "%d: lm=%5.3f,f1m=%5.3f,f2m=%5.3f,dm=%5.3f,d10m=%5.3f,c=%5.3f\n"
        ,iter,lm,f1m,f2m,dm,d10,100*cout);
      sprintf(surf_fname,"./test/lh.test%d",iter);
      MRISwrite(MRI_var->mris,surf_fname);
    }
    */
    // every 10th iteration check
    if (niter==int_smooth) {
      // demand cout < 5/10000 or iteration < 150
      if (((iter>20)&&(10000*cout<5))||(iter>150)) {
        niter--;
      };
    } else
      niter--;
  }
  fprintf(stderr,"%d iterations",iter);

  /*free memory*/
  for ( it = 0; it < mris->nvertices; it++ ) {
    for ( jt = 0; jt < 4; jt++ )
      free(dist[it][jt]);
    free(dist[it]);
  }
  free(dist);

}

