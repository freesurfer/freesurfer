/*----------------------------------------------------------
  Name: mri_label2label.c
  $Id: mri_label2label.c,v 1.3 2001/05/08 22:58:14 greve Exp $
  Author: Douglas Greve
  Purpose: Converts a label in one subject's space to a label
  in another subject's space using either talairach or spherical
  as an intermediate registration space. 

  Example 1: If you have a label from subject fred called
    broca-fred.label defined on fred's left hemispherical 
    surface and you want to convert it to sally's surface, then

    mri_label2label --srclabel broca-fred.label  --srcsubject fred 
                    --trglabel broca-sally.label --trgsubject sally
                    --regmethod surface --hemi lh

    This will map from fred to sally using sphere.reg. The registration
    surface can be changed with --surfreg.

  Example 2: You could also do the same mapping using talairach 
    space as an intermediate:

    mri_label2label --srclabel broca-fred.label  --srcsubject fred 
                    --trglabel broca-sally.label --trgsubject sally
                    --regmethod volume

    Note that no hemisphere is specified with -regmethod.
  -----------------------------------------------------------*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include "icosahedron.h"
#include "fio.h"

#include "MRIio.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mri.h"
#include "label.h"

#include "registerio.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_label2label.c,v 1.3 2001/05/08 22:58:14 greve Exp $";
char *Progname = NULL;

char  *srclabelfile = NULL;
LABEL *srclabel     = NULL;
char  *srcsubject   = NULL;
char  *trglabelfile = NULL;
LABEL *trglabel     = NULL;
char  *trgsubject   = NULL;
char  *trgsurface   = "white";

char *regmethod  = NULL;
char *hemi       = NULL;
char *surfreg = "sphere.reg";

int srcicoorder = -1;
int trgicoorder = -1;

MRI_SURFACE *SrcSurfReg;
MRI_SURFACE *TrgSurf;
MRI_SURFACE *TrgSurfReg;
MATRIX *SrcVolReg;
MATRIX *TrgVolReg;
MATRIX *InvTrgVolReg;
MATRIX *Src2TrgVolReg;

float IcoRadius = 100.0;

int debug = 0;

char *SUBJECTS_DIR = NULL;
char *MRI_DIR = NULL;
FILE *fp;

char tmpstr[2000];

/*-------------------------------------------------*/
int main(int argc, char **argv)
{
  MATRIX *xyzSrc, *xyzTrg;
  MHT *TrgHash;
  VERTEX *srcvtx, *trgvtx;
  int n,err,srcvtxno,trgvtxno,allzero;
  float dmin;

  printf("\n");

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  /*--- Get environment variables ------*/
  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR==NULL){
    fprintf(stderr,"ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  MRI_DIR = getenv("MRI_DIR") ;
  if(MRI_DIR==NULL){
    fprintf(stderr,"ERROR: MRI_DIR not defined in environment\n");
    exit(1);
  }

  /*--- Load in Source Label ------*/
  printf("Loading source label.\n");
  srclabel = LabelRead(NULL, srclabelfile);
  if(srclabel == NULL){
    fprintf(stderr,"ERROR reading %s\n",srclabelfile);
    exit(1);
  }
  printf("Found %d points in source label.\n",srclabel->n_points);
  fflush(stdout); fflush(stderr);

  /* -- Allocate the Target Label ---*/
  trglabel = LabelAlloc(srclabel->n_points,trgsubject,trglabelfile);
  trglabel->n_points = srclabel->n_points;

  /* Set up vectors */
  xyzSrc = MatrixAlloc(4,1,MATRIX_REAL);
  xyzSrc->rptr[3+1][0+1] = 1.0;
  xyzTrg = MatrixAlloc(4,1,MATRIX_REAL);

  /*--------------------- VOLUMETRIC MAPPING --------------------------*/
  if(!strcmp(regmethod,"volume")){

    printf("Starting volumetric mapping\n");

    /*** Load the Src2Tal registration ***/
    if(strcmp(srcsubject,"talairach")){
      sprintf(tmpstr,"%s/%s/mri/transforms/talairach.xfm",
        SUBJECTS_DIR,srcsubject);
      err = regio_read_mincxfm(tmpstr, &SrcVolReg);
      if(err) {
  fprintf(stderr,"ERROR reading %s\n",tmpstr);
  exit(1);
      }
    }
    else SrcVolReg = MatrixIdentity(4,NULL);

    /*** Load the Trg2Tal registration ***/
    if(strcmp(trgsubject,"talairach")){
      sprintf(tmpstr,"%s/%s/mri/transforms/talairach.xfm",
        SUBJECTS_DIR,trgsubject);

      err = regio_read_mincxfm(tmpstr, &TrgVolReg);
      if(err) {
  fprintf(stderr,"ERROR reading %s\n",tmpstr);
  exit(1);
      }
    }
    else TrgVolReg = MatrixIdentity(4,NULL);

    /* Compte the Src-to-Trg Registration */
    InvTrgVolReg = MatrixInverse(TrgVolReg,NULL);
    Src2TrgVolReg = MatrixMultiply(InvTrgVolReg,SrcVolReg,NULL);

    /* Loop through each source label and map its xyz to target */
    for(n = 0; n < srclabel->n_points; n++){

      /* load source label xyz into a vector */
      xyzSrc->rptr[0+1][0+1] = srclabel->lv[n].x;
      xyzSrc->rptr[0+2][0+1] = srclabel->lv[n].y;
      xyzSrc->rptr[0+3][0+1] = srclabel->lv[n].z;

      /* compute xyz location in target space */
      MatrixMultiply(Src2TrgVolReg,xyzSrc,xyzTrg);

      /* unload vector into target label */
      trglabel->lv[n].x = xyzTrg->rptr[0+1][0+1];
      trglabel->lv[n].y = xyzTrg->rptr[0+2][0+1];
      trglabel->lv[n].z = xyzTrg->rptr[0+3][0+1];
      trglabel->lv[n].stat = srclabel->lv[n].stat;

      /*printf("%3d  %6.4f %6.4f %6.4f    %6.4f %6.4f %6.4f\n",n,
       srclabel->lv[n].x,srclabel->lv[n].y,srclabel->lv[n].z,
       trglabel->lv[n].x,trglabel->lv[n].y,trglabel->lv[n].z);*/
    }

    MatrixFree(&SrcVolReg) ;
    MatrixFree(&TrgVolReg) ;
    MatrixFree(&InvTrgVolReg) ;
    MatrixFree(&Src2TrgVolReg) ;
    MatrixFree(&xyzSrc);
    MatrixFree(&xyzTrg);

  }/* done with volumetric mapping */

  /*--------------------- SURFACE-BASED MAPPING --------------------------*/
  if(!strcmp(regmethod,"surface")){

    printf("Starting surface-based mapping\n");

    /*** Load the source registration surface ***/
    if(strcmp(srcsubject,"ico")){
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,srcsubject,hemi,surfreg); 
      printf("Reading source registration \n %s\n",tmpstr);
      SrcSurfReg = MRISread(tmpstr);
      if(SrcSurfReg == NULL){
  fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
        exit(1);
      }
    }
    else{
      printf("Reading icosahedron, order = %d, radius = %g\n",
       srcicoorder,IcoRadius);
      SrcSurfReg = ReadIcoByOrder(srcicoorder,IcoRadius);
      if(SrcSurfReg==NULL) {
  printf("ERROR reading icosahedron\n");
  exit(1);
      }
    }

    /*** Load the target registration surface ***/
    if(strcmp(trgsubject,"ico")){
      /* load target xyz surface */
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,hemi,trgsurface); 
      printf("Reading target registration \n %s\n",tmpstr);
      TrgSurf = MRISread(tmpstr);
      if(TrgSurf == NULL){
  fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
        exit(1);
      }
      /* load target registration surface */
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,hemi,surfreg); 
      printf("Reading target registration \n %s\n",tmpstr);
      TrgSurfReg = MRISread(tmpstr);
      if(TrgSurfReg == NULL){
  fprintf(stderr,"ERROR: could not read %s\n",tmpstr);
        exit(1);
      }
    }
    else{
      printf("Reading icosahedron, order = %d, radius = %g\n",
       trgicoorder,IcoRadius);
      TrgSurfReg = ReadIcoByOrder(trgicoorder,IcoRadius);
      if(TrgSurfReg==NULL) {
  printf("ERROR reading icosahedron\n");
  exit(1);
      }
      TrgSurf = TrgSurfReg;
    }

    printf("Building target hash (res=16).\n");
    TrgHash = MHTfillVertexTableRes(TrgSurfReg, NULL,CURRENT_VERTICES,16);

    /* Loop through each source label and map its xyz to target */
    allzero = 1;
    for(n = 0; n < srclabel->n_points; n++){

      /* vertex number of the source label */
      srcvtxno = srclabel->lv[n].vno;
      if(srcvtxno < 0 || srcvtxno >= SrcSurfReg->nvertices){
  fprintf(stderr,"ERROR: label %d: vno = %d, max = %d\n",n,
    srcvtxno, SrcSurfReg->nvertices);
  exit(1);
      }

      if(srcvtxno != 0) allzero = 0;

      /* source vertex */
      srcvtx = &(SrcSurfReg->vertices[srcvtxno]);

      /* closest target vertex number */
      trgvtxno = MHTfindClosestVertexNo(TrgHash,TrgSurfReg,srcvtx,&dmin);

      /* target vertex */
      trgvtx = &(TrgSurf->vertices[trgvtxno]);

      trglabel->lv[n].x = trgvtx->x;
      trglabel->lv[n].y = trgvtx->y;
      trglabel->lv[n].z = trgvtx->z;
      trglabel->lv[n].stat = srclabel->lv[n].stat;
    }

    if(allzero){
      printf("---------------------------------------------\n");
      printf("WARNING: all source vertex numbers were zero.\n");
      printf("Make sure that the source label is surface-based.\n");      
      printf("---------------------------------------------\n");
    }

    MRISfree(&SrcSurfReg);
    MHTfree(&TrgHash);
    MRISfree(&TrgSurfReg);
    if(strcmp(trgsubject,"ico")) MRISfree(&TrgSurf);

  }/*---------- done with surface-based mapping -------------*/

  fprintf(stderr,"Writing label file %s \n",trglabelfile);
  LabelWrite(trglabel,trglabelfile);

  printf("mri_label2label: Done\n\n");

  return(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc;
  pargv = argv;
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;

    else if (!strcasecmp(option, "--version")) print_version() ;

    else if (!strcasecmp(option, "--debug"))   debug = 1;

    /* -------- source inputs ------ */
    else if (!strcmp(option, "--srcsubject")){
      if(nargc < 1) argnerr(option,1);
      srcsubject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srclabel")){
      if(nargc < 1) argnerr(option,1);
      srclabelfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--srcicoorder")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&srcicoorder);
      nargsused = 1;
    }
    /* -------- target inputs ------ */
    else if (!strcmp(option, "--trgsubject")){
      if(nargc < 1) argnerr(option,1);
      trgsubject = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trglabel")){
      if(nargc < 1) argnerr(option,1);
      trglabelfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trgsurface")){
      if(nargc < 1) argnerr(option,1);
      trgsurface = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trgicoorder")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&trgicoorder);
      nargsused = 1;
    }

    else if (!strcmp(option, "--hemi")){
      if(nargc < 1) argnerr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--regmethod")){
      if(nargc < 1) argnerr(option,1);
      regmethod = pargv[0];
      if(strcmp(regmethod,"surface") && strcmp(regmethod,"volume") &&
   strcmp(regmethod,"surf") && strcmp(regmethod,"vol")){
  fprintf(stderr,"ERROR: regmethod must be surface or volume\n");
  exit(1);
      }
      if(!strcmp(regmethod,"surf")) regmethod = "surface";
      if(!strcmp(regmethod,"vol"))  regmethod = "volume";
      nargsused = 1;
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(singledash(option))
  fprintf(stderr,"       Did you really mean -%s ?\n",option);
      exit(-1);
    }
    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  fprintf(stdout, "USAGE: %s \n",Progname) ;
  fprintf(stdout, "\n");
  fprintf(stdout, "   --srclabel   input label file \n");
  fprintf(stdout, "   --srcsubject source subject\n");
  fprintf(stdout, "   --trgsubject target subject\n");
  fprintf(stdout, "   --trglabel   output label file \n");
  fprintf(stdout, "   --regmethod  registration method (surface, volume) \n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --hemi        hemisphere (lh or rh) (with surface)\n");
  fprintf(stdout, "   --srcicoorder when srcsubject=ico\n");
  fprintf(stdout, "   --trgicoorder when trgsubject=ico\n");
  fprintf(stdout, "   --trgsurf     get xyz from this surface (white)\n");
  fprintf(stdout, "   --surfreg     surface registration (sphere.reg)  \n");
  fprintf(stdout, "\n");
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"srclabel = %s\n",  srclabelfile);
  fprintf(fp,"srcsubject = %s\n",srcsubject);
  fprintf(fp,"trgsubject = %s\n",trgsubject);
  fprintf(fp,"trglabel = %s\n",  trglabelfile);
  fprintf(fp,"regmethod = %s\n",regmethod);
  fprintf(fp,"\n");
  if(!strcmp(regmethod,"surface")){
    fprintf(fp,"hemi = %s\n",hemi);
    fprintf(fp,"trgsurface = %s\n",trgsurface);
    fprintf(fp,"surfreg = %s\n",surfreg);
  }
  fprintf(fp,"srcicoorder = %d\n",srcicoorder);
  fprintf(fp,"trgicoorder = %d\n",trgicoorder);
  fprintf(fp,"\n");

  return;
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;

  printf("
  Purpose: Converts a label in one subject's space to a label
  in another subject's space using either talairach or spherical
  as an intermediate registration space. 

  Example 1: If you have a label from subject fred called
    broca-fred.label defined on fred's left hemispherical 
    surface and you want to convert it to sally's surface, then

    mri_label2label --srclabel broca-fred.label  --srcsubject fred 
                    --trglabel broca-sally.label --trgsubject sally
                    --regmethod surface --hemi lh

    This will map from fred to sally using sphere.reg. The registration
    surface can be changed with --surfreg.

  Example 2: You could also do the same mapping using talairach 
    space as an intermediate:

    mri_label2label --srclabel broca-fred.label  --srcsubject fred 
                    --trglabel broca-sally.label --trgsubject sally
                    --regmethod volume

    Note that no hemisphere is specified with --regmethod volume.

  Notes:

  1. A label can be converted to/from talairach space by specifying
     the target/source subject as 'talairach'.
  2. A label can be converted to/from the icosahedron by specifying
     the target/source subject as 'ico'. When the source or target
     subject is specified as 'ico', then the order of the icosahedron
     must be specified with --srcicoorder/--trgicoorder.
  3. When the surface registration method is used, the xyz coordinates
     in the target label file are derived from the xyz coordinates
     from the target subject's white surface. This can be changed
     using the --trgsurf option.
  4. When the volume registration method is used, the xyz coordinates
     in the target label file are computed as xyzTrg = inv(Ttrg)*Tsrc*xyzSrc
     where Tsrc is the talairach transform in 
     srcsubject/mri/transforms/talairach.xfm, and where Ttrg is the talairach 
     transform in trgsubject/mri/transforms/talairach.xfm.
\n");



  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  fprintf(stderr, "%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if(n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(srcsubject == NULL){
    fprintf(stderr,"ERROR: no source subject specified\n");
    exit(1);
  }
  if(srclabelfile == NULL){
    fprintf(stderr,"ERROR: A source label path must be supplied\n");
    exit(1);
  }
  if(trgsubject == NULL){
    fprintf(stderr,"ERROR: no target subject specified\n");
    exit(1);
  }
  if(trglabelfile == NULL){
    fprintf(stderr,"ERROR: A target label path must be supplied\n");
    exit(1);
  }

  if(regmethod == NULL){
    fprintf(stderr,"ERROR: Must specify a registration method\n");
    exit(1);
  }

  if(!strcmp(regmethod,"surface")){
    if(hemi == NULL){
      fprintf(stderr,"ERROR: no hemisphere specified\n");
      exit(1);
    }
  }
  else{ /* volume */
    if(!strcmp(srcsubject,"ico") || !strcmp(trgsubject,"ico")){
      fprintf(stderr,"ERROR: cannot use volume registration "
        "method with subject ico\n");      
      exit(1);
    }
    if(hemi != NULL){
      fprintf(stderr,"ERROR: cannot specify hemisphere with vol reg method\n");
      exit(1);
    }
  }

  if(!strcmp(srcsubject,"ico") && srcicoorder < 0){
    fprintf(stderr,"ERROR: must specify src ico order with srcsubject=ico\n");
    exit(1);
  }

  if(!strcmp(trgsubject,"ico") && trgicoorder < 0){
    fprintf(stderr,"ERROR: must specify trg ico order with trgsubject=ico\n");
    exit(1);
  }

  if(!strcmp(regmethod,"surface") && (!strcmp(srcsubject,"talairach") ||
              !strcmp(trgsubject,"talairach"))){
    fprintf(stderr,"ERROR: cannot use talairach with surface mapping\n");
    exit(1);
  }
    


  return;
}

/*---------------------------------------------------------------*/
static int singledash(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] != '-') return(1);
  return(0);
}
