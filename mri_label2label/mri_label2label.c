/*----------------------------------------------------------
  Name: mri_label2label.c
  $Id: mri_label2label.c,v 1.11 2002/12/04 23:12:21 greve Exp $
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

#include "MRIio_old.h"
#include "error.h"
#include "diag.h"
#include "mrisurf.h"
#include "mri.h"
#include "label.h"
#include "registerio.h"
#include "mri.h"
#include "mri2.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void argnerr(char *option, int n);
static void dump_options(FILE *fp);
static int  singledash(char *flag);
static int  isflag(char *flag);
static int  nth_is_arg(int nargc, char **argv, int nth);

int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_label2label.c,v 1.11 2002/12/04 23:12:21 greve Exp $";
char *Progname = NULL;

char  *srclabelfile = NULL;
LABEL *srclabel     = NULL;
LABEL *tmplabel     = NULL;
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
float hashres = 16;
int usehash = 0;

int debug = 0;

char *SUBJECTS_DIR = NULL;
char *MRI_DIR = NULL;
FILE *fp;

char tmpstr[2000];

char *srcmaskfile, *srcmaskfmt, *srcmasksign = "abs";
int srcmaskframe = 0;
float srcmaskthresh = 0.0;
MRI *SrcMask;

int useprojabs = 0, useprojfrac = 0; 
float projabs = 0.0, projfrac = 0.0;

/*-------------------------------------------------*/
int main(int argc, char **argv)
{
  MATRIX *xyzSrc, *xyzTrg;
  MHT *TrgHash;
  VERTEX *srcvtx, *trgvtx;
  int n,err,srcvtxno,trgvtxno,allzero;
  float dmin, projdist=0.0, dx, dy, dz;
  float SubjRadius, Scale;
  char fname[2000];
  

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
      printf("Rescaling ... ");
      SubjRadius = MRISaverageRadius(SrcSurfReg) ;
      Scale = IcoRadius / SubjRadius;
      MRISscaleBrain(SrcSurfReg, SrcSurfReg, Scale);
      printf(" original radius = %g\n",SubjRadius);
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
      
    /*** Load the target surfaces ***/
    if(strcmp(trgsubject,"ico")){
      /* load target xyz surface */
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,trgsubject,
        hemi,trgsurface); 
      printf("Reading target surface \n %s\n",tmpstr);
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
      printf("Rescaling ... ");
      SubjRadius = MRISaverageRadius(TrgSurfReg) ;
      Scale = IcoRadius / SubjRadius;
      MRISscaleBrain(TrgSurfReg, TrgSurfReg, Scale);
      printf(" original radius = %g\n",SubjRadius);
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
    
    if(usehash){
      printf("Building target registration hash (res=%g).\n",hashres);
      TrgHash = MHTfillVertexTableRes(TrgSurfReg, NULL,
              CURRENT_VERTICES,hashres);
    }
    if(useprojfrac){
      sprintf(fname,"%s/%s/surf/%s.thickness",SUBJECTS_DIR,srcsubject,hemi);
      printf("Reading thickness %s\n",fname);
      MRISreadCurvatureFile(TrgSurf, fname);
      printf("Done\n");
    }

    /* handle source mask */
    if(srcmaskfile != NULL){
      printf("INFO: masking label\n");
      //SrcMask = MRIloadSurfVals(srcmaskfile, srcmaskfmt, NULL,
      //      srcsubject, hemi, NULL);
      
      SrcMask = MRISloadSurfVals(srcmaskfile, srcmaskfmt, SrcSurfReg,
         NULL,NULL,NULL);
      if(SrcMask == NULL) exit(1);
      tmplabel = MaskSurfLabel(srclabel, SrcMask, 
             srcmaskthresh, srcmasksign, srcmaskframe);
      if(tmplabel == NULL) exit(1);
      LabelFree(&srclabel) ;
      srclabel = tmplabel;
      printf("Found %d points in source label after masking.\n",
       srclabel->n_points);
      if(srclabel->n_points == 0){
  printf("ERROR: no overlap between mask and label\n");
  exit(1);
      }
      trglabel->n_points = srclabel->n_points;
    }
    
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
      if(usehash){
  trgvtxno = MHTfindClosestVertexNo(TrgHash,TrgSurfReg,srcvtx,&dmin);
  if(trgvtxno < 0){
    printf("ERROR: trgvtxno = %d < 0\n",trgvtxno);
    printf("srcvtxno = %d, dmin = %g\n",srcvtxno,dmin);
    printf("srcxyz = %g, %g, %g\n",srcvtx->x,srcvtx->y,srcvtx->z);
    exit(1);
  }
      }
      else{
  trgvtxno = MRISfindClosestVertex(TrgSurfReg,srcvtx->x,srcvtx->y,
           srcvtx->z);
      }
      /* target vertex */
      trgvtx = &(TrgSurf->vertices[trgvtxno]);

      if(useprojabs || useprojfrac){
  if(useprojabs)  projdist = projabs;
  if(useprojfrac) projdist = projfrac * trgvtx->curv;
  dx = projdist*trgvtx->nx;
  dy = projdist*trgvtx->ny;
  dz = projdist*trgvtx->nz;
      }
      else {
  dx = 0.0;
  dy = 0.0;
  dz = 0.0;
      }

      trglabel->lv[n].vno = trgvtxno;
      trglabel->lv[n].x = trgvtx->x + dx;
      trglabel->lv[n].y = trgvtx->y + dy;
      trglabel->lv[n].z = trgvtx->z + dz;
      trglabel->lv[n].stat = srclabel->lv[n].stat;
    }

    if(allzero){
      printf("---------------------------------------------\n");
      printf("WARNING: all source vertex numbers were zero.\n");
      printf("Make sure that the source label is surface-based.\n");      
      printf("---------------------------------------------\n");
    }

    MRISfree(&SrcSurfReg);
    MRISfree(&TrgSurfReg);
    if(usehash) MHTfree(&TrgHash);
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
    else if (!strcasecmp(option, "--hash"))   usehash = 1;
    else if (!strcasecmp(option, "--nohash")) usehash = 0;

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
    else if (!strcmp(option, "--srcmask")){
      if(nargc < 2) argnerr(option,2);
      srcmaskfile = pargv[0];
      sscanf(pargv[1],"%f",&srcmaskthresh); 
      nargsused = 2;
      if(nth_is_arg(nargc, pargv, 2)){
  srcmaskfmt = pargv[2]; nargsused ++;
      }
    }
    else if (!strcmp(option, "--srcmaskframe")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&srcmaskframe);  nargsused = 1;
    }
    else if (!strcmp(option, "--srcmasksign")){
      if(nargc < 1) argnerr(option,1);
      srcmasksign = pargv[0]; nargsused = 1;
      if(strcmp(srcmasksign,"abs") &&
   strcmp(srcmasksign,"pos") &&
   strcmp(srcmasksign,"neg")){
  printf("ERROR: srcmasksign = %s, must be either "
         "abs, pos, or neg\n", srcmasksign);
  exit(1);
      }
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
    else if (!strcmp(option, "--surfreg")){
      if(nargc < 1) argnerr(option,1);
      surfreg = pargv[0];
      nargsused = 1;
    }
    else if (!strcmp(option, "--trgicoorder")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%d",&trgicoorder);
      nargsused = 1;
    }
    else if (!strcmp(option, "--hashres")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&hashres);
      nargsused = 1;
    }

    else if (!strcmp(option, "--projabs")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&projabs);
      useprojabs = 1;
      nargsused = 1;
    }

    else if (!strcmp(option, "--projfrac")){
      if(nargc < 1) argnerr(option,1);
      sscanf(pargv[0],"%f",&projfrac);
      useprojfrac = 1;
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
  fprintf(stdout, "   --srcmask     surfvalfile thresh <format>\n");
  fprintf(stdout, "   --srcmasksign sign (<abs>,pos,neg)\n");
  fprintf(stdout, "   --srcmaskframe 0-based frame number <0>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "   --projabs  dist project dist mm along surf normal\n");
  fprintf(stdout, "   --projfrac frac project frac of thickness along surf normal\n");
  fprintf(stdout, "\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;

  printf(
"  Purpose: Converts a label in one subject's space to a label\n"
"  in another subject's space using either talairach or spherical\n"
"  as an intermediate registration space. \n"
"\n"
"  If a source mask is used, then the input label must have been\n"
"  created from a surface (ie, the vertex numbers are valid). The \n"
"  format can be anything supported by mri_convert or curv or paint.\n"
"  Vertices in the source label that do not meet threshold in the\n"
"  mask will be removed from the label. See Example 2.\n"
"\n"
"  Example 1: If you have a label from subject fred called\n"
"    broca-fred.label defined on fred's left hemispherical \n"
"    surface and you want to convert it to sally's surface, then\n"
"\n"
"    mri_label2label --srclabel broca-fred.label  --srcsubject fred \n"
"                    --trglabel broca-sally.label --trgsubject sally\n"
"                    --regmethod surface --hemi lh\n"
"\n"
"    This will map from fred to sally using sphere.reg. The registration\n"
"    surface can be changed with --surfreg.\n"
"\n"
"  Example 2: Same as Example 1 but with a mask\n"
"\n"
"    mri_label2label --srclabel broca-fred.label  --srcsubject fred \n"
"                    --trglabel broca-sally.label --trgsubject sally\n"
"                    --regmethod surface --hemi lh\n"
"                    --srcmask  fred-omnibus-sig 2 bfloat\n"
"\n"
"    This will load the bfloat data from fred-omnibus-sig and create\n"
"    a mask by thresholding the first frame absolute values at 2.\n"
"    To change it to only the positive values of the 3rd frame, add\n"
"         --srcmasksign pos --srcmaskframe 2   \n"
"\n"
"\n"
"  Example 3: You could also do the same mapping using talairach \n"
"    space as an intermediate:\n"
"\n"
"    mri_label2label --srclabel broca-fred.label  --srcsubject fred \n"
"                    --trglabel broca-sally.label --trgsubject sally\n"
"                    --regmethod volume\n"
"\n"
"    Note that no hemisphere is specified with --regmethod volume.\n"
"\n"
"  Notes:\n"
"\n"
"  1. A label can be converted to/from talairach space by specifying\n"
"     the target/source subject as 'talairach'.\n"
"  2. A label can be converted to/from the icosahedron by specifying\n"
"     the target/source subject as 'ico'. When the source or target\n"
"     subject is specified as 'ico', then the order of the icosahedron\n"
"     must be specified with --srcicoorder/--trgicoorder.\n"
"  3. When the surface registration method is used, the xyz coordinates\n"
"     in the target label file are derived from the xyz coordinates\n"
"     from the target subject's white surface. This can be changed\n"
"     using the --trgsurf option.\n"
"  4. When the volume registration method is used, the xyz coordinates\n"
"     in the target label file are computed as xyzTrg = inv(Ttrg)*Tsrc*xyzSrc\n"
"     where Tsrc is the talairach transform in \n"
"     srcsubject/mri/transforms/talairach.xfm, and where Ttrg is the talairach \n"
"     transform in trgsubject/mri/transforms/talairach.xfm.\n"
"  5. The registration surfaces are rescaled to a radius of 100 (including \n"
"     the ico)\n"
"  6. Projections along the surface normal can be either negative or\n"
"     positive, but can only be used with surface registration method.\n"
);



  exit(1) ;
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
  if(!strcmp(srcsubject,"ico")) fprintf(fp,"srcicoorder = %d\n",srcicoorder);
  if(!strcmp(trgsubject,"ico")) fprintf(fp,"trgicoorder = %d\n",trgicoorder);
  fprintf(fp,"usehash = %d\n",usehash);

  if(srcmaskfile != NULL){
    fprintf(fp,"srcmask %s, %s \n",srcmaskfile, srcmaskfmt);
    fprintf(fp,"srcmaskthresh %g %s\n",srcmaskthresh, srcmasksign);
    fprintf(fp,"srcmaskframe %d\n",srcmaskframe);
  }
  printf("Use ProjAbs  = %d, %g\n",useprojabs,projabs);
  printf("Use ProjFrac = %d, %g\n",useprojfrac,projfrac);

  fprintf(fp,"\n");

  return;
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
/*---------------------------------------------------------------*/
static int isflag(char *flag)
{
  int len;
  len = strlen(flag);
  if(len < 2) return(0);

  if(flag[0] == '-' && flag[1] == '-') return(1);
  return(0);
}
/*---------------------------------------------------------------*/
static int nth_is_arg(int nargc, char **argv, int nth)
{
  /* Checks that nth arg exists and is not a flag */
  /* nth is 0-based */

  /* check that there are enough args for nth to exist */
  if(nargc <= nth) return(0); 

  /* check whether the nth arg is a flag */
  if(isflag(argv[nth])) return(0);

  return(1);
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

  if(useprojabs && useprojfrac){
    fprintf(stderr,"ERROR: cannot use absolute and fractional projection\n");
    exit(1);
  }

  if( (useprojabs || useprojfrac) &&  strcmp(regmethod,"surface") ){
    fprintf(stderr,"ERROR: must use surface regmethod with absolute "
      "or fractional projection\n");
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


