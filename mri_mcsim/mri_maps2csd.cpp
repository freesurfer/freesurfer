/*
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

/*
BEGINHELP --------------------------------------------------------------

ENDHELP --------------------------------------------------------------
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
double round(double x);
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <unistd.h>

#include "macros.h"
#include "utils.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "label.h"
#include "matrix.h"
#include "annotation.h"
#include "fmriutils.h"
#include "cmdargs.h"
#include "fsglm.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"
#include "randomfields.h"

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int SaveOutput(void);
int main(int argc, char *argv[]) ;

const char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

const char *subject = NULL;
const char *surfname = "white";
const char *hemi = NULL;

char *surfpath=NULL, *csdfile=NULL;
MRIS *surf;
const char *signstr=NULL;
double thresh = -1;
char tmpstr[2000], *SUBJECTS_DIR;
char *inputs[10000];
int ninputs = 0;
CSD *csd = NULL;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  int nargs;

  csd = CSDalloc();
  csd->threshsign = -2;
  csd->thresh = -1;

  nargs = handleVersionOption(argc, argv, "mri_maps2csd");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if (argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if (checkoptsonly) return(0);
  dump_options(stdout);

  // Load the target surface
  printf("Loading %s\n",surfpath);
  surf = MRISread(surfpath);
  if(!surf) return(1);

  MRI *mri;
  int nClusters,vtx,n;
  SURFCLUSTERSUM *SurfClustList;
  printf("ninputs = %d\n",ninputs);

  strcpy(csd->simtype,"perm");
  strcpy(csd->anattype,"surface");
  strcpy(csd->subject,subject);
  strcpy(csd->hemi,hemi);
  strcpy(csd->contrast,"no-con");
  csd->nreps = ninputs;
  CSDallocData(csd);

  double threshadj = 0;
  if(csd->threshsign == 0) threshadj = csd->thresh;
  else                     threshadj = csd->thresh - log10(2.0); // one-sided test
  printf("threshadj = %g\n",threshadj);

  int m = 0;
  for(n=0; n < ninputs; n++){
    mri = MRIread(inputs[n]);
    if(mri==NULL) exit(1);
    for (vtx = 0; vtx < surf->nvertices; vtx++)
      surf->vertices[vtx].val = MRIgetVoxVal(mri,vtx,0,0,0);
    SurfClustList = sclustMapSurfClusters(surf,threshadj,-1,csd->threshsign, 0,&nClusters,NULL,NULL);
    double csize = sclustMaxClusterArea(SurfClustList, nClusters);
    int cmax,rmax,smax;
    double sigmax = MRIframeMax(mri,0,NULL,csd->threshsign,&cmax,&rmax,&smax);
    printf("%3d %4d  %g %g\n",n,nClusters,csize,sigmax);
    free(SurfClustList);
    csd->nClusters[m] = nClusters;
    csd->MaxClusterSize[m] = csize;
    csd->MaxSig[m] = sigmax;
    m++;
  }
  FILE *fp = fopen(csdfile,"w");
  fprintf(fp,"# ClusterSimulationData 2\n");
  fprintf(fp,"# mri_maps2csd simulation sim\n");
  CSDprint(fp, csd);
  fclose(fp);

  printf("mri_maps2csd done\n");
  exit(0);
}
/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--csd")) {
      if(nargc < 1) CMDargNErr(option,1);
      csdfile = pargv[0];
      nargsused = 1;
    } 
    else if (!strcasecmp(option, "--s")){
      if(nargc < 3) CMDargNErr(option,3);
      subject = pargv[0];
      hemi = pargv[1];
      surfname = pargv[2];
      nargsused = 3;
    } 
    else if (!strcasecmp(option, "--surf")) {
      if(nargc < 1) CMDargNErr(option,1);
      surfpath = pargv[0];
      nargsused = 1;
    } 
    else if (!strcmp(option, "--sd")) {
      if(nargc < 1) CMDargNErr(option,1);
      setenv("SUBJECTS_DIR",pargv[0],1);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--thresh")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&csd->thresh);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--sign")) {
      if (nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&csd->threshsign);
      nargsused = 1;
    } 
    else if(!strcasecmp(option, "--i")) {
      if (nargc < 1) CMDargNErr(option,1);
      inputs[ninputs] = pargv[0];
      ninputs++;
      nargsused = 1;
    } 
    else {
      inputs[ninputs] = option;
      ninputs++;
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
  printf("%s   inputs\n",Progname) ;
  printf("\n");
  printf("   --csd csdfile\n");
  printf("   --s subjectname hemi surf\n");
  //printf("   --surf surfpath\n");
  printf("   --thresh thresh (-log10 cluster-forming thresh)\n");
  printf("   --sign +1,-1,0 : a sign causes the thresh to be adjusted\n");
  printf("   --i input (or just put them on the command line)\n");
  printf("   \n");
  printf("   --sd SUBJECTS_DIR\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  std::cout << getVersion() << std::endl;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void) {
  std::cout << getVersion() << std::endl;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void) {

  if(surfpath == NULL && subject == NULL) {
    printf("ERROR: must specify a surface\n");
    exit(1);
  }
  if(surfpath != NULL && subject != NULL) {
    printf("ERROR: cannot spec both --s and --surf\n");
    exit(1);
  }
  if(surfpath == NULL){
    sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
    surfpath = strcpyalloc(tmpstr);
  }
  else {
    subject = "unknown";
    hemi = "unknown";
  }
  if(csdfile==NULL) {
    printf("ERROR: need to specify a csd file\n");
    exit(1);
  }
  if(ninputs == 0){
    printf("ERROR: No inputs specified\n");
    exit(1);
  }
  if(csd->thresh < 0){
    printf("ERROR: must spec --thresh\n");
    exit(1);
  }
  if(csd->threshsign == -2){
    printf("ERROR: must spec --sign +1,-1,0\n");
    exit(1);
  }

  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp) {
  fprintf(fp,"\n");
  fprintf(fp,"%s\n", getVersion().c_str());
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());
  fflush(fp);
  return;
}

