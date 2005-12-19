/*
BEGINHELP

ENDHELP
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
#include "fmriutils.h"
#include "mrisurf.h"
#include "mri2.h"
#include "fio.h"
#include "version.h"
#include "annotation.h"
#include "cmdargs.h"
#include "timer.h"
#include "matfile.h"
#include "randomfields.h"
double MRISmeanInterVertexDist(MRIS *surf);


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *subject=NULL, *hemi=NULL, *SUBJECTS_DIR=NULL;
char *surfname="white";
char *surfpath=NULL;
char tmpstr[2000];

MRIS *surf;
int dof = 100;
int nitersmax = 100;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int nargs, nthiter=0;
  nthiter=0;
  MRI *mri=NULL, *var=NULL, *mri0, *delta, *deltasm=NULL, *xyz;
  double gmax, vrfmn, vrfstd, gstd, fwhm;

  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;
  cmdline = argv2cmdline(argc,argv);
  uname(&uts);
  getcwd(cwd,2000);

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;
  if(argc == 0) usage_exit();
  parse_commandline(argc, argv);
  check_options();
  if(checkoptsonly) return(0);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR == NULL){
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
  surfpath = strcpyalloc(tmpstr);

  if(debug) dump_options(stdout);

  surf = MRISread(surfpath);
  if(surf == NULL){
    printf("ERROR: could not read %s\n",surfpath);
    exit(1);
  }
  printf("Number of vertices %d\n",surf->nvertices);
  printf("Number of faces    %d\n",surf->nfaces);
  printf("Avg IterVertex     %lf\n",MRISmeanInterVertexDist(surf));

  xyz = MRIallocSequence(surf->nvertices,1,1,MRI_FLOAT,4);
  MRIcopyMRIS(xyz,surf,0,"x");
  MRIcopyMRIS(xyz,surf,1,"y");
  MRIcopyMRIS(xyz,surf,2,"z");
  MRIcopyMRIS(xyz,surf,3,"area");
  MRIwrite(xyz,"xyz.mgh");

  delta = MRIalloc(surf->nvertices,1,1,MRI_FLOAT);
  MRIsetVoxVal(delta,(int)(surf->nvertices/2),0,0,0,1);
  MRIwrite(delta,"delta.mgh");

  deltasm = MRISgaussianSmooth(surf, delta, 2, NULL, 5.0);
  //deltasm = MRISsmoothMRI(surf, delta, 2, deltasm);
  MRIwrite(deltasm,"deltasm.mgh");

  mri0 = MRIrandn(surf->nvertices,1,1,dof,0, 1, NULL);
  mri = MRIcopy(mri0,NULL);
  
  for(nthiter = 1; nthiter <= nitersmax; nthiter++){
    //MRISsmoothMRI(surf, mri, 1, mri);
    MRISgaussianSmooth(surf, mri0, nthiter, mri, 5.0);

    var = fMRIvariance(mri, dof, 0, var);
    RFglobalStats(var, NULL, &vrfmn, &vrfstd, &gmax);
    gstd = 1/(2*sqrt(vrfmn*PI));
    fwhm = gstd*sqrt(log(256.0));
    printf("%3d %lf  %lf  %lf %lf\n",nthiter,vrfmn,vrfstd,gstd,fwhm);
    
  }

  return 0;
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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--s")){
      if(nargc < 1) CMDargNErr(option,1);
      subject = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--h")){
      if(nargc < 1) CMDargNErr(option,1);
      hemi = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--surf")){
      if(nargc < 1) CMDargNErr(option,1);
      surfname = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--dof")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&dof);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--niters")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%d",&nitersmax);
      nargsused = 1;
    }
    else{
      fprintf(stderr,"ERROR: Option %s unknown\n",option);
      if(CMDsingleDash(option))
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
  printf("USAGE: %s --s --h --surf --dof --niters \n",Progname) ;
  printf("\n");
  printf("   --s subject \n");
  printf("   --h hemi \n");
  printf("   --surf surf\n");
  printf("   --dof  dof\n");
  printf("   --niters nitersmax\n");
  printf("\n");
  printf("   --debug     turn on debugging\n");
  printf("   --checkopts don't run anything, just check options and exit\n");
  printf("   --help      print out information on how to use this program\n");
  printf("   --version   print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("WARNING: this program is not yet tested!\n");
  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"%s\n",Progname);
  fprintf(fp,"FREESURFER_HOME %s\n",getenv("FREESURFER_HOME"));
  fprintf(fp,"SUBJECTS_DIR %s\n",getenv("SUBJECTS_DIR"));
  fprintf(fp,"cwd       %s\n",cwd);
  fprintf(fp,"cmdline   %s\n",cmdline);
  fprintf(fp,"timestamp %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname   %s\n",uts.sysname);
  fprintf(fp,"hostname  %s\n",uts.nodename);
  fprintf(fp,"machine   %s\n",uts.machine);
  fprintf(fp,"user      %s\n",VERuser());
  if(subject) fprintf(fp,"subject %s\n",subject);
  if(hemi)     fprintf(fp,"hemi     %s\n",hemi);
  if(surfname) fprintf(fp,"surfname %s\n",surfname);
  fprintf(fp,"dof %d\n",dof);
  fprintf(fp,"nitersmax %d\n",nitersmax);

  return;
}
/*---------------------------------------------------------------------*/
double MRISmeanInterVertexDist(MRIS *surf)
{
  int vtx, nbrvtx, nnbrs, nthnbr;
  double dx, dy, dz, x0, y0, z0, xn, yn, zn, d;
  double dnbrsum, dnbrmn, dsum;

  dsum = 0.0;
  for(vtx = 0; vtx < surf->nvertices; vtx++){
    nnbrs = surf->vertices[vtx].vnum;
    x0 = surf->vertices[vtx].x;
    y0 = surf->vertices[vtx].y;
    z0 = surf->vertices[vtx].z;
    dnbrsum=0.0;
    for(nthnbr = 0; nthnbr < nnbrs; nthnbr++){
      nbrvtx = surf->vertices[vtx].v[nthnbr];
      xn = surf->vertices[nbrvtx].x;
      yn = surf->vertices[nbrvtx].y;
      zn = surf->vertices[nbrvtx].z;
      dx = x0-xn;
      dy = y0-yn;
      dz = z0-zn;
      d = sqrt(dx*dx + dy*dy + dz*dz);
      dnbrsum += d;
    }/* end loop over neighbor */
    dnbrmn = dnbrsum/nnbrs;
    dsum += dnbrmn;
  } /* end loop over vertex */

  return(dsum/surf->nvertices);
}
