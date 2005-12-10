

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
#include "gsl/gsl_cdf.h"
#include "pdf.h"
#include "fsgdf.h"
#include "timer.h"
#include "matfile.h"
#include "volcluster.h"
#include "surfcluster.h"


static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_diff.c,v 1.3 2005/12/10 05:49:08 greve Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *InVol1File=NULL;
char *InVol2File=NULL;
char *subject, *hemi, *SUBJECTS_DIR;
double pixthresh=0, resthresh=0;
char *DiffFile=NULL;

MRI *InVol1=NULL, *InVol2=NULL;
int CheckResolution=0;

/*---------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int nargs, r, c, s, f;
  int rmax, cmax, smax, fmax;
  double diff,maxdiff;
  double val1, val2;
  FILE *fp=NULL;

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
  dump_options(stdout);

  SUBJECTS_DIR = getenv("SUBJECTS_DIR");
  if(SUBJECTS_DIR == NULL){
    printf("ERROR: SUBJECTS_DIR not defined in environment\n");
    exit(1);
  }
  if(DiffFile){
    if(fio_FileExistsReadable(DiffFile)) unlink(DiffFile);
    if(fio_FileExistsReadable(DiffFile)){
      printf("ERROR: could not delete %s\n",DiffFile);
      exit(1);
    }
  }

  InVol1 = MRIread(InVol1File);
  if(InVol1==NULL) exit(1);

  InVol2 = MRIread(InVol2File);
  if(InVol2==NULL) exit(1);
  
  if(InVol1->width   != InVol2->width  ||
     InVol1->height  != InVol2->height ||
     InVol1->depth   != InVol2->depth  ||
     InVol1->nframes != InVol2->nframes){
    printf("ERROR: volumes differ in dimension\n");
    exit(1);
  }
  
  // Compare pixel values
  c=r=s=f=0;
  val1 = MRIgetVoxVal(InVol1,c,r,s,f);
  val2 = MRIgetVoxVal(InVol2,c,r,s,f);
  maxdiff = fabs(val1-val2);
  cmax=rmax=smax=fmax=0;
  for(c=0; c < InVol1->width; c++){
    for(r=0; r < InVol1->height; r++){
      for(s=0; s < InVol1->depth; s++){
	for(f=0; f < InVol1->nframes; f++){
	  val1 = MRIgetVoxVal(InVol1,c,r,s,f);
	  val2 = MRIgetVoxVal(InVol2,c,r,s,f);
	  diff = fabs(val1-val2);
	  if(maxdiff < diff){
	    maxdiff = diff;
	    cmax = c;
	    rmax = r;
	    smax = s;
	    fmax = f;
	  }
	}
      }
    }
  }

  printf("maxdiff %g at %d %d %d %d\n",maxdiff,cmax,rmax,smax,fmax);
  if(maxdiff > pixthresh){
    printf("Volumes differ\n");
    if(DiffFile){
      fp = fopen(DiffFile,"w");
      if(fp==NULL){
	printf("ERROR: could not open %s\n",DiffFile);
	exit(1);
      }
      dump_options(fp);
      fprintf(fp,"maxdiff %g at %d %d %d %d\n",maxdiff,cmax,rmax,smax,fmax);
      fprintf(fp,"Volumes differ in pixel value\n");
      fclose(fp);
    }
    exit(101);
  }

  if(CheckResolution){
    if(InVol1->xsize   != InVol2->xsize  ||
       InVol1->ysize   != InVol2->ysize ||
       InVol1->zsize   != InVol2->zsize){
      printf("Volumes differ in resolution\n");
      if(DiffFile){
	fp = fopen(DiffFile,"w");
	if(fp==NULL){
	  printf("ERROR: could not open %s\n",DiffFile);
	  exit(1);
	}
	dump_options(fp);
	fprintf(fp,"Volumes differ in resolution\n");
	fclose(fp);
      }
      exit(102);
    }
  }


  exit(0);
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
    else if (!strcasecmp(option, "--checkopts"))   checkoptsonly = 1;
    else if (!strcasecmp(option, "--nocheckopts")) checkoptsonly = 0;

    else if (!strcasecmp(option, "--v1")){
      if(nargc < 1) CMDargNErr(option,1);
      InVol1File = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--v2")){
      if(nargc < 1) CMDargNErr(option,1);
      InVol2File = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--thresh")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&pixthresh);
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--res-thresh")){
      if(nargc < 1) CMDargNErr(option,1);
      sscanf(pargv[0],"%lf",&resthresh);
      CheckResolution=1;
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--diff-file")){
      if(nargc < 1) CMDargNErr(option,1);
      DiffFile = pargv[0];
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
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
/* --------------------------------------------- */
static void check_options(void)
{
  if(InVol1File==NULL){
    printf("ERROR: need to spec volume 1\n");
    exit(1);
  }
  if(InVol2File==NULL){
    printf("ERROR: need to spec volume 2\n");
    exit(1);
  }
  return;
}
/* ------------------------------------------------------ */
static void usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}
/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"%s\n",Progname);
  fprintf(fp,"FREESURFER_HOME %s\n",getenv("FREESURFER_HOME"));
  fprintf(fp,"cwd       %s\n",cwd);
  fprintf(fp,"cmdline   %s\n",cmdline);
  fprintf(fp,"timestamp %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname   %s\n",uts.sysname);
  fprintf(fp,"hostname  %s\n",uts.nodename);
  fprintf(fp,"machine   %s\n",uts.machine);
  fprintf(fp,"user      %s\n",VERuser());
  fprintf(fp,"v1        %s\n",InVol1File);
  fprintf(fp,"v2        %s\n",InVol2File);
  fprintf(fp,"pixthresh %lf\n",pixthresh);
  //fprintf(fp,"checkres  %d\n",CheckResolution);
  //fprintf(fp,"resthresh %lf\n",pixthresh);
  fprintf(fp,"diff-file %s\n",DiffFile);
  return;
}
/* --------------------------------------------- */
static void print_usage(void)
{
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --v1 volfile1 : first  input volume \n");
  printf("   --v2 volfile2 : second input volume \n");
  printf("   --thresh thresh : must be greater than this \n");
  printf("   --diff-file DiffFile : store diff info in this file. \n");
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
  printf("Determines whether two volumes differ in pixel data.\n");
  printf("The voxels are compared on a voxel-by-voxel basis.\n");
  printf("If the absolute difference is greater than a threshold,\n");
  printf("the volumes are considered to be different and the program\n");
  printf("exits with value 101. If a DiffFile has been specified,\n");
  printf("that file is created and information about the environment\n");
  printf("and how the volumes were different is stored in it.\n");
  printf("The DiffFile is not created if the volumes are not different.\n");
  printf("If DiffFile exists when the program is run, it is immediately\n");
  printf("deleted.\n");
  printf("\n");
  printf("By default, the threshold is 0, but can be changed with\n");
  printf("--thresh.\n");
  printf("\n");
  printf("If an error is encountered, then exits with status=1.\n");
  printf("An error includes differing in dimension\n");
  printf("\n");
  printf("If there was no error and no difference, exit status=0\n");
  printf("\n");
  printf("EXAMPLES:\n");
  printf("\n");
  printf("mri_diff --v1 vol1.mgh --v2 vol2.img\n");
  printf("mri_diff --v1 vol1.nii --v2 vol2.img --thresh 20\n");
  printf("mri_diff --v1 vol1.mgh --v2 vol2.bhdr --thresh 20\n");
  printf("   --diff DiffFile\n");
  printf("\n");
  printf("EXAMPLE DiffFile:\n");
  printf("\n");
  printf("$Id: mri_diff.c,v 1.3 2005/12/10 05:49:08 greve Exp $\n");
  printf("mri_diff\n"
         "FREESURFER_HOME /space/cadet/2/users/cadet/freesurfer\n"
         "cwd       /space/cadet/2/users/cadet/dev/mri_diff\n"
         "cmdline   ./mri_diff --v1 /home/cadet/projects/swf/skb/skb01.v1/bold/010/f.bhdr --v2 /home/cadet/projects/swf/skb/skb01.v1/bold/010/fmc.bhdr --diff-file mydifffile.dat --thresh 10 \n"
         "timestamp 2005/12/10-04:34:08-GMT\n"
         "sysname   Linux\n"
         "hostname  sloth\n"
         "machine   x86_64\n"
         "user      cadet\n"
         "v1        /space/cadet/projects/swf/skb/skb01.v1/bold/010/f.bhdr\n"
         "v2        /space/cadet/projects/swf/skb/skb01.v1/bold/010/fmc.bhdr\n"
         "pixthresh 10.000000\n"
         "diff-file mydifffile.dat\n"
         "maxdiff 254 at 42 49 26 95\n"
	   "Volumes differ in pixel value\n");
  printf("\n");
  exit(1) ;
}
