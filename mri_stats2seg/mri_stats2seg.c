
/*--------------------------------------------------------------
Example Usage:

set subject = fsr-tst

mri_segstats \
  --in  $SUBJECTS_DIR/$subject/mri/norm.mgz \
  --seg $SUBJECTS_DIR/$subject/mri/aseg.mgz \
  --ctab-default \
  --avgwfvol stats.mgh --avgwf stats.txt \
  --sum sum.txt

./mri_stats2seg --stat stats.mgh \
  --seg $SUBJECTS_DIR/$subject/mri/aseg.mgz \
  --o asegstats.mgh

tkmedit $subject norm.mgz -aux ./asegstats.mgh\
  -segmentation $SUBJECTS_DIR/$subject/mri/aseg.mgz \
      $FREESURFER_HOME/FreeSurferColorLUT.txt
--------------------------------------------------------------*/

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

static int  parse_commandline(int argc, char **argv);
static void check_options(void);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;
static void dump_options(FILE *fp);
int main(int argc, char *argv[]) ;

static char vcid[] = "$Id: mri_stats2seg.c,v 1.4 2006/05/11 21:57:35 nicks Exp $";
char *Progname = NULL;
char *cmdline, cwd[2000];
int debug=0;
int checkoptsonly=0;
struct utsname uts;

char *TempVolFile=NULL;
char *subject, *hemi, *SUBJECTS_DIR;

char *statfile=NULL;
MRI *statmri;
char *segfile=NULL;
MRI *seg;

char *outfile=NULL;
MRI *out;


/*---------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  int nargs,r,c,s,f,segid;
  MRIS *mris;
  MRI *mri;
  double val;

  if(0){
  mris = MRISread("/space/greve/1/users/greve/subjects/fsr-tst/surf/lh.white");
  MRISreadAnnotation(mris,"/space/greve/1/users/greve/subjects/fsr-tst/label/lh.aparc.annot");
  mri = MRISannotIndex2Seg(mris);
  MRIwrite(mri,"lh.aparc.mgh");
  exit(1);
  }

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

  seg = MRIread(segfile);
  if(seg == NULL) exit(1);

  statmri = MRIread(statfile);
  if(statmri == NULL) exit(1);

  out = MRIcloneBySpace(seg,statmri->nframes);
  if(out == NULL) exit(1);

  for(c=0; c < seg->width; c++){
    for(r=0; r < seg->height; r++){
      for(s=0; s < seg->depth; s++){
	segid = MRIgetVoxVal(seg,c,r,s,0);
	//if(segid == 0) continue; 
	if(segid >= statmri->width){
	  printf("ERROR: %d %d %d segid=%d >= %d\n",c,r,s,segid,statmri->width);
	  exit(1);
	}
	for(f=0; f < statmri->nframes; f++){
	  val = MRIgetVoxVal(statmri,segid,0,0,f);;
	  MRIsetVoxVal(out,c,r,s,f,val);
	}
      }
    }
  }
  MRIwrite(out,outfile);

  printf("mri_stats2seg done\n");
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

    else if (!strcasecmp(option, "--stat")){
      if(nargc < 1) CMDargNErr(option,1);
      statfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--seg")){
      if(nargc < 1) CMDargNErr(option,1);
      segfile = pargv[0];
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--o")){
      if(nargc < 1) CMDargNErr(option,1);
      outfile = pargv[0];
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --stat mristat : stat file in an mri format\n");
  //printf("   --stat-txt stat.txt : text stat file \n");
  printf("   --seg segvol \n");
  printf("   --o out\n");
  printf("\n");
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
  if(statfile == NULL){
    printf("ERROR: need to specify a stat file\n");
    exit(1);
  }
  if(segfile == NULL){
    printf("ERROR: need to specify a seg file\n");
    exit(1);
  }
  if(outfile == NULL){
    printf("ERROR: need to specify an out file\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
  fprintf(fp,"\n");
  fprintf(fp,"%s\n",vcid);
  fprintf(fp,"cwd %s\n",cwd);
  fprintf(fp,"cmdline %s\n",cmdline);
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"user     %s\n",VERuser());

  return;
}
