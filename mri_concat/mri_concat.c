// mri_concat.c
// $Id: mri_concat.c,v 1.1 2005/07/13 03:45:28 greve Exp $

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "macros.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "error.h"
#include "diag.h"
#include "mri.h"
#include "mri2.h"
#include "version.h"

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

static char vcid[] = "$Id: mri_concat.c,v 1.1 2005/07/13 03:45:28 greve Exp $";
char *Progname = NULL;
int debug = 0;
char *inlist[100];
int ninputs = 0;
char *out = NULL;
MRI *mritmp, *mritmp0, *mriout;

/*--------------------------------------------------*/
int main(int argc, char **argv)
{
  int nargs, nthin, nframestot=0, nr=0,nc=0,ns=0, fout;
  int r,c,s,f;
  double v;

  /* rkt: check for and handle version tag */
  nargs = handle_version_option (argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1) exit (0);
  argc -= nargs;

  Progname = argv[0] ;
  argc --;
  argv++;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  if(argc == 0) usage_exit();

  parse_commandline(argc, argv);
  check_options();
  dump_options(stdout);

  for(nthin = 0; nthin < ninputs; nthin++){
    mritmp = MRIreadHeader(inlist[nthin],MRI_VOLUME_TYPE_UNKNOWN);
    if(mritmp == NULL){
      printf("ERROR: reading %s\n",inlist[nthin]);
      exit(1);
    }
    if(nthin == 0) {
      nc = mritmp->width;
      nr = mritmp->height;
      ns = mritmp->depth;
    }
    if(mritmp->width != nc || mritmp->height != nr ||
       mritmp->depth != ns){
      printf("ERROR: dimension mismatch between %s and %s\n",
	     inlist[0],inlist[nthin]);
      exit(1);
    }
    nframestot += mritmp->nframes;
    MRIfree(&mritmp);
  }
  printf("nframestot = %d\n",nframestot); 

  mriout = MRIallocSequence(nc,nr,ns,MRI_FLOAT,nframestot);
  if(mriout == NULL) exit(1);

  fout = 0;
  for(nthin = 0; nthin < ninputs; nthin++){
    mritmp = MRIread(inlist[nthin]);
    for(f=0; f < mritmp->nframes; f++){
      for(c=0; c < nc; c++){
	for(r=0; r < nr; r++){
	  for(s=0; s < ns; s++){
	    v = MRIgetVoxVal(mritmp,c,r,s,f);
	    MRIsetVoxVal(mriout,c,r,s,fout,v);
	  }
	}
      }
      fout++;
    }
    MRIfree(&mritmp);
  }

  printf("Writing to %s\n",out);
  MRIwrite(mriout,out);

  return(0);
}
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/
/*-----------------------------------------------------------------*/

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
    else if ( !strcmp(option, "--i") ) {
      if(nargc < 1) argnerr(option,1);
      inlist[ninputs] = pargv[0];
      ninputs ++;
      nargsused = 1;
    }
    else if ( !strcmp(option, "--o") ) {
      if(nargc < 1) argnerr(option,1);
      out = pargv[0];
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
  printf("USAGE: %s \n",Progname) ;
  printf("\n");
  printf("   --i invol <--i invol ...> \n");
  printf("   --o out \n");
  printf("\n");
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

  printf("Concatenates input data sets.\n");

  exit(1) ;
}
/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
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
  if(ninputs == 0){
    printf("ERROR: no inputs specified\n");
    exit(1);
  }
  if(out == NULL){
    printf("ERROR: no output specified\n");
    exit(1);
  }
  return;
}

/* --------------------------------------------- */
static void dump_options(FILE *fp)
{
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
