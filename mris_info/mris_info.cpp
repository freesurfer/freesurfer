//
// mris_info.cpp
//

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>

extern "C" {
#include "fio.h"
#include "mri.h"
#include "utils.h"
#include "gcsa.h"
#include "colortab.h"
#include "transform.h"
#include "mrisurf.h"

  char *Progname = "mris_info";
}

static int  parse_commandline(int argc, char **argv);
static void print_usage(void);
static void usage_exit(void);
static void print_help(void);
static void argnerr(char *option, int n);
static void print_version(void);

// copied from mrisurf.c
#define QUAD_FILE_MAGIC_NUMBER      (-1 & 0x00ffffff)
#define TRIANGLE_FILE_MAGIC_NUMBER  (-2 & 0x00ffffff)
#define NEW_QUAD_FILE_MAGIC_NUMBER  (-3 & 0x00ffffff)

static char vcid[] = "$Id";
using namespace std;
char *surffile=NULL, *outfile=NULL;
int debug = 0;
int PrintNVertices = 0;
int PrintNFaces    = 0;
int PrintNStrips = 0;


/*------------------------------------------------------------*/
int main(int argc, char *argv[])
{
  char ext[STRLEN] ;
  vector<string> type;
  FILE *fp;
  type.push_back("MRIS_BINARY_QUADRANGLE_FILE");
  type.push_back("MRIS_ASCII_TRIANGLE_FILE");
  type.push_back("MRIS_GEO_TRIANGLE_FILE");
  type.push_back("MRIS_TRIANGULAR_SURFACE=MRIS_ICO_SURFACE");
  type.push_back("MRIS_ICO_FILE");
  type.push_back("MRIS_VTK_FILE");

  if(argc < 2) usage_exit();
  parse_commandline(argc, argv);
  if(surffile == NULL){
    printf("ERROR: must specify a surface file\n");
    exit(1);
  }

  // Check whether it's a gcs file. If so, just print ctab
  if (!stricmp(FileNameExtension(surffile, ext), "gcs")) {
    GCSA *gcsa = GCSAread(surffile) ;
    if (!gcsa)	{
      cerr << "could not open " << surffile << endl;
      return -1;
    }
    printf("GCSA file %s opened\n", surffile) ;
    if (gcsa->ct != NULL)	{
      printf("color table:\n") ;
      CTABprint(stdout, gcsa->ct) ;
    }
    return(0) ;
  }
  
  // Open as a surface
  MRIS *mris = MRISread(surffile);
  if (!mris) {
    cerr << "could not open " << surffile << endl;
    return -1;
  }

  // Open an output file to capture values
  if(outfile != NULL){
    fp = fopen(outfile,"w");
    if(fp == NULL){
      printf("ERROR: cannot open %s\n",outfile);
      exit(1);
    }
  }
  else fp = stdout;

  if(PrintNVertices) fprintf(fp,"nvertices %d\n",mris->nvertices);
  if(PrintNFaces)    fprintf(fp,"nfaces    %d\n",mris->nfaces);
  if(PrintNStrips)   fprintf(fp,"nstrips   %d\n",mris->nstrips);

  if(outfile != NULL) fclose(fp);

  
  cout << "SURFACE INFO ======================================== " << endl;
  cout << "type        : " << type[mris->type].c_str() << endl;
  if (mris->type == MRIS_BINARY_QUADRANGLE_FILE){
    FILE *fp = fopen(surffile, "rb") ;
    int magic;
    fread3(&magic, fp) ;
    if(magic == QUAD_FILE_MAGIC_NUMBER)
      cout << "              QUAD_FILE_MAGIC_NUMBER" << endl;
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER) 
      cout << "              NEW_QUAD_FILE_MAGIC_NUMBER" << endl;
    fclose(fp);
  }

  cout << "num vertices: " << mris->nvertices << endl;
  cout << "num faces   : " << mris->nfaces << endl;
  cout << "num strips  : " << mris->nstrips << endl;
  cout << "ctr         : (" << mris->xctr << ", " << mris->yctr << ", " << mris->zctr << ")" << endl;
  cout << "vertex locs : " << (mris->useRealRAS ? "scannerRAS" : "surfaceRAS") << endl;
  if (mris->lta)
    {
      cout << "talairch.xfm: " << endl;
      MatrixPrint(stdout, mris->lta->xforms[0].m_L);
      cout << "surfaceRAS to talaraiched surfaceRAS: " << endl;
      MatrixPrint(stdout, mris->SRASToTalSRAS_);
      cout << "talairached surfaceRAS to surfaceRAS: " << endl;
      MatrixPrint(stdout, mris->TalSRASToSRAS_);
    }
  vg_print(&mris->vg); 
  
  {
    int i ;
    for (i = 0 ; i < mris->ncmds ; i++)
      printf("cmd[%d]: %s\n", i, mris->cmdlines[i]) ;
  }
  MRISfree(&mris);
}

/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if(argc < 1) usage_exit();

  nargc   = argc-1;
  pargv = &argv[1];
  while(nargc > 0){

    option = pargv[0];
    if(debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--nvertices")) PrintNVertices = 1;
    else if (!strcasecmp(option, "--nfaces"))    PrintNFaces = 1;
    else if (!strcasecmp(option, "--nstrips"))   PrintNStrips = 1;
    else if ( !strcmp(option, "--o") ) {
      if(nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    }
    else{
      surffile = option;
      nargsused = 1;
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
  printf("USAGE: %s   surfacefile\n",Progname) ;
  printf("  --o outfile : save some data to outfile\n");
  printf("  --nvertices : print nverticies\n");
  printf("  --nfaces    : print nfaces\n");
  printf("\n");
  printf("  --version   : print version and exits\n");
  printf("  --help      : no clue what this does\n");
  printf("\n");
}
/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;
  printf("\n");  
  printf("%s\n", vcid) ;
  printf("\n");
  printf("Prints out information about a surface file.\n");
  printf("\n");
  exit(1);
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
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}
