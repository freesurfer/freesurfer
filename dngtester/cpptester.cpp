#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#define export

#ifdef HAVE_ITK_LIBS
#include <vnl/vnl_matrix.h>
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#endif

// all other software are all in "C"
#ifdef __cplusplus
extern "C"
{
#endif

#include "error.h"
#include "macros.h"
#include "mri.h"
#include "transform.h"
#include "resample.h"
#include "registerio.h"
#include "version.h"

#ifdef __cplusplus
}
#endif

using namespace std;

static void printUsage(void);
static bool parseCommandLine(int argc, char *argv[]);
static int parseNextCommand(int argc, char *argv[]);

static char vcid[] =
    "$Id: lta_convert.cpp,v 1.10 2016/08/09 02:11:11 zkaufman Exp $";
char *Progname = NULL;
float myrand(float f);

int main(int argc, char *argv[])
{
  cout << vcid << endl << endl;

  // Default initialization
  int nargs = handle_version_option(argc, argv, vcid, "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit(0);
  }
  argc -= nargs;
  Progname = argv[0];
  argc--;
  argv++;
  ErrorInit(NULL, NULL, NULL);

#ifndef HAVE_ITK_LIBS
  fprintf(stderr, "%s:%d dngtester needs ITK, not available\n",__FILE__,__LINE__);
  exit(1);
#else
  while(1){
    vnl_matrix<float> M(3000,3000);
    vnl_matrix<float> Q =  M.apply(&myrand);
  }
  vnl_matrix<float> N(3,2);
  //vnl_matrix<float> Q =  M.apply(&myrand);
  //Q.print(cout);
  //vnl_matrix<float> Qinv = vnl_inverse(Q);
  //vnl_matrix<float> Qinv = vnl_matrix_inverse<float>(Q);
  //Qinv.print(cout);
  //N.fill(3);
  //M.fill_diagonal(2);
  //vnl_matrix<float> S = M*N;
  //N.print(cout);
  //printf("is_zero %d\n",M.is_zero());
  //printf("is_ident %d\n",M.is_identity());
#endif
    
  // Parse command line
  if (!parseCommandLine(argc, argv))
  {
    //printUsage();
    exit(1);
  }
}

static void printUsage(void)
{
  printf("usage\n");
}

static int parseNextCommand(int argc, char *argv[])
{
  int nargs = 0;
  char *option;

  option = argv[0] + 1;                     // remove '-'
  if (option[0] == '-')
  {
    option = option + 1;  // remove second '-'
  }
  fflush(stdout);

  return (nargs);
}

static bool parseCommandLine(int argc, char *argv[])
{
  int nargs;
  int inputargs = argc;
  for (; argc > 0 && ISOPTION(*argv[0]); argc--, argv++)
  {
    nargs = parseNextCommand(argc, argv);
    argc -= nargs;
    argv += nargs;
  }

  if (inputargs == 0)
  {
    printUsage();
    exit(1);
  }

  return true;
}

float myrand(float f)
{
  return((float)drand48());
}
