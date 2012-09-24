/**
 * @file  mris_mesh_subdivide.c
 * @brief 
 *
 * 
 */
/*
 * Original Author: jonathan polimeni
 * CVS Revision Info:
 *    $Author: jonp $
 *    $Date: 2012/09/24 15:24:26 $
 *    $Revision: 1.1 $
 *
 * Copyright © 2011-2012 The General Hospital Corporation (Boston, MA) "MGH"
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

// nint
#include "macros.h"

#include "error.h"
#include "diag.h"
#include "timer.h"

#include "mri.h"
#include "mrisurf.h"
#include "gcamorph.h"

#include "registerio.h"

#include "resample.h"

// string_to_type
#include "mri_identify.h"


char *Progname = NULL;
static int  parse_commandline(int argc, char **argv);
static void print_usage(void) ;
static void usage_exit(void);
static void print_help(void) ;
static void print_version(void) ;

static void argnerr(char *, int) ;

int mris_mesh_subdivide__loop_butterfly(MRI_SURFACE *mris, int iter);


static int debug = 0;

static char vcid[] =
  "$Id: mris_mesh_subdivide.c,v 1.1 2012/09/24 15:24:26 jonp Exp $";


int iter  = 1;

static char *surf_filename = NULL;
static char *newsurf_filename = NULL;


char *basename (char* path)
{
  char *ptr = strrchr (path, '/');
  return ptr ? ptr + 1 : (char*)path;
}


//static char  *subdividemethod_string = "nearest";
static int  subdividemethod = -1;


int main(int argc, char *argv[])
{

  int          nargs = 0;
  int          err;

  MRI_SURFACE  *mris  = NULL;

  int          msec, minutes, seconds ;
  struct timeb start ;


  //  nargs = handle_version_option (argc, argv, "$Id: mris_mesh_subdivide.c,v 1.1 2012/09/24 15:24:26 jonp Exp $", "$Name:  $");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Progname = basename(argv[0]) ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  TimerStart(&start) ;

  argc--;
  argv++;

  if (argc == 0)
  {
    usage_exit();
  }

  parse_commandline(argc, argv);


  //==--------------------------------------------------------------
  // 0) read in a surface

  mris = MRISread(surf_filename) ;
  if (!mris)
  {
    ErrorExit(ERROR_NOFILE, "%s: could not read surface %s", 
              Progname, surf_filename) ;
  }

  //==-----------------------------------------------------------
  // 1) subdivide!

  if ( subdividemethod == -1 )
    {
      err = mris_mesh_subdivide__loop_butterfly(mris, iter);
    }

  //==---------------------------------------------------------
  // N) write out surface

  // TODO: add command line string to surface header (a la HIPS)
  printf("writing to %s\n", newsurf_filename);
  MRISwrite(mris, newsurf_filename);

  MRISfree(&mris);

  msec = TimerStop(&start) ;
  seconds = nint((float)msec/1000.0f) ;
  minutes = seconds / 60 ;
  seconds = seconds % 60 ;
  fprintf(stderr, "%s took %d minutes and %d seconds.\n", 
          Progname, minutes, seconds) ;

  exit(0);
  return(0);
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv)
{
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1)
  {
    usage_exit();
  }

  nargc = argc;
  pargv = argv;
  while (nargc > 0)
  {
    option = pargv[0];
    if (debug)
    {
      printf("%d %s\n",nargc,option);
    }
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))
    {
      print_help() ;
    }
    else if (!strcasecmp(option, "--version"))
    {
      print_version() ;
    }
    else if (!strcasecmp(option, "--debug"))
    {
      debug = 1;
    }

    /*--------------*/

    else if (!strcmp(option, "--surf"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      surf_filename = pargv[0];
      if( !strcmp(surf_filename, "inflated") )
      {
        printf("\nWARNING: do you really want to subdivide the "
               "*inflated* surface?\n\n");
        exit(1);
      }
      nargsused = 1;
    }
    else if (!strcasecmp(option, "--out"))
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      newsurf_filename = pargv[0];
      nargsused = 1;
    }
    else if ( !strcmp(option, "--iterations") )
    {
      if (nargc < 1)
      {
        argnerr(option,1);
      }
      iter = atoi(pargv[0]);
      nargsused = 1;
    }

    nargc -= nargsused;
    pargv += nargsused;
  }
  return(0);
}


/* --------------------------------------------- */
static void argnerr(char *option, int n)
{
  if (n==1)
  {
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  }
  else
  {
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  }
  exit(-1);
}


/* --------------------------------------------- */
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
  printf("   --out <filename>         name for output surface (if does not contain '/'\n");
  printf("                            outputs to same directory as input surface)\n");
  printf("\n");
  printf("   --help        print out information on how to use this program\n");
  printf("   --version     print out version and exit\n");
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void)
{
  print_usage() ;

  printf(
    "This program will subdivide a surface.\n"
  ) ;

  exit(1) ;
}

/* --------------------------------------------- */
static void print_version(void)
{
  printf("%s\n", vcid) ;
  exit(1) ;
}

/* --------------------------------------------- */
int mris_mesh_subdivide__loop_butterfly(MRI_SURFACE *mris, 
                                        int iter)
{
  
  if (!mris)
  {
    ErrorExit(ERROR_NOFILE, "%s: surface invalid", 
              Progname) ;
  }

  printf("VTK: %d\n", iter);

  return(0);
}
