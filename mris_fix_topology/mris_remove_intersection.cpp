/**
 * @brief removes surface intersections
 *
 */
/*
 * Original Author: Bruce Fischl
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "macros.h"
#include "error.h"
#include "tags.h"
#include "diag.h"
#include "proto.h"
#include "mrisurf.h"
#include "mri.h"
#include "macros.h"
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "surfcluster.h"
#include "mrisurf_metricProperties.h"


int main(int argc, char *argv[]) ;

static int  get_option(int argc, char *argv[]) ;
static void usage_exit(void) ;
static void print_usage(void) ;
static void print_help(void) ;
static void print_version(void) ;

const char *Progname ;
int FillHoles = 0;

int main(int argc, char *argv[])
{
  char         **av, *in_surf_fname, *out_fname ;
  int          ac, nargs, msec ;
  MRI_SURFACE  *mris ;
  Timer then ;


  std::string cmdline = getAllInfo(argc, argv, "mris_remove_intersection");

  nargs = handleVersionOption(argc, argv, "mris_remove_intersection");
  if (nargs && argc - nargs == 1)
  {
    exit (0);
  }
  argc -= nargs;

  Gdiag = DIAG_SHOW ;

  then.reset() ;
  Progname = argv[0] ;
  ErrorInit(NULL, NULL, NULL) ;
  DiagInit(NULL, NULL, NULL) ;

  ac = argc ;
  av = argv ;
  for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++)
  {
    nargs = get_option(argc, argv) ;
    argc -= nargs ;
    argv += nargs ;
  }

  if (argc < 3)
  {
    usage_exit() ;
  }

  in_surf_fname = argv[1] ;
  out_fname = argv[2] ;

  mris = MRISread(in_surf_fname) ;
  if (!mris)
    ErrorExit(ERROR_NOFILE, "%s: could not read surface file %s",
              Progname, in_surf_fname) ;

  MRISaddCommandLine(mris, cmdline) ;

  mrisMarkIntersections(mris,FillHoles);
  int n, nintersections=0;
  for(n=0; n < mris->nvertices; n++) if(mris->vertices[n].marked) nintersections++;
  printf("Found %d intersections\n",nintersections);

  // MRISsetNeighborhoodSizeAndDist(mris, 2) ;
  MRISremoveIntersections(mris,FillHoles) ;

  printf("writing corrected surface to %s\n", out_fname) ;
  MRISwrite(mris, out_fname) ;
  msec = then.milliseconds() ;
  fprintf(stderr, "intersection removal took %2.2f hours\n",
          (float)msec/(1000.0f*60.0f*60.0f));

  exit(0) ;

  return(0) ;  /* for ansi */
}

/*----------------------------------------------------------------------
  Parameters:

  Description:
  ----------------------------------------------------------------------*/
static int
get_option(int argc, char *argv[])
{
  int  nargs = 0 ;
  char *option ;

  option = argv[1] + 1 ;            /* past '-' */
  if (!stricmp(option, "-help")||!stricmp(option, "-usage")) print_help() ;
  else if (!stricmp(option, "-version")) print_version() ;
  else if (!stricmp(option, "fill-holes"))  FillHoles = 1;
  else if (!stricmp(option, "no-fill-holes"))  FillHoles = 0;
  else if (!stricmp(option, "map"))
  {
    MRIS *surf = MRISread(argv[2]);
    if(surf==NULL) exit(1);
    if(argc > 4) {
      double projdistmm = 0;
      sscanf(argv[4],"%lf",&projdistmm);
      printf("Projecting surface %g mm\n",projdistmm);
      for(int n=0; n < surf->nvertices; n++){
	VERTEX *v = &(surf->vertices[n]);
	v->x += projdistmm*v->nx;
	v->y += projdistmm*v->ny;
	v->z += projdistmm*v->nz;
      }
    }
    mrisMarkIntersections(surf,FillHoles);
    int n, nintersections=0;
    for(n=0; n < surf->nvertices; n++){
      if(surf->vertices[n].marked) nintersections++;
    }
    printf("Found %d intersections\n",nintersections);
    #if 0
    printf("Filling any holes\n");
    SURFCLUSTERSUM *SurfClustList;
    int nClusters;
    SurfClustList = sclustMapSurfClusters(surf,0.5,-1,1,0,&nClusters,NULL,NULL);
    printf("Found %d clusters\n",nClusters);
    for(n=0; n < surf->nvertices; n++){
      if(surf->vertices[n].undefval > 1) surf->vertices[n].marked = 1;
    }
    #endif
    MRI *mri = MRIcopyMRIS(NULL,surf,0,"marked");
    int err = MRIwrite(mri,argv[3]);
    exit(err);
  }
  else switch (toupper(*option))
    {
    case 'V':
      Gdiag_no = atoi(argv[2]) ;
      nargs = 1 ;
      break ;
    case '?':
    case 'H':
    case 'U':
      print_usage() ;
      exit(1) ;
      break ;
    default:
      fprintf(stderr, "unknown option %s\n", argv[1]) ;
      exit(1) ;
      break ;
    }

  return(nargs) ;
}

static void
usage_exit(void)
{
  print_usage() ;
  exit(1) ;
}

#include "mris_remove_intersection.help.xml.h"
static void
print_usage(void)
{
  outputHelpXml(mris_remove_intersection_help_xml,mris_remove_intersection_help_xml_len);
}

static void
print_help(void)
{
  print_usage() ;
  exit(1) ;
}

static void
print_version(void)
{
  fprintf(stderr, "%s\n", getVersion().c_str()) ;
  exit(1) ;
}

