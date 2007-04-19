/**
 * @file  mris_info.cpp
 * @brief Prints out information about a surface file
 *
 */
/*
 * Original Author: Yasunari Tosa
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/04/19 19:03:08 $
 *    $Revision: 1.26 $
 *
 * Copyright (C) 2004-2007,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <iostream>
#include <iomanip>
#include <cstdio>
#include <vector>
#if (__GNUC__ < 3)
#include "/usr/include/g++-3/alloc.h"
#endif
#include <string>
#include <sys/utsname.h>

extern "C" {
#include "fio.h"
#include "mri.h"
#include "utils.h"
#include "gcsa.h"
#include "colortab.h"
#include "diag.h"
#include "transform.h"
#include "mrisurf.h"
#include "mrisutils.h"
#include "version.h"
#include "proto.h"

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

static char vcid[] = 
"$Id: mris_info.cpp,v 1.26 2007/04/19 19:03:08 nicks Exp $";
using namespace std;
char *surffile=NULL, *outfile=NULL, *curvfile=NULL;
char *SUBJECTS_DIR=NULL, *subject=NULL, *hemi=NULL, *surfname=NULL;
int debug = 0;
char tmpstr[2000];
struct utsname uts;
int talairach_flag = 0 ;
MATRIX *XFM=NULL;
int rescale = 0;
double scale=0;
int diag_vno=-1;


/*------------------------------------------------------------*/
int main(int argc, char *argv[]) {
  double InterVertexDistAvg, InterVertexDistStdDev,avgvtxarea;
  char ext[STRLEN] ;
  vector<string> type;
  FILE *fp;
  type.push_back("MRIS_BINARY_QUADRANGLE_FILE");
  type.push_back("MRIS_ASCII_TRIANGLE_FILE");
  type.push_back("MRIS_GEO_TRIANGLE_FILE");
  type.push_back("MRIS_TRIANGULAR_SURFACE=MRIS_ICO_SURFACE");
  type.push_back("MRIS_ICO_FILE");
  type.push_back("MRIS_VTK_FILE");

  if (argc < 2) usage_exit();
  parse_commandline(argc, argv);
  if (surffile == NULL) {
    printf("ERROR: must specify a surface file\n");
    exit(1);
  }

  // Check whether it's a gcs file. If so, just print ctab
  if (!stricmp(FileNameExtension(surffile, ext), "gcs")) {
    GCSA *gcsa = GCSAread(surffile) ;
    if (!gcsa) {
      cerr << "could not open " << surffile << endl;
      return -1;
    }
    printf("GCSA file %s opened\n", surffile) ;
    if (gcsa->ct != NULL) {
      printf("color table:\n") ;
      CTABprintASCII(gcsa->ct,stdout) ;
    }
    return(0) ;
  }

  // Open as a surface
  MRIS *mris = MRISread(surffile);
  if (!mris) {
    cerr << "could not open " << surffile << endl;
    return -1;
  }
  MRIScomputeMetricProperties(mris) ;

  if (diag_vno >= 0)
  {
    printf("v %d: (%2.1f, %2.1f, %2.1f)\n", 
           diag_vno,
           mris->vertices[diag_vno].x, 
           mris->vertices[diag_vno].y, 
           mris->vertices[diag_vno].z) ;
    return(0);
  }

  // attempt to load curvature info, which has the side-effect of checking
  // that the number of vertices is the same (exiting with error if not)
  if (curvfile) {
    printf("\n");
    if (MRISreadCurvatureFile(mris,curvfile)) {
      printf("\n");
      exit(1);
    }
  }

  if (rescale) {
    if (mris->group_avg_surface_area == 0) {
      printf("ERROR: cannot rescale a non-group surface\n");
      exit(1);
    }
    scale = MRISrescaleMetricProperties(mris);
    //scale = sqrt((double)mris->group_avg_surface_area/mris->total_area);
    //printf("scale = %lf\n",scale);
    //MRISscale(mris,scale);
  }

  if (talairach_flag) {
    if (! subject) {
      printf("ERROR: need --s with --t\n");
      exit(1);
    }
    XFM = DevolveXFM(subject, NULL, NULL);
    if (XFM == NULL) exit(1);
    printf("Applying talairach transform\n");
    MatrixPrint(stdout,XFM);
    MRISmatrixMultiply(mris,XFM);
  }

  cout << "SURFACE INFO ======================================== " << endl;
  cout << "type        : " << type[mris->type].c_str() << endl;
  if (mris->type == MRIS_BINARY_QUADRANGLE_FILE) {
    FILE *fp = fopen(surffile, "rb") ;
    int magic;
    fread3(&magic, fp) ;
    if (magic == QUAD_FILE_MAGIC_NUMBER)
      cout << "              QUAD_FILE_MAGIC_NUMBER" << endl;
    else if (magic == NEW_QUAD_FILE_MAGIC_NUMBER)
      cout << "              NEW_QUAD_FILE_MAGIC_NUMBER" << endl;
    fclose(fp);
  }

  InterVertexDistAvg    = mris->avg_vertex_dist;
  InterVertexDistStdDev = mris->std_vertex_dist;
  avgvtxarea = mris->avg_vertex_area;

  cout << "num vertices: " << mris->nvertices << endl;
  cout << "num faces   : " << mris->nfaces << endl;
  cout << "num strips  : " << mris->nstrips << endl;
  cout << "surface area: " << mris->total_area << endl;
  printf("AvgVtxArea       %lf\n",avgvtxarea);
  printf("AvgVtxDist       %lf\n",InterVertexDistAvg);
  printf("StdVtxDist       %lf\n",InterVertexDistStdDev);

  if (mris->group_avg_surface_area > 0) {
    cout << "group avg surface area: " << mris->group_avg_surface_area << endl;
  }
  cout << "ctr         : (" << mris->xctr 
       << ", " << mris->yctr 
       << ", " << mris->zctr 
       << ")" << endl;
  cout << "vertex locs : " 
       << (mris->useRealRAS ? "scannerRAS" : "surfaceRAS") << endl;
  if (mris->lta) {
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
  if (argc > 2) {
    int vno = atoi(argv[2]) ;
    VERTEX *v ;
    v= &mris->vertices[vno] ;
    printf("mris[%d] = (%2.1f, %2.1f, %2.1f)\n", vno, v->x, v->y, v->z) ;
  }

  uname(&uts);

  // Open an output file to capture values
  if (outfile != NULL) {
    fp = fopen(outfile,"w");
    if (fp == NULL) {
      printf("ERROR: cannot open %s\n",outfile);
      exit(1);
    }
  } else fp = stdout;

  fprintf(fp,"mris_info\n");
  fprintf(fp,"creationtime %s\n",VERcurTimeStamp());
  fprintf(fp,"sysname  %s\n",uts.sysname);
  fprintf(fp,"hostname %s\n",uts.nodename);
  fprintf(fp,"machine  %s\n",uts.machine);
  fprintf(fp,"surfacefile %s\n",surffile);
  if (SUBJECTS_DIR) fprintf(fp,"SUBJECTS_DIR  %s\n",SUBJECTS_DIR);
  if (subject)      fprintf(fp,"subject       %s\n",subject);
  if (hemi)         fprintf(fp,"hemi          %s\n",hemi);
  if (surfname)     fprintf(fp,"surfname      %s\n",surfname);
  fprintf(fp,"hemicode    %d\n",mris->hemisphere);
  fprintf(fp,"talairach_flag  %d\n",talairach_flag);
  fprintf(fp,"rescale     %lf\n",scale);
  fprintf(fp,"nvertices   %d\n",mris->nvertices);
  fprintf(fp,"nfaces      %d\n",mris->nfaces);
  fprintf(fp,"total_area  %f\n",mris->total_area);
  if (mris->group_avg_surface_area > 0)
    fprintf(fp,"group_avg_surf_area  %f\n",mris->group_avg_surface_area);
  printf("group_avg_vtxarea_loaded %d\n",mris->group_avg_vtxarea_loaded);
  fprintf(fp,"avgvtxarea  %lf\n",avgvtxarea);
  fprintf(fp,"avgvtxdist  %lf\n",InterVertexDistAvg);
  fprintf(fp,"stdvtxdist  %lf\n",InterVertexDistStdDev);
  fprintf(fp,"vtx0xyz   %f %f %f\n",mris->vertices[0].x,mris->vertices[0].y,
          mris->vertices[0].z);

  if (outfile != NULL) fclose(fp);

  MRISfree(&mris);
}


/* --------------------------------------------- */
static int parse_commandline(int argc, char **argv) {
  int  nargc , nargsused;
  char **pargv, *option ;

  if (argc < 1) usage_exit();

  nargc   = argc-1;
  pargv = &argv[1];
  while (nargc > 0) {

    option = pargv[0];
    if (debug) printf("%d %s\n",nargc,option);
    nargc -= 1;
    pargv += 1;

    nargsused = 0;

    if (!strcasecmp(option, "--help"))  print_help() ;
    else if (!strcasecmp(option, "--version")) print_version() ;
    else if (!strcasecmp(option, "--debug"))   debug = 1;
    else if (!strcasecmp(option, "--tal"))   talairach_flag = 1;
    else if (!strcasecmp(option, "--t"))   talairach_flag = 1;
    else if (!strcasecmp(option, "--r"))   rescale = 1;
    else if ( !strcmp(option, "--o") ) {
      if (nargc < 1) argnerr(option,1);
      outfile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--v") ) {
      if (nargc < 1) argnerr(option,1);
      diag_vno = atoi(pargv[0]);
      nargsused = 1;
    } else if ( !strcmp(option, "--c") ) {
      if (nargc < 1) argnerr(option,1);
      curvfile = pargv[0];
      nargsused = 1;
    } else if ( !strcmp(option, "--s") ) {
      if (nargc < 3) argnerr(option,3);
      subject  = pargv[0];
      hemi     = pargv[1];
      surfname = pargv[2];
      SUBJECTS_DIR = getenv("SUBJECTS_DIR");
      if (SUBJECTS_DIR==NULL) {
        printf("ERROR: SUBJECTS_DIR not defined in environment\n");
        exit(1);
      }
      sprintf(tmpstr,"%s/%s/surf/%s.%s",SUBJECTS_DIR,subject,hemi,surfname);
      surffile = strcpyalloc(tmpstr);
      nargsused = 3;
    } else {
      if (NULL == surffile) surffile = option;
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
  printf("USAGE: %s [options] <surfacefile>\n",Progname) ;
  printf("\nOptions:\n");
  printf("  --o outfile : save some data to outfile\n");
  printf("  --s subject hemi surfname : instead of surfacefile\n");
  printf("  --t : apply talairach xfm before reporting info\n");
  printf("  --r : rescale group surface so metrics same as "
         "avg of individuals\n");
  printf("  --v vnum : print out vertex information for vertex vnum\n") ;
  printf("  --c curvfile : load the specified curvature file, which\n");
  printf("                 has the side-effect of checking that the\n");
  printf("                 number of vertices matches the surface, and\n");
  printf("                 exiting with error if not.\n");
  printf("\n");
  printf("  --version   : print version and exits\n");
  printf("  --help      : no clue what this does\n");
  printf("\n");
}

/* --------------------------------------------- */
static void print_help(void) {
  print_usage() ;
  printf("\n");
  printf("%s\n", vcid) ;
  printf("\n");
  printf("Prints out information about a surface file.\n");
  printf("\n");
  exit(1);
}

/* --------------------------------------------- */
static void argnerr(char *option, int n) {
  if (n==1)
    fprintf(stderr,"ERROR: %s flag needs %d argument\n",option,n);
  else
    fprintf(stderr,"ERROR: %s flag needs %d arguments\n",option,n);
  exit(-1);
}

/* --------------------------------------------- */
static void print_version(void) {
  printf("%s\n", vcid) ;
  exit(1) ;
}
