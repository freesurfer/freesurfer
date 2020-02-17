
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>

#include "mri.h"
#include "macros.h"
#include "error.h"
#include "diag.h"
#include "proto.h"
#include "utils.h"
#include "timer.h"
#include "version.h"
#include "label.h"
#include "cma.h"

int main(int argc, char *argv[]) ;
static int get_option(int argc, char *argv[]) ;

const char *Progname ;
static void usage_exit(int code) ;

int main(int argc, char *argv[]) {
  char   *label_name, *surf_name, *mri_name, *outfname;
  LABEL  *area ;
  int ac, nargs;
  char   **av;

  nargs = handleVersionOption(argc, argv, "mri_dist_surf_label");
    if (nargs && argc - nargs == 1)
      exit (0);
    argc -= nargs;

    Progname = argv[0] ;
    ErrorInit(NULL, NULL, NULL) ;
    DiagInit(NULL, NULL, NULL) ;

    ac = argc ;
    av = argv ;
    for ( ; argc > 1 && ISOPTION(*argv[1]) ; argc--, argv++) {
      nargs = get_option(argc, argv) ;
      argc -= nargs ;
      argv += nargs ;
    }

    fprintf(stderr, "Computing distance from input surface and label points / waypoints.\n") ;
    fprintf(stderr, "Number of arguments %d\n", argc) ;

    if (argc < 5)
      usage_exit(1) ;

  surf_name = argv[1] ; fprintf(stderr, "surface file name is %s\n", surf_name) ; 
  label_name = argv[2]  ; fprintf(stderr, "label file name is %s\n", label_name) ;
  mri_name = argv[3]  ; fprintf(stderr, "mri vol file name is %s\n", mri_name) ;
  outfname = argv[4]  ; fprintf(stderr, "output file name is %s\n", outfname) ;

  FILE *fp = fopen(outfname, "w");

  area = LabelRead(NULL, label_name) ;
  if (!area)
    ErrorExit(ERROR_NOFILE, "%s: could not read label from %s",Progname, label_name) ;

  MRIS *mris = MRISread(surf_name);
  if (!mris) {
    ErrorExit(ERROR_NOFILE, "%s: could not open %s",Progname, surf_name) ;
    return -1;
  }

  MRI *mri = MRIread(mri_name);
  if (!mri) {
    ErrorExit(ERROR_NOFILE, "%s: could not open %s",Progname, mri_name) ;
    return -1;
  }

  MRI *mri_surf_dist;
  // mri_surf_dist = MRIScomputeDistanceToSurface(mris, NULL, 1.0); // original by BF 
  // NOTE: need to use the input volume dimensions in order to compute distance accurately from surfaces converted from gii
  MRI* tmpvol = MRIallocSequence(mri->width, mri->height, mri->depth, MRI_FLOAT, mri->nframes); 
  mri_surf_dist = MRIScomputeDistanceToSurface(mris, tmpvol, 0.5); 

  double xw, yw, zw; 
  fprintf(stderr, "That many points are in the label file: %d\n", area->n_points);
  LabelToVoxel(area, mri_surf_dist, area) ;
  for (int i = 0 ;  i  < area->n_points ; i++) { // going through all points stored in the label file
    double dist = 0;
    
    xw =  area->lv[i].x ;
    yw =  area->lv[i].y ;
    zw =  area->lv[i].z ;
    fprintf(stderr, "Label voxel coordinates: (%f, %f, %f)\n", xw, yw, zw);
    MRIsampleVolume(mri_surf_dist, xw, yw, zw, &dist) ;
    fprintf(stderr, "Distance: %2.2f\n", dist);

    fprintf(fp, "%2.2f\n", dist);
  }

  strcat(outfname, ".wd.mgz");
  MRIwrite(mri_surf_dist,  outfname);
  fclose(fp);
}
 
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  //char *option ;
  return(nargs) ;
}

static void
usage_exit(int code) {
  printf("usage: %s <surface> <label file> <output>\n",  Progname) ;
  exit(code) ;
}
