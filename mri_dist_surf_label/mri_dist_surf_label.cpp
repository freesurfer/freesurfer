
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
  char   *label_name, *surf_name, *mri_name;
  LABEL  *area ;
  int ac, nargs;
  char   **av;
  float wdist; 
  MATRIX *m_vox2vox;
  VECTOR *v1, *v2;

  nargs = handle_version_option (argc, argv, "$Id: mri_label_vals.c,v 1.16 2015/08/24 18:22:05 fischl Exp $", "$Name:  $");
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

    if (argc < 3)
      usage_exit(1) ;

  surf_name = argv[1] ; fprintf(stderr, "surface file name is %s\n", surf_name) ; 
  label_name = argv[2]  ; fprintf(stderr, "label file name is %s\n", label_name) ;
  mri_name = argv[3]  ; fprintf(stderr, "mri vol file name is %s\n", mri_name) ;

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
  mri_surf_dist = MRIScomputeDistanceToSurface(mris, NULL, 0.5); 
  MRIwrite(mri_surf_dist, "./wd.mgz");
  m_vox2vox = MRIgetVoxelToVoxelXform(mri, mri_surf_dist);
  v1 = VectorAlloc(4, MATRIX_REAL);
  v2 = VectorAlloc(4, MATRIX_REAL);
  VECTOR_ELT(v1, 4) = 1.0;
  VECTOR_ELT(v2, 4) = 1.0;

  double xw, yw, zw, xv, yv, zv;

  fprintf(stderr, "That many points are in the label file: %d\n", area->n_points);
  for (int i = 0 ;  i  < area->n_points ; i++) { // going through all points stored in the label file

    xw =  area->lv[i].x ;
    yw =  area->lv[i].y ;
    zw =  area->lv[i].z ;
    fprintf(stderr, "Label coordinates: (%f, %f, %f)\n", xw, yw, zw);
    //MRIworldToVoxel(mri, xw, yw,  zw, &xv, &yv, &zv) ;
    MRIsurfaceRASToVoxel(mri, xw, yw, zw, &xv, &yv, &zv);
    //fprintf(stderr, "Mapped label coordinates: (%f, %f, %f)\n", xv, yv, zv);

    V3_X(v1) = xv;
    V3_Y(v1) = yv;
    V3_Z(v1) = zv;
    MatrixMultiply(m_vox2vox, v1, v2);
    fprintf(stderr, "Mapped label coordinates: (%f, %f, %f)\n", V3_X(v1), V3_Y(v1), V3_Z(v1));
    fprintf(stderr, "Mapped multiplied label coordinates: (%f, %f, %f)\n", V3_X(v2), V3_Y(v2), V3_Z(v2));

    if (MRIindexNotInVolume(mri_surf_dist, V3_X(v2), V3_Y(v2), V3_Z(v2))) 
      {
	fprintf(stderr, "Distance cannot be computed! Point not in dist volume\n");
	continue;
      }
    wdist = MRIgetVoxVal(mri_surf_dist, V3_X(v2), V3_Y(v2), V3_Z(v2), 0);
    fprintf(stderr, "Distance: %f\n", wdist);
  }
  LabelToVoxel(area, mri_surf_dist, area) ;
  for (int i = 0 ;  i  < area->n_points ; i++) { // going through all points stored in the label file
    double dist ;

    xw =  area->lv[i].x ;
    yw =  area->lv[i].y ;
    zw =  area->lv[i].z ;
    fprintf(stderr, "Label voxel coordinates: (%f, %f, %f)\n", xw, yw, zw);
    MRIsampleVolume(mri_surf_dist, xw, yw, zw, &dist) ;
    fprintf(stderr, "Distance: %2.2f\n", dist);
  }

}
 
static int
get_option(int argc, char *argv[]) {
  int  nargs = 0 ;
  //char *option ;

  return(nargs) ;
}

static void
usage_exit(int code) {
  printf("usage: %s <surface> <label file>\n",  Progname) ;
  exit(code) ;
}
