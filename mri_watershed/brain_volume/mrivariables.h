/*
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


//
// mrivariables.h
//

#include "mri.h"
#include "mrisurf.h"

#ifndef c_mrivariables_h
#define c_mrivariables_h

typedef struct Cell {
  unsigned char type;
  void * next;
}
Cell;

typedef struct Basincell {
  unsigned char depth;
  unsigned long size;
  unsigned long ambiguous;
}
BasinCell;

typedef struct Bound {
  unsigned char x,y,z,val;
  struct Bound *next;
}
Bound;

typedef unsigned char Coord[3];


typedef struct STRIP_PARMS {
  /*  float  fill_level;*/
  int template_deformation;
  /*to save the surfaces into the output volume*/
  int surf_dbg;
  /*to write out the brain surface into the file surfname*/
  /*the brain surface is shrank inward of h_shrk mm*/
  int brainsurf;
  /*to write out all the BEM surfaces : brain, outer and inner skull, scalp*/
  int surf;
  /*to labelize the volume into scalp, skull, csf, white and gray*/
  int label;
  /*to use the atlas validation and correction*/
  int atlas;

  /*specify T1 volume*/
  int T1;
  /*specify the center fo the brain and its radius*/
  int cx,cy,cz,rb;

  char *surfname;
  int h_shk;
  int skull_type;
  int watershed_analyze;
  int threshold_analyze;

  int seed_coord[30][4];
  int nb_seed_points/*=0*/;
  unsigned char hpf;

  int manual_params;
  int manual_CSF_MAX,manual_TRANSITION_INTENSITY,manual_GM_INTENSITY;

}
STRIP_PARMS ;


typedef struct {
  float direction[26][3];
  MRIS *mrisphere,*mris,*mris_curv,*mris_var_curv,*mris_dCOG
  ,*mris_var_dCOG;

  double xCOG,yCOG,zCOG,rad_Brain;
  double xsCOG,ysCOG,zsCOG;
  int i_global_min/*=0*/,j_global_min,k_global_min,int_global_min;
  unsigned long estimated_size/*=0*/,main_basin_size/*=0*/;
  unsigned long brain_size /*=0*/;
  unsigned long basinnumber,basinsize;

  MRI *mri_src,*mri_dst,*mri_orig;
  int width,height,depth;

  unsigned char Imax;
  int WM_INTENSITY,WM_VARIANCE,WM_HALF_MAX,WM_HALF_MIN,WM_MAX,WM_MIN;
  int CSF_intensity,CSF_HALF_MAX,CSF_MAX,CSF_MIN;
  int GM_MIN, GM_INTENSITY,TRANSITION_INTENSITY;


  unsigned long gmnumber[256];

  /*file saving - not used */
#if OUTPUT_CURVES
  FILE *fout;
#endif
#if OUTPUT_SURFACES
  FILE *fsvout,*fsfout;
#endif

  Bound *Bound1,*Bound2;

  Cell *** Basin;

  Coord** Table[256];

  unsigned char intbasin[256];
  unsigned long tabdim[256];
  unsigned long sqrdim[256];
  unsigned long count[256];

  Coord* T1Table;
  long T1nbr;

  int decision;
  float scale;

  int atlas;
  int validation;
  int verbose_mode;

}
MRI_variables;

#endif
