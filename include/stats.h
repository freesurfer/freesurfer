#ifndef STATS_H
#define STATS_H

#include "minc_volume_io.h"
#include "matrix.h"
#include "mri.h"

typedef struct
{
  float   in_plane_res ;
  float   slice_thickness ;
  float   brightness_scale ;
  MATRIX  *fmri2mri ;
  MATRIX  *mri2fmri ;
  char    name[100] ;            /* subject's name */
} fMRI_REGISTRATION, fMRI_REG ;

#define STAT_VOLUMES   2  /* avg, std, dof for avg, dof for std */

#define DEFAULT_RESOLUTION  16

/*
  on disk in .bfloat files, the volumes are stored with each slice
  in a separate file. The files contain all the mean images (one per
  time point) followed by all the standard deviation images (one per
  time point) for that slice for each event type.
  */ 

#define ALLOC_MEANS        0x0001
#define ALLOC_STDS         0x0002
#define ALLOC_DOFS         0x0004

#define MAX_EVENTS       100

#define mri_pvals        mri_avgs

typedef struct
{
  /* from .dat file */
  int       slice_width ;   /* in voxels */
  int       slice_height ;
  int       nslices ;       /* depth in voxels */

  float     mean_dofs[MAX_EVENTS] ; /* the # dof for the means */
  float     std_dofs[MAX_EVENTS] ;  /* the # dof for the std devs */

  /* from register.dat file */
  fMRI_REG  *reg ;
  
  int       voltype; /* 0 = raw, 1 = selavg, 2 = selxavg */

  /* stuff from the .dat file */
  float     tr ;
  float     timewindow ;
  float     prestim ;
  int       nevents ;     /* # of event types (conditions) */
  int       time_per_event ;/* # of time points per event (timewindow/nbins) */

  MRI       *mri_avgs[MAX_EVENTS] ;
  MRI       *mri_stds[MAX_EVENTS] ;
  MRI       *mri_avg_dofs[MAX_EVENTS] ;
  MRI       *mri_std_dofs[MAX_EVENTS] ;

  General_transform transform ;
  Transform         *linear_transform ;
  Transform         *inverse_linear_transform ;
} STAT_VOLUME, SV  ;
  
/* can't include this before structure, as mrisurf.h includes this file. */
#include "mrisurf.h"
#include "resample.h"

/* This is so applications can specify different xforms */
/* StatReadTransform() will read in this file */
#ifdef _STATS_SRC
  char *stats_talxfm = "talairach.xfm";
  int  statnorm_float2int = FLT2INT_ROUND;
  int  stats_fixxfm = 0;
  #undef _STATS_SRC
#else
  extern char *stats_talxfm;
  extern int  statnorm_float2int;
  extern int  stats_fixxfm;
#endif

SV        *StatReadVolume(char *prefix);
SV        *StatReadTalairachVolume(char *prefix, 
                                   char *xform_fname,
                                   char *subject_name) ;
fMRI_REG  *StatReadRegistration(char *fname) ;
fMRI_REG  *StatReadTalairachRegistration(char *fname, 
                                         char *subject_name) ;
int       StatWriteVolume(SV *sv, 
                          char *prefix) ;
int       StatWriteRegistration(fMRI_REG *reg, 
                                char *fname) ;
int       StatFreeRegistration(fMRI_REG **preg) ;
int       StatFree(SV **psv) ;
SV        *StatAllocVolume(SV *sv, 
                           int nevents, 
                           int width, 
                           int height,
                           int nslices, 
                           int time_points, 
                           int which) ;
SV        *StatAllocStructuralVolume(SV *sv, 
                                     float fov, 
                                     float resolution,
                                     char *name) ;
int       StatAccumulateTalairachVolume(SV *sv_tal, 
                                        SV *sv) ;
int       StatAccumulateSurfaceVolume(SV *sv_tal, 
                                      SV *sv, 
                                      MRI_SURFACE *mris) ;
int       StatReadTransform(STAT_VOLUME *sv, 
                            char *name) ;
int       StatVolumeExists(char *prefix) ;

#endif
