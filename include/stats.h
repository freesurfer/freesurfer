#ifndef STATS_H
#define STATS_H

#include "matrix.h"
#include "volume_io.h"
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
#define AVG_VOL          0
#define STD_VOL          1
#define AVD_DOF_VOL      2
#define STD_DOF_VOL      3


#define MAX_EVENTS       10

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
  
SV        *StatReadVolume(char *prefix) ;
fMRI_REG  *StatReadRegistration(char *fname) ;
int       StatWriteVolume(SV *sv, char *prefix) ;
int       StatWriteRegistration(fMRI_REG *reg, char *fname) ;

int       StatFreeRegistration(fMRI_REG **preg) ;
int       StatFree(SV **psv) ;
#if 0
SV        *StatVolumeToTalairach(SV *sv, SV *sv_tal, int resolution) ;
#endif
SV        *StatAllocVolume(SV *sv, int nevents, int width, int height,
                           int nslices, int time_points, int track_dof) ;
SV        *StatAllocTalairachVolume(SV *sv, float fov, float resolution) ;
int       StatAccumulateTalairachVolume(SV *sv_tal, SV *sv) ;

#endif
