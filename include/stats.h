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


#ifndef STATS_H
#define STATS_H

#include "minc.h"

#include "matrix.h"
#include "mri.h"

// These are structures used by mri_segstats
typedef struct
{
  int id;
  char name[1000];
  int nhits;
  float vol;
  int red, green, blue; // 0-255
  float min, max, range, mean, std, snr;
}
STATSUMENTRY;
typedef struct
{
  int nentries;
  char *ColHeaders[100];
  int ncols;
  STATSUMENTRY *entry;
  int IsSurf;
  int UseName;
  int DoIntensity;
  const char *InIntensityName;
  const char *InIntensityUnits;
  int DoSNR; 
} SEGSTAT;

// This is the table structure that can represent the output
// of asegstats2table and aparcstats2table.
typedef struct {
  char *measure;
  int nrows, ncols;
  char **colnames;
  char **rownames;
  double **data;
  char *filename;
  MRI *mri;
}
STAT_TABLE;

typedef struct
{
  float   in_plane_res ;
  float   slice_thickness ;
  float   brightness_scale ;
  MATRIX  *fmri2mri ;
  MATRIX  *mri2fmri ;
  char    name[100] ;            /* subject's name */
}
fMRI_REGISTRATION, fMRI_REG ;

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
}
STAT_VOLUME, SV  ;

/* can't include this before structure, as mrisurf.h includes this file. */
#include "mrisurf.h"
#include "resample.h"

/* This is so applications can specify different xforms */
/* StatReadTransform() will read in this file */
#ifdef _STATS_SRC
const char *stats_talxfm = "talairach.xfm";
int  statnorm_float2int = FLT2INT_ROUND;
int  stats_fixxfm = 0;
#undef _STATS_SRC
#else
extern char *stats_talxfm;
extern int  statnorm_float2int;
extern int  stats_fixxfm;
#endif

SV        *StatReadVolume(const char *prefix);
SV        *StatReadVolume2(const char *prefix);
SV        *StatReadTalairachVolume(const char *prefix,
                                   const char *xform_fname,
                                   const char *subject_name) ;
fMRI_REG  *StatReadRegistration(const char *fname) ;
fMRI_REG  *StatReadTalairachRegistration(const char *fname,
    char *subject_name) ;
int       StatWriteVolume(SV *sv,
                          const char *prefix) ;
int       StatWriteRegistration(fMRI_REG *reg,
                                const char *fname) ;
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
                                     const char *name) ;
int       StatAccumulateTalairachVolume(SV *sv_tal,
                                        SV *sv) ;
int       StatAccumulateSurfaceVolume(SV *sv_tal,
                                      SV *sv,
                                      MRI_SURFACE *mris) ;
int       StatReadTransform(STAT_VOLUME *sv,
                            const char *name) ;
int       StatVolumeExists(const char *prefix) ;

STAT_TABLE *LoadStatTable(const char *statfile);
STAT_TABLE *AllocStatTable(int nrows, int ncols);
STAT_TABLE *InitStatTableFromMRI(MRI* mri_in, const char* tablefile);
int PrintStatTable(FILE *fp, STAT_TABLE *st);
int WriteStatTable(const char *fname, STAT_TABLE *st);

typedef struct
{
  int    label ;
  int    nvoxels ;
  double volume ;
  char   *name ;
  double int_mean ;
  double int_std ;
  double int_min ;
  double int_max ;
} FS_STAT ;
typedef struct 
{
  int     nlabels ;
  FS_STAT *labels ;
} FS_STATS ;

FS_STATS *FSstatsRead(char *fname) ;

int PrintSegStat(FILE *fp, SEGSTAT *segstat);


#endif
