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


/***************************************************************
  Name:    selxavgio.h
  Author:  Douglas Greve
  Purpose: Routines for handling header files for data created by
  selxavg or selavg (selectively averaged).
 ****************************************************************/

/* data structure for the stuff in the selxavg .dat file */
typedef struct
{
  int version;
  float TR;
  float TER;
  float TimeWindow;
  float TPreStim;
  int   Nc;    /* Total number of conditions, incl fix*/
  int   Nnnc;  /* Total number of conditions, exlc fix*/
  int   Nh;    /* Number of estimtes per condition = TW/TER*/
  float DOF;
  int   *npercond;
  int   nruns;
  int   ntp;
  int   nrows;
  int   ncols;
  int   nskip;
  int   DTOrder;
  float RescaleFactor;
  float HanningRadius;
  int   BrainAirSeg;
  int   GammaFit;
  float gfDelta[10];
  float gfTau[10];
  int   NullCondId;
  float *SumXtX;
  float *hCovMtx;
  int   nNoiseAC;
  int   *CondIdMap;
}
SXADAT;

SXADAT * ld_sxadat(const char *sxadatfile);
SXADAT * ld_sxadat_from_stem(const char *volstem);
int      sv_sxadat(SXADAT *sxadat ,const char *sxadatfile);
int      sv_sxadat_by_stem(SXADAT *sxadat ,const char *volstem);
int    dump_sxadat(FILE *fp, SXADAT *sxadat ,const char *sxadatfile);
int    free_sxadat(SXADAT **sxadat);
int is_sxa_volume(const char *volstem);
float *sxa_framepower(SXADAT *sxa, int *nframes);
