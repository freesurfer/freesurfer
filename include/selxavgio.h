/***************************************************************
  Name:    selxavgio.h
  $Id: selxavgio.h,v 1.3 2004/05/26 16:40:52 greve Exp $
  Author:  Douglas Greve
  Purpose: Routines for handling header files for data created by
  selxavg or selavg (selectively averaged).
 ****************************************************************/

/* data structure for the stuff in the selxavg .dat file */
typedef struct{
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
} SXADAT;

SXADAT * ld_sxadat(char *sxadatfile);
SXADAT * ld_sxadat_from_stem(char *volstem);
int      sv_sxadat(SXADAT *sxadat ,char *sxadatfile);
int      sv_sxadat_by_stem(SXADAT *sxadat ,char *volstem);
int    dump_sxadat(FILE *fp, SXADAT *sxadat ,char *sxadatfile);
int    free_sxadat(SXADAT **sxadat);
int is_sxa_volume(char *volstem);
float *sxa_framepower(SXADAT *sxa, int *nframes);
