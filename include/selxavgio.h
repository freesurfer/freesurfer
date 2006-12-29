/**
 * @file  selxavgio.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2006/12/29 02:09:00 $
 *    $Revision: 1.4 $
 *
 * Copyright (C) 2002-2007,
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


/***************************************************************
  Name:    selxavgio.h
  $Id: selxavgio.h,v 1.4 2006/12/29 02:09:00 nicks Exp $
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

SXADAT * ld_sxadat(char *sxadatfile);
SXADAT * ld_sxadat_from_stem(char *volstem);
int      sv_sxadat(SXADAT *sxadat ,char *sxadatfile);
int      sv_sxadat_by_stem(SXADAT *sxadat ,char *volstem);
int    dump_sxadat(FILE *fp, SXADAT *sxadat ,char *sxadatfile);
int    free_sxadat(SXADAT **sxadat);
int is_sxa_volume(char *volstem);
float *sxa_framepower(SXADAT *sxa, int *nframes);
