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



#ifndef EVSCHUTILS_H
#define EVSCHUTILS_H

#ifdef X
#undef X
#endif

#define EVS_COST_UNKNOWN   0
#define EVS_COST_EFF       1
#define EVS_COST_VRFAVG    2
#define EVS_COST_VRFSTD    3
#define EVS_COST_VRFAVGSTD 4
#define EVS_COST_IDEALXTX  5
#define EVS_COST_EFF_INV   6

typedef struct
{

  int    nevents;
  float *tevent;
  int   *eventid;
  float *weight;

  int    nEvTypes;
  int   *nEvReps;
  float  *EvDur;

  int   nthsearched;
  float cb1err;
  float eff;
  float vrfavg;
  float vrfstd;
  float vrfmin;
  float vrfmax;
  float vrfrange;
  float idealxtxerr;
  float cost;

}
EVENT_SCHEDULE, EVSCH;

EVENT_SCHEDULE *EVSAlloc(int nevents, int allocweight);
int EVSfree(EVENT_SCHEDULE **ppEvSch);
int EVSPrint(FILE *fp, EVENT_SCHEDULE *EvSch);

EVENT_SCHEDULE *EVSreadPar(char *parfile);
int EVSwritePar(char *parfile, EVENT_SCHEDULE *EvSch, char **labels,
                float tPreScan, float tMax);

int EVSmaxId(EVSCH *EvSch);

EVENT_SCHEDULE *EVSsynth(int nEvTypes, int *nPer, float *tPer,
                         float tRes, float tMax, float tPreScan,
                         int nCB1Search, float tNullMin, float tNullMax);

EVENT_SCHEDULE *RandEvSch(int nevents, int ntypes, float dtmin,
                          float dtnullavg, int randweights);
EVENT_SCHEDULE *SynthEvSch(int nEvTypes, int *nPer, float *tPer,
                           float tRes, float tMax, float tPreScan);
MATRIX *EVS2FIRmtx(int EvId, EVSCH *EvSch, float tDelay, float TR,
                   int Ntps, float PSDMin, float PSDMax, float dPSD,
                   MATRIX *X);
MATRIX *EVSfirMtxAll(EVSCH *EvSch, float tDelay, float TR, int Ntps,
                     float PSDMin, float PSDMax, float dPSD);


int EVSsort(EVSCH **EvSchList, int nList);

MATRIX *EVScb1ProbMatrix(EVSCH *EvSch);
MATRIX *EVScb1IdealProbMatrix(EVSCH *EvSch);
MATRIX *EVScb1Matrix(EVSCH *EvSch);
float   EVScb1Error(EVSCH *EvSch);

EVSCH *EVSRandSequence(int nEvTypes, int *nEvReps);
int    EVSRandTiming(EVSCH *EvSch, float *EvDur,
                     float tRes, float tMax, float tPreScan);
EVSCH *EVScb1Optimize(int nEvTypes, int *nEvReps, int nSearch);
const char  *EVScostString(int CostId);
int    EVScostId(const char *CostString);

int EVSdesignMtxStats(MATRIX *Xtask, MATRIX *Xnuis, EVSCH *EvSch,
                      MATRIX *C, MATRIX *W);
float EVScost(EVSCH *EvSch, int CostId, float *params);

int *RandPerm(int N, int *v);
int  RandPermList(int N, int *v);
int RandPermListLimit0(int N, int *v, int lim, int nitersmax);

MATRIX *EVSfirXtXIdeal(int nEvTypes, int *nEvReps, float *EvDur,
                       float TR, int Ntp,
                       float PSDMin, float PSDMax, float dPSD);
int EVSrefractory(EVSCH *sch, double alpha, double T, double dtmin);

#endif //#ifndef EVSCHUTILS_H


