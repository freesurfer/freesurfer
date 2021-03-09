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

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "diag.h"
#include "error.h"

#include "matfile.h"
#include "matrix.h"
#include "proto.h"
#include "utils.h"

#include "evschutils.h"

#ifdef const
#undef const
#endif
static int EVScompare(const void *evsch1, const void *evsch2);

/*-------------------------------------------------------------*/
EVENT_SCHEDULE *EVSAlloc(int nevents, int allocweight)
{
  EVENT_SCHEDULE *EvSch;

  EvSch = (EVENT_SCHEDULE *)calloc(sizeof(EVENT_SCHEDULE), 1);
  EvSch->nevents = nevents;

  EvSch->tevent = (float *)calloc(sizeof(float), nevents);
  EvSch->eventid = (int *)calloc(sizeof(int), nevents);
  if (allocweight) EvSch->weight = (float *)calloc(sizeof(float), nevents);

  EvSch->nEvTypes = 0;
  EvSch->nEvReps = NULL;
  EvSch->EvDur = NULL;

  return (EvSch);
}
/*-------------------------------------------------------------*/
int EVSfree(EVENT_SCHEDULE **ppEvSch)
{
  EVENT_SCHEDULE *pEvSch;

  pEvSch = *ppEvSch;

  if (pEvSch->tevent != NULL) free(pEvSch->tevent);
  if (pEvSch->eventid != NULL) free(pEvSch->eventid);
  if (pEvSch->weight != NULL) free(pEvSch->weight);
  if (pEvSch->nEvReps != NULL) free(pEvSch->nEvReps);
  if (pEvSch->EvDur != NULL) free(pEvSch->EvDur);

  free(*ppEvSch);

  *ppEvSch = NULL;

  return (0);
}

/*-------------------------------------------------------------
  EVS2FIRmtx() - convert an event schedule for a given event id into
  an FIR design matrix.
  -------------------------------------------------------------*/
MATRIX *EVS2FIRmtx(
    int EvId, EVSCH *EvSch, float tDelay, float TR, int Ntps, float PSDMin, float PSDMax, float dPSD, MATRIX *X)
{
  float tMax, tmp, PSDWindow, tPSD, PSD;
  int RSR, Npsds, nthPSD, n, rA, rB;

  /* Compute number of PSDs in the window */
  PSDWindow = PSDMax - PSDMin;
  tmp = rint(PSDWindow / dPSD) - PSDWindow / dPSD;
  if (tmp > .0001) {
    printf("ERROR: EVS2FIRmtx: PSDWindow (%g) is not an integer multiple of dPSD (%g)\n", PSDWindow, dPSD);
    return (NULL);
  }
  Npsds = rint(PSDWindow / dPSD);

  /* Compute resampling rate */
  tmp = rint(TR / dPSD) - TR / dPSD;
  if (tmp > .0001) {
    printf(
        "ERROR: EVS2FIRmtx: TR (%g) is not an integer multiple "
        " of dPSD (%g)\n",
        TR,
        dPSD);
    return (NULL);
  }
  RSR = rint(TR / dPSD);

  /* Compute the time of the last row in X */
  tMax = TR * (Ntps - 1);

  /* Create X with all zeros */
  if (X == NULL)
    X = MatrixAlloc(Ntps, Npsds, MATRIX_REAL);
  else {
    if (X->rows != Ntps || X->rows != Npsds) {
      printf("ERROR: EVS2FIRmtx: dimensions of X mismatch\n");
      return (NULL);
    }
    X = MatrixZero(Ntps, Npsds, X);
  }

  /* Fill-in the non-zero entries of X for each event */
  /* nthPSD will be the column in X */
  for (nthPSD = 0; nthPSD < Npsds; nthPSD++) {
    PSD = nthPSD * dPSD + PSDMin;

    for (n = 0; n < EvSch->nevents; n++) {
      if (EvSch->eventid[n] != EvId) continue;

      tPSD = EvSch->tevent[n] + tDelay + PSD;
      if (tPSD < 0.0) continue;
      if (tPSD > tMax) break;

      /* Could compute the time of the closest row, then skip if
fabs(tRow-tPSD) > dPSD/2 ????? This would be easily
   extended to the cases where the data were not aquired
uniformly in time. */

      /* rA would be the row of X if the rows of X were separated by dPSD */
      rA = (int)rint(tPSD / dPSD);

      /* If rA does not fall into a row of X, skip */
      if ((rA % RSR) != 0) continue;

      /* rB is the row of X */
      rB = rA / RSR;
      if (EvSch->weight != NULL)
        X->rptr[rB + 1][nthPSD + 1] = EvSch->weight[n];
      else
        X->rptr[rB + 1][nthPSD + 1] = 1;
    }
  }

  return (X);
}
/*-----------------------------------------------------------------------
  EVSfirMtxAll() - convert an event schedule to a design matrix for
  all events types. The matrices for the individual event types are
  horizontally concatenated.
  ---------------------------------------------------------------------*/
MATRIX *EVSfirMtxAll(EVSCH *EvSch, float tDelay, float TR, int Ntps, float PSDMin, float PSDMax, float dPSD)
{
  int ev;
  MATRIX *Xevfir = NULL, *Xtmp = NULL, *Xfir = NULL;

  for (ev = 1; ev <= EvSch->nEvTypes; ev++) {
    Xevfir = EVS2FIRmtx(ev, EvSch, 0, TR, Ntps, PSDMin, PSDMax, dPSD, Xevfir);
    if (Xevfir == NULL) {
      printf("ERROR: EVSfirMtxAll: Cannot compute Xevfir\n");
      return (NULL);
    }
    Xtmp = MatrixHorCat(Xfir, Xevfir, NULL);
    if (Xfir != NULL) MatrixFree(&Xfir);
    Xfir = Xtmp;

    MatrixFree(&Xevfir);
    Xevfir = NULL;
  }

  return (Xfir);
}
/*-------------------------------------------------------------*/
int EVSPrint(FILE *fp, EVENT_SCHEDULE *EvSch)
{
  int n;
  for (n = 0; n < EvSch->nevents; n++) {
    // fprintf(fp,"%3d %8.4f  %3d ",n,EvSch->tevent[n],EvSch->eventid[n]);
    fprintf(fp, "%8.4f  %3d ", EvSch->tevent[n], EvSch->eventid[n]);
    fprintf(fp, "%6.2f ", EvSch->EvDur[EvSch->eventid[n] - 1]);
    if (EvSch->weight != 0) fprintf(fp, "  %8.4f", EvSch->weight[n]);
    fprintf(fp, "\n");
  }
  return (0);
}
/*-----------------------------------------------------------------
  EVSwritePar() - write the event schedule out in paradigm file
  format. Places where the fabs(time) is tested for < some small
  number is so that when time is zero, it is not printed out as
  '-0.00'. This is PURELY for cosmetic reasons since it is usually
  the first thing that users see.
  -----------------------------------------------------------------*/
int EVSwritePar(char *parfile, EVENT_SCHEDULE *EvSch, char **labels, float tPreScan, float tMax)
{
  FILE *fp;
  float NullDur = 0.0, tNull = 0.0, t = 0.0;
  int n, id;

  fp = fopen(parfile, "w");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", parfile);
    return (1);
  }
  if (EvSch->tevent[0] > -tPreScan) {
    NullDur = EvSch->tevent[0] - (-tPreScan);
    if (fabs(tPreScan) < .0000000001)
      t = 0.0; /* purely for cosmetic reasons */
    else
      t = -tPreScan;
    fprintf(fp, "%8.4f  %3d   %6.3f   %5.4f    %10s \n", t, 0, NullDur, 1.0, "NULL");
  }

  for (n = 0; n < EvSch->nevents; n++) {
    id = EvSch->eventid[n];
    if (id == 0) continue;

    /* Print the time, id, duration, label of current event */
    if (fabs(EvSch->tevent[n]) < .0000000001)
      t = 0.0;
    else
      t = EvSch->tevent[n];
    fprintf(fp, "%8.4f  %3d   %6.3f   ", t, EvSch->eventid[n], EvSch->EvDur[id - 1]);
    if (EvSch->weight)
      fprintf(fp, "%5.4f    ", EvSch->weight[n]);
    else
      fprintf(fp, "%5.4f    ", 1.0);
    fprintf(fp, "%10s\n", labels[id - 1]);

    /* Compute duration of null event following this event (may be zero)*/
    if (n == EvSch->nevents - 1)
      NullDur = tMax - (EvSch->tevent[n] + EvSch->EvDur[id - 1]);
    else
      NullDur = EvSch->tevent[n + 1] - (EvSch->tevent[n] + EvSch->EvDur[id - 1]);

    /* If Null time is not zero, print time, id, duration, label for NULL */
    if (NullDur > .0001) {
      tNull = EvSch->tevent[n] + EvSch->EvDur[id - 1];
      if (fabs(tNull) < .0000000001)
        t = 0.0;
      else
        t = tNull;
      fprintf(fp, "%8.4f  %3d   %6.3f   %5.4f    %10s \n", t, 0, NullDur, 1.0, "NULL");
    }
  }
  fclose(fp);

  return (0);
}
/*-----------------------------------------------------------------
  EVSsynth() - synthesize an event schedule with the given number
  of event types, the given number of repetitions per event type
  and duration of each event type spread out over the given interval
  from -tPreScan to +tMax. Event onset times are constrained to
  fall on integer multiples of tRes. First-order counter-balancing
  is pre-optimized by searching over nCB1Search sequences.
-----------------------------------------------------------------*/
EVENT_SCHEDULE *EVSsynth(int nEvTypes,
                         int *nPer,
                         float *tPer,
                         float tRes,
                         float tMax,
                         float tPreScan,
                         int nCB1Search,
                         float tNullMin,
                         float tNullMax)
{
  int id, m, n, nevents, nSlotsNull, nSlotsTot, *EvSeq, nNullMax;
  float tStimTot, t, tScanTot, tNullTot;
  EVENT_SCHEDULE *EvSch;

  /* Compute the total amount of stimulation time */
  tStimTot = 0.0;
  nevents = 0;
  for (n = 0; n < nEvTypes; n++) {
    nevents += nPer[n];
    tStimTot += nPer[n] * (tPer[n] + tNullMin);
    // printf("%2d %d %g %g %g\n",n,nPer[n],tPer[n],tNullMin,tStimTot);
  }

  /* Compute the total amount of scan time and compare to stim time */
  tScanTot = tMax + tPreScan;
  if (tStimTot > tScanTot) {
    printf("ERROR: stimulation time %g (including tNullMin %g) exceeds scan time %g\n",
           tStimTot,
           nevents * tNullMin,
           tScanTot);
    return (NULL);
  }

  /* Synthesize the sequence */
  if (nCB1Search > 1)
    EvSch = EVScb1Optimize(nEvTypes, nPer, nCB1Search);
  else
    EvSch = EVScb1Optimize(nEvTypes, nPer, 1);

  /* The code  below is for synthesizing the timing */
  /* Compute the total amount of null time */
  tNullTot = tScanTot - tStimTot;

  /* Comute number of slots allocated for null */
  nSlotsNull = (int)floor(tNullTot / tRes);

  /* Compute number of slots at given temp res, include slots for null */
  nSlotsTot = nevents + nSlotsNull;

  /* Compute maximum number of back-to-back null slots */
  if (tNullMax > 0)
    nNullMax = (int)(floor((tNullMax - tNullMin) / tRes));
  else
    nNullMax = nSlotsTot;

  if (Gdiag_no > 0) {
    printf("EVSsynth(): tNullMin=%g \n", tNullMin);
    printf("EVSsynth(): tScanTot=%g, tStimTot=%g, tNullTot=%g \n", tScanTot, tStimTot, tNullTot);
    printf("EVSsynth(): nNullMax=%d, nSlotsTot=%d, nevents=%d\n", nNullMax, nSlotsTot, nevents);
  }

  /* Create a non-random sequence of 0s and 1s, 1 = Non-Null */
  EvSeq = (int *)calloc(sizeof(int), nSlotsTot);
  for (n = 0; n < nevents; n++) EvSeq[n] = 1;
  /* Randomize the sequence of nulls and non-nulls*/
  m = RandPermListLimit0(nSlotsTot, EvSeq, nNullMax, 1000000);
  if (m < 0) {
    printf("ERROR: could not enforce tNullMax=%g (ntries=1000000)\n", tNullMax);
    printf("You will need to reduce the number of time points\n");
    printf("or increase the number of presentations.\n");
    return (NULL);
  }
  // Assure that first event is non-null. Swap with first non-null
  if (EvSeq[0] != 1) {
    for (n = 0; n < nSlotsTot; n++) {
      if (EvSeq[n]) {
        EvSeq[0] = 1;
        EvSeq[n] = 0;
        break;
      }
    }
  }

  /* Compute the the timing */
  m = 0;
  t = -tPreScan;
  for (n = 0; n < nSlotsTot; n++) {
    // printf("%6.2f %d\n",t,EvSeq[n]);
    if (EvSeq[n] != 0) {
      id = EvSch->eventid[m];
      EvSch->tevent[m] = t;
      t += tPer[id - 1];
      t += tNullMin;
      m++;
    }
    else
      t += tRes;
  }

  /* Alloc and fill the number the event-type duration */
  EvSch->EvDur = (float *)calloc(sizeof(float), nEvTypes);
  for (n = 0; n < nEvTypes; n++) EvSch->EvDur[n] = tPer[n];

  free(EvSeq);
  return (EvSch);
}
/*-----------------------------------------------------------
  EVScb1Optimize() - search for an event sequence which will
  minimize the first-order counter-balancing error.
-----------------------------------------------------------*/
EVSCH *EVScb1Optimize(int nEvTypes, int *nEvReps, int nSearch)
{
  EVSCH *EvSch, *EvSchBest = NULL;
  int n;

  /* Loop over the number of search iterations */
  for (n = 0; n < nSearch; n++) {
    /* Get a random sequence of events */
    EvSch = EVSRandSequence(nEvTypes, nEvReps);

    /* Compute the CB1 cost */
    EVScb1Error(EvSch);

    /* Keep the best */
    if (n == 0)
      EvSchBest = EvSch;
    else {
      if (EvSchBest->cb1err > EvSch->cb1err) {
        EVSfree(&EvSchBest);
        EvSchBest = EvSch;
      }
      else
        EVSfree(&EvSch);
    }
  }

  return (EvSchBest);
}
/*-----------------------------------------------------------
  EVSRandSequence() - create a random sequence of events for the
  given number of event types and number of repetitions per
  event.
  -----------------------------------------------------------*/
EVSCH *EVSRandSequence(int nEvTypes, int *nEvReps)
{
  EVSCH *EvSch;
  int nevents, m, n, nthev;

  nevents = 0;
  for (n = 0; n < nEvTypes; n++) nevents += nEvReps[n];

  EvSch = EVSAlloc(nevents, 0);
  EvSch->nEvTypes = nEvTypes;

  /* Alloc and fill the number of reps (but not the duration) */
  EvSch->nEvReps = (int *)calloc(sizeof(int), nEvTypes);
  for (n = 0; n < nEvTypes; n++) EvSch->nEvReps[n] = nEvReps[n];

  /* Fill the event list with the proper number of each type */
  nthev = 0;
  for (n = 0; n < nEvTypes; n++) {
    for (m = 0; m < nEvReps[n]; m++) {
      EvSch->eventid[nthev] = n + 1;
      nthev++;
    }
  }

  /* Randomly permute the event list to randomize */
  RandPermList(nevents, EvSch->eventid);

  return (EvSch);
}
/*------------------------------------------
  EVSmaxId() - return the maximum id.
  ------------------------------------------*/
int EVSmaxId(EVSCH *EvSch)
{
  int n, id, condmax;

  condmax = -1;
  for (n = 0; n < EvSch->nevents; n++) {
    id = EvSch->eventid[n];
    if (id > condmax) condmax = id;
  }

  return (condmax);
}
/*-------------------------------------------------------------
  EVSreadPar() - read an event schedule stored in paradigm
  file format.
-------------------------------------------------------------*/
EVENT_SCHEDULE *EVSreadPar(char *parfile)
{
  FILE *fp;
  EVENT_SCHEDULE *EvSch;
  int nevents, id, ev;
  char tmpstring[2001];
  float tev;

  fp = fopen(parfile, "r");
  if (fp == NULL) {
    printf("ERROR: cannot open %s\n", parfile);
    return (NULL);
  }

  /* Count the number of non-null events */
  nevents = 0;
  while (fgets(tmpstring, 2000, fp) != NULL) {
    sscanf(tmpstring, "%f %d", &tev, &id);
    if (id != 0) nevents++;
  }
  fclose(fp);
  // printf("INFO: found %d events\n",nevents);

  EvSch = EVSAlloc(nevents, 0);

  /* Close/Open to rewind */
  fp = fopen(parfile, "r");
  nevents = 0;
  while (fgets(tmpstring, 2000, fp) != NULL) {
    sscanf(tmpstring, "%f %d", &tev, &id);
    if (id != 0) {
      EvSch->tevent[nevents] = tev;
      EvSch->eventid[nevents] = id;
      nevents++;
    }
  }
  fclose(fp);

  EvSch->nEvTypes = EVSmaxId(EvSch); /* holes? */

  /* Fill in the number of repetitions */
  EvSch->nEvReps = (int *)calloc(sizeof(int), EvSch->nEvTypes);
  for (ev = 0; ev < nevents; ev++) {
    id = EvSch->eventid[ev];
    EvSch->nEvReps[id - 1]++;
  }

  return (EvSch);
}
/*------------------------------------------------------------
  RandPerm() - returns a list of randomly permuted integers
  between 0 and N-1. Should be the same as matlab's.
  ------------------------------------------------------------*/
int *RandPerm(int N, int *v)
{
  int tmp, n, n2;

  if (v == NULL) v = (int *)calloc(sizeof(int), N);
  for (n = 0; n < N; n++) v[n] = n;

  for (n = 0; n < N; n++) {
    n2 = (int)floor(drand48() * N);
    tmp = v[n];
    v[n] = v[n2];
    v[n2] = tmp;
  }

  return (v);
}
/*------------------------------------------------------------------
 RandPermListLimit0() - randomly permute members of a list but limit
 the maximum number of times member0 can appear back-to-back. For
 example, if lim=5, then v[n]...v[n+5] could not be all be 0. This is
 done by randomly permuting the sequence until a legal one is found.
 Returns  0 upon successfully finding a legal sequence. Returns -1
 if the maximum number of iterations is exceeded.
 -----------------------------------------------------------------*/
int RandPermListLimit0(int N, int *v, int lim, int nitersmax)
{
  int niters, runlenmax, runlen, n;

  niters = 0;
  while (1) {
    RandPermList(N, v);  // permute the list
    // Count the max run length of items whose val is 0
    runlenmax = 0;
    runlen = 0;
    for (n = 0; n < N; n++) {
      if (v[n] == 0)
        runlen++;
      else {
        if (runlenmax < runlen) runlenmax = runlen;
        runlen = 0;
      }
    }
    // Termination conditions. Require first stim to be non-null,
    // I can't remember why.
    if (runlenmax <= lim && v[0] != 0) return (0);
    if (niters > nitersmax) return (-1);
    niters++;
  }

  return (0);  // should never get here
}

/*---------------------------------------------------------------
  RandPermList() - randomly permutes members of the given list.
  ------------------------------------------------------------*/
int RandPermList(int N, int *v)
{
  int *p, *vp, n;

  if (v == NULL) {
    printf("ERROR: RandPermVect: vector is null\n");
    return (1);
  }

  p = RandPerm(N, NULL);
  vp = (int *)calloc(sizeof(int), N);
  for (n = 0; n < N; n++) vp[n] = v[n];
  for (n = 0; n < N; n++) v[n] = vp[p[n]];

  free(p);
  free(vp);

  return (0);
}
/*----------------------------------------------------------
  EVSsort() - sort a list of events schedules based on the
  value in the cost element. Sorts so that the maximum
  cost is first (is descending order).
  ---------------------------------------------------------*/
int EVSsort(EVSCH **EvSchList, int nList)
{
  qsort(EvSchList, nList, sizeof(EVSCH **), EVScompare);
  return (0);
}
/*--------------------------------------------------------------
  EVScompare() - compare function suitable for the qsort() routine.
  Compares based on he cost element so that the maximum cost is first
  (is descending order).
  --------------------------------------------------------------*/
static int EVScompare(const void *evsch1, const void *evsch2)
{
  EVSCH *EvSch1, *EvSch2;

  EvSch1 = *((EVSCH **)evsch1);
  EvSch2 = *((EVSCH **)evsch2);

  // printf("cost: %g %g\n",EvSch1->cost,EvSch2->cost);

  /* Sort so that the highest cost is first */
  if (EvSch1->cost > EvSch2->cost) return (-1);
  if (EvSch1->cost < EvSch2->cost) return (+1);
  return (0);
}
/*--------------------------------------------------------------
  EVScb1Matrix() - computes the first-order counter-balance
  matrix. The number in element i,j is the number of times that
  condition i was followed by condition j. This does not include
  condition 0.
  --------------------------------------------------------------*/
MATRIX *EVScb1Matrix(EVSCH *EvSch)
{
  MATRIX *Ncb1;
  int nthev, i, j;

  Ncb1 = MatrixZero(EvSch->nEvTypes, EvSch->nEvTypes, NULL);

  for (nthev = 0; nthev < EvSch->nevents - 1; nthev++) {
    i = EvSch->eventid[nthev];
    j = EvSch->eventid[nthev + 1];
    if (i == 0 || j == 0) continue;
    Ncb1->rptr[i][j]++;
  }

  return (Ncb1);
}
/*--------------------------------------------------------------
  EVScb1ProbMatrix() - computes an estimate of the first-order
  counter-balance PROBABILITY matrix. The value in element i,j is the
  probability (or rate) that condition i was followed by condition
  j. This does not include condition 0. In theory, the rows of
  Pcb1 should sum to 1.0, but this will not be the case for
  which ever condition the sequence ends on.
  --------------------------------------------------------------*/
MATRIX *EVScb1ProbMatrix(EVSCH *EvSch)
{
  MATRIX *Ncb1, *Pcb1 = NULL;
  int i, j;

  Ncb1 = EVScb1Matrix(EvSch);
  Pcb1 = MatrixZero(EvSch->nEvTypes, EvSch->nEvTypes, NULL);

  for (i = 1; i <= EvSch->nEvTypes; i++) {
    for (j = 1; j <= EvSch->nEvTypes; j++) {
      Pcb1->rptr[i][j] = Ncb1->rptr[i][j] / EvSch->nEvReps[i - 1];
      /* Note: i-1 needed in nEvReps because its normal C array */
    }
  }
  MatrixFree(&Ncb1);

  return (Pcb1);
}
/*--------------------------------------------------------------
  EVScb1IdealProbMatrix() - computes the ideal first-order
  counter-balance PROBABILITY matrix. The value in element i,j is the
  ideal probability that condition i should be followed by condition
  j, which is just equal to the probability of condition j. This does
  not include condition 0.
  --------------------------------------------------------------*/
MATRIX *EVScb1IdealProbMatrix(EVSCH *EvSch)
{
  MATRIX *IdealPcb1;
  int i, j;

  IdealPcb1 = MatrixZero(EvSch->nEvTypes, EvSch->nEvTypes, NULL);

  for (i = 1; i <= EvSch->nEvTypes; i++) {
    for (j = 1; j <= EvSch->nEvTypes; j++) {
      // printf("i=%d, j=%d, IP = %g, nEVj=%d, nev=%d\n",i,j,
      //     IdealPcb1->rptr[i][j],EvSch->nEvReps[j-1],EvSch->nevents);
      IdealPcb1->rptr[i][j] = (float)EvSch->nEvReps[j - 1] / EvSch->nevents;
      /* Note: j-1 needed in nEvReps because its normal C array */
    }
  }

  return (IdealPcb1);
}
/*------------------------------------------------------------
  EVScb1Error() - computes a measure of how far a given sequence is
  from having ideal first-order counter-balancing. It is computed as
  the relative difference between the ideal counter-balance
  probability matrix and the actual matrix, averaged over all the
  elements in the matrix.
  ------------------------------------------------------------*/
float EVScb1Error(EVSCH *EvSch)
{
  MATRIX *IdealPcb1, *Pcb1;
  int i, j;
  float cb1err;

  IdealPcb1 = EVScb1IdealProbMatrix(EvSch);
  Pcb1 = EVScb1ProbMatrix(EvSch);

  cb1err = 0;
  for (i = 1; i <= EvSch->nEvTypes; i++) {
    for (j = 1; j <= EvSch->nEvTypes; j++) {
      cb1err += fabs(IdealPcb1->rptr[i][j] - Pcb1->rptr[i][j]) / IdealPcb1->rptr[i][j];
    }
  }
  cb1err /= (EvSch->nEvTypes * EvSch->nEvTypes);

  MatrixFree(&IdealPcb1);
  MatrixFree(&Pcb1);

  EvSch->cb1err = cb1err;

  return (cb1err);
}
/*---------------------------------------------------------------
  EVScostId() - converts from a string indicating a cost function
  to a numeric id that can be used in a case statement.
  ----------------------------------------------------------------*/
int EVScostId(const char *CostString)
{
  if (!strcmp("eff", CostString)) return (EVS_COST_EFF);
  if (!strcmp("effinv", CostString)) return (EVS_COST_EFF_INV);
  if (!strcmp("vrfavg", CostString)) return (EVS_COST_VRFAVG);
  if (!strcmp("vrfstd", CostString)) return (EVS_COST_VRFSTD);
  if (!strcmp("vrfavgstd", CostString)) return (EVS_COST_VRFAVGSTD);
  if (!strcmp("idealxtx", CostString)) return (EVS_COST_IDEALXTX);

  return (EVS_COST_UNKNOWN);
}
/*---------------------------------------------------------------
  EVScostString() - converts from a numeric id (used in a case
  statement) indicating a cost function to a human-readable
  string.
  ----------------------------------------------------------------*/
const char *EVScostString(int CostId)
{
  switch (CostId) {
    case EVS_COST_UNKNOWN:
      return ("unknown");
      break;
    case EVS_COST_EFF:
      return ("eff");
      break;
    case EVS_COST_EFF_INV:
      return ("effinv");
      break;
    case EVS_COST_VRFAVG:
      return ("vrfavg");
      break;
    case EVS_COST_VRFSTD:
      return ("vrfstd");
      break;
    case EVS_COST_VRFAVGSTD:
      return ("vrfavgstd");
      break;
    case EVS_COST_IDEALXTX:
      return ("idealxtx");
      break;
  }

  return (NULL);
}
/*------------------------------------------------------------------
  EVScost() - compute the cost of the event schedule based on
  the cost function (specified with CostId).
  ------------------------------------------------------------------*/
float EVScost(EVSCH *EvSch, int CostId, float *params)
{
  switch (CostId) {
    case EVS_COST_EFF:
      /* efficiency */
      EvSch->cost = EvSch->eff;
      break;
    case EVS_COST_EFF_INV:
      /* 1/efficiency -- get worst schedule*/
      EvSch->cost = 1. / (EvSch->eff + .000000000001);
      break;
    case EVS_COST_VRFAVG:
      /* variance reduction factor averaged across task estimates */
      EvSch->cost = EvSch->vrfavg;
      break;
    case EVS_COST_VRFSTD:
      /* variance reduction factor stddev across task estimates */
      EvSch->cost = -EvSch->vrfstd;
      break;
    case EVS_COST_VRFAVGSTD:
      /* Weighted combo of VRF average and std */
      EvSch->cost = EvSch->vrfavg - params[0] * EvSch->vrfstd;
      break;
    case EVS_COST_IDEALXTX:
      /* Weighted combo of VRF average and std */
      EvSch->cost = -EvSch->idealxtxerr;
      break;
    default:
      printf("ERROR: CostId %d unrecoginized\n", CostId);
      return (-1000000.0);
      break;
  }

  return (EvSch->cost);
}
/*--------------------------------------------------------------------
  EVSdesignMtxStats() - computes statistics about design relevant for
  optimization. stats should have at least 6 elements. Returns 1 is
  the design is singular, 0 otherwise.  ERROR: This needs to be
  modified to compute the VRF differently when there are different
  numbers of stimuli in each event type.x
  -------------------------------------------------------------------*/
int EVSdesignMtxStats(MATRIX *Xtask, MATRIX *Xnuis, EVSCH *EvSch, MATRIX *C, MATRIX *W)
{
  MATRIX *X = NULL, *Xt = NULL, *XtX = NULL;
  MATRIX *iXtX = NULL, *VRF = NULL, *Ct = NULL, *CiXtX = NULL, *CiXtXCt = NULL;
  int r, m, nTaskAvgs, nAvgs, Cfree, J;
  // int nNuisAvgs;
  float diagsum;
  double dtmp = 0;
  double dtmp1 = 0;
  double dtmp2 = 0;

  X = MatrixHorCat(Xtask, Xnuis, NULL);
  nTaskAvgs = Xtask->cols;
  // if (Xnuis != NULL)
  //   nNuisAvgs = Xnuis->cols;
  // else
  //   nNuisAvgs = 0;
  nAvgs = X->cols;

  if (W != NULL) X = MatrixMultiply(W, X, NULL);

  Xt = MatrixTranspose(X, Xt);
  XtX = MatrixMultiply(Xt, X, XtX);

  /* Compute the Inverse */
  iXtX = MatrixInverse(XtX, NULL);

  Cfree = 0;
  if (C == NULL) {
    C = MatrixZero(nTaskAvgs, nAvgs, NULL);
    for (m = 0; m < nTaskAvgs; m++) C->rptr[m + 1][m + 1] = 1;
    Cfree = 1;
  }
  J = C->rows;
  Ct = MatrixTranspose(C, NULL);

  /* Make sure that it was actually inverted */
  if (iXtX != NULL) {
    r = 0;

    CiXtX = MatrixMultiply(C, iXtX, NULL);
    CiXtXCt = MatrixMultiply(CiXtX, Ct, NULL);

    VRF = MatrixAlloc(J, 1, MATRIX_REAL);

    diagsum = 0.0;
    for (m = 0; m < J; m++) { /* exctract diag */
      diagsum += CiXtXCt->rptr[m + 1][m + 1];
      VRF->rptr[m + 1][1] = 1.0 / (CiXtXCt->rptr[m + 1][m + 1]);
    }
    EvSch->eff = 1.0 / diagsum;
    if (J > 1)
      EvSch->vrfstd = VectorStdDev(VRF, &dtmp);
    else
      EvSch->vrfstd = 0.0;
    EvSch->vrfavg = dtmp;
    EvSch->vrfrange = VectorRange(VRF, &dtmp1, &dtmp2);
    EvSch->vrfmin = dtmp1;
    EvSch->vrfmax = dtmp2;

    MatrixFree(&iXtX);
    MatrixFree(&CiXtX);
    MatrixFree(&CiXtXCt);
    MatrixFree(&VRF);
  }
  else
    r = 1;

  MatrixFree(&X);
  MatrixFree(&Xt);
  MatrixFree(&XtX);
  if (Cfree) MatrixFree(&C);
  MatrixFree(&Ct);

  return (r);
}
/*-------------------------------------------------------------
  EVSfirXtXIdeal() - computes the ideal XtX matrix for the FIR
  signal model. The XtX matrix can be broken down into nEvTypes-
  by-nEvTypes blocks, each nPSD-by-nPSD. The value at the nth row
  and mth col within the block at i,j within the block matrix
  represents the number of times one would expect to see Event
  Type i followed (n-m) dPSDs by Event Type j. If m > n, then
  j actually preceeds i.

  The block at row=i and col=j within the block matrix represents
  EVT j (col) being followed by EVT i (row).

  Within block i,j, the nth col corresponds to condition j being
  shifted by n, and the mth row corresponds to condition i being
  shifted by m, with respect to the start of the scan. In other
  words, element m,n represents condition i appearing (n-m) dPSDs
  after condition j.

  For example, for two conditions in block i=2,j=1, the number at the
  first row (m=1) and third column (n=3) represents the expected
  number of times that EVT1 (j=1) should be followed by EVT2 (i=2) by
  3-1=2 dPSDs.
  -------------------------------------------------------------*/
MATRIX *EVSfirXtXIdeal(
    int nEvTypes, int *nEvReps, float *EvDur, float TR, int Ntp, float PSDMin, float PSDMax, float dPSD)
{
  MATRIX *XtX = NULL;
  int Npsd, Navgs, npsd1, npsd2, nevt1, nevt2;
  int n, nslots, nshift;
  int r, c;
  int *EvDur_dPSD;
  float vx;

  Npsd = (int)(rint((PSDMax - PSDMin) / dPSD));
  Navgs = nEvTypes * Npsd;

  EvDur_dPSD = (int *)calloc(sizeof(int), nEvTypes);
  for (n = 0; n < nEvTypes; n++) EvDur_dPSD[n] = (int)(rint(EvDur[n] / dPSD));

  nslots = (int)(rint(TR * Ntp / dPSD));

  XtX = MatrixAlloc(Navgs, Navgs, MATRIX_REAL);

  /* Loop through the Event Type Block matrix */
  for (nevt1 = 0; nevt1 < nEvTypes; nevt1++) {   /* block rows (j)*/
    for (nevt2 = 0; nevt2 < nEvTypes; nevt2++) { /* block cols (i)*/

      for (npsd1 = 0; npsd1 < Npsd; npsd1++) {   /* rows (m)*/
        for (npsd2 = 0; npsd2 < Npsd; npsd2++) { /* cols (n)*/

          /* nshift is the number of dPSDs that EVT1 follows EVT2*/
          nshift = npsd2 - npsd1; /* n-m */
          r = nevt1 * Npsd + npsd1 + 1;
          c = nevt2 * Npsd + npsd2 + 1;

          if (nevt1 == nevt2) {
            /* diagonal block */
            if (nshift == 0)
              vx = nEvReps[nevt1];
            else if (abs(nshift) < EvDur_dPSD[nevt1])
              vx = 0.0;
            else
              vx = nEvReps[nevt1] * nEvReps[nevt1] / nslots;
          }
          else {
            /* off-diagonal block */
            if (nshift == 0 || (nshift > 0 && nshift < EvDur_dPSD[nevt2]) ||
                (nshift < 0 && -nshift < EvDur_dPSD[nevt1]))
              vx = 0.0;
            else
              vx = (float)nEvReps[nevt1] * nEvReps[nevt2] / nslots;
          }
          XtX->rptr[r][c] = vx;
        }
      }
    }
  }

  free(EvDur_dPSD);
  return (XtX);
}
/*--------------------------------------------------------
  EVSrefractory() - compute the weight for each event based
  on the amount of time since the end of the last event.
  This penalizes based on a refractory model.
  --------------------------------------------------------*/
int EVSrefractory(EVSCH *sch, double alpha, double T, double dtmin)
{
  int n, idprev;
  double tonprev, durprev, tonev, dt = 0.0, toffprev;

  if (sch->weight == NULL) sch->weight = (float *)calloc(sizeof(float), sch->nevents);

  sch->weight[0] = 1;
  for (n = 1; n < sch->nevents; n++) {
    idprev = sch->eventid[n - 1];
    tonprev = sch->tevent[n - 1];
    durprev = sch->EvDur[idprev - 1];
    toffprev = tonprev + durprev;
    tonev = sch->tevent[n];
    dt = tonev - toffprev + dtmin;
    if (dt < 0) dt = 0;
    sch->weight[n] = 1 - alpha * exp(-dt / T);
    // printf("%2d %d %g %g %g  %g %g %g\n",n,idprev,
    //   tonprev,durprev,toffprev,tonev,dt,sch->weight[n]);
  }

  // EVSPrint(stdout,sch);
  // exit(1);

  return (0);
}

#if 0
/*-----------------------------------------------------------*/
int EVSRandTiming(EVSCH *EvSch, float *EvDur,
                  float tRes, float tMax, float tPreScan)
{
  float t,tStimTot,tScanTot,tNullTot,tNullAvg,tNullSum;
  float tNull0, tNull;
  int id,n;

  tStimTot = 0.0;
  EvSch->EvDur = (float *) calloc(sizeof(float),EvSch->nEvTypes);
  for (n=0; n < EvSch->nEvTypes; n++)
  {
    tStimTot += (EvSch->nEvReps[n]*EvDur[n]);
    EvSch->EvDur[n] = EvDur[n];
  }

  /* Make sure the stimulation time does not exceed the scanning time */
  tScanTot = tMax + tPreScan;
  if (tStimTot > tScanTot)
  {
    printf("ERROR: too much stim time (%g/%g)\n",tStimTot,tScanTot);
    return(1);
  }

  /* Prep for the loop */
  tNullTot = tScanTot-tStimTot;
  tNullAvg = tNullTot/EvSch->nevents;
  tNullSum = 0.0;
  t = -tPreScan;

  /* Loop through all the events */
  for (n=0;n < EvSch->nevents; n++)
  {

    id = EvSch->eventid[n];

    /* Set the event time to the current time */
    EvSch->tevent[n] = t;

    /* Add the event duration to the current time */
    t += EvSch->EvDur[id-1];

    /* Add a random amount of null to the current time */
    if (tNullSum < tNullTot)
    {
      tNull0 = 2*drand48()*tNullAvg; /* will average tNullAvg */
      if (tRes > 0) tNull = tRes * rint(tNull0/tRes);
      else         tNull = tNull0;
      if (tNull + tNullSum > tNullTot)
      {
        /* Not much left, give it all to this null */
        tNull0 = tNullTot-tNullSum;
        if (tRes > 0) tNull = tRes * rint(tNull0/tRes);
        else         tNull = tNull0;
      }
      t += tNull;
      tNullSum += tNull;
    }
  }

  return(0);
}

/*-----------------------------------------------------------------*/
EVENT_SCHEDULE *SynthEvSch(int nEvTypes, int *nPer, float *tPer,
                           float tRes, float tMax, float tPreScan)
{
  int r,m,n,nevents, nslots, *EvSeq;
  float tStimTot,t,tScan;
  EVENT_SCHEDULE *EvSch;

  tStimTot = 0.0;
  nevents = 0;
  for (n=0;n<nEvTypes;n++)
  {
    nevents += nPer[n];
    tStimTot += (nPer[n]*tPer[n]);
  }

  tScan = tMax + tPreScan;
  if (tStimTot > tScan)
  {
    printf("ERROR: too much stim time (%g/%g)\n",tStimTot,tScan);
    return(NULL);
  }

  nslots = (int)floor(tScan/tRes);

  EvSeq = (int *) calloc(sizeof(int),nslots);

  /* Create a non-random sequence */
  r = 0;
  for (n=0;n<nEvTypes;n++)
  {
    for (m=0; m < nPer[n]; m++)
    {
      EvSeq[r] = n+1;
      r++;
    }
  }

  /* Randomize the sequence */
  RandPermList(nslots,EvSeq);

  EvSch = EVSAlloc(nevents,0);

  /* Randomize the timing */
  m = 0;
  t = -tPreScan;
  for (n=0; n < nslots; n++)
  {
    if (EvSeq[n] != 0)
    {
      EvSch->eventid[m] = EvSeq[n];
      EvSch->tevent[m]  = t;
      m++;
      t += tPer[EvSeq[n]-1];
    }
    else t += tRes;
  }

  free(EvSeq);

  EvSch->nEvTypes = nEvTypes;
  EvSch->nEvReps = (int *) calloc(sizeof(int),nEvTypes);
  EvSch->EvDur   = (float *) calloc(sizeof(float),nEvTypes);
  for (n=0;n<nEvTypes;n++)
  {
    EvSch->nEvReps[n] = nPer[n];
    EvSch->EvDur[n]   = tPer[n];
  }

  return(EvSch);
}
/*-------------------------------------------------------------*/
EVENT_SCHEDULE *RandEvSch(int nevents, int ntypes, float dtmin,
                          float dtnullavg, int randweights)
{
  EVENT_SCHEDULE *EvSch;
  float t;
  int n;

  EvSch = EVSAlloc(nevents,randweights);

  t = 0;
  for (n=0; n < EvSch->nevents; n++)
  {
    EvSch->tevent[n] = t;
    EvSch->eventid[n] = (int) floor(drand48()*ntypes) + 1;
    if (EvSch->weight != 0) EvSch->weight[n] = drand48();
    t += dtmin + drand48()*dtnullavg;
  }
  return(EvSch);
}

#endif
