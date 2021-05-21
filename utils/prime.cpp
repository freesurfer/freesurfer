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

#include "prime.h"

static int CompareFactors(const void *pf1, const void *pf2);

/* -------------------------------------------------------
   GetPrimes() - returns a list of prime numbers less than
   or equal to Nmax. The number of primes is returned
   in *Nprimes;
 ------------------------------------------------------- */
int *GetPrimes(int Nmax, int *Nprimes)
{
  int *numlist, n, m, k, k2, kmax;
  int *primes;

  if (Nmax < 1) {
    printf("ERROR: GetPrimes: cannot compute primes for %d\n", Nmax);
    return (NULL);
  }

  numlist = (int *)calloc(Nmax + 1, sizeof(int));

  for (n = 0; n <= Nmax; n++) numlist[n] = n;
  numlist[1] = 0;

  kmax = (int)sqrt(Nmax);
  for (k = 2; k <= kmax; k++) {
    if (numlist[k] != 0) {
      k2 = k * k;
      for (m = k2; m <= Nmax; m += k) numlist[m] = 0;
    }
  }

  *Nprimes = 0;
  for (n = 0; n <= Nmax; n++)
    if (numlist[n] != 0) (*Nprimes)++;
  // printf("INFO: found %d primes <=  %d\n",*Nprimes,Nmax);

  primes = (int *)calloc(*Nprimes, sizeof(int));

  m = 0;
  for (n = 0; n <= Nmax; n++) {
    if (numlist[n] != 0) {
      primes[m] = numlist[n];
      // printf("%3d %5d\n",m+1,primes[m]);
      m++;
    }
  }

  free(numlist);

  return (primes);
}
/* --------------------------------------------------------------
   GetPrimeFactors() - returns a list of prime factors of N. The
   number of prime factors is returned in *Nfactors. The factors
   are sorted from lowest to highest.
   ------------------------------------------------------- */
int *GetPrimeFactors(int N, int *Nfactors)
{
  int *pfactors, n, m, NN;
  int *allprimes, nallprimes;
  int changed;

  if (N < 1) {
    printf("ERROR: cannot compute prime factors of %d\n", N);
    return (NULL);
  }

  /* Get list of all primes below sqrt(N) */
  allprimes = GetPrimes((int)floor(sqrt(N)), &nallprimes);

  /* Go through each prime to see if it is a factor of N */
  *Nfactors = 0;
  changed = 1;
  NN = N;
  while (changed) {
    changed = 0;
    for (n = 0; n < nallprimes; n++) {
      if ((NN % allprimes[n]) == 0) {
        (*Nfactors)++;
        NN /= allprimes[n];
        changed = 1;
      }
    }
  }
  if (NN != 1) (*Nfactors)++;
  // printf("INFO: found %d prime factors for %d\n",*Nfactors,N);

  /* Now go back through and record what the factors are */
  pfactors = (int *)calloc(*Nfactors, sizeof(int));
  m = 0;
  changed = 1;
  NN = N;
  while (changed) {
    changed = 0;
    for (n = 0; n < nallprimes; n++) {
      if ((NN % allprimes[n]) == 0) {
        pfactors[m] = allprimes[n];
        NN /= allprimes[n];
        m++;
        changed = 1;
      }
    }
  }
  if (NN != 1) pfactors[m] = NN;

  free(allprimes);

  qsort(pfactors, *Nfactors, sizeof(int), CompareFactors);

  // printf("N = %d\n",N);
  // for(m=0;m<*Nfactors;m++) printf("%2d %d\n",m,pfactors[m]);
  // exit(0);

  return (pfactors);
}
/* --------------------------------------------------------------
   CompareFactors() - only for sorting.
   ------------------------------------------------------------- */
static int CompareFactors(const void *pf1, const void *pf2)
{
  int f1, f2;

  f1 = *((int *)pf1);
  f2 = *((int *)pf2);

  if (f1 < f2) return (-1);
  if (f1 > f2) return (+1);
  return (0);
}
/* --------------------------------------------------------------
   IsPrime() - returns 1 if N is prime, zero otherwise
   ------------------------------------------------------------- */
int IsPrime(int N)
{
  int *primes, nprimes, r;

  if (N < 1) {
    printf("ERROR: cannot compute primes for %d\n", N);
    return (0);
  }

  /* Get list of all primes below sqrt(N) */
  primes = GetPrimes(N, &nprimes);

  r = 0;
  if (primes[nprimes - 1] == N) r = 1;

  free(primes);
  return (r);
}
/* --------------------------------------------------------------
   GetMaxPrimeFactor() - returns the maximum prime factor of N.
   ------------------------------------------------------------- */
int GetMaxPrimeFactor(int N)
{
  int *pfactors, nfactors;
  int maxfactor;

  pfactors = GetPrimeFactors(N, &nfactors);

  maxfactor = pfactors[nfactors - 1];
  free(pfactors);
  return (maxfactor);
}
/* --------------------------------------------------------------
   GetClosestPrimeFactor() - returns the prime factor of N
   closest to P.
   ------------------------------------------------------------- */
int GetClosestPrimeFactor(int N, int P)
{
  int *pfactors, nfactors;
  int n, d, dmin, nmin, fclosest;

  pfactors = GetPrimeFactors(N, &nfactors);

  nmin = 0;
  dmin = abs(pfactors[nmin] - P);
  for (n = 0; n < nfactors; n++) {
    d = abs(pfactors[n] - P);
    if (dmin > d) {
      dmin = d;
      nmin = n;
    }
  }
  fclosest = pfactors[nmin];

  free(pfactors);

  return (fclosest);
}
/* --------------------------------------------------------------
   GetClosestPrimeFactorLess() - returns the prime factor of N
   closest to and less than P.
   ------------------------------------------------------------- */
int GetClosestPrimeFactorLess(int N, int P)
{
  int *pfactors, nfactors;
  int n, d, dmin, nmin, fclosest;

  pfactors = GetPrimeFactors(N, &nfactors);

  nmin = 0;
  dmin = abs(pfactors[nmin] - P);
  for (n = 0; n < nfactors; n++) {
    d = abs(pfactors[n] - P);
    if (dmin > d && pfactors[n] < P) {
      dmin = d;
      nmin = n;
    }
  }
  fclosest = pfactors[nmin];

  free(pfactors);

  return (fclosest);
}
#if 0
/* --------------------------------------------------------------
   GetClosestFactor() - returns the factor of N closest to, and
   greater than, P. Note that factor may not be prime.
   ------------------------------------------------------------- */
int GetClosestFactor(int N, int P)
{}
#endif
