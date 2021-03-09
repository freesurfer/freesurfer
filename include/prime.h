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



#ifndef PRIMES_H
#define PRIMES_H

int *GetPrimes(int Nmax, int *Nprimes);
int *GetPrimeFactors(int N, int *Nfactors);
int GetMaxPrimeFactor(int N);
int GetClosestPrimeFactor(int N, int P);
int GetClosestPrimeFactorLess(int N, int P);
int IsPrime(int N);


#endif /*PRIMES_H*/

