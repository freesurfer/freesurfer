/**
 * @file  prime.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:10 $
 *    $Revision: 1.4 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

