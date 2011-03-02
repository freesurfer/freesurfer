/**
 * @file  cdflib.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:09 $
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


#ifndef CDFLIB_H
#define CDFLIB_H

double algdiv(double*,double*);
double alngam(double*);
double alnrel(double*);
double apser(double*,double*,double*,double*);
double basym(double*,double*,double*,double*);
double bcorr(double*,double*);
double betaln(double*,double*);
double bfrac(double*,double*,double*,double*,double*,double*);
void bgrat(double*,double*,double*,double*,double*,double*,int*i);
double bpser(double*,double*,double*,double*);
void bratio(double*,double*,double*,double*,double*,double*,int*);
double brcmp1(int*,double*,double*,double*,double*);
double brcomp(double*,double*,double*,double*);
double bup(double*,double*,double*,double*,int*,double*);
void cdfbet(int*,double*,double*,double*,double*,double*,double*,
            int*,double*);
void cdfbin(int*,double*,double*,double*,double*,double*,double*,
            int*,double*);
void cdfchi(int*,double*,double*,double*,double*,int*,double*);
void cdfchn(int*,double*,double*,double*,double*,double*,int*,double*);
void cdff(int*,double*,double*,double*,double*,double*,int*,double*);
void cdffnc(int*,double*,double*,double*,double*,double*,double*,
            int*s,double*);
void cdfgam(int*,double*,double*,double*,double*,double*,int*,double*);
void cdfnbn(int*,double*,double*,double*,double*,double*,double*,
            int*,double*);
void cdfnor(int*,double*,double*,double*,double*,double*,int*,double*);
void cdfpoi(int*,double*,double*,double*,double*,int*,double*);
void cdft(int*,double*,double*,double*,double*,int*,double*);
void cdftnc(int*,double*,double*,double*,double*,double*,int*,double*);
void cumbet(double*,double*,double*,double*,double*,double*);
void cumbin(double*,double*,double*,double*,double*,double*);
void cumchi(double*,double*,double*,double*);
void cumchn(double*,double*,double*,double*,double*);
void cumf(double*,double*,double*,double*,double*);
void cumfnc(double*,double*,double*,double*,double*,double*);
void cumgam(double*,double*,double*,double*);
void cumnbn(double*,double*,double*,double*,double*,double*);
void cumnor(double*,double*,double*);
void cumpoi(double*,double*,double*,double*);
void cumt(double*,double*,double*,double*);
void cumtnc(double*,double*,double*,double*,double*);
double devlpl(double [],int*,double*);
double dinvnr(double *p,double *q);
void dinvr(int*,double*,double*,unsigned long*,unsigned long*);
void dstinv(double*,double*,double*,double*,double*,double*,
            double*);
double dt1(double*,double*,double*);
void dzror(int*,double*,double*,double*,double *,
           unsigned long*,unsigned long*);
void dstzr(double *zxlo,double *zxhi,double *zabstl,double *zreltl);
double erf1(double*);
double erfc1(int*,double*);
double esum(int*,double*);
double exparg(int*);
double fpser(double*,double*,double*,double*);
double gam1(double*);
void gaminv(double*,double*,double*,double*,double*,int*);
double gamln(double*);
double gamln1(double*);
double Xgamm(double*);
void grat1(double*,double*,double*,double*,double*,double*);
void gratio(double*,double*,double*,double*,int*);
double gsumln(double*,double*);
double psi(double*);
double rcomp(double*,double*);
double rexp(double*);
double rlog(double*);
double rlog1(double*);
double spmpar(int*);
double stvaln(double*);
double fifdint(double);
double fifdmax1(double,double);
double fifdmin1(double,double);
double fifdsign(double,double);
long fifidint(double);
long fifmod(long,long);
void ftnstop(char*);
extern int ipmpar(int*);

#if 0
static void E0000(int,int*,double*,double*,unsigned long*,
                  unsigned long*,double*,double*,double*,
                  double*,double*,double*,double*);
static void E0001(int,int*,double*,double*,double*,double*,
                  unsigned long*,unsigned long*,double*,double*,
                  double*,double*);
#endif

double sigf(float F, int df1, int df2) ;

#endif
