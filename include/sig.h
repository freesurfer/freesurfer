#ifndef SIG_H
#define SIG_H

double fdr2vwth(double *p, int np, double fdr);
double vwth2fdr(double *p, int np, double vwth);
int doublecompar(const void *v1, const void *v2);
double sigt(double t, int df) ;
float sigchisq(double chisq, int df) ;

#endif
