#ifndef SIG_H
#define SIG_H

double sigt(double t, int df) ;
float betacf(float a, float b, float x) ;
float gammln(float xx) ;
float betai(float a, float b, float x) ;
float sigchisq(double chisq, int df) ;
void gser(float *gamser,float a,float x, float *gln) ;
void gcf(float *gammcf, float a, float x, float *gln) ;
float gammp(float a,float x) ;

double fdrthreshold(double *p, int np, double fdr);
int doublecompar(const void *v1, const void *v2);

#endif
