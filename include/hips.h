/*
 *       FILE NAME:   hips.h
 *
 *       DESCRIPTION: hips2 prototypes
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        2/5/96
 *
*/

#ifndef HIPS_H
#define HIPS_H

/*-----------------------------------------------------
                    INCLUDE FILES
-------------------------------------------------------*/
#include <hipl_format.h>

int h_mulscale(struct header *hdi, struct header *hdo, Pixelval *b) ;
int h_divscale(struct header *hdi, struct header *hdo, Pixelval *b) ;
int h_add(struct header *hdi1, struct header *hdi2, struct header *hdo) ;
int h_linscale(struct header *hdi, struct header * hdo, float b, float c) ;
int h_div(struct header *hdi1, struct header * hdi2, struct header *hdo) ;
int h_tof(struct header *hdi, struct header *hdo) ;
int h_toc(struct header *hdi, struct header *hdo) ;
int h_reduce(struct header *hdi, struct header *hdo, float xf, float yf) ;
int h_flipquad(struct header *hdi, struct header *hdo) ;
int h_rot180(struct header *hdi, struct header *hdo) ;
int h_minmax(struct header *hd, Pixelval *minval, Pixelval *maxval,int nzflag);
struct header *h_logmap(struct header *phd,float a,int R,int rows,int cols,
      int Fx,int Fy,float disk_size) ;
struct header *h_logmap_f(struct header *phd,float a,int R,int rows,int cols,
      int Fx,int Fy,float disk_size) ;
int h_discedge(struct header *hdi, struct header *hdo,int size, float varcrit);
int h_canny(IMAGE *Isrc, IMAGE *Idst, double sigma, int mask_size,double lfrac,double hfrac,int dothin);
int h_morphdil(struct header *hdi, struct header *hde, struct header *hdo, int centerr, int centerc,
               int gray) ;

#endif
