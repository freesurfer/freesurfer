#ifndef GENERAL_H
#define GENERAL_H

void order(int *small, int *big) ;

#if 0
/* putting these prototypes in causes the canny edge detecter to fail, so 
   leave them out for now */
short nearestshort(float x) ;
void  fporder(float *big, float *small) ;
float fpclip(float x, float bound1, float bound2) ;
int nearestint(float x) ;
#endif

#endif
