#ifndef IUTILS_H
#define IUTILS_H

#if 1
/*
  these prototypes cause the canny edge detector to fail - I don't know why, 
  and I'm not going to spend the time to find out.
*/
void crop_image(unsigned char *imageptr, int *cols,int *rows, int cropcornerax,
                int cropcorneray, int cropcornerbx, int cropcornerby) ;
#endif
void get_histogram_threshold(int hgram[], int histsize, int pixelmax, 
                             int pixelmin, float fraction, 
                             int zflag, float *ht, float *lt) ;
void histogram(short *theimage, int xsize, int ysize, int pixelmax, 
               int pixelmin, int hgram[], int histsize) ;

#endif
