#ifndef _MOSAIC_H
#define _MOSAIC_H

int VolSS2MosSS(int cvol, int rvol, int svol,
		int ncvol, int nrvol, 
		int ncmos, int nrmos,
		int *cmos, int *rmos,
		int *OutOfBounds);

int MosSS2VolSS(int cmos,  int rmos, 
		int ncmos, int nrmos,
		int ncvol, int nrvol, int nsvol,
		int *cvol, int *rvol, int *svol,
		int *OutOfBounds);

int CheckMosaic(void);

#endif
