/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
 *
 * Copyright 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007
 * Washington University, Mallinckrodt Institute of Radiology.
 * All Rights Reserved.
 *
 * This software may not be reproduced, copied, or distributed without 
 * written permission of Washington University. For further information 
 * contact A. Z. Snyder.
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <endianio.h>

void flipx (float *imgf, int *pnx, int* pny, int *pnz) {
	float	*vector;
	int	ix, iy, iz, vecdim, index;

	vecdim = *pnx;
	if (!(vector = (float *) malloc (vecdim * sizeof (float)))) errm ("flipx");

	for (iz = 0; iz < *pnz; iz++) {
	for (iy = 0; iy < *pny; iy++) {
		for (ix = 0; ix < *pnx; ix++) {
			index = ix + *pnx*(iy + *pny*iz);
			vector[ix] = imgf[index];
		}
		for (ix = 0; ix < *pnx; ix++) {
			index = ix + *pnx*(iy + *pny*iz);
			imgf[index] = vector[*pnx - 1 - ix];
		}
	}}
	free (vector);
}

void flipy (float *imgf, int *pnx, int* pny, int *pnz) {
	float	*vector;
	int	ix, iy, iz, vecdim, index;

	vecdim = *pny;
	if (!(vector = (float *) malloc (vecdim * sizeof (float)))) errm ("flipy");

	for (iz = 0; iz < *pnz; iz++) {
	for (ix = 0; ix < *pnx; ix++) {
		for (iy = 0; iy < *pny; iy++) {
			index = ix + *pnx*(iy + *pny*iz);
			vector[iy] = imgf[index];
		}
		for (iy = 0; iy < *pny; iy++) {
			index = ix + *pnx*(iy + *pny*iz);
			imgf[index] = vector[*pny - 1 - iy];
		}
	}}
	free (vector);
}

void flipz (float *imgf, int *pnx, int* pny, int *pnz) {
	float	*vector;
	int	ix, iy, iz, vecdim, index;

	vecdim = *pnz;
	if (!(vector = (float *) malloc (vecdim * sizeof (float)))) errm ("flipz");

	for (iy = 0; iy < *pny; iy++) {
	for (ix = 0; ix < *pnx; ix++) {
		for (iz = 0; iz < *pnz; iz++) {
			index = ix + *pnx*(iy + *pny*iz);
			vector[iz] = imgf[index];
		}
		for (iz = 0; iz < *pnz; iz++) {
			index = ix + *pnx*(iy + *pny*iz);
			imgf[index] = vector[*pnz - 1 - iz];
		}
	}}
	free (vector);
}
