/*$Header: /space/repo/1/dev/dev/talairach_avi/cflip.c,v 1.1 2007/05/04 22:33:59 nicks Exp $*/
/*$Log: cflip.c,v $
/*Revision 1.1  2007/05/04 22:33:59  nicks
/*new talairach alignment utility, using Avi Snyders registration tools
/*
 * Revision 1.3  2007/04/25  05:00:57  avi
 * properly include standard includes
 * remove local errm()
 *
 * Revision 1.2  2004/11/18  20:50:57  rsachs
 * Fixed a logic error.
 *
 * Revision 1.1  2004/02/19  01:04:08  avi
 * Initial revision
 **/

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
