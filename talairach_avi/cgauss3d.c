/**
 * @file  cgauss3d.c
 *
 */
/*
 * Original Author: Avi Z. Snyder, Washington University
 * 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2007/05/05 00:00:06 $
 *    $Revision: 1.2 $
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

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>
#include <endianio.h>
#include <librms.h>

void gauss3d (float *image, int *pnx, int *pny, int *pnz, float *cmppix, float *fhalf) {
/*********************************************************************/
/*	3D Gaussian filter								*/
/*	Author: Avi Snyder 05-Jan-93 => C 070422 AZS					*/
/*											*/
/*	inputs: real*4 image(nx,ny,nz)	! image to be filtered				*/
/* 		real*4 cmppix(3)	! x, y, z voxel dimensions in cm		*/
/*		real*4 fhalf		! half amplitude frequency in cycles/cm		*/
/*											*/
/*	variables modified: image	! input image overwritten			*/
/*											*/
/*	subroutines called:	FFT, REALT in FORTRAN source fftsun.f or fftsol.f	*/
/* 				FFT algorithm by Richard C. Singleton			*/
/*											*/
/*	restrictions: mixed radix FFT accepts any nx, ny, nz but 			*/
/*				nx must be divisible by 2 				*/
/**********************************************************************/
	float		*a, *b;
	double		q, factor, f2, fx, fy, fz;
	int		i, k, n1, n2, n2ny, jndex;
	int		ix, iy, iz, nx, ny, nz, kx, ky, kz;
 	int		one = 1, negone = -1;

	nx = *pnx; ny = *pny; nz = *pnz;
	if (nx % 2) {
		fprintf (stderr, "gauss3d: nx not a multiple of 2\n");
		exit (-1);
	}
	n1 = nx/2;
	n2 = n1 + 1;
	n2ny = n2*ny;
	if (!(a = (float *) malloc (n2*ny*nz * sizeof (float)))) errm ("gauss3d");
	if (!(b = (float *) malloc (n2*ny*nz * sizeof (float)))) errm ("gauss3d");

	k = jndex = 0;
	for (iz = 0; iz < nz; iz++) {
	for (iy = 0; iy < ny; iy++) {
		for (i = ix = 0; ix < nx; ix += 2, i++) {
			(a + k)[i] = image[jndex++];
			(b + k)[i] = image[jndex++];
		}
		fft_   (a + k, b + k, &one, &n1, &one, &negone);
		realt_ (a + k, b + k, &one, &n1, &one, &negone);
		k += n2;
	}}
	fft_ (a, b, &nz,  &ny, &n2,   &negone);
	fft_ (a, b, &one, &nz, &n2ny, &negone);

	q = -log(2.)/((*fhalf)*(*fhalf));
	i = 0;
	for (iz = 0; iz < nz; iz++) {
		kz = (iz <= nz/2) ? iz : nz - iz;
		fz = kz/(nz*cmppix[2]);
	for (iy = 0; iy < ny; iy++) {
		ky = (iy <= ny/2) ? iy : ny - iy;
		fy = ky/(ny*cmppix[1]);
	for (ix = 0; ix < n2; ix++) {
		kx = (ix <= nx/2) ? ix : nx - ix;
		fx = kx/(nx*cmppix[0]);
		f2 = fx*fx + fy*fy + fz*fz;
		factor = exp (q*f2);
		a[i] *= factor;
 		b[i] *= factor;
		i++;
	}}}

	fft_ (a, b, &one, &nz, &n2ny, &one);
	fft_ (a, b, &nz,  &ny, &n2,   &one);
	k = jndex = 0;
	for (iz = 0; iz < nz; iz++) {
	for (iy = 0; iy < ny; iy++) {
		realt_ (a + k, b + k, &one, &n1, &one, &one);
		fft_   (a + k, b + k, &one, &n1, &one, &one);
		for (i = ix = 0; ix < nx; ix += 2, i++) {
			image[jndex++] = (a + k)[i];
			image[jndex++] = (b + k)[i];
		}
		k += n2;
	}}

	free (a);
	free (b);
}
