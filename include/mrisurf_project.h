#pragma once
/*
 * @file utilities dealing with the topology 
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ © copyright-2014,2018 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mrisurf_metricProperties.h"
void mrisSphericalProjectXYZ(float xs, float ys, float zs, float *xd, float *yd, float *zd);

void mrisSphericalProjection(MRIS *mris);

void MRISMP_projectOntoSphere(MRIS_MP* mris, double r);

void mrisAssignFaces(MRIS* mris, MHT *mht, int which_vertices);
