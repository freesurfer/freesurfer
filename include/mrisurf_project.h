#pragma once
/*
 *
 */
/*
 * surfaces Author: Bruce Fischl, extracted from mrisurf.c by Bevin Brett
 *
 * $ Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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

void mrisSphericalProjection(MRIS   * mris);
void mrisAssignFaces(MRIS* mris, MHT *mht, int which_vertices);


MRIS* MRISprojectOntoSphere   (MRIS* mris_src, MRIS* mris_dst,  double r) ;     // for now, mris_src and mris_dst must be the same pointer
MRIS* MRISprojectOntoEllipsoid(MRIS* mris_src, MRIS* mris_dst, float a, float b, float c) ;


// Ones that are supported by MRIS and MRIS_MP
//
MRIS*    MRISprojectOntoSphere(MRIS*    mris, double r);
MRIS_MP* MRISprojectOntoSphere(MRIS_MP* mris, double r);
