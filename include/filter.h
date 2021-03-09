/*
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


/*
 *       FILE NAME:   filter.h
 *
 *       DESCRIPTION: image processing filter prototypes
 *
 *       AUTHOR:      Bruce Fischl
 *       DATE:        6/18/96
 *
*/

#ifndef FILTER_H
#define FILTER_H

/* filter constants must be defined before inclusion of image.h */
#define FILTER_NONE                    0
#define FILTER_MEDIAN                  1
#define FILTER_EXP_SUM                 2
#define FILTER_GAUSSIAN                3
#define FILTER_NITSHI                  4
#define FILTER_GFA                     5
#define FILTER_DIFFUSE                 6
#define FILTER_EXP_LAPLACIAN           7
#define FILTER_LAPLACIAN               8
#define FILTER_EXPONENTIAL             9
#define FILTER_DIFFUSE_CURV            10
#define FILTER_DIFFUSE_GRAD            FILTER_DIFFUSE
#define FILTER_OFFSET_SCALE            11
#define FILTER_SIGMA                   12
#define FILTER_DIFFUSE_HV              13
#define FILTER_MEAN                    14
#define FILTER_SOBEL                   15
#define FILTER_ERODE                   16
#define FILTER_DILATE                  17
#define FILTER_OPEN                    18
#define FILTER_CLOSE                   19
#define FILTER_MINMAX                  20
#define FILTER_DIRECTION               21
#define FILTER_ERODE6                  22
#define FILTER_DILATE6                 23
#define FILTER_OPEN6                   24
#define FILTER_CLOSE6                  25
#define FILTER_PLANE_OF_LEAST_VARIANCE 26
#define FILTER_POLV_MEDIAN             27
#define FILTER_POLV_STD                28
#define FILTER_POLV_MEAN               29
#define FILTER_POLV_NORMAL_CURVATURE   30
#define FILTER_POLV_CURVATURE          31
#define FILTER_CENTRAL_PLANE_OF_LEAST_VARIANCE 32
#define FILTER_CPOLV_MEDIAN            33
#define FILTER_CPOLV_STD               34
#define FILTER_CPOLV_MEAN              35
#define FILTER_CPOLV_NORMAL_CURVATURE  36
#define FILTER_CPOLV_CURVATURE         37
#define FILTER_CPOLV_ORDER             38
#define FILTER_POLV_ORDER              39
#define FILTER_POLV_ZSCORE             40
#define FILTER_CPOLV_ZSCORE            41
#define FILTER_WM_STRAND_SIZE          42
#define FILTER_FILL                    43
#define FILTER_LOG                     44
#define FILTER_DOG                     45
#define FILTER_EXTERNAL                46
#define FILTER_MEAN_MASKED             47
#define FILTER_HISTO_MATCH             48

#define FILTER_OFFSET            0x0100

#define FILTER_MEDIAN_OFFSET     (FILTER_MEDIAN | FILTER_OFFSET)
#define FILTER_MEAN_OFFSET       (FILTER_MEAN | FILTER_OFFSET)
#define FILTER_GAUSSIAN_OFFSET   (FILTER_GAUSSIAN | FILTER_OFFSET)

#define FILTER_OFFSET_MEDIAN     FILTER_MEDIAN_OFFSET
#define FILTER_OFFSET_MEAN       FILTER_MEAN_OFFSET
#define FILTER_OFFSET_GAUSSIAN   FILTER_GAUSSIAN_OFFSET


#include "image.h"

IMAGE    *ImageNitShiFilter(IMAGE *Isrc, IMAGE *Ix, IMAGE *Iy, int wsize,
                            double sigma, IMAGE *Idst) ;
IMAGE    *ImageGreyErode(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE    *ImageGreyDilate(IMAGE *Isrc, IMAGE *Idst) ;
IMAGE    *ImageLOGFilter(IMAGE *Isrc, float sigma, IMAGE *Idst) ;
IMAGE    *ImageDOGFilter(IMAGE *Isrc, float psigma, float nsigma,IMAGE *Idst) ;

#endif
