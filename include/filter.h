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
#define FILTER_NONE              0
#define FILTER_MEDIAN            1
#define FILTER_EXP_SUM           2
#define FILTER_GAUSSIAN          3
#define FILTER_NITSHI            4
#define FILTER_GFA               5
#define FILTER_DIFFUSE           6
#define FILTER_EXP_LAPLACIAN     7
#define FILTER_LAPLACIAN         8
#define FILTER_EXPONENTIAL       9
#define FILTER_DIFFUSE_CURV      10
#define FILTER_DIFFUSE_GRAD      FILTER_DIFFUSE
#define FILTER_OFFSET_SCALE      11
#define FILTER_SIGMA             12
#define FILTER_DIFFUSE_HV        13
#define FILTER_MEAN              14
#define FILTER_SOBEL             15
#define FILTER_ERODE             16
#define FILTER_DILATE            17
#define FILTER_OPEN              18
#define FILTER_CLOSE             19
#define FILTER_MINMAX            20
#define FILTER_DIRECTION         21
#define FILTER_ERODE6            22
#define FILTER_DILATE6           23
#define FILTER_OPEN6             24
#define FILTER_CLOSE6            25
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

#endif
