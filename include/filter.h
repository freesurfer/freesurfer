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
#define FILTER_NONE          0
#define FILTER_MEDIAN        1
#define FILTER_MEDIAN_OFFSET 2
#define FILTER_EXP_SUM       3
#define FILTER_GAUSSIAN      4
#define FILTER_NITSHI        5
#define FILTER_GFA           6
#define FILTER_DIFFUSE       7
#define FILTER_EXP_LAPLACIAN 8
#define FILTER_LAPLACIAN     9
#define FILTER_EXPONENTIAL   10
#define FILTER_DIFFUSE_CURV  11
#define FILTER_DIFFUSE_GRAD  FILTER_DIFFUSE
#define FILTER_OFFSET_SCALE  12
#define FILTER_SIGMA         13
#define FILTER_DIFFUSE_HV    14
#define FILTER_MEAN          15
#define FILTER_OFFSET        0x0100

#include "image.h"

IMAGE    *ImageNitShiFilter(IMAGE *Isrc, IMAGE *Ix, IMAGE *Iy, int wsize, 
                           double sigma, IMAGE *Idst) ;

#endif
