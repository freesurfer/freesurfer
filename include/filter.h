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

#include "image.h"

IMAGE    *ImageNitShiFilter(IMAGE *Isrc, IMAGE *Ix, IMAGE *Iy, int wsize, 
                           double sigma, IMAGE *Idst) ;
IMAGE    *ImageBlur(IMAGE *Isrc, double sigma, IMAGE *Idst) ;

#define FILTER_NONE          0
#define FILTER_MEDIAN        1
#define FILTER_EXP_LAPLACIAN 2
#define FILTER_EXP_SUM       3
#define FILTER_GAUSSIAN      4
#define FILTER_NITSHI        5
#define FILTER_EXPONENTIAL   6

#endif
