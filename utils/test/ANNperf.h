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


//----------------------------------------------------------------------
// File:  ANNperf.h
// Programmer: Sunil Arya and David Mount
// Last modified: 03/04/98 (Release 0.1)
// Description: Include file for ANN performance stats
//
// Some of the code for statistics gathering has been adapted
// from the SmplStat.h package in the g++ library.
//----------------------------------------------------------------------
// Copyright (c) 1997-1998 University of Maryland and Sunil Arya and David
// Mount.  All Rights Reserved.
//
// This software and related documentation is part of the
// Approximate Nearest Neighbor Library (ANN).
//
// Permission to use, copy, and distribute this software and its
// documentation is hereby granted free of charge, provided that
// (1) it is not a component of a commercial product, and
// (2) this notice appears in all copies of the software and
//     related documentation.
//
// The University of Maryland (U.M.) and the authors make no representations
// about the suitability or fitness of this software for any purpose.  It is
// provided "as is" without express or implied warranty.
//----------------------------------------------------------------------

#ifndef ANNperf_H
#define ANNperf_H

//----------------------------------------------------------------------
//  basic includes
//----------------------------------------------------------------------
#include <stdlib.h>   // standard libaries
#include <stdio.h>   // standard I/O
#include <iostream.h>   // C++ I/O streams
#include <math.h>   // math routines
#include <values.h>

//----------------------------------------------------------------------
// kd-tree stats object
// This object is used for collecting information about a kd-tree
// or bd-tree.
//----------------------------------------------------------------------

class ANNkdStats
{   // stats on kd-tree
public:
  int dim;    // dimension of space
  int n_pts;    // no. of points
  int bkt_size;   // bucket size
  int n_lf;    // no. of leaves
  int n_tl;    // no. of trivial leaves (no points)
  int n_spl;    // no. of splitting nodes
  int n_shr;    // no. of shrinking nodes (for bd-trees)
  int depth;    // depth of tree
  float sum_ar;   // sum of leaf aspect ratios
  float avg_ar;   // average leaf aspect ratio
//
  // reset stats
  void reset(int d=0, int n=0, int bs=0)
  {
    dim = d;
    n_pts = n;
    bkt_size = bs;
    n_lf = n_tl = n_spl = n_shr = depth = 0;
    sum_ar = avg_ar = 0.0;
  }

  ANNkdStats()   // basic constructor
  { reset();
  }

  void merge(const ANNkdStats &st); // merge stats from child
};

//----------------------------------------------------------------------
//    SampStat
//
// A sample stat collects numeric (double) samples and
// returns some simple statistics.  Its mail functions
// are:
//
//  reset()  Reset to no samples.
//  += x  Include sample x.
//  samples() Return number of samples.
//  mean()  Return mean of samples.
//  stdDev() Return standard deviation
//  min()  Return minimum of samples.
//  max()  Return maximum of samples.
//----------------------------------------------------------------------
class SampStat
{
  int n;    // number of samples
  double sum;    // sum
  double sum2;   // sum of squares
  double minVal, maxVal;  // min and max
public :
  void reset()   // reset everything
  {  n = 0;
    sum = sum2 = 0;
    minVal = MAXDOUBLE;
    maxVal = -MAXDOUBLE;
  }

  SampStat()
  {
    reset();
  }  // constructor

  void operator+=(double x)  // add sample
  {
    n++;
    sum += x;
    sum2 += x*x;
    if (x < minVal) minVal = x;
    if (x > maxVal) maxVal = x;
  }

  int  samples()
  {
    return n;
  } // number of samples

  double mean()
  {
    return sum/n;
  } // mean

  // standard deviation
  double stdDev()
  {
    return sqrt((sum2 - (sum*sum)/n)/(n-1));
  }

  double min()
  {
    return minVal;
  } // minimum
  double max()
  {
    return maxVal;
  } // maximum
};

//----------------------------------------------------------------------
//      Operation count updates
//----------------------------------------------------------------------

#ifdef PERF
#define FLOP(n)  {N_float_ops += (n);}
#define LEAF(n)  {N_visit_lfs += (n);}
#define SPL(n)   {N_visit_spl += (n);}
#define SHR(n)   {N_visit_shr += (n);}
#define PTS(n)   {N_visit_pts += (n);}
#define COORD(n) {N_coord_hts += (n);}
#else
#define FLOP(n)
#define LEAF(n)
#define SPL(n)
#define SHR(n)
#define PTS(n)
#define COORD(n)
#endif

//----------------------------------------------------------------------
//  Performance statistics
// The following data and routines are used for computing
// performance statistics for nearest neighbor searching.
// Because these routines can slow the code down, they can be
// activated and deactiviated by defining the PERF variable,
// by compiling with the option: -DPERF
//----------------------------------------------------------------------

//----------------------------------------------------------------------
//  Global counters for performance measurement
//
// visit_lfs  The number of leaf nodes visited in the
//    tree.
//
// visit_spl  The number of splitting nodes visited in the
//    tree.
//
// visit_shr  The number of shrinking nodes visited in the
//    tree.
//
// visit_pts  The number of points visited in all the
//    leaf nodes visited.  Equivalently, this
//    is the number of points for which distance
//    calculations are performed.
//
// coord_hts  The number of times a coordinate of a
//    data point is accessed.  This is generally
//    less than visit_pts*d if partial distance
//    calculation is used.  This count is low
//    in the sense that if a coordinate is hit
//    many times in the same routine we may
//    count it only once.
//
// float_ops  The number of floating point operations.
//    This includes all operations in the heap
//    as well as distance calculations to boxes.
//
// average_err  The average error of each query (the
//    error of the reported point to the true
//    nearest neighbor).  For k nearest neighbors
//    the error is computed k times.
//
// rank_err  The rank error of each query (the difference
//    in the rank of the reported point and its
//    true rank).
//
// data_pts  The number of data points.  This is not
//    a counter, but used in stats computation.
//----------------------------------------------------------------------

extern int N_data_pts; // number of data points
extern int N_visit_lfs; // number of leaf nodes visited
extern int N_visit_spl; // number of splitting nodes visited
extern int N_visit_shr; // number of shrinking nodes visited
extern int N_visit_pts; // visited points for one query
extern int N_coord_hts; // coordinate hits for one query
extern int N_float_ops; // floating ops for one query
extern SampStat visit_lfs; // stats on leaf nodes visits
extern SampStat visit_spl; // stats on splitting nodes visits
extern SampStat visit_shr; // stats on shrinking nodes visits
extern SampStat visit_nds; // stats on total nodes visits
extern SampStat visit_pts; // stats on points visited
extern SampStat coord_hts; // stats on coordinate hits
extern SampStat float_ops; // stats on floating ops
extern SampStat average_err; // average error
extern SampStat rank_err; // rank error

//----------------------------------------------------------------------
//  Declaration of externally accessible routines for statistics
//----------------------------------------------------------------------

void reset_stats(int data_size); // reset stats for a set of queries

void reset_counts();   // reset counts for one queries

void update_stats();   // update stats with current counts

void print_stats(ANNbool validate); // print statistics for a run

#endif
