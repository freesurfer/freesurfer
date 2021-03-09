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
// File:  ANN.h
// Programmer: Sunil Arya and David Mount
// Last modified: 03/04/98 (Release 0.1)
// Description: Basic include file for approximate nearest
//   neighbor searching.
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

//----------------------------------------------------------------------
// ANN - approximate nearest neighbor searching
// ANN is a library for approximate nearest neighbor searching,
// based on the use of standard and priority search in kd-trees
// and balanced box-decomposition (bbd) trees.  Here are some
// references:
//
// kd-trees:
//          Friedman, Bentley, and Finkel, ``An algorithm for finding
//  best matches in logarithmic expected time,'' ACM
//  Transactions on Mathematical Software, 3(3):209-226, 1977.
//
// priority search in kd-trees:
//          Arya and Mount, ``Algorithms for fast vector quantization,''
//  Proc. of DCC '93: Data Compression Conference, eds. J. A.
//  Storer and M. Cohn, IEEE Press, 1993, 381-390.
//
// approximate nearest neighbor search and bbd trees:
//     Arya, Mount, Netanyahu, Silverman, and Wu, ``An optimal
//  algorithm for approximate nearest neighbor searching,''
//  5th Ann. ACM-SIAM Symposium on Discrete Algorithms,
//  1994, 573-582.
//----------------------------------------------------------------------

#ifndef ANN_H
#define ANN_H

//----------------------------------------------------------------------
//  basic includes
//----------------------------------------------------------------------
#include <stdlib.h>   // standard libs
#include <stdio.h>   // standard I/O (for NULL)
#include <iostream>   // I/O streams
#include <math.h>   // math includes
#include <values.h>   // special values

#define ANNversion "0.1"  // ANN version number

using namespace std;

//----------------------------------------------------------------------
//  ANNbool
// This is a simple boolean type.  Although ANSI C++ is supposed
// to support the type bool, many compilers do not have it.
//----------------------------------------------------------------------

// ANN boolean type (non ANSI C++)
enum ANNbool {ANNfalse = 0, ANNtrue = 1};

//----------------------------------------------------------------------
//  Basic Types:  ANNcoord, ANNdist, ANNidx
// ANNcoord and ANNdist are the types used for representing
// point coordinates and distances.  They can be modified by the
// user, with some care.  It is assumed that they are both numeric
// types, and that ANNdist is generally of an equal or higher type
// from ANNcoord.  A variable of type ANNdist should be large
// enough to store the sum of squared components of a variable
// of type ANNcoord for the number of dimensions needed in the
// application.  For example, the following combinations are
// legal:
//
//  ANNcoord ANNdist
//  --------- -------------------------------
//  short  short, int, long, float, double
//  int  int, long, float, double
//  long  long, float, double
//  float  float, double
//  double  double
//
// It is the user's responsibility to make sure that overflow does
// not occur in distance calculation.
//
// The code assumes that there is an infinite distance, ANN_DIST_INF
// (as large as any legal distance).  Possible values are given below:
//
//     Examples:
//     ANNdist:  double, float, long, int, short
//     ANN_DIST_INF: MAXDOUBLE, MAXFLOAT, MAXLONG, MAXINT, MAXSHORT
//
//
// ANNidx is a point index.  When the data structure is built,
// the points are given as an array.  Nearest neighbor results are
// returned as an index into this array.  To make it clearer when
// this is happening, we define the integer type ANNidx.
//
//----------------------------------------------------------------------

typedef double ANNcoord;  // coordinate data type
typedef double ANNdist;  // distance data type
typedef int ANNidx;   // point index

// largest possible distance
const ANNdist ANN_DIST_INF =  MAXDOUBLE;

//----------------------------------------------------------------------
// Self match?
// In some applications, the nearest neighbor of a point is not
// allowed to be the point itself.  This occurs, for example, when
// computing all nearest neighbors in a set.  By setting the
// parameter ANN_ALLOW_SELF_MATCH to ANNfalse, the nearest neighbor
// is the closest point whose distance from the query point is
// strictly positive.
//----------------------------------------------------------------------

const ANNbool ANN_ALLOW_SELF_MATCH = ANNtrue;

//----------------------------------------------------------------------
//  Norms and metrics:
// ANN supports any Minkowski norm for defining distance.  In
// particular, for any p >= 1, the L_p Minkowski norm defines the
// length of a d-vector (v0, v1, ..., v(d-1)) to be
//
//  (|v0|^p + |v1|^p + ... + |v(d-1)|^p)^(1/p),
//
// (where ^ denotes exponentiation, and |.| denotes absolute
// value).  The distance between two points is defined to be
// the norm of the vector joining them.  Some common distance
// metrics include
//
//  Euclidean metric p = 2
//  Manhattan metric p = 1
//  Max metric  p = infinity
//
// In the case of the max metric, the norm is computed by
// taking the maxima of the absolute values of the components.
// ANN is highly "coordinate-based" and does not support general
// distances functions (e.g. those obeying just the triangle
// inequality).  It also does not support distance functions
// based on inner-products.
//
// For the purpose of computing nearest neighbors, it is not
// necessary to compute the final power (1/p).  Thus the only
// component that is used by the program is |v(i)|^p.
//
// ANN parameterizes the distance computation through the following
// macros.  (Macros are used rather than procedures for efficiency.)
// Recall that the distance between two points is given by the length
// of the vector joining them, and the length or norm of a vector v
// is given by formula:
//
//  |v| = ROOT(POW(v0) # POW(v1) # ... # POW(v(d-1)))
//
// where ROOT, POW are unary functions and # is an associative and
// commutative binary operator satisfying:
//
//     ** POW: coord  --> dist
//     ** #: dist x dist --> dist
//     ** ROOT: dist (>0) --> double
//
// For early termination in distance calculation (partial distance
// calculation) we assume that POW and # together are monotonically
// increasing on sequences of arguments, meaning that for all
// v0..vk and y:
//
// POW(v0) #...# POW(vk) <= (POW(v0) #...# POW(vk)) # POW(y).
//
// Due to the use of incremental distance calculations in the code
// for searching k-d trees, we assume that there is an incremental
// update function DIFF(x,y) for #, such that if:
//
//      s = x0 # ... # xi # ... # xk
//
// then if s' is s with xi replaced by y, that is,
//
//      s' = x0 # ... # y # ... # xk
//
// can be computed by:
//
//      s' = s # DIFF(xi,y).
//
// Thus, if # is + then DIFF(xi,y) is (yi-x).  For the L_infinity
// norm we make use of the fact that in the program this function
// is only invoked when y > xi, and hence DIFF(xi,y)=y.
//
// Finally, for approximate nearest neighbor queries we assume
// that POW and ROOT are related such that
//
//      v*ROOT(x) = ROOT(POW(v)*x)
//
// Here are the values for the various Minkowski norms:
//
// L_p: p even:    p odd:
//  ------------------------- ------------------------
//  POW(v)  = v^p  POW(v)  = |v|^p
//  ROOT(x)  = x^(1/p) ROOT(x)  = x^(1/p)
//  #  = +  #  = +
//  DIFF(x,y) = y - x  DIFF(x,y) = y - x
//
// L_inf:
//  POW(v)  = |v|
//  ROOT(x)  = x
//  #  = max
//  DIFF(x,y)   = y
//
// By default the Euclidean norm is assumed.  To change the norm,
// uncomment the appropriate set of macros below.
//
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// Use the following for the Euclidean norm
//----------------------------------------------------------------------
#define ANN_POW(v)  ((v)*(v))
#define ANN_ROOT(x)  sqrt(x)
#define ANN_SUM(x,y)  ((x) + (y))
#define ANN_DIFF(x,y)  ((y) - (x))

//----------------------------------------------------------------------
// Use the following for the L_1 (Manhattan) norm
//----------------------------------------------------------------------
// #define ANN_POW(v)  fabs(v)
// #define ANN_ROOT(x)  (x)
// #define ANN_SUM(x,y)  ((x) + (y))
// #define ANN_DIFF(x,y) ((y) - (x))

//----------------------------------------------------------------------
// Use the following for a general L_p norm
//----------------------------------------------------------------------
// #define ANN_POW(v)  pow(fabs(v),p)
// #define ANN_ROOT(x)  pow(fabs(x),1/p)
// #define ANN_SUM(x,y)  ((x) + (y))
// #define ANN_DIFF(x,y) ((y) - (x))

//----------------------------------------------------------------------
// Use the following for the L_infinity (Max) norm
//----------------------------------------------------------------------
// #define ANN_POW(v)  fabs(v)
// #define ANN_ROOT(x)  (x)
// #define ANN_SUM(x,y)  ((x) > (y) ? (x) : (y))
// #define ANN_DIFF(x,y) (y)

//----------------------------------------------------------------------
//  Array types
//
//  ANNpoint:
// A point is represented as a (dimensionless) vector of
// coordinates, that is, as a pointer to ANNcoord.  It is the
// user's responsibility to be sure that each such vector has
// been allocated with enough components.  Because only
// pointers are stored, the values should not be altered
// through the lifetime of the nearest neighbor data structure.
//  ANNpointArray is a dimensionless array of ANNpoint.
//  ANNdistArray is a dimensionless array of ANNdist.
//  ANNidxArray is a dimensionless array of ANNidx.  This is used for
// storing buckets of points in the search trees, and for returning
// the results of k nearest neighbor queries.
//----------------------------------------------------------------------

typedef ANNcoord *ANNpoint;  // a point
typedef ANNpoint *ANNpointArray; // an array of points
typedef ANNdist  *ANNdistArray;  // an array of distances
typedef ANNidx  *ANNidxArray;  // an array of point indices

//----------------------------------------------------------------------
//  Point operations:
//
// annDist() computes the (squared) distance between a pair
// of points.  Distance calculations for queries are NOT
// performed using this routine (for reasons of efficiency).
//
// Because points (somewhat like strings in C) are stored
// as pointers.  Consequently, creating and destroying
// copies of points may require storage allocation.  These
// procedures do this.
//
// annAllocPt() and annDeallocPt() allocate a deallocate
// storage for a single point, and return a pointer to it.
// The argument to AllocPt() is used to initialize all
// components.
//
// annAllocPts() allocates an array of points as well a place
// to store their coordinates, and initializes the points to
// point to their respective coordinates.  It allocates point
// storage in a contiguous block large enough to store all the
// points.  It performs no initialization.
//
//  annDeallocPt() deallocates a point allocated by annAllocPt().
//  annDeallocPts() deallocates points allocated by annAllocPts().
//
// annCopyPt() allocates space and makes a copy of a given point.
//----------------------------------------------------------------------

ANNdist annDist(
  int   dim,  // dimension of space
  ANNpoint  p,  // points
  ANNpoint  q);

ANNpoint annAllocPt(
  int   dim,  // dimension
  ANNcoord  c = 0);  // coordinate value (all equal)

ANNpointArray annAllocPts(
  int   n,  // number of points
  int   dim);  // dimension

void annDeallocPt(
  ANNpoint  &p);  // deallocate 1 point

void annDeallocPts(
  ANNpointArray &pa);  // point array

ANNpoint annCopyPt(   // copy point
  int   dim,  // dimension
  ANNpoint  source); // point to copy

//----------------------------------------------------------------------
//  Generic approximate nearest neighbor search structure.
// ANN supports a few different data structures for
// approximate nearest neighbor searching.  All these
// data structures at a minimum support single and k-nearest
// neighbor queries described here.  The nearest neighbor
// query returns an integer identifier and the distance
// to the nearest neighbor(s).
//----------------------------------------------------------------------
class ANNpointSet
{
public:
  virtual ~ANNpointSet()
  {}  // virtual distroyer

  virtual void annkSearch(  // approx k near neighbor search
    ANNpoint q,  // query point
    int  k,  // number of near neighbors to return
    ANNidxArray nn_idx,  // nearest neighbor array (returned)
    ANNdistArray dd,  // dist to near neighbors (returned)
    double  eps=0.0  // error bound
  ) = 0;    // pure virtual (defined elsewhere)
};

//----------------------------------------------------------------------
//  Brute-force nearest neighbor search:
// For the sake of comparison and validation we include a
// trivial brute-force nearest neighbor algorithm.  It is
// simple but inefficient.  The value of epsilon is ignored
// when performing nearest neighbor calculations, since all
// are done exactly.  This data structure is very slow, and
// should not be used for most applications.  It was implemented
// by the authors to help validate the more complex data
// structures below.
//
// The nearest neighbor index returned is the index in the
// array pa[] which is passed to the constructor.
//----------------------------------------------------------------------

class ANNbruteForce: public ANNpointSet
{
  int   dim;  // dimension
  int   n_pts;  // number of points
  ANNpointArray pts;  // point array
public:
  ANNbruteForce(   // constructor from point array
    ANNpointArray pa,  // point array
    int n,    // number of points
    int dd);   // dimension

  ~ANNbruteForce();   // destructor

  virtual void annkSearch(  // approx k near neighbor search
    ANNpoint q,  // query point
    int  k,  // number of near neighbors to return
    ANNidxArray nn_idx,  // nearest neighbor array (returned)
    ANNdistArray dd,  // dist to near neighbors (returned)
    double  eps=0.0); // error bound
};

//----------------------------------------------------------------------
//  kd-tree:
// The basic search data structure supported by ann is a kd-tree.
// The tree basically consists of a root pointer.  We also store
// the dimension of the space (since it is needed for many routines
// that access the structure).  The number of data points and the
// bucket size are stored, but they are mostly information items,
// and do not affect the data structure's function.  We also store
// the bounding box for the point set.
//
// kd-trees support two searching algorithms, standard search
// (which searches nodes in tree traversal order) and priority
// search (which searches nodes in order of increasing distance
// from the query point).  For many distributions the standard
// search seems to work just fine, but priority search is safer
// for worst-case performance.
//
// The nearest neighbor index returned is the index in the
// array pa[] which is passed to the constructor.
//
// There are two methods provided for printing the tree.  Print()
// is used to produce a "human-readable" display of the tree, with
// indenation, which is handy for debugging.  Dump() produces a
// format that is suitable reading by a program.  Finally the
// method getStats() collects statistics information on the tree
// (its size, height, etc.)  See ANNperf.h for information on the
// stats structure it returns.
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// Some types and objects used by kd-tree functions
// See kd_tree.h and kd_tree.cc for definitions
//----------------------------------------------------------------------
class ANNkdStats;   // stats on kd-tree
class ANNkd_node;   // generic node in a kd-tree
typedef ANNkd_node *ANNkd_ptr; // pointer to a kd-tree node

//----------------------------------------------------------------------
// kd-tree splitting rules
// kd-trees supports a collection of different splitting rules.
// In addition to the standard kd-tree splitting rule proposed
// by Friedman, Bentley, and Finkel, we have introduced a
// number of other splitting rules, which seem to perform
// as well or better (for the distributions we have tested).
//
// The splitting methods given below allow the user to tailor
// the data structure to the particular data set.  They are
// are described in greater details in the kd_split.cc source
// file.  The method ANN_KD_SUGGEST is the method chosen (rather
// subjectively) by the implementors as the one giving the
// fastest performance, and is the default splitting method.
//----------------------------------------------------------------------

enum ANNsplitRule
{
  ANN_KD_STD,   // the optimized kd-splitting rule
  ANN_KD_MIDPT,   // midpoint split
  ANN_KD_FAIR,   // fair split
  ANN_KD_SL_MIDPT,  // sliding midpoint splitting method
  ANN_KD_SL_FAIR,   // sliding fair split method
  ANN_KD_SUGGEST
}
;  // the authors' suggestion for best

class ANNkd_tree: public ANNpointSet
{
protected:
  int   dim;  // dimension of space
  int   n_pts;  // number of points in tree
  int   bkt_size; // bucket size
  ANNpointArray pts;  // the points
  ANNidxArray  pidx;  // point indices (to pts)
  ANNkd_ptr  root;  // root of kd-tree
  ANNpoint  bnd_box_lo; // bounding box low point
  ANNpoint  bnd_box_hi; // bounding box high point

  void SkeletonTree(   // construct skeleton tree
    int n,    // number of points
    int dd,    // dimension
    int bs);   // bucket size

public:
  ANNkd_tree(    // build skeleton tree
    int  n,  // number of points
    int  dd,  // dimension
    int  bs = 1); // bucket size

  ANNkd_tree(    // build from point array
    ANNpointArray pa,  // point array
    int  n,  // number of points
    int  dd,  // dimension
    int  bs = 1,  // bucket size
    ANNsplitRule split = ANN_KD_SUGGEST); // splitting method

  ~ANNkd_tree();   // tree destructor

  virtual void annkSearch(  // approx k near neighbor search
    ANNpoint q,  // query point
    int  k,  // number of near neighbors to return
    ANNidxArray nn_idx,  // nearest neighbor array (returned)
    ANNdistArray dd,  // dist to near neighbors (returned)
    double  eps=0.0); // error bound

  virtual void annkPriSearch(  // priority k near neighbor search
    ANNpoint q,  // query point
    int  k,  // number of near neighbors to return
    ANNidxArray nn_idx,  // nearest neighbor array (returned)
    ANNdistArray dd,  // dist to near neighbors (returned)
    double  eps=0.0); // error bound

  virtual void getStats(  // compute tree statistics
    ANNkdStats &st);  // the statistics (returned)

  virtual void Print(   // print the tree (for debugging)
    ANNbool with_pts,  // print points as well?
    ostream &out);   // output stream

  virtual void Dump(   // dump entire tree
    ANNbool with_pts,  // print points as well?
    ostream &out);   // output stream
};

//----------------------------------------------------------------------
//  Box decomposition tree (bd-tree)
// The bd-tree is inherited from a kd-tree.  The main difference
// in the bd-tree and the kd-tree is a new type of internal node
// called a shrinking node (in the kd-tree there is only one type
// of internal node, a splitting node).  The shrinking node
// makes it possible to generate balanced trees in which the
// cells have bounded aspect ratio.
//
// As with splitting rules, there are a number of different
// shrinking rules.  The shrinking rule ANN_BD_NONE does no
// shrinking (and hence produces a kd-tree tree).  The rule
// ANN_BD_SUGGEST uses the implementors favorite rule.
//----------------------------------------------------------------------

enum ANNshrinkRule
{
  ANN_BD_NONE,   // no shrinking at all (just kd-tree)
  ANN_BD_SIMPLE,   // simple splitting
  ANN_BD_CENTROID,  // centroid splitting
  ANN_BD_SUGGEST
}
;  // the authors' suggested choice

class ANNbd_tree: public ANNkd_tree
{
public:
  ANNbd_tree(    // build skeleton tree
    int  n,  // number of points
    int  dd,  // dimension
    int  bs = 1)  // bucket size
      :
      ANNkd_tree(n, dd, bs)
  {} // build base kd-tree

  ANNbd_tree(    // build from point array
    ANNpointArray pa,  // point array
    int  n,  // number of points
    int  dd,  // dimension
    int  bs = 1,  // bucket size
    ANNsplitRule split  = ANN_KD_SUGGEST, // splitting rule
    ANNshrinkRule shrink = ANN_BD_SUGGEST); // shrinking rule
};

//----------------------------------------------------------------------
//  Other functions
//----------------------------------------------------------------------

void annMaxPtsVisit(   // limit max. pts to visit in search
  int   maxPts); // the limit

#endif
