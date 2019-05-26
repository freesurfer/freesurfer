#pragma once
/**
 * @file  mrishash.h
 * @brief Implements a hash table mechanism to speed comparing vertices
 *
 * The purpose of MRI hash tables is to vastly accelerate algorithms which 
 * need to compare vertices with one another or to a point.  See: 
 * http://wideman-one.com/gw/brain/fs/2007/mrishash/mrishash_100_overview.htm
 */
/*
 * Original Author: Graham Wideman, based on code by Bruce Fischl
 * CVS Revision Info:
 *    $Author: zkaufman $
 *    $Date: 2015/03/18 17:04:00 $
 *    $Revision: 1.26 $
 *
 * Copyright © 2011 The General Hospital Corporation (Boston, MA) "MGH"
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

#include "mrisurf_aaa.h"
#include "error.h"

// define the following to get a single inclusion of non-renamed functions
// Bevin thinks this is no longer used anywhere
//
#define MRISHASH_VANILLA_FUNCS


// Version
//
int MHT_gw_version(void);           // version of that unit


// Test functions
//
void MHTfindReportCounts(int * BucketsChecked, 
                         int * BucketsPresent, 
                         int * VtxNumByMHT);

// Tuning
//
void MHT_maybeParallel_begin();     // Note: Can be nested!
void MHT_maybeParallel_end();


// Support multiple representations
//
#define MHT_VIRTUAL                 
#define MHT_ABSTRACT                
#define MHT_STATIC_MEMBER           
#define MHT_FUNCTION(NAME)          MHT##NAME
#define MHT_CONST_THIS_PARAMETER    MRIS_HASH_TABLE const *mht,
#define MHT_CONST_THIS              
#define MHT_THIS_PARAMETER_NOCOMMA  MRIS_HASH_TABLE *mht 
#define MHT_MRIS_PARAMETER_NOCOMMA  MRIS *mris
#define MHT_MRIS_PARAMETER  MHT_MRIS_PARAMETER_NOCOMMA ,
#define MHT_THIS_PARAMETER  MHT_THIS_PARAMETER_NOCOMMA ,
#include "mrishash_traditional_functions.h"
