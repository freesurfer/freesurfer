/**
 * @file  gcamorphtestutils.cpp
 * @brief Utilities to help with testing GCAmorph routines
 *
 * Reference:
 * "How to Stay Sane while Programming: Collected Wisdom from Broadmoor"
 */
/*
 * Original Author: Richard Edgar
 * CVS Revision Info:
 *    $Author: rge21 $
 *    $Date: 2010/03/18 15:09:41 $
 *    $Revision: 1.17 $
 *
 * Copyright (C) 2002-2008,
 * The General Hospital Corporation (Boston, MA). 
 * All rights reserved.
 *
 * Distribution, usage and copying of this software is covered under the
 * terms found in the License Agreement file named 'COPYING' found in the
 * FreeSurfer source code root directory, and duplicated here:
 * https://surfer.nmr.mgh.harvard.edu/fswiki/FreeSurferOpenSourceLicense
 *
 * General inquiries: freesurfer@nmr.mgh.harvard.edu
 * Bug reports: analysis-bugs@nmr.mgh.harvard.edu
 *
 */

#include <vector>
#include <iostream>
using namespace std;

#include "chronometer.hpp"

#include "gcamorphtestutils.hpp"


#include "gcamorphtestutils.h"


// ======================================================================


GCAMforCMPutils::GCAMforCMPutils( void ) : varTypeMap() {
    
  // Sanity check
  if( this->nVars != 8 ) {
    cerr << __FUNCTION__ << ": Invalid nVars!" << endl;
    exit( EXIT_FAILURE );
  }
  if( this->nDims != 3 ) {
    cerr << __FUNCTION__ << ": Invalid nDims!" << endl;
    exit( EXIT_FAILURE );
  }
  
  varTypeMap[ "rx" ] = NC_DOUBLE;
  varTypeMap[ "ry" ] = NC_DOUBLE;
  varTypeMap[ "rz" ] = NC_DOUBLE;

  varTypeMap[ "origArea" ] = NC_FLOAT;
  varTypeMap[ "area" ] = NC_FLOAT;
  varTypeMap[ "area1" ] = NC_FLOAT;
  varTypeMap[ "area2" ] = NC_FLOAT;
  
  varTypeMap[ "invalid" ] = NC_CHAR;
  
  // And another sanity check
  if( varTypeMap.size() != this->nVars ) {
    cerr << __FUNCTION__ << ": Incorrect entries in varTypeMap" << endl;
    exit( EXIT_FAILURE );
  }
}


// -------------------


void GCAMforCMPutils::Write( const GCAM* src, string fName ) const {

  // Construct the filename, using fName passed by value
  fName += ".nc";
  
  std::cout << __FUNCTION__ << ": Writing file " << fName << " ... ";
  std::cout.flush();
  
  // Reference for the file
  int ncid;
  
  // Open the file
  NC_SAFE_CALL( nc_create( fName.c_str(), NC_CLOBBER, &ncid ) );
  
  
  // Set up the dimensions
  int dimIDs[nDims];
  NC_SAFE_CALL( nc_def_dim( ncid, "x", src->width, &dimIDs[this->iX] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "y", src->height, &dimIDs[this->iY] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "z", src->depth, &dimIDs[this->iZ] ) );
  
  // Set up the variable IDs
  map<string,int> varIDmap;
  map<string,nc_type>::const_iterator myIt;
  
  for( myIt = this->varTypeMap.begin();
       myIt != this->varTypeMap.end();
       myIt++ ) {
    NC_SAFE_CALL( nc_def_var( ncid,
			      myIt->first.c_str(), // Name of the variable
			      myIt->second,        // Type of the variable
			      nDims, dimIDs,
			      &varIDmap[ myIt->first ] ) );
  }
  
  // Sanity check
  if( varTypeMap.size() != varIDmap.size() ) {
    cerr << __FUNCTION__ << ": Failed to create varIDmap correctly" << endl;
      exit( EXIT_FAILURE );
  }
  
  // Make the end of the 'definition' region
  NC_SAFE_CALL( nc_enddef( ncid ) );
  
  // Pack into contiguous arrays
  const size_t nElems = src->width * src->height * src->depth;
  
  vector<double> x( nElems ), y( nElems ), z( nElems );
  vector<float> area(nElems ), area1( nElems ), area2( nElems );
  vector<float> origArea( nElems );
  vector<char> invalid( nElems );

  // Ugly loop to do the writing element by element
  for( int i=0; i<src->width; i++ ) {
    for( int j=0; j<src->height; j++ ) {
      for( int k=0; k<src->depth; k++ ) {
	const GCA_MORPH_NODE& gcamn = src->nodes[i][j][k];
	
	const size_t i1d = k + ( src->depth * ( j + ( src->height * i ) ) );

	x.at(i1d) = gcamn.x;
	y.at(i1d) = gcamn.y;
	z.at(i1d) = gcamn.z;

	area.at(i1d) = gcamn.area;
	area1.at(i1d) = gcamn.area1;
	area2.at(i1d) = gcamn.area2;
	origArea.at(i1d) = gcamn.orig_area;

	invalid.at(i1d) = gcamn.invalid;
      }
    }
  }

  // We use 'find' to get an exception if the name doesn't exist

  NC_SAFE_CALL( nc_put_var_double( ncid,
				   varIDmap.find( "rx" )->second,
				   &x[0] ) );
  NC_SAFE_CALL( nc_put_var_double( ncid,
				   varIDmap.find( "ry" )->second,
				   &y[0] ) );
  NC_SAFE_CALL( nc_put_var_double( ncid,
				   varIDmap.find( "rz" )->second,
				   &z[0] ) );

  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "origArea" )->second,
				  &origArea[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "area" )->second,
				  &area[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "area1" )->second,
				  &area1[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "area2" )->second,
				  &area2[0] ) );

  NC_SAFE_CALL( nc_put_var_text( ncid,
				 varIDmap.find( "invalid" )->second,
				 &invalid[0] ) );


  
  
  // Close the file
  NC_SAFE_CALL( nc_close( ncid ) );

  cout << "complete" << endl;
}



// ----------------

void GCAMforCMPutils::Read( GCAM** dst, string fName ) const {
  
  // Make sure input pointer is NULL
  if( *dst != NULL ) {
    cerr << __FUNCTION__
	   << ": dst pointer not NULL!"
	 << endl;
    exit( EXIT_FAILURE );
  }
  
  // Construct the filename, using fName passed by value
  fName += ".nc";
  
  std::cout << __FUNCTION__ << ": Reading file " << fName << " ... ";
  std::cout.flush();
  
  // Reference for the file
  int ncid;
  
  // Open the file
  NC_SAFE_CALL( nc_open( fName.c_str(), NC_NOWRITE, &ncid ) );
  
  // Query number of variables and dimensions
  int nDimFile, nVarFile;
  NC_SAFE_CALL( nc_inq_ndims( ncid, &nDimFile ) );
  NC_SAFE_CALL( nc_inq_nvars( ncid, &nVarFile ) );
  
  if( nDimFile != static_cast<int>(this->nDims) ) {
    cerr << "Invalid number of dimensions " << nDimFile << endl;
    exit( EXIT_FAILURE );
  }
  
  if( nVarFile != static_cast<int>(this->nVars) ) {
    std::cerr << "Invalid number of variables " << nVarFile << endl;
  }
  
  int dimIDs[3];
  size_t dimLen[3];
  
  NC_SAFE_CALL( nc_inq_dimid( ncid, "x", &dimIDs[this->iX] ) );
  NC_SAFE_CALL( nc_inq_dimid( ncid, "y", &dimIDs[this->iY] ) );
  NC_SAFE_CALL( nc_inq_dimid( ncid, "z", &dimIDs[this->iZ] ) );
  
  NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[this->iX], &dimLen[this->iX] ) );
  NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[this->iY], &dimLen[this->iY] ) );
  NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[this->iZ], &dimLen[this->iZ] ) );
  
  // Allocate the target
  *dst = GCAMalloc( dimLen[this->iX], dimLen[this->iY], dimLen[this->iZ] );
  if( *dst == NULL ) {
    std::cerr << __FUNCTION__
	      << ": GCAMalloc failed!" << std::endl;
      exit( EXIT_FAILURE );
  }
  
  // Get hold of the variable IDs
  map<string,int> varIDmap;
  map<string,nc_type>::const_iterator myIt;
  
  for( myIt = this->varTypeMap.begin();
       myIt != this->varTypeMap.end();
	 myIt++ ) {
    NC_SAFE_CALL( nc_inq_varid( ncid,
				myIt->first.c_str(), // Name of the variable
				&varIDmap[ myIt->first ] ) );
  }
  
  // Read into contiguous arrays
  const size_t nElems = (*dst)->width * (*dst)->height * (*dst)->depth;
  
  vector<double> x( nElems ), y( nElems ), z( nElems );
  vector<float> area(nElems ), area1( nElems ), area2( nElems );
  vector<float> origArea( nElems );
  vector<char> invalid( nElems );

  // We use 'find' to get an exception if the name doesn't exist

  NC_SAFE_CALL( nc_get_var_double( ncid,
				   varIDmap.find( "rx" )->second,
				   &x[0] ) );
  NC_SAFE_CALL( nc_get_var_double( ncid,
				   varIDmap.find( "ry" )->second,
				   &y[0] ) );
  NC_SAFE_CALL( nc_get_var_double( ncid,
				   varIDmap.find( "rz" )->second,
				   &z[0] ) );

  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "origArea" )->second,
				  &origArea[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "area" )->second,
				  &area[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "area1" )->second,
				  &area1[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "area2" )->second,
				  &area2[0] ) );

  NC_SAFE_CALL( nc_get_var_text( ncid,
				 varIDmap.find( "invalid" )->second,
				 &invalid[0] ) );


  // Unpack into the GCA_MORPH structure
  for( int i=0; i<(*dst)->width; i++ ) {
    for( int j=0; j<(*dst)->height; j++ ) {
      for( int k=0; k<(*dst)->depth; k++ ) {
	// Save some typing
	GCA_MORPH_NODE& gcamn = (*dst)->nodes[i][j][k];
	// Compute 1D index
	const size_t i1d = k + ( (*dst)->depth * ( j + ( (*dst)->height * i ) ) );

	gcamn.x = x.at(i1d);
	gcamn.y = y.at(i1d);
	gcamn.z = z.at(i1d);

	gcamn.orig_area = origArea.at(i1d);
	gcamn.area = area.at(i1d);
	gcamn.area1 = area1.at(i1d);
	gcamn.area2 = area2.at(i1d);

	gcamn.invalid = invalid.at(i1d);
      }
    }
  }
  
  
  NC_SAFE_CALL( nc_close( ncid ) );
  
  cout << "complete" << endl;
}




// ======================================================================


GCAMorphUtils::GCAMorphUtils( void ) : varTypeMap() {
    
  // Sanity check
  if( this->nVars != 14 ) {
    cerr << __FUNCTION__ << ": Invalid nVars!" << endl;
    exit( EXIT_FAILURE );
  }
  if( this->nDims != 3 ) {
    cerr << __FUNCTION__ << ": Invalid nDims!" << endl;
    exit( EXIT_FAILURE );
  }
  
  varTypeMap[ "rx" ] = NC_DOUBLE;
  varTypeMap[ "ry" ] = NC_DOUBLE;
  varTypeMap[ "rz" ] = NC_DOUBLE;

  varTypeMap[ "origArea" ] = NC_FLOAT;
  varTypeMap[ "origArea1" ] = NC_FLOAT;
  varTypeMap[ "origArea2" ] = NC_FLOAT;
  varTypeMap[ "area" ] = NC_FLOAT;
  varTypeMap[ "area1" ] = NC_FLOAT;
  varTypeMap[ "area2" ] = NC_FLOAT;
  
  varTypeMap[ "invalid" ] = NC_CHAR;
  varTypeMap[ "label" ] = NC_INT;
  varTypeMap[ "status" ] = NC_INT;

  varTypeMap[ "mean" ] = NC_FLOAT;
  varTypeMap[ "variance" ] = NC_FLOAT;
  
  
  // And another sanity check
  if( varTypeMap.size() != this->nVars ) {
    cerr << __FUNCTION__ << ": Incorrect entries in varTypeMap" << endl;
    exit( EXIT_FAILURE );
  }
}




void GCAMorphUtils::Write( const GCAM* src, string fName ) const {
  
  // Construct the filename, using fName passed by value
  fName += ".nc";
  
  std::cout << __FUNCTION__ << ": Writing file " << fName << " ... ";
  std::cout.flush();
  
  // Sanity check the GCAM
  if( src->ninputs != 1 ) {
    std::cerr << __FUNCTION__
	      << ": Must have ninputs=1" << std::endl;
  }

  // Reference for the file
  int ncid;
  
  // Open the file
  NC_SAFE_CALL( nc_create( fName.c_str(), NC_CLOBBER, &ncid ) );
  
  
  // Set up the dimensions
  int dimIDs[nDims];
  NC_SAFE_CALL( nc_def_dim( ncid, "x", src->width, &dimIDs[this->iX] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "y", src->height, &dimIDs[this->iY] ) );
  NC_SAFE_CALL( nc_def_dim( ncid, "z", src->depth, &dimIDs[this->iZ] ) );
  
  // Set up the variable IDs
  map<string,int> varIDmap;
  map<string,nc_type>::const_iterator myIt;
  
  for( myIt = this->varTypeMap.begin();
       myIt != this->varTypeMap.end();
       myIt++ ) {
    NC_SAFE_CALL( nc_def_var( ncid,
			      myIt->first.c_str(), // Name of the variable
			      myIt->second,        // Type of the variable
			      nDims, dimIDs,
			      &varIDmap[ myIt->first ] ) );
  }
  
  // Sanity check
  if( varTypeMap.size() != varIDmap.size() ) {
    cerr << __FUNCTION__ << ": Failed to create varIDmap correctly" << endl;
      exit( EXIT_FAILURE );
  }
  
  // Make the end of the 'definition' region
  NC_SAFE_CALL( nc_enddef( ncid ) );
  
  // Pack into contiguous arrays
  const size_t nElems = src->width * src->height * src->depth;
  
  vector<double> x( nElems ), y( nElems ), z( nElems );
  vector<float> area(nElems ), area1( nElems ), area2( nElems );
  vector<float> origArea( nElems ), origArea1( nElems ), origArea2( nElems );
  vector<char> invalid( nElems );
  vector<int> label( nElems ), status( nElems );
  vector<float> mean( nElems ), variance( nElems );

  // Ugly loop to do the writing element by element
  for( int i=0; i<src->width; i++ ) {
    for( int j=0; j<src->height; j++ ) {
      for( int k=0; k<src->depth; k++ ) {
	const GCA_MORPH_NODE& gcamn = src->nodes[i][j][k];
	
	const size_t i1d = k + ( src->depth * ( j + ( src->height * i ) ) );

	x.at(i1d) = gcamn.x;
	y.at(i1d) = gcamn.y;
	z.at(i1d) = gcamn.z;

	area.at(i1d) = gcamn.area;
	area1.at(i1d) = gcamn.area1;
	area2.at(i1d) = gcamn.area2;

	origArea.at(i1d) = gcamn.orig_area;
	origArea1.at(i1d) = gcamn.orig_area1;
	origArea2.at(i1d) = gcamn.orig_area2;

	invalid.at(i1d) = gcamn.invalid;
	label.at(i1d) = gcamn.label;
	status.at(i1d) = gcamn.status;

	if( gcamn.gc != NULL ) {
	  mean.at(i1d) = gcamn.gc->means[0];
	  variance.at(i1d) = gcamn.gc->covars[0];
	} else {
	  mean.at(i1d) = -1;
	  variance.at(i1d) = -1;
	}
      }
    }
  }

  // We use 'find' to get an exception if the name doesn't exist

  NC_SAFE_CALL( nc_put_var_double( ncid,
				   varIDmap.find( "rx" )->second,
				   &x[0] ) );
  NC_SAFE_CALL( nc_put_var_double( ncid,
				   varIDmap.find( "ry" )->second,
				   &y[0] ) );
  NC_SAFE_CALL( nc_put_var_double( ncid,
				   varIDmap.find( "rz" )->second,
				   &z[0] ) );

  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "origArea" )->second,
				  &origArea[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "origArea1" )->second,
				  &origArea1[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "origArea2" )->second,
				  &origArea2[0] ) );

  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "area" )->second,
				  &area[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "area1" )->second,
				  &area1[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "area2" )->second,
				  &area2[0] ) );

  NC_SAFE_CALL( nc_put_var_text( ncid,
				 varIDmap.find( "invalid" )->second,
				 &invalid[0] ) );

  NC_SAFE_CALL( nc_put_var_int( ncid,
				varIDmap.find( "label" )->second,
				&label[0] ) );
  NC_SAFE_CALL( nc_put_var_int( ncid,
				varIDmap.find( "status" )->second,
				&status[0] ) );

  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "mean" )->second,
				  &mean[0] ) );
  NC_SAFE_CALL( nc_put_var_float( ncid,
				  varIDmap.find( "variance" )->second,
				  &variance[0] ) );

  
  
  // Close the file
  NC_SAFE_CALL( nc_close( ncid ) );

  cout << "complete" << endl;
}



// ---------------

void GCAMorphUtils::Read( GCAM** dst, string fName ) const {

  // Make sure input pointer is NULL
  if( *dst != NULL ) {
    cerr << __FUNCTION__
	   << ": dst pointer not NULL!"
	 << endl;
    exit( EXIT_FAILURE );
  }
  
  // Construct the filename, using fName passed by value
  fName += ".nc";
  
  std::cout << __FUNCTION__ << ": Reading file " << fName << " ... ";
  std::cout.flush();
  
  // Reference for the file
  int ncid;

  // Open the file
  NC_SAFE_CALL( nc_open( fName.c_str(), NC_NOWRITE, &ncid ) );
  

  // Query number of variables and dimensions
  int nDimFile, nVarFile;
  NC_SAFE_CALL( nc_inq_ndims( ncid, &nDimFile ) );
  NC_SAFE_CALL( nc_inq_nvars( ncid, &nVarFile ) );
  
  if( nDimFile != static_cast<int>(this->nDims) ) {
    cerr << "Invalid number of dimensions " << nDimFile << endl;
    exit( EXIT_FAILURE );
  }
  
  if( nVarFile != static_cast<int>(this->nVars) ) {
    std::cerr << "Invalid number of variables " << nVarFile << endl;
  }
  
  int dimIDs[3];
  size_t dimLen[3];
  
  NC_SAFE_CALL( nc_inq_dimid( ncid, "x", &dimIDs[this->iX] ) );
  NC_SAFE_CALL( nc_inq_dimid( ncid, "y", &dimIDs[this->iY] ) );
  NC_SAFE_CALL( nc_inq_dimid( ncid, "z", &dimIDs[this->iZ] ) );
  
  NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[this->iX], &dimLen[this->iX] ) );
  NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[this->iY], &dimLen[this->iY] ) );
  NC_SAFE_CALL( nc_inq_dimlen( ncid, dimIDs[this->iZ], &dimLen[this->iZ] ) );
  
  // Allocate the target
  *dst = GCAMalloc( dimLen[this->iX], dimLen[this->iY], dimLen[this->iZ] );
  if( *dst == NULL ) {
    std::cerr << __FUNCTION__
	      << ": GCAMalloc failed!" << std::endl;
      exit( EXIT_FAILURE );
  }

  // Set this to have only one input
  (*dst)->ninputs = 1;
  
  // Get hold of the variable IDs
  map<string,int> varIDmap;
  map<string,nc_type>::const_iterator myIt;
  
  for( myIt = this->varTypeMap.begin();
       myIt != this->varTypeMap.end();
	 myIt++ ) {
    NC_SAFE_CALL( nc_inq_varid( ncid,
				myIt->first.c_str(), // Name of the variable
				&varIDmap[ myIt->first ] ) );
  }
  
  // Read into contiguous arrays
  const size_t nElems = (*dst)->width * (*dst)->height * (*dst)->depth;

  vector<double> x( nElems ), y( nElems ), z( nElems );
  vector<float> area(nElems ), area1( nElems ), area2( nElems );
  vector<float> origArea( nElems ), origArea1( nElems ), origArea2( nElems );
  vector<char> invalid( nElems );
  vector<int> label( nElems ), status( nElems );
  vector<float> mean( nElems ), variance( nElems );

  // We use 'find' to get an exception if the name doesn't exist

  NC_SAFE_CALL( nc_get_var_double( ncid,
				   varIDmap.find( "rx" )->second,
				   &x[0] ) );
  NC_SAFE_CALL( nc_get_var_double( ncid,
				   varIDmap.find( "ry" )->second,
				   &y[0] ) );
  NC_SAFE_CALL( nc_get_var_double( ncid,
				   varIDmap.find( "rz" )->second,
				   &z[0] ) );

  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "origArea" )->second,
				  &origArea[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "origArea1" )->second,
				  &origArea1[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "origArea2" )->second,
				  &origArea2[0] ) );

  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "area" )->second,
				  &area[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "area1" )->second,
				  &area1[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "area2" )->second,
				  &area2[0] ) );

  NC_SAFE_CALL( nc_get_var_text( ncid,
				 varIDmap.find( "invalid" )->second,
				 &invalid[0] ) );

  NC_SAFE_CALL( nc_get_var_int( ncid,
				varIDmap.find( "label" )->second,
				&label[0] ) );
  NC_SAFE_CALL( nc_get_var_int( ncid,
				varIDmap.find( "status" )->second,
				&status[0] ) );

  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "mean" )->second,
				  &mean[0] ) );
  NC_SAFE_CALL( nc_get_var_float( ncid,
				  varIDmap.find( "variance" )->second,
				  &variance[0] ) );

  // Unpack into the GCA_MORPH structure
  for( int i=0; i<(*dst)->width; i++ ) {
    for( int j=0; j<(*dst)->height; j++ ) {
      for( int k=0; k<(*dst)->depth; k++ ) {
	// Save some typing
	GCA_MORPH_NODE& gcamn = (*dst)->nodes[i][j][k];
	// Compute 1D index
	const size_t i1d = k + ( (*dst)->depth * ( j + ( (*dst)->height * i ) ) );

	gcamn.x = x.at(i1d);
	gcamn.y = y.at(i1d);
	gcamn.z = z.at(i1d);

	gcamn.orig_area = origArea.at(i1d);
	gcamn.orig_area1 = origArea1.at(i1d);
	gcamn.orig_area2 = origArea2.at(i1d);

	gcamn.area = area.at(i1d);
	gcamn.area1 = area1.at(i1d);
	gcamn.area2 = area2.at(i1d);

	gcamn.invalid = invalid.at(i1d);
	gcamn.label = label.at(i1d);
	gcamn.status = status.at(i1d);
	
	// Mean and variance require special handling
	if( variance.at(i1d) >= 0 ) {
	  gcamn.gc = alloc_gcs( 1, 0, 1 );
	  gcamn.gc->means[0] = mean.at(i1d);
	  gcamn.gc->covars[0] = variance.at(i1d);
	} else {
	  gcamn.gc = NULL;
	}
      }
    }
  }

   
  NC_SAFE_CALL( nc_close( ncid ) );
  
  cout << "complete" << endl;
}



// #########################################################################

static GCAMforCMPutils myCMPutils;
static GCAMorphUtils myGCAMutils;

// ======================================================================

void WriteGCAMforMetricProperties( const GCAM* src, const char* fName ) {
  /*!
    Writes the parts of a GCAM required for gcamComputeMetricProperties
    to a NetCDF file.
    Currently uses the C API, since the C++ one doesn't seem to be available
  */

  myCMPutils.Write( src, fName );
  
}



// ======================================================================

void ReadGCAMforMetricProperties( GCAM** dst, const char* fName ) {
  /*!
    Complement of WriteGCAMforMetricProperties
  */

  myCMPutils.Read( dst, fName );
}


// ======================================================================

void WriteGCAMoneInput( const GCAM* src, const char* fName ) {
  myGCAMutils.Write( src, fName );
}
