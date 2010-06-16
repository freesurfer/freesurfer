#include <cstdlib>
#include <cstdio>
#include <cmath>

#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <limits>
using namespace std;

#include <boost/program_options.hpp>
namespace bpo = boost::program_options;

#include "gcamorph.h"
#include "gcamorphtestutils.hpp"

#include "chronometer.hpp"
#include "kahan.hpp"

#ifdef FS_CUDA
#include "devicemanagement.h"
#endif



// ==========================================================

const int MATCH_RETURN = EXIT_SUCCESS;
const int DIFF_RETURN = 100;

enum compModes{ errorL2, maxAbsDifference };

const compModes compModeDefault = errorL2;
const float matchTolDefault = 1e-6;
const float diffTolDefault = 1e-4;


string varName;
string f1Name;
string f2Name;

compModes compMode = compModeDefault;
float matchTol;
float diffTol;

const char* Progname = "gcamCompare";

// ==========================================================

void ReadCommandLine( int ac, char* av[] ) {

  try {
    bpo::options_description generic("Generic options");
    generic.add_options()
      ("help", "Produce help message" )
      ("version,v", "Show source version" )
      ;

    bpo::options_description comparison("Comparision options");

    stringstream mss, dss;
    mss << "Tolerance for match (returns " << MATCH_RETURN << ")";
    dss << "Tolerance for 'diff' (returns " << DIFF_RETURN << ")";
    comparison.add_options()
      ("errL2", "Compute error in L2 norm" )
      ("maxAbsDiff", "Compute maximum absolute difference" )
      ("match,m", bpo::value<float>(&matchTol)->default_value(matchTolDefault),  mss.str().c_str() )
      ("diff,d", bpo::value<float>(&diffTol)->default_value(diffTolDefault), dss.str().c_str() )
      ;

    bpo::options_description hidden("Hidden options");
    hidden.add_options()
      ("var", bpo::value<string>(&varName), "Variable to compare" )
      ("file1", bpo::value<string>(&f1Name), "First file to compare" )
      ("file2", bpo::value<string>(&f2Name), "Second file to compare" )
      ;

    bpo::positional_options_description p;
    p.add( "var", 1 );
    p.add( "file1", 1 );
    p.add( "file2", 1 );

    bpo::options_description cmdLine;
    cmdLine.add(generic).add(comparison).add(hidden);

    bpo::options_description visible("Allowed options");
    visible.add(generic).add(comparison);


    bpo::variables_map vm;
    bpo::store( bpo::command_line_parser(ac, av).
		options(cmdLine).positional(p).run(), vm );
    bpo::notify( vm );

    if( vm.count( "errL2" ) ) {
      compMode = errorL2;
    } else if ( vm.count( "maxAbsDiff" ) ) {
      compMode = maxAbsDifference;
    }

    if( vm.count( "help" ) ) {
      cout << "Usage :" << endl;
      cout << "  gcamCompare <varName> <file1> <file2>" << endl;
      cout << visible << endl;
      exit( EXIT_SUCCESS );
    }
    if( vm.count( "version" ) ) {
      cout << __FILE__ << endl
	   << "$Author: rge21 $\n"
	   << "$Date: 2010/06/16 19:57:40 $\n"
	   << "$Revision: 1.2 $\n" 
	   << endl;
      exit( EXIT_SUCCESS );
    }
  }
  catch( exception& e ) {
    cerr << "Error: " << e.what() << endl;
    exit( EXIT_FAILURE );
  }
  catch( ... ) {
    cerr << "Unknown exception" << endl;
    exit( EXIT_FAILURE );
  }
}

// ===========================================================

//! Gets information about the given variable
void GetVariableInfo( const string file, const string var,
		      nc_type& vType, vector<size_t>& vDims ) {
  

  int ncid;

  NC_SAFE_CALL( nc_open( file.c_str(), NC_NOWRITE, &ncid ) );

  int vid;
  int ncerr = nc_inq_varid( ncid, var.c_str(), &vid );
  if( NC_NOERR != ncerr ) {
    cerr << "Variable " << var << " not found in file " << file << endl;
    exit( EXIT_FAILURE );
  }

  // Find out the type
  NC_SAFE_CALL( nc_inq_vartype( ncid, vid, &vType ) );

  // Get dimension information
  int ndims;
  NC_SAFE_CALL( nc_inq_varndims( ncid, vid, &ndims ) );
  if( ndims != 3 ) {
    cerr << "Variable " << var
	 << " in file " << file
	 << " has " << ndims << " dimensions" << endl;
    exit( EXIT_FAILURE );
  }

  
  int dimids[3];
  NC_SAFE_CALL( nc_inq_vardimid( ncid, vid, dimids ) ); 
  for( unsigned int i=0; i<3; i++ ) {
    NC_SAFE_CALL( nc_inq_dimlen( ncid, dimids[i], &vDims.at(i) ) );
  }

  NC_SAFE_CALL( nc_close( ncid ) );
}


// ===========================================================

//! Validates the two files against each other
void Validate( const string f1, const string f2, const string var,
	       nc_type& vType ) {

  nc_type f1type, f2type;
  vector<size_t> f1dims(3), f2dims(3);

  GetVariableInfo( f1, var, f1type, f1dims );
  GetVariableInfo( f2, var, f2type, f2dims );
  
  if( f1type != f2type ) {
    cerr << "Variable has different type in each file!" << endl;
    exit( EXIT_FAILURE );
  }

  for( unsigned int i=0; i<3; i++ ) {
    if( f1dims.at(i) != f2dims.at(i) ) {
      cerr << "Variables differ in dimension " << i << endl;
      exit( EXIT_FAILURE );
    }
  }

  // Assign output
  vType = f1type;

  cout << "Variable type is " << vType << endl;
  cout << "Dimensions are "
       << "( "
       << f1dims.at(0) << ", "
       << f1dims.at(1) << ", "
       << f1dims.at(2)
       << ")" << endl;
}


// ===========================================================

//! Wrapper around nc_get_var_<type>
template<typename T>
void ReadArray( const int ncid, const int vid,
		vector<T>& var ) {
  cerr << __PRETTY_FUNCTION__ << ": Unsupported type" << endl;
  exit( EXIT_FAILURE );
}

template<>
void ReadArray<double>( const int ncid, const int vid,
			vector<double>& var ) {
  NC_SAFE_CALL( nc_get_var_double( ncid, vid, &var[0] ) );
}

template<>
void ReadArray<float>( const int ncid, const int vid,
			vector<float>& var ) {
  NC_SAFE_CALL( nc_get_var_float( ncid, vid, &var[0] ) );
}

template<>
void ReadArray<char>( const int ncid, const int vid,
		      vector<char>& var ) {
  NC_SAFE_CALL( nc_get_var_text( ncid, vid, &var[0] ) );
}

template<>
void ReadArray<int>( const int ncid, const int vid,
		     vector<int>& var ) {
  NC_SAFE_CALL( nc_get_var_int( ncid, vid, &var[0] ) );
}


//! Read a variable from the file
template<typename T>
void ReadVariable( const string fName, const string varName,
		   vector<T>& var ) {
  /*!
    Read in the named variable from the given file.
    Make sure the template parameter matches
    the variable type!
  */

  // Open the file
  int ncid;
  NC_SAFE_CALL( nc_open( fName.c_str(), NC_NOWRITE, &ncid ) );

  // Get the variable ID
  int vid;
  NC_SAFE_CALL( nc_inq_varid( ncid, varName.c_str(), &vid ) );

  // Get the dimensions
  int dimids[3];
  size_t dimLens[3];
  NC_SAFE_CALL( nc_inq_vardimid( ncid, vid, dimids ) ); 
  for( unsigned int i=0; i<3; i++ ) {
    NC_SAFE_CALL( nc_inq_dimlen( ncid, dimids[i], &dimLens[i] ) );
  }

  // Resize the vector
  var.resize( dimLens[0] * dimLens[1] * dimLens[2] );

  // Read the data
  ReadArray( ncid, vid, var );

  // Close the file
  NC_SAFE_CALL( nc_close( ncid ) );
}




// ===========================================================

//! Compute the L2 norm of two vectors
template<typename T>
T ComputeL2err( const vector<T>& v1, const vector<T>& v2 ) {
  
  KahanSum<T> diff, ref;
  T l2err;

  // Check sizes
  if( v1.size() != v2.size() ) {
    cerr << __FUNCTION__ << ": Sizes not equal!" << endl;
    exit( EXIT_FAILURE );
  }

  // Perform the sums
  for( unsigned int i=0; i<v1.size(); i++ ) {
    T tmp = v1[i] - v2[i];
    diff += (tmp*tmp);
    ref += (v1[i] * v1[i]);
  }

  if( ref.GetSum() != 0 ) {
    l2err = sqrt( diff.GetSum() / ref.GetSum() );
  } else if( diff.GetSum() == 0 ) {
    l2err = 0;
  } else {
    l2err = numeric_limits<T>::max();
  }
  return( l2err );
}


template<>
char ComputeL2err( const vector<char>& v1, const vector<char>& v2 ) {
  // Check sizes
  if( v1.size() != v2.size() ) {
    cerr << __FUNCTION__ << ": Sizes not equal!" << endl;
    exit( EXIT_FAILURE );
  }
  cerr << __FUNCTION__ << ": L2 error not supported for char" << endl;
  exit( EXIT_FAILURE );
}



template<>
int ComputeL2err( const vector<int>& v1, const vector<int>& v2 ) {
  // Check sizes
  if( v1.size() != v2.size() ) {
    cerr << __FUNCTION__ << ": Sizes not equal!" << endl;
    exit( EXIT_FAILURE );
  }
  cerr << __FUNCTION__ << ": L2 error not supported for int" << endl;
  exit( EXIT_FAILURE );
}


//! Compute maximum absolute difference
template<typename T>
T ComputeMaxAbsDiff( const vector<T>& v1, const vector<T>& v2 ) {

   // Check sizes
  if( v1.size() != v2.size() ) {
    cerr << __FUNCTION__ << ": Sizes not equal!" << endl;
    exit( EXIT_FAILURE );
  }

  T myMaxDiff = 0;

  for( unsigned int i=0; i<v1.size(); i++ ) {
    T tmp;
    if( v1[i] > v2[i] ) {
      tmp = v1[i] - v2[i];
    } else {
      tmp = v2[i] - v1[i];
    }
    if( tmp > myMaxDiff ) {
      myMaxDiff = tmp;
    }
  }

  return( static_cast<T>(myMaxDiff) );
}

// ===========================================================

//! Wrapper to let us print out the value of a char, not the character
template<typename T>
string GetString( const T val ) {
  stringstream vs;

  vs << val;

  return( vs.str() );
}

template<>
string GetString<char>( const char val ) {
  return( GetString<int>( val ) );
}



template<typename T>
int CompareArrays( const string f1, const string f2,
		   const string var ) {

  int retVal = EXIT_FAILURE;

  vector<T> d1, d2;

  T cmpVal = 0;

  ReadVariable( f1, var, d1 );
  ReadVariable( f2, var, d2 );

  switch( compMode ) {
  case errorL2:
    cmpVal = ComputeL2err( d1, d2 );
    cout << "Error in L2 norm was " << cmpVal << endl;
    break;

  case maxAbsDifference:
    cmpVal = ComputeMaxAbsDiff( d1, d2 );
    cout << "Max. Abs. Diff was " << GetString(cmpVal) << endl;
    break;

  default:
    cerr << "Unrecognised comparison mode " << compMode << endl;
    exit( EXIT_FAILURE );
  }

  if( cmpVal <= matchTol ) {
    retVal = MATCH_RETURN;
  } else if( cmpVal <= diffTol ) {
    retVal = DIFF_RETURN;
  } else {
    retVal = EXIT_FAILURE;
  }

  return( retVal );
}



// ===========================================================

int main( int argc, char *argv[] ) {

  // Get command line arguments
  ReadCommandLine( argc, argv );

  cout << "var = " << varName << endl;
  cout << "file1 = " << f1Name << endl;
  cout << "file2 = " << f2Name << endl;

  // Check up on basic information
  nc_type vType;

  Validate( f1Name, f2Name, varName, vType );

  int retVal;
  switch( vType ) {
  case NC_DOUBLE:
    retVal = CompareArrays<double>( f1Name, f2Name, varName );
    break;

  case NC_FLOAT:
    retVal = CompareArrays<float>( f1Name, f2Name, varName );
    break;

  case NC_CHAR:
    retVal = CompareArrays<char>( f1Name, f2Name, varName );
    break;

  case NC_INT:
    retVal = CompareArrays<int>( f1Name, f2Name, varName );
    break;

  default:
    cerr << "Unsupported variable type " << vType << endl;
    exit( EXIT_FAILURE );
  }


  return( retVal );
}
