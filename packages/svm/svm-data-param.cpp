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


#include <stdlib.h>
#include <string.h>

#include "svm-data-param.h"

using namespace std;

void DataParam::initFileNames(int count) {
  _fileCount = count;
  fileName = new char*[_fileCount];
  for ( int i = 0; i < count; i++ )
    fileName[i] = new char[300];
}


bool DataParam::parse(const char* const* argv, int argc, int& k) {
  mode = NO_FILE;
  binary = false;
  while ( k < argc ) {
    const char* option = argv[k++];
    if ( option[0] != '-') {
      k--;
      break;
    }
    if ( !strcmp(option,"-binary") ) {
      binary = true;
      break;
    } else if ( !strcmp(option,"-file") ) {
      if ( mode != NO_FILE ) {
        cerr << "DataParam error: too many input data options specified.\n";
        mode = NO_FILE;
        return false;
      }
      initFileNames(1);
      sprintf(fileName[0], "%s", argv[k++]);
      mode = ONE;
    } else if ( !strcmp(option,"-2files") ) {
      if ( mode != NO_FILE ) {
        cerr << "DataParam error: too many input data options specified.\n";
        mode = NO_FILE;
        return false;
      }
      initFileNames(2);
      sprintf(fileName[0], "%s", argv[k++]);
      sprintf(fileName[1], "%s", argv[k++]);
      mode = TWO;
    } else if ( !strcmp(option,"-list") ) {
      if ( mode != NO_FILE ) {
        cerr << "DataParam error: too many input data options specified.\n";
        mode = NO_FILE;
        return false;
      }
      initFileNames(_posCount+_negCount);
      sprintf(fileNameFormat, "%s", argv[k++]);
      int i = 0;
      for ( ; i < _fileCount && k < argc; i++ )
        sprintf(fileName[i], "%s", argv[k++]);
      if ( i < _fileCount ) {
        cerr << "DataParam error: too few names specified with -list option.\n";
        mode = NO_FILE;
        return false;
      }
      mode = MANY;
    } else
      break;
  }

  if ( mode == NO_FILE ) {
    cerr << "DataParam Error: no data input option was specified. \n";
    return false;
  }
  return true;
}


ostream& DataParam::printUsage(ostream& s) {
  s << " (-file data_file | "
  << "-2files positive_data_file negative_data_file | "
  << "-list file_name_format name1 name2 ...) [-binary] ";

  return s;
}

ostream& DataParam::printUsageHelp(ostream& s) {
  s << "\nData input parameters: \n"
  << " -file data_file : load the entire data set from a single file\n"
  << " -2files positive_data_file negative_data_file : load the positive examples from\n"
  << "\t positive_data_file and the negative examples from negative_data_file\n"
  << " -list file_name_format name1 name2 ... : load each training example \n"
  << "\t from a separate file, create file names from the list using \n"
  << "\t file_name_format (sprintf)\n"
  << "IMPORTANT: the examples are assumed to be arranged so that all positive \n"
  << "\t examples precede all negative examples\n\t in the data files.\n\n"
  << " -binary : optional, if specified, the data are assumed to be in a binary fromat,\n"
  << "           default - text.\n";

  return s;
}


bool DataParam::readData(SvmRealMatrix& data) {
  switch (mode) {
  case ONE:
    return read(data,fileName[0],binary);
  case TWO: {
    SvmRealMatrix data1(_posCount,data.cols());
    SvmRealMatrix data2(_negCount,data.cols());
    if ( !read(data1,fileName[0],binary) )
      return false;
    if ( !read(data2,fileName[1],binary) )
      return false;
    for ( int i = 0; i < _posCount; i++ )
      data.setRow(i,data1[i]);
    for ( int i = 0; i < _negCount; i++ )
      data.setRow(i+_posCount,data2[i]);
    return true;
  }
  case MANY: {
    SvmRealVector vec(data.cols());
    for ( int i = 0; i < _posCount+_negCount; i++ ) {
      if ( !read(vec,fileName[i],binary) )
        return false;
      data.setRow(i,vec);
    }
    return true;
  }
  default:
    cerr << "DataParam Error: no data input option was specified. \n";
  }

  return false;
}





