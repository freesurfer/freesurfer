/**
 * @file  svm-data-param.h
 * @brief REPLACE_WITH_ONE_LINE_SHORT_DESCRIPTION
 *
 * REPLACE_WITH_LONG_DESCRIPTION_OR_REFERENCE
 */
/*
 * Original Author: REPLACE_WITH_FULL_NAME_OF_CREATING_AUTHOR 
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.4 $
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


////SVM-LIB////////////////////////////////////////////////////////////////
//
// Name: SvmParam
//
// This file defines all the parameters fro the svm optimization.
//
//  Polina Golland polina@ai.mit.edu
//
///////////////////////////////////////////////////////////////////////////



#ifndef __SVM_DATA_PARAM_H__
#define __SVM_DATA_PARAM_H__

#include "svm-io-format.h"
#include "svm-vector-types.h"

class DataParam {
public:
  static const int NO_FILE = 0;
  static const int ONE     = 1;
  static const int TWO     = 2;
  static const int MANY    = 3;

  int mode;
  bool binary;

private:
  int _fileCount;
  char **fileName;

  int _posCount, _negCount;
  char fileNameFormat[300];

  void initFileNames(int count);


public:

  // Constructors
  DataParam() : mode(NO_FILE), binary(false),
      _fileCount(0), fileName(NULL), _posCount(0), _negCount(0) {}

  ~DataParam() {
    if ( _fileCount > 0 ) {
      for ( int i = 0; i < _fileCount; i++ )
        delete [] fileName[i];
      delete [] fileName;
    }
  }


  void init (int posCount, int negCount) {
    _posCount = posCount;
    _negCount = negCount;
  }


  // I/O
  bool parse(const char* const* argv, int argc, int& startIndex);

  bool readData(SvmRealMatrix& data);

  static std::ostream& printUsage(std::ostream& s);
  static std::ostream& printUsageHelp(std::ostream& s);
};


#endif // __SVM_DATA_PARAM_H__




