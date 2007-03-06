/**
 * @file  vtkFSSurfaceScalarsReader.h
 * @brief Reads a scalar file and outputs a float array
 *
 * This is the equivalent of calling MRISreadValues() except that it
 * doesn't actually require an MRIS. It attemps to read a few
 * different scalar file types, and outputs a vtkFloatArray if
 * successful.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/03/06 15:25:30 $
 *    $Revision: 1.2 $
 *
 * Copyright (C) 2002-2007,
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

#ifndef vtkFSSurfaceScalarsReader_h
#define vtkFSSurfaceScalarsReader_h

#include <string>
#include "vtkPolyDataAlgorithm.h"

class vtkFloatArray;

class vtkFSSurfaceScalarsReader : public vtkPolyDataAlgorithm {

 public:
  
  static vtkFSSurfaceScalarsReader* New ();
  vtkTypeRevisionMacro( vtkFSSurfaceScalarsReader, vtkPolyDataAlgorithm );
  void PrintSelf ( ostream& os, vtkIndent indent );

  void SetFileName ( const char* ifn );
  const char* GetFileName () const;

  vtkSetMacro(NumberOfValues,int);
  vtkGetMacro(NumberOfValues,int);

 protected:

  vtkFSSurfaceScalarsReader ();
  virtual ~vtkFSSurfaceScalarsReader ();

  // Description:
  // Read the file, create the array of floats, and set our output.
  virtual int RequestData ( vtkInformation*,
			    vtkInformationVector**,
			    vtkInformationVector* iOutputVector );


  static const int NEW_SCALAR_MAGIC_NUMBER = 16777215;

  std::string FileName;
  int NumberOfValues;

 private:
  vtkFSSurfaceScalarsReader ( const vtkFSSurfaceScalarsReader& );
  void operator= ( const vtkFSSurfaceScalarsReader& );
};

#endif
