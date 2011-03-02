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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:56 $
 *    $Revision: 1.3 $
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
