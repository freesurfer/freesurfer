/**
 * @brief A VTK table that reads Freesurfer LUT files
 *
 * This is a vtkLookupTable subclass that can read the Freesurfer LUT
 * file format.
 */
/*
 * Original Author: Kevin Teich
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

#ifndef vtkFreesurferLookupTable_h
#define vtkFreesurferLookupTable_h

#include <string>
//#include "vtkRenderingCoreModule.h" // For export macro
#include "vtkLookupTable.h"

#include "colortab.h"

class vtkFreesurferLookupTable : public vtkLookupTable {

 public:
  
  static vtkFreesurferLookupTable* New();
  vtkTypeMacro( vtkFreesurferLookupTable, vtkLookupTable );

  // Description:
  // Clears and sets its internal entries by reading a Freesurfer LUT
  // file.
  void BuildFromCTAB ( COLOR_TABLE* iCtab, bool bClearZero = true );

 protected:

  vtkFreesurferLookupTable();
  virtual ~vtkFreesurferLookupTable();

private:
  vtkFreesurferLookupTable( const vtkFreesurferLookupTable& );
  void operator=( const vtkFreesurferLookupTable& );

};

#endif
