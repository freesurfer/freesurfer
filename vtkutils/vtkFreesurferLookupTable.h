/**
 * @file  vtkFreesurferLookupTable.h
 * @brief A VTK table that reads Freesurfer LUT files
 *
 * This is a vtkLookupTable subclass that can read the Freesurfer LUT
 * file format.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/18 19:11:38 $
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

#ifndef vtkFreesurferLookupTable_h
#define vtkFreesurferLookupTable_h

#include <string>

#include "vtkLookupTable.h"

extern "C" {
#include "colortab.h"
}

class vtkFreesurferLookupTable : public vtkLookupTable {

 public:
  
  static vtkFreesurferLookupTable* New();
  vtkTypeRevisionMacro( vtkFreesurferLookupTable, vtkLookupTable );

  // Description:
  // Clears and sets its internal entries by reading a Freesurfer LUT
  // file.
  void BuildFromCTAB ( COLOR_TABLE* iCtab );

 protected:

  vtkFreesurferLookupTable();
  virtual ~vtkFreesurferLookupTable();

private:
  vtkFreesurferLookupTable( const vtkFreesurferLookupTable& );
  void operator=( const vtkFreesurferLookupTable& );

};

#endif
