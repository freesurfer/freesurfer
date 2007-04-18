/**
 * @file  vtkFreesurferLookupTable.cxx
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

#include "vtkFreesurferLookupTable.h"
#include "vtkObjectFactory.h"

using namespace std;

vtkStandardNewMacro( vtkFreesurferLookupTable );
vtkCxxRevisionMacro( vtkFreesurferLookupTable, "$Revision: 1.2 $" );

vtkFreesurferLookupTable::vtkFreesurferLookupTable () {

  // Set up with all black.
  this->SetNumberOfTableValues( 1 );
  this->SetTableRange( 0, 1 );
  this->SetValueRange( 0, 0 );
  this->Build();
}

vtkFreesurferLookupTable::~vtkFreesurferLookupTable () {}

void
vtkFreesurferLookupTable::BuildFromCTAB ( COLOR_TABLE* iCtab ) {

  // Go through and make entries for each valid entry we got.
  int cEntries;
  CTABgetNumberOfTotalEntries( iCtab, &cEntries );

  // Build our VTK table with the correct number of entries.
  this->SetNumberOfTableValues( cEntries );
  this->SetTableRange( 0, cEntries );
  this->Build();

  // Zero is always clear.
  this->SetTableValue( 0, 0, 0, 0, 0 );

  // Populate those entries.
  for ( int nEntry = 1; nEntry < cEntries; nEntry++ ) {
    int bValid;
    CTABisEntryValid( iCtab, nEntry, &bValid );
    if ( bValid ) {
      float red, green, blue, alpha;
      CTABrgbaAtIndexf( iCtab, nEntry, &red, &green, &blue, &alpha );
      this->SetTableValue( nEntry, red, green, blue, alpha );
    }
  }
}
