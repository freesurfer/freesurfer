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

#include "vtkFreesurferLookupTable.h"
#include "vtkObjectFactory.h"

using namespace std;

vtkStandardNewMacro( vtkFreesurferLookupTable );

vtkFreesurferLookupTable::vtkFreesurferLookupTable () {

  // Set up with all black.
  this->SetNumberOfTableValues( 1 );
  this->SetTableRange( 0, 1 );
  this->SetValueRange( 0, 0 );
  this->Build();
}

vtkFreesurferLookupTable::~vtkFreesurferLookupTable () {}

void
vtkFreesurferLookupTable::BuildFromCTAB ( COLOR_TABLE* iCtab, bool bClearZero ) {

  // Go through and make entries for each valid entry we got.
  int cEntries;
  CTABgetNumberOfTotalEntries( iCtab, &cEntries );

  // Build our VTK table with the correct number of entries.
  this->SetNumberOfTableValues( cEntries );
  this->SetTableRange( 0, cEntries );
  this->Build();

  // Set zeros. Clear if requested (default).
  this->SetTableValue( 0, 0, 0, 0, bClearZero?0:1 );

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
