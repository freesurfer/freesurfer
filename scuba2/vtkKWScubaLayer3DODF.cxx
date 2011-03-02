/**
 * @file  vtkKWScubaLayer3DODF.cxx
 * @brief Displays 3D glyphs
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author$
 *    $Date$
 *    $Revision$
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

#include "vtkKWScubaLayer3DODF.h"

#include <string>
#include <stdexcept>

#include "ScubaCollectionProperties.h"
#include "ScubaCollectionPropertiesODF.h"
#include "ScubaInfoItem.h"
#include "vtkODFGlyph.h"
#include "vtkFSVolumeSource.h"
#include "vtkObjectFactory.h"

using namespace std;

vtkStandardNewMacro( vtkKWScubaLayer3DODF );
vtkCxxRevisionMacro( vtkKWScubaLayer3DODF, "$Revision$" );

vtkKWScubaLayer3DODF::vtkKWScubaLayer3DODF () :
  mODFProperties( NULL ) {
    
  for( int n = 0; n < 3; n++ ) {
    mGlyphs[ n ] = NULL;
    mWorldCenter[ n ] = 0.0;
    mWorldSize[ n ] = 0.0;
  }
    
}

vtkKWScubaLayer3DODF::~vtkKWScubaLayer3DODF () {
  
  for( int n=0; n<3; n++ ) {

    if ( mGlyphs[ n ] )
      mGlyphs[ n ]->Delete();
    
  }
  
}

void
vtkKWScubaLayer3DODF::SetODFProperties ( ScubaCollectionPropertiesODF const* iProperties ) {
  mODFProperties = iProperties;
}


void
vtkKWScubaLayer3DODF::Create () {

  // Bail if we don't have our source yet.
  if( NULL == mODFProperties )
    throw runtime_error( "vtkKWScubaLayer3DODF::Create: No source" );
    
}

void
vtkKWScubaLayer3DODF::AddControls ( vtkKWWidget* iPanel ) {
}

void
vtkKWScubaLayer3DODF::RemoveControls () {
}

void
vtkKWScubaLayer3DODF::DoListenToMessage ( string isMessage, void* iData ) {
}

void
vtkKWScubaLayer3DODF::GetRASBounds ( float ioBounds[6] ) const {

  if ( mODFProperties ) {

    mODFProperties->GetODFVolumeSource()->GetRASBounds( ioBounds );

  } else {

    for ( int nBound = 0; nBound < 6; nBound++ )
      ioBounds[nBound] = 0;

  }

}

void
vtkKWScubaLayer3DODF::Get3DRASZIncrementHint ( float ioHint[3]) const {

  ioHint[0] = 1.0;
  ioHint[1] = 1.0;
  ioHint[2] = 1.0;
}

void
vtkKWScubaLayer3DODF::GetInfoItems ( float iRAS[3],
                                      list<ScubaInfoItem>& ilInfo ) const {

  ScubaInfoItem info;

  info.Clear();
  info.SetLabel( "" );
  info.SetValue( "" );
  info.SetShortenHint( false );

  // Return it.
  ilInfo.push_back( info );  
  
}
