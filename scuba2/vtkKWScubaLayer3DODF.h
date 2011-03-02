/**
 * @file  vtkKWScubaLayer3DODF.h
 * @brief A vtkKWScubaLayer that displays ODF glyphs
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

#ifndef vtkKWScubaLayer3DODF_h
#define vtkKWScubaLayer3DODF_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesODF;
//ETX
class vtkODFGlyph;

class vtkKWScubaLayer3DODF : public vtkKWScubaLayer {

public:

  static vtkKWScubaLayer3DODF* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayer3DODF, vtkKWScubaLayer );

  // Description:
  // Sets the surface in the layer.
  //BTX
  void SetODFProperties ( ScubaCollectionPropertiesODF const* iProperties );
  //ETX

  // Description:
  // Listens for:
  //BTX
  void DoListenToMessage ( std::string const isMessage,
         void* const iData );
  //ETX

  // Description:
  // Load the surface and set up.
  void Create ();

  // Returns the bounds of the volume.
  virtual void GetRASBounds ( float ioBounds[6] ) const;

  // Description:
  // Returns the shortest edge distance.
  virtual void Get3DRASZIncrementHint ( float ioHint[3]) const;

  // Description:
  // Return our index and value info.
  //BTX
  void GetInfoItems ( float iRAS[3], std::list<ScubaInfoItem>& ilInfo ) const;
  //ETX  

protected:

  vtkKWScubaLayer3DODF ();
  virtual ~vtkKWScubaLayer3DODF ();

  // Description:
  // Populate a UI page with controls for this layer.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();
    
  //BTX
  ScubaCollectionPropertiesODF const* mODFProperties;

  // Pipeline -------------------------------------------------------------
  vtkODFGlyph* mGlyphs[ 3 ];
  // ----------------------------------------------------------------------

  float mWorldCenter[3];
  float mWorldSize[3];
  
  //ETX

private:
};

#endif
