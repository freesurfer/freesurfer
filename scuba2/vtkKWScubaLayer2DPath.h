/**
 * @file  vtkKWScubaLayer2DPath.h
 * @brief A vtkKWScubaLayer that displays paths.
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
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
#ifndef vtkKWScubaLayer2DPath_h
#define vtkKWScubaLayer2DPath_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesPath;
//ETX

class vtkActor;
class vtkPolyDataMapper;
class vtkPolyDataNormals;
class vtkStripper;
class vtkPlane;

class vtkKWScubaLayer2DPath : public vtkKWScubaLayer {

  public:

    static vtkKWScubaLayer2DPath* New ();
    vtkTypeRevisionMacro( vtkKWScubaLayer2DPath, vtkKWScubaLayer );
  
    // Description:
    // Sets the surface in the layer.
    //BTX
    void SetPathProperties ( ScubaCollectionPropertiesPath const* iProperties );
    //ETX

    // Description:
    // Load the surface and set up.
    void Create ();

    // Returns the bounds of the surface.
    virtual void GetRASBounds ( float ioBounds[6] ) const;
  
    // Description:
    // Returns the shortest edge distance.
    virtual void Get2DRASZIncrementHint ( float ioHint[3]) const;

    // Description:
    // Listens for:
    //BTX
    void DoListenToMessage ( std::string const isMessage,
           void* const iData );
    //ETX
  
  protected:

    vtkKWScubaLayer2DPath ();
    virtual ~vtkKWScubaLayer2DPath ();

  private:
  
    void UpdatePathColor ();

    //BTX
    ScubaCollectionPropertiesPath const* mPathProperties;
  
    vtkPlane* mSlicePlane;
//    vtkPolyDataMapper* mFastMapper;
    vtkActor* mActor;
    //ETX

    vtkStripper *mTriangleStripper;
    vtkPolyDataNormals* mNormals;
    vtkPolyDataMapper* mMapper;

    // Surface info.
    float mRASCenter[3];
  
};

#endif

