/**
 * @file  vtkKWScubaLayer2DPath.h
 * @brief A vtkKWScubaLayer that displays paths.
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/10/16 16:01:50 $
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

