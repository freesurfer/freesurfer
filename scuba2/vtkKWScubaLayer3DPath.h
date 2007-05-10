/**
 * @file  vtkKWScubaLayer3DPath.h
 * @brief A vtkKWScubaLayer that displays paths.
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: dsjen $
 *    $Date: 2007/05/10 18:48:29 $
 *    $Revision: 1.3 $
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
#ifndef vtkKWScubaLayer3DPath_h
#define vtkKWScubaLayer3DPath_h

#include <list>
#include "vtkKWScubaLayer.h"
#include "ScubaInfoItem.h"

//BTX
class ScubaCollectionPropertiesPath;
//ETX

class vtkActor;
class vtkActorCollection;
class vtkPolyDataMapper;
class vtkContourFilter;
class vtkPolyDataNormals;
class vtkStripper;
class vtkTubeFilter;
class vtkKWScaleWithEntry;
class vtkKWCheckButton;

class vtkKWScubaLayer3DPath : public vtkKWScubaLayer {

  public:

    static vtkKWScubaLayer3DPath* New ();
    vtkTypeRevisionMacro( vtkKWScubaLayer3DPath, vtkKWScubaLayer );
  
    // Description:
    // Sets the surface in the layer.
    //BTX
    void SetPathProperties ( ScubaCollectionPropertiesPath* const iProperties );
    //ETX

    // Description:
    // Load the surface and set up.
    void Create ();

    // Description:
    // Listens for: OpacityChanged, PathColorChanged
    //BTX
    void DoListenToMessage ( std::string const isMessage,
           void* const iData );
    //ETX
    
    // Description:
    // Callback for when the path representation changes
    void SetPathMode( int iMode );    
    
    
    // Description:
    // Callback for when the tube radius changes
    void SetTubeRadius( float iRadius );
    
    // Description:
    // Callback for changing the radius of the tube by the sample data
    void SetScaleTubeRadius( int ibIsScaling );

    // Description:
    // Callback for coloring the tube by the sample data
    void SetColoredTube( int ibIsColored );    
      
  protected:

    vtkKWScubaLayer3DPath ();
    virtual ~vtkKWScubaLayer3DPath ();
    
    // Description:
    // Populate a UI page with controls for this layer.
    virtual void AddControls ( vtkKWWidget* iPanel );
    virtual void RemoveControls ();

  private:
  
    void CreatePoint( vtkActor* pointActor, const int iPointIndex, const double iRadius );
    
    // Description:
    // Create the tube representation
    void CreateTube();
    
    // Description:
    // Create the multiple path tube representations.
    void CreateInitialTubes();
    
    // Description:
    // Updates the opacity of the path.
    void UpdateOpacity();
    
    // Description:
    // Updates color of the path
    void UpdatePathColor();
    
    // Description:
    // Updates the representation of the path.
    void UpdatePathMode();
        
    //BTX
    ScubaCollectionPropertiesPath const* mPathProperties;
  
    vtkActor* mPathActor;
    vtkActor* mStartPointActor;
    vtkActor* mEndPointActor;
    
    vtkActorCollection* mPathPointsCollection;
    
    vtkTubeFilter *mTubeFilter;
    vtkPolyDataMapper* mTubeMapper;
    vtkActor* mTubeActor;
    
    vtkKWScaleWithEntry *mScaleTubeRadius;    

    vtkKWCheckButton *mChkBtnScaleTubeRadius;
    vtkKWCheckButton *mChkBtnColorTubeRadius;
    //ETX

    vtkStripper *mTriangleStripper;
    vtkPolyDataNormals* mNormals;
    vtkPolyDataMapper* mMapper;

    //BTX
    // Tube representation type
    enum {
      TUBE_MODE = 0,
      THRESHOLD_MODE = 1
    };
    //ETX
    
    // can be tube or threshold
    int mPathMode;
    
    double mTubeRadius;
    
    bool mbIsTubeRadiusScaled;
    bool mbIsTubeRadiusColored;
        
};

#endif

