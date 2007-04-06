/**
 * @file  vtkKWScubaLayerCollectionPath.h
 * @brief Implementation for path viewers.
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:05 $
 *    $Revision: 1.1 $
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

#ifndef vtkKWScubaLayerCollectionPath_h
#define vtkKWScubaLayerCollectionPath_h

#include "vtkKWScubaLayerCollection.h"
#include "ScubaCollectionPropertiesPath.h"

class vtkFSVolumeSource;
class vtkSimplePointsReader;
class vtkContourFilter;
class vtkKWChangeColorButton;

class vtkKWScubaLayerCollectionPath : public vtkKWScubaLayerCollection
                                     //BTX
                                     , public ScubaCollectionPropertiesPath
                                     //ETX
{

  public:

    static vtkKWScubaLayerCollectionPath* New ();
    vtkTypeRevisionMacro( vtkKWScubaLayerCollectionPath, 
        vtkKWScubaLayerCollection );
        
    void SetVolumeFileName( const char* ifnVolume );

    vtkFSVolumeSource* GetPathVolumeSource () const;
    
    vtkSimplePointsReader* GetPathPointsSource() const;
    
    virtual vtkPolyData* GetMesh() const;

    // Description:
    // Populate a UI page with common controls for this layer type.
    virtual void AddControls ( vtkKWWidget* iPanel );
    virtual void RemoveControls ();

    void SetPathColor( const double r, const double g, const double b );

    virtual void GetPathColor( double &r, double &g, double &b ) const;

    // Description:
    // Get a pointer to the color that cooresponds to the underlying scalar value.
    virtual void GetPointColor( const int iPointIndex, double &r, double &g, double &b ) const;
    
    // Description.
    // Set the mode that we should threshold
    void SetPathThresholdMode ( int iMode );
    
    // Description.
    // Get the actual threshold value that the surface rendering of the path
    // should be generated at.
    double GetPathThreshold () const;
    
      
  protected:

    vtkKWScubaLayerCollectionPath ();
    ~vtkKWScubaLayerCollectionPath ();
    
    //BTX
    virtual vtkKWScubaLayer*
      MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode );
    //ETX

    void LoadPathVolumeFromFileName ();
    
    void MakeMesh();

    //BTX
    std::string GetFullFileName ( const char* sShortFileName ) const;
    //ETX

  private:

    vtkFSVolumeSource* mPathVolumeSource;
    vtkSimplePointsReader* mSimplePointsReader;
    vtkSimplePointsReader* mSamplesReader;
  
    //BTX
    std::string mfnPathVolume;
    //ETX

    vtkContourFilter* mContourFilter;

    vtkKWChangeColorButton *mPathColorButton;
    static const double DEFAULT_COLOR[ 3 ];
    double mPathColor[ 3 ];
    
    //BTX
    enum {
      PATH_THRESHOLD_LOW = 0,
      PATH_THRESHOLD_HIGH = 1,
    };
    //ETX
    
    double mPathThreshold;
    

};

#endif
