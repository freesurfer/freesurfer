/**
 * @file  vtkKWScubaLayerCollectionPath.h
 * @brief Implementation for path viewers.
 *
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.7 $
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

    //BTX
    const std::vector< std::vector< double* >* > * GetInitialPaths() const;
    //ETX
    
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
    
    // Description:
    // Get the scalar value at a point.
    virtual double GetPointSampleValue( const int iPointIndex ) const;
    
    // Description:
    // Get the number of samples that were read.
    virtual int GetNumberOfSamples() const;

    // Description:
    // Get the probablity at each point.
    virtual double GetPointProbabilityValue( const int iPointIndex ) const;
    
    // Description:
    // Get the number of probabilities that were read.
    virtual int GetNumberOfProbabilities() const;
    
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
    
    //BTX
    // Description:
    // Reads in the sample values.
    void ReadScalars( const char* fnSamples, std::vector< double >* ioScalars );
    //ETX

    // Description:
    // Reads the intial path points and parses them into multiple pathways.
    void ReadInitialPathPoints( const char* fnInitialPathPoints );
    
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
    
    //BTX
    
    // scalar values sampled by poistats
    std::vector< double > mSamples;

    // probabilities that the path goes through the voxel
    std::vector< double > mProbabilities;
    
    // this is a vector of pathways, each consisting of 3D points
    std::vector< std::vector< double* >* > mInitialPaths;
    //ETX

};

#endif
