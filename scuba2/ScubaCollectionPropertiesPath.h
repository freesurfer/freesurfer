/**
 * @file  ScubaCollectionPropertiesPath.h
 * @brief The common properties available to path layers
 *
 * An abstract interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 * 
 * This is only for accessing data.
 */
/*
 * Original Author: Dennis Jen
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.8 $
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

#ifndef ScubaCollectionPropertiesPath_h
#define ScubaCollectionPropertiesPath_h

#include <vector>

class vtkFSVolumeSource;
class vtkSimplePointsReader;
class vtkPolyData;

class ScubaCollectionPropertiesPath {

 public:

  // Description:
  // Get a pointer to the source volume.
  virtual vtkFSVolumeSource* GetPathVolumeSource () const = 0;

  virtual vtkSimplePointsReader* GetPathPointsSource() const = 0;
  
  virtual const std::vector< std::vector< double* >* > * GetInitialPaths() const = 0;
    
  // Description:
  // Get a pointer to the surface representation of the pathway.
  virtual vtkPolyData* GetMesh() const = 0;
  
  virtual void GetPathColor( double &r, double &g, double &b ) const = 0;
  
  // Description:
  // Get a pointer to the color that cooresponds to the underlying scalar value.
  virtual void GetPointColor( const int iPointIndex, double &r, double &g, double &b ) const = 0;
  
  virtual double GetPointSampleValue( const int iPointIndex ) const = 0;
  
  virtual int GetNumberOfSamples() const = 0;

  virtual double GetPointProbabilityValue( const int iPointIndex ) const = 0;
  
  virtual int GetNumberOfProbabilities() const = 0;

 protected:

  ScubaCollectionPropertiesPath () {};
  virtual ~ScubaCollectionPropertiesPath () {};
};

#endif
