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
 *    $Author: kteich $
 *    $Date: 2007/04/06 22:23:04 $
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

#ifndef ScubaCollectionPropertiesPath_h
#define ScubaCollectionPropertiesPath_h

class vtkFSVolumeSource;
class vtkSimplePointsReader;
class vtkPolyData;

class ScubaCollectionPropertiesPath {

 public:

  // Description:
  // Get a pointer to the source volume.
  virtual vtkFSVolumeSource* GetPathVolumeSource () const = 0;

  virtual vtkSimplePointsReader* GetPathPointsSource() const = 0;
  
  // Description:
  // Get a pointer to the surface representation of the pathway.
  virtual vtkPolyData* GetMesh() const = 0;
  
  virtual void GetPathColor( double &r, double &g, double &b ) const = 0;
  
  // Description:
  // Get a pointer to the color that cooresponds to the underlying scalar value.
  virtual void GetPointColor( const int iPointIndex, double &r, double &g, double &b ) const = 0;

};

#endif
