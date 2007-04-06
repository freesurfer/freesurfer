/**
 * @file  ScubaCollectionPropertiesMRIS.h
 * @brief The common properties available to MRIS layers
 *
 * An abstract interface implemented by a collection. Layers will get
 * a pointer to an object of this type so they can get access to
 * shared layer settings.
 */
/*
 * Original Author: Kevin Teich
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

#ifndef ScubaCollectionPropertiesMRIS_h
#define ScubaCollectionPropertiesMRIS_h

class vtkFSSurfaceSource;
class vtkAlgorithmOutput;
class vtkPolyData;

class ScubaCollectionPropertiesMRIS {

 public:

  // Description:
  // Get a pointer to the surface poly data output.
  virtual vtkFSSurfaceSource* GetSource () const = 0;

  // Description:
  // These are the output ports for fast mode (decimated) and normal
  // mode (full detail) surface meshes.
  virtual vtkAlgorithmOutput* GetNormalModeOutputPort () const = 0;
  virtual vtkAlgorithmOutput* GetFastModeOutputPort () const = 0;
  virtual vtkPolyData* GetNormalModeOutput () const = 0;
  virtual vtkPolyData* GetFastModeOutput () const = 0;
  
};

#endif
