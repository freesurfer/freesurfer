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
 *    $Author: nicks $
 *    $Date: 2011/03/02 00:04:40 $
 *    $Revision: 1.4 $
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

#ifndef ScubaCollectionPropertiesMRIS_h
#define ScubaCollectionPropertiesMRIS_h

class vtkAlgorithmOutput;
class vtkFloatArray;
class vtkFSSurfaceSource;
class vtkPolyData;
class vtkScalarsToColors;

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

  // Description:
  // Get any scalars that should be displayed.
  virtual vtkFloatArray* GetScalarsValues () const = 0;

  // Description:
  // Get the color map for scalars.
  virtual vtkScalarsToColors* GetScalarsColors () const = 0;
  
 protected:

  ScubaCollectionPropertiesMRIS () {};
  virtual ~ScubaCollectionPropertiesMRIS () {};
};

#endif
