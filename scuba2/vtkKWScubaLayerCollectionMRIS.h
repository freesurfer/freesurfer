/**
 * @file  vtkKWScubaLayerCollectionMRIS.h
 * @brief Implementation for MRIS viewers.
 *
 * In 2D, the MRIS is viewed as the intersection of the surface with a
 * slice. In 3D, the MRIS is viewed as a fully 3D mesh surface.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: kteich $
 *    $Date: 2007/07/25 19:53:47 $
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

#ifndef vtkKWScubaLayerCollectionMRIS_h
#define vtkKWScubaLayerCollectionMRIS_h

#include "vtkKWScubaLayerCollection.h"
#include "ScubaCollectionPropertiesMRIS.h"
#include "vtkSmartPointer.h"

class vtkAlgorithmOutput;
class vtkDecimatePro;
class vtkFloatArray;
class vtkFreesurferLookupTable;
class vtkFSSurfaceSource;
class vtkKWLabel;
class vtkTransformPolyDataFilter;

class vtkKWScubaLayerCollectionMRIS : public vtkKWScubaLayerCollection
                    //BTX
                    , public ScubaCollectionPropertiesMRIS
		    //ETX
{

 public:

  static vtkKWScubaLayerCollectionMRIS* New ();
  vtkTypeRevisionMacro( vtkKWScubaLayerCollectionMRIS, 
			vtkKWScubaLayerCollection );

  // Description:
  // Set the surface file name.
  void SetSurfaceFileName ( const char* ifnVolume );

  // Description:
  // Set the name of the annotation scalar file and the color table
  // file, if needed. If the color table file name is NULL, it will
  // try to find the color table from the scalar file.
  void SetAnnotationAndColorTableFileNames ( const char* ifnAnnotation,
					     const char* ifnColorTable );

  // Description:
  // Populate a UI page with common controls for this layer type.
  virtual void AddControls ( vtkKWWidget* iPanel );
  virtual void RemoveControls ();

  // Description:
  // Show a dialog to choose the annotation file to open.
  void LoadAnnotationFromDlog ();

  // Layers will get the source objects with these calls. These
  // impelement ScubaCollectionPropertiesMRIS.
  vtkFSSurfaceSource* GetSource () const;
  vtkAlgorithmOutput* GetNormalModeOutputPort () const;
  vtkAlgorithmOutput* GetFastModeOutputPort () const;
  vtkPolyData* GetNormalModeOutput () const;
  vtkPolyData* GetFastModeOutput () const;
  vtkFloatArray* GetScalarsValues () const;
  vtkScalarsToColors* GetScalarsColors () const;

 protected:

  vtkKWScubaLayerCollectionMRIS ();
  ~vtkKWScubaLayerCollectionMRIS ();

  //BTX
  virtual vtkKWScubaLayer*
    MakeLayerForDisplayMode ( vtkKWScubaView::DisplayMode iMode );
  //ETX

  void LoadSurfaceFromFileName ();

  void LoadAnnotationAndColorTableFromFileNames ();

 private:

  //BTX
  vtkSmartPointer<vtkFSSurfaceSource> mSource;
  vtkSmartPointer<vtkDecimatePro> mDecimator;
  std::string mfnSurface;
  vtkSmartPointer<vtkFloatArray> mAnnotationScalars;
  vtkSmartPointer<vtkFreesurferLookupTable> mAnnotationTable;
  std::string mfnAnnotationScalars;
  std::string mfnAnnotationTable;

  vtkSmartPointer<vtkKWLabel> mLabelAnnotationFileName;
  //ETX

};

#endif
