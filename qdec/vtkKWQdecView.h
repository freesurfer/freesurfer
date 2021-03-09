/**
 * @brief The VTK pipeline and Renderer
 *
 * Builds the VTK pipeline to display the surface and scalars. Also
 * manages the callback for clicking on a vertex and passing the point
 * along to the Fsgdf plotting code.
 */
/*
 * Original Author: Kevin Teich
 *
 * Copyright © 2021 The General Hospital Corporation (Boston, MA) "MGH"
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


#ifndef vtkKWQdecView_h
#define vtkKWQdecView_h

#include <map>
#include <string>
#include "vtkKWRenderWidget.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkSmartPointer.h"
#include "FsgdfPlot.h"

//BTX
class QdecVertexAnnotationLookup;
//ETX
class vtkActor;
class vtkCellArray;
class vtkCursor3D;
class vtkDecimatePro;
class vtkFSSurfaceScalarsReader;
class vtkFSSurfaceSource;
class vtkFloatArray;
class vtkInflatePolyData;
class vtkOutlineFilter;
class vtkPoints;
class vtkPolyData;
class vtkPolyDataMapper;
class vtkScalarBarActor;
class vtkScalarsToColors;
class vtkTransformPolyDataFilter;

class vtkKWQdecView : public vtkKWRenderWidget {

public:
  
  static vtkKWQdecView* New ();
  vtkTypeRevisionMacro( vtkKWQdecView, vtkKWRenderWidget );
  
  // Set up the interactor callback.
  void CreateWidget ();

  // This will reset the camera so it can view the entire surface from
  // a lateral angle, and will define the 'original location' that
  // RestoreView will use.
  void ResetView ();

  // Reset the view the the original location.
  void RestoreView ( const char* isHemi );

  // Dolly the camera.
  void ZoomBy ( float iFactor );

  // Get/set state of display of the cursor.
  bool GetShowCursor ();
  void SetShowCursor ( bool ibShow );

  // When a surface vertex is selected by the user, a
  // QDECEvents::UserSelectedVertex message will be sent, and the call
  // data will be a pointer to this struct about the vertex clicked.
  //BTX
  struct SurfaceVertexInformation {
    int mnVertexIndex;
  };
  //ETX

  // The surface to display. Setting a surface will delete an existing
  // one. This will also create all the pipeline objects. Nothing can
  // be drawn if there is no surface set first.
  void SetSurface ( vtkFSSurfaceSource* iSurfaceSource );

  // Set a scalar array to draw on the surface.
  void SetSurfaceScalars ( vtkFloatArray* iScalars );

  // The color table for the drawn scalars.
  void SetSurfaceScalarsColors ( vtkScalarsToColors* iScalarColors );

  // This array will be used to get values to print in the annotation
  // when a point is selected.
  void SetSurfaceLookupScalars ( vtkFloatArray* iScalars );

  // This is the overlay color scalars and color map. They are
  // assigned to a translucent surface that's just a little bit
  // inflated, floating over the main surface geometry.
  void SetSurfaceOverlayScalarsAndColors ( vtkFloatArray* iScalars,
					   vtkScalarsToColors* iColors );

  // Get and set the opacity of the overlay surface.
  void SetSurfaceOverlayOpacity ( double iOpacity );
  double GetSurfaceOverlayOpacity ();

  // This color table will be used for the ScalarBar legend.
  void SetSurfaceLegendColors ( vtkScalarsToColors* iLegendColors );

  // This is poly data to draw as an ROI. It should be a subset of the
  // surface data, so we inflate it a bit.
  void SetROI ( vtkPolyData* iROIPolyData );

  // Set a message to be displayed in the top left corner of the view.
  void SetAnnotationMessage ( const char* isMessage );

  // Show or hide the scalar bar actor and annotation actor.
  void SetShowLegend ( int ibShow );
  void SetShowAnnotation ( int ibShow );

  // Sets a lookup object so the view can get an annotation string for
  // each vertex clicked.
  void SetVertexAnnotationLookup ( QdecVertexAnnotationLookup* iLookup );

  // This is the list of points containing the loop of the current
  // selection. The window will use it as input to vtkSelectPolyData
  // for ROI operations.
  //BTX
  vtkPoints* GetSurfaceSelectionPoints ();
  //ETX

  // Rotate commands. These do a cool little animation of rotating the
  // camera arond the object.
  void AnimateCameraElevateNegative ();
  void AnimateCameraElevatePositive ();
  void AnimateCameraAzimuthNegative ();
  void AnimateCameraAzimuthPositive ();
  void AnimateCameraRollNegative ();
  void AnimateCameraRollPositive ();

  //BTX
  static double const kAnimationDegrees;
  static double const kcAnimateSteps;
  static double const kAnimateTimeInSeconds;
  //ETX

  // Select the numbered vertex.
  void SelectSurfaceVertex ( int inVertex );

  void SetVertexPlot( FsgdfPlot* iPlot ) { this->mVertexPlot = iPlot; };

 protected:

  vtkKWQdecView ();
  virtual ~vtkKWQdecView ();

  // Set up the scalar bar actor.
  void InitializeScalarBar ();
  void DeleteScalarBar ();

  // Set up the ROI actor.
  void InitializeROI ();

  // Reset (clear) the current path, add a vertex to the current path,
  // and rebuild the path line. You should rebuild the line after
  // adding a batch of points.
  void ResetPath ();
  void AddVertexToPath ( int inVertex );
  void RebuildPathLine ();

  //BTX
  std::string mfnSurface;

  // Pipeline objects.
  vtkSmartPointer<vtkFSSurfaceSource> mCurrentSurfaceSource;
  vtkSmartPointer<vtkPolyDataMapper> mSurfaceMapper;
  vtkSmartPointer<vtkActor> mSurfaceActor;
  vtkSmartPointer<vtkScalarBarActor> mScalarBar;
  vtkSmartPointer<vtkCursor3D> mCursor;
  vtkSmartPointer<vtkInflatePolyData> mSurfaceOverlayInflator;
  vtkSmartPointer<vtkPolyDataMapper> mSurfaceOverlayMapper;
  vtkSmartPointer<vtkActor> mSurfaceOverlayActor;
  vtkSmartPointer<vtkPoints> mSurfaceSelectionPoints; 
  vtkSmartPointer<vtkCellArray> mSurfaceSelectionLines;
  vtkSmartPointer<vtkPolyData> mSurfaceSelectionPolyData;
  vtkSmartPointer<vtkPolyData> mROIPolyData;
  vtkSmartPointer<vtkPolyDataMapper> mROIMapper;
  vtkSmartPointer<vtkActor> mROIActor;
  
  FsgdfPlot* mVertexPlot;

  // Because the scalars on display are actually a composite of the
  // curvature and the overlay, we need to have separate arrays for
  // the values we're actually displaying and the values we want to
  // use as the value legend and scalar bar legend.
  vtkSmartPointer<vtkFloatArray> mCurrentScalars;            // These have
  vtkSmartPointer<vtkScalarsToColors> mCurrentScalarsColors; // curv + overlay

  vtkSmartPointer<vtkFloatArray> mCurrentLookupScalars;      // These have
  vtkSmartPointer<vtkScalarsToColors> mCurrentLegendColors;  // just overlay

  // vtkSmartPointer<These scalars and colors are for the overlay surface.
  vtkSmartPointer<vtkFloatArray> mCurrentOverlayScalars;
  vtkSmartPointer<vtkScalarsToColors> mCurrentOverlayColors;

  // Our position to go to when we hit RestoreView.
  double mDefaultPosition[3];
  double mDefaultFocalPoint[3];
  double mDefaultViewUp[3];
  float  mDefaultZoom;

  // The currently selected vertex at the cursor. We use this to sync
  // the cursor when switching surfaces.
  int mnCursorVertexIndex;
  
  // The opacity of the overlay surface.
  double mOverlayOpacity;

  // Pointer to the annotation lookup.
  QdecVertexAnnotationLookup* mAnnotationLookup;

  // For keeping track of the lines we're drawing.
  bool mbInSelection;
  int mnFirstVertexInPath;
  int mnLastVertexInPath;

  // Our InteractorStyle. It's based on the trackball camera so it
  // does all that stuff, but we also listen to some mouse downs for
  // tool events, and even suppress normal dragging bevahior for the
  // drawing tool.
  friend class ViewInteractor;
  class ViewInteractor : public vtkInteractorStyleTrackballCamera {
  public:
    static ViewInteractor* New ();
    void SetView( vtkKWQdecView* iView );
    void OnLeftButtonDown();
    void OnLeftButtonUp();
    void OnMouseMove();
    int GetVertexAtPicker ();
  protected:
    ViewInteractor();
    vtkSmartPointer<vtkKWQdecView> mView;
    int mnButtonDown;
  };
  //ETX
};

#endif
