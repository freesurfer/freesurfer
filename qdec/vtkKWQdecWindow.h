/**
 * @file  vtkKWQdecWindow.h
 * @brief Main QDEC logic
 *
 * Loads in all types of data and manages them. Main logic for running
 * analysis and handling results. Composes display objects and sends
 * it to the main view. Holds settings for color scales and currently
 * displayed objects.
 */
/*
 * Original Author: Kevin Teich
 * CVS Revision Info:
 *    $Author: nicks $
 *    $Date: 2008/01/17 02:41:09 $
 *    $Revision: 1.5 $
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


#ifndef vtkKWQdecWindow_h
#define vtkKWQdecWindow_h

#include <map>
#include <string>
#include <vector>

extern "C" {
#include "colortab.h"
}

#include "vtkKWWindow.h"
#include "ProgressUpdateGUI.h"
#include "QdecProject.h"
#include "QdecVertexAnnotationLookup.h"
#include "vtkKWQdecView.h"
#include "vtkSmartPointer.h"

class vtkFSSurfaceSource;
class vtkFSSurfaceLabelSource;
class vtkFloatArray;
class vtkIntArray;
class vtkKWBltGraph;
class vtkKWCheckButton;
class vtkKWEntry;
class vtkKWHistogram;
class vtkKWListBox;
class vtkKWMenu;
class vtkKWMenuButton;
class vtkKWMultiColumnList;
class vtkKWPushButton;
class vtkKWRadioButtonSet;
class vtkKWRGBATransferFunctionEditor;
class vtkKWScale;
class vtkRGBATransferFunction;
class vtkScalarsToColors;

class vtkKWQdecWindow : public vtkKWWindow 
                        //BTX
                        , public ProgressUpdateGUI 
			, public QdecVertexAnnotationLookup
                        //ETX
{

 public:

  static vtkKWQdecWindow* New ();
  vtkTypeRevisionMacro( vtkKWQdecWindow, vtkKWWindow );

  // Set whether or not to use the histogram editor. Should be called
  // before Create.
  void SetUseHistogramEditor ( bool ibUse );

  // Override to create our interior.
  virtual void CreateWidget ();

  // Finish setting up after the KWWidgets Create calls are done.
  void FinishCreating ();

  // Implement ProgressUpdateGUI virtuals
  void BeginActionWithProgress( const char* isTitle );
  void UpdateProgressMessage( const char* isMessage );
  void UpdateProgressPercent( float iPercent );
  void EndActionWithProgress();

  // The functions that end in FromDlog are just UI wrappers that show
  // a dialog box, prompt the user to select a file or directory, and
  // then pass the result to the corresponding function without the
  // FromDlog suffix.
  void LoadDataTableFromDlog ();
  void LoadProjectFileFromDlog ();
  void LoadSurfaceFromDlog ();        //
  void LoadGDFFromDlog ();            // These are ifdef'd out
  void LoadSurfaceScalarsFromDlog (); // and only used for
  void LoadCurvatureFromDlog ();      // debugging.
  void LoadAnnotationFromDlog ();
  void LoadLabelFromDlog ();
  void SaveProjectFileFromDlog ();
  void SaveScatterPlotPostscriptFromDlog ();
  void SaveTIFFImageFromDlog ();
  void SaveGDFPostscriptFromDlog ();
  void SaveLabelFromDlog ();
  void MapLabelFromDlog ();

  // These prompt a numeric value to specify a number of steps by
  // which to smooth the values.
  void SmoothCurvatureScalarsFromDlog ();
  void SmoothSurfaceScalarsFromDlog ();

  // Load the data table and update the subjects and design page with
  // the loaded data.
  void LoadDataTable ( const char* ifnDataTable );

  // Load a project file, which will contain all the data and info we
  // need to load analyzed data.
  void LoadProjectFile ( const char* ifnProject );

  // Save analyzed data to a project file.
  void SaveProjectFile ( const char* ifnProject );

  // Load analyzed data from info contained in results object. This includes:
  // load the inflated surface, load the curvature, load the aparc.annot 
  // annotation, load any contrast files, load the stddev and regression
  // coefficients scalars, and load the GDF.
  void LoadAnalyzedData ( QdecGlmFitResults* iGlmResults );

  // Load a surface and associate it with a label. Multiple surfaces
  // may be loaded and identified with their label. Displays the
  // surface in the view.
  void LoadSurface ( const char* ifnSurface, const char* isLabel=NULL );

  // Load a GDF file using Fsgdf_Read. The ID will be passed to the
  // view so that it can use the Fsgdf_* plotting functions in Tcl
  // land to plot data.
  void LoadGDFFile ( const char* ifnGDFFile );

  // Load a surface scalar file and return the entry index associated
  // with it. Optionally takes a label name to associate with this label file.
  // Optionally takes a frame number if reading a multi-frame file. 
  // This does not change the active scalars file.
  int  LoadSurfaceScalars ( const char* ifnScalars, 
                            const char* isLabel=NULL,
                            int inFrame=0 );

  // Load a curvature file and display it on the surface in the view.
  void LoadSurfaceCurvatureScalars ( const char* ifnScalars );

  // Load an annotation file. This will create the lookup map used by
  // the view so that it can get a string to display whenever a vertex
  // is clicked, and also load the annotation as the surface overlay
  // and set it in the view.
  void LoadAnnotation ( const char* ifnScalars );

  // Load a surface scalar file and associated color map and display
  // it in the view as a translucent overlay over the
  // surface. Normally we'll get a lookup map in the scalar (.annot)
  // file itself, otherwise we need a color table file name.
  void LoadSurfaceOverlayScalars ( const char* ifnScalars,
				   const char* ifnColors=NULL );

  // Use a vtkFSSurfaceLabelSource to load a label and pass the
  // polydata to the view.
  void LoadLabel ( const char* ifnLabel );
  
  // Save the current label.
  void SaveLabel ( const char* ifnLabel );

  // Save the label to all subjects, mapped to each subject's space.
  void MapLabelToSubjects ( const char* ifnLabel );

  // Save the view's contents as a TIFF.
  void SaveTIFFImage ( const char* ifnTIFF, int iMagnificationLevel = 1 );

  // Set the current scalars index and call
  // ComposeSurfaceScalarsAndShow().
  void SetCurrentSurfaceScalars ( int inEntry );

  // Delete all the scalars.
  void ClearSurfaceScalars ();

  // Delete the curvature.
  void ClearCurvature ();

  // These call the similarly named functions in the view.
  void RestoreView ();
  void ZoomBy ( float iFactor );
  void ZoomIn ();
  void ZoomOut ();
  void ShowCursor ( int ibShow );
  void SetShowCursorFromMenu ();

  // Look up the surface associated with this label and if present,
  // set it in the view.
  void SetCurrentSurface ( const char* isLabel );

  // Set the opacity in the overlay. Passes this to the view and sets
  // the slider bar.
  void SetSurfaceOverlayOpacity ( double iOpacity );

  // Callback from the table, will call SetCurrentSurfaceScalars.
  void SetCurrentSurfaceScalarsFromTableSelection ();
  
  // Callback from the factors list boxes, these will call
  // ManageFactorListBoxSelections and pass the proper list box and
  // either maDiscreteFactorSelection or maContinuousFactorSelection.
  void DiscreteFactorsListBoxCallback ();
  void ContinuousFactorsListBoxCallback ();

  // Makes sure that a list box only has two items selected, and that
  // they are reflected properly within the given array.
  void ManageFactorListBoxSelections ( vtkKWListBox* iListBox,
                                       int iaSelections[2] );

  // Called when the user selects a continuous factor for the scatter plot.
  void ScatterPlotListBoxCallback ();

  // Take the Design form input, extract a design, and run glm.
  void AnalyzeDesign ();

  // Callbacks from entries.
  void SetSubjectsDir ( const char* isSubjectsDir );
  void SetAverageSubject ( const char* isAverageSubject );
  void SetDesignName ( const char* isDesignName );
  void SetShowCurvature ( int ibShow );
  void SetDrawCurvatureGreenRed ( int ibDraw );
  void SetSurfaceScalarsColorMin ( const char* isMin );
  void SetSurfaceScalarsColorMid ( const char* isMid );
  void SetSurfaceScalarsColorMax ( const char* isMax );
  void SetSurfaceScalarsColorReverse ( int ibReverse );
  void SetSurfaceScalarsColorShowPositive ( int ibShow );
  void SetSurfaceScalarsColorShowNegative ( int ibShow );

  // Called by the scalars window when an entry needs to be
  // completed. We make a label widget from the cell text.
  void CreateScalarTableEntry ( const char* iTable, int iRow, int iColumn,
				const char* iWidget );

  // Called when the scalar color table editor changes. This sets the
  // min/max/mid values if they are different.
  void SurfaceScalarColorsEditorChanged ();

  // Sets the color scale values, updates the entries and editor, and
  // calls ComposeSurfaceScalarsAndShow().
  void SetSurfaceScalarsColorMin ( double iMin );
  void SetSurfaceScalarsColorMid ( double iMid );
  void SetSurfaceScalarsColorMax ( double iMax );
  void SetSurfaceScalarsColors ( double iMin, double iMid, double iMax );
  void SetSurfaceScalarsColorOffset ( double iOffset );

  // Calculate the min/mid/max using mSurfaceScalarsColorsFDRRate.
  void SetSurfaceScalarsColorsUsingFDR ();

  // Callback from the FDR rate entry, just sets our value.
  void SetSurfaceScalarsColorsFDRRate ( const char* isValue );

  // Calls the view's function.
  void SelectSurfaceVertex ( int inVertex );
  
  // Add or remove the selection from the ROI. This takes the list of
  // points from view->GetSurfaceSelectionPoints(), uses them as input
  // to a vtkSelectPolyData, and uses the resulting scalars to build a
  // list of vertices to add or remove from the
  // vtkFSSurfaceLabelSource.
  void AddSelectionToROI ();
  void RemoveSelectionFromROI ();
  void ClearROI ();

  // Graph the average of the ROI in the GDF plot.
  void GraphAverageROIInGDF ();
  
  // Smooth the curvature or current scalar values.
  void SmoothCurvatureScalars ( int icSteps );
  void SmoothSurfaceScalars ( int icSteps );

  // This is called by vtkQdecNotebookObserver when a page is raised.
  void NotebookPageRaised ( const char* isTitle );

  // Callback for the notebook.
  static void NotebookRaisePageCallback
    ( vtkObject*, unsigned long, void*, void* );


  // This is called by vtkQdecBltGraphObserver when an element is
  // moused over or not moused over.
  void ScatterPlotGraphMouseoverEnterElement ( const char* isElement );
  void ScatterPlotGraphMouseoverExitElement ();

  // This is called by vtkQdecBltGraphObserver when an element is
  // right-clicked and a contextual menu is about to pop up.
  void ScatterPlotGraphSetUpContextualMenu ( const char* isElement,
						vtkKWMenu* iMenu );

  // Callbacks for the continuous element graph.
  static void ScatterPlotGraphMouseoverEnterElementCallback
    ( vtkObject*, unsigned long, void*, void* );
  static void ScatterPlotGraphMouseoverExitElementCallback
    ( vtkObject*, unsigned long, void*, void* );
  static void ScatterPlotGraphContextualMenuOpeningCallback
    ( vtkObject*, unsigned long, void*, void* );
  

  // This is called by the contextual menus associated with continuous
  // factor plot elements. It affects whether or not a subject is used
  // in the QdecGlmDesign.
  void SetExcludeSubjectID ( const char* isElement, int ibExclude );
  void SetExcludeSubjectGT ( double inExcludeGT );
  void SetExcludeSubjectGT ( const char* isExcludeGT );
  void SetExcludeSubjectLT ( double inExcludeLT );
  void SetExcludeSubjectLT ( const char* isExcludeLT );
  void ClearAllExcludedSubjects ( );

  // Implements QdecVertexAnnotationLookup, which the view will use to
  // get an annotation string for a vertex number. We use this to
  // return a region string.
  const char* GetAnnotationForVertex ( int inVertex );


 protected:

  vtkKWQdecWindow ();
  virtual ~vtkKWQdecWindow ();

  // The names of our notebook pages.
  static const char* ksSubjectsPanelName;
  static const char* ksDesignPanelName;
  static const char* ksDisplayPanelName;
  const char* mCurrentNotebookPanelName;

  // Compose the scalars and color table and display it in the view.
  void ComposeSurfaceScalarsAndShow ();

  // Enable or disable buttons and menu items based on program state.
  void UpdateCommandStatus ();

  // Update our Subjects tab menus with data from the data table.
  void UpdateSubjectsPage ();
  void UpdateNumberOfSubjects ();

  // Update our Design tab menus with data from the project.
  void UpdateDesignPage ();

  // Update our Display tab with the current scalars and other view
  // settings. Will also pack the annotation and surface frame with
  // options if that data is loaded.
  void UpdateDisplayPage ();

  // Update the colors editor with the current surface scalars or a
  // placeholder one, and arrange the color points for min/mid/max.
  void UpdateSurfaceScalarsColorsEditor();

  // Update the scatter plot of continuous factors in the subject
  // panel.
  void UpdateScatterPlot ();

  //BTX
  // Use MRISreadValuesIntoArray to load in a file and initialize a
  // vtkFloatArray.
  vtkFloatArray*
    NewScalarsFromSurfaceScalarsFile ( char const* ifn, int inIndex=0 );

  // Our callback for UserSelectedVertex events. Calls the member
  // function with the SurfaceVertexInformation structure.
  static void UserSelectedVertexCallback ( vtkObject* iCaller, 
					   unsigned long iEventId,
					   void* iClientData,
					   void* iCallData );
  
  // Called by our callback when a UserSelectedVertex event is fired.
  void UserSelectedVertex ( vtkKWQdecView::SurfaceVertexInformation& iInfo );


  // These functions make our toolbar creation code a little
  // nicer. They create a button of the specified type, set its
  // balloon text, command, and icon if any of those are not NULL, and
  // adds it to the toolbar. It returns a pointer to the button so it
  // can be enabled/disabled later.
  vtkKWPushButton* MakeToolbarButton ( vtkKWToolbar* iToolbar,
				       const char* isBalloonHelpText,
				       vtkObject* iCommandObject,
				       const char* isCommand,
				       const char* isIconKey );

  vtkKWCheckButton* MakeToolbarCheckButton ( vtkKWToolbar* iToolbar,
					     const char* isBalloonHelpText,
					     vtkObject* iCommandObject,
					     const char* isCommand,
					     const char* isIconKey );

  // Inserts a 5 pixel wide spacer into the toolbar.
  enum { kToolbarSpacerWidth = 5 };
  void AddSpacerToToolbar ( vtkKWToolbar* iToolbar, 
			    int iWidth = kToolbarSpacerWidth );

  // Mini class for creating a menu item and enabling/disabling
  // it. Parameters are self explanatory.
  class MenuItem {
  public:
    MenuItem ();

    // Creates a simple menu command.
    void MakeCommand ( vtkKWMenu* iMenu, int inItem, const char* isText,
		       vtkObject* iCommandObject, const char* isCommand, 
		       const char* isAccelerator, const char* isIconKey );

    // Create a check box item.
    void MakeCheckButton ( vtkKWMenu* iMenu, int inItem, const char* isText,
			   vtkObject* iCommandObject, const char* isCommand, 
			   const char* isAccelerator, const char* isIconKey );
    
    // Disable or enable the menu item.
    void SetStateToDisabled ();
    void SetStateToNormal ();

    // Set/get the check value if it's a check button.
    void SetSelectedState ( int ibOn );
    int GetSelectedState ();
    
    // Get the index of the menu item.
    int GetIndex () const;

  protected:
    vtkSmartPointer<vtkKWMenu> mMenu;
    int mnItem ;
  };

  // Return a list of vertex indices that are inside the closed path
  // made by the input point list.
  void FindVerticesInsidePath ( vtkPoints* iPath, vtkIntArray* ioVertices );

  // Whether or not to use the histogram editor. This seems to cause
  // problems on 64 bit machines.
  bool mbUseHistogramEditor;

  // Our main view, visible when the Display panel is active.
  vtkSmartPointer<vtkKWQdecView> mView;

  // Our graph view, visible when the Subjects panel is active.
  vtkSmartPointer<vtkKWBltGraph> mGraph;

  // The menu items associated with each command.
  MenuItem* mMenuLoadDataTable;
  MenuItem* mMenuLoadProjectFile;
  MenuItem* mMenuLoadLabel;
  MenuItem* mMenuLoadAnnotation;
  MenuItem* mMenuSaveProjectFile;
  MenuItem* mMenuSaveScatterPlotPostscript;
  MenuItem* mMenuSaveTIFF;
  MenuItem* mMenuSaveGDFPostscript;
  MenuItem* mMenuSaveLabel;
  MenuItem* mMenuMapLabel;
  MenuItem* mMenuClearCurvature;
  MenuItem* mMenuClearSurfaceScalars;
  MenuItem* mMenuRestoreView;
  MenuItem* mMenuZoomOut;
  MenuItem* mMenuZoomIn;
  MenuItem* mMenuShowCursor;
  MenuItem* mMenuAddSelectionToROI;
  MenuItem* mMenuRemoveSelectionFromROI;
  MenuItem* mMenuClearROI;
  MenuItem* mMenuSmoothCurvatureScalars;
  MenuItem* mMenuSmoothSurfaceScalars;
  MenuItem* mMenuGraphAverageROI;
  
  // The toolbar button associated with each command.
  vtkSmartPointer<vtkKWPushButton>  mBtnLoadDataTable;
  vtkSmartPointer<vtkKWPushButton>  mBtnLoadProjectFile;
  vtkSmartPointer<vtkKWPushButton>  mBtnLoadLabel;
  vtkSmartPointer<vtkKWPushButton>  mBtnSaveTIFF;
  vtkSmartPointer<vtkKWPushButton>  mBtnSaveLabel;
  vtkSmartPointer<vtkKWPushButton>  mBtnRestoreView;
  vtkSmartPointer<vtkKWPushButton>  mBtnZoomOut;
  vtkSmartPointer<vtkKWPushButton>  mBtnZoomIn;
  vtkSmartPointer<vtkKWPushButton>  mBtnCameraElevatePositive;
  vtkSmartPointer<vtkKWPushButton>  mBtnCameraElevateNegative;
  vtkSmartPointer<vtkKWPushButton>  mBtnCameraAzimuthPositive;
  vtkSmartPointer<vtkKWPushButton>  mBtnCameraAzimuthNegative;
  vtkSmartPointer<vtkKWPushButton>  mBtnCameraRollPositive;
  vtkSmartPointer<vtkKWPushButton>  mBtnCameraRollNegative;
  vtkSmartPointer<vtkKWCheckButton> mBtnShowCursor;
  vtkSmartPointer<vtkKWEntry>       mEntrySelectVertex;
  vtkSmartPointer<vtkKWPushButton>  mBtnAddSelectionToROI;
  vtkSmartPointer<vtkKWPushButton>  mBtnRemoveSelectionFromROI;
  
  // GUI panel.
  vtkSmartPointer<vtkKWUserInterfacePanel> mPanel;

  // Widgets in the Subjects panel.
  vtkSmartPointer<vtkKWEntry>       mEntrySubjectsDir;
  vtkSmartPointer<vtkKWEntry>       mEntryAverageSubject;
  vtkSmartPointer<vtkKWEntry>       mEntryDataTable;
  vtkSmartPointer<vtkKWEntry>       mEntryNumberOfSubjects;
  vtkSmartPointer<vtkKWListBox>     mListScatterPlot;
  vtkSmartPointer<vtkKWEntry>       mEntryExcludeFactor;
  vtkSmartPointer<vtkKWEntry>       mEntryExcludeSubjectGT;
  vtkSmartPointer<vtkKWEntry>       mEntryExcludeSubjectLT;

  // Widgets in the Design panel.
  vtkSmartPointer<vtkKWListBox>     mListDiscreteFactors;
  vtkSmartPointer<vtkKWListBox>     mListContinuousFactors;
  vtkSmartPointer<vtkKWMenuButton>  mMenuMeasure;
  vtkSmartPointer<vtkKWMenuButton>  mMenuHemisphere;
  vtkSmartPointer<vtkKWMenuButton>  mMenuSmoothness;
  vtkSmartPointer<vtkKWEntry>       mEntryDesignName;

  // Widgets for the Display panel.
  vtkSmartPointer<vtkKWFrame>           mFrameSurface;
  vtkSmartPointer<vtkKWRadioButtonSet>  mRadBtnSetSurface;
  vtkSmartPointer<vtkKWFrame>           mFrameCurvature;
  vtkSmartPointer<vtkKWCheckButton>     mCheckShowCurvature;
  vtkSmartPointer<vtkKWCheckButton>     mCheckDrawCurvatureGreenRed;
  vtkSmartPointer<vtkKWFrame>           mFrameOverlay;
  vtkSmartPointer<vtkKWScale>           mScaleOverlay;
  vtkSmartPointer<vtkKWFrame>           mFrameSurfaceScalars;
  vtkSmartPointer<vtkKWMultiColumnList> mTableSurfaceScalars;
  vtkSmartPointer<vtkKWCheckButton>     mCheckSurfaceScalarsColorReverse;
  vtkSmartPointer<vtkKWCheckButton>     mCheckSurfaceScalarsColorShowPositive;
  vtkSmartPointer<vtkKWCheckButton>     mCheckSurfaceScalarsColorShowNegative;
  vtkSmartPointer<vtkKWRGBATransferFunctionEditor> 
                                        mEditorSurfaceScalarColors;
  vtkSmartPointer<vtkKWHistogram>       mHistogramSurfaceScalarColors;
  vtkSmartPointer<vtkKWEntry>           mEntrySurfaceScalarsColorMin;
  vtkSmartPointer<vtkKWEntry>           mEntrySurfaceScalarsColorMid;
  vtkSmartPointer<vtkKWEntry>           mEntrySurfaceScalarsColorMax;
  vtkSmartPointer<vtkKWEntry>           mEntrySurfaceScalarsColorOffset;

  // Use these variables to see if we've packed the mFrame*
  // placeholder frames declared above. This lets us dynamically fill
  // stuff in the Display pane.
  bool mbFrameSurfaceInited;
  bool mbFrameCurvatureInited;
  bool mbFrameOverlayInited;
  bool mbFrameSurfaceScalarsInited;

  // Keep track of the order in which the items in the factor
  // listboxes were selected, so we can unselect the first one if a
  // third one is selected.
  int maDiscreteFactorSelection[2];
  int maContinuousFactorSelection[2];
  
  // The factor selected to plot in the scatter plot from
  // the Subjects panel.
  int mScatterPlotSelection;

  // The struct for a scalars object.
  typedef struct {
    int mnEntry;
    vtkSmartPointer<vtkFloatArray> mValues;
    std::string mfnSource;
    std::string msLabel;
  } SurfaceScalar;

  // Data objects.
  QdecProject* mQdecProject;
  int mcVertices;
  std::map<std::string,vtkSmartPointer<vtkFSSurfaceSource> > maSurfaceSource;
  std::string msCurrentSurfaceSource;
  int mGDFID;
  int mbGDFLoaded;
  vtkSmartPointer<vtkFloatArray> mCurvatureScalars;
  std::map<int,SurfaceScalar> maSurfaceScalars;
  int mnCurrentSurfaceScalars;
  int* maAnnotationIndicies;
  COLOR_TABLE* mAnnotationTable;
  vtkSmartPointer<vtkFloatArray> mOverlayScalars;
  vtkSmartPointer<vtkScalarsToColors> mOverlayColors;
  vtkSmartPointer<vtkFSSurfaceLabelSource> mROISource;
  
  // Indices into the scalars colors editor for our color points.
  int mnNegativeMinValue, mnNegativeMidValue, mnNegativeMaxValue;
  int mnPositiveMinValue, mnPositiveMidValue, mnPositiveMaxValue;

  // For scalars color scale drawing.
  bool mbShowCurvature;
  double mSurfaceScalarsColorMin;
  double mSurfaceScalarsColorMid;
  double mSurfaceScalarsColorMax;
  double mSurfaceScalarsColorOffset;
  bool mbSurfaceScalarsColorReverse;
  bool mbSurfaceScalarsColorShowPositive;
  bool mbSurfaceScalarsColorShowNegative;

  // If true, this will draw the curvature with green/red when there
  // is no scalar up. Otherwise, it will be drawn in binary gray, as
  // normal.
  bool mbDrawCurvatureGreenRedIfNoScalars;

  // For calculating the FDR rate.
  double mSurfaceScalarsColorsFDRRate;

  // A label for the overlay scale bar.
  std::string msOverlayDescription;

  // Initially false, until the view has been initialized by being
  // reset.
  bool mbViewInitialized;

  //ETX
};

#endif
