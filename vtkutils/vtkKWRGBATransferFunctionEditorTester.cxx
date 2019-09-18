#include "vtkKWRGBATransferFunctionEditorTester.h"

#include "mrisurf.h"

#include "vtkObjectFactory.h"
#include "vtkKWWindowBase.h"
#include "vtkKWFrame.h"
#include "vtkKWRGBATransferFunctionEditor.h"
#include "vtkRGBATransferFunction.h"
#include "vtkFloatArray.h"
#include "vtkKWHistogram.h"

extern int Rgbatransferfunctioneditortesterlib_SafeInit ( Tcl_Interp* );

const char* Progname = "vtkKWRGBATransferFunctionEditorTester";

int
main ( int iArgc, char** iArgv ) {

  Tcl_Interp* interp = vtkKWApplication::InitializeTcl( iArgc, iArgv, &cerr );
  if ( !interp ) {
    cerr << "Error initializing Tcl." << endl;
    return 1;
  }

  Rgbatransferfunctioneditortesterlib_SafeInit( interp );

  vtkKWRGBATransferFunctionEditorTester* app = 
    vtkKWRGBATransferFunctionEditorTester::New();

  app->Start( iArgc, iArgv );
  int rApp = app->GetExitStatus();

  app->Delete();

  return rApp;
}

vtkStandardNewMacro( vtkKWRGBATransferFunctionEditorTester );
vtkCxxRevisionMacro( vtkKWRGBATransferFunctionEditorTester, "$Revision: 1.3 $" );

void
vtkKWRGBATransferFunctionEditorTester::Start ( int iArgc, char** iArgv ) {
  
  vtkKWWindowBase* window = vtkKWWindowBase::New();
  window->SetApplication( this );
  this->AddWindow( window );
  window->Create();
  window->SetSize( 400, 300 );

  mEditor = vtkKWRGBATransferFunctionEditor::New();
  mEditor->SetParent( window->GetViewFrame() );
  mEditor->Create();
  
  mEditor->ExpandCanvasWidthOn();
  //      editor->SetCanvasWidth( 450 );
  mEditor->SetCanvasHeight( 150 );
  mEditor->SetLabelText("Editor");
  mEditor->SetRangeLabelPositionToTop();
  mEditor->ColorSpaceOptionMenuVisibilityOff();
  mEditor->ValueEntriesVisibilityOff();
  mEditor->SetPointPositionInValueRangeToTop();
  mEditor->SetPointStyleToCursorDown();
  mEditor->FunctionLineVisibilityOff();
  mEditor->PointGuidelineVisibilityOn();
  mEditor->PointIndexVisibilityOff();
  mEditor->SelectedPointIndexVisibilityOn();
  mEditor->MidPointEntryVisibilityOff();
  mEditor->SharpnessEntryVisibilityOff();
  mEditor->SetLabelPositionToTop();
  mEditor->ColorRampVisibilityOff();
  mEditor->ParameterTicksVisibilityOn();
  mEditor->ComputeValueTicksFromHistogramOn();
  mEditor->SetParameterTicksFormat("%-#6.0f");
  mEditor->SetFunctionChangingCommand( this, "EditorChangedFunction" );
  mEditor->SetFunctionChangedCommand( this, "EditorChangedFunction" );

  vtkRGBATransferFunction* colors = vtkRGBATransferFunction::New();
  mEditor->SetRGBATransferFunction( colors );


  const int cValues = 163842;
  float* aValues = NULL;
  int eRead = MRISreadValuesIntoArray( "vtkKWRGBATransferFunctionEditorTester-scalars.mgh", cValues, &aValues );
  if( 0 != eRead ) {
    if( aValues ) free( aValues );
    cerr << "Could not read scalar file" << endl;
    this->SetExitStatus( 1 );
    return;
  }

  vtkFloatArray* values = vtkFloatArray::New();

  values->Allocate( cValues );
  values->SetNumberOfComponents( 1 );

  for( int nValue = 0; nValue < cValues; nValue ++ )
    values->InsertNextValue( aValues[nValue] );

  free( aValues );

  vtkKWHistogram* histogram = vtkKWHistogram::New();
  histogram->BuildHistogram( values, 0 );
  mEditor->SetHistogram( histogram );

  double range[2];
  values->GetRange( range );
  mEditor->SetWholeParameterRange( range[0], range[1] );
  mEditor->SetVisibleParameterRangeToWholeParameterRange();

#if 0
  colors->AddRGBAPoint( -3, 0, 1, 1, 1 );
  colors->AddRGBAPoint( -2, 0, 0, 1, 1 );
  colors->AddRGBAPoint( -1, 0, 0, 1, 1 );
  colors->AddRGBAPoint( -0.5, 0.5, 0.5, 0.5, 1 );
  colors->AddRGBAPoint(  0.5, 0.5, 0.5, 0.5, 1 );
  colors->AddRGBAPoint(  1, 1, 0, 0, 1 );
  colors->AddRGBAPoint(  2, 1, 0, 0, 1 );
  colors->AddRGBAPoint(  3, 1, 1, 0, 1 );
  colors->Build();
  
  mEditor->SetPointCountMinimum( 8 );
  mEditor->SetPointCountMaximum( 8 );
  mEditor->SetPointSymmetry( 0, 7 ); // max
  mEditor->SetPointSymmetry( 1, 6 ); // mid
  mEditor->SetPointSymmetry( 2, 5 ); // min
  mEditor->SetPointSticky( 2, 3 ); // -min to gray
  mEditor->SetPointSticky( 4, 5 ); // gray to min
#else
  colors->AddRGBAPoint( -3, 0, 1, 1, 1 );
  colors->AddRGBAPoint( -2, 0, 0, 1, 1 );
  colors->AddRGBAPoint( -1, 0, 0, 1, 1 );
  colors->AddRGBAPoint(  1, 1, 0, 0, 1 );
  colors->AddRGBAPoint(  2, 1, 0, 0, 1 );
  colors->AddRGBAPoint(  3, 1, 1, 0, 1 );
  colors->Build();
  
  mEditor->SetPointCountMinimum( 6 );
  mEditor->SetPointCountMaximum( 6 );
  mEditor->SetPointSymmetry( 0, 5 ); // max
  mEditor->SetPointSymmetry( 1, 4 ); // mid
  mEditor->SetPointSymmetry( 2, 3 ); // min
#endif

  this->Script( "pack %s -side top -fill x -expand y",
		mEditor->GetWidgetName() );

  window->Display();

  this->Superclass::Start( iArgc, iArgv );
}

void
vtkKWRGBATransferFunctionEditorTester::EditorChangedFunction () {

}
