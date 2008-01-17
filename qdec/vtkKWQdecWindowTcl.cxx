// tcl wrapper for vtkKWQdecWindow object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkKWQdecWindow.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>

ClientData vtkKWQdecWindowNewCommand()
{
  vtkKWQdecWindow *temp = vtkKWQdecWindow::New();
  return ((ClientData)temp);
}

int vtkKWWindowCppCommand(vtkKWWindow *op, Tcl_Interp *interp,
             int argc, char *argv[]);
int VTKTCL_EXPORT vtkKWQdecWindowCppCommand(vtkKWQdecWindow *op, Tcl_Interp *interp,
             int argc, char *argv[]);

int VTKTCL_EXPORT vtkKWQdecWindowCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[])
{
  if ((argc == 2)&&(!strcmp("Delete",argv[1]))&& !vtkTclInDelete(interp))
    {
    Tcl_DeleteCommand(interp,argv[0]);
    return TCL_OK;
    }
   return vtkKWQdecWindowCppCommand((vtkKWQdecWindow *)(((vtkTclCommandArgStruct *)cd)->Pointer),interp, argc, argv);
}

int VTKTCL_EXPORT vtkKWQdecWindowCppCommand(vtkKWQdecWindow *op, Tcl_Interp *interp,
             int argc, char *argv[])
{
  int    tempi;
  double tempd;
  static char temps[80];
  int    error;

  error = 0; error = error;
  tempi = 0; tempi = tempi;
  tempd = 0; tempd = tempd;
  temps[0] = 0; temps[0] = temps[0];

  if (argc < 2)
    {
    Tcl_SetResult(interp, (char *) "Could not find requested method.", TCL_VOLATILE);
    return TCL_ERROR;
    }
  if (!interp)
    {
    if (!strcmp("DoTypecasting",argv[0]))
      {
      if (!strcmp("vtkKWQdecWindow",argv[1]))
        {
        argv[2] = (char *)((void *)op);
        return TCL_OK;
        }
      if (vtkKWWindowCppCommand((vtkKWWindow *)op,interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,(char *) "vtkKWWindow", TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkKWQdecWindow  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecWindow");
    return TCL_OK;
    }
  if ((!strcmp("GetClassName",argv[1]))&&(argc == 2))
    {
    const char    *temp20;
    temp20 = (op)->GetClassName();
    if (temp20)
      {
      Tcl_SetResult(interp, (char*)temp20, TCL_VOLATILE);
      }
    else
      {
      Tcl_ResetResult(interp);
      }
    return TCL_OK;
    }
  if ((!strcmp("IsA",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    int      temp20;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    temp20 = (op)->IsA(temp0);
    char tempResult[1024];
    sprintf(tempResult,"%i",temp20);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
    return TCL_OK;
    }
    }
  if ((!strcmp("NewInstance",argv[1]))&&(argc == 2))
    {
    vtkKWQdecWindow  *temp20;
    temp20 = (op)->NewInstance();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecWindow");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkKWQdecWindow  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkObject",interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecWindow");
    return TCL_OK;
    }
    }
  if ((!strcmp("SetUseHistogramEditor",argv[1]))&&(argc == 3))
    {
    bool   temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi ? true : false;
    if (!error)
    {
    op->SetUseHistogramEditor(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("CreateWidget",argv[1]))&&(argc == 2))
    {
    op->CreateWidget();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("FinishCreating",argv[1]))&&(argc == 2))
    {
    op->FinishCreating();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("BeginActionWithProgress",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->BeginActionWithProgress(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("UpdateProgressMessage",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->UpdateProgressMessage(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("UpdateProgressPercent",argv[1]))&&(argc == 3))
    {
    float    temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->UpdateProgressPercent(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("EndActionWithProgress",argv[1]))&&(argc == 2))
    {
    op->EndActionWithProgress();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadDataTableFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadDataTableFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadProjectFileFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadProjectFileFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadSurfaceFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadSurfaceFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadGDFFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadGDFFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadSurfaceScalarsFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadSurfaceScalarsFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadCurvatureFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadCurvatureFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadAnnotationFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadAnnotationFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadLabelFromDlog",argv[1]))&&(argc == 2))
    {
    op->LoadLabelFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SaveProjectFileFromDlog",argv[1]))&&(argc == 2))
    {
    op->SaveProjectFileFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SaveScatterPlotPostscriptFromDlog",argv[1]))&&(argc == 2))
    {
    op->SaveScatterPlotPostscriptFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SaveTIFFImageFromDlog",argv[1]))&&(argc == 2))
    {
    op->SaveTIFFImageFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SaveGDFPostscriptFromDlog",argv[1]))&&(argc == 2))
    {
    op->SaveGDFPostscriptFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SaveLabelFromDlog",argv[1]))&&(argc == 2))
    {
    op->SaveLabelFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("MapLabelFromDlog",argv[1]))&&(argc == 2))
    {
    op->MapLabelFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SmoothCurvatureScalarsFromDlog",argv[1]))&&(argc == 2))
    {
    op->SmoothCurvatureScalarsFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SmoothSurfaceScalarsFromDlog",argv[1]))&&(argc == 2))
    {
    op->SmoothSurfaceScalarsFromDlog();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("LoadDataTable",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadDataTable(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadProjectFile",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadProjectFile(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SaveProjectFile",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SaveProjectFile(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadSurface",argv[1]))&&(argc == 4))
    {
    char    *temp0;
    char    *temp1;
    error = 0;

    temp0 = argv[2];
    temp1 = argv[3];
    if (!error)
    {
    op->LoadSurface(temp0,temp1);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadGDFFile",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadGDFFile(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadSurfaceScalars",argv[1]))&&(argc == 5))
    {
    char    *temp0;
    char    *temp1;
    int      temp2;
    int      temp20;
    error = 0;

    temp0 = argv[2];
    temp1 = argv[3];
    if (Tcl_GetInt(interp,argv[4],&tempi) != TCL_OK) error = 1;
    temp2 = tempi;
    if (!error)
    {
    temp20 = (op)->LoadSurfaceScalars(temp0,temp1,temp2);
    char tempResult[1024];
    sprintf(tempResult,"%i",temp20);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadSurfaceCurvatureScalars",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadSurfaceCurvatureScalars(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadAnnotation",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadAnnotation(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadSurfaceOverlayScalars",argv[1]))&&(argc == 4))
    {
    char    *temp0;
    char    *temp1;
    error = 0;

    temp0 = argv[2];
    temp1 = argv[3];
    if (!error)
    {
    op->LoadSurfaceOverlayScalars(temp0,temp1);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("LoadLabel",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadLabel(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SaveLabel",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SaveLabel(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("MapLabelToSubjects",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->MapLabelToSubjects(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SaveTIFFImage",argv[1]))&&(argc == 4))
    {
    char    *temp0;
    int      temp1;
    error = 0;

    temp0 = argv[2];
    if (Tcl_GetInt(interp,argv[3],&tempi) != TCL_OK) error = 1;
    temp1 = tempi;
    if (!error)
    {
    op->SaveTIFFImage(temp0,temp1);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetCurrentSurfaceScalars",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetCurrentSurfaceScalars(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("ClearSurfaceScalars",argv[1]))&&(argc == 2))
    {
    op->ClearSurfaceScalars();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ClearCurvature",argv[1]))&&(argc == 2))
    {
    op->ClearCurvature();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("RestoreView",argv[1]))&&(argc == 2))
    {
    op->RestoreView();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ZoomBy",argv[1]))&&(argc == 3))
    {
    float    temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->ZoomBy(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("ZoomIn",argv[1]))&&(argc == 2))
    {
    op->ZoomIn();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ZoomOut",argv[1]))&&(argc == 2))
    {
    op->ZoomOut();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ShowCursor",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->ShowCursor(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetShowCursorFromMenu",argv[1]))&&(argc == 2))
    {
    op->SetShowCursorFromMenu();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SetCurrentSurface",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetCurrentSurface(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceOverlayOpacity",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetSurfaceOverlayOpacity(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetCurrentSurfaceScalarsFromTableSelection",argv[1]))&&(argc == 2))
    {
    op->SetCurrentSurfaceScalarsFromTableSelection();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("DiscreteFactorsListBoxCallback",argv[1]))&&(argc == 2))
    {
    op->DiscreteFactorsListBoxCallback();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ContinuousFactorsListBoxCallback",argv[1]))&&(argc == 2))
    {
    op->ContinuousFactorsListBoxCallback();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ScatterPlotListBoxCallback",argv[1]))&&(argc == 2))
    {
    op->ScatterPlotListBoxCallback();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("AnalyzeDesign",argv[1]))&&(argc == 2))
    {
    op->AnalyzeDesign();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SetSubjectsDir",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetSubjectsDir(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetAverageSubject",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetAverageSubject(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetDesignName",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetDesignName(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetShowCurvature",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetShowCurvature(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetDrawCurvatureGreenRed",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetDrawCurvatureGreenRed(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorMin",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetSurfaceScalarsColorMin(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorMid",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetSurfaceScalarsColorMid(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorMax",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetSurfaceScalarsColorMax(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorReverse",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetSurfaceScalarsColorReverse(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorShowPositive",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetSurfaceScalarsColorShowPositive(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorShowNegative",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetSurfaceScalarsColorShowNegative(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("CreateScalarTableEntry",argv[1]))&&(argc == 6))
    {
    char    *temp0;
    int      temp1;
    int      temp2;
    char    *temp3;
    error = 0;

    temp0 = argv[2];
    if (Tcl_GetInt(interp,argv[3],&tempi) != TCL_OK) error = 1;
    temp1 = tempi;
    if (Tcl_GetInt(interp,argv[4],&tempi) != TCL_OK) error = 1;
    temp2 = tempi;
    temp3 = argv[5];
    if (!error)
    {
    op->CreateScalarTableEntry(temp0,temp1,temp2,temp3);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SurfaceScalarColorsEditorChanged",argv[1]))&&(argc == 2))
    {
    op->SurfaceScalarColorsEditorChanged();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SetSurfaceScalarsColorMin",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetSurfaceScalarsColorMin(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorMid",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetSurfaceScalarsColorMid(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorMax",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetSurfaceScalarsColorMax(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColors",argv[1]))&&(argc == 5))
    {
    double   temp0;
    double   temp1;
    double   temp2;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (Tcl_GetDouble(interp,argv[3],&tempd) != TCL_OK) error = 1;
    temp1 = tempd;
    if (Tcl_GetDouble(interp,argv[4],&tempd) != TCL_OK) error = 1;
    temp2 = tempd;
    if (!error)
    {
    op->SetSurfaceScalarsColors(temp0,temp1,temp2);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorOffset",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetSurfaceScalarsColorOffset(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColorsUsingFDR",argv[1]))&&(argc == 2))
    {
    op->SetSurfaceScalarsColorsUsingFDR();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SetSurfaceScalarsColorsFDRRate",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetSurfaceScalarsColorsFDRRate(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SelectSurfaceVertex",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SelectSurfaceVertex(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("AddSelectionToROI",argv[1]))&&(argc == 2))
    {
    op->AddSelectionToROI();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("RemoveSelectionFromROI",argv[1]))&&(argc == 2))
    {
    op->RemoveSelectionFromROI();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ClearROI",argv[1]))&&(argc == 2))
    {
    op->ClearROI();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("GraphAverageROIInGDF",argv[1]))&&(argc == 2))
    {
    op->GraphAverageROIInGDF();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SmoothCurvatureScalars",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SmoothCurvatureScalars(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SmoothSurfaceScalars",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SmoothSurfaceScalars(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("NotebookPageRaised",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->NotebookPageRaised(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("ScatterPlotGraphMouseoverEnterElement",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->ScatterPlotGraphMouseoverEnterElement(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("ScatterPlotGraphMouseoverExitElement",argv[1]))&&(argc == 2))
    {
    op->ScatterPlotGraphMouseoverExitElement();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ScatterPlotGraphSetUpContextualMenu",argv[1]))&&(argc == 4))
    {
    char    *temp0;
    vtkKWMenu  *temp1;
    error = 0;

    temp0 = argv[2];
    temp1 = (vtkKWMenu *)(vtkTclGetPointerFromObject(argv[3],(char *) "vtkKWMenu",interp,error));
    if (!error)
    {
    op->ScatterPlotGraphSetUpContextualMenu(temp0,temp1);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetExcludeSubjectID",argv[1]))&&(argc == 4))
    {
    char    *temp0;
    int      temp1;
    error = 0;

    temp0 = argv[2];
    if (Tcl_GetInt(interp,argv[3],&tempi) != TCL_OK) error = 1;
    temp1 = tempi;
    if (!error)
    {
    op->SetExcludeSubjectID(temp0,temp1);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetExcludeSubjectGT",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetExcludeSubjectGT(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetExcludeSubjectGT",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetExcludeSubjectGT(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetExcludeSubjectLT",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetExcludeSubjectLT(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetExcludeSubjectLT",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetExcludeSubjectLT(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("ClearAllExcludedSubjects",argv[1]))&&(argc == 2))
    {
    op->ClearAllExcludedSubjects();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("GetAnnotationForVertex",argv[1]))&&(argc == 3))
    {
    int      temp0;
    const char    *temp20;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    temp20 = (op)->GetAnnotationForVertex(temp0);
    if (temp20)
      {
      Tcl_SetResult(interp, (char*)temp20, TCL_VOLATILE);
      }
    else
      {
      Tcl_ResetResult(interp);
      }
    return TCL_OK;
    }
    }

  if (!strcmp("ListInstances",argv[1]))
    {
    vtkTclListInstances(interp,(ClientData)vtkKWQdecWindowCommand);
    return TCL_OK;
    }

  if (!strcmp("ListMethods",argv[1]))
    {
    vtkKWWindowCppCommand(op,interp,argc,argv);
    Tcl_AppendResult(interp,"Methods from vtkKWQdecWindow:\n",NULL);
    Tcl_AppendResult(interp,"  GetSuperClassName\n",NULL);
    Tcl_AppendResult(interp,"  New\n",NULL);
    Tcl_AppendResult(interp,"  GetClassName\n",NULL);
    Tcl_AppendResult(interp,"  IsA\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  NewInstance\n",NULL);
    Tcl_AppendResult(interp,"  SafeDownCast\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetUseHistogramEditor\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  CreateWidget\n",NULL);
    Tcl_AppendResult(interp,"  FinishCreating\n",NULL);
    Tcl_AppendResult(interp,"  BeginActionWithProgress\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  UpdateProgressMessage\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  UpdateProgressPercent\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  EndActionWithProgress\n",NULL);
    Tcl_AppendResult(interp,"  LoadDataTableFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadProjectFileFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadGDFFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceScalarsFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadCurvatureFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadAnnotationFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadLabelFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveProjectFileFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveScatterPlotPostscriptFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveTIFFImageFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveGDFPostscriptFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveLabelFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  MapLabelFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SmoothCurvatureScalarsFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SmoothSurfaceScalarsFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  LoadDataTable\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadProjectFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SaveProjectFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurface\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  LoadGDFFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceScalars\t with 3 args\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceCurvatureScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadAnnotation\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceOverlayScalars\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  LoadLabel\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SaveLabel\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  MapLabelToSubjects\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SaveTIFFImage\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentSurfaceScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ClearSurfaceScalars\n",NULL);
    Tcl_AppendResult(interp,"  ClearCurvature\n",NULL);
    Tcl_AppendResult(interp,"  RestoreView\n",NULL);
    Tcl_AppendResult(interp,"  ZoomBy\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ZoomIn\n",NULL);
    Tcl_AppendResult(interp,"  ZoomOut\n",NULL);
    Tcl_AppendResult(interp,"  ShowCursor\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetShowCursorFromMenu\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentSurface\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceOverlayOpacity\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentSurfaceScalarsFromTableSelection\n",NULL);
    Tcl_AppendResult(interp,"  DiscreteFactorsListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  ContinuousFactorsListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  ScatterPlotListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  AnalyzeDesign\n",NULL);
    Tcl_AppendResult(interp,"  SetSubjectsDir\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetAverageSubject\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetDesignName\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetShowCurvature\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetDrawCurvatureGreenRed\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMin\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMid\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMax\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorReverse\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorShowPositive\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorShowNegative\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  CreateScalarTableEntry\t with 4 args\n",NULL);
    Tcl_AppendResult(interp,"  SurfaceScalarColorsEditorChanged\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMin\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMid\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMax\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColors\t with 3 args\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorOffset\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorsUsingFDR\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorsFDRRate\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SelectSurfaceVertex\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  AddSelectionToROI\n",NULL);
    Tcl_AppendResult(interp,"  RemoveSelectionFromROI\n",NULL);
    Tcl_AppendResult(interp,"  ClearROI\n",NULL);
    Tcl_AppendResult(interp,"  GraphAverageROIInGDF\n",NULL);
    Tcl_AppendResult(interp,"  SmoothCurvatureScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SmoothSurfaceScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  NotebookPageRaised\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ScatterPlotGraphMouseoverEnterElement\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ScatterPlotGraphMouseoverExitElement\n",NULL);
    Tcl_AppendResult(interp,"  ScatterPlotGraphSetUpContextualMenu\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  SetExcludeSubjectID\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  SetExcludeSubjectGT\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetExcludeSubjectGT\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetExcludeSubjectLT\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetExcludeSubjectLT\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ClearAllExcludedSubjects\n",NULL);
    Tcl_AppendResult(interp,"  GetAnnotationForVertex\t with 1 arg\n",NULL);
    return TCL_OK;
    }

  if (vtkKWWindowCppCommand((vtkKWWindow *)op,interp,argc,argv) == TCL_OK)
    {
    return TCL_OK;
    }

  if ((argc >= 2)&&(!strstr(interp->result,"Object named:")))
    {
    char temps2[256];
    sprintf(temps2,"Object named: %s, could not find requested method: %s\nor the method was called with incorrect arguments.\n",argv[0],argv[1]);
    Tcl_AppendResult(interp,temps2,NULL);
    }
    }
  catch (vtkstd::exception &e)
    {
    Tcl_AppendResult(interp, "Uncaught exception: ",  e.what(), "\n", NULL);
    return TCL_ERROR;
    }
  return TCL_ERROR;
}
