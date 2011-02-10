// tcl wrapper for vtkKWQdecWindow object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkKWQdecWindow.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>
#include <vtksys/ios/sstream>

ClientData vtkKWQdecWindowNewCommand()
{
  vtkKWQdecWindow *temp = vtkKWQdecWindow::New();
  return static_cast<ClientData>(temp);
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
   return vtkKWQdecWindowCppCommand(static_cast<vtkKWQdecWindow *>(static_cast<vtkTclCommandArgStruct *>(cd)->Pointer),interp, argc, argv);
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
    Tcl_SetResult(interp,const_cast<char *>("Could not find requested method."), TCL_VOLATILE);
    return TCL_ERROR;
    }
  if (!interp)
    {
    if (!strcmp("DoTypecasting",argv[0]))
      {
      if (!strcmp("vtkKWQdecWindow",argv[1]))
        {
        argv[2] = static_cast<char *>(static_cast<void *>(op));
        return TCL_OK;
        }
      if (vtkKWWindowCppCommand(static_cast<vtkKWWindow *>(op),interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,const_cast<char *>("vtkKWWindow"), TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkKWQdecWindow  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecWindow");
    return TCL_OK;
    }
  if ((!strcmp("GetClassName",argv[1]))&&(argc == 2))
    {
    const char    *temp20;
    temp20 = (op)->GetClassName();
    if (temp20)
      {
      Tcl_SetResult(interp, const_cast<char *>(temp20), TCL_VOLATILE);
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
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecWindow");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkKWQdecWindow  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkObject"),interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecWindow");
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
  if ((!strcmp("FinishCreating",argv[1]))&&(argc == 2))
    {
    op->FinishCreating();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SetCurrentMeasure",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetCurrentMeasure(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetCurrentSurfaceMeasure",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetCurrentSurfaceMeasure(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetDesignMatrixType",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetDesignMatrixType(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
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
  if ((!strcmp("LoadSurfaceScalars",argv[1]))&&(argc == 6))
    {
    char    *temp0;
    char    *temp1;
    char    *temp2;
    int      temp3;
    int      temp20;
    error = 0;

    temp0 = argv[2];
    temp1 = argv[3];
    temp2 = argv[4];
    if (Tcl_GetInt(interp,argv[5],&tempi) != TCL_OK) error = 1;
    temp3 = tempi;
    if (!error)
    {
    temp20 = (op)->LoadSurfaceScalars(temp0,temp1,temp2,temp3);
    char tempResult[1024];
    sprintf(tempResult,"%i",temp20);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
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
  if ((!strcmp("CreateWidget",argv[1]))&&(argc == 2))
    {
    op->CreateWidget();
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
  if ((!strcmp("SaveDataTableFromDlog",argv[1]))&&(argc == 2))
    {
    op->SaveDataTableFromDlog();
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
  if ((!strcmp("QuickSnapsTIFF",argv[1]))&&(argc == 2))
    {
    op->QuickSnapsTIFF();
    Tcl_ResetResult(interp);
    return TCL_OK;
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
  if ((!strcmp("UnloadSurfaceScalars",argv[1]))&&(argc == 2))
    {
    op->UnloadSurfaceScalars();
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
  if ((!strcmp("ShowCurvature",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->ShowCurvature(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetShowCurvatureFromMenu",argv[1]))&&(argc == 2))
    {
    op->SetShowCurvatureFromMenu();
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
  if ((!strcmp("NuisanceFactorsListBoxCallback",argv[1]))&&(argc == 2))
    {
    op->NuisanceFactorsListBoxCallback();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("NuisanceFactorsListBoxCallback",argv[1]))&&(argc == 2))
    {
    op->NuisanceFactorsListBoxCallback();
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
  if ((!strcmp("GetShowCurvature",argv[1]))&&(argc == 2))
    {
    bool   temp20;
    temp20 = (op)->GetShowCurvature();
    char tempResult[1024];
    sprintf(tempResult,"%i",(int)temp20);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
    return TCL_OK;
    }
  if ((!strcmp("SetShowCurvature",argv[1]))&&(argc == 3))
    {
    bool   temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi ? true : false;
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
  if ((!strcmp("GenerateClusterStats",argv[1]))&&(argc == 2))
    {
    op->GenerateClusterStats();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("GotoNextCluster",argv[1]))&&(argc == 2))
    {
    op->GotoNextCluster();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("GotoPrevCluster",argv[1]))&&(argc == 2))
    {
    op->GotoPrevCluster();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("GotoCluster",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->GotoCluster(temp0);
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
  if ((!strcmp("RunSimulation",argv[1]))&&(argc == 2))
    {
    op->RunSimulation();
    Tcl_ResetResult(interp);
    return TCL_OK;
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
    temp1 = (vtkKWMenu *)(vtkTclGetPointerFromObject(argv[3],const_cast<char *>("vtkKWMenu"),interp,error));
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
  if ((!strcmp("SetExcludeSubjectET",argv[1]))&&(argc == 3))
    {
    double   temp0;
    error = 0;

    if (Tcl_GetDouble(interp,argv[2],&tempd) != TCL_OK) error = 1;
    temp0 = tempd;
    if (!error)
    {
    op->SetExcludeSubjectET(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetExcludeSubjectET",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetExcludeSubjectET(temp0);
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
  if ((!strcmp("ResetStatsImportFrame",argv[1]))&&(argc == 2))
    {
    op->ResetStatsImportFrame();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("GenerateStatsDataTables",argv[1]))&&(argc == 2))
    {
    op->GenerateStatsDataTables();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("SetStatsImportItem",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetStatsImportItem(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("AddStatsToDataTable",argv[1]))&&(argc == 2))
    {
    op->AddStatsToDataTable();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("RemoveFactorFromDataTable",argv[1]))&&(argc == 2))
    {
    op->RemoveFactorFromDataTable();
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
      Tcl_SetResult(interp, const_cast<char *>(temp20), TCL_VOLATILE);
      }
    else
      {
      Tcl_ResetResult(interp);
      }
    return TCL_OK;
    }
    }
  if ((!strcmp("ComposeSurfaceScalarsAndShow",argv[1]))&&(argc == 2))
    {
    op->ComposeSurfaceScalarsAndShow();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }

  if (!strcmp("ListInstances",argv[1]))
    {
    vtkTclListInstances(interp,(ClientData)(vtkKWQdecWindowCommand));
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
    Tcl_AppendResult(interp,"  FinishCreating\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentMeasure\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentSurfaceMeasure\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetDesignMatrixType\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadDataTable\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadProjectFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurface\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  LoadGDFFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadAnnotation\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceScalars\t with 4 args\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceOverlayScalars\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceCurvatureScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadLabel\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  CreateWidget\n",NULL);
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
    Tcl_AppendResult(interp,"  SaveDataTableFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveProjectFileFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveScatterPlotPostscriptFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveTIFFImageFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveGDFPostscriptFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveLabelFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  MapLabelFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SmoothCurvatureScalarsFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SmoothSurfaceScalarsFromDlog\n",NULL);
    Tcl_AppendResult(interp,"  SaveProjectFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SaveLabel\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  MapLabelToSubjects\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SaveTIFFImage\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  QuickSnapsTIFF\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentSurfaceScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ClearSurfaceScalars\n",NULL);
    Tcl_AppendResult(interp,"  UnloadSurfaceScalars\n",NULL);
    Tcl_AppendResult(interp,"  ClearCurvature\n",NULL);
    Tcl_AppendResult(interp,"  RestoreView\n",NULL);
    Tcl_AppendResult(interp,"  ZoomBy\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ZoomIn\n",NULL);
    Tcl_AppendResult(interp,"  ZoomOut\n",NULL);
    Tcl_AppendResult(interp,"  ShowCursor\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetShowCursorFromMenu\n",NULL);
    Tcl_AppendResult(interp,"  ShowCurvature\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetShowCurvatureFromMenu\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentSurface\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceOverlayOpacity\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetCurrentSurfaceScalarsFromTableSelection\n",NULL);
    Tcl_AppendResult(interp,"  DiscreteFactorsListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  ContinuousFactorsListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  NuisanceFactorsListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  NuisanceFactorsListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  ScatterPlotListBoxCallback\n",NULL);
    Tcl_AppendResult(interp,"  AnalyzeDesign\n",NULL);
    Tcl_AppendResult(interp,"  SetSubjectsDir\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetAverageSubject\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetDesignName\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  GetShowCurvature\n",NULL);
    Tcl_AppendResult(interp,"  SetShowCurvature\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetDrawCurvatureGreenRed\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMin\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMid\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMax\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorReverse\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorShowPositive\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorShowNegative\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  GenerateClusterStats\n",NULL);
    Tcl_AppendResult(interp,"  GotoNextCluster\n",NULL);
    Tcl_AppendResult(interp,"  GotoPrevCluster\n",NULL);
    Tcl_AppendResult(interp,"  GotoCluster\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  CreateScalarTableEntry\t with 4 args\n",NULL);
    Tcl_AppendResult(interp,"  SurfaceScalarColorsEditorChanged\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMin\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMid\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorMax\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColors\t with 3 args\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorOffset\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorsUsingFDR\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColorsFDRRate\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  RunSimulation\n",NULL);
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
    Tcl_AppendResult(interp,"  SetExcludeSubjectET\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetExcludeSubjectET\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ClearAllExcludedSubjects\n",NULL);
    Tcl_AppendResult(interp,"  ResetStatsImportFrame\n",NULL);
    Tcl_AppendResult(interp,"  GenerateStatsDataTables\n",NULL);
    Tcl_AppendResult(interp,"  SetStatsImportItem\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  AddStatsToDataTable\n",NULL);
    Tcl_AppendResult(interp,"  RemoveFactorFromDataTable\n",NULL);
    Tcl_AppendResult(interp,"  GetAnnotationForVertex\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ComposeSurfaceScalarsAndShow\n",NULL);
    return TCL_OK;
    }

  if (!strcmp("DescribeMethods",argv[1]))
    {
    if(argc>3) {
      Tcl_SetResult ( interp, const_cast<char*>("Wrong number of arguments: object DescribeMethods <MethodName>"), TCL_VOLATILE ); 
      return TCL_ERROR;
 }
    if(argc==2) {

  Tcl_DString dString, dStringParent;

  Tcl_DStringInit ( &dString );

  Tcl_DStringInit ( &dStringParent );
    vtkKWWindowCppCommand(op,interp,argc,argv);
    Tcl_DStringGetResult ( interp, &dStringParent );
    Tcl_DStringAppend ( &dString, Tcl_DStringValue ( &dStringParent ), -1 );
    Tcl_DStringAppendElement ( &dString, "New" );
    Tcl_DStringAppendElement ( &dString, "GetClassName" );
    Tcl_DStringAppendElement ( &dString, "IsA" );
    Tcl_DStringAppendElement ( &dString, "NewInstance" );
    Tcl_DStringAppendElement ( &dString, "SafeDownCast" );
    Tcl_DStringAppendElement ( &dString, "SetUseHistogramEditor" );
    Tcl_DStringAppendElement ( &dString, "FinishCreating" );
    Tcl_DStringAppendElement ( &dString, "SetCurrentMeasure" );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurfaceMeasure" );
    Tcl_DStringAppendElement ( &dString, "SetDesignMatrixType" );
    Tcl_DStringAppendElement ( &dString, "LoadDataTable" );
    Tcl_DStringAppendElement ( &dString, "LoadProjectFile" );
    Tcl_DStringAppendElement ( &dString, "LoadSurface" );
    Tcl_DStringAppendElement ( &dString, "LoadGDFFile" );
    Tcl_DStringAppendElement ( &dString, "LoadAnnotation" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceScalars" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceOverlayScalars" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceCurvatureScalars" );
    Tcl_DStringAppendElement ( &dString, "LoadLabel" );
    Tcl_DStringAppendElement ( &dString, "CreateWidget" );
    Tcl_DStringAppendElement ( &dString, "BeginActionWithProgress" );
    Tcl_DStringAppendElement ( &dString, "UpdateProgressMessage" );
    Tcl_DStringAppendElement ( &dString, "UpdateProgressPercent" );
    Tcl_DStringAppendElement ( &dString, "EndActionWithProgress" );
    Tcl_DStringAppendElement ( &dString, "LoadDataTableFromDlog" );
    Tcl_DStringAppendElement ( &dString, "LoadProjectFileFromDlog" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceFromDlog" );
    Tcl_DStringAppendElement ( &dString, "LoadGDFFromDlog" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceScalarsFromDlog" );
    Tcl_DStringAppendElement ( &dString, "LoadCurvatureFromDlog" );
    Tcl_DStringAppendElement ( &dString, "LoadAnnotationFromDlog" );
    Tcl_DStringAppendElement ( &dString, "LoadLabelFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SaveDataTableFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SaveProjectFileFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SaveScatterPlotPostscriptFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SaveTIFFImageFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SaveGDFPostscriptFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SaveLabelFromDlog" );
    Tcl_DStringAppendElement ( &dString, "MapLabelFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SmoothCurvatureScalarsFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SmoothSurfaceScalarsFromDlog" );
    Tcl_DStringAppendElement ( &dString, "SaveProjectFile" );
    Tcl_DStringAppendElement ( &dString, "SaveLabel" );
    Tcl_DStringAppendElement ( &dString, "MapLabelToSubjects" );
    Tcl_DStringAppendElement ( &dString, "SaveTIFFImage" );
    Tcl_DStringAppendElement ( &dString, "QuickSnapsTIFF" );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurfaceScalars" );
    Tcl_DStringAppendElement ( &dString, "ClearSurfaceScalars" );
    Tcl_DStringAppendElement ( &dString, "UnloadSurfaceScalars" );
    Tcl_DStringAppendElement ( &dString, "ClearCurvature" );
    Tcl_DStringAppendElement ( &dString, "RestoreView" );
    Tcl_DStringAppendElement ( &dString, "ZoomBy" );
    Tcl_DStringAppendElement ( &dString, "ZoomIn" );
    Tcl_DStringAppendElement ( &dString, "ZoomOut" );
    Tcl_DStringAppendElement ( &dString, "ShowCursor" );
    Tcl_DStringAppendElement ( &dString, "SetShowCursorFromMenu" );
    Tcl_DStringAppendElement ( &dString, "ShowCurvature" );
    Tcl_DStringAppendElement ( &dString, "SetShowCurvatureFromMenu" );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurface" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceOverlayOpacity" );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurfaceScalarsFromTableSelection" );
    Tcl_DStringAppendElement ( &dString, "DiscreteFactorsListBoxCallback" );
    Tcl_DStringAppendElement ( &dString, "ContinuousFactorsListBoxCallback" );
    Tcl_DStringAppendElement ( &dString, "NuisanceFactorsListBoxCallback" );
    Tcl_DStringAppendElement ( &dString, "NuisanceFactorsListBoxCallback" );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotListBoxCallback" );
    Tcl_DStringAppendElement ( &dString, "AnalyzeDesign" );
    Tcl_DStringAppendElement ( &dString, "SetSubjectsDir" );
    Tcl_DStringAppendElement ( &dString, "SetAverageSubject" );
    Tcl_DStringAppendElement ( &dString, "SetDesignName" );
    Tcl_DStringAppendElement ( &dString, "GetShowCurvature" );
    Tcl_DStringAppendElement ( &dString, "SetShowCurvature" );
    Tcl_DStringAppendElement ( &dString, "SetDrawCurvatureGreenRed" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMin" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMid" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMax" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorReverse" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorShowPositive" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorShowNegative" );
    Tcl_DStringAppendElement ( &dString, "GenerateClusterStats" );
    Tcl_DStringAppendElement ( &dString, "GotoNextCluster" );
    Tcl_DStringAppendElement ( &dString, "GotoPrevCluster" );
    Tcl_DStringAppendElement ( &dString, "GotoCluster" );
    Tcl_DStringAppendElement ( &dString, "CreateScalarTableEntry" );
    Tcl_DStringAppendElement ( &dString, "SurfaceScalarColorsEditorChanged" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMin" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMid" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMax" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColors" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorOffset" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorsUsingFDR" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorsFDRRate" );
    Tcl_DStringAppendElement ( &dString, "RunSimulation" );
    Tcl_DStringAppendElement ( &dString, "SelectSurfaceVertex" );
    Tcl_DStringAppendElement ( &dString, "AddSelectionToROI" );
    Tcl_DStringAppendElement ( &dString, "RemoveSelectionFromROI" );
    Tcl_DStringAppendElement ( &dString, "ClearROI" );
    Tcl_DStringAppendElement ( &dString, "GraphAverageROIInGDF" );
    Tcl_DStringAppendElement ( &dString, "SmoothCurvatureScalars" );
    Tcl_DStringAppendElement ( &dString, "SmoothSurfaceScalars" );
    Tcl_DStringAppendElement ( &dString, "NotebookPageRaised" );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotGraphMouseoverEnterElement" );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotGraphMouseoverExitElement" );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotGraphSetUpContextualMenu" );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectID" );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectGT" );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectGT" );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectLT" );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectLT" );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectET" );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectET" );
    Tcl_DStringAppendElement ( &dString, "ClearAllExcludedSubjects" );
    Tcl_DStringAppendElement ( &dString, "ResetStatsImportFrame" );
    Tcl_DStringAppendElement ( &dString, "GenerateStatsDataTables" );
    Tcl_DStringAppendElement ( &dString, "SetStatsImportItem" );
    Tcl_DStringAppendElement ( &dString, "AddStatsToDataTable" );
    Tcl_DStringAppendElement ( &dString, "RemoveFactorFromDataTable" );
    Tcl_DStringAppendElement ( &dString, "GetAnnotationForVertex" );
    Tcl_DStringAppendElement ( &dString, "ComposeSurfaceScalarsAndShow" );
  Tcl_DStringResult ( interp, &dString );
  Tcl_DStringFree ( &dString );
  Tcl_DStringFree ( &dStringParent );
    return TCL_OK;
    }
    if(argc==3) {
      Tcl_DString dString;
      int SuperClassStatus;
    SuperClassStatus = vtkKWWindowCppCommand(op,interp,argc,argv);
    if ( SuperClassStatus == TCL_OK ) { return TCL_OK; }
    /* Starting function: New */
    if ( strcmp ( argv[2], "New" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "New" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for New */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "extern static vtkKWQdecWindow *New ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for New */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GetClassName */
    if ( strcmp ( argv[2], "GetClassName" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GetClassName" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GetClassName */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "const char *GetClassName ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GetClassName */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: IsA */
    if ( strcmp ( argv[2], "IsA" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "IsA" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for IsA */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "int IsA (const char *name);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for IsA */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: NewInstance */
    if ( strcmp ( argv[2], "NewInstance" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "NewInstance" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for NewInstance */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow *NewInstance ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for NewInstance */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SafeDownCast */
    if ( strcmp ( argv[2], "SafeDownCast" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SafeDownCast" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkObject" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SafeDownCast */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow *SafeDownCast (vtkObject* o);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SafeDownCast */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetUseHistogramEditor */
    if ( strcmp ( argv[2], "SetUseHistogramEditor" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetUseHistogramEditor" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "bool" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetUseHistogramEditor */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetUseHistogramEditor (bool ibUse);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetUseHistogramEditor */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: FinishCreating */
    if ( strcmp ( argv[2], "FinishCreating" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "FinishCreating" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for FinishCreating */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void FinishCreating ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for FinishCreating */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetCurrentMeasure */
    if ( strcmp ( argv[2], "SetCurrentMeasure" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetCurrentMeasure" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetCurrentMeasure */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetCurrentMeasure (const char *isMeasure);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetCurrentMeasure */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetCurrentSurfaceMeasure */
    if ( strcmp ( argv[2], "SetCurrentSurfaceMeasure" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurfaceMeasure" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetCurrentSurfaceMeasure */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetCurrentSurfaceMeasure (const char *isMeasure);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetCurrentSurfaceMeasure */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetDesignMatrixType */
    if ( strcmp ( argv[2], "SetDesignMatrixType" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetDesignMatrixType" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetDesignMatrixType */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetDesignMatrixType (const char *isType);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetDesignMatrixType */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadDataTable */
    if ( strcmp ( argv[2], "LoadDataTable" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadDataTable" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadDataTable */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadDataTable (const char *ifnDataTable);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadDataTable */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadProjectFile */
    if ( strcmp ( argv[2], "LoadProjectFile" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadProjectFile" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadProjectFile */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadProjectFile (const char *ifnProject);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadProjectFile */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadSurface */
    if ( strcmp ( argv[2], "LoadSurface" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadSurface" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurface */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadSurface (const char *ifnSurface, const char *isLabelNULL);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadSurface */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadGDFFile */
    if ( strcmp ( argv[2], "LoadGDFFile" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadGDFFile" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadGDFFile */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadGDFFile (const char *ifnGDFFile);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadGDFFile */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadAnnotation */
    if ( strcmp ( argv[2], "LoadAnnotation" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadAnnotation" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadAnnotation */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadAnnotation (const char *ifnScalars);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadAnnotation */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadSurfaceScalars */
    if ( strcmp ( argv[2], "LoadSurfaceScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurfaceScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "int LoadSurfaceScalars (const char *ifnScalars, const char *isLabelNULL, const char *isLabel2NULL, int inFrame);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadSurfaceScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadSurfaceOverlayScalars */
    if ( strcmp ( argv[2], "LoadSurfaceOverlayScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceOverlayScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurfaceOverlayScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadSurfaceOverlayScalars (const char *ifnScalars, const char *ifnColorsNULL);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadSurfaceOverlayScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadSurfaceCurvatureScalars */
    if ( strcmp ( argv[2], "LoadSurfaceCurvatureScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceCurvatureScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurfaceCurvatureScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadSurfaceCurvatureScalars (const char *ifnScalars);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadSurfaceCurvatureScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadLabel */
    if ( strcmp ( argv[2], "LoadLabel" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadLabel" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadLabel */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadLabel (const char *ifnLabel);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadLabel */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: CreateWidget */
    if ( strcmp ( argv[2], "CreateWidget" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "CreateWidget" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for CreateWidget */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "virtual void CreateWidget ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for CreateWidget */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: BeginActionWithProgress */
    if ( strcmp ( argv[2], "BeginActionWithProgress" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "BeginActionWithProgress" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for BeginActionWithProgress */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void BeginActionWithProgress (const char *isTitle);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for BeginActionWithProgress */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: UpdateProgressMessage */
    if ( strcmp ( argv[2], "UpdateProgressMessage" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "UpdateProgressMessage" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for UpdateProgressMessage */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void UpdateProgressMessage (const char *isMessage);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for UpdateProgressMessage */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: UpdateProgressPercent */
    if ( strcmp ( argv[2], "UpdateProgressPercent" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "UpdateProgressPercent" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for UpdateProgressPercent */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void UpdateProgressPercent (float iPercent);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for UpdateProgressPercent */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: EndActionWithProgress */
    if ( strcmp ( argv[2], "EndActionWithProgress" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "EndActionWithProgress" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for EndActionWithProgress */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void EndActionWithProgress ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for EndActionWithProgress */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadDataTableFromDlog */
    if ( strcmp ( argv[2], "LoadDataTableFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadDataTableFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadDataTableFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadDataTableFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadDataTableFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadProjectFileFromDlog */
    if ( strcmp ( argv[2], "LoadProjectFileFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadProjectFileFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadProjectFileFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadProjectFileFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadProjectFileFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadSurfaceFromDlog */
    if ( strcmp ( argv[2], "LoadSurfaceFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurfaceFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadSurfaceFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadSurfaceFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadGDFFromDlog */
    if ( strcmp ( argv[2], "LoadGDFFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadGDFFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadGDFFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadGDFFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadGDFFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadSurfaceScalarsFromDlog */
    if ( strcmp ( argv[2], "LoadSurfaceScalarsFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceScalarsFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurfaceScalarsFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadSurfaceScalarsFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadSurfaceScalarsFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadCurvatureFromDlog */
    if ( strcmp ( argv[2], "LoadCurvatureFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadCurvatureFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadCurvatureFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadCurvatureFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadCurvatureFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadAnnotationFromDlog */
    if ( strcmp ( argv[2], "LoadAnnotationFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadAnnotationFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadAnnotationFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadAnnotationFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadAnnotationFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: LoadLabelFromDlog */
    if ( strcmp ( argv[2], "LoadLabelFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "LoadLabelFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadLabelFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadLabelFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for LoadLabelFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveDataTableFromDlog */
    if ( strcmp ( argv[2], "SaveDataTableFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveDataTableFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveDataTableFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveDataTableFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveDataTableFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveProjectFileFromDlog */
    if ( strcmp ( argv[2], "SaveProjectFileFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveProjectFileFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveProjectFileFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveProjectFileFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveProjectFileFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveScatterPlotPostscriptFromDlog */
    if ( strcmp ( argv[2], "SaveScatterPlotPostscriptFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveScatterPlotPostscriptFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveScatterPlotPostscriptFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveScatterPlotPostscriptFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveScatterPlotPostscriptFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveTIFFImageFromDlog */
    if ( strcmp ( argv[2], "SaveTIFFImageFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveTIFFImageFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveTIFFImageFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveTIFFImageFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveTIFFImageFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveGDFPostscriptFromDlog */
    if ( strcmp ( argv[2], "SaveGDFPostscriptFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveGDFPostscriptFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveGDFPostscriptFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveGDFPostscriptFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveGDFPostscriptFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveLabelFromDlog */
    if ( strcmp ( argv[2], "SaveLabelFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveLabelFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveLabelFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveLabelFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveLabelFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: MapLabelFromDlog */
    if ( strcmp ( argv[2], "MapLabelFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "MapLabelFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for MapLabelFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void MapLabelFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for MapLabelFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SmoothCurvatureScalarsFromDlog */
    if ( strcmp ( argv[2], "SmoothCurvatureScalarsFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SmoothCurvatureScalarsFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SmoothCurvatureScalarsFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SmoothCurvatureScalarsFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SmoothCurvatureScalarsFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SmoothSurfaceScalarsFromDlog */
    if ( strcmp ( argv[2], "SmoothSurfaceScalarsFromDlog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SmoothSurfaceScalarsFromDlog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SmoothSurfaceScalarsFromDlog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SmoothSurfaceScalarsFromDlog ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SmoothSurfaceScalarsFromDlog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveProjectFile */
    if ( strcmp ( argv[2], "SaveProjectFile" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveProjectFile" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveProjectFile */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveProjectFile (const char *ifnProject);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveProjectFile */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveLabel */
    if ( strcmp ( argv[2], "SaveLabel" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveLabel" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveLabel */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveLabel (const char *ifnLabel);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveLabel */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: MapLabelToSubjects */
    if ( strcmp ( argv[2], "MapLabelToSubjects" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "MapLabelToSubjects" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for MapLabelToSubjects */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void MapLabelToSubjects (const char *ifnLabel);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for MapLabelToSubjects */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SaveTIFFImage */
    if ( strcmp ( argv[2], "SaveTIFFImage" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SaveTIFFImage" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SaveTIFFImage */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SaveTIFFImage (const char *ifnTIFF, int iMagnificationLevel);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SaveTIFFImage */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: QuickSnapsTIFF */
    if ( strcmp ( argv[2], "QuickSnapsTIFF" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "QuickSnapsTIFF" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for QuickSnapsTIFF */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void QuickSnapsTIFF ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for QuickSnapsTIFF */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetCurrentSurfaceScalars */
    if ( strcmp ( argv[2], "SetCurrentSurfaceScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurfaceScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetCurrentSurfaceScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetCurrentSurfaceScalars (int inEntry);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetCurrentSurfaceScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ClearSurfaceScalars */
    if ( strcmp ( argv[2], "ClearSurfaceScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ClearSurfaceScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ClearSurfaceScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ClearSurfaceScalars ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ClearSurfaceScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: UnloadSurfaceScalars */
    if ( strcmp ( argv[2], "UnloadSurfaceScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "UnloadSurfaceScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for UnloadSurfaceScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void UnloadSurfaceScalars ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for UnloadSurfaceScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ClearCurvature */
    if ( strcmp ( argv[2], "ClearCurvature" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ClearCurvature" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ClearCurvature */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ClearCurvature ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ClearCurvature */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: RestoreView */
    if ( strcmp ( argv[2], "RestoreView" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "RestoreView" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for RestoreView */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void RestoreView ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for RestoreView */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ZoomBy */
    if ( strcmp ( argv[2], "ZoomBy" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ZoomBy" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ZoomBy */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ZoomBy (float iFactor);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ZoomBy */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ZoomIn */
    if ( strcmp ( argv[2], "ZoomIn" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ZoomIn" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ZoomIn */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ZoomIn ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ZoomIn */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ZoomOut */
    if ( strcmp ( argv[2], "ZoomOut" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ZoomOut" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ZoomOut */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ZoomOut ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ZoomOut */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ShowCursor */
    if ( strcmp ( argv[2], "ShowCursor" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ShowCursor" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ShowCursor */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ShowCursor (int ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ShowCursor */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetShowCursorFromMenu */
    if ( strcmp ( argv[2], "SetShowCursorFromMenu" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetShowCursorFromMenu" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetShowCursorFromMenu */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetShowCursorFromMenu ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetShowCursorFromMenu */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ShowCurvature */
    if ( strcmp ( argv[2], "ShowCurvature" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ShowCurvature" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ShowCurvature */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ShowCurvature (int ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ShowCurvature */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetShowCurvatureFromMenu */
    if ( strcmp ( argv[2], "SetShowCurvatureFromMenu" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetShowCurvatureFromMenu" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetShowCurvatureFromMenu */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetShowCurvatureFromMenu ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetShowCurvatureFromMenu */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetCurrentSurface */
    if ( strcmp ( argv[2], "SetCurrentSurface" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurface" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetCurrentSurface */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetCurrentSurface (const char *isLabel);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetCurrentSurface */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceOverlayOpacity */
    if ( strcmp ( argv[2], "SetSurfaceOverlayOpacity" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceOverlayOpacity" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceOverlayOpacity */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceOverlayOpacity (double iOpacity);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceOverlayOpacity */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetCurrentSurfaceScalarsFromTableSelection */
    if ( strcmp ( argv[2], "SetCurrentSurfaceScalarsFromTableSelection" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetCurrentSurfaceScalarsFromTableSelection" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetCurrentSurfaceScalarsFromTableSelection */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetCurrentSurfaceScalarsFromTableSelection ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetCurrentSurfaceScalarsFromTableSelection */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: DiscreteFactorsListBoxCallback */
    if ( strcmp ( argv[2], "DiscreteFactorsListBoxCallback" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "DiscreteFactorsListBoxCallback" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for DiscreteFactorsListBoxCallback */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void DiscreteFactorsListBoxCallback ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for DiscreteFactorsListBoxCallback */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ContinuousFactorsListBoxCallback */
    if ( strcmp ( argv[2], "ContinuousFactorsListBoxCallback" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ContinuousFactorsListBoxCallback" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ContinuousFactorsListBoxCallback */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ContinuousFactorsListBoxCallback ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ContinuousFactorsListBoxCallback */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: NuisanceFactorsListBoxCallback */
    if ( strcmp ( argv[2], "NuisanceFactorsListBoxCallback" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "NuisanceFactorsListBoxCallback" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for NuisanceFactorsListBoxCallback */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void NuisanceFactorsListBoxCallback ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for NuisanceFactorsListBoxCallback */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: NuisanceFactorsListBoxCallback */
    if ( strcmp ( argv[2], "NuisanceFactorsListBoxCallback" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "NuisanceFactorsListBoxCallback" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for NuisanceFactorsListBoxCallback */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void NuisanceFactorsListBoxCallback ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for NuisanceFactorsListBoxCallback */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ScatterPlotListBoxCallback */
    if ( strcmp ( argv[2], "ScatterPlotListBoxCallback" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotListBoxCallback" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ScatterPlotListBoxCallback */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ScatterPlotListBoxCallback ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ScatterPlotListBoxCallback */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AnalyzeDesign */
    if ( strcmp ( argv[2], "AnalyzeDesign" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AnalyzeDesign" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AnalyzeDesign */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AnalyzeDesign ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for AnalyzeDesign */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSubjectsDir */
    if ( strcmp ( argv[2], "SetSubjectsDir" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSubjectsDir" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSubjectsDir */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSubjectsDir (const char *isSubjectsDir);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSubjectsDir */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetAverageSubject */
    if ( strcmp ( argv[2], "SetAverageSubject" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetAverageSubject" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetAverageSubject */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetAverageSubject (const char *isAverageSubject);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetAverageSubject */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetDesignName */
    if ( strcmp ( argv[2], "SetDesignName" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetDesignName" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetDesignName */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetDesignName (const char *isDesignName);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetDesignName */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GetShowCurvature */
    if ( strcmp ( argv[2], "GetShowCurvature" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GetShowCurvature" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GetShowCurvature */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "bool GetShowCurvature ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GetShowCurvature */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetShowCurvature */
    if ( strcmp ( argv[2], "SetShowCurvature" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetShowCurvature" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "bool" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetShowCurvature */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetShowCurvature (bool ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetShowCurvature */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetDrawCurvatureGreenRed */
    if ( strcmp ( argv[2], "SetDrawCurvatureGreenRed" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetDrawCurvatureGreenRed" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetDrawCurvatureGreenRed */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetDrawCurvatureGreenRed (int ibDraw);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetDrawCurvatureGreenRed */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorMin */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorMin" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMin" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorMin */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorMin (const char *isMin);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorMin */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorMid */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorMid" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMid" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorMid */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorMid (const char *isMid);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorMid */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorMax */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorMax" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMax" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorMax */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorMax (const char *isMax);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorMax */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorReverse */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorReverse" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorReverse" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorReverse */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorReverse (int ibReverse);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorReverse */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorShowPositive */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorShowPositive" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorShowPositive" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorShowPositive */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorShowPositive (int ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorShowPositive */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorShowNegative */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorShowNegative" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorShowNegative" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorShowNegative */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorShowNegative (int ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorShowNegative */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GenerateClusterStats */
    if ( strcmp ( argv[2], "GenerateClusterStats" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GenerateClusterStats" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GenerateClusterStats */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void GenerateClusterStats ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GenerateClusterStats */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GotoNextCluster */
    if ( strcmp ( argv[2], "GotoNextCluster" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GotoNextCluster" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GotoNextCluster */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void GotoNextCluster ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GotoNextCluster */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GotoPrevCluster */
    if ( strcmp ( argv[2], "GotoPrevCluster" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GotoPrevCluster" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GotoPrevCluster */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void GotoPrevCluster ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GotoPrevCluster */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GotoCluster */
    if ( strcmp ( argv[2], "GotoCluster" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GotoCluster" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GotoCluster */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void GotoCluster (int iCurrentCluster);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GotoCluster */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: CreateScalarTableEntry */
    if ( strcmp ( argv[2], "CreateScalarTableEntry" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "CreateScalarTableEntry" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for CreateScalarTableEntry */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void CreateScalarTableEntry (const char *iTable, int iRow, int iColumn, const char *iWidget);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for CreateScalarTableEntry */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SurfaceScalarColorsEditorChanged */
    if ( strcmp ( argv[2], "SurfaceScalarColorsEditorChanged" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SurfaceScalarColorsEditorChanged" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SurfaceScalarColorsEditorChanged */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SurfaceScalarColorsEditorChanged ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SurfaceScalarColorsEditorChanged */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorMin */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorMin" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMin" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorMin */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorMin (double iMin);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorMin */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorMid */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorMid" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMid" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorMid */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorMid (double iMid);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorMid */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorMax */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorMax" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorMax" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorMax */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorMax (double iMax);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorMax */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColors */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColors" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColors" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColors */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColors (double iMin, double iMid, double iMax);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColors */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorOffset */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorOffset" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorOffset" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorOffset */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorOffset (double iOffset);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorOffset */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorsUsingFDR */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorsUsingFDR" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorsUsingFDR" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorsUsingFDR */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorsUsingFDR ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorsUsingFDR */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalarsColorsFDRRate */
    if ( strcmp ( argv[2], "SetSurfaceScalarsColorsFDRRate" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColorsFDRRate" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColorsFDRRate */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColorsFDRRate (const char *isValue);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetSurfaceScalarsColorsFDRRate */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: RunSimulation */
    if ( strcmp ( argv[2], "RunSimulation" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "RunSimulation" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for RunSimulation */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void RunSimulation ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for RunSimulation */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SelectSurfaceVertex */
    if ( strcmp ( argv[2], "SelectSurfaceVertex" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SelectSurfaceVertex" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SelectSurfaceVertex */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SelectSurfaceVertex (int inVertex);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SelectSurfaceVertex */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AddSelectionToROI */
    if ( strcmp ( argv[2], "AddSelectionToROI" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AddSelectionToROI" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AddSelectionToROI */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AddSelectionToROI ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for AddSelectionToROI */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: RemoveSelectionFromROI */
    if ( strcmp ( argv[2], "RemoveSelectionFromROI" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "RemoveSelectionFromROI" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for RemoveSelectionFromROI */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void RemoveSelectionFromROI ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for RemoveSelectionFromROI */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ClearROI */
    if ( strcmp ( argv[2], "ClearROI" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ClearROI" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ClearROI */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ClearROI ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ClearROI */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GraphAverageROIInGDF */
    if ( strcmp ( argv[2], "GraphAverageROIInGDF" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GraphAverageROIInGDF" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GraphAverageROIInGDF */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void GraphAverageROIInGDF ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GraphAverageROIInGDF */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SmoothCurvatureScalars */
    if ( strcmp ( argv[2], "SmoothCurvatureScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SmoothCurvatureScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SmoothCurvatureScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SmoothCurvatureScalars (int icSteps);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SmoothCurvatureScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SmoothSurfaceScalars */
    if ( strcmp ( argv[2], "SmoothSurfaceScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SmoothSurfaceScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SmoothSurfaceScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SmoothSurfaceScalars (int icSteps);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SmoothSurfaceScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: NotebookPageRaised */
    if ( strcmp ( argv[2], "NotebookPageRaised" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "NotebookPageRaised" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for NotebookPageRaised */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void NotebookPageRaised (const char *isTitle);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for NotebookPageRaised */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ScatterPlotGraphMouseoverEnterElement */
    if ( strcmp ( argv[2], "ScatterPlotGraphMouseoverEnterElement" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotGraphMouseoverEnterElement" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ScatterPlotGraphMouseoverEnterElement */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ScatterPlotGraphMouseoverEnterElement (const char *isElement);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ScatterPlotGraphMouseoverEnterElement */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ScatterPlotGraphMouseoverExitElement */
    if ( strcmp ( argv[2], "ScatterPlotGraphMouseoverExitElement" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotGraphMouseoverExitElement" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ScatterPlotGraphMouseoverExitElement */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ScatterPlotGraphMouseoverExitElement ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ScatterPlotGraphMouseoverExitElement */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ScatterPlotGraphSetUpContextualMenu */
    if ( strcmp ( argv[2], "ScatterPlotGraphSetUpContextualMenu" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ScatterPlotGraphSetUpContextualMenu" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "vtkKWMenu" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ScatterPlotGraphSetUpContextualMenu */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ScatterPlotGraphSetUpContextualMenu (const char *isElement, vtkKWMenu *iMenu);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ScatterPlotGraphSetUpContextualMenu */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetExcludeSubjectID */
    if ( strcmp ( argv[2], "SetExcludeSubjectID" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectID" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetExcludeSubjectID */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetExcludeSubjectID (const char *isElement, int ibExclude);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetExcludeSubjectID */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetExcludeSubjectGT */
    if ( strcmp ( argv[2], "SetExcludeSubjectGT" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectGT" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetExcludeSubjectGT */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetExcludeSubjectGT (double inExcludeGT);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetExcludeSubjectGT */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetExcludeSubjectGT */
    if ( strcmp ( argv[2], "SetExcludeSubjectGT" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectGT" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetExcludeSubjectGT */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetExcludeSubjectGT (const char *isExcludeGT);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetExcludeSubjectGT */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetExcludeSubjectLT */
    if ( strcmp ( argv[2], "SetExcludeSubjectLT" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectLT" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetExcludeSubjectLT */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetExcludeSubjectLT (double inExcludeLT);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetExcludeSubjectLT */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetExcludeSubjectLT */
    if ( strcmp ( argv[2], "SetExcludeSubjectLT" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectLT" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetExcludeSubjectLT */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetExcludeSubjectLT (const char *isExcludeLT);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetExcludeSubjectLT */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetExcludeSubjectET */
    if ( strcmp ( argv[2], "SetExcludeSubjectET" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectET" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "float" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetExcludeSubjectET */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetExcludeSubjectET (double inExcludeET);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetExcludeSubjectET */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetExcludeSubjectET */
    if ( strcmp ( argv[2], "SetExcludeSubjectET" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetExcludeSubjectET" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetExcludeSubjectET */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetExcludeSubjectET (const char *isExcludeET);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetExcludeSubjectET */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ClearAllExcludedSubjects */
    if ( strcmp ( argv[2], "ClearAllExcludedSubjects" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ClearAllExcludedSubjects" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ClearAllExcludedSubjects */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ClearAllExcludedSubjects ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ClearAllExcludedSubjects */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ResetStatsImportFrame */
    if ( strcmp ( argv[2], "ResetStatsImportFrame" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ResetStatsImportFrame" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ResetStatsImportFrame */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ResetStatsImportFrame ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ResetStatsImportFrame */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GenerateStatsDataTables */
    if ( strcmp ( argv[2], "GenerateStatsDataTables" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GenerateStatsDataTables" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GenerateStatsDataTables */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void GenerateStatsDataTables ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GenerateStatsDataTables */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetStatsImportItem */
    if ( strcmp ( argv[2], "SetStatsImportItem" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetStatsImportItem" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetStatsImportItem */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetStatsImportItem (const char *isStatsImportItem);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for SetStatsImportItem */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AddStatsToDataTable */
    if ( strcmp ( argv[2], "AddStatsToDataTable" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AddStatsToDataTable" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AddStatsToDataTable */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AddStatsToDataTable ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for AddStatsToDataTable */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: RemoveFactorFromDataTable */
    if ( strcmp ( argv[2], "RemoveFactorFromDataTable" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "RemoveFactorFromDataTable" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for RemoveFactorFromDataTable */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void RemoveFactorFromDataTable ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for RemoveFactorFromDataTable */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GetAnnotationForVertex */
    if ( strcmp ( argv[2], "GetAnnotationForVertex" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GetAnnotationForVertex" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GetAnnotationForVertex */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "const char *GetAnnotationForVertex (int inVertex);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for GetAnnotationForVertex */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ComposeSurfaceScalarsAndShow */
    if ( strcmp ( argv[2], "ComposeSurfaceScalarsAndShow" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ComposeSurfaceScalarsAndShow" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ComposeSurfaceScalarsAndShow */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ComposeSurfaceScalarsAndShow ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecWindow" );
    /* Closing for ComposeSurfaceScalarsAndShow */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
   Tcl_SetResult ( interp, const_cast<char*>("Could not find method"), TCL_VOLATILE ); 
   return TCL_ERROR;
   }
 }

  if (vtkKWWindowCppCommand(static_cast<vtkKWWindow *>(op),interp,argc,argv) == TCL_OK)
    {
    return TCL_OK;
    }
    }
  catch (vtkstd::exception &e)
    {
    Tcl_AppendResult(interp, "Uncaught exception: ",  e.what(), "\n", NULL);
    return TCL_ERROR;
    }
  return TCL_ERROR;
}
