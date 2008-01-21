// tcl wrapper for vtkKWQdecView object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkKWQdecView.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>

ClientData vtkKWQdecViewNewCommand()
{
  vtkKWQdecView *temp = vtkKWQdecView::New();
  return ((ClientData)temp);
}

int vtkKWRenderWidgetCppCommand(vtkKWRenderWidget *op, Tcl_Interp *interp,
             int argc, char *argv[]);
int VTKTCL_EXPORT vtkKWQdecViewCppCommand(vtkKWQdecView *op, Tcl_Interp *interp,
             int argc, char *argv[]);

int VTKTCL_EXPORT vtkKWQdecViewCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[])
{
  if ((argc == 2)&&(!strcmp("Delete",argv[1]))&& !vtkTclInDelete(interp))
    {
    Tcl_DeleteCommand(interp,argv[0]);
    return TCL_OK;
    }
   return vtkKWQdecViewCppCommand((vtkKWQdecView *)(((vtkTclCommandArgStruct *)cd)->Pointer),interp, argc, argv);
}

int VTKTCL_EXPORT vtkKWQdecViewCppCommand(vtkKWQdecView *op, Tcl_Interp *interp,
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
      if (!strcmp("vtkKWQdecView",argv[1]))
        {
        argv[2] = (char *)((void *)op);
        return TCL_OK;
        }
      if (vtkKWRenderWidgetCppCommand((vtkKWRenderWidget *)op,interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,(char *) "vtkKWRenderWidget", TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkKWQdecView  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecView");
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
    vtkKWQdecView  *temp20;
    temp20 = (op)->NewInstance();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecView");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkKWQdecView  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkObject",interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecView");
    return TCL_OK;
    }
    }
  if ((!strcmp("CreateWidget",argv[1]))&&(argc == 2))
    {
    op->CreateWidget();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("ResetView",argv[1]))&&(argc == 2))
    {
    op->ResetView();
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
  if ((!strcmp("GetShowCursor",argv[1]))&&(argc == 2))
    {
    bool   temp20;
    temp20 = (op)->GetShowCursor();
    char tempResult[1024];
    sprintf(tempResult,"%i",(int)temp20);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
    return TCL_OK;
    }
  if ((!strcmp("SetShowCursor",argv[1]))&&(argc == 3))
    {
    bool   temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi ? true : false;
    if (!error)
    {
    op->SetShowCursor(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurface",argv[1]))&&(argc == 3))
    {
    vtkFSSurfaceSource  *temp0;
    error = 0;

    temp0 = (vtkFSSurfaceSource *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkFSSurfaceSource",interp,error));
    if (!error)
    {
    op->SetSurface(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalars",argv[1]))&&(argc == 3))
    {
    vtkFloatArray  *temp0;
    error = 0;

    temp0 = (vtkFloatArray *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkFloatArray",interp,error));
    if (!error)
    {
    op->SetSurfaceScalars(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceScalarsColors",argv[1]))&&(argc == 3))
    {
    vtkScalarsToColors  *temp0;
    error = 0;

    temp0 = (vtkScalarsToColors *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkScalarsToColors",interp,error));
    if (!error)
    {
    op->SetSurfaceScalarsColors(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceLookupScalars",argv[1]))&&(argc == 3))
    {
    vtkFloatArray  *temp0;
    error = 0;

    temp0 = (vtkFloatArray *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkFloatArray",interp,error));
    if (!error)
    {
    op->SetSurfaceLookupScalars(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetSurfaceOverlayScalarsAndColors",argv[1]))&&(argc == 4))
    {
    vtkFloatArray  *temp0;
    vtkScalarsToColors  *temp1;
    error = 0;

    temp0 = (vtkFloatArray *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkFloatArray",interp,error));
    temp1 = (vtkScalarsToColors *)(vtkTclGetPointerFromObject(argv[3],(char *) "vtkScalarsToColors",interp,error));
    if (!error)
    {
    op->SetSurfaceOverlayScalarsAndColors(temp0,temp1);
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
  if ((!strcmp("GetSurfaceOverlayOpacity",argv[1]))&&(argc == 2))
    {
    double   temp20;
    temp20 = (op)->GetSurfaceOverlayOpacity();
    char tempResult[1024];
    sprintf(tempResult,"%g",temp20);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
    return TCL_OK;
    }
  if ((!strcmp("SetSurfaceLegendColors",argv[1]))&&(argc == 3))
    {
    vtkScalarsToColors  *temp0;
    error = 0;

    temp0 = (vtkScalarsToColors *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkScalarsToColors",interp,error));
    if (!error)
    {
    op->SetSurfaceLegendColors(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetROI",argv[1]))&&(argc == 3))
    {
    vtkPolyData  *temp0;
    error = 0;

    temp0 = (vtkPolyData *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkPolyData",interp,error));
    if (!error)
    {
    op->SetROI(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetAnnotationMessage",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->SetAnnotationMessage(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetShowLegend",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetShowLegend(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("SetShowAnnotation",argv[1]))&&(argc == 3))
    {
    int      temp0;
    error = 0;

    if (Tcl_GetInt(interp,argv[2],&tempi) != TCL_OK) error = 1;
    temp0 = tempi;
    if (!error)
    {
    op->SetShowAnnotation(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("AnimateCameraElevateNegative",argv[1]))&&(argc == 2))
    {
    op->AnimateCameraElevateNegative();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("AnimateCameraElevatePositive",argv[1]))&&(argc == 2))
    {
    op->AnimateCameraElevatePositive();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("AnimateCameraAzimuthNegative",argv[1]))&&(argc == 2))
    {
    op->AnimateCameraAzimuthNegative();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("AnimateCameraAzimuthPositive",argv[1]))&&(argc == 2))
    {
    op->AnimateCameraAzimuthPositive();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("AnimateCameraRollNegative",argv[1]))&&(argc == 2))
    {
    op->AnimateCameraRollNegative();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("AnimateCameraRollPositive",argv[1]))&&(argc == 2))
    {
    op->AnimateCameraRollPositive();
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

  if (!strcmp("ListInstances",argv[1]))
    {
    vtkTclListInstances(interp,(ClientData)vtkKWQdecViewCommand);
    return TCL_OK;
    }

  if (!strcmp("ListMethods",argv[1]))
    {
    vtkKWRenderWidgetCppCommand(op,interp,argc,argv);
    Tcl_AppendResult(interp,"Methods from vtkKWQdecView:\n",NULL);
    Tcl_AppendResult(interp,"  GetSuperClassName\n",NULL);
    Tcl_AppendResult(interp,"  New\n",NULL);
    Tcl_AppendResult(interp,"  GetClassName\n",NULL);
    Tcl_AppendResult(interp,"  IsA\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  NewInstance\n",NULL);
    Tcl_AppendResult(interp,"  SafeDownCast\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  CreateWidget\n",NULL);
    Tcl_AppendResult(interp,"  ResetView\n",NULL);
    Tcl_AppendResult(interp,"  RestoreView\n",NULL);
    Tcl_AppendResult(interp,"  ZoomBy\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  GetShowCursor\n",NULL);
    Tcl_AppendResult(interp,"  SetShowCursor\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurface\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceScalarsColors\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceLookupScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceOverlayScalarsAndColors\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceOverlayOpacity\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  GetSurfaceOverlayOpacity\n",NULL);
    Tcl_AppendResult(interp,"  SetSurfaceLegendColors\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetROI\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetAnnotationMessage\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetShowLegend\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetShowAnnotation\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  AnimateCameraElevateNegative\n",NULL);
    Tcl_AppendResult(interp,"  AnimateCameraElevatePositive\n",NULL);
    Tcl_AppendResult(interp,"  AnimateCameraAzimuthNegative\n",NULL);
    Tcl_AppendResult(interp,"  AnimateCameraAzimuthPositive\n",NULL);
    Tcl_AppendResult(interp,"  AnimateCameraRollNegative\n",NULL);
    Tcl_AppendResult(interp,"  AnimateCameraRollPositive\n",NULL);
    Tcl_AppendResult(interp,"  SelectSurfaceVertex\t with 1 arg\n",NULL);
    return TCL_OK;
    }

  if (vtkKWRenderWidgetCppCommand((vtkKWRenderWidget *)op,interp,argc,argv) == TCL_OK)
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
