// tcl wrapper for vtkKWQdecView object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkKWQdecView.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>
#include <vtksys/ios/sstream>

ClientData vtkKWQdecViewNewCommand()
{
  vtkKWQdecView *temp = vtkKWQdecView::New();
  return static_cast<ClientData>(temp);
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
   return vtkKWQdecViewCppCommand(static_cast<vtkKWQdecView *>(static_cast<vtkTclCommandArgStruct *>(cd)->Pointer),interp, argc, argv);
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
    Tcl_SetResult(interp,const_cast<char *>("Could not find requested method."), TCL_VOLATILE);
    return TCL_ERROR;
    }
  if (!interp)
    {
    if (!strcmp("DoTypecasting",argv[0]))
      {
      if (!strcmp("vtkKWQdecView",argv[1]))
        {
        argv[2] = static_cast<char *>(static_cast<void *>(op));
        return TCL_OK;
        }
      if (vtkKWRenderWidgetCppCommand(static_cast<vtkKWRenderWidget *>(op),interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,const_cast<char *>("vtkKWRenderWidget"), TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkKWQdecView  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecView");
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
    vtkKWQdecView  *temp20;
    temp20 = (op)->NewInstance();
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecView");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkKWQdecView  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkObject"),interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecView");
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
  if ((!strcmp("RestoreView",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->RestoreView(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
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

    temp0 = (vtkFSSurfaceSource *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkFSSurfaceSource"),interp,error));
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

    temp0 = (vtkFloatArray *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkFloatArray"),interp,error));
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

    temp0 = (vtkScalarsToColors *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkScalarsToColors"),interp,error));
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

    temp0 = (vtkFloatArray *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkFloatArray"),interp,error));
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

    temp0 = (vtkFloatArray *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkFloatArray"),interp,error));
    temp1 = (vtkScalarsToColors *)(vtkTclGetPointerFromObject(argv[3],const_cast<char *>("vtkScalarsToColors"),interp,error));
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
    Tcl_PrintDouble(interp,temp20,tempResult);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
    return TCL_OK;
    }
  if ((!strcmp("SetSurfaceLegendColors",argv[1]))&&(argc == 3))
    {
    vtkScalarsToColors  *temp0;
    error = 0;

    temp0 = (vtkScalarsToColors *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkScalarsToColors"),interp,error));
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

    temp0 = (vtkPolyData *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkPolyData"),interp,error));
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
    vtkTclListInstances(interp,(ClientData)(vtkKWQdecViewCommand));
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
    Tcl_AppendResult(interp,"  RestoreView\t with 1 arg\n",NULL);
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
    vtkKWRenderWidgetCppCommand(op,interp,argc,argv);
    Tcl_DStringGetResult ( interp, &dStringParent );
    Tcl_DStringAppend ( &dString, Tcl_DStringValue ( &dStringParent ), -1 );
    Tcl_DStringAppendElement ( &dString, "New" );
    Tcl_DStringAppendElement ( &dString, "GetClassName" );
    Tcl_DStringAppendElement ( &dString, "IsA" );
    Tcl_DStringAppendElement ( &dString, "NewInstance" );
    Tcl_DStringAppendElement ( &dString, "SafeDownCast" );
    Tcl_DStringAppendElement ( &dString, "CreateWidget" );
    Tcl_DStringAppendElement ( &dString, "ResetView" );
    Tcl_DStringAppendElement ( &dString, "RestoreView" );
    Tcl_DStringAppendElement ( &dString, "ZoomBy" );
    Tcl_DStringAppendElement ( &dString, "GetShowCursor" );
    Tcl_DStringAppendElement ( &dString, "SetShowCursor" );
    Tcl_DStringAppendElement ( &dString, "SetSurface" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalars" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalarsColors" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceLookupScalars" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceOverlayScalarsAndColors" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceOverlayOpacity" );
    Tcl_DStringAppendElement ( &dString, "GetSurfaceOverlayOpacity" );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceLegendColors" );
    Tcl_DStringAppendElement ( &dString, "SetROI" );
    Tcl_DStringAppendElement ( &dString, "SetAnnotationMessage" );
    Tcl_DStringAppendElement ( &dString, "SetShowLegend" );
    Tcl_DStringAppendElement ( &dString, "SetShowAnnotation" );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraElevateNegative" );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraElevatePositive" );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraAzimuthNegative" );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraAzimuthPositive" );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraRollNegative" );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraRollPositive" );
    Tcl_DStringAppendElement ( &dString, "SelectSurfaceVertex" );
  Tcl_DStringResult ( interp, &dString );
  Tcl_DStringFree ( &dString );
  Tcl_DStringFree ( &dStringParent );
    return TCL_OK;
    }
    if(argc==3) {
      Tcl_DString dString;
      int SuperClassStatus;
    SuperClassStatus = vtkKWRenderWidgetCppCommand(op,interp,argc,argv);
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
    Tcl_DStringAppendElement ( &dString, "static vtkKWQdecView *New ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView *NewInstance ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView *SafeDownCast (vtkObject* o);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SafeDownCast */

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
    Tcl_DStringAppendElement ( &dString, "void CreateWidget ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for CreateWidget */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ResetView */
    if ( strcmp ( argv[2], "ResetView" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ResetView" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ResetView */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void ResetView ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for ResetView */

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
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for RestoreView */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void RestoreView (const char *isHemi);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for ZoomBy */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GetShowCursor */
    if ( strcmp ( argv[2], "GetShowCursor" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GetShowCursor" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GetShowCursor */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "bool GetShowCursor ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for GetShowCursor */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetShowCursor */
    if ( strcmp ( argv[2], "SetShowCursor" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetShowCursor" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "bool" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetShowCursor */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetShowCursor (bool ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetShowCursor */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurface */
    if ( strcmp ( argv[2], "SetSurface" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurface" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkFSSurfaceSource" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurface */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurface (vtkFSSurfaceSource *iSurfaceSource);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetSurface */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceScalars */
    if ( strcmp ( argv[2], "SetSurfaceScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkFloatArray" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalars (vtkFloatArray *iScalars);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetSurfaceScalars */

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
    Tcl_DStringAppendElement ( &dString, "vtkScalarsToColors" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceScalarsColors */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceScalarsColors (vtkScalarsToColors *iScalarColors);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetSurfaceScalarsColors */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceLookupScalars */
    if ( strcmp ( argv[2], "SetSurfaceLookupScalars" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceLookupScalars" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkFloatArray" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceLookupScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceLookupScalars (vtkFloatArray *iScalars);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetSurfaceLookupScalars */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceOverlayScalarsAndColors */
    if ( strcmp ( argv[2], "SetSurfaceOverlayScalarsAndColors" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceOverlayScalarsAndColors" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkFloatArray" );
    Tcl_DStringAppendElement ( &dString, "vtkScalarsToColors" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceOverlayScalarsAndColors */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceOverlayScalarsAndColors (vtkFloatArray *iScalars, vtkScalarsToColors *iColors);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetSurfaceOverlayScalarsAndColors */

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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetSurfaceOverlayOpacity */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: GetSurfaceOverlayOpacity */
    if ( strcmp ( argv[2], "GetSurfaceOverlayOpacity" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "GetSurfaceOverlayOpacity" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for GetSurfaceOverlayOpacity */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "double GetSurfaceOverlayOpacity ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for GetSurfaceOverlayOpacity */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetSurfaceLegendColors */
    if ( strcmp ( argv[2], "SetSurfaceLegendColors" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetSurfaceLegendColors" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkScalarsToColors" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetSurfaceLegendColors */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetSurfaceLegendColors (vtkScalarsToColors *iLegendColors);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetSurfaceLegendColors */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetROI */
    if ( strcmp ( argv[2], "SetROI" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetROI" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkPolyData" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetROI */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetROI (vtkPolyData *iROIPolyData);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetROI */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetAnnotationMessage */
    if ( strcmp ( argv[2], "SetAnnotationMessage" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetAnnotationMessage" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetAnnotationMessage */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetAnnotationMessage (const char *isMessage);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetAnnotationMessage */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetShowLegend */
    if ( strcmp ( argv[2], "SetShowLegend" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetShowLegend" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetShowLegend */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetShowLegend (int ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetShowLegend */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetShowAnnotation */
    if ( strcmp ( argv[2], "SetShowAnnotation" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetShowAnnotation" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "int" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetShowAnnotation */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void SetShowAnnotation (int ibShow);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SetShowAnnotation */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AnimateCameraElevateNegative */
    if ( strcmp ( argv[2], "AnimateCameraElevateNegative" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraElevateNegative" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AnimateCameraElevateNegative */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AnimateCameraElevateNegative ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for AnimateCameraElevateNegative */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AnimateCameraElevatePositive */
    if ( strcmp ( argv[2], "AnimateCameraElevatePositive" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraElevatePositive" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AnimateCameraElevatePositive */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AnimateCameraElevatePositive ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for AnimateCameraElevatePositive */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AnimateCameraAzimuthNegative */
    if ( strcmp ( argv[2], "AnimateCameraAzimuthNegative" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraAzimuthNegative" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AnimateCameraAzimuthNegative */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AnimateCameraAzimuthNegative ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for AnimateCameraAzimuthNegative */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AnimateCameraAzimuthPositive */
    if ( strcmp ( argv[2], "AnimateCameraAzimuthPositive" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraAzimuthPositive" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AnimateCameraAzimuthPositive */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AnimateCameraAzimuthPositive ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for AnimateCameraAzimuthPositive */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AnimateCameraRollNegative */
    if ( strcmp ( argv[2], "AnimateCameraRollNegative" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraRollNegative" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AnimateCameraRollNegative */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AnimateCameraRollNegative ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for AnimateCameraRollNegative */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: AnimateCameraRollPositive */
    if ( strcmp ( argv[2], "AnimateCameraRollPositive" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "AnimateCameraRollPositive" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for AnimateCameraRollPositive */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void AnimateCameraRollPositive ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for AnimateCameraRollPositive */

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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecView" );
    /* Closing for SelectSurfaceVertex */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
   Tcl_SetResult ( interp, const_cast<char*>("Could not find method"), TCL_VOLATILE ); 
   return TCL_ERROR;
   }
 }

  if (vtkKWRenderWidgetCppCommand(static_cast<vtkKWRenderWidget *>(op),interp,argc,argv) == TCL_OK)
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
