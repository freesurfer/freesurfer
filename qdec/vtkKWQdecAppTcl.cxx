// tcl wrapper for vtkKWQdecApp object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkKWQdecApp.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>
#include <vtksys/ios/sstream>

ClientData vtkKWQdecAppNewCommand()
{
  vtkKWQdecApp *temp = vtkKWQdecApp::New();
  return static_cast<ClientData>(temp);
}

int vtkKWApplicationCppCommand(vtkKWApplication *op, Tcl_Interp *interp,
             int argc, char *argv[]);
int VTKTCL_EXPORT vtkKWQdecAppCppCommand(vtkKWQdecApp *op, Tcl_Interp *interp,
             int argc, char *argv[]);

int VTKTCL_EXPORT vtkKWQdecAppCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[])
{
  if ((argc == 2)&&(!strcmp("Delete",argv[1]))&& !vtkTclInDelete(interp))
    {
    Tcl_DeleteCommand(interp,argv[0]);
    return TCL_OK;
    }
   return vtkKWQdecAppCppCommand(static_cast<vtkKWQdecApp *>(static_cast<vtkTclCommandArgStruct *>(cd)->Pointer),interp, argc, argv);
}

int VTKTCL_EXPORT vtkKWQdecAppCppCommand(vtkKWQdecApp *op, Tcl_Interp *interp,
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
      if (!strcmp("vtkKWQdecApp",argv[1]))
        {
        argv[2] = static_cast<char *>(static_cast<void *>(op));
        return TCL_OK;
        }
      if (vtkKWApplicationCppCommand(static_cast<vtkKWApplication *>(op),interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,const_cast<char *>("vtkKWApplication"), TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkKWQdecApp  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecApp");
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
    vtkKWQdecApp  *temp20;
    temp20 = (op)->NewInstance();
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecApp");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkKWQdecApp  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkObject"),interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkKWQdecApp");
    return TCL_OK;
    }
    }
  if ((!strcmp("Exit",argv[1]))&&(argc == 2))
    {
    int      temp20;
    temp20 = (op)->Exit();
    char tempResult[1024];
    sprintf(tempResult,"%i",temp20);
    Tcl_SetResult(interp, tempResult, TCL_VOLATILE);
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
  if ((!strcmp("LoadSurface",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadSurface(temp0);
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
  if ((!strcmp("LoadSurfaceScalars",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->LoadSurfaceScalars(temp0);
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
  if ((!strcmp("ErrorMessage",argv[1]))&&(argc == 3))
    {
    char    *temp0;
    error = 0;

    temp0 = argv[2];
    if (!error)
    {
    op->ErrorMessage(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("DisplayHelpDialog",argv[1]))&&(argc == 3))
    {
    vtkKWTopLevel  *temp0;
    error = 0;

    temp0 = (vtkKWTopLevel *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkKWTopLevel"),interp,error));
    if (!error)
    {
    op->DisplayHelpDialog(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }

  if (!strcmp("ListInstances",argv[1]))
    {
    vtkTclListInstances(interp,(ClientData)(vtkKWQdecAppCommand));
    return TCL_OK;
    }

  if (!strcmp("ListMethods",argv[1]))
    {
    vtkKWApplicationCppCommand(op,interp,argc,argv);
    Tcl_AppendResult(interp,"Methods from vtkKWQdecApp:\n",NULL);
    Tcl_AppendResult(interp,"  GetSuperClassName\n",NULL);
    Tcl_AppendResult(interp,"  New\n",NULL);
    Tcl_AppendResult(interp,"  GetClassName\n",NULL);
    Tcl_AppendResult(interp,"  IsA\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  NewInstance\n",NULL);
    Tcl_AppendResult(interp,"  SafeDownCast\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  Exit\n",NULL);
    Tcl_AppendResult(interp,"  LoadDataTable\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadProjectFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurface\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadGDFFile\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceCurvatureScalars\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadAnnotation\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  LoadSurfaceOverlayScalars\t with 2 args\n",NULL);
    Tcl_AppendResult(interp,"  LoadLabel\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetAverageSubject\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  ErrorMessage\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  DisplayHelpDialog\t with 1 arg\n",NULL);
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
    vtkKWApplicationCppCommand(op,interp,argc,argv);
    Tcl_DStringGetResult ( interp, &dStringParent );
    Tcl_DStringAppend ( &dString, Tcl_DStringValue ( &dStringParent ), -1 );
    Tcl_DStringAppendElement ( &dString, "New" );
    Tcl_DStringAppendElement ( &dString, "GetClassName" );
    Tcl_DStringAppendElement ( &dString, "IsA" );
    Tcl_DStringAppendElement ( &dString, "NewInstance" );
    Tcl_DStringAppendElement ( &dString, "SafeDownCast" );
    Tcl_DStringAppendElement ( &dString, "Exit" );
    Tcl_DStringAppendElement ( &dString, "LoadDataTable" );
    Tcl_DStringAppendElement ( &dString, "LoadProjectFile" );
    Tcl_DStringAppendElement ( &dString, "LoadSurface" );
    Tcl_DStringAppendElement ( &dString, "LoadGDFFile" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceScalars" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceCurvatureScalars" );
    Tcl_DStringAppendElement ( &dString, "LoadAnnotation" );
    Tcl_DStringAppendElement ( &dString, "LoadSurfaceOverlayScalars" );
    Tcl_DStringAppendElement ( &dString, "LoadLabel" );
    Tcl_DStringAppendElement ( &dString, "SetAverageSubject" );
    Tcl_DStringAppendElement ( &dString, "ErrorMessage" );
    Tcl_DStringAppendElement ( &dString, "DisplayHelpDialog" );
  Tcl_DStringResult ( interp, &dString );
  Tcl_DStringFree ( &dString );
  Tcl_DStringFree ( &dStringParent );
    return TCL_OK;
    }
    if(argc==3) {
      Tcl_DString dString;
      int SuperClassStatus;
    SuperClassStatus = vtkKWApplicationCppCommand(op,interp,argc,argv);
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
    Tcl_DStringAppendElement ( &dString, "static vtkKWQdecApp *New ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp *NewInstance ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp *SafeDownCast (vtkObject* o);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for SafeDownCast */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: Exit */
    if ( strcmp ( argv[2], "Exit" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "Exit" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for Exit */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "virtual int Exit ();" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for Exit */

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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
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
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurface */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadSurface (const char *ifnSurface);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
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
    Tcl_DStringAppendElement ( &dString, "void LoadGDFFile (const char *ifnGDF);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for LoadGDFFile */

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
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for LoadSurfaceScalars */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void LoadSurfaceScalars (const char *ifnScalars);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for LoadSurfaceScalars */

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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for LoadSurfaceCurvatureScalars */

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
    Tcl_DStringAppendElement ( &dString, "void LoadAnnotation (const char *ifnAnnotation);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for LoadAnnotation */

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
    Tcl_DStringAppendElement ( &dString, "void LoadSurfaceOverlayScalars (const char *ifnScalars, const char *ifnColors);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for LoadSurfaceOverlayScalars */

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
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for LoadLabel */

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
    Tcl_DStringAppendElement ( &dString, "void SetAverageSubject (const char *isAvgSubj);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for SetAverageSubject */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: ErrorMessage */
    if ( strcmp ( argv[2], "ErrorMessage" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "ErrorMessage" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "string" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for ErrorMessage */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "virtual void ErrorMessage (const char *isMessage);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for ErrorMessage */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: DisplayHelpDialog */
    if ( strcmp ( argv[2], "DisplayHelpDialog" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "DisplayHelpDialog" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkKWTopLevel" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for DisplayHelpDialog */
    Tcl_DStringAppendElement ( &dString, "" );
    Tcl_DStringAppendElement ( &dString, "void DisplayHelpDialog (vtkKWTopLevel *iTop);" );
    Tcl_DStringAppendElement ( &dString, "vtkKWQdecApp" );
    /* Closing for DisplayHelpDialog */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
   Tcl_SetResult ( interp, const_cast<char*>("Could not find method"), TCL_VOLATILE ); 
   return TCL_ERROR;
   }
 }

  if (vtkKWApplicationCppCommand(static_cast<vtkKWApplication *>(op),interp,argc,argv) == TCL_OK)
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
