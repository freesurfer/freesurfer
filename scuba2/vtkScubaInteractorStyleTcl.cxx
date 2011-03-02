// tcl wrapper for vtkScubaInteractorStyle object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkScubaInteractorStyle.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>
#include <vtksys/ios/sstream>

ClientData vtkScubaInteractorStyleNewCommand()
{
  vtkScubaInteractorStyle *temp = vtkScubaInteractorStyle::New();
  return static_cast<ClientData>(temp);
}

int vtkInteractorStyleCppCommand(vtkInteractorStyle *op, Tcl_Interp *interp,
             int argc, char *argv[]);
int VTKTCL_EXPORT vtkScubaInteractorStyleCppCommand(vtkScubaInteractorStyle *op, Tcl_Interp *interp,
             int argc, char *argv[]);

int VTKTCL_EXPORT vtkScubaInteractorStyleCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[])
{
  if ((argc == 2)&&(!strcmp("Delete",argv[1]))&& !vtkTclInDelete(interp))
    {
    Tcl_DeleteCommand(interp,argv[0]);
    return TCL_OK;
    }
   return vtkScubaInteractorStyleCppCommand(static_cast<vtkScubaInteractorStyle *>(static_cast<vtkTclCommandArgStruct *>(cd)->Pointer),interp, argc, argv);
}

int VTKTCL_EXPORT vtkScubaInteractorStyleCppCommand(vtkScubaInteractorStyle *op, Tcl_Interp *interp,
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
      if (!strcmp("vtkScubaInteractorStyle",argv[1]))
        {
        argv[2] = static_cast<char *>(static_cast<void *>(op));
        return TCL_OK;
        }
      if (vtkInteractorStyleCppCommand(static_cast<vtkInteractorStyle *>(op),interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,const_cast<char *>("vtkInteractorStyle"), TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkScubaInteractorStyle  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkScubaInteractorStyle");
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
    vtkScubaInteractorStyle  *temp20;
    temp20 = (op)->NewInstance();
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkScubaInteractorStyle");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkScubaInteractorStyle  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkObject"),interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)(temp20),"vtkScubaInteractorStyle");
    return TCL_OK;
    }
    }
  if ((!strcmp("SetWindow",argv[1]))&&(argc == 3))
    {
    vtkKWScubaWindow  *temp0;
    error = 0;

    temp0 = (vtkKWScubaWindow *)(vtkTclGetPointerFromObject(argv[2],const_cast<char *>("vtkKWScubaWindow"),interp,error));
    if (!error)
    {
    op->SetWindow(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }
  if ((!strcmp("OnMouseMove",argv[1]))&&(argc == 2))
    {
    op->OnMouseMove();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnLeftButtonDown",argv[1]))&&(argc == 2))
    {
    op->OnLeftButtonDown();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnLeftButtonUp",argv[1]))&&(argc == 2))
    {
    op->OnLeftButtonUp();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnMiddleButtonDown",argv[1]))&&(argc == 2))
    {
    op->OnMiddleButtonDown();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnMiddleButtonUp",argv[1]))&&(argc == 2))
    {
    op->OnMiddleButtonUp();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnRightButtonDown",argv[1]))&&(argc == 2))
    {
    op->OnRightButtonDown();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnRightButtonUp",argv[1]))&&(argc == 2))
    {
    op->OnRightButtonUp();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnKeyDown",argv[1]))&&(argc == 2))
    {
    op->OnKeyDown();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnKeyUp",argv[1]))&&(argc == 2))
    {
    op->OnKeyUp();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnEnter",argv[1]))&&(argc == 2))
    {
    op->OnEnter();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
  if ((!strcmp("OnLeave",argv[1]))&&(argc == 2))
    {
    op->OnLeave();
    Tcl_ResetResult(interp);
    return TCL_OK;
    }

  if (!strcmp("ListInstances",argv[1]))
    {
    vtkTclListInstances(interp,(ClientData)(vtkScubaInteractorStyleCommand));
    return TCL_OK;
    }

  if (!strcmp("ListMethods",argv[1]))
    {
    vtkInteractorStyleCppCommand(op,interp,argc,argv);
    Tcl_AppendResult(interp,"Methods from vtkScubaInteractorStyle:\n",NULL);
    Tcl_AppendResult(interp,"  GetSuperClassName\n",NULL);
    Tcl_AppendResult(interp,"  New\n",NULL);
    Tcl_AppendResult(interp,"  GetClassName\n",NULL);
    Tcl_AppendResult(interp,"  IsA\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  NewInstance\n",NULL);
    Tcl_AppendResult(interp,"  SafeDownCast\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  SetWindow\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  OnMouseMove\n",NULL);
    Tcl_AppendResult(interp,"  OnLeftButtonDown\n",NULL);
    Tcl_AppendResult(interp,"  OnLeftButtonUp\n",NULL);
    Tcl_AppendResult(interp,"  OnMiddleButtonDown\n",NULL);
    Tcl_AppendResult(interp,"  OnMiddleButtonUp\n",NULL);
    Tcl_AppendResult(interp,"  OnRightButtonDown\n",NULL);
    Tcl_AppendResult(interp,"  OnRightButtonUp\n",NULL);
    Tcl_AppendResult(interp,"  OnKeyDown\n",NULL);
    Tcl_AppendResult(interp,"  OnKeyUp\n",NULL);
    Tcl_AppendResult(interp,"  OnEnter\n",NULL);
    Tcl_AppendResult(interp,"  OnLeave\n",NULL);
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
    vtkInteractorStyleCppCommand(op,interp,argc,argv);
    Tcl_DStringGetResult ( interp, &dStringParent );
    Tcl_DStringAppend ( &dString, Tcl_DStringValue ( &dStringParent ), -1 );
    Tcl_DStringAppendElement ( &dString, "New" );
    Tcl_DStringAppendElement ( &dString, "GetClassName" );
    Tcl_DStringAppendElement ( &dString, "IsA" );
    Tcl_DStringAppendElement ( &dString, "NewInstance" );
    Tcl_DStringAppendElement ( &dString, "SafeDownCast" );
    Tcl_DStringAppendElement ( &dString, "SetWindow" );
    Tcl_DStringAppendElement ( &dString, "OnMouseMove" );
    Tcl_DStringAppendElement ( &dString, "OnLeftButtonDown" );
    Tcl_DStringAppendElement ( &dString, "OnLeftButtonUp" );
    Tcl_DStringAppendElement ( &dString, "OnMiddleButtonDown" );
    Tcl_DStringAppendElement ( &dString, "OnMiddleButtonUp" );
    Tcl_DStringAppendElement ( &dString, "OnRightButtonDown" );
    Tcl_DStringAppendElement ( &dString, "OnRightButtonUp" );
    Tcl_DStringAppendElement ( &dString, "OnKeyDown" );
    Tcl_DStringAppendElement ( &dString, "OnKeyUp" );
    Tcl_DStringAppendElement ( &dString, "OnEnter" );
    Tcl_DStringAppendElement ( &dString, "OnLeave" );
  Tcl_DStringResult ( interp, &dString );
  Tcl_DStringFree ( &dString );
  Tcl_DStringFree ( &dStringParent );
    return TCL_OK;
    }
    if(argc==3) {
      Tcl_DString dString;
      int SuperClassStatus;
    SuperClassStatus = vtkInteractorStyleCppCommand(op,interp,argc,argv);
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
    Tcl_DStringAppendElement ( &dString, "static vtkScubaInteractorStyle *New ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
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
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
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
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
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
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle *NewInstance ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
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
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle *SafeDownCast (vtkObject* o);" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for SafeDownCast */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: SetWindow */
    if ( strcmp ( argv[2], "SetWindow" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "SetWindow" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringAppendElement ( &dString, "vtkKWScubaWindow" );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for SetWindow */
    Tcl_DStringAppendElement ( &dString, " Set the window. When an event happens, the window will be passed\n along to the view with the event call. It will also tell the\n window when an event is done being handled.\n" );
    Tcl_DStringAppendElement ( &dString, "void SetWindow (vtkKWScubaWindow *iWindow);" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for SetWindow */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnMouseMove */
    if ( strcmp ( argv[2], "OnMouseMove" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnMouseMove" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnMouseMove */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnMouseMove ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnMouseMove */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnLeftButtonDown */
    if ( strcmp ( argv[2], "OnLeftButtonDown" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnLeftButtonDown" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnLeftButtonDown */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnLeftButtonDown ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnLeftButtonDown */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnLeftButtonUp */
    if ( strcmp ( argv[2], "OnLeftButtonUp" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnLeftButtonUp" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnLeftButtonUp */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnLeftButtonUp ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnLeftButtonUp */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnMiddleButtonDown */
    if ( strcmp ( argv[2], "OnMiddleButtonDown" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnMiddleButtonDown" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnMiddleButtonDown */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnMiddleButtonDown ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnMiddleButtonDown */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnMiddleButtonUp */
    if ( strcmp ( argv[2], "OnMiddleButtonUp" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnMiddleButtonUp" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnMiddleButtonUp */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnMiddleButtonUp ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnMiddleButtonUp */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnRightButtonDown */
    if ( strcmp ( argv[2], "OnRightButtonDown" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnRightButtonDown" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnRightButtonDown */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnRightButtonDown ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnRightButtonDown */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnRightButtonUp */
    if ( strcmp ( argv[2], "OnRightButtonUp" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnRightButtonUp" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnRightButtonUp */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnRightButtonUp ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnRightButtonUp */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnKeyDown */
    if ( strcmp ( argv[2], "OnKeyDown" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnKeyDown" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnKeyDown */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnKeyDown ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnKeyDown */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnKeyUp */
    if ( strcmp ( argv[2], "OnKeyUp" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnKeyUp" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnKeyUp */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnKeyUp ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnKeyUp */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnEnter */
    if ( strcmp ( argv[2], "OnEnter" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnEnter" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnEnter */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnEnter ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnEnter */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
    /* Starting function: OnLeave */
    if ( strcmp ( argv[2], "OnLeave" ) == 0 ) {
    Tcl_DStringInit ( &dString );
    Tcl_DStringAppendElement ( &dString, "OnLeave" );
    /* Arguments */
    Tcl_DStringStartSublist ( &dString );
    Tcl_DStringEndSublist ( &dString );
    /* Documentation for OnLeave */
    Tcl_DStringAppendElement ( &dString, " Implement vtkInteractorStyle callbacks so we can pass events our\n way.\n" );
    Tcl_DStringAppendElement ( &dString, "virtual void OnLeave ();" );
    Tcl_DStringAppendElement ( &dString, "vtkScubaInteractorStyle" );
    /* Closing for OnLeave */

    Tcl_DStringResult ( interp, &dString );
    Tcl_DStringFree ( &dString );
    return TCL_OK;
    }
   Tcl_SetResult ( interp, const_cast<char*>("Could not find method"), TCL_VOLATILE ); 
   return TCL_ERROR;
   }
 }

  if (vtkInteractorStyleCppCommand(static_cast<vtkInteractorStyle *>(op),interp,argc,argv) == TCL_OK)
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
