// tcl wrapper for vtkScubaInteractorStyle object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkScubaInteractorStyle.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>

ClientData vtkScubaInteractorStyleNewCommand()
{
  vtkScubaInteractorStyle *temp = vtkScubaInteractorStyle::New();
  return ((ClientData)temp);
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
   return vtkScubaInteractorStyleCppCommand((vtkScubaInteractorStyle *)(((vtkTclCommandArgStruct *)cd)->Pointer),interp, argc, argv);
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
    Tcl_SetResult(interp, (char *) "Could not find requested method.", TCL_VOLATILE);
    return TCL_ERROR;
    }
  if (!interp)
    {
    if (!strcmp("DoTypecasting",argv[0]))
      {
      if (!strcmp("vtkScubaInteractorStyle",argv[1]))
        {
        argv[2] = (char *)((void *)op);
        return TCL_OK;
        }
      if (vtkInteractorStyleCppCommand((vtkInteractorStyle *)op,interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,(char *) "vtkInteractorStyle", TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkScubaInteractorStyle  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkScubaInteractorStyle");
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
    vtkScubaInteractorStyle  *temp20;
    temp20 = (op)->NewInstance();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkScubaInteractorStyle");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkScubaInteractorStyle  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkObject",interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkScubaInteractorStyle");
    return TCL_OK;
    }
    }
  if ((!strcmp("SetWindow",argv[1]))&&(argc == 3))
    {
    vtkKWScubaWindow  *temp0;
    error = 0;

    temp0 = (vtkKWScubaWindow *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkKWScubaWindow",interp,error));
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
    vtkTclListInstances(interp,(ClientData)vtkScubaInteractorStyleCommand);
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

  if (vtkInteractorStyleCppCommand((vtkInteractorStyle *)op,interp,argc,argv) == TCL_OK)
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
