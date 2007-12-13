// tcl wrapper for vtkKWQdecApp object
//
#define VTK_WRAPPING_CXX
#define VTK_STREAMS_FWD_ONLY
#include "vtkSystemIncludes.h"
#include "vtkKWQdecApp.h"

#include "vtkTclUtil.h"
#include <vtkstd/stdexcept>

ClientData vtkKWQdecAppNewCommand()
{
  vtkKWQdecApp *temp = vtkKWQdecApp::New();
  return ((ClientData)temp);
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
   return vtkKWQdecAppCppCommand((vtkKWQdecApp *)(((vtkTclCommandArgStruct *)cd)->Pointer),interp, argc, argv);
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
    Tcl_SetResult(interp, (char *) "Could not find requested method.", TCL_VOLATILE);
    return TCL_ERROR;
    }
  if (!interp)
    {
    if (!strcmp("DoTypecasting",argv[0]))
      {
      if (!strcmp("vtkKWQdecApp",argv[1]))
        {
        argv[2] = (char *)((void *)op);
        return TCL_OK;
        }
      if (vtkKWApplicationCppCommand((vtkKWApplication *)op,interp,argc,argv) == TCL_OK)
        {
        return TCL_OK;
        }
      }
    return TCL_ERROR;
    }

  if (!strcmp("GetSuperClassName",argv[1]))
    {
    Tcl_SetResult(interp,(char *) "vtkKWApplication", TCL_VOLATILE);
    return TCL_OK;
    }

  try
    {
  if ((!strcmp("New",argv[1]))&&(argc == 2))
    {
    vtkKWQdecApp  *temp20;
    temp20 = (op)->New();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecApp");
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
    vtkKWQdecApp  *temp20;
    temp20 = (op)->NewInstance();
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecApp");
    return TCL_OK;
    }
  if ((!strcmp("SafeDownCast",argv[1]))&&(argc == 3))
    {
    vtkObject  *temp0;
    vtkKWQdecApp  *temp20;
    error = 0;

    temp0 = (vtkObject *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkObject",interp,error));
    if (!error)
    {
    temp20 = (op)->SafeDownCast(temp0);
      vtkTclGetObjectFromPointer(interp,(void *)temp20,"vtkKWQdecApp");
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

    temp0 = (vtkKWTopLevel *)(vtkTclGetPointerFromObject(argv[2],(char *) "vtkKWTopLevel",interp,error));
    if (!error)
    {
    op->DisplayHelpDialog(temp0);
    Tcl_ResetResult(interp);
    return TCL_OK;
    }
    }

  if (!strcmp("ListInstances",argv[1]))
    {
    vtkTclListInstances(interp,(ClientData)vtkKWQdecAppCommand);
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
    Tcl_AppendResult(interp,"  ErrorMessage\t with 1 arg\n",NULL);
    Tcl_AppendResult(interp,"  DisplayHelpDialog\t with 1 arg\n",NULL);
    return TCL_OK;
    }

  if (vtkKWApplicationCppCommand((vtkKWApplication *)op,interp,argc,argv) == TCL_OK)
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
