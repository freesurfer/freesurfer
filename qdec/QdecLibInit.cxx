#include "vtkTclUtil.h"
#include "vtkVersion.h"
#define VTK_TCL_TO_STRING(x) VTK_TCL_TO_STRING0(x)
#define VTK_TCL_TO_STRING0(x) #x
extern "C"
{
#if (TCL_MAJOR_VERSION == 8) && (TCL_MINOR_VERSION >= 4) && (TCL_RELEASE_LEVEL >= TCL_FINAL_RELEASE)
  typedef int (*vtkTclCommandType)(ClientData, Tcl_Interp *,int, CONST84 char *[]);
#else
  typedef int (*vtkTclCommandType)(ClientData, Tcl_Interp *,int, char *[]);
#endif
}

int vtkKWQdecAppCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWQdecAppNewCommand();
int vtkKWQdecWindowCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWQdecWindowNewCommand();
int vtkKWQdecViewCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWQdecViewNewCommand();
int vtkKWRGBATransferFunctionEditorCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWRGBATransferFunctionEditorNewCommand();
int vtkKWBltGraphCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWBltGraphNewCommand();

extern Tcl_HashTable vtkInstanceLookup;
extern Tcl_HashTable vtkPointerLookup;
extern Tcl_HashTable vtkCommandLookup;
extern void vtkTclDeleteObjectFromHash(void *);
extern void vtkTclListInstances(Tcl_Interp *interp, ClientData arg);


extern "C" {int VTK_EXPORT Qdeclib_SafeInit(Tcl_Interp *interp);}

extern "C" {int VTK_EXPORT Qdeclib_Init(Tcl_Interp *interp);}

extern void vtkTclGenericDeleteObject(ClientData cd);


int VTK_EXPORT Qdeclib_SafeInit(Tcl_Interp *interp)
{
  return Qdeclib_Init(interp);
}


int VTK_EXPORT Qdeclib_Init(Tcl_Interp *interp)
{

  vtkTclCreateNew(interp,(char *) "vtkKWQdecApp", vtkKWQdecAppNewCommand,
                  vtkKWQdecAppCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWQdecWindow", vtkKWQdecWindowNewCommand,
                  vtkKWQdecWindowCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWQdecView", vtkKWQdecViewNewCommand,
                  vtkKWQdecViewCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWRGBATransferFunctionEditor", vtkKWRGBATransferFunctionEditorNewCommand,
                  vtkKWRGBATransferFunctionEditorCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWBltGraph", vtkKWBltGraphNewCommand,
                  vtkKWBltGraphCommand);
  char pkgName[]="QdecLib";
  char pkgVers[]=VTK_TCL_TO_STRING(VTK_MAJOR_VERSION) "." VTK_TCL_TO_STRING(VTK_MINOR_VERSION);
  Tcl_PkgProvide(interp, pkgName, pkgVers);
  return TCL_OK;
}
