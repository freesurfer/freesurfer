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

int vtkKWScubaAppCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaAppNewCommand();
int vtkKWScubaApplicationSettingsInterfaceCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaApplicationSettingsInterfaceNewCommand();
int vtkKWScubaWindowCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaWindowNewCommand();
int vtkKWScubaViewCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaViewNewCommand();
int vtkKWScubaLayerCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayerNewCommand();
int vtkKWScubaLayerCollectionCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayerCollectionNewCommand();
int vtkKWScubaLayerCollectionMRICommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayerCollectionMRINewCommand();
int vtkKWScubaLayerCollectionMRISCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayerCollectionMRISNewCommand();
int vtkKWScubaLayerCollectionDTICommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayerCollectionDTINewCommand();
int vtkKWScubaLayerCollectionPathCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayerCollectionPathNewCommand();
int vtkKWScubaLayer2DMRICommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer2DMRINewCommand();
int vtkKWScubaLayer3DMRICommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer3DMRINewCommand();
int vtkKWScubaLayer2DMRISCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer2DMRISNewCommand();
int vtkKWScubaLayer3DMRISCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer3DMRISNewCommand();
int vtkKWScubaLayer2DDTICommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer2DDTINewCommand();
int vtkKWScubaLayer3DDTICommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer3DDTINewCommand();
int vtkKWScubaLayer2DPathCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer2DPathNewCommand();
int vtkKWScubaLayer3DPathCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaLayer3DPathNewCommand();
int vtkKWScubaToolCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaToolNewCommand();
int vtkKWScubaToolNavigateCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaToolNavigateNewCommand();
int vtkKWScubaToolEdit2DMRICommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWScubaToolEdit2DMRINewCommand();
int vtkKWRGBATransferFunctionEditorCommand(ClientData cd, Tcl_Interp *interp,
             int argc, char *argv[]);
ClientData vtkKWRGBATransferFunctionEditorNewCommand();

extern Tcl_HashTable vtkInstanceLookup;
extern Tcl_HashTable vtkPointerLookup;
extern Tcl_HashTable vtkCommandLookup;
extern void vtkTclDeleteObjectFromHash(void *);
extern void vtkTclListInstances(Tcl_Interp *interp, ClientData arg);


extern "C" {int VTK_EXPORT Scubalib_SafeInit(Tcl_Interp *interp);}

extern "C" {int VTK_EXPORT Scubalib_Init(Tcl_Interp *interp);}

extern void vtkTclGenericDeleteObject(ClientData cd);


int VTK_EXPORT Scubalib_SafeInit(Tcl_Interp *interp)
{
  return Scubalib_Init(interp);
}


int VTK_EXPORT Scubalib_Init(Tcl_Interp *interp)
{

  vtkTclCreateNew(interp,(char *) "vtkKWScubaApp", vtkKWScubaAppNewCommand,
                  vtkKWScubaAppCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaApplicationSettingsInterface", vtkKWScubaApplicationSettingsInterfaceNewCommand,
                  vtkKWScubaApplicationSettingsInterfaceCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaWindow", vtkKWScubaWindowNewCommand,
                  vtkKWScubaWindowCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaView", vtkKWScubaViewNewCommand,
                  vtkKWScubaViewCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer", vtkKWScubaLayerNewCommand,
                  vtkKWScubaLayerCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayerCollection", vtkKWScubaLayerCollectionNewCommand,
                  vtkKWScubaLayerCollectionCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayerCollectionMRI", vtkKWScubaLayerCollectionMRINewCommand,
                  vtkKWScubaLayerCollectionMRICommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayerCollectionMRIS", vtkKWScubaLayerCollectionMRISNewCommand,
                  vtkKWScubaLayerCollectionMRISCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayerCollectionDTI", vtkKWScubaLayerCollectionDTINewCommand,
                  vtkKWScubaLayerCollectionDTICommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayerCollectionPath", vtkKWScubaLayerCollectionPathNewCommand,
                  vtkKWScubaLayerCollectionPathCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer2DMRI", vtkKWScubaLayer2DMRINewCommand,
                  vtkKWScubaLayer2DMRICommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer3DMRI", vtkKWScubaLayer3DMRINewCommand,
                  vtkKWScubaLayer3DMRICommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer2DMRIS", vtkKWScubaLayer2DMRISNewCommand,
                  vtkKWScubaLayer2DMRISCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer3DMRIS", vtkKWScubaLayer3DMRISNewCommand,
                  vtkKWScubaLayer3DMRISCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer2DDTI", vtkKWScubaLayer2DDTINewCommand,
                  vtkKWScubaLayer2DDTICommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer3DDTI", vtkKWScubaLayer3DDTINewCommand,
                  vtkKWScubaLayer3DDTICommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer2DPath", vtkKWScubaLayer2DPathNewCommand,
                  vtkKWScubaLayer2DPathCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaLayer3DPath", vtkKWScubaLayer3DPathNewCommand,
                  vtkKWScubaLayer3DPathCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaTool", vtkKWScubaToolNewCommand,
                  vtkKWScubaToolCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaToolNavigate", vtkKWScubaToolNavigateNewCommand,
                  vtkKWScubaToolNavigateCommand);
  vtkTclCreateNew(interp,(char *) "vtkKWScubaToolEdit2DMRI", vtkKWScubaToolEdit2DMRINewCommand,
                  vtkKWScubaToolEdit2DMRICommand);
  vtkTclCreateNew(interp,(char *) "vtkKWRGBATransferFunctionEditor", vtkKWRGBATransferFunctionEditorNewCommand,
                  vtkKWRGBATransferFunctionEditorCommand);
  char pkgName[]="ScubaLib";
  char pkgVers[]=VTK_TCL_TO_STRING(VTK_MAJOR_VERSION) "." VTK_TCL_TO_STRING(VTK_MINOR_VERSION);
  Tcl_PkgProvide(interp, pkgName, pkgVers);
  return TCL_OK;
}
