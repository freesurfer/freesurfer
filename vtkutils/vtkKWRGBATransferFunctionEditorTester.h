#ifndef testvtkKWRGBATransferFunctionEditor_h
#define testvtkKWRGBATransferFunctionEditor_h

#include "vtkKWApplication.h"

class vtkKWRGBATransferFunctionEditor;

class vtkKWRGBATransferFunctionEditorTester : public vtkKWApplication {

 public:

  static vtkKWRGBATransferFunctionEditorTester* New ();
  vtkTypeRevisionMacro( vtkKWRGBATransferFunctionEditorTester, vtkKWApplication );

  virtual void Start ( int argc, char* argv[] );

  void EditorChangedFunction ();

 protected:

  vtkKWRGBATransferFunctionEditorTester () {}
  ~vtkKWRGBATransferFunctionEditorTester () {}
  
  vtkKWRGBATransferFunctionEditor* mEditor;
};

#endif
