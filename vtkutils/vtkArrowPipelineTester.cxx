#include "vtkArrowPipeline.h"

#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"

int main ( int argc, char** argv ) {

  // Make a renderer and a window.
  vtkRenderer* renderer = vtkRenderer::New();
  vtkRenderWindow* window = vtkRenderWindow::New();
  window->AddRenderer( renderer );

  // Make an interactor.
  vtkRenderWindowInteractor* interactor = vtkRenderWindowInteractor::New();
  interactor->SetRenderWindow( window );

  // Make an arrow and add it.
  vtkArrowPipeline* arrow = vtkArrowPipeline::New();
  double point[3] = { 0, 0, 0 };
  arrow->SetStartPoint( point );
  point[0] = point[1] = point[2] = 10;
  arrow->SetEndPoint( point );

  renderer->AddActor( arrow->GetActor() );

  // Let's do this.
  window->Render();
  interactor->Start();

  exit( 0 );
}
