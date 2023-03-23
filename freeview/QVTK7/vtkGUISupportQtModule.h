#ifndef VTKGUISUPPORTQT_EXPORT_H
#define VTKGUISUPPORTQT_EXPORT_H

#define VTKGUISUPPORTQT_EXPORT
#define VTKGUISUPPORTQT_NO_EXPORT

#include "vtkRenderingCoreModule.h"

#include "qsystemdetection.h"
#if defined(Q_OS_OSX) && defined(ARM64)
#  include "vtkRenderingOpenGLModule.h"
#elif Q_OS_OSX
#  include "vtkRenderingOpenGL2Module.h"
#else
#  include "vtkRenderingOpenGLModule.h"
#endif

#endif
