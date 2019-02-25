
#ifndef VTKIOIMAGE_EXPORT_H
#define VTKIOIMAGE_EXPORT_H

#ifdef VTKIOIMAGE_STATIC_DEFINE
#  define VTKIOIMAGE_EXPORT
#  define VTKIOIMAGE_NO_EXPORT
#else
#  ifndef VTKIOIMAGE_EXPORT
#    ifdef vtkIOImage_EXPORTS
        /* We are building this library */
#      define VTKIOIMAGE_EXPORT __attribute__((visibility("default")))
#    else
        /* We are using this library */
#      define VTKIOIMAGE_EXPORT __attribute__((visibility("default")))
#    endif
#  endif

#  ifndef VTKIOIMAGE_NO_EXPORT
#    define VTKIOIMAGE_NO_EXPORT __attribute__((visibility("hidden")))
#  endif
#endif

#ifndef VTKIOIMAGE_DEPRECATED
#  define VTKIOIMAGE_DEPRECATED __attribute__ ((__deprecated__))
#  define VTKIOIMAGE_DEPRECATED_EXPORT VTKIOIMAGE_EXPORT __attribute__ ((__deprecated__))
#  define VTKIOIMAGE_DEPRECATED_NO_EXPORT VTKIOIMAGE_NO_EXPORT __attribute__ ((__deprecated__))
#endif

#define DEFINE_NO_DEPRECATED 0
#if DEFINE_NO_DEPRECATED
# define VTKIOIMAGE_NO_DEPRECATED
#endif

/* AutoInit implementations.  */
#if defined(vtkIOImage_INCLUDE)
# include vtkIOImage_INCLUDE
#endif
#if defined(vtkIOImage_AUTOINIT)
# include "vtkAutoInit.h"
VTK_AUTOINIT(vtkIOImage)
#endif

#endif
