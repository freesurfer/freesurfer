project(itkutils)

SET (BSPLINES_SRCS
	bspline/itkBSplineKernelFunction.txx 
	bspline/itkBSplineKernelFunction.h 
	bspline/itkBSplineDataPointSetToImageFilter.txx 
    bspline/itkBSplineDataPointSetToImageFilter.h 
	bspline/itkPointSetToImageFilter.txx 
	bspline/itkPointSetToImageFilter.h
)

SET (UI_SRCS
  ui/CommandParser.cxx 
  ui/FreeSurferExecutable.cxx
 )
 
 SET (itkutils_SRC 
 ${UI_SRCS}
 ${BSPLINES_SRCS}
 )

add_library (itkutils ${itkutils_SRC})

