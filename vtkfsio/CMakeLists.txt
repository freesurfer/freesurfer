project(vtkfsio)
include_directories(${vtkfsio_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

SET(vtkfsio_SRCS
	vtkFSSurfaceLabelSource.cxx
	vtkFSSurfaceLabelSource.h
	vtkFSSurfaceWriter.cxx
	vtkFSSurfaceWriter.h
	vtkFSSurfaceScalarsReader.cxx
	vtkFSSurfaceScalarsReader.h
	vtkFSSurfaceSource.cxx
	vtkFSSurfaceSource.h
	vtkFSVolumeSource.cxx
	vtkFSVolumeSource.h
)

SET (MVTKFSIO_TEST_LIBS
vtkfsio
${FS_LIBS}
)

add_library(vtkfsio ${vtkfsio_SRCS})

add_executable (vtkFSSurfaceWriterTester vtkFSSurfaceWriterTester.cxx)
target_link_libraries(vtkFSSurfaceWriterTester ${MVTKFSIO_TEST_LIBS})

