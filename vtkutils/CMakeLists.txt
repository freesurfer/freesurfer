project(vtkutils)
include_directories(${vtkutils_SOURCE_DIR}
${INCLUDE_DIR_TOP} 
${VXL_INCLUDES} 
${MINC_INCLUDE_DIRS}) 

execute_process ( 
COMMAND ${VTK_RUNTIME_LIBRARY_DIRS}/vtkWrapTcl vtkKWRGBATransferFunctionEditor.h ${VTK_WRAP_HINTS} 1 vtkKWRGBATransferFunctionEditorTcl.cxx
WORKING_DIRECTORY ${vtkutils_SOURCE_DIR} 
)

execute_process ( 
COMMAND ${VTK_RUNTIME_LIBRARY_DIRS}/vtkWrapTcl vtkKWBltGraph.h ${VTK_WRAP_HINTS} 1 vtkKWBltGraphTcl.cxx
WORKING_DIRECTORY ${vtkutils_SOURCE_DIR} 
)

SET(vtkutils_SRCS
	vtkKWProgressDialog.cxx
	vtkKWRGBATransferFunctionEditor.cxx
	vtkKWRGBATransferFunctionEditorTcl.cxx
	vtkKWBltGraph.cxx vtkKWBltGraphTcl.cxx
	vtkFreesurferLookupTable.cxx
	vtkArrowPipeline.cxx
	vtkRGBATransferFunction.cxx
	vtkFDTensorGlyph.cxx
	vtkODFGlyph.cxx
	vtkInflatePolyData.cxx
	IconLoader.cxx
)
add_library (vtkutils ${vtkutils_SRCS})

SET (IconLoaderTest_SOURCES
	IconLoaderTest.cxx 
)
execute_process ( 
COMMAND ${VTK_RUNTIME_LIBRARY_DIRS}/vtkWrapTcl vtkKWRGBATransferFunctionEditorTester.h ${VTK_WRAP_HINTS} 1 vtkKWRGBATransferFunctionEditorTesterTcl.cxx
WORKING_DIRECTORY ${vtkutils_SOURCE_DIR} 
)
execute_process ( 
COMMAND ${VTK_RUNTIME_LIBRARY_DIRS}/vtkWrapTclInit RgbatransferfunctioneditortesterLibInit.data RgbatransferfunctioneditortesterLibInit.cxx
WORKING_DIRECTORY ${vtkutils_SOURCE_DIR} 
)

SET (vtkKWRGBATransferFunctionEditorTester_SRCS
	vtkKWRGBATransferFunctionEditorTester.cxx
	vtkKWRGBATransferFunctionEditorTesterTcl.cxx
	RgbatransferfunctioneditortesterLibInit.cxx
)

SET (MVTKUTILS_TEST_LIBS
vtkfsio
vtkutils
${FS_LIBS}
)

add_executable(IconLoaderTest ${IconLoaderTest_SOURCES})
target_link_libraries(IconLoaderTest ${MVTKUTILS_TEST_LIBS})

add_executable(vtkKWRGBATransferFunctionEditorTester ${vtkKWRGBATransferFunctionEditorTester_SRCS})
target_link_libraries(vtkKWRGBATransferFunctionEditorTester ${MVTKUTILS_TEST_LIBS})

add_executable(vtkArrowPipelineTester vtkArrowPipelineTester.cxx)
target_link_libraries(vtkArrowPipelineTester ${MVTKUTILS_TEST_LIBS})

