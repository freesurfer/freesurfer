# FLTK's fluid generates C++ files in the build directory
INCLUDE_DIRECTORIES(${GEMS_SOURCE_DIR}/GUI)
INCLUDE_DIRECTORIES(${GEMS_BINARY_DIR}/GUI)


# Add libary
ADD_LIBRARY( kvlGEMSGUI
             vtkFlRenderWindowInteractor.cxx
             kvlImageViewer.cxx
           )
TARGET_LINK_LIBRARIES( kvlGEMSGUI
                       ${ITK_LIBRARIES}
                       ${VTK_LIBRARIES}
                       ${FLTK_LIBRARIES} )


# Add an executable
FLTK_WRAP_UI( kvlEstimateAtlasParameters
              kvlAtlasParameterEstimationConsoleGUI.fl )

ADD_EXECUTABLE( kvlEstimateAtlasParameters
                ${kvlEstimateAtlasParameters_FLTK_UI_SRCS}
                kvlEstimateAtlasParameters.cxx
                kvlAtlasParameterEstimationConsole.cxx )
TARGET_LINK_LIBRARIES( kvlEstimateAtlasParameters
                       kvlGEMSCommon
                       kvlGEMSGUI
                      )

# Add an executable
FLTK_WRAP_UI( kvlViewMeshCollectionWithGUI
              kvlAtlasMeshViewingConsoleGUI.fl )

ADD_EXECUTABLE( kvlViewMeshCollectionWithGUI
                ${kvlViewMeshCollectionWithGUI_FLTK_UI_SRCS}
                kvlViewMeshCollectionWithGUI.cxx
                kvlAtlasMeshViewingConsole.cxx )
TARGET_LINK_LIBRARIES( kvlViewMeshCollectionWithGUI
                        kvlGEMSCommon
                        kvlGEMSGUI
                        )


