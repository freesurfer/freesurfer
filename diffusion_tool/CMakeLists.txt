project(diffusion_tool)

if(BUILD_GUIS AND Qt5_FOUND AND VTK_FOUND AND ITK_FOUND AND PETSC_FOUND AND TARGET vtkutils)

	include_directories(
		${CMAKE_CURRENT_SOURCE_DIR}
		${QT_INCLUDE_DIRS}
		)

	set(CMAKE_AUTOUIC ON)
	set(CMAKE_AUTOMOC ON)
	set(CMAKE_AUTORCC ON)

	set(SOURCES
		main.cpp
		configurationfileform.cpp
		)

	add_executable(diffusion_tool ${SOURCES})
	target_link_libraries(diffusion_tool
		Qt5::Core
		Qt5::Widgets
		)
	
	endif()
