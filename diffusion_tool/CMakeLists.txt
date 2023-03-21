project(diffusion_tool)

if(BUILD_GUIS)

	include_directories(${CMAKE_CURRENT_SOURCE_DIR})

	set(CMAKE_AUTOUIC ON)
	set(CMAKE_AUTOMOC ON)

	add_executable(diffusion_tool main.cpp configurationfileform.cpp)

  	if(Qt6_DIR)
    		target_link_libraries(diffusion_tool Qt6::Core Qt6::Widgets)
  	elseif(Qt5_DIR)
    		target_link_libraries(diffusion_tool Qt5::Core Qt5::Widgets)
  	endif()

	
endif()

