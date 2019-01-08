project(diffusion_tool)

if(BUILD_GUIS AND Qt5_FOUND)

	include_directories(${CMAKE_CURRENT_SOURCE_DIR})

	set(CMAKE_AUTOUIC ON)
	set(CMAKE_AUTOMOC ON)

	add_executable(diffusion_tool main.cpp configurationfileform.cpp)
	target_link_libraries(diffusion_tool Qt5::Core Qt5::Widgets)
	
endif()

