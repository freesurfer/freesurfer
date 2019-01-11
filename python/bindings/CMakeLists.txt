include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# temporary fix so that -g doesn't cause linker errors
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

set(FS_PYTHON_DIR ${CMAKE_SOURCE_DIR}/python/freesurfer)

add_subdirectories(algorithm)
