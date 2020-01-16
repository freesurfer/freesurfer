include_directories(${CMAKE_CURRENT_SOURCE_DIR}/utils ${FS_INCLUDE_DIRS})  # TODO remove this??

# temporary fix so that -g doesn't produce linker errors
set(CMAKE_CXX_FLAGS_DEBUG "")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "")

add_subdirectory(utils)
add_subdirectory(gems)
