project(mri_gradunwarp)

include_directories(${FS_INCLUDE_DIRS})

set(TARGET_EXE mri_gradunwarp)

add_executable(${TARGET_EXE} 
    mri_gradunwarp.cpp)

target_link_libraries(${TARGET_EXE} utils)

install(TARGETS ${TARGET_EXE} DESTINATION bin)
