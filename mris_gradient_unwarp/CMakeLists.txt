project(mris_gradient_unwarp)

include_directories(${FS_INCLUDE_DIRS})

set(TARGET_EXE mris_gradient_unwarp)

add_executable(${TARGET_EXE} 
    mris_gradient_unwarp.cpp)

target_link_libraries(${TARGET_EXE} utils)

install(TARGETS ${TARGET_EXE} DESTINATION bin)
