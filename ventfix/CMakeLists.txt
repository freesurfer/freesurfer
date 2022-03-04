project(ventfix)

include_directories(${FS_INCLUDE_DIRS})

set(TARGET_EXE ventfix)

add_executable(${TARGET_EXE} 
    ventfix_main.cpp)

target_link_libraries(${TARGET_EXE} utils)

install(TARGETS ${TARGET_EXE} DESTINATION bin)
