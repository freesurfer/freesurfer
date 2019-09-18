project(mris_ms_refine)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_ms_refine mris_ms_refine.cpp)
target_link_libraries(mris_ms_refine utils)

install(TARGETS mris_ms_refine DESTINATION bin)
