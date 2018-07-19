project(mris_ms_refine)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_ms_refine mris_ms_refine.c)
target_link_libraries(mris_ms_refine utils)

install(TARGETS mris_ms_refine DESTINATION bin)
