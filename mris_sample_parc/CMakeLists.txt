project(mris_sample_parc)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_sample_parc mris_sample_parc.c)
target_link_libraries(mris_sample_parc utils)

install(TARGETS mris_sample_parc DESTINATION bin)
