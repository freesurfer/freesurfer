project(mris_sample_parc)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_sample_parc mris_sample_parc.cpp)
target_link_libraries(mris_sample_parc utils)

install(TARGETS mris_sample_parc DESTINATION bin)
