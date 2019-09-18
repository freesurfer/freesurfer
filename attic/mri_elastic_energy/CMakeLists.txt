project(mri_elastic_energy)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_elastic_energy mri_elastic_energy.cpp)
target_link_libraries(mri_elastic_energy utils)

install(TARGETS mri_elastic_energy DESTINATION bin)
