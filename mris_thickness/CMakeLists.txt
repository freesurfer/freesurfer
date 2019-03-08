project(mris_thickness)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_thickness mris_thickness.cpp)
target_link_libraries(mris_thickness utils)
install(TARGETS mris_thickness DESTINATION bin)

add_executable(mris_intensity_profile mris_intensity_profile.cpp)
target_link_libraries(mris_intensity_profile utils)
install(TARGETS mris_intensity_profile DESTINATION bin)

add_executable(mris_gradient mris_gradient.cpp)
target_link_libraries(mris_gradient utils)
install(TARGETS mris_gradient DESTINATION bin)
