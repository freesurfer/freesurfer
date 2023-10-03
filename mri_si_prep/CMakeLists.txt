project(mri_si_prep)

# compile and install rules
include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/gems)
add_executable(mri_si_prep mri_si_prep.cpp)
add_help(mri_si_prep mri_si_prep.help.xml)
target_link_libraries(mri_si_prep utils kvlGEMSCommon)
install(TARGETS mri_si_prep DESTINATION bin)
