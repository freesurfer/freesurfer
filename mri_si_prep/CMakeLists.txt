project(mri_si_prep)

# compile and install rules
if(NOT INTEGRATE_SAMSEG)
  include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/gems)
else()
  include_directories(${FS_INCLUDE_DIRS} ${STANDALONE_SAMSEG_PATH}/gems)
endif()
add_executable(mri_si_prep mri_si_prep.cpp)
add_help(mri_si_prep mri_si_prep.help.xml)
target_link_libraries(mri_si_prep utils kvlGEMSCommon)
install(TARGETS mri_si_prep DESTINATION bin)
