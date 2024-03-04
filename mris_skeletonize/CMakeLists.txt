project(mris_skeletonize)

# compile and install rules
if(NOT INTEGRATE_SAMSEG)
  include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/gems)
else()
  include_directories(${FS_INCLUDE_DIRS} ${STANDALONE_SAMSEG_PATH}/gems)
endif()
add_executable(mris_skeletonize mris_skeletonize.cpp)
add_help(mris_skeletonize mris_skeletonize.help.xml)
target_link_libraries(mris_skeletonize utils kvlGEMSCommon)
install(TARGETS mris_skeletonize DESTINATION bin)
