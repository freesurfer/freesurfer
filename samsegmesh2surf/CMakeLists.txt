project(samsegmesh2surf)

# compile and install rules
if(NOT INTEGRATE_SAMSEG)
  include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/gems)
else()
  include_directories(${FS_INCLUDE_DIRS} ${STANDALONE_SAMSEG_PATH}/gems)
endif()  
add_executable(samsegmesh2surf samsegmesh2surf.cpp)
add_help(samsegmesh2surf samsegmesh2surf.help.xml)
target_link_libraries(samsegmesh2surf utils kvlGEMSCommon)
install(TARGETS samsegmesh2surf DESTINATION bin)
