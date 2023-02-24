project(samsegmesh2surf)

# compile and install rules
include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/gems)
add_executable(samsegmesh2surf samsegmesh2surf.cpp)
target_link_libraries(samsegmesh2surf kvlGEMSCommon utils)
install(TARGETS samsegmesh2surf DESTINATION bin)
