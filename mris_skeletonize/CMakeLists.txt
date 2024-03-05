project(mris_skeletonize)

# compile and install rules
include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_skeletonize mris_skeletonize.cpp)
add_help(mris_skeletonize mris_skeletonize.help.xml)
target_link_libraries(mris_skeletonize utils)
install(TARGETS mris_skeletonize DESTINATION bin)
