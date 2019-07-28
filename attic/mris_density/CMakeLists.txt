project(mris_density)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_density mris_density.cpp)
target_link_libraries(mris_density utils)

install(TARGETS mris_density DESTINATION bin)
