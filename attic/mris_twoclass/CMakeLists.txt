project(mris_twoclass)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_twoclass mris_twoclass.cpp)
target_link_libraries(mris_twoclass utils)

install(TARGETS mris_twoclass DESTINATION bin)
