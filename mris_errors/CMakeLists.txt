project(mris_errors)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_errors mris_errors.cpp)
target_link_libraries(mris_errors utils)

install(TARGETS mris_errors DESTINATION bin)
