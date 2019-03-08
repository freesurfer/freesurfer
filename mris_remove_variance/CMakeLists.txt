project(mris_remove_variance)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_remove_variance mris_remove_variance.cpp)
target_link_libraries(mris_remove_variance utils)

install(TARGETS mris_remove_variance DESTINATION bin)
