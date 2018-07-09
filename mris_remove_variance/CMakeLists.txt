project(mris_remove_variance)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_remove_variance mris_remove_variance.c)
target_link_libraries(mris_remove_variance utils)

install(TARGETS mris_remove_variance DESTINATION bin)
