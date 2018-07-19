project(mris_errors)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_errors mris_errors.c)
target_link_libraries(mris_errors utils)

install(TARGETS mris_errors DESTINATION bin)
