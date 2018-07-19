project(mris_compute_acorr)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_compute_acorr mris_compute_acorr.c)
target_link_libraries(mris_compute_acorr utils)

install(TARGETS mris_compute_acorr DESTINATION bin)
