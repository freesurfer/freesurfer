project(mris_compute_acorr)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_compute_acorr mris_compute_acorr.cpp)
target_link_libraries(mris_compute_acorr utils)

install(TARGETS mris_compute_acorr DESTINATION bin)
