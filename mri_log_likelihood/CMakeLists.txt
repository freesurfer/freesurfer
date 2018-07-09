project(mri_log_likelihood)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_log_likelihood mri_log_likelihood.c)
target_link_libraries(mri_log_likelihood utils)

install(TARGETS mri_log_likelihood DESTINATION bin)
