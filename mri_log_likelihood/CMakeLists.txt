project(mri_log_likelihood)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_log_likelihood mri_log_likelihood.cpp)
target_link_libraries(mri_log_likelihood utils)

install(TARGETS mri_log_likelihood DESTINATION bin)
