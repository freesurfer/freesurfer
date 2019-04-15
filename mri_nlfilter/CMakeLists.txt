project(mri_nlfilter)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_nlfilter mri_nlfilter.cpp)
target_link_libraries(mri_nlfilter utils)

install(TARGETS mri_nlfilter DESTINATION bin)
