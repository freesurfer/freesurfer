project(mri_fieldsign)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_fieldsign mri_fieldsign.cpp)
target_link_libraries(mri_fieldsign utils)

install(TARGETS mri_fieldsign DESTINATION bin)
