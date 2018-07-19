project(mri_fieldsign)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_fieldsign mri_fieldsign.c)
target_link_libraries(mri_fieldsign utils)

install(TARGETS mri_fieldsign DESTINATION bin)
