project(mri_otl)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom ${CMAKE_SOURCE_DIR}/packages/dcm2niix)

add_executable(list_otl_labels list_otl_labels.cpp)
target_link_libraries(list_otl_labels utils)
install(TARGETS list_otl_labels DESTINATION bin)
