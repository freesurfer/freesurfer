project(swi_processing)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom ${CMAKE_SOURCE_DIR}/packages/dcm2niix)

add_executable(swi_preprocess swi_preprocess.cpp cmd_line_interface.cpp)
target_link_libraries(swi_preprocess utils)
install(TARGETS swi_preprocess DESTINATION bin)

add_executable(swi_process swi_process.cpp cmd_line_interface.cpp)
target_link_libraries(swi_process utils)
install(TARGETS swi_process DESTINATION bin)

