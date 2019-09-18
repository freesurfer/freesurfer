project(dmri_tensoreig)

include_directories(${FS_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/packages/dicom)

add_executable(dmri_tensoreig dmri_tensoreig.cpp)
target_link_libraries(dmri_tensoreig utils)

install(TARGETS dmri_tensoreig DESTINATION bin)
