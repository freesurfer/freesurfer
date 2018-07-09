project(mri_probedicom)

include_directories(${CMAKE_SOURCE_DIR}/include ${GLUT_INCLUDE_DIR}/GL)

add_executable(mri_probedicom mri_probedicom.c)
target_link_libraries(mri_probedicom utils)

install(TARGETS mri_probedicom DESTINATION bin)
