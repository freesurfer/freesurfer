project(mri_parselabel)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_parselabel mri_parselabel.cpp)
target_link_libraries(mri_parselabel utils)

install(TARGETS mri_parselabel DESTINATION bin)
