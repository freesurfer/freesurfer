project(freesurfer)
cmake_minimum_required(VERSION 3.2)

include_directories(${CMAKE_SOURCE_DIR}/include
                    ${CMAKE_SOURCE_DIR}/include/dicom
                    ${CMAKE_SOURCE_DIR}/xml2
                    /usr/pubsw/packages/mni/1.5/include)

# create utils library
add_subdirectory(utils)