project(mris_fbirn_annot)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_fbirn_annot mris_fbirn_annot.c)
target_link_libraries(mris_fbirn_annot utils)

install(TARGETS mris_fbirn_annot DESTINATION bin)
