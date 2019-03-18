project(mris_fbirn_annot)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_fbirn_annot mris_fbirn_annot.cpp)
target_link_libraries(mris_fbirn_annot utils)

install(TARGETS mris_fbirn_annot DESTINATION bin)
