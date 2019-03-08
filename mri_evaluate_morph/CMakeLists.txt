project(mri_evaluate_morph)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mri_evaluate_morph mri_evaluate_morph.cpp)
target_link_libraries(mri_evaluate_morph utils)

install(TARGETS mri_evaluate_morph DESTINATION bin)
