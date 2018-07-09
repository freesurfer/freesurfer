project(mri_evaluate_morph)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mri_evaluate_morph mri_evaluate_morph.c)
target_link_libraries(mri_evaluate_morph utils)

install(TARGETS mri_evaluate_morph DESTINATION bin)
