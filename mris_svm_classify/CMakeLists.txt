project(mris_svm_classify)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mris_svm_classify mris_svm_classify.c)
target_link_libraries(mris_svm_classify utils svm)

install(TARGETS mris_svm_classify DESTINATION bin)
