project(mris_svm_classify)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_svm_classify mris_svm_classify.cpp)
target_link_libraries(mris_svm_classify utils svm)

install(TARGETS mris_svm_classify DESTINATION bin)
