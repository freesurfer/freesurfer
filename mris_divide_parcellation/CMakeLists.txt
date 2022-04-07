project(mris_divide_parcellation)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_divide_parcellation mris_divide_parcellation.cpp)
target_link_libraries(mris_divide_parcellation utils)

# temporarily remove until fixed
# add_test_script(NAME mris_divide_parcellation_test SCRIPT test.sh)

install(TARGETS mris_divide_parcellation DESTINATION bin)
