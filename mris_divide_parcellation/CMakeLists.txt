project(mris_divide_parcellation)

include_directories(${FS_INCLUDE_DIRS})

add_executable(mris_divide_parcellation mris_divide_parcellation.c)
target_link_libraries(mris_divide_parcellation utils)

add_test_script(NAME mris_divide_parcellation_test SCRIPT test.sh DEPENDS mris_divide_parcellation)

install(TARGETS mris_divide_parcellation DESTINATION bin)
