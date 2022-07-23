project(mri_em_register)

include_directories(${FS_INCLUDE_DIRS})

set(SOURCES
  mri_em_register.cpp
  findtranslation.cpp
  emregisterutils.cpp
)

add_test_script(NAME mri_em_register_test SCRIPT test.sh)

add_executable(mri_em_register ${SOURCES})
add_help(mri_em_register mri_em_register.help.xml)
target_link_libraries(mri_em_register utils)
install(TARGETS mri_em_register DESTINATION bin)
