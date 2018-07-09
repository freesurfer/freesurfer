project(mri_em_register)

include_directories(${CMAKE_SOURCE_DIR}/include)

set(SOURCES
  mri_em_register.c
  findtranslation.cpp
  emregisterutils.cpp
)

add_executable(mri_em_register ${SOURCES})
add_help(mri_em_register mri_em_register.help.xml)
target_link_libraries(mri_em_register utils)

install(TARGETS mri_em_register DESTINATION bin)
