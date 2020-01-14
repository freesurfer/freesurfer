project(AntsN4BiasFieldCorrectionFs)

include_directories(${FS_INCLUDE_DIRS})

add_compile_options(-Wno-self-assign-field)

add_executable(AntsN4BiasFieldCorrectionFs AntsN4BiasFieldCorrectionFs.cpp)
add_help(AntsN4BiasFieldCorrectionFs AntsN4BiasFieldCorrectionFs.help.xml)
target_link_libraries(AntsN4BiasFieldCorrectionFs utils)

install(TARGETS AntsN4BiasFieldCorrectionFs DESTINATION bin)
