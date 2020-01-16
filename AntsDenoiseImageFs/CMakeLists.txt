project(AntsDenoiseImageFs)

include_directories(${FS_INCLUDE_DIRS})

add_compile_options(-Wno-self-assign-field)

add_executable(AntsDenoiseImageFs AntsDenoiseImageFs.cpp)
add_help(AntsDenoiseImageFs AntsDenoiseImageFs.help.xml)
target_link_libraries(AntsDenoiseImageFs utils)

install(TARGETS AntsDenoiseImageFs DESTINATION bin)
