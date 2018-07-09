project(mkxsubjreg)

include_directories(${CMAKE_SOURCE_DIR}/include)

add_executable(mkxsubjreg mkxsubjreg.c)
target_link_libraries(mkxsubjreg utils)

install(TARGETS mkxsubjreg DESTINATION bin)
