# This file defines a few custom cmake utility functions to
# simplify the freesurfer build process


# add_subdirectories(<subdirs>)
# Simple utility to add multiple subdirectories at once
function(add_subdirectories)
  foreach(SUBDIR ${ARGN})
    add_subdirectory(${SUBDIR})
  endforeach()
endfunction()


# install_symlinks(<files> DESTINATION <dir> TYPE <type>)
# Unfortunately, cmake does not follow symlinks when installing files,
# so this is a wrapper for the install() function that will actually follow
# links. The TYPE argument should specify the install type that would normally
# be provided to the install function (i.e. PROGRAMS, FILES, etc.)
function(install_symlinks)
  cmake_parse_arguments(INSTALL "" "DESTINATION;TYPE" "" ${ARGN})
  string(TOUPPER ${INSTALL_TYPE} INSTALL_TYPE)
  foreach(arg ${INSTALL_UNPARSED_ARGUMENTS})
    get_filename_component(BASENAME ${arg} NAME)
    get_filename_component(ABS_PATH ${arg} REALPATH)
    install(${INSTALL_TYPE} ${ABS_PATH} DESTINATION ${INSTALL_DESTINATION} RENAME ${BASENAME})
  endforeach()
endfunction()


# install_tarball(<tarball> DESTINATION <dir>)
# Use this to extract a tarball to the provided install destination
function(install_tarball)
  cmake_parse_arguments(INSTALL "" "DESTINATION" "" ${ARGN})
  get_filename_component(TARBALL ${INSTALL_UNPARSED_ARGUMENTS} ABSOLUTE)
  install(CODE "
    message(STATUS \"Extracting: ${TARBALL} to ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION}\")
    file(MAKE_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION})
    execute_process(
      COMMAND bash -c \"tar -xzf ${TARBALL}\"
      WORKING_DIRECTORY ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION}
      RESULT_VARIABLE retcode
    )
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not extract ${INSTALL_UNPARSED_ARGUMENTS} - perhaps it hasn't been downloaded from the annex?\")
    endif()"
  )
endfunction()


# install_wrapped(<binary> DESTINATION <dir>)
# Some freesurfer binaries require a wrapper script before being executed.
# This function will install the program as binary.bin in the provied
# install destination and create the wrapper in the bin dir. By default,
# the script will contain:
#
# #!/bin/bash
# source $FREESURFER_HOME/sources.sh
# binary.bin "$@"
# 
# but this can be replaced with a custom string using the WRAPCODE argument.
# The .bin extension can also be replaced with the EXT argument
function(install_wrapped)
  cmake_parse_arguments(INSTALL "" "TARGETS;PROGRAMS;DESTINATION;WRAPCODE;EXT" "" ${ARGN})
  if(INSTALL_TARGETS)
    set(BINARY "${INSTALL_TARGETS}")
  elseif(INSTALL_PROGRAMS)
    set(BINARY "${INSTALL_PROGRAMS}")
  endif()
  set(WRAPCODE "#!/bin/bash\nsource $FREESURFER_HOME/sources.sh\n${BINARY}.bin \\\"$@\\\"")
  if(INSTALL_WRAPCODE)
    set(WRAPCODE ${INSTALL_WRAPCODE})
  endif()
  set(EXT "bin")
  if(INSTALL_EXT)
    set(EXT ${INSTALL_EXT})
  endif()
  if(INSTALL_TARGETS)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${BINARY} DESTINATION ${INSTALL_DESTINATION} RENAME ${BINARY}.${EXT})
  elseif(INSTALL_PROGRAMS)
    install(PROGRAMS ${BINARY} DESTINATION ${INSTALL_DESTINATION} RENAME ${BINARY}.${EXT})
  endif()
  install(CODE "
    message(STATUS \"Wrapping: ${CMAKE_INSTALL_PREFIX}/${INSTALL_DESTINATION}/${BINARY}.${EXT}\")
    execute_process(
      COMMAND bash -c \"echo -e '${WRAPCODE}' > ${CMAKE_INSTALL_PREFIX}/bin/${BINARY} &&
                        chmod +x ${CMAKE_INSTALL_PREFIX}/bin/${BINARY}\"
      RESULT_VARIABLE retcode
    )
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not wrap ${BINARY}\")
    endif()"
  )
endfunction()


# add_help(<binary> <xml>)
# Link an xml helptext to a target binary. This will create a target dependency on
# the help file and will run xxd to create the xml header during the build
function(add_help BINARY HELPTEXT)
  add_custom_command(COMMAND cd ${CMAKE_CURRENT_SOURCE_DIR} &&
                             xxd -i ${HELPTEXT} ${CMAKE_CURRENT_BINARY_DIR}/${HELPTEXT}.h
                     OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${HELPTEXT}.h
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${HELPTEXT})
  include_directories(${CMAKE_CURRENT_BINARY_DIR})
  target_sources(${BINARY} PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/${HELPTEXT}.h)
  install(FILES ${HELPTEXT} DESTINATION docs/xml)
endfunction(add_help)


# install_append_help(<script> <xml> <destination>)
# Some shell scripts also utilize a help.xml file by appending the helptext to the end of it.
# This function will setup this build command and dependency
function(install_append_help SCRIPT HELPTEXT DESTINATION)
  install(CODE "
    message(STATUS \"Installing (with helptext): ${CMAKE_INSTALL_PREFIX}/${DESTINATION}/${SCRIPT}\")
    execute_process(
      COMMAND bash -c \"
        ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR}/utils --target fsPrintHelp &&
        cp -f ${CMAKE_CURRENT_SOURCE_DIR}/${SCRIPT} ${CMAKE_INSTALL_PREFIX}/${DESTINATION} &&
        ${CMAKE_BINARY_DIR}/utils/fsPrintHelp ${CMAKE_CURRENT_SOURCE_DIR}/${HELPTEXT} >> ${CMAKE_INSTALL_PREFIX}/${DESTINATION}/${SCRIPT}\"
      OUTPUT_QUIET
      RESULT_VARIABLE retcode
    )
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not append help text to ${SCRIPT}\")
    endif()"
  )
  install(FILES ${HELPTEXT} DESTINATION docs/xml)
endfunction()


# symlink(<target> <linkname>)
# Creates a symlink at install time
function(symlink TARGET LINKNAME)
  install(CODE "
    message(STATUS \"Symlinking: ${LINKNAME} to ${TARGET}\")
    execute_process(COMMAND bash -c \"mkdir -p $(dirname ${LINKNAME}) && rm -f ${LINKNAME} && ln -s ${TARGET} ${LINKNAME}\" RESULT_VARIABLE retcode)
    if(NOT \${retcode} STREQUAL 0)
      message(FATAL_ERROR \"Could not create symlink to ${TARGET}\")
    endif()"
  )
endfunction()


# add_test_script(NAME <testname> SCRIPT <testscript> DEPENDS <targets>)
# Adds a script to the test framework. If the script calls a freesurfer binary, it
# should be specified with DEPENDS to guarantee it gets built beforehand
function(add_test_script)
  cmake_parse_arguments(TEST "" "NAME;SCRIPT" "DEPENDS" ${ARGN})
  foreach(TARGET ${TEST_DEPENDS})
    set(TEST_CMD "${TEST_CMD} ${CMAKE_COMMAND} --build ${CMAKE_CURRENT_BINARY_DIR} --target ${TARGET} &&")
  endforeach()
  add_test(${TEST_NAME} bash -c "${TEST_CMD} ${CMAKE_CURRENT_SOURCE_DIR}/${TEST_SCRIPT}")
endfunction()
