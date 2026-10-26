project(samseg)

if (FSPYTHON_INSTALL_TREE)
   set(pydir "${FSPYTHON_INSTALL_PREFIX}/python")
else()
   set(pydir "${CMAKE_INSTALL_PREFIX}/python")
endif()

install_configured(samseg samseg-long DESTINATION bin)

## should work on intel and silicon Mac (12/29/2025)
# if(NOT APPLE_ARM64)
## turn off for darwin_x86_64 (05/14/2026)
if(NOT APPLE)
   add_test_script(NAME samseg_test SCRIPT test.sh)
elseif(APPLE_ARM64)
   add_test_script(NAME samseg_test SCRIPT test.sh)
endif()
# endif()

