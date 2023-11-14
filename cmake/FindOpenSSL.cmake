# OpenSSL Find Module

if(APPLE)
   set(SSL_SHARED_LIB libcrypto.a)
   # newer homebrew installations, e.g., darwin_arm64
   set(SSL_INC_1 /opt/homebrew/opt/openssl/include)
   set(SSL_LIB_1 /opt/homebrew/opt/openssl/lib)
   # older homebrew installations, e.g., darwin_x86_64
   set(SSL_INC_2 /usr/local/opt/openssl/include)
   set(SSL_LIB_2 /usr/local/opt/openssl/lib)
else()
   set(SSL_SHARED_LIB libcrypto.so)
   # Ubuntu Linux
   set(SSL_LIB_1 /usr/lib/x86_64-linux-gnu)
   set(SSL_INC_1 /usr/include)
   # Redhat/CentOS
   set(SSL_LIB_2 /lib64)
   set(SSL_INC_2 /usr/include)
endif()

# find the include dir - search 1st for newer homebrew installation (which only exists on darwin_arm64 machine)
find_path(OpenSSL_INCLUDE_DIR HINTS ${SSL_INC_1} ${SSL_INC_2} NAMES openssl/conf.h PATH_SUFFIXES openssl)
#set(OpenSSL_INCLUDE_DIR "${OpenSSL_INCLUDE_DIR}/..")

# find the lib dir - search 1st for Ubuntu path which will not be found on CentOS/RedHat (as the reverse is not true)
find_path(OpenSSL_LIB_DIR HINTS ${SSL_LIB_1} ${SSL_LIB_2} NAMES ${SSL_SHARED_LIB})

message(STATUS "OpenSSL_INCLUDE_DIR=${OpenSSL_INCLUDE_DIR}")
message(STATUS "OpenSSL_LIB_DIR=${OpenSSL_LIB_DIR}")

