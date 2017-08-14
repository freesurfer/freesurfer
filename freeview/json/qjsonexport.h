
#ifndef QJSONEXPORT_H
#define QJSONEXPORT_H


#ifdef Q_JSONRPC_DLL
#  ifdef Q_BUILD_JSONRPC // build qjsonrpc Dll
#    define Q_JSONRPC_EXPORT Q_DECL_EXPORT
#   else // use qjsonrpc as Dll
#     define Q_JSONRPC_EXPORT Q_DECL_IMPORT
#   endif
#else
#  define Q_JSONRPC_EXPORT
#endif

#endif
