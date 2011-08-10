// Licensed under MIT license; see license.txt.

#ifdef USE_PYTHON
#include <sbl/other/Scripting.h>
#include <sbl/core/Command.h>
#include <sbl/core/PathConfig.h>
#include <python.h>
namespace sbl {


//-------------------------------------------
// SBLC PYTHON MODULE DEFINITION
//-------------------------------------------


// wrapper for disp(), warning(), and fatalError()
static PyObject *sblc_disp( PyObject *self, PyObject *args ) {
    int type;
    int indent;
    const char *message;
    if (!PyArg_ParseTuple( args, "iis", &type, &indent, &message ))
        return NULL; // fix(later): better error handling
    if (type == 0)
        disp( indent, message );
    else if (type == 1)
        warning( message );
    else if (type == 2)
        fatalError( message );
    int ret = 0;
    return Py_BuildValue( "i", ret );
}


// wrapper for execCommand()
static PyObject *sblc_execCommand( PyObject *self, PyObject *args ) {
    const char *command;
    if (!PyArg_ParseTuple( args, "s", &command ))
        return NULL; // fix(later): better error handling
    execCommand( command, false );
    int ret = 0;
    return Py_BuildValue( "i", ret );
}


// wrapper for checkCommandEvents()
static PyObject *sblc_checkCommandEvents( PyObject *self, PyObject *args ) {
    int ret = checkCommandEvents() ? 1 : 0;
    return Py_BuildValue( "i", ret );
}


// wrapper for dataPath()
static PyObject *sblc_dataPath( PyObject *self, PyObject *args ) {
    return Py_BuildValue( "s", dataPath().c_str() );
}


// our python module methods
static PyMethodDef sblcMethods[] = {
    { "disp", sblc_disp, METH_VARARGS, "Display a message." },
    { "execCommand", sblc_execCommand, METH_VARARGS, "Execute an SBL command." },
    { "checkCommandEvents", sblc_checkCommandEvents, METH_NOARGS, "Returns non-zero if user has cancelled command." },
    { "dataPath", sblc_dataPath, METH_NOARGS, "Returns the main config's dataPath." },
    { NULL, NULL, 0, NULL }
};


// define our python module
struct PyModuleDef sblcModuleDef = {
    PyModuleDef_HEAD_INIT,
    "sblc",
    NULL,
    -1,
    sblcMethods,
    NULL, NULL, NULL, NULL
};


// factory for our python module
PyMODINIT_FUNC sblcInitModule() {
    return PyModule_Create( &sblcModuleDef );
}


//-------------------------------------------
// OUR INTERFACE TO PYTHON
//-------------------------------------------


// terminate python embedding
void cleanUpPython() {
    Py_Finalize();
}


// initialize python embedding
void initPython() {
    static bool s_initDone = false;
    if (s_initDone == false) {
        PyImport_AppendInittab( "sblc", sblcInitModule );
        Py_Initialize();
        registerCleanUp( cleanUpPython );
        s_initDone = true;
    }
}


// run a python script (embedded within the current program)
void runPythonScript( Config &conf ) {

    // get command parameters
    conf.enableDualPass();
    String fileName = addDataPath( conf.readString( "scriptName" ) );
    if (conf.initialPass())
        return;

    // load and run the script
    initPython();
    FILE *file = fopen( fileName.c_str(), "r" );
    if (file) {
        PyRun_SimpleFile( file, fileName.c_str() );
        fclose( file );
    }
}


// a command for testing python interaction
void testCommand( Config &conf ) {

    // get command parameters
    conf.enableDualPass();
    int intParam = conf.readInt( "intParam" );
    String stringParam = conf.readString( "stringParam" );
    double doubleParam = conf.readDouble( "doubleParam" );
    bool boolParam = conf.readBool( "boolParam" );
    if (conf.initialPass())
        return;

    // display parameters
    disp( 1, "intParam: %d", intParam );
    disp( 1, "stringParam: %s", stringParam );
    disp( 1, "doubleParam: %f", doubleParam );
    disp( 1, "boolParam: %d", boolParam );
}


//-------------------------------------------
// INIT / CLEAN-UP
//-------------------------------------------


// register commands, etc. defined in this module
void initScripting() {
    registerCommand( "runpy", runPythonScript );
#ifdef REGISTER_TEST_COMMANDS
    registerCommand( "testcmd", testCommand );
#endif
}


} // end namespace sbl
#endif // USE_PYTHON
