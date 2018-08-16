// This file is the boilerplate that allows the optimized prmtop reader to be
// accessed from Python

#include <Python.h>

// Support versions of Python older than 2.5 that didn't define Py_ssize_t
#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#   define PY_SSIZE_T_MAX INT_MAX
#   define PY_SSIZE_T_MIN INT_MIN
#endif

// A set of macros for use with Py3
#include "CompatibilityMacros.h"

// Optimized readparm
#include "readparm.h"

void SetItem_PyDict_AndDecref(PyObject *dict, const char* key, PyObject *value) {
    PyDict_SetItemString(dict, key, value);
    Py_DECREF(value);
}

static PyObject* rdparm(PyObject *self, PyObject *args) {

    char *filename;

    if (!PyArg_ParseTuple(args, "s", &filename))
        return NULL;

    std::string fname(filename);
    std::string error_message;

    // The data we are parsing from the prmtop

    ParmDataMap parmData;
    ParmStringMap parmComments, unkParmData;
    ParmFormatMap parmFormats;
    std::vector<std::string> flagList;
    std::string version;
    ExitStatus retval;

    Py_BEGIN_ALLOW_THREADS

    retval = readparm(fname, flagList, parmData, parmComments,
                      unkParmData, parmFormats, version);

    Py_END_ALLOW_THREADS

    if (retval == NOOPEN) {
        error_message = "Could not open " + fname + " for reading";
        PyErr_SetString(PyExc_IOError, error_message.c_str());
        return NULL;
    }

    if (retval == NOVERSION) {
        error_message = "Could not find %VERSION in " + fname;
        PyErr_SetString(PyExc_TypeError, error_message.c_str());
        return NULL;
    }

    if (retval == EMPTY) {
        error_message = fname + " was empty";
        PyErr_SetString(PyExc_ValueError, error_message.c_str());
        return NULL;
    }

    if (retval == ERR) {
        error_message = "Prmtop parsing error parsing " + fname;
        PyErr_SetString(PyExc_RuntimeError, error_message.c_str());
        return NULL;
    }

    // If we got here, the parsing must have been OK. Create the parm_data,
    // formats, and comments dicts to pass back to Python
    PyObject *parm_data = PyDict_New();
    PyObject *comments = PyDict_New();
    PyObject *formats = PyDict_New();

    PyObject *unknown_flags = PyList_New((Py_ssize_t) unkParmData.size());
    PyObject *flag_list = PyList_New((Py_ssize_t) flagList.size());

    Py_ssize_t unkFlagNum = 0;

    for (size_t i = 0; i < flagList.size(); i++) {
        std::string flag = flagList[i];
        Py_ssize_t listSize;
        PyObject *list;
        if (parmFormats[flag].dataType == UNKNOWN)
            listSize = (Py_ssize_t) unkParmData[flag].size();
        else
            listSize = (Py_ssize_t) parmData[flag].size();
        list = PyList_New(listSize);
        // Now see what type this is and fill the list up accordingly
        switch (parmFormats[flag].dataType) {
            case INTEGER:
                for (Py_ssize_t j = 0; j < listSize; j++) {
                    long val = (long) parmData[flag][(size_t)j].i;
                    PyList_SET_ITEM(list, j, PyInt_FromLong(val));
                }
                break;
            case FLOAT:
                for (Py_ssize_t j = 0; j < listSize; j++) {
                    double val = parmData[flag][(size_t)j].f;
                    PyList_SET_ITEM(list, j, PyFloat_FromDouble(val));
                }
                break;
            case HOLLERITH:
                for (Py_ssize_t j = 0; j < listSize; j++) {
                    PyList_SET_ITEM(list, j,
                            PyString_FromString(parmData[flag][(size_t)j].c));
                }
                break;
            case UNKNOWN:
                for (Py_ssize_t j = 0; j < listSize; j++) {
                    std::string line = unkParmData[flag][(size_t) j];
                    PyList_SET_ITEM(list, j, PyString_FromString(line.c_str()));
                }
                // Add this to the list of unknown flags
                PyList_SET_ITEM(unknown_flags, unkFlagNum++,
                                PyString_FromString(flag.c_str()));
                break;
            default:
                // Should not be here
                PyErr_SetString(PyExc_RuntimeError, "This should be unreachable");
                return NULL;
        }
        SetItem_PyDict_AndDecref(parm_data, flag.c_str(), list);

        // Now comments
        if (parmComments.count(flag) == 0) {
            SetItem_PyDict_AndDecref(comments, flag.c_str(), PyList_New(0));
        } else {
            int ncom = parmComments[flag].size();
            PyObject *commentList = PyList_New(ncom);
            for (Py_ssize_t j = 0; j < ncom; j++) {
                std::string line = parmComments[flag][(size_t)j];
                PyList_SET_ITEM(commentList, j,
                                PyString_FromString(line.c_str()));
            }
            SetItem_PyDict_AndDecref(comments, flag.c_str(), commentList);
        }

        // Now formats
        PyObject *fmt = PyString_FromString(parmFormats[flag].fmt.c_str());
        SetItem_PyDict_AndDecref(formats, flag.c_str(), fmt);

        // Now flag list
        PyList_SET_ITEM(flag_list, (Py_ssize_t)i,
                        PyString_FromString(flag.c_str()));
    }

    PyObject *ret = PyTuple_New(6);
    PyTuple_SET_ITEM(ret, 0, parm_data);
    PyTuple_SET_ITEM(ret, 1, comments);
    PyTuple_SET_ITEM(ret, 2, formats);
    PyTuple_SET_ITEM(ret, 3, unknown_flags);
    PyTuple_SET_ITEM(ret, 4, flag_list);
    PyTuple_SET_ITEM(ret, 5, PyString_FromString(version.c_str()));

    return ret;
}

static PyMethodDef
optrdparmMethods[] = {
    { "rdparm", (PyCFunction) rdparm, METH_VARARGS,
            "Optimized prmtop file reading library written in C++"},
    { NULL },
};

#if PY_MAJOR_VERSION >= 3
static struct PyModuleDef moduledef = {
    PyModuleDef_HEAD_INIT,
    "_rdparm",                          // m_name
    "Optimized prmtop reading routine", // m_doc
    -1,                                 // m_size
    optrdparmMethods,                   // m_methods
    NULL,
    NULL,
    NULL,
    NULL,
};
#endif

#if PY_MAJOR_VERSION >= 3
PyMODINIT_FUNC
PyInit__rdparm(void) {
#else
PyMODINIT_FUNC
init_rdparm(void) {
#endif
    PyObject *m;

#if PY_MAJOR_VERSION >= 3
    m = PyModule_Create(&moduledef);
    return m;
#else
    m = Py_InitModule3("_rdparm", optrdparmMethods,
            "Optimized prmtop file reading library written in C++");
			
	(void)m; // silence unused variable warning
#endif
}
