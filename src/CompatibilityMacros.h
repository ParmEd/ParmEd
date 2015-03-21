#if PY_MAJOR_VERSION >= 3

// All of the Python3-specific macros go here

int PyObject_IS_STRING(PyObject *chk) {
    return PyUnicode_Check(chk) && PyUnicode_READY(chk) >= 0;
}


// String handling is entirely unicode now
#   define PyString_Size PyUnicode_GET_LENGTH
#   define PyString_AsString PyUnicode_AsUTF8
#   define PyString_FromString PyUnicode_FromString

#   define PY_DESTROY_TYPE Py_TYPE(self)->tp_free((PyObject *)self)

// PyInt -> PyLong
#   define PyInt_FromLong PyLong_FromLong
#   define PyInt_AsLong PyLong_AsLong

#else /* Python 2 */

#   define PyObject_IS_STRING PyString_Check
#   define PY_DESTROY_TYPE self->ob_type->tp_free((PyObject *)self)

#endif
