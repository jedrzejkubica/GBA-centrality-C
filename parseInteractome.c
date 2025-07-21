#include </usr/include/python3.9/Python.h>
#include </usr/include/python3.9/pyconfig-64.h>

#include "adjacency.h"

/*
    this module calls python scripts for parsing the interactome
*/

void *parseInteractome(char *interactomeFile) {
    printf("parsing %s\n", interactomeFile);

    Py_Initialize();

    PyRun_SimpleString("import sys; sys.path.append(\"Interactome/\")");

    PyObject *pName = PyUnicode_FromString("interactome");
    PyObject *pModule = PyImport_Import(pName);
    Py_DECREF(pName);  // free allocated memory

    if (pModule != NULL) {
        PyObject *pFunction = PyObject_GetAttrString(pModule, "parse_interactome");

        if (pFunction && PyCallable_Check(pFunction)) {
            /*
                parse_interactome() expects a filename (str),
                pArgs must be a touple
            */
            PyObject *pArgs = PyTuple_Pack(1, PyUnicode_FromString(interactomeFile));

            PyObject *pOutput = PyObject_CallObject(pFunction, pArgs);
            Py_DECREF(pArgs);

            if (pOutput != NULL && PyTuple_Check(pOutput)) {
                PyObject *pNodes = PyTuple_GetItem(pOutput, 0);
                PyObject *pMatrix = PyTuple_GetItem(pOutput, 1);
                Py_DECREF(pMatrix);

                Py_ssize_t num_nodes = (int) PyLong_AsLong(pNodes);
                printf("Number of nodes: %zd\n", num_nodes);

                Py_DECREF(pOutput);
            } else {
                PyErr_Print();
            }
            Py_XDECREF(pFunction);
        } else {
            if (PyErr_Occurred()) PyErr_Print();
        }
        Py_DECREF(pModule);
    } else {
        PyErr_Print();
    }
    Py_Finalize();

    return 0;
}