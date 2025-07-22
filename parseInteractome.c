#include </usr/include/python3.9/Python.h>
#include </usr/include/python3.9/pyconfig-64.h>

#include "adjacency.h"
#include "mem.h"

/*
    this module calls python scripts for parsing the interactome
*/


/*
    the python script should include:
    args:
    - interactome_file: str, path to SIF file (3 tab-separated columns: ENSG1 pp ENSG2)

    returns:
    - nbNodes: int
    - weights: list of floats, flattened adjacency matrix,
      weights[i*nbCols + j] contains the weight from node i to node j,
      note: weights must be 0/1 for unweighted, or in [0, 1] for weighted
*/
adjacencyMatrix *parseInteractome(char *interactomeFile) {
    printf("parsing %s\n", interactomeFile);

    adjacencyMatrix *interactome = mallocOrDie(sizeof(adjacencyMatrix), "E: OOM for interactome");

    Py_Initialize();

    PyRun_SimpleString("import sys; sys.path.append(\".\")");

    PyObject *pName = PyUnicode_FromString("parseInteractome");
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
                PyObject *pWeights = PyTuple_GetItem(pOutput, 1);

                unsigned int nbNodes = (unsigned int) PyLong_AsLong(pNodes);
                interactome->nbCols = nbNodes;

                interactome->weights = mallocOrDie(nbNodes * nbNodes * sizeof(float), "E: OOM for interactome weights");

                for (size_t i = 0; i < nbNodes * nbNodes; i++) {
                    PyObject *w = PyList_GetItem(pWeights, i);

                    if (!PyFloat_Check(w)) {
                        fprintf(stderr, "List item %zd is not a float\n", i);
                    }

                    double tmp = PyFloat_AsDouble(w);
                    float val = (float)tmp;
                    interactome->weights[i] = val;
                }
                
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

    return interactome;
}