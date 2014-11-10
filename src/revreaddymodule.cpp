/* revreaddymodule.cpp
 * author: Christoph Froehner
 *
 * This is the Python C-Extension 'revreaddy':
 * reversible reaction diffusion dynamics simulation
 *
 */

#include <Python.h>
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <stdio.h>
#include "revreaddymodule.h"
#include "Random.h"
#include "Particle.h"

using namespace std;

static PyObject* revreaddy_start(PyObject * self, PyObject * args)
{
	int lengthOfSeq;
	if (!PyArg_ParseTuple(args,"i", &lengthOfSeq) )
	{
		return NULL;
	}
	
	int nd = 1;
	npy_intp dim[1] = {(npy_intp) lengthOfSeq};
	
	double * cSeq;
	PyArrayObject * pySeq;
	
	pySeq = (PyArrayObject*) PyArray_SimpleNew(nd, dim, NPY_DOUBLE);
	cSeq = pyVectorToCArrayPointer(pySeq);
	
	// ----------- Do stuff here ------------ //

	Random random("ranlxs0");
	double vectorA[3];
	Particle h2o;


	
	// --------------------------------------//

    return PyArray_Return(pySeq);
}

double * pyVectorToCArrayPointer(PyArrayObject *arrayin)
{
	return (double *) arrayin->data;  // pointer to arrayin data as double
}

static char revreaddy_docs[] = "start( ): Will start the simulation.\n";

static PyMethodDef revreaddy_funcs[] =
{
    {"start", (PyCFunction)revreaddy_start, METH_VARARGS, revreaddy_docs},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initrevreaddy()
{
    PyObject *m = Py_InitModule3("revreaddy", revreaddy_funcs, "reversible reaction diffusion dynamics simulation");
    if (m == NULL)
    	return;
    
    /* Load numpy functionality */
    import_array();
}
