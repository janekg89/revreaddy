#include <Python.h>
//#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/arrayobject.h>
#include <gsl/gsl_rng.h> 	// GNU random number generators
#include <gsl/gsl_randist.h>// GNU random number distributions
#include <gsl/gsl_math.h>	// GNU math libraries
#include <stdio.h>
#include "revreaddymodule.h"

static PyObject* revreaddy_start(PyObject * self, PyObject * args)
{
	int lengthOfSeq;
	/* save args to according memory */
	if (!PyArg_ParseTuple(args,"i", &lengthOfSeq) )
	{
		return NULL;
	}
	
	int nd = 1;
	npy_intp dim[1] = {(npy_intp) lengthOfSeq};
	
	/* cSeq will hold array in c style and pySeq is wrapper around it. */
	double *cSeq ;
	PyArrayObject * pySeq;
	
	/* Let Python allocate memory for the array so it won't be deallocated when returning to python. */
	pySeq = (PyArrayObject*) PyArray_SimpleNew(nd, dim, NPY_DOUBLE);
	
	/* cSeq will be the pointer to the cArray data. */
	cSeq = pyvector_to_Carrayptrs(pySeq);
	
	/* generateSeq will fill the array with random numbers. */
	gsl_rng * randGen = initializeRandomNumberGenerator();
	generateSeq(cSeq, lengthOfSeq, randGen);
	
	/* deallocate memory for random number generator */
	gsl_rng_free(randGen);
	
    return PyArray_Return(pySeq);
}

double *pyvector_to_Carrayptrs(PyArrayObject *arrayin)
{
	// int n=arrayin->dimensions[0];
	return (double *) arrayin->data;  // pointer to arrayin data as double
}

void printDoubleArray(double *arr, int length)
{
	printf("[ ");
	for (int i = 0; i < length; i++) {
		printf("%f ",arr[i]);
	}
	printf("]\n");
}

void printIntArray(int *arr, int length)
{
	printf("[ ");
	for (int i = 0; i < length; i++) {
		printf("%i ",arr[i]);
	}
	printf("]\n");
}

gsl_rng * initializeRandomNumberGenerator()
{
	const gsl_rng_type * randGenType;
	gsl_rng * randGen;
	unsigned long int seed;

	gsl_rng_env_setup();
	
	/* choose type of random number generator e.g. ranlxs0, mt19937 */
	randGenType = gsl_rng_mt19937;
	
	/* randGen is the handle to generate random numbers from */
	randGen = gsl_rng_alloc(randGenType);

	seed = random_seed();
	gsl_rng_set(randGen, seed);
	printf("Random number generator of type '%s' has been intialized and seeded.\n", gsl_rng_name(randGen));
	return randGen;
}

unsigned long int random_seed()
{
	unsigned long int seed;
	FILE *devurandom;
	devurandom = fopen("/dev/urandom","r");
	fread(&seed, sizeof(seed), 1, devurandom);
	fclose(devurandom);
	printf("Got seed %lu from /dev/urandom.\n",seed);
	return(seed);
}

void generateSeq(double* seq, int length, gsl_rng * randGen)
{
	for (int i = 0; i < length; i++) {
		//seq[i] = gsl_rng_uniform(randGen);
		seq[i] = gsl_ran_gaussian_ziggurat (randGen, 1.0);
	}
}

static char revreaddy_docs[] = "start( ): Will start the simulation.\n";

static PyMethodDef revreaddy_funcs[] =
{
    {"start", (PyCFunction)revreaddy_start, METH_VARARGS, revreaddy_docs},
    {NULL, NULL, 0, NULL}
};

PyMODINIT_FUNC initrevreaddy()
{
    PyObject *m = Py_InitModule3("revreaddy", revreaddy_funcs, "interacting particle reaction diffusion simulation");
    if (m == NULL)
    	return;
    
    /* Load numpy functionality */
    import_array();
}
