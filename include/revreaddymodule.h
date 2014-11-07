double *pyvector_to_Carrayptrs(PyArrayObject *arrayin);
void printDoubleArray(double *arr, int length);
void printIntArray(int *arr, int length);
gsl_rng * initializeRandomNumberGenerator();
unsigned long int random_seed();
void generateSeq(double* seq, int length, gsl_rng * randGen);
