/* utils.h
 * author: Christoph Froehner
 *
 * several functions needed everywhere.
 */

#ifndef __UTILS_H_INCLUDED__
#define __UTILS_H_INCLUDED__

inline void printArray(double *arr, int length)
{
	printf("[ ");
	for (int i = 0; i < length; i++)
	{
		printf("%f ", arr[i]);
	}
	printf("]\n");
}

inline double squaredDistance(double * arr1, double * arr2)
{
	double result;
	result = ( arr1[0] - arr2[0] ) * ( arr1[0] - arr2[0] );
	result += ( arr1[1] - arr2[1] ) * ( arr1[1] - arr2[1] );
	result += ( arr1[2] - arr2[2] ) * ( arr1[2] - arr2[2] );
	return result;
}

#endif
