/* main.cpp
 *
 * This is a pure C++ test environment for
 * the revreaddy python extension.
 *
 */

#include <Random.h>
#include <Particle.h>
#include <ActiveParticles.h>

void printArray(double *arr, int length)
{
	printf("[ ");
	for (int i = 0; i < length; i++)
	{
		printf("%f ", arr[i]);
	}
	printf("]\n");
}

int main()
{
	Particle * p1 = new Particle();
	Particle * p2 = new Particle();
	Random * random = new Random("ranlxs0");
	double arr1[3];
	double arr2[3];
	random->normal3D(arr1);
	random->normal3D(arr2);

	ActiveParticles * ap = new ActiveParticles();
	ap->addParticle(p1);
	ap->addParticle(p2);
	printArray(ap->container[0]->position, 3);
	printArray(ap->container[1]->position, 3);
	printf("0");
	ap->container.clear();
	printf("1");
	delete ap;
	printf("2");
	delete p1;
	printf("3");
	delete p2;
	printf("4");
	delete random;
	printf("5");
	return 0;
}
