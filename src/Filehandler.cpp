/* Filehandler.cpp
 * author: Christoph Froehner
 */

#include "Filehandler.h"

void Filehandler::writeSingleParticleTrajectoryXyz(std::vector< std::array<double, 3> > trajectory)
{
	std::ofstream file;
	file.open ("trajectory.xyz");
	file << "1\nsingle particle\n";
	for (auto&& position : trajectory)
	{
		file << "A\t";
		file << position[0] << "\t";
		file << position[1] << "\t";
		file << position[2] << "\n";
	}
	file.close();
}
