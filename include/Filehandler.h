/* Filehandler.h
 * author: Christoph Froehner
 *
 * This class is responsible for writing data to files, e.g.
 * trajectories to xyz files.
 */

#include <iostream>
#include <fstream>
#include <vector>
#include <array>

class Filehandler
{
	public:
		void writeSingleParticleTrajectoryXyz(std::vector< std::array<double, 3> > trajectory);

};
