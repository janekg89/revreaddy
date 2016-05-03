/* Increments.h
 * Calculate and save the squared displacements of a specific particle type. Track via uniqueIds.
 * The output is a TxN matrix that contains the norm of the displacements. */

#ifndef REVREADDY_INCREMENTS_H
#define REVREADDY_INCREMENTS_H
#include "Observable.h"
#include "../Config.h"
#include "../World.h"
#include <string>
#include <vector>
#include <map>

class Increments : public Observable {
public:
    Increments(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename, unsigned int pTypeId);
    void setup(World* world, Config* config);
    void configure(World* world, Config* config);
    void record(World* world, double t);
    void writeToH5();
    void writeToDat();
private:
    double boxsize;
    unsigned particleTypeId;
    unsigned long numberOfParticles;
    std::vector<unsigned long long> uniqueIds;
    // full information on the position of a particle
    struct PositionTuple {
        std::vector<double> position;
        std::vector<long> boxCoordinates;
    };
    // holds the positions of particle one timestep before
    std::vector<PositionTuple> lastPositions;
    // the increments that are extracted during one record()
    struct Snapshot {
        std::vector<double> increments;
        std::vector<bool> exists;
        double time;
    };
    // the trajectory of increments
    std::vector<Snapshot> incrementsTrajectory;

    std::map<unsigned long, bool> stillExists;

    void writeToNewH5();

    void bufferTimeDependentData(boost::multi_array<double, 2> &increments, boost::multi_array<int, 2> &exists,
                                 std::vector<double> &times);

    void appendToH5();
};

#endif //REVREADDY_INCREMENTS_H
