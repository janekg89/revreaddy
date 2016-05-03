//
// Created by chris on 03.05.16.
//

#include "Increments.h"

Increments::Increments(unsigned long inRecPeriod, unsigned long inClearPeriod, std::string inFilename, unsigned int pTypeId) {
    if (inClearPeriod < inRecPeriod) {
        throw Exception("The clearing-period of Increments must be larger or equal to the recording-period.");
    }
    recPeriod = inRecPeriod;
    clearedAutomatically = true;
    clearPeriod = inClearPeriod;
    filename = inFilename;
    observableTypeName = "Increments";
    particleTypeId = pTypeId;
    isSetup = false;
}

void Increments::setup(World *world, Config *config) {
    /* add uniqueIds of all particles of type 'particleTypeId' to uniqueIds */
    for (auto&& p : world->particles) {
        if (p.typeId == particleTypeId) {
            uniqueIds.push_back( p.uniqueId );
            PositionTuple pt;
            pt.position = p.position;
            pt.boxCoordinates = p.boxCoordinates;
            lastPositions.push_back(pt);
        }
    }
    for (auto i=0; i<uniqueIds.size(); ++i) {
        stillExists[i] = true;
    }
    numberOfParticles = uniqueIds.size();
}

void Increments::configure(World *world, Config *config) {
    /* keep track of boxsize that might have changed. */
    boxsize = config->boxsize;
}

void Increments::record(World *world, double t) {
    Snapshot snapshot;
    snapshot.time = t;
    snapshot.exists.resize(numberOfParticles);
    snapshot.increments.resize(numberOfParticles);
    std::vector<double> displacement = {0.,0.,0.};
    std::vector<long> relativeBoxCoordinates = {0,0,0};
    double squaredDisplacement;
    for (auto i=0; i<uniqueIds.size(); ++i) {
        if ( stillExists[i] ) {
            auto index = this->findParticleIndex(world->particles, uniqueIds[i]);
            if (index == -1) { //particle not found
                stillExists[i] = false;
                snapshot.exists[i] = false;
                snapshot.increments[i] = 0.;
            } else { // particle does exist
                snapshot.exists[i] = true;
                // find increment of particle with uniqueId[i]/index
                relativeBoxCoordinates[0] = world->particles[index].boxCoordinates[0]
                                            - lastPositions[i].boxCoordinates[0];
                relativeBoxCoordinates[1] = world->particles[index].boxCoordinates[1]
                                            - lastPositions[i].boxCoordinates[1];
                relativeBoxCoordinates[2] = world->particles[index].boxCoordinates[2]
                                            - lastPositions[i].boxCoordinates[2];
                displacement[0] = (double) relativeBoxCoordinates[0] * this->boxsize
                        + world->particles[index].position[0]
                        - lastPositions[i].position[0];
                displacement[1] = (double) relativeBoxCoordinates[1] * this->boxsize
                        + world->particles[index].position[1]
                        - lastPositions[i].position[1];
                displacement[2] = (double) relativeBoxCoordinates[2] * this->boxsize
                        + world->particles[index].position[2]
                        - lastPositions[i].position[2];
                squaredDisplacement = displacement[0]*displacement[0]
                                    + displacement[1]*displacement[1]
                                    + displacement[2]*displacement[2];
                snapshot.increments[i] = squaredDisplacement;
                // set lastPositions
                this->lastPositions[i].position = world->particles[index].position;
                this->lastPositions[i].boxCoordinates = world->particles[index].boxCoordinates;
            }
        } else { // particle is known not to exist
            snapshot.exists[i] = false;
            snapshot.increments[i] = 0.;
        }
    }
    // append to increments Traj
    incrementsTrajectory.push_back(snapshot);
}

void Increments::writeToH5() {
    // if file exists, append to increments, exists and time in first dimension
    boost::filesystem::path path( this->filename.c_str() );
    if ( boost::filesystem::exists(path) ) {
        this->appendToH5();
    }
        // if file does not exist, create it and write datasets to i
        // chunksize of file depends on clearPeriod. one chunk of position is [clearPeriod, N, 3]
    else {
        this->writeToNewH5();
    }
}

void Increments::writeToDat() {
    this->writeToH5();
}

void Increments::writeToNewH5() {
    H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
    auto T = incrementsTrajectory.size();
    auto N = numberOfParticles;
    // buffer time dependent data here
    boost::multi_array<double,2> increments(boost::extents[T][N]);
    boost::multi_array<int,2> exists(boost::extents[T][N]);
    std::vector<double> times;
    bufferTimeDependentData(increments, exists, times);
    // write the time dependent data to new datasets
    createExtendibleDataset(file, "increments", increments);
    createExtendibleDataset(file, "exists", exists);
    createExtendibleDataset(file, "times", times);
    // write time-independent information
    createExtendibleDataset(file, "uniqueIds", this->uniqueIds);
    // finally clear the buffer
    this->incrementsTrajectory.clear();
}

void Increments::bufferTimeDependentData(boost::multi_array<double, 2> &increments, boost::multi_array<int, 2> &exists,
                                         std::vector<double> &times) {
    for (auto i=0; i<incrementsTrajectory.size(); ++i) {
        times.push_back( incrementsTrajectory[i].time );
        for (auto j=0; j<numberOfParticles; ++j) {
            increments[i][j] = incrementsTrajectory[i].increments[j];
            /* hdf5 supports no boolean. Therefore: 0 -false, 1 - true */
            if (incrementsTrajectory[i].exists[j]) {
                exists[i][j] = 1;
            } else {
                exists[i][j] = 0;
            }

        }
    }
}

void Increments::appendToH5() {
    H5::H5File file(this->filename.c_str(), H5F_ACC_RDWR);
    auto T = incrementsTrajectory.size();
    auto N = numberOfParticles;
    // buffer time dependent data
    boost::multi_array<double,2> increments(boost::extents[T][N]);
    boost::multi_array<int,2> exists(boost::extents[T][N]);
    std::vector<double> times;
    bufferTimeDependentData(increments, exists, times);
    // write the extended dataset
    extendDataset(file, "increments", increments);
    extendDataset(file, "exists", exists);
    extendDataset(file, "times", times);
    // clear the observable buffer
    this->incrementsTrajectory.clear();
}

















