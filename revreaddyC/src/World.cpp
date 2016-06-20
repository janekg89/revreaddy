/* World.cpp */

#include "World.h"
#include "simulationImpls/generateIncrements.h"

World::World() {
    this->cumulativeRuntime = 0.;
    this->energy = 0.;
    this->oldEnergy = 0.;
    this->acceptProbDiffusion = 1.;
    this->acceptProbReactions = 1.;
    this->acceptionsDiffusion = 0;
    this->rejectionsDiffusion = 0;
    this->acceptionsReactions = 0;
    this->rejectionsReactions = 0;
    this->uniqueIdCounter = 0;
}

void World::addParticle(std::vector<double> initPos, unsigned particleTypeId) {
    try {
        if (initPos.size() != 3) {
            throw Exception("Particles' initial position has dimension mismatch.");
        }
    } catch (Exception &e) {
        LOG_ERROR(e.what())
        LOG_INFO("Particle is not created.")
        return;
    }
    Particle particle;
    particle.position = initPos;
    particle.typeId = particleTypeId;
    particle.uniqueId = this->uniqueIdCounter;
    this->uniqueIdCounter += 1;
    this->particles.push_back(particle);
}

void World::removeParticle(unsigned long index) {
    if (index >= this->particles.size()) {
        throw Exception("There is no particle of the given index.");
    }
    this->particles.erase(this->particles.begin() + index);
}

unsigned long World::getNumberOfParticles() { return this->particles.size(); }

std::vector<double> World::getPosition(unsigned long index) {
    if (this->particles.size() <= index) {
        throw Exception("There is no particle of the given index.");
    }
    return this->particles[index].position;
}

void World::setPosition(unsigned long index, std::vector<double> newPos) {
    try {
        if (newPos.size() != 3) { throw Exception("New position has dimension mismatch."); }
        this->particles[index].position[0] = newPos[0];
        this->particles[index].position[1] = newPos[1];
        this->particles[index].position[2] = newPos[2];
    } catch (Exception &e) {
        LOG_ERROR(e.what())
        LOG_INFO("Particle position is not changed.")
    }
}

unsigned World::getTypeId(unsigned long index) { return this->particles[index].typeId; }

void World::setTypeId(unsigned long index, unsigned typeId) { this->particles[index].typeId = typeId; }

unsigned long long World::getUniqueId(unsigned long index) { return this->particles[index].uniqueId; }

void World::deleteAllParticles() {
    this->particles.erase(
            this->particles.begin(),
            this->particles.begin() + this->particles.size()
    );
}

/* Fractional stuff */
void World::addParticleAndIncrements(std::vector<double> initPos, unsigned particleTypeId, Random *random,
                                     unsigned long maxTime, double &timestep, double &diffConst, double &alpha) {
    addParticle(initPos, particleTypeId);
    auto uid = this->particles.back().uniqueId;
    /*how to get variables*/
    increments[uid] = janek::generateIncrements(maxTime, diffConst, timestep, alpha, random);
    incrementsIndex[uid] = 0;
}

void World::removeParticleAndIncrements(unsigned long index) {
    auto uid = this->particles[index].uniqueId;
    removeParticle(index);
    increments.erase(uid);
    incrementsIndex.erase(uid);
}