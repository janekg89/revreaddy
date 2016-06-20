//
// Created by chris on 22.04.16.
//

#include "Enzymatic.h"


Enzymatic::Enzymatic(
        std::string inName, std::vector<unsigned> inForwardTypes, std::vector<unsigned> inBackwardTypes,
        double inForwardRate, double inBackwardRate, double inReactionDistance) {
    this->name = inName;
    this->forwardTypes = inForwardTypes;
    this->backwardTypes = inBackwardTypes;
    this->forwardRate = inForwardRate;
    this->backwardRate = inBackwardRate;
    this->reactionDistance = inReactionDistance;
    this->type = "Enzymatic";
}

double Enzymatic::performForward(std::vector<unsigned long> particleIndices, double timestep, World *world,
                                 Random *random) {
    double forwardProb = forwardRate * timestep;
    double u = random->uniform();
    if (u < forwardProb) {
        // TODO this is wrong
        /* reaction occurs, particle i is the substrate and j is the catalyst.
         * particle i's type is changed from forwardTypes[0] to backwardTypes[0] */
        unsigned long indexI = particleIndices[0];
        std::vector<double> pos = world->particles[indexI].position;
        world->removeParticleAndIncrements(indexI);
        world->addParticleAndIncrements(pos, backwardTypes[0], random);
        //world->setTypeId(indexI, backwardTypes[0]);
    }
    return 1.;
}

double Enzymatic::performBackward(std::vector<unsigned long> particleIndices, double timestep, World *world,
                                  Random *random) {
    double backwardProb = backwardRate * timestep;
    double u = random->uniform();
    if (u < backwardProb) {
        /* reaction occurs, particle i's type is changed to forwardTypes[0] */
        unsigned long indexI = particleIndices[0];
        std::vector<double> pos = world->particles[indexI].position;
        world->removeParticleAndIncrements(indexI);
        world->addParticleAndIncrements(pos, forwardTypes[0], random);
    }
    return 1.;
}






