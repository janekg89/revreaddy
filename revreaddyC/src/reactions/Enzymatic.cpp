//
// Created by chris on 22.04.16.
//

#include "Enzymatic.h"

Enzymatic::Enzymatic(
        std::string inName, std::vector<unsigned> inForwardTypes, std::vector<unsigned> inBackwardTypes,
        double inForwardRate, double inBackwardRate, double inReactionDistance) {
    /* forwardTypes = [A, Catalyst] and backwardTypes = [B, Catalyst] */
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
        /* find out which particles is the substrate */
        unsigned long substrateIndex;
        if (world->particles[particleIndices[0]].typeId==forwardTypes[0]) {
            substrateIndex = particleIndices[0];
        } else if (world->particles[particleIndices[1]].typeId==forwardTypes[0]) {
            substrateIndex = particleIndices[1];
        } else {
            throw Exception("Exceptional state in performForward of Enzymatic. "
                                    "Neither of participants has the right type");
        }
        /* reaction occurs. substrate's type is changed from forwardTypes[0] to backwardTypes[0] */
        auto pos = world->particles[substrateIndex].position;
        if (world->useFractional) {
            auto alpha = world->alpha;
            world->removeParticleAndIncrements(substrateIndex);
            world->addParticleAndIncrements(pos, backwardTypes[0], random, maxTime, timestep, diffB, alpha);
        }
        else {
            world->setTypeId(substrateIndex, backwardTypes[0]);
        }
    }
    return 1.;
}

double Enzymatic::performBackward(std::vector<unsigned long> particleIndices, double timestep, World *world,
                                  Random *random) {
    double backwardProb = backwardRate * timestep;
    double u = random->uniform();
    if (u < backwardProb) {        /* find out which particles is the substrate */
        unsigned long substrateIndex;
        if (world->particles[particleIndices[0]].typeId==backwardTypes[0]) {
            substrateIndex = particleIndices[0];
        } else if (world->particles[particleIndices[1]].typeId==backwardTypes[0]) {
            substrateIndex = particleIndices[1];
        } else {
            throw Exception("Exceptional state in performForward of Enzymatic. "
                                    "Neither of participants has the right type");
        }
        /* reaction occurs, particle i's type is changed to forwardTypes[0] */
        // TODO
        world->setTypeId(substrateIndex, forwardTypes[0]);
    }
    return 1.;
}

void Enzymatic::configure(unsigned long inMaxTime, double inDiffA, double inDiffB, double inDiffC) {
    this->maxTime = inMaxTime;
    this->diffA = inDiffA;
    this->diffB = inDiffB;
    this->diffC = inDiffC;
}








