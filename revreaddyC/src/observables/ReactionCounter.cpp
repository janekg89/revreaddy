//
// Created by janekg89 on 18.08.16.
//

#include "ReactionCounter.h"

ReactionCounter::ReactionCounter(unsigned inRecPeriod, unsigned inClearPeriod, std::string inFilename, std::string reactionName) {
    this->recPeriod = inRecPeriod;
    this->clearPeriod = inClearPeriod;
    this->clearedAutomatically = false;
    this->filename = inFilename;
    this->reactionName = reactionName;
    observableTypeName = "ReactionCounter";
    isSetup = false;
}

void ReactionCounter::record(World * world, double t) {
    this->forwardCounter.push_back(world->forwardReactionCounter[reactionName]);
    this->backwardCounter.push_back(world->backwardReactionCounter[reactionName]);
    this->times.push_back(t);
}

void ReactionCounter::writeToH5() {
    H5::H5File file(this->filename.c_str(), H5F_ACC_TRUNC);
    createExtendibleDataset(file, "forwardCounter", this->forwardCounter);
    createExtendibleDataset(file, "backwardCounter", this->backwardCounter);
    createExtendibleDataset(file, "times", this->times);
}

void ReactionCounter::writeToDat() {
    ReactionCounter::writeToH5();
}

