//
// Created by janekg89 on 18.08.16.
//

#ifndef REVREADDY_REACTIONCOUNTER_H
#define REVREADDY_REACTIONCOUNTER_H

#include "Observable.h"
#include "World.h"
#include <string>
#include <vector>

class ReactionCounter : public Observable {
public:
    ReactionCounter(unsigned inRecPeriod, unsigned inClearPeriod, std::string inFilename, std::string reactionName);
    void record(World * world, double t);
    void writeToH5();
    void writeToDat();
private:
    std::vector<unsigned int> forwardCounter;
    std::vector<unsigned int> backwardCounter;
    std::vector<double> times;
    std::string reactionName;
};

#endif //REVREADDY_REACTIONCOUNTER_H
