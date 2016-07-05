//
// Created by chris on 22.04.16.
//

#ifndef REVREADDY_ENZYMATIC_H
#define REVREADDY_ENZYMATIC_H

#include "Reaction.h"

class Enzymatic : public Reaction {
public:
    Enzymatic(std::string inName, std::vector<unsigned> inForwardTypes, std::vector<unsigned> inBackwardTypes, double inForwardRate, double inBackwardRate, double inReactionDistance);

    double performForward(std::vector<unsigned long> particleIndices, double timestep, World * world, Random * random);
    double performBackward(std::vector<unsigned long> particleIndices, double timestep, World * world, Random * random);
    unsigned long maxTime;
    double diffA;
    double diffB;
    double diffC;
    void configure(unsigned long inMaxTime, double inDiffA, double inDiffB, double inDiffC);
};

#endif //REVREADDY_ENZYMATIC_H
