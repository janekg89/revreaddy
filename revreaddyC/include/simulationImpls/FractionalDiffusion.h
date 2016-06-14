/* FractionalDiffusion.h */
#ifndef __FRACTIONALDIFFUSION_H_INCLUDED__
#define __FRACTIONALDIFFUSION_H_INCLUDED__

class FractionalDiffusion : public SimulationImpl {
public:
	FractionalDiffusion(World * inWorld, Config * inConfig);
	void run(const unsigned long maxTime);
	void propagateDiffusion();
};
#endif//__FRACTIONALDIFFUSION_H_INCLUDED__