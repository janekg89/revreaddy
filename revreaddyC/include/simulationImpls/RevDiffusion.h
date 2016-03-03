/* RevDiffusion.h - child of SimulationImpl that performs the diffusion with
 * the reversible scheme (acceptance step) and the reactions without. Thus
 * only the run() method is altered compared to the default implementation.*/
 #ifndef __REVDIFFUSION_H_INCLUDED__
 #define __REVDIFFUSION_H_INCLUDED__
 #include "SimulationImpl.h"

 class RevDiffusion : public SimulationImpl {
 public:
 	RevDiffusion(World * inWorld, Config * inConfig);
 	~RevDiffusion();
 	void run(const unsigned long maxTime);
 };
 #endif //__REVDIFFUSION_H_INCLUDED__