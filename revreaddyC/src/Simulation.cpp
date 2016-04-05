/* Simulation.cpp */
#include "Simulation.h"

Simulation::Simulation(World * inWorld, Config * inConfig, std::string whichImpl) {
	LOG_TRACE("Enter Simulation Constructor.")
	if (whichImpl == "RevDiffusion") { 
		this->impl = new RevDiffusion(inWorld, inConfig); 
	}
	else if (whichImpl == "RevReactions") {
		this->impl = new RevReactions(inWorld, inConfig);
	}
	else if (whichImpl == "RevDiffusionReactions") { }
	else { this->impl = new SimulationImpl(inWorld, inConfig); }
	LOG_TRACE("Leave Simulation Constructor.")
}

Simulation::~Simulation() {
	LOG_TRACE("Enter Simulation Destructor.")
	delete this->impl;
	LOG_TRACE("Leave Simulation Destructor.")
}

void Simulation::run(const unsigned long maxTime) { impl->run(maxTime); }
bool Simulation::getUseNeighborlist() { return impl->useNeighborlist; }
void Simulation::setUseNeighborlist(const bool inUseNeighborlist) {impl->useNeighborlist = inUseNeighborlist; }

void Simulation::writeAllObservablesToFile() { impl->writeAllObservablesToFile(); }
std::string Simulation::showObservables() { return impl->showObservables(); }
void Simulation::deleteAllObservables() { impl->deleteAllObservables(); }
void Simulation::new_Trajectory(unsigned long recPeriod, std::string filename) { impl->new_Trajectory(recPeriod, filename); }
void Simulation::new_RadialDistribution(unsigned long recPeriod, std::string filename, std::vector<double> ranges, std::vector< std::array<unsigned,2> > considered) {
	impl->new_RadialDistribution(recPeriod, filename, ranges, considered);
}
void Simulation::new_MeanSquaredDisplacement(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
	impl->new_MeanSquaredDisplacement(recPeriod, filename, particleTypeId);
}
void Simulation::new_ProbabilityDensity(unsigned long recPeriod, std::string filename, unsigned pTypeId, std::vector<double> range, unsigned int coord) {
	impl->new_ProbabilityDensity(recPeriod, filename, pTypeId, range, coord);
}
void Simulation::new_Energy(unsigned long recPeriod, std::string filename) { impl->new_Energy(recPeriod, filename); }
void Simulation::new_Acceptance(unsigned long recPeriod, std::string filename, bool reactionsOrDiffusion) {
	impl->new_Acceptance(recPeriod, filename, reactionsOrDiffusion);
}
void Simulation::new_ParticleNumbers(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
	impl->new_ParticleNumbers(recPeriod, filename, particleTypeId);
}