/* Simulation.cpp */
#include "Simulation.h"

Simulation::Simulation(World * inWorld, Config * inConfig, std::string whichImpl) {
	if (whichImpl == "revDynamics") { }
	else if (whichImpl == "revReactions") { }
	else if (whichImpl == "revDynamicsReactions") { }
	else { this->impl = new SimulationImpl(inWorld, inConfig); }
}

Simulation::~Simulation() { delete this->impl; }

void Simulation::run(const unsigned long maxTime) { impl->run(maxTime); }
bool Simulation::getUseNeighborlist() { return impl->useNeighborlist; }
void Simulation::setUseNeighborlist(const bool inUseNeighborlist) {impl->useNeighborlist = inUseNeighborlist; }

void Simulation::writeAllObservablesToFile() { impl->writeAllObservablesToFile(); }
void Simulation::writeLastObservableToFile() { impl->writeLastObservableToFile(); }
std::string Simulation::showObservables() { return impl-> showObservables(); }
void Simulation::deleteAllObservables() { impl->deleteAllObservables(); }
void Simulation::deleteLastObservable() { impl->deleteLastObservable(); }
void Simulation::new_Trajectory(unsigned long recPeriod, std::string filename) { impl->new_Trajectory(recPeriod, filename); }
void Simulation::new_RadialDistribution(unsigned long recPeriod, std::string filename, std::vector<double> ranges, std::vector< std::vector<unsigned> > considered) {
	impl->new_RadialDistribution(recPeriod, filename, ranges, considered);
}
void Simulation::new_MeanSquaredDisplacement(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
	impl->new_MeanSquaredDisplacement(recPeriod, filename, particleTypeId);
}
void Simulation::new_ProbabilityDensity(unsigned long recPeriod, std::string filename, unsigned pTypeId, std::vector<double> range, unsigned int coord) {
	impl->new_ProbabilityDensity(recPeriod, filename, pTypeId, range, coord);
}
void Simulation::new_Energy(unsigned long recPeriod, std::string filename) { impl->new_Energy(recPeriod, filename); }
void Simulation::new_Acceptance(unsigned long recPeriod, std::string filename, bool reactionsOrDynamics) {
	impl->new_Acceptance(recPeriod, filename, reactionsOrDynamics);
}
void Simulation::new_ParticleNumbers(unsigned long recPeriod, std::string filename,	unsigned particleTypeId) {
	impl->new_ParticleNumbers(recPeriod, filename, particleTypeId);
}