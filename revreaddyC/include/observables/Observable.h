/* Observable.h
 * 
 * Implement different observables objects, which will be 
 * recorded during simulation. 
 */

#ifndef __OBSERVABLE_H_INCLUDED__
#define __OBSERVABLE_H_INCLUDED__
#include <vector>
#include <string>
#include <fstream>
#include <hdf5.h> // hdf5 c api
#include <H5Cpp.h> // hdf5 c++ api
#include <boost/filesystem.hpp>
#include <boost/multi_array.hpp>
#include "World.h"
#include "Config.h"
#include "Particle.h"
#include "BinaryFile.h"

class Observable {
public:
	// decide if the observable should be written during runtime or manually
	bool clearedAutomatically;
	// some observables require setup that may only happen once in contrast to configure
	bool isSetup;
	// number of timesteps between two recordings
	unsigned long int recPeriod; 
	// number of recordings between two writeBufferToFile
	// not yet used
	unsigned long int clearPeriod; 
	std::string filename;
	std::string observableTypeName;
	std::vector<unsigned long long> observedParticleIds;
	/* configure sets up the observables' parameters that could have
	 * changed between construction of the observable and the start
	 * of the simulation */
	virtual void configure(World * world, Config * config);
	/* setup cares about the initialization direclty befire the run. In
	 * contrast to configure, this must only happen once, e.g. saving
	 * uniqueIds of considered particles. */
	virtual void setup(World * world, Config * config);
	/* record data from World and write to the buffer */
	virtual void record(World * world, double t);
	/* write to file and clear the buffer */
	void writeToFile();
	virtual void writeToH5();
	virtual void writeToDat();
	long int findParticleIndex(std::vector<Particle>& particles, unsigned long long id);
	/* The following templates create and extend datasets which are extendible in the first dimension */
	template<typename T, unsigned int rank>
	void createExtendibleDatasetFromBoostArray(H5::H5File& file, const std::string& name, const boost::multi_array<T, rank>& arr);
	template<typename T, unsigned int rank>
	void extendDatasetFromBoostArray(H5::H5File& file, const std::string& name, const boost::multi_array<T, rank>& arr);

	Observable();
	virtual ~Observable();		
};

template<typename T, unsigned int rank>
void createExtendibleDatasetFromBoostArray(H5::H5File& file, const std::string& name, const boost::multi_array<T, rank>& arr) {
	hsize_t dims[rank], maxdims[rank], chunkdims[rank];
	for (auto i=0; i<rank; ++i) {
		dims[i] = arr.shape()[i];
		maxdims[i] = arr.shape()[i];
		chunkdims[i] = arr.shape()[i];
	}
	maxdims[0] = H5S_UNLIMITED; // extendible in first dimension
	chunkdims[0] = 100; // chunkdims are 100 x dims[1] x dims[2] x ...
	H5::DataSpace dspace(rank, dims, maxdims);
	H5::DSetCreatPropList plist;
	plist.setChunk(rank, chunkdims);
	createAndWrite(file, name, arr, dspace, plist);
}

/* Create and write an array of doubles. */
template<unsigned int rank> void createAndWrite(
	H5::H5File& file,
	const std::string& name,
	const boost::multi_array<double, rank>& arr,
	const H5::DataSpace& dspace,
	const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT
) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_DOUBLE);
} 

/* Create and write an array of bools. */
template<unsigned int rank> void createAndWrite(
	H5::H5File& file,
	const std::string& name,
	const boost::multi_array<bool, rank>& arr,
	const H5::DataSpace& dspace,
	const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT
) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_HBOOL, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_HBOOL);
}

/* Create and write an array of unsigned ints. */
template<unsigned int rank> void createAndWrite(
	H5::H5File& file,
	const std::string& name,
	const boost::multi_array<unsigned int, rank>& arr,
	const H5::DataSpace& dspace,
	const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT
) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_UINT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_UINT);
}

/* Create and write an array of unsigned long longs. */
template<unsigned int rank> void createAndWrite(
	H5::H5File& file,
	const std::string& name,
	const boost::multi_array<unsigned long long, rank>& arr,
	const H5::DataSpace& dspace,
	const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT
) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_ULLONG, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_ULLONG);
}

#endif // __OBSERVABLE_H_INCLUDED__
