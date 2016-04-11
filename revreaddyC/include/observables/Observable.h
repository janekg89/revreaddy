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
	Observable();
	virtual ~Observable();

	// decide if the observable should be written during runtime or manually
	bool clearedAutomatically;
	// some observables require setup that may only happen once in contrast to configure
	bool isSetup;
	// number of timesteps between two record() calls
	unsigned long int recPeriod; 
	// number of timesteps between two writeToFile() calls
	unsigned long int clearPeriod; 
	std::string filename;
	std::string observableTypeName;
	std::vector<unsigned long long> observedParticleIds;
	/* configure sets up the observables' parameters that could have
	 * changed between construction of the observable and the start
	 * of the simulation */
	virtual void configure(World * world, Config * config);
	/* setup cares about the initialization directly before the run. In
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

	// HDF5 stuff
	/* The following templates create and extend datasets which are extendible in the first dimension */
	template<typename T, unsigned long rank>
	void createExtendibleDataset(H5::H5File& file,  std::string name,  boost::multi_array<T, rank>& arr);
	template<typename T>
	void createExtendibleDataset(H5::H5File& file,  std::string name, std::vector<T>& arr);
	template<unsigned long rank>
	void createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<double, rank>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	template<unsigned long rank>
	void createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<int, rank>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	template<unsigned long rank>
	void createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<unsigned int, rank>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	template<unsigned long rank>
	void createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<unsigned long long, rank>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	void createAndWrite(H5::H5File& file,  std::string name, std::vector<double>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	void createAndWrite(H5::H5File& file,  std::string name, std::vector<int>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	void createAndWrite(H5::H5File& file,  std::string name, std::vector<unsigned int>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	void createAndWrite(H5::H5File& file,  std::string name, std::vector<unsigned long long>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist = H5::DSetCreatPropList::DEFAULT);
	template<typename T, unsigned long rank>
	void extendDataset(H5::H5File& file,  std::string name,  boost::multi_array<T, rank>& arr);
	template<typename T>
	void extendDataset(H5::H5File& file,  std::string name,  std::vector<T>& arr);
	template<unsigned long rank>
	void writeToExtended(boost::multi_array<double, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);
	template<unsigned long rank>
	void writeToExtended(boost::multi_array<int, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);
	template<unsigned long rank>
	void writeToExtended(boost::multi_array<unsigned int, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);
	template<unsigned long rank>
	void writeToExtended(boost::multi_array<unsigned long long, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);
	void writeToExtended(std::vector<double>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);
	void writeToExtended(std::vector<int>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);
	void writeToExtended(std::vector<unsigned int>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);
	void writeToExtended(std::vector<unsigned long long>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace);

};

/* hdf5 wrappers for boost multi arrays */
template<typename T, unsigned long rank>
void Observable::createExtendibleDataset(H5::H5File& file,  std::string name,  boost::multi_array<T, rank>& arr) {
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

/* Create and write a boost array of doubles. */
template<unsigned long rank> void Observable::createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<double, rank>& arr, H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_DOUBLE);
} 

/* Create and write a boost array of bools. */
template<unsigned long rank> void Observable::createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<int, rank>& arr, H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_INT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_INT);
}

/* Create and write a boost array of unsigned ints. */
template<unsigned long rank> void Observable::createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<unsigned int, rank>& arr, H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_UINT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_UINT);
}

/* Create and write a boost array of unsigned long longs. */
template<unsigned long rank> void Observable::createAndWrite(H5::H5File& file,  std::string name,  boost::multi_array<unsigned long long, rank>& arr, H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_ULLONG, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_ULLONG);
}

/* hdf5 wrappers for std vector */
template<typename T>
void Observable::createExtendibleDataset(H5::H5File& file,  std::string name, std::vector<T>& arr) {
	hsize_t dims[1], maxdims[1], chunkdims[1];
	dims[0] = arr.size(); 
	maxdims[0] = H5S_UNLIMITED;
	chunkdims[0] = 100;
	H5::DataSpace dspace(1, dims, maxdims);
	H5::DSetCreatPropList plist;
	plist.setChunk(1, chunkdims);
	createAndWrite(file, name, arr, dspace, plist);
}

/* Create and write a std vector of doubles */
void inline Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<double>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_DOUBLE);
}

/* Create and write a std vector of bools */
void inline Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<int>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_INT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_INT);
}

/* Create and write a std vector of unsigned ints */
void inline Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<unsigned int>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_UINT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_UINT);
}

/* Create and write a std vector of unsigned long longs */
void inline Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<unsigned long long>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_ULLONG, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_ULLONG);
}

/* Extend boost array */
template<typename T, unsigned long rank>
void Observable::extendDataset(H5::H5File& file,  std::string name,  boost::multi_array<T, rank>& arr) {
	H5::DataSet dset = file.openDataSet(name.c_str());
	H5::DataSpace fspace(dset.getSpace());
	// old dims
	hsize_t oldDims[rank];
	fspace.getSimpleExtentDims(oldDims); // oldDims is result array
	// dims of extending data (arr)
	hsize_t extDims[rank];
	for (auto i=0; i<rank; ++i) {
		extDims[i] = arr.shape()[i];
	}
	// check if extDims is compatible with oldDims. I.e. extDims[i] == oldDims[i], except for i=0
	for (auto i=1; i<rank; ++i) {
		if (! (extDims[i] == oldDims[i]) ) {
			throw Exception("Observable::extendDataset. The dimensions of the existing data on file are not compatible with the new data to be written.");
		}
	}
	// new dims are same as old dims except for the first dimension. This gets extended by new data
	hsize_t newDims[rank];
	for (auto i=1; i<rank; ++i) {
		newDims[i] = oldDims[i];
	}
	newDims[0] = oldDims[0] + extDims[0];
	dset.extend(newDims);
	// update filespace
	fspace = H5::DataSpace(dset.getSpace());
	// select offset and according hyperslabs
	hsize_t offset[rank];
	for (auto i=0; i<rank; ++i) {
		offset[i] = 0;
	}
	offset[0] = oldDims[0];
	fspace.selectHyperslab(H5S_SELECT_SET, extDims, offset);
	// define memspaces from which data is read
	H5::DataSpace mspace(rank, extDims);
	writeToExtended(arr, dset, mspace, fspace);
}

template<unsigned long rank>
void Observable::writeToExtended(boost::multi_array<double, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_DOUBLE, mspace, fspace);
}

template<unsigned long rank>
void Observable::writeToExtended(boost::multi_array<int, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_INT, mspace, fspace);
}

template<unsigned long rank>
void Observable::writeToExtended(boost::multi_array<unsigned int, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_UINT, mspace, fspace);
}

template<unsigned long rank>
void Observable::writeToExtended(boost::multi_array<unsigned long long, rank>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_ULLONG, mspace, fspace);
}

/* Extend std vector. */
template<typename T>
void Observable::extendDataset(H5::H5File& file,  std::string name,  std::vector<T>& arr) {
	H5::DataSet dset = file.openDataSet(name.c_str());
	H5::DataSpace fspace(dset.getSpace());
	// old dims
	hsize_t oldDims[1];
	fspace.getSimpleExtentDims(oldDims); // oldDims is result array
	// dims of extending data (arr)
	hsize_t extDims[1];
	extDims[0] = arr.size();
	hsize_t newDims[1];
	newDims[0] = oldDims[0] + extDims[0];
	dset.extend(newDims);
	// update filespace
	fspace = H5::DataSpace(dset.getSpace());
	// select offset and according hyperslabs
	hsize_t offset[1];
	offset[0] = oldDims[0];
	fspace.selectHyperslab(H5S_SELECT_SET, extDims, offset);
	// define memspaces from which data is read
	H5::DataSpace mspace(1, extDims);
	writeToExtended(arr, dset, mspace, fspace);
}

void inline Observable::writeToExtended(std::vector<double>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_DOUBLE, mspace, fspace);
}

void inline Observable::writeToExtended(std::vector<int>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_INT, mspace, fspace);
}

void inline Observable::writeToExtended(std::vector<unsigned int>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_UINT, mspace, fspace);
}

void inline Observable::writeToExtended(std::vector<unsigned long long>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_ULLONG, mspace, fspace);
}

#endif // __OBSERVABLE_H_INCLUDED__
