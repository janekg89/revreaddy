/* Observable.cpp */

#include "Observable.h"

Observable::Observable() {}

Observable::~Observable() {}

void Observable::configure(World * world, Config * config) {}
void Observable::record(World * world, double t) {}
void Observable::setup(World * world, Config * config) {}

void Observable::writeToFile() {
	// first determine the file extension
	unsigned int lastDot = this->filename.find_last_of(".");
	std::string extension = this->filename.substr(lastDot);
	if ( (extension == ".h5") || (extension == ".hdf5") ) {
		writeToH5();
	}
	else if ( (extension == ".dat") || (extension == ".txt") ) {
		writeToDat();
	}
	else {
		// write to hdf5 format per default
		writeToH5();
	}
}

void Observable::writeToH5() {
	LOG_TRACE("Enter Observable::writeToH5")
}

void Observable::writeToDat() {
	LOG_TRACE("Enter Observable::writeToDat")
}

// return the particle index within particles of particle with id id
long int Observable::findParticleIndex(std::vector<Particle>& particles, unsigned long long id) { 
	signed long max = particles.size() - 1;
	signed long min = 0;
	signed long mid = 0;
	while (max >= min) {
		mid = min + (max - min) / 2;
		if (particles[mid].uniqueId == id) {return mid;}
		else if (particles[mid].uniqueId < id) {min = mid + 1;}
		else {max = mid - 1;}
	}
	// std::cout << "id was not found\n"; //this has to be catched
	return -1;
}

// IO methods that are not templated and should not be inlined
/* Create and write a std vector of doubles */
void  Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<double>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_DOUBLE, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_DOUBLE);
}

/* Create and write a std vector of bools */
void  Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<int>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_INT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_INT);
}

/* Create and write a std vector of unsigned ints */
void  Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<unsigned int>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_UINT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_UINT);
}

/* Create and write a std vector of unsigned longs */
void  Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<unsigned long>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_UINT, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_ULONG);
}

/* Create and write a std vector of unsigned long longs */
void  Observable::createAndWrite(H5::H5File& file,  std::string name, std::vector<unsigned long long>& arr,  H5::DataSpace& dspace, const H5::DSetCreatPropList& create_plist) {
	H5::DataSet dset = file.createDataSet(name.c_str(), H5::PredType::NATIVE_ULLONG, dspace, create_plist);
	dset.write(arr.data(), H5::PredType::NATIVE_ULLONG);
}

/* Write a std vector of doubles to an extended dset*/
void  Observable::writeToExtended(std::vector<double>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_DOUBLE, mspace, fspace);
}

/* Write a std vector of ints to an extended dset*/
void  Observable::writeToExtended(std::vector<int>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_INT, mspace, fspace);
}

/* Write a std vector of uints to an extended dset*/
void  Observable::writeToExtended(std::vector<unsigned int>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_UINT, mspace, fspace);
}

/* Write a std vector of ulongs to an extended dset*/
void  Observable::writeToExtended(std::vector<unsigned long>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_ULONG, mspace, fspace);
}

/* Write a std vector of unsigned long longs to an extended dset*/
void  Observable::writeToExtended(std::vector<unsigned long long>& arr, H5::DataSet& dset, H5::DataSpace& mspace, H5::DataSpace& fspace) {
	dset.write(arr.data(), H5::PredType::NATIVE_ULLONG, mspace, fspace);
}

bool Observable::shallBeRecorded(unsigned long timeIndex) {
	return timeIndex % this->recPeriod == 0;
}
