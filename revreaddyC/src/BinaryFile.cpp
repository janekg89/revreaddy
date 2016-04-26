/* BinaryFile.cpp */

#include "BinaryFile.h"

BinaryFile::BinaryFile(std::string name) {
	this->fileId = H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	// if this returns a negative value, file could not be created
	if (this->fileId < 0) {
		throw Exception("Binary file could not be created");
	}
}

BinaryFile::~BinaryFile() {
	this->status = H5Fclose(this->fileId);
	if (this->status < 0) {
		throw Exception("Binary file could not be closed");
	}
}

void BinaryFile::addDatasetDouble(
	std::string name,
	std::vector<double>& vec)
{
	hsize_t dims[1] = { vec.size() };
	double * data = vec.data();
	this->status = H5LTmake_dataset_double(
		this->fileId,
		name.c_str(),
		1, // rank
		dims,
		data);
}

void BinaryFile::addDatasetUnsignedInt(
	std::string name,
	std::vector<unsigned int>& vec)
{
	hsize_t dims[1] = { vec.size() };
	unsigned int * data = vec.data();
	this->status = H5LTmake_dataset(
		this->fileId,
		name.c_str(),
		1,
		dims,
		H5T_NATIVE_UINT,
		data);
}

void BinaryFile::addDatasetUnsignedLong(
	std::string name,
	std::vector<unsigned long>& vec)
{
	hsize_t dims[1] = { vec.size() };
	unsigned long * data = vec.data();
	this->status = H5LTmake_dataset(
		this->fileId,
		name.c_str(),
		1,
		dims,
		H5T_NATIVE_ULONG,
		data);
}