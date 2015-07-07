/* BinaryFile.h 
 * 
 * Wrap H5LT functions further for use of std vectors */

#ifndef __BINARYFILE_H_INCLUDED__
#define __BINARYFILE_H_INCLUDED__
#include <vector>
#include <string>
#include <iostream>
#include "hdf5.h"
#include "hdf5_hl.h"

class BinaryFile
{
	private:
		hid_t fileId;
		herr_t status;
	public:
		BinaryFile(std::string name); // open new file
		~BinaryFile(); // close the file and remove this object
		// create a new dataset in the file
		void addDatasetDouble(
			std::string name,
			std::vector<double>& vec);
		void addDatasetUnsignedInt(
			std::string name,
			std::vector<unsigned int>& vec);
};

#endif // __BINARYFILE_H_INCLUDED__
