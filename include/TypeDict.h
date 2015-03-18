/* TypeDict.h
 * This structure contains information on every particle type. */

#ifndef __TYPEDICT_H_INCLUDED__
#define __TYPEDICT_H_INCLUDED__
#include <vector>
#include <string>

class TypeDict
{
	private:
		unsigned int numberOfTypes;
	public: 
		std::vector<string> names;
		std::vector<double> radii;
		std::vector<double> diffusionConstants;

		void newType(string name, double radius, double diffusionConst);
		unsigned int getNumberOfTypes();

		TypeDict();
		~TypeDict();
};

#endif // __TYPEDICT_H_INCLUDED__
