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
		std::vector<std::string> names;
		std::vector<double> radii;
		std::vector<double> diffusionConstants;
		std::vector<double> reactionRadii;

		void newType(
			std::string name,
			double radius,
			double diffusionConst,
			double reactionRadius);
		unsigned int getNumberOfTypes();

		TypeDict();
		~TypeDict();
};

#endif // __TYPEDICT_H_INCLUDED__
