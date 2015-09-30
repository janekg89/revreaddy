/* Neighborlist.h 
 * A structure that contains a list of particles for
 * every "box" with index [x][y][z] */

#ifndef __NEIGHBORLIST_H_INCLUDED__
#define __NEIGHBORLIST_H_INCLUDED__
#include <vector>

struct Neighborlist
{
	public:
	Neighborlist(unsigned inNumberBoxes)
	{
		this->numberBoxes = inNumberBoxes;
		this->neighbors = std::vector<
			std::vector<
				std::vector<
					std::vector<unsigned long>
				>
			>
		> (numberBoxes,
			std::vector<
				std::vector<
					std::vector<unsigned long>
				>
			> (numberBoxes,
				std::vector<
					std::vector<unsigned long>
				> (numberBoxes,
					std::vector<unsigned long>(0)
				)
			)
		);
	}
	~Neighborlist(){this->neighbors.clear();}

	unsigned numberBoxes;
	std::vector< // x
		std::vector< // y
			std::vector< // z
				std::vector<unsigned long> // actual list
			>
		>
	> neighbors;

	void addIndex(unsigned x, unsigned y, unsigned z, unsigned long index);
	unsigned long getIndex(unsigned x, unsigned y, unsigned z, unsigned long position);	
	unsigned long getSize(unsigned x, unsigned y, unsigned z);
};

inline void Neighborlist::addIndex(
	unsigned x,
	unsigned y,
	unsigned z,
	unsigned long index)
{
	this->neighbors[x][y][z].push_back(index);
}

inline unsigned long Neighborlist::getIndex(
	unsigned x,
	unsigned y,
	unsigned z,
	unsigned long position)
{
	return this->neighbors[x][y][z].operator[](position);
}

inline unsigned long Neighborlist::getSize(
	unsigned x,
	unsigned y,
	unsigned z)
{
	return this->neighbors[x][y][z].size();
}

#endif // __NEIGHBORLIST_H_INCLUDED__