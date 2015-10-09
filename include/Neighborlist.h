/* Neighborlist.h 
 * A structure that contains a list of particles for
 * every "box" with index [x][y][z] */

#ifndef __NEIGHBORLIST_H_INCLUDED__
#define __NEIGHBORLIST_H_INCLUDED__
#include <vector>

#define getBox(x,y,z) x + this->numberBoxes * y + this->numberBoxes * this->numberBoxes * z

struct Neighborlist
{
public:
	Neighborlist(unsigned inNumberBoxes)
	{
		this->numberBoxes = inNumberBoxes;
		this->totalNumberBoxes = numberBoxes * numberBoxes * numberBoxes;
		/*
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
		);*/
		this->neighbors = std::vector< std::vector<unsigned long> > (
			totalNumberBoxes, std::vector<unsigned long>(0)
		);
		this->sizes = std::vector<unsigned long> (totalNumberBoxes);
		std::fill(sizes.begin(),sizes.end(),0);
	}
	~Neighborlist(){this->neighbors.clear();}

	unsigned numberBoxes;
	unsigned totalNumberBoxes;
	/*
	std::vector< // x
		std::vector< // y
			std::vector< // z
				std::vector<unsigned long> // actual list
			>
		>
	> neighbors;
	*/
	std::vector< std::vector<unsigned long> > neighbors;
	std::vector< unsigned long> sizes;

	void addIndex(unsigned x, unsigned y, unsigned z, unsigned long index);
	unsigned long getIndex(unsigned x, unsigned y, unsigned z, unsigned long position);	
	unsigned long getSize(unsigned x, unsigned y, unsigned z);
	// empty all inner containers, but don't deallocate
	void clear();
private:
	//unsigned getBox(unsigned x, unsigned y, unsigned z);
};

inline void Neighborlist::addIndex(
	unsigned x,
	unsigned y,
	unsigned z,
	unsigned long index)
{
	//this->neighbors[x][y][z].push_back(index);
	this->neighbors[getBox(x,y,z)].push_back(index);	
	this->sizes[getBox(x,y,z)] += 1;
}

inline unsigned long Neighborlist::getIndex(
	unsigned x,
	unsigned y,
	unsigned z,
	unsigned long position)
{
	//return this->neighbors[x][y][z].operator[](position);
	return this->neighbors[getBox(x,y,z)].operator[](position);
}

inline unsigned long Neighborlist::getSize(
	unsigned x,
	unsigned y,
	unsigned z)
{
	//return this->neighbors[x][y][z].size();
	//return this->neighbors[getBox(x,y,z)].size();
	return this->sizes[getBox(x,y,z)];
}

inline void Neighborlist::clear()
{
	/*
	for (unsigned long x=0; x<this->numberBoxes; x++) {
		for (unsigned long y=0; y<this->numberBoxes; y++) {
			for (unsigned long z=0; z<this->numberBoxes; z++) {
				this->neighbors[x][y][z].clear();
			}
		}
	}
	*/
	for (unsigned i=0; i<this->totalNumberBoxes; i++) {
		this->neighbors[i].clear();
	}
	std::fill(sizes.begin(),sizes.end(),0);
}

/*
inline unsigned Neighborlist::getBox(
	unsigned x,
	unsigned y,
	unsigned z)
{
	return x + this->numberBoxes * y + this->numberBoxes * this->numberBoxes * z;
}
*/
#endif // __NEIGHBORLIST_H_INCLUDED__