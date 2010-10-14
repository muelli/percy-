#ifndef __SUBSET_H__
#define __SUBSET_H__

// Operations on subsets of the set of servers.  These subsets are
// stored as sorted vectors of unsigned short.

#include <vector>

using namespace std;

// Append k random elements of G[offset..end] to I
void random_subset(const vector<unsigned short> &G,
    vector<unsigned short> &I, unsigned short k, unsigned short offset = 0);

// Compute the intersection of the two *strictly sorted* vectors v1 and v2
vector<unsigned short> intersect(const vector<unsigned short> &v1,
	const vector<unsigned short> &v2);

#endif
