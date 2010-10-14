#include <NTL/ZZ.h>  // For RandomBnd()
#include "subset.h"

NTL_CLIENT

// Operations on subsets of the set of servers.  These subsets are
// stored as sorted vectors of unsigned short.

// Append k random elements of G[offset..end] to I
void random_subset(const vector<unsigned short> &G,
	vector<unsigned short> &I, unsigned short k, unsigned short offset)
{
    if (k == 0) return;

    // How many elements are there to choose from?
    unsigned short n = G.size() - offset;

    if (n == 0) return;

    // With probability k/n, choose G[offset]
    if (RandomBnd(n) < k) {
	I.push_back(G[offset]);
	random_subset(G, I, k-1, offset+1);
    } else {
	random_subset(G, I, k, offset+1);
    }
}

// Compute the intersection of the two *strictly sorted* vectors v1 and v2
vector<unsigned short> intersect(const vector<unsigned short> &v1,
	const vector<unsigned short> &v2)
{
    vector<unsigned short> intersection;
    vector<unsigned short>::const_iterator v1p, v2p;
    v1p = v1.begin();
    v2p = v2.begin();
    while (v1p != v1.end() && v2p != v2.end()) {
	if (*v1p == *v2p) {
	    // This is an element of the intersection
	    intersection.push_back(*v1p);
	    ++v1p;
	    ++v2p;
	} else if (*v1p > *v2p) {
	    ++v2p;
	} else {  // *v1p < *v2p
	    ++v1p;
	}
    }
    return intersection;
}
