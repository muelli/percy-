#ifndef __SUBSET_ITER_H__
#define __SUBSET_ITER_H__

#include <vector>

// Produce all subsets of a given vector<unsigned short> of a given size
//
// Usage:
//    vector<unsigned short> G;
//    ...
//    subset_iterator iter(G, k);
//    while(!iter.atend()) {
//	do_something_with(*iter);
//	++iter;
//    }

class subset_iterator {
private:
    std::vector<unsigned short> source;
    std::vector<unsigned short> indices;
    std::vector<unsigned short> output;
    unsigned short subset_size;
    bool finished;
public:
    subset_iterator(const std::vector<unsigned short>& source,
	    unsigned short subset_size);
    subset_iterator& operator++();
    const std::vector<unsigned short>& operator*() const;
    bool atend() const;
};

#endif
