//  Percy++
//  Copyright 2007 Ian Goldberg <iang@cs.uwaterloo.ca>
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of version 2 of the GNU General Public License as
//  published by the Free Software Foundation.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  There is a copy of the GNU General Public License in the COPYING file
//  packaged with this plugin; if you cannot find it, write to the Free
//  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
//  02111-1307  USA

// This file is an implementation of the Roth-Ruckenstein algorithm
// for finding roots of bivariate polynomials, as explained in section
// IX of:  R. J. McEliece. The Guruswami-Sudan Decoding Algorithm for
// Reed-Solomon Codes. IPN Progress Report 42-153, May 15, 2003.
// http://citeseer.ist.psu.edu/mceliece03guruswamisudan.html

#include <sstream>
#include <vector>
#include "rr_roots.h"

NTL_CLIENT

#ifdef TEST_RR

// Test the algorithm on the Example 15 from the McElice paper.
int main()
{
    ZZ modulus;
    modulus = 19;
    ZZ p1;
    p1 = 1;
    ZZ_p::init(modulus);

    stringstream ss(stringstream::in | stringstream::out);
    ZZ_pXY P;

    ss << "[[4 12 5 11 8 13] [14 14 9 16 8] [14 13 1] [2 11 1] [17]] ";
    ss >> P;

    vec_ZZ_pX roots = rr_findroots(P, 1);
    cout << "\nRoots found = " << roots << "\n";

    return 0;
}
#endif
