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

#include <fstream>
#include <vector>
#include <map>
#include <math.h>
#include <mat_ZZ_p.h>
#include <ZZ_pX.h>
#include "findroots.h"
#include "recover.h"

NTL_CLIENT

// Append k random elements of G[offset..end] to I
static void random_subset(const vector<unsigned short> &G,
	vector<unsigned short> &I, unsigned short k, unsigned short offset = 0)
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

// Return a new string consisting of s followed by the
// bytes_per_word-byte representation of wz
static string append(const string &s, const ZZ_p &wz,
	unsigned int bytes_per_word)
{
    unsigned char w[bytes_per_word];
    BytesFromZZ(w, rep(wz), bytes_per_word);
    string r = s;
    r.append((char *)w, bytes_per_word);
    return r;
}

vector<PercyResult> EasyRecover(unsigned int bytes_per_word, unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const vec_ZZ_p &values, const vec_ZZ_p &indices)
{
    vector<PercyResult> Hprime;
    vector<PercyResult>::const_iterator Hiter;

    for (Hiter = H.begin(); Hiter != H.end(); ++Hiter) {
	// Pick a random subset I of G, of size t+1
	vector<unsigned short> I;
	vector<unsigned short>::const_iterator Iiter;
	random_subset(Hiter->G, I, t+1);

	// Use Lagrange interpolation to find the unique polynomial phi
	// of degree t which matches the points indexed by I
	vec_ZZ_p I_indices, I_values;
	I_indices.SetLength(t+1);
	I_values.SetLength(t+1);
	unsigned short i = 0;
	for (Iiter = I.begin(); Iiter != I.end(); ++i, ++Iiter) {
	    I_indices[i] = indices[*Iiter];
	    I_values[i] = values[*Iiter];
	}
	ZZ_pX phi;
	interpolate(phi, I_indices, I_values);

	// Find the secret determined by phi
	ZZ_p wz;
	eval(wz, phi, ZZ_p::zero());

	// Count the number of points in G that agree, and that
	// disagree, with phi
	unsigned short numagree = 0;
	unsigned short numdisagree = 0;
	vector<unsigned short> vecagree;
	vector<unsigned short>::const_iterator Giter;
	for (Giter = Hiter->G.begin(); Giter != Hiter->G.end(); ++Giter) {
	    ZZ_p phival;
	    eval(phival, phi, indices[*Giter]);
	    if (phival == values[*Giter]) {
		++numagree;
		vecagree.push_back(*Giter);
	    } else {
		++numdisagree;
	    }
	}

	// If at least h agreed, and less than h-t disagreed, then phi
	// can be the *only* polynomial of degree t matching at least
	// h points.
	if (numagree >= h && numdisagree < h-t) {
	    PercyResult n(vecagree, append(Hiter->sigma, wz, bytes_per_word));
	    Hprime.push_back(n);
	} else {
	    // This either isn't the right polynomial, or there may be
	    // more than one.  Abort, and we'll use HardRecover.
	    Hprime.clear();
	    return Hprime;
	}
    }
    return Hprime;
}

// Compute the intersection of the two *strictly sorted* vectors v1 and v2
static vector<unsigned short> intersect(const vector<unsigned short> &v1,
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

static ZZ_p C(unsigned int n, unsigned int k)
{
    ZZ_p num, dem;
    num = 1;
    dem = 1;
    unsigned int i;

    for (i=0;i<k;++i) {
	num *= (n-i);
	dem *= (k-i);
    }

    return num/dem;
}

// Return the index of the first non-zero entry in this row, or -1 if
// it's all zeros.
static long first_non_zero(const vec_ZZ_p& R)
{
    long len = R.length();
    for (long i=0; i<len; ++i) {
	if (!IsZero(R[i])) {
	    return i;
	}
    }
    return -1;
}

static vec_ZZ_p solvemat(mat_ZZ_p& M) {
    // Do Gaussian elimination on M.
    gauss(M);

    // Build a non-zero solution.  First find the last row that isn't
    // all zeros.
    long numrows = M.NumRows();
    long numcols = M[0].length();
    long lastnonzerorow = -1;
    for (long i=numrows-1; i>=0; --i) {
	if (!IsZero(M[i])) {
	    lastnonzerorow = i;
	    break;
	}
    }

    // Initialize the solution vector to all 0's
    vec_ZZ_p soln;
    soln.SetLength(numcols);

    // We want a non-zero solution, so find the leftmost column that
    // does not contain a leading 1.
    long lastlead = -1;
    int found = 0;
    for (long i=0; i<numrows; ++i) {
	long leadingidx = first_non_zero(M[i]);
	if (leadingidx < 0 || leadingidx > lastlead + 1) {
	    soln[lastlead+1] = 1;
	    found = 1;
	    break;
	}
	lastlead = leadingidx;
    }
    if (!found && lastlead < numcols - 1) {
	soln[lastlead+1] = 1;
    }

    // Do the back substitution
    for (long i=lastnonzerorow; i>=0; --i) {
	// Find the first non-zero entry
	long leadingidx = first_non_zero(M[i]);
	ZZ_p newval;
	InnerProduct(newval, M[i], soln);
	soln[leadingidx] = -newval / M[i][leadingidx];
    }

    return soln;
}

struct RecoveryPoly {
    RecoveryPoly(vector<unsigned short> G, ZZ_pX phi) : G(G), phi(phi) {}
    vector<unsigned short> G;
    ZZ_pX phi;
};

// Find all polynomials of degree at most k that agree with the given
// (index,share) pairs at at least t of the n points.  The notation is
// from Venkatesan Guruswami and Madhu Sudan, "Improved Decoding of
// Reed-Solomon and Algebraic-Geometry Codes".
static vector<RecoveryPoly> findpolys(unsigned int k, unsigned int t,
	const vector<unsigned short> &goodservers,
	const vec_ZZ_p& indices, const vec_ZZ_p& shares, ZZ p1, ZZ p2)
{
    unsigned int n = goodservers.size();

    // Compute the r and l parameters
    unsigned int r = 1 + (unsigned int)(floor( (k*n + sqrt(k*k*n*n+4*(t*t-k*n)))/(2*(t*t-k*n)) ));
    unsigned int l = r*t - 1;

    typedef pair<unsigned int, unsigned int> pairint;
    map<pairint, unsigned int> Qmap;
    map<unsigned int, pairint> Qinv;

    // Make a mapping from <xdegree, ydegree> to a single index.
    // Also make the inverse mapping.
    unsigned int c = 0;
    unsigned int i,j,j1,j2,j1p,j2p;
    for(j=0; j <= l/k; ++j) {
	for (i=0; i <= l-j*k; ++i) {
	    pairint ij(i,j);
	    Qmap[ij] = c;
	    Qinv[c] = ij;
	    ++c;
	}
    }

    // We're going to build a large matrix A and find a non-zero
    // solution to A.v = 0.
    mat_ZZ_p A;

    A.SetDims(n * r * (r+1) / 2, c);

    std::cerr << "Generating " << n * r * (r+1) / 2 << " x " << c << " matrix...\n";

    // This is the part you need to read the paper for
    unsigned int Arow = 0;
    for (i=0; i<n; ++i) {
	for (j1=0; j1<r; ++j1) {
	    // std::cerr << i << ", " << j1 << " / " << n << ", " << r << "\n";
	    for (j2=0; j2 < r - j1; ++j2) {
		for (j2p = j2; j2p <= l/k; ++j2p) {
		    for (j1p = j1; j1p <= l - k*j2p; ++j1p) {
			pairint j1pj2p(j1p, j2p);
			A[Arow][Qmap[j1pj2p]] = C(j1p,j1) * C(j2p,j2) *
			    power(indices[goodservers[i]], j1p-j1) *
			    power(shares[goodservers[i]], j2p-j2);
		    }
		}
		++Arow;
	    }
	}
    }

    std::cerr << "Solving matrix...\n";

    vec_ZZ_p soln = solvemat(A);

    // The soln vector now consists of coefficients for a bivariate
    // polynomial P(x,y).
    vec_ZZ_pX P;  // This really should be ZZ_pXY
    P.SetLength(l/k + 1);   // The y-degree of P

    for(unsigned int i=0; i<c; ++i) {
	if (soln[i] != 0) {
	    SetCoeff(P[Qinv[i].second], Qinv[i].first, soln[i]);
	}
    }

    std::cerr << "Factoring resulting polynomial...\n";

    // It turns out that any polynomial phi(x) that we should be
    // returning (since phi is of degree at most k and agrees with the
    // input data on at least t points) is such that (y - phi(x)) is a
    // factor of P(x,y).  So we first find all roots for y of P(x,y)
    // which are polynomials of degree at most k.
    vec_ZZ_pX roots = findroots(P, k, p1, p2);

    // For each of these roots, check how many input points it agrees
    // with.  If it's at least t, add it to the list of polys to return.
    vector<RecoveryPoly> polys;
    unsigned int numroots = roots.length();
    for (unsigned int i=0; i<numroots; ++i) {
	if (deg(roots[i]) > (long)k) continue;
	vector<unsigned short>::const_iterator gooditer;
	vector<unsigned short> vecagree;
	unsigned short numagree = 0;
	for (gooditer = goodservers.begin(); gooditer != goodservers.end();
		++gooditer) {
	    ZZ_p phival;
	    eval(phival, roots[i], indices[*gooditer]);
	    if (phival == shares[*gooditer]) {
		++numagree;
		vecagree.push_back(*gooditer);
	    }
	}
	if (numagree >= t) {
	    RecoveryPoly n(vecagree, roots[i]);
	    polys.push_back(n);
	}
    }

    return polys;
}

vector<PercyResult> HardRecover(unsigned int bytes_per_word, unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const vector<unsigned short> &goodservers,
	const vec_ZZ_p &values, const vec_ZZ_p &indices, ZZ p1, ZZ p2)
{
    vector<PercyResult> Hprime;

    // Find all polynomials of degree at most t that match at least h of
    // the (index, value) pairs, mod p1*p2  (p1 might be 1).
    vector<RecoveryPoly> polys = findpolys(t, h, goodservers, indices, values,
	    p1, p2);
    vector<RecoveryPoly>::const_iterator Piter;

    for (Piter = polys.begin(); Piter != polys.end(); ++Piter) {
	// Find the secret determined by this poly
	ZZ_p wz;
	eval(wz, Piter->phi, ZZ_p::zero());

	// Find any elements (G,sigma) of H with |Piter->G \cap G| \ge h
	vector<PercyResult>::const_iterator Hiter;
	for (Hiter = H.begin(); Hiter != H.end(); ++Hiter) {
	    vector<unsigned short> intersection =
		intersect(Piter->G, Hiter->G);
	    if (intersection.size() >= h) {
		PercyResult n(intersection,
			append(Hiter->sigma, wz, bytes_per_word));
		Hprime.push_back(n);
	    }
	}

    }

    return Hprime;
}
