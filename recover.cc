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
#include "recover.h"

NTL_CLIENT

// Return a new string consisting of s followed by the
// bytes_per_word-byte representation of wz
template<>
string GSDecoder_ZZ_p::append(const string &s, const ZZ_p &wz,
	unsigned int bytes_per_word)
{
    unsigned char w[bytes_per_word];
    BytesFromZZ(w, rep(wz), bytes_per_word);
    string r = s;
    r.append((char *)w, bytes_per_word);
    return r;
}

// Return a new string consisting of s followed by the
// bytes_per_word-byte representation of wz
template<>
string GSDecoder_GF2E::append(const string &s, const GF2E &wz,
	unsigned int bytes_per_word)
{
    unsigned char w[bytes_per_word];
    BytesFromGF2X(w, rep(wz), bytes_per_word);
    string r = s;
    r.append((char *)w, bytes_per_word);
    return r;
}

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

vector<PercyResult> EasyRecover_GF28(unsigned short t,
	unsigned short h, const vector<PercyResult> &H,
	const GF28_Element *values, const GF28_Element *indices)
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
	GF28_Element I_indices[t+1], I_values[t+1];
	unsigned short i = 0;
	for (Iiter = I.begin(); Iiter != I.end(); ++i, ++Iiter) {
	    I_indices[i] = indices[*Iiter];
	    I_values[i] = values[*Iiter];
	}
	GF28_Element wz;
	wz = interpolate_GF28(I_indices, I_values, t+1, 0);

	// Count the number of points in G that agree, and that
	// disagree, with phi
	unsigned short numagree = 0;
	unsigned short numdisagree = 0;
	vector<unsigned short> vecagree;
	vector<unsigned short>::const_iterator Giter;
	for (Giter = Hiter->G.begin(); Giter != Hiter->G.end(); ++Giter) {
	    if (interpolate_GF28(I_indices, I_values, t+1, indices[*Giter])
			== values[*Giter]) {
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
	    string news(Hiter->sigma);
	    news.append((char*)&wz, 1);
	    PercyResult n(vecagree, news);
	    Hprime.push_back(n);
	} else {
	    // This either isn't the right polynomial, or there may be
	    // more than one.  Abort.
	    Hprime.clear();
	    return Hprime;
	}
    }
    return Hprime;
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
	    PercyResult n(vecagree, GSDecoder_ZZ_p::append(Hiter->sigma,
			wz, bytes_per_word));
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

#if 0
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

// Construct a bivariate polynomial P(x,y) such that, for any polynomial
// f(x) of degree at most k that agrees with at least t of the given
// points, (y-f(x)) is a factor of P(x,y).  This version is the naive
// algorithm from Venkatesan Guruswami and Madhu Sudan, "Improved
// Decoding of Reed-Solomon and Algebraic-Geometry Codes".
static ZZ_pXY interpolate_naive(unsigned int k, unsigned int t,
	const vector<unsigned short> &goodservers,
	const vec_ZZ_p& indices, const vec_ZZ_p& shares)
{
    unsigned int n = goodservers.size();

    // Compute the r and l parameters
    unsigned int r = 1 + (unsigned int)(floor( (k*n + sqrt(k*k*n*n+4*(t*t-k*n)))/(2*(t*t-k*n)) ));
    unsigned int l = r*t - 1;

    std::cerr << "Constructing (1," << k << ")-degree " << (l/k)*k << " polynomial...\n";

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

    Ccache_t Ccache;

    // This is the part you need to read the paper for
    unsigned int Arow = 0;
    for (i=0; i<n; ++i) {
	for (j1=0; j1<r; ++j1) {
	    // std::cerr << i << ", " << j1 << " / " << n << ", " << r << "\n";
	    for (j2=0; j2 < r - j1; ++j2) {
		for (j2p = j2; j2p <= l/k; ++j2p) {
		    for (j1p = j1; j1p <= l - k*j2p; ++j1p) {
			pairint j1pj2p(j1p, j2p);
			A[Arow][Qmap[j1pj2p]] = C(Ccache, j1p,j1) *
			    C(Ccache, j2p,j2) *
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
    ZZ_pXY P;
    P.SetMaxLength(l/k + 1);   // The y-degree of P

    for(unsigned int i=0; i<c; ++i) {
	if (soln[i] != 0) {
	    SetCoeff(P, Qinv[i].second, Qinv[i].first, soln[i]);
	}
    }

    return P;
}
#endif

#ifdef TEST_RECOVER

#include <sstream>
main()
{
    ZZ modulus, one;
    modulus = 10007;
    one = 1;
    ZZ_p::init(modulus);

    stringstream ss(stringstream::in | stringstream::out);

    vector<unsigned short> goodservers;
    goodservers.push_back(1);
    goodservers.push_back(2);
    goodservers.push_back(3);
    goodservers.push_back(4);
    goodservers.push_back(5);
    goodservers.push_back(6);
    vec_ZZ_p indices, shares;
    ss << "[ 1 2 3 4 5 6 ]";
    ss >> indices;
    ss << "[ 4 3 7507 1 2504 7 ]";
    ss >> shares;
    vector<RecoveryPoly> ret = findpolys(3, 5, goodservers, indices, shares,
	    one, modulus);
    for (vector<RecoveryPoly>::const_iterator iter = ret.begin(); iter != ret.end(); ++iter) {
	cout << "{ " << iter->phi << ", [ ";
	for (vector<unsigned short>::const_iterator gi = iter->G.begin(); gi != iter->G.end(); ++gi) {
	    cout << *gi << " ";
	}
	cout << "] }\n";
    }
}

#endif
