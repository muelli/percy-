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
//
//  The ZZ_pXY.cc and ZZ_pXY.h files are adapted from files in version
//  5.4 of the NTL package by Victor Shoup <victor@shoup.net>, with
//  modifications by M. Jason Hinek <mjhinek@alumni.uwaterloo.ca> and
//  Ian Goldberg <iang@cs.uwaterloo.ca>.

#ifndef NTL_ZZ_pXY__H
#define NTL_ZZ_pXY__H

#include <ZZ_pX.h>
#include <ZZ_pXFactoring.h>

NTL_OPEN_NNS




/************************************************************

                         ZZ_pXY

The class ZZ_pXY implements bivariate polynomial arithmetic modulo p.
Polynomials are represented as vec_ZZ_pX's.
If f is a ZZ_pXY, then f.rep is a vec_ZZ_pX.
The zero polynomial is represented as a zero length vector.
Otherwise. f.rep[0] is the constant-term, and f.rep[f.rep.length()-1]
is the leading coefficient, which is always non-zero.
The member f.rep is public, so the vector representation is fully
accessible.
Use the member function normalize() to strip leading zeros.

**************************************************************/


class ZZ_pXY {

public:

typedef vec_ZZ_pX VectorBaseType; 


vec_ZZ_pX rep;


/***************************************************************

          Constructors, Destructors, and Assignment

****************************************************************/


ZZ_pXY()
//  initial value 0

   { }


ZZ_pXY(INIT_SIZE_TYPE, long n) { rep.SetMaxLength(n); }

ZZ_pXY(const ZZ_pXY& a) : rep(a.rep) { }
// initial value is a


ZZ_pXY& operator=(const ZZ_pXY& a) 
   { rep = a.rep; return *this; }

~ZZ_pXY() { }

void normalize();
// strip leading zeros

void SetMaxLength(long n) 
// pre-allocate space for n coefficients.
// Value is unchanged

   { rep.SetMaxLength(n); }


void kill() 
// free space held by this polynomial.  Value becomes 0.

   { rep.kill(); }

static const ZZ_pXY& zero();


ZZ_pXY(ZZ_pXY& x, INIT_TRANS_TYPE) : rep(x.rep, INIT_TRANS) { }

inline ZZ_pXY(long i, const ZZ_pX& c);
inline ZZ_pXY(long i, const ZZ_p& c);
inline ZZ_pXY(long i, long c);

ZZ_pXY& operator=(long a);
ZZ_pXY& operator=(const ZZ_p& a);
ZZ_pXY& operator=(const ZZ_pX& a);


};



/********************************************************************

                           input and output

I/O format:

   [a_0 a_1 ... a_n],

represents the polynomial a_0 + a_1*Y + ... + a_n*Y^n.

On output, all coefficients will be polynomials of the type ZZ_pX.
amd a_n not zero (the zero polynomial is [ ]).
On input, the coefficients are arbitrary polynomials which are
then reduced modulo p, and leading zeros stripped.

*********************************************************************/


NTL_SNS istream& operator>>(NTL_SNS istream& s, ZZ_pXY& x);
NTL_SNS ostream& operator<<(NTL_SNS ostream& s, const ZZ_pXY& a);




/**********************************************************

                   Some utility routines

***********************************************************/



inline long deg(const ZZ_pXY& a) { return a.rep.length() - 1; }
// degree of a polynomial in Z_p[X][Y] (ie, degree wrt y).
// note that the zero polynomial has degree -1.

inline long degY(const ZZ_pXY& a) { return a.rep.length() - 1; }
// degree of a polynomial in Z_p[X][Y] (ie, degree wrt y) 
// note that the zero polynomial has degree -1.

long degX(const ZZ_pXY& a) ;
// degree of a polynomial in Z_p[Y][X] instead of Z_p[X][Y] (ie, degree wrt x)
// note that the zero polynomial has degree -1.



const ZZ_pX& coeff(const ZZ_pXY& a, long i);
// zero if i not in range

void GetCoeff(ZZ_pX& x, const ZZ_pXY& a, long i);
// x = a[i], or zero if i not in range

const ZZ_pX& LeadCoeff(const ZZ_pXY& a);
// zero if a == 0

const ZZ_pX& ConstTerm(const ZZ_pXY& a);
// zero if a == 0

void SetCoeff(ZZ_pXY& x, long i, long j, const ZZ_p& a);
// x[i][j] = a, error is raised if i < 0

void SetCoeff(ZZ_pXY& x, long i, const ZZ_pX& a);
// x[i] = a, error is raised if i < 0

void SetCoeff(ZZ_pXY& x, long i, long a);

void SetCoeff(ZZ_pXY& x, long i);
// x[i] = 1, error is raised if i < 0

inline ZZ_pXY::ZZ_pXY(long i, const ZZ_pX& a)
   { SetCoeff(*this, i, a); } 

inline ZZ_pXY::ZZ_pXY(long i, long a)
   { SetCoeff(*this, i, a); } 

void SetY(ZZ_pXY& y);
// y is set to the monomial Y

long IsY(const ZZ_pXY& a);
// test if a = Y

inline void clear(ZZ_pXY& x) 
// x = 0

   { x.rep.SetLength(0); }

inline void set(ZZ_pXY& x)
// x = 1

   { x.rep.SetLength(1); set(x.rep[0]); }

inline void swap(ZZ_pXY& x, ZZ_pXY& y)
// swap x & y (only pointers are swapped)

   { swap(x.rep, y.rep); }

void random(ZZ_pXY& x, long n);
inline ZZ_pXY random_ZZ_pXY(long n)
   { ZZ_pXY x; random(x, n); NTL_OPT_RETURN(ZZ_pXY, x); }
// generate a random polynomial of degree < n 

void trunc(ZZ_pXY& x, const ZZ_pXY& a, long m);
// x = a % Y^m

inline ZZ_pXY trunc(const ZZ_pXY& a, long m)
   { ZZ_pXY x; trunc(x, a, m); NTL_OPT_RETURN(ZZ_pXY, x); }

void RightShift(ZZ_pXY& x, const ZZ_pXY& a, long n);
// x = a/Y^n

inline ZZ_pXY RightShift(const ZZ_pXY& a, long n)
   { ZZ_pXY x; RightShift(x, a, n); NTL_OPT_RETURN(ZZ_pXY, x); }

void LeftShift(ZZ_pXY& x, const ZZ_pXY& a, long n);
// x = a*Y^n

inline ZZ_pXY LeftShift(const ZZ_pXY& a, long n)
   { ZZ_pXY x; LeftShift(x, a, n); NTL_OPT_RETURN(ZZ_pXY, x); }

#ifndef NTL_TRANSITION

inline ZZ_pXY operator>>(const ZZ_pXY& a, long n)
   { ZZ_pXY x; RightShift(x, a, n); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator<<(const ZZ_pXY& a, long n)
   { ZZ_pXY x; LeftShift(x, a, n); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY& operator<<=(ZZ_pXY& x, long n)
   { LeftShift(x, x, n); return x; }

inline ZZ_pXY& operator>>=(ZZ_pXY& x, long n)
   { RightShift(x, x, n); return x; }

#endif



void diff(ZZ_pXY& x, const ZZ_pXY& a);
// x = derivative of a wrt Y

inline ZZ_pXY diff(const ZZ_pXY& a)
   { ZZ_pXY x; diff(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }


void MakeMonic(ZZ_pXY& x);

void reverse(ZZ_pXY& c, const ZZ_pXY& a, long hi);

inline ZZ_pXY reverse(const ZZ_pXY& a, long hi)
   { ZZ_pXY x; reverse(x, a, hi); NTL_OPT_RETURN(ZZ_pXY, x); }

inline void reverse(ZZ_pXY& c, const ZZ_pXY& a)
{  reverse(c, a, deg(a)); }

inline ZZ_pXY reverse(const ZZ_pXY& a)
   { ZZ_pXY x; reverse(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }


/*******************************************************************

                        conversion routines

********************************************************************/



void conv(ZZ_pXY& x, long a);
void conv(ZZ_pXY& x, const ZZ& a);
void conv(ZZ_pXY& x, const ZZ_p& a);
void conv(ZZ_pXY& x, const ZZ_pX& a);
void conv(ZZ_pXY& x, const vec_ZZ_pX& a);

inline ZZ_pXY to_ZZ_pXY(long a)
   { ZZ_pXY x; conv(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY to_ZZ_pXY(const ZZ& a)
   { ZZ_pXY x; conv(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY to_ZZ_pXY(const ZZ_p& a)
   { ZZ_pXY x; conv(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY to_ZZ_pXY(const ZZ_pX& a)
   { ZZ_pXY x; conv(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY to_ZZ_pXY(const vec_ZZ_pX& a)
   { ZZ_pXY x; conv(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }



/*************************************************************

                        Comparison

**************************************************************/

long IsZero(const ZZ_pXY& a); 

long IsOne(const ZZ_pXY& a);

inline long operator==(const ZZ_pXY& a, const ZZ_pXY& b)
{
   return a.rep == b.rep;
}

inline long operator!=(const ZZ_pXY& a, const ZZ_pXY& b)
{
   return !(a == b);
}

long operator==(const ZZ_pXY& a, long b);
long operator==(const ZZ_pXY& a, const ZZ_p& b);
long operator==(const ZZ_pXY& a, const ZZ_pX& b);

inline long operator==(long a, const ZZ_pXY& b) { return b == a; }
inline long operator==(const ZZ_p& a, const ZZ_pXY& b) { return b == a; }
inline long operator==(const ZZ_pX& a, const ZZ_pXY& b) { return b == a; }

inline long operator!=(const ZZ_pXY& a, long b) { return !(a == b); }
inline long operator!=(const ZZ_pXY& a, const ZZ_p& b) { return !(a == b); }
inline long operator!=(const ZZ_pXY& a, const ZZ_pX& b) { return !(a == b); }

inline long operator!=(long a, const ZZ_pXY& b) { return !(a == b); }
inline long operator!=(const ZZ_p& a, const ZZ_pXY& b) { return !(a == b); }
inline long operator!=(const ZZ_pX& a, const ZZ_pXY& b) { return !(a == b); }


/***************************************************************

                         Addition

****************************************************************/

//
// x = a + b
//
void add(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b);
void add(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pX& b);   
void add(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_p& b);
void add(ZZ_pXY& x, const ZZ_pXY& a, long b);

inline void add(ZZ_pXY& x, const ZZ_pX& a, const ZZ_pXY& b) { add(x, b, a); }
inline void add(ZZ_pXY& x, const ZZ_p&  a, const ZZ_pXY& b) { add(x, b, a); }
inline void add(ZZ_pXY& x, long         a, const ZZ_pXY& b) { add(x, b, a); }


//
// x = a - b
//
void sub(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b);
void sub(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pX& b); 
void sub(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_p& b); 
void sub(ZZ_pXY& x, const ZZ_pXY& a, long b);

void sub(ZZ_pXY& x, const ZZ_pX& a, const ZZ_pXY& b);
void sub(ZZ_pXY& x, const ZZ_p&  a, const ZZ_pXY& b);
void sub(ZZ_pXY& x, long         a, const ZZ_pXY& b);


//
// x = -a
//
void negate(ZZ_pXY& x, const ZZ_pXY& a);


//
// operator versions for add, sub and negate
//
inline ZZ_pXY operator+(const ZZ_pXY& a, const ZZ_pXY& b)
   { ZZ_pXY x; add(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator+(const ZZ_pXY& a, const ZZ_pX& b)
   { ZZ_pXY x; add(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator+(const ZZ_pXY& a, const ZZ_p& b)
   { ZZ_pXY x; add(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator+(const ZZ_pXY& a, long b)
   { ZZ_pXY x; add(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator+(const ZZ_pX& a, const ZZ_pXY& b)
   { ZZ_pXY x; add(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator+(const ZZ_p&  a, const ZZ_pXY& b)
   { ZZ_pXY x; add(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator+(long a, const ZZ_pXY& b)
   { ZZ_pXY x; add(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }


inline ZZ_pXY operator-(const ZZ_pXY& a, const ZZ_pXY& b)
   { ZZ_pXY x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator-(const ZZ_pXY& a, const ZZ_pX& b)
   { ZZ_pXY x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator-(const ZZ_pXY& a, const ZZ_p& b)
   { ZZ_pXY x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator-(const ZZ_pXY& a, long b)
   { ZZ_pXY x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }


inline ZZ_pXY operator-(const ZZ_pX& a, const ZZ_pXY& b)
   { ZZ_pXY x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator-(const ZZ_p& a, const ZZ_pXY& b)
   { ZZ_pXY x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator-(long a, const ZZ_pXY& b)
   { ZZ_pXY x; sub(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }


inline ZZ_pXY& operator+=(ZZ_pXY& x, const ZZ_pXY& b)
   { add(x, x, b); return x; }

inline ZZ_pXY& operator+=(ZZ_pXY& x, const ZZ_pX& b)
   { add(x, x, b); return x; }

inline ZZ_pXY& operator+=(ZZ_pXY& x, const ZZ_p& b)
   { add(x, x, b); return x; }

inline ZZ_pXY& operator+=(ZZ_pXY& x, long b)
   { add(x, x, b); return x; }


inline ZZ_pXY& operator-=(ZZ_pXY& x, const ZZ_pXY& b)
   { sub(x, x, b); return x; }

inline ZZ_pXY& operator-=(ZZ_pXY& x, const ZZ_pX& b)
   { sub(x, x, b); return x; }

inline ZZ_pXY& operator-=(ZZ_pXY& x, const ZZ_p& b)
   { sub(x, x, b); return x; }

inline ZZ_pXY& operator-=(ZZ_pXY& x, long b)
   { sub(x, x, b); return x; }


inline ZZ_pXY operator-(const ZZ_pXY& a) 
   { ZZ_pXY x; negate(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY& operator++(ZZ_pXY& x) 
   { add(x, x, 1); return x; }

inline void operator++(ZZ_pXY& x, int) 
   { add(x, x, 1); }

inline ZZ_pXY& operator--(ZZ_pXY& x) 
   { sub(x, x, 1); return x; }

inline void operator--(ZZ_pXY& x, int) 
   { sub(x, x, 1); }



/*************************************************************

                      Computing y-roots

**************************************************************/


ZZ_pXY mapit(const ZZ_pXY& Q, const ZZ_p& a); 
//
// returns Q(x,xy+a)/X^m, for largest m>=0 such that division is exact 
//

vec_ZZ_p findroots(const ZZ_pX &Q);
//
// returns all y-roots of Q(x,y)
// 


long minDegX(const ZZ_pX& P);
long minX(const ZZ_pXY& Q);
ZZ comb(long n, long k);

void multiply(ZZ_pXY& x, const ZZ_pXY& Q, long a);  // double & add
void multiply(ZZ_pX& x, const ZZ_pX& Q, long a);    // double & add
void multiply(ZZ_pX& x, const ZZ_pX& Q, ZZ a);       
/*****************************************************************

                        Multiplication

******************************************************************/

//
// x = a * b  (and associated operations)
//

void mul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b); // x = a * b
void mul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pX& b);
void mul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_p& b); 
void mul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ& b); 
void mul(ZZ_pXY& x, const ZZ_pXY& a, long b);

void sqr(ZZ_pXY& x, const ZZ_pXY& a); // x = a^2
void sqr(ZZ_pXY& x, const ZZ_pX& a);
void sqr(ZZ_pXY& x, const ZZ_p& a);
void sqr(ZZ_pXY& x, long a);

void power(ZZ_pXY& x, const ZZ_pXY& a, long e); // x = a^e (e >= 0)



// inline versions

inline ZZ_pXY sqr(const ZZ_pXY& a) 
   { ZZ_pXY x; sqr(x, a); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY power(const ZZ_pXY& a, long e)
   { ZZ_pXY x; power(x, a, e); NTL_OPT_RETURN(ZZ_pXY, x); }




// other variations (of input arguments)

inline void mul(ZZ_pXY& x, const ZZ_pX& a, const ZZ_pXY& b) 
   { mul(x, b, a); }

inline void mul(ZZ_pXY& x, const ZZ_p& a, const ZZ_pXY& b) 
   { mul(x, b, a); }

inline void mul(ZZ_pXY& x, long a, const ZZ_pXY& b) 
   { mul(x, b, a); }


// operator versions

inline ZZ_pXY operator*(const ZZ_pXY& a, const ZZ_pXY& b) 
   { ZZ_pXY x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator*(const ZZ_pXY& a, const ZZ_pX& b) 
   { ZZ_pXY x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator*(const ZZ_pXY& a, const ZZ_p& b)
   { ZZ_pXY x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator*(const ZZ_pXY& a, long b)
   { ZZ_pXY x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }


inline ZZ_pXY operator*(const ZZ_pX& a, const ZZ_pXY& b)
   { ZZ_pXY x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator*(const ZZ_p& a, const ZZ_pXY& b)
   { ZZ_pXY x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }

inline ZZ_pXY operator*(long a, const ZZ_pXY& b)
   { ZZ_pXY x; mul(x, a, b); NTL_OPT_RETURN(ZZ_pXY, x); }



inline ZZ_pXY& operator*=(ZZ_pXY& x, const ZZ_pXY& b)
   { mul(x, x, b); return x; }

inline ZZ_pXY& operator*=(ZZ_pXY& x, const ZZ_pX& b)
   { mul(x, x, b); return x; }

inline ZZ_pXY& operator*=(ZZ_pXY& x, const ZZ_p& b)
   { mul(x, x, b); return x; }

inline ZZ_pXY& operator*=(ZZ_pXY& x, long b)
   { mul(x, x, b); return x; }






//
// Classical (plain) n^2 algorithms
//

void PlainMul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b);

NTL_CLOSE_NNS


#endif
