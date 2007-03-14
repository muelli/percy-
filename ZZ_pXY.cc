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

#include "ZZ_pXY.h"

NTL_START_IMPL

const ZZ_pXY& ZZ_pXY::zero()
{
   static ZZ_pXY z;
   return z;
}


ZZ_pXY& ZZ_pXY::operator=(long a)
{
   conv(*this, a);
   return *this;
}


ZZ_pXY& ZZ_pXY::operator=(const ZZ_p& a)
{
   conv(*this, a);
   return *this;
}

ZZ_pXY& ZZ_pXY::operator=(const ZZ_pX& a)
{
   conv(*this, a);
   return *this;
}


istream& operator>>(istream& s, ZZ_pXY& x)
{
   s >> x.rep;
   x.normalize();
   return s;
}

ostream& operator<<(ostream& s, const ZZ_pXY& a)
{
   return s << a.rep;
}


void ZZ_pXY::normalize()
{
   long n;
   ZZ_pX* p;

   // First normalize the coefficients
   n = rep.length();
   if (n == 0) return;
   p = rep.elts() + n;
   while (n > 0) {
       --p;
       p->normalize();
       --n;
   }
   
   n = rep.length();
   p = rep.elts() + n;
   while (n > 0 && IsZero(*--p)) {
      n--;
   }
   rep.SetLength(n);
}


long IsZero(const ZZ_pXY& a)
{
   return a.rep.length() == 0;
}


long IsOne(const ZZ_pXY& a)
{
    return a.rep.length() == 1 && IsOne(a.rep[0]);
}

void GetCoeff(ZZ_pX& x, const ZZ_pXY& a, long i)
{
   if (i < 0 || i > deg(a))
      clear(x);
   else
      x = a.rep[i];
}

void SetCoeff(ZZ_pXY& x, long i, long j, const ZZ_p& a)
{
   long k, m;

   if (i < 0) 
      Error("SetCoeff: negative index");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m) {
      /* careful: a may alias a coefficient of x */

#if 0
      long alloc = x.rep.allocated();

      if (alloc > 0 && i >= alloc) {
         ZZ_pXTemp aa_tmp;  ZZ_pX& aa = aa_tmp.val();
         aa = a;
         x.rep.SetLength(i+1);
         x.rep[i] = aa;
      }
      else {
#endif
         x.rep.SetLength(i+1);
         SetCoeff(x.rep[i], j, a);
#if 0
      }
#endif

      for (k = m+1; k < i; k++)
         clear(x.rep[k]);
   }
   else
	SetCoeff(x.rep[i], j, a);

   x.normalize();
}

void SetCoeff(ZZ_pXY& x, long i, const ZZ_pX& a)
{
   long j, m;

   if (i < 0) 
      Error("SetCoeff: negative index");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m) {
      /* careful: a may alias a coefficient of x */

#if 0
      long alloc = x.rep.allocated();

      if (alloc > 0 && i >= alloc) {
         ZZ_pXTemp aa_tmp;  ZZ_pX& aa = aa_tmp.val();
         aa = a;
         x.rep.SetLength(i+1);
         x.rep[i] = aa;
      }
      else {
#endif
         x.rep.SetLength(i+1);
         x.rep[i] = a;
#if 0
      }
#endif

      for (j = m+1; j < i; j++)
         clear(x.rep[j]);
   }
   else
      x.rep[i] = a;

   x.normalize();
}

void SetCoeff(ZZ_pXY& x, long i, long a)
{
   if (a == 1) 
      SetCoeff(x, i);
   else {
      ZZ_pX T;
      conv(T, a);
      SetCoeff(x, i, T);
   }
}

void SetCoeff(ZZ_pXY& x, long i)
{
   long j, m;

   if (i < 0) 
      Error("coefficient index out of range");

   if (NTL_OVERFLOW(i, 1, 0))
      Error("overflow in SetCoeff");

   m = deg(x);

   if (i > m) {
      x.rep.SetLength(i+1);
      for (j = m+1; j < i; j++)
         clear(x.rep[j]);
   }
   set(x.rep[i]);
   x.normalize();
}


void SetY(ZZ_pXY& y)
{
   clear(y);
   SetCoeff(y, 1);
}


long IsY(const ZZ_pXY& a)
{
   return deg(a) == 1 && IsOne(LeadCoeff(a)) && IsZero(ConstTerm(a));
}
      
      

const ZZ_pX& coeff(const ZZ_pXY& a, long i)
{
   if (i < 0 || i > deg(a))
      return ZZ_pX::zero();
   else
      return a.rep[i];
}


const ZZ_pX& LeadCoeff(const ZZ_pXY& a)
{
   if (IsZero(a))
      return ZZ_pX::zero();
   else
      return a.rep[deg(a)];
}

const ZZ_pX& ConstTerm(const ZZ_pXY& a)
{
   if (IsZero(a))
      return ZZ_pX::zero();
   else
      return a.rep[0];
}


void conv(ZZ_pXY& x, const ZZ_pX& a)
{
   if (IsZero(a))
      x.rep.SetLength(0);
   else {
      x.rep.SetLength(1);
      x.rep[0] = a;

      // note: if a aliases x.rep[i], i > 0, this code
      //       will still work, since is is assumed that
      //       SetLength(1) will not relocate or destroy x.rep[i]
   }
}

void conv(ZZ_pXY& x, long a)
{
   if (a == 0)
      clear(x);
   else if (a == 1)
      set(x);
   else {
      ZZ_pX T;
      conv(T, a);
      conv(x, T);
   }
}

void conv(ZZ_pXY& x, const ZZ& a)
{
   if (IsZero(a))
      clear(x);
   else {
      ZZ_pX T;
      conv(T, a);
      conv(x, T);
   }
}

void conv(ZZ_pXY& x, const ZZ_p& a)
{
   if (IsZero(a))
      clear(x);
   else {
      ZZ_pX T;
      conv(T, a);
      conv(x, T);
   }
}

void conv(ZZ_pXY& x, const vec_ZZ_pX& a)
{
   x.rep = a;
   x.normalize();
}


long degX(const ZZ_pXY& a)
{
  //
  // degree of polynomial in Z_p[Y][X] instead of Z_p[X][Y] 
  // 
  //
  
  long degY = a.rep.length() - 1;
  long degX = -1;

  if(degY > -1){
    degX = deg(a.rep[0]);
    for(long i = 1; i <= degY; i++){
      if( deg(a.rep[i]) > degX )
	degX = deg(a.rep[i]);
    }
    
  }
  
  return degX;
       
}



/***************************************************************

                         Addition

****************************************************************/

//
// x = a + b
//

void add(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b)
{
   long degYa = deg(a);
   long degYb = deg(b);
   long minab = min(degYa, degYb);
   long maxab = max(degYa, degYb);
   x.rep.SetLength(maxab+1);

   long i;
   const ZZ_pX *ap, *bp; 
   ZZ_pX* xp;

   for (i = minab+1, ap = a.rep.elts(), bp = b.rep.elts(), xp = x.rep.elts();
        i; i--, ap++, bp++, xp++)
      add(*xp, (*ap), (*bp));

   if (degYa > minab && &x != &a)
      for (i = degYa-minab; i; i--, xp++, ap++)
         *xp = *ap;
   else if (degYb > minab && &x != &b)
      for (i = degYb-minab; i; i--, xp++, bp++)
         *xp = *bp;
   else
      x.normalize();
}



void add(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pX& b)
{
  long n = a.rep.length(); 
  if (n == 0) {
    conv(x, b);
  }
  else if (&x == &a) {
    add(x.rep[0], a.rep[0], b);
    x.normalize();
  }
  else if (x.rep.MaxLength() == 0) {
    x = a;
    add(x.rep[0], a.rep[0], b);
    x.normalize();
  }
  else {
    // ugly...b could alias a coeff of x
    
    ZZ_pX *xp = x.rep.elts();
    add(xp[0], a.rep[0], b);
    x.rep.SetLength(n);
    xp = x.rep.elts();
    const ZZ_pX *ap = a.rep.elts();
    long i;
    for (i = 1; i < n; i++)
      xp[i] = ap[i];
    x.normalize();
  }
}

void add(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_p& b)
{
   if (a.rep.length() == 0) {
      conv(x, b);
   }
   else {
      if (&x != &a) x = a;
      add(x.rep[0], x.rep[0], b);
      x.normalize();
   }
}


void add(ZZ_pXY& x, const ZZ_pXY& a, long b)
{
   if (a.rep.length() == 0) {
      conv(x, b);
   }
   else {
      if (&x != &a) x = a;
      add(x.rep[0], x.rep[0], b);
      x.normalize();
   }
}


//
// x = a - b
//

void sub(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b)
{
   long da = deg(a);
   long db = deg(b);
   long minab = min(da, db);
   long maxab = max(da, db);
   x.rep.SetLength(maxab+1);

   long i;
   const ZZ_pX *ap, *bp; 
   ZZ_pX* xp;

   for (i = minab+1, ap = a.rep.elts(), bp = b.rep.elts(), xp = x.rep.elts();
        i; i--, ap++, bp++, xp++)
      sub(*xp, (*ap), (*bp));

   if (da > minab && &x != &a)
      for (i = da-minab; i; i--, xp++, ap++)
         *xp = *ap;
   else if (db > minab)
      for (i = db-minab; i; i--, xp++, bp++)
         negate(*xp, *bp);
   else
      x.normalize();

}

void sub(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pX& b)
{
   long n = a.rep.length();
   if (n == 0) {
      conv(x, b);
      negate(x, x);
   }
   else if (&x == &a) {
      sub(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else if (x.rep.MaxLength() == 0) {
      x = a;
      sub(x.rep[0], a.rep[0], b);
      x.normalize();
   }
   else {
      // ugly...b could alias a coeff of x

      ZZ_pX *xp = x.rep.elts();
      sub(xp[0], a.rep[0], b);
      x.rep.SetLength(n);
      xp = x.rep.elts();
      const ZZ_pX *ap = a.rep.elts();
      long i;
      for (i = 1; i < n; i++)
         xp[i] = ap[i];
      x.normalize();
   }
}

void sub(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_p b)
{
   if (IsZero(b)) {
      x = a;
      return;
   }

   if (a.rep.length() == 0) {
      x.rep.SetLength(1);
      x.rep[0] = b;
      negate(x.rep[0], x.rep[0]);
   }
   else {
      if (&x != &a) x = a;
      sub(x.rep[0], x.rep[0], b);
   }
   x.normalize();
}

void sub(ZZ_pXY& x, const ZZ_pXY& a, long b)
{
   if (b == 0) {
      x = a;
      return;
   }

   if (a.rep.length() == 0) {
      x.rep.SetLength(1);
      x.rep[0] = b;
      negate(x.rep[0], x.rep[0]);
   }
   else {
      if (&x != &a) x = a;
      sub(x.rep[0], x.rep[0], b);
   }
   x.normalize();
}


// other variants
void sub(ZZ_pXY& x, const ZZ_pX& a, const ZZ_pXY& b)
{
  ZZ_pXY temp;
  conv(temp,a);
  
  negate(x, b);
  add(x, x, temp);
}


void sub(ZZ_pXY& x, const ZZ_p& a, const ZZ_pXY& b)
{
  ZZ_pXY temp;
  conv(temp,a);
  
  negate(x, b);
  add(x, x, temp);
}

void sub(ZZ_pXY& x, long a, const ZZ_pXY& b)
{
  ZZ_pXY temp;
  conv(temp,a);
  
  negate(x, b);
  add(x, x, temp);
}



/* void sub(ZZ_pXY& x, long a, const ZZ_pXY& b) */
/* { */
/*    ZZ_pTemp TT; ZZ_p& T = TT.val();  */
/*    T = a; */

/*    negate(x, b); */
/*    add(x, x, T); */
/* } */


//
// x = -a
//

void negate(ZZ_pXY& x, const ZZ_pXY& a)
{
   long n = a.rep.length();
   x.rep.SetLength(n);

   const ZZ_pX* ap = a.rep.elts();
   ZZ_pX* xp = x.rep.elts();
   long i;

   for (i = n; i; i--, ap++, xp++)
      negate((*xp), (*ap));

}




/*****************************************************************

                        Computing y-roots

******************************************************************/

long minDegX(const ZZ_pX& P, long max){
  //
  // return m >= 0 such that X^m divides P(X)
  //

  if( IsZero(P) ) return max;
  
  long min = 0;
  while( IsZero(P.rep[min]) ){
    min++;
  }
  return min;
}


long minX(const ZZ_pXY& Q){
  //
  // return largest m >= 0 such that X^m divides Q(x,y)
  //
  
  if( IsZero(Q) ) return 0;

  long max = degX(Q);
  long minX = minDegX(Q.rep[0],max); 
  long degreeY = deg(Q);
  long currentMinX;

  for( long i = 1; i <= degreeY; i++){
    currentMinX = minDegX(Q.rep[i], max);
    if( currentMinX < minX)
      minX = currentMinX;
  }
  
  return minX;
}

ZZ comb(long n, long k){
  //
  // compute binomial coefficients
  //
  // comb(n,k) = n! / k!(n-k)!
  //

  if( k==0 ) return to_ZZ(1);
  if( k==1 ) return to_ZZ(n);

  ZZ ONE = to_ZZ(1);
  ZZ N = to_ZZ(n);
  ZZ K = to_ZZ(k);
  ZZ c = to_ZZ(1);
  if( k < n-k ){
    for(ZZ i = N; i >= N-K+ONE; i--){
      c = ((c*i) / (N-i+ONE));
    }
  }else{
    for(ZZ i = N; i >= K+ONE; i--){
      c = ((c*i) / (N-i+ONE));
    }
  } // if-then-else
  
  return c;
}


void multiply(ZZ_pXY& x, const ZZ_pXY& Q, long a)
{
  //
  // x = aQ (ie, Q added to itself a times)
  // 
  
  // do the easy stuff 
  if (a < 1) {
    Error("multiply: multiplier smaller than 1");
  }
  
  if ( IsOne(Q) || a == 1) {  // a1 = 1, 1Q = Q
    x = Q;
    return;
  }
  
  if (a == 2){
    add(x,Q,Q);
    return;
  }
  
  long dQ = deg(Q);
  
/*   if (dQ == 0) { */
/*     x = multiply(ConstTerm(Q), a); */
/*     return; */
/*   } */
  
  ZZ_pXY res;
  res.SetMaxLength(dQ + 1);
  res = 0;
  
  long k = NumBits(a);
  long i;
  
  for (i = k - 1; i >= 0; i--) {
    add(res,res,res); // double
    if (bit(a, i))
      add(res,res,Q); // add
  }
  
  x = res;
}

void multiply(ZZ_pX& x, const ZZ_pX& Q, long a)
{
  //
  // x = aQ (ie, Q added to itself a times)
  // 
  
  // do the easy stuff 
  if (a < 1) {
    Error("multiply: multiplier smaller than 1");
  }
  
  if ( IsOne(Q) || a == 1) {  // a1 = 1, 1Q = Q
    x = Q;
    return;
  }
  
  if (a == 2){
    add(x,Q,Q);
    return;
  }
  
  long dQ = deg(Q);
  
/*   if (dQ == 0) { */
/*     x = multiply(ConstTerm(Q), a); */
/*     return; */
/*   } */
  
  ZZ_pX res;
  res.SetMaxLength(dQ + 1);
  res = 0;
  
  long k = NumBits(a);
  long i;
  
  for (i = k - 1; i >= 0; i--) {
    add(res,res,res); // double
    if (bit(a, i))
      add(res,res,Q); // add
  }
  
  x = res;
}


void multiply(ZZ_pX& x, const ZZ_pX& Q, ZZ a)
{
  //
  // x = aQ (ie, Q added to itself a times)
  // 
  
  // do the easy stuff 
  if( a < to_ZZ(1) ) {
    Error("multiply: multiplier smaller than 1");
  }
  
  if( IsOne(Q) || IsOne(a) ) {  // a1 = 1, 1Q = Q
    x = Q;
    return;
  }
  
  if (a == to_ZZ(2)){
    add(x,Q,Q);
    return;
  }
  
  long dQ = deg(Q);
  
/*   if (dQ == 0) { */
/*     x = multiply(ConstTerm(Q), a); */
/*     return; */
/*   } */
  
  ZZ_pX res;
  res.SetMaxLength(dQ + 1);
  res = 0;
  
  long k = NumBits(a);
  long i;
  
  for (i = k - 1; i >= 0; i--) {
    add(res,res,res); // double
    if (bit(a, i))
      add(res,res,Q); // add
  }
  
  x = res;
}



ZZ_pX shiftX(const ZZ_pX& P, long n){
  // a =  X^n * P(X) 
  // n >= 0 for now

  long deg = P.rep.length()-1;
  
  ZZ_pX x;
  x.rep.SetLength( deg+1+n );

  for(long i = 0; i < n; i++){
    x.rep[i] = 0;
  }
  for(long i = 0; i <= deg; i++){
    x.rep[i+n] = P.rep[i];
  }
  return x;
}

ZZ_pXY backShiftX(const ZZ_pXY& P, long n){
  // returns  P(X) / X^n

  long deg = P.rep.length() - 1; 
  long degX;

  ZZ_pXY x;
  x.rep.SetLength( deg+1 );

  for(long i = 0; i <= deg; i++){ // work on each coefficient of Y
    if( IsZero(P.rep[i]) ){
      x.rep[i] = P.rep[i];
    }else{
      degX = (P.rep[i].rep.length()-1)-n;
      x.rep[i].rep.SetLength( degX + 1);
      for(long j = 0; j <= degX; j++){
	x.rep[i].rep[j] = P.rep[i].rep[j+n];
      }
    }
  }
  x.normalize();
  return x;
}


void mul(ZZ_pX& x, const ZZ_pX& a, const ZZ& b) 
{  
  ZZ_pTemp TT;  ZZ_p& T = TT.val();  
  conv(T, b);  
  mul(x, a, T);  
}  



ZZ_pXY mapit(const ZZ_pXY& Q, const ZZ_p& a){

  // 
  // if Q(x,y) has no y terms (ie deg(Q) < 1) 
  // there is nothing to do here
  //

  if( deg(Q) < 1 ) return Q;  


  //
  // Q(x,y) has some y terms, so here comes the work
  //

  ZZ_pXY temp;                    // we build up the new 
                                  // polynomial in temp

  long degY = deg(Q) ;
  temp.rep.SetLength( degY + 1 );      

  ZZ_pX coef1, coef2;             // we build up each coefficient with coef

  ZZ_p aa;

  //
  // first we let y --> xy + a
  //

  //std::cout << "--- mapit ---" << endl;

  for(long i = 0; i <= degY; i++){
    coef1 = 0;
    
    for(long j = i; j <= degY; j++){
      coef2 = 0;

      power(aa,a,j-i);
      
      mul(coef2, Q.rep[j], aa);
      //std::cout << "." ;
      mul(coef2, coef2, comb(j,i));
      // std::cout << "." ;
      add(coef1,coef1,coef2);
      //std::cout << "." << endl;
      //      std::cout << ".." << coef2 << " " << aa << endl;

    } // for j

    //std::cout << "-shift " << coef1 << endl;
    coef1 = shiftX(coef1, i);
    //    std::cout << "." << coef1 << endl;
    temp.rep[i] = coef1;

  } // for i
  
  //std::cout << "start back-shift" << endl;
  //std::cout << temp << " -- " << minX(temp) << endl;
  temp.normalize();
  temp = backShiftX(temp, minX(temp));
  //std::cout << "end back-shift" << temp << endl;



  return temp;
}



vec_ZZ_p findroots(const ZZ_pX &P){
  //
  // returns all roots of polynomial P(x)
  //   
  

  ZZ_pX Q;             //  
  Q = P;               // Q = monic version of P
  MakeMonic(Q);        //

  ZZ_pX      factor;   // one square free factor

  vec_ZZ_pX  factors;  // all factors (of factor) from Berlekamp
  long      nfactors;  //

  vec_ZZ_p  roots;     // list of roots of Q
  long     nroots = 0; //

  ZZ_p root;           // a single root of Q

  
  // first do square free decomposition 
  vec_pair_ZZ_pX_long sqrfreeQ = SquareFreeDecomp(Q);
  //  std::cout << "sqrfree = " << sqrfreeQ << endl;


  // work on each factor from square free decomposition
  for(long i = 0; i < sqrfreeQ.length(); i++){
    factor = sqrfreeQ[i].a;
    //    std::cout << "factor = " << sqrfreeQ[i] << endl;

    factors = SFBerlekamp(factor);
    nfactors = factors.length();

    //    std::cout << ".." << factors << " " << nfactors << endl;

    for(long j = 0; j < nfactors; j++){
      //    std::cout << factors[j] << deg( factors[j] ) << endl;
      if( deg(factors[j]) == 1 ){
	root = FindRoot( factors[j] );
	//	std::cout << "root " << root << endl; 
	append(roots, root); 
	//	std::cout << roots[nroots] << endl; 
	nroots++; 
      }
      
    }
    
  }
  
  return roots;

} 
























/*****************************************************************

                        Multiplication

******************************************************************/


//
// z = a * b
// 

void mul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b)
{
  long da = deg(a);
  long db = deg(b);
  
  // do the easy cases first

  if (da < 0 || db < 0){  // a or b is zero
    clear(x); 
    return; 
  }   

  if (da == 0 && db == 0){  // a and b are ZZ_pX
    x.rep.SetLength(1);     // x is ZZ_pX too then
    mul( x.rep[0], a.rep[0], b.rep[0] );
    return; 
  }

  if (da == 0){ // a is ZZ_pX and b is ZZ_pXY
    mul(x,b,a.rep[0]);
    return;
  }

  if (db == 0){ // b is ZZ_pX and a is ZZ_pXY
    mul(x,a,b.rep[0]);
    return;
  }
  
  // da,db >= 1 
  
  PlainMul(x,a,b);
}



void mul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pX& b)
{
  long da = deg(a);
  long db = deg(b);
  
  // do the easy cases first

  if (da < 0 || db < 0){  // a or b is zero
    clear(x); 
    return; 
  }   

  if (da == 0){  // a and b are both ZZ_pX
    x.rep.SetLength(da+1);
    mul( x.rep[0], a.rep[0], b );
    return; 
  }  
  
  if (db == 0){  // b is ZZ_p
    mul(x, a, b.rep[0]); 
    return; 
  }  
  
  // da >= 1 (simply multiply each coefficient of a by b)
  
  x = a;
  for(long i = 0; i <= da; i++){
    mul(x.rep[i],a.rep[i],b);
  }
  x.normalize();

}


void mul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_p& b)
{
   if (IsZero(b)) {
      clear(x);
      return;
   }

   if (IsOne(b)) {
      x = a;
      return;
   }

   ZZ_pTemp TT; ZZ_p& t = TT.val();

   long i, da;

   const ZZ_pX *ap;
   ZZ_pX* xp;

   t = b;

   da = deg(a);
   x.rep.SetLength(da+1);
   ap = a.rep.elts();
   xp = x.rep.elts();

   for (i = 0; i <= da; i++) 
      mul(xp[i], ap[i], t);

   x.normalize();
}



void mul(ZZ_pXY& x, const ZZ_pXY& a, long b)
{
   ZZ_pTemp TT;  ZZ_p& T = TT.val();
   conv(T, b);
   mul(x, a, T);
}




//
// x = a^2
//

void sqr(ZZ_pXY& x, const ZZ_pXY& a)
{
  mul(x,a,a);  // I'll take the easy way out for now!
}

void sqr(ZZ_pXY& x, const ZZ_pX& a)
{
  ZZ_pX temp;
  sqr(temp,a);  
  conv(x,temp);
}

void sqr(ZZ_pXY& x, const ZZ_p& a)
{
  ZZ_p temp;
  sqr(temp,a);  
  conv(x,temp);
}

void sqr(ZZ_pXY& x, long a)
{
  ZZ_p ap = to_ZZ_p(a);
  ZZ_p temp;
  sqr(temp, ap);
  conv(x,temp);
}


//
// x = a^e 
//

void power(ZZ_pXY& x, const ZZ_pXY& a, long e)
{
  // do the easy stuff 
  if (e < 0) {
    Error("power: negative exponent");
  }
  
  if (e == 0) {  // a^0 = 1
    x = 1;
    return;
  }
  
  if (a == 0 || a == 1 || e == 1) { 
    x = a;
    return;
  }
  
  if (e == 2){
    sqr(x,a);
    return;
  }

   long da = deg(a);

   if (da == 0) {
      x = power(ConstTerm(a), e);
      return;
   }

   if (da > (NTL_MAX_LONG-1)/e)
      Error("overflow in power");

   ZZ_pXY res;
   res.SetMaxLength(da*e + 1);
   res = 1;
   
   long k = NumBits(e);
   long i;

   for (i = k - 1; i >= 0; i--) {
      sqr(res, res);
      if (bit(e, i))
         mul(res, res, a);
   }

   x = res;
}



//
//  Classical n^2 algorithms 
//


void PlainMul(ZZ_pXY& x, const ZZ_pXY& a, const ZZ_pXY& b)
{
  long da = deg(a);
  long db = deg(b);
  
  // do the easy cases first
  
  if (da < 0 || db < 0){  // a or b is zero
    clear(x); 
    return; 
  }   
  
  if (da == 0 && db == 0){  // a and b are ZZ_pX
    x.rep.SetLength(1);     // x is ZZ_pX too then
    mul( x.rep[0], a.rep[0], b.rep[0] );
    return; 
  }
  
  if (da == 0){ // a is ZZ_pX and b is ZZ_pXY
    mul(x,b,a.rep[0]);
    return;
  }
  
  if (db == 0){ // b is ZZ_pX and a is ZZ_pXY
    mul(x,a,b.rep[0]);
    return;
  }
  
  // da,db >= 1 
  
  long d = da+db;  // degree of new polynomial
  
  ZZ_pXY la, lb;    // create new copies of a and b
  la = a;           // just in case &x == &a or &x == &b
  lb = b;           // 

  x.rep.SetLength(d+1);
  
  long i, j, jmin, jmax;
  static ZZ_pX t, accum;
  
  for (i = 0; i <= d; i++) {
    jmin = max(0, i-db);
    jmax = min(da, i);
    clear(accum);
    for (j = jmin; j <= jmax; j++) {
      mul(t, la.rep[j], lb.rep[i-j]);
      //mul(t, rep(ap[j]), rep(bp[i-j]));
      add(accum, accum, t);
    }
    SetCoeff( x, i, accum);
  }
  x.normalize();
}


long operator==(const ZZ_pXY& a, long b)
{
   if (b == 0)
      return IsZero(a);

   if (b == 1)
      return IsOne(a);

   long da = deg(a);

   if (da > 0)
      return 0;

   ZZ_p bb;
   bb = b;

   if (da < 0)
      return IsZero(bb);

   return a.rep[0] == bb;
}

long operator==(const ZZ_pXY& a, const ZZ_p& b)
{
   if (IsZero(b))
      return IsZero(a);

   long da = deg(a);

   if (da != 0)
      return 0;

   return a.rep[0] == b;
}

long operator==(const ZZ_pXY& a, const ZZ_pX& b)
{
   if (IsZero(b))
      return IsZero(a);

   long da = deg(a);

   if (da != 0)
      return 0;

   return a.rep[0] == b;
}

NTL_END_IMPL
