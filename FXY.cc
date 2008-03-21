#include "FXY.h"

NTL_OPEN_NNS

ZZ comb(long n, long k)
{
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

NTL_CLOSE_NNS
