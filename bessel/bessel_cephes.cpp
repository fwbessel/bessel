/*
   This is bessel function evaluation and integral library based on cephes.
  
   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of the <organization> nor the
      names of its contributors may be used to endorse or promote products
      derived from this software without specific prior written permission.

  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
  ARE DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
  THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
  
*/


#include <cassert>
#include <cmath>
#include <iostream>
#include "bessel_cephes.h"

extern "C"
{
  extern double cephes_special_functions_j0(double x);
  extern double cephes_special_functions_j1(double x);
  extern double cephes_special_functions_jv(double n, double x );
  extern double cephes_special_functions_struve ( double v, double x );
  extern double cephes_special_functions_gamma ( double x );
}



namespace cephes 
{
// precomputed data file 
#include "bessel_zeros_cephes.h"
#include "bessel_jp_zeros_cephes.h"
#include "bessel_jpp_zeros_cephes.h"
#include "gamma_half.h"
#include "polynomial_bessel_integral.h"

double bessel_j0(double x)
{
  return cephes_special_functions_j0(x);
}


double bessel_j1(double x)
{
  return cephes_special_functions_j1(x);
}


double bessel_jn(int n, double x)
{
  return cephes_special_functions_jv(double(n), x);
}


void bessel_jn_array(const int nmax, const double x, double *result_array)
{
  double Jnp1;
  double Jn;
  double Jnm1;
  int n;
  
  if( x < 1e-6 && x > -1e-6 )
  {
    result_array[0] = 1.0;
    for(n=1; n<=nmax; n++)
      result_array[n] = 0.0;
    return;
  }
  
  Jnp1 = bessel_jn(nmax+1, x);
  Jn   = bessel_jn(nmax, x);

  for(n=nmax; n>=0; n--) 
  {
    result_array[n] = Jn;
    Jnm1 = -Jnp1 + 2.0*n/x * Jn;
    Jnp1 = Jn;
    Jn   = Jnm1;
  }
}

void bessel_jn_array(const int nmin, const int nmax, const double x, double *result_array)
{
  double Jnp1;
  double Jn;
  double Jnm1;
  int n;
  
  Jnp1 = bessel_jn(nmax+1, x);
  Jn   = bessel_jn(nmax, x);

  for(n=nmax; n>=nmin; n--) 
  {
    result_array[n] = Jn;
    Jnm1 = -Jnp1 + 2.0*n/x * Jn;
    Jnp1 = Jn;
    Jn   = Jnm1;
  }
}


void bessel_jn_array(const int nmax, const double x,  double Jn,  double Jnp, double *result_array)
{
  double Jnm1;
  int n;
  
  for(n=nmax; n>=0; n--) 
  {
    result_array[n] = Jn;
    Jnm1 = -Jnp + 2.0*n/x * Jn;
    Jnp  = Jn;
    Jn   = Jnm1;
  }
}


double bessel_jnc(const int n, const double x)
{
  return bessel_jn(n, x)/x;
}


double bessel_jn_derivative(const int n, const double x)
{
  if(n==0) return -bessel_j1(x);
  return 0.5*(bessel_jn(n-1,x) - bessel_jn(n+1,x));
}


double bessel_jn_second_order_derivative(const int n, const double x)
{
  if(n==0) return -bessel_jn_derivative(1, x);
  return 0.5*(bessel_jn_derivative(n-1,x) - bessel_jn_derivative(n+1,x));
}


double bessel_root(const int n, const int s)
{
  return bessel_root_cephes[n][s-1];
}


double bessel_derivative_root(const int n, const int s)
{
  return bessel_jp_root_cephes[n][s-1];
}


double bessel_second_order_derivative_root(const int n, const int s)
{ 
  return bessel_jpp_root_cephes[n][s-1]; 
}


double bessel_j0_integral(double x)
{
  return x*cephes_special_functions_j0(x)+0.5*M_PI*x*(cephes_special_functions_j1(x)*cephes_special_functions_struve(0,x) - cephes_special_functions_j0(x)*cephes_special_functions_struve(1,x));
}


double bessel_jn_integral(int n, double x)
{
  if( n == 0 ) return bessel_j0_integral(x);
  
  int nn = n/2;
  if( n%2 == 1)
  {
    double res = 1 - cephes_special_functions_j0(x);
    for(int k=1; k<=nn; k++)
      res -= 2* bessel_jn(2*k,x);
    return res;
  }
  
  if( n%2 == 0)
  {
    double res = bessel_j0_integral(x);
    for(int k=0; k<nn; k++)
      res -= 2* bessel_jn(2*k+1,x);
    return res;
  }
}



double bessel_jn_tp_integral(int n, int p, double x)
{
  // p > n and p-n is odd

  if( p > n && (p-n)%2 == 1)
  {
    if( p == n+1 )
      return std::pow(x, 1+n) * bessel_jn(n+1,x);

    // FIXME
  }
  
  // otherwise
#if 0  
  //NOTE: cephes_special_functions_gamma(x) x can be as large as 171.6
  double xx = std::pow(x, p);
  double gg = cephes_special_functions_gamma(0.5*(n + p + 1))/cephes_special_functions_gamma(0.5*(n - p + 1));
  
  double ss = 0.0;
  int k=0;
  while(n+2*k+1 < x + 30.0)
  {
    double kk = cephes_special_functions_gamma(0.5*n + 0.5 + k - 0.5*p)/cephes_special_functions_gamma(0.5*n + 0.5 + k + 0.5*p + 1);
    double zz = (n+2*k+1)*kk*bessel_jn(n+2*k+1,x);
    k++;
    ss += zz;
    if( std::abs(zz) < 1e-6*std::abs(ss)) break;
  }
  return xx*gg*ss;
#endif
  
#if 1  
  //NOTE: gamma_half[x] x in [0,1,2..343]
  double xx = std::pow(x, p);
  double gg = gamma_half[(n + p + 1)]/gamma_half[(n - p + 1)];
  
  double ss = 0.0;
  int k=0;
  
  int nmax = int(x + 30.0);
  double *Jn_array = new double [1+nmax];
  bessel_jn_array(n+1, nmax, x, Jn_array);
  
  while(n+2*k+1 < nmax)
  {
    double kk = gamma_half[(n + 1 + 2*k - p)]/gamma_half[(n + 1 + 2*k + p + 2)];
    double zz = (n+2*k+1)*kk*Jn_array[n+2*k+1];
    k++;
    ss += zz;
    if( std::abs(zz) < 1e-6*std::abs(ss)) break;
  }
  
  delete [] Jn_array;
  
  return xx*gg*ss;
#endif  
  
}

#if 0
// has large numerical error if x is large (x > 2.0)
double bessel_jn_tp_integral(int n, int p, double x)
{
  double xx = std::pow(x, n+p+1)/std::pow(2.0, n);
  double ss = 0.0;
  int k=0;
  while(k<100)
  {
    double sss = std::pow(-1, k)*std::pow(0.5*x, 2*k) / cephes_special_functions_gamma(k+1) / cephes_special_functions_gamma(n+k+1)  / (n+p+1+2*k);
    k++;
    ss += sss;
    if( std::abs(sss) < 1e-6*std::abs(ss) ) break;
  }
  
  return xx*ss;
}
#endif


double bessel_j0_t1_integral(double r, double x)
{
  return 1.0/(r*r)*(r*x*bessel_j1(r*x));  
}


double bessel_j2_t1_integral(double r, double x)
{
  return 1.0/(r*r)*(-2*bessel_j0(r*x)-r*x*bessel_j1(r*x));    
}


double bessel_jn_t1_integral_intp(int n, double r, double x)
{
  double y = r*x/(2*M_PI);

  if( n>=33 ||  y > 15.0 ) 
    return bessel_jn_tp_integral(n, 1, r, x);
  
  
  const double dy = 0.005; //
  
  int intp_index_1 = floor(y/dy);
  int intp_index_2 = ceil(y/dy);

  double fr_1 = jn_t1_integral[n][intp_index_1];
  double fr_2 = jn_t1_integral[n][intp_index_2];
    
  double fr = 0.0;
  if(intp_index_1 == intp_index_2) 
  {
    fr = fr_1;
  }
  else 
  {
    double intp_r1 =  intp_index_1*dy;
    double intp_r2 =  intp_index_2*dy;
    fr  = fr_1 + (fr_2-fr_1)/(intp_r2-intp_r1)*(y - intp_r1);
  }

  return (x*x) * fr;  
}



void bessel_jn_t1_integral_intp(int n_max, double r, double x, double *result_array)
{
  double y = r*x/(2*M_PI);

  if( y > 15.0 ) return;
  
  n_max = std::min(n_max, 32);
  
  const double dy = 0.005; //
  
  int intp_index_1 = floor(y/dy);
  int intp_index_2 = ceil(y/dy);

  double intp_r1 =  intp_index_1*dy;
  double intp_r2 =  intp_index_2*dy;
    
  for(int n=0; n<=n_max; n++)
  {
    double fr_1 = (x*x)*jn_t1_integral[n][intp_index_1];
    double fr_2 = (x*x)*jn_t1_integral[n][intp_index_2];

    if(intp_index_1 == intp_index_2) 
    {
      result_array[n] = fr_1;
    }
    else 
    {
      result_array[n]  = fr_1 + (fr_2-fr_1)/(intp_r2-intp_r1)*(y - intp_r1);
    }
  }
}



double bessel_jn_tp_integral(int n, int p, double r, double x)
{
  return 1.0/std::pow(r, p+1) * bessel_jn_tp_integral(n, p, r*x);
}





double factorial(double n)
{
  return cephes_special_functions_gamma ( n+1 );
}

}



