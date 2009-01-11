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


#ifndef __bessel_cephes_h__
#define __bessel_cephes_h__


// use bessel functions from CEPHES MATHEMATICAL FUNCTION LIBRARY
// http://www.netlib.org/cephes/
namespace cephes 
{

/**
 * Bessel function of order zero
 *
 * @Returns Bessel function of order zero of the argument.
*/
double bessel_j0(const double x);


/**
 * Bessel function of order one
 *
 * @Returns Bessel function of order one of the argument.
*/
double bessel_j1(const double x);


/**
 * Bessel function of integer order
 *
 * @Returns Bessel function of order n, where n is a (possibly negative) integer.
 */
double bessel_jn(const int n, const double x);


/**
 * calculate Bessel function of order [0, nmax] at x, result_array should have the length as [0, nmax]
 */
void bessel_jn_array(const int nmax, const double x, double *result_array);


/**
 * calculate Bessel function of order [nmin, nmax] at x, result_array should have the length as [0, nmax]
 */
void bessel_jn_array(const int nmin, const int nmax,  const double x, double *result_array);


/**
 * calculate Bessel function of order [0, nmax] at x, given Jn and Jn+1, result_array should have the length as [0, nmax]
 */
void bessel_jn_array(const int nmax, const double x,  double Jn,  double Jnp, double *result_array);


/**
 * Jn(x)/x
 */
double bessel_jnc(const int n, const double x);


/**
 * derivative of integer order Bessel 
 */
double bessel_jn_derivative(const int n, const double x);

/**
 * second order derivative of integer order Bessel 
 */
double bessel_jn_second_order_derivative(const int n, const double x);


/**
 * get the sth positive root of bessel function of order n
 * this is the precomputed value, for n in [0,64], s[1,131]
 */
double bessel_root(const int n, const int s);


/**
 * get the value of J'n(root(s)), when s is the sth root of bessel Jn
 * this is the precomputed value, for n in [0,64], s[1,131]
 */
double bessel_derivative_root(const int n, const int s);

/**
 * get the value of J''n(root(s)), when s is the sth root of bessel Jn
 * this is the precomputed value, for n in [0,64], s[1,131]
 */
double bessel_second_order_derivative_root(const int n, const int s);


/**
 * compute \int_0^x j0(t) dt
 */
double bessel_j0_integral(double x);


/**
 * compute \int_0^x jn(t) dt, n>=0
 */
double bessel_jn_integral(int n, double x);

/**
 * compute \int_0^x j0(r t) t dt
 */
double bessel_j0_t1_integral(double r, double x);

/**
 * compute \int_0^x j2(r t) t dt
 */
double bessel_j2_t1_integral(double r, double x);

/**
 * compute \int_0^x jn(r t) t dt from pre-computed data
 */
double bessel_jn_t1_integral_intp(int n, double r, double x);

/**
 * compute \int_0^x jn(r t) t dt from pre-computed data
 * n in the range of [0, n_max], result_array should be pre-allcated and set to 0.0
 */
void bessel_jn_t1_integral_intp(int n_max, double r, double x, double *result_array);

/**
 * compute \int_0^x jn(t) t^p dt
 * if we want to calculate \int_0^x jn(a t) t^p dt, we can do this as 1/a^{p+1}  \int_0^ax jn(t) t^p dt
 */
double bessel_jn_tp_integral(int n, int p, double x);


/**
 * compute \int_0^x jn(r t) t^p dt
 */
double bessel_jn_tp_integral(int n, int p, double r, double x);


/**
 * calculate factorial by gamma function
 */
double factorial(double n);
}




#endif


//#include "bessel_cephes.h"

