/*
 *Copyright: Copyright (c) 2019
 *Created on 2019-5-21
 *Author:zhengaohong@zgheye.cc
 *Version 1.0.1
*/

#include <cmath>
#include <cfloat>
#include "unitily.h"
#include "defines.h"

namespace zgh {


bool equal(int x, int y) {
  return x == y;
}

bool equal(double x, double y) {
  if (x == y) {
    return true;
  }
  double abs_diff = fabs(x - y);
  x = fabs(x);
  y = fabs(y);
  double abs_max = x > y ? x : y;
  if (abs_max < DBL_MIN) {
    abs_max = DBL_MIN;
  }
  return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}

bool equal(float x, float y) {
  if (x == y) {
    return true;
  }
  float abs_diff = fabs(x - y);
  x = fabs(x);
  y = fabs(y);
  float abs_max = x > y ? x : y;
  if (abs_max < FLT_MIN) {
    abs_max = FLT_MIN;
  }
  return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * FLT_EPSILON);
}

double angle2rad(double angle) {
  return angle / 180.0 * PI;
}

double rad2angle(double rad) {
  return rad / PI * 180.0;
}

double angle_diff(double a, double b) {
  a -= b;
  while (a <= -PI) {
    a += 2.0 * PI;
  }
  while (a >= PI) {
    a -= 2.0 * PI;
  }
  if (a < 0.0) {
    a = -a;
  }
  return a;
}

double angle_diff_signed(double a, double b) {
  a -= b;
  while (a <= -PI) {
    a += 2.0 * PI;
  }
  while (a > PI) {
    a -= 2.0 * PI;
  }
  return a;
}

bool inRect(int h, int w, int x, int y) {
  return 0 <= x && x < h && 0 <= y && y < w;
}


double rotateAngle(double start_angle, double end_angle, int polarity) {
	double coverage;
	if (polarity == SAME_POL) {
		coverage = end_angle - start_angle;
	} else {
		coverage = start_angle - end_angle;
	}
	if (coverage < 0) {
    coverage += 2 * PI;
  }
	return coverage;
}

double sqr(double x) {
  return x * x;
}

int sqr(int x) {
  return x * x;
}

/*----------------------------------------------------------------------------*/
/** Computes the natural logarithm of the absolute value of
    the gamma function of x using Windschitl method.
    See http://www.rskey.org/gamma.htm

    The formula used is
    @f[
        \Gamma(x) = \sqrt{\frac{2\pi}{x}} \left( \frac{x}{e}
                    \sqrt{ x\sinh(1/x) + \frac{1}{810x^6} } \right)^x
    @f]
    so
    @f[
        \log\Gamma(x) = 0.5\log(2\pi) + (x-0.5)\log(x) - x
                      + 0.5x\log\left( x\sinh(1/x) + \frac{1}{810x^6} \right).
    @f]
    This formula is a good approximation when x > 15.
 */

double log_gamma_windschitl(double x) {
  return 0.918938533204673 + (x - 0.5) * log(x) - x
         + 0.5 * x * log(x * sinh(1.0 / x) + 1.0 / (810.0 * pow(x, 6.0)));
}

/*----------------------------------------------------------------------------*/
/*----------------------------- NFA computation ------------------------------*/
/*----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/
/** Computes the natural logarithm of the absolute value of
    the gamma function of x using the Lanczos approximation.
    See http://www.rskey.org/gamma.htm

    The formula used is
    @f[
      \Gamma(x) = \frac{ \sum_{n=0}^{N} q_n x^n }{ \Pi_{n=0}^{N} (x+n) }
                  (x+5.5)^{x+0.5} e^{-(x+5.5)}
    @f]
    so
    @f[
      \log\Gamma(x) = \log\left( \sum_{n=0}^{N} q_n x^n \right)
                      + (x+0.5) \log(x+5.5) - (x+5.5) - \sum_{n=0}^{N} \log(x+n)
    @f]
    and
      q0 = 75122.6331530,
      q1 = 80916.6278952,
      q2 = 36308.2951477,
      q3 = 8687.24529705,
      q4 = 1168.92649479,
      q5 = 83.8676043424,
      q6 = 2.50662827511.
 */

double log_gamma_lanczos(double x) {
  static double q[7] = {75122.6331530, 80916.6278952, 36308.2951477,
                 8687.24529705, 1168.92649479, 83.8676043424,
                 2.50662827511};
  double a = (x + 0.5) * log(x + 5.5) - (x + 5.5);
  double b = 0.0;
  for (int n = 0; n < 7; ++n) {
    a -= log(x + (double)n);
    b += q[n] * pow(x, (double)n);
  }
  return a + log(b);
}

double nfa(int n, int k, double p, double logNT) {
  static double inv[TABSIZE]; //  table to keep computed inverse values
  double tolerance = 0.1; //  an error of 10% in the result is accepted
  double log1term, term, bin_term, mult_term, bin_tail, err, p_term;

  /* trivial cases */
  if (n == 0 || k == 0) {
    return -logNT;
  }
  if (n == k) {
    return -logNT - (double)n * log10(p);
  }

  /* probability term */
  p_term = p / (1.0 - p);

  log1term = log_gamma((double)n + 1.0) - log_gamma((double)k + 1.0)
           - log_gamma((double)(n - k) + 1.0)
           + (double)k * log(p) + (double)(n - k) * log(1.0 - p);
  term = exp(log1term);

  /* in some cases no more computations are needed */
  if (equal(term, 0.0)) {     /* the first term is almost zero */
    if((double)k > (double)n * p) {       /* at begin or end of the tail?  */
      return -log1term / M_LN10 - logNT;  /* end: use just the first term  */
    } else {
      return -logNT;                      /* begin: the tail is roughly 1  */
    }
  }

  /* compute more terms if needed */
  bin_tail = term;
  for (int i = k + 1; i <= n; ++i) {
    bin_term = (double)(n - i + 1) * (i < TABSIZE ?
               (inv[i] != 0.0 ? inv[i] : (inv[i] = 1.0 / (double)i)) :
                1.0 / (double)i);

    mult_term = bin_term * p_term;
    term *= mult_term;
    bin_tail += term;
    if (bin_term < 1.0) {
      err = term * ((1.0 - pow(mult_term, (double)(n - i + 1))) /
                         (1.0 - mult_term) - 1.0);

      if (err < tolerance * fabs(-log10(bin_tail) - logNT) * bin_tail) {
        break;
      }
    }
  }
  double nfavalue = -log10(bin_tail) - logNT;
  return nfavalue;
}

} // namespace zgh
