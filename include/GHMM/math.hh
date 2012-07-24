// Copyright (c) 2005 The Walter and Eliza Hall Institute
// 
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
// ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef GHMM_MATH_HH_INCLUDED
#define GHMM_MATH_HH_INCLUDED

#include <GHMM/ref.hh>
#include <math.h>
#include <stdlib.h>

#include <numeric>
#include <algorithm>

namespace MATH {
  const static double LOG_ZERO = -__builtin_inf(); // -1e15;
  const static double LOG_INF  = +__builtin_inf(); // +1e15;

  static inline double logAdd(double x, double y) {
    double r;
    if (x < y) {
      r = y + log1p(exp(x - y));
    } else {
      r = x + log1p(exp(y - x));
    }
    assert(r != LOG_ZERO || (x == LOG_ZERO && y == LOG_ZERO));
    return r;
  }

  static inline double logClip(double x) {
    return std::max(std::min(log(x), LOG_INF), LOG_ZERO);
  }

  template<typename T>
  static inline T clamp(const T val, const T val_min, const T val_max) {
    return std::min(std::max(val, val_min), val_max);
  }

  class DPDF : public virtual RefObj {
  protected:
    int min_d, max_d;
    double *distrib;
    double *log_distrib;

    void reallocDistrib(int l) {
      if (distrib) delete [] distrib;
      if (log_distrib) delete [] log_distrib;
      distrib = new double[l];
      log_distrib = new double[l];
    }
    void updateLogDistrib() {
      for (int i = 0; i < max_d - min_d; i++) {
        log_distrib[i] = MATH::logClip(distrib[i]);
      }
    }
    void updateDistrib() {
      for (int i = 0; i < max_d - min_d; i++) {
        distrib[i] = exp(log_distrib[i]);
      }
    }

  public:
    typedef Ref<DPDF> Ptr;

    DPDF &operator=(const DPDF &d) {
      if (this != &d) {
        reallocDistrib(d.max_d - d.min_d);
        std::copy(d.distrib, d.distrib + d.max_d - d.min_d, distrib);
        std::copy(d.log_distrib, d.log_distrib + d.max_d - d.min_d, log_distrib);
        min_d = d.min_d;
        max_d = d.max_d;
      }
      return *this;
    }

    DPDF(const DPDF &d) : RefObj(), min_d(0), max_d(0), distrib(NULL), log_distrib(NULL) {
      *this = d;
    }

    DPDF() : RefObj(), min_d(0), max_d(0), distrib(NULL), log_distrib(NULL) {
    }

    virtual ~DPDF() {
      if (distrib) delete [] distrib;
      if (log_distrib) delete [] log_distrib;
    }

    bool normalize() {
      double freq_sum = std::accumulate(distrib, distrib + max_d - min_d, 0.0);
      if (freq_sum <= 0.0) return false;
      for (int i = 0; i < max_d - min_d; i++) {
        distrib[i] /= freq_sum;
      }
      updateLogDistrib();
      return true;
    }

    template<typename T>
    bool setDistrib(int a, int b, T freqs, bool norm = true) {
      if (b <= a) return false;

      min_d = a;
      max_d = b;
      reallocDistrib(b - a);
      std::copy(freqs, freqs + b - a, log_distrib);
      if (norm) {
        normalize();
      } else {
        updateLogDistrib();
      }
      return true;
    }

    template<typename T>
    bool setDistrib(int l, T freqs) {
      return setDistrib(0, l, freqs);
    }

    template<typename T>
    bool setLogDistrib(int a, int b, T log_freqs) {
      if (b <= a) return false;

      min_d = a;
      max_d = b;
      reallocDistrib(b - a);
      std::copy(log_freqs, log_freqs + b - a, log_distrib);
      updateDistrib();
    }

    template<typename T>
    bool setLogDistrib(int l, T log_freqs) {
      return setLogDistrib(0, l, log_freqs);
    }

    bool setDistrib(int a, int b, double d) {
      if (b <= a) return false;
      min_d = a;
      max_d = b;
      reallocDistrib(b - a);
      std::fill(distrib, distrib + b - a, d);
      std::fill(log_distrib, log_distrib + b - a, MATH::logClip(d));
      return true;
    }

    bool setLogDistrib(int a, int b, double d) {
      if (b <= a) return false;
      min_d = a;
      max_d = b;
      reallocDistrib(b - a);
      std::fill(distrib, distrib + b - a, exp(d));
      std::fill(log_distrib, log_distrib + b - a, d);
      return true;
    }

    double p(int l) const {
      if (l < min_d || l >= max_d) return 0.0;
      return distrib[l - min_d];
    }
    bool setp(int i, double d) {
      if (i < min_d || i >= max_d) return false;
      distrib[i - min_d] = d;
      log_distrib[i - min_d] = MATH::logClip(d);
      return true;
    }

    double logp(int l) const {
      if (l < min_d || l >= max_d) return MATH::LOG_ZERO;
      return log_distrib[l - min_d];
    }
    bool setlogp(int i, double d) {
      d = std::max(std::min(d, MATH::LOG_INF), MATH::LOG_ZERO);
      if (i < min_d || i >= max_d) return false;
      log_distrib[i - min_d] = d;
      distrib[i - min_d] = exp(d);
      return true;
    }

    int distribMin() const {
      return min_d;
    }

    int distribMax() const {
      return max_d;
    }
    int randZ() const {
      double r = random() / double(RAND_MAX);
      for (int i = 0; i < max_d - min_d; i++) {
        r -= distrib[i];
        if (r <= 0.0) return i + min_d;
      }
      return max_d - 1;
    }

    std::ostream &dump(std::ostream &o) const {
      o << min_d << " [ ";
      for (int i = min_d; i < max_d; i++) {
        o << p(i) << " ";
      }
      o << "] " << max_d << std::endl;
      return o;
    }
  };

  template<int base, int power> struct __pow__          { enum { val = __pow__<base, power >> 1>::val * __pow__<base, power - (power >> 1)>::val }; };
  template<int base>            struct __pow__<base, 1> { enum { val = base }; };
  template<int base>            struct __pow__<base, 0> { enum { val = 1 }; };

  template<int n> struct __fact__    { enum { val = __fact__<n - 1>::val * n }; };
  template<>      struct __fact__<0> { enum { val = 1 }; };

  static inline double gaussianCDF(double z) {
    return .5 * (1.0 + erf(z / M_SQRT2));
  }

  static inline double gaussianCDF(double mu, double sigma, double z) {
    return gaussianCDF((z - mu) / sigma);
  }

  struct GaussianKernel {
    double sigma;
    GaussianKernel(double s) : sigma(s) {
    }
    double operator()(double mu, double z) {
      return gaussianCDF(mu, sigma, z);
    }
  };

  template<typename Kernel, typename InputIterator>
  DPDF::Ptr smooth(Kernel kernel, InputIterator begin, InputIterator end, int distrib_min = -1, int distrib_max = -1) {
    DPDF::Ptr result;
    int x1, x2;
    x1 = distrib_min >= 0 ? distrib_min : *std::min_element(begin, end);
    x2 = distrib_max >= 0 ? distrib_max : *std::max_element(begin, end);
    if (x1 >= x2) return NULL;
    result = new MATH::DPDF();
    result->setDistrib(x1, x2 + 1, 0.0);
  
    for (int x = x1; x <= x2; x++) {
      double s = 0.0;
      for (InputIterator i = begin; i != end; i++) {
        s += kernel(*i, x + 0.5) - kernel(*i, x - 0.5);
      }
      result->setp(x, s);
    }

    result->normalize();

    return result;
  }
}

#endif
