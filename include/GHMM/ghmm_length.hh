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
#ifndef GHMM_GHMM_LENGTH_HH_INCLUDED
#define GHMM_GHMM_LENGTH_HH_INCLUDED

#include <GHMM/ref.hh>
#include <GHMM/math.hh>

namespace GHMM {
  namespace LENGTH {
    class Base : public virtual RefObj {
    public:
      typedef Ref<Base> Ptr;
      Base &operator=(const Base &) { return *this; }
      Base(const Base &) {}
      Base() {}
      virtual ~Base() {}
    };

    class Uniform : public LENGTH::Base, public DISTRIBUTION::BoundedContinuous {
    protected:
      double d_min, d_max;
      double v_pdf;

    public:
      typedef Ref<Uniform> Ptr;

      Uniform &operator=(const Uniform &u) {
        if (this != &u) {
          LENGTH::Base::operator=(u);
          DISTRIBUTION::BoundedContinuous::operator=(u);
          d_min = u.d_min;
          d_max = u.d_max;
          v_pdf = u.v_pdf;
        }
        return *this;
      }

      Uniform(const Uniform &u) : LENGTH::Base(), DISTRIBUTION::BoundedContinuous(), d_min(0.0), d_max(1.0), v_pdf(1.0) {
        *this = u;
      }
      Uniform() : LENGTH::Base(), DISTRIBUTION::BoundedContinuous(), d_min(0.0), d_max(1.0), v_pdf(1.0) {
      }
      Uniform(double a, double b) : LENGTH::Base(), DISTRIBUTION::BoundedContinuous(), d_min(a), d_max(b), v_pdf(1.0 / (b - a)) {
      }
      ~Uniform() {
      }

      virtual double lbound() const {
        return d_min;
      }
      virtual double rbound() const {
        return d_max;
      }
      virtual double pdf(double z) const {
        if (z >= d_max) return 0.0;
        if (z <= d_min) return 0.0;
        return v_pdf;
      }
      virtual double cdf(double z) const {
        if (z >= d_max) return 1.0;
        if (z <= d_min) return 0.0;
        return (z - d_min) * v_pdf;
      }
    };

    class Geometric : public LENGTH::Base, public DISTRIBUTION::Continuous {
    protected:
      double p_self;

    public:
      typedef Ref<Geometric> Ptr;
      Geometric &operator=(const Geometric &g) {
        if (this != &g) {
          LENGTH::Base::operator=(g);
          DISTRIBUTION::Continuous::operator=(g);
          p_self = g.p_self;
        }
        return *this;
      }

      Geometric() : LENGTH::Base(), DISTRIBUTION::Continuous(), p_self(0) {
      }
      Geometric(const Geometric &g) : LENGTH::Base(), DISTRIBUTION::Continuous(), p_self(0) {
        *this = g;
      }
      Geometric(double mean) : LENGTH::Base(), DISTRIBUTION::Continuous(), p_self(0) {
        setMean(mean);
      }

      void setPSelf(double p) {
        p_self = p;
      }
      void setMean(double mean) {
        p_self = mean / (1 + mean);
      }

      double pLength(int l) {
        return pow(p_self, l);
      }

      double logpLength(int l) const {
        return log(p_self) * l;
      }

      virtual double pdf(double z) const {
        return 0.0; // XXX
      }
      virtual double cdf(double z) const {
        return 0.0; // XXX
      }

      int minLength() const {
        return 1;
      }
      int maxLength() const {
        return 2;
      }
      int randLength() const {
        return 1;
      }

      virtual ~Geometric() {
      }
    };

    class Discrete : public LENGTH::Base, public DISTRIBUTION::Discrete, protected MATH::DPDF {
    protected:

    public:
      typedef Ref<Discrete> Ptr;

      Discrete &operator=(const Discrete &d) {
        if (this != &d) {
          LENGTH::Base::operator=(d);
          DISTRIBUTION::Discrete::operator=(d);
          MATH::DPDF::operator=(d);
        }
        return *this;
      }

      Discrete(const Discrete &d) : LENGTH::Base(), DISTRIBUTION::Discrete(), MATH::DPDF() {
        *this = d;
      }

      Discrete() : LENGTH::Base(), DISTRIBUTION::Discrete(), MATH::DPDF() {
      }

      Discrete(const MATH::DPDF::Ptr &dpdf) : LENGTH::Base(), DISTRIBUTION::Discrete(), MATH::DPDF(*dpdf) {
      }

      virtual ~Discrete() {
      }

      double pLength(int l) const {
        return p(l);
      }

      double logpLength(int l) const {
        return logp(l);
      }

      bool setLengthDistrib(int a, int b, double *d) {
        return setDistrib(a, b, d);
      }

      int minLength() const {
        return distribMin();
      }
      int maxLength() const {
        return distribMax();
      }
      int randLength() const {
        return MATH::DPDF::randZ();
      }

      std::ostream &dump(std::ostream &o) {
        return MATH::DPDF::dump(o);
      }
    };

    class Fixed : public LENGTH::Base, public DISTRIBUTION::Discrete {
    protected:
      int length;
    public:
      typedef Ref<Fixed> Ptr;
      
      Fixed &operator=(const Fixed &f) {
        if (this != &f) {
          LENGTH::Base::operator=(f);
          DISTRIBUTION::Discrete::operator=(f);
          length = f.length;
        }
        return *this;
      }
      Fixed(const Fixed &f) : LENGTH::Base(), DISTRIBUTION::Discrete(), length(1) {
        *this = f;
      }
      Fixed() : LENGTH::Base(), DISTRIBUTION::Discrete(), length(1) {
      }
      Fixed(int l) : LENGTH::Base(), DISTRIBUTION::Discrete(), length(l) {
      }

      double pLength(int l) const {
        return length == l ? 1.0 : 0.0;
      }

      double logpLength(int l) const {
        return length == l ? 0.0 : MATH::LOG_ZERO;
      }

      bool setLength(int l) {
        length = l;
        return true;
      }

      int minLength() const {
        return length;
      }
      int maxLength() const {
        return length + 1;
      }
      int randLength() const {
        return length;
      }

      virtual ~Fixed() {
      }
    };
  }
}

#endif
