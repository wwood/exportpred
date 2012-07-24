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
#ifndef GHMM_GHMM_DISTRIBUTION_HH_INCLUDED
#define GHMM_GHMM_DISTRIBUTION_HH_INCLUDED

#include <GHMM/ref.hh>
#include <GHMM/math.hh>

namespace GHMM {
  namespace DISTRIBUTION {
    class Base : public virtual RefObj {
    public:
      typedef Ref<Base> Ptr;
      Base &operator=(const Base &) { return *this; }
      Base(const Base &) {}
      Base() {}
      virtual ~Base() {}
    };

    class Discrete : public Base {
    };

    class Continuous : public Base {
    public:
      virtual double pdf(double z) const =0;
      virtual double cdf(double z) const =0;
    };

    class BoundedContinuous : public Continuous {
    public:
      virtual double lbound() const =0;
      virtual double rbound() const =0;

      virtual MATH::DPDF::Ptr discretise() {
        int x1 = (int)floor(lbound() + 0.5);
        int x2 = (int)floor(rbound() + 0.5);
        MATH::DPDF::Ptr result = new MATH::DPDF();
        result->setDistrib(x1, x2, 0.0);
        double c = 0.0;
        for (int x = x1; x < x2; x++) {
          double nc = cdf(x + 0.5);
          result->setp(x, nc - c);
          c = nc;
        }
        return result;
      }
    };
  }
}

#endif
