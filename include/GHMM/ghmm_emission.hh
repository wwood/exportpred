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
#ifndef GHMM_GHMM_EMISSION_HH_INCLUDED
#define GHMM_GHMM_EMISSION_HH_INCLUDED

#include <GHMM/ghmm_global.hh>

#include <GHMM/ref.hh>
#include <GHMM/math.hh>

#include <vector>

namespace GHMM {
  namespace EMISSION {
    class Base : public virtual RefObj {
    public:
      typedef Ref<Base> Ptr;
      Base &operator=(const Base &) { return *this; }
      Base(const Base &) {}
      Base() {}
      virtual ~Base() {}
    };

    class PositionSpecific : public EMISSION::Base {
    protected:
      std::vector<MATH::DPDF::Ptr> pssm;

    public:
      typedef Ref<PositionSpecific> Ptr;

      template<int DIR>
      class Generator {
        Generator();
        Generator(const Generator &);
        Generator &operator=(const Generator &);

      protected:
        const std::vector<MATH::DPDF::Ptr> &pssm;
        const int *seq;
        bool emitted;

      public:
        Generator(const PositionSpecific *p, const int *s, int d_min, int d_max) : pssm(p->pssm), seq(s), emitted(true) {
          emitted = !((int)pssm.size() >= d_min && (int)pssm.size() < d_max);
        }
        bool gen(int &d, double &logp) {
          if (!emitted) {
            d = pssm.size();
            const int *p = DIR == +1 ? seq : seq - d + 1;
            logp = 0.0;
            for (int i = 0; i < d; i++) {
              logp += pssm[i]->logp(*p++);
            }
            emitted = true;
            return true;
          }
          return false;
        }
      };

      PositionSpecific &operator=(const PositionSpecific &ps) {
        if (this != &ps) {
          pssm = ps.pssm;
        }
        return *this;
      }
      PositionSpecific() : EMISSION::Base(), pssm() {
      }
      PositionSpecific(const PositionSpecific &ps) : EMISSION::Base(), pssm() {
        *this = ps;
      }
      PositionSpecific(const std::vector<MATH::DPDF::Ptr> &p) : EMISSION::Base(), pssm() {
        setEmissionDistrib(p);
      }
      bool setEmissionDistrib(const std::vector<MATH::DPDF::Ptr> &p) {
        pssm = p;
        return true;
      }

      void randSequence(std::vector<int> &result, int d) const {
        assert(d == (int)pssm.size());
        for (int i = 0; i < d; i++) {
          result.push_back(pssm[i]->randZ());
        }
      }

      virtual ~PositionSpecific() {
      }
    };

    class Stateless : public EMISSION::Base, protected MATH::DPDF {
    protected:

    public:
      typedef Ref<Stateless> Ptr;

      template<int DIR>
      class Generator {
        Generator();
        Generator(const Generator &);
        Generator &operator=(const Generator &);

      protected:
        const Stateless *emit;
        const int *seq;
        int cur;
        int end;
        double accum;

      public:
        Generator(const Stateless *e, const int *s, int d_min, int d_max) : emit(e), seq(s), cur(0), end(d_max - 1), accum(0.0) {
          if (d_min == d_max) {
            cur = end;
          } else {
            while (cur < d_min - 1) {
              accum += emit->logp(seq[DIR * cur]);
              cur++;
            }
          }
        }
        bool gen(int &d, double &logp) {
          if (cur < end) {
            accum += emit->logp(seq[DIR * cur]);
            cur++;
            d = cur;
            logp = accum;
            return true;
          }
          return false;
        }
      };

      Stateless &operator=(const Stateless &d) {
        if (this != &d) {
          EMISSION::Base::operator=(d);
          MATH::DPDF::operator=(d);
        }
        return *this;
      }
      Stateless(const Stateless &d) : EMISSION::Base(), MATH::DPDF(d) {
        *this = d;
      }
      Stateless() : EMISSION::Base(), MATH::DPDF() {
      }
      Stateless(const MATH::DPDF::Ptr &dpdf) : EMISSION::Base(), MATH::DPDF(*dpdf) {
      }

      bool setEmissionDistrib(int l, double *d) {
        return setDistrib(0, l, d);
      }
      bool setEmissionDistrib(int l, double d) {
        return setDistrib(0, l, d);
      }
      bool setEmissionProb(int i, double d) {
        return setp(i, d);
      }
      bool normalizeEmissionDistrib() {
        return normalize();
      }

      void randSequence(std::vector<int> &result, int d) const {
        for (int i = 0; i < d; i++) {
          result.push_back(randZ());
        }
      }

      virtual ~Stateless() {
      }
    };
  }
}

#endif
