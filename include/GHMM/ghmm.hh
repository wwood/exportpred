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
#ifndef GHMM_GHMM_HH_INCLUDED
#define GHMM_GHMM_HH_INCLUDED

#include <GHMM/ghmm_global.hh>

#include <GHMM/ref.hh>
#include <GHMM/math.hh>

#include <GHMM/ghmm_distribution.hh>
#include <GHMM/ghmm_length.hh>
#include <GHMM/ghmm_emission.hh>

#include <vector>
#include <map>
#include <string>
#include <numeric>
#include <iostream>

namespace GHMM {
  class Parse;
  class Model;
  
  class StateBase : public virtual RefObj {
  public:
    typedef Ref<StateBase> Ptr;
    StateBase &operator=(const StateBase &) { return *this; }
    StateBase(const StateBase &) {}
    StateBase() {}
    virtual ~StateBase() {}

    virtual void preParse(int, const int *) {
    }

    virtual void nextToken(int, const int *) {
    }

    virtual void      delta(int j, const Model &model, const Parse &parse, int max_len, double &delta, int &prev_state, int &state_length) const = 0;
    virtual void      alpha(int j, const Model &model, const Parse &parse, int max_len, double &alpha) const = 0;
    virtual void alphaDelta(int j, const Model &model, const Parse &parse, int max_len, double &alpha, double &delta, int &prev_state, int &state_length) const = 0;
    virtual void       beta(int j, const Model &model, const Parse &parse, int max_len, double &beta) const = 0;

    virtual int generate(std::vector<int> &result) const = 0;
  };

  class Traceback {
  private:
    Traceback();
    Traceback &operator=(const Traceback &);
    Traceback(const Traceback &);

    friend class RefCounter<Traceback>;

    mutable int __refcount;

    inline int refcount() const {
      return __refcount;
    }
    void incref() const {
      __refcount++;
    }
    void decref() const {
      !--__refcount;
      if (!__refcount) {
        delete this;
      }
    }

    ~Traceback() {
    }
  public:

    typedef Ref<Traceback> Ptr;

    Ptr prev;
    unsigned state;
    unsigned length;

    Traceback(const Ptr &p, int s, int l) : __refcount(0), prev(p), state(s), length(l) {
    }
  };

  class Model : public virtual RefObj {
    Model();
    Model(const Model &);
    Model &operator=(const Model &);

  protected:
    std::vector<std::string> state_names;
    std::map<std::string, int> state_name_map;
    std::vector<std::vector<int> > pred_states;
    std::vector<std::vector<int> > succ_states;
    std::vector<StateBase::Ptr> states;
    double *state_trans;
    double *state_log_trans;
    int state_count;

  public:
    static const std::string BEGIN;
    static const std::string END;

    typedef Ref<Model> Ptr;

    Model(const std::vector<std::pair<std::string, StateBase::Ptr> > &, const std::map<std::pair<int, int>, double> &);
    virtual ~Model();

    int stateNumber(const std::string &s) const {
      std::map<std::string, int>::const_iterator i = state_name_map.find(s);
      if (i != state_name_map.end()) {
        return (*i).second;
      } else {
        return -1;
      }
    }
    const StateBase::Ptr &state(int n) const {
      return states[n];
    }
    const std::string &stateName(int n) const {
      return state_names[n];
    }
    int stateCount() const {
      return states.size();
    }
    const std::vector<int> predStates(int n) const {
      return pred_states[n];
    }
    const std::vector<int> succStates(int n) const {
      return succ_states[n];
    }

    double &logp(int s, int t) {
      return state_log_trans[s * state_count + t];
    }
    const double &logp(int s, int t) const {
      return state_log_trans[s * state_count + t];
    }

    double &p(int s, int t) {
      return state_log_trans[s * state_count + t];
    }
    const double &p(int s, int t) const {
      return state_trans[s * state_count + t];
    }

    int randomTransition(int s) const {
      double *tp = state_trans + s * state_count;
      double r = random() / double(RAND_MAX);
      for (int i = 0; i < state_count; i++) {
        r -= tp[i];
        if (r <= 0.0) {
          return i;
        }
      }
      return state_count - 1;
    }

    std::vector<int> generate() const {
      int state = 0;
      std::vector<int> result;
      while (1) {
        state = randomTransition(state);
        if (state == state_count - 1) {
          break;
        }

        int d = states[state]->generate(result);
        DEBUG(5,
              std::cerr << "(" << state << "," << d << ")";);
      }
      DEBUG(5,
            std::cerr << std::endl;);
      return result;
    }
  };

  class Parse : public virtual RefObj {
    Parse(const Parse &);
    Parse &operator=(const Parse &);

  protected:
    double *a, *b, *d;
    Traceback::Ptr *p;
    int *s;
    int parse_length;
    int state_count;
    int offset;

  public:
    typedef Ref<Parse> Ptr;

    Parse() : a(NULL), b(NULL), d(NULL), p(NULL), s(NULL), parse_length(0), state_count(0), offset(0) {
    }

    ~Parse() {
      if (a) delete [] a;
      if (b) delete [] b;
      if (d) delete [] d;
      if (p) delete [] p;
      if (s) delete [] s;
    }

    int idx(int state, int pos) const {
      return (pos + offset) * state_count + state;
    }

    double &alpha(int state, int pos)                   { return a[idx(state, pos)]; }
    double &beta(int state, int pos)                    { return b[idx(state, pos)]; }
    double &delta(int state, int pos)                   { return d[idx(state, pos)]; }
    Traceback::Ptr &psi(int state, int pos)             { return p[idx(state, pos)]; }
    int &seq(int pos)                                   { return s[pos + offset];    }

    const double &alpha(int state, int pos) const       { return a[idx(state, pos)]; }
    const double &beta(int state, int pos) const        { return b[idx(state, pos)]; }
    const double &delta(int state, int pos) const       { return d[idx(state, pos)]; }
    const Traceback::Ptr &psi(int state, int pos) const { return p[idx(state, pos)]; }
    const int &seq(int pos) const                       { return s[pos + offset];    }

    void traceback() {
      std::cerr << "final result:" << std::endl;
      for (int state = 1; state < state_count - 1; state++) {
        std::cerr << "  state:" << state << " delta:" << delta(state, 0) << " alpha:" << alpha(state, 0) << "(" << exp(alpha(state, 0)) << ")" << std::endl;
        Traceback::Ptr p = psi(state, 0);
        while (p != NULL) {
          std::cerr << "    " << p->state << "," << p->length << std::endl;
          p = p->prev;
        }
      }
    }

    Traceback::Ptr linkState(int state, int length, int prev) {
      Traceback::Ptr result(NULL);
      if (state != prev || length != 0) {
        Traceback::Ptr &p(psi(prev, -length));
        if (state == prev) {
          result = new Traceback(p->prev, state, p->length + length);
        } else {
          result = new Traceback(p, state, length);
        }
      }
      return result;
    }

    template<typename RandomAccessIterator>
    void parse(const Model::Ptr &model, RandomAccessIterator begin, RandomAccessIterator end) {
      RandomAccessIterator pos;
      const Model &modelRef(*model);
      parse_length = (end - begin) + 1;
      state_count = model->stateCount();

      if (a) delete [] a;
      if (b) delete [] b;
      if (d) delete [] d;
      if (p) delete [] p;
      if (s) delete [] s;

      a = new double[parse_length * state_count];
      b = new double[parse_length * state_count];
      d = new double[parse_length * state_count];
      p = new Traceback::Ptr[parse_length * state_count];
      s = new int[parse_length];

      offset = 0;

      delta(0, 0) = 0.0;
      for (int i = 1; i < state_count; i++) delta(i, 0) = MATH::LOG_ZERO;
      for (int i = 1; i < parse_length; i++) delta(0, i) = MATH::LOG_ZERO;

      alpha(0, 0) = 0.0;
      for (int i = 1; i < state_count; i++) alpha(i, 0) = MATH::LOG_ZERO;
      for (int i = 1; i < parse_length; i++) alpha(0, i) = MATH::LOG_ZERO;

      for (pos = begin; pos != end;) {
        ++offset;
        seq(0) = *pos++;

        DEBUG(8,
              std::cerr << std::endl << std::endl << "POS: " << pos - begin - 1 << std::endl;
              std::cerr << "offset=" << offset << " ch=" << seq(0) << std::endl;);

        for (int j = 1; j < state_count - 1; ++j) {
          const StateBase *js(model->state(j).ptr());
          double alpha_j, delta_j;
          int prev_state_j, state_length_j;

          js->alphaDelta(j, modelRef, *this, pos - begin, alpha_j, delta_j, prev_state_j, state_length_j);
          DEBUG(9,
                std::cerr << "delta_j=" << delta_j << " prev_state_j=" << prev_state_j << " state_length_j=" << state_length_j << std::endl;);
          psi(j, 0) = linkState(j, state_length_j, prev_state_j);
          delta(j, 0) = delta_j;
          alpha(j, 0) = alpha_j;
        }
        DEBUG(8,
              std::cerr << std::flush;
              for (int j = 1; j < state_count - 1; ++j) {
                Traceback::Ptr tb = psi(j, 0);
                if (tb == NULL) {
                  fprintf(stderr, "(null) ");
                } else {
                  fprintf(stderr, "%8.4f(%d:%4d) ", delta(j, 0), psi(j, 0)->state, psi(j,0)->length);
                }
              }
              fprintf(stderr, "\n");
              fflush(stderr);
              std::cerr << std::endl;);
      }

#if 0
      for (int i = 0; i < state_count - 1; i++) beta(i, 0) = MATH::LOG_ZERO;
      beta(state_count - 1, 0) = 0.0;

      for (pos = end; pos != begin;) {
        --offset;
        seq(0) = *--pos;
        DEBUG(8,
              std::cerr << std::endl << std::endl << "POS: " << pos - begin << std::endl;
              std::cerr << "offset=" << offset << " ch=" << *pos << std::endl;);

        for (int j = 1; j < state_count - 1; ++j) {
          const StateBase *js(model->state(j).ptr());
          double beta_j;

          js->beta(j, modelRef, *this, end - pos, beta_j);
          DEBUG(9,
                std::cerr << "beta_j=" << beta_j << std::endl;);
          beta(j, 0) = beta_j;
        }
      }

      offset = end - begin;
#endif

      DEBUG(3,
            traceback(););
    }
  };

  class ModelBuilder : public virtual RefObj {
    ModelBuilder(const ModelBuilder &);
    ModelBuilder &operator=(const ModelBuilder &);

  protected:
    std::vector<std::pair<std::string, StateBase::Ptr> > states;
    std::map<std::string, int> state_name_map;
    std::map<std::pair<int, int>, double> state_trans_map;

    int stateNum(const std::string &s);
  public:
    typedef Ref<ModelBuilder> Ptr;

    ModelBuilder();
    virtual ~ModelBuilder();

    virtual void addState(const std::string &, const StateBase::Ptr &);
    virtual void removeState(const std::string &);
    virtual void addStateTransition(const std::string &, const std::string &, double);
    virtual Model::Ptr make();
  };

  template<typename Distrib, typename Emitter>
  class State : public StateBase, public Distrib, public Emitter {
  public:
    typedef Ref<State> Ptr;

    State() : StateBase(), Distrib(), Emitter() {
    }
    State(const Distrib &d, const Emitter &e) : StateBase(), Distrib(d), Emitter(e) {
    }
    State(const Distrib &d) : StateBase(), Distrib(d), Emitter() {
    }
    State(const Emitter &e) : StateBase(), Distrib(), Emitter(e) {
    }

    State(const State &s) : StateBase(s), Distrib(s), Emitter(s) {
    }

    State &operator=(const State &s) {
      if (this != &s) {
        StateBase::operator=(s);
        Distrib::operator=(s);
        Emitter::operator=(s);
      }
      return *this;
    }
    void setEmissionDistrib(const Emitter &e) {
      Emitter::operator=(e);
    }
    void setLengthDistrib(const Distrib &d) {
      Distrib::operator=(d);
    }
    virtual void      delta(int j, const Model &model, const Parse &parse, int max_len, double &delta, int &prev_state, int &state_length) const;
    virtual void      alpha(int j, const Model &model, const Parse &parse, int max_len, double &alpha) const;
    virtual void alphaDelta(int j, const Model &model, const Parse &parse, int max_len, double &alpha, double &delta, int &prev_state, int &state_length) const;
    virtual void       beta(int j, const Model &model, const Parse &parse, int max_len, double &alpha) const;

    virtual int generate(std::vector<int> &result) const {
      int d;
      Emitter::randSequence(result, d = Distrib::randLength());
      return d;
    }
  };

  template<typename Distrib, typename Emitter>
  void State<Distrib, Emitter>::delta(int j, const Model &model, const Parse &parse, int max_len, double &delta, int &prev_state, int &state_length) const {
    const std::vector<int> &pred(model.predStates(j));

    int lmin = std::min(max_len + 1, State<Distrib, Emitter>::minLength());
    int lmax = std::min(max_len + 1, State<Distrib, Emitter>::maxLength());
    int d;
    double dprob, eprob, sprob;

    delta = MATH::LOG_ZERO;
    prev_state = j;
    state_length = 0;

    DEBUG(9,
          std::cerr << "state:" << j << std::endl;);

    typename Emitter::template Generator<-1> g(static_cast<const Emitter *>(this), &(parse.seq(0)), lmin, lmax);

    while (g.gen(d, eprob)) {
      dprob = State<Distrib, Emitter>::logpLength(d);

      DEBUG(10,
            std::cerr << "  d=" << d << " eprob=" << eprob << " dprob=" << dprob << std::endl;);

      sprob = eprob + dprob;

      for (int _i = pred.size() - 1; _i >= 0; --_i) {
        int i = pred[_i];
        double dp = sprob + model.logp(i, j) + parse.delta(i, -d);
        DEBUG(11,
              std::cerr << " model.logp(" << i << "," << j << ") =" << model.logp(i, j)
                        << " delta(" << i << "," << -d << ")=" << parse.delta(i, -d)
                        << "    --> dp=" << dp
                        << std::endl;);
        if (dp >= delta) {
          delta = dp;
          prev_state = i;
          state_length = d;
        }
      }
    }

    if (delta == MATH::LOG_ZERO) {
      prev_state = j;
      state_length = 0;
    }

    DEBUG(9,
          std::cerr << "  delta=" << delta
                    << " prev_state=" << prev_state
                    << " state_length=" << state_length
                    << std::endl;);
  }

  template<typename Distrib, typename Emitter>
  void State<Distrib, Emitter>::alpha(int j, const Model &model, const Parse &parse, int max_len, double &alpha) const {
    const std::vector<int> &pred(model.predStates(j));

    int lmin = std::min(max_len + 1, State<Distrib, Emitter>::minLength());
    int lmax = std::min(max_len + 1, State<Distrib, Emitter>::maxLength());
    int d;
    double dprob, eprob, sprob;

    alpha = MATH::LOG_ZERO;

    DEBUG(9,
          std::cerr << "state:" << j << std::endl;);

    typename Emitter::template Generator<-1> g(static_cast<const Emitter *>(this), &(parse.seq(0)), lmin, lmax);

    while (g.gen(d, eprob)) {
      dprob = State<Distrib, Emitter>::logpLength(d);

      DEBUG(10,
            std::cerr << "  d=" << d << " eprob=" << eprob << " dprob=" << dprob << std::endl;);

      sprob = eprob + dprob;

      if (sprob != MATH::LOG_ZERO) {
        double temp = MATH::LOG_ZERO;
        for (int _i = pred.size() - 1; _i >= 0; --_i) {
          int i = pred[_i];
          double ap = model.logp(i, j) + parse.alpha(i, -d);

          DEBUG(11,
                std::cerr << " model.logp(" << i << "," << j << ") =" << model.logp(i, j)
                          << " alpha(" << i << "," << -d << ")=" << parse.alpha(i, -d)
                          << "    --> ap=" << ap
                          << std::endl;);


          if (ap != MATH::LOG_ZERO) {
            DEBUG(12,
                  std::cerr << "logAdd(" << alpha << "," << ap << ") = " << MATH::logAdd(ap, alpha) << std::endl;);
            temp = MATH::logAdd(ap, temp);
          }
        }
        if (temp != MATH::LOG_ZERO) {
          alpha = MATH::logAdd(temp + sprob, alpha);
        }
      }
    }

    DEBUG(9,
          std::cerr << " alpha=" << alpha << std::endl;);
  }

  template<typename Distrib, typename Emitter>
  void State<Distrib, Emitter>::alphaDelta(int j, const Model &model, const Parse &parse, int max_len, double &alpha, double &delta, int &prev_state, int &state_length) const {
    const std::vector<int> &pred(model.predStates(j));

    int lmin = std::min(max_len + 1, State<Distrib, Emitter>::minLength());
    int lmax = std::min(max_len + 1, State<Distrib, Emitter>::maxLength());
    int d;
    double dprob, eprob, sprob;

    alpha = MATH::LOG_ZERO;
    delta = MATH::LOG_ZERO;
    prev_state = j;
    state_length = 0;

    DEBUG(9,
          std::cerr << "state:" << j << " (" << model.stateName(j) << ")" << std::endl;);

    typename Emitter::template Generator<-1> g(static_cast<const Emitter *>(this), &(parse.seq(0)), lmin, lmax);

    while (g.gen(d, eprob)) {
      dprob = State<Distrib, Emitter>::logpLength(d);
      DEBUG(10,
            std::cerr << "  d=" << d << " eprob=" << eprob << " dprob=" << dprob << std::endl;);
      sprob = eprob + dprob;
      if (sprob == MATH::LOG_ZERO) continue;

      for (int _i = pred.size() - 1; _i >= 0; --_i) {
        int i = pred[_i];
        double dp = sprob + model.logp(i, j) + parse.delta(i, -d);
        DEBUG(11,
              std::cerr << " model.logp(" << i << "," << j << ") =" << model.logp(i, j)
                        << " delta(" << i << "," << -d << ")=" << parse.delta(i, -d)
                        << "    --> dp=" << dp
                        << std::endl;);
        if (dp >= delta) {
          delta = dp;
          prev_state = i;
          state_length = d;
        }
      }

      if (sprob != MATH::LOG_ZERO) {
        double temp = MATH::LOG_ZERO;
        for (int _i = pred.size() - 1; _i >= 0; --_i) {
          int i = pred[_i];
          double ap = model.logp(i, j) + parse.alpha(i, -d);

          DEBUG(11,
                std::cerr << " model.logp(" << i << "," << j << ") =" << model.logp(i, j)
                          << " alpha(" << i << "," << -d << ")=" << parse.alpha(i, -d)
                          << "    --> ap=" << ap
                          << std::endl;);


          if (ap != MATH::LOG_ZERO) {
            DEBUG(12,
                  std::cerr << "logAdd(" << alpha << "," << ap << ") = " << MATH::logAdd(ap, alpha) << std::endl;);
            temp = MATH::logAdd(ap, temp);
          }
        }
        if (temp != MATH::LOG_ZERO) {
          alpha = MATH::logAdd(temp + sprob, alpha);
        }
      }
    }

    if (delta <= MATH::LOG_ZERO) {
      delta = MATH::LOG_ZERO;
      prev_state = j;
      state_length = 0;
    }

    DEBUG(9,
          std::cerr << "  delta=" << delta
                    << " alpha=" << alpha
                    << " prev_state=" << prev_state
                    << " state_length=" << state_length
                    << std::endl;);
  }

  template<typename Distrib, typename Emitter>
  void State<Distrib, Emitter>::beta(int j, const Model &model, const Parse &parse, int max_len, double &beta) const {
    const std::vector<int> &succ(model.succStates(j));

    int lmin = std::min(max_len + 1, State<Distrib, Emitter>::minLength());
    int lmax = std::min(max_len + 1, State<Distrib, Emitter>::maxLength());
    int d;
    double dprob, eprob, sprob;

    beta = MATH::LOG_ZERO;

    DEBUG(9,
          std::cerr << "state:" << j << std::endl;);

    typename Emitter::template Generator<+1> g(static_cast<const Emitter *>(this), &(parse.seq(0)), lmin, lmax);

    while (g.gen(d, eprob)) {
      dprob = State<Distrib, Emitter>::logpLength(d);

      DEBUG(10,
            std::cerr << "  d=" << d << " eprob=" << eprob << " dprob=" << dprob << std::endl;);

      sprob = eprob + dprob;

      if (sprob != MATH::LOG_ZERO) {
        double temp = MATH::LOG_ZERO;
        for (int _i = succ.size() - 1; _i >= 0; --_i) {
          int i = succ[_i];
          double bp = model.logp(j, i) + parse.beta(i, +d);

          DEBUG(11,
                std::cerr << " model.logp(" << i << "," << j << ") =" << model.logp(i, j)
                          << " beta(" << i << "," << -d << ")=" << parse.beta(i, -d)
                          << "    --> bp=" << bp
                          << std::endl;);

          if (bp != MATH::LOG_ZERO) {
            DEBUG(12,
                  std::cerr << "logAdd(" << beta << "," << bp << ") = " << MATH::logAdd(bp, beta) << std::endl;);
            temp = MATH::logAdd(bp, temp);
          }
        }
        if (temp != MATH::LOG_ZERO) {
          beta = MATH::logAdd(temp + sprob, beta);
        }
      }
    }

    DEBUG(9,
          std::cerr << " beta=" << beta << std::endl;);
  }

}

#include <GHMM/ghmm_util.hh>

#endif
