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
#ifndef GHMM_GHMM_UTIL_HH_INCLUDED
#define GHMM_GHMM_UTIL_HH_INCLUDED

#include <string>
#include <map>
#include <vector>
#include <sstream>
#include <GHMM/string_funcs.hh>

namespace GHMM {
  namespace UTIL {
    class Alphabet : public virtual RefObj {
      std::map<std::string, int> token_map;
      std::vector<std::string> tokens;

    public:
      static const std::string NO_TOKEN;

      typedef Ref<Alphabet> Ptr;

      Alphabet() : token_map() {
      }
      Alphabet &operator=(const Alphabet &a) {
        if (this != &a) {
          token_map = a.token_map;
        }
        return *this;
      }
      Alphabet(const Alphabet &a) : token_map() {
        *this = a;
      }
      ~Alphabet() {
      }
      int addToken(const std::string &t) {
        std::map<std::string, int>::iterator i = token_map.find(t);
        if (i != token_map.end()) {
          return (*i).second;
        }
        int idx = tokens.size();
        token_map[t] = idx;
        tokens.push_back(t);
        return idx;
      }
      void addCharTokenRange(char a, char b) {
        while (a <= b) {
          addToken(std::string(1, a++));
        }
      }
      int size() {
        return token_map.size();
      }
      const std::string &token(int idx) {
        if (idx < 0 || idx >= (int)tokens.size()) return NO_TOKEN;
        return tokens[idx];
      }
      int index(const std::string &t) {
        std::map<std::string, int>::iterator i = token_map.find(t);
        if (i != token_map.end()) {
          return (*i).second;
        }
        return -1;
      }
    };

    class EmissionDistributionParser : public virtual RefObj {
      EmissionDistributionParser();

    protected:
      Alphabet::Ptr alphabet;

    public:
      typedef Ref<EmissionDistributionParser> Ptr;

      EmissionDistributionParser(const Alphabet::Ptr &a) : alphabet(a) {
      }
      EmissionDistributionParser &operator=(const EmissionDistributionParser &e) {
        if (this != &e) {
          alphabet = e.alphabet;
        }
        return *this;
      }
      EmissionDistributionParser(const EmissionDistributionParser &e) : alphabet(NULL) {
        *this = e;
      }
      virtual ~EmissionDistributionParser() {
      }

      virtual MATH::DPDF::Ptr parse(std::istream &in) {
        MATH::DPDF::Ptr result = new MATH::DPDF();
        result->setDistrib(alphabet->size(), 0.0);
        while(!in.eof()) {
          std::string tok;
          double freq;
          while (1) {
            int ch = in.get();
            if (ch == EOF) return NULL;
            if (ch == ':') break;
            tok += ch;
          }
          tok = string::strip(tok);
          in >> freq;
          int idx = alphabet->index(tok);
          if (idx != -1) {
            result->setp(idx, freq);
          }
          if (in.fail()) return NULL;
        }
        result->normalize();
        return result;
      }
      virtual MATH::DPDF::Ptr parse(const std::string &s) {
        std::istringstream in(s);
        return parse(in);
      }
      virtual MATH::DPDF::Ptr parse(const char *s) {
        std::istringstream in(s);
        return parse(in);
      }
    };

    GHMM::StateBase::Ptr makeState(GHMM::LENGTH::Base::Ptr length, GHMM::EMISSION::Base::Ptr emission);
    GHMM::StateBase::Ptr makeMotifState(const std::vector<MATH::DPDF::Ptr> &pssm);
  }
}

#endif
