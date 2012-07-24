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
#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <fstream>

#include <GHMM/string_funcs.hh>
#include <GHMM/ghmm.hh>

std::pair<std::string, std::string> makeSignalPModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);
std::pair<std::string, std::string> makeSSModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);

int main(int argc, char **argv) {
  GHMM::UTIL::Alphabet::Ptr alphabet = new GHMM::UTIL::Alphabet();
  alphabet->addCharTokenRange('A','Z');
  GHMM::ModelBuilder mb;

  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);
  GHMM::StateBase::Ptr a_met = GHMM::UTIL::makeState(NULL, new GHMM::EMISSION::Stateless(ep->parse("M: 1")));
  mb.addState("met", a_met);
  mb.addStateTransition(GHMM::Model::BEGIN, "met", 1);

  std::pair<std::string, std::string> ss_states;

#ifdef SIGNALP_MODEL
  ss_states = makeSignalPModel(mb, alphabet);
#else
  ss_states = makeSSModel(mb, alphabet);
#endif

  mb.addStateTransition("met", ss_states.first, 1);
  mb.addStateTransition(ss_states.second, GHMM::Model::END, 1);

  GHMM::Model::Ptr model = mb.make();

  // XXX: handle this with autoconf
#if defined(__APPLE__)
  srandomdev();
#else
  srandom(time(NULL));
#endif

  for (int i = 0, l = strtoul(argv[1], NULL, 10); i < l; i++) {
    std::vector<int> temp = model->generate();
    for (int j = 0; j < (int)temp.size(); j++) {
      std::cout << (char)('A' + temp[j]);
    }
    std::cout << std::endl;
  }
}
