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

#include <GHMM/string_funcs.hh>
#include <GHMM/ghmm.hh>

int main(int argc, char **argv) {
  GHMM::UTIL::Alphabet::Ptr alphabet = new GHMM::UTIL::Alphabet();
  alphabet->addToken("heads");
  alphabet->addToken("tails");
  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);

  GHMM::EMISSION::Base::Ptr heads, tails;
  GHMM::LENGTH::Base::Ptr heads_length, tails_length;

  heads = new GHMM::EMISSION::Stateless(ep->parse("heads: 9 tails: 1"));
  tails = new GHMM::EMISSION::Stateless(ep->parse("tails: 9 heads: 1"));

  heads_length = new GHMM::LENGTH::Geometric(10);
  tails_length = new GHMM::LENGTH::Fixed(10);

  GHMM::StateBase::Ptr heads_state = GHMM::UTIL::makeState(heads_length, heads);
  GHMM::StateBase::Ptr tails_state = GHMM::UTIL::makeState(tails_length, tails);

  GHMM::ModelBuilder mb;
  mb.addState("heads", heads_state);
  mb.addState("tails", tails_state);

  mb.addStateTransition(GHMM::Model::BEGIN, "heads",          1);
  mb.addStateTransition(GHMM::Model::BEGIN, "tails",          1);
  mb.addStateTransition("heads",            "tails",          9);
  mb.addStateTransition("heads",            GHMM::Model::END, 1);
  mb.addStateTransition("tails",            "heads",          9);
  mb.addStateTransition("tails",            GHMM::Model::END, 1);

  GHMM::Model::Ptr model = mb.make();

  int a = strtoul(argv[1], NULL, 10);
  int b = strtoul(argv[2], NULL, 10);
  int c = strtoul(argv[3], NULL, 10);

  for (int i = a; i <= b; i += c) {
    srandom(i);

    std::vector<int> out = model->generate();
    for (int i = 0; i < (int)out.size(); i++) {
      std::cout << alphabet->token(out[i]) << " ";
    }
    std::cout << std::endl;

    GHMM::Parse::Ptr parse = new GHMM::Parse();

    parse->parse(model, out.begin(), out.end());
  }
}
