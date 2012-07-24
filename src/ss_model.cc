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
#include "predict_pexel.hh"

#if VERSION == 1
static int leader_raw_distrib[] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4,
  4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 8, 8, 8, 9, 9, 10, 10, 11,
  11, 11, 12, 12, 12, 12, 13, 13, 14, 14, 15, 15, 15, 15, 16, 16, 16,
  17, 17, 19, 19, 20, 20, 20, 20, 23, 24, 25, 27, 27, 27, 28, 28, 28,
  29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 32, 33, 35, 37, 37, 38,
  39, 42, 42, 43, 43, 43, 44, 44, 47, 47, 48, 48, 48, 49, 49, 49, 49,
  49, 49, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52,
  52, 52, 52, 52, 53, 53, 53, 53, 53, 55, 55, 56, 56, 59, 63, 78, 92,
  118, 143
};
#endif

#if VERSION == 2
static int leader_raw_distrib[] = {
  1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 4, 4,
  5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  6, 6, 7, 7, 7, 7, 7, 7, 7, 7, 8, 8, 9, 10, 10, 11, 11, 12, 12, 12,
  12, 13, 13, 13, 13, 14, 15, 16, 16, 16, 17, 18, 19, 20, 20, 21, 23,
  24, 25, 26, 26, 28, 28, 28, 28, 28, 28, 29, 29, 29, 29, 29, 30, 35,
  37, 37, 39, 39, 40, 42, 43, 43, 43, 43, 44, 44, 45, 46, 47, 47, 48,
  48, 48, 48, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 50, 51,
  51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 52, 53, 53, 54,
  54, 54, 57, 58, 63
};

static int hydrophobic_raw_distrib[] = {
  10, 10, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 12, 12, 12,
  12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13,
  13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14,
  14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14,
  14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16,
  16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 17, 17, 17, 17,
  17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19,
  19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 21, 22, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 25,
  25
};

#endif

std::pair<std::string, std::string> makeSSModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet) {
  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);

  GHMM::LENGTH::Discrete::Ptr leader_length;
  GHMM::EMISSION::Base::Ptr ss_leader, hydrophobic_distrib;

#if VERSION == 1
  GHMM::LENGTH::Uniform::Ptr hydrophobic_length;
#endif

#if VERSION == 2
  GHMM::LENGTH::Discrete::Ptr hydrophobic_length;
#endif

  leader_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(2.0),
                                                            leader_raw_distrib,
                                                            ENDOF(leader_raw_distrib),
                                                            1,
                                                            80));

#if VERSION == 1
  hydrophobic_length = b_hydrophobic_length = new GHMM::LENGTH::Uniform(9.5, 25.5);

  ss_leader = new GHMM::EMISSION::Stateless(ep->parse("\n\
             A:  47 C: 105 D:  96 E: 145\n\
             F: 223 G:  85 H:  58 I: 259\n\
             K: 548 L: 213 M: 197 N: 448\n\
             P:  46 Q:  62 R: 166 S: 331\n\
             T: 118 V:  94 W:  21 Y: 236"));

  hydrophobic_distrib = b_hydrophobic_distrib = new GHMM::EMISSION::Stateless(ep->parse("\n\
             A:  34 C:  75 D:   5 E:  10\n\
             F: 289 G:  57 H:  12 I: 351\n\
             K:  26 L: 446 M:  30 N:  57\n\
             P:  18 Q:  12 R:  10 S: 106\n\
             T:  77 V: 172 W:  14 Y: 104"));
#endif

#if VERSION == 2
  hydrophobic_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(1.0),
                                                                 hydrophobic_raw_distrib,
                                                                 ENDOF(hydrophobic_raw_distrib),
                                                                 10,
                                                                 25));

  ss_leader = new GHMM::EMISSION::Stateless(ep->parse("\n\
             W:  20 P:  44 A:  52 Q:  73 H:  74 G:  90\n\
             D: 103 V: 109 M: 111 C: 116 T: 127 E: 156\n\
             R: 191 L: 239 F: 248 Y: 274 I: 285 S: 366\n\
             N: 516 K: 639"));

  hydrophobic_distrib = new GHMM::EMISSION::Stateless(ep->parse("\n\
             D:   4 E:   6 Q:   6 R:   6 H:  11 K:  24\n\
             P:  24 W:  24 M:  27 A:  37 N:  40 G:  70\n\
             T:  87 C: 106 S: 113 Y: 132 V: 199 F: 396\n\
             I: 508 L: 559"));
#endif

  GHMM::StateBase::Ptr leader = GHMM::UTIL::makeState(leader_length, ss_leader);
  GHMM::StateBase::Ptr hydrophobic = GHMM::UTIL::makeState(hydrophobic_length, hydrophobic_distrib);

  mb.addState("a-leader",      leader);
  mb.addState("a-hydrophobic", hydrophobic);

  mb.addStateTransition("a-leader", "a-hydrophobic",  1);

  return std::make_pair(std::string("a-leader"), std::string("a-hydrophobic"));

}

