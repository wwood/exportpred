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

#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <fstream>

#include <GHMM/string_funcs.hh>

#include <getopt.h>

#define BLOCK_SIZE 1024

template<typename OutputIterator>
static inline void addSeq(OutputIterator &out,
                          const std::string &name,
                          const std::vector<char *> &seq,
                          int seq_len) {
  std::string s;
  s.reserve(seq_len);

  for (int i = 0; i < seq_len; i += BLOCK_SIZE) {
    s.append(seq[i / BLOCK_SIZE], std::min(BLOCK_SIZE, seq_len - i));
  }
  assert(s.size() ==(unsigned int)seq_len);
  *out++ = std::make_pair(name, s);
}

static inline void addChar(std::vector<char *> &seq,
                           int &seq_len,
                           int &seq_size,
                           char c) {
  if (!isspace(c)) {
    while (seq_len >= seq_size) {
      seq.push_back(new char[BLOCK_SIZE]);
      seq_size += BLOCK_SIZE;
    }
    seq[seq_len / BLOCK_SIZE][seq_len % BLOCK_SIZE] = c;
    seq_len++;
  }
}

template<typename OutputIterator>
void readFasta(std::istream &in,
               OutputIterator out,
               std::string (*name_xform)(const std::string &)) {
  std::string name = "";
  int c;
  std::vector<char *> seq;

  std::string l;
  int seq_len = 0;
  int seq_size = 0;
  int state = 0;

  while ((c = in.get()) != EOF) {
    switch (state) {
    case 0: {
      // looking for name
      if (c == '>') {
        state = 1;
        name = "";
      }
      break;
    }
    case 1: {
      // name
      if (c == '\n') {
        state = 2;
        if (name_xform) name = name_xform(name);
      } else {
        name.push_back(c);
      }
      break;
    }
    case 2: {
      // sequence, first char of column
      if (c == '>') {
        if (seq_len) addSeq(out, name, seq, seq_len);
        seq_len = 0;
        state = 1;
        name = "";
      } else if (c != '\n') {
        addChar(seq, seq_len, seq_size, c);
        state = 3;
      }
      break;
    }
    case 3: {
      // sequence, not first char of column
      if (c == '\n') {
        state = 2;
      } else {
        addChar(seq, seq_len, seq_size, c);
      }
      break;
    }
    }
  }

  if (seq_len) addSeq(out, name, seq, seq_len);

  for (std::vector<char *>::iterator i = seq.begin(); i != seq.end(); ++i) {
    delete [] *i;
  }
}



#if VERSION == 1

static int a_spacer_raw_distrib[] = {
  1, 2, 3, 6, 8, 8, 9, 10, 10, 11, 11, 12, 13, 13, 13, 13, 14, 14,
  14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18,
  18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
  19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21,
  21, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 24, 24,
  24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 25, 26, 26,
  26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27,
  27, 28, 29, 29, 29, 29, 29, 30, 31, 31, 31, 31, 31, 32, 32, 33, 35,
  37, 41, 43, 48, 53, 61, 89, 110, 129, 166, 182, 202, 209, 324, 495
};

#endif

#if VERSION == 2

static int a_spacer_raw_distrib[] = {
  9, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 16,
  16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
  18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
  19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21,
  21, 21, 21, 21, 21, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23,
  24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25,
  25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
  26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 29, 31, 31, 31, 31, 33,
  33
};

static int b_hydrophobic_raw_distrib[] = {
  22, 22, 23, 23, 23, 23, 23, 23, 25, 25, 25, 25, 25, 25, 25, 25, 25,
  25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27,
  27, 27, 27, 27, 28, 29, 29, 29, 29, 29, 29, 29, 32, 37
};

#endif

static int b_leader_raw_distrib[] = {
  9, 10, 11, 11, 11, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 13,
  13, 13, 14, 14, 14, 14, 14, 14, 15, 15, 15, 16, 16, 16, 16, 16, 17,
  17, 17, 17, 17, 17, 17, 17, 18, 19, 19, 19, 19, 19, 19, 20, 20, 20,
  21, 21, 21, 22, 23, 24, 25
};

GHMM::Model::Ptr makePEXELmodel() {
  GHMM::UTIL::Alphabet::Ptr alphabet = new GHMM::UTIL::Alphabet();
  alphabet->addCharTokenRange('A','Z');
  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);

  GHMM::EMISSION::Base::Ptr background, b_hydrophobic_distrib, met;

#ifdef HALDAR_MOTIF
  std::vector<MATH::DPDF::Ptr> RLE(11, MATH::DPDF::Ptr());
#else
  std::vector<MATH::DPDF::Ptr> RLE(7, MATH::DPDF::Ptr());
#endif

  std::vector<MATH::DPDF::Ptr> KLD(7, MATH::DPDF::Ptr());

  GHMM::LENGTH::Discrete::Ptr a_leader_length, a_spacer_length, b_leader_length;

  GHMM::LENGTH::Geometric::Ptr a_tail_length, b_spacer_length, b_tail_length, c_tail_length;

#if VERSION == 1
  GHMM::LENGTH::Uniform::Ptr a_hydrophobic_length, b_hydrophobic_length;
#endif

#if VERSION == 2
  GHMM::LENGTH::Discrete::Ptr a_hydrophobic_length, b_hydrophobic_length;
#endif

  a_spacer_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(1.0),
                                                            a_spacer_raw_distrib,
                                                            ENDOF(a_spacer_raw_distrib),
                                                            1,
                                                            60));

  b_leader_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(1.0),
                                                            b_leader_raw_distrib,
                                                            ENDOF(b_leader_raw_distrib),
                                                            -1,
                                                            30));
  background = new GHMM::EMISSION::Stateless(ep->parse("\n\
             A:  78883 C:  71359 D: 260979 E: 288230\n\
             F: 175488 G: 114068 H:  97688 I: 373389\n\
             K: 473828 L: 304967 M:  88773 N: 581084\n\
             P:  80295 Q: 111860 R: 106760 S: 256676\n\
             T: 164816 V: 154088 W:  19966 Y: 230362"));


#if VERSION == 1
  b_hydrophobic_length = new GHMM::LENGTH::Uniform(9.5, 25.5);

  b_hydrophobic_distrib = new GHMM::EMISSION::Stateless(ep->parse("\n\
             A:  34 C:  75 D:   5 E:  10\n\
             F: 289 G:  57 H:  12 I: 351\n\
             K:  26 L: 446 M:  30 N:  57\n\
             P:  18 Q:  12 R:  10 S: 106\n\
             T:  77 V: 172 W:  14 Y: 104"));
#endif

#if VERSION == 2
  b_hydrophobic_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(1.0),
                                                                 b_hydrophobic_raw_distrib,
                                                                 ENDOF(b_hydrophobic_raw_distrib),
                                                                 20,
                                                                 35));
  
  b_hydrophobic_distrib = new GHMM::EMISSION::Stateless(ep->parse("\n\
             D:   4 E:   6 Q:   6 R:   6 H:  11 K:  24\n\
             P:  24 W:  24 M:  27 A:  37 N:  40 G:  70\n\
             T:  87 C: 106 S: 113 Y: 132 V: 199 F: 396\n\
             I: 508 L: 559"));
#endif

  met = new GHMM::EMISSION::Stateless(ep->parse("M: 1"));

#ifdef HALDAR_MOTIF

  RLE[ 0] = ep->parse("\
    A: 0.016425 C: 0.000540 D: 0.014277 E: 0.027577 F: 0.001233 G: 0.022156 H: 0.002540\n\
    I: 0.203500 K: 0.104854 L: 0.070546 M: 0.001576 N: 0.062469 P: 0.030850 Q: 0.015745\n\
    R: 0.235895 S: 0.093319 T: 0.053231 V: 0.041329 W: 0.000384 Y: 0.001552");
  RLE[ 1] = ep->parse("\
    A: 0.035110 C: 0.019491 D: 0.032513 E: 0.054565 F: 0.090540 G: 0.012293 H: 0.069190\n\
    I: 0.077420 K: 0.092437 L: 0.032299 M: 0.011190 N: 0.145192 P: 0.020841 Q: 0.052481\n\
    R: 0.005220 S: 0.063105 T: 0.014458 V: 0.022206 W: 0.000765 Y: 0.148686");
  RLE[ 2] = ep->parse("\
    A: 0.006767 C: 0.010168 D: 0.023768 E: 0.027373 F: 0.087736 G: 0.022124 H: 0.012096\n\
    I: 0.069334 K: 0.152723 L: 0.080357 M: 0.020839 N: 0.100823 P: 0.030838 Q: 0.015603\n\
    R: 0.053420 S: 0.218135 T: 0.014804 V: 0.031933 W: 0.000390 Y: 0.020770");
  RLE[ 3] = ep->parse("\
    A: 0.000059 C: 0.000034 D: 0.000070 E: 0.000039 F: 0.000041 G: 0.000146 H: 0.000044\n\
    I: 0.000023 K: 0.000048 L: 0.000068 M: 0.000012 N: 0.000046 P: 0.000085 Q: 0.000041\n\
    R: 0.999034 S: 0.000061 T: 0.000053 V: 0.000033 W: 0.000024 Y: 0.000038");
  RLE[ 4] = ep->parse("\
    A: 0.004727 C: 0.011067 D: 0.000614 E: 0.001160 F: 0.032319 G: 0.001406 H: 0.000859\n\
    I: 0.400985 K: 0.001370 L: 0.185451 M: 0.032715 N: 0.001157 P: 0.000921 Q: 0.001221\n\
    R: 0.001350 S: 0.179883 T: 0.062928 V: 0.057912 W: 0.010406 Y: 0.011552");
  RLE[ 5] = ep->parse("\
    A: 0.002464 C: 0.000796 D: 0.000626 E: 0.001260 F: 0.012464 G: 0.000814 H: 0.000721\n\
    I: 0.033532 K: 0.001374 L: 0.896824 M: 0.031588 N: 0.000832 P: 0.001065 Q: 0.001607\n\
    R: 0.001635 S: 0.001584 T: 0.002179 V: 0.006832 W: 0.000534 Y: 0.001270");
  RLE[ 6] = ep->parse("\
    A: 0.418578 C: 0.035794 D: 0.009649 E: 0.002415 F: 0.026626 G: 0.046949 H: 0.017247\n\
    I: 0.012737 K: 0.018596 L: 0.022344 M: 0.002607 N: 0.001777 P: 0.001035 Q: 0.001720\n\
    R: 0.010255 S: 0.215815 T: 0.086576 V: 0.034855 W: 0.000586 Y: 0.033838");
  RLE[ 7] = ep->parse("\
    A: 0.012139 C: 0.000031 D: 0.075676 E: 0.688634 F: 0.000098 G: 0.010772 H: 0.000417\n\
    I: 0.000363 K: 0.002178 L: 0.000631 M: 0.000277 N: 0.000983 P: 0.000688 Q: 0.163090\n\
    R: 0.000773 S: 0.031310 T: 0.001004 V: 0.010732 W: 0.000034 Y: 0.000168");
  RLE[ 8] = ep->parse("\
    A: 0.034662 C: 0.030320 D: 0.001834 E: 0.022856 F: 0.100029 G: 0.060467 H: 0.011140\n\
    I: 0.035027 K: 0.081740 L: 0.065428 M: 0.012476 N: 0.109695 P: 0.079388 Q: 0.031994\n\
    R: 0.002601 S: 0.072368 T: 0.092184 V: 0.075345 W: 0.010280 Y: 0.070167");
  RLE[ 9] = ep->parse("\
    A: 0.026047 C: 0.010128 D: 0.043112 E: 0.142809 F: 0.020407 G: 0.022169 H: 0.021753\n\
    I: 0.040243 K: 0.095314 L: 0.080063 M: 0.011154 N: 0.158453 P: 0.002068 Q: 0.121355\n\
    R: 0.024760 S: 0.102939 T: 0.005256 V: 0.012464 W: 0.000383 Y: 0.059126");
  RLE[10] = ep->parse("\
    A: 0.054834 C: 0.000533 D: 0.129472 E: 0.066039 F: 0.049196 G: 0.022168 H: 0.021752\n\
    I: 0.030652 K: 0.172076 L: 0.041685 M: 0.001559 N: 0.100878 P: 0.030855 Q: 0.006203\n\
    R: 0.015161 S: 0.141321 T: 0.053234 V: 0.060448 W: 0.000383 Y: 0.001552");

#else

#if VERSION == 1
  RLE[0] = ep->parse("A: 1 C: 6 D: 2 E: 2 F: 11 G: 7 H: 3 I: 15 K: 23 L: 6 M: 1 N: 15 P: 1 Q: 3 R: 6 S: 38 T: 5 V: 3 W: 1 Y: 6");
  RLE[1] = ep->parse("K: 3 R: 152");
  // RLE[2] = ep->parse("I: 28 K: 11 L: 15 N: 32 S: 26 T: 9");
  RLE[2] = ep->parse("A: 2 C: 6 E: 1 F: 4 H: 3 I: 28 K: 11 L: 15 M: 2 N: 32 Q: 3 R: 3 S: 26 T: 9 V: 6 W: 1 Y: 3");
  // RLE[3] = ep->parse("L: 151");
  RLE[3] = ep->parse("F: 1 I: 2 L: 151 N: 1");
  // RLE[4] = ep->parse("A: 35 S: 48 T: 18 Y: 15 C: 9");
  RLE[4] = ep->parse("A: 35 C: 9 E: 2 F: 2 G: 5 H: 1 I: 2 K: 1 L: 3 N: 7 S: 48 T: 18 V: 7 Y: 15");
  // RLE[5] = ep->parse("D: 11 E: 109 Q: 21");
  RLE[5] = ep->parse("C: 2 D: 11 E: 109 G: 2 H: 1 K: 1 Q: 21 S: 4 T: 3 Y: 1");
  RLE[6] = ep->parse("A: 1 C: 4 D: 1 E: 6 F: 5 G: 6 H: 7 I: 5 K: 10 L: 14 M: 3 N: 15 P: 9 Q: 3 R: 5 S: 10 T: 17 V: 18 Y: 16");
#endif

#if VERSION == 2
  RLE[0] = ep->parse("M: 1 P: 1 W: 1 A: 2 D: 2 E: 2 H: 3 Q: 3 V: 4 T: 5 C: 6 Y: 6 G: 7 L: 7 R: 7 F: 13 I: 15 N: 15 K: 26 S: 40");
  RLE[1] = ep->parse("K: 7 R: 159");
  RLE[2] = ep->parse("W: 1 A: 2 E: 2 M: 2 H: 3 R: 3 Q: 4 Y: 4 F: 5 V: 6 C: 7 T: 9 K: 12 L: 16 S: 28 I: 29 N: 33");
  RLE[3] = ep->parse("N: 1 F: 2 I: 2 L: 161");
  RLE[4] = ep->parse("E: 1 H: 1 F: 2 I: 2 K: 3 L: 3 G: 4 N: 8 V: 8 C: 10 Y: 17 T: 18 A: 37 S: 52");
  RLE[5] = ep->parse("H: 1 K: 1 Y: 1 C: 2 G: 3 T: 3 S: 5 D: 15 Q: 21 E: 114");
  RLE[6] = ep->parse("D: 1 A: 3 M: 3 Q: 3 C: 4 F: 5 R: 5 G: 6 I: 6 E: 7 H: 8 P: 9 K: 10 S: 11 L: 15 N: 17 T: 17 Y: 17 V: 19");
#endif

#endif

  KLD[0] = ep->parse("A: 69 L: 1 P: 1 V: 8");
  KLD[1] = ep->parse("K: 58 R: 20");
  KLD[2] = ep->parse("D: 12 E: 10 G: 1 H: 35 N: 19 Y: 2");
  KLD[3] = ep->parse("A: 5 I: 4 L: 29 M: 11 V: 30");
  KLD[4] = ep->parse("F: 17 L: 62");
  KLD[5] = ep->parse("D: 56 E: 23");
  KLD[6] = ep->parse("D: 3 E: 13 G: 3 I: 1 K: 9 L: 1 M: 7 N: 7 Q: 2 R: 20 S: 13");

  a_tail_length = new GHMM::LENGTH::Geometric(364);

#if VERSION == 1
  b_spacer_length = new GHMM::LENGTH::Geometric(350);
  b_tail_length = new GHMM::LENGTH::Geometric(1845);
#else
  b_spacer_length = new GHMM::LENGTH::Geometric(1693);
  b_tail_length = new GHMM::LENGTH::Geometric(437);
#endif

  c_tail_length = new GHMM::LENGTH::Geometric(755);

  GHMM::StateBase::Ptr a_met = GHMM::UTIL::makeState(NULL, met);
  GHMM::StateBase::Ptr a_spacer = GHMM::UTIL::makeState(a_spacer_length, background);
  GHMM::StateBase::Ptr a_RLE = GHMM::UTIL::makeMotifState(RLE);
  GHMM::StateBase::Ptr a_tail = GHMM::UTIL::makeState(a_tail_length, background);

  GHMM::StateBase::Ptr b_met = GHMM::UTIL::makeState(NULL, met);
  GHMM::StateBase::Ptr b_leader = GHMM::UTIL::makeState(b_leader_length, background);
  GHMM::StateBase::Ptr b_KLD = GHMM::UTIL::makeMotifState(KLD);
  GHMM::StateBase::Ptr b_spacer = GHMM::UTIL::makeState(b_spacer_length, background);
  GHMM::StateBase::Ptr b_hydrophobic = GHMM::UTIL::makeState(b_hydrophobic_length, b_hydrophobic_distrib);
  GHMM::StateBase::Ptr b_tail = GHMM::UTIL::makeState(b_tail_length, background);

  GHMM::StateBase::Ptr c_met = GHMM::UTIL::makeState(NULL, met);
  GHMM::StateBase::Ptr c_tail = GHMM::UTIL::makeState(c_tail_length, background);

  GHMM::ModelBuilder mb;

  std::pair<std::string, std::string> makeSignalPModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);
  std::pair<std::string, std::string> makeSSModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);

#ifdef SIGNALP_MODEL
  std::pair<std::string, std::string> ss_states = makeSignalPModel(mb, alphabet);
#else
  std::pair<std::string, std::string> ss_states = makeSSModel(mb, alphabet);
#endif

#ifdef RLE_PATTERN
  mb.addState("a-met",         a_met);
  mb.addState("a-spacer",      a_spacer);
  mb.addState("a-RLE",         a_RLE);
  mb.addState("a-tail",        a_tail);
#endif

#ifdef KLD_PATTERN
  mb.addState("b-met",         b_met);
  mb.addState("b-leader",      b_leader);
  mb.addState("b-KLD",         b_KLD);
  mb.addState("b-spacer",      b_spacer);
  mb.addState("b-hydrophobic", b_hydrophobic);
  mb.addState("b-tail",        b_tail);
#endif

  mb.addState("c-met",         c_met);
  mb.addState("c-tail",        c_tail);

#ifdef RLE_PATTERN
  mb.addStateTransition(GHMM::Model::BEGIN, "a-met", 400);
#endif

#ifdef KLD_PATTERN
  mb.addStateTransition(GHMM::Model::BEGIN, "b-met", 100);
#endif

  mb.addStateTransition(GHMM::Model::BEGIN, "c-met", 4909);

#ifdef RLE_PATTERN
  mb.addStateTransition("a-met",            ss_states.first,  1);
  mb.addStateTransition(ss_states.second,   "a-spacer",       1);
  mb.addStateTransition("a-spacer",         "a-RLE",          1);
  mb.addStateTransition("a-RLE",            "a-tail",         1);
  mb.addStateTransition("a-tail",           GHMM::Model::END, 1);

  mb.addState("d-tail", c_tail);
  mb.addStateTransition(ss_states.second,   "d-tail",         0.01);
  mb.addStateTransition("d-tail",           GHMM::Model::END, 1);
#endif

#ifdef KLD_PATTERN
  mb.addStateTransition("b-met",            "b-leader",       1);
  mb.addStateTransition("b-leader",         "b-KLD",          1);
  mb.addStateTransition("b-KLD",            "b-spacer",       1);
  mb.addStateTransition("b-spacer",         "b-hydrophobic",  1);
  mb.addStateTransition("b-hydrophobic",    "b-tail",         1);
  mb.addStateTransition("b-tail",           GHMM::Model::END, 1);
#endif

  mb.addStateTransition("c-met",            "c-tail",         1);
  mb.addStateTransition("c-tail",           GHMM::Model::END, 1);

  return mb.make();
}

static const struct option options[] = {
  { "input",             required_argument,          0,            'i' },
  { "output",            required_argument,          0,            'o' },
  { "RLE-threshold",     required_argument,          0,            'R' },
  { "KLD-threshold",     required_argument,          0,            'K' },
  { "no-RLD",            no_argument,                0,            'r' },
  { "no-KLD",            no_argument,                0,            'k' },
  { 0,                   0,                          0,            0   }
};

std::string genParse(const std::string &sequence, const GHMM::Model::Ptr &model, GHMM::Traceback::Ptr tbp) {
  std::ostringstream out;
  std::vector<std::string> parse;
  std::string::size_type pos = sequence.size();
  while (tbp != NULL) {
    std::ostringstream o;
    pos -= tbp->length;
    o << "[" << model->stateName(tbp->state) << ":" << sequence.substr(pos, tbp->length) << "]";
    parse.push_back(o.str());
    tbp = tbp->prev;
  }
  std::copy(parse.rbegin(), parse.rend(), std::ostream_iterator<std::string>(out, ""));
  return out.str();
}

void usage(const char *progname) {
  std::cout << "Usage: " << progname << " [arguments]" << std::endl;
  std::cout << "\
\n\
--input=file            -i file        read sequences from file (-:stdin)\n\
--output=file           -o file        write results to file (-:stdout)\n\
--RLE-threshold=float   -R float       RLE threshold for positive prediction\n\
                                       (default: 4.3)\n\
--KLD-threshold=float   -K float       KLD threshold for positive prediction\n\
                                       (default: 0.0)\n\
--no-RLE                -r             turn off RLE prediction\n\
--no-KLD                -k             turn off KLD prediction\n\
\n\
";
}

int main(int argc, char **argv) {
  double RLE_threshold = 4.3;
  double KLD_threshold = 0.0;
  bool do_RLE = true;
  bool do_KLD = false;

  std::list<std::pair<std::string, std::string> > seq_list;
  std::string output = "-";

  int ch;

  while ((ch = getopt_long(argc, argv, "i:o:R:K:hkr", options, NULL)) != -1) {
    switch (ch) {
    case 'i': {
      if (!strcmp(optarg, "-")) {
        std::ifstream in(optarg);
        int old_sz = seq_list.size();
        readFasta(std::cin, std::inserter(seq_list, seq_list.begin()), NULL);
        std::cerr << seq_list.size() - old_sz << " sequences read from stdin" << std::endl;
      } else {
        std::ifstream in(optarg);
        int old_sz = seq_list.size();
        readFasta(in, std::inserter(seq_list, seq_list.begin()), NULL);
        std::cerr << seq_list.size() - old_sz << " sequences read from " << optarg << std::endl;
      }
      break;
    }
    case 'o': {
      output = optarg;
      break;
    }
    case 'R': {
      RLE_threshold = strtod(optarg, NULL);
      break;
    }
    case 'K': {
      KLD_threshold = strtod(optarg, NULL);
      break;
    }
    case 'r': {
      do_RLE = false;
      break;
    }
    case 'k': {
      do_KLD = false;
      break;
    }
    case 'h':
    case '?': {
      usage(argv[0]);
      exit(0);
    }
    }
  }

  GHMM::Model::Ptr model = makePEXELmodel();

  std::list<std::pair<std::string, std::string> >::iterator i, e;

  std::vector<std::pair<double, std::string> > rle_out, kld_out;

  for (i = seq_list.begin(), e = seq_list.end(); i != e; ++i) {
    std::string &name((*i).first);
    std::string &sequence((*i).second);
    int *seq_raw = new int[sequence.size()];

    for (int j = 0; j < (int)sequence.size(); j++) {
      if (isalpha(sequence[j])) {
        seq_raw[j] = toupper(sequence[j]) - 'A';
      } else {
        seq_raw[j] = 'X' - 'A';
      }
    }

    GHMM::Parse::Ptr parse = new GHMM::Parse();
    parse->parse(model, seq_raw, seq_raw + sequence.size());

    double alpha_rle, alpha_kld, alpha_bkg;
    alpha_rle = parse->alpha(model->stateNumber("a-tail"), 0);
    alpha_kld = parse->alpha(model->stateNumber("b-tail"), 0);
    alpha_bkg = parse->alpha(model->stateNumber("c-tail"), 0);
#if 0
    std::cerr << name
              << " alpha_rle:" << alpha_rle
              << " alpha_kld:" << alpha_kld
              << " alpha_bkg:" << alpha_bkg
              << " alpha_ssonly=" << parse->alpha(model->stateNumber("d-tail"), 0) << std::endl;
#endif

    if (alpha_rle - alpha_bkg > RLE_threshold) {
      std::ostringstream out;
      out << name << "\t"
          << "RLE" << "\t"
          << alpha_rle - alpha_bkg << "\t"
          << genParse(sequence, model, parse->psi(model->stateNumber("a-tail"), 0));
      rle_out.push_back(std::make_pair(alpha_rle - alpha_bkg, out.str()));
    }

    if (alpha_kld - alpha_bkg > KLD_threshold) {
      std::ostringstream out;
      out << name << "\t"
          << "KLD" << "\t"
          << alpha_kld - alpha_bkg << "\t"
          << genParse(sequence, model, parse->psi(model->stateNumber("b-tail"), 0));
      kld_out.push_back(std::make_pair(alpha_kld - alpha_bkg, out.str()));
    }
  }

  std::sort(rle_out.begin(), rle_out.end());
  std::sort(kld_out.begin(), kld_out.end());

  if (output == "-") {
    for (std::vector<std::pair<double, std::string> >::reverse_iterator i = rle_out.rbegin(); i != rle_out.rend(); ++i) {
      std::cout << (*i).second << std::endl;
    }
    for (std::vector<std::pair<double, std::string> >::reverse_iterator i = kld_out.rbegin(); i != kld_out.rend(); ++i) {
      std::cout << (*i).second << std::endl;
    }
  } else {
    std::ofstream out;
    out.open(output.c_str());

    for (std::vector<std::pair<double, std::string> >::reverse_iterator i = rle_out.rbegin(); i != rle_out.rend(); ++i) {
      out << (*i).second << std::endl;
    }
    for (std::vector<std::pair<double, std::string> >::reverse_iterator i = kld_out.rbegin(); i != kld_out.rend(); ++i) {
      out << (*i).second << std::endl;
    }
  }
}
