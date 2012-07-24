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

std::pair<std::string, std::string> makeSignalPModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet) {
  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);

  GHMM::EMISSION::Base::Ptr a_emit, h1_emit, cloop_emit, c6_emit, c5_emit, c4_emit, c3_emit, c2_emit, c1_emit, cut_emit;

  a_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
        A:0.0669576 C:0.0289179 D:0.0330147 E:0.06290 68:0.629 F:0.040006:0.353\n\
        G:0.0611961 H:0.022283  I:0.0269273 K:0.0721042 L:0.101807\n\
        M:0.0147688 N:0.0299838 P:0.0741557 Q:0.0418068 R:0.119932\n\
        S:0.0804025 T:0.0474921 V:0.0344917 W:0.0186883 Y:0.0221576"));

  h1_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.118119 C:0.03861 D:0.000595187 E:0.00177443 F:0.0767409\n\
    G:0.0306102 H:0.00313159 I:0.0665444 K:0.00070983 L:0.395109\n\
    M:0.0156903 N:0.00237627 P:0.0130463 Q:0.00613188 R:0.00187656\n\
    S:0.0526021 T:0.0373178 V:0.113108 W:0.0174949 Y:0.00841131"));

  cloop_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.147336 C:0.0414663 D:0.0179464 E:0.0251696 F:0.0542468\n\
    G:0.0665154 H:0.0304854 I:0.0162634 K:0.00849114 L:0.193351\n\
    M:0.0415414 N:0.0210277 P:0.0387329 Q:0.0384388 R:0.0213533\n\
    S:0.118563 T:0.0556446 V:0.0349755 W:0.0172798 Y:0.011172"));

  c6_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.153224 C:0.0246704 D:0.0116804 E:0.0201758 F:0.0563868\n\
    G:0.04778 H:0.0171105 I:0.0672533 K:0.00937792 L:0.112884\n\
    M:0.0147248 N:0.0121292 P:0.162384 Q:0.0300841 R:0.0137437\n\
    S:0.0833226 T:0.0327591 V:0.0940115 W:0.0265675 Y:0.00973044"));

  c5_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.142938 C:0.0248836 D:0.0254452 E:0.0290891 F:0.0117752\n\
    G:0.119398 H:0.0334432 I:0.0136058 K:0.0243885 L:0.0643162\n\
    M:0.00506387 N:0.0335328 P:0.101393 Q:0.0582125 R:0.0347008\n\
    S:0.124018 T:0.0840858 V:0.0478986 W:0.0151097 Y:0.006702"));

  c4_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.0927543 C:0.0364718 D:0.0134743 E:0.0330211 F:0.0257115\n\
    G:0.137211 H:0.0146256 I:0.0396322 K:0.0200229 L:0.12227\n\
    M:0.0124529 N:0.01005 P:0.0869224 Q:0.0445039 R:0.0373149\n\
    S:0.0985152 T:0.070977 V:0.0759606 W:0.0129497 Y:0.015159"));

  c3_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.270401 C:0.079943 D:0.00424344 E:0.00106659 F:0.00318821\n\
    G:0.0769233 H:0.00214294 I:0.0343464 K:0.00429208 L:0.0592223\n\
    M:0.00526526 N:0.00427151 P:0.0031593 Q:0.00316105 R:0.00850882\n\
    S:0.129776 T:0.108417 V:0.199535 W:0.00213739"));

  c2_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.0774475 C:0.0232888 D:0.0383256 E:0.0701399 F:0.0318771\n\
    G:0.0365087 H:0.0502663 I:0.0191553 K:0.0137405 L:0.170265\n\
    M:0.0191325 N:0.0353537 P:0.0116314 Q:0.0801956 R:0.0561667\n\
    S:0.11924 T:0.0505315 V:0.0359988 W:0.030043 Y:0.0306924"));

  c1_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.508548 C:0.0520235 G:0.191738 L:0.0236176 P:0.0297181\n\
    Q:0.0161297 S:0.132299 T:0.0459261"));

  cut_emit = new GHMM::EMISSION::Stateless(ep->parse("\n\
    A:0.14242   C:0.0202267 D:0.0728223 E:0.0848579 F:0.0339855\n\
    G:0.0545347 H:0.0234416 I:0.0306897 K:0.0472419 L:0.0771698\n\
    M:0.0128499 N:0.0274769 P:0.00634449 Q:0.111559 R:0.0426487\n\
    S:0.0832117 T:0.0493584 V:0.0492212 W:0.00964942 Y:0.0202909"));

  MATH::DPDF::Ptr h1_dpdf = new MATH::DPDF();
  h1_dpdf->setDistrib(6, 19, (double[]){
                        6.54977e-16,
                          3.7347e-12,
                          1.71429e-05,
                          0.00274341,
                          0.000951411,
                          0.00769368,
                          0.231647,
                          0.220435,
                          0.155337,
                          0.157697,
                          0.222918,
                          0.000560533,
                          1.92736e-08
                          });

  GHMM::LENGTH::Discrete::Ptr h1_length = new GHMM::LENGTH::Discrete(h1_dpdf);

  GHMM::StateBase::Ptr sp_a1 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(3), a_emit);
  GHMM::StateBase::Ptr sp_a2 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Geometric(15.61295), a_emit);
  GHMM::StateBase::Ptr sp_h1 =    GHMM::UTIL::makeState(h1_length, h1_emit);
  GHMM::StateBase::Ptr sp_cloop = GHMM::UTIL::makeState(new GHMM::LENGTH::Geometric(4.155), cloop_emit);
  GHMM::StateBase::Ptr sp_c9 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), cloop_emit);
  GHMM::StateBase::Ptr sp_c8 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), cloop_emit);
  GHMM::StateBase::Ptr sp_c7 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), cloop_emit);
  GHMM::StateBase::Ptr sp_c6 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), c6_emit);
  GHMM::StateBase::Ptr sp_c5 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), c5_emit);
  GHMM::StateBase::Ptr sp_c4 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), c4_emit);
  GHMM::StateBase::Ptr sp_c3 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), c3_emit);
  GHMM::StateBase::Ptr sp_c2 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), c2_emit);
  GHMM::StateBase::Ptr sp_c1 =    GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), c1_emit);
  GHMM::StateBase::Ptr sp_cut =   GHMM::UTIL::makeState(new GHMM::LENGTH::Fixed(1), cut_emit);

  mb.addState("a1",            sp_a1);
  mb.addState("a2",            sp_a2);
  mb.addState("h1",            sp_h1);
  mb.addState("cloop",         sp_cloop);
  mb.addState("c9",            sp_c9);
  mb.addState("c8",            sp_c8);
  mb.addState("c7",            sp_c7);
  mb.addState("c6",            sp_c6);
  mb.addState("c5",            sp_c5);
  mb.addState("c4",            sp_c4);
  mb.addState("c3",            sp_c3);
  mb.addState("c2",            sp_c2);
  mb.addState("c1",            sp_c1);
  mb.addState("cut",           sp_cut);

  mb.addStateTransition("a1",            "a2",           0.0217286);
  mb.addStateTransition("a2",            "h1",           0.0217286);
  mb.addStateTransition("h1",            "cloop",        0.0217286);
  mb.addStateTransition("h1",            "c9",           0.057086);
  mb.addStateTransition("h1",            "c8",           0.151762);
  mb.addStateTransition("h1",            "c7",           0.0578779);
  mb.addStateTransition("h1",            "c6",           0.0568417);
  mb.addStateTransition("h1",            "c5",           0.573097);
  mb.addStateTransition("h1",            "c4",           0.0416052);
  mb.addStateTransition("h1",            "c3",           0.0400018);

  mb.addStateTransition("cloop",         "c9",           1);
  mb.addStateTransition("c9",            "c8",           1);
  mb.addStateTransition("c8",            "c7",           1);
  mb.addStateTransition("c7",            "c6",           1);
  mb.addStateTransition("c6",            "c5",           1);
  mb.addStateTransition("c5",            "c4",           1);
  mb.addStateTransition("c4",            "c3",           1);
  mb.addStateTransition("c3",            "c2",           1);
  mb.addStateTransition("c2",            "c1",           1);
  mb.addStateTransition("c1",            "cut",          1);

  // need to add non-emitting states.
  return std::make_pair(std::string("a1"), std::string("cut"));
}
