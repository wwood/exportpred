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
#include <GHMM/ghmm.hh>

const std::string GHMM::UTIL::Alphabet::NO_TOKEN;

#define IFTYPE(klass, p) if (klass *ptr = dynamic_cast<klass *>(p.ptr()))

template<typename T, typename U>
static GHMM::StateBase::Ptr makeState_2(const T &length,
                                 const U &emission) {
  return new GHMM::State<T, U>(length, emission);
}

template<typename T>
static GHMM::StateBase::Ptr makeState_1(const T &length,
                                 GHMM::EMISSION::Base::Ptr emission) {
  IFTYPE(GHMM::EMISSION::Stateless, emission)
    return makeState_2(length, (GHMM::EMISSION::Stateless &)*ptr);
  IFTYPE(GHMM::EMISSION::PositionSpecific, emission)
    return makeState_2(length, (GHMM::EMISSION::PositionSpecific &)*ptr);
  return NULL;
}

GHMM::StateBase::Ptr GHMM::UTIL::makeState(GHMM::LENGTH::Base::Ptr length, GHMM::EMISSION::Base::Ptr emission) {
  if (length == NULL) {
    length = new GHMM::LENGTH::Fixed(1);
  }
  IFTYPE(GHMM::LENGTH::Uniform, length) {
    length = new GHMM::LENGTH::Discrete(ptr->discretise());
  }
  IFTYPE(GHMM::LENGTH::Discrete, length)
    return makeState_1((GHMM::LENGTH::Discrete &)*ptr, emission);
  IFTYPE(GHMM::LENGTH::Fixed, length)
    return makeState_1((GHMM::LENGTH::Fixed &)*ptr, emission);
  IFTYPE(GHMM::LENGTH::Geometric, length)
    return makeState_1((GHMM::LENGTH::Geometric &)*ptr, emission);
  return NULL;
}

GHMM::StateBase::Ptr GHMM::UTIL::makeMotifState(const std::vector<MATH::DPDF::Ptr> &pssm) {
  return new GHMM::State<GHMM::LENGTH::Fixed, GHMM::EMISSION::PositionSpecific>(GHMM::LENGTH::Fixed(pssm.size()), GHMM::EMISSION::PositionSpecific(pssm));
}
