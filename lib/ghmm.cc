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

const static int C_BEGIN = -1;
const static int C_END = -2;

using namespace GHMM;

const std::string Model::BEGIN = "__BEGIN__";
const std::string Model::END = "__END__";

Model::Model(const std::vector<std::pair<std::string, StateBase::Ptr> > &in_states,
             const std::map<std::pair<int, int>, double> &in_state_trans_map) :
  RefObj(), state_names(), state_name_map(), pred_states(), succ_states(), states(), state_trans(NULL), state_count(0) {

  std::vector<bool> reachable(in_states.size(), false);
  {
    int l = in_states.size() + 2;
    int *a = new int[l * l];

    std::fill(a, a + l * l, 0);
    for (std::map<std::pair<int, int>, double>::const_iterator i = in_state_trans_map.begin(), e = in_state_trans_map.end(); i != e; ++i) {
      int s = (*i).first.first;
      int t = (*i).first.second;
      if (s < 0) s += l;
      if (t < 0) t += l;
      a[s * l + t] = 1;
    }
    for (int k = 0; k < l; k++) {
      for (int i = 0; i < l; i++) {
        for (int j = 0; j < l; j++) {
          a[i * l + j] = a[i * l + j] | (a[i * l + k] & a[k * l + j]);
        }
      }
    }
    for (int i = 0; i < (int)in_states.size(); i++) {
      if (a[(l + C_BEGIN) * l + i] && a[(i * l) + (l + C_END)]) {
        reachable[i] = true;
      }
    }
    delete [] a;
  }

  std::map<int, int> state_num_remap;

  states.push_back(StateBase::Ptr());
  state_names.push_back(BEGIN);
  state_num_remap[C_BEGIN] = 0;

  for (unsigned i = 0; i < in_states.size(); i++) {
    if (reachable[i] && in_states[i].second != NULL) {
      // std::cerr << "REMAP " << i << " : " << in_states[i].first << " -> " << states.size() << std::endl;
      state_num_remap[i] = states.size();
      state_name_map[in_states[i].first] = states.size();
      state_names.push_back(in_states[i].first);
      states.push_back(in_states[i].second);
    } else {
      state_num_remap[i] = -1;
    }
  }

  state_num_remap[C_END] = states.size();
  states.push_back(StateBase::Ptr());
  state_names.push_back(END);

  state_count = states.size();

  pred_states.resize(state_count);
  succ_states.resize(state_count);

  state_trans = new double[state_count * state_count];
  state_log_trans = new double[state_count * state_count];

  std::fill(state_trans, state_trans + state_count * state_count, 0.0);

  for (std::map<std::pair<int, int>, double>::const_iterator i = in_state_trans_map.begin(), e = in_state_trans_map.end(); i != e; ++i) {
    int s = state_num_remap[(*i).first.first];
    int t = state_num_remap[(*i).first.second];
    double p = (*i).second;
    // std::cerr << "s=" << s << " -> t=" << t << " : " << p << std::endl;
    state_trans[s * state_count + t] = p;
  }

  for (int s = 0; s < state_count - 1; s++) {
    int idx = s * state_count;
    double sum = std::accumulate(state_trans + idx, state_trans + idx + state_count, 0.0);
    // std::cerr << "s=" << s << " sum=" << sum << std::endl;
    assert(sum > 0.0);

    for (int t = 0; t < state_count; t++) {
      state_trans[idx + t] /= sum;
    }


    if (LENGTH::Geometric *gp = dynamic_cast<LENGTH::Geometric *>(&*states[s])) {
      double p_self = gp->pLength(1);
      double not_p_self = 1.0 - p_self;
      for (int t = 0; t < s; t++) {
        state_trans[s * state_count + t] *= not_p_self;
      }
      state_trans[s * state_count + s] = p_self;
      for (int t = s + 1; t < state_count; t++) {
        state_trans[s * state_count + t] *= not_p_self;
      }
    }

    for (int t = 0; t < state_count; t++) {
      if (state_trans[idx + t] != 0.0) {
        pred_states[t].push_back(s);
        succ_states[s].push_back(t);
      }
    }
  }

  for (int i = 0; i < state_count * state_count; i++) {
    state_log_trans[i] = MATH::logClip(state_trans[i]);
  }

//   for (int i = 0; i < state_count; i++) {
//     for (int j = 0; j < state_count; j++) {
//       fprintf(stderr, "%9.7f ", state_trans[i * state_count + j]);
//     }
//     fprintf(stderr, "\n");
//   }
}

Model::~Model() {
  if (state_trans) delete [] state_trans;
  if (state_log_trans) delete [] state_log_trans;
}

ModelBuilder::ModelBuilder() :
  RefObj(), states(), state_name_map(), state_trans_map() {
  state_name_map[Model::BEGIN] = C_BEGIN;
  state_name_map[Model::END] = C_END;
}

ModelBuilder::~ModelBuilder() {
}

void ModelBuilder::addState(const std::string &state_name, const StateBase::Ptr &state) {
  int i = stateNum(state_name);
  if (i < 0) return;
  states[i] = std::make_pair(state_name, state);
}

void ModelBuilder::removeState(const std::string &state_name) {
  int s = stateNum(state_name);
  if (s < 0) return;

  // remove state.
  states[s].second = NULL;

  // remove transitions.
  std::map<std::pair<int, int>, double>::iterator i = state_trans_map.begin(), j, e = state_trans_map.end();
  while (i != e) {
    if ((*i).first.first == s or (*i).first.second == s) {
      j = i;
      ++j;
      state_trans_map.erase(i);
    } else {
      ++i;
    }
  }
}

int ModelBuilder::stateNum(const std::string &s) {
  std::map<std::string, int>::iterator i = state_name_map.find(s);
  if (i != state_name_map.end()) return i->second;
  state_name_map[s] = states.size();
  states.push_back(std::make_pair(s, StateBase::Ptr()));
  return states.size() - 1;
}

void ModelBuilder::addStateTransition(const std::string &src_state, const std::string &tgt_state, double freq) {
  state_trans_map[std::make_pair(stateNum(src_state), stateNum(tgt_state))] = freq;
}

Model::Ptr ModelBuilder::make() {
  return new Model(states, state_trans_map);
}
