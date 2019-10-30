#pragma once

#include <omp.h>
#include <boost/algorithm/string/replace.hpp>

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

constexpr int timing_level = SLAPH_TIMING_LEVEL;

struct TimingNode;
struct TimingEdge;

struct TimingEdge {
  TimingEdge(TimingNode *const source, TimingNode *const destination)
      : source(source), destination(destination), start(omp_get_wtime()) {}

  TimingNode *source;
  TimingNode *destination;

  double start;

  double cumtime = 0;
  int calls = 0;
};

struct TimingNode {
  TimingNode(std::string const &function, std::string const &info = "")
      : start(omp_get_wtime()), function(function), info(info) {}

  std::vector<TimingEdge *> edges;

  double start;

  double cumtime = 0.0;
  double selftime = 0.0;
  int calls = 0;
  std::string function;
  std::string info;
};

class TimingGraph {
 public:
  static TimingGraph &instance() {
    static TimingGraph instance;
    return instance;
  }

  void push(std::string const &function, std::string const &info = "");
  void pop();
  void serialize(std::ostream &ofs);
  void finalize();

 private:
  TimingGraph(){};

  TimingGraph(TimingGraph const &) = delete;

  std::vector<TimingEdge> edges_;
  std::vector<TimingNode> nodes_;

  std::vector<TimingEdge *> edge_stack_;
  std::vector<TimingNode *> node_stack_;
};

template <int level>
class TimingScope {
 public:
  TimingScope(std::string const &function, std::string const &info = "") {
#pragma omp master
    if (level <= timing_level) {
      TimingGraph::instance().push(function, info);
    }
  }

  ~TimingScope() {
#pragma omp master
    if (level <= timing_level) {
      TimingGraph::instance().pop();
    }
  }
};
