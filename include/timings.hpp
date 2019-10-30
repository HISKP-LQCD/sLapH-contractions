#pragma once

#include <omp.h>
#include <boost/algorithm/string/replace.hpp>

#include <fstream>
#include <iomanip>
#include <limits>
#include <sstream>
#include <vector>

constexpr int timing_level = SLAPH_TIMING_LEVEL;

class TimingsStreamSingleton {
 public:
  TimingsStreamSingleton() {
    if (timing_level > 0) {
      for (int i = 0; i < omp_get_max_threads(); ++i) {
        std::ostringstream oss;
        oss << "timings-thread-" << i << ".xml";
        streams_.emplace_back(oss.str());
      }

      for (auto &stream : streams_) {
        stream << std::setprecision(std::numeric_limits<double>::max_digits10);
        stream << "<timings>\n";
      }
    }
  }

  ~TimingsStreamSingleton() {
    if (timing_level > 0) {
      for (auto &stream : streams_) {
        stream << "</timings>\n";
      }
    }
  }

  static TimingsStreamSingleton &instance() {
    static TimingsStreamSingleton instance;
    return instance;
  }

  std::ofstream *get() {
    if (timing_level > 0) {
      int const tid = omp_get_thread_num();
      return &streams_.at(tid);
    } else {
      return nullptr;
    }
  }

 private:
  std::vector<std::ofstream> streams_;
};

template <int level>
class TimingScope {
 public:
  TimingScope(std::string const &function, std::string const &info = "") {
    if (level <= timing_level) {
      start_ = omp_get_wtime();
      stream_ = TimingsStreamSingleton::instance().get();
      std::string f = function;
      f = boost::replace_all_copy(f, "<", "&lt;");
      f = boost::replace_all_copy(f, ">", "&gt;");
      (*stream_) << "<call function='" << f << "' info='" << info << "'>\n";
    }
  }

  ~TimingScope() {
    if (level <= timing_level) {
      auto const end = omp_get_wtime();
      auto const duration = end - start_;
      (*stream_) << "<total time='" << duration << "' />\n</call>\n";
    }
  }

 private:
  double start_;
  std::ofstream *stream_ = nullptr;
};

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
  TimingGraph() {
    // TimingNode root("root");
    // nodes_.push_back(root);
  }

  void push(std::string const &function, std::string const &info = "");
  void pop();
  void serialize(std::ostream &ofs);


 private:
  std::vector<TimingEdge> edges_;
  std::vector<TimingNode> nodes_;

  std::vector<TimingEdge *> edge_stack_;
  std::vector<TimingNode *> node_stack_;
};
