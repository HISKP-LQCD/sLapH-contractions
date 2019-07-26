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
