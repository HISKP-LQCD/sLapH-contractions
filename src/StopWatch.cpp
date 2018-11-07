#include <StopWatch.h>

#include "typedefs.h"

#include <omp.h>

#include <algorithm>
#include <cassert>
#include <iomanip>
#include <iostream>

StopWatch::StopWatch(char const *const message, int const threads)
    : starts_(threads == -1 ? omp_get_max_threads() : threads, 0.0),
      ends_(starts_.size(), 0.0),
      used_(starts_.size(), false),
      message_(message) {
  std::cout << "[Start    ] " << message << std::endl;
}

void StopWatch::start() {
  store_to(starts_);
}

void StopWatch::stop() {
  assert(summands() > 0 && "This stopwatch has not been started yet, cannot stop.");
  store_to(ends_);
}

double StopWatch::mean() const {
#pragma omp barrier
  double sum = 0.0;
  for (int i = 0; i < ssize(starts_); ++i) {
    sum += ends_[i] - starts_[i];
  }
  return sum / summands();
}

void StopWatch::print() const {
  std::cout << "[   Finish] " << message_ << ". Duration: " << std::setprecision(10)
            << mean() << " seconds. Threads: " << summands() << std::endl;
}

void StopWatch::store_to(std::vector<double> &vec) {
  int const t_id = omp_get_thread_num();
  assert(t_id < ssize(vec) &&
         "This particular stopwatch does not support this many threads.");
  double const wtime = omp_get_wtime();
  vec.at(t_id) = wtime;
  used_[t_id] = true;
}

int StopWatch::summands() const {
  return std::count(std::begin(used_), std::end(used_), true);
}
