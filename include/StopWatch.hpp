#pragma once

#include <vector>

/**
   Usage example:

       StopWatch sw;

       #pragma omp parallel
       {
           sw.start();

           // Heavy work

           sw.stop();
       }

       double wtime = sw.mean();
*/
class StopWatch {
 public:
  StopWatch(char const *const message, int const threads = -1);

  void start();
  void stop();

  double mean() const;

  void print() const;

 private:
  void store_to(std::vector<double> &vec);

  int summands() const;

  std::vector<double> starts_;
  std::vector<double> ends_;
  std::vector<double> used_;
  char const *const message_;
};
