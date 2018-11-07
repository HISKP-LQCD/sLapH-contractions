#pragma once

#include <omp.h>

#include <iostream>
#include <string>

static inline void lt_print(const std::string msg, const double & local_timer){
  #pragma omp critical(cout)
  {
    std::cout << msg << " Thread " << omp_get_thread_num() << 
    " Timing " << local_timer << " seconds" << std::endl;
  }
}

#define LT_DECLARE \
  double local_timer;

#define LT_START \
  local_timer = omp_get_wtime();

#define LT_STOP \
	local_timer = omp_get_wtime() - local_timer;

#define LT_PRINT(msg) \
  _Pragma("omp critical(cout)") \
  { \
    std::cout << msg << " Thread " << omp_get_thread_num() << \
    " Timing " << local_timer << " seconds" << std::endl; \
  }

//lt_print((msg), local_timer);

#ifdef SLAPH_TIMING_DIAGRAMS
#define LT_DIAGRAMS_DECLARE    LT_DECLARE
#define LT_DIAGRAMS_START      LT_START
#define LT_DIAGRAMS_STOP       LT_STOP
#define LT_DIAGRAMS_PRINT(msg) LT_PRINT(msg)
#else
#define LT_DIAGRAMS_DECLARE   
#define LT_DIAGRAMS_START     
#define LT_DIAGRAMS_STOP      
#define LT_DIAGRAMS_PRINT(msg)
#endif

#ifdef SLAPH_TIMING_CORRELATOR
#define LT_CORRELATOR_DECLARE    LT_DECLARE
#define LT_CORRELATOR_START      LT_START
#define LT_CORRELATOR_STOP       LT_STOP
#define LT_CORRELATOR_PRINT(msg) LT_PRINT(msg)
#else
#define LT_CORRELATOR_DECLARE   
#define LT_CORRELATOR_START     
#define LT_CORRELATOR_STOP      
#define LT_CORRELATOR_PRINT(msg)
#endif


#ifdef SLAPH_TIMING_FINE
#define LT_FINE_DECLARE    LT_DECLARE
#define LT_FINE_START      LT_START
#define LT_FINE_STOP       LT_STOP
#define LT_FINE_PRINT(msg) LT_PRINT(msg)
#else
#define LT_FINE_DECLARE   
#define LT_FINE_START     
#define LT_FINE_STOP      
#define LT_FINE_PRINT(msg)
#endif

#ifdef SLAPH_TIMING_ULTRA_FINE
#define LT_ULTRA_FINE_DECLARE    LT_DECLARE
#define LT_ULTRA_FINE_START      LT_START
#define LT_ULTRA_FINE_STOP       LT_STOP
#define LT_ULTRA_FINE_PRINT(msg) LT_PRINT(msg)
#else
#define LT_ULTRA_FINE_DECLARE   
#define LT_ULTRA_FINE_START     
#define LT_ULTRA_FINE_STOP      
#define LT_ULTRA_FINE_PRINT(msg)
#endif

