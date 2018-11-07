/*    
*    Copyright 2018 Bartosz Kostrzewa
*    
*    This file is part of sLapH-contractions.
*
*    sLapH-contractions is free software: you can redistribute it and/or modify
*    it under the terms of the GNU General Public License as published by
*    the Free Software Foundation, either version 3 of the License, or
*    (at your option) any later version.
*
*    sLapH-contractions is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*    GNU General Public License for more details.
*
*    You should have received a copy of the GNU General Public License
*    along with sLapH-contractions.  If not, see <http://www.gnu.org/licenses/>.
*/

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

