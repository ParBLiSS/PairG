/**
 * @file    utils.hpp
 * @brief   functions for common use cases
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PSGL_UTILS_HPP
#define PSGL_UTILS_HPP

#include <random>
#include <iterator>
#include <algorithm>
#include <thread>
#include <chrono>
#include <fstream>
#include <omp.h>
#include <immintrin.h>
#include <cassert>
#include <stddef.h>

namespace psgl
{
  namespace random
  {

    /**
     * @brief       select a random element from C++ container 
     *              given a random generator
     */
    template<typename Iter, typename RandomGenerator>
      Iter select(Iter start, Iter end, RandomGenerator& g) 
      {
        std::uniform_int_distribution<> dis(0, std::distance(start, end) - 1);
        std::advance(start, dis(g));
        return start;
      }

    /**
     * @brief       select a random element from C++ container
     * @details     borrowed from stackoverflow answer by Christopher Smith
     *              https://stackoverflow.com/questions/6942273
     */
    template<typename Iter>
      Iter select(Iter start, Iter end) 
      {
        //use a fixed seed to get deterministic graph topology
        static std::mt19937 gen(41);
        return select(start, end, gen);
      }
  }

  /**
   * @brief     check if file is accessible
   */
  bool fileExists(const std::string &filename)
  {
    std::ifstream infile(filename);
    return infile.good();
  }
}

#endif
