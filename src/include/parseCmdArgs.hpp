/**
 * @file    parseCmdArgs.hpp
 * @brief   command line parsing
 * @author  Chirag Jain <cjain7@gatech.edu>
 */

#ifndef PAIRG_PARSE_CMD_HPP 
#define PAIRG_PARSE_CMD_HPP

//External includes
#include "clipp/include/clipp.h"

namespace pairg
{
  /**
   * @brief     struct to hold input parameters
   */
  struct Parameters
  {
    std::string graphfile;      //variation graph file
    std::string gmode;          //variation graph input format

    int d_low;                  //lower bound on path length
    int d_up;                   //upper bound on path length
  };

  /**
   * @brief                   parse the cmd line options
   * @param[in]   argc
   * @param[in]   argv
   * @param[out]  param       parameters are saved here
   **/
  void parseandSave(int argc, char** argv, pairg::Parameters &param)
  {
    auto cli = 
      (
       clipp::required("-r") & clipp::value("file", param.graphfile).doc("variation graph file"),
       clipp::required("-m") & 
            (clipp::required("vg").set(param.gmode) | 
            clipp::required("txt").set(param.gmode)).doc("variation graph format"),
       clipp::required("-l") & clipp::value("d1", param.d_low).doc("lower bound on path length"),
       clipp::required("-u") & clipp::value("d2", param.d_up).doc("upper bound on path length")
      );

    if(!clipp::parse(argc, argv, cli)) 
    {
      //print help page
      clipp::operator<<(std::cout, clipp::make_man_page(cli, argv[0])) << std::endl;
      exit(1);
    }
  }
}

#endif
