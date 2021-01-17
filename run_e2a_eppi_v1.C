#ifndef rune2aeppi
#define rune2aeppi
#include "e2a_eppi_v1.h"
// apapadop
#include "FilterData.h"
#include "FilterGiBUU.h"

#include <iostream>

using namespace std;


int main(int argc, char **argv)
{
  
  //this could be adjusted to take the path to data as an argument
  if( argc < 5 ){
    //cout<<"Please specify the absolute path to data, the target (3He, 56Fe, C12, 4He), the beam energy (1161, 2261 or 4461) and what you want to do (Create Plots = 0, Filter = 1)"<<endl;
    cout<<"Please specify the absolute path to data, the target (3He, 56Fe, C12, 4He), the beam energy (1161, 2261 or 4461) and what you want to do (Create Plots = 0, FilterData = 1, Filter GiBUU = 2)"<<endl;
    cout << "Path to data should be to a folder containing subfolders of the form \'e2a_<target>_<beam_energy>_v1\'" << endl;
    cout<<"================= Usage ==============="<<endl;
    //cout<<"./run_e2a_eppi_v1.cc <path_to_data> <target> <beam_energy> {0/1} (Create Plots = 0, Filter = 1)"<<endl;
    cout<<"./run_e2a_eppi_v1.cc <path_to_data> <target> <beam_energy> {0/1} (Create Plots = 0, Filter Data = 1, Filter GiBUU = 2)"<<endl;
    exit(1);
  }
  
  /*
  //this could be adjusted to take the path to data as an argument
  if( argc < 4 ){
    cout<<"Please specify the target (3He, 56Fe, C12, 4He), the beam energy (1161, 2261 or 4461) and what you want to do (Create Plots = 0, Filter = 1)"<<endl;
    cout<<"================= Usage ==============="<<endl;
    cout<<"./run_e2a_eppi_v1.cc <target> <beam_energy> {0/1} (Create Plots = 0, Filter = 1)"<<endl;
    exit(1);
  }
  
  std::string target	= argv[1];
  std::string beam_en = argv[2];
  int choice = atoi(argv[3]);
  */

  std::string path	= argv[1];
  std::string target	= argv[2];
  std::string beam_en = argv[3];
  int choice = atoi(argv[4]);
  
  //if (choice != 1 && choice != 0) {
  //  std::cout << "Unknown option for parameter 3. It should be either 0 or 1. The given value is " << choice << std::endl;
  //  return 0;
  //}

   if (choice != 2 && choice != 1 && choice != 0) {
     std::cout << "Unknown option for parameter 3. It should be 0, 1, or 2. The given value is " << choice << std::endl;
     return 0;
   }

  //creates e2a_eppi object and filter data object  
  e2a_eppi_v1	t(path,target,beam_en);      //This runs the analysis in ROOT, producing a bunch of plots
  FilterData	filterD(path,target,beam_en); //This performs basic PID then writes data to a genie compatible tree
  FilterGiBUU	filterG(path,target,beam_en); //This converts a GiBUU tree to a genie compatible tree


  if (choice == 0) {
    cout << "Creating Plots" << endl;
    t.Loop();
    cout << "Event Loop Complete. Post Loop Plot Analysis" << endl;
    t.Analyse();
  }
  
  if (choice == 1) {
    cout << "Filtering sample to genie compatabile tree" << endl;
    filterD.Loop();
  }

  if (choice == 2) {
    cout << "Filtering GiBUU sample to genie compatabile tree" << endl;
    cout << "This code is under construction, and currently does nothing. Exiting." << endl;
  //  filterG.Loop();
  }
  
  
  std::cout << "Loop complete" << std::endl;
  
  
  return 0;
}
#endif

