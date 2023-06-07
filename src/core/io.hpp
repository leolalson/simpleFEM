#pragma once
#include<iostream>
#include<fstream>
#include<vector>
#include <sstream>
#include <map>
#include <limits>
#include <iomanip>

#include "fem.hpp"

class iovtk{
  public:
  size_t numNodes;
  size_t numElements;

  iovtk(){};
  ~iovtk(){};
  iovtk(std::vector<std::vector<double>> nodeData, std::vector<std::vector<size_t>> elementData, element* elem);
  std::vector<std::vector<double>> nodes;
  std::vector<std::vector<size_t>> elements;
  std::map<std::string, std::vector<size_t>> data;
  std::vector<std::string> dataID;
  std::map<std::string, int> skiprows{{"CELL_TYPES", 0}, {"CELL_DATA", 2}};
  std::map<int, int> cellTypeMap{
                                  {1, 1},
                                  {2, 2},
                                  {3, 3},
                                  {4, 4},
                                  {5, 5}, //triangle 
                                  {9, 9}  //quad
                                };
  static iovtk readvtk(std::string fileName);
  std::vector<std::vector<size_t>> get_topology(int cellType);
  std::vector<size_t> get_data(std::string cellID, int cellType);
  void writeVtk(std::string fileName, const double solution[]);


  private:
  std::vector<std::vector<double>> read_nodes(std::string fileName);
  std::vector<std::vector<size_t>> read_elements(std::string fileName);
  std::vector<size_t> readScalarData(std::string fileName, std::string dataName);

};