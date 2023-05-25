#include<iostream>
#include<fstream>
#include<vector>
#include <sstream>
#include <map>
#include <limits>

class iovtk{
  public:
  std::string fileName;
  int dataItems;

  std::vector<std::vector<double>> nodes;
  std::vector<std::vector<size_t>> elements;
  std::vector<std::vector<size_t>> data;
  std::vector<std::string> scalarDataID;
  std::map<std::string, int> skiprows{{"CELL_TYPES", 0}, {"CELL_DATA", 2}};
  static iovtk readvtk(std::string fileName);


  private:
  std::vector<std::vector<double>> read_nodes(std::string fileName);
  std::vector<std::vector<size_t>> read_elements(std::string fileName);
  std::vector<size_t> readScalarData(std::string fileName, std::string dataName);
};