#include<iostream>
#include<fstream>
#include<vector>
#include <sstream>
#include <eigen3/Eigen/Core>
#include <map>

#include "core/utils.hpp"
#include "core/io.hpp"

class mesh{
  public:
  mesh();
  mesh(size_t dim, size_t elemType);
  ~mesh();
  std::vector<std::vector<uint32_t>> elements;
  std::vector<std::vector<double>> nodes;
  uint32_t numNodes;
  uint32_t numElements;
  size_t dim;
  size_t elemType; 
};

mesh::mesh(){};

mesh::~mesh(){};

mesh::mesh(size_t dim, size_t elemType) : dim(dim), elemType(elemType){};

int main(){

  std::string meshFileName = "rectangle_4NodeQuad.vtk";
  int dim = 1; //2D case
  int elemType = 0; //triangle
  std::vector<int> a = {1,2};

  utils::printVector(a, "a");
  std::cout << "Lalson: " << std::endl;
  iovtk vtkobj = iovtk::readvtk(meshFileName);
  utils::printVector2D(vtkobj.nodes, "nodes");
  utils::printVector2D(vtkobj.elements, "elements");
  utils::printVector(vtkobj.data[0], " 1");
  utils::printVector(vtkobj.data[1], " 2");
  return 0;
}