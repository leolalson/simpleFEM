#pragma once

#include<vector>

#include "io.hpp"
#include "fem.hpp"
#include<eigen3/Eigen/Core>

class domain{
  public:
  domain(){};
  domain(element* elem, int dim);
  ~domain(){};
  int dim;
  Eigen::MatrixXi data;
  Eigen::VectorXi tags;
  size_t size = 0;
  element* elem;

  void getDomainData(iovtk ioobj);
  std::map<size_t, size_t> getDomainTags(iovtk ioobj);

  domain getSubdomain(std::vector<size_t> subDomainTags);
};


class points{
  public:
  size_t size = 0;
  Eigen::MatrixXd data;
  points(){};
  ~points(){};
  points(Eigen::MatrixXd pointsData) : data(pointsData){ this->size = pointsData.rows();};
};

class mesh{
  public:
  mesh(){};
  mesh(element* elem);
  ~mesh(){};
  domain topology;
  points nodes;
  element* elem;

  void getMeshData(iovtk vtkobj); 
};


