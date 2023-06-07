#pragma once

#include<vector>
#include <sstream>
#include <map>
#include <eigen3/Eigen/Core>
#include <iostream>
#include <eigen3/Eigen/Dense>

  

class element{
  public:
  element();
  
  ~element();
  uint32_t numNodes;
  size_t dim;
  //size_t elemType;
  std::map<size_t, size_t> elemType;
  double thick;
  int dof;
  Eigen::MatrixXd gaussPts2D;
  Eigen::VectorXd gaussWts2D;
  Eigen::MatrixXd gaussPts1D;
  Eigen::VectorXd gaussWts1D;
  Eigen::MatrixXf K;

  static std::map<std::string, element * (*) ()> elementMap;

  virtual Eigen::MatrixXd Kmatrix(Eigen::MatrixXd coords, Eigen::MatrixXd D){
    Eigen::MatrixXd K_local(dof, dof);
    K_local.setZero();
    return K_local;
  };

  virtual Eigen::VectorXd fvector(Eigen::MatrixXd Xe, Eigen::VectorXd F){
    Eigen::VectorXd f_local(4);
    f_local.setZero();
    return f_local;
  }

  size_t get_elemType(size_t dim){
    return elemType[dim];
  }
};

class tri3 : public element{
  public:
  tri3();
  ~tri3();
  static element * create() { return new tri3(); }

  Eigen::MatrixXd Kmatrix(Eigen::MatrixXd coords, Eigen::MatrixXd D);
  Eigen::VectorXd fvector(Eigen::MatrixXd Xe, Eigen::VectorXd F);

};



